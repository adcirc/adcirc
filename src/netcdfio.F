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
      !=================================================================
      MODULE NETCDFIO

      USE GLOBAL
      USE SIZES
      USE HARM
      USE NodalAttributes
      USE WIND
      USE VERSION  ! uncomment for version 48
!--------------------------------------------
!     Force explicit declaration of all variables
!--------------------------------------------
      IMPLICIT NONE
!     -------------------
!     netCDF header file
!     -------------------
!     include '/ifs1/home/cforbes/usr/local/netcdf/include/netcdf.inc'  ! Topsail
!     include '/usr/local/netcdf/include/netcdf.inc'                    ! Mac
      include 'netcdf.inc'                                              ! Jason
!     -------------
!     Dimension ids
!     -------------
      INTEGER  timenc_dim
!     -----------------
!     Dimension lengths
!     -----------------
      INTEGER timenc_len
!     ------------
!     Variable ids
!     ------------
      INTEGER  timenc_id
!     ---------------------------------------------
!     Rank (number of dimensions) for each variable
!     ---------------------------------------------
      INTEGER  timenc_rank
!     ----------------------------
!     Parameters for each variable
!     ----------------------------
      PARAMETER (timenc_rank = 1)
!     ---------------
!     Variable shapes
!     ---------------
      INTEGER  timenc_dims(timenc_rank)
!     --------------
!     Data variables
!     --------------
      REAL(SZ), ALLOCATABLE ::  timenc(:) 
      REAL(SZ), ALLOCATABLE ::  eta22(:,:)                              ! global surface elevation
      REAL(SZ), ALLOCATABLE ::  eta11(:,:)                                ! boundary array
      REAL(SZ), ALLOCATABLE ::  etadd(:,:)                                ! boundary array
      REAL(SZ), ALLOCATABLE ::  uu22(:,:), vv22(:,:)                    ! global surface speed
      REAL(SZ), ALLOCATABLE ::  xnc(:), ync(:)                          ! coordinates
      REAL(SZ), ALLOCATABLE ::  xenc(:), yenc(:)                        ! coordinates station velocities
      REAL(SZ), ALLOCATABLE ::  xvnc(:), yvnc(:)                        ! coordinates station velocities
      REAL(SZ), ALLOCATABLE ::  pr22(:,:)                               ! pressure

      REAL(SZ), ALLOCATABLE ::  xpnc(:), ypnc(:)                        ! coordinates station pressure
      REAL(SZ), ALLOCATABLE ::  xwnc(:), ywnc(:)                        ! coordinates station wind speed


      REAL(SZ), ALLOCATABLE ::  WNCX22(:,:), WNCY22(:,:)                ! wind speed
      REAL(SZ), ALLOCATABLE ::  ET00NC(:,:)                             ! station surface elevation
      REAL(SZ), ALLOCATABLE ::  UU00nc(:,:), VV00nc(:,:)                ! station surface speed
      REAL(SZ), ALLOCATABLE ::  RMP00nc(:,:)                            ! station pressure
      REAL(SZ), ALLOCATABLE ::  RMU00nc(:,:), RMV00nc(:,:)              ! station wind speed

      DOUBLE PRECISION, PARAMETER ::  doubleval(1)=-99999.D0

      INTEGER, ALLOCATABLE ::  nbvvnc(:,:)                              ! boundary array
      INTEGER, ALLOCATABLE ::  nbounc(:)                                ! boundary array
      INTEGER, ALLOCATABLE ::  nvellnc(:)                               ! boundary array
      INTEGER, ALLOCATABLE ::  nopenc(:)                                ! boundary array
      INTEGER, ALLOCATABLE ::  nvdllnc(:)                               ! boundary array
      
      INTEGER, ALLOCATABLE ::  nodecodenc(:)                            ! boundary array
      INTEGER, ALLOCATABLE ::  noffnc(:)                                ! boundary array
      
      INTEGER nvelnc
      INTEGER netanc
!     -----------------
!     ADCIRC input variables
!     -----------------
      INTEGER fileunit
      LOGICAL fileexist  
      INTEGER irtype
      INTEGER ncmode
      
      CHARACTER*50, PARAMETER :: VAR61_LN='station sea surface elevation'
      CHARACTER*25, PARAMETER :: VAR61_ST='station_surface_elevation'

      CHARACTER*29, PARAMETER :: VAR62_LN='station velocity'
      CHARACTER*25, PARAMETER :: VAR62_ST='station_velocity'

      CHARACTER*21, PARAMETER :: VAR63_LN='sea surface elevation'
      CHARACTER*17, PARAMETER :: VAR63_ST='surface_elevation'

      CHARACTER*36, PARAMETER :: VAR71_LN='station atmospheric surface pressure'
      CHARACTER*28, PARAMETER :: VAR71_ST='station_atmospheric_pressure'
      CHARACTER*15, PARAMETER :: VAR71_UN='meters of water'

      CHARACTER*28, PARAMETER :: VAR73_LN='atmospheric surface pressure'
      CHARACTER*20, PARAMETER :: VAR73_ST='atmospheric_pressure'
      CHARACTER*15, PARAMETER :: VAR73_UN='meters of water'
      
      CHARACTER*14 , PARAMETER :: VAR_NOND='nondimensional'
      CHARACTER*4  , PARAMETER :: VAR_SECS='secs'
      CHARACTER*6  , PARAMETER :: VAR_METS='meters'
      CHARACTER*4  , PARAMETER :: VAR_FEET='feet'
      CHARACTER*2  , PARAMETER :: VAR_UP='up'
      CHARACTER*4  , PARAMETER :: VAR_DW='down'
      CHARACTER*5  , PARAMETER :: VAR_RG='right'
      CHARACTER*4  , PARAMETER :: VAR_EAST='east'
      CHARACTER*5  , PARAMETER :: VAR_NORTH='north'
      CHARACTER*13 , PARAMETER :: VAR_MPS='meters/second'
      CHARACTER*11 , PARAMETER :: VAR_FPS='feet/second'
      CHARACTER*19 , PARAMETER :: VAR_MPS_2='meters**2/second**2'
      CHARACTER*17 , PARAMETER :: VAR_FPS_2='feet**2/second**2'
      CHARACTER*23 , PARAMETER :: VAR_TUN='seconds since base_date'
      CHARACTER*12 , PARAMETER :: VAR_DEGE='degrees_east'
      CHARACTER*13 , PARAMETER :: VAR_DEGN='degrees_north'

         CONTAINS

!=======================================================================
!     Subroutine writenc63
!     Writes a netCDF file equivalent to ADCIRC's fort.63 ascii file
!
!     Dimensions:
!     node  : number of nodes in the horizontal grid
!     nele  : number of elements in the horizontal grid
!     nope  : number of elevation specified boundary forcing segments
!     neta  : total number of elevation specified boundary nodes
!     nbou  : number of normal flow specified boundary segments
!     nvel  : total number of normal flow specified boundary nodes
!     time  : number of time records
!
!     Variables :
!     time      : time
!     zeta      : Sea surface elevation 
!     x         : x-coordinate or longitude
!     y         : y-coordinate or latitude
!     element   : elements in grid
!     nvdll     : number of nodes in elevation boundary segment
!     nbdv      : node numbers on elevation specified boundary segment
!     nvell     : number of nodes in normal flow specified boundary segment
!     nbvv      : node numbers on normal flow boundary segment
!     ibtype    : boundary type
!     depth     : depth in each node
!     mindepth  : minimum depth in grid
!     maxdepth  : maximum depth in grid
!     minzetat  : minimum sea surface elevation for each output time
!     maxzetat  : maximum sea surface elevation for each output time
!     minzetan  : minimum sea surface elevation at each node
!     maxzetan  : maximum sea surface elevation at each node
!     minznodet : node of minimum sea surface elevation at each output time
!     maxznodet : node of maximum sea surface elevation at each output time
!     minztimen : time of minimum sea surface elevation at each node
!     maxztimen : time of maximum sea surface elevation at each node
!
!     Attributes:
!     model             : ADCIRC
!     version           : model version number
!     grid_type         : triangular
!     rundes            : description of simulation
!     agrid             : grid ID 
!     title             : description of the dataset
!     institution       : where the original data was produced
!     source            : method of production of the original data
!     history           : provides an audit trail for modifications to the original data
!     references        : published/web-based references that describe the data
!     comments          : miscellaneous information about the data
!     contact           : user's institution/ e-mail address
!     creation_date     : date the data was created
!     modification_date : date the data was modified
!     Parameters from fort.15 input file 
!
!     Revision history:
!
!     Date            Programmer                  Description of change
!     ----            ----------                  ----------------------
!     03/30/07  Cristina Forbes, PSGS @ UNC-IMS   Wrote original code
!     03/30/08  Cristina Forbes, PSGS @ UNC-IMS   Modified code for 
!                                                 globalio & hotstart 
!                                                 from binary files
!=======================================================================
      subroutine writenc63 (nout,fileunit,timenc_len, node_len,timenc,
     &                      zeta,nt, NDATA, rundes,runid, agrid, irtype,
     &                      dtdp, ihot, ics, nolibf,  nolifa,
     &                      nolica, nolicat, nwp,
     &                      ncor, ntip, nws, nramp, tau0,
     &                      statim, reftim, rnday, dramp,
     &                      a00, b00, c00, h0,
     &                      slam0, sfea0, ffactor,esl, cori,
     &                      ntif, nbfr, elem_len, nface_len, nm,
     &                      x,y,depth,
     &                      nope_len, neta_len, nbdv, nopenc,netanc_len,
     &                      nvdllnc,
     &                      nbou_len, nvel_len, nbvv,nbounc,nvelnc_len,
     &                      nvellnc,ibtypenc,
     &                      ncmode,
     &                      title, institution, source, history,
     &                      references, comments, host, convention,
     &                      contact,
     &                      base_date)

!--------------------------------------------
!     Force explicit declaration of all variables
!--------------------------------------------
      IMPLICIT NONE
!     -------------------
!     netCDF header file
!     -------------------
!      include 'netcdf.inc'
!     -------------------
!     Error status return
!     -------------------
      INTEGER  iret
!     ---------
!     netCDF id
!     ---------
      INTEGER  ncid
!     -------------
!     Dimension ids
!     -------------
      INTEGER  node_dim
      INTEGER  timenc_dim
      INTEGER  elem_dim
      INTEGER  nface_dim
      INTEGER  nope_dim
      INTEGER  neta_dim
      INTEGER  nbou_dim
      INTEGER  nvel_dim
      INTEGER  nvelnc_dim
      
      INTEGER  netanc_dim
!     -----------------
!     Dimension lengths
!     -----------------
      INTEGER  node_len
      INTEGER  timenc_len
      INTEGER  elem_len
      INTEGER  nface_len
      INTEGER  nope_len
      INTEGER  neta_len
      INTEGER  nbou_len
      INTEGER  nvel_len
      INTEGER  nvelnc_len
      
      INTEGER  netanc_len
!     ------------
!     Variable ids
!     ------------
      INTEGER  timenc_id
      INTEGER  zeta_id
      INTEGER  minzetat_id
      INTEGER  maxzetat_id
      INTEGER  minzetan_id
      INTEGER  maxzetan_id
      INTEGER  minznodet_id
      INTEGER  maxznodet_id
      INTEGER  minztimen_id
      INTEGER  maxztimen_id
      INTEGER  X_id
      INTEGER  Y_id
      INTEGER  DEPTH_id
      INTEGER  ELE_id
      INTEGER  mindepth_id
      INTEGER  maxdepth_id
      INTEGER  nbdv_id
      INTEGER  nbvv_id
      
      INTEGER  nvdll_id
      INTEGER  nbounc_id
      INTEGER  ibtypenc_id
      INTEGER  nvellnc_id
      
      INTEGER  nvdllnc_id
      INTEGER  nopenc_id

!     ---------------------------------------------
!     Rank (number of dimensions) for each variable
!     ---------------------------------------------
      INTEGER  timenc_rank
      INTEGER  zeta_rank
      INTEGER  minzetat_rank
      INTEGER  maxzetat_rank
      INTEGER  minzetan_rank
      INTEGER  maxzetan_rank
      INTEGER  minznodet_rank
      INTEGER  maxznodet_rank
      INTEGER  minztimen_rank
      INTEGER  maxztimen_rank
      INTEGER  depth_rank
      INTEGER  ele_rank
      INTEGER  mindepth_rank
      INTEGER  maxdepth_rank
      INTEGER  x_rank
      INTEGER  y_rank

      INTEGER  nbounc_rank
      INTEGER  nvdll_rank
      INTEGER  ibtypenc_rank
      INTEGER  nvellnc_rank

      INTEGER  nvdllnc_rank
      INTEGER  nopenc_rank

      INTEGER  nbdv_rank
      INTEGER  nbvv_rank
      
      PARAMETER (timenc_rank = 1)
      PARAMETER (zeta_rank = 2)
      PARAMETER (minzetat_rank = 1)
      PARAMETER (maxzetat_rank = 1)
      PARAMETER (minzetan_rank = 1)
      PARAMETER (maxzetan_rank = 1)
      PARAMETER (minznodet_rank = 1)
      PARAMETER (maxznodet_rank = 1)
      PARAMETER (minztimen_rank = 1)
      PARAMETER (maxztimen_rank = 1)
      PARAMETER (depth_rank = 1)
      PARAMETER (ele_rank = 2)
      PARAMETER (mindepth_rank = 0)
      PARAMETER (maxdepth_rank = 0)
      PARAMETER (x_rank = 1)
      PARAMETER (y_rank = 1)

      PARAMETER (nvdll_rank = 1)
      PARAMETER (nbounc_rank = 1)
      PARAMETER (ibtypenc_rank = 1)
      PARAMETER (nvellnc_rank = 1)

      PARAMETER (nvdllnc_rank=1)
      PARAMETER (nopenc_rank=1)

      PARAMETER (nbdv_rank = 2)
      PARAMETER (nbvv_rank = 2)
!     ---------------
!     Variable shapes
!     ---------------
      INTEGER  timenc_dims(timenc_rank)
      INTEGER  zeta_dims(zeta_rank)
      INTEGER  minzetat_dims(minzetat_rank)
      INTEGER  maxzetat_dims(maxzetat_rank)
      INTEGER  minzetan_dims(minzetan_rank)
      INTEGER  maxzetan_dims(maxzetan_rank)
      INTEGER  minznodet_dims(minznodet_rank)
      INTEGER  maxznodet_dims(maxznodet_rank)
      INTEGER  minztimen_dims(minztimen_rank)
      INTEGER  maxztimen_dims(maxztimen_rank)
      INTEGER  x_dims(x_rank)
      INTEGER  y_dims(y_rank)
      INTEGER  depth_dims(depth_rank)
      INTEGER  ele_dims(ele_rank)

      INTEGER  nvdll_dims(nvdll_rank)
      INTEGER  nbounc_dims(nbounc_rank)
      INTEGER  ibtypenc_dims(ibtypenc_rank)
      INTEGER  nvellnc_dims(nvellnc_rank)

      INTEGER  nvdllnc_dims(nvdllnc_rank)
      INTEGER  nopenc_dims(nopenc_rank)

      INTEGER  nbdv_dims(nbdv_rank)
      INTEGER  nbvv_dims(nbvv_rank)
!     --------------
!     Data variables
!     --------------
      CHARACTER*5 FILEFORT
      CHARACTER*10 FILENAME
      CHARACTER*3 FILEEXT
      INTEGER  t,ntif,nbfr,nt
      INTEGER  irtype, ihot, ics, nolibf, nolifa, nolica, nolicat, nwp
      INTEGER  ntip, nws, nramp,ncor
      INTEGER   nout
      INTEGER   element (nface_len , elem_len)
      INTEGER   nbdv (nope_len, neta_len)
      INTEGER   nbvv (nbou_len, nvel_len)

      INTEGER   nvellnc(nbou_len)
      INTEGER   nbounc (nvel_len)
      INTEGER   ibtype(nbou_len)
      INTEGER   ibtypenc(nbou_len)

      INTEGER   nvdllnc (nope_len)
      INTEGER   nopenc (neta_len)

      INTEGER   nm (elem_len, nface_len)
      INTEGER   fileunit
      REAL*8    timenc(timenc_len)
      REAL(sz)  x(node_len)
      REAL(sz)  y(node_len)
      REAL(sz)  depth(node_len)
      REAL(sz)  zeta (node_len, timenc_len)
      REAL(sz), ALLOCATABLE, SAVE ::  minzetat(:)
      REAL(sz), ALLOCATABLE, SAVE ::  maxzetat(:)
      REAL(sz), ALLOCATABLE, SAVE ::  minzetan(:)
      REAL(sz), ALLOCATABLE, SAVE ::  maxzetan(:)
      INTEGER , ALLOCATABLE, SAVE ::  minznodet(:)
      INTEGER , ALLOCATABLE, SAVE ::  maxznodet(:)
      REAL(sz), ALLOCATABLE, SAVE ::  minztimen(:)
      REAL(sz), ALLOCATABLE, SAVE ::  maxztimen(:)
      real*8  mindepth
      real*8  maxdepth
      REAL*8  dtdp, statim,reftim, rnday, dramp, tau0
      REAL*8  a00, b00, c00,h0
      REAL(sz) slam0, sfea0
      REAL(sz) slam0deg, sfea0deg
      REAL*8 ffactor,esl, cori
      CHARACTER(80) rundes, runid, agrid
      INTEGER, SAVE  :: FIRST = 1 
      INTEGER NDATA
      INTEGER, SAVE  :: NREC = 0
      INTEGER ncmode
!     -----------------
!     Attribute vectors
!     -----------------
      INTEGER intval(1),i,j
      INTEGER counti(1), starti(1), n
      INTEGER count(2), start(2)
      REAL    realval(1)
!     -----------------
!     Attribute general
!     -----------------
      CHARACTER(*) :: title, institution, source, history
      CHARACTER(*) :: references, comments, host, convention, contact
      CHARACTER(*) :: base_date
!     -----------------
!      date_string variables for time attribute
!     -----------------
      character date_string*40
      character now_date*8
      character big_ben*10
      character zone*5
      integer values(8)
      real*8 xmin, xmax
!     ---------------------------------
!      inquiring and reading variables
!     ---------------------------------

      INTEGER NDIMS, NVARS, NGATTS, UNLIMDIMID, RECID,NRECS
      CHARACTER*4 recname

!     -----------------
!     Save variables
!     -----------------

      save node_dim, elem_dim, nface_dim
      save nope_dim, neta_dim, nbou_dim, nvel_dim
      save  zeta_id, X_id, Y_id
      save minzetat_id, maxzetat_id, minzetan_id, maxzetan_id
      save minznodet_id, maxznodet_id, minztimen_id, maxztimen_id
      save DEPTH_id,ELE_id, mindepth_id, maxdepth_id, nbdv_id, nbvv_id
      save timenc_dim, timenc_id
      save mindepth, maxdepth
      
      save nvelnc_dim
      save nvdll_id, nbounc_id, ibtypenc_id, nvellnc_id
      save netanc_dim
      save nopenc_id,nvdllnc_id

!     ---------------------------
!     Opening fort.63 netCDF file
!     ---------------------------
      WRITE(6,*) "-------------------------------------------------"
      
      IF(FIRST.EQ.1.and.ncmode.ne.0) THEN
         ALLOCATE ( minzetan (node_len) )
         ALLOCATE ( maxzetan (node_len) )
         ALLOCATE ( minztimen (node_len) )
         ALLOCATE ( maxztimen (node_len) )
        
         do n=1,node_len
            minzetan  (n) =  9999.d0
            maxzetan  (n) = -9999.d0
            minztimen (n) =  9999.d0
            maxztimen (n) = -9999.d0
         end do
      ENDIF


!      IF(ncmode.eq.2) FIRST=0
!     --------------------------------
!     Switch order in array for NETCDF
!     --------------------------------
      do i=1, elem_len
         do j=1, nface_len
            element(j,i) = nm (i,j)
         end do
      end do
!     -----------------
!     Enter define mode
!     -----------------
         FILEFORT = 'fort.'
         FILEEXT = '.nc'
         WRITE(FILENAME,'(A5,I2,A3)')FILEFORT,FILEUNIT,FILEEXT
         WRITE(6,*) "Opening ",FILENAME,  " netCDF file ... "
                  
       IF((FIRST.EQ.1.AND.IHOT.EQ.0.AND.NCMODE.NE.2).OR.
     &    (FIRST.EQ.1.AND.IHOT.NE.0.AND.NCMODE.NE.2.AND.
     &     nout.EQ.-3))THEN
     
           iret = nf_create(FILENAME, NF_CLOBBER, ncid)
           CALL check_err(iret)
           FIRST=0
         
!     --------------------------------
!     Convert to degress for NETCDF I/O
!     --------------------------------
           SLAM0DEG=SLAM0
           SFEA0DEG=SFEA0
           IF (NCMODE.EQ.1) THEN
               SLAM0DEG=SLAM0*RAD2DEG
               SFEA0DEG=SFEA0*RAD2DEG
           ENDIF

!     =================
!     Define dimensions
!     =================
      iret = nf_def_dim(ncid, 'node', node_len, node_dim)
      CALL check_err(iret)
      iret = nf_def_dim(ncid, 'nele', elem_len, elem_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid, 'nface', 3, nface_dim)
      call check_err(iret)
      if(nope_len.ne.0) then
         iret = nf_def_dim(ncid, 'nope', nope_len, nope_dim)
         call check_err(iret)
      endif
      if(neta_len.ne.0) then
         iret = nf_def_dim(ncid, 'neta', neta_len, neta_dim)
         call check_err(iret)
      endif
      if(nbou_len.ne.0) then
         iret = nf_def_dim(ncid, 'nbou', nbou_len, nbou_dim)
         call check_err(iret)
      endif
      if(nvel_len.ne.0) then
         iret = nf_def_dim(ncid, 'nvel', nvel_len, nvel_dim)
         call check_err(iret)
      endif
      
      iret = nf_def_dim(ncid, 'time', nf_unlimited, timenc_dim)
      CALL check_err(iret)

!     ================
!     Define variables
!     ================
!     ----------------
!     Define time
!     ----------------

      timenc_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'time',NF_DOUBLE,timenc_rank, timenc_dims, 
     &                  timenc_id)
      CALL check_err(iret)
     
!     -------------------
!     Define X
!     -------------------
      X_dims(1) = node_dim
      iret = nf_def_var(ncid, 'x', NF_DOUBLE, X_rank, X_dims, X_id)
      CALL check_err(iret)
!     -------------------
!     Define Y coordinate
!     --------------------
      Y_dims(1) = node_dim
      iret = nf_def_var(ncid, 'y', NF_DOUBLE, Y_rank, Y_dims, Y_id)
      CALL check_err(iret)
!     -------------------
!     Define elements
!     --------------------
      ELE_dims(1) = nface_dim
      ELE_dims(2) = elem_dim
      iret = nf_def_var(ncid, 'element',NF_INT,ELE_rank,ELE_dims,ELE_id)
      CALL check_err(iret)
!     -----------------------------------------------------------------
!     Define elevation specified boundary forcing segments information
!     -----------------------------------------------------------------
      if(nope_len.ne.0.and.neta_len.ne.0) then
        nvdllnc_dims(1) = nope_dim      
        iret = nf_def_var(ncid, 'nvdll',NF_INT, nvdllnc_rank, 
     &                  nvdllnc_dims, nvdllnc_id)
        CALL check_err(iret)

        nopenc_dims(1) = neta_dim      
        iret = nf_def_var(ncid, 'nbdv',NF_INT, nopenc_rank, nopenc_dims, 
     &                  nopenc_id)
        CALL check_err(iret)

      endif

!     ----------------------------------------
!     Define normal flow boundary information
!     ----------------------------------------
      if(nbou_len.ne.0.and.nvel_len.ne.0) then

        nvellnc_dims(1) = nbou_dim      
        iret = nf_def_var(ncid, 'nvell',NF_INT, nvellnc_rank, 
     &                  nvellnc_dims, nvellnc_id)
        CALL check_err(iret)
        
        ibtypenc_dims(1) = nbou_dim      
        iret = nf_def_var(ncid, 'ibtype',NF_INT, ibtypenc_rank, 
     &                  ibtypenc_dims, ibtypenc_id)
        CALL check_err(iret)
        
        nbounc_dims(1) = nvel_dim      
        iret = nf_def_var(ncid, 'nbvv',NF_INT, nbounc_rank, nbounc_dims, 
     &                  nbounc_id)
        CALL check_err(iret)
     
       endif
!     -------------------
!     Define Z coordinate
!     --------------------
      DEPTH_dims(1) = node_dim
      iret = nf_def_var(ncid, 'depth',NF_DOUBLE, DEPTH_rank, DEPTH_dims, 
     &                  DEPTH_id)
      CALL check_err(iret)
!     --------------------
!     Define minimum depth
!     --------------------
      iret = nf_def_var(ncid, 'mindepth',NF_DOUBLE, mindepth_rank, 0, 
     &                  mindepth_id)
      CALL check_err(iret)
!     --------------------
!     Define maximum depth
!     --------------------
      iret = nf_def_var(ncid, 'maxdepth', NF_DOUBLE, maxdepth_rank, 0, 
     &                  maxdepth_id)
      CALL check_err(iret)
      
!     ----------------
!     Define zeta
!     ----------------
      zeta_dims(1) = node_dim
      zeta_dims(2) = timenc_dim
      
      iret = nf_def_var(ncid, 'zeta', NF_DOUBLE, zeta_rank, zeta_dims, 
     &                  ZETA_id)
      
      CALL check_err(iret)
!     ----------------------------
!     Define minimum zeta at times
!     ----------------------------
      minzetat_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'minzetat', NF_DOUBLE, minzetat_rank, 
     &                  minzetat_dims, minzetat_id)
      CALL check_err(iret)
!     ----------------------------
!     Define maximum zeta at times
!     ----------------------------
      maxzetat_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'maxzetat', NF_DOUBLE, maxzetat_rank, 
     &                  maxzetat_dims, maxzetat_id)
      CALL check_err(iret)

!     ----------------------------
!     Define minimum zeta at nodes
!     ----------------------------
      minzetan_dims(1) = node_dim
       iret = nf_def_var(ncid, 'minzetan', NF_DOUBLE, minzetan_rank, 
     &                  minzetan_dims, minzetan_id)
      CALL check_err(iret)
!     ----------------------------
!     Define maximum zeta at nodes
!     ----------------------------
      maxzetan_dims(1) = node_dim
      iret = nf_def_var(ncid, 'maxzetan', NF_DOUBLE, maxzetan_rank, 
     &                  maxzetan_dims, maxzetan_id)
      CALL check_err(iret)

!     ----------------------------
!     Define minimum zeta nodes at times
!     ----------------------------
      minznodet_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'minznodet', NF_INT, minznodet_rank, 
     &                  minznodet_dims, minznodet_id)
      CALL check_err(iret)
!     ----------------------------
!     Define maximum zeta nodes at times
!     ----------------------------
      maxznodet_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'maxznodet', NF_INT, maxznodet_rank, 
     &                  maxznodet_dims, maxznodet_id)
      CALL check_err(iret)
!     ----------------------------
!     Define minimum zeta nodes at nodes
!     ----------------------------
      minztimen_dims(1) = node_dim
      iret = nf_def_var(ncid, 'minztimen', NF_DOUBLE, minztimen_rank, 
     &                  minztimen_dims, minztimen_id)
      CALL check_err(iret)
!     ----------------------------
!     Define maximum zeta nodes at nodes
!     ----------------------------
      maxztimen_dims(1) = node_dim
      iret = nf_def_var(ncid, 'maxztimen', NF_DOUBLE, maxztimen_rank, 
     &                  maxztimen_dims, maxztimen_id)
      CALL check_err(iret)
       
!     -----------------
!     Assign attributes
!     -----------------
!     ----------------------------
!     Define longitude attributes
!     ----------------------------
      if(ics.eq.2) then
      
      	iret = nf_put_att_text(ncid,X_id,'long_name',9,'longitude')
      	       CALL check_err(iret)
      	iret = nf_put_att_text(ncid, X_id, 'standard_name', 9, 
     &                      'longitude')
      	       CALL check_err(iret)
      	iret = nf_put_att_text(ncid, X_id, 'units', 12, VAR_DEGE)
      	       CALL check_err(iret)
      	iret = nf_put_att_text(ncid, X_id, 'positive', 4, VAR_EAST)
      	       CALL check_err(iret)
      	iret = nf_put_att_double(ncid, X_id, 'missing_values',NF_DOUBLE, 
     &                       1, doubleval)
      	       CALL check_err(iret)
      	iret = nf_put_att_double(ncid, X_id, '_FillValue', NF_DOUBLE, 1, 
     &                       doubleval)
      	       CALL check_err(iret)      
      
!     ----------------------------
!     Define latitude attributes
!     ----------------------------
      	iret = nf_put_att_text(ncid, Y_id, 'long_name', 8, 'latitude')
      	       CALL check_err(iret)
      	iret = nf_put_att_text(ncid, Y_id, 'standard_name', 8, 
     &                       'latitude')
      	       CALL check_err(iret)
      	iret = nf_put_att_text(ncid, Y_id, 'units', 13, VAR_DEGN)
      	       CALL check_err(iret)
      	iret = nf_put_att_text(ncid, Y_id, 'positive', 5, VAR_NORTH)
      	       CALL check_err(iret)
       
      	iret = nf_put_att_double(ncid, Y_id, 'missing_values', 
     &       NF_DOUBLE, 1, doubleval)
      	       CALL check_err(iret)
      	iret = nf_put_att_double(ncid, Y_id, '_FillValue', NF_DOUBLE, 1, 
     &      doubleval)
      	       CALL check_err(iret)
 
       else
      
!     ------------------------------
!     Define x-coordinate attributes
!     ------------------------------
      	iret = nf_put_att_text(ncid, X_id, 'long_name', 22, 
     &       'Cartesian coordinate x')
      	       CALL check_err(iret)
      	iret = nf_put_att_text(ncid, X_id, 'standard_name',12, 
     &                      'x-coordinate')
      	       CALL check_err(iret)
      	       
!    -------------------------
!     determine variable units
!    -------------------------
      	if (G.LT.11.D0) then
            iret = nf_put_att_text(ncid, X_id, 'units', 6, VAR_METS) 
      	    CALL check_err(iret)
      	else
      	    iret = nf_put_att_text(ncid, X_id, 'units', 4, VAR_FEET) 
    	    CALL check_err(iret)
      	endif      	
        
      	iret = nf_put_att_text(ncid, X_id, 'positive', 5,  VAR_RG)
      	CALL check_err(iret)
      	iret = nf_put_att_double(ncid, X_id, 'missing_values',NF_DOUBLE,
     &                        1, doubleval)
      	CALL check_err(iret)
      	iret = nf_put_att_double(ncid, X_id, '_FillValue', NF_DOUBLE, 1, 
     &                       doubleval)
      	CALL check_err(iret)
!     ------------------------------
!     Define y-coordinate attributes
!     ------------------------------     
      	iret = nf_put_att_text(ncid, Y_id, 'long_name', 22, 
     &                        'Cartesian coordinate y')
      	CALL check_err(iret)
      	iret = nf_put_att_text(ncid, Y_id, 'standard_name', 12, 
     &                       'y-coordinate')
      	CALL check_err(iret)
      	if (G.LT.11.D0) then       
     	  iret = nf_put_att_text(ncid, Y_id, 'units', 6,  VAR_METS)
      	  CALL check_err(iret)
      	else 
     	  iret = nf_put_att_text(ncid, Y_id, 'units', 4,  VAR_FEET)
      	  CALL check_err(iret)      	
      	endif   
      	
      	iret = nf_put_att_text(ncid, Y_id, 'positive', 2, VAR_UP)
      	CALL check_err(iret)
      	iret = nf_put_att_double(ncid, Y_id, 'missing_values',NF_DOUBLE, 
     &                           1, doubleval)
      	CALL check_err(iret)
      	iret = nf_put_att_double(ncid, Y_id, '_FillValue', NF_DOUBLE, 1, 
     &                           doubleval)
      	CALL check_err(iret)
      
      endif
!     ------------------------------
!     Define bathymetry attributes
!     ------------------------------
      iret = nf_put_att_text(ncid, DEPTH_id,'long_name',10,'bathymetry')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, DEPTH_id, 'standard_name', 5,'depth')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, DEPTH_id, 'units', 6, VAR_METS)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, DEPTH_id, 'missing_values', 
     &                      NF_DOUBLE, 1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, DEPTH_id, '_FillValue', NF_DOUBLE,
     &                        1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, DEPTH_id, 'positive', 4, VAR_DW)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, ELE_id, 'long_name', 7, 'element')
      CALL check_err(iret)
      intval(1) = -1
      iret = nf_put_att_text(ncid, ELE_id, 'units', 14, 
     &                      VAR_NOND)
      CALL check_err(iret)
      
!     ------------------------------------
!     Define minimum bathymetry attributes
!     ------------------------------------
      iret = nf_put_att_text(ncid, mindepth_id, 'long_name', 18, 
     &                       'minimum bathymetry')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, mindepth_id, 'standard_name', 13,
     &                'minimum_depth')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, mindepth_id, 'units', 6, VAR_METS)
      CALL check_err(iret)
      
!     ------------------------------------
!     Define maximum bathymetry attributes
!     ------------------------------------
      iret = nf_put_att_text(ncid, maxdepth_id, 'long_name', 18, 
     &                       'maximum bathymetry')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, maxdepth_id, 'standard_name', 13,
     &                'maximum_depth')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, maxdepth_id, 'units', 6, VAR_METS)
      CALL check_err(iret)
      
!     -----------------------------------------------
!     Define minimum sea surface elevation attributes
!     -----------------------------------------------
      iret = nf_put_att_text(ncid, minzetat_id, 'long_name', 50, 
     &      'minimum sea surface elevation for each output time')
      CALL check_err(iret)
      if (G.LT.11.D0) then
         iret = nf_put_att_text(ncid, minzetat_id, 'units', 6, VAR_METS)
         CALL check_err(iret)
      else
         iret = nf_put_att_text(ncid, minzetat_id, 'units', 6, VAR_FEET)
         CALL check_err(iret)
      endif
      
      iret = nf_put_att_text(ncid, maxzetat_id, 'long_name', 50, 
     &       'maximum sea surface elevation for each output time')
      CALL check_err(iret)
      if (G.LT.11.D0) then
         iret = nf_put_att_text(ncid, maxzetat_id, 'units', 6, VAR_METS)
         CALL check_err(iret)
      else
         iret = nf_put_att_text(ncid, maxzetat_id, 'units', 4, VAR_FEET)
         CALL check_err(iret)
      endif
!     -----------------------------------------------
!     Define maximum sea surface elevation attributes
!     -----------------------------------------------
      iret = nf_put_att_text(ncid, minzetan_id, 'long_name', 42, 
     &      'minimum sea surface elevation at each node')
      CALL check_err(iret)
      if (G.LT.11.D0) then
         iret = nf_put_att_text(ncid, minzetan_id, 'units', 6, VAR_METS)
         CALL check_err(iret)
      else
         iret = nf_put_att_text(ncid, minzetan_id, 'units', 4, VAR_FEET)
         CALL check_err(iret)
      endif      
      iret = nf_put_att_text(ncid, maxzetan_id, 'long_name', 42, 
     &       'maximum sea surface elevation at each node')
      CALL check_err(iret)
      if (G.LT.11.D0) then
         iret = nf_put_att_text(ncid, maxzetan_id, 'units', 6, VAR_METS)
         CALL check_err(iret)
      else      
         iret = nf_put_att_text(ncid, maxzetan_id, 'units', 6, VAR_FEET)
         CALL check_err(iret)
      endif      
!     -------------------------------------------------------------
!     Define max and min in time sea surface elevation attributes 
!     -------------------------------------------------------------
      iret = nf_put_att_text(ncid, minznodet_id, 'long_name', 57, 
     &    'node of minimum sea surface elevation at each output time')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, minznodet_id, 'units',14, 
     &                      VAR_NOND)
      CALL check_err(iret)
      
      iret = nf_put_att_text(ncid, maxznodet_id, 'long_name', 57, 
     &    'node of maximum sea surface elevation at each output time')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, maxznodet_id, 'units', 14, 
     &                      VAR_NOND)
      CALL check_err(iret)
!     -------------------------------------------------------------
!     Define max and min in nodes sea surface elevation attributes
!     -------------------------------------------------------------
      iret = nf_put_att_text(ncid, minztimen_id, 'long_name', 50, 
     &     'time of minimum sea surface elevation at each node')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, minztimen_id, 'units', 4, VAR_SECS)
      CALL check_err(iret)
      
      iret = nf_put_att_text(ncid, maxztimen_id, 'long_name', 50, 
     &     'time of maximum sea surface elevation at each node')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, maxztimen_id, 'units', 4, VAR_SECS)
      CALL check_err(iret)
      

      iret = nf_put_att_text(ncid, nopenc_id, 'long_name', 53, 
     &       'node numbers on elevation specified boundary segments')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, nopenc_id, 'units', 14, 
     &                      VAR_NOND)
      CALL check_err(iret)
      
      iret = nf_put_att_text(ncid, nvdllnc_id, 'long_name', 60, 
     &'number of nodes in each elevation specified boundary segment')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, nvdllnc_id, 'units', 14, 
     &                      VAR_NOND)
      CALL check_err(iret)
!     ---------------------------------------------------------------
!     Define node numbers on normal flow boundary segments attributes
!     ----------------------------------------------------------------  

      iret = nf_put_att_text(ncid, nbounc_id, 'long_name', 59, 
     & 'node numbers in each normal flow specified boundary segment')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, nbounc_id, 'units', 14, 
     &                      VAR_NOND)
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, ibtypenc_id, 'long_name', 13, 
     &       'boundary type')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, ibtypenc_id, 'units', 14, 
     &                      VAR_NOND)
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, nvellnc_id, 'long_name', 62, 
     &'number of nodes in each normal flow specified boundary segment')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, nvellnc_id, 'units', 14, 
     &                      VAR_NOND)
      CALL check_err(iret)

!    -----------------------
!     Define time attributes
!    -----------------------
      iret = nf_put_att_text(ncid,timenc_id,'long_name',10,'model time')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid,timenc_id, 'standard_name', 4, 'time')
      CALL check_err(iret)
!      iret = nf_put_att_text(ncid, timenc_id, 'units', 7, 'seconds')
!      CALL check_err(iret)
      iret = nf_put_att_text(ncid, timenc_id, 'units', 23,VAR_TUN)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, timenc_id, 'base_date',
     .      (LEN(TRIM(base_date(1:INDEX(base_date,"!")-1)))), base_date)
      CALL check_err(iret)

!    --------------------------------------------
!     Define sea surface elevation attributes
!    --------------------------------------------
      
      iret = nf_put_att_text(ncid, zeta_id,'long_name',21, VAR63_LN)
     &                   
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, zeta_id, 'standard_name', 17, 
     &                      VAR63_ST)
      CALL check_err(iret)
      
      if (G.LT.11.D0) then
          iret = nf_put_att_text(ncid, zeta_id, 'units', 6, VAR_METS)
          CALL check_err(iret)
      else
         iret = nf_put_att_text(ncid, zeta_id, 'units', 4, VAR_FEET)
         CALL check_err(iret)
      endif      
      
      iret = nf_put_att_DOUBLE(ncid,zeta_id,'missing_values',NF_DOUBLE,
     &1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_DOUBLE(ncid, zeta_id, '_FillValue', NF_DOUBLE, 
     &                      1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, zeta_id,'positive', 2, VAR_UP)
      CALL check_err(iret)

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
      intval(1) = 1
      iret = nf_put_att_int(ncid, NCGLOBAL, 'irtype', NF_INT, 1, irtype)
      CALL check_err(iret)
      intval(1) = 0
      iret = nf_put_att_int(ncid, NCGLOBAL, 'ihot', NF_INT, 1, ihot)
      CALL check_err(iret)
      intval(1) = 1
      iret = nf_put_att_int(ncid, NCGLOBAL, 'ics', NF_INT, 1, ics)
      CALL check_err(iret)
      intval(1) = 1
      iret = nf_put_att_int(ncid, NCGLOBAL, 'nolibf', NF_INT, 1, nolibf)
      CALL check_err(iret)
      intval(1) = 1
      iret = nf_put_att_int(ncid, NCGLOBAL, 'nolifa', NF_INT, 1, nolifa)
      CALL check_err(iret)
      intval(1) = 1
      iret = nf_put_att_int(ncid, NCGLOBAL, 'nolica', NF_INT, 1, nolica)
      CALL check_err(iret)
      intval(1) = 1
      iret = nf_put_att_int(ncid, NCGLOBAL, 'nolicat', NF_INT, 1, 
     &                      nolicat)
      CALL check_err(iret)
      intval(1) = 0
      iret = nf_put_att_int(ncid, NCGLOBAL, 'nwp', NF_INT, 1, nwp)
      CALL check_err(iret)
      intval(1) = 0
      iret = nf_put_att_int(ncid, NCGLOBAL, 'ncor', NF_INT, 1, ncor)
      CALL check_err(iret)
      intval(1) = 0
      iret = nf_put_att_int(ncid, NCGLOBAL, 'ntip', NF_INT, 1, ntip)
      CALL check_err(iret)
      intval(1) = 0
      iret = nf_put_att_int(ncid, NCGLOBAL, 'nws', NF_INT, 1, nws)
      CALL check_err(iret)
      intval(1) = 1
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
      intval(1) = 5
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
      iret = nf_put_att_double(ncid, NCGLOBAL, 'ffactor', NF_DOUBLE, 1, 
     &                         ffactor)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, NCGLOBAL, 'esl', NF_DOUBLE, 1, 
     &                         esl)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, NCGLOBAL, 'cori', NF_DOUBLE, 1, 
     &                         cori)
      CALL check_err(iret)
      iret = nf_put_att_int(ncid, NCGLOBAL, 'ntif', NF_INT, 1, ntif)
      CALL check_err(iret)
      iret = nf_put_att_int(ncid, NCGLOBAL, 'nbfr', NF_INT, 1, nbfr)
      CALL check_err(iret)
!     ======================
!     Leave define mode
!     ======================
      iret = nf_enddef(ncid)
      CALL check_err(iret)
!     ======================
!     Store static variables
!     ======================
!     -------
!     Store X
!     -------
      iret = nf_put_var_double(ncid, X_id, X)
      CALL check_err(iret)
!     -------
!     Store Y
!     -------
      iret = nf_put_var_double(ncid, Y_id, Y)
      CALL check_err(iret)
!     ------------
!     Store depth
!     ------------
      iret = nf_put_var_double(ncid, DEPTH_id, depth)
      CALL check_err(iret)
!     --------------
!     Store elements
!     --------------
      count(1)=nface_len
      count(2)=elem_len
      start(1)=1
      start(2)=1
      iret=nf_put_vara_int(ncid,ele_id,start,count,element)
      call check_err(iret)
!     ------------------------------------
!     Store elevation boundary information
!     -------------------------------------
      iret = nf_put_var_int(ncid, nopenc_id, nopenc)
      CALL check_err(iret)
      iret = nf_put_var_int(ncid, nvdllnc_id, nvdllnc)
      CALL check_err(iret)
!     --------------------------------------
!     Store normal flow boundary information
!     --------------------------------------
      iret = nf_put_var_int(ncid, nbounc_id, nbounc)
      CALL check_err(iret)
      iret = nf_put_var_int(ncid, ibtypenc_id, ibtypenc)
      CALL check_err(iret)
      iret = nf_put_var_int(ncid, nvellnc_id, nvellnc)
      CALL check_err(iret)
!     --------------------
!     Store maximum depth
!     --------------------
      maxdepth = maxval (depth)
      iret = nf_put_var_double(ncid, maxdepth_id, maxdepth)
      CALL check_err(iret)
!     --------------------
!     Store minimum depth
!     --------------------
      mindepth = minval (depth)
      iret = nf_put_var_double(ncid, mindepth_id, mindepth)
      CALL check_err(iret)
!     =====================
!     Open old netcdf file
!     =====================
      ELSE
          iret = NF_OPEN(FILENAME, NF_WRITE, ncid)
          CALL check_err(iret)   
!          FIRST=0
!     ------------------
!     Inquire variables
!     ------------------
          iret=nf_inq_unlimdim(ncid,recid)
          CALL check_err(iret)   
          iret=nf_inq_dimlen(ncid,recid,nrec)
          CALL check_err(iret)   
      ENDIF
!     ------------------
!     Update Global Attributes
!     ------------------
      IF(IHOT.NE.0.AND.FIRST.EQ.1) THEN
         iret = NF_REDEF (ncid)
         CALL check_err(iret)
         iret = nf_put_att_int(ncid, NCGLOBAL, 'ihot', NF_INT, 1, ihot)
         CALL check_err(iret)
         iret = nf_put_att_int(ncid, NCGLOBAL, 'nramp', NF_INT, 1,nramp)
         CALL check_err(iret)
         iret = nf_put_att_double(ncid, NCGLOBAL, 'tau0', NF_DOUBLE, 1, 
     &                         tau0)
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
         iret = nf_put_att_double(ncid, NCGLOBAL, 'ffactor', NF_DOUBLE, 1, 
     &                         ffactor)
         CALL check_err(iret)
         iret = nf_put_att_double(ncid, NCGLOBAL, 'esl', NF_DOUBLE, 1, 
     &                         esl)
         CALL check_err(iret)
         CALL date_and_time(now_date,big_ben,zone,values)
         WRITE(date_string,71) values(1),values(2),values(3),
     &         values(5),values(6),values(7),   (values(4))/60
         iret = nf_put_att_text(ncid, NF_GLOBAL,'modification_date',
     &          len_trim(date_string),date_string)
         CALL check_err(iret)
         iret = NF_ENDDEF (ncid)
         CALL check_err(iret)
       ENDIF

      IF(ncmode.eq.0) go to 9999
      IF(FIRST.EQ.1.and.
     & ((ncmode.ne.0.and.IHOT.NE.0).or.(ncmode.eq.2.and.IHOT.eq.0).or.
     &  (ncmode.eq.2.and.IHOT.ne.0.and.nout.eq.-3)))THEN
         iret=nf_inq_varid(ncid,"time"     ,timenc_id)
         iret=nf_inq_varid(ncid,"zeta"     ,zeta_id)
         iret=nf_inq_varid(ncid,"minzetat" ,minzetat_id)
         iret=nf_inq_varid(ncid,"maxzetat" ,maxzetat_id)
         iret=nf_inq_varid(ncid,"minzetan" ,minzetan_id)
         iret=nf_inq_varid(ncid,"maxzetan" ,maxzetan_id)
         iret=nf_inq_varid(ncid,"minznodet",minznodet_id)
         iret=nf_inq_varid(ncid,"maxznodet",maxznodet_id)
         iret=nf_inq_varid(ncid,"minztimen",minztimen_id)
         iret=nf_inq_varid(ncid,"maxztimen",maxztimen_id)
      ENDIF
      IF(IHOT.NE.0.AND.FIRST.EQ.1.and.ncmode.ne.0.and.nout.eq.3)THEN
         iret = nf_get_var_double(ncid, minzetan_id, minzetan)
         CALL check_err(iret)
         iret = nf_get_var_double(ncid, maxzetan_id, maxzetan)
         CALL check_err(iret)
         iret = nf_get_var_double(ncid, minztimen_id, minztimen)
         CALL check_err(iret)
         iret = nf_get_var_double(ncid, maxztimen_id, maxztimen)
         CALL check_err(iret)
      ENDIF
      
      FIRST=0
      NREC = NREC + 1
      
!     ----------------------------------------
!     Allocating memory for min and max arrays
!     ----------------------------------------
      
      IF(.NOT. ALLOCATED(minzetat ))ALLOCATE( minzetat(timenc_len))
      IF(.NOT. ALLOCATED(maxzetat ))ALLOCATE( maxzetat(timenc_len))
      IF(.NOT. ALLOCATED(minzetan ))ALLOCATE( minzetan(node_len))
      IF(.NOT. ALLOCATED(maxzetan ))ALLOCATE( maxzetan(node_len))

      IF(.NOT. ALLOCATED(minznodet))ALLOCATE( minznodet(timenc_len))
      IF(.NOT. ALLOCATED(maxznodet))ALLOCATE( maxznodet(timenc_len))
      IF(.NOT. ALLOCATED(minztimen))ALLOCATE( minztimen(  node_len))
      IF(.NOT. ALLOCATED(maxztimen))ALLOCATE( maxztimen(  node_len))
      
!     ----------
!     Store time
!     ----------
      starti(1)=nrec
      counti(1)=timenc_len
      iret = nf_put_vara_DOUBLE(ncid, timenc_id, starti,counti, timenc)
      CALL check_err(iret)
!     WRITE(6,*)"store time "

!     ---------------------------
!     Store zeta
!     ---------------------------
      count(1)=node_len
      count(2)=timenc_len
      start(1)=1
      start(2)=nrec
      iret = nf_put_vara_DOUBLE(ncid, zeta_id, start, count, zeta)
      CALL check_err(iret)

!     -------------------------------------
!     Calculating zeta minimum in each node
!     -------------------------------------
      do n=1,node_len
 	 if(zeta (n,1).lt.minzetan(n) ) then
 	    minzetan (n) = zeta (n,1)
	    minztimen(n) = timenc(1)
 	 endif
      end do
!     --------------------------------------
!     Calculating zeta minimum for each time
!     --------------------------------------
      minzetat (1) = 9999.d0
      do n=1,node_len
 	 if(zeta (n,1).gt.-9999.d0.and.zeta(n,1).lt. minzetat (1) )then 
 	    minzetat (1) = zeta (n,1)
 	    minznodet(1) = n
 	 endif
      end do

      iret=nf_put_vara_DOUBLE(ncid, minzetat_id,(/nrec/),(/1/),minzetat)
      CALL check_err(iret)
      iret=nf_put_vara_INT(ncid, minznodet_id,(/nrec/),(/1/),minznodet)
      CALL check_err(iret)

!     -------------------------------------
!     Calculating zeta maximum in each node
!     -------------------------------------
      do n=1,node_len
         if(zeta (n,1).gt.maxzetan(n) ) then
            maxzetan (n) = zeta (n,1)
            maxztimen(n) = timenc(1)
         endif
      end do
!     -------------------------------------
!     Calculating zeta maximum for each time
!     -------------------------------------
      maxzetat (1) = -9999.d0
      do n=1,node_len
 	 if(zeta(n,1).gt.maxzetat(1))then
 	    maxzetat (1) = zeta (n,1)
 	    maxznodet(1) = n
 	 endif
      end do

      iret=nf_put_vara_DOUBLE(ncid,maxzetat_id,(/nrec/),(/1/),maxzetat)
      CALL check_err(iret)
      iret=nf_put_vara_INT(ncid,maxznodet_id,(/nrec/),(/1/),maxznodet)
      CALL check_err(iret)

!     ---------------------------
!     Closing fort.X3 netCDF file
!     ---------------------------
      IF(NREC.EQ.NDATA) THEN
         iret = nf_put_var_double(ncid, minzetan_id, minzetan)
         CALL check_err(iret)
         iret = nf_put_var_double(ncid, maxzetan_id, maxzetan)
         CALL check_err(iret)
         iret = nf_put_var_double(ncid, minztimen_id, minztimen)
         CALL check_err(iret)
         iret = nf_put_var_double(ncid, maxztimen_id, maxztimen)
         CALL check_err(iret)
      ENDIF
      
      WRITE(6,*) "Closing ", FILENAME, " netCDF file ... "
      WRITE(6,*) "-----------------------------------------------------"
      
9999  iret = nf_close(ncid)
      CALL check_err(iret)
      WRITE(6,*)'end subroutine writenc63'
     
      end subroutine writenc63
!=======================================================================
!     -----------------------
!     Subroutine check errors
!     -----------------------
      subroutine check_err(iret)
      INTEGER iret
      if (iret .ne. NF_NOERR) then
         print *, nf_strerror(iret)
         stop
      endif
      end subroutine check_err
!=======================================================================
!     Subroutine writenc64
!     Writes a netCDF file equivalent to ADCIRC's fort.64 ascii file
!
!     Dimensions:
!     node  : number of nodes in the horizontal grid
!     nele  : number of elements in the horizontal grid
!     nope  : number of elevation specified boundary forcing segments
!     neta  : total number of elevation specified boundary nodes
!     nbou  : number of normal flow specified boundary segments
!     nvel  : total number of normal flow specified boundary nodes
!     time  : number of time records
!
!     Variables :
!     time      : time
!     u-vel     : u velocity in x or east/west direction 
!     v-vel     : v velocity in y or north/south direction
!     x         : x-coordinate or longitude
!     y         : y-coordinate or latitude
!     element   : elements in grid
!     nvdll     : number of nodes in elevation boundary segment
!     nbdv      : node numbers on elevation specified boundary segment
!     nvell     : number of nodes in normal flow specified boundary segment
!     nbvv      : node numbers on normal flow boundary segment
!     ibtype    : boundary type
!     depth     : depth in each node
!     mindepth  : minimum depth in grid
!     maxdepth  : maximum depth in grid
!     minuvelt  : minimum vertically averaged u at each time
!     maxuvelt  : maximum vertically averaged u at each time
!     minvvelt  : minimum vertically averaged v at each time
!     maxvvelt  : maximum vertically averaged v at each time
!     minuveln  : minimum vertically averaged u at each node
!     maxuveln  : maximum vertically averaged u at each node
!     minvveln  : minimum vertically averaged v at each node
!     maxvveln  : maximum vertically averaged v at each node
!     minunodet : node where minimum u occurs at each output time
!     maxunodet : node where maximum u occurs at each output time
!     minvnodet : node where minimum v occurs at each output time
!     maxvnodet : node where maximum v occurs at each output time
!     minutimen : time when minimum u occurs at each node
!     maxutimen : time when maximum u occurs at each node
!     minvtimen : time when minimum v occurs at each node
!     maxvtimen : time when maximum v occurs at each node
!
!     Attributes:
!     model             : ADCIRC
!     version           : model version number
!     grid_type         : triangular
!     rundes            : description of simulation
!     agrid             : grid ID 
!     title             : description of the dataset
!     institution       : where the original data was produced
!     source            : method of production of the original data
!     history           : provides an audit trail for modifications to the original data
!     references        : published or web-based references that describe the data or methods used to produce it
!     comments          : miscellaneous information about the data
!     contact           : user's institution/ e-mail address
!     creation_date     : date the data was created
!     modification_date : date the data was modified
!     Parameters from fort.15 input file 
!
!     Revision history:
!
!     Date            Programmer                  Description of change
!     ----            ----------                  ----------------------
!     03/30/07  Cristina Forbes, PSGS @ UNC-IMS   Wrote original code
!     03/30/08  Cristina Forbes, PSGS @ UNC-IMS   Modified code for 
!                                                 globalio & hotstart 
!                                                 from binary files
!=======================================================================
      subroutine writenc64 (nout,fileunit,timenc_len,node_len,timenc,
     &                      nt, NDATA, u,v,
     &                      rundes, runid, agrid, irtype,
     &                      dtdp, ihot, ics, nolibf, nolifa,
     &                      nolica, nolicat, nwp, 
     &                      ncor, ntip, nws, nramp, tau0,
     &                      statim, reftim, rnday, dramp,
     &                      a00, b00, c00, h0,
     &                      slam0, sfea0, ffactor,esl, cori,
     &                      ntif, nbfr, elem_len, nface_len, nm,
     &                      x,y,depth,
     &                      nope_len, neta_len, nbdv,nopenc,netanc_len,
     &                      nvdllnc, 
     &                      nbou_len, nvel_len, nbvv,nbounc,nvelnc_len,
     &                      nvellnc,ibtypenc,
     &                      ncmode,
     &                      title, institution, source, history,
     &                      references, comments, host, convention,
     &                      contact,
     &                      base_date)


!--------------------------------------------
!     Force explicit declaration of all variables
!--------------------------------------------
      IMPLICIT NONE
!     -------------------
!     Error status return
!     -------------------
      INTEGER  iret
!     ---------
!     netCDF id
!     ---------
      INTEGER  ncid
!     -------------
!     Dimension ids
!     -------------
      INTEGER  node_dim
      INTEGER  timenc_dim
      INTEGER  elem_dim
      INTEGER  nface_dim
      INTEGER  nope_dim
      INTEGER  neta_dim
      INTEGER  nbou_dim
      INTEGER  nvel_dim
      INTEGER  nvelnc_dim
      INTEGER  netanc_dim
!     -----------------
!     Dimension lengths
!     -----------------
      INTEGER  node_len
      INTEGER  timenc_len
      INTEGER  elem_len
      INTEGER  nface_len
      INTEGER  nope_len
      INTEGER  neta_len
      INTEGER  nbou_len
      INTEGER  nvel_len
      INTEGER  nvelnc_len
      INTEGER  netanc_len
!     ------------
!     Variable ids
!     ------------
      INTEGER  timenc_id
      INTEGER  X_id
      INTEGER  Y_id
      INTEGER  DEPTH_id
      INTEGER  ELE_id
      INTEGER  mindepth_id
      INTEGER  maxdepth_id
      INTEGER  nbdv_id
      INTEGER  nbvv_id
      INTEGER  u_id
      INTEGER  v_id
      INTEGER  minut_id
      INTEGER  maxut_id
      INTEGER  minvt_id
      INTEGER  maxvt_id
      INTEGER  minun_id
      INTEGER  maxun_id
      INTEGER  minvn_id
      INTEGER  maxvn_id
      
      INTEGER  minunodet_id
      INTEGER  maxunodet_id
      INTEGER  minutimen_id
      INTEGER  maxutimen_id
      INTEGER  minvnodet_id
      INTEGER  maxvnodet_id
      INTEGER  minvtimen_id
      INTEGER  maxvtimen_id

      INTEGER  nvdll_id
      INTEGER  nbounc_id
      INTEGER  ibtypenc_id
      INTEGER  nvellnc_id
      INTEGER  nvdllnc_id
      INTEGER  nopenc_id
      
!     ---------------------------------------------
!     Rank (number of dimensions) for each variable
!     ---------------------------------------------
      INTEGER  timenc_rank
      INTEGER  depth_rank
      INTEGER  ele_rank
      INTEGER  mindepth_rank
      INTEGER  maxdepth_rank
      INTEGER  x_rank
      INTEGER  y_rank
      INTEGER  nvdll_rank
      INTEGER  nbdv_rank
      INTEGER  nbvv_rank
      INTEGER  u_rank
      INTEGER  v_rank
      INTEGER  minut_rank
      INTEGER  maxut_rank
      INTEGER  minvt_rank
      INTEGER  maxvt_rank
      INTEGER  minun_rank
      INTEGER  maxun_rank
      INTEGER  minvn_rank
      INTEGER  maxvn_rank

      INTEGER  nbounc_rank
      INTEGER  ibtypenc_rank
      INTEGER  nvellnc_rank
      INTEGER  nvdllnc_rank
      INTEGER  nopenc_rank

      INTEGER  minunodet_rank
      INTEGER  maxunodet_rank
      INTEGER  minutimen_rank
      INTEGER  maxutimen_rank
      
      INTEGER  minvnodet_rank
      INTEGER  maxvnodet_rank
      INTEGER  minvtimen_rank
      INTEGER  maxvtimen_rank

      PARAMETER (timenc_rank = 1)
      PARAMETER (depth_rank = 1)
      PARAMETER (ele_rank = 2)
      PARAMETER (mindepth_rank = 0)
      PARAMETER (maxdepth_rank = 0)
      PARAMETER (x_rank = 1)
      PARAMETER (y_rank = 1)
      PARAMETER (nbdv_rank = 2)
      PARAMETER (nbvv_rank = 2)
      PARAMETER (u_rank = 2)
      PARAMETER (v_rank = 2)
      PARAMETER (minut_rank = 1)
      PARAMETER (maxut_rank = 1)
      PARAMETER (minvt_rank = 1)
      PARAMETER (maxvt_rank = 1)
      PARAMETER (minun_rank = 1)
      PARAMETER (maxun_rank = 1)
      PARAMETER (minvn_rank = 1)
      PARAMETER (maxvn_rank = 1)

      PARAMETER (minunodet_rank = 1)
      PARAMETER (maxunodet_rank = 1)
      PARAMETER (minutimen_rank = 1)
      PARAMETER (maxutimen_rank = 1)

      PARAMETER (minvnodet_rank = 1)
      PARAMETER (maxvnodet_rank = 1)
      PARAMETER (minvtimen_rank = 1)
      PARAMETER (maxvtimen_rank = 1)

      PARAMETER (nvdll_rank = 1)
      PARAMETER (nbounc_rank = 1)
      PARAMETER (ibtypenc_rank = 1)
      PARAMETER (nvellnc_rank = 1)
      PARAMETER (nvdllnc_rank=1)
      PARAMETER (nopenc_rank=1)

!     ---------------
!     Variable shapes
!     ---------------
      INTEGER  timenc_dims(timenc_rank)
      INTEGER  x_dims(x_rank)
      INTEGER  y_dims(y_rank)
      INTEGER  depth_dims(depth_rank)
      INTEGER  ele_dims(ele_rank)
      INTEGER  nbdv_dims(nbdv_rank)
      INTEGER  nbvv_dims(nbvv_rank)
      INTEGER  u_dims(u_rank)
      INTEGER  v_dims(v_rank)
      INTEGER  minut_dims(minut_rank)
      INTEGER  maxut_dims(maxut_rank)
      INTEGER  minvt_dims(minvt_rank)
      INTEGER  maxvt_dims(maxvt_rank)
      INTEGER  minun_dims(minun_rank)
      INTEGER  maxun_dims(maxun_rank)
      INTEGER  minvn_dims(minvn_rank)
      INTEGER  maxvn_dims(maxvn_rank)

      INTEGER  nvdll_dims(nvdll_rank)
      INTEGER  nbounc_dims(nbounc_rank)
      INTEGER  ibtypenc_dims(ibtypenc_rank)
      INTEGER  nvellnc_dims(nvellnc_rank)
      INTEGER  nvdllnc_dims(nvdllnc_rank)
      INTEGER  nopenc_dims(nopenc_rank)

      INTEGER  minunodet_dims(minunodet_rank)
      INTEGER  maxunodet_dims(maxunodet_rank)
      INTEGER  minutimen_dims(minutimen_rank)
      INTEGER  maxutimen_dims(maxutimen_rank)

      INTEGER  minvnodet_dims(minvnodet_rank)
      INTEGER  maxvnodet_dims(maxvnodet_rank)
      INTEGER  minvtimen_dims(minvtimen_rank)
      INTEGER  maxvtimen_dims(maxvtimen_rank)

!     --------------
!     Data variables
!     --------------
      CHARACTER*5 FILEFORT
      CHARACTER*10 FILENAME
      CHARACTER*3 FILEEXT
      INTEGER   t,ntif,nbfr,nt
      INTEGER   irtype, ihot, ics, nolibf, nolifa, nolica, nolicat, nwp
      INTEGER   ntip, nws, nramp,ncor
      INTEGER   nout
      INTEGER   element (nface_len , elem_len)
      INTEGER   nbdv (nope_len, neta_len)
      INTEGER   nbvv (nbou_len, nvel_len)
      
      INTEGER   nvellnc(nbou_len)
      INTEGER   nbounc (nvel_len)
      INTEGER   ibtype(nbou_len)
      INTEGER   ibtypenc(nbou_len)
      INTEGER   nvdllnc (nope_len)
      INTEGER   nopenc (neta_len)
      
      INTEGER   nm (elem_len, nface_len)
      INTEGER   fileunit
      REAL*8    timenc(timenc_len)
      REAL(sz)  x(node_len)
      REAL(sz)  y(node_len)
      REAL(sz)  depth(node_len)
      REAL(sz)  u(node_len, timenc_len)
      REAL(sz)  v(node_len, timenc_len)
      REAL(sz), ALLOCATABLE, SAVE ::   minut(:)
      REAL(sz), ALLOCATABLE, SAVE ::   maxut(:)
      REAL(sz), ALLOCATABLE, SAVE ::   minvt(:)
      REAL(sz), ALLOCATABLE, SAVE ::   maxvt(:)
      REAL(sz), ALLOCATABLE, SAVE ::   minun(:)
      REAL(sz), ALLOCATABLE, SAVE ::   maxun(:)
      REAL(sz), ALLOCATABLE, SAVE ::   minvn(:)
      REAL(sz), ALLOCATABLE, SAVE ::   maxvn(:)

      INTEGER , ALLOCATABLE, SAVE ::  minunodet(:)
      INTEGER , ALLOCATABLE, SAVE ::  maxunodet(:)
      REAL(sz), ALLOCATABLE, SAVE ::  minutimen64(:)
      REAL(sz), ALLOCATABLE, SAVE ::  maxutimen64(:)
      
      INTEGER , ALLOCATABLE, SAVE ::  minvnodet(:)
      INTEGER , ALLOCATABLE, SAVE ::  maxvnodet(:)
      REAL(sz), ALLOCATABLE, SAVE ::  minvtimen64(:)
      REAL(sz), ALLOCATABLE, SAVE ::  maxvtimen64(:)

      REAL  minull
      REAL  maxull
      REAL  minvll
      REAL  maxvll
      real*8  mindepth
      real*8  maxdepth
      REAL*8 dtdp, statim,reftim, rnday, dramp, tau0
      REAL*8 a00, b00, c00,h0
      REAL(sz) slam0, sfea0, ffactor,esl, cori
      REAL(sz) slam0deg, sfea0deg
      CHARACTER(80) rundes, runid, agrid
      INTEGER, SAVE  :: FIRSTV = 1 
      INTEGER NDATA
      INTEGER, SAVE  :: NREC = 0
      INTEGER ncmode

!     -----------------
!     Attribute vectors
!     -----------------
      INTEGER intval(1),i,j
      INTEGER counti(1), starti(1), n
      INTEGER count(2), start(2)
      REAL    realval(1)
!     -----------------
!     Attribute general
!     -----------------
      CHARACTER(*) :: title, institution, source, history
      CHARACTER(*) :: references, comments, host, convention, contact
      CHARACTER(*) :: base_date

!     -----------------
!      date_string variables for time attribute
!     -----------------
      character date_string*40
      character now_date*8
      character big_ben*10
      character zone*5
      integer values(8)
!     ---------------------------------
!      inquiring and reading variables
!     ---------------------------------

      INTEGER NDIMS, NVARS, NGATTS, UNLIMDIMID, RECID,NRECS
      CHARACTER*4 recname
      
!     -------------------------
!     Save dims, ids, variables
!     -------------------------
      save node_dim,  elem_dim, nface_dim
      save nope_dim, neta_dim, nbou_dim, nvel_dim
      save u_id, v_id, X_id, Y_id
      save minut_id, maxut_id, minun_id, maxun_id
      save minvt_id, maxvt_id, minvn_id, maxvn_id

      save minunodet_id, maxunodet_id, minutimen_id, maxutimen_id
      save minvnodet_id, maxvnodet_id, minvtimen_id, maxvtimen_id

      save DEPTH_id,ELE_id, mindepth_id, maxdepth_id, nbdv_id, nbvv_id
      save timenc_dim, timenc_id
      save mindepth, maxdepth

      save nvelnc_dim
      save nvdll_id, nbounc_id, ibtypenc_id, nvellnc_id
      save netanc_dim
      save nopenc_id,nvdllnc_id

!     ---------------------------
!     Opening fort.64 netCDF file
!     ---------------------------
      WRITE(6,*) "-----------------------------------------------------"
       
      IF(FIRSTV.EQ.1.and.ncmode.NE.0) THEN
         
         ALLOCATE ( minun (node_len) )
         ALLOCATE ( maxun (node_len) )
         ALLOCATE ( minvn (node_len) )
         ALLOCATE ( maxvn (node_len) )
         
         ALLOCATE ( minutimen64 (node_len) )
         ALLOCATE ( maxutimen64 (node_len) )
         ALLOCATE ( minvtimen64 (node_len) )
         ALLOCATE ( maxvtimen64 (node_len) )

         do n=1,node_len
            minun (n) =  9999.d0
            maxun (n) = -9999.d0
            minvn (n) =  9999.d0
            maxvn (n) = -9999.d0
            minutimen64 (n) =  9999.d0
            maxutimen64 (n) = -9999.d0
            minvtimen64 (n) =  9999.d0
            maxvtimen64 (n) = -9999.d0
         end do
      ENDIF
       
!     --------------------------------
!     Switch order in array for NETCDF
!     --------------------------------
      do i=1, elem_len
         do j=1, nface_len
            element(j,i) = nm (i,j)
         end do
      end do
      
!     -----------------
!     Enter define mode
!     -----------------
         FILEFORT = 'fort.'
         FILEEXT = '.nc'
         WRITE(FILENAME,'(A5,I2,A3)')FILEFORT,FILEUNIT,FILEEXT
         WRITE(6,*) "Opening ",FILENAME,  " netCDF file ... "
         
       IF((FIRSTV.EQ.1.AND.IHOT.EQ.0.AND.NCMODE.NE.2).OR.
     &    (FIRSTV.EQ.1.AND.IHOT.NE.0.AND.NCMODE.NE.2.AND.
     &     nout.EQ.-3))THEN

           iret = nf_create(FILENAME, NF_CLOBBER, ncid)
           CALL check_err(iret)
           FIRSTV=0
         
!     --------------------------------
!     Convert to degress for NETCDF I/O
!     --------------------------------
           SLAM0DEG=SLAM0
           SFEA0DEG=SFEA0
           IF (NCMODE.EQ.1) THEN
               SLAM0DEG=SLAM0*RAD2DEG
               SFEA0DEG=SFEA0*RAD2DEG
           ENDIF
!     -----------------
!     Define dimensions
!     -----------------
      iret = nf_def_dim(ncid, 'node', node_len, node_dim)
      CALL check_err(iret)
      iret = nf_def_dim(ncid, 'nele', elem_len, elem_dim)
      CALL check_err(iret)
      iret = nf_def_dim(ncid, 'nface', 3, nface_dim)
      CALL check_err(iret)
      if(nope_len.ne.0) then
         iret = nf_def_dim(ncid, 'nope', nope_len, nope_dim)
         call check_err(iret)
      endif
      if(neta_len.ne.0) then
         iret = nf_def_dim(ncid, 'neta', neta_len, neta_dim)
         call check_err(iret)
      endif
      if(nbou_len.ne.0) then
         iret = nf_def_dim(ncid, 'nbou', nbou_len, nbou_dim)
         call check_err(iret)
      endif
      if(nvel_len.ne.0) then
         iret = nf_def_dim(ncid, 'nvel', nvel_len, nvel_dim)
         call check_err(iret)
      endif
      
      iret = nf_def_dim(ncid, 'time', nf_unlimited, timenc_dim)
      CALL check_err(iret)
      
!     ================
!     Define variables
!     ================
!     ----------------
!     Define time
!     ----------------
      timenc_dims(1) = timenc_dim
      iret = nf_def_var(ncid,'time',NF_DOUBLE, timenc_rank, timenc_dims, 
     &                  timenc_id)
      CALL check_err(iret)
      
!     -------------------
!     Define x
!     -------------------
      x_dims(1) = node_dim
      iret = nf_def_var(ncid, 'x', NF_DOUBLE, x_rank, x_dims, x_id)
      CALL check_err(iret)
!     -------------------
!     Define Y coordinate
!     --------------------
      Y_dims(1) = node_dim
      iret = nf_def_var(ncid, 'y', NF_DOUBLE, Y_rank, Y_dims, Y_id)
      CALL check_err(iret)
!     -------------------
!     Define elements
!     --------------------
      ELE_dims(1) = nface_dim
      ELE_dims(2) = elem_dim
      iret = nf_def_var(ncid, 'element',NF_INT,ELE_rank,ELE_dims,ELE_id)
      CALL check_err(iret)
!     -----------------------------------------------------------------
!     Define elevation specified boundary forcing segments information
!     -----------------------------------------------------------------
      if(nope_len.ne.0.and.neta_len.ne.0) then
        nvdllnc_dims(1) = nope_dim      
        iret = nf_def_var(ncid, 'nvdll',NF_INT, nvdllnc_rank, 
     &                  nvdllnc_dims, nvdllnc_id)
        CALL check_err(iret)

        nopenc_dims(1) = neta_dim      
        iret = nf_def_var(ncid, 'nbbv',NF_INT, nopenc_rank, nopenc_dims, 
     &                  nopenc_id)
        CALL check_err(iret)

      endif

!     ----------------------------------------------------------
!     Define normal flow specified boundary segments information
!     ----------------------------------------------------------
      if(nbou_len.ne.0.and.nvel_len.ne.0) then

        nvellnc_dims(1) = nbou_dim      
        iret = nf_def_var(ncid, 'nvell',NF_INT, nvellnc_rank, 
     &                  nvellnc_dims, nvellnc_id)
        CALL check_err(iret)
        
        ibtypenc_dims(1) = nbou_dim      
        iret = nf_def_var(ncid, 'ibtype',NF_INT, ibtypenc_rank, 
     &                  ibtypenc_dims, ibtypenc_id)
        CALL check_err(iret)
        nbounc_dims(1) = nvel_dim      
        iret = nf_def_var(ncid, 'nbvv',NF_INT, nbounc_rank, nbounc_dims, 
     &                  nbounc_id)
        CALL check_err(iret)
     
       endif

!     -------------------
!     Define Z coordinate
!     --------------------
      DEPTH_dims(1) = node_dim
      iret = nf_def_var(ncid, 'depth', NF_DOUBLE,DEPTH_rank, DEPTH_dims, 
     &                  DEPTH_id)
      CALL check_err(iret)
!     --------------------
!     Define minimum depth
!     --------------------
      iret = nf_def_var(ncid, 'mindepth', NF_DOUBLE, mindepth_rank, 0, 
     &                  mindepth_id)
      CALL check_err(iret)
!     --------------------
!     Define maximum depth
!     --------------------
      iret = nf_def_var(ncid, 'maxdepth', NF_DOUBLE, maxdepth_rank, 0, 
     &                  maxdepth_id)
      CALL check_err(iret)

!     ----------------
!     Define u vel
!     ----------------
      u_dims(1) = node_dim
      u_dims(2) = timenc_dim
      iret = nf_def_var(ncid, 'u-vel', NF_DOUBLE, u_rank, u_dims, 
     &                  u_id)
      CALL check_err(iret)

!     ----------------
!     Define v vel
!     ----------------
      v_dims(1) = node_dim
      v_dims(2) = timenc_dim
      iret = nf_def_var(ncid, 'v-vel', NF_DOUBLE, v_rank, v_dims, 
     &                  v_id)
      CALL check_err(iret)

!      WRITE(6,*)"Define MIN u/v at times"
!     -------------------------
!     Define minimum u at times
!     -------------------------
      minut_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'minuvelt', NF_DOUBLE, minut_rank, 
     &                  minut_dims, minut_id)
      CALL check_err(iret)
!     -------------------------
!     Define maximum u at times
!     -------------------------
      maxut_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'maxuvelt', NF_DOUBLE, maxut_rank, 
     &                  maxut_dims, maxut_id)
      CALL check_err(iret)

!     -------------------------
!     Define minimum v at times
!     -------------------------
      minvt_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'minvvelt', NF_DOUBLE, minvt_rank, 
     &                  minvt_dims, minvt_id)
      CALL check_err(iret)
      
!     -------------------------
!     Define maximum v at times
!     -------------------------
      maxvt_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'maxvvelt', NF_DOUBLE, maxvt_rank, 
     &                  maxvt_dims, maxvt_id)
      CALL check_err(iret)

!     -------------------------
!     Define minimum u at times
!     -------------------------
      minun_dims(1) = node_dim
      iret = nf_def_var(ncid, 'minuveln', NF_DOUBLE, minun_rank, 
     &                  minun_dims, minun_id)
      CALL check_err(iret)
!     -------------------------
!     Define maximum u at nodes
!     -------------------------
      maxun_dims(1) = node_dim
      iret = nf_def_var(ncid, 'maxuveln', NF_DOUBLE, maxun_rank, 
     &                  maxun_dims, maxun_id)
      CALL check_err(iret)

!     -------------------------
!     Define minimum v at nodes
!     -------------------------
      minvn_dims(1) = node_dim
      iret = nf_def_var(ncid, 'minvveln', NF_DOUBLE, minvn_rank, 
     &                  minvn_dims, minvn_id)
      CALL check_err(iret)
!     -------------------------
!     Define maximum v at nodes
!     -------------------------
      maxvn_dims(1) = node_dim
      iret = nf_def_var(ncid, 'maxvveln', NF_DOUBLE, maxvn_rank, 
     &                  maxvn_dims, maxvn_id)
      CALL check_err(iret)
!     ----------------------------
!     Define minimum u nodes at times
!     ----------------------------
      minunodet_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'minunodet', NF_INT, minunodet_rank, 
     &                  minunodet_dims, minunodet_id)
      CALL check_err(iret)
!     ----------------------------
!     Define maximum u nodes at times
!     ----------------------------
      maxunodet_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'maxunodet', NF_INT, maxunodet_rank, 
     &                  maxunodet_dims, maxunodet_id)
      CALL check_err(iret)
!     ----------------------------
!     Define minimum u nodes at nodes
!     ----------------------------
      minutimen_dims(1) = node_dim
      iret = nf_def_var(ncid, 'minutimen', NF_DOUBLE, minutimen_rank, 
     &                  minutimen_dims, minutimen_id)
      CALL check_err(iret)
!     ----------------------------
!     Define maximum u nodes at nodes
!     ----------------------------
      maxutimen_dims(1) = node_dim
      iret = nf_def_var(ncid, 'maxutimen', NF_DOUBLE, maxutimen_rank, 
     &                  maxutimen_dims, maxutimen_id)
      CALL check_err(iret)
            
!     ----------------------------
!     Define minimum v nodes at times
!     ----------------------------
      minvnodet_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'minvnodet', NF_INT, minvnodet_rank, 
     &                  minvnodet_dims, minvnodet_id)
      CALL check_err(iret)
!     ----------------------------
!     Define maximum v nodes at times
!     ----------------------------
      maxvnodet_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'maxvnodet', NF_INT, maxvnodet_rank, 
     &                  maxvnodet_dims, maxvnodet_id)
      CALL check_err(iret)
!     ----------------------------
!     Define minimum v nodes at nodes
!     ----------------------------
      minvtimen_dims(1) = node_dim
      iret = nf_def_var(ncid, 'minvtimen', NF_DOUBLE, minvtimen_rank, 
     &                  minvtimen_dims, minvtimen_id)
      CALL check_err(iret)
!     ----------------------------
!     Define maximum v nodes at nodes
!     ----------------------------
      maxvtimen_dims(1) = node_dim
      iret = nf_def_var(ncid, 'maxvtimen', NF_DOUBLE, maxvtimen_rank, 
     &                  maxvtimen_dims, maxvtimen_id)
      CALL check_err(iret)

!     -----------------
!     Assign attributes
!     -----------------

      if(ics.eq.2) then
      iret = nf_put_att_text(ncid, X_id, 'long_name', 9, 'longitude')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, X_id, 'standard_name', 9, 
     &                      'longitude')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, X_id, 'units', 12, VAR_DEGE)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, X_id, 'positive', 4, VAR_EAST)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, X_id, 'missing_values', NF_DOUBLE,
     &                        1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, X_id, '_FillValue', NF_DOUBLE, 1, 
     &                       doubleval)
      CALL check_err(iret)      
      
      iret = nf_put_att_text(ncid, Y_id, 'long_name', 8, 'latitude')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, Y_id, 'standard_name', 8, 
     &                       'latitude')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, Y_id, 'units', 13, VAR_DEGN)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, Y_id, 'positive', 5, VAR_NORTH)
      CALL check_err(iret)
       
      iret = nf_put_att_double(ncid, Y_id, 'missing_values', NF_DOUBLE, 
     &       1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, Y_id, '_FillValue', NF_DOUBLE, 1, 
     &       doubleval)
      CALL check_err(iret)

      else
      
      iret = nf_put_att_text(ncid, X_id, 'long_name', 22, 
     &       'Cartesian coordinate x')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, X_id, 'standard_name', 12, 
     &                      'x_coordinate')
      CALL check_err(iret)
      
      if (G.LT.11.D0) then
         iret = nf_put_att_text(ncid, X_id, 'units', 6,  VAR_METS)
         CALL check_err(iret)
      else
         iret = nf_put_att_text(ncid, X_id, 'units', 4,  VAR_FEET)
         CALL check_err(iret)
      endif
      
      iret = nf_put_att_text(ncid, X_id, 'positive', 5,  VAR_RG)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, X_id, 'missing_values', NF_DOUBLE, 
     &                       1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, X_id, '_FillValue', NF_DOUBLE, 1, 
     &                       doubleval)
      CALL check_err(iret)
!     ------------------------------
!     Define y-coordinate attributes
!     ------------------------------     
      
      iret = nf_put_att_text(ncid, Y_id, 'long_name', 22, 
     &       'Cartesian coordinate y')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, Y_id, 'standard_name', 12, 
     &                       'y_coordinate')
      CALL check_err(iret)
      if (G.LT.11.D0) then
         iret = nf_put_att_text(ncid, Y_id, 'units', 6, VAR_METS)
         CALL check_err(iret)
      else
         iret = nf_put_att_text(ncid, Y_id, 'units', 4, VAR_FEET)
         CALL check_err(iret)
      endif
      iret = nf_put_att_text(ncid, Y_id, 'positive', 2, VAR_UP)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, Y_id, 'missing_values', 
     &       NF_DOUBLE, 1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, Y_id, '_FillValue', NF_DOUBLE, 1, 
     &       doubleval)
      CALL check_err(iret)
      
      endif

!     ------------------------------
!     Define bathymetry attributes
!     ------------------------------

      iret = nf_put_att_text(ncid, DEPTH_id,'long_name',10,'bathymetry')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, DEPTH_id, 'standard_name', 5,'depth')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, DEPTH_id, 'units', 6, VAR_METS)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, DEPTH_id, 'missing_values', 
     &                      NF_DOUBLE,  1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid,DEPTH_id, '_FillValue',NF_DOUBLE, 1, 
     &                       doubleval)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, DEPTH_id, 'positive', 4, VAR_DW)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, ELE_id, 'long_name', 7, 'element')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, ELE_id, 'units', 14, VAR_NOND)
      intval(1) = -1
      
!    ---------------------------
!     minimum/maximum bathymetry
!    ---------------------------

      iret = nf_put_att_text(ncid, mindepth_id, 'long_name', 18, 
     &                       'minimum bathymetry')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, mindepth_id, 'standard_name', 13,
     &                'minimum_depth')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, mindepth_id, 'units', 6, VAR_METS)
      CALL check_err(iret)
      
      iret = nf_put_att_text(ncid, maxdepth_id, 'long_name', 18, 
     &                       'maximum bathymetry')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, maxdepth_id, 'standard_name', 13,
     &           'maximum_depth')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, maxdepth_id, 'units', 6, VAR_METS)
      CALL check_err(iret)
      
!    ----------------------------------------------------
!     minimum/maximum e/w and n/s velocities at each time
!    ----------------------------------------------------

!    --------------------------------------------------
!     determine variable names according to coordinates
!    ---------------------------------------------------
      if (ics.eq.2) then
         iret = nf_put_att_text(ncid, minut_id, 'long_name', 53, 
     &   'minimum vertically averaged e/w velocity at each time')
      else
          iret = nf_put_att_text(ncid, minut_id, 'long_name', 64, 
     &'minimum vertically averaged velocity in x-direction at each time
     &')
      endif
      CALL check_err(iret)
      
      iret = nf_put_att_text(ncid, minut_id, 'standard_name', 18,
     &                'minimum_u_velocity')
      CALL check_err(iret)
!    -------------------------
!     determine variable units
!    -------------------------
      if (G.LT.11.D0) then
          iret = nf_put_att_text(ncid, minut_id,'units',13, VAR_MPS)
      else
          iret = nf_put_att_text(ncid, minut_id,'units',11, VAR_FPS)
      endif
      CALL check_err(iret)
      
!    --------------------------------------------------
!     determine variable names according to coordinates
!    ---------------------------------------------------
      if (ics.eq.2) then
          iret = nf_put_att_text(ncid, maxut_id, 'long_name', 53, 
     &    'maximum vertically averaged e/w velocity at each time')
      else
          iret = nf_put_att_text(ncid, maxut_id, 'long_name', 64, 
     &     'maximum vertically averaged velocity in x-direction at each 
     &time')
      endif
      CALL check_err(iret)
       
      iret = nf_put_att_text(ncid, maxut_id, 'standard_name', 18,
     &                'maximum_u_velocity')
      CALL check_err(iret)
!    -------------------------
!     determine variable units
!    -------------------------
      if (G.LT.11.D0) then
          iret = nf_put_att_text(ncid, maxut_id,'units',13, VAR_MPS)
      else
          iret = nf_put_att_text(ncid, maxut_id,'units',11, VAR_FPS)
      endif
      CALL check_err(iret)
      
!    --------------------------------------------------
!     determine variable names according to coordinates
!    ---------------------------------------------------
      if (ics.eq.2) then
          iret = nf_put_att_text(ncid, minvt_id, 'long_name', 53, 
     &    'minimum vertically averaged n/s velocity at each time')
      else
          iret = nf_put_att_text(ncid, minvt_id, 'long_name', 64, 
     &'minimum vertically averaged velocity in y-direction at each time
     &')
      endif
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, minvt_id, 'standard_name', 18,
     &                'minimum_v_velocity')
      CALL check_err(iret)
!    -------------------------
!     determine variable units
!    -------------------------
      if (G.LT.11.D0) then
          iret = nf_put_att_text(ncid, minvt_id,'units',13, VAR_MPS)
      else
          iret = nf_put_att_text(ncid, minvt_id,'units',11, VAR_FPS)
      endif
      CALL check_err(iret)
      
!    --------------------------------------------------
!     determine variable names according to coordinates
!    ---------------------------------------------------
      if (ics.eq.2) then
          iret = nf_put_att_text(ncid, maxvt_id, 'long_name', 53, 
     &          'maximum vertically averaged n/s velocity at each time')
      else
            iret = nf_put_att_text(ncid, maxvt_id, 'long_name', 64, 
     &'maximum vertically averaged velocity in y-direction at each time
     &')
      endif
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, maxvt_id, 'standard_name', 18,
     &                'maximum_v_velocity')
      CALL check_err(iret)
!    -------------------------
!     determine variable units
!    -------------------------
      if (G.LT.11.D0) then
          iret = nf_put_att_text(ncid, maxvt_id,'units',13, VAR_MPS)
      else
          iret = nf_put_att_text(ncid, maxvt_id,'units',11, VAR_FPS)
      endif
      CALL check_err(iret)

!    ----------------------------------------------------
!     minimum/maximum e/w and n/s velocities at each node
!    ----------------------------------------------------
!    --------------------------------------------------
!     determine variable names according to coordinates
!    ---------------------------------------------------
      if (ics.eq.2) then
         iret = nf_put_att_text(ncid, minun_id, 'long_name', 53, 
     &    'minimum vertically averaged e/w velocity at each node')
      else
          iret = nf_put_att_text(ncid, minun_id, 'long_name', 64, 
     &'minimum vertically averaged velocity in x-direction at each node
     &')
      endif
      CALL check_err(iret)      
      iret = nf_put_att_text(ncid, minun_id, 'standard_name', 18,
     &                'minimum_u_velocity')
      CALL check_err(iret)
!    -------------------------
!     determine variable units
!    -------------------------
      if (G.LT.11.D0) then
          iret = nf_put_att_text(ncid, minun_id,'units',13, VAR_MPS)
      else
          iret = nf_put_att_text(ncid, minun_id,'units',11, VAR_FPS)
      endif
      CALL check_err(iret)
      
!    --------------------------------------------------
!     determine variable names according to coordinates
!    ---------------------------------------------------
      if (ics.eq.2) then
          iret = nf_put_att_text(ncid, maxun_id, 'long_name', 53, 
     &    'maximum vertically averaged e/w velocity at each node')
      else      
          iret = nf_put_att_text(ncid, maxun_id, 'long_name', 64, 
     &'maximum vertically averaged velocity in x-direction at each node
     &')
      endif
      CALL check_err(iret)
       
      iret = nf_put_att_text(ncid, maxun_id, 'standard_name', 18,
     &                'maximum_u_velocity')
      CALL check_err(iret)
!    -------------------------
!     determine variable units
!    -------------------------
      if (G.LT.11.D0) then
          iret = nf_put_att_text(ncid, maxun_id,'units',13, VAR_MPS)
      else
          iret = nf_put_att_text(ncid, maxun_id,'units',11, VAR_FPS)
      endif
      CALL check_err(iret)
      
!    --------------------------------------------------
!     determine variable name according to coordinates
!    ---------------------------------------------------
      if (ics.eq.2) then
          iret = nf_put_att_text(ncid, minvn_id, 'long_name', 53, 
     &    'minimum vertically averaged n/s velocity at each node')
      else
          iret = nf_put_att_text(ncid, minvn_id, 'long_name', 64, 
     &'minimum vertically averaged velocity in y-direction at each node
     &')
      endif
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, minvn_id, 'standard_name', 18,
     &                'minimum_v_velocity')
      CALL check_err(iret)
!    -------------------------
!     determine variable units
!    -------------------------
      if (G.LT.11.D0) then
         iret = nf_put_att_text(ncid, minvn_id,'units',13, VAR_MPS)
      else
         iret = nf_put_att_text(ncid, minvn_id,'units',11, VAR_FPS)
      endif
      CALL check_err(iret)
      
!    --------------------------------------------------
!     determine variable name according to coordinates
!    ---------------------------------------------------
      if (ics.eq.2) then
          iret = nf_put_att_text(ncid, maxvn_id, 'long_name', 53, 
     &          'maximum vertically averaged n/s velocity at each node')
      else      
          iret = nf_put_att_text(ncid, maxvn_id, 'long_name', 64, 
     &'maximum vertically averaged velocity in y-direction at each node
     &')
      endif
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, maxvn_id, 'standard_name', 18,
     &                'maximum_v_velocity')
      CALL check_err(iret)
!    -------------------------
!     determine variable units
!    -------------------------
      if (G.LT.11.D0) then
          iret = nf_put_att_text(ncid, maxvn_id,'units',13,  VAR_MPS)
      else
          iret = nf_put_att_text(ncid, maxvn_id,'units',11, VAR_FPS)
      endif
      CALL check_err(iret)

!    ----------------------------------------------------
!     minimum/maximum e/w and n/s velocities at each time
!    ----------------------------------------------------

!    --------------------------------------------------
!     determine variable names according to coordinates
!    ---------------------------------------------------

      if (ics.eq.2) then
         iret = nf_put_att_text(ncid, minunodet_id, 'long_name', 58, 
     &   'node where minimum e/w velocity occurs at each output time')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, minunodet_id, 'standard_name', 28,
     &                'node_of_minimum_e/w_velocity')
      CALL check_err(iret)
      else
          iret = nf_put_att_text(ncid, minunodet_id, 'long_name', 69, 
     &'node where minimum velocity in x-direction occurs at each output 
     &time')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, minunodet_id, 'standard_name', 26,
     &                'node_of_minimum_u_velocity')
      CALL check_err(iret)
      endif
      
      iret = nf_put_att_text(ncid, minunodet_id,'units',14, VAR_NOND)
      CALL check_err(iret)
      
!    --------------------------------------------------
!     determine variable name according to coordinates
!    ---------------------------------------------------
      if (ics.eq.2) then
         iret = nf_put_att_text(ncid, maxunodet_id, 'long_name', 58, 
     &    'node where maximum e/w velocity occurs at each output time')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, maxunodet_id, 'standard_name', 28,
     &                'node_of_maximum_e/w_velocity')
      CALL check_err(iret)
      else
         iret = nf_put_att_text(ncid, maxunodet_id, 'long_name', 69, 
     &'node where maximum velocity in x-direction occurs at each output 
     &time')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, maxunodet_id, 'standard_name', 26,
     &                'node_of_maximum_u_velocity')
      CALL check_err(iret)
      endif
       
      iret = nf_put_att_text(ncid, maxunodet_id,'units',14,  VAR_NOND)
      CALL check_err(iret)
      
!    --------------------------------------------------
!     determine variable name according to coordinates
!    ---------------------------------------------------
      if (ics.eq.2) then
         iret = nf_put_att_text(ncid, minvnodet_id, 'long_name', 58, 
     &   'node where minimum n/s velocity occurs at each output time')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, minvnodet_id, 'standard_name', 28,
     &                'node_of_minimum_n/s_velocity')
      CALL check_err(iret)
      else
          iret = nf_put_att_text(ncid, minvnodet_id, 'long_name', 69, 
     &'node where minimum velocity in y-direction occurs at each output 
     &time')
       CALL check_err(iret)
      iret = nf_put_att_text(ncid, minvnodet_id, 'standard_name', 26,
     &                'node_of_minimum_v_velocity')
      CALL check_err(iret)
      endif
      iret = nf_put_att_text(ncid, minvnodet_id,'units',14,  VAR_NOND)
      CALL check_err(iret)
!    --------------------------------------------------
!     determine variable name according to coordinates
!    ---------------------------------------------------
      if (ics.eq.2) then
         iret = nf_put_att_text(ncid, maxvnodet_id, 'long_name', 58, 
     &    'node where maximum n/s velocity occurs at each output time')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, maxvnodet_id, 'standard_name', 28,
     &                'node_of_maximum_n/s_velocity')
      CALL check_err(iret)
      else
         iret = nf_put_att_text(ncid, maxvnodet_id, 'long_name', 69, 
     &'node where maximum velocity in y-direction occurs at each output 
     &time')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, maxvnodet_id, 'standard_name', 26,
     &                'node_of_maximum_v_velocity')
      CALL check_err(iret)
      endif
      iret = nf_put_att_text(ncid, maxvnodet_id,'units',14, VAR_NOND)
      CALL check_err(iret)
!    ----------------------------------------------------
!     minimum/maximum e/w and n/s velocities at each node
!    ----------------------------------------------------
!    --------------------------------------------------
!     determine variable name according to coordinates
!    ---------------------------------------------------
      if (ics.eq.2) then
         iret = nf_put_att_text(ncid, minutimen_id, 'long_name', 50, 
     &   'time when minimum e/w velocity occurs at each node')
      CALL check_err(iret)

      
      iret = nf_put_att_text(ncid, minutimen_id, 'standard_name', 28,
     &                'time_of_minimum_e/w_velocity')
      CALL check_err(iret)
      else

         iret = nf_put_att_text(ncid, minutimen_id, 'long_name', 61, 
     &  'time when minimum velocity in x-direction occurs at each node')
      endif
      CALL check_err(iret)

      
      iret = nf_put_att_text(ncid, minutimen_id, 'standard_name', 26,
     &                'time_of_minimum_u_velocity')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, minutimen_id,'units',4, VAR_SECS)
      CALL check_err(iret)
      
!    --------------------------------------------------
!     determine variable name according to coordinates
!    ---------------------------------------------------
      if (ics.eq.2) then
          iret = nf_put_att_text(ncid, maxutimen_id, 'long_name', 50, 
     &    'time when maximum e/w velocity occurs at each node')
      CALL check_err(iret)
       
      iret = nf_put_att_text(ncid, maxutimen_id, 'standard_name', 28,
     &                'time_of_maximum_e/w_velocity')
      CALL check_err(iret)
      else

          iret = nf_put_att_text(ncid, maxutimen_id, 'long_name', 61, 
     &  'time when maximum velocity in x-direction occurs at each node')
      CALL check_err(iret)
       
      iret = nf_put_att_text(ncid, maxutimen_id, 'standard_name', 26,
     &                'time_of_maximum_u_velocity')
      CALL check_err(iret)
      endif
      iret = nf_put_att_text(ncid, maxutimen_id,'units',4, VAR_SECS)
      CALL check_err(iret)
      
!    --------------------------------------------------
!     determine variable name according to coordinates
!    ---------------------------------------------------
      if (ics.eq.2) then
         iret = nf_put_att_text(ncid, minvtimen_id, 'long_name', 50, 
     &   'time when minimum n/s velocity occurs at each node')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, minvtimen_id, 'standard_name', 28,
     &                'time_of_minimum_n/s_velocity')
      CALL check_err(iret)
      else

         iret = nf_put_att_text(ncid, minvtimen_id, 'long_name', 61, 
     &  'time when minimum velocity in y-direction occurs at each node')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, minvtimen_id, 'standard_name', 26,
     &                'time_of_minimum_v_velocity')
      CALL check_err(iret)
      endif
      iret = nf_put_att_text(ncid, minvtimen_id,'units',4, VAR_SECS)
      CALL check_err(iret)
      
!    --------------------------------------------------
!     determine variable name according to coordinates
!    ---------------------------------------------------
      if (ics.eq.2) then
          iret = nf_put_att_text(ncid, maxvtimen_id, 'long_name', 50, 
     &      'time when maximum n/s velocity occurs at each_node')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, maxvtimen_id, 'standard_name', 28,
     &                'time_of_maximum_n/s_velocity')
      CALL check_err(iret)
       else

         iret = nf_put_att_text(ncid, maxvtimen_id, 'long_name', 61, 
     &  'time when maximum velocity in y-direction occurs at each node')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, maxvtimen_id, 'standard_name', 26,
     &                'time_of_maximum_v_velocity')
      CALL check_err(iret)
      endif
      iret = nf_put_att_text(ncid, maxvtimen_id,'units',4, VAR_SECS)
      CALL check_err(iret)

!     -----------------------------------------------------------------------
!     Define node numbers on elevation specified boundary segments attributes
!     -----------------------------------------------------------------------      
      iret = nf_put_att_text(ncid, nopenc_id, 'long_name', 53, 
     &       'node numbers on elevation specified boundary segments')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, nopenc_id, 'units', 14, 
     &                      VAR_NOND)
      CALL check_err(iret)
      
      iret = nf_put_att_text(ncid, nvdllnc_id, 'long_name', 60, 
     &'number of nodes in each elevation specified boundary segment')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, nvdllnc_id, 'units', 14, 
     &                      VAR_NOND)
      CALL check_err(iret)
!     ---------------------------------------------------------------
!     Define node numbers on normal flow boundary segments attributes
!     ----------------------------------------------------------------  

      iret = nf_put_att_text(ncid, nbounc_id, 'long_name', 59, 
     & 'node numbers in each normal flow specified boundary segment')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, nbounc_id, 'units', 14, 
     &                      VAR_NOND)
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, ibtypenc_id, 'long_name', 13, 
     &       'boundary type')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, ibtypenc_id, 'units', 14, 
     &                      VAR_NOND)
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, nvellnc_id, 'long_name', 62, 
     &'number of nodes in each normal flow specified boundary segment')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, nvellnc_id, 'units', 14, 
     &                      VAR_NOND)
      CALL check_err(iret)

!    -----------------------
!     Define time attributes
!    -----------------------
      iret = nf_put_att_text(ncid,timenc_id,'long_name',10,'model time')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, timenc_id, 'standard_name', 4,'time')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, timenc_id, 'units', 23,VAR_TUN)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, timenc_id, 'base_date',
     . (LEN(TRIM(base_date(1:INDEX(base_date,"!")-1)))), base_date)
      CALL check_err(iret)

!    -----------
!     Define u velocity
!    -----------

!    --------------------------------------------------
!     determine variable name according to coordinates
!    ---------------------------------------------------
      if (ics.eq.2) then
          iret = nf_put_att_text(ncid, u_id, 'long_name', 32, 
     &                       'vertically averaged e/w velocity')
       else
         iret = nf_put_att_text(ncid, u_id, 'long_name', 43, 
     &            'vertically averaged velocity in x-direction')
      endif
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, u_id, 'standard_name', 10, 
     &                      'u_velocity')
      CALL check_err(iret)
!    -------------------------
!     determine variable units
!    -------------------------
      if (G.LT.11.D0) then
          iret = nf_put_att_text(ncid, u_id, 'units',13, VAR_MPS)
      else
          iret = nf_put_att_text(ncid, u_id, 'units',11, VAR_FPS)
      endif
      CALL check_err(iret)
!    --------------------------------------------------
!     determine variable name according to coordinates
!    ---------------------------------------------------
      if (ics.eq.2) then
          iret = nf_put_att_text(ncid, u_id, 'positive', 4, VAR_EAST)
      else
          iret = nf_put_att_text(ncid, u_id, 'positive', 5, VAR_RG)
      endif
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, u_id, 'missing_values', 
     &       NF_DOUBLE, 1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, u_id, '_FillValue', NF_DOUBLE, 1, 
     &                       doubleval)
      CALL check_err(iret)
      
      iret = nf_put_att_double(ncid, u_id, 'dry_Value', NF_DOUBLE, 1, 
     &                       doubleval)
      CALL check_err(iret)
!    -----------
!     Define v velocity
!    -----------
!    --------------------------------------------------
!     determine variable name according to coordinates
!    ---------------------------------------------------
      if (ics.eq.2) then
          iret = nf_put_att_text(ncid, v_id, 'long_name', 32, 
     &           'vertically averaged n/s velocity')
      else
          iret = nf_put_att_text(ncid, v_id, 'long_name', 43, 
     &           'vertically averaged velocity in y-direction')
      endif
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, v_id, 'standard_name', 10, 
     &                      'v_velocity')
      CALL check_err(iret)
!    -------------------------
!     determine variable units
!    -------------------------
      if (G.LT.11.D0) then
         iret = nf_put_att_text(ncid, v_id, 'units',13, VAR_MPS)
      else
         iret = nf_put_att_text(ncid, v_id, 'units',11, VAR_FPS)
      endif
      CALL check_err(iret)

!    --------------------------------------------------
!     determine variable name according to coordinates
!    ---------------------------------------------------
      if (ics.eq.2) then
          iret = nf_put_att_text(ncid, v_id, 'positive', 5, VAR_NORTH)
      else
          iret = nf_put_att_text(ncid, v_id, 'positive', 2, VAR_UP)
      endif
      CALL check_err(iret)

      iret = nf_put_att_DOUBLE(ncid, v_id, 'missing_values', 
     &       NF_DOUBLE, 1, doubleval)
      CALL check_err(iret)
!      doubleval(1) = -99999.
      iret = nf_put_att_DOUBLE(ncid, v_id, '_FillValue', NF_DOUBLE, 1, 
     &                       doubleval)
      CALL check_err(iret)
!      doubleval(1) = -99999.
      iret = nf_put_att_DOUBLE(ncid, v_id, 'dry_Value', NF_DOUBLE, 1, 
     &                       doubleval)
      CALL check_err(iret)

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
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'description', 
     &       LEN(TRIM(rundes(1:INDEX(rundes,"!")-1))), rundes)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'grid',
     &       LEN(TRIM(agrid)), agrid)       
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'title',
     &       LEN(TRIM(title(1:INDEX(title,"!")-1))), 
     &       title)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'institution', 
     &       LEN(TRIM(institution(1:INDEX(institution,"!")-1))), 
     &       institution)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'source', 
     &       LEN(TRIM(source(1:INDEX(source,"!")-1))),
     &       source)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'history', 
     &       LEN(TRIM(history(1:INDEX(history,"!")-1))),
     &       history)
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
     &       LEN(TRIM(host(1:INDEX(host,"!")-1))), 
     &       host)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'convention', 
     &       LEN(TRIM(convention(1:INDEX(convention,"!")-1))), 
     &       convention)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'contact', 
     &       LEN(TRIM(contact(1:INDEX(contact,"!")-1))), 
     &       contact)
      CALL check_err(iret)
 
 
      CALL date_and_Time(now_date,big_ben,zone,values)
      write(date_string,71) values(1),values(2),values(3)
     & ,values(5),values(6),values(7),   (values(4))/60
 71   format(I4,'-',I2.2,'-',i2.2,' ',i2,':',i2.2,':',i2.2,' '
     &  ,i3.2,':00')
      iret = nf_put_att_text(ncid, NF_GLOBAL,'creation_date'
     & ,len_trim(date_string),date_string)
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
      iret = nf_put_att_int(ncid, NCGLOBAL, 'irtype', NF_INT, 1, irtype)
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
      iret = nf_put_att_double(ncid, NCGLOBAL, 'ffactor', NF_DOUBLE, 1, 
     &                         ffactor)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, NCGLOBAL, 'esl', NF_DOUBLE, 1, 
     &                         esl)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, NCGLOBAL, 'cori', NF_DOUBLE, 1, 
     &                         cori)
      CALL check_err(iret)
      iret = nf_put_att_int(ncid, NCGLOBAL, 'ntif', NF_INT, 1, ntif)
      CALL check_err(iret)
      iret = nf_put_att_int(ncid, NCGLOBAL, 'nbfr', NF_INT, 1, nbfr)
      CALL check_err(iret)
     
!     ======================
!     Leave define mode
!     ======================
      iret = nf_enddef(ncid)
      CALL check_err(iret)
!     ======================
!     Store static variables
!     ======================
!     -------
!     Store X
!     -------
      iret = nf_put_var_double(ncid, X_id, X)
      CALL check_err(iret)
!     -------
!     Store Y
!     -------
      iret = nf_put_var_double(ncid, Y_id, Y)
      CALL check_err(iret)
!     ------------
!     Store depth
!     ------------
      iret = nf_put_var_double(ncid, DEPTH_id, depth)
      CALL check_err(iret)
!     --------------
!     Store elements
!     --------------
      count(1)=nface_len
      count(2)=elem_len
      start(1)=1
      start(2)=1
      iret=nf_put_vara_int(ncid,ele_id,start,count,element)
      CALL check_err(iret)
!     ------------------------------------
!     Store elevation boundary information
!     -------------------------------------
      iret = nf_put_var_int(ncid, nopenc_id, nopenc)
      CALL check_err(iret)
      iret = nf_put_var_int(ncid, nvdllnc_id, nvdllnc)
      CALL check_err(iret)
!     --------------------------------------
!     Store normal flow boundary information
!     --------------------------------------
      iret = nf_put_var_int(ncid, nbounc_id, nbounc)
      CALL check_err(iret)
      iret = nf_put_var_int(ncid, ibtypenc_id, ibtypenc)
      CALL check_err(iret)
      iret = nf_put_var_int(ncid, nvellnc_id, nvellnc)
      CALL check_err(iret)      
!     --------------------
!     Store maximum depth
!     --------------------
      maxdepth = maxval (depth)
      iret = nf_put_var_double(ncid, maxdepth_id, maxdepth)
      CALL check_err(iret)
!     --------------------
!     Store minimum depth
!     --------------------
      mindepth = minval (depth)
      iret = nf_put_var_double(ncid, mindepth_id, mindepth)
      CALL check_err(iret)
      
!     =====================
!     Open old netcdf file
!     =====================
      ELSE
          iret = NF_OPEN(FILENAME, NF_WRITE, ncid)
          CALL check_err(iret)   
!     ------------------
!     Inquire variables
!     ------------------
          iret=nf_inq_unlimdim(ncid,recid)
          CALL check_err(iret)   
          iret=nf_inq_dimlen(ncid,recid,nrec)
          CALL check_err(iret)   
      ENDIF
!     ------------------
!     Update Global Attributes
!     ------------------
      IF(IHOT.NE.0.AND.FIRSTV.EQ.1) THEN
         iret = NF_REDEF (ncid)
         CALL check_err(iret)
         iret = nf_put_att_int(ncid, NCGLOBAL, 'ihot', NF_INT, 1, ihot)
         CALL check_err(iret)
         iret = nf_put_att_int(ncid, NCGLOBAL, 'nramp', NF_INT, 1,nramp)
         CALL check_err(iret)
         iret = nf_put_att_double(ncid, NCGLOBAL, 'tau0', NF_DOUBLE, 1, 
     &                         tau0)
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
         iret = nf_put_att_double(ncid, NCGLOBAL, 'ffactor', NF_DOUBLE, 1, 
     &                         ffactor)
         CALL check_err(iret)
         iret = nf_put_att_double(ncid, NCGLOBAL, 'esl', NF_DOUBLE, 1, 
     &                         esl)
         CALL check_err(iret)
         CALL date_and_time(now_date,big_ben,zone,values)
         WRITE(date_string,71) values(1),values(2),values(3),
     &         values(5),values(6),values(7),   (values(4))/60
         iret = nf_put_att_text(ncid, NF_GLOBAL,'modification_date',
     &          len_trim(date_string),date_string)
         CALL check_err(iret)
         iret = NF_ENDDEF (ncid)
         CALL check_err(iret)
       ENDIF

      IF(ncmode.eq.0) go to 9999
      IF(FIRSTV.EQ.1.and.
     & ((ncmode.ne.0.and.IHOT.NE.0).or.(ncmode.eq.2.and.IHOT.eq.0).or.
     &  (ncmode.eq.2.and.IHOT.ne.0.and.nout.eq.-3)))THEN
         iret=nf_inq_varid(ncid,"time",timenc_id)
         iret=nf_inq_varid(ncid,"u-vel",u_id)
         iret=nf_inq_varid(ncid,"v-vel",v_id)
         iret=nf_inq_varid(ncid,"minuvelt",minut_id)
         iret=nf_inq_varid(ncid,"maxuvelt",maxut_id)
         iret=nf_inq_varid(ncid,"minvvelt",minvt_id)
         iret=nf_inq_varid(ncid,"maxvvelt",maxvt_id)
         iret=nf_inq_varid(ncid,"minuveln",minun_id)
         iret=nf_inq_varid(ncid,"maxuveln",maxun_id)
         iret=nf_inq_varid(ncid,"minvveln",minvn_id)
         iret=nf_inq_varid(ncid,"maxvveln",maxvn_id)
         iret=nf_inq_varid(ncid,"minunodet",minunodet_id)
         iret=nf_inq_varid(ncid,"maxunodet",maxunodet_id)
         iret=nf_inq_varid(ncid,"minvnodet",minvnodet_id)
         iret=nf_inq_varid(ncid,"maxvnodet",maxvnodet_id)
         iret=nf_inq_varid(ncid,"minutimen",minutimen_id)
         iret=nf_inq_varid(ncid,"maxutimen",maxutimen_id)
         iret=nf_inq_varid(ncid,"minvtimen",minvtimen_id)
         iret=nf_inq_varid(ncid,"maxvtimen",maxvtimen_id)
      ENDIF
      IF(IHOT.NE.0.AND.FIRSTV.EQ.1.and.ncmode.ne.0.and.nout.eq.3)THEN
         iret = nf_get_var_double(ncid, minun_id, minun)
         CALL check_err(iret)
         iret = nf_get_var_double(ncid, maxun_id, maxun)
         CALL check_err(iret)
         iret = nf_get_var_double(ncid, minvn_id, minvn)
         CALL check_err(iret)
         iret = nf_get_var_double(ncid, maxvn_id, maxvn)
         CALL check_err(iret)
         
         iret = nf_get_var_double(ncid, minutimen_id, minutimen64)
         CALL check_err(iret)
         iret = nf_get_var_double(ncid, maxutimen_id, maxutimen64)
         CALL check_err(iret)
         iret = nf_get_var_double(ncid, minvtimen_id, minvtimen64)
         CALL check_err(iret)
         iret = nf_get_var_double(ncid, maxvtimen_id, maxvtimen64)
         CALL check_err(iret)
      ENDIF
      
      FIRSTV=0
      NREC = NREC + 1
!     ----------
!     Store time
!     ----------
      starti(1)=nrec
      counti(1)=timenc_len
      iret = nf_put_vara_DOUBLE(ncid, timenc_id, starti,counti, timenc)
      CALL check_err(iret)
!     ---------------------------
!     Store u
!     ---------------------------
      count(1)=node_len
      count(2)=timenc_len
      start(1)=1
      start(2)=nrec
      iret = nf_put_vara_DOUBLE(ncid, u_id, start, count, u)
      CALL check_err(iret)
!     ---------------------------
!     Store v
!     ---------------------------
      count(1)=node_len
      count(2)=timenc_len
      start(1)=1
      start(2)=nrec
      iret = nf_put_vara_DOUBLE(ncid, v_id, start, count, v)
      CALL check_err(iret)
      
!     ----------------------------------------
!     Allocating memory for min and max arrays
!     ----------------------------------------
      IF(.NOT. ALLOCATED(minut))ALLOCATE( minut(timenc_len))
      IF(.NOT. ALLOCATED(maxut))ALLOCATE( maxut(timenc_len))
      IF(.NOT. ALLOCATED(minun))ALLOCATE( minun(node_len  ))
      IF(.NOT. ALLOCATED(maxun))ALLOCATE( maxun(node_len  ))

      IF(.NOT. ALLOCATED(minvt))ALLOCATE( minvt(timenc_len))
      IF(.NOT. ALLOCATED(maxvt))ALLOCATE( maxvt(timenc_len))
      IF(.NOT. ALLOCATED(minvn))ALLOCATE( minvn(node_len  ))
      IF(.NOT. ALLOCATED(maxvn))ALLOCATE( maxvn(node_len  ))

      IF(.NOT. ALLOCATED(minunodet  ))ALLOCATE( minunodet(timenc_len))
      IF(.NOT. ALLOCATED(maxunodet  ))ALLOCATE( maxunodet(timenc_len))
      IF(.NOT. ALLOCATED(minutimen64))ALLOCATE( minutimen64(node_len))
      IF(.NOT. ALLOCATED(maxutimen64))ALLOCATE( maxutimen64(node_len))

      IF(.NOT. ALLOCATED(minvnodet  ))ALLOCATE( minvnodet(timenc_len))
      IF(.NOT. ALLOCATED(maxvnodet  ))ALLOCATE( maxvnodet(timenc_len))
      IF(.NOT. ALLOCATED(minvtimen64))ALLOCATE( minvtimen64(node_len))
      IF(.NOT. ALLOCATED(maxvtimen64))ALLOCATE( maxvtimen64(node_len))

!     ---------------------------
!     Calculating u minimum
!     ---------------------------
!     ----------------------------------
!     Calculating u minimum in each node
!     ----------------------------------
      do n=1,node_len
         if(u (n, 1).lt.minun(n) ) then
 	    minun (n)       = u (n,1)
 	    minutimen64 (n) = timenc(1)
 	 endif
      end do

!     ----------------------------------
!     Calculating u minimum for each time
!     ----------------------------------
      minut (1) = 99999.
      do n=1,node_len
 	 if(u (n,1).gt.-99999.and.u(n,1).lt. minut (1) ) then
            minut (1)    = u (n,1)
            minunodet(1) = n
         endif
      end do

      iret = nf_put_vara_DOUBLE(ncid, minut_id,(/nrec/),(/1/),minut)
      CALL check_err(iret)
      iret=nf_put_vara_INT(ncid, minunodet_id,(/nrec/),(/1/),minunodet)
      CALL check_err(iret)
!     ---------------------
!     Calculating u maximum
!     ---------------------
!     ----------------------------------
!     Calculating u maximum in each node
!     ----------------------------------
      do n=1,node_len
 	 if(u (n,1).gt.maxun(n) ) then
            maxun (n)       = u (n,1)
            maxutimen64 (n) = timenc(1)
 	 endif
      end do
!     -----------------------------------
!     Calculating u maximum for each time
!     -----------------------------------
      maxut (1) = -9999.
      do n=1,node_len
 	 if(u(n,1).gt.maxut(1))  then
            maxut (1)    = u (n,1)
            maxunodet(1) = n 		    
 	 endif
      end do
      iret = nf_put_vara_DOUBLE(ncid,maxut_id,(/nrec/),(/1/),maxut)
      CALL check_err(iret)
      iret=nf_put_vara_INT(ncid,maxunodet_id,(/nrec/),(/1/),maxunodet)
      CALL check_err(iret)

!     ---------------------------
!     Calculating v minimum
!     ---------------------------
!     ----------------------------------
!     Calculating v minimum in each node
!     ----------------------------------
      do n=1,node_len
 	 if(v (n,1).lt.minvn(n) ) then
            minvn (n)       = v (n,1)
            minvtimen64 (n) = timenc(1)
 	 endif
      end do
!     -----------------------------------
!     Calculating v minimum for each time
!     -----------------------------------
      minvt (1) = 99999.
      do n=1,node_len
         if(v (n,1).gt.-99999.and.v(n,1).lt. minvt (1) ) then
            minvt (1)    = v (n,1)
            minvnodet(1) = n
         endif
      end do

      iret = nf_put_vara_DOUBLE(ncid, minvt_id,(/nrec/),(/1/),minvt)
      CALL check_err(iret)
      iret=nf_put_vara_INT(ncid, minvnodet_id,(/nrec/),(/1/),minvnodet)
      CALL check_err(iret)
!     ---------------------------
!     Calculating v maximum
!     ---------------------------
!     ----------------------------------
!     Calculating v maximum in each node
!     ----------------------------------
      do n=1,node_len
         if(v (n,1).gt.maxvn(n) ) then
            maxvn (n)      = v (n,1)
            maxvtimen64(n) = timenc(1)
         endif
      end do
!     -----------------------------------
!     Calculating v maximum for each time
!     -----------------------------------
      maxvt (1) = -9999.
      do n=1,node_len
         if(v(n,1).gt.maxvt(1)) then
            maxvt (1) = v (n,1)
            maxvnodet(1) = n 		    
         endif
      end do
      iret = nf_put_vara_DOUBLE(ncid,maxvt_id,(/nrec/),(/1/),maxvt)
      CALL check_err(iret)
      iret=nf_put_vara_INT(ncid,maxvnodet_id,(/nrec/),(/1/),maxvnodet)
      CALL check_err(iret)

!     ---------------------------
!     Closing fort.X4 netCDF file
!     ---------------------------
      IF(NREC.EQ.NDATA) THEN
         iret = nf_put_var_double(ncid, minun_id, minun)
         CALL check_err(iret)
         iret = nf_put_var_double(ncid, maxun_id, maxun)
         CALL check_err(iret)
         iret = nf_put_var_double(ncid, minvn_id, minvn)
         CALL check_err(iret)
         iret = nf_put_var_double(ncid, maxvn_id, maxvn)
         CALL check_err(iret)
         
         iret = nf_put_var_double(ncid, minutimen_id, minutimen64)
         CALL check_err(iret)
         iret = nf_put_var_double(ncid, maxutimen_id, maxutimen64)
         CALL check_err(iret)
         iret = nf_put_var_double(ncid, minvtimen_id, minvtimen64)
         CALL check_err(iret)
         iret = nf_put_var_double(ncid, maxvtimen_id, maxvtimen64)
         CALL check_err(iret)
      ENDIF

!     ---------------------------
!     Closing netCDF file
!     ---------------------------
      WRITE(6,*) "Closing ", FILENAME, " netCDF file ... "
      WRITE(6,*) "-----------------------------------------------------"
      
9999  iret = nf_close(ncid)
      CALL check_err(iret)
      WRITE(6,*)'end subroutine writenc64'

      end subroutine writenc64     

!=======================================================================
!     Subroutine writenc73
!     Writes a netCDF file equivalent to ADCIRC's fort.73 ascii file
!
!     Dimensions:
!     node  : number of nodes in the horizontal grid
!     nele  : number of elements in the horizontal grid
!     nope  : number of elevation specified boundary forcing segments
!     neta  : total number of elevation specified boundary nodes
!     nbou  : number of normal flow specified boundary segments
!     nvel  : total number of normal flow specified boundary nodes
!     time  : number of time records
!
!     Variables :
!     time      : time
!     zeta      : pressure 
!     x         : x-coordinate or longitude
!     y         : y-coordinate or latitude
!     element   : elements in grid
!     nvdll     : number of nodes in elevation boundary segment
!     nbdv      : node numbers on elevation specified boundary segment
!     nvell     : number of nodes in normal flow specified boundary segment
!     nbvv      : node numbers on normal flow boundary segment
!     ibtype    : boundary type
!     depth     : depth in each node
!     mindepth  : minimum depth in grid
!     maxdepth  : maximum depth in grid
!     minprest  : minimum pressure for each output time
!     maxprest  : maximum pressure for each output time
!     minpresn  : minimum pressure at each node
!     maxpresn  : maximum pressure at each node
!     minpnodet : node of minimum pressure at each output time
!     maxpnodet : node of maximum pressure at each output time
!     minptimen : time of minimum pressure at each node
!     maxptimen : time of maximum pressure at each node
!
!     Attributes:
!     model             : ADCIRC
!     version           : model version number
!     grid_type         : triangular
!     rundes            : description of simulation
!     agrid             : grid ID 
!     title             : description of the dataset
!     institution       : where the original data was produced
!     source            : method of production of the original data
!     history           : provides audit trail for modifications to the original data
!     references        : published/web-based references that describe data
!     comments          : miscellaneous information about the data
!     contact           : user's institution/ e-mail address
!     creation_date     : date the data was created
!     modification_date : date the data was modified
!     Parameters from fort.15 input file 
!
!     Revision history:
!
!     Date            Programmer                  Description of change
!     ----            ----------                  ----------------------
!     03/30/07  Cristina Forbes, PSGS @ UNC-IMS   Wrote original code
!     03/30/08  Cristina Forbes, PSGS @ UNC-IMS   Modified code for 
!                                                 globalio & hotstart 
!                                                 from binary files
!=======================================================================
      subroutine writenc73 (nout,fileunit,timenc_len, node_len,timenc,
     &                      zeta,nt, NDATA, rundes,runid, agrid, irtype,
     &                      dtdp, ihot, ics, nolibf,  nolifa,
     &                      nolica, nolicat, nwp,
     &                      ncor, ntip, nws, nramp, tau0,
     &                      statim, reftim, rnday, dramp,
     &                      a00, b00, c00, h0,
     &                      slam0, sfea0, ffactor,esl, cori,
     &                      ntif, nbfr, elem_len, nface_len, nm,
     &                      x,y,depth,
     &                      nope_len, neta_len, nbdv, nopenc,netanc_len,
     &                      nvdllnc,
     &                      nbou_len, nvel_len, nbvv,nbounc,nvelnc_len,
     &                      nvellnc,ibtypenc,
     &                      ncmode,
     &                      title, institution, source, history,
     &                      references, comments, host, convention,
     &                      contact,
     &                      base_date)

!--------------------------------------------
!     Force explicit declaration of all variables
!--------------------------------------------
      IMPLICIT NONE
!     -------------------
!     netCDF header file
!     -------------------
!      include 'netcdf.inc'
!     -------------------
!     Error status return
!     -------------------
      INTEGER  iret
!     ---------
!     netCDF id
!     ---------
      INTEGER  ncid
!     -------------
!     Dimension ids
!     -------------
      INTEGER  node_dim
      INTEGER  timenc_dim
      INTEGER  elem_dim
      INTEGER  nface_dim
      INTEGER  nope_dim
      INTEGER  neta_dim
      INTEGER  nbou_dim
      INTEGER  nvel_dim
      INTEGER  slen_dim
      INTEGER  scode_dim
      INTEGER  nvelnc_dim
      INTEGER  netanc_dim
!     -----------------
!     Dimension lengths
!     -----------------
      INTEGER  node_len
      INTEGER  timenc_len
      INTEGER  elem_len
      INTEGER  nface_len
      INTEGER  nope_len
      INTEGER  neta_len
      INTEGER  nbou_len
      INTEGER  nvel_len
      INTEGER  nvelnc_len
      INTEGER  netanc_len
!     ------------
!     Variable ids
!     ------------
      INTEGER  timenc_id
      INTEGER  zeta_id
      INTEGER  minzetat_id
      INTEGER  maxzetat_id
      INTEGER  minzetan_id
      INTEGER  maxzetan_id
      INTEGER  minznodet_id
      INTEGER  maxznodet_id
      INTEGER  minztimen_id
      INTEGER  maxztimen_id
      INTEGER  station_id
      INTEGER  scode_id
      INTEGER  X_id
      INTEGER  Y_id
      INTEGER  DEPTH_id
      INTEGER  ELE_id
      INTEGER  mindepth_id
      INTEGER  maxdepth_id
      INTEGER  nbdv_id
      INTEGER  nbvv_id

      INTEGER  nvdll_id
      INTEGER  nbounc_id
      INTEGER  ibtypenc_id
      INTEGER  nvellnc_id
      INTEGER  nvdllnc_id
      INTEGER  nopenc_id
!     ---------------------------------------------
!     Rank (number of dimensions) for each variable
!     ---------------------------------------------
      INTEGER  timenc_rank
      INTEGER  zeta_rank
      INTEGER  station_rank
      INTEGER  scode_rank
      INTEGER  minzetat_rank
      INTEGER  maxzetat_rank
      INTEGER  minzetan_rank
      INTEGER  maxzetan_rank
      INTEGER  minznodet_rank
      INTEGER  maxznodet_rank
      INTEGER  minztimen_rank
      INTEGER  maxztimen_rank
      INTEGER  depth_rank
      INTEGER  ele_rank
      INTEGER  mindepth_rank
      INTEGER  maxdepth_rank
      INTEGER  x_rank
      INTEGER  y_rank

      INTEGER  nbounc_rank
      INTEGER  nvdll_rank
      INTEGER  ibtypenc_rank
      INTEGER  nvellnc_rank
      INTEGER  nvdllnc_rank
      INTEGER  nopenc_rank
      
      INTEGER  nbdv_rank
      INTEGER  nbvv_rank
      PARAMETER (timenc_rank = 1)
      PARAMETER (station_rank = 2)
      PARAMETER (scode_rank = 2)
      PARAMETER (zeta_rank = 2)
      PARAMETER (minzetat_rank = 1)
      PARAMETER (maxzetat_rank = 1)
      PARAMETER (minzetan_rank = 1)
      PARAMETER (maxzetan_rank = 1)
      PARAMETER (minznodet_rank = 1)
      PARAMETER (maxznodet_rank = 1)
      PARAMETER (minztimen_rank = 1)
      PARAMETER (maxztimen_rank = 1)
      PARAMETER (depth_rank = 1)
      PARAMETER (ele_rank = 2)
      PARAMETER (mindepth_rank = 0)
      PARAMETER (maxdepth_rank = 0)
      PARAMETER (x_rank = 1)
      PARAMETER (y_rank = 1)

      PARAMETER (nvdll_rank = 1)
      PARAMETER (nbounc_rank = 1)
      PARAMETER (ibtypenc_rank = 1)
      PARAMETER (nvellnc_rank = 1)

      PARAMETER (nvdllnc_rank=1)
      PARAMETER (nopenc_rank=1)

      PARAMETER (nbdv_rank = 2)
      PARAMETER (nbvv_rank = 2)
!     ---------------
!     Variable shapes
!     ---------------
      INTEGER  timenc_dims(timenc_rank)
      INTEGER  zeta_dims(zeta_rank)
      INTEGER  station_dims  (station_rank)
      INTEGER  scode_dims    (scode_rank)
      INTEGER  minzetat_dims(minzetat_rank)
      INTEGER  maxzetat_dims(maxzetat_rank)
      INTEGER  minzetan_dims(minzetan_rank)
      INTEGER  maxzetan_dims(maxzetan_rank)
      INTEGER  minznodet_dims(minznodet_rank)
      INTEGER  maxznodet_dims(maxznodet_rank)
      INTEGER  minztimen_dims(minztimen_rank)
      INTEGER  maxztimen_dims(maxztimen_rank)
      INTEGER  x_dims(x_rank)
      INTEGER  y_dims(y_rank)
      INTEGER  depth_dims(depth_rank)
      INTEGER  ele_dims(ele_rank)

      INTEGER  nvdll_dims(nvdll_rank)
      INTEGER  nbounc_dims(nbounc_rank)
      INTEGER  ibtypenc_dims(ibtypenc_rank)
      INTEGER  nvellnc_dims(nvellnc_rank)
      INTEGER  nvdllnc_dims(nvdllnc_rank)
      INTEGER  nopenc_dims(nopenc_rank)

      INTEGER  nbdv_dims(nbdv_rank)
      INTEGER  nbvv_dims(nbvv_rank)
!     --------------
!     Data variables
!     --------------
      CHARACTER*5 FILEFORT
      CHARACTER*10 FILENAME
      CHARACTER*3 FILEEXT
      INTEGER   t,ntif,nbfr,nt
      INTEGER   irtype, ihot, ics, nolibf, nolifa, nolica, nolicat, nwp
      INTEGER   ntip, nws, nramp,ncor
      INTEGER   nvel, neta, nout
      INTEGER   element (nface_len , elem_len)
      INTEGER   nbdv (nope_len, neta_len)
      INTEGER   nbvv (nbou_len, nvel_len)

      INTEGER   nvellnc(nbou_len)
      INTEGER   nbounc (nvel_len)
      INTEGER   ibtype(nbou_len)
      INTEGER   ibtypenc(nbou_len)
      INTEGER   nvdllnc (nope_len)
      INTEGER   nopenc (neta_len)

      INTEGER   nm (elem_len, nface_len)
      INTEGER   fileunit
      REAL*8    timenc(timenc_len)
      REAL(sz)  x(node_len)
      REAL(sz)  y(node_len)
      REAL(sz)  depth(node_len)
      REAL(sz)  zeta (node_len, timenc_len)
      REAL(sz), ALLOCATABLE, SAVE ::  minzetat(:)
      REAL(sz), ALLOCATABLE, SAVE ::  maxzetat(:)
      REAL(sz), ALLOCATABLE, SAVE ::  minzetan(:)
      REAL(sz), ALLOCATABLE, SAVE ::  maxzetan(:)
      INTEGER , ALLOCATABLE, SAVE ::  minznodet(:)
      INTEGER , ALLOCATABLE, SAVE ::  maxznodet(:)
      REAL(sz), ALLOCATABLE, SAVE ::  minztimen73(:)
      REAL(sz), ALLOCATABLE, SAVE ::  maxztimen73(:)
      real*8  mindepth
      real*8  maxdepth
      REAL*8  dtdp, statim,reftim, rnday, dramp, tau0
      REAL*8  a00, b00, c00,h0
      REAL(sz) slam0, sfea0
      REAL(sz) slam0deg, sfea0deg
      REAL*8 ffactor,esl, cori
      CHARACTER(80) rundes, runid, agrid
      INTEGER, SAVE  :: FIRST73 = 1 
      INTEGER NDATA
      INTEGER, SAVE  :: NREC = 0
      INTEGER ncmode
!     -----------------
!     Attribute vectors
!     -----------------
      INTEGER intval(1),i,j
      INTEGER counti(1), starti(1), n
      INTEGER count(2), start(2)
      REAL    realval(1)
!      DOUBLE PRECISION  doubleval(1)
!     -----------------
!     Attribute general
!     -----------------
      CHARACTER(*) :: title, institution, source, history
      CHARACTER(*) :: references, comments, host, convention, contact
      CHARACTER(*) :: base_date
!     -----------------
!      date_string variables for time attribute
!     -----------------
      character date_string*40
      character now_date*8
      character big_ben*10
      character zone*5
      integer   values(8)
      real*8    xmin, xmax

!     ---------------------------------
!      inquiring and reading variables
!     ---------------------------------

      INTEGER NDIMS, NVARS, NGATTS, UNLIMDIMID, RECID,NRECS
      CHARACTER*4 recname

!     -------------------------
!     Save dims, ids, variables
!     -------------------------
      save node_dim, elem_dim, nface_dim
      save nope_dim, neta_dim, nbou_dim, nvel_dim
      save zeta_id, X_id, Y_id
      save minzetat_id, maxzetat_id, minzetan_id, maxzetan_id
      save minznodet_id, maxznodet_id, minztimen_id, maxztimen_id
      save DEPTH_id,ELE_id, mindepth_id, maxdepth_id, nbdv_id, nbvv_id
      save timenc_dim, timenc_id
      save station_id, scode_id
      save mindepth, maxdepth

      save nvelnc_dim
      save nvdll_id, nbounc_id, ibtypenc_id, nvellnc_id
      save netanc_dim
      save nopenc_id,nvdllnc_id

!     ---------------------------
!     Opening fort.63 netCDF file
!     ---------------------------
      WRITE(6,*) "-----------------------------------------------------"
!      NREC = NREC + 1
      
      IF(FIRST73.EQ.1.and.ncmode.NE.0) THEN
         
         ALLOCATE ( minzetan    (node_len) )
         ALLOCATE ( maxzetan    (node_len) )
         ALLOCATE ( minztimen73 (node_len) )
         ALLOCATE ( maxztimen73 (node_len) )
        
         do n=1,node_len
            minzetan    (n) =  9999.d0
            maxzetan    (n) = -9999.d0
            minztimen73 (n) =  9999.d0
            maxztimen73 (n) = -9999.d0
         end do
      ENDIF

!     --------------------------------
!     Switch order in array for NETCDF
!     --------------------------------
      do i=1, elem_len
         do j=1, nface_len
            element(j,i) = nm (i,j)
         end do
      end do
!     -----------------
!     Enter define mode
!     -----------------
         FILEFORT = 'fort.'
         FILEEXT = '.nc'
         WRITE(FILENAME,'(A5,I2,A3)')FILEFORT,FILEUNIT,FILEEXT
         WRITE(6,*) "Opening ",FILENAME,  " netCDF file ... "
                  
       IF((FIRST73.EQ.1.AND.IHOT.EQ.0.AND.NCMODE.NE.2).OR.
     &    (FIRST73.EQ.1.AND.IHOT.NE.0.AND.NCMODE.NE.2.AND.
     &     nout.EQ.-3))THEN

           iret = nf_create(FILENAME, NF_CLOBBER, ncid)
           CALL check_err(iret)
           FIRST73=0
         
!     --------------------------------
!     Convert to degress for NETCDF I/O
!     --------------------------------
           SLAM0DEG=SLAM0
           SFEA0DEG=SFEA0
           IF (NCMODE.EQ.1) THEN
               SLAM0DEG=SLAM0*RAD2DEG
               SFEA0DEG=SFEA0*RAD2DEG
           ENDIF

!     =================
!     Define dimensions
!     =================
      iret = nf_def_dim(ncid, 'node', node_len, node_dim)
      CALL check_err(iret)
      iret = nf_def_dim(ncid, 'nele', elem_len, elem_dim)
      CALL check_err(iret)
      iret = nf_def_dim(ncid, 'nface', 3, nface_dim)
      CALL check_err(iret)
      if(nope_len.ne.0) then
         iret = nf_def_dim(ncid, 'nope', nope_len, nope_dim)
         call check_err(iret)
      endif
      if(neta_len.ne.0) then
         iret = nf_def_dim(ncid, 'neta', neta_len, neta_dim)
         call check_err(iret)
      endif
      if(nbou_len.ne.0) then
         iret = nf_def_dim(ncid, 'nbou', nbou_len, nbou_dim)
         call check_err(iret)
      endif
      if(nvel_len.ne.0) then
         iret = nf_def_dim(ncid, 'nvel', nvel_len, nvel_dim)
         call check_err(iret)
      endif
      iret = nf_def_dim(ncid, 'time', nf_unlimited, timenc_dim)
      CALL check_err(iret)
!     ================
!     Define variables
!     ================
!     ----------------
!     Define time
!     ----------------

      timenc_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'time',NF_DOUBLE,timenc_rank, timenc_dims, 
     &                  timenc_id)
      CALL check_err(iret)
     
!     -------------------
!     Define X
!     -------------------
      X_dims(1) = node_dim
      iret = nf_def_var(ncid, 'x', NF_DOUBLE, X_rank, X_dims, X_id)
      CALL check_err(iret)
!     -------------------
!     Define Y coordinate
!     --------------------
      Y_dims(1) = node_dim
      iret = nf_def_var(ncid, 'y', NF_DOUBLE, Y_rank, Y_dims, Y_id)
      CALL check_err(iret)
!     -------------------
!     Define elements
!     --------------------
      ELE_dims(1) = nface_dim
      ELE_dims(2) = elem_dim
      iret = nf_def_var(ncid, 'element',NF_INT,ELE_rank,ELE_dims,ELE_id)
      CALL check_err(iret)
      
!     -----------------------------------------------------------------
!     Define elevation specified boundary forcing segments information
!     -----------------------------------------------------------------
      if(nope_len.ne.0.and.neta_len.ne.0) then
        nvdllnc_dims(1) = nope_dim      
        iret = nf_def_var(ncid, 'nvdll',NF_INT, nvdllnc_rank, 
     &                  nvdllnc_dims, nvdllnc_id)
        CALL check_err(iret)

        nopenc_dims(1) = neta_dim      
        iret = nf_def_var(ncid, 'nbdv',NF_INT, nopenc_rank, nopenc_dims, 
     &                  nopenc_id)
        CALL check_err(iret)

      endif
!     ----------------------------------------------------------
!     Define normal flow specified boundary segments information
!     ----------------------------------------------------------
      if(nbou_len.ne.0.and.nvel_len.ne.0) then

        nvellnc_dims(1) = nbou_dim      
        iret = nf_def_var(ncid, 'nvell',NF_INT, nvellnc_rank, 
     &                  nvellnc_dims, nvellnc_id)
        CALL check_err(iret)
        
        ibtypenc_dims(1) = nbou_dim      
        iret = nf_def_var(ncid, 'ibtype',NF_INT, ibtypenc_rank, 
     &                  ibtypenc_dims, ibtypenc_id)
        CALL check_err(iret)
     
         nbounc_dims(1) = nvel_dim      
        iret = nf_def_var(ncid, 'nbvv',NF_INT, nbounc_rank, nbounc_dims, 
     &                  nbounc_id)
        CALL check_err(iret)
      endif
!     -------------------
!     Define Z coordinate
!     --------------------
      DEPTH_dims(1) = node_dim
      iret = nf_def_var(ncid, 'depth',NF_DOUBLE, DEPTH_rank, DEPTH_dims, 
     &                  DEPTH_id)
      CALL check_err(iret)
!     --------------------
!     Define minimum depth
!     --------------------
      iret = nf_def_var(ncid, 'mindepth',NF_DOUBLE, mindepth_rank, 0, 
     &                  mindepth_id)
      CALL check_err(iret)
!     --------------------
!     Define maximum depth
!     --------------------
      iret = nf_def_var(ncid, 'maxdepth', NF_DOUBLE, maxdepth_rank, 0, 
     &                  maxdepth_id)
      CALL check_err(iret)
      
!     ----------------
!     Define zeta
!     ----------------
      zeta_dims(1) = node_dim
      zeta_dims(2) = timenc_dim     
      iret = nf_def_var(ncid,'pressure',NF_DOUBLE, zeta_rank, zeta_dims, 
     &                  zeta_id)            
      CALL check_err(iret)
!     ----------------------------
!     Define minimum zeta at times
!     ----------------------------
      minzetat_dims(1) = timenc_dim

      iret = nf_def_var(ncid, 'minprest', NF_DOUBLE, minzetat_rank, 
     &                  minzetat_dims, minzetat_id)
      CALL check_err(iret)
!     ----------------------------
!     Define maximum zeta at times
!     ----------------------------
      maxzetat_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'maxprest', NF_DOUBLE, maxzetat_rank, 
     &                  maxzetat_dims, maxzetat_id)
      CALL check_err(iret)
!     ----------------------------
!     Define minimum zeta at nodes
!     ----------------------------
      minzetan_dims(1) = node_dim
       iret = nf_def_var(ncid, 'minpresn', NF_DOUBLE, minzetan_rank, 
     &                  minzetan_dims, minzetan_id)
      CALL check_err(iret)
!     ----------------------------
!     Define maximum zeta at nodes
!     ----------------------------
      maxzetan_dims(1) = node_dim
      iret = nf_def_var(ncid, 'maxpresn', NF_DOUBLE, maxzetan_rank, 
     &                  maxzetan_dims, maxzetan_id)
      CALL check_err(iret)
      
!     ----------------------------
!     Define minimum zeta nodes at times
!     ----------------------------
      minznodet_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'minpnodet', NF_INT, minznodet_rank, 
     &                  minznodet_dims, minznodet_id)
      CALL check_err(iret)
!     ----------------------------
!     Define maximum zeta nodes at times
!     ----------------------------
      maxznodet_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'maxpnodet', NF_INT, maxznodet_rank, 
     &                  maxznodet_dims, maxznodet_id)
      CALL check_err(iret)
!     ----------------------------
!     Define minimum zeta nodes at nodes
!     ----------------------------
      minztimen_dims(1) = node_dim
      iret = nf_def_var(ncid, 'minptimen', NF_DOUBLE, minztimen_rank, 
     &                  minztimen_dims, minztimen_id)
      CALL check_err(iret)
!     ----------------------------
!     Define maximum zeta nodes at nodes
!     ----------------------------
      maxztimen_dims(1) = node_dim
      iret = nf_def_var(ncid, 'maxptimen', NF_DOUBLE, maxztimen_rank, 
     &                  maxztimen_dims, maxztimen_id)
      CALL check_err(iret)
!     -----------------
!     Assign attributes
!     -----------------      
!     ----------------------------
!     Define longitude attributes
!     ----------------------------
      if(ics.eq.2) then
      
      	iret = nf_put_att_text(ncid,X_id,'long_name',9,'longitude')
      	       CALL check_err(iret)
      	iret = nf_put_att_text(ncid, X_id, 'standard_name', 9, 
     &                      'longitude')
      	       CALL check_err(iret)
      	iret = nf_put_att_text(ncid, X_id, 'units', 12, VAR_DEGE)
      	       CALL check_err(iret)
      	iret = nf_put_att_text(ncid, X_id, 'positive', 4, VAR_EAST)
      	       CALL check_err(iret)

      	iret = nf_put_att_double(ncid, X_id, 'missing_values',NF_DOUBLE, 
     &                       1, doubleval)
      	       CALL check_err(iret)

      	iret = nf_put_att_double(ncid, X_id, '_FillValue', NF_DOUBLE, 1, 
     &                       doubleval)
      	       CALL check_err(iret)            
!     ----------------------------
!     Define latitude attributes
!     ----------------------------
      	iret = nf_put_att_text(ncid, Y_id, 'long_name', 8, 'latitude')
      	       CALL check_err(iret)
      	iret = nf_put_att_text(ncid, Y_id, 'standard_name', 8, 
     &                       'latitude')
      	       CALL check_err(iret)
      	iret = nf_put_att_text(ncid, Y_id, 'units', 13, VAR_DEGN)
      	       CALL check_err(iret)
      	iret = nf_put_att_text(ncid, Y_id, 'positive', 5, VAR_NORTH)
      	       CALL check_err(iret)
       
 
      	iret = nf_put_att_double(ncid, Y_id, 'missing_values', 
     &       NF_DOUBLE, 1, doubleval)
      	       CALL check_err(iret)
 
      	iret = nf_put_att_double(ncid, Y_id, '_FillValue', NF_DOUBLE, 1, 
     &      doubleval)
      	       CALL check_err(iret)
       else
      
!     ------------------------------
!     Define x-coordinate attributes
!     ------------------------------
      iret = nf_put_att_text(ncid, X_id, 'long_name', 22, 
     &                      'Cartesian coordinate x')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, X_id, 'standard_name',12, 
     &                       'x-coordinate')
      CALL check_err(iret)
      if (G.LT.11.D0) then
          iret = nf_put_att_text(ncid, X_id, 'units', 6,  VAR_METS)
      	  CALL check_err(iret)
      else
          iret = nf_put_att_text(ncid, X_id, 'units', 4,  VAR_FEET)
      	  CALL check_err(iret)
      endif
      iret = nf_put_att_text(ncid, X_id, 'positive', 5,  VAR_RG)
      CALL check_err(iret)

      iret = nf_put_att_double(ncid, X_id, 'missing_values',NF_DOUBLE,
     &                         1, doubleval)
      CALL check_err(iret)

      iret = nf_put_att_double(ncid, X_id, '_FillValue', NF_DOUBLE, 1, 
     &                        doubleval)
      CALL check_err(iret)
!     ------------------------------
!     Define y-coordinate attributes
!     ------------------------------     
      iret = nf_put_att_text(ncid, Y_id, 'long_name', 22, 
     &                        'Cartesian coordinate y')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, Y_id, 'standard_name', 12, 
     &                       'y-coordinate')
      CALL check_err(iret)
      if (G.LT.11.D0) then
	  iret = nf_put_att_text(ncid, Y_id, 'units', 6, VAR_METS)
      else
     	  iret = nf_put_att_text(ncid, Y_id, 'units', 4, VAR_FEET)
      endif
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, Y_id, 'positive', 2, VAR_UP)
      CALL check_err(iret)
 
      iret = nf_put_att_double(ncid, Y_id, 'missing_values',NF_DOUBLE, 
     &                         1, doubleval)
      	CALL check_err(iret)
 
      	iret = nf_put_att_double(ncid, Y_id, '_FillValue', NF_DOUBLE, 1, 
     &         doubleval)
      	CALL check_err(iret)      
      endif
!     ------------------------------
!     Define bathymetry attributes
!     ------------------------------
      iret = nf_put_att_text(ncid, DEPTH_id,'long_name',10,'bathymetry')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, DEPTH_id, 'standard_name', 5,'depth')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, DEPTH_id, 'units', 6, VAR_METS)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, DEPTH_id, 'missing_values', 
     &                      NF_DOUBLE, 1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, DEPTH_id, '_FillValue', NF_DOUBLE,
     &                        1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, DEPTH_id, 'positive', 4, VAR_DW)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, ELE_id, 'long_name', 7, 'element')
      CALL check_err(iret)
      intval(1) = -1
      iret = nf_put_att_text(ncid, ELE_id, 'units', 14, 
     &                      VAR_NOND)
      CALL check_err(iret)
      
!     ------------------------------------
!     Define minimum bathymetry attributes
!     ------------------------------------
      iret = nf_put_att_text(ncid, mindepth_id, 'long_name', 18, 
     &                       'minimum bathymetry')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, mindepth_id, 'standard_name', 13,
     &                'minimum_epth')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, mindepth_id, 'units', 6, VAR_METS)
      CALL check_err(iret)
      
!     ------------------------------------
!     Define maximum bathymetry attributes
!     ------------------------------------
      iret = nf_put_att_text(ncid, maxdepth_id, 'long_name', 18, 
     &                       'maximum bathymetry')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, maxdepth_id, 'standard_name', 13,
     &                'maximum_depth')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, maxdepth_id, 'units', 6, VAR_METS)
      CALL check_err(iret)
      
!     -------------------------------------------------------------
!     Define max and min in time sea surface pressure attributes 
!     -------------------------------------------------------------
      
      iret = nf_put_att_text(ncid, minzetat_id, 'long_name', 50, 
     &      'minimum sea surface pressure for each output time')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, minzetat_id, 'units', 15, VAR73_UN)
      CALL check_err(iret)
      
      iret = nf_put_att_text(ncid, maxzetat_id, 'long_name', 50, 
     &       'maximum sea surface pressure for each output time')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, maxzetat_id, 'units',15,VAR73_UN) 
       CALL check_err(iret)
      
!     -----------------------------------------------
!     Define maximum sea surface pressure attributes
!     -----------------------------------------------
      iret = nf_put_att_text(ncid, minzetan_id, 'long_name', 42, 
     &      'minimum sea surface pressure at each node')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, minzetan_id, 'units', 15, VAR73_UN)
      CALL check_err(iret)
      
      iret = nf_put_att_text(ncid, maxzetan_id, 'long_name', 42, 
     &       'maximum sea surface pressure at each node')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, maxzetan_id, 'units', 15, VAR73_UN)
      CALL check_err(iret)

!     -------------------------------------------------------------
!     Define max and min in time sea surface elevation attributes 
!     -------------------------------------------------------------
      iret = nf_put_att_text(ncid, minznodet_id, 'long_name', 64, 
     &     'node of minimum atmospheric surface pressure at each output
     & time')
             CALL check_err(iret)
      iret = nf_put_att_text(ncid, minznodet_id, 'units',14, VAR_NOND)
             CALL check_err(iret)
      
      iret = nf_put_att_text(ncid, maxznodet_id, 'long_name', 64, 
     &     'node of maximum atmospheric surface pressure at each output
     & time')
             CALL check_err(iret)
      iret = nf_put_att_text(ncid, maxznodet_id, 'units', 14, VAR_NOND)
             CALL check_err(iret)
!     -------------------------------------------------------------
!     Define max and min in nodes sea surface elevation attributes
!     -------------------------------------------------------------
      iret = nf_put_att_text(ncid, minztimen_id, 'long_name', 57, 
     &     'time of minimum atmospheric surface pressure at each node')
            CALL check_err(iret)
      iret = nf_put_att_text(ncid, minztimen_id, 'units', 4, VAR_SECS)
            CALL check_err(iret)
      
      iret = nf_put_att_text(ncid, maxztimen_id, 'long_name', 58, 
     &     'time of maximum atmospheric surface pressure at each node')
            CALL check_err(iret)
      iret = nf_put_att_text(ncid, maxztimen_id, 'units', 4, VAR_SECS)
            CALL check_err(iret)

!     -----------------------------------------------------------------------
!     Define node numbers on elevation specified boundary segments attributes
!     -----------------------------------------------------------------------      
      iret = nf_put_att_text(ncid, nopenc_id, 'long_name', 53, 
     &       'node numbers on elevation specified boundary segments')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, nopenc_id, 'units', 14, 
     &                      VAR_NOND)
      CALL check_err(iret)
      
      iret = nf_put_att_text(ncid, nvdllnc_id, 'long_name', 60, 
     &'number of nodes in each elevation specified boundary segment')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, nvdllnc_id, 'units', 14, 
     &                      VAR_NOND)
      CALL check_err(iret)
!     ---------------------------------------------------------------
!     Define node numbers on normal flow boundary segments attributes
!     ----------------------------------------------------------------  

      iret = nf_put_att_text(ncid, nbounc_id, 'long_name', 59, 
     & 'node numbers in each normal flow specified boundary segment')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, nbounc_id, 'units', 14, 
     &                      VAR_NOND)
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, ibtypenc_id, 'long_name', 13, 
     &       'boundary type')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, ibtypenc_id, 'units', 14, 
     &                      VAR_NOND)
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, nvellnc_id, 'long_name', 62, 
     &'number of nodes in each normal flow specified boundary segment')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, nvellnc_id, 'units', 14, 
     &                      VAR_NOND)
      CALL check_err(iret)

!    -----------------------
!     Define time attributes
!    -----------------------
      iret = nf_put_att_text(ncid,timenc_id,'long_name',10,'model time')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid,timenc_id, 'standard_name', 4, 'time')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, timenc_id, 'units', 23, 
     .       'seconds since base_date')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, timenc_id, 'base_date',
     .       (LEN(TRIM(base_date(1:INDEX(base_date,"!")-1)))), base_date)
      CALL check_err(iret)

!    --------------------------------------------
!     Define sea surface pressure attributes
!    --------------------------------------------

      iret = nf_put_att_text(ncid, zeta_id,'long_name',21, VAR73_LN)
     &                   
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, zeta_id, 'standard_name', 17, 
     &                      VAR73_ST)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, zeta_id, 'units', 15, VAR73_UN)
      CALL check_err(iret)
      iret = nf_put_att_DOUBLE(ncid,zeta_id,'missing_values',NF_DOUBLE,
     &1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_DOUBLE(ncid, zeta_id, '_FillValue', NF_DOUBLE, 
     &                      1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, zeta_id,'positive', 2, VAR_UP)
      CALL check_err(iret)
      
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
      intval(1) = 1
      iret = nf_put_att_int(ncid, NCGLOBAL, 'irtype', NF_INT, 1, irtype)
      CALL check_err(iret)
      intval(1) = 0
      iret = nf_put_att_int(ncid, NCGLOBAL, 'ihot', NF_INT, 1, ihot)
      CALL check_err(iret)
      intval(1) = 1
      iret = nf_put_att_int(ncid, NCGLOBAL, 'ics', NF_INT, 1, ics)
      CALL check_err(iret)
      intval(1) = 1
      iret = nf_put_att_int(ncid, NCGLOBAL, 'nolibf', NF_INT, 1, nolibf)
      CALL check_err(iret)
      intval(1) = 1
      iret = nf_put_att_int(ncid, NCGLOBAL, 'nolifa', NF_INT, 1, nolifa)
      CALL check_err(iret)
      intval(1) = 1
      iret = nf_put_att_int(ncid, NCGLOBAL, 'nolica', NF_INT, 1, nolica)
      CALL check_err(iret)
      intval(1) = 1
      iret = nf_put_att_int(ncid, NCGLOBAL, 'nolicat', NF_INT, 1, 
     &                      nolicat)
      CALL check_err(iret)
      intval(1) = 0
      iret = nf_put_att_int(ncid, NCGLOBAL, 'nwp', NF_INT, 1, nwp)
      CALL check_err(iret)
      intval(1) = 0
      iret = nf_put_att_int(ncid, NCGLOBAL, 'ncor', NF_INT, 1, ncor)
      CALL check_err(iret)
      intval(1) = 0
      iret = nf_put_att_int(ncid, NCGLOBAL, 'ntip', NF_INT, 1, ntip)
      CALL check_err(iret)
      intval(1) = 0
      iret = nf_put_att_int(ncid, NCGLOBAL, 'nws', NF_INT, 1, nws)
      CALL check_err(iret)
      intval(1) = 1
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
      iret = nf_put_att_double(ncid, NCGLOBAL, 'ffactor', NF_DOUBLE, 1, 
     &                         ffactor)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, NCGLOBAL, 'esl', NF_DOUBLE, 1, 
     &                         esl)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, NCGLOBAL, 'cori', NF_DOUBLE, 1, 
     &                         cori)
      CALL check_err(iret)
      intval(1) = 0
      iret = nf_put_att_int(ncid, NCGLOBAL, 'ntif', NF_INT, 1, ntif)
      CALL check_err(iret)
      intval(1) = 1
      iret = nf_put_att_int(ncid, NCGLOBAL, 'nbfr', NF_INT, 1, nbfr)
      CALL check_err(iret)
      
!     ======================
!     Leave define mode
!     ======================
      iret = nf_enddef(ncid)
      CALL check_err(iret)
!     ======================
!     Store static variables
!     ======================
!     -------
!     Store X
!     -------
      iret = nf_put_var_double(ncid, X_id, X)
      CALL check_err(iret)
!     -------
!     Store Y
!     -------
      iret = nf_put_var_double(ncid, Y_id, Y)
      CALL check_err(iret)
!     ------------
!     Store depth
!     ------------
      iret = nf_put_var_double(ncid, DEPTH_id, depth)
      CALL check_err(iret)
!     --------------
!     Store elements
!     --------------
      count(1)=nface_len
      count(2)=elem_len
      start(1)=1
      start(2)=1
      iret=nf_put_vara_int(ncid,ele_id,start,count,element)
      CALL check_err(iret)
!     ------------------------------------
!     Store elevation boundary information
!     -------------------------------------
      iret = nf_put_var_int(ncid, nopenc_id, nopenc)
      CALL check_err(iret)
      iret = nf_put_var_int(ncid, nvdllnc_id, nvdllnc)
      CALL check_err(iret)
!     --------------------------------------
!     Store normal flow boundary information
!     --------------------------------------
      iret = nf_put_var_int(ncid, nbounc_id, nbounc)
      CALL check_err(iret)
      iret = nf_put_var_int(ncid, ibtypenc_id, ibtypenc)
      CALL check_err(iret)
      iret = nf_put_var_int(ncid, nvellnc_id, nvellnc)
      CALL check_err(iret)
!     --------------------
!     Store maximum depth
!     --------------------
      maxdepth = maxval (depth)
      iret = nf_put_var_double(ncid, maxdepth_id, maxdepth)
      CALL check_err(iret)
!     --------------------
!     Store minimum depth
!     --------------------
      mindepth = minval (depth)
      iret = nf_put_var_double(ncid, mindepth_id, mindepth)
      CALL check_err(iret)

!     =====================
!     Open old netcdf file
!     =====================
      ELSE
          iret = NF_OPEN(FILENAME, NF_WRITE, ncid)
          CALL check_err(iret)   
!     ------------------
!     Inquire variables
!     ------------------
          iret=nf_inq_unlimdim(ncid,recid)
          CALL check_err(iret)   
          iret=nf_inq_dimlen(ncid,recid,nrec)
          CALL check_err(iret)   
      ENDIF
!     ------------------
!     Update Global Attributes
!     ------------------
      IF(IHOT.NE.0.AND.FIRST73.EQ.1) THEN
         iret = NF_REDEF (ncid)
         CALL check_err(iret)
         iret = nf_put_att_int(ncid, NCGLOBAL, 'ihot', NF_INT, 1, ihot)
         CALL check_err(iret)
         iret = nf_put_att_int(ncid, NCGLOBAL, 'nramp', NF_INT, 1,nramp)
         CALL check_err(iret)
         iret = nf_put_att_double(ncid, NCGLOBAL, 'tau0', NF_DOUBLE, 1, 
     &                         tau0)
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
         iret = nf_put_att_double(ncid, NCGLOBAL, 'ffactor', NF_DOUBLE, 1, 
     &                         ffactor)
         CALL check_err(iret)
         iret = nf_put_att_double(ncid, NCGLOBAL, 'esl', NF_DOUBLE, 1, 
     &                         esl)
         CALL check_err(iret)
         CALL date_and_time(now_date,big_ben,zone,values)
         WRITE(date_string,71) values(1),values(2),values(3),
     &         values(5),values(6),values(7),   (values(4))/60
         iret = nf_put_att_text(ncid, NF_GLOBAL,'modification_date',
     &          len_trim(date_string),date_string)
         CALL check_err(iret)
         iret = NF_ENDDEF (ncid)
         CALL check_err(iret)
!        WRITE(6,*)"finish global att"
      ENDIF  

      IF(ncmode.eq.0) go to 9999
      IF(FIRST73.EQ.1.and.
     & ((ncmode.ne.0.and.IHOT.NE.0).or.(ncmode.eq.2.and.IHOT.eq.0).or.
     &  (ncmode.eq.2.and.IHOT.ne.0.and.nout.eq.-3)))THEN
         iret=nf_inq_varid(ncid,"time"     ,timenc_id)
         iret=nf_inq_varid(ncid,"pressure" ,zeta_id)
         iret=nf_inq_varid(ncid,"minprest" ,minzetat_id)
         iret=nf_inq_varid(ncid,"maxprest" ,maxzetat_id)
         iret=nf_inq_varid(ncid,"minpresn" ,minzetan_id)
         iret=nf_inq_varid(ncid,"maxpresn" ,maxzetan_id)
         iret=nf_inq_varid(ncid,"minpnodet",minznodet_id)
         iret=nf_inq_varid(ncid,"maxpnodet",maxznodet_id)
         iret=nf_inq_varid(ncid,"minptimen",minztimen_id)
         iret=nf_inq_varid(ncid,"maxptimen",maxztimen_id)
      ENDIF
      IF(IHOT.NE.0.AND.FIRST73.EQ.1.and.ncmode.ne.0.and.nout.eq.3)THEN
         iret = nf_get_var_double(ncid, minzetan_id, minzetan)
         CALL check_err(iret)
         iret = nf_get_var_double(ncid, maxzetan_id, maxzetan)
         CALL check_err(iret)
         iret = nf_get_var_double(ncid, minztimen_id, minztimen73)
         CALL check_err(iret)
         iret = nf_get_var_double(ncid, maxztimen_id, maxztimen73)
         CALL check_err(iret)
      ENDIF
      
      FIRST73=0
      NREC = NREC + 1
      
!     ----------
!     Store time
!     ----------
      starti(1)=nrec
      counti(1)=timenc_len
      iret = nf_put_vara_DOUBLE(ncid, timenc_id, starti,counti, timenc)
      CALL check_err(iret)
!     ---------------------------
!     Store zeta
!     ---------------------------
      count(1)=node_len
      count(2)=timenc_len
      start(1)=1
      start(2)=nrec
      iret = nf_put_vara_DOUBLE(ncid, zeta_id, start, count, zeta)
      CALL check_err(iret)
!     ----------------------------------------
!     Allocating memory for min and max arrays
!     ----------------------------------------
      
      IF(.NOT. ALLOCATED(minzetat))ALLOCATE( minzetat(timenc_len))
      IF(.NOT. ALLOCATED(maxzetat))ALLOCATE( maxzetat(timenc_len))
      IF(.NOT. ALLOCATED(minzetan))ALLOCATE( minzetan(node_len))
      IF(.NOT. ALLOCATED(maxzetan))ALLOCATE( maxzetan(node_len))

      IF(.NOT. ALLOCATED(minznodet  ))ALLOCATE( minznodet  (timenc_len))
      IF(.NOT. ALLOCATED(maxznodet  ))ALLOCATE( maxznodet  (timenc_len))
      IF(.NOT. ALLOCATED(minztimen73))ALLOCATE( minztimen73(  node_len))
      IF(.NOT. ALLOCATED(maxztimen73))ALLOCATE( maxztimen73(  node_len))

!     -------------------------------------
!     Calculating zeta minimum in each node
!     -------------------------------------
      do n=1,node_len
         if(zeta (n,1).lt.minzetan(n) ) then
 	    minzetan    (n) = zeta (n,1)
 	    minztimen73 (n) = timenc(1)
         endif
      end do

!     ---------------------------------------
!     Calculating zeta minimum for each time
!     ---------------------------------------
      minzetat (1) = 9999.d0
      do n=1,node_len
         if(zeta (n,1).gt.-9999.d0.and.zeta(n,1).lt. minzetat (1) ) then 
           minzetat (1) = zeta (n,1)
           minznodet(1) = n
        endif
      end do

      iret=nf_put_vara_DOUBLE(ncid, minzetat_id,(/nrec/),(/1/),minzetat)
           CALL check_err(iret)
      iret=nf_put_vara_INT(ncid, minznodet_id,(/nrec/),(/1/),minznodet)
           CALL check_err(iret)

!     ---------------------------------------
!     Calculating zeta maximum in each node
!     ---------------------------------------
      do n=1,node_len
 	 if(zeta (n,1).gt.maxzetan(n) ) then
 	    maxzetan (n)    = zeta (n,1)
            maxztimen73 (n) = timenc (1)
         endif
      end do
!     ---------------------------------------
!     Calculating zeta maximum for each time
!     ---------------------------------------
      maxzetat (1) = -9999.d0
      do n=1,node_len
         if(zeta(n,1).gt.maxzetat(1))  then
            maxzetat (1) = zeta (n,1)
            maxznodet(1) = n
         endif
      end do

      iret=nf_put_vara_DOUBLE(ncid,maxzetat_id,(/nrec/),(/1/),maxzetat)
      CALL check_err(iret)
      iret=nf_put_vara_INT(ncid,maxznodet_id,(/nrec/),(/1/),maxznodet)
      CALL check_err(iret)

!     ---------------------------
!     Closing fort.X3 netCDF file
!     ---------------------------
      IF(NREC.EQ.NDATA) THEN
         iret = nf_put_var_double(ncid, minzetan_id, minzetan)
         CALL check_err(iret)
         iret = nf_put_var_double(ncid, maxzetan_id, maxzetan)
         CALL check_err(iret)
         iret = nf_put_var_double(ncid, minztimen_id, minztimen73)
         CALL check_err(iret)
         iret = nf_put_var_double(ncid, maxztimen_id, maxztimen73)
         CALL check_err(iret)
      ENDIF
      
      WRITE(6,*) "Closing ", FILENAME, " netCDF file ... "
      WRITE(6,*) "-----------------------------------------------------"
      
9999  iret = nf_close(ncid)
      CALL check_err(iret)
      WRITE(6,*)'end subroutine writenc73'
     
      end subroutine writenc73
!=======================================================================
!     Subroutine writenc74
!     Writes a netCDF file equivalent to ADCIRC's fort.74 ascii file
!
!     Dimensions:
!     node  : number of nodes in the horizontal grid
!     nele  : number of elements in the horizontal grid
!     nope  : number of elevation specified boundary forcing segments
!     neta  : total number of elevation specified boundary nodes
!     nbou  : number of normal flow specified boundary segments
!     nvel  : total number of normal flow specified boundary nodes
!     time  : number of time records
!
!     Variables :
!     time       : time
!     windx      : u velocity in x or east/west direction 
!     windy      : v velocity in y or north/south direction
!     x          : x-coordinate or longitude
!     y          : y-coordinate or latitude
!     element    : elements in grid
!     nvdll      : number of nodes in elevation boundary segment
!     nbdv       : node numbers on elevation specified boundary segment
!     nvell      : number of nodes in normal flow specified boundary segment
!     nbvv       : node numbers on normal flow boundary segment
!     ibtype     : boundary type
!     depth      : depth in each node
!     mindepth   : minimum depth in grid
!     maxdepth   : maximum depth in grid
!     minwxt     : minimum vertically averaged u at each time
!     maxwxt     : maximum vertically averaged u at each time
!     minwyt     : minimum vertically averaged v at each time
!     maxwyt     : maximum vertically averaged v at each time
!     minwxn     : minimum vertically averaged u at each node
!     maxwxn     : maximum vertically averaged u at each node
!     minwyn     : minimum vertically averaged v at each node
!     maxwyn     : maximum vertically averaged v at each node
!     minwxnodet : node where minimum u occurs at each output time
!     maxwxnodet : node where maximum u occurs at each output time
!     minwynodet : node where minimum v occurs at each output time
!     maxwynodet : node where maximum v occurs at each output time
!     minwxtimen : time when minimum u occurs at each node
!     maxwxtimen : time when maximum u occurs at each node
!     minwytimen : time when minimum v occurs at each node
!     maxwytimen : time when maximum v occurs at each node
!
!     Attributes:
!     model             : ADCIRC
!     version           : model version number
!     grid_type         : triangular
!     rundes            : description of simulation
!     agrid             : grid ID 
!     title             : description of the dataset
!     institution       : where the original data was produced
!     source            : method of production of the original data
!     history           : provides an audit trail for modifications to the original data
!     references        : published or web-based references that describe the data or methods used to produce it
!     comments          : miscellaneous information about the data
!     contact           : user's institution/ e-mail address
!     creation_date     : date the data was created
!     modification_date : date the data was modified
!     Parameters from fort.15 input file 
!
!     Revision history:
!
!     Date            Programmer                  Description of change
!     ----            ----------                  ----------------------
!     03/30/07  Cristina Forbes, PSGS @ UNC-IMS   Wrote original code
!     03/30/08  Cristina Forbes, PSGS @ UNC-IMS   Modified code for 
!                                                 globalio & hotstart 
!                                                 from binary files
!=======================================================================
      subroutine writenc74 (nout,fileunit,timenc_len,node_len,timenc,
     &                      nt, NDATA, u,v,
     &                      rundes, runid, agrid, irtype,
     &                      dtdp, ihot, ics, nolibf, nolifa,
     &                      nolica, nolicat, nwp, 
     &                      ncor, ntip, nws, nramp, tau0,
     &                      statim, reftim, rnday, dramp,
     &                      a00, b00, c00, h0,
     &                      slam0, sfea0, ffactor,esl, cori,
     &                      ntif, nbfr, elem_len, nface_len, nm,
     &                      x,y,depth,
     &                      nope_len, neta_len, nbdv, nopenc,netanc_len,
     &                      nvdllnc,
     &                      nbou_len, nvel_len, nbvv,nbounc,nvelnc_len,
     &                      nvellnc,ibtypenc,
     &                      ncmode,
     &                      title, institution, source, history,
     &                      references, comments, host, convention,
     &                      contact,
     &                      base_date)


!--------------------------------------------
!     Force explicit declaration of all variables
!--------------------------------------------
      IMPLICIT NONE
!     -------------------
!     netCDF header file
!     -------------------
 !     include 'netcdf.inc'
!     -------------------
!     Error status return
!     -------------------
      INTEGER  iret
!     ---------
!     netCDF id
!     ---------
      INTEGER  ncid
!     -------------
!     Dimension ids
!     -------------
      INTEGER  node_dim
      INTEGER  timenc_dim
      INTEGER  elem_dim
      INTEGER  nface_dim
      INTEGER  nope_dim
      INTEGER  neta_dim
      INTEGER  nbou_dim
      INTEGER  nvel_dim
      INTEGER  nvelnc_dim
      INTEGER  netanc_dim
!     -----------------
!     Dimension lengths
!     -----------------
      INTEGER  node_len
      INTEGER  timenc_len
      INTEGER  elem_len
      INTEGER  nface_len
      INTEGER  nope_len
      INTEGER  neta_len
      INTEGER  nbou_len
      INTEGER  nvel_len
      INTEGER  nvelnc_len
      INTEGER  netanc_len
!     ------------
!     Variable ids
!     ------------
      INTEGER  timenc_id
      INTEGER  X_id
      INTEGER  Y_id
      INTEGER  DEPTH_id
      INTEGER  ELE_id
      INTEGER  mindepth_id
      INTEGER  maxdepth_id
      INTEGER  nbdv_id
      INTEGER  nbvv_id
      INTEGER  u_id
      INTEGER  v_id
      INTEGER  minut_id
      INTEGER  maxut_id
      INTEGER  minvt_id
      INTEGER  maxvt_id
      INTEGER  minun_id
      INTEGER  maxun_id
      INTEGER  minvn_id
      INTEGER  maxvn_id
      
      INTEGER  minunodet_id
      INTEGER  maxunodet_id
      INTEGER  minutimen_id
      INTEGER  maxutimen_id
      INTEGER  minvnodet_id
      INTEGER  maxvnodet_id
      INTEGER  minvtimen_id
      INTEGER  maxvtimen_id
      
      INTEGER  nvdll_id
      INTEGER  nbounc_id
      INTEGER  ibtypenc_id
      INTEGER  nvellnc_id
      INTEGER  nvdllnc_id
      INTEGER  nopenc_id
!     ---------------------------------------------
!     Rank (number of dimensions) for each variable
!     ---------------------------------------------
      INTEGER  timenc_rank
      INTEGER  depth_rank
      INTEGER  ele_rank
      INTEGER  mindepth_rank
      INTEGER  maxdepth_rank
      INTEGER  x_rank
      INTEGER  y_rank

      INTEGER  nbounc_rank
      INTEGER  nvdll_rank
      INTEGER  ibtypenc_rank
      INTEGER  nvellnc_rank
      INTEGER  nvdllnc_rank
      INTEGER  nopenc_rank

      INTEGER  nbdv_rank
      INTEGER  nbvv_rank
      INTEGER  u_rank
      INTEGER  v_rank
      INTEGER  minut_rank
      INTEGER  maxut_rank
      INTEGER  minvt_rank
      INTEGER  maxvt_rank
      INTEGER  minun_rank
      INTEGER  maxun_rank
      INTEGER  minvn_rank
      INTEGER  maxvn_rank

      INTEGER  minunodet_rank
      INTEGER  maxunodet_rank
      INTEGER  minutimen_rank
      INTEGER  maxutimen_rank
      
      INTEGER  minvnodet_rank
      INTEGER  maxvnodet_rank
      INTEGER  minvtimen_rank
      INTEGER  maxvtimen_rank

      PARAMETER (timenc_rank = 1)
      PARAMETER (depth_rank = 1)
      PARAMETER (ele_rank = 2)
      PARAMETER (mindepth_rank = 0)
      PARAMETER (maxdepth_rank = 0)
      PARAMETER (x_rank = 1)
      PARAMETER (y_rank = 1)
      
      PARAMETER (nvdll_rank = 1)
      PARAMETER (nbounc_rank = 1)
      PARAMETER (ibtypenc_rank = 1)
      PARAMETER (nvellnc_rank = 1)
      PARAMETER (nvdllnc_rank=1)
      PARAMETER (nopenc_rank=1)

      PARAMETER (nbdv_rank = 2)
      PARAMETER (nbvv_rank = 2)
      PARAMETER (u_rank = 2)
      PARAMETER (v_rank = 2)
      PARAMETER (minut_rank = 1)
      PARAMETER (maxut_rank = 1)
      PARAMETER (minvt_rank = 1)
      PARAMETER (maxvt_rank = 1)
      PARAMETER (minun_rank = 1)
      PARAMETER (maxun_rank = 1)
      PARAMETER (minvn_rank = 1)
      PARAMETER (maxvn_rank = 1)

      PARAMETER (minunodet_rank = 1)
      PARAMETER (maxunodet_rank = 1)
      PARAMETER (minutimen_rank = 1)
      PARAMETER (maxutimen_rank = 1)

      PARAMETER (minvnodet_rank = 1)
      PARAMETER (maxvnodet_rank = 1)
      PARAMETER (minvtimen_rank = 1)
      PARAMETER (maxvtimen_rank = 1)

!     ---------------
!     Variable shapes
!     ---------------
      INTEGER  timenc_dims(timenc_rank)
      INTEGER  x_dims(x_rank)
      INTEGER  y_dims(y_rank)
      INTEGER  depth_dims(depth_rank)
      INTEGER  ele_dims(ele_rank)

      INTEGER  nvdll_dims(nvdll_rank)
      INTEGER  nbounc_dims(nbounc_rank)
      INTEGER  ibtypenc_dims(ibtypenc_rank)
      INTEGER  nvellnc_dims(nvellnc_rank)
      INTEGER  nvdllnc_dims(nvdllnc_rank)
      INTEGER  nopenc_dims(nopenc_rank)

      INTEGER  nbdv_dims(nbdv_rank)
      INTEGER  nbvv_dims(nbvv_rank)
      INTEGER  u_dims(u_rank)
      INTEGER  v_dims(v_rank)
      INTEGER  minut_dims(minut_rank)
      INTEGER  maxut_dims(maxut_rank)
      INTEGER  minvt_dims(minvt_rank)
      INTEGER  maxvt_dims(maxvt_rank)
      INTEGER  minun_dims(minun_rank)
      INTEGER  maxun_dims(maxun_rank)
      INTEGER  minvn_dims(minvn_rank)
      INTEGER  maxvn_dims(maxvn_rank)
      
      INTEGER  minunodet_dims(minunodet_rank)
      INTEGER  maxunodet_dims(maxunodet_rank)
      INTEGER  minutimen_dims(minutimen_rank)
      INTEGER  maxutimen_dims(maxutimen_rank)

      INTEGER  minvnodet_dims(minvnodet_rank)
      INTEGER  maxvnodet_dims(maxvnodet_rank)
      INTEGER  minvtimen_dims(minvtimen_rank)
      INTEGER  maxvtimen_dims(maxvtimen_rank)
!     --------------
!     Data variables
!     --------------
      INTEGER t,ntif,nbfr,nt
      INTEGER irtype, ihot, ics, nolibf, nolifa, nolica, nolicat, nwp
      INTEGER ntip, nws, nramp,ncor
      INTEGER nvel, neta, nout
      INTEGER element (nface_len , elem_len)
      INTEGER nbdv (nope_len, neta_len)
      INTEGER nbvv (nbou_len, nvel_len)
      
      INTEGER   nvellnc(nbou_len)
      INTEGER   nbounc (nvel_len)
      INTEGER   ibtype(nbou_len)
      INTEGER   ibtypenc(nbou_len)
      INTEGER   nvdllnc (nope_len)
      INTEGER   nopenc (neta_len)

      INTEGER nm (elem_len, nface_len)
      INTEGER fileunit
      CHARACTER*5 FILEFORT
      CHARACTER*10 FILENAME
      CHARACTER*3 FILEEXT
      REAL*8    timenc(timenc_len)
      REAL(sz)  x(node_len)
      REAL(sz)  y(node_len)
      REAL(sz)  depth(node_len)
      REAL(sz)  u(node_len, timenc_len)
      REAL(sz)  v(node_len, timenc_len)
      REAL(sz), ALLOCATABLE, SAVE ::   minut(:)
      REAL(sz), ALLOCATABLE, SAVE ::   maxut(:)
      REAL(sz), ALLOCATABLE, SAVE ::   minvt(:)
      REAL(sz), ALLOCATABLE, SAVE ::   maxvt(:)
      REAL(sz), ALLOCATABLE, SAVE ::   minun(:)
      REAL(sz), ALLOCATABLE, SAVE ::   maxun(:)
      REAL(sz), ALLOCATABLE, SAVE ::   minvn(:)
      REAL(sz), ALLOCATABLE, SAVE ::   maxvn(:)

      INTEGER , ALLOCATABLE, SAVE ::  minunodet(:)
      INTEGER , ALLOCATABLE, SAVE ::  maxunodet(:)
      REAL(sz), ALLOCATABLE, SAVE ::  minutimen74(:)
      REAL(sz), ALLOCATABLE, SAVE ::  maxutimen74(:)
      
      INTEGER , ALLOCATABLE, SAVE ::  minvnodet(:)
      INTEGER , ALLOCATABLE, SAVE ::  maxvnodet(:)
      REAL(sz), ALLOCATABLE, SAVE ::  minvtimen74(:)
      REAL(sz), ALLOCATABLE, SAVE ::  maxvtimen74(:)

      REAL  minull
      REAL  maxull
      REAL  minvll
      REAL  maxvll
      real*8  mindepth
      real*8  maxdepth
      REAL*8 dtdp, statim,reftim, rnday, dramp, tau0
      REAL*8 a00, b00, c00,h0
      REAL(sz) slam0, sfea0, ffactor,esl, cori
      REAL(sz) slam0deg, sfea0deg
      CHARACTER(80) rundes, runid, agrid
      INTEGER, SAVE  :: FIRST74 = 1 
      INTEGER NDATA
      INTEGER, SAVE  :: NREC = 0
      INTEGER ncmode

!     -----------------
!     Attribute vectors
!     -----------------
      INTEGER intval(1),i,j
      INTEGER counti(1), starti(1), n
      INTEGER count(2), start(2)
      REAL    realval(1)
!     -----------------
!     Attribute general
!     -----------------
      CHARACTER(*) :: title, institution, source, history
      CHARACTER(*) :: references, comments, host, convention, contact
      CHARACTER(*) :: base_date

!     -----------------
!      date_string variables for time attribute
!     -----------------
      character date_string*40
      character now_date*8
      character big_ben*10
      character zone*5
      integer   values(8)
!     ---------------------------------
!      inquiring and reading variables
!     ---------------------------------

      INTEGER NDIMS, NVARS, NGATTS, UNLIMDIMID, RECID,NRECS
      CHARACTER*4 recname

!     -------------------------
!     Save dims, ids, variables
!     -------------------------
      save node_dim,  elem_dim, nface_dim
      save nope_dim, neta_dim, nbou_dim, nvel_dim
      save u_id, v_id, X_id, Y_id
      save minut_id, maxut_id, minun_id, maxun_id
      save minvt_id, maxvt_id, minvn_id, maxvn_id
      
      save minunodet_id, maxunodet_id, minutimen_id, maxutimen_id
      save minvnodet_id, maxvnodet_id, minvtimen_id, maxvtimen_id

      save DEPTH_id,ELE_id, mindepth_id, maxdepth_id, nbdv_id, nbvv_id
      save timenc_dim, timenc_id
      save mindepth, maxdepth

      save nvelnc_dim
      save nvdll_id, nbounc_id, ibtypenc_id, nvellnc_id
      save netanc_dim
      save nopenc_id,nvdllnc_id

!     ---------------------------
!     Opening fort.74 netCDF file
!     ---------------------------
      WRITE(6,*) "-----------------------------------------------------"
      IF(FIRST74.EQ.1.and.ncmode.NE.0) THEN
         
         ALLOCATE ( minun (node_len) )
         ALLOCATE ( maxun (node_len) )
         ALLOCATE ( minvn (node_len) )
         ALLOCATE ( maxvn (node_len) )
         
         ALLOCATE ( minutimen74 (node_len) )
         ALLOCATE ( maxutimen74 (node_len) )
         ALLOCATE ( minvtimen74 (node_len) )
         ALLOCATE ( maxvtimen74 (node_len) )

         do n=1,node_len
            minun (n) =  9999.d0
            maxun (n) = -9999.d0
            minvn (n) =  9999.d0
            maxvn (n) = -9999.d0
            minutimen74 (n) =  9999.d0
            maxutimen74 (n) = -9999.d0
            minvtimen74 (n) =  9999.d0
            maxvtimen74 (n) = -9999.d0
         end do
      ENDIF
!     --------------------------------
!     Switch order in array for NETCDF
!     --------------------------------
      do i=1, elem_len
         do j=1, nface_len
            element(j,i) = nm (i,j)
         end do
      end do
      
!     -----------------
!     Enter define mode
!     -----------------
         FILEFORT = 'fort.'
         FILEEXT = '.nc'
         WRITE(FILENAME,'(A5,I2,A3)')FILEFORT,FILEUNIT,FILEEXT
         WRITE(6,*) "Opening ",FILENAME,  " netCDF file ... "

       IF((FIRST74.EQ.1.AND.IHOT.EQ.0.AND.NCMODE.NE.2).OR.
     &    (FIRST74.EQ.1.AND.IHOT.NE.0.AND.NCMODE.NE.2.AND.
     &     nout.EQ.-3))THEN

           iret = nf_create(FILENAME, NF_CLOBBER, ncid)
           CALL check_err(iret)
           FIRST74=0

!     --------------------------------
!     Convert to degress for NETCDF I/O
!     --------------------------------
           SLAM0DEG=SLAM0
           SFEA0DEG=SFEA0
           IF (NCMODE.EQ.1) THEN
               SLAM0DEG=SLAM0*RAD2DEG
               SFEA0DEG=SFEA0*RAD2DEG
           ENDIF
!     -----------------
!     Define dimensions
!     -----------------
      iret = nf_def_dim(ncid, 'node', node_len, node_dim)
      CALL check_err(iret)
      iret = nf_def_dim(ncid, 'nele', elem_len, elem_dim)
      CALL check_err(iret)
      iret = nf_def_dim(ncid, 'nface', 3, nface_dim)
      CALL check_err(iret)
      if(nope_len.ne.0) then
         iret = nf_def_dim(ncid, 'nope', nope_len, nope_dim)
         call check_err(iret)
      endif
      if(neta_len.ne.0) then
         iret = nf_def_dim(ncid, 'neta', neta_len, neta_dim)
         call check_err(iret)
      endif
      if(nbou_len.ne.0) then
         iret = nf_def_dim(ncid, 'nbou', nbou_len, nbou_dim)
         call check_err(iret)
      endif
      if(nvel_len.ne.0) then
         iret = nf_def_dim(ncid, 'nvel', nvel_len, nvel_dim)
         call check_err(iret)
      endif
      
      iret = nf_def_dim(ncid, 'time', nf_unlimited, timenc_dim)
      CALL check_err(iret)

!     ================
!     Define variables
!     ================
!     ----------------
!     Define time
!     ----------------
      timenc_dims(1) = timenc_dim
      iret = nf_def_var(ncid,'time',NF_DOUBLE, timenc_rank, timenc_dims, 
     &                  timenc_id)
      CALL check_err(iret)
      
!     -------------------
!     Define x
!     -------------------
      x_dims(1) = node_dim
      iret = nf_def_var(ncid, 'x', NF_DOUBLE, x_rank, x_dims, x_id)
      CALL check_err(iret)
!     -------------------
!     Define Y coordinate
!     --------------------
      Y_dims(1) = node_dim
      iret = nf_def_var(ncid, 'y', NF_DOUBLE, Y_rank, Y_dims, Y_id)
      CALL check_err(iret)
!     -------------------
!     Define elements
!     --------------------
      ELE_dims(1) = nface_dim
      ELE_dims(2) = elem_dim
      iret = nf_def_var(ncid, 'element',NF_INT,ELE_rank,ELE_dims,ELE_id)
      CALL check_err(iret)

!     -----------------------------------------------------------------
!     Define elevation specified boundary forcing segments information
!     -----------------------------------------------------------------
      if(nope_len.ne.0.and.neta_len.ne.0) then
        nvdllnc_dims(1) = nope_dim      
        iret = nf_def_var(ncid, 'nvdll',NF_INT, nvdllnc_rank, 
     &                  nvdllnc_dims, nvdllnc_id)
        CALL check_err(iret)

        nopenc_dims(1) = neta_dim      
        iret = nf_def_var(ncid, 'nbdv',NF_INT, nopenc_rank, nopenc_dims, 
     &                  nopenc_id)
        CALL check_err(iret)

      endif
!     ----------------------------------------------------------
!     Define normal flow specified boundary segments information
!     ----------------------------------------------------------
      if(nbou_len.ne.0.and.nvel_len.ne.0) then

        nvellnc_dims(1) = nbou_dim      
        iret = nf_def_var(ncid, 'nvell',NF_INT, nvellnc_rank, 
     &                  nvellnc_dims, nvellnc_id)
        CALL check_err(iret)
        
        ibtypenc_dims(1) = nbou_dim      
        iret = nf_def_var(ncid, 'ibtype',NF_INT, ibtypenc_rank, 
     &                  ibtypenc_dims, ibtypenc_id)
        CALL check_err(iret)
        nbounc_dims(1) = nvel_dim      
        iret = nf_def_var(ncid, 'nbvv',NF_INT, nbounc_rank, nbounc_dims, 
     &                  nbounc_id)
        CALL check_err(iret)
     
       endif
!     -------------------
!     Define Z coordinate
!     --------------------
      DEPTH_dims(1) = node_dim
      iret = nf_def_var(ncid, 'depth', NF_DOUBLE,DEPTH_rank, DEPTH_dims, 
     &                  DEPTH_id)
      CALL check_err(iret)
!     --------------------
!     Define minimum depth
!     --------------------
      iret = nf_def_var(ncid, 'mindepth', NF_DOUBLE, mindepth_rank, 0, 
     &                  mindepth_id)
      CALL check_err(iret)
!     --------------------
!     Define maximum depth
!     --------------------
      iret = nf_def_var(ncid, 'maxdepth', NF_DOUBLE, maxdepth_rank, 0, 
     &                  maxdepth_id)
      CALL check_err(iret)

!     ----------------
!     Define u vel
!     ----------------
      u_dims(1) = node_dim
      u_dims(2) = timenc_dim
      iret = nf_def_var(ncid, 'windx', NF_DOUBLE, u_rank, u_dims, 
     &                  u_id)
      CALL check_err(iret)
!     ----------------
!     Define v vel
!     ----------------
      v_dims(1) = node_dim
      v_dims(2) = timenc_dim
      iret = nf_def_var(ncid, 'windy', NF_DOUBLE, v_rank, v_dims, 
     &                  v_id)
      CALL check_err(iret)
!     -------------------------
!     Define minimum u at times
!     -------------------------
      minut_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'minwxt', NF_DOUBLE, minut_rank, 
     &                  minut_dims, minut_id)
      CALL check_err(iret)
!     -------------------------
!     Define maximum u at times
!     -------------------------
      maxut_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'maxwxt', NF_DOUBLE, maxut_rank, 
     &                  maxut_dims, maxut_id)
      CALL check_err(iret)
!     -------------------------
!     Define minimum v at times
!     -------------------------
      minvt_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'minwyt', NF_DOUBLE, minvt_rank, 
     &                  minvt_dims, minvt_id)
      CALL check_err(iret)
!     -------------------------
!     Define maximum v at times
!     -------------------------
      maxvt_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'maxwyt', NF_DOUBLE, maxvt_rank, 
     &                  maxvt_dims, maxvt_id)
      CALL check_err(iret)
!     -------------------------
!     Define minimum u at times
!     -------------------------
      minun_dims(1) = node_dim
      iret = nf_def_var(ncid, 'minwxn', NF_DOUBLE, minun_rank, 
     &                  minun_dims, minun_id)
      CALL check_err(iret)
!     -------------------------
!     Define maximum u at nodes
!     -------------------------
      maxun_dims(1) = node_dim
      iret = nf_def_var(ncid, 'maxwxn', NF_DOUBLE, maxun_rank, 
     &                  maxun_dims, maxun_id)
      CALL check_err(iret)
!     -------------------------
!     Define minimum v at nodes
!     -------------------------
      minvn_dims(1) = node_dim
      iret = nf_def_var(ncid, 'minwyn', NF_DOUBLE, minvn_rank, 
     &                  minvn_dims, minvn_id)
      CALL check_err(iret)
!     -------------------------
!     Define maximum v at nodes
!     -------------------------
      maxvn_dims(1) = node_dim
      iret = nf_def_var(ncid, 'maxwyn', NF_DOUBLE, maxvn_rank, 
     &                  maxvn_dims, maxvn_id)
      CALL check_err(iret)
!     ----------------------------
!     Define minimum u nodes at times
!     ----------------------------
      minunodet_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'minwxnodet', NF_INT, minunodet_rank, 
     &                  minunodet_dims, minunodet_id)
      CALL check_err(iret)
!     ----------------------------
!     Define maximum u nodes at times
!     ----------------------------
      maxunodet_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'maxwxnodet', NF_INT, maxunodet_rank, 
     &                  maxunodet_dims, maxunodet_id)
      CALL check_err(iret)
!     ----------------------------
!     Define minimum u nodes at nodes
!     ----------------------------
      minutimen_dims(1) = node_dim
      iret = nf_def_var(ncid, 'minwxtimen', NF_DOUBLE, minutimen_rank, 
     &                  minutimen_dims, minutimen_id)
      CALL check_err(iret)
!     ----------------------------
!     Define maximum u nodes at nodes
!     ----------------------------
      maxutimen_dims(1) = node_dim
      iret = nf_def_var(ncid, 'maxwxtimen', NF_DOUBLE, maxutimen_rank, 
     &                  maxutimen_dims, maxutimen_id)
      CALL check_err(iret)

            
!     ----------------------------
!     Define minimum v nodes at times
!     ----------------------------
      minvnodet_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'minwynodet', NF_INT, minvnodet_rank, 
     &                  minvnodet_dims, minvnodet_id)
      CALL check_err(iret)
!     ----------------------------
!     Define maximum v nodes at times
!     ----------------------------
      maxvnodet_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'maxwynodet', NF_INT, maxvnodet_rank, 
     &                  maxvnodet_dims, maxvnodet_id)
      CALL check_err(iret)
!     ----------------------------
!     Define minimum v nodes at nodes
!     ----------------------------
      minvtimen_dims(1) = node_dim
      iret = nf_def_var(ncid, 'minwytimen', NF_DOUBLE, minvtimen_rank, 
     &                  minvtimen_dims, minvtimen_id)
      CALL check_err(iret)
!     ----------------------------
!     Define maximum v nodes at nodes
!     ----------------------------
      maxvtimen_dims(1) = node_dim
      iret = nf_def_var(ncid, 'maxwytimen', NF_DOUBLE, maxvtimen_rank, 
     &                  maxvtimen_dims, maxvtimen_id)
      CALL check_err(iret)      
!     -----------------
!     Assign attributes
!     -----------------

      if(ics.eq.2) then
      
       iret = nf_put_att_text(ncid, X_id, 'long_name', 9, 'longitude')
       CALL check_err(iret)
       iret = nf_put_att_text(ncid, X_id, 'standard_name', 9, 
     &                       'longitude')
       CALL check_err(iret)
       iret = nf_put_att_text(ncid, X_id, 'units', 12, VAR_DEGE)
       CALL check_err(iret)
       iret = nf_put_att_text(ncid, X_id, 'positive', 4, VAR_EAST)
       CALL check_err(iret)
       iret = nf_put_att_double(ncid, X_id, 'missing_values', NF_DOUBLE,
     &                          1, doubleval)
       CALL check_err(iret)
       iret = nf_put_att_double(ncid, X_id, '_FillValue', NF_DOUBLE, 1, 
     &                          doubleval)
       CALL check_err(iret)      
      
       iret = nf_put_att_text(ncid, Y_id, 'long_name', 8, 'latitude')
       CALL check_err(iret)
       iret = nf_put_att_text(ncid, Y_id, 'standard_name', 8, 
     &                       'latitude')
       CALL check_err(iret)
       iret = nf_put_att_text(ncid, Y_id, 'units', 13,  VAR_DEGN)
       CALL check_err(iret)
       iret = nf_put_att_text(ncid, Y_id, 'positive', 5, VAR_NORTH)
       CALL check_err(iret)
       
       iret = nf_put_att_double(ncid, Y_id, 'missing_values', NF_DOUBLE, 
     &       1, doubleval)
       CALL check_err(iret)
       iret = nf_put_att_double(ncid, Y_id, '_FillValue', NF_DOUBLE, 1, 
     &       doubleval)
       CALL check_err(iret)

      else
      
       iret = nf_put_att_text(ncid, X_id, 'long_name', 22, 
     &       'Cartesian coordinate x')
       CALL check_err(iret)
       iret = nf_put_att_text(ncid, X_id, 'standard_name', 12, 
     &                      'x_coordinate')
       CALL check_err(iret)
!    -------------------------
!     determine variable units
!    -------------------------
      if (G.LT.11.D0) then
         iret = nf_put_att_text(ncid, X_id, 'units', 6,  VAR_METS) 
      else
         iret = nf_put_att_text(ncid, X_id, 'units', 4,  VAR_FEET) 
      endif
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, X_id, 'positive', 5, VAR_RG)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, X_id, 'missing_values', NF_DOUBLE, 
     &                       1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, X_id, '_FillValue', NF_DOUBLE, 1, 
     &                       doubleval)
      CALL check_err(iret)
!     ------------------------------
!     Define y-coordinate attributes
!     ------------------------------     
      
      iret = nf_put_att_text(ncid, Y_id, 'long_name', 22, 
     &       'Cartesian coordinate y')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, Y_id, 'standard_name', 12, 
     &                       'y_coordinate')
      CALL check_err(iret)
!    -------------------------
!     determine variable units
!    -------------------------
      if (G.LT.11.D0) then
         iret = nf_put_att_text(ncid, Y_id, 'units', 6, VAR_METS)
      else
         iret = nf_put_att_text(ncid, Y_id, 'units', 4, VAR_FEET) 
      endif
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, Y_id, 'positive', 2, VAR_UP)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, Y_id, 'missing_values', 
     &       NF_DOUBLE, 1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, Y_id, '_FillValue', NF_DOUBLE, 1, 
     &       doubleval)
      CALL check_err(iret)
      
      endif

!     ------------------------------
!     Define bathymetry attributes
!     ------------------------------

      iret = nf_put_att_text(ncid, DEPTH_id,'long_name',10,'bathymetry')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, DEPTH_id, 'standard_name', 5,'depth')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, DEPTH_id, 'units', 6, VAR_METS)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, DEPTH_id, 'missing_values', 
     &                      NF_DOUBLE,  1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid,DEPTH_id, '_FillValue',NF_DOUBLE, 1, 
     &                       doubleval)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, DEPTH_id, 'positive', 4, VAR_DW)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, ELE_id, 'long_name', 7, 'element')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, ELE_id, 'units', 14,  VAR_NOND)
      intval(1) = -1
!    ---------------------------
!     minimum/maximum bathymetry
!    ---------------------------
      iret = nf_put_att_text(ncid, mindepth_id, 'long_name', 18, 
     &                       'minimum bathymetry')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, mindepth_id, 'standard_name', 13,
     &                'minimum_depth')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, mindepth_id, 'units', 6, VAR_METS)
      CALL check_err(iret)
      
      iret = nf_put_att_text(ncid, maxdepth_id, 'long_name', 18, 
     &                       'maximum bathymetry')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, maxdepth_id, 'standard_name', 13,
     &           'maximum_depth')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, maxdepth_id, 'units', 6, VAR_METS)
      CALL check_err(iret)
      
!    ----------------------------------------------------------
!     determine variable names according to coordinates and nws
!    ----------------------------------------------------------
      if (ics.eq.2) then
        if (abs(nws).eq.1.or.abs(nws).eq.2) then
             iret = nf_put_att_text(ncid, minut_id, 'long_name', 36, 
     &       'minimum e/w wind stress at each time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minut_id, 'standard_name', 23,
     &       'minimum_e/w_wind_stress')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvt_id, 'long_name', 36, 
     &       'minimum n/s wind stress at each time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvt_id, 'standard_name', 23,
     &       'minimum_n/s_wind_stress')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxut_id, 'long_name', 36, 
     &       'maximum e/w wind stress at each time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxut_id, 'standard_name', 23,
     &       'maximum_e/w_wind_stress')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvt_id, 'long_name', 36, 
     &       'maximum n/s wind stress at each time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvt_id, 'standard_name', 23,
     &       'maximum_n/s_wind_stress')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minun_id, 'long_name', 36, 
     &       'minimum e/w wind stress at each node')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minun_id, 'standard_name', 23,
     &       'minimum_e/w_wind_stress')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxun_id, 'long_name', 36, 
     &        'maximum e/w wind stress at each node')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxun_id, 'standard_name', 23,
     &       'maximum_e/w_wind_stress')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvn_id, 'long_name', 36, 
     &       'minimum n/s wind stress at each node')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvn_id, 'standard_name', 23,
     &       'minimum_n/s_wind_stress')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvn_id, 'long_name', 36, 
     &       'maximum n/s wind stress at each node')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvn_id, 'standard_name', 23,
     &       'maximum_n/s_wind_stress')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, minunodet_id, 'long_name', 
     &        51, 'node of minimum e/w wind stress at each output time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minunodet_id, 'standard_name',
     &        31,'node_of_minimum_e/w_wind_stress')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, maxunodet_id, 'long_name', 
     &        51, 'node of maximum e/w wind stress at each output time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxunodet_id, 'standard_name',
     &        31,'node_of_maximum_e/w_wind_stress')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, minvnodet_id, 'long_name', 
     &        51, 'node of minimum n/s wind stress at each output time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvnodet_id, 'standard_name',
     &        31,'node_of_minimum_n/s_wind_stress')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, maxvnodet_id, 'long_name', 
     &        51, 'node of maximum n/s wind stress at each output time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvnodet_id, 'standard_name',
     &        31,'node_of_maximum_n/s_wind_stress')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, minutimen_id, 'long_name', 
     &        44, 'time of minimum e/w wind stress at each node')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minutimen_id, 'standard_name',
     &        31,'time_of_minimum_e/w_wind_stress')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, maxutimen_id, 'long_name', 
     &        44, 'time of maximum e/w wind stress at each node')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxutimen_id, 'standard_name',
     &        31,'time_of_maximum_e/w_wind_stress')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, minvtimen_id, 'long_name', 
     &        44, 'time of minimum n/s wind stress at each node')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvtimen_id, 'standard_name',
     &        31,'time_of_minimum_n/s_wind_stress')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, maxvtimen_id, 'long_name', 
     &        44, 'time of maximum n/s wind stress at each node')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvtimen_id, 'standard_name',
     &        31,'time_of_maximum_n/s_wind_stress')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, u_id, 'long_name', 15, 
     &                       'e/w wind stress')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, u_id, 'standard_name', 15, 
     &                      'e/w_wind_stress')
             CALL check_err(iret)
             iret = nf_put_att_text(ncid, u_id, 'positive', 4, VAR_EAST)

              iret = nf_put_att_text(ncid, v_id, 'long_name', 15, 
     &                       'n/s wind stress')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, v_id, 'standard_name', 15, 
     &                      'n/s_wind_stress')
             CALL check_err(iret)
            iret = nf_put_att_text(ncid, v_id, 'positive', 5, VAR_NORTH)



        else
              iret = nf_put_att_text(ncid, minut_id, 'long_name', 38, 
     &       'minimum e/w wind velocity at each time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minut_id, 'standard_name', 25,
     &       'minimum_e/w_wind_velocity')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, minvt_id, 'long_name', 38, 
     &       'minimum n/s wind velocity at each time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvt_id, 'standard_name', 25,
     &       'minimum_n/s_wind_velocity')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, maxut_id, 'long_name', 38, 
     &       'maximum e/w wind velocity at each time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxut_id, 'standard_name', 25,
     &       'maximum_e/w_wind_velocity')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, maxvt_id, 'long_name', 38, 
     &       'maximum n/s wind velocity at each time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvt_id, 'standard_name', 25,
     &       'maximum_n/s_wind_velocity')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, minun_id, 'long_name', 38, 
     &       'minimum e/w wind velocity at each node')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minun_id, 'standard_name', 25,
     &       'minimum_e/w_wind_velocity')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxun_id, 'long_name', 38, 
     &       'maximum e/w wind velocity at each node')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxun_id, 'standard_name', 25,
     &       'maximum_e/w_wind_velocity')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvn_id, 'long_name', 38, 
     &       'minimum n/s wind velocity at each node')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvn_id, 'standard_name', 25,
     &       'minimum_n/s_wind_velocity')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvn_id, 'long_name', 38, 
     &       'maximum n/s wind velocity at each node')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvn_id,'standard_name', 25,
     &       'maximum_n/s_wind_velocity')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, minunodet_id, 'long_name', 
     &       53,'node of minimum e/w wind velocity at each output time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minunodet_id, 'standard_name',
     &        33,'node_of_minimum_e/w_wind_velocity')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, maxunodet_id, 'long_name', 
     &       53,'node of maximum e/w wind velocity at each output time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxunodet_id, 'standard_name',
     &        33,'node_of_maximum_e/w_wind_velocity')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, minvnodet_id, 'long_name', 
     &       53,'node of minimum n/s wind velocity at each output time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvnodet_id, 'standard_name',
     &        33,'node_of_minimum_n/s_wind_velocity')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, maxvnodet_id, 'long_name', 
     &       53,'node of maximum n/s wind velocity at each output time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvnodet_id, 'standard_name',
     &        33,'node_of_maximum_n/s_wind_velocity')
              CALL check_err(iret)

              iret = nf_put_att_text(ncid, minutimen_id, 'long_name', 
     &       46,'time of minimum e/w wind velocity at each node')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minutimen_id, 'standard_name',
     &        33,'time_of_minimum_e/w_wind_velocity')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, maxutimen_id, 'long_name', 
     &       46,'time of maximum e/w wind velocity at each node')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxutimen_id, 'standard_name',
     &        33,'time_of_maximum_e/w_wind_velocity')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, minvtimen_id, 'long_name', 
     &       46,'time of minimum n/s wind velocity at each node')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvtimen_id, 'standard_name',
     &        33,'time_of_minimum_n/s_wind_velocity')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, maxvtimen_id, 'long_name', 
     &       46,'time of maximum n/s wind velocity at each node')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvtimen_id, 'standard_name',
     &        33,'time_of_maximum_n/s_wind_velocity')
              CALL check_err(iret)

              iret = nf_put_att_text(ncid, u_id, 'long_name', 17, 
     &                       'e/w wind velocity')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, u_id, 'standard_name', 17, 
     &                      'e/w_wind_velocity')
             CALL check_err(iret)
             iret = nf_put_att_text(ncid, u_id, 'positive', 4, VAR_EAST)

              iret = nf_put_att_text(ncid, v_id, 'long_name', 17, 
     &                       'n/s wind velocity')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, v_id, 'standard_name', 17, 
     &                      'n/s_wind_velocity')
             CALL check_err(iret)
            iret = nf_put_att_text(ncid, v_id, 'positive', 5, VAR_NORTH)


          endif
       else
         if (abs(nws).eq.1.or.abs(nws).eq.2) then
             iret = nf_put_att_text(ncid, minut_id, 'long_name', 47, 
     &       'minimum wind stress in x-direction at each time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minut_id, 'standard_name', 35,
     &       'minimum_wind_stress_in_x_direction')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvt_id, 'long_name', 47, 
     &       'minimum wind stress in y-direction at each time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvt_id, 'standard_name', 35,
     &       'minimum_wind_stress_in_y_direction')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxut_id, 'long_name', 47, 
     &       'maximum wind stress in x-direction at each time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxut_id, 'standard_name', 35,
     &       'maximum_wind_stress_in_x_direction')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvt_id, 'long_name', 47, 
     &       'maximum wind stress in y-direction at each time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvt_id, 'standard_name', 35,
     &       'maximum_wind_stress_in_y_direction')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minun_id, 'long_name', 47, 
     &       'minimum wind stress in x-direction at each node')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minun_id, 'standard_name', 35,
     &       'minimum_wind_stress_in_x_direction')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxun_id, 'long_name', 47, 
     &       'maximum wind stress in x-direction at each node')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxun_id, 'standard_name', 35,
     &       'maximum_wind_stress_in_x_direction')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvn_id, 'long_name', 47, 
     &       'minimum wind stress in y-direction at each node')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvn_id, 'standard_name', 35,
     &       'minimum_wind_stress_in_y_direction')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvn_id, 'long_name', 47, 
     &       'maximum wind stress in y-direction at each node')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvn_id, 'standard_name', 35,
     &       'maximum_wind_stress_in_y_direction')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minunodet_id, 'long_name', 62, 
     &       'node of minimum wind stress in x-direction at each output 
     &time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minunodet_id, 'standard_name', 
     &       43,'node_of_minimum_wind_stress_in_x_direction')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxunodet_id, 'long_name', 62, 
     &       'node of maximum wind stress in x-direction at each output 
     &time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxunodet_id, 'standard_name', 
     &       43,'node_of_maximum_wind_stress_in_x_direction')
              CALL check_err(iret)

             iret = nf_put_att_text(ncid, minvnodet_id, 'long_name', 62, 
     &       'node of minimum wind stress in y-direction at each output 
     &time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvnodet_id, 'standard_name', 
     &       43,'node_of_minimum_wind_stress_in_y_direction')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvnodet_id, 'long_name', 62, 
     &       'node of maximum wind stress in y-direction at each output 
     &time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvnodet_id, 'standard_name', 
     &       43,'node_of_maximum_wind_stress_in_y_direction')
              CALL check_err(iret)

             iret = nf_put_att_text(ncid, minutimen_id, 'long_name', 55, 
     &       'time of minimum wind stress in x-direction at each node')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minutimen_id, 'standard_name', 
     &       43,'time_of_minimum_wind_stress_in_x_direction')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxutimen_id, 'long_name', 55, 
     &       'time of maximum wind stress in x-direction at each node')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxutimen_id, 'standard_name', 
     &       43,'time_of_maximum_wind_stress_in_x_direction')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvtimen_id, 'long_name', 55, 
     &       'time of minimum wind stress in y-direction at each node')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvtimen_id, 'standard_name', 
     &       43,'time_of_minimum_wind_stress_in_y_direction')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvtimen_id, 'long_name', 55, 
     &       'time of maximum wind stress in y-direction at each node')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvtimen_id, 'standard_name', 
     &       43,'time_of_maximum_wind_stress_in_y_direction')
              CALL check_err(iret)

              iret = nf_put_att_text(ncid, u_id, 'long_name', 26, 
     &                       'wind stress in x-direction')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, u_id, 'standard_name', 26, 
     &                      'wind_stress_in_x_direction')
             CALL check_err(iret)
             iret = nf_put_att_text(ncid, u_id, 'positive', 5, VAR_RG)

              iret = nf_put_att_text(ncid, v_id, 'long_name', 26, 
     &                       'wind stress in y-direction')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, v_id, 'standard_name', 26, 
     &                      'wind_stress_in_y_direction')
             CALL check_err(iret)
             iret = nf_put_att_text(ncid, v_id, 'positive', 2, VAR_UP)

              

         else
             iret = nf_put_att_text(ncid, minut_id, 'long_name', 49, 
     &       'minimum wind velocity in x-direction at each time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minut_id, 'standard_name', 36,
     &       'minimum_wind_velocity_in_x_direction')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvt_id, 'long_name', 49, 
     &       'minimum wind velocity in y-direction at each time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvt_id, 'standard_name', 36,
     &       'minimum_wind_velocity_in_y_direction')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxut_id, 'long_name', 49, 
     &       'maximum wind velocity in x-direction at each time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxut_id, 'standard_name', 36,
     &       'maximum_wind_velocity_in_x_direction')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvt_id, 'long_name', 49, 
     &       'maximum wind velocity in y-direction at each time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvt_id, 'standard_name', 36,
     &       'maximum_wind_velocity_in_y_direction')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minun_id, 'long_name', 49, 
     &       'minimum wind velocity in x-direction at each node')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minun_id, 'standard_name', 36,
     &       'minimum_wind_velocity_in_x_direction')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxun_id, 'long_name', 49, 
     &       'maximum wind velocity in x-direction at each node')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxun_id, 'standard_name', 36,
     &       'maximum_wind_velocity_in_x_direction')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvn_id, 'long_name', 49, 
     &       'minimum wind velocity in y-direction at each node')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvn_id, 'standard_name', 36,
     &       'minimum_wind_velocity_in_y_direction')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvn_id, 'long_name', 49, 
     &       'maximum wind velocity in y-direction at each node')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvn_id, 'standard_name', 36,
     &       'maximum_wind_velocity_in_y_direction')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minunodet_id, 'long_name', 64, 
     &            'node of minimum wind velocity in x-direction at each 
     &output time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minunodet_id, 'standard_name', 
     &       54,'node_of_minimum_wind_velocity_in_x_direction')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxunodet_id, 'long_name', 64, 
     &            'node of maximum wind velocity in x-direction at each 
     &output time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxunodet_id, 'standard_name', 
     &       54,'node_of_maximum_wind_velocity_in_x_direction')
              CALL check_err(iret)

             iret = nf_put_att_text(ncid, minvnodet_id, 'long_name', 64, 
     &            'node of minimum wind velocity in y-direction at each 
     &output time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvnodet_id, 'standard_name', 
     &       54,'node_of_minimum_wind_velocity_in_y_direction')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvnodet_id, 'long_name', 64, 
     &            'node of maximum wind velocity in y-direction at each 
     &output time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvnodet_id, 'standard_name', 
     &       54,'node_of_maximum_wind_velocity_in_y_direction')
              CALL check_err(iret)

             iret = nf_put_att_text(ncid, minutimen_id, 'long_name', 57, 
     &            'time of minimum wind velocity in x-direction at each 
     &node')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minutimen_id, 'standard_name', 
     &       54,'time_of_minimum_wind_velocity_in_x_direction')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxutimen_id, 'long_name', 57, 
     &            'time of maximum wind velocity in x-direction at each 
     &node')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxutimen_id, 'standard_name', 
     &       54,'time_of_maximum_wind_velocity_in_x_direction')
              CALL check_err(iret)

             iret = nf_put_att_text(ncid, minvtimen_id, 'long_name', 57, 
     &            'time of minimum wind velocity in y-direction at each 
     &node')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvtimen_id, 'standard_name', 
     &       54,'time_of_minimum_wind_velocity_in_y_direction')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvtimen_id, 'long_name', 57, 
     &            'time of maximum wind velocity in y-direction at each 
     &node')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvtimen_id, 'standard_name', 
     &       54,'time_of_maximum_wind_velocity_in_y_direction')
              CALL check_err(iret)

              iret = nf_put_att_text(ncid, u_id, 'long_name', 28, 
     &                       'wind velocity in x-direction')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, u_id, 'standard_name', 28, 
     &                      'wind_velocity_in_x_direction')
             CALL check_err(iret)
             iret = nf_put_att_text(ncid, u_id, 'positive', 5, VAR_RG)

              iret = nf_put_att_text(ncid, v_id, 'long_name', 28, 
     &                       'wind velocity in y-direction')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, v_id, 'standard_name', 28, 
     &                      'wind_velocity_in_y_direction')
             CALL check_err(iret)
             iret = nf_put_att_text(ncid, v_id, 'positive', 2, VAR_UP)
             CALL check_err(iret)

         endif
      endif

!    -------------------------
!     determine variable units
!    -------------------------
      if (G.LT.11.D0) then
         if (abs(nws).gt.2.and.abs(nws).lt.100) then
            iret = nf_put_att_text(ncid, minut_id,'units',13, VAR_MPS)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, minvt_id,'units',13, VAR_MPS)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, maxut_id,'units',13, VAR_MPS)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, maxvt_id,'units',13, VAR_MPS)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, minun_id,'units',13, VAR_MPS)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, maxun_id,'units',13, VAR_MPS)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, minvn_id,'units',13, VAR_MPS)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, maxvn_id,'units',13, VAR_MPS)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, u_id, 'units',13, VAR_MPS)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, v_id, 'units',13, VAR_MPS)
            CALL check_err(iret)
         else
            iret = nf_put_att_text(ncid, minut_id,'units',19, VAR_MPS_2)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, minvt_id,'units',19, VAR_MPS_2)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, maxut_id,'units',19, VAR_MPS_2)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, maxvt_id,'units',19, VAR_MPS_2)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, minun_id,'units',19, VAR_MPS_2)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, maxun_id,'units',19, VAR_MPS_2)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, minvn_id,'units',19, VAR_MPS_2)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, maxvn_id,'units',19, VAR_MPS_2)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, u_id, 'units',19, VAR_MPS_2)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, v_id, 'units',19, VAR_MPS_2)
            CALL check_err(iret)         
         endif
      else
         if (abs(nws).gt.2.and.abs(nws).lt.100) then
            iret = nf_put_att_text(ncid, minut_id,'units',11, VAR_FPS)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, minvt_id,'units',11, VAR_FPS)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, maxut_id,'units',11, VAR_FPS)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, maxvt_id,'units',11, VAR_FPS)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, minun_id,'units',11, VAR_FPS)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, maxun_id,'units',11, VAR_FPS)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, minvn_id,'units',11, VAR_FPS)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, maxvn_id,'units',11, VAR_FPS)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, u_id, 'units',11, VAR_FPS)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, v_id, 'units',11, VAR_FPS)
            CALL check_err(iret)
         else
            iret = nf_put_att_text(ncid, minut_id,'units',17, VAR_FPS_2)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, minvt_id,'units',17, VAR_FPS_2)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, maxut_id,'units',17, VAR_FPS_2)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, maxvt_id,'units',17, VAR_FPS_2)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, minun_id,'units',17, VAR_FPS_2)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, maxun_id,'units',17, VAR_FPS_2)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, minvn_id,'units',17, VAR_FPS_2)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, maxvn_id,'units',17, VAR_FPS_2)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, u_id, 'units',17, VAR_FPS_2)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, v_id, 'units',17, VAR_FPS_2)
            CALL check_err(iret)
         endif
      endif
      CALL check_err(iret)
      
!    -------------------------------------
!     other variables units attributes
!    ------------------------------------      

      iret = nf_put_att_text(ncid, minunodet_id,'units',14,  VAR_NOND)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, maxunodet_id,'units',14,  VAR_NOND)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, minvnodet_id,'units',14,  VAR_NOND)
      CALL check_err(iret)      
      iret = nf_put_att_text(ncid, maxvnodet_id,'units',14,  VAR_NOND)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, minutimen_id,'units',4, VAR_SECS)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, maxutimen_id,'units',4, VAR_SECS)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, minvtimen_id,'units',4, VAR_SECS)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, maxvtimen_id,'units',4, VAR_SECS)
      CALL check_err(iret)

!    -------------------------------------
!     u and v velocities other attributes
!    ------------------------------------      

      CALL check_err(iret)
      iret = nf_put_att_double(ncid, u_id, 'missing_values', 
     &       NF_DOUBLE, 1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, u_id, '_FillValue', NF_DOUBLE, 1, 
     &                       doubleval)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, u_id, 'dry_Value', NF_DOUBLE, 1, 
     &                       doubleval)
      CALL check_err(iret)
      iret = nf_put_att_DOUBLE(ncid, v_id, 'missing_values', 
     &       NF_DOUBLE, 1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_DOUBLE(ncid, v_id, '_FillValue', NF_DOUBLE, 1, 
     &                       doubleval)
      CALL check_err(iret)
      iret = nf_put_att_DOUBLE(ncid, v_id, 'dry_Value', NF_DOUBLE, 1, 
     &                       doubleval)
      CALL check_err(iret)
    
!     -----------------------------------------------------------------------
!     Define node numbers on elevation specified boundary segments attributes
!     -----------------------------------------------------------------------      
      iret = nf_put_att_text(ncid, nopenc_id, 'long_name', 53, 
     &       'node numbers on elevation specified boundary segments')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, nopenc_id, 'units', 14, 
     &                      VAR_NOND)
      CALL check_err(iret)
      
      iret = nf_put_att_text(ncid, nvdllnc_id, 'long_name', 60, 
     &'number of nodes in each elevation specified boundary segment')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, nvdllnc_id, 'units', 14, 
     &                      VAR_NOND)
      CALL check_err(iret)
!     ---------------------------------------------------------------
!     Define node numbers on normal flow boundary segments attributes
!     ----------------------------------------------------------------  

      iret = nf_put_att_text(ncid, nbounc_id, 'long_name', 59, 
     & 'node numbers in each normal flow specified boundary segment')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, nbounc_id, 'units', 14, 
     &                      VAR_NOND)
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, ibtypenc_id, 'long_name', 13, 
     &       'boundary type')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, ibtypenc_id, 'units', 14, 
     &                      VAR_NOND)
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, nvellnc_id, 'long_name', 62, 
     &'number of nodes in each normal flow specified boundary segment')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, nvellnc_id, 'units', 14, 
     &                      VAR_NOND)
      CALL check_err(iret)

!    -----------------------
!     Define time attributes
!    -----------------------
      iret = nf_put_att_text(ncid,timenc_id,'long_name',10,'model time')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, timenc_id, 'standard_name', 4,'time')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, timenc_id, 'units', 23,VAR_TUN)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, timenc_id, 'base_date',
     . (LEN(TRIM(base_date(1:INDEX(base_date,"!")-1)))), base_date)
      CALL check_err(iret)

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
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'description', 
     &       LEN(TRIM(rundes(1:INDEX(rundes,"!")-1))), rundes)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'grid',
     &       LEN(TRIM(agrid)), agrid)       
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'title',
     &       LEN(TRIM(title(1:INDEX(title,"!")-1))), 
     &       title)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'institution', 
     &       LEN(TRIM(institution(1:INDEX(institution,"!")-1))), 
     &       institution)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'source', 
     &       LEN(TRIM(source(1:INDEX(source,"!")-1))),
     &       source)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'history', 
     &       LEN(TRIM(history(1:INDEX(history,"!")-1))),
     &       history)
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
     &       LEN(TRIM(host(1:INDEX(host,"!")-1))), 
     &       host)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'convention', 
     &       LEN(TRIM(convention(1:INDEX(convention,"!")-1))), 
     &       convention)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'contact', 
     &       LEN(TRIM(contact(1:INDEX(contact,"!")-1))), 
     &       contact)
      CALL check_err(iret)
 
 
      CALL date_and_Time(now_date,big_ben,zone,values)
      write(date_string,71) values(1),values(2),values(3)
     & ,values(5),values(6),values(7),   (values(4))/60
 71   format(I4,'-',I2.2,'-',i2.2,' ',i2,':',i2.2,':',i2.2,' '
     &  ,i3.2,':00')
      iret = nf_put_att_text(ncid, NF_GLOBAL,'creation_date'
     & ,len_trim(date_string),date_string)
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
      iret = nf_put_att_int(ncid, NCGLOBAL, 'irtype', NF_INT, 1, irtype)
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
      intval(1) = 5
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
      iret = nf_put_att_double(ncid, NCGLOBAL, 'ffactor', NF_DOUBLE, 1, 
     &                         ffactor)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, NCGLOBAL, 'esl', NF_DOUBLE, 1, 
     &                         esl)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, NCGLOBAL, 'cori', NF_DOUBLE, 1, 
     &                         cori)
      CALL check_err(iret)
      iret = nf_put_att_int(ncid, NCGLOBAL, 'ntif', NF_INT, 1, ntif)
      CALL check_err(iret)
      iret = nf_put_att_int(ncid, NCGLOBAL, 'nbfr', NF_INT, 1, nbfr)
      CALL check_err(iret)
     
!     ======================
!     Leave define mode
!     ======================
      iret = nf_enddef(ncid)
      CALL check_err(iret)
!     ======================
!     Store static variables
!     ======================
!     -------
!     Store X
!     -------
      iret = nf_put_var_double(ncid, X_id, X)
      CALL check_err(iret)
!     -------
!     Store Y
!     -------
      iret = nf_put_var_double(ncid, Y_id, Y)
      CALL check_err(iret)
!     ------------
!     Store depth
!     ------------
      iret = nf_put_var_double(ncid, DEPTH_id, depth)
      CALL check_err(iret)
!     --------------
!     Store elements
!     --------------
      count(1)=nface_len
      count(2)=elem_len
      start(1)=1
      start(2)=1
      iret=nf_put_vara_int(ncid,ele_id,start,count,element)
      CALL check_err(iret)
!     ------------------------------------
!     Store elevation boundary information
!     -------------------------------------
      iret = nf_put_var_int(ncid, nopenc_id, nopenc)
      CALL check_err(iret)
      iret = nf_put_var_int(ncid, nvdllnc_id, nvdllnc)
      CALL check_err(iret)
!     --------------------------------------
!     Store normal flow boundary information
!     --------------------------------------
      iret = nf_put_var_int(ncid, nbounc_id, nbounc)
      CALL check_err(iret)
      iret = nf_put_var_int(ncid, ibtypenc_id, ibtypenc)
      CALL check_err(iret)
      iret = nf_put_var_int(ncid, nvellnc_id, nvellnc)
      CALL check_err(iret)      
!     --------------------
!     Store maximum depth
!     --------------------
      maxdepth = maxval (depth)
      iret = nf_put_var_double(ncid, maxdepth_id, maxdepth)
      CALL check_err(iret)
!     --------------------
!     Store minimum depth
!     --------------------
      mindepth = minval (depth)
      iret = nf_put_var_double(ncid, mindepth_id, mindepth)
      CALL check_err(iret)

!     ---------------------
!     Open old netcdf file
!     ---------------------
      ELSE
          iret = NF_OPEN(FILENAME, NF_WRITE, ncid)
          CALL check_err(iret)   
!     ------------------
!     Inquire variables
!     ------------------
          iret=nf_inq_unlimdim(ncid,recid)
          CALL check_err(iret)   
          iret=nf_inq_dimlen(ncid,recid,nrec)
          CALL check_err(iret)   
      ENDIF
!     ------------------
!     Update Global Attributes
!     ------------------
      IF(IHOT.NE.0.AND.FIRST74.EQ.1) THEN
         iret = NF_REDEF (ncid)
         CALL check_err(iret)
         iret = nf_put_att_int(ncid, NCGLOBAL, 'ihot', NF_INT, 1, ihot)
         CALL check_err(iret)
         iret = nf_put_att_int(ncid, NCGLOBAL, 'nramp', NF_INT, 1,nramp)
         CALL check_err(iret)
         iret = nf_put_att_double(ncid, NCGLOBAL, 'tau0', NF_DOUBLE, 1, 
     &                         tau0)
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
         iret = nf_put_att_double(ncid, NCGLOBAL, 'ffactor', NF_DOUBLE,  
     &                         1, ffactor)
         CALL check_err(iret)
         iret = nf_put_att_double(ncid, NCGLOBAL, 'esl', NF_DOUBLE, 1, 
     &                         esl)
         CALL check_err(iret)
         CALL date_and_time(now_date,big_ben,zone,values)
         WRITE(date_string,71) values(1),values(2),values(3),
     &         values(5),values(6),values(7),   (values(4))/60
         iret = nf_put_att_text(ncid, NF_GLOBAL,'modification_date',
     &          len_trim(date_string),date_string)
         CALL check_err(iret)
         iret = NF_ENDDEF (ncid)
         CALL check_err(iret)
      ENDIF
      
      IF(ncmode.eq.0) go to 9999
      IF(FIRST74.EQ.1.and.
     & ((ncmode.ne.0.and.IHOT.NE.0).or.(ncmode.eq.2.and.IHOT.eq.0).or.
     &  (ncmode.eq.2.and.IHOT.ne.0.and.nout.eq.-3)))THEN
         iret=nf_inq_varid(ncid,"time",timenc_id)
         iret=nf_inq_varid(ncid,"windx",u_id)
         iret=nf_inq_varid(ncid,"windy",v_id)
         iret=nf_inq_varid(ncid,"minwxt",minut_id)
         iret=nf_inq_varid(ncid,"maxwxt",maxut_id)
         iret=nf_inq_varid(ncid,"minwyt",minvt_id)
         iret=nf_inq_varid(ncid,"maxwyt",maxvt_id)
         iret=nf_inq_varid(ncid,"minwxn",minun_id)
         iret=nf_inq_varid(ncid,"maxwxn",maxun_id)
         iret=nf_inq_varid(ncid,"minwyn",minvn_id)
         iret=nf_inq_varid(ncid,"maxwyn",maxvn_id)
         iret=nf_inq_varid(ncid,"minwxnodet",minunodet_id)
         iret=nf_inq_varid(ncid,"maxwxnodet",maxunodet_id)
         iret=nf_inq_varid(ncid,"minwynodet",minvnodet_id)
         iret=nf_inq_varid(ncid,"maxwynodet",maxvnodet_id)
         iret=nf_inq_varid(ncid,"minwxtimen",minutimen_id)
         iret=nf_inq_varid(ncid,"maxwxtimen",maxutimen_id)
         iret=nf_inq_varid(ncid,"minwytimen",minvtimen_id)
         iret=nf_inq_varid(ncid,"maxwytimen",maxvtimen_id)
      ENDIF
      IF(IHOT.NE.0.AND.FIRST74.EQ.1.and.ncmode.ne.0.and.nout.eq.3)THEN
         iret = nf_get_var_double(ncid, minun_id, minun)
         CALL check_err(iret)
         iret = nf_get_var_double(ncid, maxun_id, maxun)
         CALL check_err(iret)
         iret = nf_get_var_double(ncid, minvn_id, minvn)
         CALL check_err(iret)
         iret = nf_get_var_double(ncid, maxvn_id, maxvn)
         CALL check_err(iret)

         iret = nf_get_var_double(ncid, minutimen_id, minutimen74)
         CALL check_err(iret)
         iret = nf_get_var_double(ncid, maxutimen_id, maxutimen74)
         CALL check_err(iret)
         iret = nf_get_var_double(ncid, minvtimen_id, minvtimen74)
         CALL check_err(iret)
         iret = nf_get_var_double(ncid, maxvtimen_id, maxvtimen74)
         CALL check_err(iret)
      ENDIF
      
      FIRST74=0
      NREC = NREC + 1

!     ----------
!     Store time
!     ----------
      starti(1)=nrec
      counti(1)=timenc_len
      iret = nf_put_vara_DOUBLE(ncid, timenc_id, starti,counti, timenc)
      CALL check_err(iret)
      
!     ---------------------------
!     Store u
!     ---------------------------
      count(1)=node_len
      count(2)=timenc_len
      start(1)=1
      start(2)=nrec
      iret = nf_put_vara_DOUBLE(ncid, u_id, start, count, u)
      CALL check_err(iret)
!     ---------------------------
!     Store v
!     ---------------------------
      count(1)=node_len
      count(2)=timenc_len
      start(1)=1
      start(2)=nrec
      iret = nf_put_vara_DOUBLE(ncid, v_id, start, count, v)
      CALL check_err(iret)
      
!     ----------------------------------------
!     Allocating memory for min and max arrays
!     ----------------------------------------
      IF(.NOT. ALLOCATED(minut))ALLOCATE( minut(timenc_len))
      IF(.NOT. ALLOCATED(maxut))ALLOCATE( maxut(timenc_len))
      IF(.NOT. ALLOCATED(minun))ALLOCATE( minun(node_len  ))
      IF(.NOT. ALLOCATED(maxun))ALLOCATE( maxun(node_len  ))

      IF(.NOT. ALLOCATED(minvt))ALLOCATE( minvt(timenc_len))
      IF(.NOT. ALLOCATED(maxvt))ALLOCATE( maxvt(timenc_len))
      IF(.NOT. ALLOCATED(minvn))ALLOCATE( minvn(node_len  ))
      IF(.NOT. ALLOCATED(maxvn))ALLOCATE( maxvn(node_len  ))

      IF(.NOT. ALLOCATED(minunodet))  ALLOCATE( minunodet(timenc_len))
      IF(.NOT. ALLOCATED(maxunodet))  ALLOCATE( maxunodet(timenc_len))
      IF(.NOT. ALLOCATED(minutimen74))ALLOCATE( minutimen74(node_len))
      IF(.NOT. ALLOCATED(maxutimen74))ALLOCATE( maxutimen74(node_len))

      IF(.NOT. ALLOCATED(minvnodet))  ALLOCATE( minvnodet(timenc_len))
      IF(.NOT. ALLOCATED(maxvnodet))  ALLOCATE( maxvnodet(timenc_len))
      IF(.NOT. ALLOCATED(minvtimen74))ALLOCATE( minvtimen74(node_len))
      IF(.NOT. ALLOCATED(maxvtimen74))ALLOCATE( maxvtimen74(node_len))

!     ---------------------
!     Calculating u minimum
!     ---------------------
!     ----------------------------------
!     Calculating u minimum in each node
!     ----------------------------------
      do n=1,node_len
         if(u (n, 1).lt.minun(n) ) then
           minun (n)       = u (n,1)
           minutimen74 (n) = timenc(1)
         endif
      end do

!     -----------------------------------
!     Calculating u minimum for each time
!     -----------------------------------
      minut (1) = 99999.
      do n=1,node_len
         if(u (n,1).gt.-99999.and.u(n,1).lt. minut (1) ) then
            minut (1)    = u (n,1)
            minunodet(1) = n
         endif
      end do

      iret = nf_put_vara_DOUBLE(ncid, minut_id,(/nrec/),(/1/),minut)
      CALL check_err(iret)
      iret=nf_put_vara_INT(ncid, minunodet_id,(/nrec/),(/1/),minunodet)
      CALL check_err(iret)
!     ----------------------------------
!     Calculating u maximum at each node
!     ----------------------------------
      do n=1,node_len
 	 if(u (n,1).gt.maxun(n) ) then
           maxun (n)       = u (n,1)
           maxutimen74 (n) = timenc(1)
	 endif
      end do
!     -----------------------------------
!     Calculating u maximum for each time
!     ------------------------------------
      maxut (1) = -9999.
      do n=1,node_len
 	 if(u(n,1).gt.maxut(1))  then
            maxut (1)    = u (n,1)
            maxunodet(1) = n 		    
         endif
      end do
      iret = nf_put_vara_DOUBLE(ncid,maxut_id,(/nrec/),(/1/),maxut)
      CALL check_err(iret)
      iret=nf_put_vara_INT(ncid,maxunodet_id,(/nrec/),(/1/),maxunodet)
      CALL check_err(iret)
!     -------------------------------
!     Calculating v minimum each node
!     -------------------------------
      do n=1,node_len
 	 if(v (n,1).lt.minvn(n) ) then
 	    minvn (n)       = v (n,1)
            minvtimen74 (n) = timenc(1)
 	 endif
      end do
!     -----------------------------------
!     Calculating v minimum for each time
!     -----------------------------------
      minvt (1) = 99999.
      do n=1,node_len
         if(v (n,1).gt.-99999.and.v(n,1).lt. minvt (1) )  then
            minvt (1)    = v (n,1)
            minvnodet(1) = n
         endif
      end do

      iret = nf_put_vara_DOUBLE(ncid, minvt_id,(/nrec/),(/1/),minvt)
      CALL check_err(iret)
      iret=nf_put_vara_INT(ncid, minvnodet_id,(/nrec/),(/1/),minvnodet)
      CALL check_err(iret)
!     ---------------------------
!     Calculating v maximum
!     ---------------------------
!     ------------
!     Calculating v maximum in each node
!     ------------
      do n=1,node_len
         if(v (n,1).gt.maxvn(n) ) then
 	    maxvn (n)      = v (n,1)
 	    maxvtimen74(n) = timenc(1)
         endif
      end do
!     --------------
!     Calculating v maximum for each time
!     --------------
      maxvt (1) = -9999.
      do n=1,node_len
 	 if(v(n,1).gt.maxvt(1)) then
 	    maxvt (1)    = v (n,1)
 	    maxvnodet(1) = n 		    
 	 endif
      end do
      iret = nf_put_vara_DOUBLE(ncid,maxvt_id,(/nrec/),(/1/),maxvt)
      CALL check_err(iret)
      iret=nf_put_vara_INT(ncid,maxvnodet_id,(/nrec/),(/1/),maxvnodet)
      CALL check_err(iret)

!     ---------------------------
!     Closing fort.X4 netCDF file
!     ---------------------------
      IF(NREC.EQ.NDATA) THEN
         iret = nf_put_var_double(ncid, minun_id, minun)
         CALL check_err(iret)
         iret = nf_put_var_double(ncid, maxun_id, maxun)
         CALL check_err(iret)
         iret = nf_put_var_double(ncid, minvn_id, minvn)
         CALL check_err(iret)
         iret = nf_put_var_double(ncid, maxvn_id, maxvn)
         CALL check_err(iret)

         iret = nf_put_var_double(ncid, minutimen_id, minutimen74)
         CALL check_err(iret)
         iret = nf_put_var_double(ncid, maxutimen_id, maxutimen74)
         CALL check_err(iret)
         iret = nf_put_var_double(ncid, minvtimen_id, minvtimen74)
         CALL check_err(iret)
         iret = nf_put_var_double(ncid, maxvtimen_id, maxvtimen74)
         CALL check_err(iret)

      ENDIF

!     ---------------------------
!     Closing netCDF file
!     ---------------------------
      WRITE(6,*) "Closing ", FILENAME, " netCDF file ... "
      WRITE(6,*) "-----------------------------------------------------"
      
9999  iret = nf_close(ncid)
      CALL check_err(iret)
      WRITE(6,*)'end subroutine writenc74'

      end subroutine writenc74
      
      
!=======================================================================
!     Subroutine writenc61
!     Writes a netCDF file equivalent to ADCIRC's fort.61 ascii file
!
!     Dimensions:
!     station  : number of stations
!     namelen  : length of station's name
!     time     : number of time records
!
!     Variables :
!     time      : time
!     station_name : station's name
!     zeta      : Sea surface elevation 
!     x         : x-coordinate or longitude
!     y         : y-coordinate or latitude
!     minzetat  : minimum sea surface elevation for each output time
!     maxzetat  : maximum sea surface elevation for each output time
!     minzetan  : minimum sea surface elevation at each station
!     maxzetan  : maximum sea surface elevation at each station
!     minznodet : node of minimum sea surface elevation at each output time
!     maxznodet : node of maximum sea surface elevation at each output time
!     minztimen : time of minimum sea surface elevation at each station
!     maxztimen : time of maximum sea surface elevation at each station
!
!     Attributes:
!     model             : ADCIRC
!     version           : model version number
!     grid_type         : triangular
!     rundes            : description of simulation
!     agrid             : grid ID 
!     title             : description of the dataset
!     institution       : where the original data was produced
!     source            : method of production of the original data
!     history           : provides an audit trail for modifications to the original data
!     references        : published/web-based references that describe the data 
!     comments          : miscellaneous information about the data
!     contact           : user's institution/ e-mail address
!     creation_date     : date the data was created
!     modification_date : date the data was modified
!     Parameters from fort.15 input file 
!
!     Revision history:
!
!     Date            Programmer                  Description of change
!     ----            ----------                  ----------------------
!     03/30/07  Cristina Forbes, PSGS @ UNC-IMS   Wrote original code
!     03/30/08  Cristina Forbes, PSGS @ UNC-IMS   Modified code for 
!                                                 globalio & hotstart 
!                                                 from binary files
!=======================================================================
      subroutine writenc61 (nout,fileunit,timenc_len, node_len,timenc,
     &                      zeta,nt, NDATA, rundes,runid, agrid, irtype,
     &                      dtdp, ihot, ics, nolibf,  nolifa,
     &                      nolica, nolicat, nwp,
     &                      ncor, ntip, nws, nramp, tau0,
     &                      statim, reftim, rnday, dramp,
     &                      a00, b00, c00, h0,
     &                      slam0, sfea0, ffactor,esl, cori,
     &                      ntif, nbfr, elem_len, nface_len, nm,
     &                      x,y,depth,
     &                      nope_len, neta_len, nbdv, 
     &                      nbou_len, nvel_len, nbvv,
     &                      ncmode,
     &                      title, institution, source, history,
     &                      references, comments, host, convention,
     &                      contact,
     &                      base_date)

!--------------------------------------------
!     Force explicit declaration of all variables
!--------------------------------------------
      IMPLICIT NONE
!     -------------------
!     netCDF header file
!     -------------------
!      include 'netcdf.inc'
!     -------------------
!     Error status return
!     -------------------
      INTEGER  iret
!     ---------
!     netCDF id
!     ---------
      INTEGER  ncid
!     -------------
!     Dimension ids
!     -------------
      INTEGER  node_dim
      INTEGER  timenc_dim
      INTEGER  elem_dim
      INTEGER  nface_dim
      INTEGER  nope_dim
      INTEGER  neta_dim
      INTEGER  nbou_dim
      INTEGER  nvel_dim
      INTEGER  slen_dim
      INTEGER  scode_dim
      
!     -----------------
!     Dimension lengths
!     -----------------
      INTEGER  node_len
      INTEGER  timenc_len
      INTEGER  elem_len
      INTEGER  nface_len
      INTEGER  nope_len
      INTEGER  neta_len
      INTEGER  nbou_len
      INTEGER  nvel_len
!     ------------
!     Variable ids
!     ------------
      INTEGER  timenc_id
      INTEGER  zeta_id
      INTEGER  minzetat_id
      INTEGER  maxzetat_id
      INTEGER  minzetan_id
      INTEGER  maxzetan_id
      INTEGER  minznodet_id
      INTEGER  maxznodet_id
      INTEGER  minztimen_id
      INTEGER  maxztimen_id
      INTEGER  station_id
      INTEGER  scode_id
      INTEGER  X_id
      INTEGER  Y_id
      INTEGER  DEPTH_id
      INTEGER  ELE_id
      INTEGER  mindepth_id
      INTEGER  maxdepth_id
      INTEGER  nbdv_id
      INTEGER  nbvv_id
!     ---------------------------------------------
!     Rank (number of dimensions) for each variable
!     ---------------------------------------------
      INTEGER  timenc_rank
      INTEGER  zeta_rank
      INTEGER  station_rank
      INTEGER  scode_rank
      INTEGER  minzetat_rank
      INTEGER  maxzetat_rank
      INTEGER  minzetan_rank
      INTEGER  maxzetan_rank
      INTEGER  minznodet_rank
      INTEGER  maxznodet_rank
      INTEGER  minztimen_rank
      INTEGER  maxztimen_rank
      INTEGER  depth_rank
      INTEGER  ele_rank
      INTEGER  mindepth_rank
      INTEGER  maxdepth_rank
      INTEGER  x_rank
      INTEGER  y_rank
      INTEGER  nvdll_rank
      INTEGER  nbdv_rank
      INTEGER  nbvv_rank
      PARAMETER (timenc_rank = 1)
      PARAMETER (station_rank = 2)
      PARAMETER (scode_rank = 2)
      PARAMETER (zeta_rank = 2)
      PARAMETER (minzetat_rank = 1)
      PARAMETER (maxzetat_rank = 1)
      PARAMETER (minzetan_rank = 1)
      PARAMETER (maxzetan_rank = 1)
      PARAMETER (minznodet_rank = 1)
      PARAMETER (maxznodet_rank = 1)
      PARAMETER (minztimen_rank = 1)
      PARAMETER (maxztimen_rank = 1)
      PARAMETER (depth_rank = 1)
      PARAMETER (ele_rank = 2)
      PARAMETER (mindepth_rank = 0)
      PARAMETER (maxdepth_rank = 0)
      PARAMETER (x_rank = 1)
      PARAMETER (y_rank = 1)
      PARAMETER (nvdll_rank = 1)
      PARAMETER (nbdv_rank = 2)
      PARAMETER (nbvv_rank = 2)
!     ---------------
!     Variable shapes
!     ---------------
      INTEGER  timenc_dims   (timenc_rank)
      INTEGER  zeta_dims     (zeta_rank)
      INTEGER  station_dims  (station_rank)
      INTEGER  scode_dims    (scode_rank)
      INTEGER  minzetat_dims (minzetat_rank)
      INTEGER  maxzetat_dims (maxzetat_rank)
      INTEGER  minzetan_dims (minzetan_rank)
      INTEGER  maxzetan_dims (maxzetan_rank)
      INTEGER  minznodet_dims(minznodet_rank)
      INTEGER  maxznodet_dims(maxznodet_rank)
      INTEGER  minztimen_dims(minztimen_rank)
      INTEGER  maxztimen_dims(maxztimen_rank)
      INTEGER  x_dims        (x_rank)
      INTEGER  y_dims        (y_rank)
      INTEGER  depth_dims    (depth_rank)
      INTEGER  ele_dims      (ele_rank)
      INTEGER  nvdll_dims    (nvdll_rank)
      INTEGER  nbdv_dims     (nbdv_rank)
      INTEGER  nbvv_dims     (nbvv_rank)
!     --------------
!     Data variables
!     --------------
      CHARACTER*5 FILEFORT
      CHARACTER*10 FILENAME
      CHARACTER*3 FILEEXT
      INTEGER   t,ntif,nbfr,nt
      INTEGER   irtype, ihot, ics, nolibf, nolifa, nolica, nolicat, nwp
      INTEGER   ntip, nws, nramp,ncor
      INTEGER   nvel, neta, nout
      INTEGER   element (nface_len , elem_len)
      INTEGER   nbdv (nope_len, neta_len)
      INTEGER   nbvv (nbou_len, nvel_len)
      INTEGER   nm (elem_len, nface_len)
      INTEGER   fileunit
      REAL*8    timenc(timenc_len)
      REAL(sz)  x(node_len)
      REAL(sz)  y(node_len)
      REAL(sz)  depth(node_len)
      REAL(sz)  zeta (node_len, timenc_len)
      REAL(sz), ALLOCATABLE, SAVE ::  minzetat(:)
      REAL(sz), ALLOCATABLE, SAVE ::  maxzetat(:)
      REAL(sz), ALLOCATABLE, SAVE ::  minzetan61(:)
      REAL(sz), ALLOCATABLE, SAVE ::  maxzetan61(:)
      INTEGER , ALLOCATABLE, SAVE ::  minznodet(:)
      INTEGER , ALLOCATABLE, SAVE ::  maxznodet(:)
      REAL(sz), ALLOCATABLE, SAVE ::  minztimen61(:)
      REAL(sz), ALLOCATABLE, SAVE ::  maxztimen61(:)
      real*8    mindepth
      real*8    maxdepth
      REAL*8    dtdp, statim,reftim, rnday, dramp, tau0
      REAL*8    a00, b00, c00,h0
      REAL*8    ffactor,esl, cori
      REAL*8    nele
      REAL(sz)  slam0, sfea0
      REAL(sz)  slam0deg, sfea0deg
      CHARACTER(80) rundes, runid, agrid
      INTEGER, SAVE  :: FIRST61 = 1 
      INTEGER NDATA
      INTEGER, SAVE  :: NREC = 0
      INTEGER ncmode
!     -----------------
!     Attribute vectors
!     -----------------
      INTEGER  intval(1),i,j
      INTEGER  counti(1), starti(1), n
      INTEGER  count(2), start(2)
      REAL     realval(1)
!      DOUBLE PRECISION  doubleval(1)
!     -----------------
!     Attribute general
!     -----------------
      CHARACTER(*) :: title, institution, source, history
      CHARACTER(*) :: references, comments, host, convention, contact
      CHARACTER(*) :: base_date
!     -----------------
!     date_string variables for time attribute
!     -----------------
      character date_string*40
      character now_date*8
      character big_ben*10
      character zone*5
      integer values(8)
      real*8 xmin, xmax
!     ---------------------------------
!      inquiring and reading variables
!     ---------------------------------

      INTEGER NDIMS, NVARS, NGATTS, UNLIMDIMID, RECID,NRECS
      CHARACTER*4 recname

!     -------------------------
!     Save dims, ids, variables
!     -------------------------
      save node_dim, elem_dim, nface_dim
      save nope_dim, neta_dim, nbou_dim, nvel_dim
      save zeta_id, X_id, Y_id
      save minzetat_id, maxzetat_id, minzetan_id, maxzetan_id
      save minznodet_id, maxznodet_id, minztimen_id, maxztimen_id
      save DEPTH_id,ELE_id, mindepth_id, maxdepth_id, nbdv_id, nbvv_id
      save timenc_dim, timenc_id
      save station_id, scode_id
      save mindepth, maxdepth
!     ---------------------------
!     Opening fort.61 netCDF file
!     ---------------------------
      WRITE(6,*) "-----------------------------------------------------"

      IF(FIRST61.EQ.1.and.ncmode.NE.0) THEN
         
         ALLOCATE ( minzetan61  (node_len) )
         ALLOCATE ( maxzetan61  (node_len) )
         ALLOCATE ( minztimen61 (node_len) )
         ALLOCATE ( maxztimen61 (node_len) )
         
         DO n=1,node_len
            minzetan61  (n) =  9999.d0
            maxzetan61  (n) = -9999.d0
            minztimen61 (n) =  9999.d0
            maxztimen61 (n) = -9999.d0
         END DO         
      ENDIF
      
!     --------------------------------
!     Switch order in array for NETCDF
!     --------------------------------
      do i=1, elem_len
         do j=1, nface_len
            element(j,i) = nm (i,j)
         end do
      end do
!     -----------------
!     Enter define mode
!     -----------------
         FILEFORT = 'fort.'
         FILEEXT = '.nc'
!        WRITE(6,'(A5,I2,A3)')FILEFORT,FILEUNIT,FILEEXT
         WRITE(FILENAME,'(A5,I2,A3)')FILEFORT,FILEUNIT,FILEEXT
         WRITE(6,*) "Opening ",FILENAME,  " netCDF file ... "
                  
       IF((FIRST61.EQ.1.AND.IHOT.EQ.0.AND.NCMODE.NE.2).OR.
     &    (FIRST61.EQ.1.AND.IHOT.NE.0.AND.NCMODE.NE.2.AND.
     &     nout.EQ.-3))THEN

           iret = nf_create(FILENAME, NF_CLOBBER, ncid)
           CALL check_err(iret)
           FIRST61=0
         
!     --------------------------------
!     Convert to degress for NETCDF I/O
!     --------------------------------
           SLAM0DEG=SLAM0
           SFEA0DEG=SFEA0
           IF (NCMODE.EQ.1) THEN
               SLAM0DEG=SLAM0*RAD2DEG
               SFEA0DEG=SFEA0*RAD2DEG
           ENDIF

!     =================
!     Define dimensions
!     =================
      iret = nf_def_dim(ncid, 'station', node_len, node_dim)
             CALL check_err(iret)
      iret = nf_def_dim(ncid, 'namelen', SNAMLEN, slen_dim)
             CALL check_err(iret)
      iret = nf_def_dim(ncid, 'time', nf_unlimited, timenc_dim)
             CALL check_err(iret)
!     ================
!     Define variables
!     ================
!     ----------------
!     Define time
!     ----------------
      timenc_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'time',NF_DOUBLE,timenc_rank, timenc_dims, 
     &                  timenc_id)
      CALL check_err(iret)
!     -------------------
!     Define X
!     -------------------
      X_dims(1) = node_dim
      iret = nf_def_var(ncid, 'x', NF_DOUBLE, X_rank, X_dims, X_id)
      CALL check_err(iret)
!     -------------------
!     Define Y coordinate
!     --------------------
      Y_dims(1) = node_dim
      iret = nf_def_var(ncid, 'y', NF_DOUBLE, Y_rank, Y_dims, Y_id)
      CALL check_err(iret)
!     ---------------------
!     Define stations name
!     --------------------
      station_dims(1) = slen_dim
      station_dims(2) = node_dim
      iret = nf_def_var(ncid, 'station_name', NF_CHAR, station_rank, 
     1      station_dims, station_id)
      CALL check_err(iret)
      
!     ----------------
!     Define zeta
!     ----------------
      zeta_dims(1) = node_dim
      zeta_dims(2) = timenc_dim
      
      iret = nf_def_var(ncid, 'zeta', NF_DOUBLE, zeta_rank, zeta_dims, 
     &                  ZETA_id)
      
      CALL check_err(iret)
!     ----------------------------
!     Define minimum zeta at times
!     ----------------------------
      minzetat_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'minzetat', NF_DOUBLE, minzetat_rank, 
     &                  minzetat_dims, minzetat_id)
      CALL check_err(iret)
!     ----------------------------
!     Define maximum zeta at times
!     ----------------------------
      maxzetat_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'maxzetat', NF_DOUBLE, maxzetat_rank, 
     &                  maxzetat_dims, maxzetat_id)
      CALL check_err(iret)
!     ----------------------------
!     Define minimum zeta at nodes
!     ----------------------------
      minzetan_dims(1) = node_dim
      iret = nf_def_var(ncid, 'minzetan', NF_DOUBLE, minzetan_rank, 
     &                  minzetan_dims, minzetan_id)
      CALL check_err(iret)
!     ----------------------------
!     Define maximum zeta at nodes
!     ----------------------------
      maxzetan_dims(1) = node_dim
      iret = nf_def_var(ncid, 'maxzetan', NF_DOUBLE, maxzetan_rank, 
     &                  maxzetan_dims, maxzetan_id)
      CALL check_err(iret)
!     ----------------------------
!     Define minimum zeta nodes at times
!     ----------------------------
      minznodet_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'minznodet', NF_INT, minznodet_rank, 
     &                  minznodet_dims, minznodet_id)
      CALL check_err(iret)
!     ----------------------------
!     Define maximum zeta nodes at times
!     ----------------------------
      maxznodet_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'maxznodet', NF_INT, maxznodet_rank, 
     &                  maxznodet_dims, maxznodet_id)
      CALL check_err(iret)
!     ----------------------------
!     Define minimum zeta nodes at nodes
!     ----------------------------
      minztimen_dims(1) = node_dim
      iret = nf_def_var(ncid, 'minztimen', NF_DOUBLE, minztimen_rank, 
     &                  minztimen_dims, minztimen_id)
      CALL check_err(iret)
!     ----------------------------
!     Define maximum zeta nodes at nodes
!     ----------------------------
      maxztimen_dims(1) = node_dim
      iret = nf_def_var(ncid, 'maxztimen', NF_DOUBLE, maxztimen_rank, 
     &                  maxztimen_dims, maxztimen_id)
      CALL check_err(iret)

!     ============================
!     Assign attributes
!     ============================
!     ----------------------------
!     Define longitude attributes
!     ----------------------------
      if(ics.eq.2) then
      
      	iret = nf_put_att_text(ncid,X_id,'long_name',9,'longitude')
      	       CALL check_err(iret)
      	iret = nf_put_att_text(ncid, X_id, 'standard_name', 9, 
     &                      'longitude')
      	       CALL check_err(iret)
      	iret = nf_put_att_text(ncid, X_id, 'units', 12,   VAR_DEGE)
      	       CALL check_err(iret)
      	iret = nf_put_att_text(ncid, X_id, 'positive', 4, VAR_EAST)
      	       CALL check_err(iret)

      	iret = nf_put_att_double(ncid, X_id, 'missing_values',NF_DOUBLE, 
     &                       1, doubleval)
      	       CALL check_err(iret)

      	iret = nf_put_att_double(ncid, X_id, '_FillValue', NF_DOUBLE, 1, 
     &                       doubleval)
      	       CALL check_err(iret)      
      
!     ----------------------------
!     Define latitude attributes
!     ----------------------------
      	iret = nf_put_att_text(ncid, Y_id, 'long_name', 8, 'latitude')
      	       CALL check_err(iret)
      	iret = nf_put_att_text(ncid, Y_id, 'standard_name', 8, 
     &                       'latitude')
      	       CALL check_err(iret)
      	iret = nf_put_att_text(ncid, Y_id, 'units', 13,  VAR_DEGN)
      	       CALL check_err(iret)
      	iret = nf_put_att_text(ncid, Y_id, 'positive', 5, VAR_NORTH)
      	       CALL check_err(iret)
       
 
      	iret = nf_put_att_double(ncid, Y_id, 'missing_values', 
     &         NF_DOUBLE, 1, doubleval)
      	       CALL check_err(iret)
 
      	iret = nf_put_att_double(ncid, Y_id, '_FillValue', NF_DOUBLE, 1, 
     &         doubleval)
      	       CALL check_err(iret)

       else
      
!     ------------------------------
!     Define x-coordinate attributes
!     ------------------------------
      	iret = nf_put_att_text(ncid, X_id, 'long_name', 22, 
     &       'Cartesian coordinate x')
      	       CALL check_err(iret)
      	iret = nf_put_att_text(ncid, X_id, 'standard_name',12, 
     &                      'x-coordinate')
      	       CALL check_err(iret)
      	if (G.LT.11.D0) then
               iret = nf_put_att_text(ncid, X_id, 'units', 6, VAR_METS)  
      	       CALL check_err(iret)
        else
      	       iret = nf_put_att_text(ncid, X_id, 'units', 4, VAR_FEET)  
      	       CALL check_err(iret)
        endif


      	iret = nf_put_att_text(ncid, X_id, 'positive', 5, VAR_RG)
      	       CALL check_err(iret)

      	iret = nf_put_att_double(ncid, X_id, 'missing_values',NF_DOUBLE,
     &                        1, doubleval)
      	       CALL check_err(iret)

      	iret = nf_put_att_double(ncid, X_id, '_FillValue', NF_DOUBLE, 1, 
     &                       doubleval)
      	       CALL check_err(iret)
!     ------------------------------
!     Define y-coordinate attributes
!     ------------------------------     
      	iret = nf_put_att_text(ncid, Y_id, 'long_name', 22, 
     &       'Cartesian coordinate y')
      	       CALL check_err(iret)
      	iret = nf_put_att_text(ncid, Y_id, 'standard_name', 12, 
     &                       'y-coordinate')
      	       CALL check_err(iret)
      	if (G.LT.11.D0) then
     	    iret = nf_put_att_text(ncid, Y_id, 'units', 6, VAR_METS)
      	    CALL check_err(iret)
      	else
     	    iret = nf_put_att_text(ncid, Y_id, 'units', 4, VAR_FEET)
      	    CALL check_err(iret)
      	endif
      	iret = nf_put_att_text(ncid, Y_id, 'positive', 5, VAR_UP)
      	       CALL check_err(iret)
 
      	iret = nf_put_att_double(ncid, Y_id, 'missing_values',NF_DOUBLE, 
     &       1, doubleval)
      	       CALL check_err(iret)
 
      	iret = nf_put_att_double(ncid, Y_id, '_FillValue', NF_DOUBLE, 1, 
     &       doubleval)
      	       CALL check_err(iret)
      
      endif
!     -------------------------------------------------------------
!     Define max and min in time sea surface elevation attributes 
!     -------------------------------------------------------------
      iret = nf_put_att_text(ncid, minzetat_id, 'long_name', 58, 
     &     'minimum station sea surface elevation for each output time')
             CALL check_err(iret)
      if (G.LT.11.D0) then
         iret = nf_put_att_text(ncid, minzetat_id, 'units', 6, VAR_METS)
         CALL check_err(iret)
      else
         iret = nf_put_att_text(ncid, minzetat_id, 'units', 4, VAR_FEET)
         CALL check_err(iret)
      endif
      
      iret = nf_put_att_text(ncid, maxzetat_id, 'long_name', 58, 
     &     'maximum station sea surface elevation for each output time')
      	     CALL check_err(iret)
      if (G.LT.11.D0) then
         iret = nf_put_att_text(ncid, maxzetat_id, 'units', 6, VAR_METS)
         CALL check_err(iret)
      else
         iret = nf_put_att_text(ncid, maxzetat_id, 'units', 4, VAR_FEET)
         CALL check_err(iret)
      endif
!     -------------------------------------------------------------
!     Define max and min in nodes sea surface elevation attributes
!     -------------------------------------------------------------
      iret = nf_put_att_text(ncid, minzetan_id, 'long_name', 45, 
     &      'minimum sea surface elevation at each station')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, maxzetan_id, 'long_name', 45, 
     &       'maximum sea surface elevation at each station')
      CALL check_err(iret)
      if (G.LT.11.D0) then
         iret = nf_put_att_text(ncid, minzetan_id, 'units', 6, VAR_METS)
         CALL check_err(iret)
         iret = nf_put_att_text(ncid, maxzetan_id, 'units', 6, VAR_METS)
         CALL check_err(iret)
      else
         iret = nf_put_att_text(ncid, minzetan_id, 'units', 4, VAR_FEET)
         CALL check_err(iret)
         iret = nf_put_att_text(ncid, maxzetan_id, 'units', 6, VAR_FEET)
         CALL check_err(iret)
      endif

!     -------------------------------------------------------------
!     Define max and min in time sea surface elevation attributes 
!     -------------------------------------------------------------
      iret = nf_put_att_text(ncid, minznodet_id, 'long_name', 79, 
     &       'station where minimum sea surface elevation occurs at each 
     &output time')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, minznodet_id, 'units',14, VAR_NOND)
      CALL check_err(iret)
      
      iret = nf_put_att_text(ncid, maxznodet_id, 'long_name', 79, 
     &       'station where maximum sea surface elevation occurs at each 
     &output time')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, maxznodet_id, 'units', 14, VAR_NOND)
      CALL check_err(iret)
!     -------------------------------------------------------------
!     Define max and min in nodes sea surface elevation attributes
!     -------------------------------------------------------------
      iret = nf_put_att_text(ncid, minztimen_id, 'long_name', 62, 
     & 'time when minimum sea surface elevation occurs at each station')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, minztimen_id, 'units', 4, VAR_SECS)
      CALL check_err(iret)
      
      iret = nf_put_att_text(ncid, maxztimen_id, 'long_name', 62, 
     & 'time when maximum sea surface elevation occurs at each station')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, maxztimen_id, 'units', 4, VAR_SECS)
      CALL check_err(iret)

!     -------------------------------------------------------------
!     Define station names and code attributes
!     -------------------------------------------------------------

      iret = nf_put_att_text(ncid,station_id,'long_name',12,
     & 'station name')
      CALL check_err(iret)
      
!    -----------------------
!     Define time attributes
!    -----------------------
      iret = nf_put_att_text(ncid,timenc_id,'long_name',10,'model time')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid,timenc_id, 'standard_name', 4, 'time')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, timenc_id, 'units', 23, 
     .       'seconds since base_date')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, timenc_id, 'base_date',
     .       (LEN(TRIM(base_date(1:INDEX(base_date,"!")-1)))), base_date)
      CALL check_err(iret)

!    --------------------------------------------
!     Define sea surface elevation attributes
!    --------------------------------------------
      
      iret = nf_put_att_text(ncid, zeta_id,'long_name',29, VAR61_LN)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, zeta_id, 'standard_name', 25, 
     &                      VAR61_ST)
      CALL check_err(iret)
      if (G.LT.11.D0) then
          iret = nf_put_att_text(ncid, zeta_id, 'units', 6, VAR_METS)
          CALL check_err(iret)
      else
          iret = nf_put_att_text(ncid, zeta_id, 'units', 6, VAR_FEET)
          CALL check_err(iret)
      endif
      iret = nf_put_att_DOUBLE(ncid,zeta_id,'missing_values',NF_DOUBLE,
     &1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_DOUBLE(ncid, zeta_id, '_FillValue', NF_DOUBLE, 
     &                      1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, zeta_id,'positive', 2, VAR_UP)
      CALL check_err(iret)
 
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
      intval(1) = 1
      iret = nf_put_att_int(ncid, NCGLOBAL, 'irtype', NF_INT, 1, irtype)
      CALL check_err(iret)
      intval(1) = 0
      iret = nf_put_att_int(ncid, NCGLOBAL, 'ihot', NF_INT, 1, ihot)
      CALL check_err(iret)
      intval(1) = 1
      iret = nf_put_att_int(ncid, NCGLOBAL, 'ics', NF_INT, 1, ics)
      CALL check_err(iret)
      intval(1) = 1
      iret = nf_put_att_int(ncid, NCGLOBAL, 'nolibf', NF_INT, 1, nolibf)
      CALL check_err(iret)
      intval(1) = 1
      iret = nf_put_att_int(ncid, NCGLOBAL, 'nolifa', NF_INT, 1, nolifa)
      CALL check_err(iret)
      intval(1) = 1
      iret = nf_put_att_int(ncid, NCGLOBAL, 'nolica', NF_INT, 1, nolica)
      CALL check_err(iret)
      intval(1) = 1
      iret = nf_put_att_int(ncid, NCGLOBAL, 'nolicat', NF_INT, 1, 
     &                      nolicat)
      CALL check_err(iret)
      intval(1) = 0
      iret = nf_put_att_int(ncid, NCGLOBAL, 'nwp', NF_INT, 1, nwp)
      CALL check_err(iret)
      intval(1) = 0
      iret = nf_put_att_int(ncid, NCGLOBAL, 'ncor', NF_INT, 1, ncor)
      CALL check_err(iret)
      intval(1) = 0
      iret = nf_put_att_int(ncid, NCGLOBAL, 'ntip', NF_INT, 1, ntip)
      CALL check_err(iret)
      intval(1) = 0
      iret = nf_put_att_int(ncid, NCGLOBAL, 'nws', NF_INT, 1, nws)
      CALL check_err(iret)
      intval(1) = 1
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
      intval(1) = 5
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
      iret = nf_put_att_double(ncid, NCGLOBAL, 'ffactor', NF_DOUBLE, 1, 
     &                         ffactor)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, NCGLOBAL, 'esl', NF_DOUBLE, 1, 
     &                         esl)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, NCGLOBAL, 'cori', NF_DOUBLE, 1, 
     &                         cori)
      CALL check_err(iret)
      intval(1) = 0
      iret = nf_put_att_int(ncid, NCGLOBAL, 'ntif', NF_INT, 1, ntif)
      CALL check_err(iret)
      intval(1) = 1
      iret = nf_put_att_int(ncid, NCGLOBAL, 'nbfr', NF_INT, 1, nbfr)
      CALL check_err(iret)
!     ======================
!     Leave define mode
!     ======================
      iret = nf_enddef(ncid)
      CALL check_err(iret)
!     ======================
!     Store static variables
!     ======================
!     -------
!     Store X
!     -------
      iret = nf_put_var_double(ncid, X_id, X)
      CALL check_err(iret)
!     -------
!     Store Y
!     -------
      iret = nf_put_var_double(ncid, Y_id, Y)
      CALL check_err(iret)
!     --------------------
!     Store station name
!     --------------------
 
      do i=1,node_len
         start(1)=1
         start(2)=i
         count(1)= LEN(TRIM(STATNAME(i,SNAMLEN)))
         count(2)=1
         iret = nf_put_vara_text(ncid,station_id,start,count,
     &   STATNAME(i,SNAMLEN) )
         CALL check_err(iret)
      end do
      
!     =====================
!     Open old netcdf file
!     =====================
      ELSE
          iret = NF_OPEN(FILENAME, NF_WRITE, ncid)
          CALL check_err(iret)   
!     ------------------
!     Inquire variables
!     ------------------
          iret=nf_inq_unlimdim(ncid,recid)
          CALL check_err(iret)   
          iret=nf_inq_dimlen(ncid,recid,nrec)
          CALL check_err(iret)   
      ENDIF
!     ------------------
!     Update Global Attributes
!     ------------------
      IF(IHOT.NE.0.AND.FIRST61.EQ.1) THEN
         iret = NF_REDEF (ncid)
         CALL check_err(iret)
         iret = nf_put_att_int(ncid, NCGLOBAL, 'ihot', NF_INT, 1, ihot)
         CALL check_err(iret)
         iret = nf_put_att_int(ncid, NCGLOBAL, 'nramp', NF_INT, 1,nramp)
         CALL check_err(iret)
         iret = nf_put_att_double(ncid, NCGLOBAL, 'tau0', NF_DOUBLE, 1, 
     &                         tau0)
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
         iret = nf_put_att_double(ncid, NCGLOBAL, 'ffactor', NF_DOUBLE,
     &                          1, ffactor)
         CALL check_err(iret)
         iret = nf_put_att_double(ncid, NCGLOBAL, 'esl', NF_DOUBLE, 1, 
     &                         esl)
         CALL check_err(iret)
         CALL date_and_time(now_date,big_ben,zone,values)
         WRITE(date_string,71) values(1),values(2),values(3),
     &         values(5),values(6),values(7),   (values(4))/60
         iret = nf_put_att_text(ncid, NF_GLOBAL,'modification_date',
     &          len_trim(date_string),date_string)
         CALL check_err(iret)
         iret = NF_ENDDEF (ncid)
         CALL check_err(iret)
       ENDIF
      IF(ncmode.eq.0) go to 9999
      IF(FIRST61.EQ.1.and.
     & ((ncmode.ne.0.and.IHOT.NE.0).or.(ncmode.eq.2.and.IHOT.eq.0).or.
     &  (ncmode.eq.2.and.IHOT.ne.0.and.nout.eq.-3)))THEN
         iret=nf_inq_varid(ncid,"time",timenc_id)
         iret=nf_inq_varid(ncid,"zeta",zeta_id)
         iret=nf_inq_varid(ncid,"minzetat",minzetat_id)
         iret=nf_inq_varid(ncid,"maxzetat",maxzetat_id)
         iret=nf_inq_varid(ncid,"minzetan",minzetan_id)
         iret=nf_inq_varid(ncid,"maxzetan",maxzetan_id)
         iret=nf_inq_varid(ncid,"minznodet",minznodet_id)
         iret=nf_inq_varid(ncid,"maxznodet",maxznodet_id)
         iret=nf_inq_varid(ncid,"minztimen",minztimen_id)
         iret=nf_inq_varid(ncid,"maxztimen",maxztimen_id)
      ENDIF
      IF(IHOT.NE.0.AND.FIRST61.EQ.1.and.ncmode.ne.0.and.nout.eq.3)THEN
         iret = nf_get_var_double(ncid, minzetan_id, minzetan61)
         CALL check_err(iret)
         iret = nf_get_var_double(ncid, maxzetan_id, maxzetan61)
         CALL check_err(iret)
         iret = nf_get_var_double(ncid, minztimen_id, minztimen61)
         CALL check_err(iret)
         iret = nf_get_var_double(ncid, maxztimen_id, maxztimen61)
         CALL check_err(iret)
      ENDIF
      
      FIRST61=0
      NREC = NREC + 1

!     ----------
!     Store time
!     ----------
      starti(1)=nrec
      counti(1)=timenc_len
      iret = nf_put_vara_DOUBLE(ncid, timenc_id, starti,counti, timenc)
      CALL check_err(iret)
!     ---------------------------
!     Store zeta
!     ---------------------------
      count(1)=node_len
      count(2)=timenc_len
      start(1)=1
      start(2)=nrec
      iret = nf_put_vara_DOUBLE(ncid, zeta_id, start, count, zeta)
      CALL check_err(iret)
!     ----------------------------------------
!     Allocating memory for min and max arrays
!     ----------------------------------------
      IF(.NOT. ALLOCATED(minzetat))  ALLOCATE( minzetat(timenc_len))
      IF(.NOT. ALLOCATED(maxzetat))  ALLOCATE( maxzetat(timenc_len))
      IF(.NOT. ALLOCATED(minzetan61))ALLOCATE( minzetan61(node_len))
      IF(.NOT. ALLOCATED(maxzetan61))ALLOCATE( maxzetan61(node_len))

      IF(.NOT. ALLOCATED(minznodet))  ALLOCATE( minznodet(timenc_len))
      IF(.NOT. ALLOCATED(maxznodet))  ALLOCATE( maxznodet(timenc_len))
      IF(.NOT. ALLOCATED(minztimen61))ALLOCATE( minztimen61(node_len))
      IF(.NOT. ALLOCATED(maxztimen61))ALLOCATE( maxztimen61(node_len))

!     --------------------------------------
!     Calculating zeta minimum  in each node
!     --------------------------------------
      do n=1,node_len
         if(zeta (n,1).lt.minzetan61(n) )then
            minzetan61  (n) = zeta (n,1)
 	    minztimen61 (n) = timenc(1)
         endif
      end do
!     --------------------------------------
!     Calculating zeta minimum for each time
!     --------------------------------------
      minzetat (1) = 9999.d0
      do n=1,node_len
 	 if(zeta (n,1).gt.-9999.d0.and.zeta(n,1).lt. minzetat (1) )then 
            minzetat (1) = zeta (n,1)
            minznodet(1) = n
         endif
      end do

      iret=nf_put_vara_DOUBLE(ncid, minzetat_id,(/nrec/),(/1/),minzetat)
      CALL check_err(iret)
      iret=nf_put_vara_INT(ncid, minznodet_id,(/nrec/),(/1/),minznodet)
      CALL check_err(iret)
!     -------------------------------------
!     Calculating zeta maximum in each node
!     -------------------------------------
      do n=1,node_len
         if(zeta (n,1).gt.maxzetan61(n) ) then
           maxzetan61  (n) = zeta (n,1)
           maxztimen61 (n) = timenc(1)
         endif
      end do
!     -------------------------------------
!     Calculating zeta maximum for each time
!     -------------------------------------
      maxzetat (1) = -9999.d0
      do n=1,node_len
         if(zeta(n,1).gt.maxzetat(1))then
 	    maxzetat (1) = zeta (n,1)
 	    maxznodet(1) = n
         endif
      end do

      iret=nf_put_vara_DOUBLE(ncid,maxzetat_id,(/nrec/),(/1/),maxzetat)
      CALL check_err(iret)
      iret=nf_put_vara_INT(ncid,maxznodet_id,(/nrec/),(/1/),maxznodet)
      CALL check_err(iret)

!     ---------------------------
!     Closing fort.X3 netCDF file
!     ---------------------------
      IF(NREC.EQ.NDATA) THEN
         iret = nf_put_var_double(ncid, minzetan_id, minzetan61)
         CALL check_err(iret)
         iret = nf_put_var_double(ncid, maxzetan_id, maxzetan61)
         CALL check_err(iret)
         iret = nf_put_var_double(ncid, minztimen_id, minztimen61)
         CALL check_err(iret)
         iret = nf_put_var_double(ncid, maxztimen_id, maxztimen61)
         CALL check_err(iret)
      ENDIF
      
      WRITE(6,*) "Closing ", FILENAME, " netCDF file ... "
      WRITE(6,*) "-----------------------------------------------------"
      
9999  iret = nf_close(ncid)
      CALL check_err(iret)
      WRITE(6,*)'end subroutine writenc61'
     
      end subroutine writenc61      

!=======================================================================
!     Subroutine writenc62
!     Writes a netCDF file equivalent to ADCIRC's fort.62 ascii file
!
!     Dimensions:
!     station  : number of stations
!     namelen  : length of station's name
!     time     : number of time records
!
!     Variables :
!     time      : time
!     station_name : station's name
!     u-vel     : u velocity in x or east/west direction 
!     v-vel     : v velocity in y or north/south direction
!     x         : x-coordinate or longitude
!     y         : y-coordinate or latitude
!     minuvelt  : minimum vertically averaged u at each output time
!     maxuvelt  : maximum vertically averaged u at each output time
!     minvvelt  : minimum vertically averaged v at each output time
!     maxvvelt  : maximum vertically averaged v at each output time
!     minuveln  : minimum vertically averaged u at each station
!     maxuveln  : maximum vertically averaged u at each station
!     minvveln  : minimum vertically averaged v at each station
!     maxvveln  : maximum vertically averaged v at each station
!     minunodet : node where minimum u occurs at each output time
!     maxunodet : node where maximum u occurs at each output time
!     minvnodet : node where minimum v occurs at each output time
!     maxvnodet : node where maximum v occurs at each output time
!     minutimen : time when minimum u occurs at each station
!     maxutimen : time when maximum u occurs at each station
!     minvtimen : time when minimum v occurs at each station
!     maxvtimen : time when maximum v occurs at each station
!
!     Attributes:
!     model             : ADCIRC
!     version           : model version number
!     grid_type         : triangular
!     rundes            : description of simulation
!     agrid             : grid ID 
!     title             : description of the dataset
!     institution       : where the original data was produced
!     source            : method of production of the original data
!     history           : provides an audit trail for modifications to the original data
!     references        : published or web-based references that describe the data or methods used to produce it
!     comments          : miscellaneous information about the data
!     contact           : user's institution/ e-mail address
!     creation_date     : date the data was created
!     modification_date : date the data was modified
!     Parameters from fort.15 input file 
!
!     Revision history:
!
!     Date            Programmer                  Description of change
!     ----            ----------                  ----------------------
!     03/30/07  Cristina Forbes, PSGS @ UNC-IMS   Wrote original code
!     03/30/08  Cristina Forbes, PSGS @ UNC-IMS   Modified code for 
!                                                 globalio & hotstart 
!                                                 from binary files
!=======================================================================
      subroutine writenc62 (nout,fileunit,timenc_len,node_len,timenc,
     &                      nt, NDATA, u,v,
     &                      rundes, runid, agrid, irtype,
     &                      dtdp, ihot, ics, nolibf, nolifa,
     &                      nolica, nolicat, nwp, 
     &                      ncor, ntip, nws, nramp, tau0,
     &                      statim, reftim, rnday, dramp,
     &                      a00, b00, c00, h0,
     &                      slam0, sfea0, ffactor,esl, cori,
     &                      ntif, nbfr, elem_len, nface_len, nm,
     &                      x,y,depth,
     &                      nope_len, neta_len, nbdv, 
     &                      nbou_len, nvel_len, nbvv,
     &                      ncmode,
     &                      title, institution, source, history,
     &                      references, comments, host, convention,
     &                      contact,
     &                      base_date)


!--------------------------------------------
!     Force explicit declaration of all variables
!--------------------------------------------
      IMPLICIT NONE
!     -------------------
!     netCDF header file
!     -------------------
 !     include 'netcdf.inc'
!     -------------------
!     Error status return
!     -------------------
      INTEGER  iret
!     ---------
!     netCDF id
!     ---------
      INTEGER  ncid
!     -------------
!     Dimension ids
!     -------------
      INTEGER  node_dim
      INTEGER  timenc_dim
      INTEGER  elem_dim
      INTEGER  nface_dim
      INTEGER  nope_dim
      INTEGER  neta_dim
      INTEGER  nbou_dim
      INTEGER  nvel_dim
      INTEGER  slen_dim
      INTEGER  scode_dim
!     -----------------
!     Dimension lengths
!     -----------------
      INTEGER  node_len
      INTEGER  timenc_len
      INTEGER  elem_len
      INTEGER  nface_len
      INTEGER  nope_len
      INTEGER  neta_len
      INTEGER  nbou_len
      INTEGER  nvel_len
!     ------------
!     Variable ids
!     ------------
      INTEGER  timenc_id
      INTEGER  X_id
      INTEGER  Y_id
      INTEGER  DEPTH_id
      INTEGER  ELE_id
      INTEGER  mindepth_id
      INTEGER  maxdepth_id
      INTEGER  nbdv_id
      INTEGER  nbvv_id
      INTEGER  u_id
      INTEGER  v_id
      INTEGER  minut_id
      INTEGER  maxut_id
      INTEGER  minvt_id
      INTEGER  maxvt_id
      INTEGER  minun_id
      INTEGER  maxun_id
      INTEGER  minvn_id
      INTEGER  maxvn_id
      
      INTEGER  minunodet_id
      INTEGER  maxunodet_id
      INTEGER  minutimen_id
      INTEGER  maxutimen_id
      INTEGER  minvnodet_id
      INTEGER  maxvnodet_id
      INTEGER  minvtimen_id
      INTEGER  maxvtimen_id
      INTEGER  station_id
      INTEGER  scode_id
!     ---------------------------------------------
!     Rank (number of dimensions) for each variable
!     ---------------------------------------------
      INTEGER  timenc_rank
      INTEGER  station_rank
      INTEGER  scode_rank
      INTEGER  depth_rank
      INTEGER  ele_rank
      INTEGER  mindepth_rank
      INTEGER  maxdepth_rank
      INTEGER  x_rank
      INTEGER  y_rank
      INTEGER  nvdll_rank
      INTEGER  nbdv_rank
      INTEGER  nbvv_rank
      INTEGER  u_rank
      INTEGER  v_rank
      INTEGER  minut_rank
      INTEGER  maxut_rank
      INTEGER  minvt_rank
      INTEGER  maxvt_rank
      INTEGER  minun_rank
      INTEGER  maxun_rank
      INTEGER  minvn_rank
      INTEGER  maxvn_rank

      INTEGER  minunodet_rank
      INTEGER  maxunodet_rank
      INTEGER  minutimen_rank
      INTEGER  maxutimen_rank
      
      INTEGER  minvnodet_rank
      INTEGER  maxvnodet_rank
      INTEGER  minvtimen_rank
      INTEGER  maxvtimen_rank

      PARAMETER (timenc_rank = 1)
      PARAMETER (scode_rank = 2)
      PARAMETER (station_rank = 2)
      PARAMETER (depth_rank = 1)
      PARAMETER (ele_rank = 2)
      PARAMETER (mindepth_rank = 0)
      PARAMETER (maxdepth_rank = 0)
      PARAMETER (x_rank = 1)
      PARAMETER (y_rank = 1)
      PARAMETER (nvdll_rank = 1)
      PARAMETER (nbdv_rank = 2)
      PARAMETER (nbvv_rank = 2)
      PARAMETER (u_rank = 2)
      PARAMETER (v_rank = 2)
      PARAMETER (minut_rank = 1)
      PARAMETER (maxut_rank = 1)
      PARAMETER (minvt_rank = 1)
      PARAMETER (maxvt_rank = 1)
      PARAMETER (minun_rank = 1)
      PARAMETER (maxun_rank = 1)
      PARAMETER (minvn_rank = 1)
      PARAMETER (maxvn_rank = 1)
      
      PARAMETER (minunodet_rank = 1)
      PARAMETER (maxunodet_rank = 1)
      PARAMETER (minutimen_rank = 1)
      PARAMETER (maxutimen_rank = 1)

      PARAMETER (minvnodet_rank = 1)
      PARAMETER (maxvnodet_rank = 1)
      PARAMETER (minvtimen_rank = 1)
      PARAMETER (maxvtimen_rank = 1)
!     ---------------
!     Variable shapes
!     ---------------
      INTEGER  timenc_dims(timenc_rank)
      INTEGER  station_dims  (station_rank)
      INTEGER  x_dims(x_rank)
      INTEGER  y_dims(y_rank)
      INTEGER  depth_dims(depth_rank)
      INTEGER  ele_dims(ele_rank)
      INTEGER  nvdll_dims(nvdll_rank)
      INTEGER  nbdv_dims(nbdv_rank)
      INTEGER  nbvv_dims(nbvv_rank)
      INTEGER  u_dims(u_rank)
      INTEGER  v_dims(v_rank)
      INTEGER  minut_dims(minut_rank)
      INTEGER  maxut_dims(maxut_rank)
      INTEGER  minvt_dims(minvt_rank)
      INTEGER  maxvt_dims(maxvt_rank)
      INTEGER  minun_dims(minun_rank)
      INTEGER  maxun_dims(maxun_rank)
      INTEGER  minvn_dims(minvn_rank)
      INTEGER  maxvn_dims(maxvn_rank)

      INTEGER  minunodet_dims(minunodet_rank)
      INTEGER  maxunodet_dims(maxunodet_rank)
      INTEGER  minutimen_dims(minutimen_rank)
      INTEGER  maxutimen_dims(maxutimen_rank)

      INTEGER  minvnodet_dims(minvnodet_rank)
      INTEGER  maxvnodet_dims(maxvnodet_rank)
      INTEGER  minvtimen_dims(minvtimen_rank)
      INTEGER  maxvtimen_dims(maxvtimen_rank)

!     --------------
!     Data variables
!     --------------
      CHARACTER*5 FILEFORT
      CHARACTER*10 FILENAME
      CHARACTER*3 FILEEXT
      INTEGER   t,ntif,nbfr,nt
      INTEGER   irtype, ihot, ics, nolibf, nolifa, nolica, nolicat, nwp
      INTEGER   ntip, nws, nramp,ncor
      INTEGER   nvel, neta, nout
      INTEGER   element (nface_len , elem_len)
      INTEGER   nbdv (nope_len, neta_len)
      INTEGER   nbvv (nbou_len, nvel_len)
      INTEGER   nm (elem_len, nface_len)
      INTEGER   fileunit
      REAL*8    timenc(timenc_len)
      REAL(sz)  x(node_len)
      REAL(sz)  y(node_len)
      REAL(sz)  depth(node_len)
      REAL(sz)  u(node_len, timenc_len)
      REAL(sz)  v(node_len, timenc_len)
      REAL(sz), ALLOCATABLE, SAVE ::   minut(:)
      REAL(sz), ALLOCATABLE, SAVE ::   maxut(:)
      REAL(sz), ALLOCATABLE, SAVE ::   minvt(:)
      REAL(sz), ALLOCATABLE, SAVE ::   maxvt(:)
      REAL(sz), ALLOCATABLE, SAVE ::   minun62(:)
      REAL(sz), ALLOCATABLE, SAVE ::   maxun62(:)
      REAL(sz), ALLOCATABLE, SAVE ::   minvn62(:)
      REAL(sz), ALLOCATABLE, SAVE ::   maxvn62(:)

      INTEGER , ALLOCATABLE, SAVE ::  minunodet(:)
      INTEGER , ALLOCATABLE, SAVE ::  maxunodet(:)
      REAL(sz), ALLOCATABLE, SAVE ::  minutimen62(:)
      REAL(sz), ALLOCATABLE, SAVE ::  maxutimen62(:)
      
      INTEGER , ALLOCATABLE, SAVE ::  minvnodet(:)
      INTEGER , ALLOCATABLE, SAVE ::  maxvnodet(:)
      REAL(sz), ALLOCATABLE, SAVE ::  minvtimen62(:)
      REAL(sz), ALLOCATABLE, SAVE ::  maxvtimen62(:)

      REAL  minull
      REAL  maxull
      REAL  minvll
      REAL  maxvll
      real*8  mindepth
      real*8  maxdepth
      REAL*8 dtdp, statim,reftim, rnday, dramp, tau0
      REAL*8 a00, b00, c00,h0
      REAL(sz) slam0, sfea0, ffactor,esl, cori
      REAL(sz) slam0deg, sfea0deg
      CHARACTER(80) rundes, runid, agrid
      INTEGER, SAVE  :: FIRST62 = 1 
      INTEGER NDATA
      INTEGER, SAVE  :: NREC = 0
      INTEGER ncmode

!     -----------------
!     Attribute vectors
!     -----------------
      INTEGER  intval(1),i,j
      REAL  realval(1)
      INTEGER counti(1), starti(1), n
      INTEGER count(2), start(2)
!      DOUBLE PRECISION  doubleval(1)
!     -----------------
!     Attribute general
!     -----------------
      CHARACTER(*) :: title, institution, source, history
      CHARACTER(*) :: references, comments, host, convention, contact
      CHARACTER(*) :: base_date

!     -----------------
!      date_string variables for time attribute
!     -----------------
      character date_string*40
      character now_date*8
      character big_ben*10
      character zone*5
      integer values(8)
!     ---------------------------------
!      inquiring and reading variables
!     ---------------------------------

      INTEGER NDIMS, NVARS, NGATTS, UNLIMDIMID, RECID,NRECS
      CHARACTER*4 recname

!     -------------------------
!     Save dims, ids, variables
!     -------------------------
      save node_dim,  elem_dim, nface_dim
      save nope_dim, neta_dim, nbou_dim, nvel_dim
      save u_id, v_id, X_id, Y_id
      save minut_id, maxut_id, minun_id, maxun_id
      save minvt_id, maxvt_id, minvn_id, maxvn_id

      save minunodet_id, maxunodet_id, minutimen_id, maxutimen_id
      save minvnodet_id, maxvnodet_id, minvtimen_id, maxvtimen_id

      save DEPTH_id,ELE_id, mindepth_id, maxdepth_id, nbdv_id, nbvv_id
      save timenc_dim, timenc_id
      save station_id, scode_id
      save mindepth, maxdepth

!     ---------------------------
!     Opening fort.62 netCDF file
!     ---------------------------
      WRITE(6,*) "-----------------------------------------------------"

       
      IF(FIRST62.EQ.1.and.ncmode.NE.0) THEN
         ALLOCATE ( minun62 (node_len) )
         ALLOCATE ( maxun62 (node_len) )
         ALLOCATE ( minvn62 (node_len) )
         ALLOCATE ( maxvn62 (node_len) )
         
         ALLOCATE ( minutimen62 (node_len) )
         ALLOCATE ( maxutimen62 (node_len) )
         ALLOCATE ( minvtimen62 (node_len) )
         ALLOCATE ( maxvtimen62 (node_len) )
         
         DO n=1,node_len
            minun62 (n) =  9999.d0
            maxun62 (n) = -9999.d0
            minvn62 (n) =  9999.d0
            maxvn62 (n) = -9999.d0

            minutimen62 (n) =  9999.d0
            maxutimen62 (n) = -9999.d0
            minvtimen62 (n) =  9999.d0
            maxvtimen62 (n) = -9999.d0
         END DO
       
      ENDIF
        
!     --------------------------------
!     Switch order in array for NETCDF
!     --------------------------------
      do i=1, elem_len
         do j=1, nface_len
            element(j,i) = nm (i,j)
         end do
      end do
      
!     -----------------
!     Enter define mode
!     -----------------
         FILEFORT = 'fort.'
         FILEEXT = '.nc'
         WRITE(FILENAME,'(A5,I2,A3)')FILEFORT,FILEUNIT,FILEEXT
         WRITE(6,*) "Opening ",FILENAME,  " netCDF file ... "
         
       IF((FIRST62.EQ.1.AND.IHOT.EQ.0.AND.NCMODE.NE.2).OR.
     &    (FIRST62.EQ.1.AND.IHOT.NE.0.AND.NCMODE.NE.2.AND.
     &     nout.EQ.-3))THEN

          iret = nf_create(FILENAME, NF_CLOBBER, ncid)
          CALL check_err(iret)
          FIRST62=0
!     --------------------------------
!     Convert to degress for NETCDF I/O
!     --------------------------------
          SLAM0DEG=SLAM0
          SFEA0DEG=SFEA0
          IF (NCMODE.EQ.1) THEN
              SLAM0DEG=SLAM0*RAD2DEG
              SFEA0DEG=SFEA0*RAD2DEG
          ENDIF
!     -----------------
!     Define dimensions
!     -----------------
      iret = nf_def_dim(ncid, 'station', node_len, node_dim)
      CALL check_err(iret)
      iret = nf_def_dim(ncid, 'namelen', SNAMLEN, slen_dim)
             CALL check_err(iret)
      iret = nf_def_dim(ncid, 'time', nf_unlimited, timenc_dim)
      CALL check_err(iret)

!     ================
!     Define variables
!     ================
!     ----------------
!     Define time
!     ----------------
      timenc_dims(1) = timenc_dim
      iret = nf_def_var(ncid,'time',NF_DOUBLE, timenc_rank, timenc_dims, 
     &                  timenc_id)
      CALL check_err(iret)
      
!     -------------------
!     Define x
!     -------------------
      x_dims(1) = node_dim
      iret = nf_def_var(ncid, 'x', NF_DOUBLE, x_rank, x_dims, x_id)
      CALL check_err(iret)
!     -------------------
!     Define Y coordinate
!     --------------------
      Y_dims(1) = node_dim
      iret = nf_def_var(ncid, 'y', NF_DOUBLE, Y_rank, Y_dims, Y_id)
      CALL check_err(iret)
!     ---------------------
!     Define stations name
!     --------------------
      station_dims(1) = slen_dim
      station_dims(2) = node_dim
      iret = nf_def_var(ncid, 'station_name', NF_CHAR, station_rank, 
     1      station_dims, station_id)
      CALL check_err(iret)

!     ----------------
!     Define u vel
!     ----------------
      u_dims(1) = node_dim
      u_dims(2) = timenc_dim
      iret = nf_def_var(ncid, 'u-vel', NF_DOUBLE, u_rank, u_dims, 
     &                  u_id)
      CALL check_err(iret)

!     ----------------
!     Define v vel
!     ----------------
      v_dims(1) = node_dim
      v_dims(2) = timenc_dim
      iret = nf_def_var(ncid, 'v-vel', NF_DOUBLE, v_rank, v_dims, 
     &                  v_id)
      CALL check_err(iret)

!     -------------------------
!     Define minimum u at times
!     -------------------------
      minut_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'minuvelt', NF_DOUBLE, minut_rank, 
     &                  minut_dims, minut_id)
      CALL check_err(iret)
!     -------------------------
!     Define maximum u at times
!     -------------------------
      maxut_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'maxuvelt', NF_DOUBLE, maxut_rank, 
     &                  maxut_dims, maxut_id)
      CALL check_err(iret)

!     -------------------------
!     Define minimum v at times
!     -------------------------
      minvt_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'minvvelt', NF_DOUBLE, minvt_rank, 
     &                  minvt_dims, minvt_id)
      CALL check_err(iret)
      
!     -------------------------
!     Define maximum v at times
!     -------------------------
      maxvt_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'maxvvelt', NF_DOUBLE, maxvt_rank, 
     &                  maxvt_dims, maxvt_id)
      CALL check_err(iret)

!     -------------------------
!     Define minimum u at times
!     -------------------------
      minun_dims(1) = node_dim
      iret = nf_def_var(ncid, 'minuveln', NF_DOUBLE, minun_rank, 
     &                  minun_dims, minun_id)
      CALL check_err(iret)
!     -------------------------
!     Define maximum u at nodes
!     -------------------------
      maxun_dims(1) = node_dim
      iret = nf_def_var(ncid, 'maxuveln', NF_DOUBLE, maxun_rank, 
     &                  maxun_dims, maxun_id)
      CALL check_err(iret)
!     -------------------------
!     Define minimum v at nodes
!     -------------------------
      minvn_dims(1) = node_dim
      iret = nf_def_var(ncid, 'minvveln', NF_DOUBLE, minvn_rank, 
     &                  minvn_dims, minvn_id)
      CALL check_err(iret)
!     -------------------------
!     Define maximum v at nodes
!     -------------------------
      maxvn_dims(1) = node_dim
      iret = nf_def_var(ncid, 'maxvveln', NF_DOUBLE, maxvn_rank, 
     &                  maxvn_dims, maxvn_id)
      CALL check_err(iret)
      
!     ----------------------------
!     Define minimum u nodes at times
!     ----------------------------
      minunodet_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'minunodet', NF_INT, minunodet_rank, 
     &                  minunodet_dims, minunodet_id)
      CALL check_err(iret)
!     ----------------------------
!     Define maximum u nodes at times
!     ----------------------------
      maxunodet_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'maxunodet', NF_INT, maxunodet_rank, 
     &                  maxunodet_dims, maxunodet_id)
      CALL check_err(iret)
!     ----------------------------
!     Define minimum u nodes at nodes
!     ----------------------------
      minutimen_dims(1) = node_dim
      iret = nf_def_var(ncid, 'minutimen', NF_DOUBLE, minutimen_rank, 
     &                  minutimen_dims, minutimen_id)
      CALL check_err(iret)
!     ----------------------------
!     Define maximum u nodes at nodes
!     ----------------------------
      maxutimen_dims(1) = node_dim
      iret = nf_def_var(ncid, 'maxutimen', NF_DOUBLE, maxutimen_rank, 
     &                  maxutimen_dims, maxutimen_id)
      CALL check_err(iret)
            
!     ----------------------------
!     Define minimum v nodes at times
!     ----------------------------
      minvnodet_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'minvnodet', NF_INT, minvnodet_rank, 
     &                  minvnodet_dims, minvnodet_id)
      CALL check_err(iret)
!     ----------------------------
!     Define maximum v nodes at times
!     ----------------------------
      maxvnodet_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'maxvnodet', NF_INT, maxvnodet_rank, 
     &                  maxvnodet_dims, maxvnodet_id)
      CALL check_err(iret)
!     ----------------------------
!     Define minimum v nodes at nodes
!     ----------------------------
      minvtimen_dims(1) = node_dim
      iret = nf_def_var(ncid, 'minvtimen', NF_DOUBLE, minvtimen_rank, 
     &                  minvtimen_dims, minvtimen_id)
      CALL check_err(iret)
!     ----------------------------
!     Define maximum v nodes at nodes
!     ----------------------------
      maxvtimen_dims(1) = node_dim
      iret = nf_def_var(ncid, 'maxvtimen', NF_DOUBLE, maxvtimen_rank, 
     &                  maxvtimen_dims, maxvtimen_id)
      CALL check_err(iret)
      
!     -----------------
!     Assign attributes
!     -----------------

      if(ics.eq.2) then
      
      iret = nf_put_att_text(ncid, X_id, 'long_name', 9, 'longitude')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, X_id, 'standard_name', 9, 
     &                      'longitude')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, X_id, 'units', 12,  VAR_DEGE)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, X_id, 'positive', 4, VAR_EAST)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, X_id, 'missing_values', NF_DOUBLE,
     &                        1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, X_id, '_FillValue', NF_DOUBLE, 1, 
     &                       doubleval)
      CALL check_err(iret)      
      
      iret = nf_put_att_text(ncid, Y_id, 'long_name', 8, 'latitude')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, Y_id, 'standard_name', 8, 
     &                       'latitude')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, Y_id, 'units', 13,  VAR_DEGN)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, Y_id, 'positive', 5, VAR_NORTH)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, Y_id, 'missing_values', NF_DOUBLE, 
     &       1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, Y_id, '_FillValue', NF_DOUBLE, 1, 
     &       doubleval)
      CALL check_err(iret)

       else
      
      iret = nf_put_att_text(ncid, X_id, 'long_name', 22, 
     &       'Cartesian coordinate x')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, X_id, 'standard_name', 12, 
     &                      'x_coordinate')
      CALL check_err(iret)
      if (G.LT.11.D0) then
         iret = nf_put_att_text(ncid, X_id, 'units', 6, VAR_METS)
          CALL check_err(iret)
      else
         iret = nf_put_att_text(ncid, X_id, 'units', 4, VAR_FEET)
          CALL check_err(iret)
      endif
      iret = nf_put_att_text(ncid, X_id, 'positive', 5, VAR_RG)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, X_id, 'missing_values', NF_DOUBLE, 
     &                       1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, X_id, '_FillValue', NF_DOUBLE, 1, 
     &                       doubleval)
      CALL check_err(iret)
!     ------------------------------
!     Define y-coordinate attributes
!     ------------------------------     
      
      iret = nf_put_att_text(ncid, Y_id, 'long_name', 22, 
     &       'Cartesian coordinate y')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, Y_id, 'standard_name', 12, 
     &                       'y_coordinate')
      CALL check_err(iret)
      if (G.LT.11.D0) then
         iret = nf_put_att_text(ncid, Y_id, 'units', 6, VAR_METS)
         CALL check_err(iret)
      else
         iret = nf_put_att_text(ncid, Y_id, 'units', 4, VAR_FEET)
         CALL check_err(iret)
      endif
      iret = nf_put_att_text(ncid, Y_id, 'positive', 2, VAR_UP)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, Y_id, 'missing_values', 
     &       NF_DOUBLE, 1, doubleval)
      CALL check_err(iret)
!      doubleval(1) = -99999
      iret = nf_put_att_double(ncid, Y_id, '_FillValue', NF_DOUBLE, 1, 
     &       doubleval)
      CALL check_err(iret)
      
      endif
      
!    ----------------------------------------------------
!     minimum/maximum e/w and n/s velocities at each time
!    ----------------------------------------------------

      if (ics.eq.2) then
      iret = nf_put_att_text(ncid, minut_id, 'long_name', 68, 
     & 'minimum station vertically averaged e/w velocity at each output 
     &time')
      CALL check_err(iret)

      else

      iret = nf_put_att_text(ncid, minut_id, 'long_name', 71, 
     &     'minimum vertically averaged velocity in x-direction at each 
     &output time')
      CALL check_err(iret)
      endif
      
      iret = nf_put_att_text(ncid, minut_id, 'standard_name', 26,
     &                'minimum_station_u_velocity')
      CALL check_err(iret)

      if (G.LT.11.D0) then
          iret = nf_put_att_text(ncid, minut_id,'units',13, VAR_MPS)
          CALL check_err(iret)
      else
          iret = nf_put_att_text(ncid, minut_id,'units',11, VAR_FPS)
          CALL check_err(iret)
      endif
      
      if (ics.eq.2) then
          iret = nf_put_att_text(ncid, maxut_id, 'long_name', 68, 
     & 'maximum station vertically averaged e/w velocity at each output 
     &time')
          CALL check_err(iret)

       else

      iret = nf_put_att_text(ncid, maxut_id, 'long_name', 78, 
     &     'maximum vertically averaged velocity in x-direction at each 
     &output time')
      CALL check_err(iret)
      endif
       
      iret = nf_put_att_text(ncid, maxut_id, 'standard_name', 26,
     &                'maximum_station_u_velocity')
      CALL check_err(iret)
      if (G.LT.11.D0) then
         iret = nf_put_att_text(ncid, maxut_id,'units',13, VAR_MPS)
      else
         iret = nf_put_att_text(ncid, maxut_id,'units',11, VAR_FPS)
      endif
      CALL check_err(iret)
      
      if (ics.eq.2) then
         iret = nf_put_att_text(ncid, minvt_id, 'long_name', 69, 
     &  'minimum station vertically averaged n/s velocity at each output 
     & time')
      else

      iret = nf_put_att_text(ncid, minvt_id, 'long_name', 79, 
     &  'minimum station vertically averaged velocity in y-direction at 
     &each output time')
      endif
     
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, minvt_id, 'standard_name', 26,
     &                'minimum_station_v_velocity')
      CALL check_err(iret)
      if (G.LT.11.D0) then
          iret = nf_put_att_text(ncid, minvt_id,'units',13, VAR_MPS)
      else
          iret = nf_put_att_text(ncid, minvt_id,'units',11, VAR_FPS)
      endif
      CALL check_err(iret)
      
      if (ics.eq.2) then
         iret = nf_put_att_text(ncid, maxvt_id, 'long_name', 68, 
     &        'maximum station vertically averaged n/s velocity at each 
     &output time')
      else
         iret = nf_put_att_text(ncid, maxvt_id, 'long_name', 79, 
     &  'maximum station vertically averaged velocity in y-direction at 
     &each output time')
      endif
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, maxvt_id, 'standard_name', 26,
     &                'maximum_station_v_velocity')
      CALL check_err(iret)
      if (G.LT.11.D0) then
          iret = nf_put_att_text(ncid, maxvt_id,'units',13, VAR_MPS)
      else
          iret = nf_put_att_text(ncid, maxvt_id,'units',11, VAR_FPS)
      endif
      CALL check_err(iret)

!    ----------------------------------------------------
!     minimum/maximum e/w and n/s velocities at each node
!    ----------------------------------------------------

      if (ics.eq.2) then
      
          iret = nf_put_att_text(ncid, minun_id, 'long_name', 56, 
     &    'minimum vertically averaged e/w velocity at each station')
          CALL check_err(iret)

       else

          iret = nf_put_att_text(ncid, minun_id, 'long_name', 67, 
     &     'minimum vertically averaged velocity in x-direction at each 
     &station')
          CALL check_err(iret)
          
      endif
      
      iret = nf_put_att_text(ncid, minun_id, 'standard_name', 26,
     &                'minimum_station_u_velocity')
      CALL check_err(iret)
      if (G.LT.11.D0) then
          iret = nf_put_att_text(ncid, minun_id,'units',13, VAR_MPS)
      else
          iret = nf_put_att_text(ncid, minun_id,'units',11, VAR_FPS)
      endif
      CALL check_err(iret)
      
      if (ics.eq.2) then
      iret = nf_put_att_text(ncid, maxun_id, 'long_name', 56, 
     & 'maximum vertically averaged e/w velocity at each station')
      CALL check_err(iret)

       else

      iret = nf_put_att_text(ncid, maxun_id, 'long_name', 67, 
     &     'maximum vertically averaged velocity in x-direction at each 
     &station')
      CALL check_err(iret)
      endif
       
      iret = nf_put_att_text(ncid, maxun_id, 'standard_name', 26,
     &                'maximum_station_u_velocity')
      CALL check_err(iret)
      
      if (G.LT.11.D0) then
          iret = nf_put_att_text(ncid, maxun_id,'units',13, VAR_MPS)
      else
          iret = nf_put_att_text(ncid, maxun_id,'units',11, VAR_FPS)
      endif
      CALL check_err(iret)
      
      if (ics.eq.2) then
      iret = nf_put_att_text(ncid, minvn_id, 'long_name', 56, 
     &      'minimum vertically averaged n/s velocity at each station')

      else

      iret = nf_put_att_text(ncid, minvn_id, 'long_name', 66, 
     &     'minimum vertically averaged velocity in y-direction at each 
     &station')
      endif
     
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, minvn_id, 'standard_name', 26,
     &                'minimum_station_v_velocity')
      CALL check_err(iret)
      if (G.LT.11.D0) then
          iret = nf_put_att_text(ncid, minvn_id,'units',13, VAR_MPS)
      else
          iret = nf_put_att_text(ncid, minvn_id,'units',11, VAR_FPS)
      endif
      CALL check_err(iret)
      
      if (ics.eq.2) then
      iret = nf_put_att_text(ncid, maxvn_id, 'long_name', 56, 
     &      'maximum vertically averaged n/s velocity at each station')
      CALL check_err(iret)

       else

      iret = nf_put_att_text(ncid, maxvn_id, 'long_name', 66, 
     &     'maximum vertically averaged velocity in y-direction at each 
     &station')
      CALL check_err(iret)
      endif

      iret = nf_put_att_text(ncid, maxvn_id, 'standard_name', 26,
     &                'maximum_station_v_velocity')
      CALL check_err(iret)
!    --------------------------
!     determine variable units
!    --------------------------
      if (G.LT.11.D0) then
          iret = nf_put_att_text(ncid, maxvn_id,'units',13, VAR_MPS)
      else
          iret = nf_put_att_text(ncid, maxvn_id,'units',11, VAR_FPS)
      endif
      CALL check_err(iret)

!    ----------------------------------------------------
!     minimum/maximum e/w and n/s velocities at each time
!    ----------------------------------------------------
      
      if (ics.eq.2) then
          iret = nf_put_att_text(ncid, minunodet_id, 'long_name', 52, 
     &    'station of minimum e/w velocity at each output time')
           CALL check_err(iret)

       else

          iret = nf_put_att_text(ncid, minunodet_id, 'long_name', 62, 
     &    'station of minimum velocity in x-direction at each output 
     &time')
          CALL check_err(iret)
      endif
      
      iret = nf_put_att_text(ncid, minunodet_id, 'standard_name', 31,
     &                'station_of_minimum_e/w_velocity')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, minunodet_id,'units',14,VAR_NOND)
      CALL check_err(iret)

      if (ics.eq.2) then
          iret = nf_put_att_text(ncid, maxunodet_id, 'long_name', 52, 
     &    'station of maximum e/w velocity at each output time')
          CALL check_err(iret)
       else

           iret = nf_put_att_text(ncid, maxunodet_id, 'long_name', 62, 
     &     'station of maximum velocity in x-direction at each output 
     &time')
           CALL check_err(iret)
      endif
       
      iret = nf_put_att_text(ncid, maxunodet_id, 'standard_name', 31,
     &                'station_of_maximum_e/w_velocity')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, maxunodet_id,'units',14, VAR_NOND)
      CALL check_err(iret)
      
      if (ics.eq.2) then
          iret = nf_put_att_text(ncid, minvnodet_id, 'long_name', 52, 
     &    'station of minimum n/s velocity at each output time')

       else

          iret = nf_put_att_text(ncid, minvnodet_id, 'long_name', 62, 
     &    'station of minimum velocity in y-direction at each output 
     &time')
      endif
      CALL check_err(iret)
      
      iret = nf_put_att_text(ncid, minvnodet_id, 'standard_name', 31,
     &                'station_of_minimum_n/s_velocity')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, minvnodet_id,'units',14, VAR_NOND)
      CALL check_err(iret)
      
      if (ics.eq.2) then
          iret = nf_put_att_text(ncid, maxvnodet_id, 'long_name', 52, 
     &    'station of maximum n/s velocity at each output time')
 
       else

          iret = nf_put_att_text(ncid, maxvnodet_id, 'long_name', 62, 
     &    'station of maximum velocity in y-direction at each output 
     &time')
      endif
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, maxvnodet_id, 'standard_name', 31,
     &                'station_of_maximum_n/s_velocity')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, maxvnodet_id,'units',14,VAR_NOND)
      CALL check_err(iret)

!    ----------------------------------------------------
!     minimum/maximum e/w and n/s velocities at each node
!    ----------------------------------------------------
      
!    ----------------------
!     determine coordinates
!    ----------------------
      if (ics.eq.2) then
         iret = nf_put_att_text(ncid, minutimen_id, 'long_name', 44, 
     &   'time of minimum e/w velocity at each station')

      else

          iret = nf_put_att_text(ncid, minutimen_id, 'long_name', 55, 
     &   'time of minimum velocity in x-direction at each station')

      endif
      CALL check_err(iret)
    
      iret = nf_put_att_text(ncid, minutimen_id, 'standard_name', 28,
     &                'time_of_minimum_e/w_velocity')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, minutimen_id,'units',4, VAR_SECS)
      CALL check_err(iret)
      
!    ----------------------
!     determine coordinates
!    ----------------------
      if (ics.eq.2) then
          iret = nf_put_att_text(ncid, maxutimen_id, 'long_name', 44, 
     &   'time of maximum e/w velocity at each station')

       else

          iret = nf_put_att_text(ncid, maxutimen_id, 'long_name', 45, 
     &    'time of maximum velocity in x-direction at each station')
      endif
      CALL check_err(iret)
       
      iret = nf_put_att_text(ncid, maxutimen_id, 'standard_name', 28,
     &                'time_of_maximum_e/w_velocity')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, maxutimen_id,'units',4, VAR_SECS)
      CALL check_err(iret)

!    ----------------------
!     determine coordinates
!    ----------------------
      if (ics.eq.2) then
          iret = nf_put_att_text(ncid, minvtimen_id, 'long_name', 44, 
     &    'time of minimum n/s velocity at each station')
       else
          iret = nf_put_att_text(ncid, minvtimen_id, 'long_name', 55, 
     &    'time of minimum velocity in y-direction at each station')
      endif
     
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, minvtimen_id, 'standard_name', 28,
     &                'time_of_minimum_n/s_velocity')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, minvtimen_id,'units',4, VAR_SECS)
      CALL check_err(iret)
     
      if (ics.eq.2) then
      iret = nf_put_att_text(ncid, maxvtimen_id, 'long_name', 44, 
     &    'time of maximum n/s velocity at each station')
      CALL check_err(iret)

       else

      iret = nf_put_att_text(ncid, maxvtimen_id, 'long_name', 55, 
     &'time of maximum velocity in y-direction at each station')
      CALL check_err(iret)
      endif

      iret = nf_put_att_text(ncid, maxvtimen_id, 'standard_name', 28,
     &                'time_of_maximum_n/s_velocity')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, maxvtimen_id,'units',4, VAR_SECS)
      CALL check_err(iret)
      
!     -------------------------------------------------------------
!     Define station names and code attributes
!     -------------------------------------------------------------
      iret = nf_put_att_text(ncid,station_id,'long_name',12,
     & 'station name')
      CALL check_err(iret)

!    -----------------------
!     Define time attributes
!    -----------------------
      iret = nf_put_att_text(ncid,timenc_id,'long_name',10,'model time')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, timenc_id, 'standard_name', 4,'time')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, timenc_id, 'units', 23, VAR_TUN)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, timenc_id, 'base_date',
     . (LEN(TRIM(base_date(1:INDEX(base_date,"!")-1)))), base_date)
      CALL check_err(iret)

!    -----------
!     Define u velocity
!    -----------
!    ----------------------
!     determine coordinates
!    ----------------------

      if (ics.eq.2) then
          iret = nf_put_att_text(ncid, u_id, 'long_name', 38, 
     &    'station vertically averaged e/w velocity')
      else

      iret = nf_put_att_text(ncid, u_id, 'long_name', 49, 
     &       'station vertically averaged velocity in x-direction')
      endif
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, u_id, 'standard_name', 18, 
     &       'station_u_velocity')
      CALL check_err(iret)
!    ----------------------
!     determine variable units
!    ----------------------
      if (G.LT.11.D0) then
          iret = nf_put_att_text(ncid, u_id, 'units',13, VAR_MPS)
      else
          iret = nf_put_att_text(ncid, u_id, 'units',11, VAR_FPS)
      endif
      CALL check_err(iret)
      
!    ----------------------
!     determine coordinates
!    ----------------------
      if (ics.eq.2) then
          iret = nf_put_att_text(ncid, u_id, 'positive', 4, VAR_EAST)
      else
          iret = nf_put_att_text(ncid, u_id, 'positive', 5, VAR_RG)
      endif
      CALL check_err(iret)
      
      iret = nf_put_att_double(ncid, u_id, 'missing_values', 
     &       NF_DOUBLE, 1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, u_id, '_FillValue', NF_DOUBLE, 1, 
     &                       doubleval)
      CALL check_err(iret)
      
      iret = nf_put_att_double(ncid, u_id, 'dry_Value', NF_DOUBLE, 1, 
     &                       doubleval)
      CALL check_err(iret)
!    -----------
!     Define v velocity
!    -----------
!    ----------------------
!     determine coordinates
!    ----------------------

      if (ics.eq.2) then
         iret = nf_put_att_text(ncid, v_id, 'long_name', 38, 
     &         'station vertically averaged n/s velocity')

       else

          iret = nf_put_att_text(ncid, v_id, 'long_name', 49, 
     &          'station vertically averaged velocity in y-direction')
      endif
      CALL check_err(iret)


      iret = nf_put_att_text(ncid, v_id, 'standard_name', 18, 
     &                      'station_v_velocity')
      CALL check_err(iret)
!    ----------------------
!     determine variable units
!    ----------------------
      if (G.LT.11.D0) then
          iret = nf_put_att_text(ncid, v_id, 'units',13, VAR_MPS)
      else
          iret = nf_put_att_text(ncid, v_id, 'units',11, VAR_FPS)
      endif
      CALL check_err(iret)

!    ----------------------
!     determine coordinates
!    ----------------------
      if (ics.eq.2) then
         iret = nf_put_att_text(ncid, v_id, 'positive', 5, VAR_NORTH)
       else
         iret = nf_put_att_text(ncid, v_id, 'positive', 2, VAR_UP)
      endif
      CALL check_err(iret)

      iret = nf_put_att_DOUBLE(ncid, v_id, 'missing_values', 
     &       NF_DOUBLE, 1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_DOUBLE(ncid, v_id, '_FillValue', NF_DOUBLE, 1, 
     &                       doubleval)
      CALL check_err(iret)
      iret = nf_put_att_DOUBLE(ncid, v_id, 'dry_Value', NF_DOUBLE, 1, 
     &                       doubleval)
      CALL check_err(iret)

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
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'description', 
     &       LEN(TRIM(rundes(1:INDEX(rundes,"!")-1))), rundes)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'grid',
     &       LEN(TRIM(agrid(1:INDEX(agrid,"!")-1))), 
     &       agrid)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'title',
     &       LEN(TRIM(title(1:INDEX(title,"!")-1))), 
     &       title)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'institution', 
     &       LEN(TRIM(institution(1:INDEX(institution,"!")-1))), 
     &       institution)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'source', 
     &       LEN(TRIM(source(1:INDEX(source,"!")-1))),
     &       source)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'history', 
     &       LEN(TRIM(history(1:INDEX(history,"!")-1))),
     &       history)
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
     &       LEN(TRIM(host(1:INDEX(host,"!")-1))), 
     &       host)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'convention', 
     &       LEN(TRIM(convention(1:INDEX(convention,"!")-1))), 
     &       convention)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'contact', 
     &       LEN(TRIM(contact(1:INDEX(contact,"!")-1))), 
     &       contact)
      CALL check_err(iret)
 
 
      CALL date_and_Time(now_date,big_ben,zone,values)
      write(date_string,71) values(1),values(2),values(3)
     & ,values(5),values(6),values(7),   (values(4))/60
 71   format(I4,'-',I2.2,'-',i2.2,' ',i2,':',i2.2,':',i2.2,' '
     &  ,i3.2,':00')
      iret = nf_put_att_text(ncid, NF_GLOBAL,'creation_date'
     & ,len_trim(date_string),date_string)
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
      iret = nf_put_att_int(ncid, NCGLOBAL, 'irtype', NF_INT, 1, irtype)
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
      intval(1) = 5
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
      iret = nf_put_att_double(ncid, NCGLOBAL, 'ffactor', NF_DOUBLE, 1, 
     &                         ffactor)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, NCGLOBAL, 'esl', NF_DOUBLE, 1, 
     &                         esl)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, NCGLOBAL, 'cori', NF_DOUBLE, 1, 
     &                         cori)
      CALL check_err(iret)
      iret = nf_put_att_int(ncid, NCGLOBAL, 'ntif', NF_INT, 1, ntif)
      CALL check_err(iret)
      iret = nf_put_att_int(ncid, NCGLOBAL, 'nbfr', NF_INT, 1, nbfr)
      CALL check_err(iret)
     
!     ======================
!     Leave define mode
!     ======================
      iret = nf_enddef(ncid)
      CALL check_err(iret)
      
!     ======================
!     Store static variables
!     ======================
!     -------
!     Store X
!     -------
      iret = nf_put_var_double(ncid, X_id, X)
      CALL check_err(iret)
!     -------
!     Store Y
!     -------
      iret = nf_put_var_double(ncid, Y_id, Y)
      CALL check_err(iret)
!     --------------------
!     Store station name
!     --------------------
      do i=1,node_len
         start(1)=1
         start(2)=i
         count(1)= LEN(TRIM(STATNAMEV(i,SNAMLEN)))
         count(2)=1
         iret = nf_put_vara_text(ncid,station_id,start,count,
     &   STATNAMEV(i,SNAMLEN) )
         CALL check_err(iret)
      end do

!     =====================
!     Open old netcdf file
!     =====================
      ELSE
          iret = NF_OPEN(FILENAME, NF_WRITE, ncid)
          CALL check_err(iret)   
!     ------------------
!     Inquire variables
!     ------------------
          iret=nf_inq_unlimdim(ncid,recid)
          CALL check_err(iret)   
          iret=nf_inq_dimlen(ncid,recid,nrec)
          CALL check_err(iret)   
      ENDIF
!     ------------------
!     Update Global Attributes
!     ------------------
      IF(IHOT.NE.0.AND.FIRST62.EQ.1) THEN
         iret = NF_REDEF (ncid)
         CALL check_err(iret)
         iret = nf_put_att_int(ncid, NCGLOBAL, 'ihot', NF_INT, 1, ihot)
         CALL check_err(iret)
         iret = nf_put_att_int(ncid, NCGLOBAL, 'nramp', NF_INT, 1,nramp)
         CALL check_err(iret)
         iret = nf_put_att_double(ncid, NCGLOBAL, 'tau0', NF_DOUBLE, 1, 
     &                         tau0)
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
         iret = nf_put_att_double(ncid, NCGLOBAL, 'ffactor', NF_DOUBLE, 1, 
     &                         ffactor)
         CALL check_err(iret)
         iret = nf_put_att_double(ncid, NCGLOBAL, 'esl', NF_DOUBLE, 1, 
     &                         esl)
         CALL check_err(iret)
         CALL date_and_time(now_date,big_ben,zone,values)
         WRITE(date_string,71) values(1),values(2),values(3),
     &         values(5),values(6),values(7),   (values(4))/60
         iret = nf_put_att_text(ncid, NF_GLOBAL,'modification_date',
     &          len_trim(date_string),date_string)
         CALL check_err(iret)
         iret = NF_ENDDEF (ncid)
         CALL check_err(iret)
!        WRITE(6,*)"finish global att"
       ENDIF
      IF(ncmode.eq.0) go to 9999
      IF(FIRST62.EQ.1.and.
     & ((ncmode.ne.0.and.IHOT.NE.0).or.(ncmode.eq.2.and.IHOT.eq.0).or.
     &  (ncmode.eq.2.and.IHOT.ne.0.and.nout.eq.-3)))THEN
         iret=nf_inq_varid(ncid,"time",timenc_id)
         iret=nf_inq_varid(ncid,"u-vel",u_id)
         iret=nf_inq_varid(ncid,"v-vel",v_id)
         iret=nf_inq_varid(ncid,"minuvelt",minut_id)
         iret=nf_inq_varid(ncid,"maxuvelt",maxut_id)
         iret=nf_inq_varid(ncid,"minvvelt",minvt_id)
         iret=nf_inq_varid(ncid,"maxvvelt",maxvt_id)
         iret=nf_inq_varid(ncid,"minuveln",minun_id)
         iret=nf_inq_varid(ncid,"maxuveln",maxun_id)
         iret=nf_inq_varid(ncid,"minvveln",minvn_id)
         iret=nf_inq_varid(ncid,"maxvveln",maxvn_id)
         iret=nf_inq_varid(ncid,"minunodet",minunodet_id)
         iret=nf_inq_varid(ncid,"maxunodet",maxunodet_id)
         iret=nf_inq_varid(ncid,"minvnodet",minvnodet_id)
         iret=nf_inq_varid(ncid,"maxvnodet",maxvnodet_id)
         iret=nf_inq_varid(ncid,"minutimen",minutimen_id)
         iret=nf_inq_varid(ncid,"maxutimen",maxutimen_id)
         iret=nf_inq_varid(ncid,"minvtimen",minvtimen_id)
         iret=nf_inq_varid(ncid,"maxvtimen",maxvtimen_id)
      ENDIF
!     ----------------------------------------
!     Allocating memory for min and max arrays
!     ----------------------------------------
      IF(.NOT. ALLOCATED(minut)  )ALLOCATE( minut(timenc_len))
      IF(.NOT. ALLOCATED(maxut)  )ALLOCATE( maxut(timenc_len))
      IF(.NOT. ALLOCATED(minun62))ALLOCATE( minun62(node_len))
      IF(.NOT. ALLOCATED(maxun62))ALLOCATE( maxun62(node_len))

      IF(.NOT. ALLOCATED(minvt)  )ALLOCATE( minvt(timenc_len))
      IF(.NOT. ALLOCATED(maxvt)  )ALLOCATE( maxvt(timenc_len))
      IF(.NOT. ALLOCATED(minvn62))ALLOCATE( minvn62(node_len))
      IF(.NOT. ALLOCATED(maxvn62))ALLOCATE( maxvn62(node_len))


      IF(.NOT. ALLOCATED(minunodet)  )ALLOCATE( minunodet  (timenc_len))
      IF(.NOT. ALLOCATED(maxunodet)  )ALLOCATE( maxunodet  (timenc_len))
      IF(.NOT. ALLOCATED(minutimen62))ALLOCATE( minutimen62(  node_len))
      IF(.NOT. ALLOCATED(maxutimen62))ALLOCATE( maxutimen62(  node_len))

      IF(.NOT. ALLOCATED(minvnodet)  )ALLOCATE( minvnodet  (timenc_len))
      IF(.NOT. ALLOCATED(maxvnodet)  )ALLOCATE( maxvnodet  (timenc_len))
      IF(.NOT. ALLOCATED(minvtimen62))ALLOCATE( minvtimen62(  node_len))
      IF(.NOT. ALLOCATED(maxvtimen62))ALLOCATE( maxvtimen62(  node_len))

      IF(IHOT.NE.0.AND.FIRST62.EQ.1.and.ncmode.ne.0.and.nout.eq.3)THEN
         iret = nf_get_var_double(ncid, minun_id, minun62)
         CALL check_err(iret)
         iret = nf_get_var_double(ncid, maxun_id, maxun62)
         CALL check_err(iret)
         iret = nf_get_var_double(ncid, minvn_id, minvn62)
         CALL check_err(iret)
         iret = nf_get_var_double(ncid, maxvn_id, maxvn62)
         CALL check_err(iret)
         
         iret = nf_get_var_double(ncid, minutimen_id, minutimen62)
         CALL check_err(iret)
         iret = nf_get_var_double(ncid, maxutimen_id, maxutimen62)
         CALL check_err(iret)
         iret = nf_get_var_double(ncid, minvtimen_id, minvtimen62)
         CALL check_err(iret)
         iret = nf_get_var_double(ncid, maxvtimen_id, maxvtimen62)
         CALL check_err(iret)
      ENDIF
      
      FIRST62=0
      NREC = NREC + 1
      
!     ----------
!     Store time
!     ----------
      starti(1)=nrec
      counti(1)=timenc_len
      iret = nf_put_vara_DOUBLE(ncid, timenc_id, starti,counti, timenc)
      CALL check_err(iret)
!     ---------------------------
!     Store u
!     ---------------------------
      count(1)=node_len
      count(2)=timenc_len
      start(1)=1
      start(2)=nrec
      iret = nf_put_vara_DOUBLE(ncid, u_id, start, count, u)
      CALL check_err(iret)
!     ---------------------------
!     Store v
!     ---------------------------
      count(1)=node_len
      count(2)=timenc_len
      start(1)=1
      start(2)=nrec
      iret = nf_put_vara_DOUBLE(ncid, v_id, start, count, v)
      CALL check_err(iret)

!     ------------------------------------
!     Calculating u minimum for each node
!     ------------------------------------
      do n=1,node_len
 	 if(u (n, 1).lt.minun62(n) ) then
 	    minun62     (n) = u (n,1)
 	    minutimen62 (n) = timenc(1)
 	 endif
      end do
!     ------------------------------------
!     Calculating u minimum for each time
!     ------------------------------------
      minut (1) = 99999.
      do n=1,node_len
         if(u (n,1).gt.-99999.and.u(n,1).lt. minut (1) )  then
            minut (1)    = u (n,1)
            minunodet(1) = n
         endif
      end do

      iret = nf_put_vara_DOUBLE(ncid, minut_id,(/nrec/),(/1/),minut)
      CALL check_err(iret)
      iret=nf_put_vara_INT(ncid, minunodet_id,(/nrec/),(/1/),minunodet)
      CALL check_err(iret)
      
!     ------------------------------------
!     Calculating u maximum for each node
!     ------------------------------------
      do n=1,node_len
 	 if(u (n,1).gt.maxun62(n) ) then
 	    maxun62     (n) = u (n,1)
            maxutimen62 (n) = timenc(1)
 	 endif
      end do
!     ------------------------------------
!     Calculating u maximum for each time
!     ------------------------------------
      maxut (1) = -9999.
      do n=1,node_len
 	 if(u(n,1).gt.maxut(1))  then
 	    maxut (1)    = u (n,1)
 	    maxunodet(1) = n 		    
         endif
      end do
      iret = nf_put_vara_DOUBLE(ncid,maxut_id,(/nrec/),(/1/),maxut)
      CALL check_err(iret)
      iret=nf_put_vara_INT(ncid,maxunodet_id,(/nrec/),(/1/),maxunodet)
      CALL check_err(iret)
!     ---------------------------------
!     Calculating v minimum each node
!     --------------------------------
      do n=1,node_len
         if(v (n,1).lt.minvn62(n) ) then
 	    minvn62     (n) = v (n,1)
	    minvtimen62 (n) = timenc(1)
         endif
      end do
!     ------------------------------------
!     Calculating v minimum for each time
!     ------------------------------------
      minvt (1) = 99999.
      do n=1,node_len
         if(v (n,1).gt.-99999.and.v(n,1).lt. minvt (1) ) then
            minvt (1)    = v (n,1)
            minvnodet(1) = n
         endif
      end do

      iret = nf_put_vara_DOUBLE(ncid, minvt_id,(/nrec/),(/1/),minvt)
      CALL check_err(iret)
      iret=nf_put_vara_INT(ncid, minvnodet_id,(/nrec/),(/1/),minvnodet)
      CALL check_err(iret)
!     -----------------------------------
!     Calculating v maximum in each node
!     -----------------------------------
      do n=1,node_len
         if(v (n,1).gt.maxvn62(n) ) then
            maxvn62    (n) = v (n,1)
            maxvtimen62(n) = timenc(1)
         endif
      end do
!     ------------------------------------
!     Calculating v maximum for each time
!     ------------------------------------
      maxvt (1) = -9999.
      do n=1,node_len
         if(v(n,1).gt.maxvt(1)) then
 	    maxvt (1)    = v (n,1)
 	    maxvnodet(1) = n 		    
         endif
      end do
      iret = nf_put_vara_DOUBLE(ncid,maxvt_id,(/nrec/),(/1/),maxvt)
      CALL check_err(iret)
      iret=nf_put_vara_INT(ncid,maxvnodet_id,(/nrec/),(/1/),maxvnodet)
      CALL check_err(iret)


!     ---------------------------
!     Closing fort.X2 netCDF file
!     ---------------------------
      IF(NREC.EQ.NDATA) THEN
         iret = nf_put_var_double(ncid, minun_id, minun62)
         CALL check_err(iret)
         iret = nf_put_var_double(ncid, maxun_id, maxun62)
         CALL check_err(iret)
         iret = nf_put_var_double(ncid, minvn_id, minvn62)
         CALL check_err(iret)
         iret = nf_put_var_double(ncid, maxvn_id, maxvn62)
         CALL check_err(iret)
         
         iret = nf_put_var_double(ncid, minutimen_id, minutimen62)
         CALL check_err(iret)
         iret = nf_put_var_double(ncid, maxutimen_id, maxutimen62)
         CALL check_err(iret)
         iret = nf_put_var_double(ncid, minvtimen_id, minvtimen62)
         CALL check_err(iret)
         iret = nf_put_var_double(ncid, maxvtimen_id, maxvtimen62)
         CALL check_err(iret)

      ENDIF

!     ---------------------------
!     Closing netCDF file
!     ---------------------------
      WRITE(6,*) "Closing ", FILENAME, " netCDF file ... "
      WRITE(6,*) "-----------------------------------------------------"
      
9999  iret = nf_close(ncid)
      CALL check_err(iret)
      WRITE(6,*)'end subroutine writenc62'

      end subroutine writenc62     

!=======================================================================
!     Subroutine writenc71
!     Writes a netCDF file equivalent to ADCIRC's fort.71 ascii file
!
!     station  : number of stations
!     namelen  : length of station's name
!     time     : number of time records
!
!     Variables :
!     time      : time
!     station_name : station's name
!     zeta      : pressure 
!     x         : x-coordinate or longitude
!     y         : y-coordinate or latitude
!     element   : elements in grid
!     nvdll     : number of nodes in elevation boundary segment
!     nbdv      : node numbers on elevation specified boundary segment
!     nvell     : number of nodes in normal flow specified boundary segment
!     nbvv      : node numbers on normal flow boundary segment
!     ibtype    : boundary type
!     depth     : depth in each node
!     mindepth  : minimum depth in grid
!     maxdepth  : maximum depth in grid
!     minprest  : minimum pressure for each output time
!     maxprest  : maximum pressure for each output time
!     minpresn  : minimum pressure at each station
!     maxpresn  : maximum pressure at each station
!     minpnodet : node of minimum pressure at each output time
!     maxpnodet : node of maximum pressure at each output time
!     minptimen : time of minimum pressure at each station
!     maxptimen : time of maximum pressure at each station
!
!     Attributes:
!     model             : ADCIRC
!     version           : model version number
!     grid_type         : triangular
!     rundes            : description of simulation
!     agrid             : grid ID 
!     title             : description of the dataset
!     institution       : where the original data was produced
!     source            : method of production of the original data
!     history           : provides audit trail for modifications to the original data
!     references        : published/web-based references that describe data
!     comments          : miscellaneous information about the data
!     contact           : user's institution/ e-mail address
!     creation_date     : date the data was created
!     modification_date : date the data was modified
!     Parameters from fort.15 input file 
!
!     Revision history:
!
!     Date            Programmer                  Description of change
!     ----            ----------                  ----------------------
!     03/30/07  Cristina Forbes, PSGS @ UNC-IMS   Wrote original code
!     03/30/08  Cristina Forbes, PSGS @ UNC-IMS   Modified code for 
!                                                 globalio & hotstart 
!                                                 from binary files
!=======================================================================
      subroutine writenc71 (nout,fileunit,timenc_len, node_len,timenc,
     &                      zeta,nt, NDATA, rundes,runid, agrid, irtype,
     &                      dtdp, ihot, ics, nolibf,  nolifa,
     &                      nolica, nolicat, nwp,
     &                      ncor, ntip, nws, nramp, tau0,
     &                      statim, reftim, rnday, dramp,
     &                      a00, b00, c00, h0,
     &                      slam0, sfea0, ffactor,esl, cori,
     &                      ntif, nbfr, elem_len, nface_len, nm,
     &                      x,y,depth,
     &                      nope_len, neta_len, nbdv, 
     &                      nbou_len, nvel_len, nbvv,
     &                      ncmode,
     &                      title, institution, source, history,
     &                      references, comments, host, convention,
     &                      contact,
     &                      base_date)

!--------------------------------------------
!     Force explicit declaration of all variables
!--------------------------------------------
      IMPLICIT NONE
!     -------------------
!     netCDF header file
!     -------------------
!      include 'netcdf.inc'
!     -------------------
!     Error status return
!     -------------------
      INTEGER  iret
!     ---------
!     netCDF id
!     ---------
      INTEGER  ncid
!     -------------
!     Dimension ids
!     -------------
      INTEGER  node_dim
      INTEGER  timenc_dim
      INTEGER  elem_dim
      INTEGER  nface_dim
      INTEGER  nope_dim
      INTEGER  neta_dim
      INTEGER  nbou_dim
      INTEGER  nvel_dim
      INTEGER  slen_dim
      INTEGER  scode_dim
!     -----------------
!     Dimension lengths
!     -----------------
      INTEGER  node_len
      INTEGER  timenc_len
      INTEGER  elem_len
      INTEGER  nface_len
      INTEGER  nope_len
      INTEGER  neta_len
      INTEGER  nbou_len
      INTEGER  nvel_len
!     ------------
!     Variable ids
!     ------------
      INTEGER  timenc_id
      INTEGER  zeta_id
      INTEGER  minzetat_id
      INTEGER  maxzetat_id
      INTEGER  minzetan_id
      INTEGER  maxzetan_id
      INTEGER  minznodet_id
      INTEGER  maxznodet_id
      INTEGER  minztimen_id
      INTEGER  maxztimen_id
      INTEGER  station_id
      INTEGER  scode_id
      INTEGER  X_id
      INTEGER  Y_id
      INTEGER  DEPTH_id
      INTEGER  ELE_id
      INTEGER  mindepth_id
      INTEGER  maxdepth_id
      INTEGER  nbdv_id
      INTEGER  nbvv_id
!     ---------------------------------------------
!     Rank (number of dimensions) for each variable
!     ---------------------------------------------
      INTEGER  timenc_rank
      INTEGER  zeta_rank
      INTEGER  station_rank
      INTEGER  scode_rank
      INTEGER  minzetat_rank
      INTEGER  maxzetat_rank
      INTEGER  minzetan_rank
      INTEGER  maxzetan_rank
      INTEGER  minznodet_rank
      INTEGER  maxznodet_rank
      INTEGER  minztimen_rank
      INTEGER  maxztimen_rank
      INTEGER  depth_rank
      INTEGER  ele_rank
      INTEGER  mindepth_rank
      INTEGER  maxdepth_rank
      INTEGER  x_rank
      INTEGER  y_rank
      INTEGER  nvdll_rank
      INTEGER  nbdv_rank
      INTEGER  nbvv_rank
      PARAMETER (timenc_rank = 1)
      PARAMETER (station_rank = 2)
      PARAMETER (scode_rank = 2)
      PARAMETER (zeta_rank = 2)
      PARAMETER (minzetat_rank = 1)
      PARAMETER (maxzetat_rank = 1)
      PARAMETER (minzetan_rank = 1)
      PARAMETER (maxzetan_rank = 1)
      PARAMETER (minznodet_rank = 1)
      PARAMETER (maxznodet_rank = 1)
      PARAMETER (minztimen_rank = 1)
      PARAMETER (maxztimen_rank = 1)
      PARAMETER (depth_rank = 1)
      PARAMETER (ele_rank = 2)
      PARAMETER (mindepth_rank = 0)
      PARAMETER (maxdepth_rank = 0)
      PARAMETER (x_rank = 1)
      PARAMETER (y_rank = 1)
      PARAMETER (nvdll_rank = 1)
      PARAMETER (nbdv_rank = 2)
      PARAMETER (nbvv_rank = 2)
!     ---------------
!     Variable shapes
!     ---------------
      INTEGER  timenc_dims(timenc_rank)
      INTEGER  zeta_dims(zeta_rank)
      INTEGER  station_dims  (station_rank)
      INTEGER  scode_dims    (scode_rank)
      INTEGER  minzetat_dims(minzetat_rank)
      INTEGER  maxzetat_dims(maxzetat_rank)
      INTEGER  minzetan_dims(minzetan_rank)
      INTEGER  maxzetan_dims(maxzetan_rank)
      INTEGER  minznodet_dims(minznodet_rank)
      INTEGER  maxznodet_dims(maxznodet_rank)
      INTEGER  minztimen_dims(minztimen_rank)
      INTEGER  maxztimen_dims(maxztimen_rank)
      INTEGER  x_dims(x_rank)
      INTEGER  y_dims(y_rank)
      INTEGER  depth_dims(depth_rank)
      INTEGER  ele_dims(ele_rank)
      INTEGER  nvdll_dims(nvdll_rank)
      INTEGER  nbdv_dims(nbdv_rank)
      INTEGER  nbvv_dims(nbvv_rank)
!     --------------
!     Data variables
!     --------------
      CHARACTER*5 FILEFORT
      CHARACTER*10 FILENAME
      CHARACTER*3 FILEEXT
      INTEGER   t,ntif,nbfr,nt
      INTEGER   irtype, ihot, ics, nolibf, nolifa, nolica, nolicat, nwp
      INTEGER   ntip, nws, nramp,ncor
      INTEGER   nvel, neta, nout
      INTEGER   element (nface_len , elem_len)
      INTEGER   nbdv (nope_len, neta_len)
      INTEGER   nbvv (nbou_len, nvel_len)
      INTEGER   nm (elem_len, nface_len)
      INTEGER   fileunit
      REAL*8    timenc(timenc_len)
      REAL(sz)  x(node_len)
      REAL(sz)  y(node_len)
      REAL(sz)  depth(node_len)
      REAL(sz)  zeta (node_len, timenc_len)
      REAL(sz), ALLOCATABLE, SAVE ::  minzetat(:)
      REAL(sz), ALLOCATABLE, SAVE ::  maxzetat(:)
      REAL(sz), ALLOCATABLE, SAVE ::  minzetan71(:)
      REAL(sz), ALLOCATABLE, SAVE ::  maxzetan71(:)
      INTEGER , ALLOCATABLE, SAVE ::  minznodet(:)
      INTEGER , ALLOCATABLE, SAVE ::  maxznodet(:)
      REAL(sz), ALLOCATABLE, SAVE ::  minztimen71(:)
      REAL(sz), ALLOCATABLE, SAVE ::  maxztimen71(:)
      real*8    mindepth
      real*8    maxdepth
      REAL*8    dtdp, statim,reftim, rnday, dramp, tau0
      REAL*8    a00, b00, c00,h0
      REAL*8    ffactor,esl, cori
      REAL(sz)  slam0, sfea0
      REAL(sz)  slam0deg, sfea0deg
      CHARACTER(80) rundes, runid, agrid
      INTEGER, SAVE  :: FIRST71 = 1 
      INTEGER NDATA
      INTEGER, SAVE  :: NREC = 0
      INTEGER ncmode
!     -----------------
!     Attribute vectors
!     -----------------
      INTEGER intval(1),i,j
      INTEGER counti(1), starti(1), n
      INTEGER count(2), start(2)
      REAL    realval(1)
!     -----------------
!     Attribute general
!     -----------------
      CHARACTER(*) :: title, institution, source, history
      CHARACTER(*) :: references, comments, host, convention, contact
      CHARACTER(*) :: base_date
!     -----------------
!      date_string variables for time attribute
!     -----------------
      character date_string*40
      character now_date*8
      character big_ben*10
      character zone*5
      integer   values(8)
      real*8    xmin, xmax
!     ---------------------------------
!      inquiring and reading variables
!     ---------------------------------

      INTEGER NDIMS, NVARS, NGATTS, UNLIMDIMID, RECID,NRECS
      CHARACTER*4 recname

!     -------------------------
!     Save dims, ids, variables
!     -------------------------
      save node_dim, elem_dim, nface_dim
      save nope_dim, neta_dim, nbou_dim, nvel_dim
      save zeta_id, X_id, Y_id
      save minzetat_id, maxzetat_id, minzetan_id, maxzetan_id
      save minznodet_id, maxznodet_id, minztimen_id, maxztimen_id
      save DEPTH_id,ELE_id, mindepth_id, maxdepth_id, nbdv_id, nbvv_id
      save timenc_dim, timenc_id
      save station_id, scode_id
      save mindepth, maxdepth
!     ---------------------------
!     Opening fort.63 netCDF file
!     ---------------------------
      WRITE(6,*) "-----------------------------------------------------"

!      NREC = NREC + 1

      IF(FIRST71.EQ.1.and.ncmode.NE.0) THEN
         
         ALLOCATE ( minzetan71  (node_len) )
         ALLOCATE ( maxzetan71  (node_len) )
         ALLOCATE ( minztimen71 (node_len) )
         ALLOCATE ( maxztimen71 (node_len) )
         
         do n=1,node_len
            minzetan71  (n) =  9999.d0
            maxzetan71  (n) = -9999.d0
            minztimen71 (n) =  9999.d0
            maxztimen71 (n) = -9999.d0
         end do
       ENDIF

!     --------------------------------
!     Switch order in array for NETCDF
!     --------------------------------
      do i=1, elem_len
      do j=1, nface_len
         element(j,i) = nm (i,j)
      end do
      end do
!     -----------------
!     Enter define mode
!     -----------------
         FILEFORT = 'fort.'
         FILEEXT = '.nc'
!        WRITE(6,'(A5,I2,A3)')FILEFORT,FILEUNIT,FILEEXT
         WRITE(FILENAME,'(A5,I2,A3)')FILEFORT,FILEUNIT,FILEEXT
         WRITE(6,*) "Opening ",FILENAME,  " netCDF file ... "
                  
       IF((FIRST71.EQ.1.AND.IHOT.EQ.0.AND.NCMODE.NE.2).OR.
     &    (FIRST71.EQ.1.AND.IHOT.NE.0.AND.NCMODE.NE.2.AND.
     &     nout.EQ.-3))THEN
     
           iret = nf_create(FILENAME, NF_CLOBBER, ncid)
           CALL check_err(iret)
           FIRST71=0
         
!     --------------------------------
!     Convert to degress for NETCDF I/O
!     --------------------------------
           SLAM0DEG=SLAM0
           SFEA0DEG=SFEA0
           IF (NCMODE.EQ.1) THEN
               SLAM0DEG=SLAM0*RAD2DEG
               SFEA0DEG=SFEA0*RAD2DEG
           ENDIF

!     =================
!     Define dimensions
!     =================
      iret = nf_def_dim(ncid, 'station', node_len, node_dim)
             CALL check_err(iret)
      iret = nf_def_dim(ncid, 'namelen', SNAMLEN, slen_dim)
             CALL check_err(iret)
      iret = nf_def_dim(ncid, 'time', nf_unlimited, timenc_dim)
             CALL check_err(iret)

!     ================
!     Define variables
!     ================
!     ----------------
!     Define time
!     ----------------

      timenc_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'time',NF_DOUBLE,timenc_rank, timenc_dims, 
     &                  timenc_id)
      CALL check_err(iret)
     
!     -------------------
!     Define X
!     -------------------
      X_dims(1) = node_dim
      iret = nf_def_var(ncid, 'x', NF_DOUBLE, X_rank, X_dims, X_id)
      CALL check_err(iret)
!     -------------------
!     Define Y coordinate
!     --------------------
      Y_dims(1) = node_dim
      iret = nf_def_var(ncid, 'y', NF_DOUBLE, Y_rank, Y_dims, Y_id)
      CALL check_err(iret)
!     ---------------------
!     Define stations name
!     --------------------
      station_dims(1) = slen_dim
      station_dims(2) = node_dim
      iret = nf_def_var(ncid, 'station_name', NF_CHAR, station_rank, 
     1      station_dims, station_id)
      CALL check_err(iret)
      
!     ----------------
!     Define zeta
!     ----------------
      zeta_dims(1) = node_dim
      zeta_dims(2) = timenc_dim
      
      iret = nf_def_var(ncid, 'pressure', NF_DOUBLE, zeta_rank, zeta_dims, 
     &                  ZETA_id)
      
      CALL check_err(iret)
!     ----------------------------
!     Define minimum zeta at times
!     ----------------------------
      minzetat_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'minprest', NF_DOUBLE, minzetat_rank, 
     &                  minzetat_dims, minzetat_id)
      CALL check_err(iret)
!     ----------------------------
!     Define maximum zeta at times
!     ----------------------------
      maxzetat_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'maxprest', NF_DOUBLE, maxzetat_rank, 
     &                  maxzetat_dims, maxzetat_id)
      CALL check_err(iret)

!     ----------------------------
!     Define minimum zeta at nodes
!     ----------------------------
      minzetan_dims(1) = node_dim
       iret = nf_def_var(ncid, 'minpresn', NF_DOUBLE, minzetan_rank, 
     &                  minzetan_dims, minzetan_id)
      CALL check_err(iret)
!     ----------------------------
!     Define maximum zeta at nodes
!     ----------------------------
      maxzetan_dims(1) = node_dim
      iret = nf_def_var(ncid, 'maxpresn', NF_DOUBLE, maxzetan_rank, 
     &                  maxzetan_dims, maxzetan_id)
      CALL check_err(iret)
      
!     ----------------------------
!     Define minimum zeta nodes at times
!     ----------------------------
      minznodet_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'minpnodet', NF_INT, minznodet_rank, 
     &                  minznodet_dims, minznodet_id)
      CALL check_err(iret)
!     ----------------------------
!     Define maximum zeta nodes at times
!     ----------------------------
      maxznodet_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'maxpnodet', NF_INT, maxznodet_rank, 
     &                  maxznodet_dims, maxznodet_id)
      CALL check_err(iret)
!     ----------------------------
!     Define minimum zeta nodes at nodes
!     ----------------------------
      minztimen_dims(1) = node_dim
      iret = nf_def_var(ncid, 'minptimen', NF_DOUBLE, minztimen_rank, 
     &                  minztimen_dims, minztimen_id)
      CALL check_err(iret)
!     ----------------------------
!     Define maximum zeta nodes at nodes
!     ----------------------------
      maxztimen_dims(1) = node_dim
      iret = nf_def_var(ncid, 'maxptimen', NF_DOUBLE, maxztimen_rank, 
     &                  maxztimen_dims, maxztimen_id)
      CALL check_err(iret)       
!     -----------------
!     Assign attributes
!     -----------------
!     ----------------------------
!     Define longitude attributes
!     ----------------------------
!     ----------------------------
!     setermine coordinates
!     ----------------------------
      if(ics.eq.2) then
      
      	iret = nf_put_att_text(ncid,X_id,'long_name',9,'longitude')
      	       CALL check_err(iret)
      	iret = nf_put_att_text(ncid, X_id, 'standard_name', 9, 
     &                      'longitude')
      	       CALL check_err(iret)
      	iret = nf_put_att_text(ncid, X_id, 'units', 12, VAR_DEGE)
      	       CALL check_err(iret)
      	iret = nf_put_att_text(ncid, X_id, 'positive', 4, VAR_EAST)
      	       CALL check_err(iret)

      	iret = nf_put_att_double(ncid, X_id, 'missing_values',NF_DOUBLE, 
     &                       1, doubleval)
      	       CALL check_err(iret)

      	iret = nf_put_att_double(ncid, X_id, '_FillValue', NF_DOUBLE, 1, 
     &                       doubleval)
      	       CALL check_err(iret)      
      
!     ----------------------------
!     Define latitude attributes
!     ----------------------------
      	iret = nf_put_att_text(ncid, Y_id, 'long_name', 8, 'latitude')
      	       CALL check_err(iret)
      	iret = nf_put_att_text(ncid, Y_id, 'standard_name', 8, 
     &                       'latitude')
      	       CALL check_err(iret)
      	iret = nf_put_att_text(ncid, Y_id, 'units', 13, VAR_DEGN)
      	       CALL check_err(iret)
      	iret = nf_put_att_text(ncid, Y_id, 'positive', 5, VAR_NORTH)
      	       CALL check_err(iret)
       
 
      	iret = nf_put_att_double(ncid, Y_id, 'missing_values', 
     &         NF_DOUBLE, 1, doubleval)
      	       CALL check_err(iret)
 
      	iret = nf_put_att_double(ncid, Y_id, '_FillValue', NF_DOUBLE, 1, 
     &         doubleval)
      	       CALL check_err(iret)

       else
      
!     ------------------------------
!     Define x-coordinate attributes
!     ------------------------------
      	iret = nf_put_att_text(ncid, X_id, 'long_name', 22, 
     &       'Cartesian coordinate x')
      	       CALL check_err(iret)
      	iret = nf_put_att_text(ncid, X_id, 'standard_name',12, 
     &                      'x-coordinate')
      	       CALL check_err(iret)
        if (G.LT.11.D0) then
      	    iret = nf_put_att_text(ncid, X_id, 'units', 6,VAR_METS)
      	else
      	    iret = nf_put_att_text(ncid, X_id, 'units', 6,VAR_FEET)
      	endif
      	CALL check_err(iret)
      	iret = nf_put_att_text(ncid, X_id, 'positive', 5, VAR_RG)
      	CALL check_err(iret)

      	iret = nf_put_att_double(ncid, X_id, 'missing_values',NF_DOUBLE,
     &                        1, doubleval)
      	CALL check_err(iret)

      	iret = nf_put_att_double(ncid, X_id, '_FillValue', NF_DOUBLE, 1, 
     &                       doubleval)
      	CALL check_err(iret)
!     ------------------------------
!     Define y-coordinate attributes
!     ------------------------------     
      	iret = nf_put_att_text(ncid, Y_id, 'long_name', 22, 
     &       'Cartesian coordinate y')
      	CALL check_err(iret)
      	iret = nf_put_att_text(ncid, Y_id, 'standard_name', 12, 
     &                       'y-coordinate')
      	CALL check_err(iret)
!     ----------------------
!     Determine coordinates
!     ---------------------    
        if (G.LT.11.D0) then
     	   iret = nf_put_att_text(ncid, Y_id, 'units', 6,VAR_METS)
     	else
     	   iret = nf_put_att_text(ncid, Y_id, 'units', 6,VAR_FEET)
     	endif
      	CALL check_err(iret)
      	iret = nf_put_att_text(ncid, Y_id, 'positive', 5, VAR_UP)
      	CALL check_err(iret)
 
      	iret = nf_put_att_double(ncid, Y_id, 'missing_values',NF_DOUBLE, 
     &       1, doubleval)
      	CALL check_err(iret)
 
      	iret = nf_put_att_double(ncid, Y_id, '_FillValue', NF_DOUBLE, 1, 
     &       doubleval)
      	       CALL check_err(iret)
      
      endif
      
!     -------------------------------------------------------------
!     Define max and min in time sea surface pressure attributes 
!     -------------------------------------------------------------
      iret = nf_put_att_text(ncid, minzetat_id, 'long_name', 65, 
     &             'minimum station atmospheric surface pressure at each 
     &output time')
            CALL check_err(iret)
      iret = nf_put_att_text(ncid, minzetat_id, 'units', 15,  'meters of
     & water')
            CALL check_err(iret)
      
      iret = nf_put_att_text(ncid, maxzetat_id, 'long_name', 65, 
     &             'maximum station atmospheric surface pressure at each 
     &output time')
            CALL check_err(iret)
      iret = nf_put_att_text(ncid, maxzetat_id, 'units', 15, VAR71_UN)
            CALL check_err(iret)
!     -------------------------------------------------------------
!     Define max and min in nodes sea surface pressure attributes
!     -------------------------------------------------------------
      iret = nf_put_att_text(ncid, minzetan_id, 'long_name', 52, 
     &      'minimum atmospheric surface pressure at each station')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, minzetan_id, 'units', 15,  VAR71_UN)
      CALL check_err(iret)
      
      iret = nf_put_att_text(ncid, maxzetan_id, 'long_name', 52, 
     &      'maximum atmospheric surface pressure at each station')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, maxzetan_id, 'units', 15,  VAR71_UN)
      CALL check_err(iret)
!     -------------------------------------------------------------
!     Define max and min in time sea surface elevation attributes 
!     -------------------------------------------------------------
      iret = nf_put_att_text(ncid, minznodet_id, 'long_name', 68, 
     &  'station of minimum atmospheric surface pressure at each output
     & time')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, minznodet_id, 'units',14, VAR_NOND)
      CALL check_err(iret)
      
      iret = nf_put_att_text(ncid, maxznodet_id, 'long_name', 68, 
     &  'station of maximum atmospheric surface pressure at each output
     & time')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, maxznodet_id, 'units', 14,VAR_NOND) 
      CALL check_err(iret)
!     -------------------------------------------------------------
!     Define max and min in nodes sea surface elevation attributes
!     -------------------------------------------------------------
      iret = nf_put_att_text(ncid, minztimen_id, 'long_name', 60, 
     &   'time of minimum atmospheric surface pressure at each station')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, minztimen_id, 'units', 4, VAR_SECS)
      CALL check_err(iret)
      
      iret = nf_put_att_text(ncid, maxztimen_id, 'long_name', 60, 
     &  'time of maximum atmospheric surface pressure at each station')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, maxztimen_id, 'units', 4, VAR_SECS)
      CALL check_err(iret)
      
!     -------------------------------------------------------------
!     Define station names and code attributes
!     -------------------------------------------------------------

      iret = nf_put_att_text(ncid,station_id,'long_name',12,
     & 'station name')
      CALL check_err(iret)
!    -----------------------
!     Define time attributes
!    -----------------------
      iret = nf_put_att_text(ncid,timenc_id,'long_name',10,'model time')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid,timenc_id, 'standard_name', 4, 'time')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, timenc_id, 'units', 23, VAR_TUN)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, timenc_id, 'base_date',
     .       (LEN(TRIM(base_date(1:INDEX(base_date,"!")-1)))), base_date)
      CALL check_err(iret)

!    --------------------------------------------
!     Define sea surface pressure attributes
!    --------------------------------------------
      
      iret = nf_put_att_text(ncid, zeta_id,'long_name',36, VAR71_LN)
     &                   
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, zeta_id, 'standard_name', 28, 
     &                      VAR71_ST)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, zeta_id, 'units', 15, VAR71_UN)
      CALL check_err(iret)
!      doubleval(1) = -99999.
      iret = nf_put_att_DOUBLE(ncid,zeta_id,'missing_values',NF_DOUBLE,
     &1, doubleval)
      CALL check_err(iret)
!      doubleval(1) = -99999.
      iret = nf_put_att_DOUBLE(ncid, zeta_id, '_FillValue', NF_DOUBLE, 
     &                      1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, zeta_id,'positive', 2, VAR_UP)
      CALL check_err(iret)

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
      intval(1) = 1
      iret = nf_put_att_int(ncid, NCGLOBAL, 'irtype', NF_INT, 1, irtype)
      CALL check_err(iret)
      intval(1) = 0
      iret = nf_put_att_int(ncid, NCGLOBAL, 'ihot', NF_INT, 1, ihot)
      CALL check_err(iret)
      intval(1) = 1
      iret = nf_put_att_int(ncid, NCGLOBAL, 'ics', NF_INT, 1, ics)
      CALL check_err(iret)
      intval(1) = 1
      iret = nf_put_att_int(ncid, NCGLOBAL, 'nolibf', NF_INT, 1, nolibf)
      CALL check_err(iret)
      intval(1) = 1
      iret = nf_put_att_int(ncid, NCGLOBAL, 'nolifa', NF_INT, 1, nolifa)
      CALL check_err(iret)
      intval(1) = 1
      iret = nf_put_att_int(ncid, NCGLOBAL, 'nolica', NF_INT, 1, nolica)
      CALL check_err(iret)
      intval(1) = 1
      iret = nf_put_att_int(ncid, NCGLOBAL, 'nolicat', NF_INT, 1, 
     &                      nolicat)
      CALL check_err(iret)
      intval(1) = 0
      iret = nf_put_att_int(ncid, NCGLOBAL, 'nwp', NF_INT, 1, nwp)
      CALL check_err(iret)
      intval(1) = 0
      iret = nf_put_att_int(ncid, NCGLOBAL, 'ncor', NF_INT, 1, ncor)
      CALL check_err(iret)
      intval(1) = 0
      iret = nf_put_att_int(ncid, NCGLOBAL, 'ntip', NF_INT, 1, ntip)
      CALL check_err(iret)
      intval(1) = 0
      iret = nf_put_att_int(ncid, NCGLOBAL, 'nws', NF_INT, 1, nws)
      CALL check_err(iret)
      intval(1) = 1
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
      intval(1) = 5
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
      iret = nf_put_att_double(ncid, NCGLOBAL, 'ffactor', NF_DOUBLE, 1, 
     &                         ffactor)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, NCGLOBAL, 'esl', NF_DOUBLE, 1, 
     &                         esl)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, NCGLOBAL, 'cori', NF_DOUBLE, 1, 
     &                         cori)
      CALL check_err(iret)
      iret = nf_put_att_int(ncid, NCGLOBAL, 'ntif', NF_INT, 1, ntif)
      CALL check_err(iret)
      iret = nf_put_att_int(ncid, NCGLOBAL, 'nbfr', NF_INT, 1, nbfr)
      CALL check_err(iret)
      
!     ======================
!     Leave define mode
!     ======================
      iret = nf_enddef(ncid)
      CALL check_err(iret)
!     ======================
!     Store static variables
!     ======================
!     -------
!     Store Y
!     -------
      iret = nf_put_var_double(ncid, Y_id, Y)
      CALL check_err(iret)
!     -------
!     Store X
!     -------
      iret = nf_put_var_double(ncid, X_id, X)
      CALL check_err(iret)

!     --------------------
!     Store station name
!     --------------------
 
      do i=1,node_len
         start(1)=1
         start(2)=i
         count(1)= LEN(TRIM(STATNAMEM(i,SNAMLEN)))
         count(2)=1
         iret = nf_put_vara_text(ncid,station_id,start,count,
     &   STATNAMEM(i,SNAMLEN) )
         CALL check_err(iret)
      end do
!     =====================
!     Open old netcdf file
!     =====================
      ELSE
          iret = NF_OPEN(FILENAME, NF_WRITE, ncid)
          CALL check_err(iret)   
!     ------------------
!     Inquire variables
!     ------------------
          iret=nf_inq_unlimdim(ncid,recid)
          CALL check_err(iret)   
          iret=nf_inq_dimlen(ncid,recid,nrec)
          CALL check_err(iret)   
      ENDIF
!     ------------------
!     Update Global Attributes
!     ------------------
      IF(IHOT.NE.0.AND.FIRST71.EQ.1) THEN
         iret = NF_REDEF (ncid)
         CALL check_err(iret)
         iret = nf_put_att_int(ncid, NCGLOBAL, 'ihot', NF_INT, 1, ihot)
         CALL check_err(iret)
         iret = nf_put_att_int(ncid, NCGLOBAL, 'nramp', NF_INT, 1,nramp)
         CALL check_err(iret)
         iret = nf_put_att_double(ncid, NCGLOBAL, 'tau0', NF_DOUBLE, 1, 
     &                         tau0)
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
         iret = nf_put_att_double(ncid, NCGLOBAL, 'ffactor', NF_DOUBLE, 1, 
     &                         ffactor)
         CALL check_err(iret)
         iret = nf_put_att_double(ncid, NCGLOBAL, 'esl', NF_DOUBLE, 1, 
     &                         esl)
         CALL check_err(iret)
         CALL date_and_time(now_date,big_ben,zone,values)
         WRITE(date_string,71) values(1),values(2),values(3),
     &         values(5),values(6),values(7),   (values(4))/60
         iret = nf_put_att_text(ncid, NF_GLOBAL,'modification_date',
     &          len_trim(date_string),date_string)
         CALL check_err(iret)
         iret = NF_ENDDEF (ncid)
         CALL check_err(iret)
       ENDIF
      IF(ncmode.eq.0) go to 9999
      IF(FIRST71.EQ.1.and.
     & ((ncmode.ne.0.and.IHOT.NE.0).or.(ncmode.eq.2.and.IHOT.eq.0).or.
     &  (ncmode.eq.2.and.IHOT.ne.0.and.nout.eq.-3)))THEN
         iret=nf_inq_varid(ncid,"time"     ,timenc_id)
         iret=nf_inq_varid(ncid,"pressure" ,zeta_id)
         iret=nf_inq_varid(ncid,"minprest" ,minzetat_id)
         iret=nf_inq_varid(ncid,"maxprest" ,maxzetat_id)
         iret=nf_inq_varid(ncid,"minpresn" ,minzetan_id)
         iret=nf_inq_varid(ncid,"maxpresn" ,maxzetan_id)
         iret=nf_inq_varid(ncid,"minpnodet",minznodet_id)
         iret=nf_inq_varid(ncid,"maxpnodet",maxznodet_id)
         iret=nf_inq_varid(ncid,"minptimen",minztimen_id)
         iret=nf_inq_varid(ncid,"maxptimen",maxztimen_id)
      ENDIF
      IF(IHOT.NE.0.AND.FIRST71.EQ.1.and.ncmode.ne.0.and.nout.eq.3)THEN
         iret = nf_get_var_double(ncid, minzetan_id, minzetan71)
         CALL check_err(iret)
         iret = nf_get_var_double(ncid, maxzetan_id, maxzetan71)
         CALL check_err(iret)
         iret = nf_get_var_double(ncid, minztimen_id, minztimen71)
         CALL check_err(iret)
         iret = nf_get_var_double(ncid, maxztimen_id, maxztimen71)
         CALL check_err(iret)
      ENDIF
      
      FIRST71=0
      NREC = NREC + 1
      
!     ----------
!     Store time
!     ----------

      starti(1)=nrec
      counti(1)=timenc_len
      iret = nf_put_vara_DOUBLE(ncid, timenc_id, starti,counti, timenc)
      CALL check_err(iret)
!     ---------------------------
!     Store zeta
!     ---------------------------
      count(1)=node_len
      count(2)=timenc_len
      start(1)=1
      start(2)=nrec
      iret = nf_put_vara_DOUBLE(ncid, zeta_id, start, count, zeta)
      CALL check_err(iret)
!     ----------------------------------------
!     Allocating memory for min and max arrays
!     ----------------------------------------
      IF(.NOT. ALLOCATED(minzetat))ALLOCATE( minzetat(timenc_len))
      IF(.NOT. ALLOCATED(maxzetat))ALLOCATE( maxzetat(timenc_len))
      IF(.NOT. ALLOCATED(minzetan71))ALLOCATE( minzetan71(node_len))
      IF(.NOT. ALLOCATED(maxzetan71))ALLOCATE( maxzetan71(node_len))

      IF(.NOT. ALLOCATED(minznodet))ALLOCATE( minznodet(timenc_len))
      IF(.NOT. ALLOCATED(maxznodet))ALLOCATE( maxznodet(timenc_len))
      IF(.NOT. ALLOCATED(minztimen71))ALLOCATE( minztimen71(  node_len))
      IF(.NOT. ALLOCATED(maxztimen71))ALLOCATE( maxztimen71(  node_len))

!     --------------------------------------
!     Calculating zeta minimum in each node
!     --------------------------------------
      do n=1,node_len
 	 if(zeta (n,1).lt.minzetan71(n) ) then
 	    minzetan71  (n) = zeta (n,1)
            minztimen71 (n) = timenc(1)
 	 endif
      end do
!     --------------------------------------
!     Calculating zeta minimum for each time
!     --------------------------------------
      minzetat (1) = 9999.d0
      do n=1,node_len
         if(zeta (n,1).gt.-9999.d0.and.zeta(n,1).lt. minzetat (1))then 
           minzetat (1) = zeta (n,1)
           minznodet(1) = n
        endif
      end do

      iret=nf_put_vara_DOUBLE(ncid, minzetat_id,(/nrec/),(/1/),minzetat)
           CALL check_err(iret)
      iret=nf_put_vara_INT(ncid, minznodet_id,(/nrec/),(/1/),minznodet)
           CALL check_err(iret)

!     --------------------------------------
!     Calculating zeta maximum in each node
!     --------------------------------------
      do n=1,node_len
         if(zeta (n,1).gt.maxzetan71(n) ) then
            maxzetan71 (n) = zeta (n,1)
            maxztimen71(n) = timenc(1)
          endif
      end do
!     --------------------------------------
!     Calculating zeta maximum for each time
!     --------------------------------------
      maxzetat (1) = -9999.d0
      do n=1,node_len
         if(zeta(n,1).gt.maxzetat(1))  then
 	    maxzetat (1) = zeta (n,1)
 	    maxznodet(1) = n
 	 endif
      end do

      iret=nf_put_vara_DOUBLE(ncid,maxzetat_id,(/nrec/),(/1/),maxzetat)
      CALL check_err(iret)
      iret=nf_put_vara_INT(ncid,maxznodet_id,(/nrec/),(/1/),maxznodet)
      CALL check_err(iret)

!     ---------------------------
!     Closing fort.X3 netCDF file
!     ---------------------------
      IF(NREC.EQ.NDATA) THEN
         iret = nf_put_var_double(ncid, minzetan_id, minzetan71)
         CALL check_err(iret)
         iret = nf_put_var_double(ncid, maxzetan_id, maxzetan71)
         CALL check_err(iret)
         iret = nf_put_var_double(ncid, minztimen_id, minztimen71)
         CALL check_err(iret)
         iret = nf_put_var_double(ncid, maxztimen_id, maxztimen71)
         CALL check_err(iret)
      ENDIF
      
      WRITE(6,*) "Closing ", FILENAME, " netCDF file ... "
      WRITE(6,*) "-----------------------------------------------------"
      
9999  iret = nf_close(ncid)
      CALL check_err(iret)
      WRITE(6,*)'end subroutine writenc71'
     
      end subroutine writenc71      


!=======================================================================
!     Subroutine writenc72
!     Writes a netCDF file equivalent to ADCIRC's fort.72 ascii file
!
!     Dimensions:
!     station  : number of stations
!     namelen  : length of station's name
!     time     : number of time records
!
!     Variables :
!     time      : time
!     station_name : station's name
!     windx      : u velocity in x or east/west direction 
!     windy      : v velocity in y or north/south direction
!     x          : x-coordinate or longitude
!     y          : y-coordinate or latitude
!     minwxt     : minimum vertically averaged u at each output time
!     maxwxt     : maximum vertically averaged u at each output time
!     minwyt     : minimum vertically averaged v at each output time
!     maxwyt     : maximum vertically averaged v at each output time
!     minwxn     : minimum vertically averaged u at each node
!     maxwxn     : maximum vertically averaged u at each node
!     minwyn     : minimum vertically averaged v at each node
!     maxwyn     : maximum vertically averaged v at each node
!     minwxnodet : node where minimum u occurs at each output time
!     maxwxnodet : node where maximum u occurs at each output time
!     minwynodet : node where minimum v occurs at each output time
!     maxwynodet : node where maximum v occurs at each output time
!     minwxtimen : time when minimum u occurs at each node
!     maxwxtimen : time when maximum u occurs at each node
!     minwytimen : time when minimum v occurs at each node
!     maxwytimen : time when maximum v occurs at each node
!
!     Attributes:
!     model             : ADCIRC
!     version           : model version number
!     grid_type         : triangular
!     rundes            : description of simulation
!     agrid             : grid ID 
!     title             : description of the dataset
!     institution       : where the original data was produced
!     source            : method of production of the original data
!     history           : provides an audit trail for modifications to the original data
!     references        : published or web-based references that describe the data or methods used to produce it
!     comments          : miscellaneous information about the data
!     contact           : user's institution/ e-mail address
!     creation_date     : date the data was created
!     modification_date : date the data was modified
!     Parameters from fort.15 input file 
!
!     Revision history:
!
!     Date            Programmer                  Description of change
!     ----            ----------                  ----------------------
!     03/30/07  Cristina Forbes, PSGS @ UNC-IMS   Wrote original code
!     03/30/08  Cristina Forbes, PSGS @ UNC-IMS   Modified code for 
!                                                 globalio & hotstart 
!                                                 from binary files
!=======================================================================
      subroutine writenc72 (nout,fileunit,timenc_len,node_len,timenc,
     &                      nt, NDATA, u,v,
     &                      rundes, runid, agrid, irtype,
     &                      dtdp, ihot, ics, nolibf, nolifa,
     &                      nolica, nolicat, nwp, 
     &                      ncor, ntip, nws, nramp, tau0,
     &                      statim, reftim, rnday, dramp,
     &                      a00, b00, c00, h0,
     &                      slam0, sfea0, ffactor,esl, cori,
     &                      ntif, nbfr, elem_len, nface_len, nm,
     &                      x,y,depth,
     &                      nope_len, neta_len, nbdv, 
     &                      nbou_len, nvel_len, nbvv,
     &                      ncmode,
     &                      title, institution, source, history,
     &                      references, comments, host, convention,
     &                      contact,
     &                      base_date)


!--------------------------------------------
!     Force explicit declaration of all variables
!--------------------------------------------
      IMPLICIT NONE
!     -------------------
!     netCDF header file
!     -------------------
 !     include 'netcdf.inc'
!     -------------------
!     Error status return
!     -------------------
      INTEGER  iret
!     ---------
!     netCDF id
!     ---------
      INTEGER  ncid
!     -------------
!     Dimension ids
!     -------------
      INTEGER  node_dim
      INTEGER  timenc_dim
      INTEGER  elem_dim
      INTEGER  nface_dim
      INTEGER  nope_dim
      INTEGER  neta_dim
      INTEGER  nbou_dim
      INTEGER  nvel_dim
      INTEGER  slen_dim
      INTEGER  scode_dim
!     -----------------
!     Dimension lengths
!     -----------------
      INTEGER  node_len
      INTEGER  timenc_len
      INTEGER  elem_len
      INTEGER  nface_len
      INTEGER  nope_len
      INTEGER  neta_len
      INTEGER  nbou_len
      INTEGER  nvel_len
!     ------------
!     Variable ids
!     ------------
      INTEGER  timenc_id
      INTEGER  X_id
      INTEGER  Y_id
      INTEGER  DEPTH_id
      INTEGER  ELE_id
      INTEGER  mindepth_id
      INTEGER  maxdepth_id
      INTEGER  nbdv_id
      INTEGER  nbvv_id
      INTEGER  u_id
      INTEGER  v_id
      INTEGER  minut_id
      INTEGER  maxut_id
      INTEGER  minvt_id
      INTEGER  maxvt_id
      INTEGER  minun_id
      INTEGER  maxun_id
      INTEGER  minvn_id
      INTEGER  maxvn_id
      
      INTEGER  minunodet_id
      INTEGER  maxunodet_id
      INTEGER  minutimen_id
      INTEGER  maxutimen_id
      INTEGER  minvnodet_id
      INTEGER  maxvnodet_id
      INTEGER  minvtimen_id
      INTEGER  maxvtimen_id
      INTEGER  station_id
      INTEGER  scode_id
!     ---------------------------------------------
!     Rank (number of dimensions) for each variable
!     ---------------------------------------------
      INTEGER  timenc_rank
      INTEGER  station_rank
      INTEGER  scode_rank
      INTEGER  depth_rank
      INTEGER  ele_rank
      INTEGER  mindepth_rank
      INTEGER  maxdepth_rank
      INTEGER  x_rank
      INTEGER  y_rank
      INTEGER  nvdll_rank
      INTEGER  nbdv_rank
      INTEGER  nbvv_rank
      INTEGER  u_rank
      INTEGER  v_rank
      INTEGER  minut_rank
      INTEGER  maxut_rank
      INTEGER  minvt_rank
      INTEGER  maxvt_rank
      INTEGER  minun_rank
      INTEGER  maxun_rank
      INTEGER  minvn_rank
      INTEGER  maxvn_rank

      INTEGER  minunodet_rank
      INTEGER  maxunodet_rank
      INTEGER  minutimen_rank
      INTEGER  maxutimen_rank
      
      INTEGER  minvnodet_rank
      INTEGER  maxvnodet_rank
      INTEGER  minvtimen_rank
      INTEGER  maxvtimen_rank

      PARAMETER (timenc_rank = 1)
      PARAMETER (scode_rank = 2)
      PARAMETER (station_rank = 2)
      PARAMETER (depth_rank = 1)
      PARAMETER (ele_rank = 2)
      PARAMETER (mindepth_rank = 0)
      PARAMETER (maxdepth_rank = 0)
      PARAMETER (x_rank = 1)
      PARAMETER (y_rank = 1)
      PARAMETER (nvdll_rank = 1)
      PARAMETER (nbdv_rank = 2)
      PARAMETER (nbvv_rank = 2)
      PARAMETER (u_rank = 2)
      PARAMETER (v_rank = 2)
      PARAMETER (minut_rank = 1)
      PARAMETER (maxut_rank = 1)
      PARAMETER (minvt_rank = 1)
      PARAMETER (maxvt_rank = 1)
      PARAMETER (minun_rank = 1)
      PARAMETER (maxun_rank = 1)
      PARAMETER (minvn_rank = 1)
      PARAMETER (maxvn_rank = 1)

      PARAMETER (minunodet_rank = 1)
      PARAMETER (maxunodet_rank = 1)
      PARAMETER (minutimen_rank = 1)
      PARAMETER (maxutimen_rank = 1)

      PARAMETER (minvnodet_rank = 1)
      PARAMETER (maxvnodet_rank = 1)
      PARAMETER (minvtimen_rank = 1)
      PARAMETER (maxvtimen_rank = 1)

!     ---------------
!     Variable shapes
!     ---------------
      INTEGER  timenc_dims(timenc_rank)
      INTEGER  station_dims  (station_rank)
      INTEGER  x_dims(x_rank)
      INTEGER  y_dims(y_rank)
      INTEGER  depth_dims(depth_rank)
      INTEGER  ele_dims(ele_rank)
      INTEGER  nvdll_dims(nvdll_rank)
      INTEGER  nbdv_dims(nbdv_rank)
      INTEGER  nbvv_dims(nbvv_rank)
      INTEGER  u_dims(u_rank)
      INTEGER  v_dims(v_rank)
      INTEGER  minut_dims(minut_rank)
      INTEGER  maxut_dims(maxut_rank)
      INTEGER  minvt_dims(minvt_rank)
      INTEGER  maxvt_dims(maxvt_rank)
      INTEGER  minun_dims(minun_rank)
      INTEGER  maxun_dims(maxun_rank)
      INTEGER  minvn_dims(minvn_rank)
      INTEGER  maxvn_dims(maxvn_rank)

      INTEGER  minunodet_dims(minunodet_rank)
      INTEGER  maxunodet_dims(maxunodet_rank)
      INTEGER  minutimen_dims(minutimen_rank)
      INTEGER  maxutimen_dims(maxutimen_rank)

      INTEGER  minvnodet_dims(minvnodet_rank)
      INTEGER  maxvnodet_dims(maxvnodet_rank)
      INTEGER  minvtimen_dims(minvtimen_rank)
      INTEGER  maxvtimen_dims(maxvtimen_rank)

!     --------------
!     Data variables
!     --------------
      CHARACTER*5 FILEFORT
      CHARACTER*10 FILENAME
      CHARACTER*3 FILEEXT
      INTEGER   t,ntif,nbfr,nt
      INTEGER   irtype, ihot, ics, nolibf, nolifa, nolica, nolicat, nwp
      INTEGER   ntip, nws, nramp,ncor
      INTEGER   nvel, neta, nout
      INTEGER   element (nface_len , elem_len)
      INTEGER   nbdv (nope_len, neta_len)
      INTEGER   nbvv (nbou_len, nvel_len)
      INTEGER   nm (elem_len, nface_len)
      INTEGER   fileunit
      REAL*8    timenc(timenc_len)
      REAL(sz)  x(node_len)
      REAL(sz)  y(node_len)
      REAL(sz)  depth(node_len)
      REAL(sz)  u(node_len, timenc_len)
      REAL(sz)  v(node_len, timenc_len)
      REAL(sz), ALLOCATABLE, SAVE ::   minut(:)
      REAL(sz), ALLOCATABLE, SAVE ::   maxut(:)
      REAL(sz), ALLOCATABLE, SAVE ::   minvt(:)
      REAL(sz), ALLOCATABLE, SAVE ::   maxvt(:)
      REAL(sz), ALLOCATABLE, SAVE ::   minun72(:)
      REAL(sz), ALLOCATABLE, SAVE ::   maxun72(:)
      REAL(sz), ALLOCATABLE, SAVE ::   minvn72(:)
      REAL(sz), ALLOCATABLE, SAVE ::   maxvn72(:)
      
      INTEGER , ALLOCATABLE, SAVE ::  minunodet(:)
      INTEGER , ALLOCATABLE, SAVE ::  maxunodet(:)
      REAL(sz), ALLOCATABLE, SAVE ::  minutimen72(:)
      REAL(sz), ALLOCATABLE, SAVE ::  maxutimen72(:)
      
      INTEGER , ALLOCATABLE, SAVE ::  minvnodet(:)
      INTEGER , ALLOCATABLE, SAVE ::  maxvnodet(:)
      REAL(sz), ALLOCATABLE, SAVE ::  minvtimen72(:)
      REAL(sz), ALLOCATABLE, SAVE ::  maxvtimen72(:)

      REAL  minull
      REAL  maxull
      REAL  minvll
      REAL  maxvll
      real*8  mindepth
      real*8  maxdepth
      REAL*8 dtdp, statim,reftim, rnday, dramp, tau0
      REAL*8 a00, b00, c00,h0
      REAL(sz) slam0, sfea0, ffactor,esl, cori
      REAL(sz) slam0deg, sfea0deg
      CHARACTER(80) rundes, runid, agrid
      INTEGER, SAVE  :: FIRST72 = 1 
      INTEGER NDATA
      INTEGER, SAVE  :: NREC = 0
      INTEGER ncmode

!     -----------------
!     Attribute vectors
!     -----------------
      INTEGER intval(1),i,j
      INTEGER counti(1), starti(1), n
      INTEGER count(2), start(2)
      REAL    realval(1)
!     -----------------
!     Attribute general
!     -----------------
      CHARACTER(*) :: title, institution, source, history
      CHARACTER(*) :: references, comments, host, convention, contact
      CHARACTER(*) :: base_date

!     -----------------
!      date_string variables for time attribute
!     -----------------
      character date_string*40
      character now_date*8
      character big_ben*10
      character zone*5
      integer   values(8)
!     ---------------------------------
!      inquiring and reading variables
!     ---------------------------------

      INTEGER NDIMS, NVARS, NGATTS, UNLIMDIMID, RECID,NRECS
      CHARACTER*4 recname

!     -------------------------
!     Save dims, ids, variables
!     -------------------------
      save node_dim,  elem_dim, nface_dim
      save nope_dim, neta_dim, nbou_dim, nvel_dim
      save u_id, v_id, X_id, Y_id
      save minut_id, maxut_id, minun_id, maxun_id
      save minvt_id, maxvt_id, minvn_id, maxvn_id

      save minunodet_id, maxunodet_id, minutimen_id, maxutimen_id
      save minvnodet_id, maxvnodet_id, minvtimen_id, maxvtimen_id

      save DEPTH_id,ELE_id, mindepth_id, maxdepth_id, nbdv_id, nbvv_id
      save timenc_dim, timenc_id
      save station_id, scode_id
      save mindepth, maxdepth

!     ---------------------------
!     Opening fort.64 netCDF file
!     ---------------------------
      WRITE(6,*) "-----------------------------------------------------"

         IF(FIRST72.EQ.1.and.ncmode.NE.0) THEN
         
         ALLOCATE ( minun72 (node_len) )
         ALLOCATE ( maxun72 (node_len) )
         ALLOCATE ( minvn72 (node_len) )
         ALLOCATE ( maxvn72 (node_len) )
         
         ALLOCATE ( minutimen72 (node_len) )
         ALLOCATE ( maxutimen72 (node_len) )
         ALLOCATE ( minvtimen72 (node_len) )
         ALLOCATE ( maxvtimen72 (node_len) )
         
         do n=1,node_len
            minun72 (n) =  9999.d0
            maxun72 (n) = -9999.d0
            minvn72 (n) =  9999.d0
            maxvn72 (n) = -9999.d0

            minutimen72 (n) =  9999.d0
            maxutimen72 (n) = -9999.d0
            minvtimen72 (n) =  9999.d0
            maxvtimen72 (n) = -9999.d0

         end do
      ENDIF

!     --------------------------------
!     Switch order in array for NETCDF
!     --------------------------------
      do i=1, elem_len
         do j=1, nface_len
            element(j,i) = nm (i,j)
         end do
      end do
      
!     -----------------
!     Enter define mode
!     -----------------
         FILEFORT = 'fort.'
         FILEEXT = '.nc'
         WRITE(FILENAME,'(A5,I2,A3)')FILEFORT,FILEUNIT,FILEEXT
         WRITE(6,*) "Opening ",FILENAME,  " netCDF file ... "
         
       IF((FIRST72.EQ.1.AND.IHOT.EQ.0.AND.NCMODE.NE.2).OR.
     &    (FIRST72.EQ.1.AND.IHOT.NE.0.AND.NCMODE.NE.2.AND.
     &     nout.EQ.-3))THEN

           iret = nf_create(FILENAME, NF_CLOBBER, ncid)
           CALL check_err(iret)
           FIRST72=0
!     --------------------------------
!     Convert to degress for NETCDF I/O
!     --------------------------------
           SLAM0DEG=SLAM0
           SFEA0DEG=SFEA0
           IF (NCMODE.EQ.1) THEN
               SLAM0DEG=SLAM0*RAD2DEG
               SFEA0DEG=SFEA0*RAD2DEG
           ENDIF
!     -----------------
!     Define dimensions
!     -----------------
      iret = nf_def_dim(ncid, 'station', node_len, node_dim)
      CALL check_err(iret)
      iret = nf_def_dim(ncid, 'namelen', SNAMLEN, slen_dim)
      CALL check_err(iret)
      iret = nf_def_dim(ncid, 'time', nf_unlimited, timenc_dim)
      CALL check_err(iret)

!     ================
!     Define variables
!     ================
!     ----------------
!     Define time
!     ----------------
      timenc_dims(1) = timenc_dim
      iret = nf_def_var(ncid,'time',NF_DOUBLE, timenc_rank, timenc_dims, 
     &                  timenc_id)
      CALL check_err(iret)
      
!     -------------------
!     Define x
!     -------------------
      x_dims(1) = node_dim
      iret = nf_def_var(ncid, 'x', NF_DOUBLE, x_rank, x_dims, x_id)
      CALL check_err(iret)
!     -------------------
!     Define Y coordinate
!     --------------------
      Y_dims(1) = node_dim
      iret = nf_def_var(ncid, 'y', NF_DOUBLE, Y_rank, Y_dims, Y_id)
      CALL check_err(iret)
!     ---------------------
!     Define stations name
!     --------------------
      station_dims(1) = slen_dim
      station_dims(2) = node_dim
      iret = nf_def_var(ncid, 'station_name', NF_CHAR, station_rank, 
     1      station_dims, station_id)
      CALL check_err(iret)

!     ----------------
!     Define u vel
!     ----------------
      u_dims(1) = node_dim
      u_dims(2) = timenc_dim
      iret = nf_def_var(ncid, 'windx', NF_DOUBLE, u_rank, u_dims, 
     &                  u_id)
      CALL check_err(iret)

!     ----------------
!     Define v vel
!     ----------------
      v_dims(1) = node_dim
      v_dims(2) = timenc_dim
      iret = nf_def_var(ncid, 'windy', NF_DOUBLE, v_rank, v_dims, 
     &                  v_id)
      CALL check_err(iret)

!     -------------------------
!     Define minimum u at times
!     -------------------------
      minut_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'minuvelt', NF_DOUBLE, minut_rank, 
     &                  minut_dims, minut_id)
      CALL check_err(iret)
!     -------------------------
!     Define maximum u at times
!     -------------------------
      maxut_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'maxuvelt', NF_DOUBLE, maxut_rank, 
     &                  maxut_dims, maxut_id)
      CALL check_err(iret)

!     -------------------------
!     Define minimum v at times
!     -------------------------
      minvt_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'minvvelt', NF_DOUBLE, minvt_rank, 
     &                  minvt_dims, minvt_id)
      CALL check_err(iret)
      
!     -------------------------
!     Define maximum v at times
!     -------------------------
      maxvt_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'maxvvelt', NF_DOUBLE, maxvt_rank, 
     &                  maxvt_dims, maxvt_id)
      CALL check_err(iret)

!     -------------------------
!     Define minimum u at nodes
!     -------------------------
      minun_dims(1) = node_dim
      iret = nf_def_var(ncid, 'minuveln', NF_DOUBLE, minun_rank, 
     &                  minun_dims, minun_id)
      CALL check_err(iret)
!     -------------------------
!     Define maximum u at nodes
!     -------------------------
      maxun_dims(1) = node_dim
      iret = nf_def_var(ncid, 'maxuveln', NF_DOUBLE, maxun_rank, 
     &                  maxun_dims, maxun_id)
      CALL check_err(iret)

!     -------------------------
!     Define minimum v at nodes
!     -------------------------
      minvn_dims(1) = node_dim
      iret = nf_def_var(ncid, 'minvveln', NF_DOUBLE, minvn_rank, 
     &                  minvn_dims, minvn_id)
      CALL check_err(iret)
!     -------------------------
!     Define maximum v at nodes
!     -------------------------
      maxvn_dims(1) = node_dim
      iret = nf_def_var(ncid, 'maxvveln', NF_DOUBLE, maxvn_rank, 
     &                  maxvn_dims, maxvn_id)
      CALL check_err(iret)            
!     ----------------------------
!     Define minimum u nodes at times
!     ----------------------------
      minunodet_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'minunodet', NF_INT, minunodet_rank, 
     &                  minunodet_dims, minunodet_id)
      CALL check_err(iret)
!     ----------------------------
!     Define maximum u nodes at times
!     ----------------------------
      maxunodet_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'maxunodet', NF_INT, maxunodet_rank, 
     &                  maxunodet_dims, maxunodet_id)
      CALL check_err(iret)
!     ----------------------------
!     Define minimum u nodes at nodes
!     ----------------------------
      minutimen_dims(1) = node_dim
      iret = nf_def_var(ncid, 'minutimen', NF_DOUBLE, minutimen_rank, 
     &                  minutimen_dims, minutimen_id)
      CALL check_err(iret)
!     ----------------------------
!     Define maximum u nodes at nodes
!     ----------------------------
      maxutimen_dims(1) = node_dim
      iret = nf_def_var(ncid, 'maxutimen', NF_DOUBLE, maxutimen_rank, 
     &                  maxutimen_dims, maxutimen_id)
      CALL check_err(iret)

            
!     ----------------------------
!     Define minimum v nodes at times
!     ----------------------------
      minvnodet_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'minvnodet', NF_INT, minvnodet_rank, 
     &                  minvnodet_dims, minvnodet_id)
      CALL check_err(iret)
!     ----------------------------
!     Define maximum v nodes at times
!     ----------------------------
      maxvnodet_dims(1) = timenc_dim
      iret = nf_def_var(ncid, 'maxvnodet', NF_INT, maxvnodet_rank, 
     &                  maxvnodet_dims, maxvnodet_id)
      CALL check_err(iret)
!     ----------------------------
!     Define minimum v nodes at nodes
!     ----------------------------
      minvtimen_dims(1) = node_dim
      iret = nf_def_var(ncid, 'minvtimen', NF_DOUBLE, minvtimen_rank, 
     &                  minvtimen_dims, minvtimen_id)
      CALL check_err(iret)
!     ----------------------------
!     Define maximum v nodes at nodes
!     ----------------------------
      maxvtimen_dims(1) = node_dim
      iret = nf_def_var(ncid, 'maxvtimen', NF_DOUBLE, maxvtimen_rank, 
     &                  maxvtimen_dims, maxvtimen_id)
      CALL check_err(iret)
      
!     -----------------
!     Assign attributes
!     -----------------

      if(ics.eq.2) then
      
      iret = nf_put_att_text(ncid, X_id, 'long_name', 9, 'longitude')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, X_id, 'standard_name', 9, 
     &                      'longitude')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, X_id, 'units', 12,  VAR_DEGE)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, X_id, 'positive', 4, VAR_EAST)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, X_id, 'missing_values', NF_DOUBLE,
     &                        1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, X_id, '_FillValue', NF_DOUBLE, 1, 
     &                       doubleval)
      CALL check_err(iret)      
      
      iret = nf_put_att_text(ncid, Y_id, 'long_name', 8, 'latitude')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, Y_id, 'standard_name', 8, 
     &                       'latitude')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, Y_id, 'units', 13,  VAR_DEGN)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, Y_id, 'positive', 5, VAR_NORTH)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, Y_id, 'missing_values', NF_DOUBLE, 
     &       1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, Y_id, '_FillValue', NF_DOUBLE, 1, 
     &       doubleval)
      CALL check_err(iret)

       else

      
      iret = nf_put_att_text(ncid, X_id, 'long_name', 22, 
     &       'Cartesian coordinate x')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, X_id, 'standard_name', 12, 
     &                      'x_coordinate')
      CALL check_err(iret)
!    ----------------------
!     determine variable units
!    ----------------------
      if (G.LT.11.D0) then
          iret = nf_put_att_text(ncid, X_id, 'units', 6, VAR_METS)
      else
          iret = nf_put_att_text(ncid, X_id, 'units', 6, VAR_FEET)
      endif
      CALL check_err(iret)
      
      iret = nf_put_att_text(ncid, X_id, 'positive', 5, VAR_RG)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, X_id, 'missing_values', NF_DOUBLE, 
     &                       1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, X_id, '_FillValue', NF_DOUBLE, 1, 
     &                       doubleval)
      CALL check_err(iret)
!     ------------------------------
!     Define y-coordinate attributes
!     ------------------------------     
      
      iret = nf_put_att_text(ncid, Y_id, 'long_name', 22, 
     &       'Cartesian coordinate y')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, Y_id, 'standard_name', 12, 
     &                       'y_coordinate')
      CALL check_err(iret)
!    ----------------------
!     determine variable units
!    ----------------------
      if (G.LT.11.D0) then
         iret = nf_put_att_text(ncid, Y_id, 'units', 6,  VAR_METS)
      else
         iret = nf_put_att_text(ncid, Y_id, 'units', 6,  VAR_FEET)
      endif
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, Y_id, 'positive', 2, VAR_UP)
      CALL check_err(iret)

      iret = nf_put_att_double(ncid, Y_id, 'missing_values', 
     &       NF_DOUBLE, 1, doubleval)
      CALL check_err(iret)

      iret = nf_put_att_double(ncid, Y_id, '_FillValue', NF_DOUBLE, 1, 
     &       doubleval)
      CALL check_err(iret)
      
      endif
!!+++++++++++++++++++++++++++++++++++++++++++++++
    !    ----------------------------------------------------------
!     determine variable names according to coordinates and nws
!    ----------------------------------------------------------
      if (ics.eq.2) then
        if (abs(nws).eq.1.or.abs(nws).eq.2) then
             iret = nf_put_att_text(ncid, minut_id, 'long_name', 51, 
     &       'minimum station e/w wind stress at each output time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minut_id, 'standard_name', 31,
     &       'minimum_station_e/w_wind_stress')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvt_id, 'long_name', 51, 
     &       'minimum station n/s wind stress at each output time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvt_id, 'standard_name', 31,
     &       'minimum_station_n/s_wind_stress')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxut_id, 'long_name', 51, 
     &       'maximum station e/w wind stress at each output time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxut_id, 'standard_name', 31,
     &       'maximum_station_e/w_wind_stress')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvt_id, 'long_name', 51, 
     &       'maximum station n/s wind stress at each output time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvt_id, 'standard_name', 31,
     &       'maximum_station_n/s_wind_stress')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minun_id, 'long_name', 39, 
     &       'minimum e/w wind stress at each station')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minun_id, 'standard_name', 39,
     &       'minimum_e/w_wind_stress_at_each_station')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxun_id, 'long_name', 39, 
     &        'maximum e/w wind stress at each station')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxun_id, 'standard_name', 39,
     &       'maximum_e/w_wind_stress_at_each_station')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvn_id, 'long_name', 39, 
     &       'minimum station n/s wind stress at each station')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvn_id, 'standard_name', 39,
     &       'minimum_n/s_wind_stress_at_each_station')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvn_id, 'long_name', 39, 
     &       'maximum n/s wind stress at each station')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvn_id, 'standard_name', 39,
     &       'maximum_n/s_wind_stress_at_each_station')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, minunodet_id, 'long_name',54, 
     &        'station of minimum e/w wind stress at each output time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minunodet_id, 'standard_name',
     &        34,'station_of_minimum_e/w_wind_stress')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, maxunodet_id, 'long_name',54, 
     &        'station of maximum e/w wind stress at each output time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxunodet_id, 'standard_name',
     &        34,'station_of_maximum_e/w_wind_stress')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, minvnodet_id, 'long_name',54,
     &        'station of minimum n/s wind stress at each output time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvnodet_id, 'standard_name',
     &        34,'station_of_minimum_n/s_wind_stress')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, maxvnodet_id, 'long_name',51, 
     &         'station of maximum n/s wind stress at each output time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvnodet_id, 'standard_name',
     &        34,'station_of_maximum_n/s_wind_stress')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, minutimen_id, 'long_name', 
     &        47, 'time of minimum e/w wind stress at each station')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minutimen_id, 'standard_name',
     &        31,'time_of_minimum_e/w_wind_stress')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, maxutimen_id, 'long_name', 
     &        47, 'time of maximum e/w wind stress at each station')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxutimen_id, 'standard_name',
     &        31,'time_of_maximum_e/w_wind_stress')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, minvtimen_id, 'long_name', 
     &        47, 'time of minimum n/s wind stress at each station')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvtimen_id, 'standard_name',
     &        31,'time_of_minimum_n/s_wind_stress')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, maxvtimen_id, 'long_name', 
     &        47, 'time of maximum n/s wind stress at each station')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvtimen_id, 'standard_name',
     &        31,'time_of_maximum_n/s_wind_stress')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, u_id, 'long_name', 23, 
     &                       'station e/w wind stress')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, u_id, 'standard_name', 23, 
     &                      'station_e/w_wind_stress')
             CALL check_err(iret)
             iret = nf_put_att_text(ncid, u_id, 'positive', 4, VAR_EAST)

              iret = nf_put_att_text(ncid, v_id, 'long_name', 23, 
     &                       'station n/s wind stress')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, v_id, 'standard_name', 23, 
     &                      'station_n/s_wind_stress')
             CALL check_err(iret)
            iret = nf_put_att_text(ncid, v_id, 'positive', 5, VAR_NORTH)



        else
              iret = nf_put_att_text(ncid, minut_id, 'long_name', 53, 
     &       'minimum station e/w wind velocity at each output time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minut_id, 'standard_name', 53,
     &       'minimum_station_e/w_wind_velocity_at_each_output_time')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, minvt_id, 'long_name', 53, 
     &       'minimum station n/s wind velocity at each output time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvt_id, 'standard_name', 53,
     &       'minimum_station_n/s_wind_velocity_at_each_output_time')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, maxut_id, 'long_name', 53, 
     &       'maximum station e/w wind velocity at each output time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxut_id, 'standard_name', 53,
     &       'maximum_station_e/w_wind_velocity_at_each_output_time')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, maxvt_id, 'long_name', 53, 
     &       'maximum station n/s wind velocity at each output time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvt_id, 'standard_name', 53,
     &       'maximum_station_n/s_wind_velocity_at_each_output_time')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, minun_id, 'long_name', 41, 
     &       'minimum e/w wind velocity at each station')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minun_id, 'standard_name', 41,
     &       'minimum_e/w_wind_velocity_at_each_station')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxun_id, 'long_name', 41, 
     &       'maximum e/w wind velocity at each station')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxun_id, 'standard_name', 41,
     &       'maximum_e/w_wind_velocity_at_each_station')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvn_id, 'long_name', 41, 
     &       'minimum n/s wind velocity at each station')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvn_id, 'standard_name', 41,
     &       'minimum_n/s_wind_velocity_at_each_station')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvn_id, 'long_name', 41, 
     &       'maximum n/s wind velocity at each station')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvn_id,'standard_name', 41,
     &       'maximum_n/s_wind_velocity_at_each_station')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, minunodet_id, 'long_name',56, 
     &       'station of minimum e/w wind velocity at each output time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minunodet_id, 'standard_name',
     &        56,
     &       'station_of_minimum_e/w_wind_velocity_at_each_output_time')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, maxunodet_id, 'long_name',56, 
     &        'station of maximum e/w wind velocity at each output time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxunodet_id, 'standard_name',
     &        56,
     &      'station_of_maximum_e/w_wind_velocity_at_each_output_time')
             CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvnodet_id, 'long_name',56, 
     &       'station of minimum n/s wind velocity at each output time')
             CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvnodet_id, 'standard_name',
     &       56,
     &       'station_of_minimum_n/s_wind_velocity_at_each_output_time')
             CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvnodet_id, 'long_name',56, 
     &      'station of maximum n/s wind velocity at each output time')
             CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvnodet_id, 'standard_name',
     &       56,
     &       'station_of_maximum_n/s_wind_velocity_at_each_output_time')
             CALL check_err(iret)

             iret = nf_put_att_text(ncid, minutimen_id, 'long_name', 
     &       49,'time of minimum e/w wind velocity at each station')
             CALL check_err(iret)
             iret = nf_put_att_text(ncid, minutimen_id, 'standard_name',
     &       49,'time_of_minimum_e/w_wind_velocity_at_each_station')
             CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxutimen_id, 'long_name', 
     &       49,'time of maximum e/w wind velocity at each station')
             CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxutimen_id, 'standard_name',
     &       49,'time_of_maximum_e/w_wind_velocity_at_each_station')
             CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvtimen_id, 'long_name', 
     &       49,'time of minimum n/s wind velocity at each station')
             CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvtimen_id, 'standard_name',
     &       49,'time_of_minimum_n/s_wind_velocity_at_each_station')
             CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvtimen_id, 'long_name', 
     &       49,'time of maximum n/s wind velocity at each station')
             CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvtimen_id, 'standard_name',
     &       49,'time_of_maximum_n/s_wind_velocity_at_each_station')
             CALL check_err(iret)

              iret = nf_put_att_text(ncid, u_id, 'long_name', 25, 
     &                       'station e/w wind velocity')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, u_id, 'standard_name', 25, 
     &                      'station_e/w_wind_velocity')
             CALL check_err(iret)
             iret = nf_put_att_text(ncid, u_id, 'positive', 4, VAR_EAST)

              iret = nf_put_att_text(ncid, v_id, 'long_name', 25, 
     &                       'station n/s wind velocity')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, v_id, 'standard_name', 25, 
     &                      'station_n/s_wind_velocity')
             CALL check_err(iret)
             iret = nf_put_att_text(ncid, v_id, 'positive', 5,VAR_NORTH)


          endif
       else
         if (abs(nws).eq.1.or.abs(nws).eq.2) then
             iret = nf_put_att_text(ncid, minut_id, 'long_name', 62, 
     & 'minimum station wind stress in x-direction at each output time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minut_id, 'standard_name', 62,
     & 'minimum_station_wind_stress_in_x_direction_at_each_output_time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvt_id, 'long_name', 62, 
     & 'minimum station wind stress in y-direction at each output time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvt_id, 'standard_name', 62,
     & 'minimum_station_wind_stress_in_y_direction_at_each_output_time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxut_id, 'long_name', 62, 
     & 'maximum station wind stress in x-direction at each output time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxut_id, 'standard_name', 62,
     & 'maximum_station_wind_stress_in_x_direction_at_each_output_time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvt_id, 'long_name', 62, 
     & 'maximum station wind stress in y-direction at each output time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvt_id, 'standard_name', 62,
     & 'maximum_station_wind_stress_in_y_direction_at_each_output_time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minun_id, 'long_name', 50, 
     &       'minimum wind stress in x-direction at each station')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minun_id, 'standard_name', 50,
     &       'minimum_wind_stress_in_x_direction_at_each_station')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxun_id, 'long_name', 50, 
     &       'maximum wind stress in x-direction at each station')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxun_id, 'standard_name', 50,
     &       'maximum_wind_stress_in_x_direction_at_each_station')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvn_id, 'long_name', 50, 
     &       'minimum wind stress in y-direction at each station')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvn_id, 'standard_name', 50,
     &       'minimum_wind_stress_in_y_direction_at_each_station')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvn_id, 'long_name', 50, 
     &       'maximum wind stress in y-direction at each station')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvn_id, 'standard_name', 50,
     &       'maximum_wind_stress_in_y_directionat_each_station')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minunodet_id, 'long_name', 65, 
     &           'station of minimum wind stress in x-direction at each  
     &output time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minunodet_id, 'standard_name', 
     &       65, 'station_of_minimum_wind_stress_in_x_direction_at_each_
     &output_time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxunodet_id, 'long_name', 65, 
     &           'station of maximum wind stress in x-direction at each 
     &output time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxunodet_id, 'standard_name', 
     &       65, 'station_of_maximum_wind_stress_in_x_direction_at_each_
     &output_time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvnodet_id, 'long_name', 62, 
     &           'station of minimum wind stress in y-direction at each  
     &output time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvnodet_id, 'standard_name', 
     &       62, 'station_of_minimum_wind_stress_in_y_direction_at_each_
     &output_time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvnodet_id, 'long_name', 62, 
     &           'station of maximum wind stress in y-direction at each 
     &output time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvnodet_id, 'standard_name', 
     &       62, 'station_of_maximum_wind_stress_in_y_direction_at_each_
     &output_time')
              CALL check_err(iret)

            iret = nf_put_att_text(ncid, minutimen_id, 'long_name', 58, 
     &     'time of minimum wind stress in x-direction at each station')
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, minutimen_id, 'standard_name', 
     &       58,
     &     'time_of_minimum_wind_stress_in_x_direction_at_each_station')
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, maxutimen_id, 'long_name', 58, 
     &     'time of maximum wind stress in x-direction at each station')
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, maxutimen_id, 'standard_name', 
     &       58,
     &     'time_of_maximum_wind_stress_in_x_direction_at_each_station')
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, minvtimen_id, 'long_name', 58, 
     &     'time of minimum wind stress in y-direction at each station')
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, minvtimen_id, 'standard_name', 
     &      58,
     &     'time_of_minimum_wind_stress_in_y_direction_at_each_station')
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, maxvtimen_id, 'long_name', 58, 
     &     'time of maximum wind stress in y-direction at each station')
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, maxvtimen_id, 'standard_name', 
     &      58,
     &     'time_of_maximum_wind_stress_in_y_direction_at_each_station')
            CALL check_err(iret)

            iret = nf_put_att_text(ncid, u_id, 'long_name', 33, 
     &                       'station wind stress in x-direction')
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, u_id, 'standard_name', 33, 
     &                      'station_wind_stress_in_x_direction')
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, u_id, 'positive', 5, VAR_RG)

            iret = nf_put_att_text(ncid, v_id, 'long_name', 33, 
     &                       'station wind stress in y-direction')
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, v_id, 'standard_name', 33, 
     &                      'station_wind_stress_in_y_direction')
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, v_id, 'positive', 2, VAR_UP)
             
         else
             iret = nf_put_att_text(ncid, minut_id, 'long_name', 64, 
     &            'minimum station wind velocity in x-direction at each 
     &output time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minut_id, 'standard_name', 64,
     &            'minimum_station_wind_velocity_in_x_direction_at_each_
     &output_time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvt_id, 'long_name', 64, 
     &            'minimum station wind velocity in y-direction at each 
     &output time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvt_id, 'standard_name', 64,
     &            'minimum_station_wind_velocity_in_y_direction_at_each_
     &output_time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxut_id, 'long_name', 64, 
     &            'maximum station wind velocity in x-direction at each 
     &output time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxut_id, 'standard_name', 64,
     &            'maximum_station_wind_velocity_in_x_direction_at_each_
     &output_time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvt_id, 'long_name', 64, 
     &            'maximum station wind velocity in y-direction at each 
     &output time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvt_id, 'standard_name', 64,
     &            'maximum_station_wind_velocity_in_y_direction_at_each_
     &output_time')
              CALL check_err(iret)
              
             iret = nf_put_att_text(ncid, minun_id, 'long_name', 52, 
     &       'minimum wind velocity in x-direction at each station')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minun_id, 'standard_name', 52,
     &       'minimum_wind_velocity_in_x_direction_at_each_station')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxun_id, 'long_name', 52, 
     &       'maximum wind velocity in x-direction at each station')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxun_id, 'standard_name', 52,
     &       'maximum_wind_velocity_in_x_direction_at_each_station')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvn_id, 'long_name', 52, 
     &       'minimum wind velocity in y-direction at each station')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvn_id, 'standard_name', 52,
     &       'minimum_wind_velocity_in_y_direction_at_each_station')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvn_id, 'long_name', 52, 
     &       'maximum wind velocity in y-direction at each station')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvn_id, 'standard_name', 52,
     &       'maximum_wind_velocity_in_y_direction_at_each_station')
              CALL check_err(iret)
              
             iret = nf_put_att_text(ncid, minunodet_id, 'long_name', 67, 
     &         'station of minimum wind velocity in x-direction at each 
     &output time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minunodet_id, 'standard_name', 
     &      67,'station_of_minimum_wind_velocity_in_x_direction_at_each_
     &output_time')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, maxunodet_id, 'long_name', 67, 
     &         'station of maximum wind velocity in x-direction at each 
     &output time')
             CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxunodet_id, 'standard_name', 
     &      67,'station_of_maximum_wind_velocity_in_x_direction_at_each_
     &output_time')
              CALL check_err(iret)

             iret = nf_put_att_text(ncid, minvnodet_id, 'long_name', 67, 
     &         'station of minimum wind velocity in y-direction at each 
     &output time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvnodet_id, 'standard_name', 
     &      67,'station_of_minimum_wind_velocity_in_y_direction_at_each_
     &output_time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvnodet_id, 'long_name', 67, 
     &         'station of maximum wind velocity in y-direction at each 
     &output time')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvnodet_id, 'standard_name', 
     &      67,'station_of_maximum_wind_velocity_in_y_direction_at_each_
     &output_time')
              CALL check_err(iret)

             iret = nf_put_att_text(ncid, minutimen_id, 'long_name', 60, 
     &            'time of minimum wind velocity in x-direction at each 
     &station')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minutimen_id, 'standard_name', 
     &         60,'time_of_minimum_wind_velocity_in_x_direction_at_each_
     &station')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxutimen_id, 'long_name', 60, 
     &            'time of maximum wind velocity in x-direction at each 
     &station')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxutimen_id, 'standard_name', 
     &         60,'time_of_maximum_wind_velocity_in_x_direction_at_each_
     &station')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvtimen_id, 'long_name', 60, 
     &            'time of minimum wind velocity in y-direction at each 
     &station')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, minvtimen_id, 'standard_name', 
     &       60,  'time_of_minimum_wind_velocity_in_y_direction_at_each_
     &station')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvtimen_id, 'long_name', 60, 
     &            'time of maximum wind velocity in y-direction at each 
     &station')
              CALL check_err(iret)
             iret = nf_put_att_text(ncid, maxvtimen_id, 'standard_name', 
     &       60,  'time_of_maximum_wind_velocity_in_y_direction_at_each_
     &station')
              CALL check_err(iret)

              iret = nf_put_att_text(ncid, u_id, 'long_name', 36, 
     &                       'station wind velocity in x-direction')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, u_id, 'standard_name', 36, 
     &                      'station_wind_velocity_in_x_direction')
             CALL check_err(iret)
             iret = nf_put_att_text(ncid, u_id, 'positive', 5, VAR_RG)

              iret = nf_put_att_text(ncid, v_id, 'long_name', 36, 
     &                       'station_wind velocity in y-direction')
              CALL check_err(iret)
              iret = nf_put_att_text(ncid, v_id, 'standard_name', 36, 
     &                      'station_wind_velocity_in_y_direction')
             CALL check_err(iret)
             iret = nf_put_att_text(ncid, v_id, 'positive', 2, VAR_UP)
             CALL check_err(iret)

         endif
      endif

!    -------------------------
!     determine variable units
!    -------------------------
      if (G.LT.11.D0) then
         if (abs(nws).gt.2.and.abs(nws).lt.100) then
            iret = nf_put_att_text(ncid, minut_id,'units',13, VAR_MPS)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, minvt_id,'units',13, VAR_MPS)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, maxut_id,'units',13, VAR_MPS)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, maxvt_id,'units',13, VAR_MPS)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, minun_id,'units',13, VAR_MPS)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, maxun_id,'units',13, VAR_MPS)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, minvn_id,'units',13, VAR_MPS)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, maxvn_id,'units',13, VAR_MPS)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, u_id, 'units',13, VAR_MPS)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, v_id, 'units',13, VAR_MPS)
            CALL check_err(iret)
         else
            iret = nf_put_att_text(ncid, minut_id,'units',19, VAR_MPS_2)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, minvt_id,'units',19, VAR_MPS_2)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, maxut_id,'units',19, VAR_MPS_2)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, maxvt_id,'units',19, VAR_MPS_2)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, minun_id,'units',19, VAR_MPS_2)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, maxun_id,'units',19, VAR_MPS_2)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, minvn_id,'units',19, VAR_MPS_2)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, maxvn_id,'units',19, VAR_MPS_2)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, u_id, 'units',19, VAR_MPS_2)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, v_id, 'units',19, VAR_MPS_2)
            CALL check_err(iret)         
         endif
      else
         if (abs(nws).gt.2.and.abs(nws).lt.100) then
            iret = nf_put_att_text(ncid, minut_id,'units',11, VAR_FPS)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, minvt_id,'units',11, VAR_FPS)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, maxut_id,'units',11, VAR_FPS)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, maxvt_id,'units',11, VAR_FPS)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, minun_id,'units',11, VAR_FPS)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, maxun_id,'units',11, VAR_FPS)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, minvn_id,'units',11, VAR_FPS)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, maxvn_id,'units',11, VAR_FPS)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, u_id, 'units',11, VAR_FPS)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, v_id, 'units',11, VAR_FPS)
            CALL check_err(iret)
         else
            iret = nf_put_att_text(ncid, minut_id,'units',17, VAR_FPS_2)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, minvt_id,'units',17, VAR_FPS_2)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, maxut_id,'units',17, VAR_FPS_2)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, maxvt_id,'units',17, VAR_FPS_2)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, minun_id,'units',17, VAR_FPS_2)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, maxun_id,'units',17, VAR_FPS_2)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, minvn_id,'units',17, VAR_FPS_2)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, maxvn_id,'units',17, VAR_FPS_2)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, u_id, 'units',17, VAR_FPS_2)
            CALL check_err(iret)
            iret = nf_put_att_text(ncid, v_id, 'units',17, VAR_FPS_2)
            CALL check_err(iret)
         endif
      endif
      CALL check_err(iret)
      
!    -------------------------------------
!     other variables units attributes
!    ------------------------------------      

      iret = nf_put_att_text(ncid, minunodet_id,'units',14,  VAR_NOND)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, maxunodet_id,'units',14,  VAR_NOND)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, minvnodet_id,'units',14,  VAR_NOND)
      CALL check_err(iret)      
      iret = nf_put_att_text(ncid, maxvnodet_id,'units',14,  VAR_NOND)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, minutimen_id,'units',4, VAR_SECS)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, maxutimen_id,'units',4, VAR_SECS)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, minvtimen_id,'units',4, VAR_SECS)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, maxvtimen_id,'units',4, VAR_SECS)
      CALL check_err(iret)

!    -------------------------------------
!     u and v velocities other attributes
!    ------------------------------------      

      CALL check_err(iret)
      iret = nf_put_att_double(ncid, u_id, 'missing_values', 
     &       NF_DOUBLE, 1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, u_id, '_FillValue', NF_DOUBLE, 1, 
     &                       doubleval)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, u_id, 'dry_Value', NF_DOUBLE, 1, 
     &                       doubleval)
      CALL check_err(iret)
      iret = nf_put_att_DOUBLE(ncid, v_id, 'missing_values', 
     &       NF_DOUBLE, 1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_DOUBLE(ncid, v_id, '_FillValue', NF_DOUBLE, 1, 
     &                       doubleval)
      CALL check_err(iret)
      iret = nf_put_att_DOUBLE(ncid, v_id, 'dry_Value', NF_DOUBLE, 1, 
     &                       doubleval)
      CALL check_err(iret)
    
  
!!+++++++++++++++++++++++++++++++++++++++++++++++
!     -------------------------------------------------------------
!     Define station names and code attributes
!     -------------------------------------------------------------
      iret = nf_put_att_text(ncid,station_id,'long_name',12,
     & 'station name')
      CALL check_err(iret)


!    -----------------------
!     Define time attributes
!    -----------------------
      iret = nf_put_att_text(ncid,timenc_id,'long_name',10,'model time')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, timenc_id, 'standard_name', 4,'time')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, timenc_id, 'units', 23,VAR_TUN)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, timenc_id, 'base_date',
     . (LEN(TRIM(base_date(1:INDEX(base_date,"!")-1)))), base_date)
      CALL check_err(iret)
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
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'description', 
     &       LEN(TRIM(rundes(1:INDEX(rundes,"!")-1))), rundes)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'grid',
     &       LEN(TRIM(agrid)), agrid)       
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'title',
     &       LEN(TRIM(title(1:INDEX(title,"!")-1))), 
     &       title)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'institution', 
     &       LEN(TRIM(institution(1:INDEX(institution,"!")-1))), 
     &       institution)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'source', 
     &       LEN(TRIM(source(1:INDEX(source,"!")-1))),
     &       source)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'history', 
     &       LEN(TRIM(history(1:INDEX(history,"!")-1))),
     &       history)
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
     &       LEN(TRIM(host(1:INDEX(host,"!")-1))), 
     &       host)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'convention', 
     &       LEN(TRIM(convention(1:INDEX(convention,"!")-1))), 
     &       convention)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'contact', 
     &       LEN(TRIM(contact(1:INDEX(contact,"!")-1))), 
     &       contact)
      CALL check_err(iret)
 
 
      CALL date_and_Time(now_date,big_ben,zone,values)
      write(date_string,71) values(1),values(2),values(3)
     & ,values(5),values(6),values(7),   (values(4))/60
 71   format(I4,'-',I2.2,'-',i2.2,' ',i2,':',i2.2,':',i2.2,' '
     &  ,i3.2,':00')
      iret = nf_put_att_text(ncid, NF_GLOBAL,'creation_date'
     & ,len_trim(date_string),date_string)
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
      iret = nf_put_att_int(ncid, NCGLOBAL, 'irtype', NF_INT, 1, irtype)
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
      iret = nf_put_att_double(ncid, NCGLOBAL, 'ffactor', NF_DOUBLE, 1, 
     &                         ffactor)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, NCGLOBAL, 'esl', NF_DOUBLE, 1, 
     &                         esl)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, NCGLOBAL, 'cori', NF_DOUBLE, 1, 
     &                         cori)
      CALL check_err(iret)
      iret = nf_put_att_int(ncid, NCGLOBAL, 'ntif', NF_INT, 1, ntif)
      CALL check_err(iret)
      iret = nf_put_att_int(ncid, NCGLOBAL, 'nbfr', NF_INT, 1, nbfr)
      CALL check_err(iret)
     
!     ======================
!     Leave define mode
!     ======================
      iret = nf_enddef(ncid)
      CALL check_err(iret)
!     ======================
!     Store static variables
!     ======================
!     -------
!     Store X
!     -------
      iret = nf_put_var_double(ncid, X_id, X)
      CALL check_err(iret)
!     -------
!     Store Y
!     -------
      iret = nf_put_var_double(ncid, Y_id, Y)
      CALL check_err(iret)
!     --------------------
!     Store station name
!     --------------------
      do i=1,node_len
         start(1)=1
         start(2)=i
         count(1)= LEN(TRIM(STATNAMEM(i,SNAMLEN)))
         count(2)=1
         iret = nf_put_vara_text(ncid,station_id,start,count,
     &   STATNAMEM(i,SNAMLEN) )
         CALL check_err(iret)
      end do

!     =====================
!     Open old netcdf file
!     =====================
      ELSE
          iret = NF_OPEN(FILENAME, NF_WRITE, ncid)
          CALL check_err(iret)   
!     ------------------
!     Inquire variables
!     ------------------
          iret=nf_inq_unlimdim(ncid,recid)
          CALL check_err(iret)   
          iret=nf_inq_dimlen(ncid,recid,nrec)
          CALL check_err(iret)   
      ENDIF
!     ------------------
!     Update Global Attributes
!     ------------------
      IF(IHOT.NE.0.AND.FIRST72.EQ.1) THEN
         iret = NF_REDEF (ncid)
         CALL check_err(iret)
         iret = nf_put_att_int(ncid, NCGLOBAL, 'ihot', NF_INT, 1, ihot)
         CALL check_err(iret)
         iret = nf_put_att_int(ncid, NCGLOBAL, 'nramp', NF_INT, 1,nramp)
         CALL check_err(iret)
         iret = nf_put_att_double(ncid, NCGLOBAL, 'tau0', NF_DOUBLE, 1, 
     &                         tau0)
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
         iret = nf_put_att_double(ncid, NCGLOBAL, 'ffactor', NF_DOUBLE, 1, 
     &                         ffactor)
         CALL check_err(iret)
         iret = nf_put_att_double(ncid, NCGLOBAL, 'esl', NF_DOUBLE, 1, 
     &                         esl)
         CALL check_err(iret)
         CALL date_and_time(now_date,big_ben,zone,values)
         WRITE(date_string,71) values(1),values(2),values(3),
     &         values(5),values(6),values(7),   (values(4))/60
         iret = nf_put_att_text(ncid, NF_GLOBAL,'modification_date',
     &          len_trim(date_string),date_string)
         CALL check_err(iret)
         iret = NF_ENDDEF (ncid)
         CALL check_err(iret)
       ENDIF
      IF(ncmode.eq.0) go to 9999
      IF(FIRST72.EQ.1.and.
     & ((ncmode.ne.0.and.IHOT.NE.0).or.(ncmode.eq.2.and.IHOT.eq.0).or.
     &  (ncmode.eq.2.and.IHOT.ne.0.and.nout.eq.-3)))THEN
         iret=nf_inq_varid(ncid,"time",timenc_id)
         iret=nf_inq_varid(ncid,"windx",u_id)
         iret=nf_inq_varid(ncid,"windy",v_id)
         iret=nf_inq_varid(ncid,"minuvelt",minut_id)
         iret=nf_inq_varid(ncid,"maxuvelt",maxut_id)
         iret=nf_inq_varid(ncid,"minvvelt",minvt_id)
         iret=nf_inq_varid(ncid,"maxvvelt",maxvt_id)
         iret=nf_inq_varid(ncid,"minuveln",minun_id)
         iret=nf_inq_varid(ncid,"maxuveln",maxun_id)
         iret=nf_inq_varid(ncid,"minvveln",minvn_id)
         iret=nf_inq_varid(ncid,"maxvveln",maxvn_id)
         iret=nf_inq_varid(ncid,"minunodet",minunodet_id)
         iret=nf_inq_varid(ncid,"maxunodet",maxunodet_id)
         iret=nf_inq_varid(ncid,"minvnodet",minvnodet_id)
         iret=nf_inq_varid(ncid,"maxvnodet",maxvnodet_id)
         iret=nf_inq_varid(ncid,"minutimen",minutimen_id)
         iret=nf_inq_varid(ncid,"maxutimen",maxutimen_id)
         iret=nf_inq_varid(ncid,"minvtimen",minvtimen_id)
         iret=nf_inq_varid(ncid,"maxvtimen",maxvtimen_id)
      ENDIF
      IF(IHOT.NE.0.AND.FIRST72.EQ.1.and.ncmode.ne.0.and.nout.eq.3)THEN
         iret = nf_get_var_double(ncid, minun_id, minun72)
         CALL check_err(iret)
         iret = nf_get_var_double(ncid, maxun_id, maxun72)
         CALL check_err(iret)
         iret = nf_get_var_double(ncid, minvn_id, minvn72)
         CALL check_err(iret)
         iret = nf_get_var_double(ncid, maxvn_id, maxvn72)
         CALL check_err(iret)
         
         iret = nf_get_var_double(ncid, minutimen_id, minutimen72)
         CALL check_err(iret)
         iret = nf_get_var_double(ncid, maxutimen_id, maxutimen72)
         CALL check_err(iret)
         iret = nf_get_var_double(ncid, minvtimen_id, minvtimen72)
         CALL check_err(iret)
         iret = nf_get_var_double(ncid, maxvtimen_id, maxvtimen72)
         CALL check_err(iret)
      ENDIF
      
      FIRST72=0
      NREC = NREC + 1
      
!     ----------
!     Store time
!     ----------
      starti(1)=nrec
      counti(1)=timenc_len
      iret = nf_put_vara_DOUBLE(ncid, timenc_id, starti,counti, timenc)
      CALL check_err(iret)
!     ---------------------------
!     Store u
!     ---------------------------
      count(1)=node_len
      count(2)=timenc_len
      start(1)=1
      start(2)=nrec
      iret = nf_put_vara_DOUBLE(ncid, u_id, start, count, u)
      CALL check_err(iret)
!     ---------------------------
!     Store v
!     ---------------------------
      count(1)=node_len
      count(2)=timenc_len
      start(1)=1
      start(2)=nrec
      iret = nf_put_vara_DOUBLE(ncid, v_id, start, count, v)
      CALL check_err(iret)
      
!     ----------------------------------------
!     Allocating memory for min and max arrays
!     ----------------------------------------
      IF(.NOT. ALLOCATED(minut))      ALLOCATE( minut(timenc_len))
      IF(.NOT. ALLOCATED(maxut))      ALLOCATE( maxut(timenc_len))
      IF(.NOT. ALLOCATED(minun72))    ALLOCATE( minun72(node_len))
      IF(.NOT. ALLOCATED(maxun72))    ALLOCATE( maxun72(node_len))

      IF(.NOT. ALLOCATED(minvt))      ALLOCATE( minvt(timenc_len))
      IF(.NOT. ALLOCATED(maxvt))      ALLOCATE( maxvt(timenc_len))
      IF(.NOT. ALLOCATED(minvn72))    ALLOCATE( minvn72(node_len))
      IF(.NOT. ALLOCATED(maxvn72))    ALLOCATE( maxvn72(node_len))

      IF(.NOT. ALLOCATED(minunodet))  ALLOCATE( minunodet(timenc_len))
      IF(.NOT. ALLOCATED(maxunodet))  ALLOCATE( maxunodet(timenc_len))
      IF(.NOT. ALLOCATED(minutimen72))ALLOCATE( minutimen72(node_len))
      IF(.NOT. ALLOCATED(maxutimen72))ALLOCATE( maxutimen72(node_len))

      IF(.NOT. ALLOCATED(minvnodet))  ALLOCATE( minvnodet(timenc_len))
      IF(.NOT. ALLOCATED(maxvnodet))  ALLOCATE( maxvnodet(timenc_len))
      IF(.NOT. ALLOCATED(minvtimen72))ALLOCATE( minvtimen72(node_len))
      IF(.NOT. ALLOCATED(maxvtimen72))ALLOCATE( maxvtimen72(node_len))
!     -----------------------------------
!     Calculating u minimum in each node
!     -----------------------------------
      do n=1,node_len
         if(u (n, 1).lt.minun72(n) ) then
            minun72 (n)     = u (n,1)
 	    minutimen72 (n) = timenc(1)
         endif
      end do
!     -----------------------------------
!     Calculating u minimum for each time
!     -----------------------------------
      minut (1) = 99999.
      do n=1,node_len
 	 if(u (n,1).gt.-99999.and.u(n,1).lt. minut (1) ) then
           minut (1)    = u (n,1)
           minunodet(1) = n
         endif
      end do

      iret = nf_put_vara_DOUBLE(ncid, minut_id,(/nrec/),(/1/),minut)
      CALL check_err(iret)
      iret=nf_put_vara_INT(ncid, minunodet_id,(/nrec/),(/1/),minunodet)
      CALL check_err(iret)
      
!     -----------------------------------
!     Calculating u maximum in each node
!     -----------------------------------
      do n=1,node_len
         if(u (n,1).gt.maxun72(n) ) then
           maxun72 (n)     = u (n,1)
           maxutimen72 (n) = timenc(1)
         endif
      end do
!     -----------------------------------
!     Calculating u maximum for each time
!     -----------------------------------
      maxut (1) = -9999.
      do n=1,node_len
         if(u(n,1).gt.maxut(1))  then
           maxut (1)    = u (n,1)
           maxunodet(1) = n 		    
         endif
      end do
      iret = nf_put_vara_DOUBLE(ncid,maxut_id,(/nrec/),(/1/),maxut)
      CALL check_err(iret)
      iret=nf_put_vara_INT(ncid,maxunodet_id,(/nrec/),(/1/),maxunodet)
      CALL check_err(iret)
!     --------------------------------
!     Calculating v minimum each node
!     --------------------------------
      do n=1,node_len
         if(v (n,1).lt.minvn72(n) ) then
           minvn72     (n) = v (n,1)
           minvtimen72 (n) = timenc(1)
         endif
      end do
!     -----------------------------------
!     Calculating v minimum for each time
!     ------------------------------------
      minvt (1) = 99999.
      do n=1,node_len
         if(v (n,1).gt.-99999.and.v(n,1).lt. minvt (1) ) then
            minvt (1) = v (n,1)
            minvnodet(1) = n
         endif
      end do

      iret = nf_put_vara_DOUBLE(ncid, minvt_id,(/nrec/),(/1/),minvt)
      CALL check_err(iret)
      iret=nf_put_vara_INT(ncid, minvnodet_id,(/nrec/),(/1/),minvnodet)
      CALL check_err(iret)
!     ----------------------------------
!     Calculating v maximum in each node
!     ----------------------------------
      do n=1,node_len
         if(v (n,1).gt.maxvn72(n) ) then
 		   maxvn72 (n) = v (n,1)
 		   maxvtimen72(n) = timenc(1)
         endif
      end do
!     -----------------------------------
!     Calculating v maximum for each time
!     -----------------------------------
      maxvt (1) = -9999.
      do n=1,node_len
         if(v(n,1).gt.maxvt(1)) then
           maxvt (1)    = v (n,1)
           maxvnodet(1) = n 		    
         endif
      end do
      iret = nf_put_vara_DOUBLE(ncid,maxvt_id,(/nrec/),(/1/),maxvt)
      CALL check_err(iret)
      iret=nf_put_vara_INT(ncid,maxvnodet_id,(/nrec/),(/1/),maxvnodet)
      CALL check_err(iret)

!     ---------------------------
!     Closing fort.X2 netCDF file
!     ---------------------------
      IF(NREC.EQ.NDATA) THEN
         iret = nf_put_var_double(ncid, minun_id, minun72)
         CALL check_err(iret)
         iret = nf_put_var_double(ncid, maxun_id, maxun72)
         CALL check_err(iret)
         iret = nf_put_var_double(ncid, minvn_id, minvn72)
         CALL check_err(iret)
         iret = nf_put_var_double(ncid, maxvn_id, maxvn72)
         CALL check_err(iret)
         
         iret = nf_put_var_double(ncid, minutimen_id, minutimen72)
         CALL check_err(iret)
         iret = nf_put_var_double(ncid, maxutimen_id, maxutimen72)
         CALL check_err(iret)
         iret = nf_put_var_double(ncid, minvtimen_id, minvtimen72)
         CALL check_err(iret)
         iret = nf_put_var_double(ncid, maxvtimen_id, maxvtimen72)
         CALL check_err(iret)
         
      ENDIF

!     ---------------------------
!     Closing netCDF file
!     ---------------------------
      WRITE(6,*) "Closing ", FILENAME, " netCDF file ... "
      WRITE(6,*) "-----------------------------------------------------"
      
9999  iret = nf_close(ncid)
      CALL check_err(iret)
      WRITE(6,*)'end subroutine writenc72'

      end subroutine writenc72     

!=======================================================================
!     Subroutine writenc67
!     Writes a netCDF file equivalent to ADCIRC's fort.64 ascii file
!
!     Dimensions:
!     node  : number of nodes in the horizontal grid
!     nele  : number of elements in the horizontal grid
!     nope  : number of elevation specified boundary forcing segments
!     neta  : total number of elevation specified boundary nodes
!     nbou  : number of normal flow specified boundary segments
!     nvel  : total number of normal flow specified boundary nodes
!     time  : number of time records
!
!     Variables :
!     time      : time
!     u-vel     : u velocity in x or east/west direction 
!     v-vel     : v velocity in y or north/south direction
!     x         : x-coordinate or longitude
!     y         : y-coordinate or latitude
!     element   : elements in grid
!     nvdll     : number of nodes in elevation boundary segment
!     nbdv      : node numbers on elevation specified boundary segment
!     nvell     : number of nodes in normal flow specified boundary segment
!     nbvv      : node numbers on normal flow boundary segment
!     ibtype    : boundary type
!     depth     : depth in each node
!     mindepth  : minimum depth in grid
!     maxdepth  : maximum depth in grid
!     minuvelt  : minimum vertically averaged u at each time
!     maxuvelt  : maximum vertically averaged u at each time
!     minvvelt  : minimum vertically averaged v at each time
!     maxvvelt  : maximum vertically averaged v at each time
!     minuveln  : minimum vertically averaged u at each node
!     maxuveln  : maximum vertically averaged u at each node
!     minvveln  : minimum vertically averaged v at each node
!     maxvveln  : maximum vertically averaged v at each node
!     minunodet : node where minimum u occurs at each output time
!     maxunodet : node where maximum u occurs at each output time
!     minvnodet : node where minimum v occurs at each output time
!     maxvnodet : node where maximum v occurs at each output time
!     minutimen : time when minimum u occurs at each node
!     maxutimen : time when maximum u occurs at each node
!     minvtimen : time when minimum v occurs at each node
!     maxvtimen : time when maximum v occurs at each node
!
!     Attributes:
!     model             : ADCIRC
!     version           : model version number
!     grid_type         : triangular
!     rundes            : description of simulation
!     agrid             : grid ID 
!     title             : description of the dataset
!     institution       : where the original data was produced
!     source            : method of production of the original data
!     history           : provides an audit trail for modifications to the original data
!     references        : published or web-based references that describe the data or methods used to produce it
!     comments          : miscellaneous information about the data
!     contact           : user's institution/ e-mail address
!     creation_date     : date the data was created
!     modification_date : date the data was modified
!     Parameters from fort.15 input file 
!
!     Revision history:
!
!     Date            Programmer                  Description of change
!     ----            ----------                  ----------------------
!     03/30/07  Cristina Forbes, PSGS @ UNC-IMS   Wrote original code
!     03/30/08  Cristina Forbes, PSGS @ UNC-IMS   Modified code for 
!                                                 globalio & hotstart 
!                                                 from binary files
!=======================================================================
      subroutine writenc67 (nout,fileunit,timenc_len,node_len,timenc,
     &                      nt, NDATA, u,v,
     &                      zeta1, zeta2, zetad, nodecodenc,noffnc,
     &                      imhs, iths,iestp,nscoue,ivstp,nscouv, icstp, 
     &                      nscouc, ipstp, iwstp, nscoum, igep, nscouge,
     &                      igvp, nscougv, igcp, nscougc, igpp, igwp,
     &                      nscougw,
     &                      rundes, runid, agrid, irtype,
     &                      dtdp, ihot, ics, nolibf, nolifa,
     &                      nolica, nolicat, nwp, 
     &                      ncor, ntip, nws, nramp, tau0,
     &                      statim, reftim, rnday, dramp,
     &                      a00, b00, c00, h0,
     &                      slam0, sfea0, ffactor,esl, cori,
     &                      ntif, nbfr, elem_len, nface_len, nm,
     &                      x,y,depth,
     &                      nope_len, neta_len, nbdv,nopenc,netanc_len,
     &                      nvdllnc, 
     &                      nbou_len, nvel_len, nbvv,nbounc,nvelnc_len,
     &                      nvellnc,ibtypenc,
     &                      ncmode,
     &                      title, institution, source, history,
     &                      references, comments, host, convention,
     &                      contact,
     &                      base_date)


!--------------------------------------------
!     Force explicit declaration of all variables
!--------------------------------------------
      IMPLICIT NONE
!     -------------------
!     Error status return
!     -------------------
      INTEGER  iret
!     ---------
!     netCDF id
!     ---------
      INTEGER  ncid
!     -------------
!     Dimension ids
!     -------------
      INTEGER  node_dim
      INTEGER  timenc_dim
      INTEGER  elem_dim
      INTEGER  nface_dim
      INTEGER  nope_dim
      INTEGER  neta_dim
      INTEGER  nbou_dim
      INTEGER  nvel_dim
      INTEGER  nvelnc_dim
      INTEGER  netanc_dim
!     -----------------
!     Dimension lengths
!     -----------------
      INTEGER  node_len
      INTEGER  timenc_len
      INTEGER  elem_len
      INTEGER  nface_len
      INTEGER  nope_len
      INTEGER  neta_len
      INTEGER  nbou_len
      INTEGER  nvel_len
      INTEGER  nvelnc_len
      INTEGER  netanc_len
!     ------------
!     Variable ids
!     ------------
      INTEGER  timenc_id
      INTEGER  X_id
      INTEGER  Y_id
      INTEGER  DEPTH_id
      INTEGER  ELE_id
      INTEGER  mindepth_id
      INTEGER  maxdepth_id
      INTEGER  nbdv_id
      INTEGER  nbvv_id
      
      INTEGER  u_id
      INTEGER  v_id
      INTEGER  zeta1_id
      INTEGER  zeta2_id
      INTEGER  zetad_id
      INTEGER  nodecodenc_id
      INTEGER  noffnc_id

      INTEGER  nvdll_id
      INTEGER  nbounc_id
      INTEGER  ibtypenc_id
      INTEGER  nvellnc_id
      INTEGER  nvdllnc_id
      INTEGER  nopenc_id
      
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
!     ---------------------------------------------
!     Rank (number of dimensions) for each variable
!     ---------------------------------------------
      INTEGER  timenc_rank
      INTEGER  depth_rank
      INTEGER  ele_rank
      INTEGER  mindepth_rank
      INTEGER  maxdepth_rank
      INTEGER  x_rank
      INTEGER  y_rank
      INTEGER  nvdll_rank
      INTEGER  nbdv_rank
      INTEGER  nbvv_rank
      INTEGER  u_rank
      INTEGER  v_rank
      INTEGER  zeta1_rank
      INTEGER  zeta2_rank
      INTEGER  zetad_rank
      INTEGER  nodecodenc_rank
      INTEGER  noffnc_rank

      INTEGER  nbounc_rank
      INTEGER  ibtypenc_rank
      INTEGER  nvellnc_rank
      INTEGER  nvdllnc_rank
      INTEGER  nopenc_rank

      INTEGER imhs_rank
      INTEGER iths_rank
      INTEGER IESTP_rank
      INTEGER NSCOUE_rank
      INTEGER IVSTP_rank
      INTEGER NSCOUV_rank
      INTEGER ICSTP_rank
      INTEGER NSCOUC_rank
      INTEGER IPSTP_rank
      INTEGER IWSTP_rank
      INTEGER NSCOUM_rank
      INTEGER IGEP_rank
      INTEGER NSCOUGE_rank
      INTEGER IGVP_rank
      INTEGER NSCOUGV_rank
      INTEGER IGCP_rank
      INTEGER NSCOUGC_rank
      INTEGER IGPP_rank
      INTEGER IGWP_rank
      INTEGER NSCOUGW_rank

      PARAMETER (timenc_rank = 1)
      PARAMETER (depth_rank = 1)
      PARAMETER (ele_rank = 2)
      PARAMETER (mindepth_rank = 0)
      PARAMETER (maxdepth_rank = 0)
      PARAMETER (x_rank = 1)
      PARAMETER (y_rank = 1)
      PARAMETER (nbdv_rank = 2)
      PARAMETER (nbvv_rank = 2)
      PARAMETER (u_rank = 2)
      PARAMETER (v_rank = 2)
      PARAMETER (zeta1_rank = 2)
      PARAMETER (zeta2_rank = 2)
      PARAMETER (zetad_rank = 2)
      PARAMETER (nodecodenc_rank = 1)
      PARAMETER (noffnc_rank = 1)

      PARAMETER (nvdll_rank = 1)
      PARAMETER (nbounc_rank = 1)
      PARAMETER (ibtypenc_rank = 1)
      PARAMETER (nvellnc_rank = 1)
      PARAMETER (nvdllnc_rank=1)
      PARAMETER (nopenc_rank=1)

      PARAMETER (imhs_rank = 0)
      PARAMETER (iths_rank = 0)
      PARAMETER (IESTP_rank = 0)
      PARAMETER (NSCOUE_rank = 0)
      PARAMETER (IVSTP_rank = 0)
      PARAMETER (NSCOUV_rank = 0)
      PARAMETER (ICSTP_rank = 0)
      PARAMETER (NSCOUC_rank = 0)
      PARAMETER (IPSTP_rank = 0)
      PARAMETER (IWSTP_rank = 0)
      PARAMETER (NSCOUM_rank = 0)
      PARAMETER (IGEP_rank = 0)
      PARAMETER (NSCOUGE_rank = 0)
      PARAMETER (IGVP_rank = 0)
      PARAMETER (NSCOUGV_rank = 0)
      PARAMETER (IGCP_rank = 0)
      PARAMETER (NSCOUGC_rank = 0)
      PARAMETER (IGPP_rank = 0)
      PARAMETER (IGWP_rank = 0)
      PARAMETER (NSCOUGW_rank = 0)
      
!     ---------------
!     Variable shapes
!     ---------------
      INTEGER  timenc_dims(timenc_rank)
      INTEGER  x_dims(x_rank)
      INTEGER  y_dims(y_rank)
      INTEGER  depth_dims(depth_rank)
      INTEGER  ele_dims(ele_rank)
      INTEGER  nbdv_dims(nbdv_rank)
      INTEGER  nbvv_dims(nbvv_rank)
      INTEGER  u_dims(u_rank)
      INTEGER  v_dims(v_rank)
      INTEGER  zeta1_dims(zeta1_rank)
      INTEGER  zeta2_dims(zeta2_rank)
      INTEGER  zetad_dims(zetad_rank)
      INTEGER  nodecodenc_dims(nodecodenc_rank)
      INTEGER  noffnc_dims(noffnc_rank)

      INTEGER  nvdll_dims(nvdll_rank)
      INTEGER  nbounc_dims(nbounc_rank)
      INTEGER  ibtypenc_dims(ibtypenc_rank)
      INTEGER  nvellnc_dims(nvellnc_rank)
      INTEGER  nvdllnc_dims(nvdllnc_rank)
      INTEGER  nopenc_dims(nopenc_rank)

!     --------------
!     Data variables
!     --------------
      CHARACTER*5 FILEFORT
      CHARACTER*10 FILENAME
      CHARACTER*3 FILEEXT
      INTEGER   t,ntif,nbfr,nt
      INTEGER   irtype, ihot, ics, nolibf, nolifa, nolica, nolicat, nwp
      INTEGER   ntip, nws, nramp,ncor
      INTEGER   nout
      INTEGER   element (nface_len , elem_len)
      INTEGER   nbdv (nope_len, neta_len)
      INTEGER   nbvv (nbou_len, nvel_len)
      
      INTEGER   nvellnc(nbou_len)
      INTEGER   nbounc (nvel_len)
      INTEGER   ibtype(nbou_len)
      INTEGER   ibtypenc(nbou_len)
      INTEGER   nvdllnc (nope_len)
      INTEGER   nopenc (neta_len)
      
      INTEGER   nm (elem_len, nface_len)
      INTEGER   fileunit
      REAL*8    timenc(timenc_len)
      REAL(sz)  x(node_len)
      REAL(sz)  y(node_len)
      REAL(sz)  depth(node_len)
      REAL(sz)  u(node_len, timenc_len)
      REAL(sz)  v(node_len, timenc_len)
      REAL(sz)  zeta1 (node_len, timenc_len)
      REAL(sz)  zeta2 (node_len, timenc_len)
      REAL(sz)  zetad (node_len, timenc_len)
      INTEGER   nodecodenc(node_len)
      INTEGER   noffnc(elem_len)

      real*8  mindepth
      real*8  maxdepth
      REAL*8 dtdp, statim,reftim, rnday, dramp, tau0
      REAL*8 a00, b00, c00,h0
      REAL(sz) slam0, sfea0, ffactor,esl, cori
      REAL(sz) slam0deg, sfea0deg
      CHARACTER(80) rundes, runid, agrid
      INTEGER, SAVE  :: FIRST67 = 1 
      INTEGER, SAVE  :: FIRST68 = 1 
      INTEGER NDATA
      INTEGER, SAVE  :: NREC = 0
      INTEGER ncmode
      
      INTEGER imhs
      INTEGER iths
      INTEGER IESTP
      INTEGER NSCOUE
      INTEGER IVSTP
      INTEGER NSCOUV
      INTEGER ICSTP
      INTEGER NSCOUC
      INTEGER IPSTP
      INTEGER IWSTP
      INTEGER NSCOUM
      INTEGER IGEP
      INTEGER NSCOUGE
      INTEGER IGVP
      INTEGER NSCOUGV
      INTEGER IGCP
      INTEGER NSCOUGC
      INTEGER IGPP
      INTEGER IGWP
      INTEGER NSCOUGW

      

!     -----------------
!     Attribute vectors
!     -----------------
      INTEGER intval(1),i,j
      INTEGER counti(1), starti(1), n
      INTEGER count(2), start(2)
      REAL    realval(1)
!     -----------------
!     Attribute general
!     -----------------
      CHARACTER(*) :: title, institution, source, history
      CHARACTER(*) :: references, comments, host, convention, contact
      CHARACTER(*) :: base_date

!     -----------------
!      date_string variables for time attribute
!     -----------------
      character date_string*40
      character now_date*8
      character big_ben*10
      character zone*5
      integer values(8)
!     ---------------------------------
!      inquiring and reading variables
!     ---------------------------------

      INTEGER NDIMS, NVARS, NGATTS, UNLIMDIMID, RECID,NRECS
      CHARACTER*4 recname
      
!     -------------------------
!     Save dims, ids, variables
!     -------------------------
      save node_dim,  elem_dim, nface_dim
      save nope_dim, neta_dim, nbou_dim, nvel_dim
      save u_id, v_id, zeta1_id, zeta2_id, zetad_id
      save nodecodenc_id, noffnc_id
      save X_id, Y_id

      save DEPTH_id,ELE_id, mindepth_id, maxdepth_id, nbdv_id, nbvv_id
      save timenc_dim, timenc_id
      save mindepth, maxdepth

      save nvelnc_dim
      save nvdll_id, nbounc_id, ibtypenc_id, nvellnc_id
      save netanc_dim
      save nopenc_id,nvdllnc_id


      save imhs_id
      save iths_id
      save IESTP_id
      save NSCOUE_id
      save IVSTP_id
      save NSCOUV_id
      save ICSTP_id
      save NSCOUC_id
      save IPSTP_id
      save IWSTP_id
      save NSCOUM_id
      save IGEP_id
      save NSCOUGE_id
      save IGVP_id
      save NSCOUGV_id
      save IGCP_id
      save NSCOUGC_id
      save IGPP_id
      save IGWP_id
      save NSCOUGW_id


!     ---------------------------
!     Opening fort.67 netCDF file
!     ---------------------------
      WRITE(6,*) "-----------------------------------------------------"
           PRINT *, " *********  node_len = ", node_len
           PRINT *, " *********  elem_len = ", elem_len
       
!      do i=1, node_len
!            PRINT *, "i ",i, " nodecodenc(i) = ", nodecodenc(i)
!      end do
!            PRINT *, " nodecodenc(11) = ", nodecodenc(11)
!      do i=1, elem_len
!            PRINT *, "i ",i, " noffnc(i) = ", noffnc(i)
!      end do
!            PRINT *, " noffnc(11) = ", noffnc(i)

!      IF(FIRST67.EQ.1.and.ncmode.NE.0) THEN
         
!      ENDIF
       
!     --------------------------------
!     Switch order in array for NETCDF
!     --------------------------------
      do i=1, elem_len
         do j=1, nface_len
            element(j,i) = nm (i,j)
         end do
      end do
      
!     -----------------
!     Enter define mode
!     -----------------
         FILEFORT = 'fort.'
         FILEEXT = '.nc'
         WRITE(FILENAME,'(A5,I2,A3)')FILEFORT,FILEUNIT,FILEEXT
         WRITE(6,*) "Opening ",FILENAME,  " netCDF file ... "
         
!       IF((FIRST67.EQ.1.AND.IHOT.EQ.0.AND.NCMODE.NE.2).OR.
!     &    (FIRST67.EQ.1.AND.IHOT.NE.0.AND.NCMODE.NE.2.AND.
!     &     nout.EQ.-3))THEN

       IF( ((FIRST67.EQ.1.AND.IHOT.EQ.0.AND.NCMODE.NE.2).OR.
     &     (FIRST68.EQ.1.AND.IHOT.EQ.0.AND.NCMODE.NE.2)).OR.
     &    ((FIRST67.EQ.1.AND.IHOT.NE.0.AND.NCMODE.NE.2.AND.
     &     nout.EQ.-3).OR.
     &     (FIRST68.EQ.1.AND.IHOT.NE.0.AND.NCMODE.NE.2.AND.
     &     nout.EQ.-3)))THEN



           iret = nf_create(FILENAME, NF_CLOBBER, ncid)
           CALL check_err(iret)
           IF(FILEUNIT.EQ.67)FIRST67=0
           IF(FILEUNIT.EQ.68)FIRST68=0
           WRITE(6,*) "FIRST67 ",FIRST67
           WRITE(6,*) "FIRST68 ",FIRST68
         
!     --------------------------------
!     Convert to degress for NETCDF I/O
!     --------------------------------
           SLAM0DEG=SLAM0
           SFEA0DEG=SFEA0
           IF (NCMODE.EQ.1) THEN
               SLAM0DEG=SLAM0*RAD2DEG
               SFEA0DEG=SFEA0*RAD2DEG
           ENDIF
!     -----------------
!     Define dimensions
!     -----------------
         WRITE(6,*) "Define dimensions "
      iret = nf_def_dim(ncid, 'node', node_len, node_dim)
      CALL check_err(iret)
      iret = nf_def_dim(ncid, 'nele', elem_len, elem_dim)
      CALL check_err(iret)
      iret = nf_def_dim(ncid, 'nface', 3, nface_dim)
      CALL check_err(iret)
      if(nope_len.ne.0) then
         iret = nf_def_dim(ncid, 'nope', nope_len, nope_dim)
         call check_err(iret)
      endif
      if(neta_len.ne.0) then
         iret = nf_def_dim(ncid, 'neta', neta_len, neta_dim)
         call check_err(iret)
      endif
      if(nbou_len.ne.0) then
         iret = nf_def_dim(ncid, 'nbou', nbou_len, nbou_dim)
         call check_err(iret)
      endif
      if(nvel_len.ne.0) then
         iret = nf_def_dim(ncid, 'nvel', nvel_len, nvel_dim)
         call check_err(iret)
      endif
      
      iret = nf_def_dim(ncid, 'time', nf_unlimited, timenc_dim)
      CALL check_err(iret)
      
!     ================
!     Define variables
!     ================
         WRITE(6,*) "Define variables "
!     ----------------
!     Define time
!     ----------------
      timenc_dims(1) = timenc_dim
      iret = nf_def_var(ncid,'time',NF_DOUBLE, timenc_rank, timenc_dims, 
     &                  timenc_id)
      CALL check_err(iret)
      
!     -------------------
!     Define x
!     -------------------
      x_dims(1) = node_dim
      iret = nf_def_var(ncid, 'x', NF_DOUBLE, x_rank, x_dims, x_id)
      CALL check_err(iret)
!     -------------------
!     Define Y coordinate
!     --------------------
      Y_dims(1) = node_dim
      iret = nf_def_var(ncid, 'y', NF_DOUBLE, Y_rank, Y_dims, Y_id)
      CALL check_err(iret)
!     -------------------
!     Define elements
!     --------------------
      ELE_dims(1) = nface_dim
      ELE_dims(2) = elem_dim
      iret = nf_def_var(ncid, 'element',NF_INT,ELE_rank,ELE_dims,ELE_id)
      CALL check_err(iret)
!     -----------------------------------------------------------------
!     Define elevation specified boundary forcing segments information
!     -----------------------------------------------------------------
      if(nope_len.ne.0.and.neta_len.ne.0) then
        nvdllnc_dims(1) = nope_dim      
        iret = nf_def_var(ncid, 'nvdll',NF_INT, nvdllnc_rank, 
     &                  nvdllnc_dims, nvdllnc_id)
        CALL check_err(iret)

        nopenc_dims(1) = neta_dim      
        iret = nf_def_var(ncid, 'nbbv',NF_INT, nopenc_rank, nopenc_dims, 
     &                  nopenc_id)
        CALL check_err(iret)

      endif

!     ----------------------------------------------------------
!     Define normal flow specified boundary segments information
!     ----------------------------------------------------------
      if(nbou_len.ne.0.and.nvel_len.ne.0) then

        nvellnc_dims(1) = nbou_dim      
        iret = nf_def_var(ncid, 'nvell',NF_INT, nvellnc_rank, 
     &                  nvellnc_dims, nvellnc_id)
        CALL check_err(iret)
        
        ibtypenc_dims(1) = nbou_dim      
        iret = nf_def_var(ncid, 'ibtype',NF_INT, ibtypenc_rank, 
     &                  ibtypenc_dims, ibtypenc_id)
        CALL check_err(iret)
        nbounc_dims(1) = nvel_dim      
        iret = nf_def_var(ncid, 'nbvv',NF_INT, nbounc_rank, nbounc_dims, 
     &                  nbounc_id)
        CALL check_err(iret)
     
       endif

         WRITE(6,*) "Define depth "
!     -------------------
!     Define Z coordinate
!     --------------------
      DEPTH_dims(1) = node_dim
      iret = nf_def_var(ncid, 'depth', NF_DOUBLE,DEPTH_rank, DEPTH_dims, 
     &                  DEPTH_id)
      CALL check_err(iret)
!     --------------------
!     Define minimum depth
!     --------------------
      iret = nf_def_var(ncid, 'mindepth', NF_DOUBLE, mindepth_rank, 0, 
     &                  mindepth_id)
      CALL check_err(iret)
!     --------------------
!     Define maximum depth
!     --------------------
      iret = nf_def_var(ncid, 'maxdepth', NF_DOUBLE, maxdepth_rank, 0, 
     &                  maxdepth_id)
      CALL check_err(iret)

!     ----------------
!     Define zeta
!     ----------------
      zeta1_dims(1) = node_dim
      zeta1_dims(2) = timenc_dim
      
      iret = nf_def_var(ncid, 'zeta1',NF_DOUBLE, zeta1_rank, zeta1_dims, 
     &                  ZETA1_id)
      
      CALL check_err(iret)
!     ----------------
!     Define zeta
!     ----------------
      zeta2_dims(1) = node_dim
      zeta2_dims(2) = timenc_dim
      
      iret = nf_def_var(ncid, 'zeta2',NF_DOUBLE, zeta2_rank, zeta2_dims, 
     &                  ZETA2_id)
      
      CALL check_err(iret)

!     ----------------
!     Define zeta
!     ----------------
      zetad_dims(1) = node_dim
      zetad_dims(2) = timenc_dim
      
      iret = nf_def_var(ncid, 'zetad',NF_DOUBLE, zetad_rank, zetad_dims, 
     &                  ZETAd_id)
      
      CALL check_err(iret)


         WRITE(6,*) "Define vel "
!     ----------------
!     Define u vel
!     ----------------
      u_dims(1) = node_dim
      u_dims(2) = timenc_dim
      iret = nf_def_var(ncid, 'u-vel', NF_DOUBLE, u_rank, u_dims, 
     &                  u_id)
      CALL check_err(iret)

!     ----------------
!     Define v vel
!     ----------------
      v_dims(1) = node_dim
      v_dims(2) = timenc_dim
      iret = nf_def_var(ncid, 'v-vel', NF_DOUBLE, v_rank, v_dims, 
     &                  v_id)
      CALL check_err(iret)

!      WRITE(6,*)"Define MIN u/v at times"

         WRITE(6,*) "hotstart parameters "
!     --------------------
!     Define hotstart parameters
!     --------------------
      iret = nf_def_var(ncid, 'nodecode', NF_INT, nodecodenc_rank, 
     &                  node_dim, nodecodenc_id)
      CALL check_err(iret)
         WRITE(6,*) "noff: "
      iret = nf_def_var(ncid, 'noff', NF_INT, noffnc_rank, elem_dim, 
     &                  noffnc_id)
      CALL check_err(iret)
         WRITE(6,*) "imhs: "
      iret = nf_def_var(ncid, 'imhs', NF_INT, imhs_rank, 0, 
     &                  imhs_id)
      CALL check_err(iret)
      iret = nf_def_var(ncid, 'iths', NF_INT, iths_rank, 0, 
     &                  iths_id)
      CALL check_err(iret)
      iret = nf_def_var(ncid, 'iestp', NF_INT, iestp_rank, 0, 
     &                  iestp_id)
      CALL check_err(iret)
      iret = nf_def_var(ncid, 'nscoue', NF_INT, nscoue_rank, 0, 
     &                  nscoue_id)
      CALL check_err(iret)

         WRITE(6,*) "ivstp: "
      iret = nf_def_var(ncid, 'ivstp', NF_INT, ivstp_rank, 0, 
     &                  ivstp_id)
      CALL check_err(iret)
      iret = nf_def_var(ncid, 'nscouv', NF_INT, nscouv_rank, 0, 
     &                  nscouv_id)
      CALL check_err(iret)
      iret = nf_def_var(ncid, 'icstp', NF_INT, icstp_rank, 0, 
     &                  icstp_id)
      CALL check_err(iret)
      iret = nf_def_var(ncid, 'nscouc', NF_INT, nscouc_rank, 0, 
     &                  nscouc_id)
      CALL check_err(iret)

      iret = nf_def_var(ncid, 'ipstp', NF_INT, ipstp_rank, 0, 
     &                  ipstp_id)
      CALL check_err(iret)
      iret = nf_def_var(ncid, 'iwstp', NF_INT, iwstp_rank, 0, 
     &                  iwstp_id)
      CALL check_err(iret)
      iret = nf_def_var(ncid, 'nscoum', NF_INT, nscoum_rank, 0, 
     &                  nscoum_id)
      CALL check_err(iret)
      iret = nf_def_var(ncid, 'igep', NF_INT, igep_rank, 0, 
     &                  igep_id)
      CALL check_err(iret)
      
      iret = nf_def_var(ncid, 'nscouge', NF_INT, nscouge_rank, 0, 
     &                  nscouge_id)
      CALL check_err(iret)
      iret = nf_def_var(ncid, 'igvp', NF_INT, igvp_rank, 0, 
     &                  igvp_id)
      CALL check_err(iret)
      iret = nf_def_var(ncid, 'nscougv', NF_INT, nscougv_rank, 0, 
     &                  nscougv_id)
      CALL check_err(iret)
      iret = nf_def_var(ncid, 'igcp', NF_INT, igcp_rank, 0, 
     &                  igcp_id)
      CALL check_err(iret)
         WRITE(6,*) "nscougc: "

      iret = nf_def_var(ncid, 'nscougc', NF_INT, nscougc_rank, 0, 
     &                  nscougc_id)
      CALL check_err(iret)
      iret = nf_def_var(ncid, 'igpp', NF_INT, igpp_rank, 0, 
     &                  igpp_id)
      CALL check_err(iret)
      iret = nf_def_var(ncid, 'igwp', NF_INT, igwp_rank, 0, 
     &                  igwp_id)
      CALL check_err(iret)
      iret = nf_def_var(ncid, 'nscougw', NF_INT, nscougw_rank, 0, 
     &                  nscougw_id)
      CALL check_err(iret)
!     -----------------
!     Assign attributes
!     -----------------
         WRITE(6,*) "Assign attributes: "

      if(ics.eq.2) then
      iret = nf_put_att_text(ncid, X_id, 'long_name', 9, 'longitude')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, X_id, 'standard_name', 9, 
     &                      'longitude')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, X_id, 'units', 12, VAR_DEGE)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, X_id, 'positive', 4, VAR_EAST)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, X_id, 'missing_values', NF_DOUBLE,
     &                        1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, X_id, '_FillValue', NF_DOUBLE, 1, 
     &                       doubleval)
      CALL check_err(iret)      
      
      iret = nf_put_att_text(ncid, Y_id, 'long_name', 8, 'latitude')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, Y_id, 'standard_name', 8, 
     &                       'latitude')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, Y_id, 'units', 13, VAR_DEGN)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, Y_id, 'positive', 5, VAR_NORTH)
      CALL check_err(iret)
       
      iret = nf_put_att_double(ncid, Y_id, 'missing_values', NF_DOUBLE, 
     &       1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, Y_id, '_FillValue', NF_DOUBLE, 1, 
     &       doubleval)
      CALL check_err(iret)

      else
      
      iret = nf_put_att_text(ncid, X_id, 'long_name', 22, 
     &       'Cartesian coordinate x')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, X_id, 'standard_name', 12, 
     &                      'x_coordinate')
      CALL check_err(iret)
      
      if (G.LT.11.D0) then
         iret = nf_put_att_text(ncid, X_id, 'units', 6,  VAR_METS)
         CALL check_err(iret)
      else
         iret = nf_put_att_text(ncid, X_id, 'units', 4,  VAR_FEET)
         CALL check_err(iret)
      endif
      
      iret = nf_put_att_text(ncid, X_id, 'positive', 5,  VAR_RG)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, X_id, 'missing_values', NF_DOUBLE, 
     &                       1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, X_id, '_FillValue', NF_DOUBLE, 1, 
     &                       doubleval)
      CALL check_err(iret)
!     ------------------------------
!     Define y-coordinate attributes
!     ------------------------------     
      
      iret = nf_put_att_text(ncid, Y_id, 'long_name', 22, 
     &       'Cartesian coordinate y')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, Y_id, 'standard_name', 12, 
     &                       'y_coordinate')
      CALL check_err(iret)
      if (G.LT.11.D0) then
         iret = nf_put_att_text(ncid, Y_id, 'units', 6, VAR_METS)
         CALL check_err(iret)
      else
         iret = nf_put_att_text(ncid, Y_id, 'units', 4, VAR_FEET)
         CALL check_err(iret)
      endif
      iret = nf_put_att_text(ncid, Y_id, 'positive', 2, VAR_UP)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, Y_id, 'missing_values', 
     &       NF_DOUBLE, 1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, Y_id, '_FillValue', NF_DOUBLE, 1, 
     &       doubleval)
      CALL check_err(iret)
      
      endif

!     ------------------------------
!     Define bathymetry attributes
!     ------------------------------

      iret = nf_put_att_text(ncid, DEPTH_id,'long_name',10,'bathymetry')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, DEPTH_id, 'standard_name', 5,'depth')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, DEPTH_id, 'units', 6, VAR_METS)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, DEPTH_id, 'missing_values', 
     &                      NF_DOUBLE,  1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid,DEPTH_id, '_FillValue',NF_DOUBLE, 1, 
     &                       doubleval)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, DEPTH_id, 'positive', 4, VAR_DW)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, ELE_id, 'long_name', 7, 'element')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, ELE_id, 'units', 14, VAR_NOND)
      intval(1) = -1
      
!    ---------------------------
!     minimum/maximum bathymetry
!    ---------------------------

      iret = nf_put_att_text(ncid, mindepth_id, 'long_name', 18, 
     &                       'minimum bathymetry')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, mindepth_id, 'standard_name', 13,
     &                'minimum_depth')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, mindepth_id, 'units', 6, VAR_METS)
      CALL check_err(iret)
      
      iret = nf_put_att_text(ncid, maxdepth_id, 'long_name', 18, 
     &                       'maximum bathymetry')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, maxdepth_id, 'standard_name', 13,
     &           'maximum_depth')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, maxdepth_id, 'units', 6, VAR_METS)
      CALL check_err(iret)

!    ---------------------------
!     Hotstart variables
!    ---------------------------
         WRITE(6,*) "Assign nodecodenc_id: "


       iret = nf_put_att_text(ncid, nodecodenc_id, 'long_name', 82, 
     &                       'node code indicating whether the node is
     & presently wet (active) or dry (inactive)')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, nodecodenc_id, 'standard_name', 18,
     &                'wet_dry_indicator')
      CALL check_err(iret)
         WRITE(6,*) "Assign noffnc_id: "
      iret = nf_put_att_text(ncid, noffnc_id, 'long_name', 34, 
     &                       'Element based wetting/drying flag')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, noffnc_id, 'standard_name', 28,
     &                'element_wetting/drying_flag')
      CALL check_err(iret)

 
 
         WRITE(6,*) "Assign imhs_id: "

      iret = nf_put_att_text(ncid, imhs_id, 'long_name', 18, 
     &                       'model type')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, imhs_id, 'standard_name', 10,
     &                'model_type')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, iths_id, 'long_name', 59, 
     &                       'model time step number since the beginning
     & of the model run')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, iths_id, 'standard_name', 15,
     &           'model_time_step')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, iestp_id, 'long_name', 176, 
     &                       'line number (for ASCII output) or record
     & number (for binary output) of the most recent entry in the
     & elevation time series at specified elevation recording stations
     & output file')
      CALL check_err(iret)
    
      iret = nf_put_att_text(ncid, iestp_id, 'standard_name', 54,
     &      'line/record_number_of_last_entry_in_elev_rec_stations')
      CALL check_err(iret)
      
      iret = nf_put_att_text(ncid, nscoue_id, 'long_name', 150, 
     &                       'time step counter to determine when the
     & next entry will be written to the elevation time series at
     & specified elevation recording Stations output file')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, nscoue_id, 'standard_name', 15,
     &           'time_step_counter_for_next_entry_elev_rec_stations)')
      CALL check_err(iret)
      
      iret = nf_put_att_text(ncid, ivstp_id, 'long_name', 190, 
     &                       ' line number (for ASCII output) or record 
     & number (for binary output) of the most recent entry in the
     & depth-averaged velocity time series at specified velocity 
     & recording stations output file')
      CALL check_err(iret)
      
      iret = nf_put_att_text(ncid, ivstp_id, 'standard_name', 48,
     &               'record_number_of_last_entry_in_vel_rec_stations')
      CALL check_err(iret)

       iret = nf_put_att_text(ncid, nscouv_id, 'long_name', 165, 
     &                       'time step counter to determine when the
     & next entry will be written to the Depth-averaged Velocity Time
     & series at specified velocity recording stations output file.')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, nscouv_id, 'standard_name', 50,
     &           'time_step_counter_for_next_entry_vel_rec_stations')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, icstp_id, 'long_name', 190, 
     &                       'line number (for ASCII output) or record
     & number (for binary output) of the most recent entry in the scalar
     & concentration time series at specified concentration recording
     & stations output file')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, icstp_id, 'standard_name', 49,
     &               'record_number_of_last_entry_in_conc_rec_stations')
      CALL check_err(iret)
     
      iret = nf_put_att_text(ncid, nscouc_id, 'long_name', 165, 
     &                       'time step counter to determine when the
     & next entry will be written to the scalar concentration Time 
     &series at specified concentration Recording Stations output file')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, nscouc_id, 'standard_name', 51,
     &           'time_step_counter_for_next_entry_conc_rec_stations')
      CALL check_err(iret)

         WRITE(6,*) "Assign ipstp_id: "

      iret = nf_put_att_text(ncid, ipstp_id, 'long_name', 192, 
     &                       'line number (for ASCII output) or record
     & number (for binary output) of the most recent entry in the 
     &atmospheric pressure time series at specified meteorological 
     &recording stations')
      CALL check_err(iret)
     
      iret = nf_put_att_text(ncid, ipstp_id, 'standard_name', 57,
     &       'record_number_of_last_entry_of_atm_press_at_rec_stations')
      CALL check_err(iret)
     
         WRITE(6,*) "Assign iwstp_id: "
      iret = nf_put_att_text(ncid, iwstp_id, 'long_name', 190, 
     &                       'line number (for ASCII output) or record
     & number (for binary output) of the most recent entry in the Wind 
     & Velocity time series at specified meteorological recording 
     &stations')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, iwstp_id, 'standard_name', 56,
     &       'record_number_of_last_entry_of_wind_vel_at_rec_stations')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, nscoum_id, 'long_name', 244, 
     &                       'time step counter to determine when the
     & next entry will be written to the atmospheric pressure time 
     &series at specified meteorological recording stations and wind 
     &velocity time series at specified meteorological recording 
     &stations output files')
      CALL check_err(iret)

         WRITE(6,*) "Assign nscoum_id: "
      iret = nf_put_att_text(ncid, nscoum_id, 'standard_name', 60,
     &   'time_step_counter_of_atm_press_and_wind_vel_at_rec_stations')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, igep_id, 'long_name', 165, 
     &                       'line number (for ASCII output) or record 
     &number (for binary output) of the most recent entry in the 
     &Elevation time series at All Nodes in the Model Grid output file')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, igep_id, 'standard_name', 51,
     &       'record_number_of_last_entry_of_elev_at_model_nodes')
      CALL check_err(iret)

         WRITE(6,*) "Assign nscouge_id: "
      iret = nf_put_att_text(ncid, nscouge_id, 'long_name', 139, 
     &                       'time step counter to determine when the
     & next entry will be written to the  Elevation time series at All 
     &Nodes in the Model Grid output file')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, nscouge_id, 'standard_name', 41,
     &   'time_step_counter_of_elev_at_model_nodes')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, igvp_id, 'long_name', 179, 
     &                       'line number (for ASCII output) or record 
     &number (for binary output) of the most recent entry in the 
     &Depth-averaged velocity time series at All Nodes in the Model 
     &Grid output file')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, igvp_id, 'standard_name', 50,
     &       'record_number_of_last_entry_of_vel_at_model_nodes')
      CALL check_err(iret)

         WRITE(6,*) "Assign nscougv_id: "
      iret = nf_put_att_text(ncid, nscougv_id, 'long_name', 153, 
     &                       'time step counter to determine when the 
     &next entry will be written to the Depth-averaged velocity time 
     &series at All Nodes in the Model Grid output file')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, nscougv_id, 'standard_name', 40,
     &   'time_step_counter_of_vel_at_model_nodes')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, igcp_id, 'long_name', 176, 
     &                       'line number (for ASCII output) or record 
     &number (for binary output) of the most recent entry in the 
     &scalar Concentration time series at All Nodes in the Model 
     &Grid output file')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, igcp_id, 'standard_name', 51,
     &       'record_number_of_last_entry_of_conc_at_model_nodes')
      CALL check_err(iret)

         WRITE(6,*) "Assign nscougc_id: "
      iret = nf_put_att_text(ncid, nscougc_id, 'long_name', 150, 
     &                       'time step counter to determine when the 
     &next entry will be written to the scalar Concentration time 
     &series at All Nodes in the Model Grid output file')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, nscougc_id, 'standard_name', 41,
     &   'time_step_counter_of_conc_at_model_nodes')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, igpp_id, 'long_name', 176, 
     &                       'line number (for ASCII output) or record 
     &number (for binary output) of the most recent entry in the 
     &atmospheric pressure time series at All Nodes in the Model 
     &Grid output file')
      CALL check_err(iret)

         WRITE(6,*) "Assign igpp_id: "
      iret = nf_put_att_text(ncid, igpp_id, 'standard_name', 61,
     &       'record_number_of_last_entry_of_atm_press_at_model_nodes')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, igwp_id, 'long_name', 179, 
     &                       'line number (for ASCII output) or record 
     &number (for binary output) of the most recent entry in the 
     &Wind Stress or velocity time series at All Nodes in the Model 
     &Grid output file')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, igwp_id, 'standard_name', 55,
     &       'record_number_of_last_entry_of_wind_vel_at_model_nodes')
      CALL check_err(iret)

         WRITE(6,*) "Assign nscougw_id: "
      iret = nf_put_att_text(ncid, nscougw_id, 'long_name', 222, 
     &                       'time step counter to determine when the 
     &next entry will be written to the atmospheric pressure time series 
     &at All Nodes in the Model Grid and Wind Stress or velocity time 
     &series at All Nodes in the Model Grid output files')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, nscougw_id, 'standard_name', 59,
     &   'time_step_counter_of_atm_press_and_wind_vel_at_model_nodes')
      CALL check_err(iret)


!     -----------------------------------------------------------------------
!     Define node numbers on elevation specified boundary segments attributes
!     -----------------------------------------------------------------------      
      iret = nf_put_att_text(ncid, nopenc_id, 'long_name', 53, 
     &       'node numbers on elevation specified boundary segments')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, nopenc_id, 'units', 14, 
     &                      VAR_NOND)
      CALL check_err(iret)
      
      iret = nf_put_att_text(ncid, nvdllnc_id, 'long_name', 60, 
     &'number of nodes in each elevation specified boundary segment')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, nvdllnc_id, 'units', 14, 
     &                      VAR_NOND)
      CALL check_err(iret)
!     ---------------------------------------------------------------
!     Define node numbers on normal flow boundary segments attributes
!     ----------------------------------------------------------------  

      iret = nf_put_att_text(ncid, nbounc_id, 'long_name', 59, 
     & 'node numbers in each normal flow specified boundary segment')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, nbounc_id, 'units', 14, 
     &                      VAR_NOND)
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, ibtypenc_id, 'long_name', 13, 
     &       'boundary type')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, ibtypenc_id, 'units', 14, 
     &                      VAR_NOND)
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, nvellnc_id, 'long_name', 62, 
     &'number of nodes in each normal flow specified boundary segment')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, nvellnc_id, 'units', 14, 
     &                      VAR_NOND)
      CALL check_err(iret)

!    -----------------------
!     Define time attributes
!    -----------------------
      iret = nf_put_att_text(ncid,timenc_id,'long_name',10,'model time')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, timenc_id, 'standard_name', 4,'time')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, timenc_id, 'units', 23,VAR_TUN)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, timenc_id, 'base_date',
     . (LEN(TRIM(base_date(1:INDEX(base_date,"!")-1)))), base_date)
      CALL check_err(iret)

!    -----------
!     Define u velocity
!    -----------
         WRITE(6,*) "Define u velocity: "

!    --------------------------------------------------
!     determine variable name according to coordinates
!    ---------------------------------------------------
      if (ics.eq.2) then
          iret = nf_put_att_text(ncid, u_id, 'long_name', 32, 
     &                       'vertically averaged e/w velocity')
       else
         iret = nf_put_att_text(ncid, u_id, 'long_name', 43, 
     &            'vertically averaged velocity in x-direction')
      endif
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, u_id, 'standard_name', 10, 
     &                      'u_velocity')
      CALL check_err(iret)
!    -------------------------
!     determine variable units
!    -------------------------
      if (G.LT.11.D0) then
          iret = nf_put_att_text(ncid, u_id, 'units',13, VAR_MPS)
      else
          iret = nf_put_att_text(ncid, u_id, 'units',11, VAR_FPS)
      endif
      CALL check_err(iret)
!    --------------------------------------------------
!     determine variable name according to coordinates
!    ---------------------------------------------------
      if (ics.eq.2) then
          iret = nf_put_att_text(ncid, u_id, 'positive', 4, VAR_EAST)
      else
          iret = nf_put_att_text(ncid, u_id, 'positive', 5, VAR_RG)
      endif
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, u_id, 'missing_values', 
     &       NF_DOUBLE, 1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, u_id, '_FillValue', NF_DOUBLE, 1, 
     &                       doubleval)
      CALL check_err(iret)
      
      iret = nf_put_att_double(ncid, u_id, 'dry_Value', NF_DOUBLE, 1, 
     &                       doubleval)
      CALL check_err(iret)
!    -----------
!     Define v velocity
!    -----------
!    --------------------------------------------------
!     determine variable name according to coordinates
!    ---------------------------------------------------
      if (ics.eq.2) then
          iret = nf_put_att_text(ncid, v_id, 'long_name', 32, 
     &           'vertically averaged n/s velocity')
      else
          iret = nf_put_att_text(ncid, v_id, 'long_name', 43, 
     &           'vertically averaged velocity in y-direction')
      endif
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, v_id, 'standard_name', 10, 
     &                      'v_velocity')
      CALL check_err(iret)
!    -------------------------
!     determine variable units
!    -------------------------
      if (G.LT.11.D0) then
         iret = nf_put_att_text(ncid, v_id, 'units',13, VAR_MPS)
      else
         iret = nf_put_att_text(ncid, v_id, 'units',11, VAR_FPS)
      endif
      CALL check_err(iret)

!    --------------------------------------------------
!     determine variable name according to coordinates
!    ---------------------------------------------------
      if (ics.eq.2) then
          iret = nf_put_att_text(ncid, v_id, 'positive', 5, VAR_NORTH)
      else
          iret = nf_put_att_text(ncid, v_id, 'positive', 2, VAR_UP)
      endif
      CALL check_err(iret)

      iret = nf_put_att_DOUBLE(ncid, v_id, 'missing_values', 
     &       NF_DOUBLE, 1, doubleval)
      CALL check_err(iret)
!      doubleval(1) = -99999.
      iret = nf_put_att_DOUBLE(ncid, v_id, '_FillValue', NF_DOUBLE, 1, 
     &                       doubleval)
      CALL check_err(iret)
!      doubleval(1) = -99999.
      iret = nf_put_att_DOUBLE(ncid, v_id, 'dry_Value', NF_DOUBLE, 1, 
     &                       doubleval)
      CALL check_err(iret)

!    --------------------------------------------
!     Define sea surface elevation attributes
!    --------------------------------------------
      
      iret = nf_put_att_text(ncid, zeta1_id,'long_name',44, 
     &'surface elevation at the previous time step')
                        
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, zeta1_id, 'standard_name', 44, 
     &         'surface_elevation_at_previous_time step')
      CALL check_err(iret)
      
      if (G.LT.11.D0) then
          iret = nf_put_att_text(ncid, zeta1_id, 'units', 6, VAR_METS)
          CALL check_err(iret)
      else
         iret = nf_put_att_text(ncid, zeta1_id, 'units', 4, VAR_FEET)
         CALL check_err(iret)
      endif      
      
      iret = nf_put_att_DOUBLE(ncid,zeta1_id,'missing_values',NF_DOUBLE,
     &1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_DOUBLE(ncid, zeta1_id, '_FillValue', NF_DOUBLE, 
     &                      1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, zeta1_id,'positive', 2, VAR_UP)
      CALL check_err(iret)
!    --------------------------------------------
!     Define sea surface elevation attributes
!    --------------------------------------------
      
      iret = nf_put_att_text(ncid, zeta2_id,'long_name',43, 
     & 'surface elevation at the current time step')
                       
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, zeta2_id, 'standard_name', 43, 
     &         'surface_elevation_at_the_current_time_step')
      CALL check_err(iret)
      
      if (G.LT.11.D0) then
          iret = nf_put_att_text(ncid, zeta2_id, 'units', 6, VAR_METS)
          CALL check_err(iret)
      else
         iret = nf_put_att_text(ncid, zeta2_id, 'units', 4, VAR_FEET)
         CALL check_err(iret)
      endif      
      
      iret = nf_put_att_DOUBLE(ncid,zeta2_id,'missing_values',NF_DOUBLE,
     &1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_DOUBLE(ncid, zeta2_id, '_FillValue', NF_DOUBLE, 
     &                      1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, zeta2_id,'positive', 2, VAR_UP)
      CALL check_err(iret)
!    --------------------------------------------
!     Define sea surface elevation attributes
!    --------------------------------------------

      iret = nf_put_att_text(ncid, zetad_id,'long_name',37,       
     &    'elevation at flux specified boundary')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, zetad_id, 'standard_name', 37, 
     &                      'elevation_at_flux_specified_boundary')
      CALL check_err(iret)
      
      if (G.LT.11.D0) then
          iret = nf_put_att_text(ncid, zetad_id, 'units', 6, VAR_METS)
          CALL check_err(iret)
      else
         iret = nf_put_att_text(ncid, zetad_id, 'units', 4, VAR_FEET)
         CALL check_err(iret)
      endif      
      
      iret = nf_put_att_DOUBLE(ncid,zetad_id,'missing_values',NF_DOUBLE,
     &1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_DOUBLE(ncid, zetad_id, '_FillValue', NF_DOUBLE, 
     &                      1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, zetad_id,'positive', 2, VAR_UP)
      CALL check_err(iret)
!     -----------------
!     Global attributes
!     -----------------
         WRITE(6,*) "Global attributes: "

      iret = nf_put_att_text(ncid, NF_GLOBAL, 'model', 6, 'ADCIRC')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'version', 
     &       LEN(TRIM(ADC_VERSION)), ADC_VERSION)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'grid_type', 10, 
     &                      'Triangular')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'description', 
     &       LEN(TRIM(rundes(1:INDEX(rundes,"!")-1))), rundes)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'grid',
     &       LEN(TRIM(agrid)), agrid)       
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'title',
     &       LEN(TRIM(title(1:INDEX(title,"!")-1))), 
     &       title)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'institution', 
     &       LEN(TRIM(institution(1:INDEX(institution,"!")-1))), 
     &       institution)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'source', 
     &       LEN(TRIM(source(1:INDEX(source,"!")-1))),
     &       source)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'history', 
     &       LEN(TRIM(history(1:INDEX(history,"!")-1))),
     &       history)
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
     &       LEN(TRIM(host(1:INDEX(host,"!")-1))), 
     &       host)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'convention', 
     &       LEN(TRIM(convention(1:INDEX(convention,"!")-1))), 
     &       convention)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'contact', 
     &       LEN(TRIM(contact(1:INDEX(contact,"!")-1))), 
     &       contact)
      CALL check_err(iret)
 
 
      CALL date_and_time(now_date,big_ben,zone,values)
      write(date_string,71) values(1),values(2),values(3)
     & ,values(5),values(6),values(7),   (values(4))/60
 71   format(I4,'-',I2.2,'-',i2.2,' ',i2,':',i2.2,':',i2.2,' '
     &  ,i3.2,':00')
      iret = nf_put_att_text(ncid, NF_GLOBAL,'creation_date'
     & ,len_trim(date_string),date_string)
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
      iret = nf_put_att_int(ncid, NCGLOBAL, 'irtype', NF_INT, 1, irtype)
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
      iret = nf_put_att_double(ncid, NCGLOBAL, 'ffactor', NF_DOUBLE, 1, 
     &                         ffactor)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, NCGLOBAL, 'esl', NF_DOUBLE, 1, 
     &                         esl)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, NCGLOBAL, 'cori', NF_DOUBLE, 1, 
     &                         cori)
      CALL check_err(iret)
      iret = nf_put_att_int(ncid, NCGLOBAL, 'ntif', NF_INT, 1, ntif)
      CALL check_err(iret)
      iret = nf_put_att_int(ncid, NCGLOBAL, 'nbfr', NF_INT, 1, nbfr)
      CALL check_err(iret)
     
!     ======================
!     Leave define mode
!     ======================
      iret = nf_enddef(ncid)
      CALL check_err(iret)
!     ======================
!     Store static variables
!     ======================
!     -------
!     Store X
!     -------
      iret = nf_put_var_double(ncid, X_id, X)
      CALL check_err(iret)
!     -------
!     Store Y
!     -------
      iret = nf_put_var_double(ncid, Y_id, Y)
      CALL check_err(iret)
!     ------------
!     Store depth
!     ------------
      iret = nf_put_var_double(ncid, DEPTH_id, depth)
      CALL check_err(iret)
!     --------------
!     Store elements
!     --------------
      count(1)=nface_len
      count(2)=elem_len
      start(1)=1
      start(2)=1
      iret=nf_put_vara_int(ncid,ele_id,start,count,element)
      CALL check_err(iret)
!     ------------------------------------
!     Store elevation boundary information
!     -------------------------------------
      iret = nf_put_var_int(ncid, nopenc_id, nopenc)
      CALL check_err(iret)
      iret = nf_put_var_int(ncid, nvdllnc_id, nvdllnc)
      CALL check_err(iret)
!     --------------------------------------
!     Store normal flow boundary information
!     --------------------------------------
      iret = nf_put_var_int(ncid, nbounc_id, nbounc)
      CALL check_err(iret)
      iret = nf_put_var_int(ncid, ibtypenc_id, ibtypenc)
      CALL check_err(iret)
      iret = nf_put_var_int(ncid, nvellnc_id, nvellnc)
      CALL check_err(iret)      
         WRITE(6,*) "Store maximum depth: "
!     --------------------
!     Store maximum depth
!     --------------------
      maxdepth = maxval (depth)
      iret = nf_put_var_double(ncid, maxdepth_id, maxdepth)
      CALL check_err(iret)
!     --------------------
!     Store minimum depth
!     --------------------
      mindepth = minval (depth)
      iret = nf_put_var_double(ncid, mindepth_id, mindepth)
      CALL check_err(iret)
      
!     =====================
!     Open old netcdf file
!     =====================
      ELSE
          iret = NF_OPEN(FILENAME, NF_WRITE, ncid)
          CALL check_err(iret)   
!     ------------------
!     Inquire variables
!     ------------------
          iret=nf_inq_unlimdim(ncid,recid)
          CALL check_err(iret)   
          iret=nf_inq_dimlen(ncid,recid,nrec)
          CALL check_err(iret)   
      ENDIF
!     ------------------
!     Update Global Attributes
!     ------------------
!      IF(IHOT.NE.0.AND.FIRST67.EQ.1) THEN
       IF((IHOT.NE.0.AND.FIRST67.EQ.1).OR.
     &    (IHOT.NE.0.AND.FIRST68.EQ.1)) THEN
         iret = NF_REDEF (ncid)
         CALL check_err(iret)
         iret = nf_put_att_int(ncid, NCGLOBAL, 'ihot', NF_INT, 1, ihot)
         CALL check_err(iret)
         iret = nf_put_att_int(ncid, NCGLOBAL, 'nramp', NF_INT, 1,nramp)
         CALL check_err(iret)
         iret = nf_put_att_double(ncid, NCGLOBAL, 'tau0', NF_DOUBLE, 1, 
     &                         tau0)
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
         iret = nf_put_att_double(ncid, NCGLOBAL, 'ffactor', NF_DOUBLE, 1, 
     &                         ffactor)
         CALL check_err(iret)
         iret = nf_put_att_double(ncid, NCGLOBAL, 'esl', NF_DOUBLE, 1, 
     &                         esl)
         CALL check_err(iret)
         CALL date_and_time(now_date,big_ben,zone,values)
         WRITE(date_string,71) values(1),values(2),values(3),
     &         values(5),values(6),values(7),   (values(4))/60
         iret = nf_put_att_text(ncid, NF_GLOBAL,'modification_date',
     &          len_trim(date_string),date_string)
         CALL check_err(iret)
         iret = NF_ENDDEF (ncid)
         CALL check_err(iret)
       ENDIF

      IF(ncmode.eq.0) go to 9999
      IF((FIRST67.EQ.1.and.
     & ((ncmode.ne.0.and.IHOT.NE.0).or.(ncmode.eq.2.and.IHOT.eq.0).or.
     &  (ncmode.eq.2.and.IHOT.ne.0.and.nout.eq.-3))).OR.
     & (FIRST68.EQ.1.and.
     & ((ncmode.ne.0.and.IHOT.NE.0).or.(ncmode.eq.2.and.IHOT.eq.0).or.
     &  (ncmode.eq.2.and.IHOT.ne.0.and.nout.eq.-3)))
     &  )THEN
              WRITE(6,*) "FIRST67.EQ.1. ncmode.ne.0"
              WRITE(6,*) "FIRST68.EQ.1. ncmode.ne.0"

         iret=nf_inq_varid(ncid,"time",timenc_id)
         iret=nf_inq_varid(ncid,"zeta1",zeta1_id)
         iret=nf_inq_varid(ncid,"zeta2",zeta2_id)
         iret=nf_inq_varid(ncid,"zetad",zetad_id)
         iret=nf_inq_varid(ncid,"u-vel",u_id)
         iret=nf_inq_varid(ncid,"v-vel",v_id)
         iret=nf_inq_varid(ncid,"nodecode",nodecodenc_id)
         iret=nf_inq_varid(ncid,"noff",noffnc_id)
         iret=nf_inq_varid(ncid,"imhs",imhs_id)
         iret=nf_inq_varid(ncid,"iths",iths_id)
         iret=nf_inq_varid(ncid,"iestp",iestp_id)
         iret=nf_inq_varid(ncid,"nscoue",nscoue_id)
         iret=nf_inq_varid(ncid,"ivstp",ivstp_id)
         iret=nf_inq_varid(ncid,"nscouv",nscouv_id)
         iret=nf_inq_varid(ncid,"icstp",icstp_id)
         iret=nf_inq_varid(ncid,"nscouc",nscouc_id)
         iret=nf_inq_varid(ncid,"ipstp",ipstp_id)
         iret=nf_inq_varid(ncid,"iwstp",iwstp_id)
         iret=nf_inq_varid(ncid,"nscoum",nscoum_id)
         iret=nf_inq_varid(ncid,"igep",igep_id)
         iret=nf_inq_varid(ncid,"nscouge",nscouge_id)
         iret=nf_inq_varid(ncid,"igvp",igvp_id)
         iret=nf_inq_varid(ncid,"nscougv",nscougv_id)
         iret=nf_inq_varid(ncid,"igcp",igcp_id)
         iret=nf_inq_varid(ncid,"nscougc",nscougc_id)
         iret=nf_inq_varid(ncid,"igpp",igpp_id)
         iret=nf_inq_varid(ncid,"igwp",igwp_id)
         iret=nf_inq_varid(ncid,"nscougw",nscougw_id)
      ENDIF
!      IF(IHOT.NE.0.AND.FIRST67.EQ.1.and.ncmode.ne.0.and.nout.eq.3)THEN
!         iret = nf_get_var_double(ncid, minun_id, minun)
!         CALL check_err(iret)
!      ENDIF
      
      IF(FILEUNIT.EQ.67)FIRST67=0
      IF(FILEUNIT.EQ.68)FIRST68=0
!      NREC = NREC + 1
      NREC = 1
!     ----------
!     Store time
!     ----------
      starti(1)=nrec
      counti(1)=timenc_len
      iret = nf_put_vara_DOUBLE(ncid, timenc_id, starti,counti, timenc)
      CALL check_err(iret)
!     ---------------------------
!     Store zeta1
!     ---------------------------
              WRITE(6,*) "Store zeta1"
      count(1)=node_len
      count(2)=timenc_len
      start(1)=1
      start(2)=nrec
      iret = nf_put_vara_DOUBLE(ncid, zeta1_id, start, count, zeta1)
      CALL check_err(iret)
!     ---------------------------
!     Store zeta2
!     ---------------------------
              WRITE(6,*) "Store zeta2"
      count(1)=node_len
      count(2)=timenc_len
      start(1)=1
      start(2)=nrec
      iret = nf_put_vara_DOUBLE(ncid, zeta2_id, start, count, zeta2)
      CALL check_err(iret)
!     ---------------------------
!     Store zetad
!     ---------------------------
              WRITE(6,*) "Store zetad"
      count(1)=node_len
      count(2)=timenc_len
      start(1)=1
      start(2)=nrec
      iret = nf_put_vara_DOUBLE(ncid, zetad_id, start, count, zetad)
      CALL check_err(iret)

!     ---------------------------
!     Store u
!     ---------------------------
              WRITE(6,*) "Store u"
      count(1)=node_len
      count(2)=timenc_len
      start(1)=1
      start(2)=nrec
      iret = nf_put_vara_DOUBLE(ncid, u_id, start, count, u)
      CALL check_err(iret)
!     ---------------------------
!     Store v
!     ---------------------------
      count(1)=node_len
      count(2)=timenc_len
      start(1)=1
      start(2)=nrec
      iret = nf_put_vara_DOUBLE(ncid, v_id, start, count, v)
      CALL check_err(iret)
!     ---------------
!     Store nodecode
!     ---------------
      WRITE(6,*) "Store nodecodenc"
      iret = nf_put_var_int(ncid, nodecodenc_id, nodecodenc)
      CALL check_err(iret)
!     ---------------
!     Store noff
!     ---------------
      WRITE(6,*) "Store noffnc"
      iret = nf_put_var_int(ncid, noffnc_id, noffnc)
      CALL check_err(iret)
!     ------------------------
!     Store hotstart variables
!     ------------------------
               WRITE(6,*) "Store hotstart variables"
      iret = nf_put_var_int(ncid, imhs_id, imhs)
      CALL check_err(iret)
      iret = nf_put_var_int(ncid, iths_id, iths)
      CALL check_err(iret)
      iret = nf_put_var_int(ncid, iestp_id, iestp)
      CALL check_err(iret)
      iret = nf_put_var_int(ncid, nscoue_id, nscoue)
      CALL check_err(iret)
      iret = nf_put_var_int(ncid, ivstp_id, ivstp)
      CALL check_err(iret)
      iret = nf_put_var_int(ncid, nscouv_id, nscouv)
      CALL check_err(iret)
      iret = nf_put_var_int(ncid, icstp_id, icstp)
      CALL check_err(iret)
      iret = nf_put_var_int(ncid, nscouc_id, nscouc)
      CALL check_err(iret)
      iret = nf_put_var_int(ncid, ipstp_id, ipstp)
      CALL check_err(iret)
      iret = nf_put_var_int(ncid, iwstp_id, iwstp)
      CALL check_err(iret)
      iret = nf_put_var_int(ncid, nscoum_id, nscoum)
      CALL check_err(iret)
               WRITE(6,*) "Store igep"
      iret = nf_put_var_int(ncid, igep_id, igep)
      CALL check_err(iret)
      iret = nf_put_var_int(ncid, nscouge_id, nscouge)
      CALL check_err(iret)
      iret = nf_put_var_int(ncid, igvp_id, igvp)
      CALL check_err(iret)
      iret = nf_put_var_int(ncid, nscougv_id, nscougv)
      CALL check_err(iret)
      iret = nf_put_var_int(ncid, igcp_id, igcp)
      CALL check_err(iret)
      iret = nf_put_var_int(ncid, nscougc_id, nscougc)
      CALL check_err(iret)
      iret = nf_put_var_int(ncid, igpp_id, igpp)
      CALL check_err(iret)
      iret = nf_put_var_int(ncid, igwp_id, igwp)
      CALL check_err(iret)
      iret = nf_put_var_int(ncid, nscougw_id, nscougw)
      CALL check_err(iret)
!     ---------------------------
!     Closing fort.X4 netCDF file
!     ---------------------------
!      IF(NREC.EQ.NDATA) THEN
!         iret = nf_put_var_double(ncid, minun_id, minun)
!         CALL check_err(iret)
!         iret = nf_put_var_double(ncid, maxun_id, maxun)
!         CALL check_err(iret)
!         iret = nf_put_var_double(ncid, minvn_id, minvn)
!         CALL check_err(iret)
!         iret = nf_put_var_double(ncid, maxvn_id, maxvn)
!         CALL check_err(iret)
!      ENDIF

!     ---------------------------
!     Closing netCDF file
!     ---------------------------
      WRITE(6,*) "Closing ", FILENAME, " netCDF file ... "
      WRITE(6,*) "-----------------------------------------------------"
      
9999  iret = nf_close(ncid)
      CALL check_err(iret)
      WRITE(6,*)'end subroutine writenc67'

      end subroutine writenc67     

!=======================================================================
!     Subroutine readnc67
!     Writes a netCDF file equivalent to ADCIRC's fort.64 ascii file
!
!     Dimensions:
!     node  : number of nodes in the horizontal grid
!     nele  : number of elements in the horizontal grid
!     nope  : number of elevation specified boundary forcing segments
!     neta  : total number of elevation specified boundary nodes
!     nbou  : number of normal flow specified boundary segments
!     nvel  : total number of normal flow specified boundary nodes
!     time  : number of time records
!
!     Variables :
!     time      : time
!     u-vel     : u velocity in x or east/west direction 
!     v-vel     : v velocity in y or north/south direction
!     x         : x-coordinate or longitude
!     y         : y-coordinate or latitude
!     element   : elements in grid
!     nvdll     : number of nodes in elevation boundary segment
!     nbdv      : node numbers on elevation specified boundary segment
!     nvell     : number of nodes in normal flow specified boundary segment
!     nbvv      : node numbers on normal flow boundary segment
!     ibtype    : boundary type
!     depth     : depth in each node
!     mindepth  : minimum depth in grid
!     maxdepth  : maximum depth in grid
!     minuvelt  : minimum vertically averaged u at each time
!     maxuvelt  : maximum vertically averaged u at each time
!     minvvelt  : minimum vertically averaged v at each time
!     maxvvelt  : maximum vertically averaged v at each time
!     minuveln  : minimum vertically averaged u at each node
!     maxuveln  : maximum vertically averaged u at each node
!     minvveln  : minimum vertically averaged v at each node
!     maxvveln  : maximum vertically averaged v at each node
!     minunodet : node where minimum u occurs at each output time
!     maxunodet : node where maximum u occurs at each output time
!     minvnodet : node where minimum v occurs at each output time
!     maxvnodet : node where maximum v occurs at each output time
!     minutimen : time when minimum u occurs at each node
!     maxutimen : time when maximum u occurs at each node
!     minvtimen : time when minimum v occurs at each node
!     maxvtimen : time when maximum v occurs at each node
!
!     Attributes:
!     model             : ADCIRC
!     version           : model version number
!     grid_type         : triangular
!     rundes            : description of simulation
!     agrid             : grid ID 
!     title             : description of the dataset
!     institution       : where the original data was produced
!     source            : method of production of the original data
!     history           : provides an audit trail for modifications to the original data
!     references        : published or web-based references that describe the data or methods used to produce it
!     comments          : miscellaneous information about the data
!     contact           : user's institution/ e-mail address
!     creation_date     : date the data was created
!     modification_date : date the data was modified
!     Parameters from fort.15 input file 
!
!     Revision history:
!
!     Date            Programmer                  Description of change
!     ----            ----------                  ----------------------
!     03/30/07  Cristina Forbes, PSGS @ UNC-IMS   Wrote original code
!     03/30/08  Cristina Forbes, PSGS @ UNC-IMS   Modified code for 
!                                                 globalio & hotstart 
!                                                 from binary files
!=======================================================================
      subroutine readnc67   (nout,fileunit,timenc_len,node_len,timenc,
     &                      nt, NDATA, u,v,
     &                      zeta1, zeta2, zetad, nodecodenc,noffnc,
     &                      imhs, iths,iestp,nscoue,ivstp,nscouv, icstp, 
     &                      nscouc, ipstp, iwstp, nscoum, igep, nscouge,
     &                      igvp, nscougv, igcp, nscougc, igpp, igwp,
     &                      nscougw,
     &                      rundes, runid, agrid, irtype,
     &                      dtdp, ihot, ics, nolibf, nolifa,
     &                      nolica, nolicat, nwp, 
     &                      ncor, ntip, nws, nramp, tau0,
     &                      statim, reftim, rnday, dramp,
     &                      a00, b00, c00, h0,
     &                      slam0, sfea0, ffactor,esl, cori,
     &                      ntif, nbfr, elem_len, nface_len, nm,
     &                      x,y,depth,
     &                      nope_len, neta_len, nbdv,nopenc,netanc_len,
     &                      nvdllnc, 
     &                      nbou_len, nvel_len, nbvv,nbounc,nvelnc_len,
     &                      nvellnc,ibtypenc,
     &                      ncmode,
     &                      title, institution, source, history,
     &                      references, comments, host, convention,
     &                      contact,
     &                      base_date)


!--------------------------------------------
!     Force explicit declaration of all variables
!--------------------------------------------
      IMPLICIT NONE
!     -------------------
!     Error status return
!     -------------------
      INTEGER  iret
!     ---------
!     netCDF id
!     ---------
      INTEGER  ncid
!     -------------
!     Dimension ids
!     -------------
      INTEGER  node_dim
      INTEGER  timenc_dim
      INTEGER  elem_dim
      INTEGER  nface_dim
      INTEGER  nope_dim
      INTEGER  neta_dim
      INTEGER  nbou_dim
      INTEGER  nvel_dim
      INTEGER  nvelnc_dim
      INTEGER  netanc_dim
!     -----------------
!     Dimension lengths
!     -----------------
      INTEGER  node_len
      INTEGER  timenc_len
      INTEGER  elem_len
      INTEGER  nface_len
      INTEGER  nope_len
      INTEGER  neta_len
      INTEGER  nbou_len
      INTEGER  nvel_len
      INTEGER  nvelnc_len
      INTEGER  netanc_len
!     ------------
!     Variable ids
!     ------------
      INTEGER  timenc_id
      INTEGER  X_id
      INTEGER  Y_id
      INTEGER  DEPTH_id
      INTEGER  ELE_id
      INTEGER  mindepth_id
      INTEGER  maxdepth_id
      INTEGER  nbdv_id
      INTEGER  nbvv_id
      
      INTEGER  u_id
      INTEGER  v_id
      INTEGER  zeta1_id
      INTEGER  zeta2_id
      INTEGER  zetad_id
      INTEGER  nodecodenc_id
      INTEGER  noffnc_id

      INTEGER  nvdll_id
      INTEGER  nbounc_id
      INTEGER  ibtypenc_id
      INTEGER  nvellnc_id
      INTEGER  nvdllnc_id
      INTEGER  nopenc_id
      
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
!     ---------------------------------------------
!     Rank (number of dimensions) for each variable
!     ---------------------------------------------
      INTEGER  timenc_rank
      INTEGER  depth_rank
      INTEGER  ele_rank
      INTEGER  mindepth_rank
      INTEGER  maxdepth_rank
      INTEGER  x_rank
      INTEGER  y_rank
      INTEGER  nvdll_rank
      INTEGER  nbdv_rank
      INTEGER  nbvv_rank
      INTEGER  u_rank
      INTEGER  v_rank
      INTEGER  zeta1_rank
      INTEGER  zeta2_rank
      INTEGER  zetad_rank
      INTEGER  nodecodenc_rank
      INTEGER  noffnc_rank

      INTEGER  nbounc_rank
      INTEGER  ibtypenc_rank
      INTEGER  nvellnc_rank
      INTEGER  nvdllnc_rank
      INTEGER  nopenc_rank

      INTEGER imhs_rank
      INTEGER iths_rank
      INTEGER IESTP_rank
      INTEGER NSCOUE_rank
      INTEGER IVSTP_rank
      INTEGER NSCOUV_rank
      INTEGER ICSTP_rank
      INTEGER NSCOUC_rank
      INTEGER IPSTP_rank
      INTEGER IWSTP_rank
      INTEGER NSCOUM_rank
      INTEGER IGEP_rank
      INTEGER NSCOUGE_rank
      INTEGER IGVP_rank
      INTEGER NSCOUGV_rank
      INTEGER IGCP_rank
      INTEGER NSCOUGC_rank
      INTEGER IGPP_rank
      INTEGER IGWP_rank
      INTEGER NSCOUGW_rank

      PARAMETER (timenc_rank = 1)
      PARAMETER (depth_rank = 1)
      PARAMETER (ele_rank = 2)
      PARAMETER (mindepth_rank = 0)
      PARAMETER (maxdepth_rank = 0)
      PARAMETER (x_rank = 1)
      PARAMETER (y_rank = 1)
      PARAMETER (nbdv_rank = 2)
      PARAMETER (nbvv_rank = 2)
      PARAMETER (u_rank = 2)
      PARAMETER (v_rank = 2)
      PARAMETER (zeta1_rank = 2)
      PARAMETER (zeta2_rank = 2)
      PARAMETER (zetad_rank = 2)
      PARAMETER (nodecodenc_rank = 1)
      PARAMETER (noffnc_rank = 1)

      PARAMETER (nvdll_rank = 1)
      PARAMETER (nbounc_rank = 1)
      PARAMETER (ibtypenc_rank = 1)
      PARAMETER (nvellnc_rank = 1)
      PARAMETER (nvdllnc_rank=1)
      PARAMETER (nopenc_rank=1)

      PARAMETER (imhs_rank = 0)
      PARAMETER (iths_rank = 0)
      PARAMETER (IESTP_rank = 0)
      PARAMETER (NSCOUE_rank = 0)
      PARAMETER (IVSTP_rank = 0)
      PARAMETER (NSCOUV_rank = 0)
      PARAMETER (ICSTP_rank = 0)
      PARAMETER (NSCOUC_rank = 0)
      PARAMETER (IPSTP_rank = 0)
      PARAMETER (IWSTP_rank = 0)
      PARAMETER (NSCOUM_rank = 0)
      PARAMETER (IGEP_rank = 0)
      PARAMETER (NSCOUGE_rank = 0)
      PARAMETER (IGVP_rank = 0)
      PARAMETER (NSCOUGV_rank = 0)
      PARAMETER (IGCP_rank = 0)
      PARAMETER (NSCOUGC_rank = 0)
      PARAMETER (IGPP_rank = 0)
      PARAMETER (IGWP_rank = 0)
      PARAMETER (NSCOUGW_rank = 0)
      
!     ---------------
!     Variable shapes
!     ---------------
      INTEGER  timenc_dims(timenc_rank)
      INTEGER  x_dims(x_rank)
      INTEGER  y_dims(y_rank)
      INTEGER  depth_dims(depth_rank)
      INTEGER  ele_dims(ele_rank)
      INTEGER  nbdv_dims(nbdv_rank)
      INTEGER  nbvv_dims(nbvv_rank)
      INTEGER  u_dims(u_rank)
      INTEGER  v_dims(v_rank)
      INTEGER  zeta1_dims(zeta1_rank)
      INTEGER  zeta2_dims(zeta2_rank)
      INTEGER  zetad_dims(zetad_rank)
      INTEGER  nodecodenc_dims(nodecodenc_rank)
      INTEGER  noffnc_dims(noffnc_rank)

      INTEGER  nvdll_dims(nvdll_rank)
      INTEGER  nbounc_dims(nbounc_rank)
      INTEGER  ibtypenc_dims(ibtypenc_rank)
      INTEGER  nvellnc_dims(nvellnc_rank)
      INTEGER  nvdllnc_dims(nvdllnc_rank)
      INTEGER  nopenc_dims(nopenc_rank)

!     --------------
!     Data variables
!     --------------
      CHARACTER*5 FILEFORT
      CHARACTER*10 FILENAME
      CHARACTER*3 FILEEXT
      INTEGER   t,ntif,nbfr,nt
      INTEGER   irtype, ihot, ics, nolibf, nolifa, nolica, nolicat, nwp
      INTEGER   ntip, nws, nramp,ncor
      INTEGER   nout
      INTEGER   element (nface_len , elem_len)
      INTEGER   nbdv (nope_len, neta_len)
      INTEGER   nbvv (nbou_len, nvel_len)
      
      INTEGER   nvellnc(nbou_len)
      INTEGER   nbounc (nvel_len)
      INTEGER   ibtype(nbou_len)
      INTEGER   ibtypenc(nbou_len)
      INTEGER   nvdllnc (nope_len)
      INTEGER   nopenc (neta_len)
      
      INTEGER   nm (elem_len, nface_len)
      INTEGER   fileunit
      REAL*8    timenc(timenc_len)
      REAL(sz)  x(node_len)
      REAL(sz)  y(node_len)
      REAL(sz)  depth(node_len)
      REAL(sz)  u(node_len, timenc_len)
      REAL(sz)  v(node_len, timenc_len)
      REAL(sz)  zeta1 (node_len, timenc_len)
      REAL(sz)  zeta2 (node_len, timenc_len)
      REAL(sz)  zetad (node_len, timenc_len)
      INTEGER   nodecodenc(node_len)
      INTEGER   noffnc(elem_len)

      real*8  mindepth
      real*8  maxdepth
      REAL*8 dtdp, statim,reftim, rnday, dramp, tau0
      REAL*8 a00, b00, c00,h0
      REAL(sz) slam0, sfea0, ffactor,esl, cori
      REAL(sz) slam0deg, sfea0deg
      CHARACTER(80) rundes, runid, agrid
      INTEGER, SAVE  :: FIRST67 = 1 
      INTEGER, SAVE  :: FIRST68 = 1 
      INTEGER NDATA
      INTEGER, SAVE  :: NREC = 0
      INTEGER ncmode
      
      INTEGER imhs
      INTEGER iths
      INTEGER IESTP
      INTEGER NSCOUE
      INTEGER IVSTP
      INTEGER NSCOUV
      INTEGER ICSTP
      INTEGER NSCOUC
      INTEGER IPSTP
      INTEGER IWSTP
      INTEGER NSCOUM
      INTEGER IGEP
      INTEGER NSCOUGE
      INTEGER IGVP
      INTEGER NSCOUGV
      INTEGER IGCP
      INTEGER NSCOUGC
      INTEGER IGPP
      INTEGER IGWP
      INTEGER NSCOUGW

      

!     -----------------
!     Attribute vectors
!     -----------------
      INTEGER intval(1),i,j
      INTEGER counti(1), starti(1), n
      INTEGER count(2), start(2)
      REAL    realval(1)
!     -----------------
!     Attribute general
!     -----------------
      CHARACTER(*) :: title, institution, source, history
      CHARACTER(*) :: references, comments, host, convention, contact
      CHARACTER(*) :: base_date

!     -----------------
!      date_string variables for time attribute
!     -----------------
      character date_string*40
      character now_date*8
      character big_ben*10
      character zone*5
      integer values(8)
!     ---------------------------------
!      inquiring and reading variables
!     ---------------------------------

      INTEGER NDIMS, NVARS, NGATTS, UNLIMDIMID, RECID,NRECS
      CHARACTER*4 recname
      
!     -------------------------
!     Save dims, ids, variables
!     -------------------------
      save node_dim,  elem_dim, nface_dim
      save nope_dim, neta_dim, nbou_dim, nvel_dim
      save u_id, v_id, zeta1_id, zeta2_id, zetad_id
      save nodecodenc_id, noffnc_id
      save X_id, Y_id

      save DEPTH_id,ELE_id, mindepth_id, maxdepth_id, nbdv_id, nbvv_id
      save timenc_dim, timenc_id
      save mindepth, maxdepth

      save nvelnc_dim
      save nvdll_id, nbounc_id, ibtypenc_id, nvellnc_id
      save netanc_dim
      save nopenc_id,nvdllnc_id


      save imhs_id
      save iths_id
      save IESTP_id
      save NSCOUE_id
      save IVSTP_id
      save NSCOUV_id
      save ICSTP_id
      save NSCOUC_id
      save IPSTP_id
      save IWSTP_id
      save NSCOUM_id
      save IGEP_id
      save NSCOUGE_id
      save IGVP_id
      save NSCOUGV_id
      save IGCP_id
      save NSCOUGC_id
      save IGPP_id
      save IGWP_id
      save NSCOUGW_id


!     ---------------------------
!     Opening fort.67 netCDF file
!     ---------------------------
      WRITE(6,*) "-----------------------------------------------------"
           PRINT *, " *********  node_len = ", node_len
           PRINT *, " *********  elem_len = ", elem_len
       
!      do i=1, node_len
!            PRINT *, "i ",i, " nodecodenc(i) = ", nodecodenc(i)
!      end do
!      do i=1, elem_len
!            PRINT *, "i ",i, " noffnc(i) = ", noffnc(i)
!      end do

!      IF(FIRST67.EQ.1.and.ncmode.NE.0) THEN
         
!      ENDIF
       
!     --------------------------------
!     Switch order in array for NETCDF
!     --------------------------------
       if(ncmode.ne.0) then

       do i=1, elem_len
         do j=1, nface_len
            element(j,i) = nm (i,j)
         end do
       end do
       endif 
      
!     -----------------
!     Enter define mode
!     -----------------
         FILEFORT = 'fort.'
         FILEEXT = '.nc'
         WRITE(FILENAME,'(A5,I2,A3)')FILEFORT,FILEUNIT,FILEEXT
         WRITE(6,*) "Opening ",FILENAME,  " netCDF file ... "

         PRINT *,"netcdfio: ncmode = ",ncmode," IHOT= ",IHOT
         
!       IF((FIRST67.EQ.1.AND.IHOT.EQ.0.AND.NCMODE.NE.2).OR.
!     &    (FIRST67.EQ.1.AND.IHOT.NE.0.AND.NCMODE.NE.2.AND.
!     &     nout.EQ.-3))THEN

       IF( ((FIRST67.EQ.1.AND.IHOT.EQ.0.AND.NCMODE.NE.2).OR.
     &     (FIRST68.EQ.1.AND.IHOT.EQ.0.AND.NCMODE.NE.2)).OR.
     &    ((FIRST67.EQ.1.AND.IHOT.NE.0.AND.NCMODE.NE.2.AND.
     &     nout.EQ.-3).OR.
     &     (FIRST68.EQ.1.AND.IHOT.NE.0.AND.NCMODE.NE.2.AND.
     &     nout.EQ.-3)))THEN

         PRINT *,"netcdfio: inside first if"


           iret = nf_create(FILENAME, NF_CLOBBER, ncid)
           CALL check_err(iret)
           IF(FILEUNIT.EQ.67)FIRST67=0
           IF(FILEUNIT.EQ.68)FIRST68=0
           WRITE(6,*) "FIRST67 ",FIRST67
           WRITE(6,*) "FIRST68 ",FIRST68
         
!     --------------------------------
!     Convert to degress for NETCDF I/O
!     --------------------------------
           SLAM0DEG=SLAM0
           SFEA0DEG=SFEA0
           IF (NCMODE.EQ.1) THEN
               SLAM0DEG=SLAM0*RAD2DEG
               SFEA0DEG=SFEA0*RAD2DEG
           ENDIF
!     -----------------
!     Define dimensions
!     -----------------
         WRITE(6,*) "Define dimensions "
      iret = nf_def_dim(ncid, 'node', node_len, node_dim)
      CALL check_err(iret)
      iret = nf_def_dim(ncid, 'nele', elem_len, elem_dim)
      CALL check_err(iret)
      iret = nf_def_dim(ncid, 'nface', 3, nface_dim)
      CALL check_err(iret)
      if(nope_len.ne.0) then
         iret = nf_def_dim(ncid, 'nope', nope_len, nope_dim)
         call check_err(iret)
      endif
      if(neta_len.ne.0) then
         iret = nf_def_dim(ncid, 'neta', neta_len, neta_dim)
         call check_err(iret)
      endif
      if(nbou_len.ne.0) then
         iret = nf_def_dim(ncid, 'nbou', nbou_len, nbou_dim)
         call check_err(iret)
      endif
      if(nvel_len.ne.0) then
         iret = nf_def_dim(ncid, 'nvel', nvel_len, nvel_dim)
         call check_err(iret)
      endif
      
      iret = nf_def_dim(ncid, 'time', nf_unlimited, timenc_dim)
      CALL check_err(iret)
      
!     ================
!     Define variables
!     ================
         WRITE(6,*) "Define variables "
!     ----------------
!     Define time
!     ----------------
      timenc_dims(1) = timenc_dim
      iret = nf_def_var(ncid,'time',NF_DOUBLE, timenc_rank, timenc_dims, 
     &                  timenc_id)
      CALL check_err(iret)
      
!     -------------------
!     Define x
!     -------------------
      x_dims(1) = node_dim
      iret = nf_def_var(ncid, 'x', NF_DOUBLE, x_rank, x_dims, x_id)
      CALL check_err(iret)
!     -------------------
!     Define Y coordinate
!     --------------------
      Y_dims(1) = node_dim
      iret = nf_def_var(ncid, 'y', NF_DOUBLE, Y_rank, Y_dims, Y_id)
      CALL check_err(iret)
!     -------------------
!     Define elements
!     --------------------
      ELE_dims(1) = nface_dim
      ELE_dims(2) = elem_dim
      iret = nf_def_var(ncid, 'element',NF_INT,ELE_rank,ELE_dims,ELE_id)
      CALL check_err(iret)
!     -----------------------------------------------------------------
!     Define elevation specified boundary forcing segments information
!     -----------------------------------------------------------------
      if(nope_len.ne.0.and.neta_len.ne.0) then
        nvdllnc_dims(1) = nope_dim      
        iret = nf_def_var(ncid, 'nvdll',NF_INT, nvdllnc_rank, 
     &                  nvdllnc_dims, nvdllnc_id)
        CALL check_err(iret)

        nopenc_dims(1) = neta_dim      
        iret = nf_def_var(ncid, 'nbbv',NF_INT, nopenc_rank, nopenc_dims, 
     &                  nopenc_id)
        CALL check_err(iret)

      endif

!     ----------------------------------------------------------
!     Define normal flow specified boundary segments information
!     ----------------------------------------------------------
      if(nbou_len.ne.0.and.nvel_len.ne.0) then

        nvellnc_dims(1) = nbou_dim      
        iret = nf_def_var(ncid, 'nvell',NF_INT, nvellnc_rank, 
     &                  nvellnc_dims, nvellnc_id)
        CALL check_err(iret)
        
        ibtypenc_dims(1) = nbou_dim      
        iret = nf_def_var(ncid, 'ibtype',NF_INT, ibtypenc_rank, 
     &                  ibtypenc_dims, ibtypenc_id)
        CALL check_err(iret)
        nbounc_dims(1) = nvel_dim      
        iret = nf_def_var(ncid, 'nbvv',NF_INT, nbounc_rank, nbounc_dims, 
     &                  nbounc_id)
        CALL check_err(iret)
     
       endif

         WRITE(6,*) "Define depth "
!     -------------------
!     Define Z coordinate
!     --------------------
      DEPTH_dims(1) = node_dim
      iret = nf_def_var(ncid, 'depth', NF_DOUBLE,DEPTH_rank, DEPTH_dims, 
     &                  DEPTH_id)
      CALL check_err(iret)
!     --------------------
!     Define minimum depth
!     --------------------
      iret = nf_def_var(ncid, 'mindepth', NF_DOUBLE, mindepth_rank, 0, 
     &                  mindepth_id)
      CALL check_err(iret)
!     --------------------
!     Define maximum depth
!     --------------------
      iret = nf_def_var(ncid, 'maxdepth', NF_DOUBLE, maxdepth_rank, 0, 
     &                  maxdepth_id)
      CALL check_err(iret)

!     ----------------
!     Define zeta
!     ----------------
      zeta1_dims(1) = node_dim
      zeta1_dims(2) = timenc_dim
      
      iret = nf_def_var(ncid, 'zeta1',NF_DOUBLE, zeta1_rank, zeta1_dims, 
     &                  ZETA1_id)
      
      CALL check_err(iret)
!     ----------------
!     Define zeta
!     ----------------
      zeta2_dims(1) = node_dim
      zeta2_dims(2) = timenc_dim
      
      iret = nf_def_var(ncid, 'zeta2',NF_DOUBLE, zeta2_rank, zeta2_dims, 
     &                  ZETA2_id)
      
      CALL check_err(iret)

!     ----------------
!     Define zeta
!     ----------------
      zetad_dims(1) = node_dim
      zetad_dims(2) = timenc_dim
      
      iret = nf_def_var(ncid, 'zetad',NF_DOUBLE, zetad_rank, zetad_dims, 
     &                  ZETAd_id)
      
      CALL check_err(iret)


         WRITE(6,*) "Define vel "
!     ----------------
!     Define u vel
!     ----------------
      u_dims(1) = node_dim
      u_dims(2) = timenc_dim
      iret = nf_def_var(ncid, 'u-vel', NF_DOUBLE, u_rank, u_dims, 
     &                  u_id)
      CALL check_err(iret)

!     ----------------
!     Define v vel
!     ----------------
      v_dims(1) = node_dim
      v_dims(2) = timenc_dim
      iret = nf_def_var(ncid, 'v-vel', NF_DOUBLE, v_rank, v_dims, 
     &                  v_id)
      CALL check_err(iret)

!      WRITE(6,*)"Define MIN u/v at times"

         WRITE(6,*) "hotstart parameters "
!     --------------------
!     Define hotstart parameters
!     --------------------
      iret = nf_def_var(ncid, 'nodecode', NF_INT, nodecodenc_rank, 
     &                  node_dim, nodecodenc_id)
      CALL check_err(iret)
         WRITE(6,*) "noff: "
      iret = nf_def_var(ncid, 'noff', NF_INT, noffnc_rank, elem_dim, 
     &                  noffnc_id)
      CALL check_err(iret)
         WRITE(6,*) "imhs: "
      iret = nf_def_var(ncid, 'imhs', NF_INT, imhs_rank, 0, 
     &                  imhs_id)
      CALL check_err(iret)
      iret = nf_def_var(ncid, 'iths', NF_INT, iths_rank, 0, 
     &                  iths_id)
      CALL check_err(iret)
      iret = nf_def_var(ncid, 'iestp', NF_INT, iestp_rank, 0, 
     &                  iestp_id)
      CALL check_err(iret)
      iret = nf_def_var(ncid, 'nscoue', NF_INT, nscoue_rank, 0, 
     &                  nscoue_id)
      CALL check_err(iret)

         WRITE(6,*) "ivstp: "
      iret = nf_def_var(ncid, 'ivstp', NF_INT, ivstp_rank, 0, 
     &                  ivstp_id)
      CALL check_err(iret)
      iret = nf_def_var(ncid, 'nscouv', NF_INT, nscouv_rank, 0, 
     &                  nscouv_id)
      CALL check_err(iret)
      iret = nf_def_var(ncid, 'icstp', NF_INT, icstp_rank, 0, 
     &                  icstp_id)
      CALL check_err(iret)
      iret = nf_def_var(ncid, 'nscouc', NF_INT, nscouc_rank, 0, 
     &                  nscouc_id)
      CALL check_err(iret)

      iret = nf_def_var(ncid, 'ipstp', NF_INT, ipstp_rank, 0, 
     &                  ipstp_id)
      CALL check_err(iret)
      iret = nf_def_var(ncid, 'iwstp', NF_INT, iwstp_rank, 0, 
     &                  iwstp_id)
      CALL check_err(iret)
      iret = nf_def_var(ncid, 'nscoum', NF_INT, nscoum_rank, 0, 
     &                  nscoum_id)
      CALL check_err(iret)
      iret = nf_def_var(ncid, 'igep', NF_INT, igep_rank, 0, 
     &                  igep_id)
      CALL check_err(iret)
      
      iret = nf_def_var(ncid, 'nscouge', NF_INT, nscouge_rank, 0, 
     &                  nscouge_id)
      CALL check_err(iret)
      iret = nf_def_var(ncid, 'igvp', NF_INT, igvp_rank, 0, 
     &                  igvp_id)
      CALL check_err(iret)
      iret = nf_def_var(ncid, 'nscougv', NF_INT, nscougv_rank, 0, 
     &                  nscougv_id)
      CALL check_err(iret)
      iret = nf_def_var(ncid, 'igcp', NF_INT, igcp_rank, 0, 
     &                  igcp_id)
      CALL check_err(iret)
         WRITE(6,*) "nscougc: "

      iret = nf_def_var(ncid, 'nscougc', NF_INT, nscougc_rank, 0, 
     &                  nscougc_id)
      CALL check_err(iret)
      iret = nf_def_var(ncid, 'igpp', NF_INT, igpp_rank, 0, 
     &                  igpp_id)
      CALL check_err(iret)
      iret = nf_def_var(ncid, 'igwp', NF_INT, igwp_rank, 0, 
     &                  igwp_id)
      CALL check_err(iret)
      iret = nf_def_var(ncid, 'nscougw', NF_INT, nscougw_rank, 0, 
     &                  nscougw_id)
      CALL check_err(iret)
!     -----------------
!     Assign attributes
!     -----------------
         WRITE(6,*) "Assign attributes: "

      if(ics.eq.2) then
      iret = nf_put_att_text(ncid, X_id, 'long_name', 9, 'longitude')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, X_id, 'standard_name', 9, 
     &                      'longitude')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, X_id, 'units', 12, VAR_DEGE)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, X_id, 'positive', 4, VAR_EAST)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, X_id, 'missing_values', NF_DOUBLE,
     &                        1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, X_id, '_FillValue', NF_DOUBLE, 1, 
     &                       doubleval)
      CALL check_err(iret)      
      
      iret = nf_put_att_text(ncid, Y_id, 'long_name', 8, 'latitude')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, Y_id, 'standard_name', 8, 
     &                       'latitude')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, Y_id, 'units', 13, VAR_DEGN)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, Y_id, 'positive', 5, VAR_NORTH)
      CALL check_err(iret)
       
      iret = nf_put_att_double(ncid, Y_id, 'missing_values', NF_DOUBLE, 
     &       1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, Y_id, '_FillValue', NF_DOUBLE, 1, 
     &       doubleval)
      CALL check_err(iret)

      else
      
      iret = nf_put_att_text(ncid, X_id, 'long_name', 22, 
     &       'Cartesian coordinate x')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, X_id, 'standard_name', 12, 
     &                      'x_coordinate')
      CALL check_err(iret)
      
      if (G.LT.11.D0) then
         iret = nf_put_att_text(ncid, X_id, 'units', 6,  VAR_METS)
         CALL check_err(iret)
      else
         iret = nf_put_att_text(ncid, X_id, 'units', 4,  VAR_FEET)
         CALL check_err(iret)
      endif
      
      iret = nf_put_att_text(ncid, X_id, 'positive', 5,  VAR_RG)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, X_id, 'missing_values', NF_DOUBLE, 
     &                       1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, X_id, '_FillValue', NF_DOUBLE, 1, 
     &                       doubleval)
      CALL check_err(iret)
!     ------------------------------
!     Define y-coordinate attributes
!     ------------------------------     
      
      iret = nf_put_att_text(ncid, Y_id, 'long_name', 22, 
     &       'Cartesian coordinate y')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, Y_id, 'standard_name', 12, 
     &                       'y_coordinate')
      CALL check_err(iret)
      if (G.LT.11.D0) then
         iret = nf_put_att_text(ncid, Y_id, 'units', 6, VAR_METS)
         CALL check_err(iret)
      else
         iret = nf_put_att_text(ncid, Y_id, 'units', 4, VAR_FEET)
         CALL check_err(iret)
      endif
      iret = nf_put_att_text(ncid, Y_id, 'positive', 2, VAR_UP)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, Y_id, 'missing_values', 
     &       NF_DOUBLE, 1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, Y_id, '_FillValue', NF_DOUBLE, 1, 
     &       doubleval)
      CALL check_err(iret)
      
      endif

!     ------------------------------
!     Define bathymetry attributes
!     ------------------------------

      iret = nf_put_att_text(ncid, DEPTH_id,'long_name',10,'bathymetry')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, DEPTH_id, 'standard_name', 5,'depth')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, DEPTH_id, 'units', 6, VAR_METS)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, DEPTH_id, 'missing_values', 
     &                      NF_DOUBLE,  1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid,DEPTH_id, '_FillValue',NF_DOUBLE, 1, 
     &                       doubleval)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, DEPTH_id, 'positive', 4, VAR_DW)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, ELE_id, 'long_name', 7, 'element')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, ELE_id, 'units', 14, VAR_NOND)
      intval(1) = -1
      
!    ---------------------------
!     minimum/maximum bathymetry
!    ---------------------------

      iret = nf_put_att_text(ncid, mindepth_id, 'long_name', 18, 
     &                       'minimum bathymetry')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, mindepth_id, 'standard_name', 13,
     &                'minimum_depth')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, mindepth_id, 'units', 6, VAR_METS)
      CALL check_err(iret)
      
      iret = nf_put_att_text(ncid, maxdepth_id, 'long_name', 18, 
     &                       'maximum bathymetry')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, maxdepth_id, 'standard_name', 13,
     &           'maximum_depth')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, maxdepth_id, 'units', 6, VAR_METS)
      CALL check_err(iret)

!    ---------------------------
!     Hotstart variables
!    ---------------------------
         WRITE(6,*) "Assign nodecodenc_id: "


       iret = nf_put_att_text(ncid, nodecodenc_id, 'long_name', 82, 
     &                       'node code indicating whether the node is
     & presently wet (active) or dry (inactive)')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, nodecodenc_id, 'standard_name', 18,
     &                'wet_dry_indicator')
      CALL check_err(iret)
         WRITE(6,*) "Assign noffnc_id: "
      iret = nf_put_att_text(ncid, noffnc_id, 'long_name', 34, 
     &                       'Element based wetting/drying flag')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, noffnc_id, 'standard_name', 28,
     &                'element_wetting/drying_flag')
      CALL check_err(iret)

 
 
         WRITE(6,*) "Assign imhs_id: "

      iret = nf_put_att_text(ncid, imhs_id, 'long_name', 18, 
     &                       'model type')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, imhs_id, 'standard_name', 10,
     &                'model_type')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, iths_id, 'long_name', 59, 
     &                       'model time step number since the beginning
     & of the model run')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, iths_id, 'standard_name', 15,
     &           'model_time_step')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, iestp_id, 'long_name', 176, 
     &                       'line number (for ASCII output) or record
     & number (for binary output) of the most recent entry in the
     & elevation time series at specified elevation recording stations
     & output file')
      CALL check_err(iret)
    
      iret = nf_put_att_text(ncid, iestp_id, 'standard_name', 54,
     &      'line/record_number_of_last_entry_in_elev_rec_stations')
      CALL check_err(iret)
      
      iret = nf_put_att_text(ncid, nscoue_id, 'long_name', 150, 
     &                       'time step counter to determine when the
     & next entry will be written to the elevation time series at
     & specified elevation recording Stations output file')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, nscoue_id, 'standard_name', 15,
     &           'time_step_counter_for_next_entry_elev_rec_stations)')
      CALL check_err(iret)
      
      iret = nf_put_att_text(ncid, ivstp_id, 'long_name', 190, 
     &                       ' line number (for ASCII output) or record 
     & number (for binary output) of the most recent entry in the
     & depth-averaged velocity time series at specified velocity 
     & recording stations output file')
      CALL check_err(iret)
      
      iret = nf_put_att_text(ncid, ivstp_id, 'standard_name', 48,
     &               'record_number_of_last_entry_in_vel_rec_stations')
      CALL check_err(iret)

       iret = nf_put_att_text(ncid, nscouv_id, 'long_name', 165, 
     &                       'time step counter to determine when the
     & next entry will be written to the Depth-averaged Velocity Time
     & series at specified velocity recording stations output file.')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, nscouv_id, 'standard_name', 50,
     &           'time_step_counter_for_next_entry_vel_rec_stations')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, icstp_id, 'long_name', 190, 
     &                       'line number (for ASCII output) or record
     & number (for binary output) of the most recent entry in the scalar
     & concentration time series at specified concentration recording
     & stations output file')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, icstp_id, 'standard_name', 49,
     &               'record_number_of_last_entry_in_conc_rec_stations')
      CALL check_err(iret)
     
      iret = nf_put_att_text(ncid, nscouc_id, 'long_name', 165, 
     &                       'time step counter to determine when the
     & next entry will be written to the scalar concentration Time 
     &series at specified concentration Recording Stations output file')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, nscouc_id, 'standard_name', 51,
     &           'time_step_counter_for_next_entry_conc_rec_stations')
      CALL check_err(iret)

         WRITE(6,*) "Assign ipstp_id: "

      iret = nf_put_att_text(ncid, ipstp_id, 'long_name', 192, 
     &                       'line number (for ASCII output) or record
     & number (for binary output) of the most recent entry in the 
     &atmospheric pressure time series at specified meteorological 
     &recording stations')
      CALL check_err(iret)
     
      iret = nf_put_att_text(ncid, ipstp_id, 'standard_name', 57,
     &       'record_number_of_last_entry_of_atm_press_at_rec_stations')
      CALL check_err(iret)
     
         WRITE(6,*) "Assign iwstp_id: "
      iret = nf_put_att_text(ncid, iwstp_id, 'long_name', 190, 
     &                       'line number (for ASCII output) or record
     & number (for binary output) of the most recent entry in the Wind 
     & Velocity time series at specified meteorological recording 
     &stations')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, iwstp_id, 'standard_name', 56,
     &       'record_number_of_last_entry_of_wind_vel_at_rec_stations')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, nscoum_id, 'long_name', 244, 
     &                       'time step counter to determine when the
     & next entry will be written to the atmospheric pressure time 
     &series at specified meteorological recording stations and wind 
     &velocity time series at specified meteorological recording 
     &stations output files')
      CALL check_err(iret)

         WRITE(6,*) "Assign nscoum_id: "
      iret = nf_put_att_text(ncid, nscoum_id, 'standard_name', 60,
     &   'time_step_counter_of_atm_press_and_wind_vel_at_rec_stations')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, igep_id, 'long_name', 165, 
     &                       'line number (for ASCII output) or record 
     &number (for binary output) of the most recent entry in the 
     &Elevation time series at All Nodes in the Model Grid output file')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, igep_id, 'standard_name', 51,
     &       'record_number_of_last_entry_of_elev_at_model_nodes')
      CALL check_err(iret)

         WRITE(6,*) "Assign nscouge_id: "
      iret = nf_put_att_text(ncid, nscouge_id, 'long_name', 139, 
     &                       'time step counter to determine when the
     & next entry will be written to the  Elevation time series at All 
     &Nodes in the Model Grid output file')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, nscouge_id, 'standard_name', 41,
     &   'time_step_counter_of_elev_at_model_nodes')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, igvp_id, 'long_name', 179, 
     &                       'line number (for ASCII output) or record 
     &number (for binary output) of the most recent entry in the 
     &Depth-averaged velocity time series at All Nodes in the Model 
     &Grid output file')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, igvp_id, 'standard_name', 50,
     &       'record_number_of_last_entry_of_vel_at_model_nodes')
      CALL check_err(iret)

         WRITE(6,*) "Assign nscougv_id: "
      iret = nf_put_att_text(ncid, nscougv_id, 'long_name', 153, 
     &                       'time step counter to determine when the 
     &next entry will be written to the Depth-averaged velocity time 
     &series at All Nodes in the Model Grid output file')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, nscougv_id, 'standard_name', 40,
     &   'time_step_counter_of_vel_at_model_nodes')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, igcp_id, 'long_name', 176, 
     &                       'line number (for ASCII output) or record 
     &number (for binary output) of the most recent entry in the 
     &scalar Concentration time series at All Nodes in the Model 
     &Grid output file')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, igcp_id, 'standard_name', 51,
     &       'record_number_of_last_entry_of_conc_at_model_nodes')
      CALL check_err(iret)

         WRITE(6,*) "Assign nscougc_id: "
      iret = nf_put_att_text(ncid, nscougc_id, 'long_name', 150, 
     &                       'time step counter to determine when the 
     &next entry will be written to the scalar Concentration time 
     &series at All Nodes in the Model Grid output file')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, nscougc_id, 'standard_name', 41,
     &   'time_step_counter_of_conc_at_model_nodes')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, igpp_id, 'long_name', 176, 
     &                       'line number (for ASCII output) or record 
     &number (for binary output) of the most recent entry in the 
     &atmospheric pressure time series at All Nodes in the Model 
     &Grid output file')
      CALL check_err(iret)

         WRITE(6,*) "Assign igpp_id: "
      iret = nf_put_att_text(ncid, igpp_id, 'standard_name', 61,
     &       'record_number_of_last_entry_of_atm_press_at_model_nodes')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, igwp_id, 'long_name', 179, 
     &                       'line number (for ASCII output) or record 
     &number (for binary output) of the most recent entry in the 
     &Wind Stress or velocity time series at All Nodes in the Model 
     &Grid output file')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, igwp_id, 'standard_name', 55,
     &       'record_number_of_last_entry_of_wind_vel_at_model_nodes')
      CALL check_err(iret)

         WRITE(6,*) "Assign nscougw_id: "
      iret = nf_put_att_text(ncid, nscougw_id, 'long_name', 222, 
     &                       'time step counter to determine when the 
     &next entry will be written to the atmospheric pressure time series 
     &at All Nodes in the Model Grid and Wind Stress or velocity time 
     &series at All Nodes in the Model Grid output files')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, nscougw_id, 'standard_name', 59,
     &   'time_step_counter_of_atm_press_and_wind_vel_at_model_nodes')
      CALL check_err(iret)


!     -----------------------------------------------------------------------
!     Define node numbers on elevation specified boundary segments attributes
!     -----------------------------------------------------------------------      
      iret = nf_put_att_text(ncid, nopenc_id, 'long_name', 53, 
     &       'node numbers on elevation specified boundary segments')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, nopenc_id, 'units', 14, 
     &                      VAR_NOND)
      CALL check_err(iret)
      
      iret = nf_put_att_text(ncid, nvdllnc_id, 'long_name', 60, 
     &'number of nodes in each elevation specified boundary segment')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, nvdllnc_id, 'units', 14, 
     &                      VAR_NOND)
      CALL check_err(iret)
!     ---------------------------------------------------------------
!     Define node numbers on normal flow boundary segments attributes
!     ----------------------------------------------------------------  

      iret = nf_put_att_text(ncid, nbounc_id, 'long_name', 59, 
     & 'node numbers in each normal flow specified boundary segment')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, nbounc_id, 'units', 14, 
     &                      VAR_NOND)
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, ibtypenc_id, 'long_name', 13, 
     &       'boundary type')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, ibtypenc_id, 'units', 14, 
     &                      VAR_NOND)
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, nvellnc_id, 'long_name', 62, 
     &'number of nodes in each normal flow specified boundary segment')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, nvellnc_id, 'units', 14, 
     &                      VAR_NOND)
      CALL check_err(iret)

!    -----------------------
!     Define time attributes
!    -----------------------
      iret = nf_put_att_text(ncid,timenc_id,'long_name',10,'model time')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, timenc_id, 'standard_name', 4,'time')
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, timenc_id, 'units', 23,VAR_TUN)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, timenc_id, 'base_date',
     . (LEN(TRIM(base_date(1:INDEX(base_date,"!")-1)))), base_date)
      CALL check_err(iret)

!    -----------
!     Define u velocity
!    -----------
         WRITE(6,*) "Define u velocity: "

!    --------------------------------------------------
!     determine variable name according to coordinates
!    ---------------------------------------------------
      if (ics.eq.2) then
          iret = nf_put_att_text(ncid, u_id, 'long_name', 32, 
     &                       'vertically averaged e/w velocity')
       else
         iret = nf_put_att_text(ncid, u_id, 'long_name', 43, 
     &            'vertically averaged velocity in x-direction')
      endif
      CALL check_err(iret)

      iret = nf_put_att_text(ncid, u_id, 'standard_name', 10, 
     &                      'u_velocity')
      CALL check_err(iret)
!    -------------------------
!     determine variable units
!    -------------------------
      if (G.LT.11.D0) then
          iret = nf_put_att_text(ncid, u_id, 'units',13, VAR_MPS)
      else
          iret = nf_put_att_text(ncid, u_id, 'units',11, VAR_FPS)
      endif
      CALL check_err(iret)
!    --------------------------------------------------
!     determine variable name according to coordinates
!    ---------------------------------------------------
      if (ics.eq.2) then
          iret = nf_put_att_text(ncid, u_id, 'positive', 4, VAR_EAST)
      else
          iret = nf_put_att_text(ncid, u_id, 'positive', 5, VAR_RG)
      endif
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, u_id, 'missing_values', 
     &       NF_DOUBLE, 1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, u_id, '_FillValue', NF_DOUBLE, 1, 
     &                       doubleval)
      CALL check_err(iret)
      
      iret = nf_put_att_double(ncid, u_id, 'dry_Value', NF_DOUBLE, 1, 
     &                       doubleval)
      CALL check_err(iret)
!    -----------
!     Define v velocity
!    -----------
!    --------------------------------------------------
!     determine variable name according to coordinates
!    ---------------------------------------------------
      if (ics.eq.2) then
          iret = nf_put_att_text(ncid, v_id, 'long_name', 32, 
     &           'vertically averaged n/s velocity')
      else
          iret = nf_put_att_text(ncid, v_id, 'long_name', 43, 
     &           'vertically averaged velocity in y-direction')
      endif
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, v_id, 'standard_name', 10, 
     &                      'v_velocity')
      CALL check_err(iret)
!    -------------------------
!     determine variable units
!    -------------------------
      if (G.LT.11.D0) then
         iret = nf_put_att_text(ncid, v_id, 'units',13, VAR_MPS)
      else
         iret = nf_put_att_text(ncid, v_id, 'units',11, VAR_FPS)
      endif
      CALL check_err(iret)

!    --------------------------------------------------
!     determine variable name according to coordinates
!    ---------------------------------------------------
      if (ics.eq.2) then
          iret = nf_put_att_text(ncid, v_id, 'positive', 5, VAR_NORTH)
      else
          iret = nf_put_att_text(ncid, v_id, 'positive', 2, VAR_UP)
      endif
      CALL check_err(iret)

      iret = nf_put_att_DOUBLE(ncid, v_id, 'missing_values', 
     &       NF_DOUBLE, 1, doubleval)
      CALL check_err(iret)
!      doubleval(1) = -99999.
      iret = nf_put_att_DOUBLE(ncid, v_id, '_FillValue', NF_DOUBLE, 1, 
     &                       doubleval)
      CALL check_err(iret)
!      doubleval(1) = -99999.
      iret = nf_put_att_DOUBLE(ncid, v_id, 'dry_Value', NF_DOUBLE, 1, 
     &                       doubleval)
      CALL check_err(iret)

!    --------------------------------------------
!     Define sea surface elevation attributes
!    --------------------------------------------
      
      iret = nf_put_att_text(ncid, zeta1_id,'long_name',44, 
     &'surface elevation at the previous time step')
                        
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, zeta1_id, 'standard_name', 44, 
     &         'surface_elevation_at_previous_time step')
      CALL check_err(iret)
      
      if (G.LT.11.D0) then
          iret = nf_put_att_text(ncid, zeta1_id, 'units', 6, VAR_METS)
          CALL check_err(iret)
      else
         iret = nf_put_att_text(ncid, zeta1_id, 'units', 4, VAR_FEET)
         CALL check_err(iret)
      endif      
      
      iret = nf_put_att_DOUBLE(ncid,zeta1_id,'missing_values',NF_DOUBLE,
     &1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_DOUBLE(ncid, zeta1_id, '_FillValue', NF_DOUBLE, 
     &                      1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, zeta1_id,'positive', 2, VAR_UP)
      CALL check_err(iret)
!    --------------------------------------------
!     Define sea surface elevation attributes
!    --------------------------------------------
      
      iret = nf_put_att_text(ncid, zeta2_id,'long_name',43, 
     & 'surface elevation at the current time step')
                       
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, zeta2_id, 'standard_name', 43, 
     &         'surface_elevation_at_the_current_time_step')
      CALL check_err(iret)
      
      if (G.LT.11.D0) then
          iret = nf_put_att_text(ncid, zeta2_id, 'units', 6, VAR_METS)
          CALL check_err(iret)
      else
         iret = nf_put_att_text(ncid, zeta2_id, 'units', 4, VAR_FEET)
         CALL check_err(iret)
      endif      
      
      iret = nf_put_att_DOUBLE(ncid,zeta2_id,'missing_values',NF_DOUBLE,
     &1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_DOUBLE(ncid, zeta2_id, '_FillValue', NF_DOUBLE, 
     &                      1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, zeta2_id,'positive', 2, VAR_UP)
      CALL check_err(iret)
!    --------------------------------------------
!     Define sea surface elevation attributes
!    --------------------------------------------

      iret = nf_put_att_text(ncid, zetad_id,'long_name',37,       
     &    'elevation at flux specified boundary')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, zetad_id, 'standard_name', 37, 
     &                      'elevation_at_flux_specified_boundary')
      CALL check_err(iret)
      
      if (G.LT.11.D0) then
          iret = nf_put_att_text(ncid, zetad_id, 'units', 6, VAR_METS)
          CALL check_err(iret)
      else
         iret = nf_put_att_text(ncid, zetad_id, 'units', 4, VAR_FEET)
         CALL check_err(iret)
      endif      
      
      iret = nf_put_att_DOUBLE(ncid,zetad_id,'missing_values',NF_DOUBLE,
     &1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_DOUBLE(ncid, zetad_id, '_FillValue', NF_DOUBLE, 
     &                      1, doubleval)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, zetad_id,'positive', 2, VAR_UP)
      CALL check_err(iret)
!     -----------------
!     Global attributes
!     -----------------
         WRITE(6,*) "Global attributes: "

      iret = nf_put_att_text(ncid, NF_GLOBAL, 'model', 6, 'ADCIRC')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'version', 
     &       LEN(TRIM(ADC_VERSION)), ADC_VERSION)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'grid_type', 10, 
     &                      'Triangular')
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'description', 
     &       LEN(TRIM(rundes(1:INDEX(rundes,"!")-1))), rundes)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'grid',
     &       LEN(TRIM(agrid)), agrid)       
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'title',
     &       LEN(TRIM(title(1:INDEX(title,"!")-1))), 
     &       title)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'institution', 
     &       LEN(TRIM(institution(1:INDEX(institution,"!")-1))), 
     &       institution)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'source', 
     &       LEN(TRIM(source(1:INDEX(source,"!")-1))),
     &       source)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'history', 
     &       LEN(TRIM(history(1:INDEX(history,"!")-1))),
     &       history)
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
     &       LEN(TRIM(host(1:INDEX(host,"!")-1))), 
     &       host)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'convention', 
     &       LEN(TRIM(convention(1:INDEX(convention,"!")-1))), 
     &       convention)
      CALL check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'contact', 
     &       LEN(TRIM(contact(1:INDEX(contact,"!")-1))), 
     &       contact)
      CALL check_err(iret)
 
 
      CALL date_and_time(now_date,big_ben,zone,values)
      write(date_string,71) values(1),values(2),values(3)
     & ,values(5),values(6),values(7),   (values(4))/60
 71   format(I4,'-',I2.2,'-',i2.2,' ',i2,':',i2.2,':',i2.2,' '
     &  ,i3.2,':00')
      iret = nf_put_att_text(ncid, NF_GLOBAL,'creation_date'
     & ,len_trim(date_string),date_string)
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
      iret = nf_put_att_int(ncid, NCGLOBAL, 'irtype', NF_INT, 1, irtype)
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
      iret = nf_put_att_double(ncid, NCGLOBAL, 'ffactor', NF_DOUBLE, 1, 
     &                         ffactor)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, NCGLOBAL, 'esl', NF_DOUBLE, 1, 
     &                         esl)
      CALL check_err(iret)
      iret = nf_put_att_double(ncid, NCGLOBAL, 'cori', NF_DOUBLE, 1, 
     &                         cori)
      CALL check_err(iret)
      iret = nf_put_att_int(ncid, NCGLOBAL, 'ntif', NF_INT, 1, ntif)
      CALL check_err(iret)
      iret = nf_put_att_int(ncid, NCGLOBAL, 'nbfr', NF_INT, 1, nbfr)
      CALL check_err(iret)
     
!     ======================
!     Leave define mode
!     ======================
      iret = nf_enddef(ncid)
      CALL check_err(iret)
!     ======================
!     Store static variables
!     ======================
!     -------
!     Store X
!     -------
      iret = nf_put_var_double(ncid, X_id, X)
      CALL check_err(iret)
!     -------
!     Store Y
!     -------
      iret = nf_put_var_double(ncid, Y_id, Y)
      CALL check_err(iret)
!     ------------
!     Store depth
!     ------------
      iret = nf_put_var_double(ncid, DEPTH_id, depth)
      CALL check_err(iret)
!     --------------
!     Store elements
!     --------------
      count(1)=nface_len
      count(2)=elem_len
      start(1)=1
      start(2)=1
      iret=nf_put_vara_int(ncid,ele_id,start,count,element)
      CALL check_err(iret)
!     ------------------------------------
!     Store elevation boundary information
!     -------------------------------------
      iret = nf_put_var_int(ncid, nopenc_id, nopenc)
      CALL check_err(iret)
      iret = nf_put_var_int(ncid, nvdllnc_id, nvdllnc)
      CALL check_err(iret)
!     --------------------------------------
!     Store normal flow boundary information
!     --------------------------------------
      iret = nf_put_var_int(ncid, nbounc_id, nbounc)
      CALL check_err(iret)
      iret = nf_put_var_int(ncid, ibtypenc_id, ibtypenc)
      CALL check_err(iret)
      iret = nf_put_var_int(ncid, nvellnc_id, nvellnc)
      CALL check_err(iret)      
         WRITE(6,*) "Store maximum depth: "
!     --------------------
!     Store maximum depth
!     --------------------
      maxdepth = maxval (depth)
      iret = nf_put_var_double(ncid, maxdepth_id, maxdepth)
      CALL check_err(iret)
!     --------------------
!     Store minimum depth
!     --------------------
      mindepth = minval (depth)
      iret = nf_put_var_double(ncid, mindepth_id, mindepth)
      CALL check_err(iret)
      
!     =====================
!     Open old netcdf file
!     =====================
      ELSE
          iret = NF_OPEN(FILENAME, NF_WRITE, ncid)
          CALL check_err(iret)   
!     ------------------
!     Inquire variables
!     ------------------
          iret=nf_inq_unlimdim(ncid,recid)
          CALL check_err(iret)   
          iret=nf_inq_dimlen(ncid,recid,nrec)
          CALL check_err(iret)   
      ENDIF
!     ------------------
!     Update Global Attributes
!     ------------------
!      IF(IHOT.NE.0.AND.FIRST67.EQ.1) THEN
       IF((IHOT.NE.0.AND.FIRST67.EQ.1).OR.
     &    (IHOT.NE.0.AND.FIRST68.EQ.1)) THEN
         iret = NF_REDEF (ncid)
         CALL check_err(iret)
         iret = nf_put_att_int(ncid, NCGLOBAL, 'ihot', NF_INT, 1, ihot)
         CALL check_err(iret)
         iret = nf_put_att_int(ncid, NCGLOBAL, 'nramp', NF_INT, 1,nramp)
         CALL check_err(iret)
         iret = nf_put_att_double(ncid, NCGLOBAL, 'tau0', NF_DOUBLE, 1, 
     &                         tau0)
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
         iret = nf_put_att_double(ncid, NCGLOBAL, 'ffactor', NF_DOUBLE, 1, 
     &                         ffactor)
         CALL check_err(iret)
         iret = nf_put_att_double(ncid, NCGLOBAL, 'esl', NF_DOUBLE, 1, 
     &                         esl)
         CALL check_err(iret)
         CALL date_and_time(now_date,big_ben,zone,values)
         WRITE(date_string,71) values(1),values(2),values(3),
     &         values(5),values(6),values(7),   (values(4))/60
         iret = nf_put_att_text(ncid, NF_GLOBAL,'modification_date',
     &          len_trim(date_string),date_string)
         CALL check_err(iret)
         iret = NF_ENDDEF (ncid)
         CALL check_err(iret)
       ENDIF
       PRINT *,"netcdfio: end mod"

!      IF(ncmode.eq.0) go to 9999
!      IF((FIRST67.EQ.1.and.
!     & ((ncmode.ne.0.and.IHOT.NE.0).or.(ncmode.eq.2.and.IHOT.eq.0).or.
!     &  (ncmode.eq.2.and.IHOT.ne.0.and.nout.eq.-3))).OR.
!     & (FIRST68.EQ.1.and.
!     & ((ncmode.ne.0.and.IHOT.NE.0).or.(ncmode.eq.2.and.IHOT.eq.0).or.
!     &  (ncmode.eq.2.and.IHOT.ne.0.and.nout.eq.-3)))
!     &  )THEN
      IF( (ncmode.ne.2.and.IHOT.NE.0))THEN

              WRITE(6,*) "FIRST67.EQ.1. ncmode.ne.0"
              WRITE(6,*) "FIRST68.EQ.1. ncmode.ne.0"

         iret=nf_inq_varid(ncid,"time",timenc_id)
         iret=nf_inq_varid(ncid,"zeta1",zeta1_id)
         iret=nf_inq_varid(ncid,"zeta2",zeta2_id)
         iret=nf_inq_varid(ncid,"zetad",zetad_id)
         iret=nf_inq_varid(ncid,"u-vel",u_id)
         iret=nf_inq_varid(ncid,"v-vel",v_id)
         iret=nf_inq_varid(ncid,"nodecode",nodecodenc_id)
         iret=nf_inq_varid(ncid,"noff",noffnc_id)
         iret=nf_inq_varid(ncid,"imhs",imhs_id)
         iret=nf_inq_varid(ncid,"iths",iths_id)
         iret=nf_inq_varid(ncid,"iestp",iestp_id)
         iret=nf_inq_varid(ncid,"nscoue",nscoue_id)
         iret=nf_inq_varid(ncid,"ivstp",ivstp_id)
         iret=nf_inq_varid(ncid,"nscouv",nscouv_id)
         iret=nf_inq_varid(ncid,"icstp",icstp_id)
         iret=nf_inq_varid(ncid,"nscouc",nscouc_id)
         iret=nf_inq_varid(ncid,"ipstp",ipstp_id)
         iret=nf_inq_varid(ncid,"iwstp",iwstp_id)
         iret=nf_inq_varid(ncid,"nscoum",nscoum_id)
         iret=nf_inq_varid(ncid,"igep",igep_id)
         iret=nf_inq_varid(ncid,"nscouge",nscouge_id)
         iret=nf_inq_varid(ncid,"igvp",igvp_id)
         iret=nf_inq_varid(ncid,"nscougv",nscougv_id)
         iret=nf_inq_varid(ncid,"igcp",igcp_id)
         iret=nf_inq_varid(ncid,"nscougc",nscougc_id)
         iret=nf_inq_varid(ncid,"igpp",igpp_id)
         iret=nf_inq_varid(ncid,"igwp",igwp_id)
         iret=nf_inq_varid(ncid,"nscougw",nscougw_id)
      ENDIF
!      IF(IHOT.NE.0.AND.FIRST67.EQ.1.and.ncmode.ne.0.and.nout.eq.3)THEN
!         iret = nf_get_var_double(ncid, minun_id, minun)
!         CALL check_err(iret)
!      ENDIF
      
      IF(FILEUNIT.EQ.67)FIRST67=0
      IF(FILEUNIT.EQ.68)FIRST68=0
!      NREC = NREC + 1
      NREC = 1
!     ----------
!     Store time
!     ----------
      starti(1)=nrec
      counti(1)=timenc_len
      iret = nf_get_vara_DOUBLE(ncid, timenc_id, starti,counti, timenc)
      CALL check_err(iret)
      PRINT *," NETCDFIO: TIMENC ", timenc(timenc_len)
!     ---------------------------
!     Store zeta1
!     ---------------------------
              WRITE(6,*) "Store zeta1"
      count(1)=node_len
      count(2)=timenc_len
      start(1)=1
      start(2)=nrec
      iret = nf_get_vara_DOUBLE(ncid, zeta1_id, start, count, zeta1)
      CALL check_err(iret)
      DO i=1, node_len
         WRITE(6,*) "netcdfio: zeta1",i, zeta1(i,1)
      END DO
!     ---------------------------
!     Store zeta2
!     ---------------------------
              WRITE(6,*) "Store zeta2"
      count(1)=node_len
      count(2)=timenc_len
      start(1)=1
      start(2)=nrec
      iret = nf_get_vara_DOUBLE(ncid, zeta2_id, start, count, zeta2)
      CALL check_err(iret)
      DO i=1, node_len
         WRITE(6,*) "netcdfio: zeta2",i, zeta2(i,1)
      END DO
!     ---------------------------
!     Store zetad
!     ---------------------------
              WRITE(6,*) "Store zetad"
      count(1)=node_len
      count(2)=timenc_len
      start(1)=1
      start(2)=nrec
      iret = nf_get_vara_DOUBLE(ncid, zetad_id, start, count, zetad)
      CALL check_err(iret)

!     ---------------------------
!     Store u
!     ---------------------------
              WRITE(6,*) "Store u"
      count(1)=node_len
      count(2)=timenc_len
      start(1)=1
      start(2)=nrec
      iret = nf_get_vara_DOUBLE(ncid, u_id, start, count, u)
      CALL check_err(iret)
!     ---------------------------
!     Store v
!     ---------------------------
      count(1)=node_len
      count(2)=timenc_len
      start(1)=1
      start(2)=nrec
      iret = nf_get_vara_DOUBLE(ncid, v_id, start, count, v)
      CALL check_err(iret)
!     ---------------
!     Store nodecode
!     ---------------
      WRITE(6,*) "Store nodecodenc"
      iret = nf_get_var_int(ncid, nodecodenc_id, nodecodenc)
      CALL check_err(iret)
!!      WRITE(6,*) "Store nodecodenc"
!!      start(1)=1
!!      count(1)=node_len
!!      iret = nf_get_vara_int(ncid, nodecodenc_id, start, count, 
!!     &                       nodecodenc)
!!      CALL check_err(iret)
      DO i=1, node_len
           WRITE(6,*) "netcdfio: nodecodenc",i, nodecodenc(i)
      END DO
!     ---------------
!     Store noff
!     ---------------
       WRITE(6,*) "Store noffnc"
       iret = nf_get_var_int(ncid, noffnc_id, noffnc)
       CALL check_err(iret)
       WRITE(6,*) "finish getting noffnc"
!!      WRITE(6,*) "Store noffnc"
!!      WRITE(6,*) "elem_len = ",elem_len
!!      WRITE(6,*) "noffnc_id = ",noffnc_id
!!      start(1)=1
!!      count(1)=2
!      count(1)=elem_len

!!      iret = nf_get_vara_int(ncid, noffnc_id, start, count, noffnc)
!!      CALL check_err(iret)
      
      DO i=1, elem_len
         WRITE(6,*)  "netcdfio: noff",i, noffnc(i)
      END DO
!     ------------------------
!     Store hotstart variables
!     ------------------------
               WRITE(6,*) "Store hotstart variables"
      iret = nf_get_var_int(ncid, imhs_id, imhs)
      CALL check_err(iret)
      iret = nf_get_var_int(ncid, iths_id, iths)
      CALL check_err(iret)
      iret = nf_get_var_int(ncid, iestp_id, iestp)
      CALL check_err(iret)
      iret = nf_get_var_int(ncid, nscoue_id, nscoue)
      CALL check_err(iret)
      iret = nf_get_var_int(ncid, ivstp_id, ivstp)
      CALL check_err(iret)
      iret = nf_get_var_int(ncid, nscouv_id, nscouv)
      CALL check_err(iret)
      iret = nf_get_var_int(ncid, icstp_id, icstp)
      CALL check_err(iret)
      iret = nf_get_var_int(ncid, nscouc_id, nscouc)
      CALL check_err(iret)
      iret = nf_get_var_int(ncid, ipstp_id, ipstp)
      CALL check_err(iret)
      iret = nf_get_var_int(ncid, iwstp_id, iwstp)
      CALL check_err(iret)
      iret = nf_get_var_int(ncid, nscoum_id, nscoum)
      CALL check_err(iret)
               WRITE(6,*) "Store igep"
      iret = nf_get_var_int(ncid, igep_id, igep)
               WRITE(6,*) "IGEP = ", IGEP
      CALL check_err(iret)
      iret = nf_get_var_int(ncid, nscouge_id, nscouge)
      CALL check_err(iret)
      iret = nf_get_var_int(ncid, igvp_id, igvp)
      CALL check_err(iret)
      iret = nf_get_var_int(ncid, nscougv_id, nscougv)
      CALL check_err(iret)
      iret = nf_get_var_int(ncid, igcp_id, igcp)
      CALL check_err(iret)
      iret = nf_get_var_int(ncid, nscougc_id, nscougc)
      CALL check_err(iret)
      iret = nf_get_var_int(ncid, igpp_id, igpp)
      CALL check_err(iret)
      iret = nf_get_var_int(ncid, igwp_id, igwp)
      CALL check_err(iret)
      iret = nf_get_var_int(ncid, nscougw_id, nscougw)
      CALL check_err(iret)
!     ---------------------------
!     Closing fort.X4 netCDF file
!     ---------------------------
!      IF(NREC.EQ.NDATA) THEN
!         iret = nf_put_var_double(ncid, minun_id, minun)
!         CALL check_err(iret)
!         iret = nf_put_var_double(ncid, maxun_id, maxun)
!         CALL check_err(iret)
!         iret = nf_put_var_double(ncid, minvn_id, minvn)
!         CALL check_err(iret)
!         iret = nf_put_var_double(ncid, maxvn_id, maxvn)
!         CALL check_err(iret)
!      ENDIF

!     ---------------------------
!     Closing netCDF file
!     ---------------------------
      WRITE(6,*) "Closing ", FILENAME, " netCDF file ... "
      WRITE(6,*) "-----------------------------------------------------"
      
9999  iret = nf_close(ncid)
      CALL check_err(iret)
      WRITE(6,*)'end subroutine readnc67'

      end subroutine readnc67     


      END MODULE NETCDFIO

