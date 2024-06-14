!-------------------------------------------------------------------------------!
!
! ADCIRC - The ADvanced CIRCulation model
! Copyright (C) 1994-2023 R.A. Luettich, Jr., J.J. Westerink
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
C----------------------------------------------------------------------------
C
C                           MODULE PRE_GLOBAL
C
C----------------------------------------------------------------------------
C
C                  For use with ADCPREP Version 2.2 (  09/20/2006 )
C
C                     current for ADCIRC v46.34   09/20/2006
C----------------------------------------------------------------------------

      MODULE PRE_GLOBAL
      USE ADC_CONSTANTS,ONLY: REarth,DEG2RAD,RAD2DEG
      USE GLOBAL_3DVS, ONLY: NFEN
      USE PRESIZES
      USE VERSION
      USE HASHTABLE
      IMPLICIT NONE

C     Land Boundary Nodes and Segment Definitions:

C     PARM14 section
      INTEGER NELG,NNODG
      INTEGER NOPE, NETA, NETA_MAX
      INTEGER NBOU  ! Number of Global Land Boundary Segments
      INTEGER NVEL  ! Total Number of Global Land Boundary Nodes
      INTEGER,ALLOCATABLE::NNEG(:,:)
      INTEGER,ALLOCATABLE::NVDLL(:)
      INTEGER,ALLOCATABLE::NBDV(:,:)
      INTEGER,ALLOCATABLE::NVELL(:)  ! NVELL(K) Number of Global Land Boundary
                                     ! Nodes of Segment K
      INTEGER,ALLOCATABLE::NBVV(:,:) ! NBVV(K,I) Global Node Number of
C                                    ! I-th Node on Land Boundary Segment K
      INTEGER,ALLOCATABLE::IBTYPE(:) ! Type of discharge boundary segment K

      INTEGER,ALLOCATABLE::IBTYPEE(:) ! Type of elevation boundary segment K
      INTEGER,ALLOCATABLE::IBCONNR(:,:)
      INTEGER,ALLOCATABLE::LBCODE(:) ! LBCODE(I) Boundary Type of Land Boundary
                                     ! Node I
      INTEGER :: numLBCodeValues ! number of values in LBCODE array
      INTEGER,ALLOCATABLE::WEIR(:)
      INTEGER,ALLOCATABLE::WEIRD(:)

C.... DW, periodic bc
      INTEGER:: NPERBC 
      INTEGER, ALLOCATABLE:: IPERCONN(:,:), IPERPTR(:) ! a simple sparse CRS format 
C.... DW

      TYPE(ipair),ALLOCATABLE, TARGET :: node_dict(:)
      INTEGER,ALLOCATABLE             :: labels(:)

C     METIS interface section
      INTEGER,ALLOCATABLE ::  PROC(:)
C
C     GRID14 section
      REAL(8)  ,ALLOCATABLE ::    X(:),Y(:),SLAM(:),SFEA(:)
      REAL(8) ,ALLOCATABLE ::   DP(:)
      REAL(8) SL0,SF0            ! center of cpp projection, in radians
      REAL(8), ALLOCATABLE :: SL1(:),SF1(:) ! nodal coordinates, in radians
C
C     BARRIER14 section
      REAL(8) ,ALLOCATABLE ::  BAR1(:,:),BAR2(:,:),BAR3(:,:)
      REAL(8) ,ALLOCATABLE ::  BAR4(:,:),BAR5(:,:),BAR6(:,:)

C     FLOWBC section
      INTEGER NFLBN    ! NFLBN Number of Global Flow Boundary Nodes
      INTEGER NFLBNP   ! NFLBNP(PE) Number of Flow boundary Nodes on PE
      INTEGER,ALLOCATABLE::FLBN(:)  ! FLBN(I) Global Node number of
C                                   ! Ith Flow Boundary Node
      INTEGER,ALLOCATABLE::FLBNX(:) ! FLBNX(I) Index of Ith Flow Boundary Node
      INTEGER,ALLOCATABLE::FLBNXP(:)! FLBNXP(I) Index of Ith Flow Boundary Node
      INTEGER,ALLOCATABLE::FLUX14_ARY(:)
      REAL(8),ALLOCATABLE::FLUX20_ARY(:,:)
      INTEGER EXIST_FLUX
      LOGICAL APERIODIC_FLOW_BC     ! .true. if EXIST_FLUX is .true. .AND.
                                    ! NFFR=0 in fort.15
                                    
      ! @jasonfleming: add boundary weights for perturbing the domain
      ! decomposition in subroutine metis()
      logical :: strictBoundaries ! true if subdomain boundaries must match fulldomain
      integer, allocatable :: boundaryWeights(:) ! mnp

! kmd - updates for rivers in a baroclinic simulation
      INTEGER EXIST_BC_TS
      LOGICAL APERIODIC_BC_TS
      INTEGER, ALLOCATABLE :: BCTS14_ARY(:)
! kmd - Evan's changes for rivers above MSL
      LOGICAL River_above_MSL
C
C     STRING14 section
      CHARACTER(80) AGRID,NSIZES,SIZEMSG
      CHARACTER(80) NOPEMSG,NETAMSG
      CHARACTER(80) NBOUMSG,NVELMSG
      CHARACTER(80),ALLOCATABLE   ::  NVELLMSG(:),NVDLLMSG(:)
C
C     PARM15  section
      INTEGER  NOUTE,NSPOOLE,NSTAE,NSTAE_MAX
      INTEGER  NOUTV,NSPOOLV,NSTAV,NSTAV_MAX
      INTEGER  NOUTM,NSPOOLM,NSTAM,NSTAM_MAX
      INTEGER  NOUTC,NSPOOLC,NSTAC,NSTAC_MAX
      INTEGER  NFOVER,NABOUT,NSCREEN
      INTEGER  IHOT,ICS,IM
      INTEGER  NOLIBF,NOLIFA,NOLICA,NOLICAT
      INTEGER  NWP,NCOR,NTIP,NWS,NRAMP
      INTEGER  NRS
      INTEGER  NCICE  !tcm v49.64.01 added ice concentration
      INTEGER  NTIF,NBFR,NFFR
      INTEGER  NOUTGE,NSPOOLGE
      INTEGER  NOUTGV,NSPOOLGV
      INTEGER  NOUTGC,NSPOOLGC
      INTEGER  NOUTGW,NSPOOLGW
      INTEGER  NHARFR,NHSINC
      INTEGER  NHSTAR
      INTEGER  IREFYR,IREFMO,IREFDAY,IREFHR,IREFMIN
      INTEGER  NWLAT,NWLON,NFLUXF,NFLUXIB,NFLUXIBP,NFLUXB
      INTEGER,ALLOCATABLE ::   NNSEG(:),NNSVG(:),NNSCG(:),NNSMG(:)
C
C     jgf46.00 Nodal attributes section.
      CHARACTER(len=80), ALLOCATABLE :: useNodalAttrNames(:) ! labels
C
C     STRING15 section
C
      CHARACTER(80) RUNDES,RUNID,OVERMSG,ABOUTMSG,SCREENMSG,HOTMSG
      CHARACTER(80) ICSMSG,IMMSG,IBFMSG,IFAMSG,ICAMSG,ICATMSG,NWPMSG,
     &  NCORMSG
      CHARACTER(80) NTIPMSG,NWSMSG,RAMPMSG,GMSG,TAU0MSG,DTMSG,
     &  STATMSG,REFTMSG
      CHARACTER(80) RNDAYMSG,DRAMPMSG,COEFMSG,H0MSG,SLMSG,TAUMSG,
     &  ESLMSG,CORIMSG
      CHARACTER(80)  NTIFMSG,NBFRMSG,ANGMSG,NFFRMSG,STAEMSG,
     &  NSTAEMSG,STAVMSG,RSTIMMSG
      CHARACTER(80)  NSTAVMSG,STACMSG,NSTACMSG,STAMMSG,NSTAMMSG
      CHARACTER(80)  OUTGEMSG,OUTGVMSG,OUTGCMSG,OUTGWMSG
      CHARACTER(80)  HARFRMSG,HARPARMSG,OUTHARMSG,HSTARMSG
      CHARACTER(80)  SOLVMSG,WSMSG1,WSMSG2
      CHARACTER(132),ALLOCATABLE :: STAELOC(:),STAVLOC(:),STACLOC(:)
      CHARACTER(132),ALLOCATABLE :: STAMLOC(:)
      CHARACTER*80,ALLOCATABLE :: HAFREMSG(:)
      CHARACTER*80,ALLOCATABLE :: TIPOTAG(:),BOUNTAG(:),FBOUNTAG(:)
      CHARACTER*80,ALLOCATABLE :: ALPHA1(:),ALPHA2(:),FREQMSG(:),
     &  QNMSG(:,:)
      CHARACTER*80,ALLOCATABLE :: AMIGMSG(:),EMOMSG(:,:),TPKMSG(:)

C
C     PARM15-3DVS section
      INTEGER :: ISLIP,ICG,IEVC
      INTEGER :: I3DSD,NSPO3DSD       ! NSTA3DD is in presizes
      INTEGER :: I3DSV,NSPO3DSV       ! NSTA3DV is in presizes
      INTEGER :: I3DST,NSPO3DST       ! NSTA3DT is in presizes
      REAL(8) :: TO3DSDS,TO3DSDF
      REAL(8) :: TO3DSVS,TO3DSVF
      REAL(8) :: TO3DSTS,TO3DSTF
      INTEGER,ALLOCATABLE :: ISDHOUT(:),ISVHOUT(:),ISTHOUT(:)
      INTEGER :: I3DGD,NSPO3DGD
      INTEGER :: I3DGV,NSPO3DGV
      INTEGER :: I3DGT,NSPO3DGT
      REAL(8) :: TO3DGDS,TO3DGDF
      REAL(8) :: TO3DGVS,TO3DGVF
      REAL(8) :: TO3DGTS,TO3DGTF
      REAL(8) :: KP,Z0S,Z0B,ALP1,ALP2,ALP3,EVMIN,EVCON,THETA1,THETA2
      REAL(8),ALLOCATABLE,TARGET :: X3DD(:),Y3DD(:),SX3DD(:),SY3DD(:)
      REAL(8),ALLOCATABLE,TARGET :: SL3DD(:),SF3DD(:)
      REAL(8),ALLOCATABLE,TARGET :: X3DV(:),Y3DV(:),SX3DV(:),SY3DV(:)
      REAL(8),ALLOCATABLE,TARGET :: SL3DV(:),SF3DV(:)
      REAL(8),ALLOCATABLE,TARGET :: X3DT(:),Y3DT(:),SX3DT(:),SY3DT(:)
      REAL(8),ALLOCATABLE,TARGET :: SL3DT(:),SF3DT(:)
      CHARACTER*132,ALLOCATABLE :: STA3DDLOC(:),STA3DVLOC(:),
     &   STA3DTLOC(:)
      INTEGER,ALLOCATABLE ::   NNS3DDG(:),NNS3DVG(:),NNS3DTG(:)
      REAL(8),ALLOCATABLE :: EVTot(:) !(NFEN) Vertical eddy viscosity
      INTEGER IGC                      ! vertical grid code in fort.15

C   kmd48.33bc variables for 3D baroclinic and hot start
      REAL(8) :: RBCTIMEINC, SBCTIMEINC, TBCTIMEINC    ! flag for bcs on baroclinic
      REAL(8) :: BCSTATIM, SBCSTATIM, TBCSTATIM       ! flag for bcs on baroclinic
      REAL(8) :: TTBCSTATIM, TTBCTIMEINC       ! flag for bcs on baroclinic
      INTEGER N3DSD, I3DSDRec, N3DSV, I3DSVRec, N3DST, I3DSTRec
      INTEGER N3DGD, I3DGDRec, N3DGV, I3DGVRec, N3DGT, I3DGTRec
      REAL(8),ALLOCATABLE :: DUU(:),DUV(:),DVV(:)
      REAL(8),ALLOCATABLE :: UU(:),VV(:)
      REAL(8),ALLOCATABLE :: BSX(:),BSY(:)
      REAL(8),ALLOCATABLE :: WZ(:,:), q20(:,:), l(:,:)
      REAL(8),ALLOCATABLE :: SigT(:,:), Sal(:,:), Temp(:,:)
      REAL(8),ALLOCATABLE :: RealQ(:,:), ImagQ(:,:)
      INTEGER,ALLOCATABLE :: NODECODE(:)
      INTEGER IESTP,NSCOUE,IVSTP,NSCOUV,ICSTP,NSCOUC,IPSTP,IWSTP,NSCOUM,
     &        IGEP,NSCOUGE,IGVP,NSCOUGV,IGCP,NSCOUGC,IGPP,IGWP,NSCOUGW

C     jgf45.12 Terms added for the transport parameters
C     diffusion coefficients for the salinity field
      REAL(8) :: NLSD, NVSD
C     diffusion coefficients for the temperature field
      REAL(8) :: NLTD, NVTD
c     time weighting coefficient for the vertical diffusion term in the
c     transport equation
      REAL(8) :: Alp4
C     Flag for the top boundary condition for the temperature field in
C     the transport equation
!      INTEGER :: NTF ! kmd48.33bc took out with new evaluation of boundary condition
C
C     STRING15-3DVS section
      CHARACTER(80) :: IDENMSG,SLIPMSG,Z0MSG,ALPMSG,FEMSG
      CHARACTER(80) :: EVCMSG,THETAMSG
      CHARACTER(80) :: DSDMSG,DSVMSG,DSTMSG,DGDMSG,DGVMSG,DGTMSG
      CHARACTER(80) :: NSTA3DDMSG,NSTA3DVMSG,NSTA3DTMSG

C  kmd48.33bc variables for boundary conditions
      CHARACTER(80) :: RESBCFLAGMSG, BCTIMEMSG, BCSTATMSG, TBCTIMEMSG
      CHARACTER(80) :: SPONGEDISTMSG, EqnstateMSG
C
C     ELESTAT section
      REAL(8) TOUTSE,TOUTFE,TOUTSGE,TOUTFGE
      REAL(8),ALLOCATABLE, TARGET :: XEL(:),YEL(:),SLEL(:),SFEL(:)
C
C     VELSTAT section
      REAL(8) TOUTSGV,TOUTFGV,TOUTSV,TOUTFV
      REAL(8),ALLOCATABLE, TARGET :: XEV(:),YEV(:),SLEV(:),SFEV(:)
C
C     CONSTAT section
      REAL(8) TOUTSC,TOUTFC,TOUTSGC,TOUTFGC
      REAL(8),ALLOCATABLE, TARGET :: XEC(:),YEC(:),SLEC(:),SFEC(:)
C
C     METSTAT section
      REAL(8) TOUTSM,TOUTFM,TOUTSGM,TOUTFGM
      REAL(8),ALLOCATABLE, TARGET ::  XEM(:),YEM(:),SLEM(:),SFEM(:)
      REAL(8) TOUTSGW,TOUTFGW
C
C     INPUT15 section
      INTEGER    NODEDRYMIN,NODEWETRMP
      REAL(8)   G
      REAL(8)   TAU0,RNDAY,DRAMP
      REAL(8)   A00,B00,C00,H0,VELMIN
      REAL(8)    SLAM0,SFEA0
      REAL(8)   TAU,CF,ESL,CORI,ANGINN
      REAL(8)    DTDP
      REAL(8)   ESLM,ESLC,HBREAK,FTHETA,FGAMMA
      REAL(8)    REFSEC
      REAL(8)   WLATMAX,WLONMIN,WLATINC,WLONINC
      REAL(8),ALLOCATABLE ::  FF(:),FACE(:)
      REAL(8) ,ALLOCATABLE ::  EMO(:,:),EFA(:,:)
      REAL(8),ALLOCATABLE ::  TPK(:),ETRF(:),FFT(:),FACET(:)
      REAL(8),ALLOCATABLE ::  QNAM(:,:),QNPH(:,:)
      REAL(8),ALLOCATABLE ::  FFF(:),FFACE(:)
C
C     Time varying tau0 min and max
      CHARACTER(80) :: TAU0LIMMSG
      REAL(8) Tau0FullDomainMin, Tau0FullDomainMax
C
C     INPUT15D section
      REAL(8)   DT
      REAL(8)    STATIM,REFTIM,WTIMINC
      REAL(8),ALLOCATABLE ::   AMIG(:),AMIGT(:),FAMIG(:)
      REAL(8)  CICE_TIMINC   ! tcm v49.64.01 addition for ice

c.... tcm v50.66.02 added for time varying bathymetry
      INTEGER  NDDT   !tcm v50.66.02 added for time varying bathymetry
      REAL(8)  BTIMINC ! TIME INCREMENT (SECONDS) FOR BATHYMETRY CHANGES USED FOR TIME VARYING BATHYMETRY
      REAL(8) :: BCHGTIMINC ! time increment (seconds) over which bathymetry changes during a btimeinc interval
      LOGICAL :: FOUND_TBC_NML   !flag to determine if the timebathycontrol namelist was present
      NAMELIST /TimeBathyControl/ NDDT,BTIMINC,BCHGTIMINC

c.... !tcm v52.30.01 added for forced use of Smagorinski Eddy Viscosity
      REAL(8) :: SMAG_UPPER_LIM = 100.d0  ! m^2/s (default upper limit)
      REAL(8) :: SMAG_LOWER_LIM = 1.0d-8  ! m^2/s (default lower limit)
      LOGICAL :: SMAG_COMP_FLAG = .FALSE.   !BY DEFAULT DO NOT USE SMAG
      LOGICAL :: FOUND_Smag_NML  ! flag to determine if the Smag_Control namelist was present
      NAMELIST /Smag_Control/ SMAG_UPPER_LIM,SMAG_LOWER_LIM,SMAG_COMP_FLAG


      ! tcm v50.79 added for metControl namelist
      LOGICAL :: FOUND_METCON_NML  !flag to determine if the metControl namelist was present

      ! zc - added for tvw namelist
      LOGICAL :: FOUND_TVW_NML

Casey 121019: Added multiplication factor to be used before sending winds to coupled wave models.
      LOGICAL :: FOUND_WC_NML
#ifdef ADCNETCDF
Casey 180318: Added NWS=13
      LOGICAL :: FOUND_OWINC_NML
#endif

      ! DW 2024
      LOGICAL :: found_NWS14Grib2NetCdf_Namelist = .false. 
      NAMELIST /WindGrib2NetCdf/ read_NWS14_NetCdf_using_core_0

      logical :: foundWetDryControlNamelist = .false. ! for namelist 
      logical :: outputNodeCode = .false.
      logical :: outputNOFF = .false.
      logical :: noffActive = .true.
      logical :: StatPartWetFix = .false.  ! tcm v53.01.02
      integer :: How2FixStatPartWet = 0    ! tcm v52.01.02 (Can also be 1 in v54.)

!JLW: add subgrid namelist control
      logical :: foundSubgridControlNamelist = .false.
      character(200) :: subgridFilename
      logical :: level0 = .true. 
      logical :: level1 = .false.

      ! jgf52.08.02: For inundation output control. Wanted to just
      ! use the variables of the same name from the write_output module,
      ! but this must wait until build dependencies are resolved. 
      logical :: foundInundationOutputControlNamelist = .false.
      logical :: inundationOutput = .false.
      real(8) :: inunThresh = 0.6d0

      ! jgf: added support for water level correction
      logical :: foundCorrectionControlNamelist = .false. 
      character(len=2048) :: dynamicWaterLevelCorrectionFileName='null' ! water level correction file 
      integer :: dynamicWaterLevelCorrectionSkipSnaps=0 ! number of correction data sets to skip 
      real(8) :: dynamicWaterLevelCorrectionMultiplier = 1.d0 ! multiply each correction value by this
      real(8) :: dynamicWaterLevelCorrectionRampStart = 0.d0 ! seconds after run start
      real(8) :: dynamicWaterLevelCorrectionRampEnd = 0.d0   ! seconds after run start     
      character(len=10) :: dynamicWaterLevelCorrectionRampReferenceTime = 'hotstart'

C     VEW1D channel namelist 11/06/2023 sb
      logical :: foundVEW1DChannelControlNamelist = .false. 
      logical :: activateVEW1DChannelWetPerimeter

C     Full luni-solar tidal potential namelist: 2024
      logical:: foundTidalPotentialControlNamelist = .false. 
      LOGICAL:: UseFullTIPFormula = .FALSE. 
      INTEGER:: TIPOrder = 2
      CHARACTER (LEN=24):: TIPStartDate = 'Basedate' 
      CHARACTER (LEN=24):: MoonSunPositionComputeMethod = 'JM' 
      CHARACTER (LEN=80):: MoonSunCoordFile = 'none' 
      LOGICAL:: IncludeNutation = .TRUE.  
      REAL (8):: k2value = 0.302D0 
      REAL (8):: h2value = 0.609D0    
      NAMELIST /TidalPotentialControl/ UseFullTIPFormula, TIPOrder,
     &    TIPStartDate, MoonSunPositionComputeMethod, MoonSunCoordFile,
     &    IncludeNutation, k2Value, h2Value
      
C     SOLVER  section
      INTEGER ITITER,ISLDIA,ITMAX
      REAL(8) CONVCR
C
C     jgf48.03 netCDF metadata
      CHARACTER(80) title
      CHARACTER(80) institution
      CHARACTER(80) source
      CHARACTER(80) history
      CHARACTER(80) references
      CHARACTER(80) comments
      CHARACTER(80) host
      CHARACTER(80) convention
      CHARACTER(80) contact
      CHARACTER(80) base_date
      INTEGER NT
      INTEGER NTRSPE, NTRSPV, NTRSPM, NDSETSE, NDSETSV, NDSETSM
      INTEGER NDSETSW
      INTEGER NHY
C
C
C--------------------------------------------------------------------------C
C                                                                          C
C              DATA DECOMPOSITION DECLARATIONS BEGIN HERE                  C
C                                                                          C
C--------------------------------------------------------------------------C
C
C--Local Map Variable Declarations
C
C     LOCALI section
      INTEGER NPROC,MSHAR,NWEIR
      INTEGER,ALLOCATABLE::NELP(:), NNODP(:), NNEP(:,:,:)
      INTEGER,ALLOCATABLE::NOD_RES_TOT(:),NSTACP(:), NSTAMP(:)
      INTEGER,ALLOCATABLE::NWEIRP(:),NSTAEP(:),NSTAVP(:)
      INTEGER,ALLOCATABLE::NOPEP(:),NETAP(:),NVDLLP(:)
      INTEGER,ALLOCATABLE::NBDVP(:,:)
      INTEGER,ALLOCATABLE::NBOUP(:)    ! NBOUP(PE) Number of Land Boundary
C                                      ! Segments on PE
      INTEGER,ALLOCATABLE::NVELLP(:)   ! NVELLP(K,PE) Number of Land Boundary
C                                      ! Nodes of Segment K on PE
      INTEGER,ALLOCATABLE::NVELP(:)    ! NVELP(PE) Total Number of Land
C                                      ! Boundary Nodes on PE
      INTEGER,ALLOCATABLE::NBVVP(:,:)  ! NBVVP(K,I) Local Node Number of I-th
C                                      ! Node on Land Boundary Segment K on PE
      INTEGER,ALLOCATABLE::IBTYPEP(:,:)! IBTYPEP(K,PE) Type Land Boundary
C                                      ! Segment K on PE (0=mainland,1=island)
      INTEGER,ALLOCATABLE::IBTYPEEP(:,:)! IBTYPEEP(K,PE) Type Elev Boundary
C                                      ! Segment K on PE
      INTEGER,ALLOCATABLE::LBCODEP(:,:)! LBCODEP(I,PE) Boundary Type of Land
C                                      ! Boundary Node I on PE
      INTEGER,ALLOCATABLE::NFLUXFP(:)
      INTEGER,ALLOCATABLE::IBCONNRP(:,:)

      INTEGER,ALLOCATABLE::NNSTA3DDP(:)! Number of 3D density stations on (PE)
      INTEGER,ALLOCATABLE::NNSTA3DVP(:)! Number of 3D velocity stations on (PE)
      INTEGER,ALLOCATABLE::NNSTA3DTP(:)! Number of 3D turbulence stationson(PE)
C
C     DIAGS section
      REAL(8),ALLOCATABLE ::   PROC_SV(:)
C
C--Local-to-Global Mapping Variables
C
C     LOC2G section
      INTEGER,ALLOCATABLE  ::   IMAP_NOD_LG(:,:)
      INTEGER,ALLOCATABLE  ::   IMAP_EL_LG(:,:),   EL_SHARE(:)
      INTEGER,ALLOCATABLE  ::   IMAP_STAE_LG(:,:), STAE_SHARE(:)
      INTEGER,ALLOCATABLE  ::   IMAP_STAV_LG(:,:), STAV_SHARE(:)
      INTEGER,ALLOCATABLE  ::   IMAP_STAC_LG(:,:), STAC_SHARE(:)
      INTEGER,ALLOCATABLE  ::   IMAP_STAM_LG(:,:), STAM_SHARE(:)
      INTEGER,ALLOCATABLE  ::   OBNODE_LG(:,:)
      INTEGER,ALLOCATABLE  ::   WEIRP_LG(:,:)
      INTEGER,ALLOCATABLE  ::   WEIRDP_LG(:,:)
      INTEGER,ALLOCATABLE  ::   LBINDEX_LG(:,:)
C
C     jgf45.11 Subdomain->Fulldomain element mappings for 3D recording stations
      INTEGER,ALLOCATABLE::IMAP_STA3DD_LG(:,:) ! density
      INTEGER,ALLOCATABLE::IMAP_STA3DV_LG(:,:) ! velocity
      INTEGER,ALLOCATABLE::IMAP_STA3DT_LG(:,:) ! turbulence
C
C--Global-to-Local Mapping Variables
C
C     GLOB2L section
      INTEGER,ALLOCATABLE::ITOTPROC(:)
      INTEGER,ALLOCATABLE  ::   IMAP_NOD_GL(:,:),IMAP_NOD_GL2(:,:)
C
C--Message-Passing Variables
C
C     MSGTAB section
      INTEGER,ALLOCATABLE  ::   NUM_COMM_PE(:), COMM_PE_NUM(:,:)
      INTEGER,ALLOCATABLE  ::   IRECV_TOT(:,:),IRECV(:)
      INTEGER,ALLOCATABLE  ::   ISEND_TOT(:,:),ISEND(:)

      integer, save       :: FileFmtVersion      ! File format version.

C...tcm v49.48.01 -- Added for new point in element searches
C-- KDTREE2 ENABLED SEARCH ALGORITHM VARIABLES
C
      INTEGER :: SRCHDP=100
      REAL(8), ALLOCATABLE ::  RMAX(:),BCXY(:,:)

CCCCCCCCCCCCCCCCCC DW
      LOGICAL:: prep5354flag = .FALSE.
      LOGICAL:: prep2001flag = .FALSE.
      INTEGER,ALLOCATABLE :: nodeNumSP(:)
C-------------------end of data declarations------------------------C


      CONTAINS


      function VERSION_NUMBER(Major, Minor, Rev) result(vn)
      implicit none
      integer   :: vn, major, minor, rev

      vn = ior(ior(ishft(major,20),ishft(minor,10)),rev)
      return
      end function


      function CMP_VERSION_NUMBERS(a,b) result(match)
      implicit none
      integer  a,b
      logical  :: match

      match = ishft(a,-10) == ishft(b,-10)
      return
      end function


c----------------------------------------------------------------------
      SUBROUTINE ALLOC_MAIN1()
c----------------------------------------------------------------------
      use memory_usage
      integer(8) :: nbytes = 0
C
C     Allocate space for Arrays except those dimensioned by MNPP and MNEP
C
      ALLOCATE ( PROC(MNP) )
      nbytes = nbytes + 4*mnp
      ALLOCATE ( X(MNP),Y(MNP),DP(MNP),SLAM(MNP),SFEA(MNP) )
      nbytes = nbytes + 40*mnp
      ALLOCATE ( SL1(MNP),SF1(MNP))
      nbytes = nbytes + 16*mnp
      ALLOCATE ( NNEG(3,MNE) )
      nbytes = nbytes + 12*mne
      ALLOCATE ( NVDLL(MNOPE),NBDV(MNOPE,MNETA),IBTYPEE(MNOPE) )
      nbytes = nbytes + 4*(mneta+1)*mnope + 4*mnope
      ALLOCATE ( NVELL(MNBOU),NBVV(MNBOU,0:MNVEL),IBTYPE(MNBOU) )
      nbytes = nbytes + 4*mnbou + 4*(mnvel+2)*mnbou
      ALLOCATE ( IBCONNR(MNBOU,MNVEL),LBCODE(MNVEL) )
      nbytes = nbytes + 4*(mnbou+1)*mnvel
      ALLOCATE ( WEIR(MNVEL),WEIRD(MNVEL) )
      nbytes = nbytes + 8*mnvel
      ALLOCATE ( BAR1(MNBOU,MNVEL),BAR2(MNBOU,MNVEL),BAR3(MNBOU,MNVEL) )
      nbytes = nbytes + 12*mnbou*mnvel
      if (meshHasWeirWithPipes.eqv..true.) then
         ALLOCATE ( BAR4(MNBOU,MNVEL),BAR5(MNBOU,MNVEL),BAR6(MNBOU,MNVEL) )
         nbytes = nbytes + 12*mnbou*mnvel
      endif
      ALLOCATE ( FLBN(MNVEL),FLBNX(MNVEL),FLBNXP(MNVEL) )
      nbytes = nbytes + 12*mnvel
      ALLOCATE ( NVDLLMSG(MNOPE),NVELLMSG(MNBOU+1) )
      nbytes = nbytes + 4*mnope+4*(mnbou+1)
      ALLOCATE ( NNSEG(MNSTAE),NNSVG(MNSTAV),NNSCG(MNSTAC),
     &           NNSMG(MNSTAM) )
      ALLOCATE ( STAELOC(MNSTAE),STAVLOC(MNSTAV),STACLOC(MNSTAC) )
      ALLOCATE ( STAMLOC(MNSTAM) )
      nbytes = nbytes + 8*(mnstae+mnstav+mnstam+mnstac)
      ALLOCATE ( TIPOTAG(MNTIF),BOUNTAG(MNBFR),FBOUNTAG(MNFFR) )
      nbytes = nbytes + 4*(mntif+mnbfr+mnffr)
      ALLOCATE ( ALPHA1(MNBFR),ALPHA2(MNFFR),FREQMSG(MNFFR),
     &           QNMSG(MNFFR,MNVEL) )
      nbytes = nbytes + 4*mnbfr+8*mnffr+4*mnffr*mnvel
      ALLOCATE ( AMIGMSG(MNBFR),EMOMSG(MNBFR,MNETA),TPKMSG(MNTIF) )
      nbytes = nbytes + 4*mnbfr+4*mntif+4*mnbfr*mneta
      ALLOCATE ( XEL(MNSTAE),YEL(MNSTAE),SLEL(MNSTAE),SFEL(MNSTAE) )
      nbytes = nbytes + 12*mnstae
      ALLOCATE ( XEV(MNSTAV),YEV(MNSTAV),SLEV(MNSTAV),SFEV(MNSTAV) )
      nbytes = nbytes + 12*mnstav
      ALLOCATE ( XEC(MNSTAC),YEC(MNSTAC),SLEC(MNSTAC),SFEC(MNSTAC) )
      nbytes = nbytes + 12*mnstac
      ALLOCATE ( XEM(MNSTAM),YEM(MNSTAM),SLEM(MNSTAM),SFEM(MNSTAM) )
      nbytes = nbytes + 12*mnstam
      ALLOCATE ( FF(MNBFR),FACE(MNBFR) )
      nbytes = nbytes + 8*mnbfr
      ALLOCATE ( EMO(MNBFR,MNETA),EFA(MNBFR,MNETA) )
      nbytes = nbytes + 8*mnbfr*mneta
      ALLOCATE ( TPK(MNTIF),ETRF(MNTIF),FFT(MNTIF),FACET(MNTIF) )
      nbytes = nbytes + 16*mntif
      ALLOCATE ( QNAM(MNFFR,MNVEL),QNPH(MNFFR,MNVEL) )
      nbytes = nbytes + 8*mnffr*mnvel
      ALLOCATE ( FFF(MNFFR),FFACE(MNFFR) )
      ALLOCATE ( AMIG(MNBFR),AMIGT(MNTIF), FAMIG(MNFFR) )
      nbytes = nbytes + 8*mnbfr+8*mntif+8*mnffr
      ALLOCATE ( NELP(MNPROC), NNODP(MNPROC) )
      ALLOCATE ( NOD_RES_TOT(MNPROC) )
      ALLOCATE ( NWEIRP(MNPROC) )
      ALLOCATE ( NSTAEP(MNPROC),NSTAVP(MNPROC),NSTACP(MNPROC) )
      ALLOCATE ( NSTAMP(MNPROC) )
      nbytes = nbytes + 32*mnproc
      ALLOCATE ( NOPEP(MNPROC),NETAP(MNPROC),NVDLLP(MNOPE),
     &           IBTYPEEP(MNOPE,MNPROC))
      nbytes = nbytes + 8*mnproc + 4*mnope + 4*mnope*mnproc
      ALLOCATE ( NBDVP(MNOPE,MNETA) )
      nbytes = nbytes + 4*mnope*mneta
      ALLOCATE ( NBOUP(MNPROC),NVELP(MNPROC),NVELLP(MNBOU) )
      nbytes = nbytes + 8*mnproc+4*mnbou
      ALLOCATE ( NBVVP(MNBOU,0:MNVEL), IBTYPEP(MNBOU,MNPROC) )
      nbytes = nbytes + 4*mnbou*(mnvel+1)+4*mnproc*mnbou
      ALLOCATE ( LBCODEP(MNVEL,MNPROC) )
      nbytes = nbytes + 4*mnvel*mnproc
      ALLOCATE ( ITOTPROC(MNP),NFLUXFP(MNPROC) )
      nbytes = nbytes + 4*(mnp+mnproc)
      ALLOCATE ( IBCONNRP(MNBOU,MNVEL) )
      nbytes = nbytes + 4*(mnbou*mnvel)
      ALLOCATE ( PROC_SV(MNPROC) )
      nbytes = nbytes + 4*mnproc
      ALLOCATE ( IMAP_STAE_LG(MNSTAE,MNPROC) )
      ALLOCATE ( IMAP_STAV_LG(MNSTAV,MNPROC) )
      ALLOCATE ( IMAP_STAC_LG(MNSTAC,MNPROC) )
      ALLOCATE ( IMAP_STAM_LG(MNSTAM,MNPROC) )
      nbytes = nbytes + 4*mnproc*(mnstae+mnstav+mnstac+mnstam)
      ALLOCATE ( OBNODE_LG(MNETA,MNPROC) )
      nbytes = nbytes + 4*mneta*mnproc
      ALLOCATE ( WEIRP_LG(MNVEL,MNPROC) )
      ALLOCATE ( WEIRDP_LG(MNVEL,MNPROC) )
      nbytes = nbytes + 4*mnvel*mnproc
      ALLOCATE ( LBINDEX_LG(MNBOU,MNVEL) )
      nbytes = nbytes + 4*mnbou*mnvel
      ALLOCATE ( IMAP_NOD_GL(2,MNP),IMAP_NOD_GL2(2*MNEI,MNP) )
      nbytes = nbytes + 8*mnp + 8*mnei*mnp
      ALLOCATE ( NUM_COMM_PE(MNPROC), COMM_PE_NUM(MNPROC,MNPROC) )
      nbytes = nbytes + 4*(mnproc+mnproc*mnproc)
      ALLOCATE ( IRECV_TOT(MNPROC,MNPROC),IRECV(MNP) )
      ALLOCATE ( ISEND_TOT(MNPROC,MNPROC),ISEND(MNP) )
      nbytes = nbytes + 8*(mnp + mnproc*mnproc)
      ALLOCATE ( STAE_SHARE(MNSTAE), STAV_SHARE(MNSTAV),
     &    STAM_SHARE(MNSTAM))
      ALLOCATE ( STAC_SHARE(MNSTAC))
      nbytes = nbytes + 4*(mnstae+mnstav+mnstac+mnstam)
      STAE_SHARE(:) = -1
      STAV_SHARE(:) = -1
      STAM_SHARE(:) = -1
      STAC_SHARE(:) = -1

C....tcm 44.48.01 allocate memory for fast Search Tree
      allocate (rmax(mne),bcxy(2,mne))
      nbytes = nbytes + 8*mne + 2*8*mne

      call memory_alloc(nbytes)
      print *,"from alloc_main1: "
      call memory_status()
C
      RETURN
C---------------------------------------------------------------------
      END SUBROUTINE ALLOC_MAIN1
C---------------------------------------------------------------------
C
C
C---------------------------------------------------------------------
C---------------------------------------------------------------------
      SUBROUTINE ALLOC_MAIN2()
      use memory_usage
      integer(8) :: nbytes
C
C     Allocate space for Arrays dimensioned by MNPP and MNEP
C
      nbytes = 0
      ALLOCATE ( IMAP_NOD_LG(MNPP,MNPROC),IMAP_EL_LG(MNEP,MNPROC) )
      ALLOCATE ( NNEP(3,MNEP,MNPROC) )
      nbytes = nbytes + 4*mnpp*mnproc + 4*mnep*mnproc + 12*mnep*mnproc
      ALLOCATE(EL_SHARE(NELG))
      EL_SHARE(:)   = -1
      
      nbytes = nbytes + mnep
      call memory_alloc(nbytes)
      print *,"from alloc_main2: "
      call memory_status()
C
      RETURN
C---------------------------------------------------------------------
      END SUBROUTINE ALLOC_MAIN2
C---------------------------------------------------------------------
C
C
C---------------------------------------------------------------------
C---------------------------------------------------------------------
      subroutine dealloc_main2()
      use memory_usage
      integer(8) :: nbytes
C
C     Deallocate space for Arrays dimensioned by MNPP and MNEP
C
      nbytes = 0
      deallocate(imap_nod_lg,imap_el_lg)
      deallocate(nnep)
      nbytes = nbytes + 4*mnpp*mnproc + 4*mnep*mnproc + 12*mnep*mnproc
      deallocate(el_share)
      
      nbytes = nbytes + mnep
      call memory_dealloc(nbytes)
      print *,"from dealloc_main2: "
      call memory_status()
c
      return
c---------------------------------------------------------------------
      end subroutine dealloc_main2
C---------------------------------------------------------------------
C
C
C
C--------------------------------------------------------------------------C
C                                                                          C
C                DEFINITIONS OF DOMAIN DECOMPOSITION VARIABLES             C
C                                                                          C
C--------------------------------------------------------------------------C
C
C Processing Element Definitions:
C
C   MNPROC               = Maximum Number of PEs
C   NPROC                = Actual Number of PEs - this should eventually be
C                          dropped and MNPROC used throughout preprocessor
C                          routines - RL
C   MSHAR                = Max. Number PEs assigned to any Global Node
C
C Nodal and Element Definitions:
C
C   X(I)                 = X-coordinate of Global Node I
C   Y(I)                 = Y-coordinate of Global Node I
C   DP(I)                = Bathymetry of Global Node I
C
C   NELG                 = Number of Global Elements
C   NELP(PE)             = Number of Elements Assigned to PE
C
C   NNODG                = Number of Global Nodes
C   NNODP(PE)            = Number of Nodes Assigned to PE
C
C   NOD_RES_TOT(PE)      = Number of Resident Nodes on PE
C   ITOTPROC(I)          = Number of PEs assigned to Global Node I
C
C   NNEG(3,I)            = Three Nodes of Global Element I
C   NNEP(3,I,PE)         = Three Nodes of Element I on PE
C
C   IMAP_NOD_GL(1,I)     = PE assigned to Global Node I
C   IMAP_NOD_GL(2,I)     = Local Node Number of Global Node I
C
C   IMAP_NOD_LG(I,PE)    = Global Node Number of Local Node I on PE
C   IMAP_EL_LG(I,PE)     = Global Element Number of Local Element I on PE
C
C   IMAP_NOD_GL2(2(PE-1)+1,I)  = PE assigned to Global Node I
C   IMAP_NOD_GL2(2(PE-1)+2,I)  = Local Node Number of Global Node I on PE
C
C Open Boundary Nodes and Segment Definitions:
C
C   NETA                 = Number of Global Open Boundary Nodes
C   NETAP(PE)            = Number of Open Boundary Nodes on PE
C
C   NOPE                 = Number of Global Open Boundary Segments
C   NOPEP(PE)            = Number of Open Boundary Segments on PE
C
C   NVDLL(K)             = Number of Nodes on Global Open Boundary Segment K
C   NVDLLP(K)            = Number of Nodes on Open Boundary Segment K on PE
C
C   NBDV(K,I)            = Global Node Number of I-th Node on Open Boundary
C                          Segment K
C   NBDVP(K,I)           = Local Node Number of I-th Node on Open Boundary
C                          Segment K on PE
C
C   OBNODE_LG(I,PE)      = Global Open Boundary Node Number of Local
C                          Open Boundary Node I on PE
C
Cvjp modified array to drop last dimension to save memory space
C   LBINDEX_LG(K,I,PE)   = Global Index of I-th Node on Land Boundary Segment
C                          K on PE
C
C   NWEIR                = Total Number of Global Weir Land Boundary Pairs
C   NWEIRP(PE)           = Total Number of Land Boundary Nodes on PE
C
C   WEIRP_LG(I,PE)       = Global Node Number of I-th Weir Node on PE
C   WEIRDP_LG(I,PE)      = Global Node Number of I-th Dual Weir Node on PE
C
C Elevation Station Definitions:
C
C   NSTAEP(PE)           = Number of Elevation Stations on PE
C   IMAP_STAE_LG(I,PE)   = Global Number of Local Elevation Station I on PE
C
C Velocity Station Definitions:
C
C   NSTAVP(PE)           = Number of Velocity Stations on PE
C   IMAP_STAV_LG(I,PE)   = Global Number of Local Velocity  Station I on PE
C
C Concentration Station Definitions:
C
C   NSTACP(PE)           = Number of Concentration Stations on PE
C   IMAP_STAC_LG(I,PE)   = Global Number of Local Concentration Station I on PE
C
C Meterological Station Definitions:
C
C   NSTAMP(PE)           = Number of Meterological Stations on PE
C   IMAP_STAM_LG(I,PE)   = Global Number of Local Meterological Station I on PE
C
C
C Message-Passing Definitions:
C
C   NUM_COMM_PE(PE)      = Number of PEs communicating with PE
C   COMM_PE_NUM(IPE,PE)  = IPE-th PE communicating with PE
C
C   IRECV_TOT(IPE,PE)    = Number of Nodes Received by PE from IPE
C   ISEND_TOT(IPE,PE)    = Number of Nodes Sent by PE to IPE
C
C   PROC_SV(PE)          = Surface-to-Volume Ratio on PE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      END MODULE PRE_GLOBAL
