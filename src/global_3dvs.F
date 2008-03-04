C******************************************************************************
C PADCIRC VERSION 45.12 03/17/2006                                            *
C  last changes in this file VERSION 45.11                                    *
C                                                                             *
C This module declares the 3D exclusive global variables.  It also uses some  *
C general global variables that were declared in module GLOBAL.               *
C The module also allocates the required arrays.                              *
C                                                                             *
C******************************************************************************
C
      MODULE GLOBAL_3DVS
C...
C...BRING IN NECESSARY VARIABLES FROM GLOBAL
C...
      USE SIZES,
     &     ONLY : SZ, !size of reals (4 or 8 bytes), set during compile time
     &     NBYTE,
     &     MNP,
     &     MNE,
     &     MNEI,
     &     MNVEL,
     &     MNFEN, 
     &     MYProc,
     &     INPUTDIR,
     &     LOCALDIR

      USE GLOBAL,
     &     ONLY :
     &       X, Y,                         !X(MNP),Y(MNP) = nodal coordinates
     &       DP,                           !DP(MNP) = still water depth
     &       NM,                           !NM(MNE,3) = nodal table for each element
     &       NNeigh,                       !NNeigh(MNP) = number of neighboring nodes for each node
     &       NEITAB,                       !NeiTab(MNP,MNEI) = table of neighbor nodes for each node
     &       NEITABELE,                    !NeiTabEle(MNP,MNei) = table of neighboring elements for each node 
     &       AREAS,                        !AREAS(MNE) = 2*Element Area
     &       WSX1, WSY1,                   !WSX1(MNP),WSY1(MNP) = Wind stress components at time level s  
     &       WSX2, WSY2,                   !WSX2(MNP),WSY2(MNP) = Wind stress components at time level s+1
     &       BSX => BSX1, BSY => BSY1,     !BSX(MNP),BSY(MNP) = bottom stress computed after velocity solution
     &       UU => UU2, VV => VV2,         !UU(MNP),VV(MNP) = vertically averaged velocity components 
     &       DAFluxX => QX2, DAFluxY => QY2,
     &       DUU => DUU1, DUV => DUV1, DVV => DVV1,   !DUU(MNP),DVV(MNP),DUV(MNP) = velocity dispersion terms
     &       DASigT,
     &       ETA1, ETA2,                   !ETA1(MNP),ETA2(MNP) = water surf elev at time levels s, s+1
     &       CORIF,                        !CORIF(MNP) = nodal values of Coriolis parameter
     &       BTP => MOM_LV_X,              !BTP(MNP) = barotropic pressure (incl TP & wl) at time levels s+1/2
     &       LBARRAY_POINTER,              !LBARRAY_POINTER(MNP) = pointer into array of land or flux boundary
     &       LBCODEI,                      !LBCODEI(MNVel) = array of land or flux boundary codes
     &       CSII, SIII,                   !CSII(MNVel),SIII(MNVel) = cosine and sine of normal flux boundary node
     &       QNORMSP1 => QN2,              !QNORMSP1(MNVel) = specified normal flux boundary condition at time level s+1
     &       ICS,
     &       NP,                           !number of horizontal nodes
     &       NE,                           !number of horizontal elements
     &       IFNLCT,                       !nonlinear advection flag (1=yes, 0=no) (jgf45.11)
     &       IFNLFA,                       !nonlinear finite amplitude flag (1=yes,0=no) (jgf45.11)
     &       NScreen,                      !flag to suppress or allow screen output
     &       ScreenUnit,                   !i/o for screen output
     &       IHOT,                         !flag to specify if coldstart (=0) or hotstart file to read from (63 or 64)
     &       HSS,                          !hotstart output file structure
     &       IHOTSTP,                      !hotstart file record counter
     &       NWS,
     &       IDen,                         !flag telling whether salinity, temperature and/or sigmaT are being used
     &       RhoWat0, SigT0, G,
     &       RUNDES, RDES4, RDES8,
     &       RUNID, RID4, RID8,
     &       AGRID, AID4, AID8,
     &       CBaroclinic,
     &       C3DVS,                        
     &       DelT => DTDP,                 !DelT = altername for timestep real(8)
     &       DTDP,                         !DTDP = timestep declared as real(8)
     &       SLAM0, SFEA0,                 !center point of CPP spherical projection
     &       Deg2Rad, Rad2Deg,
     &       Day2Sec, Sec2Day


C...DECLARE 3D GLOBAL ARRAYS

C     Gamma(NFEN): horizontal velocity solution in the complex form u +
C     iv for a specific node
      COMPLEX(SZ),ALLOCATABLE  :: Gamma(:)
C     Q(MNP,NFEN): horizontal velocity solution in the complex form u +
C     iv for all nodes
      COMPLEX(SZ),ALLOCATABLE  :: Q(:,:)     
C     Sigma(NFEN): stretched vertical coordinate levels (-1 to 1)
      REAL(SZ),ALLOCATABLE :: SIGMA(:)   
C     EVTot(NFEN): vertical eddy viscosity values at a node
      REAL(SZ),ALLOCATABLE ::	EVTot(:) 
C     INM(NFEN,3): Integral used in vertical FE assembly
      REAL(SZ),ALLOCATABLE :: INM(:,:)
C     LVN(NFEN): Integral used in vertical FE assembly
      REAL(SZ),ALLOCATABLE ::	LVN(:)
C     WZ(MNP,NFEN): "z" vertical velocity, all nodes
      REAL(SZ),ALLOCATABLE :: WZ(:,:)   
C     SIGT(MNP,NFEN): sigma T, all nodes
      REAL(SZ),ALLOCATABLE :: SIGT(:,:)
C     TEMP(MNP,NFEN): temperature, all nodes
      REAL(SZ),ALLOCATABLE :: Temp(:,:)  
C     Sal(MNP,NFEN): salinity, all nodes
      REAL(SZ),ALLOCATABLE :: Sal(:,:)  
C     BCP(MNP,NFEN): baroclinic pressure, integrated down from surface,
C     all nodes
      REAL(SZ),ALLOCATABLE :: BCP(:,:)
C     Q20(MNP,NFEN): Turbulent Kinetic Energy computed by Mellor-Yamada
C     L2.5 closure, all nodes
      REAL(SZ),ALLOCATABLE :: Q20(:,:)
C     L(MNP,NFEN): turbulent length scale computed by Mellor-Yamada L2.5
C     closure, all nodes
      REAL(SZ),ALLOCATABLE :: L(:,:)
C     EV(MNP,NFEN): Vertical eddy viscosity, all nodes
      REAL(SZ),ALLOCATABLE :: EV(:,:) 

C     kmd45.12 Added in variables for transport 

C     Gammatrans(NFEN): horizontal transport soln for a specific node
      REAL(SZ),ALLOCATABLE :: Gammatrans(:)
C     Tempkp1(MNP,NFEN): temperature of all nodes at future time level
      REAL(SZ),ALLOCATABLE :: Tempkp1(:,:)  
C     Salkp1(MNP,NFEN): salinity of all nodes at future time level 
      REAL(SZ),ALLOCATABLE :: Salkp1(:,:)
     
C     kmd45.12 Added in these variables for the top boundary conditions
C     associated with the temperature field
      REAL(SZ), ALLOCATABLE :: qsurfkp1(:), qsurf(:)

C     kmd45.12 Added in the variables for saving 3D velocities for two
C     time levels

C     Qkp1(MNP,NFEN): horizontal velocity in the complex form u + iv for
C     all nodes at the future time level
      COMPLEX(SZ),ALLOCATABLE  :: Qkp1(:,:)
C     WZkp1(MNP,NFEN): "z" vertical velocity at all nodes for the future
C     time level
      REAL(SZ),ALLOCATABLE :: WZkp1(:,:)

!     kmd45.12 added in for baroclinic changes
      COMPLEX(SZ), ALLOCATABLE :: BPG(:,:)      


C...  DECLARE 3D GLOBAL SCALARS

      REAL(SZ) :: GORho                !=G/RhoWat0
      REAL(SZ) :: GORhoOAMB            !=GORho/AMB
      REAL(SZ) :: KP                   !3D bottom friction coefficient
      REAL(SZ) :: EVMin, EVCon         !minimum vertical eddy viscosity and vertical eddy viscosity constant
      REAL(SZ) :: Z0S, Z0B             !surface and bottom roughnesses
      REAL(SZ) :: DTALP2, DT1MALP2     !DelT*alpha2, DelT*(1-alpha2)
      REAL(SZ) :: DTALP3, DT1MALP3     !DelT*alpha3, DelT*(1-alpha3)
      REAL(SZ) :: THETA1, THETA2       !timestepping coefficients for MY L2.5 TKE calculations
      REAL(SZ), PARAMETER ::  A= 1.D0  !Value of Sigma at the surface
      REAL(SZ), PARAMETER ::  B=-1.D0  !Value of Sigma at the bottom
      REAL(SZ), PARAMETER ::  AMB=A-B

      INTEGER  :: NFEN                 !Number of nodes in the vertical
      INTEGER  :: ISLIP                !bottom slip coefficient flag (0=no slip, 1=linear, 2=quadratic)
      INTEGER  :: IEVC                 !eddy viscosity flag determine type of EV used
      INTEGER  :: IStart3D = 0         !=0 when code starts up, incremented each time VSSOL is called, used to initialize EV

      COMPLEX(SZ) :: IDTALP1              !i*DelT*alpha1 
      COMPLEX(SZ) :: IDT1MALP1            !i*DelT*(1-alpha1)
      COMPLEX(SZ),  PARAMETER ::  I=(0.0d0, 1.0d0)   !imaginary i

C...Terms added for the transport parameters !kmd45.11 transport
      REAL(SZ) :: NLSD, NVSD           !diffusion coefficients for the salinity field
      REAL(SZ) :: NLTD, NVTD           !diffusion coefficients for the temperature field 
      REAL(SZ) :: Alp4                 !time weighting coefficient for the vertical diffusion term in the transport equation
      INTEGER :: NTF                  !Flag for the top boundary condition for the temperature field in the transport equation
      REAL(SZ) :: DTALP4, DT1MALP4     !DelT*alpha4, DelT*(1-alpha4)

C...Declare and initialize logical variables

      LOGICAL  :: C3D_BTrans     = .FALSE.      !if true, then 3D prognostic baroclinic run
      LOGICAL  :: C3D_PTrans     = .FALSE.      !if true, then 3D passive scalar transport included
      LOGICAL  :: turb_allocated = .FALSE.


C...3D Density station output variables

      INTEGER :: I3DSD        !flag for whether to have 3D station density output and format of this output
      INTEGER :: NSpo3DSD     !3D station density output interval in timesteps
      INTEGER :: NTO3DSDS     !TO3DSDS converted to timesteps after StaTime
      INTEGER :: NTO3DSDF     !TO3DSDF converted to timesteps after StaTime
      INTEGER :: NDSet3DSD    !total number of times that 3D density station data will be written
      INTEGER :: N3DSD        !counter to determine if it is time to write 3D density station data
      INTEGER :: I3DSDRec     !record counter for binary file writes
      INTEGER :: NSta3DD      !number of 3D density stations
      INTEGER :: NSta3DD_G    !jgf47.06 global number of 3D density stations
      INTEGER :: MNSta3DD     !larger of 1 or NSta3DD, used to dimension arrays
      INTEGER, ALLOCATABLE :: NE3DD(:)     !element number containing 3D density station N
      REAL(8),ALLOCATABLE :: StaI3DD1(:)   !interpolating factor used to compute output at 3D density station N
      REAL(8),ALLOCATABLE :: StaI3DD2(:)   !interpolating factor used to compute output at 3D density station N
      REAL(8),ALLOCATABLE :: StaI3DD3(:)   !interpolating factor used to compute output at 3D density station N
      REAL(SZ),ALLOCATABLE :: SigTSta(:)   !3D Sigma T station output - not a global variable but easier to allocate here 
      REAL(SZ),ALLOCATABLE :: TempSta(:)   !3D Temperature station output - not a global variable but easier to allocate here
      REAL(SZ),ALLOCATABLE :: SalSta(:)    !3D Salinity station output - not a global variable but easier to allocate here

C...3D Velocity station output variables

      INTEGER :: I3DSV        !flag for whether to have 3D station velocity output and format of this output
      INTEGER :: NSpo3DSV     !3D station velocity output interval in timesteps
      INTEGER :: NTO3DSVS     !TO3DSVS converted to timesteps after StaTime
      INTEGER :: NTO3DSVF     !TO3DSVF converted to timesteps after StaTime
      INTEGER :: NDSet3DSV    !total number of times that 3D velocity station data will be written
      INTEGER :: N3DSV        !counter to determine if it is time to write 3D velocity station data
      INTEGER :: I3DSVRec     !record counter for binary file writes
      INTEGER :: NSta3DV      !number of 3D velocity stations
      INTEGER :: NSta3DV_G    !jgf47.06 global number of 3D velocity stations
      INTEGER :: MNSta3DV     !larger of 1 or NSta3DV, used to dimension arrays
      INTEGER, ALLOCATABLE :: NE3DV(:)     !element number containing 3D velocity station N
      REAL(8),ALLOCATABLE :: StaI3DV1(:)   !interpolating factor used to compute output at 3D velocity station N
      REAL(8),ALLOCATABLE :: StaI3DV2(:)   !interpolating factor used to compute output at 3D velocity station N
      REAL(8),ALLOCATABLE :: StaI3DV3(:)   !interpolating factor used to compute output at 3D velocity station N
      REAL(SZ),ALLOCATABLE :: WZSta(:)     !vertical velocity station output - not a global variable but easier to allocate here
      COMPLEX,ALLOCATABLE :: qSta(:)       !3D velocity station output - not a global variable but easier to allocate here

C...3D Turbulence station output variables

      INTEGER :: I3DST        !flag for whether to have 3D station turbulence output and format of this output
      INTEGER :: NSpo3DST     !3D station turbulence output interval in timesteps
      INTEGER :: NTO3DSTS     !TO3DSTS converted to timesteps after StaTime
      INTEGER :: NTO3DSTF     !TO3DSTF converted to timesteps after StaTime
      INTEGER :: NDSet3DST    !total number of times that 3D turbulence station data will be written
      INTEGER :: N3DST        !counter to determine if it is time to write 3D turbulence station data
      INTEGER :: I3DSTRec     !record counter for binary file writes
      INTEGER :: NSta3DT      !number of 3D turbulence stations
      INTEGER :: NSta3DT_G    !jgf47.06 global number of 3D turbulence stations
      INTEGER :: MNSta3DT     !larger of 1 or NSta3DT, used to dimension arrays
      INTEGER, ALLOCATABLE :: NE3DT(:)     !element number containing 3D turbulence station N
      REAL(8),ALLOCATABLE :: StaI3DT1(:)   !interpolating factor used to compute output at 3D turbulence station N
      REAL(8),ALLOCATABLE :: StaI3DT2(:)   !interpolating factor used to compute output at 3D turbulence station N
      REAL(8),ALLOCATABLE :: StaI3DT3(:)   !interpolating factor used to compute output at 3D turbulence station N
      REAL(SZ),ALLOCATABLE :: q20Sta(:)    !TKE station output - not a global variable but easier to allocate here
      REAL(SZ),ALLOCATABLE :: lSta(:)      !turb length scale station output - not a global variable but easier to allocate here
      REAL(SZ),ALLOCATABLE ::	EVSta(:)     !vertical viscosity station output - not a global variable but easier to allocate here

C...3D Global density output variables

      INTEGER :: I3DGD        !flag for whether to have 3D global density output and format of this output
      INTEGER :: NSpo3DGD     !3D global density output interval in timesteps
      INTEGER :: NTO3DGDS     !TO3DSDS converted to timesteps after StaTime
      INTEGER :: NTO3DGDF     !TO3DSDF converted to timesteps after StaTime
      INTEGER :: NDSet3DGD    !total number of times that 3D global density data will be written
      INTEGER :: N3DGD        !counter to determine if it is time to write 3D global density data
      INTEGER :: I3DGDRec     !record counter for binary file writes

C...3D Global velocity output variables

      INTEGER :: I3DGV        !flag for whether to have 3D global velocity output and format of this output
      INTEGER :: NSpo3DGV     !3D global velocity output interval in timesteps
      INTEGER :: NTO3DGVS     !TO3DGVS converted to timesteps after StaTime
      INTEGER :: NTO3DGVF     !TO3DGVF converted to timesteps after StaTime
      INTEGER :: NDSet3DGV    !total number of times that 3D global velocity data will be written
      INTEGER :: N3DGV        !counter to determine if it is time to write 3D global velocity data
      INTEGER :: I3DGVRec     !record counter for binary file writes

C...3D Global turbulence output variables

      INTEGER :: I3DGT        !flag for whether to have 3D global turbulence output and format of this output
      INTEGER :: NSpo3DGT     !3D global turbulence output interval in timesteps
      INTEGER :: NTO3DGTS     !TO3DGTS converted to timesteps after StaTime
      INTEGER :: NTO3DGTF     !TO3DGTF converted to timesteps after StaTime
      INTEGER :: NDSet3DGT    !total number of times that 3D global turbulence data will be written
      INTEGER :: N3DGT        !counter to determine if it is time to write 3D global turbulence data
      INTEGER :: I3DGTRec     !record pointer for binary file writes

C-------------------end of data declarations----------------------------------C


      CONTAINS




C
C     Allocate space for arrays used in 3D VS routines
C
      SUBROUTINE ALLOC_3DVS()
      ALLOCATE( SIGMA(NFEN), EVTOT(NFEN) )
      ALLOCATE( GAMMA(NFEN), INM(NFEN,3), LVN(NFEN) )
      ALLOCATE( Q(MNP,NFEN), WZ(MNP,NFEN) )
      ALLOCATE( SIGT(MNP,NFEN), TEMP(MNP,NFEN) )
      ALLOCATE( SAL(MNP,NFEN), BCP(MNP,NFEN) )
      ALLOCATE( Q20(MNP,NFEN), L(MNP,NFEN), EV(MNP,NFEN) )
      ALLOCATE( Gammatrans(NFEN) )
      ALLOCATE( Qkp1(MNP,NFEN), WZkp1(MNP,NFEN) )
      ALLOCATE( SALkp1(MNP,NFEN), TEMPkp1(MNP,NFEN) )
      ALLOCATE( qsurfkp1(MNP), qsurf(MNP) )
      ALLOCATE( BPG(MNP,NFEN) )
      RETURN
      END SUBROUTINE ALLOC_3DVS

C
C     Allocate space for arrays used in 3D Density station output
C
      SUBROUTINE ALLOC_3DSD()
      ALLOCATE( NE3DD(MNSta3DD) ) 
      ALLOCATE( StaI3DD1(MNSta3DD), StaI3DD2(MNSta3DD),
     &                              StaI3DD3(MNSta3DD) )
      ALLOCATE( SigTSta(NFEN), SalSta(NFEN), TempSta(NFEN) )
      RETURN
      END SUBROUTINE ALLOC_3DSD      

C
C     Allocate space for arrays used in 3D Velocity station output
C
      SUBROUTINE ALLOC_3DSV()
      ALLOCATE( NE3DV(MNSta3DV) ) 
      ALLOCATE( StaI3DV1(MNSta3DV), StaI3DV2(MNSta3DV),
     &                              StaI3DV3(MNSta3DV) )
      ALLOCATE( qSta(NFEN), WZSta(NFEN) )
      RETURN
      END SUBROUTINE ALLOC_3DSV

C
C     Allocate space for arrays used in 3D Turbulence station output
C
      SUBROUTINE ALLOC_3DST()
      ALLOCATE( NE3DT(MNSta3DT) ) 
      ALLOCATE( StaI3DT1(MNSta3DT), StaI3DT2(MNSta3DT),
     &                              StaI3DT3(MNSta3DT) )
      ALLOCATE( q20Sta(NFEN), lSta(NFEN), EVSta(NFEN) )
      RETURN
      END SUBROUTINE ALLOC_3DST


      END MODULE GLOBAL_3DVS
