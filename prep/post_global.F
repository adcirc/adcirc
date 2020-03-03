C     jgf Updated for 45.07 to add handling of fort.71, fort.72, fort.73
C         and the NOFF array in output files 
C     jgf Updated for 45.06 09/07/2005 b/c of changes in hot start file format
C     jgf Updated for 45.11 02/02/2006 b/c for 3D recording stations
C     jgf Updated for 45.12 03/17/2006 
      MODULE POST_GLOBAL
      IMPLICIT NONE
C...SET NUMBER OF BYTES "SZ" IN REAL(SZ) DECLARATIONS       

#ifdef REAL4
      INTEGER, PARAMETER :: SZ = 4
#endif
#ifdef REAL8
      INTEGER, PARAMETER :: SZ = 8
#endif

      INTEGER MNPROC,MNP,MNE,MNPP,MNSTAE,MNSTAV,MNHARF,MNWLAT,MNWLON
!$omp threadprivate(MNPROC,MNP,MNE,MNPP,MNSTAE,MNSTAV,MNHARF,MNWLAT,MNWLON)
      INTEGER MNSTAC !jgf46.02 # of conc. rec. sta. in full domain
!$omp threadprivate(MNSTAC)
      INTEGER MNSTAM ! number of met. recording stations in full domain
!$omp threadprivate(MNSTAM)
      INTEGER MNEP   !jgf45.07 Max number of elements in any subdomain
!$omp threadprivate(MNEP)
C
C   PARM14   section
      INTEGER NELG,NNODG
!$omp threadprivate(NELG,NNODG)
C
C   STRING14 section
      CHARACTER*80 AGRID
!$omp threadprivate(AGRID)
C
C   PARM15  section
      INTEGER IM
!$omp threadprivate(IM)
      INTEGER NWS
!$omp threadprivate(NWS)
      INTEGER NOUTE,NSPOOLE,NSTAE
!$omp threadprivate(NOUTE,NSPOOLE,NSTAE)
      INTEGER NOUTV,NSPOOLV,NSTAV
!$omp threadprivate(NOUTV,NSPOOLV,NSTAV)
      INTEGER NOUTC,NSPOOLC,NSTAC
!$omp threadprivate(NOUTC,NSPOOLC,NSTAC)
      INTEGER NOUTM,NSPOOLM,NSTAM
!$omp threadprivate(NOUTM,NSPOOLM,NSTAM)
      INTEGER NOUTGE,NSPOOLGE
!$omp threadprivate(NOUTGE,NSPOOLGE)
      INTEGER NOUTGV,NSPOOLGV
!$omp threadprivate(NOUTGV,NSPOOLGV)
      INTEGER NOUTGC,NSPOOLGC
!$omp threadprivate(NOUTGC,NSPOOLGC)
      INTEGER NOUTGW,NSPOOLGW
!$omp threadprivate(NOUTGW,NSPOOLGW)
      INTEGER NHASE,NHASV,NHAGE,NHAGV
!$omp threadprivate(NHASE,NHASV,NHAGE,NHAGV)
C
C
C  STRING15 section
      CHARACTER*80 RUNDES,RUNID
!$omp threadprivate(RUNDES,RUNID)
C
C  PARM15-3DVS section
      INTEGER IDEN ! 3D run type
!$omp threadprivate(IDEN)
      INTEGER NFEN ! number of vertical nodes 
!$omp threadprivate(NFEN)
      INTEGER :: I3DSD,NSPO3DSD,NSTA3DD
!$omp threadprivate(I3DSD,NSPO3DSD,NSTA3DD)
      INTEGER :: I3DSV,NSPO3DSV,NSTA3DV
!$omp threadprivate(I3DSV,NSPO3DSV,NSTA3DV)
      INTEGER :: I3DST,NSPO3DST,NSTA3DT
!$omp threadprivate(I3DST,NSPO3DST,NSTA3DT)
      REAL(SZ) :: TO3DSDS,TO3DSDF
!$omp threadprivate(TO3DSDS,TO3DSDF)
      REAL(SZ) :: TO3DSVS,TO3DSVF
!$omp threadprivate(TO3DSVS,TO3DSVF)
      REAL(SZ) :: TO3DSTS,TO3DSTF
!$omp threadprivate(TO3DSTS,TO3DSTF)
      INTEGER :: I3DGD,NSPO3DGD
!$omp threadprivate(I3DGD,NSPO3DGD)
      INTEGER :: I3DGV,NSPO3DGV
!$omp threadprivate(I3DGV,NSPO3DGV)
      INTEGER :: I3DGT,NSPO3DGT
!$omp threadprivate(I3DGT,NSPO3DGT)
      REAL(SZ) :: TO3DGDS,TO3DGDF
!$omp threadprivate(TO3DGDS,TO3DGDF)
      REAL(SZ) :: TO3DGVS,TO3DGVF
!$omp threadprivate(TO3DGVS,TO3DGVF)
      REAL(SZ) :: TO3DGTS,TO3DGTF
!$omp threadprivate(TO3DGTS,TO3DGTF)
C
C  ELESTAT section
      REAL(SZ) TOUTSE,TOUTFE,TOUTSGE,TOUTFGE
!$omp threadprivate(TOUTSE,TOUTFE,TOUTSGE,TOUTFGE)
C
C  VELSTAT section
      REAL(SZ) TOUTSGV,TOUTFGV,TOUTSV,TOUTFV
!$omp threadprivate(TOUTSGV,TOUTFGV,TOUTSV,TOUTFV)
C
C  CONSTAT  section
      REAL(SZ) TOUTSC,TOUTFC,TOUTSGC,TOUTFGC
!$omp threadprivate(TOUTSC,TOUTFC,TOUTSGC,TOUTFGC)
      REAL(SZ) TOUTSGW,TOUTFGW
!$omp threadprivate(TOUTSGW,TOUTFGW)
C  Meteorological recording station section
      REAL(SZ) TOUTFM, TOUTSM
!$omp threadprivate(TOUTFM,TOUTSM)

C
C--Degress-to-Radians and Radians-to-Degrees
C
C  CONVERT section
      REAL(8) DEG2RAD,RAD2DEG,R
!$omp threadprivate(DEG2RAD,RAD2DEG,R)
C
C
C--Used to DASD file operations
C
C DASD
      INTEGER NBYTE
!$omp threadprivate(NBYTE)
C
C--------------------------------------------------------------------------C
C                                                                          C
C              DATA DECOMPOSITION DECLARATIONS BEGIN HERE                  C
C                                                                          C
C--------------------------------------------------------------------------C
C
C--Local Map Variable Declarations
C
Cvjp Note:  NSTACP  and IMAP_STAC_LG are not used yet
C
C LOCALI section
      INTEGER NPROC
!$omp threadprivate(NPROC)
      INTEGER, ALLOCATABLE ::  NNODP(:),NELP(:)
!$omp threadprivate(NNODP,NELP)
      INTEGER, ALLOCATABLE ::  NOD_RES_TOT(:)
!$omp threadprivate(NOD_RES_TOT)
      INTEGER, ALLOCATABLE ::  NSTAEP(:), NSTAVP(:), NSTACP(:)
!$omp threadprivate(NSTAEP,NSTAVP,NSTACP)
C     Number of subdomain atmospheric pressure recording stations
      INTEGER, ALLOCATABLE ::  NSTAMP(:) 
!$omp threadprivate(NSTAMP)
C
C
C--Local-to-Global Mapping Variables
C
C  LOC2G section
      INTEGER, ALLOCATABLE :: IMAP_NOD_LG(:,:)
!$omp threadprivate(IMAP_NOD_LG)
      INTEGER, ALLOCATABLE :: IMAP_EL_LG(:,:)!jgf45.07
!$omp threadprivate(IMAP_EL_LG)
      INTEGER, ALLOCATABLE :: IMAP_STAE_LG(:,:)
!$omp threadprivate(IMAP_STAE_LG)
      INTEGER, ALLOCATABLE :: IMAP_STAV_LG(:,:)
!$omp threadprivate(IMAP_STAV_LG)
      INTEGER, ALLOCATABLE :: IMAP_STAC_LG(:,:)
!$omp threadprivate(IMAP_STAC_LG)
      INTEGER, ALLOCATABLE :: IMAP_STAM_LG(:,:) ! atmospheric pressure
!$omp threadprivate(IMAP_STAM_LG)
C
C--Global-to-Local Mapping Variables
C
C  GLOB2L section
      INTEGER, ALLOCATABLE :: IMAP_NOD_GL(:,:)
!$omp threadprivate(IMAP_NOD_GL)
C
C--3DVS section
C
      INTEGER,ALLOCATABLE  ::  NNSTA3DDP(:),NNSTA3DVP(:),NNSTA3DTP(:)
!$omp threadprivate(NNSTA3DDP,NNSTA3DVP,NNSTA3DTP)
      INTEGER,ALLOCATABLE  ::  IMAP_STA3DD_LG(:,:)
!$omp threadprivate(IMAP_STA3DD_LG)
      INTEGER,ALLOCATABLE  ::  IMAP_STA3DV_LG(:,:)
!$omp threadprivate(IMAP_STA3DV_LG)
      INTEGER,ALLOCATABLE  ::  IMAP_STA3DT_LG(:,:)
!$omp threadprivate(IMAP_STA3DT_LG)
C
C--------------------------------------------------------------------------C
c   END OF DECLARATIONS
C--------------------------------------------------------------------------C

       CONTAINS

       SUBROUTINE ALLOC_MAIN1()
C                     
      ALLOCATE ( NNODP(MNPROC), NELP(MNPROC))
      ALLOCATE ( NOD_RES_TOT(MNPROC))
      ALLOCATE ( NSTAEP(MNPROC), NSTAVP(MNPROC), NSTACP(MNPROC))
      ALLOCATE ( NSTAMP(MNPROC) )
      ALLOCATE ( IMAP_NOD_LG(MNPP,MNPROC))
      ALLOCATE ( IMAP_EL_LG(MNEP,MNPROC)) !jgf45.07
      ALLOCATE ( IMAP_NOD_GL(2,MNP))

      ALLOCATE ( IMAP_STAE_LG(MNSTAE,MNPROC))
      ALLOCATE ( IMAP_STAV_LG(MNSTAV,MNPROC))
      ALLOCATE ( IMAP_STAC_LG(MNSTAC,MNPROC))
      ALLOCATE ( IMAP_STAM_LG(MNSTAM,MNPROC))
C
      END SUBROUTINE ALLOC_MAIN1
      END MODULE POST_GLOBAL

