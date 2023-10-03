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
C     jgf Updated for 45.07 to add handling of fort.71, fort.72, fort.73
C         and the NOFF array in output files 
C     jgf Updated for 45.06 09/07/2005 b/c of changes in hot start file format
C     jgf Updated for 45.11 02/02/2006 b/c for 3D recording stations
C     jgf Updated for 45.12 03/17/2006 
      MODULE POST_GLOBAL
      IMPLICIT NONE

      INTEGER MNPROC,MNP,MNE,MNPP,MNSTAE,MNSTAV,MNHARF,MNWLAT,MNWLON
      INTEGER MNSTAC !jgf46.02 # of conc. rec. sta. in full domain
      INTEGER MNSTAM ! number of met. recording stations in full domain
      INTEGER MNEP   !jgf45.07 Max number of elements in any subdomain
C
C   PARM14   section
      INTEGER NELG,NNODG
C
C   STRING14 section
      CHARACTER*80 AGRID
C
C   PARM15  section
      INTEGER IM
      INTEGER NWS
      INTEGER NOUTE,NSPOOLE,NSTAE
      INTEGER NOUTV,NSPOOLV,NSTAV
      INTEGER NOUTC,NSPOOLC,NSTAC
      INTEGER NOUTM,NSPOOLM,NSTAM
      INTEGER NOUTGE,NSPOOLGE
      INTEGER NOUTGV,NSPOOLGV
      INTEGER NOUTGC,NSPOOLGC
      INTEGER NOUTGW,NSPOOLGW
      INTEGER NHASE,NHASV,NHAGE,NHAGV
C
C
C  STRING15 section
      CHARACTER*80 RUNDES,RUNID
C
C  PARM15-3DVS section
      INTEGER IDEN ! 3D run type
      INTEGER NFEN ! number of vertical nodes 
      INTEGER :: I3DSD,NSPO3DSD,NSTA3DD
      INTEGER :: I3DSV,NSPO3DSV,NSTA3DV
      INTEGER :: I3DST,NSPO3DST,NSTA3DT
      REAL(8) :: TO3DSDS,TO3DSDF
      REAL(8) :: TO3DSVS,TO3DSVF
      REAL(8) :: TO3DSTS,TO3DSTF
      INTEGER :: I3DGD,NSPO3DGD
      INTEGER :: I3DGV,NSPO3DGV
      INTEGER :: I3DGT,NSPO3DGT
      REAL(8) :: TO3DGDS,TO3DGDF
      REAL(8) :: TO3DGVS,TO3DGVF
      REAL(8) :: TO3DGTS,TO3DGTF
C
C  ELESTAT section
      REAL(8) TOUTSE,TOUTFE,TOUTSGE,TOUTFGE
C
C  VELSTAT section
      REAL(8) TOUTSGV,TOUTFGV,TOUTSV,TOUTFV
C
C  CONSTAT  section
      REAL(8) TOUTSC,TOUTFC,TOUTSGC,TOUTFGC
      REAL(8) TOUTSGW,TOUTFGW
C  Meteorological recording station section
      REAL(8) TOUTFM, TOUTSM

C
C--Degress-to-Radians and Radians-to-Degrees
C
C  CONVERT section
      REAL(8) DEG2RAD,RAD2DEG,R
C
C
C--Used to DASD file operations
C
C DASD
      INTEGER NBYTE
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
      INTEGER, ALLOCATABLE ::  NNODP(:),NELP(:)
      INTEGER, ALLOCATABLE ::  NOD_RES_TOT(:)
      INTEGER, ALLOCATABLE ::  NSTAEP(:), NSTAVP(:), NSTACP(:)
C     Number of subdomain atmospheric pressure recording stations
      INTEGER, ALLOCATABLE ::  NSTAMP(:) 
C
C
C--Local-to-Global Mapping Variables
C
C  LOC2G section
      INTEGER, ALLOCATABLE :: IMAP_NOD_LG(:,:)
      INTEGER, ALLOCATABLE :: IMAP_EL_LG(:,:)!jgf45.07
      INTEGER, ALLOCATABLE :: IMAP_STAE_LG(:,:)
      INTEGER, ALLOCATABLE :: IMAP_STAV_LG(:,:)
      INTEGER, ALLOCATABLE :: IMAP_STAC_LG(:,:)
      INTEGER, ALLOCATABLE :: IMAP_STAM_LG(:,:) ! atmospheric pressure
C
C--Global-to-Local Mapping Variables
C
C  GLOB2L section
      INTEGER, ALLOCATABLE :: IMAP_NOD_GL(:,:)
C
C--3DVS section
C
      INTEGER,ALLOCATABLE  ::  NNSTA3DDP(:),NNSTA3DVP(:),NNSTA3DTP(:)
      INTEGER,ALLOCATABLE  ::  IMAP_STA3DD_LG(:,:)
      INTEGER,ALLOCATABLE  ::  IMAP_STA3DV_LG(:,:)
      INTEGER,ALLOCATABLE  ::  IMAP_STA3DT_LG(:,:)
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

