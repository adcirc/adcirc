!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                         M O D U L E  S U B G R I D
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! EXECUTES THE SUBGRID SUBROUTINES THROUGHOUT THE CODE !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FIRST WE CREATE ALL OF THE VARIABLES NEEDED FOR THE SUBGRID ROUTINES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      module subgrid
      use sizes, only: sz
      use global, only: nodecode, nnodecode
     &   , eta2, tk, tkm, nolifa, tkm
     &   , H0, ifnlfa, uu1, vv1, ncchange
      use mesh, only: ne, np, nm, totalArea, mju, areas, dp
      
      

! JLW: initialize the arrays for use in the subgrid
! code
      implicit none

      ! Flag for subgrid corrections to wetting and drying
      LOGICAL :: level0
      ! Flag for subgrid corrections to bottom friction and advection
      LOGICAL :: level1
      ! Number of surface elevations used in subgrid preprocessor
!      INTEGER :: numSurLevs
!      ! Number of possible phi values
      INTEGER :: numPhi
      ! Array of Phi values from 0 to 1
<<<<<<< Updated upstream
      REAL(SZ), DIMENSION(:), ALLOCATABLE :: setPhi
      ! Max depth used for calculating lookup table
      REAL(SZ):: maxLookupDepth
      ! Min depth used for calculating lookup table
      REAL(SZ) :: minLookupDepth
      ! Array of element numbers in the subgrid area
      INTEGER, DIMENSION(:), ALLOCATABLE :: subgridEleList
      ! Depths corresponding to each phi value on the element
      REAL(SZ), DIMENSION(:,:,:), ALLOCATABLE :: wetFracDepthEle
      ! Grid Averaged Total Water Depths corresponding to each phi valu
      ! ! on the element
      REAL(SZ), DIMENSION(:,:,:), ALLOCATABLE :: gridDepthEle
      ! Area of sub elements
      REAL(SZ), DIMENSION(:,:), ALLOCATABLE :: area
      ! Advection corrections corresponding to each phi value on the
      ! element
      REAL(SZ), DIMENSION(:,:,:), ALLOCATABLE :: cadvEle
=======
      REAL(8), DIMENSION(:), ALLOCATABLE :: setPhi
>>>>>>> Stashed changes
      ! averaged bottom friction coefficient without level 1 correction
      REAL(SZ), DIMENSION(:,:), ALLOCATABLE :: cfVertTab
      ! averaged bottom friction coefficient with level 1 correction
<<<<<<< Updated upstream
      REAL(SZ), DIMENSION(:,:), ALLOCATABLE :: cmfVertTab
=======
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: cmfVertTab
      ! averaged advection with level 1 correction
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: cadvVertTab
>>>>>>> Stashed changes
      ! Grid Averaged Total Water Depths corresponding to each phi valu
      ! on the vertex
      REAL(SZ), DIMENSION(:,:), ALLOCATABLE :: gridDepthVertTab
      ! Wet Averaged Total Water Depths corresponding to each phi value
      ! on the vertex
      REAL(SZ), DIMENSION(:,:), ALLOCATABLE :: wetDepthVertTab
      ! Depths corresponding to each phi value on the vertex
      REAL(SZ), DIMENSION(:,:), ALLOCATABLE :: wetFracVertTab
      ! Array of vertex numbers in the subgrid area
      INTEGER, DIMENSION(:), ALLOCATABLE :: subgridVertList
<<<<<<< Updated upstream
      !! Array of elemental wet fractions for zeta at TS k
      !REAL(SZ), DIMENSION(:,:), ALLOCATABLE :: wetFracEleETA1
      ! Array of elemental wet fractions for zeea at TS k+1
      REAL(SZ), DIMENSION(:,:), ALLOCATABLE :: wetFracEleETA2
      ! Array of elemental grid depths for zeta at k
      REAL(SZ), DIMENSION(:,:), ALLOCATABLE :: gridDepthEleETA1
      ! Array of elemental grid depths for zeta at k+1
      REAL(SZ), DIMENSION(:,:), ALLOCATABLE :: gridDepthEleETA2
      !! Array of elemental advection corrections for zeta at k
      !REAL(SZ), DIMENSION(:,:), ALLOCATABLE :: cadvEleETA1
      ! Array of elemental advection corrections for zeta at k+1
      REAL(SZ), DIMENSION(:,:), ALLOCATABLE :: cadvEleETA2
      ! Array of elemental average wet fractions for the entire element
      ! for zeta at k
      REAL(SZ), DIMENSION(:), ALLOCATABLE :: wetFracAvgETA1
      ! Array of elemental average wet fractions for the entire element
      ! for zeta at k+1
      REAL(SZ), DIMENSION(:), ALLOCATABLE :: wetFracAvgETA2
      !! Array to hold average elemental grid depths
      !REAL(SZ), DIMENSION(:), ALLOCATABLE :: HGAvgETA1
      !! Array to hold average elemental grid depths
      !REAL(SZ), DIMENSION(:), ALLOCATABLE :: HGAvgETA2
      ! Array to hold wet/dry state of element
      INTEGER, DIMENSION(:), ALLOCATABLE :: NCEleArray
      ! Array to hold vertex averged wet fractions for zeta at k+1
      REAL(SZ), DIMENSION(:), ALLOCATABLE :: wetFracVertETA2
=======
      ! Array to hold vertex averged wet fractions for zeta at k
      REAL(8), DIMENSION(:), ALLOCATABLE :: wetFracVertETA1
      ! Array to hold vertex averged wet fractions for zeta at k+1
      REAL(8), DIMENSION(:), ALLOCATABLE :: wetFracVertETA2
      ! Array to hold vertex grid water depths for zeta at k
      REAL(8), DIMENSION(:), ALLOCATABLE :: gridDepthVertETA1
>>>>>>> Stashed changes
      ! Array to hold vertex grid water depths for zeta at k+1
      REAL(SZ), DIMENSION(:), ALLOCATABLE :: gridDepthVertETA2
      ! Array to hold vertex wet water depths for zeta at k+1
      REAL(SZ), DIMENSION(:), ALLOCATABLE :: wetDepthVertETA2
      ! Array to hold vertex bottom friction coefficients for zeta at
      ! k+1
      REAL(SZ), DIMENSION(:), ALLOCATABLE :: cfVertETA2
      ! Array to hold vertex corrected bottom friction coefficients for
      ! zeta at k+1
<<<<<<< Updated upstream
      REAL(SZ), DIMENSION(:), ALLOCATABLE :: cmfVertETA2
=======
      REAL(8), DIMENSION(:), ALLOCATABLE :: cmfVertETA2
      ! Array to hold vertex corrected advection for
      ! zeta at k+1
      REAL(8), DIMENSION(:), ALLOCATABLE :: cadvVertETA2
>>>>>>> Stashed changes
      ! Minimum depth in depth array used to preprocess lookup table 
      REAL(SZ) :: minDepth
      ! Maximum depth in depth array used to preprocess lookup table 
      REAL(SZ) :: maxDepth

      contains

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!
!                  INTIALIZE VARIABLES
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! THIS SECTION OF CODE INITIALIZES THE ARRAYS NEEDED FOR THE LOOKUP
! TABLES 

      subroutine initVar

!JLW: initialize local arrays for subgrid variables
      IF(ALLOCATED(wetFracVertETA1))DEALLOCATE(wetFracVertETA1)
      ALLOCATE(wetFracVertETA1(NP))
      IF(ALLOCATED(gridDepthVertETA1))DEALLOCATE(gridDepthVertETA1)
      ALLOCATE(gridDepthVertETA1(NP))
      IF(ALLOCATED(wetFracVertETA2))DEALLOCATE(wetFracVertETA2)
      ALLOCATE(wetFracVertETA2(NP))
      IF(ALLOCATED(gridDepthVertETA2))DEALLOCATE(gridDepthVertETA2)
      ALLOCATE(gridDepthVertETA2(NP))
      IF(ALLOCATED(wetDepthVertETA2))DEALLOCATE(wetDepthVertETA2)
      ALLOCATE(wetDepthVertETA2(NP))
      IF(level1)THEN
         IF(ALLOCATED(cmfVertETA2))DEALLOCATE(cmfVertETA2)
         ALLOCATE(cmfVertETA2(NP))
         IF(ALLOCATED(cadvVertETA2))DEALLOCATE(cadvVertETA2)
         ALLOCATE(cadvVertETA2(NP))
!JLW: will need to add vertex based advection
      ELSE
         IF(ALLOCATED(cfVertETA2))DEALLOCATE(cfVertETA2)
         ALLOCATE(cfVertETA2(NP))
      ENDIF

!JLW: set intial value of all arrays to 0
      wetFracVertETA1(:) = 0.d0
      wetFracVertETA2(:) = 0.d0
      gridDepthVertETA1(:) = 0.d0
      gridDepthVertETA2(:) = 0.d0
      wetDepthVertETA2(:) = 0.d0
      cmfVertETA2(:) = 0.d0
      cadvVertETA2(:) = 1.d0
      cfVertETA2(:) = 0.d0
!JLW: will need to add vertex based advection

!----------------------------------------------------------------------
!----------------------------------------------------------------------
      end subroutine initVar
!----------------------------------------------------------------------
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!
!                  SUBROUTINE READ IN NETCDF FILE
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! READS IN THE NETCDF LOOKUP TABLE AND POPULATES THE SUBGRID ARRAYS BOTH
! FOR SERIAL AND PARALLEL ADCIRC. THIS SUBROUTINE USES A FUNCTION
! "CHECK" TO READ IN THE NETCDF FILE. THE "CHECK" SUBROUTINE IS LOCATED
! AT THE BOTTOM OF THIS FILE.

      subroutine readSubgridLookup
      use netcdf
<<<<<<< Updated upstream
      use global, only: imap_el_lg, nodes_lg
=======
      ! use global, only: imap_el_lg, nodes_lg, allMessage, ERROR
      use global, only: nodes_lg, allMessage, ERROR
>>>>>>> Stashed changes
      
      implicit none

      ! Total number of elements and vertices in mesh (can most
      ! definitely be replaced by variables already in ADCIRC) 
      !integer :: numElems, numVerts, numVertsPerEle
      integer :: numVerts
      ! Used for transfering from global to local
      ! integer :: currEle, currNode
      integer :: currNode
      ! Arrays for holding global variables
<<<<<<< Updated upstream
      REAL(SZ), DIMENSION(:,:,:), ALLOCATABLE :: globalEleLookup
      REAL(SZ), DIMENSION(:,:), ALLOCATABLE :: globalVertLookup
      REAL(SZ), DIMENSION(:,:), ALLOCATABLE :: globalarea
      INTEGER, DIMENSION(:), ALLOCATABLE :: globalSubEleList
=======
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: globalVertLookup
>>>>>>> Stashed changes
      INTEGER, DIMENSION(:), ALLOCATABLE :: globalSubVertList
      ! variables needed to read in main NETCDF lookup table
      INTEGER :: NC_ERR, IS, NC_ID, NC_VAR 
      INTEGER :: I

!JLW: adding subgrid lookup table read in
!JLW: first open subgrid lookup table and read in dimensions
      CALL CHECK(NF90_OPEN("avgVars.nc",NF90_NOWRITE,NC_ID))
      
      CALL CHECK(NF90_INQ_DIMID(NC_ID,"numNode",NC_VAR))
      CALL CHECK(NF90_INQUIRE_DIMENSION(NC_ID,NC_VAR,len=numVerts))
      
      CALL CHECK(NF90_INQ_DIMID(NC_ID,"numPhi",NC_VAR))
      CALL CHECK(NF90_INQUIRE_DIMENSION(NC_ID,NC_VAR,len=numPhi))

      IF(ALLOCATED(setPhi)) DEALLOCATE(setPhi)
      ALLOCATE(setPhi(numPhi))

      CALL Check(NF90_INQ_VARID(NC_ID,"phiSet",NC_VAR))
      CALL Check(NF90_GET_VAR(NC_ID,NC_VAR,setPhi))

!JLW: go through and read in each of the vertex averaged variables
!JLW: allocate arrays  
      IF(level1)THEN
         IF(ALLOCATED(cmfVertTab)) DEALLOCATE(cmfVertTab)
         ALLOCATE(cmfVertTab(NP,numPhi))
         IF(ALLOCATED(cadvVertTab)) DEALLOCATE(cadvVertTab)
         ALLOCATE(cadvVertTab(NP,numPhi))
      ELSE
         IF(ALLOCATED(cfVertTab)) DEALLOCATE(cfVertTab)
         ALLOCATE(cfVertTab(NP,numPhi))
      ENDIF
      IF(ALLOCATED(gridDepthVertTab))DEALLOCATE(gridDepthVertTab)
      ALLOCATE(gridDepthVertTab(NP,numPhi))
      IF(ALLOCATED(wetDepthVertTab))DEALLOCATE(wetDepthVertTab)
      ALLOCATE(wetDepthVertTab(NP,numPhi))
      IF(ALLOCATED(wetFracVertTab))DEALLOCATE(wetFracVertTab)
      ALLOCATE(wetFracVertTab(NP,numPhi))
      IF(ALLOCATED(subgridVertList))DEALLOCATE(subgridVertList)
      ALLOCATE(subgridVertList(NP))

!!!!!!!!!!!!!!!!!! bottom friction coefficient !!!!!!!!!!!!!!!!!!!!!!
      IF(ALLOCATED(globalVertLookup))DEALLOCATE(globalVertLookup)
      ALLOCATE(globalVertLookup(numPhi,numVerts))

      IF(level1)THEN
         CALL Check(NF90_INQ_VARID(NC_ID,"cmfVertex",NC_VAR))
         CALL Check(NF90_GET_VAR(NC_ID,NC_VAR,globalVertLookup))
         DO I = 1,NP
#ifdef CMPI
            currNode = abs(nodes_lg(I))
#else
            currNode = I
#endif
            cmfVertTab(I,:) = globalVertLookup(:,currNode)
         ENDDO
      ELSE
         CALL Check(NF90_INQ_VARID(NC_ID,"cfVertex",NC_VAR))
         CALL Check(NF90_GET_VAR(NC_ID,NC_VAR,globalVertLookup))
         DO I = 1,NP
#ifdef CMPI
            currNode = abs(nodes_lg(I))
#else
            currNode = I
#endif
            cfVertTab(I,:) = globalVertLookup(:,currNode)
         ENDDO
      ENDIF
!!!!!!!!!!!!!!!!!! ADVECTION CORRECTION VETEX !!!!!!!!!!!!!!!!!!
      IF(ALLOCATED(globalVertLookup))DEALLOCATE(globalVertLookup)
      ALLOCATE(globalVertLookup(numPhi,numVerts))
      IF(level1)THEN
         CALL Check(NF90_INQ_VARID(NC_ID,"cadvVertex",NC_VAR))
         CALL Check(NF90_GET_VAR(NC_ID,NC_VAR,globalVertLookup))
         DO I = 1,NP
#ifdef CMPI
            currNode = abs(nodes_lg(I))
#else
            currNode = I
#endif
            cadvVertTab(I,:) = globalVertLookup(:,currNode)
         ENDDO
      ENDIF
!!!!!!!!!!!!!!!!! grid depth vertex !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF(ALLOCATED(globalVertLookup))DEALLOCATE(globalVertLookup)
      ALLOCATE(globalVertLookup(numPhi,numVerts))

      CALL Check(NF90_INQ_VARID(NC_ID,"gridTotWatDepthVertex",NC_VAR))
      CALL Check(NF90_GET_VAR(NC_ID,NC_VAR,globalVertLookup))
      DO I = 1,NP
#ifdef CMPI
            currNode = abs(nodes_lg(I))
#else
            currNode = I
#endif
            gridDepthVertTab(I,:) = globalVertLookup(:,currNode)
      ENDDO
!!!!!!!!!!!!!!!!! wet depth vertex !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF(ALLOCATED(globalVertLookup))DEALLOCATE(globalVertLookup)
      ALLOCATE(globalVertLookup(numPhi,numVerts))

      CALL Check(NF90_INQ_VARID(NC_ID,"wetTotWatDepthVertex",NC_VAR))
      CALL Check(NF90_GET_VAR(NC_ID,NC_VAR,globalVertLookup))
      DO I = 1,NP
#ifdef CMPI
            currNode = abs(nodes_lg(I))
#else
            currNode = I
#endif
            wetDepthVertTab(I,:) = globalVertLookup(:,currNode)
      ENDDO
!!!!!!!!!!!!!!!!!!!! wet fraction depths vertex !!!!!!!!!!!!!!!!!!!!!!
      IF(ALLOCATED(globalVertLookup))DEALLOCATE(globalVertLookup)
      ALLOCATE(globalVertLookup(numPhi,numVerts))

      CALL Check(NF90_INQ_VARID(NC_ID,"wetFractionVertex",NC_VAR))
      CALL Check(NF90_GET_VAR(NC_ID,NC_VAR,globalVertLookup))
      DO I = 1,NP
#ifdef CMPI
            currNode = abs(nodes_lg(I))
#else
            currNode = I
#endif
            wetFracVertTab(I,:) = globalVertLookup(:,currNode)
      ENDDO
      IF(ALLOCATED(globalVertLookup))DEALLOCATE(globalVertLookup)
!!!!!!!!!!!!!!!!! subgrid flag for vertex !!!!!!!!!!!!!!!!!!!!!!!!
      IF(ALLOCATED(globalSubVertList))DEALLOCATE(globalSubVertList)
      ALLOCATE(globalSubVertList(numVerts))

      CALL Check(NF90_INQ_VARID(NC_ID,"binaryVertexList",NC_VAR))
      CALL Check(NF90_GET_VAR(NC_ID,NC_VAR,globalSubVertList))

      DO I = 1,NP
#ifdef CMPI
         currNode = abs(nodes_lg(I))
#else
         currNode = I
#endif
         subgridVertList(I) = globalSubVertList(currNode)
      ENDDO
      IF(ALLOCATED(globalSubVertList))DEALLOCATE(globalSubVertList)

!JLW: close netcdf lookup table file
      NC_ERR = NF90_CLOSE(NC_ID)

!----------------------------------------------------------------------
!----------------------------------------------------------------------
      end subroutine readSubgridLookup
!----------------------------------------------------------------------
!----------------------------------------------------------------------


!---------------------------------------------------------------------
!---------------------------------------------------------------------
!
!                  GET VERTEX AVERAGED VARIABLES
!
!---------------------------------------------------------------------
!---------------------------------------------------------------------
! THIS CODE IS VERY SIMILAR TO THE ELEMENTAL SUBROUTINE BUT FOR SUBRID
! VARIABLES RELATED TO THE VERTEX AREAS. FOR EACH NEW WATER SURFACE
! ELEVATION SUBGRID VERTEX QUANTITIES ARE LOOKED UP. 

      subroutine getVertLookup

      use sizes, only : sz, mne, mnp
      use mesh, only : ne, np, dp, mju, totalArea, nm, x, y, areas
     &   , NNeigh, NeiTab
      use global, only : H0

      implicit none
      integer :: ie, i, j
      integer :: numGreater
      REAL(8) :: HABSMIN

!JLW: set k = k + 1 for new timestep
      HABSMIN = 0.8d0*H0
      wetFracVertETA1 = wetFracVertETA2
      gridDepthVertETA1 = gridDepthVertETA2
!JLW: now loop through the vertices and lookup vertex averages
      DO I = 1,NP
         IF(subgridVertList(I).EQ.1)THEN
!JLW: find how many of the depths are greater than the current wet area
!fraction
            numGreater = 0 !initialize counter 
            DO j = 1,numPhi
               IF(eta2(i).LT.wetFracVertTab(i,j))THEN
                  numGreater = numGreater + 1
               ENDIF
            ENDDO
            IF(numGreater.EQ.numPhi)THEN
!JLW: this means that none of the depths in the array are greater than
!the current depth and this node is dry
               wetFracVertETA2(i) = setPhi(1)
               gridDepthVertETA2(i) = gridDepthVertTab(i,1)
               wetDepthVertETA2(i) = wetDepthVertTab(i,1)
               IF(level1)THEN
                  cmfVertETA2(i) = cmfVertTab(i,1)
                  cadvVertETA2(i) = cadvVertTab(i,1)
               ELSE
                  cfVertETA2(i) = cfVertTab(i,1)
               ENDIF
            ELSEIF(numGreater.EQ.0)THEN
!JLW: this means that none of the depths in the array are greater than
!the current depth and this node is always wet
               wetFracVertETA2(i) = setPhi(numPhi)
               gridDepthVertETA2(i) = gridDepthVertTab(i,numPhi)
     &            + (eta2(i) - wetFracVertTab(i,numPhi))
               wetDepthVertETA2(i) = wetDepthVertTab(i,numPhi)
     &            + (eta2(i) - wetFracVertTab(i,numPhi))
               IF(level1)THEN
                  cmfVertETA2(i) = cmfVertTab(i,numPhi)
                  cadvVertETA2(i) = cadvVertTab(i,numPhi)
               ELSE
                  cfVertETA2(i) = cfVertTab(i,numPhi)
               ENDIF
            ELSEIF(wetFracVertTab(i,numPhi-numGreater).NE.
!JLW: as long as the depth we are indexing is not equal to the minimum
!depth that we used to calculate the subgrid lookup table we interpolate
!between depths to find subgrid values
     &         minDepth)THEN
               wetFracVertETA2(i) = ((eta2(i)
     &            - wetFracVertTab(i,numPhi-numGreater))
     &            /(wetFracVertTab(i,numPhi-numGreater+1)
     &            - wetFracVertTab(i,numPhi-numGreater))
     &            *(setPhi(numPhi-numGreater+1)
     &            - setPhi(numPhi-numGreater))
     &            + setPhi(numPhi-numGreater))
               gridDepthVertETA2(i) = ((eta2(i)
     &            - wetFracVertTab(i,numPhi-numGreater))
     &            /(wetFracVertTab(i,numPhi-numGreater+1)
     &            - wetFracVertTab(i,numPhi-numGreater))
     &            *(gridDepthVertTab(i,numPhi-numGreater+1)
     &            - gridDepthVertTab(i,numPhi-numGreater))
     &            + gridDepthVertTab(i,numPhi-numGreater))
               wetDepthVertETA2(i) = ((eta2(i)
     &            - wetFracVertTab(i,numPhi-numGreater))
     &            /(wetFracVertTab(i,numPhi-numGreater+1)
     &            - wetFracVertTab(i,numPhi-numGreater))
     &            *(wetDepthVertTab(i,numPhi-numGreater+1)
     &            - wetDepthVertTab(i,numPhi-numGreater))
     &            + wetDepthVertTab(i,numPhi-numGreater))
               IF(level1)THEN
                  cmfVertETA2(i) = ((eta2(i)
     &               - wetFracVertTab(i,numPhi-numGreater))
     &               /(wetFracVertTab(i,numPhi-numGreater+1)
     &               - wetFracVertTab(i,numPhi-numGreater))
     &               *(cmfVertTab(i,numPhi-numGreater+1)
     &               - cmfVertTab(i,numPhi-numGreater))
     &               + cmfVertTab(i,numPhi-numGreater))
                  cadvVertETA2(i) = ((eta2(i)
     &               - wetFracVertTab(i,numPhi-numGreater))
     &               /(wetFracVertTab(i,numPhi-numGreater+1)
     &               - wetFracVertTab(i,numPhi-numGreater))
     &               *(cadvVertTab(i,numPhi-numGreater+1)
     &               - cadvVertTab(i,numPhi-numGreater))
     &               + cadvVertTab(i,numPhi-numGreater))
               ELSE
                  cfVertETA2(i) = ((eta2(i)
     &               - wetFracVertTab(i,numPhi-numGreater))
     &               /(wetFracVertTab(i,numPhi-numGreater+1)
     &               - wetFracVertTab(i,numPhi-numGreater))
     &               *(cfVertTab(i,numPhi-numGreater+1)
     &               - cfVertTab(i,numPhi-numGreater))
     &               + cfVertTab(i,numPhi-numGreater))
               ENDIF
            ELSE
!JLW: if the current water level is equal to the first depth of the
!array wet set our averaged variables to fully wet
               wetFracVertETA2(i) = setPhi(numPhi-numGreater+1)
               gridDepthVertETA2(i) = 
     &            gridDepthVertTab(i,numPhi-numGreater+1)
               wetDepthVertETA2(i) = 
     &            wetDepthVertTab(i,numPhi-numGreater+1)
               IF(level1)THEN
                  cmfVertETA2(i) = cmfVertTab(i,numPhi-numGreater+1)
                  cadvVertETA2(i) = cadvVertTab(i,numPhi-numGreater+1)
               ELSE
                  cfVertETA2(i) = cfVertTab(i,numPhi-numGreater+1)
               ENDIF
            ENDIF
         ELSE
         !JLW: if not in subgrid area set grid total depth (this was a
         !bug)
            gridDepthVertETA2(i) = eta2(i) + dp(i)
            wetDepthVertETA2(i) = eta2(i) + dp(i)
            wetFracVertETA2(i) = 1.d0
         ENDIF
      ENDDO
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      end subroutine getVertLookup
!----------------------------------------------------------------------
!----------------------------------------------------------------------


!---------------------------------------------------------------------
!---------------------------------------------------------------------
!
<<<<<<< Updated upstream
!                  WET DRY CHECK AND UPDATE
!
!---------------------------------------------------------------------
!---------------------------------------------------------------------
! AFTER THE ELEMENTAL SUBGRID ROUTINE AND THE VERTEX SUBROUTINE ARE
! CALLED AND THE SUBGRID ARRAYS ARE UPDATED, THIS WET DRY CODE IS CALLED
! TO UPDATE THE WET/DRY STATE OF THE ELEMENTS AND VERTICES. THIS
! CONSISTS OF A SIMPLE CHECK OF TOTAL WATER DEPTH H. IN ADDITION,
! UPDATED COEFFICIENTS OF FRICTION AND AREAS ARE ALSO UPDATED.

      subroutine subgridwetdry

      use NodalAttributes, only: BFCdLLimit, fric
     &   , fgamma, ftheta, manningsn, hbreak, loadManningsN
     &   , iflinbf, ifnlbf, ifhybf 
      use constants, only: g

      implicit none
      integer :: ie, i, j
      real(8) :: areaEle
      integer :: NCEle
      real(8) :: uv1
      integer :: nm1, nm2, nm3
      REAL(8) :: H1N1, H1N2, H1N3
      REAL(8) :: HGAvg2
      REAL(8) :: H2

!JLW: handle the vertex wetting and drying 
      nnodecode(:) = 0
      DO i = 1,np
         IF(subgridVertList(i).EQ.1)THEN
            H2 = gridDepthVertETA2(i)
         ELSE
            H2 = dp(i) + eta2(i) 
         ENDIF
         IF(H2.GT.H0)THEN
            nnodecode(i) = 1
         ELSE
            nnodecode(i) = 0
         ENDIF
!         IF(subgridVertList(i).EQ.1)THEN
!            IF(gridDepthVertETA2(i).GT.H0)THEN
!               nnodecode(i) = 1
!            ELSE
!               nnodecode(i) = 0
!            ENDIF
!         ELSE
!            IF(dp(i)+eta2(i).GT.H0)THEN
!               nnodecode(i) = 1
!            ELSE
!               nnodecode(i) = 0
!            ENDIF
!         ENDIF
!JLW: debating on removing -- don't think this is really necessary but
!depends on where and when the friction coefficients are updated in the
!timestep code.
         !JLW: now recaluate tk for the newly wet nodes
         IF((nodecode(i).EQ.0).AND.(nnodecode(i).EQ.1))THEN
            IF(LoadManningsN)THEN
               IF(subgridVertList(i).EQ.1)THEN
                  IF(level1)THEN
                     fric(i) = cmfVertETA2(i)
                  ELSE
                     fric(i) = cfVertETA2(i)
                  ENDIF
               ELSE
                  fric(i) = g*manningsn(i)**2.d0
     &               /((dp(i)+ifnlfa*eta2(i))**(1.d0/3.d0))
               ENDIF
               IF(fric(i).LT.BFCdLLimit)THEN
                  fric(i) = BFCdLLimit
               ENDIF
            ENDIF
            uv1 = SQRT(uu1(i)*uu1(i)+vv1(i)*vv1(i))
            IF(subgridVertList(i).EQ.1)THEN
               IF(wetDepthVertETA2(i).GT.0.0)THEN
                  tk(i) = fric(i)*
     &                 ( iflinbf +       ! linear
     &                 (uv1/wetDepthVertETA2(i)) * (ifnlbf ! nonlinear
     &                 + ifhybf*(1+(hbreak/wetDepthVertETA2(i))**ftheta)
     &                 **(fgamma/ftheta)))
               ENDIF
            ELSE
               IF(dp(i)+eta2(i).GT.0.0)THEN
                  tk(i) = fric(i)*
     &                 ( iflinbf +       ! linear
     &                 (uv1/(dp(i)+eta2(i))) * (ifnlbf ! nonlinear
     &                 + ifhybf*(1+(hbreak/(dp(i)+eta2(i)))**ftheta)
     &                 **(fgamma/ftheta)))
               ENDIF
            ENDIF
         ENDIF
      ENDDO
!JLW: update nodecode
      nodecode = nnodecode
!JLW: handle elemental wetting and drying 
      mju(:) = 0
      totalArea(:) = 0.d0
      DO ie = 1,ne
         nm1 = nm(ie,1)
         nm2 = nm(ie,2)
         nm3 = nm(ie,3)
!JLW: calculate the average total water depth
         H1N1 = gridDepthEleETA2(ie,1)
         H1N2 = gridDepthEleETA2(ie,2)
         H1N3 = gridDepthEleETA2(ie,3)
         HGAvg2 = (H1N1+H1N2+H1N3)/3.d0
         IF(subgridEleList(ie).EQ.1)THEN
!JLW: use the average total water depth to determine wet dry state of
!the element
            IF(HGAvg2.GT.H0)THEN
               NCEle = 1
            ELSE
               NCEle = 0
            ENDIF
         ELSE
            NCEle = nodecode(nm1)*nodecode(nm2)*nodecode(nm3)
         ENDIF
         NCEleArray(ie) = NCEle
         areaEle = NCEle*areas(IE)/2.d0
         mju(nm1)=mju(nm1)+NCEle
         mju(nm2)=mju(nm2)+NCEle
         mju(nm3)=mju(nm3)+NCEle
         totalArea(nm1)=totalArea(nm1)+areaEle
         totalArea(nm2)=totalArea(nm2)+areaEle
         totalArea(nm3)=totalArea(nm3)+areaEle
      ENDDO

!JLW: make sure to activate ncchange
      ncchange = 1

!----------------------------------------------------------------------
!----------------------------------------------------------------------
      end subroutine subgridwetdry
!----------------------------------------------------------------------
!----------------------------------------------------------------------


!---------------------------------------------------------------------
!---------------------------------------------------------------------
!
=======
>>>>>>> Stashed changes
!                 READ NETCDF UTILITY
!
!---------------------------------------------------------------------
!---------------------------------------------------------------------

!JLW: add subroutine for reading in netcdf
      SUBROUTINE Check(Status)
        USE netcdf
        INTEGER,INTENT(IN) :: Status
        IF(Status.NE.NF90_NOERR)THEN
            WRITE(*,'(A,A,A)') "FATAL ERRORfrom"
     &         ,TRIM(NF90_STRERROR(Status)),"."
            STOP
        ENDIF

!---------------------------------------------------------------------
!---------------------------------------------------------------------
      END SUBROUTINE
!---------------------------------------------------------------------
!---------------------------------------------------------------------



      end module subgrid
