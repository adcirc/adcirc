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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                         M O D U L E  S U B G R I D
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! EXECUTES THE SUBGRID SUBROUTINES THROUGHOUT THE CODE !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FIRST WE CREATE ALL OF THE VARIABLES NEEDED FOR THE SUBGRID ROUTINES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      module subgrid
      use global, only: nodecode, nnodecode
     &   , eta2, tk, tkm, nolifa, tkm
     &   , H0, ifnlfa, uu1, vv1, ncchange
      use mesh, only: ne, np, nm, totalArea, mju, areas, dp
      
      

! JLW: initialize the arrays for use in the subgrid
! code
      implicit none
      CHARACTER(200) :: subgridFilename
      ! Flag for subgrid corrections to wetting and drying
      LOGICAL :: level0
      ! Flag for subgrid corrections to bottom friction and advection
      LOGICAL :: level1
      ! Number of surface elevations used in subgrid preprocessor
      INTEGER :: numSurLevs
      ! Number of possible phi values
      INTEGER :: numPhi
      ! Array of Phi values from 0 to 1
      REAL(8), DIMENSION(:), ALLOCATABLE :: setPhi
      ! Max depth used for calculating lookup table
      REAL(8):: maxLookupDepth
      ! Min depth used for calculating lookup table
      REAL(8) :: minLookupDepth
      ! Array of element numbers in the subgrid area
      INTEGER, DIMENSION(:), ALLOCATABLE :: subgridEleList
      ! Depths corresponding to each phi value on the element
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: wetFracDepthEle
      ! Grid Averaged Total Water Depths corresponding to each phi valu
      ! ! on the element
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: gridDepthEle
      ! Area of sub elements
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: area
      ! Advection corrections corresponding to each phi value on the
      ! element
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: cadvEle
      ! averaged bottom friction coefficient without level 1 correction
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: cfVertTab
      ! averaged bottom friction coefficient with level 1 correction
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: cmfVertTab
      ! Grid Averaged Total Water Depths corresponding to each phi valu
      ! on the vertex
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: gridDepthVertTab
      ! Wet Averaged Total Water Depths corresponding to each phi value
      ! on the vertex
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: wetDepthVertTab
      ! Depths corresponding to each phi value on the vertex
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: wetFracVertTab
      ! Array of vertex numbers in the subgrid area
      INTEGER, DIMENSION(:), ALLOCATABLE :: subgridVertList
      !! Array of elemental wet fractions for zeta at TS k
      !REAL(8), DIMENSION(:,:), ALLOCATABLE :: wetFracEleETA1
      ! Array of elemental wet fractions for zeea at TS k+1
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: wetFracEleETA2
      ! Array of elemental grid depths for zeta at k
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: gridDepthEleETA1
      ! Array of elemental grid depths for zeta at k+1
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: gridDepthEleETA2
      !! Array of elemental advection corrections for zeta at k
      !REAL(8), DIMENSION(:,:), ALLOCATABLE :: cadvEleETA1
      ! Array of elemental advection corrections for zeta at k+1
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: cadvEleETA2
      ! Array of elemental average wet fractions for the entire element
      ! for zeta at k
      REAL(8), DIMENSION(:), ALLOCATABLE :: wetFracAvgETA1
      ! Array of elemental average wet fractions for the entire element
      ! for zeta at k+1
      REAL(8), DIMENSION(:), ALLOCATABLE :: wetFracAvgETA2
      !! Array to hold average elemental grid depths
      !REAL(8), DIMENSION(:), ALLOCATABLE :: HGAvgETA1
      !! Array to hold average elemental grid depths
      !REAL(8), DIMENSION(:), ALLOCATABLE :: HGAvgETA2
      ! Array to hold wet/dry state of element
      INTEGER, DIMENSION(:), ALLOCATABLE :: NCEleArray
      ! Array to hold vertex averged wet fractions for zeta at k+1
      REAL(8), DIMENSION(:), ALLOCATABLE :: wetFracVertETA2
      ! Array to hold vertex grid water depths for zeta at k+1
      REAL(8), DIMENSION(:), ALLOCATABLE :: gridDepthVertETA2
      ! Array to hold vertex wet water depths for zeta at k+1
      REAL(8), DIMENSION(:), ALLOCATABLE :: wetDepthVertETA2
      ! Array to hold vertex bottom friction coefficients for zeta at
      ! k+1
      REAL(8), DIMENSION(:), ALLOCATABLE :: cfVertETA2
      ! Array to hold vertex corrected bottom friction coefficients for
      ! zeta at k+1
      REAL(8), DIMENSION(:), ALLOCATABLE :: cmfVertETA2
      ! Minimum depth in depth array used to preprocess lookup table 
      REAL(8) :: minDepth
      ! Maximum depth in depth array used to preprocess lookup table 
      REAL(8) :: maxDepth

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
      IF(ALLOCATED(wetFracVertETA2))DEALLOCATE(wetFracVertETA2)
      ALLOCATE(wetFracVertETA2(NP))
      IF(ALLOCATED(gridDepthVertETA2))DEALLOCATE(gridDepthVertETA2)
      ALLOCATE(gridDepthVertETA2(NP))
      IF(ALLOCATED(wetDepthVertETA2))DEALLOCATE(wetDepthVertETA2)
      ALLOCATE(wetDepthVertETA2(NP))
      IF(level1)THEN
         IF(ALLOCATED(cmfVertETA2))DEALLOCATE(cmfVertETA2)
         ALLOCATE(cmfVertETA2(NP))
      ELSE
         IF(ALLOCATED(cfVertETA2))DEALLOCATE(cfVertETA2)
         ALLOCATE(cfVertETA2(NP))
      ENDIF

!JLW: set intial value of all arrays to 0
      wetFracVertETA2(:) = 0.d0
      gridDepthVertETA2(:) = 0.d0
      wetDepthVertETA2(:) = 0.d0
      cmfVertETA2(:) = 0.d0
      cfVertETA2(:) = 0.d0

!JLW: initialize local arrays for subgrid element variables
      !IF(ALLOCATED(wetFracEleETA1))DEALLOCATE(wetFracEleETA1)
      !ALLOCATE(wetFracEleETA1(NE,3))
      IF(ALLOCATED(wetFracEleETA2))DEALLOCATE(wetFracEleETA2)
      ALLOCATE(wetFracEleETA2(NE,3))
      IF(ALLOCATED(gridDepthEleETA1))DEALLOCATE(gridDepthEleETA1)
      ALLOCATE(gridDepthEleETA1(NE,3))
      IF(ALLOCATED(gridDepthEleETA2))DEALLOCATE(gridDepthEleETA2)
      ALLOCATE(gridDepthEleETA2(NE,3))
      IF(ALLOCATED(wetFracAvgETA1))DEALLOCATE(wetFracAvgETA1)
      ALLOCATE(wetFracAvgETA1(NE))
      IF(ALLOCATED(wetFracAvgETA2))DEALLOCATE(wetFracAvgETA2)
      ALLOCATE(wetFracAvgETA2(NE))
      !IF(ALLOCATED(HGAvgETA1))DEALLOCATE(HGAvgETA1)
      !ALLOCATE(HGAvgETA1(NE))
      !IF(ALLOCATED(HGAvgETA2))DEALLOCATE(HGAvgETA2)
      !ALLOCATE(HGAvgETA2(NE))
      IF(ALLOCATED(NCEleArray))DEALLOCATE(NCEleArray)
      ALLOCATE(NCEleArray(NE))
      IF(level1)THEN
         !IF(ALLOCATED(cadvEleETA1))DEALLOCATE(cadvEleETA1)
         !ALLOCATE(cadvEleETA1(NE,3))
         IF(ALLOCATED(cadvEleETA2))DEALLOCATE(cadvEleETA2)
         ALLOCATE(cadvEleETA2(NE,3))
      ENDIF

!JLW: set intial value of all arrays to 0
      !wetFracEleETA1(:,:) = 0.d0
      wetFracEleETA2(:,:) = 0.d0
      gridDepthEleETA1(:,:) = 0.d0
      gridDepthEleETA2(:,:) = 0.d0
      wetFracAvgETA1(:) = 0.d0
      wetFracAvgETA2(:) = 0.d0
      !HGAvgETA1(:) = 0.d0
      !HGAvgETA2(:) = 0.d0
      NCEleArray(:) = 0
      IF(level1)THEN
         !cadvEleETA1(:,:) = 0.d0
         cadvEleETA2(:,:) = 0.d0
      ENDIF

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
#ifdef ADCNETCDF          
      use netcdf
      use global, only: imap_el_lg, nodes_lg, allMessage, ERROR
      
      implicit none

      ! Total number of elements and vertices in mesh (can most
      ! definitely be replaced by variables already in ADCIRC) 
      integer :: numElems, numVerts, numVertsPerEle
      ! Used for transfering from global to local
      integer :: currEle, currNode
      ! Arrays for holding global variables
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: globalEleLookup
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: globalVertLookup
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: globalarea
      INTEGER, DIMENSION(:), ALLOCATABLE :: globalSubEleList
      INTEGER, DIMENSION(:), ALLOCATABLE :: globalSubVertList
      ! variables needed to read in main NETCDF lookup table
      INTEGER :: NC_ERR, IS, NC_ID, NC_VAR, NC_DIM_ELEM, NC_DIM_VERT
      INTEGER :: I
      LOGICAL :: file_exists

      INQUIRE(FILE=TRIM(subgridFilename), EXIST=file_exists)
      IF(.NOT.file_exists)THEN
          call allMessage(ERROR, "subgrid lookup file does not exist")
          call terminate()
          call exit(1) 
      ENDIF

!JLW: adding subgrid lookup table read in
!JLW: first open subgrid lookup table and read in dimensions
      CALL CHECK(NF90_OPEN(TRIM(subgridFilename),NF90_NOWRITE,NC_ID))
      
      CALL CHECK(NF90_INQ_DIMID(NC_ID,"numEle",NC_VAR))
      CALL CHECK(NF90_INQUIRE_DIMENSION(NC_ID,NC_VAR,len=numElems))
      
      CALL CHECK(NF90_INQ_DIMID(NC_ID,"numNode",NC_VAR))
      CALL CHECK(NF90_INQUIRE_DIMENSION(NC_ID,NC_VAR,len=numVerts))
      
      CALL CHECK(NF90_INQ_DIMID(NC_ID,"numSfcElevs",NC_VAR))
      CALL CHECK(NF90_INQUIRE_DIMENSION(NC_ID,NC_VAR,len=numSurLevs))
      
      CALL CHECK(NF90_INQ_DIMID(NC_ID,"numPhi",NC_VAR))
      CALL CHECK(NF90_INQUIRE_DIMENSION(NC_ID,NC_VAR,len=numPhi))
!JLW: go through elemental arrays first
!JLW: allocate arrays to fill
      IF(ALLOCATED(wetFracDepthEle))DEALLOCATE(wetFracDepthEle)
      ALLOCATE(wetFracDepthEle(NE,3,numPhi))
      IF(ALLOCATED(gridDepthEle))DEALLOCATE(gridDepthEle)
      ALLOCATE(gridDepthEle(NE,3,numPhi))
      IF(ALLOCATED(cadvEle))DEALLOCATE(cadvEle)
      ALLOCATE(cadvEle(NE,3,numPhi))
      IF(ALLOCATED(area)) DEALLOCATE(area)
      ALLOCATE(area(NE,3))
      IF(ALLOCATED(setPhi)) DEALLOCATE(setPhi)
      ALLOCATE(setPhi(numPhi))
      IF(ALLOCATED(subgridEleList)) DEALLOCATE(subgridEleList)
      ALLOCATE(subgridEleList(NE))

!!!!!!!!!!!!!!!! MINIMUM AND MAXIMUM PREPROCESS DEPTHS !!!!!!!!!!!!!!!
      CALL Check(NF90_INQ_VARID(NC_ID,"minSurfElev",NC_VAR))
      CALL Check(NF90_GET_VAR(NC_ID,NC_VAR,minDepth))
      CALL Check(NF90_INQ_VARID(NC_ID,"maxSurfElev",NC_VAR))
      CALL Check(NF90_GET_VAR(NC_ID,NC_VAR,maxDepth))
!!!!!!!!!!!!!!!!!!!!! READ IN SET WET FRACTIONS !!!!!!!!!!!!!!!!!!!!!!
      CALL Check(NF90_INQ_VARID(NC_ID,"phiSet",NC_VAR))
      CALL Check(NF90_GET_VAR(NC_ID,NC_VAR,setPhi))
!!!!!!!!!!!!!!!!!!!!! wet fraction depths !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF(ALLOCATED(globalEleLookup))DEALLOCATE(globalEleLookup)
      ALLOCATE(globalEleLookup(numElems,3,numPhi))
      CALL Check(NF90_INQ_VARID(NC_ID,"wetFractionDepths",NC_VAR))
      CALL Check(NF90_GET_VAR(NC_ID,NC_VAR,globalEleLookup))

      DO I = 1,NE
#ifdef CMPI
         currEle = abs(imap_el_lg(I))
#else
         currEle = I
#endif
         wetFracDepthEle(I,:,:) = globalEleLookup(currEle,:,:)
      ENDDO
!!!!!!!!!!!!!!!!!!!!!!! total water depths !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!JLW: use only one global elemental array to save space
      IF(ALLOCATED(globalEleLookup))DEALLOCATE(globalEleLookup)
      ALLOCATE(globalEleLookup(numElems,3,numPhi))

      CALL Check(NF90_INQ_VARID(NC_ID,"totWatDepth",NC_VAR))
      CALL Check(NF90_GET_VAR(NC_ID,NC_VAR,globalEleLookup))

      DO I = 1,NE
#ifdef CMPI
         currEle = abs(imap_el_lg(I))
#else
         currEle = I
#endif
         gridDepthEle(I,:,:) = globalEleLookup(currEle,:,:)
      ENDDO
!JLW: now deallocate the large array
      IF(ALLOCATED(globalEleLookup))DEALLOCATE(globalEleLookup)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! advection !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF(level1)THEN
!JLW: use only one global elemental array to save space
         IF(ALLOCATED(globalEleLookup))DEALLOCATE(globalEleLookup)
         ALLOCATE(globalEleLookup(numElems,3,numPhi))

         CALL Check(NF90_INQ_VARID(NC_ID,"cadv",NC_VAR))
         CALL Check(NF90_GET_VAR(NC_ID,NC_VAR,globalEleLookup))

         DO I = 1,NE
#ifdef CMPI
            currEle = abs(imap_el_lg(I))
#else
            currEle = I
#endif
            cadvEle(I,:,:) = globalEleLookup(currEle,:,:)
         ENDDO
!JLW: now deallocate the large array
         IF(ALLOCATED(globalEleLookup))DEALLOCATE(globalEleLookup)
      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! areas !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF(ALLOCATED(globalarea))DEALLOCATE(globalarea)
      ALLOCATE(globalarea(3,numElems))

      CALL Check(NF90_INQ_VARID(NC_ID,"area",NC_VAR))
      CALL Check(NF90_GET_VAR(NC_ID,NC_VAR,globalarea))

      DO I = 1,NE
#ifdef CMPI
         currEle = abs(imap_el_lg(I))
#else
         currEle = I
#endif
         area(I,:) = globalarea(:,currEle)
      ENDDO
!JLW: now deallocate the large array
      IF(ALLOCATED(globalarea))DEALLOCATE(globalarea)
!!!!!!!!!!!!!!!!!!!!!!!!!! subgrid region flag !!!!!!!!!!!!!!!!!!!!!!!!
      IF(ALLOCATED(globalSubEleList))DEALLOCATE(globalSubEleList)
      ALLOCATE(globalSubEleList(numElems))

      CALL Check(NF90_INQ_VARID(NC_ID,"binaryElementList",NC_VAR))
      CALL Check(NF90_GET_VAR(NC_ID,NC_VAR,globalSubEleList))

      DO I = 1,NE
#ifdef CMPI
         currEle = abs(imap_el_lg(I))
#else
         currEle = I
#endif
         subgridEleList(I) = globalSubEleList(currEle)
      ENDDO
      IF(ALLOCATED(globalSubEleList)) DEALLOCATE(globalSubEleList)

!JLW: now repeat this process with all of the vertex lookup tables  
      IF(level1)THEN
         IF(ALLOCATED(cmfVertTab)) DEALLOCATE(cmfVertTab)
         ALLOCATE(cmfVertTab(NP,numPhi))
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

#endif ADCNETCDF      

!----------------------------------------------------------------------
!----------------------------------------------------------------------
      end subroutine readSubgridLookup
!----------------------------------------------------------------------
!----------------------------------------------------------------------


!---------------------------------------------------------------------
!---------------------------------------------------------------------
!
!                  GET ELEMENTAL AVERAGED VARIABLES
!
!---------------------------------------------------------------------
!---------------------------------------------------------------------
! AFTER THE GWCE IS CALLED IN THE TIMESTEP FILE, THIS SUBROUTINE IS
! CALLED TO UPDATE THE ELEMENTAL SUBGRID VARIABLES WITH THE NEWLY
! AVAIABLE WATER SURFACE ELEVATIONS.

      subroutine getEleLookup()
      use mesh, only: ne, nm, totalarea, mju, areas, dp

      implicit none
      integer :: I,J,IE
      integer :: NM1,NM2,NM3
      integer :: numGreater

!JLW: update previous timesteps 
!      wetFracEleETA1 = wetFracEleETA2
      gridDepthEleETA1 = gridDepthEleETA2
      wetFracAvgETA1 = wetFracAvgETA2
!      HGAvgETA1 = HGAvgETA2
!      IF(level1)THEN
!         cadvEleETA1 = cadvEleETA2
!      ENDIF
      MJU(:) = 0
      TotalArea(:) = 0.d0

!JLW: now loop through the elements are get subgrid variables from the
! the most recent water levels calculated by the gwce.

      DO IE = 1,NE
         NM1 = NM(IE,1)
         NM2 = NM(IE,2)
         NM3 = NM(IE,3)
!JLW: if the element is within the subgrid area
         IF(subgridEleList(IE).EQ.1)THEN
            DO I = 1,3
!JLW: intitialize counter 
               numGreater = 0
               DO J = 1,numPhi
!JLW: check to see how many depths in the list of wet fraction depths
! are greater than the current water surface elevation
                  IF(ETA2(NM(IE,I)).LT.wetFracDepthEle(IE,I,J))THEN
                     numGreater = numGreater + 1
                  ENDIF
               ENDDO
               IF(numGreater.EQ.numPhi)THEN
!JLW: this means that none of the depths in the array are greater than
!the current depth and this node is dry
                  wetFracEleETA2(IE,I) = setPhi(1)
                  gridDepthEleETA2(IE,I) = gridDepthEle(IE,I,1)
                  IF(level1)THEN
                     cadvEleETA2(IE,I) = cadvEle(IE,I,1)
                  ENDIF
               ELSEIF(numGreater.EQ.0)THEN
!JLW: this means that none of the depths in the array are greater than
!the current depth and this node is always wet
                  wetFracEleETA2(IE,I) = setPhi(numPhi)   
                  gridDepthEleETA2(IE,I) = gridDepthEle(IE,I,numPhi)
     &               + (ETA2(NM(IE,I)) - wetFracDepthEle(IE,I,numPhi))
                  IF(level1)THEN
                     cadvEleETA2(IE,I) = cadvEle(IE,I,numPhi)
                  ENDIF
               ELSEIF(wetFracDepthEle(IE,I,numPhi-numGreater)
     &           .NE.minLookupDepth)THEN
!JLW: as long as the depth we are indexing is not equal to the minimum
!depth that we used to calculate the subgrid lookup table we interpolate
!between depths to find subgrid values
                  wetFracEleETA2(IE,I) = ((ETA2(NM(IE,I))
     &               - wetFracDepthEle(IE,I,numPhi-numGreater))
     &               /(wetFracDepthEle(IE,I,numPhi-numGreater+1)
     &               - wetFracDepthEle(IE,I,numPhi-numGreater))
     &               *(setPhi(numPhi-numGreater+1)
     &               - setPhi(numPhi-numGreater))
     &               + setPhi(numPhi-numGreater))
                  gridDepthEleETA2(IE,I) = ((ETA2(NM(IE,I))
     &               - wetFracDepthEle(IE,I,numPhi-numGreater))
     &               /(wetFracDepthEle(IE,I,numPhi-numGreater+1)
     &               - wetFracDepthEle(IE,I,numPhi-numGreater))
     &               *(gridDepthEle(IE,I,numPhi-numGreater+1)
     &               - gridDepthEle(IE,I,numPhi-numGreater))
     &               + gridDepthEle(IE,I,numPhi-numGreater))
                  IF(level1)THEN
                     cadvEleETA2(IE,I) = ((ETA2(NM(IE,I))
     &                  - wetFracDepthEle(IE,I,numPhi-numGreater))
     &                  /(wetFracDepthEle(IE,I,numPhi-numGreater+1)
     &                  - wetFracDepthEle(IE,I,numPhi-numGreater))
     &                  *(cadvEle(IE,I,numPhi-numGreater+1)
     &                  - cadvEle(IE,I,numPhi-numGreater))
     &                  + cadvEle(IE,I,numPhi-numGreater))
                  ENDIF
               ELSE
!JLW: if the current water level is equal to the first depth of the
!array wet set our averaged variables to fully wet
                  wetFracEleETA2(IE,I) = setPhi(numPhi-numGreater+1)
                  gridDepthEleETA2(IE,I) = 
     &               gridDepthEle(IE,I,numPhi-numGreater+1)
                  IF(level1)THEN
                     cadvEleETA2(IE,I) =
     &                  cadvEle(IE,I,numPhi-numGreater+1)
                  ENDIF
               ENDIF
            ENDDO
!JLW: calculate the average wet area fraction for the element
            wetFracAvgETA2(IE) = (wetFracEleETA2(IE,1)
     &         + wetFracEleETA2(IE,2)
     &         + wetFracEleETA2(IE,3))/3.d0
         ELSE
!JLW: if the element is not contained within the subgrid area
            wetFracAvgETA2(IE) = 1
            ! NEED TO MOVE THIS TO wet/dry but just for areas outside
            ! the subgrid area since I need the updated node codes
            !NCEle = NODECODE(NM1)*NODECODE(NM2)*NODECODE(NM3)
            DO I = 1,3
!JLW: calcualte the water depth from bathymetry of mesh for areas
!outside the subgrid
               gridDepthEleETA2(IE,I) = DP(NM(IE,I))
     &            + IFNLFA*ETA2(NM(IE,I))
!JLW: set cadv to 1 for areas outside the subgrid so that it does not
!contribute to the calculations 
               cadvEleETA2(IE,I) = 1.d0
            ENDDO
         ENDIF
!JLW: calculate the average total water depth for the element
!         HGAvgETA2(IE) = (gridDepthEleETA2(IE,1) 
!     &      + gridDepthEleETA2(IE,2)
!     &      + gridDepthEleETA2(IE,3))/3.d0
      ENDDO

 
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      end subroutine getEleLookup
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

      use sizes, only : mne, mnp
      use mesh, only : ne, np, dp, mju, totalArea, nm, x, y, areas
     &   , NNeigh, NeiTab

      implicit none
      integer :: ie, i, j
      integer :: numGreater
      REAL(8) :: H2

!JLW: set all temporary nodecodes to 0
      NNODECODE(:) = 0
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
               ELSE
                  cfVertETA2(i) = cfVertTab(i,numPhi-numGreater+1)
               ENDIF
            ENDIF
!!JLW: get water depth to determine nnodecode
!            H2 = gridDepthVertETA2(i)
!         ELSE
!            H2 = dp(i) + eta2(i) 
         ENDIF
!!JLW: set nnodecode based on water depth
!         IF(H2.GT.H0)THEN
!            nnodecode(i) = 1
!         ELSE
!            nnodecode(i) = 0
!         ENDIF
      ENDDO
!!JLW: update nodecode
!      nodecode = nnodecode
!!JLW: make sure to activate ncchange
!      ncchange = 1
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      end subroutine getVertLookup
!----------------------------------------------------------------------
!----------------------------------------------------------------------


!---------------------------------------------------------------------
!---------------------------------------------------------------------
!
!                  WET DRY CHECK AND UPDATE
!
!---------------------------------------------------------------------
!---------------------------------------------------------------------
! AFTER THE ELEMENTAL SUBGRID ROUTINE AND THE VERTEX SUBROUTINE ARE
! CALLED AND THE SUBGRID ARRAYS ARE UPDATED, THIS WET DRY CODE IS CALLED
! TO UPDATE THE WET/DRY STATE OF THE ELEMENTS AND VERTICES. THIS
! CONSISTS OF A SIMPLE CHECK OF TOTAL WATER DEPTH H. IN ADDITION,
! UPDATED COEFFICIENTS OF FRICTION AND AREAS ARE ALSO UPDATED.

      subroutine subgridwetdry(BFCdLLimit, fric, fgamma, ftheta, 
     &                         manningsn, hbreak,
     &                         iflinbf, ifnlbf, ifhybf)

!      use NodalAttributes, only: BFCdLLimit, fric
!     &   , fgamma, ftheta, manningsn, hbreak, loadManningsN
!     &   , iflinbf, ifnlbf, ifhybf 

      use constants, only: g
      use NODAL_ATTRIBUTE_TYPE_MODULE, ONLY: NODAL_ATTRIBUTE_VECTOR_REAL_T
      IMPLICIT NONE

      REAL(8), INTENT(IN)    :: BFCdLLimit
      REAL(8), INTENT(IN)    :: fgamma
      REAL(8), INTENT(IN)    :: ftheta
      REAL(8), INTENT(IN)    :: hbreak
      INTEGER, INTENT(IN)    :: iflinbf
      INTEGER, INTENT(IN)    :: ifnlbf
      INTEGER, INTENT(IN)    :: ifhybf
      REAL(8), INTENT(INOUT) :: fric(*)
      TYPE(NODAL_ATTRIBUTE_VECTOR_REAL_T) :: ManningsN

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
            IF(ManningsN%active)THEN
               IF(subgridVertList(i).EQ.1)THEN
                  IF(level1)THEN
                     fric(i) = cmfVertETA2(i)
                  ELSE
                     fric(i) = cfVertETA2(i)
                  ENDIF
               ELSE
                  fric(i) = g*manningsn%values(i)**2.d0
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
!                 READ NETCDF UTILITY
!
!---------------------------------------------------------------------
!---------------------------------------------------------------------

!JLW: add subroutine for reading in netcdf
#ifdef ADCNETCDF      
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
#endif
!---------------------------------------------------------------------
!---------------------------------------------------------------------

      SUBROUTINE terminate(NO_MPI_FINALIZE)
#ifdef CMPI
      USE MESSENGER
#endif
      USE GLOBAL, ONLY: INFO,allMessage,
     &                  setMessageSource,unsetMessageSource
      IMPLICIT NONE
      LOGICAL, OPTIONAL :: NO_MPI_FINALIZE
      call setMessageSource("terminate")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      call allMessage(INFO,"ADCIRC Terminating.")

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
#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.") ! should be unreachable
#endif
      call unsetMessageSource()
      END SUBROUTINE terminate

      end module subgrid
