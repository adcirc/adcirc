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
#ifdef ADCNETCDF      
      use netcdf_error, only: check_err
#endif      

! JLW: initialize the arrays for use in the subgrid
! code
      implicit none
      CHARACTER(200) :: subgridFilename
      ! Flag for subgrid corrections to wetting and drying
      LOGICAL :: level0
      ! Flag for subgrid corrections to bottom friction and advection
      LOGICAL :: level1
      ! Number of surface elevations used in subgrid preprocessor
!      INTEGER :: numSurLevs
!      ! Number of possible phi values
      INTEGER :: numPhi
      ! Array of Phi values from 0 to 1
      REAL(8), DIMENSION(:), ALLOCATABLE :: setPhi
      ! averaged bottom friction coefficient without level 1 correction
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: cfVertTab
      ! averaged bottom friction coefficient with level 1 correction
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: cmfVertTab
      ! averaged advection with level 1 correction
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: cadvVertTab
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
      ! Array to hold vertex averged wet fractions for zeta at k
      REAL(8), DIMENSION(:), ALLOCATABLE :: wetFracVertETA1
      ! Array to hold vertex averged wet fractions for zeta at k+1
      REAL(8), DIMENSION(:), ALLOCATABLE :: wetFracVertETA2
      ! Array to hold vertex grid water depths for zeta at k
      REAL(8), DIMENSION(:), ALLOCATABLE :: gridDepthVertETA1
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
      ! Array to hold vertex corrected advection for
      ! zeta at k+1
      REAL(8), DIMENSION(:), ALLOCATABLE :: cadvVertETA2
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

!JLW: set intial value of all arrays to 0
      wetFracVertETA1(:) = 0.d0
      wetFracVertETA2(:) = 0.d0
      gridDepthVertETA1(:) = 0.d0
      gridDepthVertETA2(:) = 0.d0
      wetDepthVertETA2(:) = 0.d0

      IF(level1)THEN
         IF(ALLOCATED(cmfVertETA2))DEALLOCATE(cmfVertETA2)
         ALLOCATE(cmfVertETA2(NP))
         IF(ALLOCATED(cadvVertETA2))DEALLOCATE(cadvVertETA2)
         ALLOCATE(cadvVertETA2(NP))
         cmfVertETA2(:) = 0.d0
         cadvVertETA2(:) = 1.d0
!JLW: will need to add vertex based advection
      ELSE
         IF(ALLOCATED(cfVertETA2))DEALLOCATE(cfVertETA2)
         ALLOCATE(cfVertETA2(NP))
         cfVertETA2(:) = 0.d0
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
      ! use global, only: imap_el_lg, nodes_lg, allMessage, ERROR
      use global, only: nodes_lg, allMessage, ERROR
      
      implicit none

      ! Total number of elements and vertices in mesh (can most
      ! definitely be replaced by variables already in ADCIRC) 
      !integer :: numElems, numVerts, numVertsPerEle
      integer :: numVerts
      ! Used for transfering from global to local
      ! integer :: currEle, currNode
      integer :: currNode
      ! Arrays for holding global variables
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: globalVertLookup
      INTEGER, DIMENSION(:), ALLOCATABLE :: globalSubVertList
      ! variables needed to read in main NETCDF lookup table
      INTEGER :: NC_ERR, IS, NC_ID, NC_VAR 
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
      CALL CHECK_ERR(NF90_OPEN(TRIM(subgridFilename),NF90_NOWRITE,NC_ID))
      
      CALL CHECK_ERR(NF90_INQ_DIMID(NC_ID,"numNode",NC_VAR))
      CALL CHECK_ERR(NF90_INQUIRE_DIMENSION(NC_ID,NC_VAR,len=numVerts))
      
      CALL CHECK_ERR(NF90_INQ_DIMID(NC_ID,"numPhi",NC_VAR))
      CALL CHECK_ERR(NF90_INQUIRE_DIMENSION(NC_ID,NC_VAR,len=numPhi))

      IF(ALLOCATED(setPhi)) DEALLOCATE(setPhi)
      ALLOCATE(setPhi(numPhi))

      CALL check_err(NF90_INQ_VARID(NC_ID,"phiSet",NC_VAR))
      CALL check_err(NF90_GET_VAR(NC_ID,NC_VAR,setPhi))

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
         CALL check_err(NF90_INQ_VARID(NC_ID,"cmfVertex",NC_VAR))
         CALL check_err(NF90_GET_VAR(NC_ID,NC_VAR,globalVertLookup))
         DO I = 1,NP
#ifdef CMPI
            currNode = abs(nodes_lg(I))
#else
            currNode = I
#endif
            cmfVertTab(I,:) = globalVertLookup(:,currNode)
         ENDDO
      ELSE
         CALL check_err(NF90_INQ_VARID(NC_ID,"cfVertex",NC_VAR))
         CALL check_err(NF90_GET_VAR(NC_ID,NC_VAR,globalVertLookup))
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
         CALL check_err(NF90_INQ_VARID(NC_ID,"cadvVertex",NC_VAR))
         CALL check_err(NF90_GET_VAR(NC_ID,NC_VAR,globalVertLookup))
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

      CALL check_err(NF90_INQ_VARID(NC_ID,"gridTotWatDepthVertex",NC_VAR))
      CALL check_err(NF90_GET_VAR(NC_ID,NC_VAR,globalVertLookup))
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

      CALL check_err(NF90_INQ_VARID(NC_ID,"wetTotWatDepthVertex",NC_VAR))
      CALL check_err(NF90_GET_VAR(NC_ID,NC_VAR,globalVertLookup))
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

      CALL check_err(NF90_INQ_VARID(NC_ID,"wetFractionVertex",NC_VAR))
      CALL check_err(NF90_GET_VAR(NC_ID,NC_VAR,globalVertLookup))
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

      CALL check_err(NF90_INQ_VARID(NC_ID,"binaryVertexList",NC_VAR))
      CALL check_err(NF90_GET_VAR(NC_ID,NC_VAR,globalSubVertList))

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
!                 READ NETCDF UTILITY
!
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
