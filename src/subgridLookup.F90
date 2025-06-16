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
   use global, only: nodecode, nnodecode &
                     , eta2, tk, tkm, nolifa, tkm &
                     , H0, ifnlfa, uu1, vv1, ncchange
   use mesh, only: ne, np, nm, totalArea, mju, areas, dp
#ifdef ADCNETCDF
   use netcdf_error, only: check_err
#endif

   ! JLW: initialize the arrays for use in the subgrid
   ! code
   implicit none
   character(200) :: subgridFilename = 'null'
   ! Flag for subgrid corrections to wetting and drying
   logical :: level0
   ! Flag for subgrid corrections to bottom friction and advection
   logical :: level1
   !JLW: adding logical for dphidt to be turned off - set default to off
   logical :: phi_time_derivative = .false.
   integer :: dphidt ! used to 0 out the term
   ! Number of possible phi values
   integer :: numPhi
   ! Array of Phi values from 0 to 1
   real(8), dimension(:), allocatable :: setPhi
   ! averaged bottom friction coefficient without level 1 correction
   real(8), dimension(:, :), allocatable :: cfVertTab
   ! averaged bottom friction coefficient with level 1 correction
   real(8), dimension(:, :), allocatable :: cmfVertTab
   ! averaged advection with level 1 correction
   real(8), dimension(:, :), allocatable :: cadvVertTab
   ! Grid Averaged Total Water Depths corresponding to each phi valu
   ! on the vertex
   real(8), dimension(:, :), allocatable :: gridDepthVertTab
   ! Wet Averaged Total Water Depths corresponding to each phi value
   ! on the vertex
   real(8), dimension(:, :), allocatable :: wetDepthVertTab
   ! Depths corresponding to each phi value on the vertex
   real(8), dimension(:, :), allocatable :: wetFracVertTab
   ! Array of vertex numbers in the subgrid area
   integer, dimension(:), allocatable :: subgridVertList
   ! Array to hold vertex averged wet fractions for zeta at k
   real(8), dimension(:), allocatable :: wetFracVertETA1
   ! Array to hold vertex averged wet fractions for zeta at k+1
   real(8), dimension(:), allocatable :: wetFracVertETA2
   ! Array to hold vertex grid water depths for zeta at k
   real(8), dimension(:), allocatable :: gridDepthVertETA1
   ! Array to hold vertex grid water depths for zeta at k+1
   real(8), dimension(:), allocatable :: gridDepthVertETA2
   ! Array to hold vertex wet water depths for zeta at k+1
   real(8), dimension(:), allocatable :: wetDepthVertETA2
   ! Array to hold vertex bottom friction coefficients for zeta at
   ! k+1
   real(8), dimension(:), allocatable :: cfVertETA2
   ! Array to hold vertex corrected bottom friction coefficients for
   ! zeta at k+1
   real(8), dimension(:), allocatable :: cmfVertETA2
   ! Array to hold vertex corrected advection for
   ! zeta at k+1
   real(8), dimension(:), allocatable :: cadvVertETA2
   ! Minimum depth in depth array used to preprocess lookup table
   real(8) :: minDepth
   ! Maximum depth in depth array used to preprocess lookup table
!   real(8) :: maxDepth

   private

   public :: initVar, getVertLookup, readSubgridLookup, level0, level1, dphidt, &
             phi_time_derivative, cmfVertETA2, cfVertETA2, subgridVertList, &
             wetDepthVertETA2, wetFracVertETA1, wetFracVertETA2, cadvVertETA2, &
             gridDepthVertETA1, gridDepthVertETA2, wetFracVertTab, wetDepthVertTab, &
             gridDepthVertTab, numPhi, subgridFilename

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
      if (allocated(wetFracVertETA1)) deallocate (wetFracVertETA1)
      allocate (wetFracVertETA1(NP))
      if (allocated(gridDepthVertETA1)) deallocate (gridDepthVertETA1)
      allocate (gridDepthVertETA1(NP))
      if (allocated(wetFracVertETA2)) deallocate (wetFracVertETA2)
      allocate (wetFracVertETA2(NP))
      if (allocated(gridDepthVertETA2)) deallocate (gridDepthVertETA2)
      allocate (gridDepthVertETA2(NP))
      if (allocated(wetDepthVertETA2)) deallocate (wetDepthVertETA2)
      allocate (wetDepthVertETA2(NP))

      !JLW: set intial value of all arrays to 0
      wetFracVertETA1(:) = 0.d0
      wetFracVertETA2(:) = 0.d0
      gridDepthVertETA1(:) = 0.d0
      gridDepthVertETA2(:) = 0.d0
      wetDepthVertETA2(:) = 0.d0

      if (level1) then
         if (allocated(cmfVertETA2)) deallocate (cmfVertETA2)
         allocate (cmfVertETA2(NP))
         if (allocated(cadvVertETA2)) deallocate (cadvVertETA2)
         allocate (cadvVertETA2(NP))
         cmfVertETA2(:) = 0.d0
         cadvVertETA2(:) = 1.d0
         !JLW: will need to add vertex based advection
      else
         if (allocated(cfVertETA2)) deallocate (cfVertETA2)
         allocate (cfVertETA2(NP))
         cfVertETA2(:) = 0.d0
      end if

   end subroutine initVar
   !----------------------------------------------------------------------

   !----------------------------------------------------------------------
   !
   !                  SUBROUTINE READ IN NETCDF FILE
   !
   !----------------------------------------------------------------------
   ! READS IN THE NETCDF LOOKUP TABLE AND POPULATES THE SUBGRID ARRAYS BOTH
   ! FOR SERIAL AND PARALLEL ADCIRC. THIS SUBROUTINE USES A FUNCTION
   ! "CHECK" TO READ IN THE NETCDF FILE. THE "CHECK" SUBROUTINE IS LOCATED
   ! AT THE BOTTOM OF THIS FILE.

   subroutine readSubgridLookup
#ifdef ADCNETCDF
      use netcdf, only: NF90_OPEN, NF90_NOWRITE, NF90_INQ_DIMID, &
                        NF90_INQUIRE_DIMENSION, NF90_INQ_VARID, &
                        NF90_GET_VAR, NF90_CLOSE
      use global, only: allMessage, ERROR
#ifdef CMPI
      use global, only: nodes_lg
#endif

      implicit none

      ! Total number of elements and vertices in mesh (can most
      ! definitely be replaced by variables already in ADCIRC)
      integer :: numVerts
      ! Used for transfering from global to local
      ! integer :: currEle, currNode
      integer :: currNode
      ! Arrays for holding global variables
      real(8), dimension(:, :), allocatable :: globalVertLookup
      integer, dimension(:), allocatable :: globalSubVertList
      ! variables needed to read in main NETCDF lookup table
      integer :: NC_ERR, NC_ID, NC_VAR
      integer :: I
      logical :: file_exists

      inquire (FILE=trim(subgridFilename), EXIST=file_exists)
      if (.not. file_exists) then
         call allMessage(ERROR, "subgrid lookup file does not exist")
         call terminate()
         call exit(1)
      end if

!JLW: adding subgrid lookup table read in
!JLW: first open subgrid lookup table and read in dimensions
      call CHECK_ERR(NF90_OPEN(trim(subgridFilename), NF90_NOWRITE, NC_ID))

      call CHECK_ERR(NF90_INQ_DIMID(NC_ID, "numNode", NC_VAR))
      call CHECK_ERR(NF90_INQUIRE_DIMENSION(NC_ID, NC_VAR, len=numVerts))

      call CHECK_ERR(NF90_INQ_DIMID(NC_ID, "numPhi", NC_VAR))
      call CHECK_ERR(NF90_INQUIRE_DIMENSION(NC_ID, NC_VAR, len=numPhi))

      if (allocated(setPhi)) deallocate (setPhi)
      allocate (setPhi(numPhi))

      call check_err(NF90_INQ_VARID(NC_ID, "phiSet", NC_VAR))
      call check_err(NF90_GET_VAR(NC_ID, NC_VAR, setPhi))

!JLW: go through and read in each of the vertex averaged variables
!JLW: allocate arrays
      if (level1) then
         if (allocated(cmfVertTab)) deallocate (cmfVertTab)
         allocate (cmfVertTab(NP, numPhi))
         if (allocated(cadvVertTab)) deallocate (cadvVertTab)
         allocate (cadvVertTab(NP, numPhi))
      else
         if (allocated(cfVertTab)) deallocate (cfVertTab)
         allocate (cfVertTab(NP, numPhi))
      end if
      if (allocated(gridDepthVertTab)) deallocate (gridDepthVertTab)
      allocate (gridDepthVertTab(NP, numPhi))
      if (allocated(wetDepthVertTab)) deallocate (wetDepthVertTab)
      allocate (wetDepthVertTab(NP, numPhi))
      if (allocated(wetFracVertTab)) deallocate (wetFracVertTab)
      allocate (wetFracVertTab(NP, numPhi))
      if (allocated(subgridVertList)) deallocate (subgridVertList)
      allocate (subgridVertList(NP))

!!!!!!!!!!!!!!!!!! bottom friction coefficient !!!!!!!!!!!!!!!!!!!!!!
      if (allocated(globalVertLookup)) deallocate (globalVertLookup)
      allocate (globalVertLookup(numPhi, numVerts))

      if (level1) then
         call check_err(NF90_INQ_VARID(NC_ID, "cmfVertex", NC_VAR))
         call check_err(NF90_GET_VAR(NC_ID, NC_VAR, globalVertLookup))
         do I = 1, NP
#ifdef CMPI
            currNode = abs(nodes_lg(I))
#else
            currNode = I
#endif
            cmfVertTab(I, :) = globalVertLookup(:, currNode)
         end do
      else
         call check_err(NF90_INQ_VARID(NC_ID, "cfVertex", NC_VAR))
         call check_err(NF90_GET_VAR(NC_ID, NC_VAR, globalVertLookup))
         do I = 1, NP
#ifdef CMPI
            currNode = abs(nodes_lg(I))
#else
            currNode = I
#endif
            cfVertTab(I, :) = globalVertLookup(:, currNode)
         end do
      end if
!!!!!!!!!!!!!!!!!! ADVECTION CORRECTION VETEX !!!!!!!!!!!!!!!!!!
      if (allocated(globalVertLookup)) deallocate (globalVertLookup)
      allocate (globalVertLookup(numPhi, numVerts))
      if (level1) then
         call check_err(NF90_INQ_VARID(NC_ID, "cadvVertex", NC_VAR))
         call check_err(NF90_GET_VAR(NC_ID, NC_VAR, globalVertLookup))
         do I = 1, NP
#ifdef CMPI
            currNode = abs(nodes_lg(I))
#else
            currNode = I
#endif
            cadvVertTab(I, :) = globalVertLookup(:, currNode)
         end do
      end if
!!!!!!!!!!!!!!!!! grid depth vertex !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (allocated(globalVertLookup)) deallocate (globalVertLookup)
      allocate (globalVertLookup(numPhi, numVerts))

      call check_err(NF90_INQ_VARID(NC_ID, "gridTotWatDepthVertex", NC_VAR))
      call check_err(NF90_GET_VAR(NC_ID, NC_VAR, globalVertLookup))
      do I = 1, NP
#ifdef CMPI
         currNode = abs(nodes_lg(I))
#else
         currNode = I
#endif
         gridDepthVertTab(I, :) = globalVertLookup(:, currNode)
      end do
!!!!!!!!!!!!!!!!! wet depth vertex !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (allocated(globalVertLookup)) deallocate (globalVertLookup)
      allocate (globalVertLookup(numPhi, numVerts))

      call check_err(NF90_INQ_VARID(NC_ID, "wetTotWatDepthVertex", NC_VAR))
      call check_err(NF90_GET_VAR(NC_ID, NC_VAR, globalVertLookup))
      do I = 1, NP
#ifdef CMPI
         currNode = abs(nodes_lg(I))
#else
         currNode = I
#endif
         wetDepthVertTab(I, :) = globalVertLookup(:, currNode)
      end do
!!!!!!!!!!!!!!!!!!!! wet fraction depths vertex !!!!!!!!!!!!!!!!!!!!!!
      if (allocated(globalVertLookup)) deallocate (globalVertLookup)
      allocate (globalVertLookup(numPhi, numVerts))

      call check_err(NF90_INQ_VARID(NC_ID, "wetFractionVertex", NC_VAR))
      call check_err(NF90_GET_VAR(NC_ID, NC_VAR, globalVertLookup))
      do I = 1, NP
#ifdef CMPI
         currNode = abs(nodes_lg(I))
#else
         currNode = I
#endif
         wetFracVertTab(I, :) = globalVertLookup(:, currNode)
      end do
      if (allocated(globalVertLookup)) deallocate (globalVertLookup)
!!!!!!!!!!!!!!!!! subgrid flag for vertex !!!!!!!!!!!!!!!!!!!!!!!!
      if (allocated(globalSubVertList)) deallocate (globalSubVertList)
      allocate (globalSubVertList(numVerts))

      call check_err(NF90_INQ_VARID(NC_ID, "binaryVertexList", NC_VAR))
      call check_err(NF90_GET_VAR(NC_ID, NC_VAR, globalSubVertList))

      do I = 1, NP
#ifdef CMPI
         currNode = abs(nodes_lg(I))
#else
         currNode = I
#endif
         subgridVertList(I) = globalSubVertList(currNode)
      end do
      if (allocated(globalSubVertList)) deallocate (globalSubVertList)

!JLW: close netcdf lookup table file
      NC_ERR = NF90_CLOSE(NC_ID)

#endif

   end subroutine readSubgridLookup
   !----------------------------------------------------------------------

   !---------------------------------------------------------------------
   !
   !                  GET VERTEX AVERAGED VARIABLES
   !
   !---------------------------------------------------------------------
   !---------------------------------------------------------------------
   ! THIS CODE IS VERY SIMILAR TO THE ELEMENTAL SUBROUTINE BUT FOR SUBRID
   ! VARIABLES RELATED TO THE VERTEX AREAS. FOR EACH NEW WATER SURFACE
   ! ELEVATION SUBGRID VERTEX QUANTITIES ARE LOOKED UP.

   subroutine getVertLookup()

      use mesh, only: np, dp
      use global, only: H0

      implicit none
      real(8), parameter :: eps = epsilon(1.d0)
      integer :: i, j
      integer :: numGreater
      real(8) :: HABSMIN

      !JLW: set k = k + 1 for new timestep
      HABSMIN = 0.8d0*H0
      wetFracVertETA1 = wetFracVertETA2
      gridDepthVertETA1 = gridDepthVertETA2
      !JLW: now loop through the vertices and lookup vertex averages
      do I = 1, NP
         if (subgridVertList(I) == 1) then
            !JLW: find how many of the depths are greater than the current wet area
            !fraction
            numGreater = 0 !initialize counter
            do j = 1, numPhi
               if (eta2(i) < wetFracVertTab(i, j)) then
                  numGreater = numGreater + 1
               end if
            end do
            if (numGreater == numPhi) then
               !JLW: this means that none of the depths in the array are greater than
               !the current depth and this node is dry
               wetFracVertETA2(i) = setPhi(1)
               gridDepthVertETA2(i) = gridDepthVertTab(i, 1)
               wetDepthVertETA2(i) = wetDepthVertTab(i, 1)
               if (level1) then
                  cmfVertETA2(i) = cmfVertTab(i, 1)
                  cadvVertETA2(i) = cadvVertTab(i, 1)
               else
                  cfVertETA2(i) = cfVertTab(i, 1)
               end if
            elseif (numGreater == 0) then
               !JLW: this means that none of the depths in the array are greater than
               !the current depth and this node is always wet
               wetFracVertETA2(i) = setPhi(numPhi)
               gridDepthVertETA2(i) = gridDepthVertTab(i, numPhi) &
                                      + (eta2(i) - wetFracVertTab(i, numPhi))
               wetDepthVertETA2(i) = wetDepthVertTab(i, numPhi) &
                                     + (eta2(i) - wetFracVertTab(i, numPhi))
               !JLW 20250615: making it so if a node is fully wet, use traditional
               ! manning's calculations with nodal attribute values
               if (level1) then
                  cmfVertETA2(i) = cmfVertTab(i, numPhi)
                  cadvVertETA2(i) = cadvVertTab(i, numPhi)
               else
                  cfVertETA2(i) = cfVertTab(i, numPhi)
               end if

            elseif (abs(wetFracVertTab(i, numPhi - numGreater) - minDepth) > eps) then
               wetFracVertETA2(i) = ((eta2(i) &
                                      - wetFracVertTab(i, numPhi - numGreater)) &
                                     /(wetFracVertTab(i, numPhi - numGreater + 1) &
                                       - wetFracVertTab(i, numPhi - numGreater)) &
                                     *(setPhi(numPhi - numGreater + 1) &
                                       - setPhi(numPhi - numGreater)) &
                                     + setPhi(numPhi - numGreater))
               gridDepthVertETA2(i) = ((eta2(i) &
                                        - wetFracVertTab(i, numPhi - numGreater)) &
                                       /(wetFracVertTab(i, numPhi - numGreater + 1) &
                                         - wetFracVertTab(i, numPhi - numGreater)) &
                                       *(gridDepthVertTab(i, numPhi - numGreater + 1) &
                                         - gridDepthVertTab(i, numPhi - numGreater)) &
                                       + gridDepthVertTab(i, numPhi - numGreater))
               wetDepthVertETA2(i) = ((eta2(i) &
                                       - wetFracVertTab(i, numPhi - numGreater)) &
                                      /(wetFracVertTab(i, numPhi - numGreater + 1) &
                                        - wetFracVertTab(i, numPhi - numGreater)) &
                                      *(wetDepthVertTab(i, numPhi - numGreater + 1) &
                                        - wetDepthVertTab(i, numPhi - numGreater)) &
                                      + wetDepthVertTab(i, numPhi - numGreater))
               if (level1) then
                  cmfVertETA2(i) = ((eta2(i) &
                                     - wetFracVertTab(i, numPhi - numGreater)) &
                                    /(wetFracVertTab(i, numPhi - numGreater + 1) &
                                      - wetFracVertTab(i, numPhi - numGreater)) &
                                    *(cmfVertTab(i, numPhi - numGreater + 1) &
                                      - cmfVertTab(i, numPhi - numGreater)) &
                                    + cmfVertTab(i, numPhi - numGreater))
                  cadvVertETA2(i) = ((eta2(i) &
                                      - wetFracVertTab(i, numPhi - numGreater)) &
                                     /(wetFracVertTab(i, numPhi - numGreater + 1) &
                                       - wetFracVertTab(i, numPhi - numGreater)) &
                                     *(cadvVertTab(i, numPhi - numGreater + 1) &
                                       - cadvVertTab(i, numPhi - numGreater)) &
                                     + cadvVertTab(i, numPhi - numGreater))
               else
                  cfVertETA2(i) = ((eta2(i) &
                                    - wetFracVertTab(i, numPhi - numGreater)) &
                                   /(wetFracVertTab(i, numPhi - numGreater + 1) &
                                     - wetFracVertTab(i, numPhi - numGreater)) &
                                   *(cfVertTab(i, numPhi - numGreater + 1) &
                                     - cfVertTab(i, numPhi - numGreater)) &
                                   + cfVertTab(i, numPhi - numGreater))
               end if
            else
               !JLW: if the current water level is equal to the first depth of the
               !array wet set our averaged variables to fully wet
               wetFracVertETA2(i) = setPhi(numPhi - numGreater + 1)
               gridDepthVertETA2(i) = &
                  gridDepthVertTab(i, numPhi - numGreater + 1)
               wetDepthVertETA2(i) = &
                  wetDepthVertTab(i, numPhi - numGreater + 1)
               if (level1) then
                  cmfVertETA2(i) = cmfVertTab(i, numPhi - numGreater + 1)
                  cadvVertETA2(i) = cadvVertTab(i, numPhi - numGreater + 1)
               else
                  cfVertETA2(i) = cfVertTab(i, numPhi - numGreater + 1)
               end if
            end if
         else
            !JLW: if not in subgrid area set grid total depth (this was a
            !bug)
            gridDepthVertETA2(i) = eta2(i) + dp(i)
            wetDepthVertETA2(i) = eta2(i) + dp(i)
            wetFracVertETA2(i) = 1.d0
         end if
      end do

   end subroutine getVertLookup
   !----------------------------------------------------------------------

      subroutine terminate()
#ifdef CMPI
      use MESSENGER, only: MSG_FINI, subdomainFatalError
#endif
      use GLOBAL, only: INFO, allMessage, &
                        setMessageSource, unsetMessageSource
#if defined(MESH_TRACE) || defined(ALL_TRACE)
      use GLOBAL, only: DEBUG
#endif
      implicit none
      call setMessageSource("terminate")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      call allMessage(INFO, "ADCIRC Terminating.")

#ifdef CMPI
      subdomainFatalError = .true.
      call MSG_FINI()
#endif
      call exit(1)
      !
#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.") ! should be unreachable
#endif
      call unsetMessageSource()
   end subroutine terminate

end module subgrid
