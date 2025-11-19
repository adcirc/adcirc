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
!----------------------------------------------------------------------
module GL2LOC_MAPPING
   !----------------------------------------------------------------------
   !     This module exists to create a method to more efficiently read in
   !     input data on core 0 and distribute that data to the appropriate
   !     subdomain. It was initially written for use in reading in
   !     netCDF-format global baroclinic forcing data to be used in loose,
   !     one-way coupling to GOFS. My hope is that it will be useful for
   !     others.
   !
   !     Written by: Coleman Blakely 11/2022
   !----------------------------------------------------------------------
   use SIZES, only: MNPROC, MYPROC
   use MESH, only: NP ! local number of nodes
   use GLOBAL, only: nodes_lg, np_g, COMM

   use mpi_f08, only: MPI_INTEGER, MPI_Bcast, MPI_Gather, &
                      MPI_Gatherv, MPI_Scatterv

   implicit none

   private

   ! Start of variables used
   integer, dimension(:), allocatable :: GL2LOC
   ! sending information
   real(8), dimension(:), allocatable :: SENDBUF_REAL
   ! counts of data to send and displacements
   integer, dimension(:), allocatable :: SENDCOUNTS, DISPLS
   ! recieving info
   real(8), dimension(:), allocatable :: RECVBUF_REAL
   integer :: RECVCOUNT

   public :: MAPTOLOCAL_REAL, BcastToLocal_Int, BcastToLocal_2DRealArray

contains

   !----------------------------------------------------------------------
   subroutine MAPTOLOCAL_REAL(GLOBALDATA, LOCALDATA)
      use mpi_f08, only: MPI_DOUBLE_PRECISION
      implicit none
      real(8), intent(IN) :: GLOBALDATA(:)
      real(8), intent(OUT) :: LOCALDATA(:)
      integer :: ierr
      integer :: kk

      ! build the global mapping (GL2LOC) if this is the first call
      call BUILD_GL2LOC()
      ! To avoid allocating/deallocating a bunch just use a dmy array
      ! that has size equal to the maximum number of nodes on a
      ! processor.

      if (MYPROC == 0) then
         ! Put global data into sendbufr in the correct order
         do kk = 1, sum(SENDCOUNTS)
            SENDBUF_REAL(kk) = GLOBALDATA(GL2LOC(kk))
         end do
      end if
      call MPI_SCATTERV(SENDBUF_REAL, SENDCOUNTS, DISPLS, MPI_DOUBLE_PRECISION, &
                        RECVBUF_REAL, RECVCOUNT, MPI_DOUBLE_PRECISION, &
                        0, COMM, ierr)
      LOCALDATA = RECVBUF_REAL

   end subroutine MAPTOLOCAL_REAL

   !----------------------------------------------------------------------
   subroutine BUILD_GL2LOC()
      implicit none

      logical, save :: first_call = .true.
      integer :: kk, ierr
      integer, allocatable :: local_gl2loc(:)

      ! Return if we have already built the table
      if (first_call .eqv. .false.) return
      first_call = .false.
      ! On processor 0, allocate sendcounts and displs
      if (MYPROC == 0) then
         allocate (SENDCOUNTS(MNPROC))
         allocate (DISPLS(MNPROC))
      else
         allocate (sendcounts(1))
         allocate (DISPLS(1))
      end if
      ! get recvcount on all processors
      RECVCOUNT = NP
      ! allocate recvbuf_real
      allocate (RECVBUF_REAL(RECVCOUNT))
      ! build local mapping tables
      allocate (local_gl2loc(RECVCOUNT))
      do kk = 1, RECVCOUNT
         local_gl2loc(kk) = abs(nodes_lg(kk))
      end do
      ! gather all recvcount values on proc 0
      call MPI_GATHER(RECVCOUNT, 1, MPI_INTEGER, &
                      SENDCOUNTS, 1, MPI_INTEGER, &
                      0, COMM, ierr)
      if (MYPROC == 0) then
         DISPLS(1) = 0
         do kk = 2, MNPROC
            DISPLS(kk) = sum(SENDCOUNTS(1:kk - 1))
         end do
         ! allocate gl2loc
         allocate (GL2LOC(sum(SENDCOUNTS)))
         allocate (SENDBUF_REAL(sum(SENDCOUNTS)))
      else
         allocate (SENDBUF_REAL(1))
         allocate (gl2loc(1))
      end if
      ! make total mapping
      call MPI_GATHERV(local_gl2loc, RECVCOUNT, MPI_INTEGER, &
                       GL2LOC, SENDCOUNTS, DISPLS, MPI_INTEGER, &
                       0, COMM, ierr)
      deallocate (local_gl2loc)
   end subroutine BUILD_GL2LOC
   !----------------------------------------------------------------------

   !----------------------------------------------------------------------
   subroutine BcastToLocal_Int(Val)
      !----------------------------------------------------------------------
      implicit none
      integer, intent(INOUT) :: Val
      integer :: ierr
      call MPI_BCast(Val, 1, MPI_Integer, 0, COMM, ierr)
      !----------------------------------------------------------------------
   end subroutine BcastToLocal_Int
   !----------------------------------------------------------------------

   !----------------------------------------------------------------------
   subroutine BcastToLocal_2DRealArray(Val, NX, NY)
      use mpi_f08, only: MPI_REAL
      implicit none
      real(4), intent(INOUT) :: Val(:, :)
      integer, intent(IN) :: NX, NY
      real(4), allocatable :: TMP(:)
      integer :: ii, jj, NumVals, kk, ierr

      NumVals = NX*NY
      allocate (TMP(NumVals))
      if (MYPROC == 0) then
         kk = 1
         do ii = 1, NX
            do jj = 1, NY
               TMP(kk) = Val(ii, jj)
               kk = kk + 1
            end do
         end do
      end if
      call MPI_BCAST(TMP, NumVals, MPI_REAL, 0, COMM, ierr)
      kk = NumVals
      do ii = NX, 1, -1
         do jj = NY, 1, -1
            Val(ii, jj) = TMP(kk)
            kk = kk - 1
         end do
      end do

      deallocate (TMP)
   end subroutine BcastToLocal_2DRealArray
   !----------------------------------------------------------------------
end module GL2LOC_MAPPING
!----------------------------------------------------------------------
