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
C----------------------------------------------------------------------
      MODULE GL2LOC_MAPPING 
C----------------------------------------------------------------------
C     This module exists to create a method to more efficiently read in
C     input data on core 0 and distribute that data to the appropriate
C     subdomain. It was initially written for use in reading in
C     netCDF-format global baroclinic forcing data to be used in loose,
C     one-way coupling to GOFS. My hope is that it will be useful for
C     others.
C
C     Written by: Coleman Blakely 11/2022
C----------------------------------------------------------------------
      USE SIZES, ONLY : MNPROC, MYPROC
      USE MESH, ONLY : NP ! local number of nodes
      USE GLOBAL, ONLY : nodes_lg, np_g, COMM, realtype
!#ifdef HAVE_MPI_MOD
#ifdef CMPI
      USE MPI
#endif      
      IMPLICIT NONE
      ! Start of variables used
      INTEGER,DIMENSION(:),ALLOCATABLE,PRIVATE :: GL2LOC
      ! sending information
      REAL(8),DIMENSION(:),ALLOCATABLE,PRIVATE :: SENDBUF_REAL
      ! counts of data to send and displacements
      INTEGER,DIMENSION(:),ALLOCATABLE,PRIVATE :: SENDCOUNTS, DISPLS
      ! recieving info
      REAL(8),DIMENSION(:),ALLOCATABLE,PRIVATE :: RECVBUF_REAL
      INTEGER,PRIVATE :: RECVCOUNT

C----------------------------------------------------------------------
      CONTAINS
C----------------------------------------------------------------------
      SUBROUTINE MAPTOLOCAL_REAL(GLOBALDATA,LOCALDATA)
C----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: GLOBALDATA(:)
      REAL(8),INTENT(OUT) :: LOCALDATA(:)
      INTEGER :: kk, ierr, jj, TAG
#ifdef CMPI      
      INTEGER :: RECV_REQ
      INTEGER :: RECV_STAT(MPI_STATUS_SIZE), TEST_STAT(MPI_STATUS_SIZE)
#else
      INTEGER :: loc_mess_stat
#endif     
      TAG = 24
      ! build the global mapping (GL2LOC) if this is the first call
      CALL BUILD_GL2LOC()
      ! To avoid allocating/deallocating a bunch just use a dmy array
      ! that has size equal to the maximum number of nodes on a
      ! processor.
#ifdef CMPI
      IF (MYPROC.EQ.0) THEN
         ! Put global data into sendbufr in the correct order
         DO kk = 1,sum(SENDCOUNTS)
            SENDBUF_REAL(kk) = GLOBALDATA(GL2LOC(kk))
         ENDDO
      ENDIF
      CALL MPI_SCATTERV(SENDBUF_REAL,SENDCOUNTS,DISPLS,REALTYPE,
     &                  RECVBUF_REAL,RECVCOUNT,REALTYPE,
     &                  0,COMM,ierr)
      LOCALDATA = RECVBUF_REAL
#else
      LOCALDATA = GLOBALDATA      
#endif
C----------------------------------------------------------------------
      END SUBROUTINE MAPTOLOCAL_REAL
C----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE BUILD_GL2LOC()
C----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL,SAVE :: first_call = .TRUE.
      INTEGER :: localNumNodes
      INTEGER,ALLOCATABLE :: local_gl2loc(:)
      INTEGER :: kk, ii, ierr, TAG, dmycounter
#ifdef CMPI      
      INTEGER :: stat(MPI_STATUS_SIZE)
#else
      INTEGER :: stat
#endif      
      TAG = 100
#ifdef CMPI
      ! Return if we have already built the table
      IF (first_call.EQV..FALSE.) RETURN
      first_call = .FALSE.
      ! On processor 0, allocate sendcounts and displs
      IF (MYPROC.EQ.0) THEN
         ALLOCATE( SENDCOUNTS(MNPROC) )
         ALLOCATE( DISPLS(MNPROC) )
      ELSE
         ALLOCATE( sendcounts(1) )
         ALLOCATE( DISPLS(1) )
      ENDIF
      ! get recvcount on all processors 
      RECVCOUNT = NP
      ! allocate recvbuf_real
      ALLOCATE( RECVBUF_REAL(RECVCOUNT) )
      ! build local mapping tables
      ALLOCATE( local_gl2loc(RECVCOUNT) )
      DO kk = 1,RECVCOUNT
         local_gl2loc(kk) = ABS(nodes_lg(kk))
      ENDDO
      ! gather all recvcount values on proc 0
      CALL MPI_GATHER(RECVCOUNT, 1, MPI_INTEGER,
     &                SENDCOUNTS, 1, MPI_INTEGER, 
     &                0, COMM, ierr)
      IF (MYPROC.eq.0) THEN
         DISPLS(1) = 0
         DO kk = 2,MNPROC
            DISPLS(kk) = SUM(SENDCOUNTS(1:kk-1))
         ENDDO
         ! allocate gl2loc
         ALLOCATE( GL2LOC(SUM(SENDCOUNTS)) )
         ALLOCATE( SENDBUF_REAL(SUM(SENDCOUNTS)) )
      ELSE
         ALLOCATE( SENDBUF_REAL(1) )
         ALLOCATE( gl2loc(1) )
      ENDIF
      ! make total mapping
      CALL MPI_GATHERV(local_gl2loc, RECVCOUNT, MPI_INTEGER,
     &                 GL2LOC, SENDCOUNTS, DISPLS, MPI_INTEGER,
     &                 0, COMM, ierr)
      DEALLOCATE( local_gl2loc )
      RETURN
#endif
C----------------------------------------------------------------------
      END SUBROUTINE BUILD_GL2LOC
C----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE BcastToLocal_Int(Val)
C----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER,INTENT(INOUT) :: Val
      INTEGER :: ierr

#ifdef CMPI
      CALL MPI_BCast(Val,1,MPI_Integer,0,COMM,ierr)
#else
      Val = Val
#endif
C----------------------------------------------------------------------
      END SUBROUTINE BcastToLocal_Int
C----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE BcastToLocal_2DRealArray(Val,NX, NY)
C----------------------------------------------------------------------
      IMPLICIT NONE
      REAL,INTENT(INOUT) :: Val(:,:)
      INTEGER,INTENT(IN) :: NX, NY
      REAL(8),ALLOCATABLE :: TMP(:)
      INTEGER :: ii, jj, NumVals, kk, ierr
#ifdef CMPI
      NumVals = NX*NY
      ALLOCATE( TMP(NumVals) )
      IF (MYPROC.EQ.0) THEN
         kk = 1
         DO ii = 1,NX
            DO jj = 1,NY
               TMP(kk) = Val(ii,jj)
               kk = kk + 1
            ENDDO
         ENDDO
      ENDIF
      CALL MPI_BCAST(TMP,NumVals,REALTYPE,0,COMM,ierr)
      kk = NumVals
      DO ii = NX,1,-1
         DO jj = NY,1,-1
            Val(ii,jj) = TMP(kk)
            kk = kk - 1
         ENDDO
      ENDDO
      
      DEALLOCATE( TMP )
#endif
C----------------------------------------------------------------------
      END SUBROUTINE BcastToLocal_2DRealArray
C----------------------------------------------------------------------
      ENDMODULE GL2LOC_MAPPING 
C----------------------------------------------------------------------
