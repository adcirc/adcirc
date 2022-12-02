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
      USE SIZES, ONLY : MNPROC, MYPROC, SZ
      USE MESH, ONLY : NP ! local number of nodes
      USE GLOBAL, ONLY : nodes_lg, np_g, COMM, realtype
#IFDEF HAVE_MPI_MOD
      USE MPI
#ENDIF      
      IMPLICIT NONE
      ! Define derived datatype for mapping. The format is fairly
      ! simple. GLIDX will be of size MNPROC and for each subdomain will
      ! contain the number of local nodes and the global node number of
      ! each local node. This will be stored on core 0
      TYPE :: GLIDX
         INTEGER :: NumNodes
         INTEGER,ALLOCATABLE :: GL(:)
         REAL(SZ),ALLOCATABLE :: realbuffer(:)
      END TYPE GLIDX
      ! Start of variables used
      TYPE(GLIDX),DIMENSION(:),ALLOCATABLE,PRIVATE :: GL2LOC
      ! array of send requests
      INTEGER,DIMENSION(:),ALLOCATABLE,PRIVATE :: SEND_REQ
      ! array of send statuses
      INTEGER,DIMENSION(:,:),ALLOCATABLE,PRIVATE :: SEND_STAT
      INTEGER,PRIVATE :: NPmax ! max number of points in a subdomain
      ! to prevent allocation/deallocation slowing us down define
      ! variables that persist
      REAL(SZ),ALLOCATABLE,PRIVATE :: tmp_real(:)

C----------------------------------------------------------------------
      CONTAINS
C----------------------------------------------------------------------
      SUBROUTINE MAPTOLOCAL_REAL(GLOBALDATA,LOCALDATA)
C----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(SZ),INTENT(IN) :: GLOBALDATA(:)
      REAL(SZ),INTENT(OUT) :: LOCALDATA(:)
      INTEGER :: kk, ierr, jj, TAG
#IFDEF CMPI      
      INTEGER :: RECV_REQ
      INTEGER :: RECV_STAT(MPI_STATUS_SIZE), TEST_STAT(MPI_STATUS_SIZE)
#ELSE
      INTEGER :: loc_mess_stat
#ENDIF     
      TAG = 24
      ! build the global mapping (GL2LOC) if this is the first call
!      IF (MYPROC.EQ.0) THEN
!         WRITE(16,*) 'About to build gl2loc'
!      ENDIF
      CALL BUILD_GL2LOC()
!      WRITE(16,*) 'Built Global-to-Local mapping'
      ! To avoid allocating/deallocating a bunch just use a dmy array
      ! that has size equal to the maximum number of nodes on a
      ! processor.
#IFDEF CMPI
      ! loop through each processor and build the local arrays
      IF (MYPROC.EQ.0) THEN
         ! handle processor 0 first
         DO jj = 1,GL2LOC(1)%NumNodes
            LOCALDATA(jj) = GLOBALDATA(GL2LOC(1)%GL(jj))
         ENDDO
         ! handle rest of processors
         DO kk = 2,MNPROC
            tmp_real = 0_SZ
            DO jj = 1,GL2LOC(kk)%NumNodes
               GL2LOC(kk)%realbuffer(jj) = GLOBALDATA(GL2LOC(kk)%GL(jj))
            ENDDO
            WRITE(16,*) 'About to send info to proc ',kk
            CALL MPI_ISEND(GL2LOC(kk)%realbuffer,
     &                     GL2LOC(kk)%NumNodes, REALTYPE, kk-1,
     &                     TAG, COMM, SEND_REQ(kk-1), ierr) 
            WRITE(16,*) 'Just sent info to proc ',kk
         ENDDO
         CALL MPI_WAITALL(MNPROC-1, SEND_REQ, SEND_STAT, ierr)
         WRITE(16,*) 'Done waiting on sends'
      ELSE
         WRITE(16,*) 'About to receive data from proc 0'
         CALL MPI_IRECV(LOCALDATA, NP, REALTYPE, 0, MPI_ANY_TAG, COMM,
     &                  RECV_REQ, ierr)
         WRITE(16,*) 'Just received data from proc 0'
         CALL MPI_WAIT(RECV_REQ, RECV_STAT, ierr)
      ENDIF
#ELSE
      LOCALDATA = GLOBALDATA      
#ENDIF
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
      INTEGER :: kk, ierr, TAG
#IFDEF CMPI      
      INTEGER :: stat(MPI_STATUS_SIZE)
#ELSE
      INTEGER :: stat
#ENDIF      
      TAG = 100
#IFDEF CMPI
      ! Return if we have already build the table
      IF (first_call.EQV..FALSE.) RETURN
      first_call = .FALSE.
      ! make dummy variables to put in gl2loc on each processor
      localNumNodes = NP
      NPmax = localNumNodes
      ALLOCATE( local_gl2loc(localNumNodes) )
      DO kk = 1,localNumNodes
         local_gl2loc(kk) = abs(nodes_lg(kk))
      END DO
      WRITE(16,*) 'PROC ',MYPROC,' built local_gl2loc'
      ! First, put the data from proc 0 in the proper place in GL2LOC
      IF (MYPROC.EQ.0) THEN
         ! allocate send requests for irecv/isend. Note that we only
         ! need MNPROC-1 request tags since we don't send to proc 0
         ALLOCATE( SEND_REQ(MNPROC-1) )
         ALLOCATE( SEND_STAT(MPI_STATUS_SIZE,MNPROC-1) )
         ALLOCATE( GL2LOC(MNPROC) )
         GL2LOC(1)%NumNodes = localNumNodes
         ALLOCATE( GL2LOC(1)%GL(GL2LOC(1)%NumNodes) )
         GL2LOC(1)%GL = local_gl2loc
         DEALLOCATE( local_gl2loc )
         ! Gather all the mapping info onto proc 0
         DO kk = 1,MNPROC-1
            CALL MPI_RECV(localNumNodes, 1, MPI_INTEGER, kk, 
     &                    MPI_ANY_TAG, COMM, stat, ierr)
            ALLOCATE( local_gl2loc(localNumNodes) )
            CALL MPI_RECV(local_gl2loc, localNumNodes, MPI_INTEGER, kk,
     &                    MPI_ANY_TAG, COMM, stat, ierr)
            GL2LOC(kk+1)%NumNodes = localNumNodes
            ALLOCATE( GL2LOC(kk+1)%GL(GL2LOC(kk+1)%NumNodes) )
            ALLOCATE( GL2LOC(kk+1)%realbuffer(GL2LOC(kk+1)%NumNodes) )
            GL2LOC(kk+1)%GL = local_gl2loc
            DEALLOCATE( local_gl2loc )
            ! check if NP of this is greater than NPmax
            NPmax = max(NPmax,localNumNodes)
         END DO
         WRITE(16,*) 'Received all mapping info'
      ELSE
         CALL MPI_SEND(localNumNodes, 1, MPI_INTEGER, 0, TAG, COMM,ierr)
         CALL MPI_SEND(local_gl2loc, localNumNodes, MPI_INTEGER, 0, TAG,
     &                 COMM,ierr)
         DEALLOCATE( local_gl2loc )
         WRITE(16,*) 'Sent mapping info to proc 1'
      ENDIF
#ENDIF
C----------------------------------------------------------------------
      END SUBROUTINE BUILD_GL2LOC
C----------------------------------------------------------------------
      ENDMODULE GL2LOC_MAPPING 
C----------------------------------------------------------------------
