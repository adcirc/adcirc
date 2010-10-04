C******************************************************************************
C PADCIRC VERSION 47.xx 10/13/2006                                             *
C******************************************************************************
C 
      MODULE MESSENGER
      USE SIZES
      USE GLOBAL, ONLY: C3D, COMM, NODES_LG, IMAP_EL_LG,
     &    IMAP_STAE_LG, IMAP_STAV_LG, IMAP_STAM_LG, IMAP_STAC_LG,
     &    COMM_WRITER, COMM_WRITEH, COMM_HSLEEP, WRITER_ID, SIG_VAL, SIG_WRITE, SIG_TERM,
     &    FLOAT_TYPE, REALTYPE, DBLETYPE
#ifdef HAVE_MPI_MOD
      USE MPI  
#endif
      IMPLICIT NONE
      SAVE

 
C--------------------------------------------------------------------------
C  This module supplies the MPI Message-Passing Interface for PADCIRC.
C  Uses asynchronous communication with buffer packing as performance 
C  enhancement for "cluster" architectures.
C--------------------------------------------------------------------------

 
C  Message-Passing Array space
 
      INTEGER ::  MPI_COMM_ADCIRC                     ! Local communicator
      INTEGER ::  NEIGHPROC, RDIM, IERR
      INTEGER ::  TAG = 100
      INTEGER ::  COMM_COMP                           ! COMMUNICATOR FOR COMPUTATION
      INTEGER ::  GROUP_WORLD, GROUP_COMP
      INTEGER, ALLOCATABLE ::  GROUP_WRITER(:)        ! GROUPS FOR GLOBAL FILE WRITING
      INTEGER, ALLOCATABLE ::  GROUP_WRITEH(:)        ! GROUPS FOR HOTSTART FILE WRITING !st3 100711
      INTEGER, ALLOCATABLE ::  GROUP_HSLEEP(:)        ! GROUPS FOR HOTSTART FILE WRITING !st3 100711
      LOGICAL, ALLOCATABLE :: RESNODE(:)

      INTEGER, ALLOCATABLE :: IPROC(:), NNODELOC(:), NNODSEND(:), 
     &    NNODRECV(:), IBELONGTO(:),ISENDLOC(:,:), IRECVLOC(:,:), 
     &    ISENDBUF(:,:), IRECVBUF(:,:)
 
      INTEGER, ALLOCATABLE :: REQ_I1(:), REQ_I2(:)
      INTEGER, ALLOCATABLE :: STAT_I1(:,:), STAT_I2(:,:)
      INTEGER, ALLOCATABLE :: REQ_R1(:), REQ_R2(:), REQ_R3(:)
      INTEGER, ALLOCATABLE :: STAT_R1(:,:), STAT_R2(:,:), STAT_R3(:,:)
      INTEGER, ALLOCATABLE :: REQ_R3D(:), STAT_R3D(:,:)
      INTEGER, ALLOCATABLE :: REQ_C3D(:), STAT_C3D(:,:)
      INTEGER, ALLOCATABLE :: INDX(:)
      REAL(SZ),ALLOCATABLE :: SENDBUF(:,:), RECVBUF(:,:)
C
C
      CONTAINS

C---------------------end of data declarations--------------------------------C


      SUBROUTINE MSG_INIT (MPI_COMM)
C--------------------------------------------------------------------------
C  Routine performs following steps:
C   (1)  define mpi data types to be used
C   (2)  initialize MPI, 
C   (3)  get number of processors,
C   (4)  get MPI rank of processor 
C   (5)  initialize adcirc COMM communicator MPI_COMM_WORLD
C  vjp  10/3/2006
C--------------------------------------------------------------------------
      IMPLICIT NONE

#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif

      INTEGER, OPTIONAL :: MPI_COMM
      INTEGER :: I
      INTEGER, ALLOCATABLE :: RANKS(:)
      INTEGER :: SENDBUF, RECVBUF
      INTEGER :: IRANK_SLEEP(2)        !st3 100711  for hsfile


#ifdef  REAL4
      REALTYPE = MPI_REAL
      DBLETYPE = MPI_DOUBLE_PRECISION
      float_type = REALTYPE                    ! used in globalio
#else
      REALTYPE = MPI_DOUBLE_PRECISION
      DBLETYPE = MPI_DOUBLE_PRECISION
      float_type = REALTYPE                    ! used in globalio
#endif

      IF (PRESENT(MPI_COMM)) THEN
C.......Duplicate communicator passed from outside
        CALL MPI_COMM_DUP (MPI_COMM,MPI_COMM_ADCIRC,IERR)
      ELSE
C.......Initialize MPI
        CALL MPI_INIT (IERR)
C.......Duplicate communicator
        CALL MPI_COMM_DUP (MPI_COMM_WORLD,MPI_COMM_ADCIRC,IERR)
      ENDIF

      CALL MPI_COMM_SIZE (MPI_COMM_ADCIRC,MNALLPROC,IERR)   ! Get number of procs
      CALL MPI_COMM_RANK (MPI_COMM_ADCIRC,MYPROC,IERR)   ! Get MPI rank
      
      MNPROC = MNALLPROC - MNWPROC             ! MNALLPROC = MNPROC + MNWPROC
     &                   - MNWPROH             !             + MNWPROH  !st3 100711 for hsfile


      ALLOCATE(RANKS(MNPROC+1))

C...  Create a communicator for computation
      DO I=1,MNPROC
        RANKS(I) = I-1
      ENDDO

      CALL MPI_COMM_GROUP(MPI_COMM_ADCIRC,GROUP_WORLD,IERR)
      CALL MPI_GROUP_INCL(GROUP_WORLD,MNPROC,RANKS,GROUP_COMP,IERR)
      CALL MPI_COMM_CREATE(MPI_COMM_ADCIRC,GROUP_COMP,COMM_COMP,IERR)

      WRITER_ID = 0

      IF(MNWPROC > 0) THEN
C...  Allocate memory for groups and communicators
        ALLOCATE(GROUP_WRITER(MNWPROC),COMM_WRITER(MNWPROC))

C...  Create communicators for global file writings
        DO I=1,MNWPROC
          RANKS(MNPROC+1) = MNPROC - 1 + I

          CALL MPI_GROUP_INCL(GROUP_WORLD,MNPROC+1,RANKS,
     &                        GROUP_WRITER(I),IERR)
          CALL MPI_COMM_CREATE(MPI_COMM_ADCIRC,GROUP_WRITER(I),
     &                         COMM_WRITER(I),IERR)

          IF(MYPROC == MNPROC - 1 + I) THEN
            WRITER_ID = I
          ENDIF
        ENDDO

      ENDIF
!st3 for hsfile 05.14.2010
#ifdef NETCDF
      MNWPROH = 0
#endif
#ifndef CMPI
      MNWPROH = 0
#endif
      IF(MNWPROH > 0) THEN
        IRANK_SLEEP(1) = 0
        ALLOCATE(GROUP_WRITEH(MNWPROH),COMM_WRITEH(MNWPROH))
        ALLOCATE(GROUP_HSLEEP(MNWPROH),COMM_HSLEEP(MNWPROH))
        DO I=1,MNWPROH
          RANKS(MNPROC+1) = MNPROC + MNWPROC - 1 + I

          CALL MPI_GROUP_INCL(GROUP_WORLD,MNPROC+1,RANKS,
     &                        GROUP_WRITEH(I),IERR)
          CALL MPI_COMM_CREATE(MPI_COMM_ADCIRC,GROUP_WRITEH(I),
     &                         COMM_WRITEH(I),IERR)

          IRANK_SLEEP(2) = MNPROC + MNWPROC - 1 + I
          CALL MPI_GROUP_INCL(GROUP_WORLD,2,IRANK_SLEEP,
     &                        GROUP_HSLEEP(I),IERR)
          CALL MPI_COMM_CREATE(MPI_COMM_ADCIRC,GROUP_HSLEEP(I),
     &                         COMM_HSLEEP(I),IERR)

          IF(MYPROC == MNPROC + MNWPROC - 1 + I) THEN
            WRITER_ID = I
          ENDIF
        ENDDO
      ENDIF


      DEALLOCATE(RANKS)

      COMM = COMM_COMP

      RETURN
      END SUBROUTINE MSG_INIT


      SUBROUTINE MSG_FINI (NO_MPI_FINALIZE)
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      LOGICAL, OPTIONAL :: NO_MPI_FINALIZE
      INTEGER I
C--------------------------------------------------------------------------
C  Delete MPI resources and Shutdown MPI library.
C  vjp  8/29/1999
C--------------------------------------------------------------------------
      IF(MNWPROC > 0) THEN
        IF(MYPROC.eq.0) THEN
          DO I=1,MNWPROC
            WRITE(16,*)'PROC ',MYPROC,' IS SENDING SIG_TERM TO WRITER ',
     &           I
            CALL MPI_SEND(SIG_TERM,1,MPI_INTEGER,MNPROC,
     &           TAG,COMM_WRITER(I),IERR)
          ENDDO
        ENDIF
      ENDIF

      IF(MNWPROH > 0) THEN             !st3 100711 for hsfile
        IF(MYPROC.eq.0) THEN
          DO I=1,MNWPROH
            WRITE(16,*)'PROC ',MYPROC,' IS SENDING SIG_TERM TO HSWRITER',
     &           I
            CALL MPI_BARRIER(COMM_HSLEEP(I),IERR)
            CALL MPI_SEND(SIG_TERM,1,MPI_INTEGER,MNPROC,
     &           TAG,COMM_WRITEH(I),IERR)
          ENDDO
        ENDIF
      ENDIF

      IF (PRESENT(NO_MPI_FINALIZE)) THEN
         IF (.NOT.NO_MPI_FINALIZE) THEN
            CALL MPI_FINALIZE(IERR)
            IF (MYPROC.EQ.0)  
     &      PRINT *, "MPI terminated with Status = ",IERR      
         ENDIF
      ELSE
         CALL MPI_FINALIZE(IERR)
         IF (MYPROC.EQ.0)  
     &   PRINT *, "MPI terminated with Status = ",IERR      
      ENDIF

      RETURN
      END SUBROUTINE MSG_FINI



      SUBROUTINE MSG_TABLE () 
      USE GLOBAL
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      INTEGER :: IDPROC, NLOCAL, I, J, jdumy_loc
      INTEGER :: jdumy,jdumy_G,jdumy_max, inputFileFmtVn
      CHARACTER(10) :: BlkName
      LOGICAL FileFound
C
C--------------------------------------------------------------------------
C  Routine preforms following steps:
C
C   (1) Read Message-Passing Information from file "fort.18"
C   (2) Determine resident nodes: RESNODE(I) is true  if I is resident node
C   (3) Determine ghost nodes:    RESNODE(I) is false if I is ghost node    
C   (4) Determine number of neighbor subdomains
C   (5) MPI rank of each neighbor and number of ghosts nodes to receive
C   (6) Read Message-Passing Receive List
C   (7) MPI rank of each neighbor and number of ghosts nodes to send
C   (8) Read Message-Passing Send List
C  vjp  10/13/2006
C--------------------------------------------------------------------------
 
      FileFound = .FALSE.
      INQUIRE(FILE=TRIM(INPUTDIR)//'/'//'fort.18',EXIST=FileFound)
      IF (FileFound.eqv..FALSE.) THEN
         WRITE(*,*) "ERROR: ",myProc,": The file ",
     &      TRIM(INPUTDIR)//'/'//'fort.18 was not found.'
         WRITE(*,9973)      ! execution terminated
         CALL MSG_FINI()
         STOP
      ENDIF
C
      OPEN(18,FILE=TRIM(INPUTDIR)//'/'//'fort.18',STATUS='OLD')
C
      READ(18,3020) BlkName, inputFileFmtVn
      if (Trim(BlkName) /= 'FileFmt' .and.
     $  .not. CMP_VERSION_NUMBERS(FileFmtVersion, inputFileFmtVn)) then
        write(16,*) 'File Format of Fort.18 does not match aborting'
        stop
      end if

C Read Global number of elements and Local-to_Global element map ( used by module global_io )
Casey 100209: Changed from FMT=3015.
      READ(18,'(8X,3I12)') NE_G
      ALLOCATE ( IMAP_EL_LG(MNE) )
      DO I=1,MNE
         READ(18,*) IMAP_EL_LG(I)
      ENDDO

C Read Global number of nodes and Local-to_Global node map ( used by module global_io )
      READ(18,3015) NP_G
      ALLOCATE ( NODES_LG(MNP) )
      DO I=1,MNP
         READ(18,*) NODES_LG(I)
      ENDDO

C  This information is provided for relocalizing fort.15
C  Just read past it.
      READ(18,'(8X,I8)') jdumy                          ! nfluxf for subdomain
      READ(18,'(8X,3I8)') jdumy_g, jdumy_max, jdumy_loc ! neta for subdomain
      DO I=1, jdumy_loc
         READ(18,'(I8)') jdumy                          ! obnode_lg table
      ENDDO

C Read Global indexes of Elevation Station nodes ( used by module global_io )
      READ(18,3015) NSTAE_G
      IF (NSTAE > 0) THEN
        ALLOCATE ( IMAP_STAE_LG(NSTAE) )
        DO I=1,NSTAE
           READ(18,'(I8)') IMAP_STAE_LG(I)              
        ENDDO
      ENDIF

C Read Global indexes of Velocity  Station nodes ( used by module global_io )
      READ(18,3015) NSTAV_G
      IF (NSTAV > 0) THEN
        ALLOCATE ( IMAP_STAV_LG(NSTAV) )
        DO I=1,NSTAV
           READ(18,'(I8)') IMAP_STAV_LG(I)              
        ENDDO
      ENDIF

C Read Global indexes of Meteorlogical Station nodes ( used by module global_io )
      READ(18,3015) NSTAM_G
      IF (NSTAM > 0) THEN
        ALLOCATE ( IMAP_STAM_LG(NSTAM) )
        DO I=1,NSTAM
           READ(18,'(I8)') IMAP_STAM_LG(I)              
        ENDDO
      ENDIF

C Read Global indexes of Concentration Station nodes ( used by module global_io )
      READ(18,3015) NSTAC_G
      IF (NSTAC > 0) THEN
        ALLOCATE ( IMAP_STAC_LG(NSTAC) )
        DO I=1,NSTAC
           READ(18,'(I8)') IMAP_STAC_LG(I)              
        ENDDO
      ENDIF

C---------------------------------------------------------------------------------
C--Message-Passing tables start here
C---------------------------------------------------------------------------------

      READ(18,3010) IDPROC,NLOCAL    
      ALLOCATE ( NNODELOC(NLOCAL) )
      READ(18,1130) (NNODELOC(I), I=1,NLOCAL)
      ALLOCATE ( IBELONGTO(MNP),RESNODE(MNP) )
 
      DO I=1,MNP
         IBELONGTO(I) = 0
      ENDDO
      DO I=1,NLOCAL
         IBELONGTO(NNODELOC(I)) = IDPROC + 1
      ENDDO
      DO I=1, MNP
         IF (IBELONGTO(I)-1.EQ.MYPROC) THEN
           RESNODE(I) = .TRUE.
         ELSE 
           RESNODE(I) = .FALSE.
         ENDIF
      ENDDO
 
      READ(18,3015) NEIGHPROC

      RDIM = 2*NEIGHPROC
      ALLOCATE( INDX(RDIM) )
      ALLOCATE( IPROC(NEIGHPROC),NNODRECV(NEIGHPROC) )
      ALLOCATE( IRECVLOC(MNP,NEIGHPROC) )
 
      DO J=1,NEIGHPROC
         READ(18,3010) IPROC(J),NNODRECV(J)
         READ(18,1130) (IRECVLOC(I,J), I=1,NNODRECV(J))
      ENDDO
 
      ALLOCATE( NNODSEND(NEIGHPROC) )
      ALLOCATE( ISENDLOC(MNP,NEIGHPROC) )
 
      DO J=1,NEIGHPROC
         READ(18,3010) IPROC(J),NNODSEND(J)
         READ(18,1130) (ISENDLOC(I,J), I=1,NNODSEND(J))
      ENDDO
 
      CLOSE(18)
      RETURN
1130  FORMAT(8X,9I8)
3010  FORMAT(8X,2I8)
3015  FORMAT(8X,3I8)
3020  format(a8,I8)
9973  FORMAT(/,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
      END SUBROUTINE MSG_TABLE
 


      SUBROUTINE MSG_START ()
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
C--------------------------------------------------------------------------
C  Routine allocates message-passing space
C  vjp  10/01/2006
C--------------------------------------------------------------------------
      INTEGER :: J
      ALLOCATE ( ISENDBUF(MNP,NEIGHPROC), IRECVBUF(MNP,NEIGHPROC) )
      ALLOCATE ( REQ_I1(RDIM),REQ_I2(RDIM) )
      ALLOCATE ( REQ_R1(RDIM),REQ_R2(RDIM),REQ_R3(RDIM) )
      ALLOCATE ( STAT_I1(MPI_STATUS_SIZE,RDIM),       
     &           STAT_I2(MPI_STATUS_SIZE,RDIM) )
      ALLOCATE ( STAT_R1(MPI_STATUS_SIZE,RDIM),       
     &           STAT_R2(MPI_STATUS_SIZE,RDIM),
     &           STAT_R3(MPI_STATUS_SIZE,RDIM) )
 
      IF (C3D) THEN
         ALLOCATE ( SENDBUF(2*MNP*MNFEN,NEIGHPROC) )
         ALLOCATE ( RECVBUF(2*MNP*MNFEN,NEIGHPROC) )
         ALLOCATE ( REQ_R3D(RDIM) )
         ALLOCATE ( STAT_R3D(MPI_STATUS_SIZE,RDIM) )
         ALLOCATE ( REQ_C3D(RDIM) )
         ALLOCATE ( STAT_C3D(MPI_STATUS_SIZE,RDIM) )
      ELSE
         ALLOCATE ( SENDBUF(MNP,NEIGHPROC) )
         ALLOCATE ( RECVBUF(MNP,NEIGHPROC) )
      ENDIF
      RETURN
      END SUBROUTINE MSG_START


      SUBROUTINE UPDATEI( IVEC1, IVEC2, NMSG )
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      INTEGER,   INTENT(IN) :: NMSG
      INTEGER,   INTENT(INOUT) :: IVEC1(*),IVEC2(*)
      INTEGER :: N,I,J,NCOUNT,NFINI,TOT
C--------------------------------------------------------------------------
C  Update 1 or 2 Integer Arrays's Ghost Cells using asynchronous
C  and persistent message-passing.
C
C  vjp  8/06/1999
C--------------------------------------------------------------------------
                             !..Pack 1 or 2 Messages
      DO J=1,NEIGHPROC
         NCOUNT = 0
         DO I=1,NNODSEND(J)
            NCOUNT = NCOUNT+1
            ISENDBUF(NCOUNT,J)=IVEC1(ISENDLOC(I,J))
         ENDDO
         IF (NMSG.GT.1) THEN
           DO I=1,NNODSEND(J)
              NCOUNT = NCOUNT+1
              ISENDBUF(NCOUNT,J)=IVEC2(ISENDLOC(I,J))
           ENDDO
         ENDIF
      ENDDO
                          ! Send/receive messages to/from all neighbors
      IF (NMSG.EQ.1) THEN
        DO J=1,NEIGHPROC   
           CALL MPI_IRECV( IRECVBUF(1,J), NNODRECV(J), 
     &             MPI_INTEGER,IPROC(J), TAG, COMM, REQ_I1(J),IERR)
           CALL MPI_ISEND( ISENDBUF(1,J), NNODSEND(J), 
     &       MPI_INTEGER,IPROC(J), TAG, COMM, REQ_I1(J+NEIGHPROC),IERR )
        ENDDO
      ELSE
        DO J=1,NEIGHPROC   
           CALL MPI_IRECV( IRECVBUF(1,J), 2*NNODRECV(J), 
     &             MPI_INTEGER,IPROC(J), TAG, COMM, REQ_I2(J),IERR)
           CALL MPI_ISEND( ISENDBUF(1,J), 2*NNODSEND(J), 
     &       MPI_INTEGER,IPROC(J), TAG, COMM, REQ_I2(J+NEIGHPROC),IERR )
        ENDDO
      ENDIF
                          !..Unpack Received messages as they arrive  
      IF (NMSG.EQ.1) THEN   
        TOT = 0
        DO WHILE (TOT.LT.RDIM)
           DO N=1, RDIM
              INDX(N) = 0
           ENDDO
           CALL MPI_WAITSOME( RDIM,REQ_I1,NFINI,INDX,STAT_I1,IERR )
           TOT = TOT + NFINI
           DO N=1, NFINI
              IF (INDX(N).GT.0.AND.INDX(N).LE.RDIM)  THEN
                IF (INDX(N).LE.NEIGHPROC) THEN
                  J = INDX(N)
                  NCOUNT = 0
                  DO I=1,NNODRECV(J)
                     NCOUNT = NCOUNT+1
                     IVEC1(IRECVLOC(I,J)) = IRECVBUF(NCOUNT,J)
                  ENDDO
                ENDIF
              ENDIF
           ENDDO
        ENDDO
      ELSE
        TOT = 0
        DO WHILE (TOT.LT.RDIM)
           DO N=1, RDIM
              INDX(N) = 0
           ENDDO
           CALL MPI_WAITSOME( RDIM,REQ_I2,NFINI,INDX,STAT_I2,IERR )
           TOT = TOT + NFINI
           DO N=1, NFINI
              IF (INDX(N).GT.0.AND.INDX(N).LE.RDIM)  THEN
                IF (INDX(N).LE.NEIGHPROC) THEN
                  J = INDX(N)
                  NCOUNT = 0
                  DO I=1,NNODRECV(J)
                     NCOUNT = NCOUNT+1
                     IVEC1(IRECVLOC(I,J)) = IRECVBUF(NCOUNT,J)
                  ENDDO
                  DO I=1,NNODRECV(J)
                     NCOUNT = NCOUNT+1
                     IVEC2(IRECVLOC(I,J)) = IRECVBUF(NCOUNT,J)
                  ENDDO
                ENDIF
              ENDIF
           ENDDO
        ENDDO
      ENDIF
 999  RETURN
      END SUBROUTINE UPDATEI


      SUBROUTINE UPDATER( VEC1, VEC2, VEC3, NMSG )
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      INTEGER,  INTENT(IN) ::  NMSG
      REAL(SZ), INTENT(INOUT) ::  VEC1(*),VEC2(*),VEC3(*)
      INTEGER :: N,I,J,NCOUNT,NFINI,TOT
C--------------------------------------------------------------------------
C  Update 1, 2, or 3 Integer Arrays's Ghost Cells using asynchronous
C  and persistent message-passing.
C
C  vjp  8/06/1999
C--------------------------------------------------------------------------
 
      DO J=1,NEIGHPROC
         NCOUNT = 0
         DO I=1,NNODSEND(J)
            NCOUNT = NCOUNT+1
            SENDBUF(NCOUNT,J)=VEC1(ISENDLOC(I,J))
         ENDDO
         IF (NMSG.GT.1) THEN
           DO I=1,NNODSEND(J)
              NCOUNT = NCOUNT+1
              SENDBUF(NCOUNT,J)=VEC2(ISENDLOC(I,J))
           ENDDO
         ENDIF
         IF (NMSG.GT.2) THEN
           DO I=1,NNODSEND(J)
              NCOUNT = NCOUNT+1
              SENDBUF(NCOUNT,J)=VEC3(ISENDLOC(I,J))
           ENDDO
         ENDIF
      ENDDO
 
      IF (NMSG.EQ.1) THEN
        DO J=1,NEIGHPROC  
           CALL MPI_IRECV( RECVBUF(1,J), NNODRECV(J), 
     &             REALTYPE,IPROC(J), TAG, COMM, REQ_R1(J),IERR)
            CALL MPI_ISEND( SENDBUF(1,J), NNODSEND(J), 
     &         REALTYPE,IPROC(J), TAG, COMM, REQ_R1(J+NEIGHPROC),IERR)
        ENDDO
      ELSEIF (NMSG.EQ.2) THEN
        DO J=1,NEIGHPROC  
           CALL MPI_IRECV( RECVBUF(1,J), 2*NNODRECV(J), 
     &        REALTYPE,IPROC(J), TAG, COMM, REQ_R2(J),IERR)
           CALL MPI_ISEND( SENDBUF(1,J), 2*NNODSEND(J), 
     &        REALTYPE,IPROC(J), TAG, COMM, REQ_R2(J+NEIGHPROC),IERR)
        ENDDO
      ELSE
        DO J=1,NEIGHPROC  
           CALL MPI_IRECV( RECVBUF(1,J), 3*NNODRECV(J), 
     &             REALTYPE,IPROC(J), TAG, COMM, REQ_R3(J),IERR)
           CALL MPI_ISEND( SENDBUF(1,J), 3*NNODSEND(J), 
     &        REALTYPE,IPROC(J), TAG, COMM, REQ_R3(J+NEIGHPROC),IERR)
        ENDDO
      ENDIF
 
      IF (NMSG.EQ.1) THEN   
        TOT = 0
        DO WHILE (TOT.LT.RDIM)
           DO N=1, RDIM
              INDX(N) = 0
           ENDDO
           CALL MPI_WAITSOME( RDIM,REQ_R1,NFINI,INDX,STAT_R1,IERR )
           TOT = TOT + NFINI
           DO N=1, NFINI
              IF (INDX(N).GT.0.AND.INDX(N).LE.RDIM)  THEN
                IF (INDX(N).LE.NEIGHPROC) THEN
                  J = INDX(N)
                  NCOUNT = 0
                  DO I=1,NNODRECV(J)
                     NCOUNT = NCOUNT+1
                     VEC1(IRECVLOC(I,J)) = RECVBUF(NCOUNT,J)
                  ENDDO
                ENDIF
              ENDIF
           ENDDO
        ENDDO
        GOTO 999
      ELSEIF (NMSG.EQ.2) THEN
        TOT = 0
        DO WHILE (TOT.LT.RDIM)
           DO N=1, RDIM
              INDX(N) = 0
           ENDDO
           CALL MPI_WAITSOME( RDIM,REQ_R2,NFINI,INDX,STAT_R2,IERR )
           TOT = TOT + NFINI
           DO N=1, NFINI
              IF (INDX(N).GT.0.AND.INDX(N).LE.RDIM)  THEN
                IF (INDX(N).LE.NEIGHPROC) THEN
                  J = INDX(N)
                  NCOUNT = 0
                  DO I=1,NNODRECV(J)
                     NCOUNT = NCOUNT+1
                     VEC1(IRECVLOC(I,J)) = RECVBUF(NCOUNT,J)
                  ENDDO
                  DO I=1,NNODRECV(J)
                     NCOUNT = NCOUNT+1
                     VEC2(IRECVLOC(I,J)) = RECVBUF(NCOUNT,J)
                  ENDDO
                ENDIF
              ENDIF
           ENDDO
        ENDDO
        GOTO 999
      ELSE
        TOT = 0
        DO WHILE (TOT.LT.RDIM)
           DO N=1, RDIM
              INDX(N) = 0
           ENDDO
           CALL MPI_WAITSOME( RDIM,REQ_R3,NFINI,INDX,STAT_R3,IERR )
           TOT = TOT + NFINI
cdebug     print *, myproc, tot,nfini,INDX(1),INDX(2)
           DO N=1, NFINI
              IF (INDX(N).GT.0.AND.INDX(N).LE.RDIM)  THEN
                IF (INDX(N).LE.NEIGHPROC) THEN
                  J = INDX(N)
                  NCOUNT = 0
                  DO I=1,NNODRECV(J)
                     NCOUNT = NCOUNT+1
                     VEC1(IRECVLOC(I,J)) = RECVBUF(NCOUNT,J)
                  ENDDO
                  DO I=1,NNODRECV(J)
                     NCOUNT = NCOUNT+1
                     VEC2(IRECVLOC(I,J)) = RECVBUF(NCOUNT,J)
                  ENDDO
                  DO I=1,NNODRECV(J)
                     NCOUNT = NCOUNT+1
                     VEC3(IRECVLOC(I,J)) = RECVBUF(NCOUNT,J)
                  ENDDO
                ENDIF
              ENDIF
           ENDDO
        ENDDO
        GOTO 999
      ENDIF
  
 999  RETURN
      END SUBROUTINE UPDATER


      SUBROUTINE UPDATER3D( VEC )
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      REAL(SZ), INTENT(INOUT) ::  VEC(MNP,MNFEN)
      INTEGER :: N,I,J,K,NCOUNT,NFINI,TOT
C--------------------------------------------------------------------------
C  Update 1 Three-dimensional Real Arrays's Ghost Cells using asynchronous
C  and persistent message-passing.
C
C  tjc  6/24/2002
C--------------------------------------------------------------------------
 
                             !..Pack Messages
      DO J=1,NEIGHPROC
         NCOUNT = 0
         DO I=1,NNODSEND(J)
            DO K=1,MNFEN
               NCOUNT = NCOUNT+1
               SENDBUF(NCOUNT,J)=VEC(ISENDLOC(I,J),K)
            ENDDO
         ENDDO
      ENDDO
              ! Send/receive messages to/from all neighbors
      DO J=1,NEIGHPROC  
         CALL MPI_IRECV( RECVBUF(1,J), MNFEN*NNODRECV(J), 
     &        REALTYPE,IPROC(J), TAG, COMM, REQ_R3D(J),IERR)
         CALL MPI_ISEND( SENDBUF(1,J), MNFEN*NNODSEND(J), 
     &        REALTYPE,IPROC(J), TAG, COMM, REQ_R3D(J+NEIGHPROC),IERR)
      ENDDO

              !..Unpack Received messages as they arrive     
      TOT = 0
      DO WHILE (TOT.LT.RDIM)
         DO N=1, RDIM
            INDX(N) = 0
         ENDDO
         CALL MPI_WAITSOME( RDIM,REQ_R3D,NFINI,INDX,STAT_R3D,IERR )
         TOT = TOT + NFINI
         DO N=1, NFINI
            IF (INDX(N).GT.0.AND.INDX(N).LE.RDIM)  THEN
              IF (INDX(N).LE.NEIGHPROC) THEN
                J = INDX(N)
                NCOUNT = 0
                DO I=1,NNODRECV(J)
                   DO K=1,MNFEN
                      NCOUNT = NCOUNT+1
                      VEC(IRECVLOC(I,J),K) = RECVBUF(NCOUNT,J)
                   ENDDO
                ENDDO
              ENDIF
            ENDIF
         ENDDO
      ENDDO
      RETURN
      END SUBROUTINE UPDATER3D


      SUBROUTINE UPDATEC3D( VEC )
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      COMPLEX(SZ), INTENT(INOUT) ::  VEC(MNP,MNFEN)
      INTEGER :: N,I,J,K,NCOUNT,NFINI,TOT
C--------------------------------------------------------------------------
C  Update 1 Three-dimensional Complex Arrays' Ghost Cells using asynchronous
C  and persistent message-passing.
C  tjc  6/24/2002
C--------------------------------------------------------------------------
  
                             !..Pack Messages
      DO J=1,NEIGHPROC
         NCOUNT = 0
         DO I=1,NNODSEND(J)
            DO K=1,MNFEN
               NCOUNT = NCOUNT+1
               SENDBUF(NCOUNT,J)=REAL(VEC(ISENDLOC(I,J),K))
               NCOUNT = NCOUNT+1
               SENDBUF(NCOUNT,J)=AIMAG(VEC(ISENDLOC(I,J),K))
            ENDDO
         ENDDO
      ENDDO
                     
              ! Send/receive messages to/from all neighbors
 
      DO J=1,NEIGHPROC  
         CALL MPI_IRECV( RECVBUF(1,J), 2*MNFEN*NNODRECV(J), 
     &        REALTYPE,IPROC(J), TAG, COMM, REQ_C3D(J),IERR)
C     jgf48.02 Commented out the following two lines and added the 
C     two after that to fix bug.
C         CALL MPI_SEND_INIT ( SENDBUF(1,J), 2*MNFEN*NNODSEND(J), 
C     &        REALTYPE,IPROC(J), TAG, COMM, REQ_C3D(J+NEIGHPROC),IERR)
         CALL MPI_ISEND ( SENDBUF(1,J), 2*MNFEN*NNODSEND(J), 
     &        REALTYPE,IPROC(J), TAG, COMM, REQ_C3D(J+NEIGHPROC),IERR)
      ENDDO

              !..Unpack Received messages as they arrive     
      TOT = 0
      DO WHILE (TOT.LT.RDIM)
         DO N=1, RDIM
            INDX(N) = 0
         ENDDO
         CALL MPI_WAITSOME( RDIM,REQ_C3D,NFINI,INDX,STAT_C3D,IERR )
         TOT = TOT + NFINI
         DO N=1, NFINI
            IF (INDX(N).GT.0.AND.INDX(N).LE.RDIM)  THEN
              IF (INDX(N).LE.NEIGHPROC) THEN
                J = INDX(N)
                NCOUNT = 0
                DO I=1,NNODRECV(J)
                   DO K=1,MNFEN
                      !jgf48.50 Make sure CMPLX returns the right KIND
#ifdef REAL4
                      VEC(IRECVLOC(I,J),K) = CMPLX(RECVBUF(NCOUNT+1,J),
     &                           RECVBUF(NCOUNT+2,J),4)
#else
                      VEC(IRECVLOC(I,J),K) = CMPLX(RECVBUF(NCOUNT+1,J),
     &                           RECVBUF(NCOUNT+2,J),8)
#endif
                      NCOUNT = NCOUNT+2
                   ENDDO
                ENDDO
              ENDIF
            ENDIF
         ENDDO
      ENDDO
C 
      RETURN
      END SUBROUTINE UPDATEC3D



      function psdot( n, sx, sy ) result(gsum)
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      integer n, i
      real(8) lsum,gsum
      real(sz) sx(*),sy(*)
      integer kount !jgf46.00 added
C------------------------------------------------------------------------------
C  Parallel version of SDOT for ITPACKV module
C------------------------------------------------------------------------------
 
      gsum = 0.0d0
      if (n.le.0) return
 
      lsum = 0.0D0
      do i = 1,n
         if (RESNODE(i)) then
            lsum = lsum  + sx(i)*sy(i)
         endif
      enddo
 
      kount = 1                              
      call MPI_ALLREDUCE( lsum, gsum, kount, DBLETYPE,
     &     MPI_SUM, COMM, ierr)
 
      return
      end function psdot


      subroutine ps2dots( n, sd, sdt ,dot3rray) 
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      integer n, i, idx
      real(8) lsum(2),gsum(2)
      real(sz) sd(*),sdt(*)
      real(sz) dot3rray(*)
      integer kount             !jgf46.00 added
C------------------------------------------------------------------------------
C  Parallel version of 2-SDOTs for ITPACKV module
C  jbr  6/17/00
C------------------------------------------------------------------------------
 
      dot3rray(1) = 0.0
      dot3rray(2) = 0.0
      dot3rray(3) = 0.0
      if (n.le.0) return
 
      lsum(1) = 0.0D0
      lsum(2) = 0.0D0

      do i = 1,n
         if (RESNODE(i)) lsum(1) = lsum(1)  + sd(i)*sd(i)
      enddo
      do i = 1,n
         if (RESNODE(i)) lsum(2) = lsum(2)  + sd(i)*sdt(i)
      enddo
      
      kount=2 
      call MPI_ALLREDUCE( lsum, gsum, kount, DBLETYPE,
     &     MPI_SUM, COMM, ierr)
 
      dot3rray(1) = gsum(1)
      dot3rray(2) = gsum(2)
      dot3rray(3) = 1.0
      return
      end subroutine ps2dots

      subroutine ps3dots( n, sd, sdt ,su,dot3rray)
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      integer n, i, idx
      real(8) lsum(3),gsum(3)
      real(sz) sd(*),sdt(*),su(*)
      real(sz) dot3rray(*)
      integer kount             ! jgf46.00 added
C------------------------------------------------------------------------------
C  Parallel version of 3-SDOTs for ITPACKV module
C  jbr  6/17/00
C------------------------------------------------------------------------------
 
      dot3rray(1) = 0.0
      dot3rray(2) = 0.0
      dot3rray(3) = 0.0
      if (n.le.0) return
 
      lsum(1) = 0.0D0
      lsum(2) = 0.0D0
      lsum(3) = 0.0D0

      do i = 1,n
         if (RESNODE(i)) lsum(1) = lsum(1)  + sd(i)*sd(i)
      enddo
      do i = 1,n
         if (RESNODE(i)) lsum(2) = lsum(2)  + sd(i)*sdt(i)
      enddo
      do i = 1,n
         if (RESNODE(i)) lsum(3) = lsum(3)  + su(i)*su(i)
      enddo
 
      kount=3 
      call MPI_ALLREDUCE( lsum, gsum, kount, DBLETYPE,
     &     MPI_SUM, COMM, ierr)
 
      dot3rray(1) = gsum(1)
      dot3rray(2) = gsum(2)
      dot3rray(3) = gsum(3)
      return
      end subroutine ps3dots


      SUBROUTINE ALLNODES( TOTNODES )
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      INTEGER I   
      INTEGER LNODES,TOTNODES
      INTEGER kount !jgf46.00 added
C------------------------------------------------------------------------------
C  Compute Number of nodes in entire domain.   
C------------------------------------------------------------------------------
 
      LNODES = 0
      DO I=1,MNP
         IF (RESNODE(I)) LNODES = LNODES + 1
      ENDDO
 
      kount=1  
      CALL MPI_ALLREDUCE(LNODES,TOTNODES,kount,MPI_INTEGER,MPI_SUM,
     &                   COMM,IERR)
      END SUBROUTINE ALLNODES


C------------------------------------------------------------------------------
C               S U B R O U T I N E   W E T D R Y S U M 
C------------------------------------------------------------------------------
C     jgf45.06 The GWCE left hand side must be reset whenever wetting or
C     drying has caused the grid to change. In parallel execution,
C     wetting or drying may occur in one subdomain but not another
C     during any particular time step. In this case, the NCChange flag
C     will be 1 in one subdomain but 0 in another, causing the
C     subdomains to get out of sync with each other's MPI calls. PADCIRC
C     will necessarily hang under these circumstances. This subroutine
C     sums the NCChange flags from all subdomains. If the sum comes back
C     greater than 0, they must all recompute the GWCE lhs, thus
C     preventing desynchronization of MPI communications.
C------------------------------------------------------------------------------
C
      SUBROUTINE WetDrySum( NCCHANGE )
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      INTEGER NCCHANGE    !input flag,=1 if this subdomain has wetted or dried
      INTEGER SumNCChange !sum total of all flags from all subdomains
      INTEGER kount       !jgf46.00 to avoid compiler bug on certain platforms

      SumNCChange = 0
      kount=1
      call MPI_ALLREDUCE( NCCHANGE, SumNCChange, kount, MPI_INTEGER,
     &     MPI_SUM, COMM, ierr)
      NCCHANGE = SumNCChange!resets GWCE for all subdomains if any s.d. resets
      RETURN
      END SUBROUTINE WetDrySum


C------------------------------------------------------------------------------
C               S U B R O U T I N E   W A R N  E L E V  S U M 
C------------------------------------------------------------------------------
C     jgf46.11 If the warning elevation was exceeded in one subdomain,
C     that information is propagated to the other subdomains so that the
C     velocities can be dumped to a file for debugging.
C------------------------------------------------------------------------------
      SUBROUTINE WarnElevSum( WarnElevExceeded )
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      INTEGER WarnElevExceeded !=1 if this subdomain has exceeded warning elev
      INTEGER SumWarnElevExceeded !sum total of all flags from all subdomains
      INTEGER kount       ! to avoid compiler bug on certain platforms

      SumWarnElevExceeded = 0
      kount=1
      call MPI_ALLREDUCE( WarnElevExceeded, SumWarnElevExceeded, kount,
     &     MPI_INTEGER, MPI_SUM, COMM, ierr)
      WarnElevExceeded = SumWarnElevExceeded
C------------------------------------------------------------------------------
      END SUBROUTINE WarnElevSum
C------------------------------------------------------------------------------

      SUBROUTINE MSG_IBCAST( array, n)
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      INTEGER :: n, array(:)
C------------------------------------------------------------------------------
C  Broadcast integer array from processor 0.
C  vjp  9/26/2006
C------------------------------------------------------------------------------
C
      CALL MPI_BCAST(array, n, mpi_integer,0,comm,ierr)
      RETURN
      END SUBROUTINE MSG_IBCAST


      SUBROUTINE MSG_LBCAST(array, n)
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      INTEGER :: n
      LOGICAL :: array(:)
      CALL MPI_BCAST(array, n, mpi_logical, 0, comm, ierr)
      RETURN
      END SUBROUTINE MSG_LBCAST


      SUBROUTINE MSG_CBCAST( msg, n)
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      INTEGER n
      CHARACTER(*) :: msg
C------------------------------------------------------------------------------
C  Broadcast integer array from processor 0.
C  vjp  9/26/2006
C------------------------------------------------------------------------------
C
      CALL MPI_BCAST(msg, n, mpi_character, 0, comm, ierr)
      RETURN
      END SUBROUTINE MSG_CBCAST


      SUBROUTINE MSG_RBCAST( array, n)
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      INTEGER n
      REAL(sz) array(:,:,:)
C------------------------------------------------------------------------------
C  Broadcast integer array from processor 0.
C  vjp  9/26/2006
C------------------------------------------------------------------------------
C
      CALL MPI_BCAST(array, n, realtype,0,comm,ierr)
      RETURN
      END SUBROUTINE MSG_RBCAST

Casey 090327: Implement Seizo's buffering of radiation stress broadcasts.
!=====Seizo Dividing Messenger
      SUBROUTINE MSG_RBCASTD( array, in, jn, kn)
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      INTEGER, INTENT(IN) :: in, jn, kn
      REAL(sz), INTENT(INOUT) :: array(in,jn,kn)
      INTEGER i, j, k, icount, ntimes, nlast
      REAL(sz), ALLOCATABLE :: buffer(:)
      INTEGER, ALLOCATABLE :: nbox(:)
      INTEGER, Parameter :: limit_buff=1000000
      allocate( buffer(in*jn*kn) )
!
      ntimes = in*jn*kn / limit_buff
      nlast = in*jn*kn - limit_buff * ntimes
      allocate( nbox(ntimes+1) )
      do i = 1, ntimes
         nbox(i) = limit_buff
      enddo
      if( nlast /= 0 ) then
        ntimes = ntimes + 1
        nbox(ntimes) = nlast
      endif
!
      icount = 0
      do i = 1, in
         do j = 1, jn
            do k = 1, kn
               icount = icount + 1
               buffer(icount) =  array(i,j,k)
            enddo
         enddo
      enddo
C
      do i = 1, ntimes
         icount = limit_buff*(i-1) + 1
         CALL MPI_BCAST(buffer(icount), nbox(i), realtype,0,comm,ierr)
      enddo
c
      icount = 0
      do i = 1, in
         do j = 1, jn
            do k = 1, kn
               icount = icount + 1
               array(i,j,k) = buffer(icount)
            enddo
         enddo
      enddo
      deallocate( buffer, nbox )
      RETURN
      END SUBROUTINE MSG_RBCASTD

      subroutine MSG_BARRIER()
C      use messenger
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      integer myerr
      call mpi_barrier(comm,myerr)
      end subroutine MSG_BARRIER
!=============================================================================

      END MODULE MESSENGER

      subroutine MSG_ABORT()
C      use messenger
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      integer myerr, myerrcode
      call mpi_abort(comm,myerrcode,myerr)
      end subroutine MSG_ABORT

