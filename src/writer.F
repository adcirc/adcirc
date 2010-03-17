C******************************************************************************
C MODULE WRITER
C  This module is for writer processors which receives values from compute 
C  processors and writes them out asynchronously.
C                                                              05 Feb 2007, sb
C******************************************************************************
C 
      MODULE WRITER
      USE SIZES
      USE GLOBAL, ONLY: C3D, COMM,
     &     BUFSIZE_MAX, BUF, resultBuf, OutputDataDescript_t, float_type
      USE MESSENGER, ONLY: TAG, COMM_WRITER, WRITER_ID,
     &     SIG_VAL, SIG_WRITE, SIG_TERM, 
     &     MSG_FINI, REALTYPE, DBLETYPE
#ifdef HAVE_MPI_MOD
      use mpi  
      IMPLICIT NONE
#else
      IMPLICIT NONE
#endif

 
C  variables

      INTEGER, PARAMETER :: NUM_BUF_MAX = 4
      INTEGER  :: NNOD_GLOBAL(NUM_BUF_MAX)
      INTEGER  :: NUM_ITEMS_PER_RECORD(NUM_BUF_MAX)
      INTEGER  :: ITS(NUM_BUF_MAX)
      REAL(SZ) :: DEFAULT_VALUES(NUM_BUF_MAX)
      REAL(SZ) :: TIMES(NUM_BUF_MAX)
      CHARACTER(256)  :: WFILENAMES(NUM_BUF_MAX)
      INTEGER :: BUFID, WCOMMID
      LOGICAL :: NEED_FLUSH
      
      REAL(SZ), SAVE, POINTER :: RESGBUF(:,:)
      REAL(SZ), ALLOCATABLE, TARGET :: RESGBUF1(:,:),RESGBUF2(:,:)
      REAL(SZ), ALLOCATABLE, TARGET :: RESGBUF3(:,:),RESGBUF4(:,:)

      CONTAINS

C---------------------end of data declarations--------------------------------C


      SUBROUTINE WRITER_INIT ()

      WCOMMID = 0
      NEED_FLUSH = .FALSE.
      
      RETURN
      END SUBROUTINE WRITER_INIT


      SUBROUTINE WRITER_MAIN ()

      USE GLOBAL,ONLY : ScreenUnit
      IMPLICIT NONE

#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif

      INTEGER :: SIG,IERR
      INTEGER :: STAT(MPI_STATUS_SIZE)

C...  Initializing
      BUFID = 0

C...  Main loop
      DO WHILE(.TRUE.)
        CALL MPI_RECV(SIG, 1, MPI_INTEGER, 0, TAG,
     &       COMM_WRITER(WRITER_ID),STAT,IERR)

        IF(SIG == SIG_VAL) THEN
          WRITE(ScreenUnit,*) 'WRITER PROC ',WRITER_ID,
     &         ' RECEIVED SIG_VAL'
          CALL WRITER_RECV_VALUES()

        ELSE IF(SIG == SIG_WRITE) THEN
          WRITE(ScreenUnit,*) 'WRITER PROC ',WRITER_ID,
     &         ' RECEIVED SIG_WRITE'
          CALL WRITER_WRITE()

        ELSE IF(SIG == SIG_TERM) THEN
          WRITE(ScreenUnit,*) 'WRITER PROC ',WRITER_ID,
     &         ' RECEIVED SIG_TERM'
          EXIT
        ENDIF
      ENDDO

      CALL MSG_FINI()
      STOP

      END SUBROUTINE WRITER_MAIN

 
      subroutine WRITE_GBL_FILE_THROUGH_WRITER
     &     (filename, descript, time, it, store_cmd, default)
      use global
      implicit none
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif

      external store_cmd
      integer      :: num
      integer      :: ierr
      integer      :: bufsize
      integer      :: istart, iend, it
      real(8)      :: time
      character(*) :: filename
      character(256) :: filenamebuf
      real(sz)     :: default
      integer      :: INTBUF(3)
      real(sz)     :: REBUF(2)
      integer      :: wcomm
      type (OutputDataDescript_t) :: descript


      NEED_FLUSH = .TRUE.

      CALL WRITER_GET_NEXT_WRITER_COMM(WCOMM)

#ifdef WRITER_DEBUG
      PRINT *, 'COMM_WRITE IS ', wcomm, ' (PROC ',MYPROC,')'
#endif

C...  Proc 0 does some special things
      if (myproc == 0) then

C...    Send SIG_VAL to a writer to get it ready to receive values.
        CALL MPI_SEND(SIG_VAL,1,MPI_INTEGER,MNPROC,TAG,
     &       wcomm,IERR)

C...    Send some parameters associated to this file
        INTBUF(1) = NP_G
        INTBUF(2) = descript % num_items_per_record
        INTBUF(3) = IT
        CALL MPI_SEND(INTBUF,3,MPI_INTEGER,MNPROC,TAG,
     &       WCOMM,IERR)

        REBUF(1) = default
        REBUF(2) = time
        CALL MPI_SEND(REBUF,2,REALTYPE,MNPROC,TAG,
     &       WCOMM,IERR)

        if(len_trim(filename).ge.256) then
          stop 'FILE NAME LENGTH MUST BE LESS THAN 256'
        endif

        filenamebuf = ""
        filenamebuf(1:len_trim(filename))=filename(1:len_trim(filename))

        CALL MPI_SEND(filenamebuf,256,MPI_CHARACTER,MNPROC,TAG,
     &       WCOMM,IERR)
      endif

C...  Compute the buffer size and others
      bufsize = min(BUFSIZE_MAX, 
     &    descript % num_items_per_record * descript % num_fd_records)
      num     = bufsize / descript % num_items_per_record
      iend    = num
      istart  = 1

      do while (istart < iend)

        !------------------------------------------------------------
        ! Initialize
        !------------------------------------------------------------
        buf(:)  = descript % initial_value

        call store_cmd(descript, istart, iend)

#ifdef WRITER_DEBUG
        PRINT *, 'PROC ', MYPROC, ' IS SENDING FOR ',
     &       TRIM(filename)
#endif

C...    The following mpi_reduce sends array values to rank MNPROC,
C...    which is a writer proc.
C...    Note that ranks from 0 to (MNPROC-1) are compute procs, and 
C...    rank MNPROC of writer communicators is always a writer proc.
        call mpi_reduce(buf, resultBuf, bufsize, float_type, MPI_SUM,
     &    MNPROC, WCOMM, ierr)

        istart = iend + 1
        iend   = min(istart + num - 1, descript % num_fd_records)
        num    = iend - istart + 1
      end do

      return
      end subroutine WRITE_GBL_FILE_THROUGH_WRITER

C...
C...  SUBROUTINE WRITER_STOREONE
C...  Partly adopted from globalio.F      
C...
      subroutine Writer_storeOne(descript, istart, iend)
      use global
      implicit none

      type (OutputDataDescript_t) :: descript
      integer :: i, istart, iend, iglobal
      integer :: ioffset

      ioffset = istart - 1
      do i = 1, descript % num_records_this
        iglobal = descript % imap(i)
        if (istart <= iglobal .and. iglobal <= iend) then
          buf(iglobal-ioffset) = descript % array(i)
        end if
      end do

      return
      end subroutine Writer_storeOne

C...
C...  SUBROUTINE WRITER_STORETWO
C...  Partly adopted from globalio.F      
C...
      subroutine Writer_storeTwo(descript, istart, iend)
      use global
      implicit none

      type (OutputDataDescript_t) :: descript
      integer :: i, istart, iend, iglobal
      integer :: ioffset, j

      ioffset = istart - 1
      do i = 1, descript % num_records_this
        iglobal = descript % imap(i)
        if (istart <= iglobal .and. iglobal <= iend) then
          j = 2*(iglobal - ioffset) - 1 
          buf(j)     = descript % array(i)
          buf(j + 1) = descript % array2(i)
        end if
      end do
      return
      end subroutine Writer_storeTwo

C...
C...  SUBROUTINE WRITER_STORE63
C...  Partly adopted from globalio.F      
C...
      subroutine Writer_store63(descript, istart, iend)
      use global
      implicit none

      type (OutputDataDescript_t) :: descript
      integer :: i, istart, iend, iglobal
      integer :: ioffset

      ioffset = istart - 1
      do i = 1, np
        iglobal = descript % imap(i)
        if (istart <= iglobal .and. iglobal <= iend) then
          if (nodecode(i) == 1) buf(iglobal - ioffset ) = eta2(i)
          if (nodecode(i) == 0) buf(iglobal - ioffset ) = -99999.0
        end if
      end do

      return
      end subroutine Writer_store63


      subroutine WRITER_RECV_VALUES()
      use global
      implicit none
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      integer      :: num, i, j
      integer      :: ierr
      integer      :: bufsize
      integer      :: istart, iend
      character(256) :: filenamebuf
      integer      :: intbuf(3)
      real(sz)     :: rebuf(2)
      integer      :: wcomm
      INTEGER      :: STAT(MPI_STATUS_SIZE)

C...  BUFID determines the buffer to store values in.
      BUFID = BUFID + 1

C...  Writers use their own communicator
      wcomm = comm_writer(writer_id)

C...  Receive some parameters necessary to write a file
      CALL MPI_RECV(INTBUF,3,MPI_INTEGER,0,TAG,
     &     WCOMM,STAT,IERR)
      NNOD_GLOBAL(BUFID) = INTBUF(1)
      NUM_ITEMS_PER_RECORD(BUFID) = INTBUF(2)
      ITS(BUFID) = INTBUF(3)
#ifdef WRITER_DEBUG
      print *, 'RECEIVED NNOD_GLOBAL = ',NNOD_GLOBAL(BUFID)
      print *, 'RECEIVED NUM_ITEMS_PER_REC = ',NUM_ITEMS_PER_RECORD(BUFID)
      print *, 'RECEIVED ITS = ',ITS(BUFID)
#endif

      IF((NUM_ITEMS_PER_RECORD(BUFID).NE.1).AND.
     &   (NUM_ITEMS_PER_RECORD(BUFID).NE.2)) THEN
        STOP 'WRITER: NUM_ITEMS_PER_RECORD MUST BE 1 OR 2'
      ENDIF

      CALL MPI_RECV(REBUF,2,REALTYPE,0,TAG,
     &     wcomm,STAT,IERR)
      DEFAULT_VALUES(BUFID) = REBUF(1)
      TIMES(BUFID) = REBUF(2)
#ifdef WRITER_DEBUG
      print *, 'RECEIVED DEFAULT = ',DEFAULT_VALUES(BUFID)
      print *, 'RECEIVED TIME = ',TIMES(BUFID)
#endif

      CALL MPI_RECV(filenamebuf,256,MPI_CHARACTER,0,TAG,
     &     wcomm,STAT,IERR)
#ifdef WRITER_DEBUG
      print *, 'RECEIVED FILENAME = ',TRIM(filenamebuf)
#endif
      WFILENAMES(BUFID) = ""
      WFILENAMES(BUFID)(1:LEN_TRIM(filenamebuf))
     &     = filenamebuf(1:LEN_TRIM(filenamebuf))

C...  Prepare buffer to store receiving values
      CALL WRITER_SET_RESGBUF(BUFID,.TRUE.)

C...  Compute the buffer size and others
      bufsize = min(BUFSIZE_MAX, 
     &    num_items_per_record(BUFID) * NNOD_GLOBAL(BUFID))
      num     = bufsize / num_items_per_record(BUFID)
      iend    = num
      istart  = 1

      do while (istart < iend)

#ifdef WRITER_DEBUG
c        PRINT *, 'PROC ', MYPROC, ' IS RECEIVING FOR ',
C     &       TRIM(WFILENAMES(BUFID))
#endif
        buf(:)  = 0.D0

C...    The following mpi_reduce sends array values to rank MNPROC,
C...    which is a writer proc.
C...    Note that ranks from 0 to (MNPROC-1) are compute procs, and 
C...    rank MNPROC of writer communicators is always a writer proc.
        call mpi_reduce(buf, resultBuf, bufsize, float_type, MPI_SUM,
     &    MNPROC, wcomm, ierr)

        j = 1
        do i = istart, iend
          resgbuf(1,i) = resultBuf(j)
          j = j+1

          if(NUM_ITEMS_PER_RECORD(BUFID).EQ.2) then
            resgbuf(2,i) = resultBuf(j)
            j = j+1
          endif
        end do

        istart = iend + 1
        iend   = min(istart + num - 1, NNOD_GLOBAL(BUFID))
        num    = iend - istart + 1
      end do

      end subroutine WRITER_RECV_VALUES

      SUBROUTINE WRITER_SET_RESGBUF(ID_BUF,DO_ALLOCATION)

      INTEGER ID_BUF
      LOGICAL DO_ALLOCATION

      IF(ID_BUF.EQ.1) THEN
        IF(.NOT.ALLOCATED(RESGBUF1)) THEN
          IF(DO_ALLOCATION) THEN
            ALLOCATE(RESGBUF1(2,NNOD_GLOBAL(ID_BUF)))
          ELSE
            GOTO 999
          ENDIF
        ENDIF
        RESGBUF => RESGBUF1

      ELSE IF(ID_BUF.EQ.2) THEN
        IF(.NOT.ALLOCATED(RESGBUF2)) THEN
          IF(DO_ALLOCATION) THEN
            ALLOCATE(RESGBUF2(2,NNOD_GLOBAL(ID_BUF)))
          ELSE
            GOTO 999
          ENDIF
        ENDIF
        RESGBUF => RESGBUF2

      ELSE IF(ID_BUF.EQ.3) THEN
        IF(.NOT.ALLOCATED(RESGBUF3)) THEN
          IF(DO_ALLOCATION) THEN
            ALLOCATE(RESGBUF3(2,NNOD_GLOBAL(ID_BUF)))
          ELSE
            GOTO 999
          ENDIF
        ENDIF
        RESGBUF => RESGBUF3

      ELSE IF(ID_BUF.EQ.4) THEN
        IF(.NOT.ALLOCATED(RESGBUF4)) THEN
          IF(DO_ALLOCATION) THEN
            ALLOCATE(RESGBUF4(2,NNOD_GLOBAL(ID_BUF)))
          ELSE
            GOTO 999
          ENDIF
        ENDIF
        RESGBUF => RESGBUF4

      ELSE
        STOP 'ID_BUF MUST BE 1-4.'

      ENDIF

      RETURN

 999  CONTINUE
      STOP 'WRITER.F: REQUESTED RESGBUF IS NOT ALLOCATED.'

      END SUBROUTINE WRITER_SET_RESGBUF

      SUBROUTINE FLUSH_WRITERS()
      use global
      implicit none
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      INTEGER I,IERR

      IF(.NOT.NEED_FLUSH) RETURN

      IF(MYPROC.eq.0) THEN
        DO I=1,MNWPROC
          WRITE(ScreenUnit,*) 'PROC ',MYPROC,
     &         ' IS SENDING SIG_WRITE TO WRITER ',I
          CALL MPI_SEND(SIG_WRITE,1,MPI_INTEGER,MNPROC,
     &         TAG,COMM_WRITER(I),IERR)
        ENDDO

      ENDIF
      
      NEED_FLUSH = .FALSE.

Casey 100125: Per Seizo's suggestion, comment out this line to allow
C             the cycling of writer cores.
C     WCOMMID = 0

      RETURN
      END SUBROUTINE FLUSH_WRITERS

      SUBROUTINE WRITER_WRITE()
      use global
      implicit none
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif

      INTEGER :: WCOMM
      INTEGER :: IBUF,INOD
      INTEGER :: nNonDefVal
      INTEGER, PARAMETER :: FID = 100
      REAL(SZ) :: VAL1,VAL2,DEFAULT

      WCOMM = COMM_WRITER(WRITER_ID)

#ifdef WRITER_DEBUG
      print *, 'WRITER ',WRITER_ID,' IS WRITING...'
#endif

      DO IBUF=1,BUFID

C...    Validity check
        IF((NUM_ITEMS_PER_RECORD(IBUF).NE.1).AND.
     &       (NUM_ITEMS_PER_RECORD(IBUF).NE.2)) THEN
          WRITE(ScreenUnit, *) 'WRITER.F : NUM_ITEMS_PER_NODE ',
     &         'MUST BE EITHER ONE OR TWO :', NUM_ITEMS_PER_RECORD(IBUF)
          STOP
        ENDIF

        DEFAULT = DEFAULT_VALUES(IBUF)

C...    Retrieve a RESGBUF corresponding to BUFID
        CALL WRITER_SET_RESGBUF(IBUF, .FALSE.)

C...    Open file
        OPEN(FID, FILE=WFILENAMES(IBUF), ACCESS='SEQUENTIAL',
     &       POSITION='APPEND')
        
C...    Count non-default values
        nNonDefVal = 0
C....   Scalar data
        IF(NUM_ITEMS_PER_RECORD(IBUF) == 1) THEN
          DO INOD=1,NNOD_GLOBAL(IBUF)
            VAL1 = RESGBUF(1,INOD)

            IF(VAL1 == DEFAULT) CYCLE

            nNonDefVal = nNonDefVal + 1
          ENDDO

C....   Vector data
        ELSE IF(NUM_ITEMS_PER_RECORD(IBUF) == 2) THEN
          DO INOD=1,NNOD_GLOBAL(IBUF)
            VAL1 = RESGBUF(1,INOD)
            VAL2 = RESGBUF(2,INOD)

            IF((VAL1 == DEFAULT).AND.
     &         (VAL2 == DEFAULT)) CYCLE

            nNonDefVal = nNonDefVal + 1
          ENDDO

        ENDIF
        
C...    Writes a time stamp
        write(fid, 100) TIMES(IBUF), ITS(IBUF), nNonDefVal, DEFAULT
 100    FORMAT(2x,1pE20.10,5X,I10,5X,I10,5X,1pE20.10)

C...    Write data
C....   Scalar data
        IF(NUM_ITEMS_PER_RECORD(IBUF) == 1) THEN
          DO INOD=1,NNOD_GLOBAL(IBUF)
            VAL1 = RESGBUF(1,INOD)

            IF(VAL1 == DEFAULT) CYCLE

            WRITE(FID,1000) INOD, VAL1
          ENDDO

C....   Vector data
        ELSE IF(NUM_ITEMS_PER_RECORD(IBUF) == 2) THEN
          DO INOD=1,NNOD_GLOBAL(IBUF)
            VAL1 = RESGBUF(1,INOD)
            VAL2 = RESGBUF(2,INOD)

            IF((VAL1 == DEFAULT).AND.
     &         (VAL2 == DEFAULT)) CYCLE

            WRITE(FID,2000) INOD, VAL1, VAL2
          ENDDO

        ENDIF
        
C...    Close file
        CLOSE(FID)
      ENDDO

      BUFID = 0
#ifdef WRITER_DEBUG
      print *, 'WRITER ',WRITER_ID,' FINISHED WRITING'
#endif

      RETURN
 1000 format(i8, 1x, 1pE20.10)
 2000 format(i8, 1x, 1pE20.10, 1pE20.10)

      END SUBROUTINE WRITER_WRITE

      SUBROUTINE WRITER_GET_NEXT_WRITER_COMM(NEXT_WCOMM)
      INTEGER :: NEXT_WCOMM
      INTEGER :: NEXT_WCOMMID
      
      IF(WCOMMID >= MNWPROC) THEN
        WCOMMID = 0
      ENDIF

      WCOMMID = WCOMMID + 1
      NEXT_WCOMMID = MOD(WCOMMID - 1,MNWPROC) + 1
      NEXT_WCOMM = COMM_WRITER(NEXT_WCOMMID)

      END SUBROUTINE WRITER_GET_NEXT_WRITER_COMM

      END MODULE WRITER

