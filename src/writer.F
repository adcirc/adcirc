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
     &     BUFSIZE_MAX, BUF, resultBuf, OutputDataDescript_t, float_type,
     &     DEBUG, ECHO, INFO, WARNING, ERROR, logMessage, allMessage,
     &     setMessageSource, unsetMessageSource
      USE MESSENGER, ONLY: TAG, COMM_WRITER, WRITER_ID,
     &     SIG_VAL, SIG_WRITE, SIG_TERM, SIG_PAUSE,
     &     MSG_FINI, REALTYPE, DBLETYPE
#ifdef HAVE_MPI_MOD
      use mpi
      IMPLICIT NONE
#else
      IMPLICIT NONE
#endif

C  variables

      INTEGER, PARAMETER :: NUM_BUF_MAX = 7      !st3 100708: Increase for SWAN OUTPUT
      INTEGER  :: NNOD_GLOBAL(NUM_BUF_MAX)
      INTEGER  :: NUM_ITEMS_PER_RECORD(NUM_BUF_MAX)
      INTEGER  :: ITS(NUM_BUF_MAX)
      INTEGER  :: FULLREC_SWITCH(NUM_BUF_MAX)  ! tcm v50.27 added for full record option
      INTEGER  :: NETCDF_SWITCH(NUM_BUF_MAX)   ! cobell ... for netcdf
      INTEGER  :: LUN(NUM_BUF_MAX)
      REAL(SZ) :: DEFAULT_VALUES(NUM_BUF_MAX)
      REAL(SZ) :: TIMES(NUM_BUF_MAX)
      CHARACTER(256)  :: WFILENAMES(NUM_BUF_MAX)
      INTEGER :: BUFID, WCOMMID
      LOGICAL :: NEED_FLUSH

      REAL(SZ), SAVE, POINTER :: RESGBUF(:,:)
      REAL(SZ), ALLOCATABLE, TARGET :: RESGBUF1(:,:),RESGBUF2(:,:)
      REAL(SZ), ALLOCATABLE, TARGET :: RESGBUF3(:,:),RESGBUF4(:,:)
      REAL(SZ), ALLOCATABLE, TARGET :: RESGBUF5(:,:),RESGBUF6(:,:) !st3 100708: Increase for SWAN OUTPUT
      REAL(SZ), ALLOCATABLE, TARGET :: RESGBUF7(:,:)               !st3 100708: Increase for SWAN OUTPUT
!
      INTEGER, ALLOCATABLE :: NWLOADED(:)   !st3 100708: Round Robin Chamber

      CONTAINS

C---------------------end of data declarations--------------------------------C


C-----------------------------------------------------------------------
C     S U B R O U T I N E   W R I T E R  _  I N I T
C-----------------------------------------------------------------------
C     Initializes round robin chamber; clears mpi communicator id for
C     this processor; clears flush flag.
C-----------------------------------------------------------------------
      SUBROUTINE WRITER_INIT ()
      IMPLICIT NONE
      call setMessageSource("writer_init")
#if defined(WRITER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      WCOMMID = 0
      NEED_FLUSH = .FALSE.
      ALLOCATE ( NWLOADED(1:MNWPROC) )   !st3 100708: Allocate Round Robin Chamber
      NWLOADED(1:MNWPROC) = 0            !st3 100708: Zeroc    Round Robin Chamber

#if defined(WRITER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN
C-----------------------------------------------------------------------
      END SUBROUTINE WRITER_INIT
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     S U B R O U T I N E   W R I T E R  _  M A I N
C-----------------------------------------------------------------------
C     Writer processor goes into an infinite loop where it waits for
C     either a signal to receive data, or a signal to write output.
C     When it receives a signal, it performs the specified action,
C     then goes back to waiting for another signal.
C-----------------------------------------------------------------------
      SUBROUTINE WRITER_MAIN ()
      USE GLOBAL,ONLY : ScreenUnit
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      INTEGER :: SIG,IERR
      INTEGER :: STAT(MPI_STATUS_SIZE)
C
      call setMessageSource("writer_main")
#if defined(WRITER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
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

! tcm addeded v49.52.01 SIG_PAUSE
        ELSE IF(SIG == SIG_PAUSE) THEN
           WRITE(ScreenUnit,*) 'WRITER PROC ',WRITER_ID,
     &         ' RECEIVED SIG_PAUSE'
           EXIT

        ELSE IF(SIG == SIG_TERM) THEN
          WRITE(ScreenUnit,*) 'WRITER PROC ',WRITER_ID,
     &         ' RECEIVED SIG_TERM'

          EXIT
        ENDIF
      ENDDO

! tcm commented out v49.52.01
!      CALL MSG_FINI()
!      STOP
#if defined(WRITER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE WRITER_MAIN
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E   W R I T E R  _  P A U S E
C-----------------------------------------------------------------------
C     Sends a signal to all writer processors to pause, which just
C     means that they exit their infinite loop and return.
C-----------------------------------------------------------------------
      subroutine WRITER_PAUSE()
      use global
      use messenger
      implicit none
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      integer i
C
      call setMessageSource("writer_pause")
#if defined(WRITER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
C
      DO I=1,MNWPROC
         CALL MPI_SEND( SIG_PAUSE, 1, MPI_INTEGER, MNPROC,
     &                  TAG,COMM_WRITER(I), IERR)
      ENDDO
C
C
#if defined(WRITER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      return
C-----------------------------------------------------------------------
      end subroutine WRITER_PAUSE
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     S U B R O U T I N E
C        W R I T E _ G B L _ F I L E _ T H R O U G H _ W R I T E R
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      subroutine WRITE_GBL_FILE_THROUGH_WRITER
     &     (tlun, filename, descript, TimeLoc, it, store_cmd, default_val)
      use global
      implicit none
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      external store_cmd
      integer      :: tlun !...Cobell-For Netcdf unit numbers
      integer      :: num
      integer      :: ierr
      integer      :: bufsize
      integer      :: istart, iend, it
      integer      :: fullrec   !1 = true, 0 = false   tcm v50.27 added
      integer      :: w_netcdf
      real(8)      :: TimeLoc
      character(*) :: filename
      character(256) :: filenamebuf
      real(sz)     :: default_val
      integer      :: INTBUF(6)  !tcm v50.27 increased size from 3 to 4, cobell increase to 6
      real(sz)     :: REBUF(2)
      integer      :: wcomm
      type (OutputDataDescript_t) :: descript
C
      call setMessageSource("write_gbl_file_through_writer")
#if defined(WRITER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      NEED_FLUSH = .TRUE.
      CALL WRITER_GET_NEXT_WRITER_COMM(WCOMM)

#ifdef WRITER_DEBUG
      PRINT *, 'COMM_WRITE IS ', wcomm, ' (PROC ',MYPROC,')'
#endif

C...  Proc 0 does some special things
      if (myproc == 0) then

C...     Send SIG_VAL to a writer to get it ready to receive values.
         CALL MPI_SEND(SIG_VAL,1,MPI_INTEGER,MNPROC,TAG,
     &       wcomm,IERR)

C...     Send some parameters associated to this file
C....    tcm v50.27 added for full/compact record
         FULLREC = 1  !full record output
         W_NETCDF = 0
         IF (ABS(descript % specifier) == 4) THEN
            !compact record
            FULLREC = 0
Cobell - Allow writer processors to output netcdf data            
         ELSEIF ( ABS(descript % specifier) == 3 .OR. 
     &            ABS(descript % specifier) == 5 ) THEN 
            W_NETCDF = 1
         ENDIF

C
         INTBUF(1) = NP_G
         INTBUF(2) = descript % num_items_per_record
         INTBUF(3) = IT
         INTBUF(4) = FULLREC   !tcm v50.27 added
         INTBUF(5) = W_NETCDF  !cobell added
         INTBUF(6) = tlun
         CALL MPI_SEND(INTBUF,6,MPI_INTEGER,MNPROC,TAG,
     &       WCOMM,IERR)

         REBUF(1) = default_val
         REBUF(2) = TimeLoc
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

C...     The following mpi_reduce sends array values to rank MNPROC,
C...     which is a writer proc.
C...     Note that ranks from 0 to (MNPROC-1) are compute procs, and
C...     rank MNPROC of writer communicators is always a writer proc.
         call mpi_reduce(buf, resultBuf, bufsize, float_type, MPI_SUM,
     &     MNPROC, WCOMM, ierr)
         istart = iend + 1
         iend   = min(istart + num - 1, descript % num_fd_records)
         num    = iend - istart + 1
      end do
C
#if defined(WRITER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      return
C-----------------------------------------------------------------------
      end subroutine WRITE_GBL_FILE_THROUGH_WRITER
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     S U B R O U T I N E     W R I T E R _ S T O R E   O N E
C-----------------------------------------------------------------------
C...  Partly adopted from globalio.F
C-----------------------------------------------------------------------
      subroutine Writer_storeOne(descript, istart, iend)
      use global
      implicit none

      type (OutputDataDescript_t) :: descript
      integer :: i, istart, iend, iglobal
      integer :: ioffset

      call setMessageSource("writer_storeone")
#if defined(WRITER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif


      ioffset = istart - 1
      do i = 1, descript % num_records_this
        iglobal = descript % imap(i)
        if (istart <= iglobal .and. iglobal <= iend) then
          buf(iglobal-ioffset) = descript % array(i)
        end if
      end do


#if defined(WRITER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      return
C-----------------------------------------------------------------------
      end subroutine Writer_storeOne
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     S U B R O U T I N E     W R I T E R _ S T O R E   T W O
C-----------------------------------------------------------------------
C...  Partly adopted from globalio.F
C-----------------------------------------------------------------------
      subroutine Writer_storeTwo(descript, istart, iend)
      use global
      implicit none

      type (OutputDataDescript_t) :: descript
      integer :: i, istart, iend, iglobal
      integer :: ioffset, j
C
      call setMessageSource("writer_storetwo")
#if defined(WRITER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      ioffset = istart - 1
      do i = 1, descript % num_records_this
         iglobal = descript % imap(i)
         if (istart <= iglobal .and. iglobal <= iend) then
            j = 2*(iglobal - ioffset) - 1
            buf(j)     = descript % array(i)
            buf(j + 1) = descript % array2(i)
         end if
      end do

#if defined(WRITER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      return
C-----------------------------------------------------------------------
      end subroutine Writer_storeTwo
C-----------------------------------------------------------------------



C-----------------------------------------------------------------------
C     S U B R O U T I N E     W R I T E R _ S T O R E   6 3
C-----------------------------------------------------------------------
C...  Partly adopted from globalio.F
C-----------------------------------------------------------------------
      subroutine Writer_store63(descript, istart, iend)
      use global
      implicit none

      type (OutputDataDescript_t) :: descript
      integer :: i, istart, iend, iglobal
      integer :: ioffset

      call setMessageSource("writer_store63")
#if defined(WRITER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      ioffset = istart - 1
      do i = 1, np
         iglobal = descript % imap(i)
         if (istart <= iglobal .and. iglobal <= iend) then
            if (nodecode(i) == 1) buf(iglobal - ioffset ) = eta2(i)
            if (nodecode(i) == 0) buf(iglobal - ioffset ) = -99999.0
         end if
      end do

#if defined(WRITER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      return
C-----------------------------------------------------------------------
      end subroutine Writer_store63
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     S U B R O U T I N E     W R I T E R _ R E C V _ V A L U E S
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
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
      integer      :: intbuf(6)  !tcm v50.27 increased size from 3 to 4, cobell to 6
      real(sz)     :: rebuf(2)
      integer      :: wcomm
      INTEGER      :: STAT(MPI_STATUS_SIZE)
      character(len=8) :: cdigit              !st3 100708: split file

      call setMessageSource("writer_recv_values")
#if defined(WRITER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

C...  BUFID determines the buffer to store values in.
      BUFID = BUFID + 1

C...  Writers use their own communicator
      wcomm = comm_writer(writer_id)

C...  Receive some parameters necessary to write a file
      CALL MPI_RECV(INTBUF,6,MPI_INTEGER,0,TAG,
     &     WCOMM,STAT,IERR)
      NNOD_GLOBAL(BUFID) = INTBUF(1)
      NUM_ITEMS_PER_RECORD(BUFID) = INTBUF(2)
      ITS(BUFID) = INTBUF(3)
      FULLREC_SWITCH(BUFID) = INTBUF(4)  !tcm v50.27 added
      NETCDF_SWITCH(BUFID) = INTBUF(5)   !cobell added for netcdf
      LUN(BUFID) = INTBUF(6)

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
      if( isplit /= 0 ) then                                    !st3 100708: split file
         write(cdigit,'(i8.8)') NINT(TIMES(BUFID))              !st3 split files are named
         filenamebuf = trim(adjustl(filenamebuf)) //'.'//cdigit !st3 by integer part of
      endif                                                     !st3 simulation time
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

C...     The following mpi_reduce sends array values to rank MNPROC,
C...     which is a writer proc.
C...     Note that ranks from 0 to (MNPROC-1) are compute procs, and
C...     rank MNPROC of writer communicators is always a writer proc.
         call mpi_reduce(buf, resultBuf, bufsize, float_type, MPI_SUM,
     &     MNPROC, wcomm, ierr)

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

#if defined(WRITER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      return
C-----------------------------------------------------------------------
      end subroutine WRITER_RECV_VALUES
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     S U B R O U T I N E     W R I T E R _ S E T _ R E S G B U F
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE WRITER_SET_RESGBUF(ID_BUF,DO_ALLOCATION)
      IMPLICIT NONE
      INTEGER ID_BUF
      LOGICAL DO_ALLOCATION

      call setMessageSource("writer_set_resgbuf")
#if defined(WRITER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

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

!st3 100708: Increase for SWAN OUTPUT
#ifdef CSWAN
      ELSE IF(ID_BUF.EQ.5) THEN
         IF(.NOT.ALLOCATED(RESGBUF5)) THEN
            IF(DO_ALLOCATION) THEN
               ALLOCATE(RESGBUF5(2,NNOD_GLOBAL(ID_BUF)))
            ELSE
               GOTO 999
            ENDIF
         ENDIF
         RESGBUF => RESGBUF5

      ELSE IF(ID_BUF.EQ.6) THEN
         IF(.NOT.ALLOCATED(RESGBUF6)) THEN
            IF(DO_ALLOCATION) THEN
               ALLOCATE(RESGBUF6(2,NNOD_GLOBAL(ID_BUF)))
            ELSE
               GOTO 999
            ENDIF
         ENDIF
         RESGBUF => RESGBUF6

      ELSE IF(ID_BUF.EQ.7) THEN
         IF(.NOT.ALLOCATED(RESGBUF7)) THEN
            IF(DO_ALLOCATION) THEN
               ALLOCATE(RESGBUF7(2,NNOD_GLOBAL(ID_BUF)))
            ELSE
               GOTO 999
            ENDIF
         ENDIF
         RESGBUF => RESGBUF7
#endif

      ELSE
         STOP 'ID_BUF MUST BE 1-7.'
      ENDIF

#if defined(WRITER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN

 999  CONTINUE
      STOP 'WRITER.F: REQUESTED RESGBUF IS NOT ALLOCATED.'
C-----------------------------------------------------------------------
      END SUBROUTINE WRITER_SET_RESGBUF
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     S U B R O U T I N E     F L U S H _ W R I T E R S
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE FLUSH_WRITERS()
      use global
      implicit none
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      INTEGER I,IERR

      call setMessageSource("flush_writers")
#if defined(WRITER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      IF(.NOT.NEED_FLUSH) THEN
#if defined(WRITER_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG,"Return.")
#endif
         call unsetMessageSource()
         RETURN
      ENDIF

!st3 100708: If split option is off, load the blank cartridge to all chamber.
      IF( ISPLIT == 0 ) THEN
         NWLOADED(1:MNWPROC) = 1
      ENDIF
      IF(MYPROC.eq.0) THEN
         DO I=1,MNWPROC
!st3 100708: Flash only loaded chamber
            IF( NWLOADED(I) == 1 ) then    !st3 NWLOADED()=0:unloaded, 1:loaded
               WRITE(ScreenUnit,*) 'PROC ',MYPROC,
     &           ' IS SENDING SIG_WRITE TO WRITER ',I
               CALL MPI_SEND(SIG_WRITE,1,MPI_INTEGER,MNPROC,
     &              TAG,COMM_WRITER(I),IERR)
            ENDIF
            NWLOADED(I) = 0
         ENDDO
      ENDIF

      NEED_FLUSH = .FALSE.

Casey 100125: Per Seizo's suggestion, comment out this line to allow
C             the cycling of writer cores.
C     WCOMMID = 0
#if defined(WRITER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN
C-----------------------------------------------------------------------
      END SUBROUTINE FLUSH_WRITERS
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     S U B R O U T I N E     W R I T E R _ W R I T E
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE WRITER_WRITE()
      use global
#ifdef ADCNETCDF
      use NETCDFIO,only:writeOutArrayNetCDFWriter
#endif      
      implicit none
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      INTEGER :: WCOMM,IERR
      INTEGER :: IBUF,INOD,I,J
      INTEGER :: MYRECORD,MYLUN
      INTEGER :: nNonDefVal
      INTEGER :: FULLREC  !tcm v50.27 added
      INTEGER :: W_NETCDF
      INTEGER, PARAMETER :: FID = 100
      REAL(SZ) :: VAL1,VAL2,DEFAULT_VAL

      call setMessageSource("writer_write")
#if defined(WRITER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      
      WCOMM = COMM_WRITER(WRITER_ID)

#ifdef WRITER_DEBUG
      print *, 'WRITER ',WRITER_ID,' IS WRITING...'
#endif

      DO IBUF=1,BUFID

C...     Validity check
         IF((NUM_ITEMS_PER_RECORD(IBUF).NE.1).AND.
     &       (NUM_ITEMS_PER_RECORD(IBUF).NE.2)) THEN
            WRITE(ScreenUnit, *) 'WRITER.F : NUM_ITEMS_PER_NODE ',
     &         'MUST BE EITHER ONE OR TWO :', NUM_ITEMS_PER_RECORD(IBUF)
            STOP
         ENDIF

         DEFAULT_VAL = DEFAULT_VALUES(IBUF)
         FULLREC = FULLREC_SWITCH(IBUF)
         W_NETCDF = NETCDF_SWITCH(IBUF) 

C...     Retrieve a RESGBUF corresponding to BUFID
         CALL WRITER_SET_RESGBUF(IBUF, .FALSE.)

C...     Open file
         IF(W_NETCDF.NE.1)THEN
            if(isplit == 0 ) then             !Not Split               !st3  split file 100708
            OPEN(FID, FILE=WFILENAMES(IBUF), ACCESS='SEQUENTIAL',
     &       POSITION='APPEND')
            else                              !Split                   !st3  split file 100708
                OPEN(FID, FILE=WFILENAMES(IBUF), STATUS='REPLACE',       !st3  split file 100708
     &              ACTION='WRITE')                                     !st3  split file 100708
            endif
         ENDIF   

         IF(W_NETCDF.NE.1)THEN
c... tcm v50.27 added full record writing options
            IF(FULLREC == 1) THEN
                nNonDefVal = NNOD_GLOBAL(IBUF)
            ELSE
C...        Count non-default values
                nNonDefVal = 0
C....       Scalar data
                IF(NUM_ITEMS_PER_RECORD(IBUF) == 1) THEN
                    DO INOD=1,NNOD_GLOBAL(IBUF)
                        VAL1 = RESGBUF(1,INOD)
                        IF(VAL1 == DEFAULT_VAL) CYCLE
                        nNonDefVal = nNonDefVal + 1
                    ENDDO
C....   Vector data
                ELSE IF(NUM_ITEMS_PER_RECORD(IBUF) == 2) THEN
                    DO INOD=1,NNOD_GLOBAL(IBUF)
                        VAL1 = RESGBUF(1,INOD)
                        VAL2 = RESGBUF(2,INOD)

                        IF((VAL1 == DEFAULT_VAL).AND.
     &                  (VAL2 == DEFAULT_VAL)) CYCLE

                        nNonDefVal = nNonDefVal + 1
                    ENDDO
                ENDIF
            ENDIF

C...        Writes a time stamp
c.          tcm v50.27 added for full record option
            IF (FULLREC == 0) THEN
                WRITE(FID, 100) TIMES(IBUF), ITS(IBUF), nNonDefVal, DEFAULT_VAL
            ELSE
                WRITE(FID, 110) TIMES(IBUF), ITS(IBUF)
            ENDIF
 100        FORMAT(2x,1pE20.10E3,5X,I10,5X,I10,5X,1pE20.10E3)  !changed formats to match those in globalio.F
 110        FORMAT(2x,1pE20.10E3,5X,I10)  !tcm v50.27 added for full record

C...        Write data
C....       Scalar data
            IF(NUM_ITEMS_PER_RECORD(IBUF) == 1) THEN
                DO INOD=1,NNOD_GLOBAL(IBUF)
                    VAL1 = RESGBUF(1,INOD)

                    IF((FULLREC==0).AND.(VAL1 == DEFAULT_VAL)) CYCLE

                    WRITE(FID,1000) INOD, VAL1
                ENDDO

C....   Vector data
            ELSE IF(NUM_ITEMS_PER_RECORD(IBUF) == 2) THEN
                DO INOD=1,NNOD_GLOBAL(IBUF)
                    VAL1 = RESGBUF(1,INOD)
                    VAL2 = RESGBUF(2,INOD)

                    IF((FULLREC==0).AND.((VAL1 == DEFAULT_VAL).AND.
     &                  (VAL2 == DEFAULT_VAL))) CYCLE

                    WRITE(FID,2000) INOD, VAL1, VAL2
                ENDDO

            ENDIF

C...        Close file
            CLOSE(FID)
         ELSE
Cobell ... Added to allow writer core to write NetCDF output without
C          needing access to the entire descript array which would have
C          required extra message passing steps.
#ifdef ADCNETCDF         
            IF(NUM_ITEMS_PER_RECORD(IBUF).EQ.1)THEN
                CALL writeOutArrayNetCDFWriter(lun(ibuf), times(ibuf), 
     &              its(ibuf),1,default_values(ibuf),RESGBUF(1,:) )
            ELSEIF(NUM_ITEMS_PER_RECORD(IBUF).EQ.2)THEN
                CALL writeOutArrayNetCDFWriter(lun(ibuf), times(ibuf),
     &              its(ibuf),2,default_values(ibuf),RESGBUF(1,:),
     &              RESGBUF(2,:) )
            ENDIF
#else
         WRITE(ScreenUnit,'(A)') "ERROR: NetCDF is not available. "//
     &                           "No output has been written."    
#endif         
         ENDIF
      ENDDO

      BUFID = 0
#ifdef WRITER_DEBUG
      print *, 'WRITER ',WRITER_ID,' FINISHED WRITING'
#endif

#if defined(WRITER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN
C.... tcm v50.27 changed formats to match those in globalio.F
 1000 FORMAT(2X,I8, 2X, 1PE20.10E3) !format(i8, 1x, 1pE20.10)
 2000 FORMAT(2X,I8, 2X, 1PE20.10E3, 1PE20.10E3) !format(i8, 1x, 1pE20.10, 1pE20.10)
C-----------------------------------------------------------------------
      END SUBROUTINE WRITER_WRITE
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     S U B R O U T I N E
C         W R I T E R _ G E T _ N E X T _ W R I T E R _ C O M M
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE WRITER_GET_NEXT_WRITER_COMM(NEXT_WCOMM)
      IMPLICIT NONE
      INTEGER :: NEXT_WCOMM
      INTEGER :: NEXT_WCOMMID

      call setMessageSource("writer_get_next_writer_comm")
#if defined(WRITER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      IF(WCOMMID >= MNWPROC) THEN
        WCOMMID = 0
      ENDIF

      WCOMMID = WCOMMID + 1
      NEXT_WCOMMID = MOD(WCOMMID - 1,MNWPROC) + 1
      NEXT_WCOMM = COMM_WRITER(NEXT_WCOMMID)
      NWLOADED(NEXT_WCOMMID) = 1                    !st3 100708: Load the bullet to Round Robin Chamber

#if defined(WRITER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN
C-----------------------------------------------------------------------
      END SUBROUTINE WRITER_GET_NEXT_WRITER_COMM
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      END MODULE WRITER
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
!
!
!
!!
C_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
C MODULE HSWRITER  !st3
C  This module is for writer processors for hotstart files based on S.Bunya's
C_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
C
      MODULE HSWRITER
        USE SIZES
        USE GLOBAL
        USE MESSENGER, ONLY: TAG, COMM_WRITEH, COMM_HSLEEP, WRITER_ID,
     &       SIG_WRITE, SIG_TERM, SIG_PAUSE, MSG_FINI
#ifdef HAVE_MPI_MOD
        use mpi
#endif
      IMPLICIT NONE
      INTEGER :: WCOMMID, WCOMM
      CONTAINS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!?

!------------------------------------------------------------
      SUBROUTINE HSWRITER_INIT ()
        WCOMMID = 0
      END SUBROUTINE HSWRITER_INIT
!
!------------------------------------------------------------
!
      SUBROUTINE HSWRITER_MAIN ()
        USE GLOBAL,ONLY : ScreenUnit
        IMPLICIT NONE
#ifndef HAVE_MPI_MOD
        include 'mpif.h'
#endif
        INTEGER :: SIG,IERR
        INTEGER :: STAT(MPI_STATUS_SIZE)

        call setMessageSource("HSWRITER_MAIN")
!
        WCOMMID=WRITER_ID
        WCOMM=COMM_WRITEH(WCOMMID)
C...  Main loop
      DO WHILE(.TRUE.)
        CALL MPI_BARRIER(COMM_HSLEEP(WRITER_ID),IERR) !st3
        CALL MPI_RECV(SIG, 1, MPI_INTEGER, 0, TAG,
     &       COMM_WRITEH(WRITER_ID),STAT,IERR)

        IF(SIG == SIG_WRITE) THEN
          WRITE(ScreenUnit,*) 'HSWRITER PROC ',WRITER_ID,
     &         ' RECEIVED SIG_WRITE'
          CALL writeHotstart_through_HSwriter(0.0d0, 0)

!tcm v49.52.01 added SIG_PAUSE
        ELSE IF(SIG == SIG_PAUSE) then
           WRITE(ScreenUnit,*) 'HSWRITER PROC ',WRITER_ID,
     &         ' RECEIVED SIG_PAUSE'
           EXIT

        ELSE IF(SIG == SIG_TERM) THEN
          WRITE(ScreenUnit,*) 'HSWRITER PROC ',WRITER_ID,
     &         ' RECEIVED SIG_TERM'
          EXIT
        ENDIF
      ENDDO
!tcm v49.52.01 commented out next two lines
!      CALL MSG_FINI()
!      STOP

      call  unsetMessageSource()
      
      END SUBROUTINE HSWRITER_MAIN
!
!------------------------------------------------------------
! tcm v49.52.01 added the HSWRITER_PAUSE subroutine
!
      subroutine HSWRITER_PAUSE()
      use global
      use messenger
      implicit none
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      integer i

      DO I=1,MNWPROH
         CALL MPI_SEND( SIG_PAUSE, 1, MPI_INTEGER, MNPROC,
     &                  TAG,COMM_WRITEH(I), IERR)
      ENDDO

      return
      end subroutine HSWRITER_PAUSE
!
!------------------------------------------------------------

C
C-----------------------------------------------------------------------
C     SUBROUTINE BASED ON "WRITEHOTSTART"
C-----------------------------------------------------------------------
C     This subroutine was copied from write_output.F writeHotstart
C       and fitted to writer core algorithm  !st3 05.15.2010
C
C      ******* THIS SUBROUTINE DOES NOT DEAL WITH NETCDF
C-----------------------------------------------------------------------
      SUBROUTINE writeHotstart_through_HSwriter(TimeLoc, IT)

      USE SIZES
      USE GLOBAL
      USE GLOBAL_IO, ONLY: packOne, unpackOne, packOneInt, unpackOneInt,
     &    packTwo, unpackTwo, HEADER_MAX
!     &    collectFullDomainArray, collectFullDomainIntArray
      USE MESSENGER
!
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      INTEGER :: STAT(MPI_STATUS_SIZE)
      integer :: i

      REAL(8), intent(in) :: TimeLoc
      INTEGER, intent(in) :: IT
!
      double precision :: ttime, dpack(1)  !st3
      integer          :: iit,   ipack(26)  !st3
      character(len=2048) ::     cjunk1  !st3
      character(len=20) ::       cjunk2  !st3  tcm v51.26 changed from 10 to 20
      integer :: iNHSINC  !tcm v51.26
!
      type(OutputDataDescript_t) :: descript
      INTEGER npx, nex
      CHARACTER*9 :: itstr  !tcm v51.26 added for time-stamped hot start file name adjustments

      type(OutputDataDescript_t), SAVE :: Elev1Descript
      type(OutputDataDescript_t), SAVE :: Elev2Descript
      type(OutputDataDescript_t), SAVE :: VelDescript
      type(OutputDataDescript_t), SAVE :: CH1Descript
      type(OutputDataDescript_t), SAVE :: EtaDiscDescript
      type(OutputDataDescript_t), SAVE :: NodeCodeDescript
      type(OutputDataDescript_t), SAVE :: NOFFDescript
      LOGICAL,  SAVE :: FirstCall = .true.
!

      call setMessageSource("writeHotstart_through_HSwriter")

      if( myproc <= MNPROC-1 ) then
        CALL GET_NEXT_HSWRITER_COMM(WCOMM)
      endif
!
      iit = 0
      ttime=0.0d0
!
      if( myproc==0 ) then
        CALL MPI_BARRIER(COMM_HSLEEP(WCOMMID),IERR) !st3
        CALL MPI_SEND(SIG_WRITE,1,MPI_INTEGER,MNPROC,TAG,WCOMM,IERR)
        ipack( 1) = IT
        ipack( 2) = NHSTAR
        ipack( 3) = IM
        ipack( 4) = NP_G
        ipack( 5) = NE_G
        ipack( 6) = hss%lun
        ipack( 7) = FileFmtVersion
        ipack( 8) = IESTP
        ipack( 9) = NSCOUE
        ipack(10) = IVSTP
        ipack(11) = NSCOUV
        ipack(12) = ICSTP
        ipack(13) = NSCOUC
        ipack(14) = IPSTP
        ipack(15) = IWSTP
        ipack(16) = NSCOUM
        ipack(17) = IGEP
        ipack(18) = NSCOUGE
        ipack(19) = IGVP
        ipack(20) = NSCOUGV
        ipack(21) = IGCP
        ipack(22) = NSCOUGC
        ipack(23) = IGPP
        ipack(24) = IGWP
        ipack(25) = NSCOUGW
        ipack(26) = NHSINC  !tcm v51.26

        dpack(1) = TimeLoc
        cjunk1 = trim(adjustl(HOTSTARTDIR))
        cjunk2 = trim(adjustl(hss%filename))

        CALL MPI_SEND(ipack,26,MPI_INTEGER,MNPROC,TAG,WCOMM,IERR)
        CALL MPI_SEND(dpack,1,MPI_DOUBLE_PRECISION,MNPROC,TAG,WCOMM,IERR)
        CALL MPI_SEND(WRITE_LOCAL_HOT_START_FILES,1,MPI_LOGICAL,
     &                MNPROC,TAG,WCOMM,IERR)
        CALL MPI_SEND(cjunk1,2048,MPI_CHARACTER,MNPROC,TAG,WCOMM,IERR)
        CALL MPI_SEND(cjunk2,  20,MPI_CHARACTER,MNPROC,TAG,WCOMM,IERR)  !tcm v51.26 changed to length 20
      endif

      if( WRITER_ID /= 0 ) then
        CALL MPI_RECV(ipack,26,MPI_INTEGER,0,TAG,WCOMM,STAT,IERR)
        CALL MPI_RECV(dpack,1,MPI_DOUBLE_PRECISION,0,TAG,WCOMM,STAT,IERR)
        CALL MPI_RECV(WRITE_LOCAL_HOT_START_FILES,1,MPI_LOGICAL,
     &                0, TAG,WCOMM,STAT,IERR)
        CALL MPI_RECV(cjunk1,2048,MPI_CHARACTER,0,TAG,WCOMM,STAT,IERR)
        CALL MPI_RECV(cjunk2,  20,MPI_CHARACTER,0,TAG,WCOMM,STAT,IERR)  !tcm v51.26 changed from 10 to 20
        iit            = ipack( 1)
        NHSTAR         = ipack( 2)
        IM             = ipack( 3)
        NP_G           = ipack( 4)
        NE_G           = ipack( 5)
        hss%lun        = ipack( 6)
        FileFmtVersion = ipack( 7)
        IESTP          = ipack( 8)
        NSCOUE         = ipack( 9)
        IVSTP          = ipack(10)
        NSCOUV         = ipack(11)
        ICSTP          = ipack(12)
        NSCOUC         = ipack(13)
        IPSTP          = ipack(14)
        IWSTP          = ipack(15)
        NSCOUM         = ipack(16)
        IGEP           = ipack(17)
        NSCOUGE        = ipack(18)
        IGVP           = ipack(19)
        NSCOUGV        = ipack(20)
        IGCP           = ipack(21)
        NSCOUGC        = ipack(22)
        IGPP           = ipack(23)
        IGWP           = ipack(24)
        NSCOUGW        = ipack(25)
        iNHSINC        = ipack(26) !tcm v51.26
        
        ttime  = dpack(1)
        HOTSTARTDIR =trim(adjustl(cjunk1))
        hss%filename=trim(adjustl(cjunk2))
        NP = 0
        NE = 0
      endif
      CALL MPI_BARRIER(WCOMM,IERR) !st3

      IF (WRITE_LOCAL_HOT_START_FILES) THEN
         write(6,*) 'HSWriter can not treat Local HOT START FILES'
         STOP
      ENDIF
!
      WRITE(16,*) "***********WRITE HOTSTART*********************"
      WRITE(16,*) "NHSTAR = ", NHSTAR
      WRITE(16,*) "FirstCall = ", FirstCall
      IF (FirstCall) THEN
         IF ( WRITER_ID /= 0  ) THEN
            IF (.not.ALLOCATED(ETA1_g)) THEN
              write(16,*) 'Allocating Elev1Descript'
              ALLOCATE(ETA1_g(NP_G))
            ENDIF
         ENDIF
         Elev1Descript % specifier            = NHSTAR
         Elev1Descript % initial_value        = 0.0
         Elev1Descript % num_items_per_record = 1
         Elev1Descript % num_fd_records       = NP_G
         Elev1Descript % num_records_this     = NP
         Elev1Descript % imap                 => NODES_LG
         Elev1Descript % array                => ETA1
         Elev1Descript % array_g              => ETA1_g

         IF ( WRITER_ID /= 0 ) THEN
            IF (.not.ALLOCATED(ETA2_g)) THEN
               write(16,*) 'Allocating Elev2Descript'
               ALLOCATE(ETA2_g(NP_G))
            ENDIF
         ENDIF
         Elev2Descript % specifier            = NHSTAR
         Elev2Descript % initial_value        = 0.0
         Elev2Descript % num_items_per_record = 1
         Elev2Descript % num_fd_records       = NP_G
         Elev2Descript % num_records_this     = NP
         Elev2Descript % imap                 => NODES_LG
         Elev2Descript % array                => ETA2
         Elev2Descript % array_g              => ETA2_g

         IF ( WRITER_ID /= 0 ) THEN
            IF (.not.ALLOCATED(UU2_g)) THEN
               write(16,*) 'Allocating VelDescript'
               ALLOCATE(UU2_g(NP_G))
               ALLOCATE(VV2_g(NP_G))
            ENDIF
         ENDIF
         VelDescript % specifier            = NHSTAR
         VelDescript % initial_value        = 0.0
         VelDescript % num_items_per_record = 2
         VelDescript % num_fd_records       = NP_G
         VelDescript % num_records_this     = NP
         VelDescript % imap                 => NODES_LG
         VelDescript % array                => UU2
         VelDescript % array2               => VV2
         VelDescript % array_g              => UU2_g
         VelDescript % array2_g             => VV2_g

         IF (IM.eq.10) THEN
            IF ( WRITER_ID /= 0 ) THEN
               write(16,*) 'Allocating CH1Descript'
               ALLOCATE(CH1_g(NP_G))
            ENDIF
            CH1Descript % specifier            = NHSTAR
            CH1Descript % initial_value        = 0.0
            CH1Descript % num_items_per_record = 1
            CH1Descript % num_fd_records       = NP_G
            CH1Descript % num_records_this     = NP
            CH1Descript % imap                 => NODES_LG
            CH1Descript % array                => CH1
            CH1Descript % array_g              => CH1_g
         ENDIF

         IF ( WRITER_ID /= 0 ) THEN
            write(16,*) 'Allocating EtaDiscDescript'
            ALLOCATE(EtaDisc_g(NP_G))
         ENDIF
         EtaDiscDescript % specifier            = NHSTAR
         EtaDiscDescript % initial_value        = 0.0
         EtaDiscDescript % num_items_per_record = 1
         EtaDiscDescript % num_fd_records       = NP_G
         EtaDiscDescript % num_records_this     = NP
         EtaDiscDescript % imap                 => NODES_LG
         EtaDiscDescript % array                => EtaDisc
         EtaDiscDescript % array_g              => EtaDisc_g

         IF ( WRITER_ID /= 0 ) THEN
            write(16,*) 'Allocating NodeCodeDescript'
            ALLOCATE(NodeCode_g(NP_G))
         ENDIF
         NodeCodeDescript % specifier            = NHSTAR
         NodeCodeDescript % int_initial_value        = 0
         NodeCodeDescript % num_items_per_record = 1
         NodeCodeDescript % num_fd_records       = NP_G
         NodeCodeDescript % num_records_this     = NP
         NodeCodeDescript % imap                 => NODES_LG
         NodeCodeDescript % iarray               => NODECODE
         NodeCodeDescript % iarray_g             => NODECODE_g

         IF ( WRITER_ID /= 0 ) THEN
            write(16,*) 'Allocating NOFFDescript'
            ALLOCATE(NOFF_g(NE_G))
         ENDIF
         NOFFDescript % specifier            = NHSTAR
         NOFFDescript % int_initial_value        = 0
         NOFFDescript % num_items_per_record = 1
         NOFFDescript % num_fd_records       = NE_G
         NOFFDescript % num_records_this     = NE
         NOFFDescript % imap                 => IMAP_EL_LG
         NOFFDescript % iarray               => NOFF
         NOFFDescript % iarray_g             => NOFF_g

         FirstCall = .false.
      ENDIF

C     collect up the data from subdomains if running in parallel
         CALL collectFullDomainArrayw(Elev1Descript, packOne, unpackOne)
         CALL collectFullDomainArrayw(Elev2Descript, packOne, unpackOne)
         CALL collectFullDomainArrayw(VelDescript, packTwo, unpackTwo)
         IF (IM.eq.10) THEN
            CALL collectFullDomainArrayw(CH1Descript, packOne, unpackOne)
         ENDIF
         CALL collectFullDomainArrayw(EtaDiscDescript, packOne, unpackOne)
         CALL collectFullDomainIntArrayw(NodeCodeDescript, packOneInt, unpackOneInt)
         CALL collectFullDomainIntArrayw(NOFFDescript, packOneInt, unpackOneInt)

! tcm v49.54 -- changed noff to noff_g to fix an allocation bug in hot start writer
         WHERE (NOFF_G.gt.1) NOFF_g = 1
!
         WRITE(16,*) "FINISH collectFullDomainArray"
      SELECT CASE (NHSTAR)
      !tcm v51.26 added nhstar=-1 for time-stamped hot starts
      CASE(-1,1) ! nonportable binary (1 for backwards compatibility)
         NPX = NP_G
         NEX = NE_G
         IF ((WRITER_ID /= 0)) THEN
            !tcm v51.26 mod for time-stamped hot start files nhstar=-1
            IF (NHSTAR.eq.-1) THEN
               hss % filename(:) = ' '
               hss % filename = 'fort.68_'
               itstr(:) = ' '
               WRITE(itstr,'(I9.9)') IIT  !note that IIT is the correct value and not IT for writers
               hss % filename(9:17) = itstr
            ENDIF
            OPEN(hss % lun ,FILE=TRIM(adjustl(HOTSTARTDIR))//'/'//trim(adjustl(hss%filename)),
     $      ACCESS='DIRECT',RECL=8)
            IHOTSTP=1
            WRITE(hss % lun,REC=IHOTSTP) FileFmtVersion
            IHOTSTP = IHOTSTP + 1
            WRITE(hss % lun,REC=IHOTSTP) IM   ; IHOTSTP = IHOTSTP + 1
            WRITE(hss % lun,REC=IHOTSTP) ttime; IHOTSTP = IHOTSTP + 1 !st3
            WRITE(hss % lun,REC=IHOTSTP) iit  ; IHOTSTP = IHOTSTP + 1 !st3
            WRITE(hss % lun,REC=IHOTSTP) NPX  ; IHOTSTP = IHOTSTP + 1
            WRITE(hss % lun,REC=IHOTSTP) NEX  ; IHOTSTP = IHOTSTP + 1
            WRITE(hss % lun,REC=IHOTSTP) NPX  ; IHOTSTP = IHOTSTP + 1
            WRITE(hss % lun,REC=IHOTSTP) NEX  ; IHOTSTP = IHOTSTP + 1
         ENDIF
         IF ( WRITER_ID /= 0 ) THEN
            DO I=1, Elev1Descript % num_fd_records
               WRITE(hss % lun,REC=IHOTSTP) Elev1Descript % array_g(I)
               IHOTSTP=IHOTSTP+1
            ENDDO
            DO I=1, Elev2Descript % num_fd_records
               WRITE(hss % lun,REC=IHOTSTP) Elev2Descript % array_g(I)
               IHOTSTP=IHOTSTP+1
            ENDDO
            DO I=1, EtaDiscDescript % num_fd_records
               WRITE(hss % lun,REC=IHOTSTP) EtaDiscDescript % array_g(I)
               IHOTSTP=IHOTSTP+1
            ENDDO
            DO I=1, VelDescript % num_fd_records
               WRITE(hss % lun,REC=IHOTSTP) VelDescript % array_g(I)
               IHOTSTP=IHOTSTP+1
            ENDDO
            DO I=1, VelDescript % num_fd_records
               WRITE(hss % lun,REC=IHOTSTP) VelDescript % array2_g(I)
               IHOTSTP=IHOTSTP+1
            ENDDO
            IF (IM.eq.10) THEN
               DO I=1, CH1Descript % num_fd_records
                  WRITE(hss % lun,REC=IHOTSTP) CH1Descript % array_g(I)
                  IHOTSTP=IHOTSTP+1
               ENDDO
            ENDIF
            DO I=1, NodeCodeDescript % num_fd_records
              WRITE(hss % lun,REC=IHOTSTP)NodeCodeDescript % iarray_g(I)
              IHOTSTP=IHOTSTP+1
            ENDDO
            DO I=1, NOFFDescript % num_fd_records
               WRITE(hss % lun,REC=IHOTSTP) NOFFDescript % iarray_g(I)
               IHOTSTP=IHOTSTP+1
            ENDDO
         ENDIF
C
         IF (WRITER_ID /= 0) THEN
            WRITE(hss % lun,REC=IHOTSTP) IESTP   ; IHOTSTP = IHOTSTP + 1
            WRITE(hss % lun,REC=IHOTSTP) NSCOUE  ; IHOTSTP = IHOTSTP + 1

            WRITE(hss % lun,REC=IHOTSTP) IVSTP   ; IHOTSTP = IHOTSTP + 1
            WRITE(hss % lun,REC=IHOTSTP) NSCOUV  ; IHOTSTP = IHOTSTP + 1

            WRITE(hss % lun,REC=IHOTSTP) ICSTP   ; IHOTSTP = IHOTSTP + 1
            WRITE(hss % lun,REC=IHOTSTP) NSCOUC  ; IHOTSTP = IHOTSTP + 1

            WRITE(hss % lun,REC=IHOTSTP) IPSTP   ; IHOTSTP = IHOTSTP + 1
            WRITE(hss % lun,REC=IHOTSTP) IWSTP   ; IHOTSTP = IHOTSTP + 1
            WRITE(hss % lun,REC=IHOTSTP) NSCOUM  ; IHOTSTP = IHOTSTP + 1

            WRITE(hss % lun,REC=IHOTSTP) IGEP    ; IHOTSTP = IHOTSTP + 1
            WRITE(hss % lun,REC=IHOTSTP) NSCOUGE ; IHOTSTP = IHOTSTP + 1

            WRITE(hss % lun,REC=IHOTSTP) IGVP    ; IHOTSTP = IHOTSTP + 1
            WRITE(hss % lun,REC=IHOTSTP) NSCOUGV ; IHOTSTP = IHOTSTP + 1

            WRITE(hss % lun,REC=IHOTSTP) IGCP    ; IHOTSTP = IHOTSTP + 1
            WRITE(hss % lun,REC=IHOTSTP) NSCOUGC ; IHOTSTP = IHOTSTP + 1

            WRITE(hss % lun,REC=IHOTSTP) IGPP    ; IHOTSTP = IHOTSTP + 1
            WRITE(hss % lun,REC=IHOTSTP) IGWP    ; IHOTSTP = IHOTSTP + 1
            WRITE(hss % lun,REC=IHOTSTP) NSCOUGW ; IHOTSTP = IHOTSTP + 1
         ENDIF
      CASE(2) ! ascii?
         write(16,*) 'ASCII Hotstart file not implemented yet.'

      CASE(3) ! netcdf
         write(16,*) 'NetCDF Hotstart undergoing testing.'

      CASE DEFAULT
         write(ScreenUnit,*) 'The value of NHSTAR=',NHSTAR
         write(ScreenUnit,*) 'is not supported.'
         write(ScreenUnit,*) 'Hotstart file not written.'
      END SELECT
C...
C...    CLOSE THE HOT START OUTPUT FILE
C...
      IF(ABS(NHSTAR).EQ.1) THEN                                              ! added by mcf 7/14/08 tcm v51.26 added abs(nhstar)
         if (WRITER_ID /= 0) CLOSE(hss % lun)
      ENDIF                                                             ! end NHSTAR=1

      IF(NSCREEN.NE.0.AND. WRITER_ID /= 0 ) THEN
         WRITE(ScreenUnit,24541) hss % lun, WRITER_ID, IIT, ttime  !tcm v51.26 changed IT to IIT and TimeLoc to ttime since
         WRITE(16,24541) hss % lun, WRITER_ID, IIT, ttime !IT and TimeLoc are not known by writer processor id /=0  (echoing to writer file)
      ELSE
         WRITE(16,24541) hss % lun,WRITER_ID,IT,TimeLoc
      ENDIF
24541 FORMAT(1X,'HOT START OUTPUT WRITTEN TO UNIT ',I2,'  HSWriter:',i2
     &    ' AT TIME STEP = ',I9,' TIME = ',E15.8)
      IF(hss % lun.EQ.67) THEN
C        jgf45.07 added option to stop ADCIRC after writing hot start file.
         IF (NHSTAR.EQ.67) THEN
            WRITE(ScreenUnit,*) 'NHSTAR.EQ.67; ADCIRC stopping.'
            WRITE(16,*) 'NHSTAR.EQ.67; ADCIRC stopping.'
         ENDIF
         hss % lun      = 68
         IF(NHSTAR.EQ.1) THEN
            hss % filename(:) = ' '
            hss % filename = 'fort.68'                     ! added by mcf 7/14/08
         ENDIF
         !tcm v51.26 mod for time-stamped nhstar=-1  This section is just a failsafe
         IF (NHSTAR.eq.-1) THEN
            hss % filename(:) = ' '
            hss % filename = 'fort.68_'
            itstr(:) = ' '
            WRITE(itstr,'(I9.9)') IIT+iNHSINC  !setting the name for the next hot start file
            hss % filename(9:17) = itstr
         ENDIF         
         
      ELSE
         IF (NHSTAR.EQ.68) THEN
            WRITE(ScreenUnit,*) 'NHSTAR.EQ.68; ADCIRC stopping.'
            WRITE(16,*) 'NHSTAR.EQ.68; ADCIRC stopping.'
         ENDIF
         hss % lun=67
         IF(NHSTAR.EQ.1) THEN
            hss % filename(:) = ' '
            hss % filename = 'fort.67'                     ! added by mcf 7/14/08
         ENDIF
         !tcm v51.26 mod for time-stamped nhstar=-1  This section is just a failsafe
         IF (NHSTAR.eq.-1) THEN
            hss % filename(:) = ' '
            hss % filename = 'fort.68_'
            itstr(:) = ' '
            WRITE(itstr,'(I9.9)') IIT+iNHSINC  !setting the name for the next hot start file
            hss % filename(9:17) = itstr
         ENDIF         
         
      ENDIF
      
      call unsetMessageSource()
      
      END SUBROUTINE writeHotstart_through_HSwriter
!
      SUBROUTINE GET_NEXT_HSWRITER_COMM(NEXT_WCOMM)
      include 'mpif.h'
      INTEGER,intent(inout) :: NEXT_WCOMM
      INTEGER :: NEXT_WCOMMID

        IF(WCOMMID >= MNWPROH) THEN
          WCOMMID = 0
        ENDIF

        WCOMMID = WCOMMID + 1
        NEXT_WCOMMID = MOD(WCOMMID - 1,MNWPROH) + 1
        NEXT_WCOMM = COMM_WRITEH(NEXT_WCOMMID)

      END SUBROUTINE GET_NEXT_HSWRITER_COMM
!
!     !--------------------------------------------------------------
      subroutine collectFullDomainArrayw(descript, pack_cmd, unpack_cmd)
!     !--------------------------------------------------------------
      implicit none
#ifdef CMPI
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
#endif
      type (OutputDataDescript_t) :: descript
      external pack_cmd
      external unpack_cmd
#ifdef CMPI
!     ! the subroutine used to write the file
      integer      :: ierr, status(MPI_STATUS_SIZE), request
      integer, save:: tagbase = 6000
      integer      :: iproc
      integer      :: bufsize
      integer      :: ibucket
      integer      :: istart     ! vector tuple to start with
      integer      :: iend       ! vector tuple to end on
      integer      :: tag
!     ! number of vector tuples in the buffer
      integer      :: num
      integer      :: i, j, k

      bufsize = min(BUFSIZE_MAX,
     &    descript % num_items_per_record * descript % num_fd_records)
      num     = bufsize / descript % num_items_per_record
      iend    = num
      istart  = 1

      if (tagbase == 5000) then
         tagbase = 6000
      else
         tagbase = 5000
      endif
      ibucket = 0

      do while (istart < iend)

!        !------------------------------------------------------------
!        ! Initialize
!        !------------------------------------------------------------
         buf(:)  = descript % initial_value
         ibucket = ibucket + 1
         tag     = tagbase + mod(ibucket, 8)

         call pack_cmd(descript, istart, iend)
         call mpi_reduce(buf, resultBuf, bufsize, float_type, MPI_SUM,
     &                   MNPROC, WCOMM, ierr)

         if (WRITER_ID /= 0) then
            call unpack_cmd(descript, istart, iend)
         end if
         istart = iend + 1
         iend   = min(istart + num - 1, descript % num_fd_records)
         num    = iend - istart + 1
      end do
#endif
!     !--------------------------------------------------------------
      end subroutine collectFullDomainArrayw
!     !--------------------------------------------------------------
!
!     !--------------------------------------------------------------
      subroutine collectFullDomainIntArrayw(descript,
     &                                     pack_cmd, unpack_cmd)
!     !--------------------------------------------------------------
      implicit none
#ifdef CMPI
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
#endif
      type (OutputDataDescript_t) :: descript
      external pack_cmd
      external unpack_cmd
#ifdef CMPI
!     ! the subroutine used to write the file
      integer      :: ierr, status(MPI_STATUS_SIZE), request
      integer, save:: tagbase = 6000
      integer      :: iproc
      integer      :: bufsize
      integer      :: ibucket
      integer      :: istart     ! vector tuple to start with
      integer      :: iend       ! vector tuple to end on
      integer      :: tag
!     ! number of vector tuples in the buffer
      integer      :: num
      integer      :: i, j, k

      bufsize = min(BUFSIZE_MAX,
     &    descript % num_items_per_record * descript % num_fd_records)
      num     = bufsize / descript % num_items_per_record
      iend    = num
      istart  = 1
      if (tagbase == 5000) then
         tagbase = 6000
      else
         tagbase = 5000
      endif
      ibucket = 0

      do while (istart < iend)
!      ! Initialize
         integerBuffer(:)  = descript % int_initial_value
         ibucket = ibucket + 1
         tag     = tagbase + mod(ibucket, 8)

         call pack_cmd(descript, istart, iend)
         call mpi_reduce(integerBuffer, integerResultBuffer, bufsize,
     &                 MPI_INTEGER, MPI_SUM, MNPROC, WCOMM, ierr)

         if (WRITER_ID /= 0) then
            call unpack_cmd(descript, istart, iend)
         end if
         istart = iend + 1
         iend   = min(istart + num - 1, descript % num_fd_records)
         num    = iend - istart + 1
      end do
#endif
C! CMPI
!     !--------------------------------------------------------------
      end subroutine collectFullDomainIntArrayw
!     !--------------------------------------------------------------
!
      END MODULE HSWRITER
