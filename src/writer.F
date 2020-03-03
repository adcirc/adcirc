C******************************************************************************
C MODULE WRITER
C  This module is for writer processors which receives values from compute
C  processors and writes them out asynchronously.
C                                                              05 Feb 2007, sb
C******************************************************************************
C
      module writer
      use CkLunMod, only : CK_LUN
      use sizes, only : sz, mnproc, mnwproc, OFF, ASCII, SPARSE_ASCII,
     &     NETCDF3, NETCDF4, XDMF, iSplit, globalDir, meshType, 
     &     controlType, meshFileName_g, controlFileName_g
      use global, only: c3d, comm,
     &     bufsize_max, buf, resultbuf, outputdatadescript_t, float_type,
     &     debug, echo, info, warning, error, logmessage, allmessage,
     &     setmessagesource, unsetmessagesource, initLogging, nabout,
     &     scratchMessage, 
     &     sig_val, sig_write, sig_pause, comm_writer, writer_id,
     &     sig_term, realtype, dbletype, sig_mesh
      use messenger, only: tag, msg_fini
#ifdef HAVE_MPI_MOD
      use mpi
      implicit none
#else
      implicit none
#endif
      ! tcm v50.27 increased size from 3 to 4, cobell increase to 6,
      ! jgf51.21.27 decreased size of intbuf to 5
      ! jgf51.21.29 increased size of intbuf to 6 for XDMF
      ! jgf52.08.25: increased intbuf to 7 to transmit nabout
      integer, parameter :: numIntMetadata = 7
      integer :: intBuf(numIntMetadata)  ! used to transmit pieces of descript via mpi
!$omp threadprivate(intBuf)
      integer, parameter :: numRealMetadata = 2
      real(sz) :: reBuf(numRealMetadata) ! used to transmit pieces of descript via mpi
!$omp threadprivate(reBuf)
      !
      ! maximum number of buffers that each writer processor has
      ! for storing data; each writer will only allocate the space that
      ! it actually needs, so a large value (e.g., 25) does not consume
      ! any extra or unnecessary resources
      integer, parameter :: num_buf_max = 25  !st3 100708: increase for swan output
      !
      ! data structure on the writer processor that mirrors the corresponding
      ! metadata on the compute processor; this is dimensioned by the 
      ! num_buf_max parameter
      type(OutputDataDescript_t), allocatable :: writer_descripts(:)
!$omp threadprivate(writer_descripts)
      ! 
      ! .true. if a writer processor has received data into one or more
      ! of its local buffers but has not written it yet; the nwloaded
      ! array is dimensioned equal to the number of writer processors;
      ! this is an array that is updated and used by compute processors
      ! not writer processors
      logical, allocatable :: nwloaded(:)   !st3 100708: round robin chamber
!$omp threadprivate(nwloaded)
      !
      ! buffer number on a particular writer processor; each writer 
      ! processor manages its own bufID so it knows where to store the
      ! next dataset it receives
      integer :: bufID    
!$omp threadprivate(bufID)
      !
      ! wcommID is incremented by and used only by compute processors to
      ! determine which writer processor to send data to each time
      integer :: wcommID  ! writer processor number (1 -> MNWPROC)
!$omp threadprivate(wcommID)

C---------------------end of data declarations--------------------------------C

      !--------    
      CONTAINS
      !--------

C-----------------------------------------------------------------------
C     S U B R O U T I N E   W R I T E R  _  I N I T
C-----------------------------------------------------------------------
C     jgf51.21.27: This subroutine is executed by each writer processor.
C 
C     Initializes round robin chamber; clears mpi communicator id for
C     this processor; clears flush flag.
C-----------------------------------------------------------------------
      subroutine writer_init ()
      implicit none
      integer i

      call setMessageSource("writer_init")
#if defined(WRITER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      wcommID = 0
      bufID = 0
      allocate (nwloaded(1:mnwproc)) !st3 100708: allocate round robin chamber
      nwloaded(1:mnwproc) = .false.    !st3 100708: zeroc    round robin chamber
      allocate(writer_descripts(num_buf_max))
      do i=1,num_buf_max
         writer_descripts(i) % initialized = .false.
         writer_descripts(i) % num_fd_records = 0
         writer_descripts(i) % num_items_per_record = 0
         writer_descripts(i) % specifier = 0
         writer_descripts(i) % lun = 0
         writer_descripts(i) % outputTimeStepIncrement = 0
         writer_descripts(i) % alternate_value = -99999.d0
         writer_descripts(i) % wMD % bufferInitialized = .false.
         writer_descripts(i) % wMD % bufferLoaded = .false.
         writer_descripts(i) % wMD % myTimeStep = 0
         writer_descripts(i) % wMD % myTime = 0.d0
      end do
      
#if defined(WRITER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN
C-----------------------------------------------------------------------
      END SUBROUTINE WRITER_INIT
C-----------------------------------------------------------------------
      

C-----------------------------------------------------------------------
C     S U B R O U T I N E   S E N D L A B E L I N F O T O W R I T E R S
C-----------------------------------------------------------------------
C     This subroutine sends the label information from the hashtable
C     to the writer processors so they can write the output data
C     correctly
C-----------------------------------------------------------------------
      SUBROUTINE sendLabelInfoToWriters()
C-----------------------------------------------------------------------
          use mpi
          use sizes,only: myproc,mnwproc,mnproc
          use messenger,only: mpi_comm_adcirc
          use global,only: labels_g,np_g
          implicit none
          INTEGER :: ierr
          INTEGER :: I,N
          INTEGER :: mpistat(mpi_status_size)

          IF(myproc.EQ.0)THEN
              DO I = 0,MNWPROC-1
                 CALL mpi_send(np_g,1,mpi_integer,mnproc+i,0,mpi_comm_adcirc,ierr)
                 CALL mpi_send(labels_g,np_g,mpi_integer,mnproc+i,0,mpi_comm_adcirc,ierr)
              ENDDO
          ELSEIF(myproc.GT.mnproc-1.and.myproc.LT.mnproc+mnwproc)THEN
              CALL mpi_recv(n,1,mpi_integer,0,0,mpi_comm_adcirc,mpistat,ierr)
              allocate(labels_g(1:n))
              CALL mpi_recv(labels_g,n,mpi_integer,0,0,mpi_comm_adcirc,mpistat,ierr)
          ENDIF

          RETURN

C-----------------------------------------------------------------------
      END SUBROUTINE sendLabelInfoToWriters
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     S U B R O U T I N E   W R I T E R  _  M A I N
C-----------------------------------------------------------------------
C     jgf51.21.27: This subroutine is executed by each writer processor.
C
C     Writer processor goes into an infinite loop where it waits for
C     either a signal to receive data, or a signal to write output.
C     When it receives a signal, it performs the specified action,
C     then goes back to waiting for another signal.
C-----------------------------------------------------------------------
      subroutine writer_main ()
      USE GLOBAL,ONLY : ScreenUnit
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      integer :: sig,ierr
      integer :: mpistat(mpi_status_size)
C
      call setMessageSource("writer_main")
#if defined(WRITER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      !
      ! infinite loop
      do 
         call mpi_recv(sig, 1, mpi_integer, 0, tag,
     &       comm_writer(writer_id), mpistat, ierr)

         select case(sig)
         case(sig_val) 
            write(scratchMessage,
     &       '("Writer processor ",i0," received sig_val.")') writer_id
            call logMessage(INFO,scratchMessage)
            call writer_recv_values()
         case(sig_write)
            write(scratchMessage,
     &       '("Writer processor ",i0," received sig_write.")') writer_id
            call logMessage(INFO,scratchMessage)
            call writer_write()
         case(sig_mesh)
            write(scratchMessage,
     &       '("Writer processor ",i0," received sig_mesh.")') writer_id
            call logMessage(INFO,scratchMessage)
            call writerReadMesh()   
         case(sig_pause)   ! tcm addeded v49.52.01 SIG_PAUSE
            write(scratchMessage,
     &       '("Writer processor ",i0," received sig_pause.")') writer_id
            call logMessage(INFO,scratchMessage)
            exit
         case(sig_term) 
            write(scratchMessage,
     &       '("Writer processor ",i0," received sig_term.")') writer_id
            call logMessage(INFO,scratchMessage)
            exit
         case default
            write(scratchMessage,
     &     '("Writer processor ",i0," received invalid signal ",i0,".")') 
     &       writer_id, sig
            exit
         end select
      end do

#if defined(WRITER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      end subroutine writer_main
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E   S E N D   D A T A   T O   W R I T E R
C-----------------------------------------------------------------------
C     jgf51.21.27: This subroutine is executed by all compute processors
C     to send metadata and then output data to a writer processor. 
C-----------------------------------------------------------------------
      subroutine sendDataToWriter(descript, TimeLoc, it, store_cmd)
      use sizes, only : myproc
      use global      
      implicit none
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      type (OutputDataDescript_t), intent(in) :: descript
      real(8) :: TimeLoc  ! output time in seconds
      integer :: it      ! current time step
      external store_cmd ! command to pack the buffer for mpi transmission
C   
      integer :: ierr
      integer :: num
      integer :: bufsize
      integer :: istart, iend
C
      call setMessageSource("sendDataToWriter")
#if defined(WRITER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      !
      ! jgf51.21.27: Figure out which mpi communicator to use with this
      ! writer processor.
      if (wcommID.ge.mnwproc) then
         wcommID = 0
      endif
      wcommID = wcommID + 1  
      ! st3 100708: load the bullet to round robin chamber
      nwloaded(wcommID) = .true.  
      !
      ! jgf51.21.27: proc 0 signals the writer processor to let it know 
      ! that it is time to write and then sends it some metadata
      if (myproc.eq.0) then
         write(scratchMessage,'(a,i0,a)') 
     &    'Processor 0 is sending sig_val to writer ',wcommID,'.'
         call logMessage(DEBUG,scratchMessage)
         !
         ! Send SIG_VAL to a writer to get it ready to receive values.
         call mpi_send(SIG_VAL, 1, mpi_integer, mnproc, tag,
     &      comm_writer(wcommID), ierr)

C...     Send minimal set of metadata associated with this file
C....    tcm v50.27 added for full/compact record
Cobell - Allow writer processors to output netcdf data
         ! jgf51.21.29: Added outputTimeStepIncrement since it is needed
         ! by XDMF so it will know when to finalize the output
         INTBUF(1) = descript % num_fd_records
         INTBUF(2) = descript % num_items_per_record
         INTBUF(3) = descript % specifier
         INTBUF(4) = descript % lun
         INTBUF(5) = descript % outputTimeStepIncrement
         INTBUF(6) = it
         INTBUF(7) = nabout ! jgf52.08.25: need logging level
         call mpi_send(intbuf, numIntMetaData, mpi_integer, mnproc, tag, 
     &      comm_writer(wcommID),ierr)
         REBUF(1) = descript % alternate_value
         REBUF(2) = timeLoc
         ! send dry value and time in seconds
         call mpi_send(rebuf, numRealMetaData, realtype, mnproc, tag, 
     &      comm_writer(wcommID),ierr)
         ! send file name
         call mpi_send(descript % file_name, 1024, mpi_character, 
     &      mnproc, tag, comm_writer(wcommID), ierr)
         write(scratchMessage,'(a,i0,a,a,a,i0,a)') 
     &     'Processor ',myproc,' is sending file name ',
     &     trim(descript%file_name),' to writer ',wcommID,'.'
         call logMessage(DEBUG,scratchMessage)
      endif

C...  Compute the buffer size
      bufsize = min(BUFSIZE_MAX,
     &    descript % num_items_per_record * descript % num_fd_records)
      num     = bufsize / descript % num_items_per_record
      iend    = num
      istart  = 1
      ! 
      ! now all compute processors send this data to the dedicated writer processor
      do while (istart.le.iend)

         buf(:)  = descript % initial_value
         call store_cmd(descript, istart, iend)

         write(scratchMessage,'(a,i0,a,a,a)') 'proc ',myproc,
     &    ' is sending for ',trim(descript % file_name),'.'
         call logMessage(DEBUG,scratchMessage)

C...     The following mpi_reduce sends array values to rank MNPROC,
C...     which is the writer proc in the group.
C...     Note that ranks from 0 to (MNPROC-1) are compute procs, and
C...     rank MNPROC of writer communicators is always the writer proc.
         call mpi_reduce(buf, resultBuf, bufsize, float_type, MPI_SUM,
     &      MNPROC, comm_writer(wcommID), ierr)
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
      end subroutine sendDataToWriter
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     S U B R O U T I N E     W R I T E R _ R E C V _ V A L U E S
C-----------------------------------------------------------------------
C     jgf51.21.27: This subroutine is executed by a writer processor
C     that has received a sig_val from processor 0 indicating that data
C     are about to be sent. The writer processor receives some metadata
C     before receiving the actual output data. The data are not actually
C     written in this subroutine.
C-----------------------------------------------------------------------
      subroutine writer_recv_values()
      implicit none
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      integer :: num, i, j
      integer :: ierr
      integer :: bufsize
      integer :: istart, iend
      integer :: mpistat(MPI_STATUS_SIZE)
      logical :: descriptFound ! true if we found the particular buffer we need
      character(len=8) :: cdigit  !st3 100708: split file

      call setMessageSource("writer_recv_values")
#if defined(WRITER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      !
      ! Receive minimal set of metadata required to write the file
      call mpi_recv(intbuf, numIntMetaData, mpi_integer, 0, tag, 
     &   comm_writer(writer_id), mpistat, ierr)
      !
      ! Select the writer buffer that will receive the values. 
      descriptFound = .false.
      do i=1,num_buf_max
         ! look for the descript that is already set up to write this data
         if (writer_descripts(i) % lun .eq. intbuf(4)) then
            bufID = i
            descriptFound = .true.
            exit
         endif
      end do
      ! If this writer has not written this type of output file
      ! yet, then find a buffer that has not been used yet.
      ! NetCDF and XDMF data files get initialized only once on each writer, 
      ! so on each writer, we'll reuse the same buffer each time we 
      ! write to these files. For ascii files it doesn't actually matter
      ! whether we use the same or different buffers each time, but 
      ! treating all the output types the same keeps things simple.
      if (descriptFound.eqv..false.) then
         do i=1,num_buf_max
            if ( writer_descripts(i)%lun.eq.0 ) then
               bufID = i
               descriptFound = .true.
               exit
            endif
          end do
      endif
      ! 
      ! place the metadata into the selected buffer           
      writer_descripts(bufid) % num_fd_records = intbuf(1)
      writer_descripts(bufid) % num_items_per_record = intbuf(2)
      writer_descripts(bufid) % specifier = intbuf(3) !jgf51.21.27 added
      writer_descripts(bufid) % lun = intbuf(4)
      writer_descripts(bufid) % outputTimeStepIncrement = intbuf(5)
      writer_descripts(bufid) % wMD % myTimeStep = intbuf(6)
      nabout = intbuf(7)
      !
      ! receive "dry" data value and current simulation time
      call mpi_recv(rebuf, numRealMetaData, realtype, 0, tag, 
     &   comm_writer(writer_id), mpistat, ierr)
      writer_descripts(bufid) % alternate_value = rebuf(1)
      writer_descripts(bufid) % wMD % myTime = rebuf(2)
      ! 
      ! receive file name
      writer_descripts(bufID) % file_name(:) = " "
      call mpi_recv(writer_descripts(bufID)%file_name, 1024, 
     & mpi_character, 0, tag, comm_writer(writer_id), mpistat, ierr)
      write(scratchMessage,'(a,i0,a,i0,a,a,a)') 'Writer ',writer_id,
     &   ' buffer ', bufID,' received file name ',
     &   trim(writer_descripts(bufID)%file_name)
      call logMessage(DEBUG,scratchMessage)
      !
      !st3 100708: split file
      if (isplit.eqv..true.) then                             
         !st3 split files are named by integer part of simulation time
         write(cdigit,'(i8.8)') 
     &      nint(writer_descripts(bufid) % wMD % myTime)        
         writer_descripts(bufid) % file_name = 
     &      trim(adjustl(writer_descripts(bufid) % file_name))
     &      //'.'//cdigit 
      endif                                                     
      !
      ! Prepare buffer to store received values if this is the first time
      ! for this buffer on this writer
      if ( writer_descripts(bufid) % wMD % bufferInitialized.eqv..false.) then
         allocate(writer_descripts(bufid) 
     &      % array_g(writer_descripts(bufID) % num_fd_records) )
         ! also allocate the 2nd component, even though we may not 
         ! need it for storing this particular data set
         allocate(writer_descripts(bufID) 
     &      % array2_g(writer_descripts(bufID) % num_fd_records) )
         ! don't need to allocate memory for this buffer on this
         ! writer again
         writer_descripts(bufid) % wMD % bufferInitialized = .true.
      endif

C...  Compute the writer mpi buffer size 
      bufsize = min(BUFSIZE_MAX, 
     &    writer_descripts(bufid) % num_items_per_record
     &   * writer_descripts(bufid) % num_fd_records)
      num     = bufsize / writer_descripts(bufid) % num_items_per_record
      iend    = num
      istart  = 1
      !
      ! writer processor now actually receives the data
      do while (istart.le.iend)

         buf(:)  = 0.D0

C...     The following mpi_reduce sends array values to rank MNPROC,
C...     which is a writer proc.
C...     Note that ranks from 0 to (MNPROC-1) are compute procs, and
C...     rank MNPROC of writer communicators is always a writer proc.
         call mpi_reduce(buf, resultBuf, bufsize, float_type, MPI_SUM,
     &     MNPROC, comm_writer(writer_id), ierr)

         j = 1
         do i = istart, iend
            writer_descripts(bufid) % array_g(i) = resultBuf(j)
            j = j+1

            if (writer_descripts(bufid) % num_items_per_record.EQ.2) then
               writer_descripts(bufid) % array2_g(i) = resultBuf(j)
               j = j+1
            endif
         end do

         istart = iend + 1
         iend   = min(istart + num - 1, 
     &                writer_descripts(bufid) % num_fd_records )
         num    = iend - istart + 1
      end do
      writer_descripts(bufid) % wMD % bufferLoaded = .true.

#if defined(WRITER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      return
C-----------------------------------------------------------------------
      end subroutine writer_recv_values
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     S U B R O U T I N E     F L U S H _ W R I T E R S
C-----------------------------------------------------------------------
C     jgf51.21.27: Processor 0 will now check to see which writer 
C     processors have received data but have not written it; 
C     send signals to these writer processors and make each one
C     actually write the data for all the buffers it may be
C     holding. 
C-----------------------------------------------------------------------
      SUBROUTINE FLUSH_WRITERS()
      use sizes, only : myproc
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
      !
      ! return immediately if no writer processor has received data 
      ! since the last flush
      if ( .not.any(nwloaded).eqv..true. ) then
#if defined(WRITER_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG,"Return.")
#endif
         call unsetMessageSource()
         return
      endif
      !
      ! processor 0 loops over the writer processors and tells the ones
      ! that have received data in one or more of their buffers that they
      ! need to write it all now
      if(myproc.eq.0) then
         do i=1,mnwproc
            !st3 100708: flash only loaded chamber
            if (nwloaded(i).eqv..true.) then   !st3 nwloaded()=false:unloaded, true:loaded
               write(scratchMessage,'(a,i0,a,i0,a)') 'Proc ',myproc,
     &           ' is sending sig_write to writer ',i,'.'
               call logMessage(DEBUG, scratchMessage)
               call mpi_send(sig_write,1,mpi_integer,mnproc,
     &              tag,comm_writer(i),ierr)
               nwloaded(i) = .false.
            endif
         enddo
      endif

#if defined(WRITER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      return
C-----------------------------------------------------------------------
      end subroutine flush_writers
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     S U B R O U T I N E     W R I T E R _ W R I T E
C-----------------------------------------------------------------------
C     Called when a dedicated writer processor receives a call to write
C     out the data in its buffers. 
C
C     jgf51.21.27: Modified this to use the original 
C     descript % specifier integer instead of FULL_REC and W_NETCDF 
C-----------------------------------------------------------------------
      subroutine writer_write()
      use global
#ifdef ADCNETCDF
      use NETCDFIO, only : initNetCDFOutputFile, writeOutArrayNetCDF
#endif      
#ifdef ADCXDMF
      use XDMFIO, only : initOutputXDMF, writeOutArrayXDMF, writeControlXDMF
#endif
      implicit none
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      integer :: ierr
      integer :: ibuf,i,j    ! loop counters
      integer :: nnondefval  ! non default values for sparse output
      !
      ! the following logicals are dimensioned by the number of nodes (np_g)
      logical, allocatable :: areDefault1(:) ! true where array_g equals alternate_value
      logical, allocatable :: areDefault2(:) ! true where array2_g equals alternate_value
      integer, parameter :: fid = 100
      character(len=20) :: fpos  ! file position specifier on open statement
      character(len=20) :: fstat ! file status specifier on open statement
      logical :: nerr             ! netcdf error status indicator

      call setMessageSource("writer_write")
#if defined(WRITER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      write(scratchMessage,'(a,i0,a)') 
     &    'Writer ',writer_id,' is writing...'
      call logMessage(INFO, scratchMessage)      
      !
      ! Loop over the buffers on this writer processor and write out
      ! the data that are stored in them 
      do ibuf=1, num_buf_max
         if (writer_descripts(ibuf) % wMD % bufferLoaded .eqv. .false.) then
            !write(6,'(a,i0,a,i0,a)') 'DEBUG: Writer ',writer_id,' buffer ',ibuf,' is empty.'
            cycle
         endif
         write(scratchMessage,'(a,i0,a,i0,a,a,a)') 'Writer ',writer_id,
     &       ' buffer ',ibuf,
     &       ' is full and will be written to the file "',
     &        trim(writer_descripts(ibuf) % file_name),'"'
         call logMessage(DEBUG, scratchMessage)
         !
         ! for ascii and sparse ascii, open the file appropriately
         if ( (abs(writer_descripts(ibuf) % specifier).eq.ASCII).or.
     &        (abs(writer_descripts(ibuf) % specifier).eq.SPARSE_ASCII) 
     &      ) then
            fstat = 'old'
            fpos = 'append'
            if (isplit.eqv..true.) then  !st3  split file 100708
               fstat = 'replace'
               fpos = 'asis'
            endif
            open(fid+CK_LUN, file=writer_descripts(ibuf) % file_name, 
     &             status=fstat, position=fpos, action='write')
         endif
         ! 
         select case(abs(writer_descripts(ibuf) % specifier))
         !
         !    A S C I I
         !
         ! tcm v50.27 added full record writing options
         case(ASCII)
            ! write dataset header
            write(fid+CK_LUN, 100) writer_descripts(ibuf) % wMD % myTime,
     &            writer_descripts(ibuf) % wMD % myTimeStep 
            ! write dataset
            select case(writer_descripts(ibuf)%num_items_per_record)
            case(1)
               do i=1, writer_descripts(ibuf) % num_fd_records
                  write(fid+CK_LUN,1000) labels_g(i), writer_descripts(ibuf) % array_g(i)
               end do
            case(2)
               do i=1, writer_descripts(ibuf) % num_fd_records
                  write(fid+CK_LUN,2000) labels_g(i), writer_descripts(ibuf) % array_g(i),
     &               writer_descripts(ibuf) % array2_g(i)
               end do
            case default
               write(scratchMessage,'(a,i0,a)') 'Writer cannot write ',
     &            writer_descripts(ibuf) % num_items_per_record,
     &            ' component data.'
               call allMessage(ERROR,scratchMessage)
            end select
            close(fid+CK_LUN)
         !
         !  S P A R S E   A S C I I 
         !
         case(SPARSE_ASCII)
            ! Count non-default values
            allocate(areDefault1(writer_descripts(ibuf) 
     &                                             % num_fd_records))
            areDefault1 = .true.
            where ( writer_descripts(ibuf) % array_g .ne.
     &                  writer_descripts(ibuf) % alternate_value ) 
               areDefault1 = .false.
            end where          
            !
            ! write data
            select case(writer_descripts(ibuf)%num_items_per_record)
            !
            !.... Sparse ASCII Scalar data
            case(1) 
               nNonDefVal = count(areDefault1.eqv..false.)
               ! write dataset header
               write(fid+CK_LUN, 100) writer_descripts(ibuf) % wMD % myTime,
     &               writer_descripts(ibuf) % wMD % myTimeStep, 
     &               nNonDefVal, 
     &               writer_descripts(ibuf) % alternate_value
               ! write sparse scalar dataset
               do i=1, writer_descripts(ibuf) % num_fd_records
                  if (areDefault1(i).eqv..false.) then
                     write(fid+CK_LUN,1000) labels_g(i), 
     &                     writer_descripts(ibuf) % array_g(i)
                  endif
               end do
               close(fid+CK_LUN)
               deallocate(areDefault1)
            !
            !.... Sparse ASCII Vector data
            case(2) 
               ! now we have to count the non default values on
               ! the second component as well
               allocate(areDefault2(writer_descripts(ibuf) 
     &                                             % num_fd_records))
               areDefault2 = .true.
               where ( writer_descripts(ibuf) % array2_g .ne.
     &                  writer_descripts(ibuf) % alternate_value ) 
                  areDefault2 = .false. 
               end where
               nNonDefVal = count((areDefault1.eqv..false.)
     &                          .and.(areDefault2.eqv..false.))
               ! write dataset header
               write(fid+CK_LUN, 100) writer_descripts(ibuf) % wMD % myTime,
     &               writer_descripts(ibuf) % wMD % myTimeStep, 
     &               nNonDefVal, 
     &               writer_descripts(ibuf) % alternate_value
               ! write sparse vector dataset
               do i=1, writer_descripts(ibuf) % num_fd_records
                  if ((areDefault1(i).eqv..false.).or.
     &                   (areDefault2(i).eqv..false.)) then
                     write(fid+CK_LUN,2000) labels_g(i), 
     &                     writer_descripts(ibuf) % array_g(i),
     &                     writer_descripts(ibuf) % array2_g(i)
                  endif
               end do
               close(fid+CK_LUN)
               deallocate(areDefault1, areDefault2)
            case default
               write(scratchMessage,'(a,i0,a)') 'Writer cannot write ',
     &               writer_descripts(ibuf) % num_items_per_record,
     &               ' component data.'
               call allMessage(ERROR, scratchMessage)
            end select
            close(fid+CK_LUN)
         !
         !     N E T C D F
         ! 
         case(NETCDF3, NETCDF4)
Cobell ... Added to allow writer core to write NetCDF output without
C          needing access to the entire descript array which would have
C          required extra message passing steps.
#ifdef ADCNETCDF         
            if (writer_descripts(ibuf)%initialized.eqv..false.) then
               call initNetCDFOutputFile(writer_descripts(ibuf), nerr)
               writer_descripts(ibuf)%initialized = .true.
            endif   
            call writeOutArrayNetCDF(writer_descripts(ibuf) % lun,
     &         writer_descripts(ibuf) % wMD % myTime,
     &         writer_descripts(ibuf) % wMD % myTimeStep,
     &         writer_descripts(ibuf) ) 
#else
            write(scratchMessage,'(a)') 
     &          'NetCDF is not available. No output has been written.'
            call allMessage(ERROR, scratchMessage)
#endif         
         !
         !    X D M F
         ! 
         case(XDMF)
#ifdef ADCXDMF
            if (writer_descripts(ibuf)%initialized.eqv..false.) then
               write(scratchMessage,'(a,i0,a,a)') ' Writer ',writer_id,
     &          ' is calling initOutputXDMF for file ',
     &          trim(writer_descripts(ibuf) % file_name)
               call logMessage(DEBUG, scratchMessage)
               call initOutputXDMF(writer_descripts(ibuf))
               call writeControlXDMF(
     &             writer_descripts(ibuf) % xdmfMD % xdmfFortranObj)
               writer_descripts(ibuf)%initialized = .true.
            endif   
            write(scratchMessage,'(a,i0,a,a)') 'Writer ',writer_id,
     &          ' is calling writeOutArrayXDMF for file ',
     &          trim(writer_descripts(ibuf) % file_name)
            call logMessage(DEBUG, scratchMessage)
            call writeOutArrayXDMF(
     &         writer_descripts(ibuf) % wMD % myTime,
     &         writer_descripts(ibuf) % wMD % myTimeStep,
     &         writer_descripts(ibuf) )

#else
            call allMessage(ERROR,
     &         'XDMF support was not compiled into this executable.')
#endif


         case default
            write(scratchMessage,'(a,i0,a)') 
     &         'Output specifier ',
     &         writer_descripts(ibuf) % specifier,
     &         ' is not supported by dedicated writer processors.'   
            call allMessage(ERROR,scratchMessage)
         end select
         writer_descripts(ibuf) % wMD % bufferLoaded = .false.
      end do
      ! set the buffer counter back to zero on this writer processor
      bufid = 0

      ! changed formats to match those in globalio.F
 100  FORMAT(2x,1pE20.10E3,5X,I10,5X,I10,5X,1pE20.10E3)  
      ! tcm v50.27 added for full record
 110  FORMAT(2x,1pE20.10E3,5X,I10)             

      write(scratchMessage,'(a,i0,a)') 'Writer ',writer_id,
     &    ' finished writing.' 
      call allMessage(INFO, scratchMessage)

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
C     S U B R O U T I N E   S E N D   I N I T   W R I T E R   X D M F  
C-----------------------------------------------------------------------
C     jgf51.21.27: This subroutine is executed by processor 0 to send
C     mesh file info to the writer processor so it can read the mesh
C     and use it during the writing of XDMF output.  
C-----------------------------------------------------------------------
      subroutine sendInitWriterXDMF(descript)
      use sizes, only : myproc
      use global      
      implicit none
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      type (OutputDataDescript_t), intent(in) :: descript
      character(len=80) :: fileKeyword
      character(len=80) :: fileValue
      integer :: ierr
C   
      call setMessageSource("sendInitWriterXDMF")
#if defined(WRITER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      ! jgf51.21.27: Figure out which mpi communicator to use with this
      ! writer processor.
      if (wcommID.ge.mnwproc) then
         wcommID = 0
      endif
      wcommID = wcommID + 1  
      !
      ! jgf51.21.27: proc 0 signals the writer processor to tell it that
      ! the mesh info is coming, and then sends the data
      if (myproc.eq.0) then
         write(scratchMessage,'(a,i0,a,i0,a)') 
     &  'Processor ',myproc,' is sending sig_mesh to writer ',wcommID,'.'
         call logMessage(INFO, scratchMessage)
         !
         ! Send sig_mesh to a writer to get it ready to receive metadata.
         call mpi_send(sig_mesh, 1, mpi_integer, mnproc, tag,
     &      comm_writer(wcommID), ierr)
         ! 
         ! read the fulldomain mesh file name and control file name
         ! so that the writer can load them up and use them (it only
         ! uses the control file for generating metadata) ... the 
         ! 'fulldomainInputFiles' file is created by adcprep
         call openFileForRead(30,trim(globaldir)//'/fulldomainInputFiles', ierr)
         do 
            read(30+CK_LUN,fmt='(a80)',end=10) fileKeyword
            read(30+CK_LUN,fmt='(a80)',end=10) fileValue
            select case(trim(fileKeyword))
            case('meshType')
               call getFormatInteger(fileValue, meshType)
            case('controlType')
               call getFormatInteger(fileValue, controlType)
            case('meshFileName')
               meshFileName_g = trim(globaldir) // '/' // fileValue
            case('controlFileName') 
               controlFileName_g = trim(globaldir) // '/' // fileValue
            case default
               ! ignore the keywords and values we are not interested in
            end select
         end do
10       close(30+CK_LUN)
         call mpi_send(meshFileName_g, 2048, mpi_character,
     &      mnproc, tag, comm_writer(wcommID), ierr)
         call mpi_send(controlFileName_g, 2048, mpi_character,
     &      mnproc, tag, comm_writer(wcommID), ierr)
         ! mesh file type
         call mpi_send(meshType, 1, mpi_integer, mnproc, tag,
     &      comm_writer(wcommID), ierr)
         ! control file type
         call mpi_send(controlType, 1, mpi_integer, mnproc, tag,
     &      comm_writer(wcommID), ierr)
         ! ending time step
         call mpi_send(nt, 1, mpi_integer, mnproc, tag,
     &      comm_writer(wcommID), ierr)
      endif
C
#if defined(WRITER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      return
C-----------------------------------------------------------------------
      end subroutine sendInitWriterXDMF
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     S U B R O U T I N E     W R I T E R   R E A D   M E S H 
C-----------------------------------------------------------------------
C     jgf51.21.27: This subroutine causes a dedicated writer processor
C     to read the fulldomain mesh file so that it can later be written
C     into XDMF output files. 
C-----------------------------------------------------------------------
      subroutine writerReadMesh()
      use global
      use mesh, only : readMesh, x, y, slam, sfea
      use sizes, only: MESHFILENAME,CONTROLFILENAME

#ifdef ADCXDMF
      use xdmfio, only : meshInitialized
      use control, only : readControlFile
#endif
      implicit none
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      integer      :: ierr
      integer      :: mpistat(mpi_status_size)
      character(len=8) :: cdigit              !st3 100708: split file

      call setMessageSource("writerReadMesh")
#if defined(WRITER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

#ifdef ADCXDMF
      ! 
      ! receive file name of mesh file
      call mpi_recv(meshFileName, 2048, mpi_character, 0, tag, 
     &   comm_writer(writer_id), mpistat, ierr)
      call mpi_recv(controlFileName, 2048, mpi_character, 0, tag, 
     &   comm_writer(writer_id), mpistat, ierr)
      !
      ! Receive the file type of the mesh file
      call mpi_recv(meshType, 1, mpi_integer, 0, tag, 
     &   comm_writer(writer_id), mpistat, ierr)
      call mpi_recv(controlType, 1, mpi_integer, 0, tag, 
     &   comm_writer(writer_id), mpistat, ierr)

      write(scratchMessage,'(a,i0,a,a)') 'Writer ',writer_id,
     &   ' received mesh file name ',trim(meshFileName)
      call logMessage(INFO, scratchMessage)
      write(scratchMessage,'(a,a,a)') 'Writer GLOBALDIR is ',
     &  trim(GLOBALDIR),'.'
      call logMessage(INFO,scratchMessage)
      ! ending time step
      call mpi_recv(nt, 1, mpi_integer, 0, tag,
     &      comm_writer(writer_id), mpistat, ierr)
      !
      ! now read the mesh so it can be used to create XDMF output files
      if (meshInitialized.eqv..false.) then
         call readMesh()
         ! the x and y will be written to the output filem but the
         ! readMesh routine only populates the slam and sfea arrays
         x = slam
         y = sfea
         meshInitialized = .true.
      endif
      ! read the fort.15 file for metadata     
      call readControlFile(controlFileName_g,.false.)
#endif

#if defined(WRITER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      return
C-----------------------------------------------------------------------
      end subroutine writerReadMesh
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     S U B R O U T I N E   W R I T E R  _  P A U S E
C-----------------------------------------------------------------------
C     Sends a signal to all writer processors to pause, which just
C     means that they exit their infinite loop and return. This 
C     subroutine is only called by processor 0.
C-----------------------------------------------------------------------
      subroutine writer_pause()
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
      !
      ! iterate over the writer processors and use their designated
      ! mpi communicator to signal them to pause
      do i=1,mnwproc
         call mpi_send( sig_pause, 1, mpi_integer, mnproc,
     &                  tag, comm_writer(i), ierr)
      enddo
C
C
#if defined(WRITER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      return
C-----------------------------------------------------------------------
      end subroutine writer_pause
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
      use CkLunMod, only : CK_LUN
        USE SIZES
        USE GLOBAL
        USE MESSENGER, ONLY: TAG, COMM_WRITEH, COMM_HSLEEP, WRITER_ID,
     &       SIG_TERM, MSG_FINI
#ifdef HAVE_MPI_MOD
        use mpi
#endif
      IMPLICIT NONE
      INTEGER :: WCOMMID, WCOMM
!$omp threadprivate(WCOMMID,WCOMM)
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
          WRITE(ScreenUnit+CK_LUN,*) 'HSWRITER PROC ',WRITER_ID,
     &         ' RECEIVED SIG_WRITE'
          CALL writeHotstart_through_HSwriter(0.0d0, 0)

!tcm v49.52.01 added SIG_PAUSE
        ELSE IF(SIG == SIG_PAUSE) then
           WRITE(ScreenUnit+CK_LUN,*) 'HSWRITER PROC ',WRITER_ID,
     &         ' RECEIVED SIG_PAUSE'
           EXIT

        ELSE IF(SIG == SIG_TERM) THEN
          WRITE(ScreenUnit+CK_LUN,*) 'HSWRITER PROC ',WRITER_ID,
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
      USE MESH, ONLY : NE, NP
      USE GLOBAL_IO, ONLY: packOne, unpackOne, 
     &    packTwo, unpackTwo, HEADER_MAX
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
!$omp threadprivate(Elev1Descript)
      type(OutputDataDescript_t), SAVE :: Elev2Descript
!$omp threadprivate(Elev2Descript)
      type(OutputDataDescript_t), SAVE :: VelDescript
!$omp threadprivate(VelDescript)
      type(OutputDataDescript_t), SAVE :: CH1Descript
!$omp threadprivate(CH1Descript)
      type(OutputDataDescript_t), SAVE :: EtaDiscDescript
!$omp threadprivate(EtaDiscDescript)
      type(OutputDataDescript_t), SAVE :: NodeCodeDescript
!$omp threadprivate(NodeCodeDescript)
      type(OutputDataDescript_t), SAVE :: NOFFDescript
!$omp threadprivate(NOFFDescript)
      LOGICAL,  SAVE :: FirstCall = .true.
!$omp threadprivate(FirstCall)
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
        cjunk1 = trim(adjustl(GLOBALDIR))
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
        GLOBALDIR=trim(adjustl(cjunk1))
        hss%filename=trim(adjustl(cjunk2))
        NP = 0
        NE = 0
      endif
      CALL MPI_BARRIER(WCOMM,IERR) !st3

      IF (WRITE_LOCAL_HOT_START_FILES) THEN
         write(6+CK_LUN,*) 'HSWriter can not treat Local HOT START FILES'
         STOP
      ENDIF
!
      WRITE(16+CK_LUN,*) "***********WRITE HOTSTART*********************"
      WRITE(16+CK_LUN,*) "NHSTAR = ", NHSTAR
      WRITE(16+CK_LUN,*) "FirstCall = ", FirstCall
      IF (FirstCall) THEN
         IF ( WRITER_ID /= 0  ) THEN
            IF (.not.ALLOCATED(ETA1_g)) THEN
              write(16+CK_LUN,*) 'Allocating Elev1Descript'
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
               write(16+CK_LUN,*) 'Allocating Elev2Descript'
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
               write(16+CK_LUN,*) 'Allocating VelDescript'
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
               write(16+CK_LUN,*) 'Allocating CH1Descript'
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
            write(16+CK_LUN,*) 'Allocating EtaDiscDescript'
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
            write(16+CK_LUN,*) 'Allocating NodeCodeDescript'
            ALLOCATE(NodeCode_g(NP_G))
         ENDIF
         NodeCodeDescript % specifier            = NHSTAR
         NodeCodeDescript % initial_value        = 0.d0
         NodeCodeDescript % num_items_per_record = 1
         NodeCodeDescript % num_fd_records       = NP_G
         NodeCodeDescript % num_records_this     = NP
         NodeCodeDescript % imap                 => NODES_LG
         NodeCodeDescript % iarray               => NODECODE
         NodeCodeDescript % iarray_g             => NODECODE_g

         IF ( WRITER_ID /= 0 ) THEN
            write(16+CK_LUN,*) 'Allocating NOFFDescript'
            ALLOCATE(NOFF_g(NE_G))
         ENDIF
         NOFFDescript % specifier            = NHSTAR
         NOFFDescript % initial_value        = 0.d0
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
         CALL collectFullDomainIntArrayw(NodeCodeDescript, packOne, unpackOne)
         CALL collectFullDomainIntArrayw(NOFFDescript, packOne, unpackOne)

! tcm v49.54 -- changed noff to noff_g to fix an allocation bug in hot start writer
         WHERE (NOFF_G.gt.1) NOFF_g = 1
!
         WRITE(16+CK_LUN,*) "FINISH collectFullDomainArray"
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
            OPEN(hss % lun +CK_LUN,FILE=TRIM(adjustl(GLOBALDIR))//'/'//trim(adjustl(hss%filename)),
     $      ACCESS='DIRECT',RECL=8)
            IHOTSTP=1
            WRITE(hss % lun+CK_LUN,REC=IHOTSTP) FileFmtVersion
            IHOTSTP = IHOTSTP + 1
            WRITE(hss % lun+CK_LUN,REC=IHOTSTP) IM   ; IHOTSTP = IHOTSTP + 1
            WRITE(hss % lun+CK_LUN,REC=IHOTSTP) ttime; IHOTSTP = IHOTSTP + 1 !st3
            WRITE(hss % lun+CK_LUN,REC=IHOTSTP) iit  ; IHOTSTP = IHOTSTP + 1 !st3
            WRITE(hss % lun+CK_LUN,REC=IHOTSTP) NPX  ; IHOTSTP = IHOTSTP + 1
            WRITE(hss % lun+CK_LUN,REC=IHOTSTP) NEX  ; IHOTSTP = IHOTSTP + 1
            WRITE(hss % lun+CK_LUN,REC=IHOTSTP) NPX  ; IHOTSTP = IHOTSTP + 1
            WRITE(hss % lun+CK_LUN,REC=IHOTSTP) NEX  ; IHOTSTP = IHOTSTP + 1
         ENDIF
         IF ( WRITER_ID /= 0 ) THEN
            DO I=1, Elev1Descript % num_fd_records
               WRITE(hss % lun+CK_LUN,REC=IHOTSTP) Elev1Descript % array_g(I)
               IHOTSTP=IHOTSTP+1
            ENDDO
            DO I=1, Elev2Descript % num_fd_records
               WRITE(hss % lun+CK_LUN,REC=IHOTSTP) Elev2Descript % array_g(I)
               IHOTSTP=IHOTSTP+1
            ENDDO
            DO I=1, EtaDiscDescript % num_fd_records
               WRITE(hss % lun+CK_LUN,REC=IHOTSTP) EtaDiscDescript % array_g(I)
               IHOTSTP=IHOTSTP+1
            ENDDO
            DO I=1, VelDescript % num_fd_records
               WRITE(hss % lun+CK_LUN,REC=IHOTSTP) VelDescript % array_g(I)
               IHOTSTP=IHOTSTP+1
            ENDDO
            DO I=1, VelDescript % num_fd_records
               WRITE(hss % lun+CK_LUN,REC=IHOTSTP) VelDescript % array2_g(I)
               IHOTSTP=IHOTSTP+1
            ENDDO
            IF (IM.eq.10) THEN
               DO I=1, CH1Descript % num_fd_records
                  WRITE(hss % lun+CK_LUN,REC=IHOTSTP) CH1Descript % array_g(I)
                  IHOTSTP=IHOTSTP+1
               ENDDO
            ENDIF
            DO I=1, NodeCodeDescript % num_fd_records
              WRITE(hss % lun+CK_LUN,REC=IHOTSTP)NodeCodeDescript % iarray_g(I)
              IHOTSTP=IHOTSTP+1
            ENDDO
            DO I=1, NOFFDescript % num_fd_records
               WRITE(hss % lun+CK_LUN,REC=IHOTSTP) NOFFDescript % iarray_g(I)
               IHOTSTP=IHOTSTP+1
            ENDDO
         ENDIF
C
         IF (WRITER_ID /= 0) THEN
            WRITE(hss % lun+CK_LUN,REC=IHOTSTP) IESTP   ; IHOTSTP = IHOTSTP + 1
            WRITE(hss % lun+CK_LUN,REC=IHOTSTP) NSCOUE  ; IHOTSTP = IHOTSTP + 1

            WRITE(hss % lun+CK_LUN,REC=IHOTSTP) IVSTP   ; IHOTSTP = IHOTSTP + 1
            WRITE(hss % lun+CK_LUN,REC=IHOTSTP) NSCOUV  ; IHOTSTP = IHOTSTP + 1

            WRITE(hss % lun+CK_LUN,REC=IHOTSTP) ICSTP   ; IHOTSTP = IHOTSTP + 1
            WRITE(hss % lun+CK_LUN,REC=IHOTSTP) NSCOUC  ; IHOTSTP = IHOTSTP + 1

            WRITE(hss % lun+CK_LUN,REC=IHOTSTP) IPSTP   ; IHOTSTP = IHOTSTP + 1
            WRITE(hss % lun+CK_LUN,REC=IHOTSTP) IWSTP   ; IHOTSTP = IHOTSTP + 1
            WRITE(hss % lun+CK_LUN,REC=IHOTSTP) NSCOUM  ; IHOTSTP = IHOTSTP + 1

            WRITE(hss % lun+CK_LUN,REC=IHOTSTP) IGEP    ; IHOTSTP = IHOTSTP + 1
            WRITE(hss % lun+CK_LUN,REC=IHOTSTP) NSCOUGE ; IHOTSTP = IHOTSTP + 1

            WRITE(hss % lun+CK_LUN,REC=IHOTSTP) IGVP    ; IHOTSTP = IHOTSTP + 1
            WRITE(hss % lun+CK_LUN,REC=IHOTSTP) NSCOUGV ; IHOTSTP = IHOTSTP + 1

            WRITE(hss % lun+CK_LUN,REC=IHOTSTP) IGCP    ; IHOTSTP = IHOTSTP + 1
            WRITE(hss % lun+CK_LUN,REC=IHOTSTP) NSCOUGC ; IHOTSTP = IHOTSTP + 1

            WRITE(hss % lun+CK_LUN,REC=IHOTSTP) IGPP    ; IHOTSTP = IHOTSTP + 1
            WRITE(hss % lun+CK_LUN,REC=IHOTSTP) IGWP    ; IHOTSTP = IHOTSTP + 1
            WRITE(hss % lun+CK_LUN,REC=IHOTSTP) NSCOUGW ; IHOTSTP = IHOTSTP + 1
         ENDIF
      CASE(2) ! ascii?
         write(16+CK_LUN,*) 'ASCII Hotstart file not implemented yet.'

      CASE(3) ! netcdf
         write(16+CK_LUN,*) 'NetCDF Hotstart undergoing testing.'

      CASE DEFAULT
         write(ScreenUnit+CK_LUN,*) 'The value of NHSTAR=',NHSTAR
         write(ScreenUnit+CK_LUN,*) 'is not supported.'
         write(ScreenUnit+CK_LUN,*) 'Hotstart file not written.'
      END SELECT
C...
C...    CLOSE THE HOT START OUTPUT FILE
C...
      IF(ABS(NHSTAR).EQ.1) THEN                                              ! added by mcf 7/14/08 tcm v51.26 added abs(nhstar)
         if (WRITER_ID /= 0) CLOSE(hss % lun+CK_LUN)
      ENDIF                                                             ! end NHSTAR=1

      IF(NSCREEN.NE.0.AND. WRITER_ID /= 0 ) THEN
         WRITE(ScreenUnit+CK_LUN,24541) hss % lun, WRITER_ID, IIT, ttime  !tcm v51.26 changed IT to IIT and TimeLoc to ttime since
         WRITE(16+CK_LUN,24541) hss % lun, WRITER_ID, IIT, ttime !IT and TimeLoc are not known by writer processor id /=0  (echoing to writer file)
      ELSE
         WRITE(16+CK_LUN,24541) hss % lun,WRITER_ID,IT,TimeLoc
      ENDIF
24541 FORMAT(1X,'HOT START OUTPUT WRITTEN TO UNIT ',I2,'  HSWriter:',i2
     &    ' AT TIME STEP = ',I9,' TIME = ',E15.8)
      IF(hss % lun.EQ.67) THEN
C        jgf45.07 added option to stop ADCIRC after writing hot start file.
         IF (NHSTAR.EQ.67) THEN
            WRITE(ScreenUnit+CK_LUN,*) 'NHSTAR.EQ.67; ADCIRC stopping.'
            WRITE(16+CK_LUN,*) 'NHSTAR.EQ.67; ADCIRC stopping.'
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
            WRITE(ScreenUnit+CK_LUN,*) 'NHSTAR.EQ.68; ADCIRC stopping.'
            WRITE(16+CK_LUN,*) 'NHSTAR.EQ.68; ADCIRC stopping.'
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
#ifdef CMPI      
#ifndef HAVE_MPI_MOD      
      include 'mpif.h'
#endif      
#endif
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
!$omp threadprivate(tagbase)
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
!$omp threadprivate(tagbase)
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
         integerBuffer(:)  = int(descript % initial_value)
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
!
      END MODULE HSWRITER
