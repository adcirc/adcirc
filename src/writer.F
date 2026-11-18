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
#include "logging_macros.h"
C******************************************************************************
C MODULE WRITER
C  This module is for writer processors which receives values from compute
C  processors and writes them out asynchronously.
C                                                              05 Feb 2007, sb
C******************************************************************************
C
      module writer
      use mpi_f08, only: MPI_DOUBLE_PRECISION, MPI_Comm
      use mod_sizes, only : mnproc, mnwproc, OFF, ASCII, SPARSE_ASCII,
     &     NETCDF3, NETCDF4, XDMF, iSplit, globalDir, meshType,
     &     controlType, meshFileName_g, controlFileName_g
      use global, only: c3d, comm,
     &    bufsize_max, buf, resultbuf, outputdatadescript_t,
     &     sig_val, sig_write, sig_ini, sig_pause, comm_writer,
     &     sig_term, sig_mesh, writer_id
      use mod_logging, only: debug, echo, info, warning, error,
     &     logmessage, allmessage, t_log_scope, init_log_scope,
     &     nabout
      use messenger, only: tag
      implicit none
      ! tcm v50.27 increased size from 3 to 4, cobell increase to 6,
      ! jgf51.21.27 decreased size of intbuf to 5
      ! jgf51.21.29 increased size of intbuf to 6 for XDMF
      ! jgf52.08.25: increased intbuf to 7 to transmit nabout
      integer, parameter :: numIntMetadata = 7
      integer :: intBuf(numIntMetadata)  ! used to transmit pieces of descript via mpi
      integer, parameter :: numRealMetadata = 2
      real(8) :: reBuf(numRealMetadata) ! used to transmit pieces of descript via mpi
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
      ! 
      ! .true. if a writer processor has received data into one or more
      ! of its local buffers but has not written it yet; the nwloaded
      ! array is dimensioned equal to the number of writer processors;
      ! this is an array that is updated and used by compute processors
      ! not writer processors
      logical, allocatable :: nwloaded(:)   !st3 100708: round robin chamber
      !
      ! buffer number on a particular writer processor; each writer 
      ! processor manages its own bufID so it knows where to store the
      ! next dataset it receives
      integer :: bufID    
      !
      ! wcommID is incremented by and used only by compute processors to
      ! determine which writer processor to send data to each time
      integer :: wcommID  ! writer processor number (1 -> MNWPROC)

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

      LOG_SCOPE_TRACED("writer_init", WRITER_TRACING)
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
          use mpi_f08, only: MPI_Send, MPI_Recv, MPI_INTEGER, MPI_Status
          use mod_sizes,only: myproc, mnwproc, mnproc
          use messenger,only: mpi_comm_adcirc
          use global,only: labels_g, np_g
          implicit none
          INTEGER :: ierr
          INTEGER :: I,N
          type(MPI_Status) :: mpistat

          IF(myproc.EQ.0)THEN
              DO I = 0,MNWPROC-1
                 CALL mpi_send(sig_ini,1,mpi_integer,mnproc,tag,
     &                         comm_writer(i+1))
                 CALL mpi_send(np_g,1,mpi_integer,mnproc+i,0,
     &                         mpi_comm_adcirc)
                 CALL mpi_send(labels_g,np_g,mpi_integer,mnproc+i,0,
     &                         mpi_comm_adcirc)
              ENDDO
          ELSEIF(myproc.GT.mnproc-1.and.myproc.LT.mnproc+mnwproc)THEN
              CALL mpi_recv(n,1,mpi_integer,0,0,mpi_comm_adcirc,mpistat)
              if (.not.allocated(labels_g)) allocate(labels_g(1:n))   !tcm 20240804 added test for allocated
              CALL mpi_recv(labels_g,n,mpi_integer,0,0,mpi_comm_adcirc,mpistat)
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
      use mod_logging, only : ScreenUnit
      use mpi_f08, only: MPI_Recv, MPI_INTEGER, MPI_Status
      use messenger, only: msg_barrier, mpi_comm_adcirc
      IMPLICIT NONE
      integer :: sig, ierr, waittime, nextProc, prevProc, sig_dummy
      type(MPI_Status) :: mpistat
      character(len=16) :: TimeIni
      character(len=1024) :: scratchMessage
C
      LOG_SCOPE_TRACED("writer_main", WRITER_TRACING)
      !
      ! infinite loop

      do
         call mpi_recv(sig, 1, mpi_integer, 0, tag,
     &       comm_writer(writer_id), mpistat)

         select case(sig)
         case(sig_ini) 
            write(scratchMessage,
     &       '("Writer processor ",i0," received sig_ini.")') writer_id
            call logMessage(INFO,scratchMessage)
            call sendLabelInfoToWriters()
         case(sig_val) 
            write(scratchMessage,
     &       '("Writer processor ",i0," received sig_val.")') writer_id
            call logMessage(INFO,scratchMessage)
            call writer_recv_values()
         case(sig_write)
            write(scratchMessage,
     &      '("Writer processor ",i0," received sig_write.")') writer_id
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
      use mod_sizes, only : myproc
      use global
      use mod_logging, only: allmessage, debug,
     &     logmessage, t_log_scope, init_log_scope
      use mpi_f08, only: MPI_Send, MPI_Reduce, MPI_INTEGER,
     &                   MPI_CHARACTER, MPI_SUM
      implicit none
      type (OutputDataDescript_t), intent(in) :: descript
      real(8) :: TimeLoc  ! output time in seconds
      integer :: it      ! current time step
      external store_cmd ! command to pack the buffer for mpi transmission
C
      integer :: ierr
      integer :: num
      integer :: bufsize
      integer :: istart, iend
      character(1024) :: scratchMessage
C
      LOG_SCOPE_TRACED("sendDataToWriter", WRITER_TRACING)
      !
      ! jgf51.21.27: Figure out which mpi communicator to use with this
      ! writer processor.
      if (wcommID.ge.mnwproc) wcommID = 0
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
     &      comm_writer(wcommID))

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
         INTBUF(7) = nabout%level()
         call mpi_send(intbuf, numIntMetaData, mpi_integer, mnproc, tag,
     &      comm_writer(wcommID))
         REBUF(1) = descript % alternate_value
         REBUF(2) = timeLoc
         ! send dry value and time in seconds
         call mpi_send(rebuf, numRealMetaData, MPI_DOUBLE_PRECISION, mnproc, tag,
     &      comm_writer(wcommID))
         ! send file name
         call mpi_send(descript % file_name, 1024, mpi_character,
     &      mnproc, tag, comm_writer(wcommID))
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
         call mpi_reduce(buf, resultBuf, bufsize, MPI_DOUBLE_PRECISION, MPI_SUM,
     &      MNPROC, comm_writer(wcommID))
         istart = iend + 1
         iend   = min(istart + num - 1, descript % num_fd_records)
         num    = iend - istart + 1
      end do
C
C-----------------------------------------------------------------------
      end subroutine sendDataToWriter


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
      use mpi_f08, only: MPI_Recv, MPI_Reduce, MPI_INTEGER,
     &                   MPI_CHARACTER, MPI_SUM, MPI_Status
      use mod_logging, only: allmessage, debug,
     &     logmessage, t_log_level, t_log_scope, init_log_scope
      implicit none
      integer :: num, i, j
      integer :: ierr
      integer :: bufsize
      integer :: istart, iend
      type(MPI_Status) :: mpistat
      logical :: descriptFound ! true if we found the particular buffer we need
      character(len=8) :: cdigit  !st3 100708: split file
      character(len=1024) :: scratchMessage

      LOG_SCOPE_TRACED("writer_recv_values", WRITER_TRACING)
      !
      ! Receive minimal set of metadata required to write the file
      call mpi_recv(intbuf, numIntMetaData, mpi_integer, 0, tag,
     &   comm_writer(writer_id), mpistat)
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
      nabout = t_log_level(intbuf(7))
      !
      ! receive "dry" data value and current simulation time
      call mpi_recv(rebuf, numRealMetaData, MPI_DOUBLE_PRECISION, 0, tag,
     &   comm_writer(writer_id), mpistat)
      writer_descripts(bufid) % alternate_value = rebuf(1)
      writer_descripts(bufid) % wMD % myTime = rebuf(2)
      !
      ! receive file name
      writer_descripts(bufID) % file_name(:) = " "
      call mpi_recv(writer_descripts(bufID)%file_name, 1024,
     & mpi_character, 0, tag, comm_writer(writer_id), mpistat)
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
         call mpi_reduce(buf, resultBuf, bufsize, MPI_DOUBLE_PRECISION, MPI_SUM,
     &     MNPROC, comm_writer(writer_id))

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
      use mod_sizes, only : myproc
      use global
      use mod_logging, only: allmessage, debug,
     &     logmessage, t_log_scope, init_log_scope
      use mpi_f08, only: MPI_Send, MPI_INTEGER
      implicit none
      INTEGER :: I, IERR
      character(1024) :: scratchMessage

      LOG_SCOPE_TRACED("flush_writers", WRITER_TRACING)
      !
      ! return immediately if no writer processor has received data 
      ! since the last flush
      if ( .not.any(nwloaded).eqv..true. ) then
         return
      endif
      !
      ! processor 0 loops over the writer processors and tells the ones
      ! that have received data in one or more of their buffers that they
      ! need to write it all now
      if(myproc.eq.0) then
         do i=1,mnwproc
            !st3 100708: flush only loaded chamber
            if (nwloaded(i).eqv..true.) then   !st3 nwloaded()=false:unloaded, true:loaded
               write(scratchMessage,'(a,i0,a,i0,a)') 'Proc ',myproc,
     &           ' is sending sig_write to writer ',i,'.'
               call logMessage(DEBUG, scratchMessage)
               call mpi_send(sig_write,1,mpi_integer,mnproc,
     &              tag,comm_writer(i))
               nwloaded(i) = .false.
            endif
         enddo
      endif

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
      use mod_logging, only: allmessage, debug, info,
     &     error, logmessage, t_log_scope, init_log_scope
#ifdef ADCNETCDF
      use NETCDFIO, only : initNetCDFOutputFile, writeOutArrayNetCDF
#endif
#ifdef ADCXDMF
      use XDMFIO, only : initOutputXDMF, writeOutArrayXDMF, writeControlXDMF
#endif
      implicit none
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
      character(len=1024) :: scratchMessage
      logical :: nerr             ! netcdf error status indicator

      LOG_SCOPE_TRACED("writer_write", WRITER_TRACING)
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
            open(fid, file=writer_descripts(ibuf) % file_name, 
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
            write(fid, 100) writer_descripts(ibuf) % wMD % myTime,
     &            writer_descripts(ibuf) % wMD % myTimeStep 
            ! write dataset
            select case(writer_descripts(ibuf)%num_items_per_record)
            case(1)
               do i=1, writer_descripts(ibuf) % num_fd_records
                  write(fid,1000) labels_g(i), writer_descripts(ibuf) % array_g(i)
               end do
            case(2)
               do i=1, writer_descripts(ibuf) % num_fd_records
                  write(fid,2000) labels_g(i), writer_descripts(ibuf) % array_g(i),
     &               writer_descripts(ibuf) % array2_g(i)
               end do
            case default
               write(scratchMessage,'(a,i0,a)') 'Writer cannot write ',
     &            writer_descripts(ibuf) % num_items_per_record,
     &            ' component data.'
               call allMessage(ERROR,scratchMessage)
            end select
            close(fid)
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
               write(fid, 100) writer_descripts(ibuf) % wMD % myTime,
     &               writer_descripts(ibuf) % wMD % myTimeStep, 
     &               nNonDefVal, 
     &               writer_descripts(ibuf) % alternate_value
               ! write sparse scalar dataset
               do i=1, writer_descripts(ibuf) % num_fd_records
                  if (areDefault1(i).eqv..false.) then
                     write(fid,1000) labels_g(i), 
     &                     writer_descripts(ibuf) % array_g(i)
                  endif
               end do
               close(fid)
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
               write(fid, 100) writer_descripts(ibuf) % wMD % myTime,
     &               writer_descripts(ibuf) % wMD % myTimeStep, 
     &               nNonDefVal, 
     &               writer_descripts(ibuf) % alternate_value
               ! write sparse vector dataset
               do i=1, writer_descripts(ibuf) % num_fd_records
                  if ((areDefault1(i).eqv..false.).or.
     &                   (areDefault2(i).eqv..false.)) then
                     write(fid,2000) labels_g(i), 
     &                     writer_descripts(ibuf) % array_g(i),
     &                     writer_descripts(ibuf) % array2_g(i)
                  endif
               end do
               close(fid)
               deallocate(areDefault1, areDefault2)
            case default
               write(scratchMessage,'(a,i0,a)') 'Writer cannot write ',
     &               writer_descripts(ibuf) % num_items_per_record,
     &               ' component data.'
               call allMessage(ERROR, scratchMessage)
            end select
            close(fid)
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
      use mod_sizes, only : myproc
      use global
      use mod_io, only: openFileForRead
      use mod_logging, only: allmessage, debug, info,
     &     logmessage, t_log_scope, init_log_scope
      use mpi_f08, only: MPI_Send, MPI_INTEGER, MPI_CHARACTER
      implicit none
      type (OutputDataDescript_t), intent(in) :: descript
      character(len=80) :: fileKeyword
      character(len=80) :: fileValue
      character(len=1024) :: scratchMessage
      integer :: ierr
C
      LOG_SCOPE_TRACED("sendInitWriterXDMF", WRITER_TRACING)
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
     &      comm_writer(wcommID))
         ! 
         ! read the fulldomain mesh file name and control file name
         ! so that the writer can load them up and use them (it only
         ! uses the control file for generating metadata) ... the 
         ! 'fulldomainInputFiles' file is created by adcprep
         call openFileForRead(30,trim(globaldir)//'/fulldomainInputFiles',
     &      ierr, required=.FALSE.)
         do 
            read(30,fmt='(a80)',end=10) fileKeyword
            read(30,fmt='(a80)',end=10) fileValue
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
10       close(30)
         call mpi_send(meshFileName_g, 2048, mpi_character,
     &      mnproc, tag, comm_writer(wcommID))
         call mpi_send(controlFileName_g, 2048, mpi_character,
     &      mnproc, tag, comm_writer(wcommID))
         ! mesh file type
         call mpi_send(meshType, 1, mpi_integer, mnproc, tag,
     &      comm_writer(wcommID))
         ! control file type
         call mpi_send(controlType, 1, mpi_integer, mnproc, tag,
     &      comm_writer(wcommID))
         ! ending time step
         call mpi_send(nt, 1, mpi_integer, mnproc, tag,
     &      comm_writer(wcommID))
      endif
C
C-----------------------------------------------------------------------
      end subroutine sendInitWriterXDMF


C-----------------------------------------------------------------------
C     S U B R O U T I N E     W R I T E R   R E A D   M E S H 
C-----------------------------------------------------------------------
C     jgf51.21.27: This subroutine causes a dedicated writer processor
C     to read the fulldomain mesh file so that it can later be written
C     into XDMF output files. 
C-----------------------------------------------------------------------
      subroutine writerReadMesh()
      use global
      use mod_logging, only: allmessage, debug, info,
     &     logmessage, t_log_scope, init_log_scope
      use mesh, only : readMesh, x, y, slam, sfea
      use mod_sizes, only: MESHFILENAME,CONTROLFILENAME
      use mpi_f08, only: MPI_Recv, MPI_INTEGER, MPI_CHARACTER, MPI_Status
#ifdef ADCXDMF
      use xdmfio, only : meshInitialized
      use control, only : readControlFile
#endif
      implicit none
      integer      :: ierr
      type(MPI_Status) :: mpistat
      character(len=8) :: cdigit              !st3 100708: split file
      character(len=1024) :: scratchMessage

      LOG_SCOPE_TRACED("writerReadMesh", WRITER_TRACING)

#ifdef ADCXDMF
      ! 
      ! receive file name of mesh file
      call mpi_recv(meshFileName, 2048, mpi_character, 0, tag, 
     &   comm_writer(writer_id), mpistat)
      call mpi_recv(controlFileName, 2048, mpi_character, 0, tag,
     &   comm_writer(writer_id), mpistat)
      !
      ! Receive the file type of the mesh file
      call mpi_recv(meshType, 1, mpi_integer, 0, tag,
     &   comm_writer(writer_id), mpistat)
      call mpi_recv(controlType, 1, mpi_integer, 0, tag,
     &   comm_writer(writer_id), mpistat)

      write(scratchMessage,'(a,i0,a,a)') 'Writer ',writer_id,
     &   ' received mesh file name ',trim(meshFileName)
      call logMessage(INFO, scratchMessage)
      write(scratchMessage,'(a,a,a)') 'Writer GLOBALDIR is ',
     &  trim(GLOBALDIR),'.'
      call logMessage(INFO,scratchMessage)
      ! ending time step
      call mpi_recv(nt, 1, mpi_integer, 0, tag,
     &      comm_writer(writer_id), mpistat)
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
      use mod_logging, only: allmessage, debug,
     &     t_log_scope, init_log_scope
      use messenger, only: tag, writers_active
      use mpi_f08, only: MPI_Send, MPI_INTEGER
      implicit none
      integer i
      integer ierr

      LOG_SCOPE_TRACED("writer_pause", WRITER_TRACING)

      if(.not.writers_active)return
      writers_active = .false.
      !
      ! iterate over the writer processors and use their designated
      ! mpi communicator to signal them to pause
      do i=1,mnwproc
         call mpi_send( sig_pause, 1, mpi_integer, mnproc,
     &                  tag, comm_writer(i))
      enddo
C
C
      return
C-----------------------------------------------------------------------
      end subroutine writer_pause
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      END MODULE WRITER
C-----------------------------------------------------------------------
