      !--------------------------------------------------------------
      ! Global IO module.
      !
      !  This module exists to write out global output files on 
      !  parallel machines and one file on uniproc machines.
      !--------------------------------------------------------------

      module global_io
      use SIZES
      use VERSION
      use GLOBAL
C      use GLOBAL, ONLY: NP, NODES_LG, IMAP_EL_LG,
C     &    IMAP_STAE_LG, IMAP_STAV_LG, IMAP_STAM_LG, IMAP_STAC_LG,
C     &    BUFSIZE_MAX, buf, resultBuf, float_type, OutputDataDescript_t,
C     &    COMM, nodecode
#ifdef CMPI
#ifdef HAVE_MPI_MOD
      use mpi
#endif
#endif
      contains

      !--------------------------------------------------------------
      !                  S U B R O U T I N E   
      !     C O L L E C T  F U L L  D O M A I N  A R R A Y 
      !--------------------------------------------------------------
      ! jgf48.03 Collects array data from each subdomain. 
      !--------------------------------------------------------------
      subroutine collectFullDomainArray(descript, pack_cmd, unpack_cmd)
      USE SIZES
      USE GLOBAL
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
      ! the subroutine used to write the file
      integer      :: ierr, status(MPI_STATUS_SIZE), request
      integer, save:: tagbase = 6000
      integer      :: iproc
      integer      :: bufsize
      integer      :: ibucket    
      integer      :: istart     ! vector tuple to start with
      integer      :: iend       ! vector tuple to end on
      integer      :: tag        
      ! number of vector tuples in the buffer
      integer      :: num        
      integer      :: i, j, k
C
      call setMessageSource("collectFullDomainArray")
#if defined(GLOBALIO_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      bufsize = min(BUFSIZE_MAX, 
     &    descript % num_items_per_record * descript % num_fd_records)
C     num will be less than the number of full domain records if the 
C     buffer is too small to hold all the records in the full domain,
C     in this case it is the number of records passed back to proc 0
C     on each iteration of the while loop below
      num     = bufsize / descript % num_items_per_record
      iend    = num
      istart  = 1

      if (tagbase == 5000) then
         tagbase = 6000
      else
         tagbase = 5000
      endif
      ibucket = 0

      do while (istart.le.iend)

         !------------------------------------------------------------
         ! Initialize
         !------------------------------------------------------------
         buf(:)  = descript % initial_value
         ibucket = ibucket + 1
         tag     = tagbase + mod(ibucket, 8)

C        now pack the buffer 
         call pack_cmd(descript, istart, iend)
C        now send data to processor 0
         call mpi_reduce(buf, resultBuf, bufsize, float_type, MPI_SUM,
     &                   0, COMM, ierr)
         if (myproc == 0) then
            call unpack_cmd(descript, istart, iend)
         end if
C        set new starting position to just after the
C        current ending position in full domain array
         istart = iend + 1
C        set new ending position to either the current start plus the
C        number of records that will fit in the buffer (minus 1), 
C        or to the end of the full domain array, whichever is less
         iend   = min(istart + num - 1, descript % num_fd_records)
         num    = iend - istart + 1
      end do
#endif   
C! CMPI

 1000 format(2x, i8, 2x, 1pE20.10E3, 1pE20.10E3, 1pE20.10E3, 1pE20.10E3)
 1100 FORMAT(2x,1pE20.10E3,5X,I10)
#if defined(GLOBALIO_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      !--------------------------------------------------------------
      end subroutine collectFullDomainArray
      !--------------------------------------------------------------

      !--------------------------------------------------------------
      !                  S U B R O U T I N E   
      !     C O L L E C T  F U L L  D O M A I N  I N T  A R R A Y 
      !--------------------------------------------------------------
      ! jgf48.17 Collects integer array data from each subdomain. 
      !--------------------------------------------------------------
      subroutine collectFullDomainIntArray(descript, 
     &                                    pack_cmd, unpack_cmd,variable)
      USE SIZES
      USE GLOBAL
      implicit none
#ifdef CMPI
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
#endif
      type (OutputDataDescript_t) :: descript
      external pack_cmd
      external unpack_cmd
      integer :: variable
#ifdef CMPI
      ! the subroutine used to write the file
      integer      :: ierr, status(MPI_STATUS_SIZE), request
      integer, save:: tagbase = 6000
      integer      :: iproc
      integer      :: bufsize
      integer      :: ibucket    
      integer      :: istart     ! vector tuple to start with
      integer      :: iend       ! vector tuple to end on
      integer      :: tag        
      ! number of vector tuples in the buffer
      integer      :: num        
      integer      :: i, j, k

      bufsize = min(BUFSIZE_MAX, 
     &    descript % num_items_per_record * descript % num_fd_records)
C     num will be less than the number of full domain records if the 
C     buffer is too small to hold all the records in the full domain,
C     in this case it is the number of records passed back to proc 0
C     on each iteration of the while loop below
      num     = bufsize / descript % num_items_per_record
      iend    = num
      istart  = 1
      if (tagbase == 5000) then
         tagbase = 6000
      else
         tagbase = 5000
      endif
      ibucket = 0

      do while (istart.le.iend)

         !------------------------------------------------------------
         ! Initialize
         !------------------------------------------------------------
         integerBuffer(:)  = descript % int_initial_value
         ibucket = ibucket + 1
         tag     = tagbase + mod(ibucket, 8)

C        now pack the buffer 
         call pack_cmd(descript, istart, iend)
C        now send data to processor 0
ckmd48.33bc - changed this to eliminate problems with NOFF=2.
         IF (variable.EQ.1) THEN 
           call mpi_reduce(integerBuffer, integerResultBuffer, bufsize, 
     &                 MPI_INTEGER, MPI_SUM, 0, COMM, ierr)
         ELSE IF (variable.eq.2) THEN 
           call mpi_reduce(integerBuffer, integerResultBuffer, bufsize, 
     &                 MPI_INTEGER, MPI_MAX, 0, COMM, ierr)
         END IF 

         if (myproc == 0) then
            call unpack_cmd(descript, istart, iend)
         end if
C        set new starting position to just after the
C        current ending position in full domain array
         istart = iend + 1
C        set new ending position to either the current start plus the
C        number of records that will fit in the buffer (minus 1), 
C        or to the end of the full domain array, whichever is less
         iend   = min(istart + num - 1, descript % num_fd_records)
         num    = iend - istart + 1
      end do
#endif   
C! CMPI

 1000 format(2x, i8, 2x, 1pE20.10E3, 1pE20.10E3, 1pE20.10E3, 1pE20.10E3)
 1100 FORMAT(2x,1pE20.10E3,5X,I10)
      !--------------------------------------------------------------
      end subroutine collectFullDomainIntArray
      !--------------------------------------------------------------

      !--------------------------------------------------------------
      ! S U B R O U T I N E   P A C K  O N E
      !--------------------------------------------------------------
      !  jgf48.03 Subroutine to store a single array of real numbers
      !  to a buffer 
      !--------------------------------------------------------------
      subroutine packOne(descript, istart, iend)
      USE SIZES
      USE GLOBAL
      implicit none

      type (OutputDataDescript_t) :: descript
      integer :: i, j, istart, iend, iglobal
      integer :: ioffset

#ifdef CMPI
      ioffset = istart - 1
      do i = 1, descript % num_records_this
        iglobal = descript % imap(i)
        if (istart <= iglobal .and. iglobal <= iend) then
          buf(iglobal-ioffset) = descript % array(i)
        end if
      end do
#endif

      return
      !--------------------------------------------------------------
      end subroutine packOne
      !--------------------------------------------------------------

      !--------------------------------------------------------------
      ! S U B R O U T I N E   P A C K  O N E  I N T
      !--------------------------------------------------------------
      !  jgf48.17 Subroutine to store a single array of integers
      !  to a buffer 
      !--------------------------------------------------------------
      subroutine packOneInt(descript, istart, iend)
      USE SIZES
      USE GLOBAL
      implicit none

      type (OutputDataDescript_t) :: descript
      integer :: i, j, istart, iend, iglobal
      integer :: ioffset

#ifdef CMPI
      ioffset = istart - 1
      do i = 1, descript % num_records_this
        iglobal = descript % imap(i)
        if (istart <= iglobal .and. iglobal <= iend) then
          integerBuffer(iglobal-ioffset) = descript % iarray(i)
        end if
      end do
#endif

      return
      !--------------------------------------------------------------
      end subroutine packOneInt
      !--------------------------------------------------------------

      !--------------------------------------------------------------
      !       S U B R O U T I N E   P A C K 6 3 
      !--------------------------------------------------------------
      !  jgf48.03 Subroutine similar to packOne, except dry node 
      !  elevations are set to -99999.0 during the packing process.
      !--------------------------------------------------------------
      subroutine pack63(descript, istart, iend)
      USE SIZES
      USE GLOBAL
      implicit none

      type (OutputDataDescript_t) :: descript
      integer :: i, istart, iend, iglobal
      integer :: ioffset

#ifdef CMPI
      ioffset = istart - 1
      do i = 1, np
         iglobal = descript % imap(i)
         if (istart <= iglobal .and. iglobal <= iend) then
            if (nodecode(i) == 1) then
               buf(iglobal - ioffset ) = descript % array(i) !eta2(i)
            else
               buf(iglobal - ioffset ) = -99999.0
            endif
         endif
      end do
#endif

      return
 1000 format(2x, i8, 2x, 1pE20.10E3)
      !--------------------------------------------------------------
      end subroutine pack63
      !--------------------------------------------------------------

      !--------------------------------------------------------------
      !       S U B R O U T I N E   P A C K 8 3 
      !--------------------------------------------------------------
      !  jgf48.03 Subroutine similar to packOne, except dry node 
      !  elevations are set to -99999.0 during the packing process.
      !--------------------------------------------------------------
      subroutine pack83(descript, istart, iend)
      USE SIZES
      USE GLOBAL
      implicit none

      type (OutputDataDescript_t) :: descript
      integer :: i, istart, iend, iglobal
      integer :: ioffset

#ifdef CMPI
      ioffset = istart - 1
      do i = 1, np
         iglobal = descript % imap(i)
         if (istart <= iglobal .and. iglobal <= iend) then
            if (nodecode(i) == 1) then
               buf(iglobal - ioffset ) = descript % array(i) !eta2(i)
            else
               buf(iglobal - ioffset ) = -99999.0
            endif
         endif
      end do
#endif

      return
 1000 format(2x, i8, 2x, 1pE20.10E3)
      !--------------------------------------------------------------
      end subroutine pack83
      !--------------------------------------------------------------

      !--------------------------------------------------------------
      ! S U B R O U T I N E   U N P A C K  O N E
      !--------------------------------------------------------------
      !  Subroutine to retrieve a single array of real numbers
      !  from a buffer 
      !--------------------------------------------------------------
      subroutine unpackOne(descript, istart, iend)
      USE SIZES
      USE GLOBAL
      implicit none

      type (OutputDataDescript_t) :: descript
      integer :: i, j, istart, iend, iglobal
      integer :: ioffset

       j = 1
       do i = istart, iend
          descript % array_g(i) = resultBuf(j)
          j = j + 1
      end do
      return
      !--------------------------------------------------------------
      end subroutine unpackOne
      !--------------------------------------------------------------

      !--------------------------------------------------------------
      ! S U B R O U T I N E   U N P A C K  O N E  I N T
      !--------------------------------------------------------------
      !  Subroutine to retrieve a single array of integers
      !  from a buffer 
      !--------------------------------------------------------------
      subroutine unpackOneInt(descript, istart, iend)
      USE SIZES
      USE GLOBAL
      implicit none

      type (OutputDataDescript_t) :: descript
      integer :: i, j, istart, iend, iglobal
      integer :: ioffset

       j = 1
       do i = istart, iend
          descript % iarray_g(i) = integerResultBuffer(j)
          j = j + 1
      end do
      return
      !--------------------------------------------------------------
      end subroutine unpackOneInt
      !--------------------------------------------------------------

      !--------------------------------------------------------------
      ! S U B R O U T I N E   P A C K  T W O
      !--------------------------------------------------------------
      ! Subroutine to pack two interleaved arrays of real
      ! numbers into a buffer
      !--------------------------------------------------------------
      subroutine packTwo(descript, istart, iend)
      USE SIZES
      USE GLOBAL
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
 1000 format(2x, i8, 2x, 1pE20.10E3, 1pE20.10E3)
      !--------------------------------------------------------------
      end subroutine packTwo
      !--------------------------------------------------------------

      !--------------------------------------------------------------
      ! S U B R O U T I N E  U N P A C K  T W O
      !--------------------------------------------------------------
      ! Subroutine to unpack two interleaved arrays of real
      ! numbers into a buffer
      !--------------------------------------------------------------
      subroutine unpackTwo(descript, istart, iend)
      USE SIZES
      USE GLOBAL
      implicit none

      type (OutputDataDescript_t) :: descript
      integer :: i, j, istart, iend

      j = 1
      do i = istart, iend
         descript % array_g(i) = resultBuf(j)
         descript % array2_g(i) = resultBuf(j+1)
         j = j + 2
      end do
      return
      !--------------------------------------------------------------
      end subroutine unpackTwo
      !--------------------------------------------------------------

      !--------------------------------------------------------------
      ! S U B R O U T I N E   P A C K  M  B Y  N P 
      !--------------------------------------------------------------
      ! Subroutine to pack an m x np array of real
      ! numbers into a buffer
      !--------------------------------------------------------------
      subroutine packMbyNP(descript, istart, iend)
      USE SIZES
      USE GLOBAL
      implicit none

      type (OutputDataDescript_t) :: descript
      integer :: i, istart, iend, iglobal
      integer :: ioffset, j
      integer :: k ! frequency counter
C
      call setMessageSource("packMbyNP")
#if defined(GLOBALIO_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      ioffset = istart - 1
      do i = 1, descript % num_records_this
         iglobal = descript % imap(i)
         if (istart <= iglobal .and. iglobal <= iend) then
            j = descript%num_items_per_record * (iglobal-ioffset)
     &          - (descript%num_items_per_record - 1 )
            do k = 0, ( descript % num_items_per_record - 1 )
               buf(j+k) = descript % array2D(k+1,i)
            end do
         end if
      end do
C
 1000 format(2x, i8, 2x, 1pE20.10E3, 1pE20.10E3)
#if defined(GLOBALIO_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      return
      !--------------------------------------------------------------
      end subroutine packMbyNP
      !--------------------------------------------------------------

      !--------------------------------------------------------------
      ! S U B R O U T I N E  U N P A C K  M  B Y  N P
      !--------------------------------------------------------------
      ! Subroutine to unpack an m x np array of real
      ! numbers out of a buffer
      !--------------------------------------------------------------
      subroutine unpackMbyNP(descript, istart, iend)
      USE SIZES
      USE GLOBAL
      implicit none

      type (OutputDataDescript_t) :: descript
      integer :: i, j, istart, iend
      integer k ! frequency counter

      call setMessageSource("unpackMbyNP")
#if defined(GLOBALIO_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      j = 1
      do i = istart, iend
         do k = 1, descript % num_items_per_record
            descript % array2D_g(k,i) = resultBuf(j)
            j = j + 1
         end do
      end do
#if defined(GLOBALIO_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      return
      !--------------------------------------------------------------
      end subroutine unpackMbyNP
      !--------------------------------------------------------------

      !--------------------------------------------------------------
      ! S U B R O U T I N E   P A C K  N P  B Y  M 
      !--------------------------------------------------------------
      ! Subroutine to pack an np x m array of real
      ! numbers into a buffer
      !--------------------------------------------------------------
      subroutine packNPbyM(descript, istart, iend)
      USE SIZES
      USE GLOBAL
      implicit none

      type (OutputDataDescript_t) :: descript
      integer :: i, istart, iend, iglobal
      integer :: ioffset, j
      integer :: k ! frequency counter
C
      call setMessageSource("packNPbyM")
#if defined(GLOBALIO_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      ioffset = istart - 1
      do i = 1, descript % num_records_this
         iglobal = descript % imap(i)
         if (istart <= iglobal .and. iglobal <= iend) then
            j = descript%num_items_per_record * (iglobal-ioffset)
     &          - (descript%num_items_per_record - 1 )
            do k = 0, ( descript % num_items_per_record - 1 )
               buf(j+k) = descript % array2D(i,k+1)
            end do
         end if
      end do
C
 1000 format(2x, i8, 2x, 1pE20.10E3, 1pE20.10E3)
#if defined(GLOBALIO_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      return
      !--------------------------------------------------------------
      end subroutine packNPbyM
      !--------------------------------------------------------------

      !--------------------------------------------------------------
      ! S U B R O U T I N E  U N P A C K   N P  B Y  M
      !--------------------------------------------------------------
      ! Subroutine to unpack an np x m array of real
      ! numbers out of a buffer
      !--------------------------------------------------------------
      subroutine unpackNPbyM(descript, istart, iend)
      USE SIZES
      USE GLOBAL
      implicit none

      type (OutputDataDescript_t) :: descript
      integer :: i, j, istart, iend
      integer k ! frequency counter

      call setMessageSource("unpackNPbyM")
#if defined(GLOBALIO_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      j = 1
      do i = istart, iend
         do k = 1, descript % num_items_per_record
            descript % array2D_g(i,k) = resultBuf(j)
            j = j + 1
         end do
      end do
#if defined(GLOBALIO_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      return
      !--------------------------------------------------------------
      end subroutine unpackNPbyM
      !--------------------------------------------------------------

      !--------------------------------------------------------------   
      !    S U B R O U T I N E   O P E N _ G B L _ F I L E 
      !--------------------------------------------------------------
      ! Open Global File: opens the file, opens the write_header
      !                   routine and closes the file
      !--------------------------------------------------------------
      subroutine open_gbl_file(lun, filename, size_g, size_this,
     &                          write_header)
      USE SIZES
      USE GLOBAL
      implicit none

      external write_header ! subroutine used to actually write the header
      integer      :: lun, size_this, size_g, szz
      character(*) :: filename

      szz = size_g
      if (mnproc == 1 .or. WRITE_LOCAL_FILES) szz = size_this

      if (myproc == 0 .or. WRITE_LOCAL_FILES) then
        open(lun, file=filename, status='UNKNOWN')
        call write_header(lun,szz)
        close(lun)
      end if
      return
      !--------------------------------------------------------------
      end subroutine open_gbl_file
      !--------------------------------------------------------------

      !--------------------------------------------------------------   
      !    S U B R O U T I N E   O P E N _ M I N M A X _ F I L E 
      !--------------------------------------------------------------
      ! Open Global File: opens the file, replacing the existing file,
      !                   executes the write_header routine and closes
      !                   the file.
      !--------------------------------------------------------------
      subroutine open_minmax_file(lun, filename, size_g, size_this,
     &                          write_header)
      USE SIZES
      USE GLOBAL
      implicit none

      external write_header ! subroutine used to actually write the header
      integer      :: lun, size_this, size_g, szz
      character(*) :: filename

      szz = size_g
      if (mnproc == 1 .or. WRITE_LOCAL_FILES) szz = size_this

      if (myproc == 0 .or. WRITE_LOCAL_FILES) then
        open(lun, file=filename, status='REPLACE')
        call write_header(lun,szz)
        close(lun)
      end if
      return
      !--------------------------------------------------------------
      end subroutine open_minmax_file
      !--------------------------------------------------------------


      !--------------------------------------------------------------
      !     S U B R O U T I N E   W R I T E _ G B L _ F I L E
      !--------------------------------------------------------------
      ! Write Global File : It opens the file.
      !   On uniproc machines it calls the "write_cmd"
      !   routine to write the record out.
      !
      !   On parallel machine it performs a bucket
      !   brigade by copying a block of data to proc 0
      !   which it writes out.
      !--------------------------------------------------------------

      subroutine write_gbl_file(lun, filename, descript, TimeLoc, it, 
     &                          write_cmd)
      USE SIZES
      USE GLOBAL
      implicit none
#ifdef CMPI
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
#endif
      integer lun
      character(*) :: filename
      type (OutputDataDescript_t) :: descript
      real(8)  :: TimeLoc
      integer  :: it
      ! the subroutine used to write the file
      external write_cmd
#ifdef CMPI
      integer      :: ierr, status(MPI_STATUS_SIZE), request
#endif
      integer, save:: tagbase = 6000
      integer      :: iproc
      integer      :: bufsize
      integer      :: ibucket    
      integer      :: istart     ! vector tuple to start with
      integer      :: iend       ! vector tuple to end on
      integer      :: tag        
      ! number of vector tuples in the buffer
      integer      :: num        
      integer      :: i, j, k

C      write(16,*) 'About to open file.' !jgfdebug48.03
      if (myproc == 0 .or. WRITE_LOCAL_FILES) then
         open(lun, file=filename, access='SEQUENTIAL',
     &        position='APPEND')
         write(lun, 1100) TimeLoc, it
      endif

      if (WRITE_LOCAL_FILES) then
         call write_cmd(lun, descript, 1, descript % num_records_this)
         close(lun)
         return
      end if

#ifdef CMPI
      bufsize = min(BUFSIZE_MAX, 
     &    descript % num_items_per_record * descript % num_fd_records)
C     num will be less than the number of full domain records if the 
C     buffer is too small to hold all the records in the full domain,
C     in this case it is the number of records passed back to proc 0
C     on each iteration of the while loop below
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

         !------------------------------------------------------------
         ! Initialize
         !------------------------------------------------------------
         buf(:)  = descript % initial_value
         ibucket = ibucket + 1
         tag     = tagbase + mod(ibucket, 8)
         
C         write(16,*) 'About to call write_cmd.' !jgfdebug48.03
         call write_cmd(lun, descript, istart, iend)
C         write(16,*) 'About to call mpi_reduce.' !jgfdebug48.03 
         call mpi_reduce(buf, resultBuf, bufsize, float_type, MPI_SUM,
     &                   0, COMM, ierr)
C         write(16,*) 'Completed call to mpi_reduce.' !jgfdebug48.03
         j = 1
         if (myproc == 0) then
            do i = istart, iend
!               write(16,*) i, j, k !jgfdebug48.03                      ! commented out by mcf : k has not been defined
               write(lun, 1000) i, (resultBuf(k), k = j, j +
     &                          descript % num_items_per_record - 1)
               j = j + descript % num_items_per_record
            end do
         end if

C        set new starting position to just after the
C        current ending position in full domain array
         istart = iend + 1
C        set new ending position to either the current start plus the
C        number of records that will fit in the buffer (minus 1), 
C        or to the end of the full domain array, whichever is less
         iend   = min(istart + num - 1, descript % num_fd_records)
         num    = iend - istart + 1
      end do
      if (myproc == 0) then
         close(lun)
      endif
#endif   
C! CMPI

 1000 format(2x, i8, 2x, 1pE20.10E3, 1pE20.10E3, 1pE20.10E3, 1pE20.10E3)
 1100 FORMAT(2x,1pE20.10E3,5X,I10)
      !--------------------------------------------------------------
      end subroutine write_gbl_file
      !--------------------------------------------------------------


      !--------------------------------------------------------------
      !      S U B R O U T I N E   
      !      W R I T E _ G B L _ 6 3 _ S K I P _ D R Y N O D E
      !--------------------------------------------------------------
      ! Write Global File Skip Default : sb 11/10/2006
      !      File writing subroutine for fort.63
      !      It opens the file.
      !      On uniproc machines it calls the "write_cmd"
      !      routine to write the record out.
      !
      !      On parallel machine it performs a bucket
      !      brigade by copying a block of data to proc 0
      !      which it writes out.
      !
      !      The specified default values will not be written 
      !      out.
      !--------------------------------------------------------------

      subroutine write_gbl_63_skip_drynode
     &           (lun, filename, descript, TimeLoc, it, write_cmd)
      USE SIZES
      USE GLOBAL
      implicit none
#ifdef CMPI
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
#endif

      external write_cmd
      integer      :: num, i, j, k
#ifdef CMPI
      integer      :: ierr, status(MPI_STATUS_SIZE), request
#endif
      integer, save:: tagbase = 6000
      integer      :: lun, iproc
      integer      :: bufsize, ibucket
      integer      :: istart, iend, it, tag
#ifdef CMPI
      integer      :: nLWetNodes, nGWetNodes
      integer      :: iglobal
#endif
      real(8)      :: TimeLoc
      character(*) :: filename
      type (OutputDataDescript_t) :: descript


      if(descript%num_records_this.ne.np) then
        stop 'FATAL ERROR: descript%num_records_this.ne.np'
      endif

#ifndef CMPI
C     !-----IF CMPI IS NOT DEFINED -----
      open(lun, file=filename, access='SEQUENTIAL', position='APPEND')
      write(lun, 1100) TimeLoc, it

      call write_cmd(lun, descript, 1, descript % num_records_this)
      close(lun)
      return

#else
C    !-----IF CMPI IS DEFINED -----

      if (WRITE_LOCAL_FILES) then
        open(lun, file=filename, access='SEQUENTIAL', position='APPEND')
        write(lun, 1100) TimeLoc, it
        call write_cmd(lun, descript, 1, descript % num_records_this)
        close(lun)

        return !---THERE IS A RETURN STATEMENT HERE---

      end if

      ! Count the number of wet nodes
      nLWetNodes = 0
      do i=1,np
        iglobal = descript%imap(i)
        if(iglobal.le.0) cycle ! cycle if node i is ghost

        if(nodecode(i) == 0) cycle ! cycle if node i is dry

        nLWetNodes=nLWetNodes+1
      enddo

      call mpi_reduce(nLWetNodes, nGWetNodes, 1, MPI_INTEGER,
     &     MPI_SUM, 0, comm, ierr)

      if (myproc == 0) then
        open(lun, file=filename, access='SEQUENTIAL', position='APPEND')
        write(lun, 1101) TimeLoc, it, nGWetNodes, -99999.d0
      endif

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

        !------------------------------------------------------------
        ! Initialize
        !------------------------------------------------------------
        buf(:)  = descript % initial_value
        ibucket = ibucket + 1
        tag     = tagbase + mod(ibucket, 8)

        call write_cmd(lun, descript, istart, iend)
        call mpi_reduce(buf, resultBuf, bufsize, float_type, MPI_SUM, 0,
     &    COMM, ierr)
        if (myproc == 0) then
          do i = istart, iend
            j = 1 + (i-istart)*descript % num_items_per_record

            !If values are eqaul to the default value, skip this node
            if(resultBuf(j) == -99999.0) cycle

            write(lun, 1000) i, (resultBuf(k), k = j, j +
     &        descript % num_items_per_record - 1)
          end do
        end if

        istart = iend + 1
        iend   = min(istart + num - 1, descript % num_fd_records)
        num    = iend - istart + 1
      end do
      if (myproc == 0) then
        close(lun)
      endif
#endif   
C! CMPI

 1000 format(2x, i8, 2x, 1pE20.10E3, 1pE20.10E3, 1pE20.10E3, 1pE20.10E3)
 1100 FORMAT(2x,1pE20.10E3,5X,I10)
 1101 FORMAT(2x,1pE20.10E3,5X,I10,5X,I10,5X,1pE20.10E3)      
      !--------------------------------------------------------------
      end subroutine write_gbl_63_skip_drynode
      !--------------------------------------------------------------


      !--------------------------------------------------------------
      !       S U B R O U T I N E   
      !       W R I T E _ G B L _ F I L E _ S K I P _ D E F A U L T
      !--------------------------------------------------------------
      ! Write Global File Skip Default : sb 11/10/2006
      !       It opens the file.
      !       On uniproc machines it calls the "write_cmd"
      !       routine to write the record out.
      !
      !       On parallel machine it performs a bucket
      !       brigade by copying a block of data to proc 0
      !       which it writes out.
      !
      !       The specified default values will not be written 
      !       out.
      !--------------------------------------------------------------

      subroutine write_gbl_file_skip_default
     &           (lun, filename, descript, TimeLoc, it, write_cmd, default)
      USE SIZES
      USE GLOBAL
      implicit none
#ifdef CMPI
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
#endif

      external write_cmd
      integer      :: num, i, j, k
#ifdef CMPI
      integer      :: ierr, status(MPI_STATUS_SIZE), request
#endif
      integer, save:: tagbase = 6000
      integer      :: lun, iproc
      integer      :: bufsize, ibucket
      integer      :: istart, iend, it, tag
#ifdef CMPI
      integer      :: nLNonDefValues, nGNonDefValues
      integer      :: iglobal
#endif
      real(8)      :: default
      real(8)      :: TimeLoc
      character(*) :: filename
      type (OutputDataDescript_t) :: descript

      if(descript%num_records_this.ne.np) then
        stop 'FATAL ERROR: descript%num_records_this.ne.np'
      endif

#ifndef CMPI 
C    !-----IF CMPI IS NOT DEFINED -----
      open(lun, file=filename, access='SEQUENTIAL', position='APPEND')
      write(lun, 1100) TimeLoc, it

      call write_cmd(lun, descript, 1, descript % num_records_this)
      close(lun)
      return

#else   
C     !-----IF CMPI IS DEFINED -----

      if (WRITE_LOCAL_FILES) then
        open(lun, file=filename, access='SEQUENTIAL', position='APPEND')
        write(lun, 1100) TimeLoc, it
        call write_cmd(lun, descript, 1, descript % num_records_this)
        close(lun)

        return !---THERE IS A RETURN STATEMENT HERE---

      end if

      ! Count the number of non-default values
      nLNonDefValues = 0
      do i=1,np
        iglobal = descript%imap(i)
        if(iglobal.le.0) cycle ! cycle if node i is ghost

        if(descript%num_items_per_record.eq.1) then
          if(descript%array (i).eq.default) cycle
        else if(descript%num_items_per_record.eq.2) then
          if(descript%array (i).eq.default.and.
     &       descript%array2(i).eq.default) cycle
        endif

        nLNonDefValues=nLNonDefValues+1
      enddo

      call mpi_reduce(nLNonDefValues, nGNonDefValues, 1, MPI_INTEGER,
     &     MPI_SUM, 0, comm, ierr)

      if (myproc == 0) then
        open(lun, file=filename, access='SEQUENTIAL', position='APPEND')
        write(lun, 1101) TimeLoc, it, nGNonDefValues, default
      endif

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

        !------------------------------------------------------------
        ! Initialize
        !------------------------------------------------------------
        buf(:)  = descript % initial_value
        ibucket = ibucket + 1
        tag     = tagbase + mod(ibucket, 8)

        call write_cmd(lun, descript, istart, iend)
        call mpi_reduce(buf, resultBuf, bufsize, float_type, MPI_SUM, 0,
     &    COMM, ierr)
        if (myproc == 0) then
          do i = istart, iend
            j = 1 + (i-istart)*descript % num_items_per_record

            !If values are eqaul to the default value, skip this node
            if(descript % num_items_per_record == 1) then
              if(resultBuf(j) == default) cycle
            else if(descript % num_items_per_record == 2) then
              if((resultBuf(j)   == default).and.
     &           (resultBuf(j+1) == default)) cycle
            endif

            write(lun, 1000) i, (resultBuf(k), k = j, j +
     &        descript % num_items_per_record - 1)
          end do
        end if

        istart = iend + 1
        iend   = min(istart + num - 1, descript % num_fd_records)
        num    = iend - istart + 1
      end do
      if (myproc == 0) then
        close(lun)
      endif
#endif   
C! CMPI

 1000 format(2x, i8, 2x, 1pE20.10E3, 1pE20.10E3, 1pE20.10E3, 1pE20.10E3)
 1100 FORMAT(2x,1pE20.10E3,5X,I10)
 1101 FORMAT(2x,1pE20.10E3,5X,I10,5X,I10,5X,1pE20.10E3)
      !--------------------------------------------------------------
      end subroutine write_gbl_file_skip_default
      !--------------------------------------------------------------

      !--------------------------------------------------------------
      ! These are the store, header routines. For each file (fort.61, ...) 
      ! that is being globalize, there is a pair of routines : header61 and
      ! store61 etc.  The header routine writes out the header once
      ! and the store routine copies the data (heights, velocities, etc)
      ! to the buffer in parallel mode or in serial mode it just write out
      ! the data to file.
      !--------------------------------------------------------------


      !--------------------------------------------------------------
      ! S U B R O U T I N E   S T O R E  O N E
      !--------------------------------------------------------------
      ! jgf47.05: Subroutine to store a single array of real numbers
      !   +to disk in text format if WRITE_LOCAL_FILES is .true.
      !   +to a buffer if we are running in parallel
      !--------------------------------------------------------------
      subroutine storeOne(lun, descript, istart, iend)
      USE SIZES
      USE GLOBAL
      implicit none

      type (OutputDataDescript_t) :: descript
      integer :: i, istart, iend, iglobal, lun
      integer :: ioffset

      if (WRITE_LOCAL_FILES) then
         do i = 1, descript % num_records_this
            write(lun , 1000) i, descript % array(i)
         end do
         return
      end if

#ifdef CMPI
      ioffset = istart - 1
      do i = 1, descript % num_records_this
         iglobal = descript % imap(i)
         if (istart <= iglobal .and. iglobal <= iend) then
            buf(iglobal-ioffset) = descript % array(i)
         end if
      end do
#endif

      return
 1000 format(2x, i8, 2x, 1pE20.10E3)
      !--------------------------------------------------------------
      end subroutine storeOne
      !--------------------------------------------------------------


      !--------------------------------------------------------------
      ! S U B R O U T I N E   S T O R E  T W O
      !--------------------------------------------------------------
      ! jgf47.05: Subroutine to store two interleaved arrays of real
      ! numbers
      !  +to disk in text format if WRITE_LOCAL_FILES is .true.
      !  +to a buffer if we are running in parallel
      !--------------------------------------------------------------
      subroutine storeTwo(lun, descript, istart, iend)
      USE SIZES
      USE GLOBAL
      implicit none

      type (OutputDataDescript_t) :: descript
      integer :: i, istart, iend, iglobal, lun
      integer :: ioffset, j

      if (WRITE_LOCAL_FILES) then
        do i = 1, descript % num_records_this
          write(lun , 1000) i, descript % array(i),
     $      descript % array2(i)
        end do
        return
      end if

#ifdef CMPI
      ioffset = istart - 1
      do i = 1, descript % num_records_this
        iglobal = descript % imap(i)
        if (istart <= iglobal .and. iglobal <= iend) then
          j = 2*(iglobal - ioffset) - 1 
          buf(j)     = descript % array(i)
          buf(j + 1) = descript % array2(i)
        end if
      end do
#endif
      return
 1000 format(2x, i8, 2x, 1pE20.10E3, 1pE20.10E3)
      !--------------------------------------------------------------
      end subroutine storeTwo
      !--------------------------------------------------------------


      !--------------------------------------------------------------
      !       S U B R O U T I N E   S T O R E 6 3 
      !--------------------------------------------------------------
      !
      !--------------------------------------------------------------
      subroutine store63(lun, descript, istart, iend)
      USE SIZES
      USE GLOBAL
      implicit none

      type (OutputDataDescript_t) :: descript
      integer :: i, istart, iend, iglobal, lun
      integer :: ioffset


      if (WRITE_LOCAL_FILES) then
        if (descript % ConsiderWetDry.eqv..true.) then
           do i = 1, np
             if (nodecode(i) == 1) write(lun, 1000) i, descript % array(i)  !eta2(i)
             if (nodecode(i) == 0) write(lun, 1000) i, descript % alternate_value
           end do
        else
           do i = 1,np
              write(lun , 1000) i, descript % array(i)
           end do
        endif
        return
      end if

#ifdef CMPI
      ioffset = istart - 1
      if (descript % ConsiderWetDry.eqv..true.) then
         do i = 1, np
            iglobal = descript % imap(i)
            if (istart <= iglobal .and. iglobal <= iend) then
               if (nodecode(i) == 1) buf(iglobal - ioffset ) = descript % array(i)  !eta2(i)
               if (nodecode(i) == 0) buf(iglobal - ioffset ) = descript % alternate_value
           end if
         end do
      else
         do i = 1, np
            iglobal = descript % imap(i)
            if (istart <= iglobal .and. iglobal <= iend) then
               buf(iglobal - ioffset ) = descript % array(i)
           end if
         end do
      endif
#endif

      return
 1000 format(2x, i8, 2x, 1pE20.10E3)
      !--------------------------------------------------------------
      end subroutine store63
      !--------------------------------------------------------------


      !--------------------------------------------------------------
      !     S U B R O U T I N E   H E A D E R 6 1 
      !--------------------------------------------------------------
      ! Elevation for Elevation Recording station: fort.61
      !--------------------------------------------------------------
      subroutine header61(lun, size_this)
      USE SIZES
      USE GLOBAL
      implicit none
      integer :: lun, size_this

      write(lun, 1000) rundes, runid, agrid
      write(lun, 1010) ntrspe, size_this, dtdp*nspoole, nspoole, 1, 
     &                 FileFmtVersion

 1000 FORMAT(1X,A32,2X,A24,2X,A24)
 1010 FORMAT(1X,I10,1X,I10,1X,E15.7E3,1X,I5,1X,I5,1X,
     &       'FileFmtVersion: ',I10)
      !--------------------------------------------------------------
      end subroutine header61
      !--------------------------------------------------------------

      !--------------------------------------------------------------
      !    S U B R O U T I N E   H E A D E R 6 2
      !--------------------------------------------------------------
      ! Velocity for Velocity Recording station: fort.62
      !--------------------------------------------------------------
      subroutine header62(lun, size_this)
      USE SIZES
      USE GLOBAL
      implicit none
      integer :: lun, size_this

      write(lun, 1000) rundes, runid, agrid
      write(lun, 1010) ntrspv, size_this, dtdp*nspoolv, nspoolv, 2,
     &     FileFmtVersion

 1000 FORMAT(1X,A32,2X,A24,2X,A24)
 1010 FORMAT(1X,I10,1X,I10,1X,E15.7E3,1X,I5,1X,I5,1X,
     &      'FileFmtVersion: ',I10)
      !--------------------------------------------------------------
      end subroutine header62
      !--------------------------------------------------------------

      !--------------------------------------------------------------
      !    S U B R O U T I N E   H E A D E R 6 3
      !--------------------------------------------------------------
      ! Elevation for all nodes: fort.63
      !--------------------------------------------------------------
      subroutine header63(lun, size_this)
      USE SIZES
      USE GLOBAL
      implicit none
      integer :: lun, size_this

      write(lun, 1000) rundes, runid, agrid
      write(lun, 1010) ndsetse, size_this, dtdp*nspoolge, nspoolge, 1,
     &                 FileFmtVersion

 1000 FORMAT(1X,A32,2X,A24,2X,A24)
 1010 FORMAT(1X,I10,1X,I10,1X,E15.7E3,1X,I5,1X,I5,1X,
     &       'FileFmtVersion: ',I10)
      !--------------------------------------------------------------
      end subroutine header63
      !--------------------------------------------------------------

      !--------------------------------------------------------------
      ! S U B R O U T I N E   H E A D E R _ M A X
      !--------------------------------------------------------------
      ! Maximum elevation for all nodes: fort.63
      !--------------------------------------------------------------
      subroutine header_max(lun, size_this)
      USE SIZES
      USE GLOBAL
      implicit none
      integer :: lun, size_this

      write(lun, 1000) rundes, runid, agrid
      write(lun, 1010) 1, size_this, 1.d0, 1, 1, FileFmtVersion

 1000 FORMAT(1X,A32,2X,A24,2X,A24)
 1010 FORMAT(1X,I10,1X,I10,1X,E15.7E3,1X,I5,1X,I5,1X,
     &      'FileFmtVersion: ',I10)
      !--------------------------------------------------------------
      end subroutine header_max
      !--------------------------------------------------------------

      !--------------------------------------------------------------
      !  S U B R O U T I N E   H E A D E R 6 4 
      !--------------------------------------------------------------
      ! Velocity for all nodes: fort.64
      !--------------------------------------------------------------
      subroutine header64(lun, size_this)
      USE SIZES
      USE GLOBAL
      implicit none
      integer :: lun, size_this

      write(lun, 1000) rundes, runid, agrid
      write(lun, 1010) ndsetsv, size_this, dtdp*nspoolgv, nspoolgv, 2, 
     &                 FileFmtVersion

 1000 FORMAT(1X,A32,2X,A24,2X,A24)
 1010 FORMAT(1X,I10,1X,I10,1X,E15.7E3,1X,I5,1X,I5,1X,
     &      'FileFmtVersion: ',I10)
      !--------------------------------------------------------------
      end subroutine header64
      !--------------------------------------------------------------

      !--------------------------------------------------------------
      ! S U B R O U T I N E   H E A D E R 7 1
      !--------------------------------------------------------------
      ! Pressure for Meteorological Recording station: fort.71
      !--------------------------------------------------------------
      subroutine header71(lun, size_this)
      USE SIZES
      USE GLOBAL
      implicit none
      integer :: lun, size_this

      write(lun, 1000) rundes, runid, agrid
      write(lun, 1010) ntrspm, size_this, dtdp*nspoolm, nspoolm, 1, 
     &                 FileFmtVersion

 1000 FORMAT(1X,A32,2X,A24,2X,A24)
 1010 FORMAT(1X,I10,1X,I10,1X,E15.7E3,1X,I5,1X,I5,1X,
     &       'FileFmtVersion: ',I10)
      !--------------------------------------------------------------
      end subroutine header71
      !--------------------------------------------------------------

      !--------------------------------------------------------------
      ! S U B R O U T I N E   H E A D E R 7 2
      !--------------------------------------------------------------
      ! Wind Velocity for Meteorological Recording stations : fort.72
      !--------------------------------------------------------------
      subroutine header72(lun, size_this)
      USE SIZES
      USE GLOBAL
      implicit none
      integer :: lun, size_this

      write(lun, 1000) rundes, runid, agrid
      write(lun, 1010) ntrspm, size_this, dtdp*nspoolm, nspoolm, 2,
     &                  FileFmtVersion

 1000 FORMAT(1X,A32,2X,A24,2X,A24)
 1010 FORMAT(1X,I10,1X,I10,1X,E15.7E3,1X,I5,1X,I5,1X,
     &       'FileFmtVersion: ',I10)
      !--------------------------------------------------------------
      end subroutine header72
      !--------------------------------------------------------------

      !--------------------------------------------------------------
      ! S U B R O U T I N E   H E A D E R 7 3
      !--------------------------------------------------------------
      ! Atmospheric Pressure for all nodes: fort.73
      !--------------------------------------------------------------
      subroutine header73(lun, size_this)
      USE SIZES
      USE GLOBAL
      implicit none
      integer :: lun, size_this

      write(lun, 1000) rundes, runid, agrid
      write(lun, 1010) ndsetsw, size_this, dtdp*nspoolgw, nspoolgw, 1,
     &                 FileFmtVersion

 1000 FORMAT(1X,A32,2X,A24,2X,A24)
 1010 FORMAT(1X,I10,1X,I10,1X,E15.7E3,1X,I5,1X,I5,1X,
     &       'FileFmtVersion: ',I10)
      !--------------------------------------------------------------
      end subroutine header73
      !--------------------------------------------------------------

      !--------------------------------------------------------------
      ! S U B R O U T I N E   H E A D E R 7 4      
      !--------------------------------------------------------------
      ! Wind Stress for all nodes: fort.74
      !--------------------------------------------------------------
      subroutine header74(lun, size_this)
      USE SIZES
      USE GLOBAL
      implicit none
      integer :: lun, size_this

      write(lun, 1000) rundes, runid, agrid
      write(lun, 1010) ndsetsw, size_this, dtdp*nspoolgw, nspoolgw, 2, 
     &                 FileFmtVersion

 1000 FORMAT(1X,A32,2X,A24,2X,A24)
 1010 FORMAT(1X,I10,1X,I10,1X,E15.7E3,1X,I5,1X,I5,1X,
     &       'FileFmtVersion: ',I10)
      !--------------------------------------------------------------
      end subroutine header74
      !--------------------------------------------------------------

      !--------------------------------------------------------------
      ! S U B R O U T I N E   H E A D E R 9 1
      !--------------------------------------------------------------
      ! Ice Concentration for Meteorological Recording station: fort.91
      ! tcm v49.64.01 -- added
      !--------------------------------------------------------------
      subroutine header91(lun, size_this)
      USE SIZES
      USE GLOBAL
      implicit none
      integer :: lun, size_this

      write(lun, 1000) rundes, runid, agrid
      write(lun, 1010) ntrspm, size_this, dtdp*nspoolm, nspoolm, 1, 
     &                 FileFmtVersion

 1000 FORMAT(1X,A32,2X,A24,2X,A24)
 1010 FORMAT(1X,I10,1X,I10,1X,E15.7E3,1X,I5,1X,I5,1X,
     &       'FileFmtVersion: ',I10)
      !--------------------------------------------------------------
      end subroutine header91
      !--------------------------------------------------------------

      !--------------------------------------------------------------
      ! S U B R O U T I N E   H E A D E R 9 3
      !--------------------------------------------------------------
      ! Ice Concentration Field for all nodes: fort.93
      ! tcm v49.64.01 -- added
      !--------------------------------------------------------------
      subroutine header93(lun, size_this)
      USE SIZES
      USE GLOBAL
      implicit none
      integer :: lun, size_this

      write(lun, 1000) rundes, runid, agrid
      write(lun, 1010) ndsetsw, size_this, dtdp*nspoolgw, nspoolgw, 1,
     &                 FileFmtVersion

 1000 FORMAT(1X,A32,2X,A24,2X,A24)
 1010 FORMAT(1X,I10,1X,I10,1X,E15.7E3,1X,I5,1X,I5,1X,
     &       'FileFmtVersion: ',I10)
      !--------------------------------------------------------------
      end subroutine header93
      !--------------------------------------------------------------


      !--------------------------------------------------------------
      ! S U B R O U T I N E   H E A D E R 8 1
      !--------------------------------------------------------------
      ! Concentration Recording station: fort.81
      !--------------------------------------------------------------
      subroutine header81(lun, size_this)
      USE SIZES
      USE GLOBAL
      implicit none
      integer :: lun, size_this

      write(lun, 1000) rundes, runid, agrid
      write(lun, 1010) ntrspc, size_this, dtdp*nspoolc, nspoolc, 1,
     &                  FileFmtVersion

 1000 FORMAT(1X,A32,2X,A24,2X,A24)
 1010 FORMAT(1X,I10,1X,I10,1X,E15.7E3,1X,I5,1X,I5,1X,
     &      'FileFmtVersion: ',I10)
      !--------------------------------------------------------------
      end subroutine header81
      !--------------------------------------------------------------

      !--------------------------------------------------------------
      ! S U B R O U T I N E  
      !    W R I T E  S T A T I O N  H E A D E R 
      !--------------------------------------------------------------
      ! jgf47.05: This subroutine will write the header for station
      ! files (scalars or vectors, 2D or 3D). This is intended to 
      ! eventually replace the header subroutines "headerXX" for 
      ! station files.
      !--------------------------------------------------------------
      subroutine writeStationHeader(lun, numSta, typeStr)
      USE SIZES
      USE GLOBAL
      use global_3dvs, only: 
     &     ndset3dsv, nspo3dsv, nfen,
     &     ndset3dst, nspo3dst,
     &     ndset3dsd, nspo3dsd
      implicit none
      integer lun               ! logical unit number to write to
      integer numSta            ! number of stations to write
      character(len=10) typeStr ! the type of header to write 
       
      write(lun, 1000) rundes, runid, agrid

      select case(trim(typeStr))
      case("Elev")
          write(lun, 1010) ntrspe, numSta, dtdp*nspoole, nspoole, 1, 
     &         FileFmtVersion
      case("Vel2D")
          write(lun, 1010) ntrspv, numSta, dtdp*nspoolv, nspoolv, 2,
     &         FileFmtVersion
      case("Vel3D")
          write(lun, 1020) ndset3dsv, numSta, dtdp*nspo3dsv, 
     &         nspo3dsv, nfen, 3, FileFmtVersion
      case("Turb3D")
          write(lun, 1020) ndset3dst, numSta, dtdp*nspo3dst, 
     &         nspo3dst, nfen, 3, FileFmtVersion
      case("Dens2D")
          write(ScreenUnit,*) 
     &         'ERROR: 2D density station output not implemented.'
          stop
      case("Dens3D")
          write(lun, 1020) ndset3dsd, numSta, dtdp*nspo3dsd, 
     &         nspo3dsd, nfen, 3, FileFmtVersion
      case("Press","Wind")
          write(lun, 1010) ntrspm, numSta, dtdp*nspoolm, nspoolm, 1, 
     &                 FileFmtVersion
      case("Conc2D")
          write(lun, 1010) ntrspc, numSta, dtdp*nspoolc, nspoolc, 1,
     &                  FileFmtVersion
      case("Conc3D")
          write(ScreenUnit,*) 
     &         'ERROR: 3D conc. station output not implemented.'
          stop
      case default
          write(ScreenUnit,*) 
     &         'ERROR: Station header type unrecongnized.'
          stop
      end select

 1000 FORMAT(1X,A32,2X,A24,2X,A24)
 1010 FORMAT(1X,I10,1X,I10,1X,E15.7E3,1X,I5,1X,I5,1X,
     &      'FileFmtVersion: ',I10)
 1020 FORMAT(1X,I10,1X,I10,1X,E15.7E3,1X,I5,1X,I5,1X,I5,1X,
     &      'FileFmtVersion: ',I10)
      !--------------------------------------------------------------
      end subroutine writeStationHeader
      !--------------------------------------------------------------

      !--------------------------------------------------------------
      ! S U B R O U T I N E  
      !    W R I T E  D O M A I N  H E A D E R 
      !--------------------------------------------------------------
      ! jgf47.05: This subroutine will write the header for domain 
      ! output files, i.e., not station files. It handles scalars or 
      ! vectors, 2D or 3D). This is intended to eventually replace 
      ! the header subroutines "headerXX" for domain files.
      !--------------------------------------------------------------
      subroutine writeDomainHeader(lun, filename, numNodes_g, 
     &            numNodes, typeStr)
      USE SIZES
      USE GLOBAL
      use global_3dvs, only:
     &   ndset3dgv, nspo3dgv, nfen,
     &   ndset3dgt, nspo3dgt,
     &   ndset3dgd, nspo3dgd
      implicit none
      integer lun               ! logical unit number to write to
      character(*) filename     ! name of file to create
      integer numNodes          ! number of nodes in this domain
      integer numNodes_g        ! number of nodes in the full domain
      character(len=10) typeStr ! the type of header to write 
      integer nodes             ! number of nodes to write into header
C
C     Unless we are processor 0 or we are supposed to write local
C     files, there is nothing to do but *** RETURN EARLY ***
      if ( (myproc.ne.0).and.(.not.WRITE_LOCAL_FILES) ) return 
C
C     Set the appropriate number of nodes
      if ( (myproc.eq.1).or.(WRITE_LOCAL_FILES) ) then
         nodes = numNodes
      else
         nodes = numNodes_g
      endif
C
C     Open the file
      open(lun, file=filename, status='UNKNOWN')
      write(lun, 1000) rundes, runid, agrid

      select case(trim(typeStr))
      case("Elev")
          write(lun, 1010) ndsetse, nodes, dtdp*nspoolge,
     &                 nspoolge, 1, FileFmtVersion
      case("ElevMax")
          write(lun, 1010) 1, numNodes, 1.d0, 1, 1, FileFmtVersion
      case("Vel2D")
          write(lun, 1010) ndsetsv, nodes, dtdp*nspoolgv,
     &                     nspoolgv, 2, FileFmtVersion
      case("Vel3D")
          write(lun, 1020) ndset3dgv, nodes, dtdp*nspo3dgv,
     &                     nspo3dgv, nfen, 3, FileFmtVersion 
      case("Turb3D")
          write(lun, 1020) ndset3dgt, nodes, dtdp*nspo3dgt,
     &                     nspo3dgt, nfen, 3, FileFmtVersion 
      case("Dens2D")
          write(ScreenUnit,*) 
     &         'ERROR: 2D density output not implemented.'
          stop
      case("Dens3D")
          write(lun, 1020) ndset3dgd, nodes, dtdp*nspo3dgd,
     &                     nspo3dgd, nfen, 3, FileFmtVersion 
      case("Press","Wind")
          write(lun, 1010) ndsetsw, nodes, dtdp*nspoolgw,
     &                     nspoolgw, 1, FileFmtVersion
      case("Conc2D")
          write(ScreenUnit,*) 
     &         'ERROR: 2D conc. output not implemented.'
          stop
      case("Conc3D")
          write(ScreenUnit,*) 
     &         'ERROR: 3D conc. output not implemented.'
          stop
      case("Tau0")
          write(lun, 1010) ndsetse, nodes, dtdp*nspoolge,
     &                 nspoolge, 1, FileFmtVersion
      case("sponge")
          write(lun, 1010) 1, nodes, dtdp*nspoolge,
     &                 nspoolge, 1, FileFmtVersion
      case default
          write(ScreenUnit,*) 
     &         'ERROR: Full domain header type unrecongnized.'
      end select
C
C     close the file
      close(lun)

 1000 FORMAT(1X,A32,2X,A24,2X,A24)
 1010 FORMAT(1X,I10,1X,I10,1X,E15.7E3,1X,I5,1X,I5,1X,
     &      'FileFmtVersion: ',I10)
 1020 FORMAT(1X,I10,1X,I10,1X,E15.7E3,1X,I5,1X,I5,1X,I5,1X,
     &      'FileFmtVersion: ',I10)
      !--------------------------------------------------------------
      end subroutine writeDomainHeader
      !--------------------------------------------------------------

C----------------------------------------------------------------------
C       S U B R O U T I N E   R E A D  M I N  M A X 
C----------------------------------------------------------------------
C
C     jgf48.4636 Subroutine to read in continuous min and max files,
C     if they are present. 
C
C----------------------------------------------------------------------
      SUBROUTINE readMinMax(minMaxFN,array)
      USE SIZES , ONLY : SZ, MNPROC, GLOBALDIR
      USE GLOBAL , ONLY : NSCREEN, ScreenUnit, MyProc
      IMPLICIT NONE
      CHARACTER(*) :: minMaxFN ! name of min/max file to read
      REAL(SZ), dimension(:), target :: array ! min/max array
      INTEGER node               ! node number of the min/max data
      CHARACTER(len=80) skipline ! dummy variable for min/max header data
      INTEGER rawNode            ! node number, as read from file
      REAL(SZ) rawDatum          ! array member value at node, as read from file
      INTEGER numLinesInMinMaxFile ! count the lines to report to log file
      LOGICAL tooFewMinMaxLines  ! true if couldn't read enough lines from file
      LOGICAL FileFound  ! true if the file exists
      INTEGER ErrorIO      ! zero if file opened successfully
      INTEGER I
C
      FileFound=.False.
      INQUIRE(FILE=TRIM(GLOBALDIR)//'/'//trim(minMaxFN),
     &   EXIST=FileFound)

      IF (.not.FileFound) THEN
         WRITE(16,2331) trim(minMaxFN)  ! INFO: min/max file  
         WRITE(16,1011)           ! was not found.
         WRITE(16,3333) trim(minMaxFN)  ! WARNING: values not preserved
         IF (NScreen.ne.0.and.MyProc.eq.0) THEN
            WRITE(ScreenUnit,2331) trim(minMaxFN)
            WRITE(ScreenUnit,1011)
            WRITE(ScreenUnit,3333) trim(minMaxFN)! WARNING: values not preserved
         ENDIF
         RETURN  ! early return
      ENDIF
C
C     Now open the file and use the contents to initialize the 
C     min or max array.
      WRITE(16,440) minMaxFN
 440  FORMAT(/,9X,'Min/max file ', (A),' was found. Opening file.')
      OPEN(963,FILE=TRIM(GLOBALDIR)//'/'//trim(minMaxFN),
     &   STATUS='OLD',ACCESS='SEQUENTIAL',IOSTAT=ErrorIO)
      IF (ErrorIO.GT.0) THEN
         WRITE(16,3001) trim(minMaxFN) ! ERROR: min/max file
         WRITE(16,1005)         ! exists but can't be opened
         IF (NScreen.ne.0.and.MyProc.eq.0) THEN
            WRITE(ScreenUnit,3001) trim(minMaxFN)
            WRITE(ScreenUnit,1005)
         ENDIF
         RETURN  ! early return
      ENDIF
C 
C     The file exists and is open, so read the data.
      numLinesInMinMaxFile = 0
      READ(963,*,END=3990,ERR=3990) skipline   ! 1st header line
      READ(963,*,END=3990,ERR=3990) skipline   ! 2nd header line
      READ(963,*,END=3990,ERR=3990) skipline   ! TIME and IT line
      DO I=1, size(array)
         READ(963,*,END=3990,ERR=3990) rawNode, rawDatum 
         IF ((rawNode.ge.1).and.(rawNode.le.size(array))) THEN
            node = rawNode
            array(node) = rawDatum
         ELSE 
            EXIT ! jump out of loop if node number is out of bounds
         ENDIF
         numLinesInMinMaxFile = numLinesInMinMaxFile + 1
      ENDDO
 3990 CONTINUE ! jump to here when end of file is reached
      WRITE(16,450) numLinesInMinMaxFile, trim(minMaxFN)
      IF (numLinesInMinMaxFile.lt.size(array)) THEN
         WRITE(16,*) "ERROR: Not enough data in ",trim(minMaxFN),"."
      ENDIF
 450  FORMAT(/,'Finished reading ',(I8),' lines from the ',
     &   (A12),' file.')
      CLOSE(963)
      RETURN
 1011 FORMAT('was not found.')
 1005 FORMAT('exists but cannot be opened.')
 2331 FORMAT('INFO: The full domain min/max file ', (A12))
 3001 FORMAT('ERROR: The full domain min/max file ', (A12))
 3333 FORMAT('WARNING: Values from ',(A12),
     &       ' will not reflect the solution prior to this hotstart.')
C-----------------------------------------------------------------------
      END SUBROUTINE readMinMax
C-----------------------------------------------------------------------




      !--------------------------------------------------------------
      end module global_io
      !--------------------------------------------------------------
