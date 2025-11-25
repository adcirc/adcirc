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
      program inflate
!
! Use this code to re-inflate the sparse format of the 
! special 63/64/73/74 format from FEMA/LACPR project.
!
! Two arguments required: input and output file names.
! One optional argument
!
!
      implicit none
      character*80 line
      character*256 filefrom, fileto, buf
      character*256 chunkify   ! used to break up output file into chunks
      integer numSnapsPerChunk  ! chunk size
      integer chunkNum             ! chunk number being output
      character(len=8) chunkFileName ! name of output chunk file 
      logical chunksOut
      integer snapNum  ! number of snaps written to the current chunk
      integer i,j,k,nspoolge,iter
      integer nnum,np,nodes,nsnaps
      real*8 x,y,time,dtdp,default
      real*8,allocatable :: u(:), v(:)
      integer in, out 
      integer iargc

      chunksOut = .false.

      if(COMMAND_ARGUMENT_COUNT().lt.2) then
         print *, 
     &   'Usage: inflate_xxxx <input file> <output file> chunkify'
         print *, "  <input file> can be 'STDIN'"
         print *, "  <output file> can be 'STDOUT'"
         print *, '  chunkify is optional, and will generate'
         print *, '  an output file in chunks'
         stop
      endif

      call GET_COMMAND_ARGUMENT(1, filefrom)
      call GET_COMMAND_ARGUMENT(2, fileto)
      if ( COMMAND_ARGUMENT_COUNT().eq.3 ) call GET_COMMAND_ARGUMENT(3,chunkify)

      if (chunkify.eq."chunkify") chunksOut = .true.

      if(filefrom.eq."STDIN") then
         in = 5
      else
         in = 1
         open(in,file=filefrom,action="read")
      endif

      if(fileto.eq."STDOUT") then
         out = 6
      else
         out = 2
         if (chunksOut) then
             chunkNum = 10
             chunkFileName = "chunk.10"
             open(out,file=chunkFileName,action="write")
         else
             open(out,file=fileto,action="write")
         endif
      endif

      read(in,'(a80)') line
      write(out,10) line
      read(in,*) nsnaps,nodes,dtdp,nspoolge,k
C
C     jgf47.05 Set nsnaps to the actual number of datasets in the 
C     file, so that SMS can use it.
      call countsnaps(nsnaps, in)
C
C     jgf47.07 Determine how many datasets will appear in each chunk
      if (chunksOut) then
         numSnapsPerChunk = nsnaps / 10;
      endif      
C
C     Begin writing output file
      write(out,15) nsnaps,nodes,dtdp,nspoolge,k

      allocate(u(nodes))
      allocate(v(nodes))
C
C     jgf47.03 Changed to infinite loop b/c nsnaps is not
C     a reliable indicator of the number of datasets in the
C     file, b/c it can't be updated during a hotstart when the
C     file is appended
      i = 0
      snapNum = 1 
      if(k.eq.1) then
         do 
            read(in,*,end=95) time, iter, np, default
            if (chunksOut) then
                if (snapNum.gt.numSnapsPerChunk) then
                    close(out)
                    chunkNum = chunkNum + 1
                    write(chunkFileName(7:8),2) chunkNum
                    open(out,file=chunkFileName,action="write")
                    snapNum = 1
                endif
            endif
            write(out,20) time,iter
            do j=1,nodes
               u(j) = default
            enddo
            do j=1,np
               read(in,*) nnum,x
               u(nnum) = x
            enddo
            do j=1,nodes
               write(buf,30) j, u(j)
               write(out,*) trim(adjustl(buf))
            end do
            i = i + 1
            snapNum = snapNum + 1
         enddo
      else
         do 
            read(in,*,end=95) time, iter, np, default
            if (chunksOut) then
                if (snapNum.gt.numSnapsPerChunk) then
                    close(out)
                    chunkNum = chunkNum + 1
                    write(chunkFileName(7:8),2) chunkNum
                    open(out,file=chunkFileName,action="write")
                    snapNum = 1
                endif
            endif
            write(out,20) time,iter
            do j=1,nodes
               u(j) = default
               v(j) = default
            enddo
            do j=1,np
               read(in,*) nnum,x,y
               u(nnum) = x
               v(nnum) = y
            enddo
            do j=1,nodes
               write(buf,40) j, u(j), v(j)
               write(out,*) trim(adjustl(buf))
            end do
            i = i + 1
         enddo
      endif

95    continue
      write(0,*) 'PROCESSED ', i,' datasets'
      deallocate(u,v)

      close(1)
      close(2)

  2   format(i2)
 10   format(a80) 
 15   format(1X,I10,1X,I10,1X,E15.7,1X,I5,1X,I5)
 20   format(2X,E20.10,5X,I10)
 30   format(i8,1x,e15.8)
 40   format(i8,1x,e13.6,1x,e13.6)
C     ----------------------------------------------------------------
      end program inflate
C     ----------------------------------------------------------------



C     ----------------------------------------------------------------
C         S U B R O U T I N E   C O U N T S N A P S 
C     ----------------------------------------------------------------
C     jgf47.05 Subroutine to count the number of datasets in the file. 
C     this is needed because SMS relies on this number if it is positive,
C     and this number will be wrong after a hotstart. Also, it can be
C     negative if there was a hotstart before the output is scheduled
C     to begin. 
C     ----------------------------------------------------------------
      subroutine countsnaps(nsnaps, in)
      implicit none
      integer, intent(out) :: nsnaps ! number of datasets in the file
      integer, intent(in) :: in      ! unit number of input file
      character(len=80) skipped      ! data that must be skipped
      real*8 skiptime                ! time (seconds), not needed
      integer skipiter               ! time step number, not needed
      integer np                     ! number of nodes in this dataset
      real*8 skipdefault             ! value of unwritten nodes, not needed
      integer j                      ! counter for nodes in dataset

      nsnaps = 0
      do 
         read(in,*,end=195) skiptime, skipiter, np, skipdefault
         do j=1,np
            read(in,*) skipped
         enddo
         nsnaps = nsnaps + 1
      enddo

195   write(0,*) 'Finished counting. File contains ', nsnaps, 
     &           ' datasets.'
      rewind(in)
      read(in,*) skipped
      read(in,*) skipped
      return
C     ----------------------------------------------------------------
      end subroutine countsnaps
C     ----------------------------------------------------------------
