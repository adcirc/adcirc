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
      program radedit

      implicit none
      integer :: ni, nj, in, out, nsteps, ios

      integer :: iargc, i, j, k, istep
      character(80) :: fnameIn, fnameOut, buf
      character(40) :: buf

      real(8), allocatable :: rs(:,:,:)

      real(8) :: x , y

      if (iargc() < 2) then
        print *, "Usage: radedit file num_steps"
        stop
      endif


      call getarg(1, fnameIn)
      call getarg(2, buf)

      read(buf,'(i10)') nsteps

      print *, "nsteps:", nsteps

      in = 10
      open(in, file=fnameIn)

      fnameOut = trim(fnameIn)//".new"

      out = 11
      open(out, file= fnameOut)

      
      read(in,*) ni, nj, x , y
      read(in,*) buf

      write(out,*) ni, nj, x, y
      write(out,*) buf

      allocate (rs(2, ni, nj))

      do istep = 1, nsteps

        do j = nj, 1, -1
          read(in,*, iostat=ios) ((rs(k,i,j), k= 1,2), i = 1,ni)
        end do

        do j = nj, 1, -1
          write(out,1000) ((rs(k,i,j), k= 1,2), i = 1,ni)
        end do

      end do
 1000 format(1x,1p5e15.7)
      end program
