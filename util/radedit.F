      program radedit

      implicit none
      integer :: ni, nj, in, out, nsteps, ios
!$omp threadprivate(ni,nj,in,out,nsteps,ios)

      integer :: iargc, i, j, k, istep
!$omp threadprivate(iargc,i,j,k,istep)
      character(80) :: fnameIn, fnameOut, buf
!$omp threadprivate(fnameIn,fnameOut,buf)
      character(40) :: buf
!$omp threadprivate(buf)

      real(8), allocatable :: rs(:,:,:)
!$omp threadprivate(rs)

      real(8) :: x , y
!$omp threadprivate(x,y)

      if (iargc() < 2) then
        print *, "Usage: radedit file num_steps"
        stop
      endif


      call getarg(1, fnameIn)
      call getarg(2, buf)

      read(buf,'(i10)') nsteps

      print *, "nsteps:", nsteps

      in = 10
      open(in+CK_LUN, file=fnameIn)

      fnameOut = trim(fnameIn)//".new"

      out = 11
      open(out+CK_LUN, file= fnameOut)

      
      read(in+CK_LUN,*) ni, nj, x , y
      read(in+CK_LUN,*) buf

      write(out+CK_LUN,*) ni, nj, x, y
      write(out+CK_LUN,*) buf

      allocate (rs(2, ni, nj))

      do istep = 1, nsteps

        do j = nj, 1, -1
          read(in+CK_LUN,*, iostat=ios) ((rs(k,i,j), k= 1,2), i = 1,ni)
        end do

        do j = nj, 1, -1
          write(out+CK_LUN,1000) ((rs(k,i,j), k= 1,2), i = 1,ni)
        end do

      end do
 1000 format(1x,1p5e15.7)
      end program
