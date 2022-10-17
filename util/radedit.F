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
      open(AMPI_LUN(in), file=fnameIn)
      call AMPI_LUN_REGISTER(AMPI_LUN(in), file=fnameIn)

      fnameOut = trim(fnameIn)//".new"

      out = 11
      open(AMPI_LUN(out), file= fnameOut)
      call AMPI_LUN_REGISTER(AMPI_LUN(out), file= fnameOut)

      
      read(AMPI_LUN(in),*) ni, nj, x , y
      read(AMPI_LUN(in),*) buf

      write(AMPI_LUN(out),*) ni, nj, x, y
      write(AMPI_LUN(out),*) buf

      allocate (rs(2, ni, nj))

      do istep = 1, nsteps

        do j = nj, 1, -1
          read(AMPI_LUN(in),*, iostat=ios) ((rs(k,i,j), k= 1,2), i = 1,ni)
        end do

        do j = nj, 1, -1
          write(AMPI_LUN(out),1000) ((rs(k,i,j), k= 1,2), i = 1,ni)
        end do

      end do
 1000 format(1x,1p5e15.7)
      end program
