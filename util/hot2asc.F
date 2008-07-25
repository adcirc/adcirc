      program hot2asc
      implicit none
      integer np, ne, i, irec, iargc, kk
      integer np_g, np_a, ne_g, ne_a
      integer imhs, iths, fileVersion
      integer,parameter :: IN = 10, OUT=11, NVARS=18
      character(80) :: inbuf, title
      character(7) :: binfname      
      character(11) :: ascfname      
      real(8) :: time

      
      character(8) :: vars(NVARS) =
     &  (/ "IESTP   ", "NSCOUE  ", "IVSTP   ", "NSCOUV  ", "ICSTP   ", 
     &     "NSCOUC  ",
     &     "IPSTP   ", "IWSTP   ", "NSCOUM  ", "IGEP    ", "NSCOUGE ", 
     &     "IGVP    ",
     &     "NSCOUGV ", "IGCP    ", "NSCOUGC ", "IGPP    ", "IGWP    ", 
     &     "NSCOUGW " /)


      if (iargc() < 1) then
        print *, 'Usage hot2asc hotstart_file'
        stop
      endif

      print *, "Version: "//
     $  "$Id: hot2asc.F,v 1.1.2.4 2006-10-19 21:44:16 mclay Exp $"

      call getarg(1,inbuf)
      binfname = trim(inbuf)
      ascfname = binfname // '.asc'
      print *, "binfname = ", binfname
      print *, "ascfname = ", ascfname

      open(in, file=binfname, ACCESS='DIRECT',RECL=8)
      open(out,file=ascfname, FORM='FORMATTED',  ACCESS='SEQUENTIAL')

      irec = 1 
      read(in,rec=irec) fileVersion    ; irec = irec + 1
      write(out, '(1x,''Major: '',i3,'' Minor: '', i3, '' Rev: '',i3)')
     $  ishft(fileVersion,-20), iand(1023,ishft(fileVersion,-10)),
     $  iand(1023,fileVersion)
      read(in,rec=irec) imhs    ; irec = irec + 1
      write(out,'(A,i8)') "imhs = ", imhs
      
      read(in,rec=irec) time    ; irec = irec + 1
      write(out,'(A,e25.16)') "time = ", time
 
      read(in,rec=irec) iths   ; irec = irec + 1
      write(out,'(A,i10)') "iths = ", iths

      read(in,rec=irec) np_g   ; irec = irec + 1
      write(out,'(A,i10)') "NP_G = ", np_g

      read(in,rec=irec) ne_g   ; irec = irec + 1
      write(out,'(A,i10)') "NE_G = ", ne_g

      read(in,rec=irec) np_a   ; irec = irec + 1
      write(out,'(A,i10)') "NP_A = ", np_a

      read(in,rec=irec) ne_a   ; irec = irec + 1
      write(out,'(A,i10)') "NE_A = ", ne_a

      np = np_g
      ne = ne_g
      
      call dsply(in,out,irec,"ETA1",np)
      call dsply(in,out,irec,"ETA2",np)
      call dsply(in,out,irec,"EtaDisc",np)
      call dsply(in,out,irec,"UU2",np)
      call dsply(in,out,irec,"VV2",np)

      if (imhs ==  10) call dsply(in,out,irec,"CH1",np)

      call idsply(in,out,irec,"NODECODE",np)
      call idsply(in,out,irec,"NOFF",ne)
      
      do i = 1, NVARS
        irec = irec + 1
        read(in, rec = irec) kk
        write(out, '(a,''='',i10)') vars(i), kk
      end do

      print *, "final irec = ", irec

      close(in)
      close(out)
      stop
      end

      subroutine dsply(in,out, irec, varname, size)

      implicit none
      integer :: in, out, irec, size, i
      real (8) :: x
      character(*) :: varname

      write(out, 1000) trim(varname)

      do i = 1, size
        !irec = irec + 1
        read(in, rec=irec) x ; irec = irec + 1
        write(out, 1010) trim(varname), i, x
      end do

      return
 1000 format("#---- ",a," ----")
 1010 format(a10, i8, ":", 1pe20.10)
      end subroutine


      subroutine idsply(in,out, irec, varname, size)

      implicit none
      integer :: in, out, irec, size, i, x
      character(*) :: varname

      write(out, 1000) trim(varname)

      do i = 1, size
        !irec = irec + 1
        read(in, rec=irec) x ; irec = irec + 1
        write(out, 1010) trim(varname), i, x
      end do

      return
 1000 format("#---- ",a," ----")
 1010 format(a10, i8, ":", I10)
      end subroutine
            
