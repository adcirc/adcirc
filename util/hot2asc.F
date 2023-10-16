      program hot2asc
      implicit none
      integer np, ne, i, j, n, irec, iargc, kk
      integer np_g, np_a, ne_g, ne_a
      integer imhs, iths, fileVersion
      integer nf
      integer nfreq
      integer mm
      integer icha
      integer nstae ! number of elevation stations for harmonic analysis
      integer nstav ! number of velocity station for harmonic analysis
      integer npha  ! number of nodes, from harmonic analysis section
      integer nhage ! full domain harmonic analysis on elevation (on/off)
      integer nhagv ! full domain harmonic analysis on velocity (on/off)
      integer nhase ! station harmonic analysis on elevation (on/off)
      integer nhasv ! station harmonic analysis on velocity (on/off)
      CHARACTER(8) FNAM8(2)
      integer,parameter :: IN = 10, OUT=11, NVARS=18, NHSVARS=12
      character(80) :: arg ! holds command line arguments for processing 
      character(80) :: title
      character(7) :: binfname      
      character(11) :: ascfname      
      character(2) :: output3D
      logical :: get3D 
      logical :: getHarmonic
      logical :: getReSynth  ! true if time series resynthesis data requested
      real(8) :: time
      integer :: itud
      real(8) :: realvar ! assume that reals are 8 bytes (very likely)
      logical :: fileFound  ! true if the requested file was found
      integer :: errorIO ! io status variable 
      
      character(8) :: vars(NVARS) =
     &  (/ "IESTP   ", "NSCOUE  ", "IVSTP   ", "NSCOUV  ", "ICSTP   ", 
     &     "NSCOUC  ",
     &     "IPSTP   ", "IWSTP   ", "NSCOUM  ", "IGEP    ", "NSCOUGE ", 
     &     "IGVP    ",
     &     "NSCOUGV ", "IGCP    ", "NSCOUGC ", "IGPP    ", "IGWP    ", 
     &     "NSCOUGW " /)

      character(8) :: hsvars(NHSVARS) =
     &  (/ "NZ      ", "NF      ", "MM      ", "NP      ", 
     &     "NSTAE   ",
     &     "NSTAV   ", "NHASE   ", "NHASV   ", "NHAGE   ", "NHAGV   ", 
     &     "ICALL   ",
     &     "NFREQ   " /) 

      print *, "Version: v48.xx"

      i = 1
      if (iargc() < i) then
         print *, 'Usage hot2asc hotstart_filename [harmonic] [resynth]'
         stop
      else
        call getarg(i,arg) ! filename
      endif

      binfname = trim(arg)
      ascfname = binfname // '.asc'
      print *, "binfname = ", binfname
      print *, "ascfname = ", ascfname

      fileFound = .false.
      inquire(file=binfname,exist=fileFound)
      if (.not.fileFound) then
         print *,'ERROR: A file named ',binfname,' was not found.'
         print *,'Please check the file name and try again. Stopping.'
         stop
      else 
         print *,'File named ',binfname,' was found. Opening file.'
      endif

      open(AMPI_LUN(in), file=binfname, ACCESS='DIRECT', RECL=8, IOSTAT=errorIO)
      call AMPI_LUN_REGISTER(AMPI_LUN(in), file=binfname, ACCESS='DIRECT', RECL=8, IOSTAT=errorIO)
      if ( errorIO.gt.0) then
         print *,'ERROR: The file ',binfname,' cannot be opened.'
         print *,'Stopping.'
      endif

      open(AMPI_LUN(out),file=ascfname, FORM='FORMATTED', ACCESS='SEQUENTIAL',
     &    iostat=errorIO)
      call AMPI_LUN_REGISTER(AMPI_LUN(out),file=ascfname, FORM='FORMATTED', ACCESS='SEQUENTIAL',
     &    iostat=errorIO)
      if ( errorIO.gt.0) then
         print *,'ERROR: The file ',binfname,' cannot be opened.'
         print *,'Stopping.'
      endif

      !jgf48.38 Make it possible to extract harmonic and time series 
      ! resynthesis data but user must decide whether they are present 
      ! or not ... hot start files are not self describing
      get3D = .false.
      getHarmonic = .false.
      getReSynth = .false.
 
      i = i + 1
      do j = i, iargc()
         call getarg(j,arg)
         select case (arg)
         case("harmonic","Harmonic","HARMONIC") 
            getHarmonic = .true.
            print *,'Harmonic data were requested.'
            print *,'They will be retrieved, if present.'
         case("resynth","ReSynth","RESYNTH")
            getReSynth = .true.
            print *,'Time series resynthesis data were requested.'
            print *,'They will be retrieved, if present.' 
         case default
            print *,'Command line argument not understood:'
            print *,arg
         end select 
      end do 

      irec = 1 
      read(AMPI_LUN(in),rec=irec) fileVersion    ; irec = irec + 1
      write(AMPI_LUN(out), '(1x,''Major: '',i3,'' Minor: '', i3, '' Rev: '',i3)')
     $ ishft(fileVersion,-20), iand(1023,ishft(fileVersion,-10)),
     $ iand(1023,fileVersion)
      read(AMPI_LUN(in),rec=irec) imhs    ; irec = irec + 1
      write(AMPI_LUN(out),'(A,i8)') "imhs = ", imhs

      !jgf48.38 need to support 3D data as well.
      select case (imhs)
         case(1,2,11,21,31)
            get3D = .true.
            print *,'File contains 3D data, this is not supported.'
            if (getHarmonic) then
               getHarmonic = .false.
               print *,'Harmonic data cannot currently be converted'
               print *,'when 3D data are present.'
            endif
         case default
            print *,'File does not contain 3D data.'
      end select
      
      read(AMPI_LUN(in),rec=irec) time    ; irec = irec + 1
      write(AMPI_LUN(out),'(A,e25.16)') "time = ", time
 
      read(AMPI_LUN(in),rec=irec) iths   ; irec = irec + 1
      write(AMPI_LUN(out),'(A,i10)') "iths = ", iths

      read(AMPI_LUN(in),rec=irec) np_g   ; irec = irec + 1
      write(AMPI_LUN(out),'(A,i10)') "NP_G = ", np_g

      read(AMPI_LUN(in),rec=irec) ne_g   ; irec = irec + 1
      write(AMPI_LUN(out),'(A,i10)') "NE_G = ", ne_g

      read(AMPI_LUN(in),rec=irec) np_a   ; irec = irec + 1
      write(AMPI_LUN(out),'(A,i10)') "NP_A = ", np_a

      read(AMPI_LUN(in),rec=irec) ne_a   ; irec = irec + 1
      write(AMPI_LUN(out),'(A,i10)') "NE_A = ", ne_a

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
        read(AMPI_LUN(in), rec = irec) kk
        irec = irec + 1
        write(AMPI_LUN(out), '(a,''='',i10)') vars(i), kk
      end do

      ! jgf48.38 add support for reading harmonic analysis data
      if (getHarmonic) then
         read(AMPI_LUN(in), rec = irec+1) icha
         irec = irec + 1
         do i = 1, NHSVARS
            read(AMPI_LUN(in), rec = irec+i) kk
            write(AMPI_LUN(out), '(a,''='',i10)') hsvars(i), kk
            ! save the values of certain parameters for later use
            select case(hsvars(i))
               case("NF")    ; nf = kk
               case("NFREQ") ; nfreq = kk
               case("MM")    ; mm = kk
               case("NSTAE") ; nstae = kk
               case("NSTAV") ; nstav = kk
               case("NP")    ; npha = kk
               case("NHAGE") ; nhage = kk
               case("NHAGV") ; nhagv = kk
               case("NHASE") ; nhase = kk
               case("NHASV") ; nhasv = kk
               case default
                  ! do nothing; not a variable we need to record for later
            end select            
         end do
         irec = irec + NHSVARS
         do i=1,nfreq+nf
            read(AMPI_LUN(in), rec=irec+1) fnam8(1) 
            write(AMPI_LUN(out),'(A,a)') "fnam8(1) = ", fnam8(1)  
            read(AMPI_LUN(in), rec=irec+2) fnam8(2) 
            write(AMPI_LUN(out),'(A,a)') "fnam8(2) = ", fnam8(2) 
            irec = irec + 2
            read(AMPI_LUN(in),rec=irec+1) realvar    
            write(AMPI_LUN(out),'(A,i2,A,e25.16)') "hafreq(",i,") = ", realvar
            read(AMPI_LUN(in),rec=irec+2) realvar    
            write(AMPI_LUN(out),'(A,i2,A,e25.16)') "haff(",i,")   = ", realvar
            read(AMPI_LUN(in),rec=irec+3) realvar   
            write(AMPI_LUN(out),'(A,i2,A,e25.16)') "haface(",i,") = ", realvar
            irec = irec + 3
         end do

         read(AMPI_LUN(in),rec=irec+1) realvar      
         write(AMPI_LUN(out),'(A,e25.16)') "timeud = ", realvar
         read(AMPI_LUN(in),rec=irec+2) itud         
         write(AMPI_LUN(out),'(A,i10)') "itud = ", itud
         !jgf48.38 this must increment 3 rather than 2 as you might
         ! expect, because itud is a 4 byte integer, and when the 
         ! read statement steps forward by "1" record, it is really
         ! only stepping forward by 4 bytes, so we need to step 
         ! forward another record before reading again
         irec = irec + 3

         do i=1,mm
            do j=1,mm
               read(AMPI_LUN(in),rec=irec) realvar    ; irec = irec + 1
               write(AMPI_LUN(out),'(A,i2,A,i2,A,e25.16)') 
     &               "ha(",i,",",j,") = ", realvar
            end do
         end do

         if ( nhase.eq.1 ) then
            do n=1,nstae
               do i=1,mm
                  read(AMPI_LUN(in), rec=irec) realvar 
                  write(AMPI_LUN(out),'(A,i2,A,i2,A,e25.16)') 
     &                   "staelv(",i,",",n,") = ", realvar
                  irec = irec + 1
               end do
            end do
         endif

         if ( nhasv.eq.1) then
            do n=1,nstav
               do i=1,mm
                  read(AMPI_LUN(in), rec=irec) realvar 
                  write(AMPI_LUN(out),'(A,i2,A,i2,A,e25.16)') 
     &                "staulv(",i,",",n,") = ", realvar
                  irec = irec + 1
                  read(AMPI_LUN(in), rec=irec) realvar 
                  write(AMPI_LUN(out),'(A,i2,A,i2,A,e25.16)') 
     &                "stavlv(",i,",",n,") = ", realvar
                  irec = irec + 1
               end do
            end do
         endif        

         if ( nhage.eq.1 ) then
            do n=1,npha
               do i=1,mm
                  read(AMPI_LUN(in), rec=irec) realvar 
                  write(AMPI_LUN(out),'(A,i2,A,i2,A,e25.16)') 
     &                   "gloelv(",i,",",n,") = ", realvar
                  irec = irec + 1
               end do
            end do
         endif

         if ( nhagv.eq.1) then
            do n=1,npha
               do i=1,mm
                  read(AMPI_LUN(in), rec=irec) realvar 
                  write(AMPI_LUN(out),'(A,i2,A,i2,A,e25.16)') 
     &                "gloulv(",i,",",n,") = ", realvar
                  irec = irec + 1
                  read(AMPI_LUN(in), rec=irec) realvar 
                  write(AMPI_LUN(out),'(A,i2,A,i2,A,e25.16)') 
     &                "glovlv(",i,",",n,") = ", realvar
                  irec = irec + 1
               end do
            end do
         endif        

      endif

      if (getReSynth) then
         read(AMPI_LUN(in), rec=irec) kk ; irec = irec + 1
         write(AMPI_LUN(out),'(A,i10)') "NE_A = ", kk
         if (nhage.eq.1) then
            do i = 1, npha
               read(AMPI_LUN(in), rec=irec) realvar ; irec = irec + 1
               write(AMPI_LUN(out),'(A,i10,A,e25.16)') 
     &           "ELAV(",i,") = ", realvar
               read(AMPI_LUN(in), rec=irec) realvar ; irec = irec + 1
               write(AMPI_LUN(out),'(A,i10,A,e25.16)') 
     &           "ELVA(",i,") = ", realvar
            end do
         endif   
         if ( nhagv.eq.1 ) then
            do i = 1, np
               read(AMPI_LUN(in), rec=irec) realvar ; irec = irec + 1
               write(AMPI_LUN(out),'(A,i10,A,e25.16)') 
     &           "XVELAV(",i,") = ", realvar
               read(AMPI_LUN(in), rec=irec) realvar ; irec = irec + 1
               write(AMPI_LUN(out),'(A,i10,A,e25.16)') 
     &           "YVELAV(",i,") = ", realvar
               read(AMPI_LUN(in), rec=irec) realvar ; irec = irec + 1
               write(AMPI_LUN(out),'(A,i10,A,e25.16)') 
     &           "XVELVA(",i,") = ", realvar
               read(AMPI_LUN(in), rec=irec) realvar ; irec = irec + 1
               write(AMPI_LUN(out),'(A,i10,A,e25.16)') 
     &           "YVELVA(",i,") = ", realvar
            end do
         endif
      endif 

      print *, "final irec = ", irec

      close(AMPI_LUN(in))
      close(AMPI_LUN(out))
      stop
      end

      subroutine dsply(in,out, irec, varname, size)

      implicit none
      integer :: in, out, irec, size, i
      real (8) :: x
      character(*) :: varname

      write(AMPI_LUN(out), 1000) trim(varname)

      do i = 1, size
        !irec = irec + 1
        read(AMPI_LUN(in), rec=irec) x ; irec = irec + 1
        write(AMPI_LUN(out), 1010) trim(varname), i, x
      end do

      return
 1000 format("#---- ",a," ----")
 1010 format(a10, i8, ":", 1pe20.10)
      end subroutine


      subroutine idsply(in,out, irec, varname, size)

      implicit none
      integer :: in, out, irec, size, i, x
      character(*) :: varname

      write(AMPI_LUN(out), 1000) trim(varname)

      do i = 1, size
        !irec = irec + 1
        read(AMPI_LUN(in), rec=irec) x ; irec = irec + 1
        write(AMPI_LUN(out), 1010) trim(varname), i, x
      end do

      return
 1000 format("#---- ",a," ----")
 1010 format(a10, i8, ":", I10)
      end subroutine
            
