      program hot2asc
      implicit none
      integer np, ne, i, j, n, irec, iargc, kk
!$omp threadprivate(np,ne,i,j,n,irec,iargc,kk)
      integer np_g, np_a, ne_g, ne_a
!$omp threadprivate(np_g,np_a,ne_g,ne_a)
      integer imhs, iths, fileVersion
!$omp threadprivate(imhs,iths,fileVersion)
      integer nf
!$omp threadprivate(nf)
      integer nfreq
!$omp threadprivate(nfreq)
      integer mm
!$omp threadprivate(mm)
      integer icha
!$omp threadprivate(icha)
      integer nstae ! number of elevation stations for harmonic analysis
!$omp threadprivate(nstae)
      integer nstav ! number of velocity station for harmonic analysis
!$omp threadprivate(nstav)
      integer npha  ! number of nodes, from harmonic analysis section
!$omp threadprivate(npha)
      integer nhage ! full domain harmonic analysis on elevation (on/off)
!$omp threadprivate(nhage)
      integer nhagv ! full domain harmonic analysis on velocity (on/off)
!$omp threadprivate(nhagv)
      integer nhase ! station harmonic analysis on elevation (on/off)
!$omp threadprivate(nhase)
      integer nhasv ! station harmonic analysis on velocity (on/off)
!$omp threadprivate(nhasv)
      CHARACTER(8) FNAM8(2)
!$omp threadprivate(FNAM8)
      integer,parameter :: IN = 10, OUT=11, NVARS=18, NHSVARS=12
      character(80) :: arg ! holds command line arguments for processing 
!$omp threadprivate(arg)
      character(80) :: title
!$omp threadprivate(title)
      character(7) :: binfname      
!$omp threadprivate(binfname)
      character(11) :: ascfname      
!$omp threadprivate(ascfname)
      character(2) :: output3D
!$omp threadprivate(output3D)
      logical :: get3D 
!$omp threadprivate(get3D)
      logical :: getHarmonic
!$omp threadprivate(getHarmonic)
      logical :: getReSynth  ! true if time series resynthesis data requested
!$omp threadprivate(getReSynth)
      real(8) :: time
!$omp threadprivate(time)
      integer :: itud
!$omp threadprivate(itud)
      real(8) :: realvar ! assume that reals are 8 bytes (very likely)
!$omp threadprivate(realvar)
      logical :: fileFound  ! true if the requested file was found
!$omp threadprivate(fileFound)
      integer :: errorIO ! io status variable 
!$omp threadprivate(errorIO)
      
      character(8) :: vars(NVARS) =
     &  (/ "IESTP   ", "NSCOUE  ", "IVSTP   ", "NSCOUV  ", "ICSTP   ", 
     &     "NSCOUC  ",
     &     "IPSTP   ", "IWSTP   ", "NSCOUM  ", "IGEP    ", "NSCOUGE ", 
     &     "IGVP    ",
     &     "NSCOUGV ", "IGCP    ", "NSCOUGC ", "IGPP    ", "IGWP    ", 
     &     "NSCOUGW " /)
!$omp threadprivate(vars)

      character(8) :: hsvars(NHSVARS) =
     &  (/ "NZ      ", "NF      ", "MM      ", "NP      ", 
     &     "NSTAE   ",
     &     "NSTAV   ", "NHASE   ", "NHASV   ", "NHAGE   ", "NHAGV   ", 
     &     "ICALL   ",
     &     "NFREQ   " /) 
!$omp threadprivate(hsvars)

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

      open(in+CK_LUN, file=binfname, ACCESS='DIRECT', RECL=8, IOSTAT=errorIO)
      if ( errorIO.gt.0) then
         print *,'ERROR: The file ',binfname,' cannot be opened.'
         print *,'Stopping.'
      endif

      open(out+CK_LUN,file=ascfname, FORM='FORMATTED', ACCESS='SEQUENTIAL',
     &     iostat=errorIO)
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
      read(in+CK_LUN,rec=irec) fileVersion    ; irec = irec + 1
      write(out+CK_LUN, '(1x,''Major: '',i3,'' Minor: '', i3, '' Rev: '',i3)')
     $  ishft(fileVersion,-20), iand(1023,ishft(fileVersion,-10)),
     $  iand(1023,fileVersion)
      read(in+CK_LUN,rec=irec) imhs    ; irec = irec + 1
      write(out+CK_LUN,'(A,i8)') "imhs = ", imhs

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
      
      read(in+CK_LUN,rec=irec) time    ; irec = irec + 1
      write(out+CK_LUN,'(A,e25.16)') "time = ", time
 
      read(in+CK_LUN,rec=irec) iths   ; irec = irec + 1
      write(out+CK_LUN,'(A,i10)') "iths = ", iths

      read(in+CK_LUN,rec=irec) np_g   ; irec = irec + 1
      write(out+CK_LUN,'(A,i10)') "NP_G = ", np_g

      read(in+CK_LUN,rec=irec) ne_g   ; irec = irec + 1
      write(out+CK_LUN,'(A,i10)') "NE_G = ", ne_g

      read(in+CK_LUN,rec=irec) np_a   ; irec = irec + 1
      write(out+CK_LUN,'(A,i10)') "NP_A = ", np_a

      read(in+CK_LUN,rec=irec) ne_a   ; irec = irec + 1
      write(out+CK_LUN,'(A,i10)') "NE_A = ", ne_a

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
        read(in+CK_LUN, rec = irec) kk
        irec = irec + 1
        write(out+CK_LUN, '(a,''='',i10)') vars(i), kk
      end do

      ! jgf48.38 add support for reading harmonic analysis data
      if (getHarmonic) then
         read(in+CK_LUN, rec = irec+1) icha
         irec = irec + 1
         do i = 1, NHSVARS
            read(in+CK_LUN, rec = irec+i) kk
            write(out+CK_LUN, '(a,''='',i10)') hsvars(i), kk
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
            read(in+CK_LUN, rec=irec+1) fnam8(1) 
            write(out+CK_LUN,'(A,a)') "fnam8(1) = ", fnam8(1)  
            read(in+CK_LUN, rec=irec+2) fnam8(2) 
            write(out+CK_LUN,'(A,a)') "fnam8(2) = ", fnam8(2) 
            irec = irec + 2
            read(in+CK_LUN,rec=irec+1) realvar    
            write(out+CK_LUN,'(A,i2,A,e25.16)') "hafreq(",i,") = ", realvar
            read(in+CK_LUN,rec=irec+2) realvar    
            write(out+CK_LUN,'(A,i2,A,e25.16)') "haff(",i,")   = ", realvar
            read(in+CK_LUN,rec=irec+3) realvar   
            write(out+CK_LUN,'(A,i2,A,e25.16)') "haface(",i,") = ", realvar
            irec = irec + 3
         end do

         read(in+CK_LUN,rec=irec+1) realvar      
         write(out+CK_LUN,'(A,e25.16)') "timeud = ", realvar
         read(in+CK_LUN,rec=irec+2) itud         
         write(out+CK_LUN,'(A,i10)') "itud = ", itud
         !jgf48.38 this must increment 3 rather than 2 as you might
         ! expect, because itud is a 4 byte integer, and when the 
         ! read statement steps forward by "1" record, it is really
         ! only stepping forward by 4 bytes, so we need to step 
         ! forward another record before reading again
         irec = irec + 3

         do i=1,mm
            do j=1,mm
               read(in+CK_LUN,rec=irec) realvar    ; irec = irec + 1
               write(out+CK_LUN,'(A,i2,A,i2,A,e25.16)') 
     &                "ha(",i,",",j,") = ", realvar
            end do
         end do

         if ( nhase.eq.1 ) then
            do n=1,nstae
               do i=1,mm
                  read(in+CK_LUN, rec=irec) realvar 
                  write(out+CK_LUN,'(A,i2,A,i2,A,e25.16)') 
     &                    "staelv(",i,",",n,") = ", realvar
                  irec = irec + 1
               end do
            end do
         endif

         if ( nhasv.eq.1) then
            do n=1,nstav
               do i=1,mm
                  read(in+CK_LUN, rec=irec) realvar 
                  write(out+CK_LUN,'(A,i2,A,i2,A,e25.16)') 
     &                 "staulv(",i,",",n,") = ", realvar
                  irec = irec + 1
                  read(in+CK_LUN, rec=irec) realvar 
                  write(out+CK_LUN,'(A,i2,A,i2,A,e25.16)') 
     &                 "stavlv(",i,",",n,") = ", realvar
                  irec = irec + 1
               end do
            end do
         endif        

         if ( nhage.eq.1 ) then
            do n=1,npha
               do i=1,mm
                  read(in+CK_LUN, rec=irec) realvar 
                  write(out+CK_LUN,'(A,i2,A,i2,A,e25.16)') 
     &                    "gloelv(",i,",",n,") = ", realvar
                  irec = irec + 1
               end do
            end do
         endif

         if ( nhagv.eq.1) then
            do n=1,npha
               do i=1,mm
                  read(in+CK_LUN, rec=irec) realvar 
                  write(out+CK_LUN,'(A,i2,A,i2,A,e25.16)') 
     &                 "gloulv(",i,",",n,") = ", realvar
                  irec = irec + 1
                  read(in+CK_LUN, rec=irec) realvar 
                  write(out+CK_LUN,'(A,i2,A,i2,A,e25.16)') 
     &                 "glovlv(",i,",",n,") = ", realvar
                  irec = irec + 1
               end do
            end do
         endif        

      endif

      if (getReSynth) then
         read(in+CK_LUN, rec=irec) kk ; irec = irec + 1
         write(out+CK_LUN,'(A,i10)') "NE_A = ", kk
         if (nhage.eq.1) then
            do i = 1, npha
               read(in+CK_LUN, rec=irec) realvar ; irec = irec + 1
               write(out+CK_LUN,'(A,i10,A,e25.16)') 
     &            "ELAV(",i,") = ", realvar
               read(in+CK_LUN, rec=irec) realvar ; irec = irec + 1
               write(out+CK_LUN,'(A,i10,A,e25.16)') 
     &            "ELVA(",i,") = ", realvar
            end do
         endif   
         if ( nhagv.eq.1 ) then
            do i = 1, np
               read(in+CK_LUN, rec=irec) realvar ; irec = irec + 1
               write(out+CK_LUN,'(A,i10,A,e25.16)') 
     &            "XVELAV(",i,") = ", realvar
               read(in+CK_LUN, rec=irec) realvar ; irec = irec + 1
               write(out+CK_LUN,'(A,i10,A,e25.16)') 
     &            "YVELAV(",i,") = ", realvar
               read(in+CK_LUN, rec=irec) realvar ; irec = irec + 1
               write(out+CK_LUN,'(A,i10,A,e25.16)') 
     &            "XVELVA(",i,") = ", realvar
               read(in+CK_LUN, rec=irec) realvar ; irec = irec + 1
               write(out+CK_LUN,'(A,i10,A,e25.16)') 
     &            "YVELVA(",i,") = ", realvar
            end do
         endif
      endif 

      print *, "final irec = ", irec

      close(in+CK_LUN)
      close(out+CK_LUN)
      stop
      end

      subroutine dsply(in,out, irec, varname, size)

      implicit none
      integer :: in, out, irec, size, i
      real (8) :: x
      character(*) :: varname

      write(out+CK_LUN, 1000) trim(varname)

      do i = 1, size
        !irec = irec + 1
        read(in+CK_LUN, rec=irec) x ; irec = irec + 1
        write(out+CK_LUN, 1010) trim(varname), i, x
      end do

      return
 1000 format("#---- ",a," ----")
 1010 format(a10, i8, ":", 1pe20.10)
      end subroutine


      subroutine idsply(in,out, irec, varname, size)

      implicit none
      integer :: in, out, irec, size, i, x
      character(*) :: varname

      write(out+CK_LUN, 1000) trim(varname)

      do i = 1, size
        !irec = irec + 1
        read(in+CK_LUN, rec=irec) x ; irec = irec + 1
        write(out+CK_LUN, 1010) trim(varname), i, x
      end do

      return
 1000 format("#---- ",a," ----")
 1010 format(a10, i8, ":", I10)
      end subroutine
            
