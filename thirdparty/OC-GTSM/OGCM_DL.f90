!----------------------------------------------------------------------
!-----------------------------------------------------------------------
      PROGRAM OGCM_DL
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      use datetime_module
      use MPI
      use netcdf
      use gsw_mod_toolbox, only: gsw_SA_from_SP, gsw_rho, gsw_CT_from_t&
          , gsw_Nsquared, gsw_mlp
      implicit none
!     This program downloads OGCM data (e.g. GOFS 3.1 HYCOM), and computes  
!     depth-averaged baroclinic pressure gradients, buoyancy frequencies, 
!     sea surface density, and mixed-layer depth ratio on the structured 
!     grid, outputting it to a compressed NetCDF file
!
!      -and (if required)- 
!
!     Computes depth-averaged momentum dispersion, momentum dispersion
!     as a quadratic bottom friction coefficient, and depth-averaged 
!     kinetic energy from OGCM 3D velocities

!     Use below two lines for output in single precision (saves space)
      integer  :: nf90_type = nf90_float  ! Type to store matrix arrays in NetCDF
      integer,parameter  :: sz = 4 !
!     Use below two lines for output in double precision
!      integer  :: nf90_type = nf90_double ! Type to store matrix arrays in NetCDF
!     integer,parameter  :: sz = 8 !

      character(len=80) :: BC3D_Name ! Name of the BC3D netCDF file 
      integer  :: BC3Ddir     ! orientaton of longitude 
      integer  :: NZ, NX, NY  ! dimensions of 3D BC model required
      integer  :: BC3D_DT = 3 ! Default delta t for the OGCM (hours)
      integer  :: dfl = 2     ! Compression level for NetCDF4
      integer  :: TMULT       ! BC3D_DT multiplier (to skip times)
      integer  :: is          ! Start of the x-dimension when writing out array
      integer  :: tsind = 1   ! Start time index in netCDF file
      integer  :: OutType = 1 ! Output type; = 1 baroclinicity, = 2 velocity
      character(len=16) :: TS, TE ! Start and end times
      character(len=3) :: BCServer = 'ftp' ! Server type
      real*8,parameter :: BPGlim = 0.1d0 !upper/lower bound for BPGs
      real*8,parameter :: LatUL = 89d0 !upper limit for Lat in haversine 
      real*8,allocatable,dimension(:) :: BC3D_Lon, BC3D_Lat, BC3D_Z
      real(sz),allocatable,dimension(:,:,:) :: BC3D_SP, BC3D_T, BC3D_BCP
      real(sz),allocatable,dimension(:,:) :: BC2D_NM, BC2D_NB, BC2D_BX,&
                                BC2D_SigTS, BC2D_BY, BC2D_MLD, BC2D_CD,&
                                BC2D_DispX, BC2D_DispY, BC2D_KE
      real*8,allocatable,dimension(:,:) :: DX, DY, DUU, DVV, DUV, DU,  &
                                           DV, B
!     Variables for temporal integration
      type(datetime) :: TSDT, TEDT, CurDT
      integer :: ierr, myProc, mnProc, nextProc, prevProc
      real*8,parameter :: RhoWat0 = 1d3, PI = 3.141592653589793D0,     &
                       deg2rad = PI/180d0, Rearth = 6378206.4d0,       &
                       G = 9.80665d0
      real(sz),parameter :: SigT0 = -999d0, DFV = 1d0 !DFV = 1d4 
      ! Netcdf output info
      character(len=40) :: BC2D_Name ! Name of Output File
      integer :: ncid, ncformat, timenc_id, timenc_dim_id, NX_dim_id,  &
              CD_id, NYYY_dim_id, NY_dim_id, NYY_dim_id, MLD_id, NM_id,&
              SigTS_id, BPGX_id, BPGY_id, NB_id, lon_id, lat_id, KE_id,&
              latc_id, lats_id, lonc_id, strlen_dim_id, DispX_id, &
              DispY_id
!-----------------------------------------------------------------------

!.......Initialize MPI
      call MPI_Init(ierr)
      call MPI_Comm_Size (MPI_Comm_World,mnProc,ierr)   ! Get number of procs
      call MPI_Comm_Rank (MPI_Comm_World,myProc,ierr)   ! Get MPI rank
!      if (myProc.eq.0) then
!         write(6,*) 'mnProc = ',mnProc
!      endif
!      write(6,*) 'myProc = ',myProc
!..... 
!     Get next and previous processors to inform when netCDF
!     file is free to write out into
      nextProc = myProc + 1
      if (nextProc.gt.mnProc-1) nextProc = 0    
      prevProc = myProc - 1
      if (prevProc.lt.0) prevProc = mnProc - 1    
 
!.......Read the Control File
      call Read_Input_File() 

!.......Start the download and process loop
      call OGCM_Run()
  
!.......Finish up the program 
      call MPI_Finalize(ierr)
      stop
     
!-----------------------------------------------------------------------
      contains
!-----------------------------------------------------------------------
!     S U B R O U T I N E  R E A D _ I N P U T _ F I L E
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine Read_Input_File()
      implicit none
      logical :: Fexist
!
      ! Read from standard input
      if (myProc.eq.0) then
         read(5,*) ! Read blank
         read(5,'(A16)') TS;           
         read(5,'(A16)') TE;           
         TSDT = strptime(trim(TS),"%Y-%m-%d %H:%M")
         TEDT = strptime(trim(TE),"%Y-%m-%d %H:%M")
         write(6,*) TSDT%isoformat(' ')            
         write(6,*) TEDT%isoformat(' ')            
         ! Check this a valid time
         if (.not.TSDT%isvalid().or..not.TEDT%isvalid()) then
            write(6,*) 'ERROR: Invalid datetime TS or TE: '   &
                     //'Must be formatted as yyyy-MM-dd HH:mm'
            call MPI_Abort(MPI_Comm_World,0,ierr)
         endif
         read(5,*) TMULT; print *,'TMULT = ', TMULT          
         read(5,'(A3)') BCServer; print *, 'BCServer = ', BCServer
         read(5,*) OutType
         if (OutType == 0) then
            write(6,*) 'INFO: Compute & output buoyancy frequencies'
         elseif (OutType == 1) then
            write(6,*) 'INFO: Compute & output all T and S variables'
         elseif (OutType == 2) then
            write(6,*) 'INFO: Compute & output T, S, U, V'
         else
            write(6,*) 'ERROR: Invalid OutType = ',OutType
            call MPI_Finalize(ierr); stop
         endif
         read(5,'(A40)') BC2D_Name; print *, 'BC2D_Name = ', BC2D_Name
         ncformat = nf90_hdf5
      endif
      ! Broadcast information to all processors
      call MPI_Bcast(TS,16,MPI_Character,0,MPI_Comm_World,ierr)
      call MPI_Bcast(TE,16,MPI_Character,0,MPI_Comm_World,ierr)
      TSDT = strptime(TS,"%Y-%m-%d %H:%M")
      TEDT = strptime(TE,"%Y-%m-%d %H:%M")
      call MPI_Bcast(TMULT,1,MPI_Integer,0,MPI_Comm_World,ierr)
      call MPI_Bcast(BCServer,3,MPI_Character,0,MPI_Comm_World,ierr)
      call MPI_Bcast(OutType,1,MPI_Integer,0,MPI_Comm_World,ierr)
      call MPI_Bcast(BC2D_Name,40,MPI_Character,0,MPI_Comm_World,ierr)
!
      end subroutine Read_Input_File
!
!-----------------------------------------------------------------------
!     S U B R O U T I N E  R E A D _ I N P U T _ F I L E
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine OGCM_Run()
      
      implicit none
      type(timedelta) :: EndCheck
      character(len=3):: hhh 
      integer :: IT, ii
      logical :: OKflag

      ! Set the output filename
      write(hhh,'(I3)') myProc
      BC3D_Name = 'HYCOM'//trim(adjustl(hhh))//'.nc' 
      CurDT = TSDT + timedelta(hours=myProc*BC3D_DT*TMULT)
      write(6,*) 'MyProc = ',myProc,'CurDT = ',CurDT%isoformat(' ')
      EndCheck = TEDT - CurDT; IT = 0
      do while (EndCheck%total_seconds() >= 0)
         do ii = 1,max(OutType,1)
            ! Download new OGCM NetCDF file
            call BC3DDOWNLOAD(CurDT,ii,OKflag)
            if (.not.OKflag) exit
            ! Read the OGCM NetCDF file
            call Read_BC3D_NetCDF(ii)
            ! Calculate the new BC2D terms.
            call Calc_BC2D_Terms(ii)
         end do
         ! Put new BC2D terms in netCDF output file
         call UpdateNetCDF(IT,OKflag)
         ! Update the time 
         CurDT = CurDT + timedelta(hours=mnProc*BC3D_DT*TMULT)
         EndCheck = TEDT - CurDT
         IT = IT + 1
      enddo
!
      end subroutine OGCM_Run
!
!-----------------------------------------------------------------------
!     S U B R O U T I N E  R E A D _ B C 3 D _ N E T C D F
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine Read_BC3D_NetCDF(flag)
      implicit none
      integer,intent(in) :: flag
      logical :: errL
      integer :: NC_ID
      
      ! Open NETCDF 
      call Check_err(NF90_OPEN(BC3D_Name,NF90_NOWRITE,NC_ID))
    
      ! Get the dimensions and time and spatial arrays of the data
      ! or check the newly downloaded data
      call Get_LonLatDepthTime(NC_ID)
     
      ! Read all the necessary temporal varying data
      select case(flag)
        case(1) ! Baroclinicity
          write(6,*) myProc,'Processing T & S'
          ! Practical Salinity
          BC3D_SP = read_nc_var(NC_ID,'salinity  ',1)
          ! Temperature 
          BC3D_T  = read_nc_var(NC_ID,'water_temp',1)
        case(2) ! Velocity
          write(6,*) myProc,'Processing U & V'
          ! East-West velocity
          BC3D_SP = read_nc_var(NC_ID,'water_u   ',1)
          ! North-South Velocity
          BC3D_T  = read_nc_var(NC_ID,'water_v   ',1)
      end select

      ! Close NETCDF
      call Check_err(NF90_CLOSE(NC_ID))
!
!-----------------------------------------------------------------------
      end subroutine Read_BC3D_NetCDF
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine Get_LonLatDepthTime(NC_ID)
      implicit none
      integer,intent(in) :: NC_ID
      integer  :: Temp_ID, i, j, ii 
      real*8   :: BC3D_Lon_s, LatVal 
        
      if (allocated(BC3D_Lon)) then
         ! Test the latitude, longitude and z variables 
         call Check_err(NF90_INQ_VARID(NC_ID,'lon',Temp_ID))
         call Check_err(NF90_GET_VAR(NC_ID,Temp_ID,BC3D_Lon))
      else
         ! First call
         call Check_err(NF90_INQ_DIMID(NC_ID,'lat',Temp_ID))
         call Check_err(NF90_INQUIRE_DIMENSION(NC_ID,Temp_ID,len=NY))
         call Check_err(NF90_INQ_DIMID(NC_ID,'lon',Temp_ID))
         call Check_err(NF90_INQUIRE_DIMENSION(NC_ID,Temp_ID,len=NX))
         call Check_err(NF90_INQ_DIMID(NC_ID,'depth',Temp_ID))
         call Check_err(NF90_INQUIRE_DIMENSION(NC_ID,Temp_ID,len=NZ))
         
         ! Allocate the lat lon and z arrays first 
         allocate(BC3D_Lat(NY),BC3D_Lon(NX),BC3D_Z(NZ))
         allocate(BC3D_SP(NX,NY,NZ),BC3D_T(NX,NY,NZ),BC3D_BCP(NX,NY,NZ))
         allocate(BC2D_NB(NX,NY),BC2D_NM(NX,NY),&
                  BC2D_SigTS(NX,NY),BC2D_MLD(NX,NY))
         allocate(BC2D_BX(NX,NY),BC2D_BY(NX,NY-1),DX(NX,NY),DY(NX,NY-1))
         if (OutType.eq.2) then
            allocate(BC2D_CD(NX,NY-2),BC2D_DispX(NX,NY-2),&
                     BC2D_DispY(NX,NY-2))
            allocate(DUU(NX,NY), DVV(NX,NY), BC2D_KE(NX,NY),&
                     DUV(NX,NY), DU(NX,NY), DV(NX,NY), B(NX,NY))
         endif

         ! Read the latitude, longitude and z variables 
         call Check_err(NF90_INQ_VARID(NC_ID,'lat',Temp_ID))
         call Check_err(NF90_GET_VAR(NC_ID,Temp_ID,BC3D_Lat))
         call Check_err(NF90_INQ_VARID(NC_ID,'lon',Temp_ID))
         call Check_err(NF90_GET_VAR(NC_ID,Temp_ID,BC3D_Lon))
         call Check_err(NF90_INQ_VARID(NC_ID,'depth',Temp_ID))
         call Check_err(NF90_GET_VAR(NC_ID,Temp_ID,BC3D_Z, &
                        start=[1],count=[NZ]))
         ! Pre-computing the distances on the sphere
         do j = 1,NY 
            do i = 1,NX
               LatVal  = min(LatUL,BC3D_Lat(j))
               ii = i + 1
               if (i.eq.NX) ii = 1
               DX(i,j) = haversine(BC3D_Lon(i),BC3D_Lon(ii),   &
                                   LatVal,LatVal) 
               if (j < NY) then
                  DY(i,j) = haversine(BC3D_Lon(i),BC3D_Lon(i), &
                                      BC3D_Lat(j),BC3D_Lat(j+1)) 
               endif 
            enddo
         enddo
      endif
      BC3D_Lon_s = BC3D_Lon(1)
      if (BC3D_Lon_s < 0d0) then
         ! Represents seam at 180/-180
         BC3Ddir = -1
         !write(6,*) 'Lon of ',trim(BC3D_Name),&
         !           ' is defined from -180 to 180'
         is = 1
      else
         ! Represents seam at 0/360
         BC3Ddir =  1
         is = maxloc(BC3D_Lon,dim=1,mask=BC3D_Lon < 180d0) + 1
         !write(6,*) 'Lon of ',trim(BC3D_Name),&
         !           ' is defined from 0 to 360, is = ',is
      endif
      
      end subroutine Get_LonLatDepthTime
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      function read_nc_var(NC_ID,varname,BCIT_IN) result (Var)
      implicit none

      integer,intent(in) :: NC_ID, BCIT_IN
      character(10),intent(in) :: varname
      integer :: Temp_ID, i, j, k 
      real(sz) :: FV, SF, OS
      real(sz),allocatable :: Var(:,:,:)

      call Check_err(NF90_INQ_VARID(NC_ID,trim(varname),Temp_ID))
      call Check_err(NF90_GET_ATT(NC_ID,Temp_ID,'_FillValue',FV))
      call Check_err(NF90_GET_ATT(NC_ID,Temp_ID,'scale_factor',SF))
      call Check_err(NF90_GET_ATT(NC_ID,Temp_ID,'add_offset',OS))
      allocate(Var(NX,NY,NZ)) 
      call Check_err(NF90_GET_VAR(NC_ID,Temp_ID,Var))
      ! Add on a little for buffer
      FV = FV + 1d-3
      do j = 1,NY
         do i = 1,NX
            do k = 1,NZ
               if (Var(i,j,k) > FV) Var(i,j,k) = Var(i,j,k)*SF+OS
            enddo
         enddo
      enddo
!
!-----------------------------------------------------------------------
      end function read_nc_var
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     S U B R O U T I N E   C H E C K  _  E R R
!-----------------------------------------------------------------------
!     jgf49.17.02 Checks the return value from netCDF calls; if there
!     was an error, it writes the error message to the screen
!-----------------------------------------------------------------------
      subroutine check_err(iret)
      implicit none

      integer, intent(in) :: iret
      

      if (iret .ne. nf90_noerr) then
         write(6,*) 'netcdf check_err error' 
         call MPI_Abort(MPI_Comm_World,iret,ierr)
      endif
!-----------------------------------------------------------------------
      end subroutine check_err
!-----------------------------------------------------------------------
!----------------------------------------------------------------------
      subroutine BC3DDOWNLOAD(TimeIn,flag,OKflag)
! PV - 02/16/2022
#ifdef __INTEL_COMPILER
      use ifport, only: systemqq, sleep, getlasterrorqq
#endif

      implicit none
     
      type(datetime),intent(in) :: TimeIn
      integer,intent(in) :: flag
      logical,intent(out) :: OKflag
      logical(4) :: resOK
      type(datetime)     :: TimeOff
      character(len=200) :: fileserver, vars, times, filemid, & 
                            fileend, expt, options, filegrid
      character(len=280) :: FullPath, FullPathP, line
      character(len=3)   :: hhh 
      character(len=4)   :: yyyy 
      character(len=5)   :: command 
      integer            :: iter, stat, fid, resINT, &
                            mpistat(mpi_status_size)
      integer*8          :: Fsize, FSizeTrue 
       
      ! OKflag is false by default    
      OKflag = .false.
      
      ! Add on the time based on processor number 
      TimeOff  = TimeIn - timedelta(hours=12)
      yyyy     = TimeOff%strftime("%Y")

      ! GLBv0.08: GOFS 3.1 on lat lon structured grid
	  ! CPB 1/21/2022: changed to updated GOFS 3.1 grid (GLBy0.08)
      ! Get servername
      if (BCServer == 'ftp') then
         ! FTP
         fileserver = 'ftp://ftp.hycom.org/datasets/'
         command    = 'curl ' !'wget '
         options    = ' -s --connect-timeout 30'
      elseif (BCServer == 'ncs') then
         ! NCSS
         fileserver = '"http://ncss.hycom.org/thredds/ncss/'
         command    = 'curl '
         options    = ' -s --connect-timeout 30'
      elseif (BCServer == 'opd') then
         ! OpenDap
         command    = 'ncks '
         fileserver = '"http://tds.hycom.org/thredds/dodsC/'
         options    = ' -4 -O '
      endif

      filegrid = 'GLBv0.08/'

      if (TimeOff%getYear().ge.2018) then
         ! GOFS 3.1 analysis
         filegrid = 'GLBy0.08/'
         expt = 'expt_93.0'
         filemid = '/hycom_glby_930_'
         if (BCServer == 'ftp') then
            expt = trim(expt)//'/data/hindcasts/'//yyyy
         endif
      elseif (TimeOff%getYear().eq.2017) then
         ! GOFS 3.1 analysis
         if (TimeOff%getMonth().ge.10) then
            expt = 'expt_92.9'
            filemid = '/hycom_glbv_929_'
         elseif (TimeOff%getMonth().ge.6) then
            expt = 'expt_57.7'
            filemid = '/hycom_GLBv0.08_577_'
         elseif (TimeOff%getMonth().ge.2) then
            expt = 'expt_92.8'
            filemid = '/hycom_glbv_928_'
         else
            expt = 'expt_57.2'
            filemid = '/hycom_GLBv0.08_572_'
         endif
         if (BCServer == 'ftp') then
            expt = trim(expt)//'/data/hindcasts/'//yyyy
         endif
      elseif (TimeOff%getYear().eq.2016.and.TimeOff%getMonth().ge.5) then
         ! GOFS 3.1 analysis
         filemid = '/hycom_GLBv0.08_572_'
         expt = 'expt_57.2'
         if (BCServer == 'ftp') then
            expt = trim(expt)//'/data/hindcasts/'//yyyy
         endif
      elseif (TimeOff%getYear().eq.2016.or.(TimeOff%getYear().eq.2015.and.&
              TimeOff%getYear().eq.12.and.TimeOff%getDay().eq.31)) then
         ! GOFS 3.1 analysis
         filemid = '/hycom_GLBv0.08_563_'
         expt = 'expt_56.3'
         if (BCServer == 'ftp') then
            expt = trim(expt)//'/data/hindcasts/'//yyyy
         endif
      elseif (TimeOff%getYear().ge.1994.or.&
             (TimeOff%getYear().eq.1993.and.TimeOff%getMonth().ge.10)) then
         ! GOFS 3.1 reanalysis, expt 53.X with the year
         expt = 'expt_53.X/data/'//yyyy
         if (TimeOff%getYear().lt.2000) then
            filemid = '/hycom_GLBv0.08_530_'
         elseif (TimeOff%getYear().eq.2000) then
            filemid = '/hycom_GLBv0.08_531_'
         elseif (TimeOff%getYear().lt.2004) then
            filemid = '/hycom_GLBv0.08_532_'
         elseif (TimeOff%getYear().lt.2006.and.TimeOff%getMonth().lt.7) then
            filemid = '/hycom_GLBv0.08_533_'
         elseif (TimeOff%getYear().lt.2008) then
            filemid = '/hycom_GLBv0.08_534_'
         elseif (TimeOff%getYear().lt.2010) then
            filemid = '/hycom_GLBv0.08_535_'
         elseif (TimeOff%getYear().lt.2012) then
            filemid = '/hycom_GLBv0.08_536_'
         elseif (TimeOff%getYear().lt.2014) then
            filemid = '/hycom_GLBv0.08_537_'
         elseif (TimeOff%getYear().eq.2014) then
            filemid = '/hycom_GLBv0.08_538_'
         elseif (TimeOff%getYear().lt.2015) then
            filemid = '/hycom_GLBv0.08_539_'
         endif
      else
         write(6,*) 'Error: no GOFS 3.1 data available before Oct 1993.'
      endif
      hhh = ''
      write(hhh,'(I3)') myProc
      if (BCServer == 'ftp') then
         vars    = ''
         times   = TimeOff%strftime("%Y%m%d")//'12_t0'&
                 //TimeOff%strftime("%H")
         if (expt(6:6) == '9') then
            if (flag == 1) then
               fileend = '_ts3z.nc"'
            elseif (flag == 2) then
               fileend = '_uv3z.nc"'
            endif
         else
            if (flag == 2) then
               OKflag = .true.; return
            endif 
            fileend = '.nc"'
         endif
      elseif (BCServer == 'ncs') then
         if (flag == 1) then
            vars = '?var=salinity&var=water_temp'
         elseif (flag == 2) then
            vars = '?var=water_u&var=water_v'
         endif
         times   = '&time='//TimeIn%strftime("%Y-%m-%dT%H")&
                 //'%3A00%3A00Z'
         filemid = ''
         fileend = '&vertStride=1&addLatLon=true&accept=netcdf4"'
      elseif (BCServer == 'opd') then
         if (flag == 1) then
            vars = '" -v salinity,water_temp'
         elseif (flag == 2) then
            vars = '" -v water_u,water_v'
         endif
         times   = ' -d time,'//trim(hhh)
         fileend = ''
         filemid = ''
      endif
      ! Get the final path
      Fullpath = trim(fileserver)//trim(filegrid)//trim(expt)//trim(vars) &
               //trim(filemid)//trim(times)//trim(fileend) 
      
      ! Let's get expected filesize
      if (BCServer == 'ftp') then
         FSizeTrue = 0; 
         !do while(FSizeTrue.eq.0) 
            ! Get filesize at remote location
! PV - 02/16/2022
#ifdef __INTEL_COMPILER
            resOK = systemqq(command//trim(Fullpath)//trim(options)//  &
                       ' -I -L -o filesize'//trim(adjustl(hhh))//'.txt')
#else
            resOK = .false.
            resINT = system(command//trim(Fullpath)//trim(options)//  &
                       ' -I -L -o filesize'//trim(adjustl(hhh))//'.txt')
            if (resINT == 0) resOK = .true.
#endif

            if (.not.resOK) then
! PV - 02/16/2022
#ifdef __INTEL_COMPILER
               resINT = getlasterrorqq()
#else
               resINT = ierrno()
#endif
               write(6,*) 'Error: ',resINT,myProc
            endif
            open(newunit=fid,                                     &
                 file='filesize'//trim(adjustl(hhh))//'.txt',     &
                 status='old',action='read')
               do ! read until eof
                  read(fid,'(A280)',iostat=stat) line
                  if (stat.ne.0) exit
                  if (line(1:15).eq.'Content-Length:') then
                     read(line(17:280),*) FSizeTrue
                     exit
                  endif
               enddo
            close(fid)
            if (FSizeTrue.eq.0) then
               ! Non-existant file. 
               ! Try forecast
               !TimeOff = TimeOff - timedelta(hours=24)
               !times   = TimeOff%strftime("%Y%m%d")//'12_t0'       &
               !        //TimeOff%getHour + 24
               !Fullpath = trim(fileserver)//trim(expt)//trim(vars) &
               !        //trim(filemid)//trim(times)//trim(fileend)
               !if (FullPath.eq.FullPathP) then
               !   write(6,*) myProc,                               &
               !             'Previous and current paths are equal'
               !endif
               ! Try BC3D_DT hours ahead
               !TimeOff = TimeOff + timedelta(hours=BC3D_DT)
               !times   = TimeOff%strftime("%Y%m%d")//'12_t0'       &
               !        //TimeOff%strftime("%H")
               !Fullpath = trim(fileserver)//trim(expt)//trim(vars) &
               !        //trim(filemid)//trim(times)//trim(fileend)
               write(6,*) myProc,'no data available at: ',trim(FullPath)
               return ! return error
            endif 
         !enddo
      else
         ! Expect over a 500MB
         FSizeTrue = 5e8
      endif
      
      write(6,*) myProc, 'Downloading: ',trim(FullPath)
      iter = 0; resOK = .false.
      do while (.not.resOK.and.iter < 25)
! PV - 02/16/2022
#ifdef __INTEL_COMPILER
            resOK = systemqq(command//trim(Fullpath)//trim(options)//&
                       ' -o '//trim(BC3D_Name))
#else
         resOK = .false.
         resINT = system(command//trim(Fullpath)//trim(options)//&
                         ' -o '//trim(BC3D_Name))
         if (resINT == 0) resOK = .true.
#endif

         if (resOK) then
            ! Expect over a 500MB
            inquire(file=trim(BC3D_Name),size=Fsize)
            if (Fsize < FSizeTrue) then
               write(6,*) 'FSize of '//trim(BC3D_Name)//': ',     &
                           FSize, FSizeTrue
               resOK = .false.
            endif
         endif
         if (.not.resOK) then
            ! Let's wait 5 seconds before trying again
            write(6,*) 'Problem downloading GOFS 3.1 NetCDF. '&
                    // 'Try again after 5 s'
            call sleep(5) ; iter = iter + 1; 
         endif
      enddo
      if (.not.resOK) then
         write(6,*) 'Error in downloading GOFS 3.1 NetCDF.'
         call MPI_Abort(MPI_Comm_World,0,ierr) 
      else
         write(6,*) myProc,'download successfull'
         OKflag = .true.
      endif
!----------------------------------------------------------------------
      end subroutine BC3DDOWNLOAD
!-----------------------------------------------------------------------
!     S U B R O U T I N E  C A L C _ B C 2 D _ T E R M S
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine Calc_BC2D_Terms(flag)
      implicit none
     
      integer,intent(in) :: flag 
      real*8 :: rho, rhoavg, bx_avg, by_avg, bx_z, by_z, DZ
      real*8 :: uavg, vavg, DUUX, DVVY, DUVX, DUVY, Vel, CD, CDx, CDy 
      real*8 :: SA(NZ), CT(NZ), Lat(NZ), N2(NZ-1), zmid(NZ-1)
      integer  :: i, j, k, ii, kb, kbx, kby, ip, im
      real(sz) :: FV = -3d4, FVP = -3d4 + 1d-3, Vs = 1d-3 

      select case(flag)
      
      !! FOR SALINITY & TEMPERATURE 
      case(1)
      ! Loop over all cells to get the baroclinic pressure and the
      ! buoyancy frequencies
!!$omp parallel private(Lat,kb,DZ,rho,rhoavg,kbx,kby,bx_avg,by_avg, &
!!$omp bx_z,by_z,SA,CT,N2,zmid)
!!$omp do
      do j = 1,NY
         ! Need an array of latitude values for Nsquared call
         Lat = BC3D_Lat(j)
         do i = 1,NX
            ! Initialize the free surface density
            BC2D_SigTS(i,j) = SigT0
            kb = 0; 
            do k = 1,NZ
               if (k > 1) DZ = BC3D_Z(k) - BC3D_Z(k-1)
               ! For the baroclinic pressure gradient and surface
               ! density 
               if (BC3D_SP(i,j,k) > FVP .and. BC3D_T(i,j,k) > FVP) then
                  ! Change SP to SA
                  SA(k) = gsw_SA_from_SP(max(2d0,BC3D_SP(i,j,k)),&
                          BC3D_Z(k),BC3D_Lon(i),BC3D_Lat(j))
                  ! Change T to CT
                  CT(k) = gsw_CT_from_T(SA(k),dble(BC3D_T(i,j,k)),&
                          BC3D_Z(k))
                  ! Save previous rho
                  if (k > 1) rhoavg = rho
                  ! Calculate rho
                  rho = gsw_rho(SA(k),CT(k),BC3D_Z(k))
                  if (k > 1) then
                     ! The trapezoidal rule
                     rhoavg = 0.5d0*(rhoavg + rho)
                     BC3D_BCP(i,j,k) = BC3D_BCP(i,j,k-1)    &
                                     + DZ*(rhoavg - RhoWat0)
                  else
                     BC2D_SigTS(i,j) = rho - RhoWat0
                     BC3D_BCP(i,j,k) = 0.0d0
                  endif
                  kb = k
               else
                  ! Reached the bottom
                  BC3D_BCP(i,j,k) = FV
               endif
            enddo
            ! The buoyancy frequencies
            if (kb > 1) then
               ! Get the buoyancy frequency profile
               call gsw_Nsquared(SA(1:kb),CT(1:kb),BC3D_Z(1:kb),  &
                                 Lat(1:kb),N2(1:kb-1),zmid(1:kb-1))
               ! Get bottom (seabed) value
               BC2D_NB(i,j) = sqrt(max(0.0d0,N2(kb-1)))
               ! Integrate (here N2 is defined in the middle of depth
               ! layers (zmid), and I assume it to be the average in
               ! that depth range for simplicity) 
               BC2D_NM(i,j)  = 0.0d0  
               do k = 1,kb-1
                  DZ = BC3D_Z(k+1) - BC3D_Z(k)
                  BC2D_NM(i,j) = BC2D_NM(i,j) + DZ*sqrt(max(0.d0,N2(k)))
               enddo
               ! Divide by the depth
               BC2D_NM(i,j) = BC2D_NM(i,j)/BC3D_Z(kb)
               ! Get the mixed-layer depth
               !BC2D_MLD(i,j) = min(BC3D_Z(kb),&
               !                gsw_mlp(SA(1:kb),CT(1:kb),BC3D_Z(1:kb)))
               ! Get the mixed-layer depth ratio
               BC2D_MLD(i,j) = min(DFV,&
                     gsw_mlp(SA(1:kb),CT(1:kb),BC3D_Z(1:kb))/BC3D_Z(kb))
            else
               ! Let's just set to zero value 
               BC2D_NB(i,j)  = 0.0d0
               BC2D_NM(i,j)  = 0.0d0
               BC2D_MLD(i,j) = DFV
            endif
         enddo
      enddo 
!!$omp end do
      ! Now calculate the gradients of the BCP and integrate,
      ! and calculate gradients of DUU, DVU and DVV
      ! (we do central difference about the mid-point)
!!$omp do
      do j = 1,NY 
         do i = 1,NX
            kbx = 0; kby = 0;
            do k = 1,NZ
               if (k > 1) then
                  bx_avg = bx_z
                  by_avg = by_z 
                  ! If we hit bottom
                  if (BC3D_BCP(i,j,k) < FVP) exit
                  ii = i + 1
                  if (i == NX) ii = 1
                  if (BC3D_BCP(ii,j,k) > FVP) then
                     ! x-gradient at this level
                     bx_z = ( BC3D_BCP(ii,j,k) - BC3D_BCP(i,j,k) ) &
                          / DX(i,j)
                     ! trapezoidal rule
                     bx_avg = 0.5d0*(bx_avg + bx_z)
                     BC2D_BX(i,j) = BC2D_BX(i,j) +      &
                                   (BC3D_Z(k)-BC3D_Z(k-1))*bx_avg
                     kbx = k
                  endif
                  if (j < NY) then
                     if (BC3D_BCP(i,j+1,k) > FVP) then
                        ! y-gradient at this level
                        by_z = ( BC3D_BCP(i,j+1,k) -      &
                                 BC3D_BCP(i,j,k) ) / DY(i,j)
                        ! trapezoidal rule
                        by_avg = 0.5d0*(by_avg + by_z)
                        BC2D_BY(i,j) = BC2D_BY(i,j) +     & 
                                      (BC3D_Z(k)-BC3D_Z(k-1))*by_avg
                        kby = k
                     endif
                  endif
               else
                  ! Is always zero at top
                  bx_z = 0.0d0
                  by_z = 0.0d0
                  BC2D_BX(i,j) = 0.0d0
                  if (j < NY) BC2D_BY(i,j) = 0.0d0
               endif
            enddo
            ! Get the depth-averaged value
            if (kbx > 1) then
               BC2D_BX(i,j) = G/RhoWat0*min(BPGlim,max(-BPGlim, &
                              BC2D_BX(i,j)/BC3D_Z(kbx)))

            endif
            if (kby > 1) then
               BC2D_BY(i,j) = G/RhoWat0*min(BPGlim,max(-BPGlim, &
                              BC2D_BY(i,j)/BC3D_Z(kby)))
            endif
         enddo
      enddo
!!$omp end do
!!$omp end parallel

      !! FOR VELOCITIES 
      case(2)
!!$omp parallel private(kb,DZ,uavg,vavg,ip,im,DUUX,DUVY,DUVX,DVVY,Vel,CDx,CDy)
!!$omp do
      ! Loop over all cells to get the depth-integrated momentum dispersion
      do j = 1,NY
         do i = 1,NX
            ! Initialize the dispersion values
            DUU(i,j) = 0d0; DVV(i,j) = 0d0; DUV(i,j) = 0d0;
            DU(i,j) = 0d0; DV(i,j) = 0d0; B(i,j) = 0d0; 
            BC2D_KE(i,j) = 0d0; kb = 0; 
            do k = 2,NZ
               ! For depth-integrated dispersion 
               if (BC3D_SP(i,j,k) < FVP .or. BC3D_T(i,j,k) < FVP) exit
               DZ = BC3D_Z(k) - BC3D_Z(k-1)
               ! The trapezoidal rule
               uavg = 0.5d0*(BC3D_SP(i,j,k) + BC3D_SP(i,j,k-1))
               vavg = 0.5d0*(BC3D_T(i,j,k) + BC3D_T(i,j,k-1))
               DUU(i,j) = DUU(i,j) + DZ*uavg*uavg
               DVV(i,j) = DVV(i,j) + DZ*vavg*vavg
               DUV(i,j) = DUV(i,j) + DZ*uavg*vavg
               DU(i,j)  = DU(i,j) + DZ*uavg
               DV(i,j)  = DV(i,j) + DZ*vavg
               kb = k
            enddo
            ! Now calculate dispersion values in depth-averaged form
            if (kb > 1) then
               B(i,j)   = BC3D_Z(kb)
               ! Depth-averaged velocities
               DU(i,j)  = DU(i,j)/B(i,j)
               DV(i,j)  = DV(i,j)/B(i,j)
               ! Depth-averaged kinetic energy
               BC2D_KE(i,j)  = 0.5d0*( DUU(i,j) + DVV(i,j) ) / B(i,j)
               ! Depth-averaged dispersion
               DUU(i,j) = DUU(i,j)/B(i,j) - DU(i,j)*DU(i,j) 
               DVV(i,j) = DVV(i,j)/B(i,j) - DV(i,j)*DV(i,j)
               DUV(i,j) = DUV(i,j)/B(i,j) - DU(i,j)*DV(i,j)
            endif
         enddo
      enddo 
!!$omp end do

      ! Now calculate the gradients of DUU, DUV and DVV
      ! (we do central difference about original point)
!!$omp do
      do j = 2,NY-1 
         do i = 1,NX
            ip = i + 1
            im = i - 1
            if (ip > NX) ip = 1
            if (im < 1)  im = NX
                        
            if (abs(DU(i,j)) < Vs .or. abs(DV(i,j)) < Vs .or.      &
                B(ip,j) < BC3D_Z(2) .or. B(im,j) < BC3D_Z(2) .or.  & 
                B(i,j+1) < BC3D_Z(2) .or. B(i,j-1) < BC3D_Z(2)) then
               BC2D_CD(i,j-1) = 0d0
               BC2D_DispX(i,j-1) = 0d0
               BC2D_DispY(i,j-1) = 0d0
               cycle
            endif

            DUUX = 0.5d0*( DUU(ip,j) - DUU(im,j) ) / DX(i,j)
            DUVX = 0.5d0*( DUV(ip,j) - DUV(im,j) ) / DX(i,j)
            DVVY = 0.5d0*( DVV(i,j+1)  - DVV(i,j-1) ) / DY(i,j)
            DUVY = 0.5d0*( DUV(i,j+1)  - DUV(i,j-1) ) / DY(i,j)

            BC2D_DispX(i,j-1) = DUUX + DUVY
            BC2D_DispY(i,j-1) = DVVY + DUVX

            CD = BC2D_DispX(i,j-1)*DU(i,j) + BC2D_DispY(i,j-1)*DV(i,j)
            if (CD <= 0d0) then
               BC2D_CD(i,j-1) = 0d0
            else
               Vel  = sqrt(DU(i,j)*DU(i,j) + DV(i,j)*DV(i,j)) 
               BC2D_CD(i,j-1) = CD * B(i,j) / Vel**3
            endif 
            !CDx  = B(i,j) * ( DUUX + DUVY ) / ( DU(i,j) * Vel )
            !CDy  = B(i,j) * ( DVVY + DUVX ) / ( DV(i,j) * Vel )   
            !BC2D_CD(i,j-1) = sqrt(CDx*CDx + CDy*CDy)
            !BC2D_CD(i,j-1) = max(abs(CDx),abs(CDy))
         enddo
      enddo
!!$omp end do
!!$omp end parallel
      end select
!-----------------------------------------------------------------------
      end subroutine Calc_BC2D_Terms
!-----------------------------------------------------------------------
      function haversine(deglon1,deglon2,deglat1,deglat2) result (dist)
          ! great circle distance -- adapted from Matlab 
          real*8,intent(in) :: deglat1,deglon1,deglat2,deglon2
          real*8 :: a,c,dist,dlat,dlon,lat1,lat2
 
          dlat = deg2rad*(deglat2-deglat1)
          dlon = deg2rad*(deglon2-deglon1)
          lat1 = deg2rad*(deglat1)
          lat2 = deg2rad*(deglat2)
          a = ( sin(0.5d0*dlat) )**2 + &
                cos(lat1) * cos(lat2) * ( sin(0.5d0*dlon) )**2
          c = 2d0*asin( sqrt(a) )
          dist = Rearth*c
      end function haversine
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     S U B R O U T I N E    I N I T _ N E T C D F
!-----------------------------------------------------------------------
      subroutine initNetCDF()
      implicit none
      logical :: Fexist 
      integer :: ierr, data_dims(3), TEE 
      type(timedelta) :: DT
      type(datetime)  :: TNDT
      character(len=16) :: TSS

      if (myProc.eq.0) then   
         inquire(file=BC2D_Name,exist=Fexist)
      else
         Fexist = .true.
      endif
      if (.not.Fexist) then
         ! Open file  
         call check_err(nf90_create(BC2D_Name, ncformat, ncid))
         ! Define dimensions
         call check_err(nf90_def_dim(ncid, 'time', nf90_unlimited, &
                        timenc_dim_id))
         call check_err(nf90_def_dim(ncid, 'strlen', 16, strlen_dim_id))
         call check_err(nf90_def_dim(ncid,'NX',NX,NX_dim_id))
         call check_err(nf90_def_dim(ncid,'NY',NY,NY_dim_id))
         if (OutType.gt.0) &
           call check_err(nf90_def_dim(ncid,'NYY',NY-1,NYY_dim_id))
		 ! CPB 01/21/2022: changed condition to always define NYYY_dim_id because of netcdf error
		 ! This should not be an issue since it just results in the latitude for Disp to be written
		 ! even if it is not used in the adcirc run
         if (OutType.gt.0) &
           call check_err(nf90_def_dim(ncid,'NYYY',NY-2,NYYY_dim_id))
         ! Define vars 
         data_dims = [strlen_dim_id, timenc_dim_id, 1]
         call def_var_att(ncid,'time',nf90_char, data_dims(1:2), &
                          timenc_id,'UTC datetime','YYYY-MM-DD HH:mm')
         call def_var_att(ncid,'lon',nf90_double, [NX_dim_id], &
                          lon_id,'longitude','degrees')
         call def_var_att(ncid,'lat',nf90_double, [NY_dim_id], &
                          lat_id,'latitude','degrees')
         if (OutType.gt.0) then
           call def_var_att(ncid,'lonc',nf90_double, [NX_dim_id], &
                          lonc_id,'longitude for BPGX','degrees')
           call def_var_att(ncid,'latc',nf90_double, [NYY_dim_id], &
                          latc_id,'latitude for BPGY','degrees')
           call def_var_att(ncid,'lats',nf90_double, [NYYY_dim_id], &
                          lats_id,'latitude for Disp','degrees')
           data_dims = [NX_dim_id, NY_dim_id, timenc_dim_id]
           call def_var_att(ncid,'BPGX',nf90_type, data_dims, BPGX_id, &
                        'east-west depth-averaged baroclinic pressure '&
                        //'gradient','ms^-2')
           data_dims = [NX_dim_id, NYY_dim_id, timenc_dim_id]
           call def_var_att(ncid,'BPGY',nf90_type, data_dims, BPGY_id, &
                      'north-south depth-averaged baroclinic pressure '&
                        //'gradient','ms^-2')
           data_dims = [NX_dim_id, NY_dim_id, timenc_dim_id]
           call def_var_att(ncid,'SigTS',nf90_type, data_dims,   &
                          SigTS_id, 'surface sigmat density',  &
                          'kgm^-3',SigT0)
           call def_var_att(ncid,'MLD',nf90_type, data_dims, MLD_id, &
                          'mixed-layer depth ratio','[]',DFV)
         !                 'mixed-layer depth','m',DFV)
         endif
         data_dims = [NX_dim_id, NY_dim_id, timenc_dim_id]
         call def_var_att(ncid,'NB',nf90_type, data_dims, NB_id, &
                         'buoyancy frequency at the seabed','s^-1')
         call def_var_att(ncid,'NM',nf90_type, data_dims, NM_id, &
                         'depth-averaged buoyancy frequency','s^-1')
         if (OutType.eq.2) then
            call def_var_att(ncid,'KE',nf90_type, data_dims, KE_id, &
                         'depth-averaged kinetic energy','m^2s^-2')
            data_dims = [NX_dim_id, NYYY_dim_id, timenc_dim_id]
            call def_var_att(ncid,'CDisp',nf90_type, data_dims, CD_id, &
                   'momentum dispersion quadratic bottom friction','[]')
            call def_var_att(ncid,'DispX',nf90_type, data_dims, &
                  DispX_id, 'depth-averaged x-momentum dispersion',&
                  'ms^-2')
            call def_var_att(ncid,'DispY',nf90_type, data_dims, &
                  DispY_id, 'depth-averaged y-momentum dispersion',&
                  'ms^-2')
         endif
         ! Allowing vars to deflate
         call check_err(nf90_def_var_deflate(ncid, lon_id, 1, 1, dfl))
         call check_err(nf90_def_var_deflate(ncid, lat_id, 1, 1, dfl))
         if (OutType.gt.0) then
            call check_err(nf90_def_var_deflate(ncid, lonc_id, 1,1,dfl))
            call check_err(nf90_def_var_deflate(ncid, latc_id, 1,1,dfl))
            call check_err(nf90_def_var_deflate(ncid, lats_id, 1,1,dfl))
            call check_err(nf90_def_var_deflate(ncid, BPGX_id, 1,1,dfl))
            call check_err(nf90_def_var_deflate(ncid, BPGY_id, 1,1,dfl))
            call check_err(nf90_def_var_deflate(ncid, SigTS_id,1,1,dfl))
            call check_err(nf90_def_var_deflate(ncid, MLD_id,  1,1,dfl))
         endif
         call check_err(nf90_def_var_deflate(ncid, NB_id, 1, 1, dfl))
         call check_err(nf90_def_var_deflate(ncid, NM_id, 1, 1, dfl))
         if (OutType.eq.2) then
            call check_err(nf90_def_var_deflate(ncid, KE_id,   1,1,dfl))
            call check_err(nf90_def_var_deflate(ncid, CD_id,   1,1,dfl))
            call check_err(nf90_def_var_deflate(ncid, DispX_id,1,1,dfl))
            call check_err(nf90_def_var_deflate(ncid, DispY_id,1,1,dfl))
         endif
         ! Close define mode
         call check_err(nf90_close(ncid))
         
         ! Put X, Y, Z on
         call check_err(nf90_open(BC2D_Name, nf90_write, ncid))
         if (is == 1) then
            ! Already -180/180 orientation
            call check_err(nf90_put_var(ncid, lon_id, BC3D_Lon))
            if (OutType.gt.0) &
               call check_err(nf90_put_var(ncid, lonc_id,   &
                 [0.5d0*(BC3D_Lon(1:NX-1) + BC3D_Lon(2:NX)),&
                  0.5d0*(BC3D_Lon(NX) + BC3D_Lon(1)+360d0)]))
         else
            ! Change to -180/180 orientation
            call check_err(nf90_put_var(ncid, lon_id,      &
                 [BC3D_Lon(is:NX)-360d0, BC3D_Lon(1:is-1)]))
            if (OutType.gt.0) &
               call check_err(nf90_put_var(ncid, lonc_id,             &
                 [0.5d0*(BC3D_Lon(is:NX-1) + BC3D_Lon(is+1:NX))-360d0,&
                  0.5d0*(BC3D_Lon(NX)-360d0 + BC3D_Lon(1)),           &
                  0.5d0*(BC3D_Lon(1:is-1)  + BC3D_Lon(2:is))]))
         endif
         call check_err(nf90_put_var(ncid, lat_id, BC3D_Lat))
         if (OutType.gt.0) &
            call check_err(nf90_put_var(ncid, latc_id, &
              0.5d0*(BC3D_Lat(1:NY-1) + BC3D_Lat(2:NY))))
         if (OutType.eq.2) &
            call check_err(nf90_put_var(ncid, lats_id,BC3D_Lat(2:NY-1)))
         call check_err(nf90_close(ncid))
      endif
      ! Barrier to ensure wait until netcdf is created by first
      ! processor
      call MPI_Barrier(MPI_Comm_World,ierr)     
      
      if (Fexist) then 
         ! Get the dim and var ids
         call check_err(nf90_open(BC2D_Name, nf90_nowrite, ncid))
         call check_err(nf90_inq_dimid(ncid,'time', timenc_dim_id))
         call Check_err(nf90_inquire_dimension(ncid,&
                        timenc_dim_id,len=TEE))
         call check_err(nf90_inq_varid(ncid,'time', timenc_id))
         ! Let's get the start time index
         do tsind = 1,TEE
            call check_err(nf90_get_var(ncid, timenc_id, &
                           TSS,start=[1, tsind],count=[16, 1]))
            TNDT = strptime(trim(TSS),"%Y-%m-%d %H:%M")
            DT   = TSDT - TNDT
            if (DT%total_seconds() <= 0) exit
         enddo
         if (myProc.eq.0) write(6,*) 'tsind = ',tsind
         if (OutType.gt.0) then
            call check_err(nf90_inq_varid(ncid,'BPGX', BPGX_id))
            call check_err(nf90_inq_varid(ncid,'BPGY', BPGY_id))
            call check_err(nf90_inq_varid(ncid,'SigTS', SigTS_id))
            call check_err(nf90_inq_varid(ncid,'MLD', MLD_id))
         endif
         call check_err(nf90_inq_varid(ncid,'NB', NB_id))    
         call check_err(nf90_inq_varid(ncid,'NM', NM_id))
         if (OutType.eq.2) then
            call check_err(nf90_inq_varid(ncid,'KE', KE_id))
            call check_err(nf90_inq_varid(ncid,'CDisp', CD_id))
            call check_err(nf90_inq_varid(ncid,'DispX', DispX_id))
            call check_err(nf90_inq_varid(ncid,'DispY', DispY_id))
         endif
         call check_err(nf90_close(ncid))
      endif
 
      end subroutine initNetCDF
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     S U B R O U T I N E    D E F _ V A R _ A T T
!-----------------------------------------------------------------------
      subroutine def_var_att(ncid, Sname, NFtype, dims, var_id, Lname, &
                 Units, FillValue, Factor, Offset)
      implicit none
      integer,intent(in)  :: ncid, NFtype 
      integer,intent(in),dimension(:) :: dims
      integer,intent(out) :: var_id   
      character(*),intent(in) :: Sname, Lname, Units 
      real(sz),intent(in),optional :: FillValue, Offset, Factor

      call check_err(nf90_def_var(ncid, Sname, NFtype, dims, var_id))
      call check_err(nf90_put_att(ncid, var_id, 'long_name', Lname))
      call check_err(nf90_put_att(ncid, var_id, 'units', Units))
      
      if (present(FillValue)) then
         call Check_err(nf90_put_att(ncid, var_ID, &
                        '_FillValue',FillValue))
        !call Check_err(nf90_put_att(ncid, var_ID,'scale_factor',Factor))
        !call Check_err(nf90_put_att(ncid, var_ID,'add_offset',Offset))
      endif
!
      end subroutine def_var_att       
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     S U B R O U T I N E    I N I T _ N E T C D F
!-----------------------------------------------------------------------
      subroutine UpdateNetCDF(IT,OKflag)
      implicit none
      integer, intent(in) :: IT
      logical, intent(in) :: OKflag
      integer,dimension(3) :: start, kount, kountX, kountY, kountYY
      integer  :: startT, mpistat(mpi_status_size)
    
      ! Make netcdf/get dimension and variable IDs 
      if (IT.eq.0) then 
         call initNetCDF()
      endif

      !start = [1, 1, mnProc*IT + myProc + tsind];
      kount = [NX, NY, 1];
      kountY = [NX, NY-1, 1];
      kountYY = [NX, NY-2, 1];
 
      if ( (IT > 0 .or. myProc > 0) .and. mnProc > 1) then
         ! Receiving message from prev processor to start writing
         call mpi_recv(startT,1,mpi_integer,prevProc,0,MPI_Comm_World,&
                       mpistat,ierr)
      endif

      if (OKflag) then
         if ((IT.eq.0.and.myProc.eq.0).or.mnProc.eq.1) then
            startT = tsind + IT
         else
            startT = startT + 1
         endif
         start = [1, 1, startT];
         write(6,*) 'UpdateNetCDF:',myProc,&
                    CurDT%strftime("%Y-%m-%d %H:%M")
         call check_err(nf90_open(BC2D_Name, nf90_write, ncid))
         call check_err(nf90_put_var(ncid, timenc_id,    &
                        CurDT%strftime("%Y-%m-%d %H:%M"),&
                        [1, start(3)],[16, kount(3)]))
         if (OutType.gt.0) then
            call put_var(ncid, BPGX_id,BC2D_BX, start, kount)
            call put_var(ncid, BPGY_id, BC2D_BY, start, kountY)
            call put_var(ncid, SigTS_id, BC2D_SigTS, start, kount)
            call put_var(ncid, MLD_id, BC2D_MLD, start, kount)
         endif
         call put_var(ncid, NB_id, BC2D_NB, start, kount)
         call put_var(ncid, NM_id, BC2D_NM, start, kount)
         if (OutType.eq.2) then
            call put_var(ncid, KE_id, BC2D_KE, start, kount)
            call put_var(ncid, CD_id, BC2D_CD, start, kountYY)
            call put_var(ncid, DispX_id, BC2D_DispX, start, kountYY)
            call put_var(ncid, DispY_id, BC2D_DispY, start, kountYY)
         endif
         call check_err(nf90_close(ncid))
      else
         write(6,*) 'Skipping Update of NetCDF:',&
                     myProc,CurDT%strftime("%Y-%m-%d %H:%M")
      endif

      if (mnProc > 1) then
         ! Telling next processor to start writing 
         call mpi_send(startT,1,mpi_integer,nextProc,0,&
                       MPI_Comm_World,ierr)
      endif
 
      end subroutine UpdateNetCDF
!-----------------------------------------------------------------------
!     S U B R O U T I N E    P U T _ V A R
!-----------------------------------------------------------------------
      subroutine put_var(ncid, var_id, InArray, start, kount)
      implicit none
      integer,intent(in)  :: ncid, var_id, start(3), kount(3) 
      real(sz),intent(in) :: InArray(:,:)
      real(sz),allocatable :: Temp2D(:,:)

      if (is > 1) then
         ! If 0/360, re-orientate
         allocate(Temp2D(kount(1),kount(2)))
         Temp2D(1:kount(1)-is+1,:) = InArray(is:kount(1),:)
         Temp2D(kount(1)-is:kount(1),:)  = InArray(1:is-1,:)
         ! Write out into netcdf
         call check_err(nf90_put_var(ncid, var_id, Temp2D,  &
                        start, kount))
      else
         ! Write outinto netcdf
         call check_err(nf90_put_var(ncid, var_id, InArray, &
                        start, kount))
      endif
!
      end subroutine put_var       
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      END PROGRAM OGCM_DL
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
