!>
!! @mainpage ADCIRC NUOPC Cap
!! @author Saeed Moghimi (moghimis@gmail.com)
!! @date 15/1/17 Original documentation
!------------------------------------------------------
!LOG-----------------
!
!
!


      MODULE hwrf_mod
      USE netcdf
      use esmf

      IMPLICIT NONE
      real(ESMF_KIND_R8), allocatable     :: LONS(:), LATS(:),TIMES(:)
      real(ESMF_KIND_R8), allocatable     :: UGRD10(:,:,:), VGRD10(:,:,:)
      real(ESMF_KIND_R8), allocatable     :: PRMSL (:,:,:)
      integer               :: nlat, nlon, ntime
      character (len = 280) :: FILE_NAME

      ! reading data time management info 
      integer               :: wrf_int, wrf_num, wrf_den
      character (len = 280) :: wrf_dir, wrf_nam
!-----------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------
!- Sub !!!????
!-----------------------------------------------------------------------
    SUBROUTINE read_hwrf_nc()
      IMPLICIT NONE
      character (len = *), parameter :: LAT_NAME    = "latitude"
      character (len = *), parameter :: LON_NAME    = "longitude"
      character (len = *), parameter :: REC_NAME    = "time"
      character (len = *), parameter :: PRMSL_NAME  = "PRMSL_meansealevel"
      character (len = *), parameter :: UGRD10_NAME = "UGRD_10maboveground"
      character (len = *), parameter :: VGRD10_NAME = "VGRD_10maboveground"
      integer :: ncid, LON_dimid, LAT_dimid, rec_dimid
      integer :: PRMSL_dimid, UGRD10_dimid, VGRD10_dimid

      integer, parameter :: NDIMS = 3
      integer :: LON_varid, LAT_varid, rec_varid
      integer :: PRMSL_varid, UGRD10_varid, VGRD10_varid
      integer :: start(NDIMS),count(NDIMS)
      logical :: THERE
      !integer :: dimids(NDIMS)

      integer :: lat, lon,i, iret

      print *, FILE_NAME
      INQUIRE( FILE= FILE_NAME, EXIST=THERE ) 
      if ( .not. THERE)  stop 'HWRF file does not exist !'

      ncid = 0
      ! Open the file.
      call check(  nf90_open(trim(FILE_NAME), NF90_NOWRITE, ncid))
      
      ! Get ID of unlimited dimension
      call check( nf90_inquire(ncid, unlimitedDimId = rec_dimid) )

      ! Get ID of limited dimension
      call check( nf90_inq_dimid(ncid, LON_NAME, LON_dimid) )
      call check( nf90_inq_dimid(ncid, LAT_NAME, LAT_dimid) )

      ! How many values of "lat" are there?
      call check(nf90_inquire_dimension(ncid, lon_dimid, len = nlon) )
      call check(nf90_inquire_dimension(ncid, lat_dimid, len = nlat) )

      ! What is the name of the unlimited dimension, how many records are there?
      call check(nf90_inquire_dimension(ncid, rec_dimid, len = ntime))

      ! Get the varids of the pressure and temperature netCDF variables.
      call check( nf90_inq_varid(ncid, LAT_NAME,    LAT_varid) )
      call check( nf90_inq_varid(ncid, LON_NAME,    LON_varid) )
      call check( nf90_inq_varid(ncid, REC_NAME,    rec_varid) )
      call check( nf90_inq_varid(ncid, PRMSL_NAME,  PRMSL_varid) )
      call check( nf90_inq_varid(ncid, UGRD10_NAME, UGRD10_varid) )
      call check( nf90_inq_varid(ncid, VGRD10_NAME, VGRD10_varid) )
      !print *, ntime,nlat,nlon
      !allocate vars
      if(.not. allocated(LATS))   allocate (LATS  (1:nlat))
      if(.not. allocated(LONS))   allocate (LONS  (1:nlon))
      if(.not. allocated(TIMES))  allocate (TIMES (1:ntime))
      if(.not. allocated(UGRD10)) allocate (UGRD10(nlon,nlat,ntime))
      if(.not. allocated(VGRD10)) allocate (VGRD10(nlon,nlat,ntime))
      if(.not. allocated(PRMSL))  allocate (PRMSL (nlon,nlat,ntime))
      ! read vars
      call check(nf90_get_var(ncid, LAT_varid, LATS))
      call check(nf90_get_var(ncid, LON_varid, LONS))
      call check(nf90_get_var(ncid, rec_varid, TIMES))
      !TODO: Why the order is other way???? Might change the whole forcing fields!!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< IMPORTANT <<<<<
      !TODO: plot input and output<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      start = (/ 1, 1, 1/)
      count = (/nlon,nlat,ntime/)
      print *, 'nlon > ',nlon,' nlat > ',nlat,' ntime > ',ntime
      print *, start+count
      !!!!                                                              print *,size(PRMSL(:,nlat,ntime))
      call check( nf90_get_var(ncid,PRMSL_varid,   PRMSL, start, count) )
      call check( nf90_get_var(ncid,UGRD10_varid, UGRD10, start, count) )
      call check( nf90_get_var(ncid,VGRD10_varid, VGRD10, start, count) )

      print *, FILE_NAME

      !do i=1,nlat
      !print *, 'PRMSL (:,ilat,1) > ','ilat = ',i,'    ', PRMSL(:,i,1)
      !end do
    !
    END SUBROUTINE
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine check(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop 
    end if
  end subroutine check
  
  
!-------------------------------------------------------------------------
    FUNCTION Replace_Text (s,text,rep)  RESULT(outs)
    CHARACTER(*)        :: s,text,rep
    CHARACTER(LEN(s)+300) :: outs     ! provide outs with extra 100 char len
    INTEGER             :: i, nt, nr

    outs = s ; nt = LEN_TRIM(text) ; nr = LEN_TRIM(rep)
        DO
           i = INDEX(outs,text(:nt)) ; IF (i == 0) EXIT
           outs = outs(:i-1) // rep(:nr) // outs(i+nt:)
        END DO
    END FUNCTION Replace_Text  
  
 
    subroutine update_hwrf_filename (YY, MM, DD, H)
        integer             :: YY, MM, DD, H
        CHARACTER(len=280)      :: inps     ! provide outs with extra 100 char len
        CHARACTER(len=4)        :: year
        CHARACTER(len=2)        :: mon,day
        CHARACTER(len=3)        :: hours

        ! example:  wrf_nam:  andrew04l.YYYYMMDD00.hwrfprs.d123.0p06.fHHH.grb2.nc
        inps = trim(wrf_nam)
        
        write(year,"(I4.4)") YY
        inps =  Replace_Text (inps,'YYYY',year)
      
        write(mon,"(I2.2)")  MM
        inps =  Replace_Text (inps,'MM',mon)    

        write(day,"(I2.2)")  DD
        inps =  Replace_Text (inps,'DD',day)    
        
        !past hours from start date
        write(hours,"(I3.3)") H
        inps =  Replace_Text (inps,'HHH',hours)        

        FILE_NAME =  TRIM(wrf_dir)//'/'//TRIM(inps)
        

    END subroutine update_hwrf_filename 
 
 
     subroutine read_config()

    character(ESMF_MAXPATHLEN)    :: fname ! config file name
    type(ESMF_Config)             :: cf     ! the Config itself
    integer                       :: rc

    rc = ESMF_SUCCESS
    
   !Initiate reading resource file
    cf = ESMF_ConfigCreate(rc=rc)  ! Create the empty Config
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    fname = "config.rc" ! Name the Resource File
    call ESMF_ConfigLoadFile(cf, fname, rc=rc) ! Load the Resource File
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

   ! read time hwrf data time interval info
   ! call ESMF_ConfigGetAttribute(cf, wrf_int, label="wrf_int:",default=3600,rc=rc)
   ! call ESMF_ConfigGetAttribute(cf, wrf_num, label="wrf_num:",default=0  , rc=rc)
   ! call ESMF_ConfigGetAttribute(cf, wrf_den, label="wrf_den:",default=1  , rc=rc)
   ! if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
   !   line=__LINE__, &
   !   file=__FILE__)) &
   !   return  ! bail out
      
    call ESMF_ConfigGetAttribute(cf, wrf_dir, label="wrf_dir:",default='hwrf_data'  , rc=rc)
    call ESMF_ConfigGetAttribute(cf, wrf_nam, label="wrf_nam:", &
         default='andrew04l.YYYYMMDD00.hwrfprs.d123.0p06.f0HH.grb2.nc'  , rc=rc)    
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out      
        
    call ESMF_ConfigDestroy(cf, rc=rc) ! Destroy the Config
        
    end subroutine read_config

  

END MODULE
!-----------------------------------------
