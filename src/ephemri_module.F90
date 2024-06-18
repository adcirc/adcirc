!
! By Al Cerrone,
! Jun 2024.
!
module ephemri_module

#ifdef ADCNETCDF         
  use netcdf
#endif

  use astroformod, only: km2AU 

  implicit none

  private::  interpolate, L, adjust_RA, check

contains

  subroutine HEAVENLY_OBJS_COORDS_FROM_TABLE( MoonSunCoor, julian_date_loc, MoonSunCoordFile, IERR )
    ! use netcdf
    implicit none

    real(8), intent(in) :: julian_date_loc
    character(len=*), intent(in) :: MoonSunCoordFile
    real(8):: MoonSunCoor(3, 2)
    INTEGER:: IERR

    integer :: ncid, time_varid, lunar_distance_varid, solar_distance_varid
    integer :: lunar_ra_varid, solar_ra_varid, lunar_dec_varid, solar_dec_varid
    integer :: retval, j, len_times
    integer, parameter :: NC_ERR = -1
    real(8), allocatable :: times(:), lunar_distances(:), solar_distances(:)
    real(8), allocatable :: lunar_ras(:), solar_ras(:), lunar_decs(:), solar_decs(:)
    real(8) :: julian_datetime_2000, seconds_between

    real(8) :: ratmp(4) 
    logical :: match
    integer :: time_dimid, ndimsp, dimlen

#ifndef ADCNETCDF
    WRITE(*,*) " Error:"
    WRITE(*,*) "   To include this subroutine,"
    WRITE(*,*) "    1. Compile this file witth -DADCNETCDF in compiler flags."
    WRITE(*,*) "    2. A fortran netcdf library is requried."  
    WRITE(*,*) "    3. Include the location of the netcdf fortran module file"// &
                       "(-I$(netcdf_include_path))"//" to complie this file" 
    WRITE(*,*) "    4. Link the netcdf library when building an execubale"//
                      "(-L$(netcdf_library_path) -lnetcdf -lnetcdff."
    ierr = 1 ;

    RETURN ;     
#else
    integer:: dimids(nf90_max_dims) 

    ! Open the NetCDF file
    call check(nf90_open(MoonSunCoordFile, nf90_nowrite, ncid))

    ! Get variable IDs
    call check(nf90_inq_varid(ncid, 'time', time_varid))
    call check(nf90_inq_varid(ncid, 'lunar_distance', lunar_distance_varid))
    call check(nf90_inq_varid(ncid, 'solar_distance', solar_distance_varid))
    call check(nf90_inq_varid(ncid, 'lunar_ra', lunar_ra_varid))
    call check(nf90_inq_varid(ncid, 'solar_ra', solar_ra_varid))
    call check(nf90_inq_varid(ncid, 'lunar_dec', lunar_dec_varid))
    call check(nf90_inq_varid(ncid, 'solar_dec', solar_dec_varid))

    ! Inquire about the time variable to get dimension IDs
    call check(nf90_inquire_variable(ncid, time_varid, ndims=ndimsp, dimids=dimids))
    
    ! Get the length of the time dimension
    call check(nf90_inquire_dimension(ncid, dimids(1), len=dimlen))
    len_times = dimlen
    
    ! Allocate arrays
    allocate(times(len_times), lunar_distances(len_times), solar_distances(len_times))
    allocate(lunar_ras(len_times), solar_ras(len_times), lunar_decs(len_times), solar_decs(len_times))

    ! Read variables
    call check(nf90_get_var(ncid, time_varid, times))
    call check(nf90_get_var(ncid, lunar_distance_varid, lunar_distances))
    call check(nf90_get_var(ncid, solar_distance_varid, solar_distances))
    call check(nf90_get_var(ncid, lunar_ra_varid, lunar_ras))
    call check(nf90_get_var(ncid, solar_ra_varid, solar_ras))
    call check(nf90_get_var(ncid, lunar_dec_varid, lunar_decs))
    call check(nf90_get_var(ncid, solar_dec_varid, solar_decs))

    ! Close the NetCDF file
    call check(nf90_close(ncid))

    ! Calculate seconds between the provided date and the reference date
    julian_datetime_2000 = 2451544.5d0
    seconds_between = (julian_date_loc - julian_datetime_2000) * 86400.0d0

    match = .false.

    do j = 1, len_times
      if (times(j) > seconds_between) then

        ! Lunar distance interpolation
        if (j >= 3 .and. j <= len_times - 2) then
          match = .true.
          call interpolate(times(j-2:j+1), lunar_distances(j-2:j+1), seconds_between, MoonSunCoor(3,1))
        else if (j == 2) then
          match = .true.
          call interpolate(times(1:2), lunar_distances(1:2), seconds_between, MoonSunCoor(3,1))
        else if (j == len_times) then
          match = .true.
          call interpolate(times(len_times-1:len_times), lunar_distances(len_times-1:len_times), seconds_between, MoonSunCoor(3,1))
        end if

        ! Solar distance interpolation
        if (j >= 3 .and. j <= len_times - 2) then
          match = .true.
          call interpolate(times(j-2:j+1), solar_distances(j-2:j+1), seconds_between, MoonSunCoor(3,2))
        else if (j == 2) then
          match = .true.
          call interpolate(times(1:2), solar_distances(1:2), seconds_between, MoonSunCoor(3,2))
        else if (j == len_times) then
          match = .true.
          call interpolate(times(len_times-1:len_times), solar_distances(len_times-1:len_times), seconds_between, MoonSunCoor(3,2))
        end if

        ! Lunar right ascension interpolation
        ratmp = 0.D0 ;
        if (j >= 3 .and. j <= len_times - 2) then
          match = .true.

          ratmp(1:4) = lunar_ras(j-2:j+1) ;
          call adjust_RA( ratmp(1:4) ) ;

!!          call interpolate(times(j-2:j+1), lunar_ras(j-2:j+1), seconds_between, MoonSunCoor(1,1))
          call interpolate(times(j-2:j+1), ratmp(1:4), seconds_between, MoonSunCoor(1,1))
        else if (j == 2) then
          match = .true.

          ratmp(1:2) = lunar_ras(1:2) ;
          call adjust_RA( ratmp(1:2) ) ;

!!          call interpolate(times(1:2), lunar_ras(1:2), seconds_between, MoonSunCoor(1,1))
          call interpolate(times(1:2), ratmp(1:2), seconds_between, MoonSunCoor(1,1))
        else if (j == len_times) then
          match = .true.

          ratmp(1:2) = lunar_ras(len_times-1:len_times) ;
          call adjust_RA( ratmp(1:2) ) ;

!!          call interpolate(times(len_times-1:len_times), lunar_ras(len_times-1:len_times), seconds_between, MoonSunCoor(1,1))
          call interpolate(times(len_times-1:len_times), ratmp(1:2), seconds_between, MoonSunCoor(1,1))
        end if
        if ( match ) then
            MoonSunCoor(1,1) = modulo( MoonsunCoor(1,1), 360.D0 ) ; 
        end if
       
        ! Solar right ascension interpolation
        ratmp = 0.D0 ;
        if (j >= 3 .and. j <= len_times - 2) then
          match = .true.

          ratmp(1:4) = solar_ras(j-2:j+1) ;
          call adjust_RA( ratmp(1:4) ) ;

!!          call interpolate(times(j-2:j+1), solar_ras(j-2:j+1), seconds_between, MoonSunCoor(1,2))
          call interpolate(times(j-2:j+1), ratmp(1:4), seconds_between, MoonSunCoor(1,2))
        else if (j == 2) then
          match = .true.

          ratmp(1:2) = solar_ras(1:2) ;
          call adjust_RA( ratmp(1:2) ) ;

!!          call interpolate(times(1:2), solar_ras(1:2), seconds_between, MoonSunCoor(1,2))
          call interpolate(times(1:2), ratmp(1:2), seconds_between, MoonSunCoor(1,2))
        else if (j == len_times) then
          match = .true.

          ratmp(1:2) = solar_ras(len_times-1:len_times) ;
          call adjust_RA( ratmp ) ; 
          
!!          call interpolate(times(len_times-1:len_times), solar_ras(len_times-1:len_times), seconds_between, MoonSunCoor(1,2))
          call interpolate(times(len_times-1:len_times), ratmp(1:2), seconds_between, MoonSunCoor(1,2))
        end if
        if ( match ) then
            MoonSunCoor(1,2) = modulo( MoonsunCoor(1,2), 360.D0 ) ; 
        end if

        ! Lunar declination interpolation
        if (j >= 3 .and. j <= len_times - 2) then
          match = .true.
          call interpolate(times(j-2:j+1), lunar_decs(j-2:j+1), seconds_between, MoonSunCoor(2,1))
        else if (j == 2) then
          match = .true.
          call interpolate(times(1:2), lunar_decs(1:2), seconds_between, MoonSunCoor(2,1))
        else if (j == len_times) then
          match = .true.
          call interpolate(times(len_times-1:len_times), lunar_decs(len_times-1:len_times), seconds_between, MoonSunCoor(2,1))
        end if

        ! Solar declination interpolation
        if (j >= 3 .and. j <= len_times - 2) then
          match = .true.
          call interpolate(times(j-2:j+1), solar_decs(j-2:j+1), seconds_between, MoonSunCoor(2,2))
        else if (j == 2) then
          match = .true.
          call interpolate(times(1:2), solar_decs(1:2), seconds_between, MoonSunCoor(2,2))
        else if (j == len_times) then
          match = .true.
          call interpolate(times(len_times-1:len_times), solar_decs(len_times-1:len_times), seconds_between, MoonSunCoor(2,2))
        end if

        exit
      end if
    end do

    IERR = 0 ; 
    if (.not. match) then
      ! write(*,*) 'Error: Date not within database.'
      IERR = 2 ; 
      !stop
    end if

    ! Convert solar distance to AU
    ! MoonSunCoor(3,2) = MoonSunCoor(3,2) * 6.6845871226706E-9
    MoonSunCoor(3,2) = MoonSunCoor(3,2) * km2AU ; 
#endif

  end subroutine HEAVENLY_OBJS_COORDS_FROM_TABLE

  ! check whether RA change cross the zero hr line.
  ! if so, adjust accordingly                  
  subroutine adjust_RA( RA, fadjust )
    implicit none

    LOGICAL, optional:: fadjust
    REAL (8), INTENT(INOUT):: RA(:)

    INTEGER:: I, N
    LOGICAL:: CrossZeroRA 

    N = ubound( RA, 1 ) ; 

    CrossZeroRA = .FALSE. ;
    DO I = 1, N - 1
      IF ( abs(RA(I+1) - RA(I)) > 120.D0 ) THEN
        CrossZeroRA = .TRUE. ;
      END IF
    END DO

    IF ( CrossZeroRA ) THEN
      WHERE ( RA < 180.D0 ) RA = RA + 360.D0 ;        
    END IF

    IF ( present(fadjust) ) fadjust = CrossZeroRA ;  

    RETURN ; 
  end subroutine adjust_RA


  ! Subroutine to perform Lagrange interpolation
  subroutine interpolate(x, y, xi, yi)
    real(8), intent(in) :: x(:), y(:), xi
    real(8), intent(out) :: yi
    integer :: k
    real(8) :: term
    yi = 0.0d0
    do k = 1, size(x)
      call L(k, xi, x, term)
      yi = yi + y(k) * term
    end do
  end subroutine interpolate

  ! Subroutine to calculate Lagrange polynomial term
  subroutine L(k, xi, x, term)
    integer, intent(in) :: k
    real(8), intent(in) :: xi, x(:)
    real(8), intent(out) :: term
    integer :: i
    term = 1.0d0
    do i = 1, size(x)
      if (i /= k) then
        term = term * (xi - x(i)) / (x(k) - x(i))
      end if
    end do
  end subroutine L

  ! Subroutine for error checking
  subroutine check(status)
    integer, intent(in) :: status
 
#ifdef ADCNETCDF    
    if (status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
#endif

  end subroutine check

end module ephemri_module
