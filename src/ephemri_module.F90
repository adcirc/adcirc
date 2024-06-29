!
! By Al Cerrone,
! Jun 2024.
!
module ephemri_module

#ifdef ADCNETCDF         
  use netcdf
#endif

#ifndef TIPSTANDALONE
  use global, only: RNDAY
#endif

  use astroformod, only: km2AU 

  implicit none

  INTEGER, PRIVATE:: len_times = 0 
  REAL(8), PRIVATE:: tbeg(2) = 0.D0, tend(2) = 0.D0  
  REAL(8), PRIVATE, ALLOCATABLE, DIMENSION(:):: times, lunar_distances, &
       & solar_distances, lunar_ras, solar_ras, lunar_decs, solar_decs   

  REAL(8), PRIVATE:: tcache = 2592000.D0 ; ! cache 30 days of data for a faster retrival when
                                           ! HAVENLY_OBJS_COORDS_FROM_TABLE is called multiple 
                                           ! times in a single program  

  private::  interpolate, L, adjust_RA, check, GET_RANK_SIMPLE_SEARCH, &
          GET_RANK_UNIFORM, GET_RANK_BINARY_SEARCH 


contains

  subroutine HEAVENLY_OBJS_COORDS_FROM_TABLE( MoonSunCoor, julian_date_loc, MoonSunCoordFile, IERR, UniformDT  )
    ! use netcdf
    implicit none

    real(8), intent(in) :: julian_date_loc
    character(len=*), intent(in) :: MoonSunCoordFile
    real(8):: MoonSunCoor(3, 2)
    INTEGER:: IERR
    LOGICAL, optional:: UniformDT 

    integer:: ncid, time_varid, lunar_distance_varid, solar_distance_varid
    integer:: lunar_ra_varid, solar_ra_varid, lunar_dec_varid, solar_dec_varid
    integer :: retval, j, ii
    integer, parameter :: NC_ERR = -1
!    integer::  len_times
!    real(8), allocatable :: times(:), lunar_distances(:), solar_distances(:)
!    real(8), allocatable :: lunar_ras(:), solar_ras(:), lunar_decs(:), solar_decs(:)
    real(8) :: julian_datetime_2000, seconds_between

    LOGICAL:: UniformRankSearch = .false. ; 
    INTEGER:: lenarr, iabeg, iaend
    real(8), ALLOCATABLE :: tmparr(:)

    real(8) :: ratmp(4) 
    logical :: match
    integer :: time_dimid, ndimsp, dimlen
    integer :: idarr(4) 

    logical:: first = .true., recache = .false.  

#ifndef ADCNETCDF
    WRITE(*,*) " Error:"
    WRITE(*,*) "   To include this subroutine,"
    WRITE(*,*) "    1. Compile this file witth -DADCNETCDF in compiler flags."
    WRITE(*,*) "    2. A netcdf and fortran netcdf library are requried."  
    WRITE(*,*) "    3. Include the location of the netcdf fortran module file"// &
                       "(-I$(netcdf_include_path))"//" to complie this file" 
    WRITE(*,*) "    4. Link the netcdf library when building an execubale"//
                      "(-L$(netcdf_library_path) -lnetcdf -lnetcdff."
    ierr = 1 ;

    RETURN ;     
#else
    integer:: dimids(nf90_max_dims) 

    ! Calculate seconds between the provided date and the reference date !
    julian_datetime_2000 = 2451544.5d0
    seconds_between = (julian_date_loc - julian_datetime_2000) * 86400.0d0

#ifndef TIPSTANDALONE
    ! For ADCIRC, keep a portion of ephemeride data that covers an entire run period  
    tcache = RNDAY*86400.0D0 
#endif

    IF ( present(UniformDT) ) THEN
      UniformRankSearch = UniformDT ;
    END IF   

    IF ( .not. first ) THEN
      IF ( seconds_between > times(len_times - 4) ) THEN
        recache = .true. ; 
      END IF         
    END IF

    IF ( first .OR. recache ) THEN    
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
      
      lenarr = dimlen ;
      allocate( tmparr(lenarr) ) ; tmparr = 0.D0 ;
      call check(nf90_get_var( ncid, time_varid, tmparr )) ;

      tbeg(2) = tmparr(1) ;
      tend(2) = tmparr(lenarr) ; 

      IF ( seconds_between < tbeg(2) .OR. seconds_between > tend(2) ) THEN
        IERR = 2 ;
        return ;
      END IF

      IF ( .not. UniformRankSearch ) THEN
        iabeg = GET_RANK_BINARY_SEARCH( seconds_between, tmparr, lenarr ) - 4 ;    
      ELSE
        iabeg = GET_RANK_UNIFORM( seconds_between, tmparr, lenarr) - 4 ;
      END IF        
      if ( iabeg < 1 ) iabeg = 1 ; 

      IF ( .not. UniformRankSearch ) THEN
        iaend = GET_RANK_BINARY_SEARCH( seconds_between + tcache, tmparr, lenarr ) + 4 ;
      ELSE  
        iaend = GET_RANK_UNIFORM( seconds_between + tcache, tmparr, lenarr ) + 4 ;
      END IF        
      if ( iaend > lenarr ) iaend = lenarr ; 

      tbeg(1) = tmparr(iabeg) ; 
      tend(1) = tmparr(iaend) ; 

      ! Allocate arrays !
      len_times = (iaend - iabeg) + 1 ; 
      if ( allocated( times ) ) then
        deallocate( times ) ; 
        deallocate( lunar_distances, solar_distances ) ;
        deallocate( lunar_ras, solar_ras ) ; 
        deallocate( lunar_decs, solar_decs ) ; 
      endif    
      allocate(times(len_times), lunar_distances(len_times), solar_distances(len_times))
      allocate(lunar_ras(len_times), solar_ras(len_times), lunar_decs(len_times), solar_decs(len_times))
  
      times = tmparr(iabeg:iaend) ; 
      call check( nf90_get_var( ncid, lunar_distance_varid, lunar_distances, (/ iabeg /), (/ len_times /) ))
      call check( nf90_get_var( ncid, solar_distance_varid, solar_distances, (/ iabeg /), (/ len_times /) ))
      call check( nf90_get_var( ncid, lunar_ra_varid, lunar_ras, (/ iabeg /), (/ len_times /) ))
      call check( nf90_get_var( ncid, solar_ra_varid, solar_ras, (/ iabeg /), (/ len_times /) ))
      call check( nf90_get_var( ncid, lunar_dec_varid, lunar_decs, (/ iabeg /), (/ len_times /) ))
      call check( nf90_get_var( ncid, solar_dec_varid, solar_decs, (/ iabeg /), (/ len_times /) ))
  
      ! Close the NetCDF file
      call check(nf90_close(ncid))
      deallocate(tmparr) ;

      first = .false. ;
      if ( recache ) recache = .false. ;
    END IF

    IF ( .not. UniformRankSearch ) THEN
      J = GET_RANK_BINARY_SEARCH( seconds_between, times, len_times ) ;
    ELSE
      J = GET_RANK_UNIFORM( seconds_between, times, len_times ) ; 
    END IF    

    match = .false. ;
    if ( j >= 3 .and. j <= len_times - 2 ) then
       match = .true. ;
       idarr = (/ ( j - ii, ii = -2, 1 ) /)  ;
    else if ( j == 2 ) then
       match = .true. ;
       idarr = (/ (ii, ii = 1, 4) /)  ; 
    else if ( j == len_times ) then
       match = .true. 
       idarr = (/ (len_times - 3 + ii, ii = 0, 3 ) /) ;
    end if     

    if ( match ) then
       ! Lunar distance interpolation
       call interpolate(times(idarr), lunar_distances(idarr), seconds_between, MoonSunCoor(3,1))

       ! Solar distance interpolation   
       call interpolate(times(idarr), solar_distances(idarr), seconds_between, MoonSunCoor(3,2))
      
       ! Lunar right ascension interpolation
       ratmp(1:4) = lunar_ras(idarr) ;
       call adjust_RA( ratmp(1:4) ) ;
       call interpolate( times(idarr), ratmp(1:4), seconds_between, MoonSunCoor(1,1) )
       MoonSunCoor(1,1) = modulo( MoonsunCoor(1,1), 360.D0 ) ; 

       ! Solar right ascension interpolation
       ratmp(1:4) = solar_ras(idarr) ;
       call adjust_RA( ratmp(1:4) ) ;
       call interpolate(times(idarr), ratmp(1:4), seconds_between, MoonSunCoor(1,2))
       MoonSunCoor(1,2) = modulo( MoonsunCoor(1,2), 360.D0 ) ; 

       ! Lunar declination interpolation
       call interpolate(times(idarr), lunar_decs(idarr), seconds_between, MoonSunCoor(2,1))

       ! Solar declination interpolation
       call interpolate(times(idarr), solar_decs(idarr), seconds_between, MoonSunCoor(2,2))
    end if   

    if (.not. match) then
      IERR = 2 ; 
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

  ! Return a rank j in ARR such that arr(j-1) < val <= arr(j1) !  
  FUNCTION GET_RANK_SIMPLE_SEARCH( val, arr, len ) RESULT( rank ) 
    IMPLICIT NONE

    INTEGER:: rank 
    REAL(8), INTENT(IN):: val
    REAL(8), INTENT(IN):: arr(:)
    INTEGER, INTENT(IN):: len

    INTEGER:: J

    IF ( val < arr(1) ) THEN 
       rank = 0  ;
       return ;
    END IF

    IF ( val > arr(len) ) THEN
       rank = len + 1  ;
       return ; 
    END IF

    RANK = 0 ; 
    DO J = 1, LEN
       IF ( ARR(J) >= val ) THEN
         RANK = J ;
         EXIT ;
       END IF
    END DO
    
  END FUNCTION GET_RANK_SIMPLE_SEARCH     

  ! Find a rank j in arr such that arr(j-1) < val < = arr(j)         !        
  ! in a unifrom table, i,e, arr(i+1) - arr(i) = arr(i+2) - arr(i+1) !
  FUNCTION GET_RANK_UNIFORM( val, arr, len ) RESULT( rank )
     IMPLICIT NONE

     INTEGER:: rank
     REAL(8), INTENT(IN):: val 
     REAL(8), INTENT(IN):: arr(:)
     INTEGER, INTENT(IN):: len

     REAL (8):: ds

     IF ( val < arr(1) ) THEN 
        rank = 0  ;
        return ;
     END IF
 
     IF ( val > arr(len) ) THEN
        rank = len + 1  ;
        return ; 
     END IF
 
     ds = (arr(2) - arr(1)) ; 
     rank = ceiling( (val - arr(1))/ds ) + 1 ;

  END FUNCTION GET_RANK_UNIFORM

  ! Find a rank j in arr such that arr(j-1) < val < = arr(j)  !        
  ! using a binary search                                     !
  FUNCTION GET_RANK_BINARY_SEARCH( val, arr, len ) RESULT( rank )
     IMPLICIT NONE

     INTEGER:: rank
     REAL(8), INTENT(IN):: val 
     REAL(8), INTENT(IN):: arr(:)
     INTEGER:: len

     INTEGER:: low, high, mid

     IF ( val < arr(1) ) THEN 
        rank = 0  ;
        return ;
     END IF
     IF ( val > arr(len) ) THEN
        rank = len + 1  ;
        return ; 
     END IF

     low = 1 ;
     high = len ; 
     DO
       mid = (low + high)/2 ; 
       IF ( abs(val - arr(mid)) < 1.0e-12 ) THEN
           ! Exact match is found !
           rank = mid ; 
           EXIT ;
       END IF
 
       IF ( val < arr(mid) ) THEN
           high = mid ;    
       ELSE     
           low = mid ;
       END IF
 
       IF ( high == low + 1 ) THEN
          ! found the interval to which val  !
          ! falls inbetween is found         !     
          rank = high ;
          EXIT ;
       END IF
 
       IF ( high < low ) THEN
         WRITE(*,*) "Error in GET_RANK_BINARY_SEARCH()" ; 
         rank = -1 ;
         EXIT ; 
       END IF
    END DO
                 
  END FUNCTION GET_RANK_BINARY_SEARCH       

end module ephemri_module
