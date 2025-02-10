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
!
! By Al Cerrone,
! Jun 2024.
!
module mod_ephemerides

  implicit none

  type t_ephemerides
    private
    logical :: first = .true.
    integer :: len_times
    character(len=256) :: MoonSunCoordFile = "none"
    real(8) :: tbeg(2) = (/0d0, 0d0/)
    real(8) :: tend(2) = (/0d0, 0d0/)
    real(8), allocatable :: times(:), lunar_distances(:), solar_distances(:)
    real(8), allocatable :: lunar_ras(:), solar_ras(:), lunar_decs(:), solar_decs(:)
    real(8) :: tcache = 2592000.d0 ! cache 30 days of data for a faster retrival when
    ! HAVENLY_OBJS_COORDS_FROM_TABLE is called multiple
    ! times in a single program
  contains
#ifdef ADCNETCDF
    procedure, pass(self) :: HEAVENLY_OBJS_COORDS_FROM_TABLE
    procedure, pass(self), private :: recache_data
#endif
    procedure, pass(self), private :: reallocate_arrays
    procedure, pass(self), private :: init => initialize_ephemerides
  end type t_ephemerides

  interface t_ephemerides
    module procedure createEphemerides
  end interface t_ephemerides

  private

  public :: t_ephemerides

#ifdef ADCNETCDF
  public :: HEAVENLY_OBJS_COORDS_FROM_TABLE
#endif

contains

  subroutine initialize_ephemerides(self, in_rnday, in_MoonSunCoordFile)
    implicit none
    class(t_ephemerides), intent(inout) :: self
    real(8), intent(in) :: in_rnday
    character(len=*), intent(in) :: in_MoonSunCoordFile

    self%first = .true.
    self%len_times = 0
    self%MoonSunCoordFile = in_MoonSunCoordFile

    ! Keep a portion of ephemeride data that covers an entire run period
    self%tcache = in_rnday*86400.0d0

  end subroutine initialize_ephemerides

  function createEphemerides(rnday, MoonSunCoordFile) result(ephemerides)
    implicit none
    real(8), intent(in) :: rnday
    character(len=*), intent(in) :: MoonSunCoordFile
    type(t_ephemerides) :: ephemerides

    call ephemerides%init(rnday, MoonSunCoordFile)

  end function createEphemerides

#ifdef ADCNETCDF
  subroutine HEAVENLY_OBJS_COORDS_FROM_TABLE(self, MoonSunCoor, julian_date_loc, ierr, UniformDT)
    use mod_astronomic, only: km2AU
#ifdef CMPI
    use messenger, only: msg_fini
#endif
    implicit none
    class(t_ephemerides), intent(inout) :: self
    real(8), intent(in) :: julian_date_loc
    real(8) :: MoonSunCoor(3, 2)
    integer, intent(out) :: ierr
    logical, optional :: UniformDT

    integer :: j, ii
    real(8) :: julian_datetime_2000, seconds_between
    logical :: UniformRankSearch = .false.
    real(8) :: ratmp(4)
    logical :: match
    integer :: idarr(4)

    ierr = 0

    ! Calculate seconds between the provided date and the reference date !
    julian_datetime_2000 = 2451544.5d0
    seconds_between = (julian_date_loc - julian_datetime_2000)*86400.0d0

    if (present(UniformDT)) then
      UniformRankSearch = UniformDT
    end if

    if (.not. self%first) then
      if (seconds_between > self%times(self%len_times - 4)) then
        ierr = self%recache_data(seconds_between, UniformRankSearch)
      end if
    else
      ierr = self%recache_data(seconds_between, UniformRankSearch)
    end if

    if (.not. UniformRankSearch) then
      J = GET_RANK_BINARY_SEARCH(seconds_between, self%times, self%len_times); 
    else
      J = GET_RANK_UNIFORM(seconds_between, self%times, self%len_times); 
    end if

    match = .false.; 
    if (j >= 3 .and. j <= self%len_times - 2) then
      match = .true.; 
      idarr = (/(j - ii, ii=-2, 1)/); 
    else if (j == 2) then
      match = .true.; 
      idarr = (/(ii, ii=1, 4)/); 
    else if (j == self%len_times) then
      match = .true.
      idarr = (/(self%len_times - 3 + ii, ii=0, 3)/); 
    end if

    if (match) then
      ! Lunar distance interpolation
      call interpolate(self%times(idarr), self%lunar_distances(idarr), seconds_between, MoonSunCoor(3, 1))

      ! Solar distance interpolation
      call interpolate(self%times(idarr), self%solar_distances(idarr), seconds_between, MoonSunCoor(3, 2))

      ! Lunar right ascension interpolation
      ratmp(1:4) = self%lunar_ras(idarr); 
      call adjust_RA(ratmp(1:4)); 
      call interpolate(self%times(idarr), ratmp(1:4), seconds_between, MoonSunCoor(1, 1))
      MoonSunCoor(1, 1) = modulo(MoonsunCoor(1, 1), 360.d0); 
      ! Solar right ascension interpolation
      ratmp(1:4) = self%solar_ras(idarr); 
      call adjust_RA(ratmp(1:4)); 
      call interpolate(self%times(idarr), ratmp(1:4), seconds_between, MoonSunCoor(1, 2))
      MoonSunCoor(1, 2) = modulo(MoonsunCoor(1, 2), 360.d0); 
      ! Lunar declination interpolation
      call interpolate(self%times(idarr), self%lunar_decs(idarr), seconds_between, MoonSunCoor(2, 1))

      ! Solar declination interpolation
      call interpolate(self%times(idarr), self%solar_decs(idarr), seconds_between, MoonSunCoor(2, 2))
    end if

    if (.not. match) then
      IERR = 2; 
    end if

    ! Convert solar distance to AU
    ! MoonSunCoor(3,2) = MoonSunCoor(3,2) * 6.6845871226706E-9
    MoonSunCoor(3, 2) = MoonSunCoor(3, 2)*km2AU;

  end subroutine HEAVENLY_OBJS_COORDS_FROM_TABLE

  integer function recache_data(self, seconds_between, UniformRankSearch) result(ierr)

    use netcdf
    use netcdf_error, only: check_err

    implicit none
    class(t_ephemerides), intent(inout) :: self
    real(8), intent(in) :: seconds_between
    logical, intent(in) :: UniformRankSearch

    integer :: ndimsp, dimlen
    integer :: ncid, time_varid, lunar_distance_varid, solar_distance_varid
    integer :: lunar_ra_varid, solar_ra_varid, lunar_dec_varid, solar_dec_varid
    integer :: lenarr, iabeg, iaend
    integer :: dimids(nf90_max_dims)
    real(8), allocatable :: tmparr(:)

    ierr = 0

    ! Open the NetCDF file
    call check_err(nf90_open(self%MoonSunCoordFile, nf90_nowrite, ncid))

    ! Get variable IDs
    call check_err(nf90_inq_varid(ncid, 'time', time_varid))
    call check_err(nf90_inq_varid(ncid, 'lunar_distance', lunar_distance_varid))
    call check_err(nf90_inq_varid(ncid, 'solar_distance', solar_distance_varid))
    call check_err(nf90_inq_varid(ncid, 'lunar_ra', lunar_ra_varid))
    call check_err(nf90_inq_varid(ncid, 'solar_ra', solar_ra_varid))
    call check_err(nf90_inq_varid(ncid, 'lunar_dec', lunar_dec_varid))
    call check_err(nf90_inq_varid(ncid, 'solar_dec', solar_dec_varid))

    ! Inquire about the time variable to get dimension IDs
    call check_err(nf90_inquire_variable(ncid, time_varid, ndims=ndimsp, dimids=dimids))

    ! Get the length of the time dimension
    call check_err(nf90_inquire_dimension(ncid, dimids(1), len=dimlen))

    lenarr = dimlen
    allocate (tmparr(lenarr)); tmparr = 0.d0
    call check_err(nf90_get_var(ncid, time_varid, tmparr))
    self%tbeg(2) = tmparr(1)
    self%tend(2) = tmparr(lenarr)
    if (seconds_between < self%tbeg(2) .or. seconds_between > self%tend(2)) then
      ierr = 2
      return
    end if

    if (.not. UniformRankSearch) then
      iabeg = GET_RANK_BINARY_SEARCH(seconds_between, tmparr, lenarr) - 4; 
    else
      iabeg = GET_RANK_UNIFORM(seconds_between, tmparr, lenarr) - 4; 
    end if

    if (iabeg < 1) iabeg = 1; 
    if (.not. UniformRankSearch) then
      iaend = GET_RANK_BINARY_SEARCH(seconds_between + self%tcache, tmparr, lenarr) + 4
    else
      iaend = GET_RANK_UNIFORM(seconds_between + self%tcache, tmparr, lenarr) + 4
    end if
    if (iaend > lenarr) iaend = lenarr; 
    self%tbeg(1) = tmparr(iabeg); 
    self%tend(1) = tmparr(iaend); 
    self%len_times = (iaend - iabeg) + 1; 
    call self%reallocate_arrays()

    self%times = tmparr(iabeg:iaend); 
    call check_err(nf90_get_var(ncid, lunar_distance_varid, self%lunar_distances, (/iabeg/), (/self%len_times/)))
    call check_err(nf90_get_var(ncid, solar_distance_varid, self%solar_distances, (/iabeg/), (/self%len_times/)))
    call check_err(nf90_get_var(ncid, lunar_ra_varid, self%lunar_ras, (/iabeg/), (/self%len_times/)))
    call check_err(nf90_get_var(ncid, solar_ra_varid, self%solar_ras, (/iabeg/), (/self%len_times/)))
    call check_err(nf90_get_var(ncid, lunar_dec_varid, self%lunar_decs, (/iabeg/), (/self%len_times/)))
    call check_err(nf90_get_var(ncid, solar_dec_varid, self%solar_decs, (/iabeg/), (/self%len_times/)))

    ! Close the NetCDF file
    call check_err(nf90_close(ncid))
    deallocate (tmparr)
    self%first = .false.

  end function recache_data
#endif

  subroutine reallocate_arrays(self)
    implicit none
    class(t_ephemerides), intent(inout) :: self

    ! Allocate arrays !
    if (allocated(self%times)) then
      deallocate (self%times)
      deallocate (self%lunar_distances)
      deallocate (self%solar_distances)
      deallocate (self%lunar_ras, self%solar_ras)
      deallocate (self%lunar_decs, self%solar_decs)
    end if
    allocate (self%times(self%len_times))
    allocate (self%lunar_distances(self%len_times))
    allocate (self%solar_distances(self%len_times))
    allocate (self%lunar_ras(self%len_times))
    allocate (self%solar_ras(self%len_times))
    allocate (self%lunar_decs(self%len_times))
    allocate (self%solar_decs(self%len_times))
  end subroutine reallocate_arrays

  ! check whether RA change cross the zero hr line.
  ! if so, adjust accordingly
  subroutine adjust_RA(RA, fadjust)
    implicit none

    logical, optional :: fadjust
    real(8), intent(INOUT) :: RA(:)

    integer :: I, N
    logical :: CrossZeroRA

    N = ubound(RA, 1); 
    CrossZeroRA = .false.; 
    do I = 1, N - 1
      if (abs(RA(I + 1) - RA(I)) > 120.d0) then
        CrossZeroRA = .true.; 
      end if
    end do

    if (CrossZeroRA) then
      where (RA < 180.d0) RA = RA + 360.d0; 
    end if

    if (present(fadjust)) fadjust = CrossZeroRA; 
    return; 
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
      yi = yi + y(k)*term
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
        term = term*(xi - x(i))/(x(k) - x(i))
      end if
    end do
  end subroutine L

  ! Return a rank j in ARR such that arr(j-1) < val <= arr(j1) !
!  function GET_RANK_SIMPLE_SEARCH(val, arr, len) result(rank)
!    implicit none
!
!    integer :: rank
!    real(8), intent(IN) :: val
!    real(8), intent(IN) :: arr(:)
!    integer, intent(IN) :: len
!
!    integer :: J
!
!    if (val < arr(1)) then
!      rank = 0;
!      return;
!    end if
!
!    if (val > arr(len)) then
!      rank = len + 1;
!      return;
!    end if
!
!    RANK = 0;
!    do J = 1, LEN
!      if (ARR(J) >= val) then
!        RANK = J;
!        exit;
!      end if
!    end do
!
!  end function GET_RANK_SIMPLE_SEARCH

  ! Find a rank j in arr such that arr(j-1) < val < = arr(j)         !
  ! in a unifrom table, i,e, arr(i+1) - arr(i) = arr(i+2) - arr(i+1) !
  function GET_RANK_UNIFORM(val, arr, len) result(rank)
    implicit none

    integer :: rank
    real(8), intent(IN) :: val
    real(8), intent(IN) :: arr(:)
    integer, intent(IN) :: len

    real(8) :: ds

    if (val < arr(1)) then
      rank = 0; 
      return; 
    end if

    if (val > arr(len)) then
      rank = len + 1; 
      return; 
    end if

    ds = (arr(2) - arr(1)); 
    rank = ceiling((val - arr(1))/ds) + 1; 
  end function GET_RANK_UNIFORM

  ! Find a rank j in arr such that arr(j-1) < val < = arr(j)  !
  ! using a binary search                                     !
  function GET_RANK_BINARY_SEARCH(val, arr, len) result(rank)
    implicit none

    integer :: rank
    real(8), intent(IN) :: val
    real(8), intent(IN) :: arr(:)
    integer :: len

    integer :: low, high, mid

    if (val < arr(1)) then
      rank = 0; 
      return; 
    end if
    if (val > arr(len)) then
      rank = len + 1; 
      return; 
    end if

    low = 1; 
    high = len; 
    do
      mid = (low + high)/2; 
      if (abs(val - arr(mid)) < 1.0e-12) then
        ! Exact match is found !
        rank = mid; 
        exit; 
      end if

      if (val < arr(mid)) then
        high = mid; 
      else
        low = mid; 
      end if

      if (high == low + 1) then
        ! found the interval to which val  !
        ! falls inbetween is found         !
        rank = high; 
        exit; 
      end if

      if (high < low) then
        write (*, *) "Error in GET_RANK_BINARY_SEARCH()"; 
        rank = -1; 
        exit; 
      end if
    end do

  end function GET_RANK_BINARY_SEARCH

end module mod_ephemerides
