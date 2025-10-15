!-------------------------------------------------------------------------------!
!
! ADCIRC - The ADvanced CIRCulation model
! Copyright (C) 1994-2025 R.A. Luettich, Jr., J.J. Westerink
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
module mod_sun_coordinates

   use mod_astronomic, only: t_astronomic_values

   implicit none

   interface SUN_COORDINATES
      module procedure SUN_COORDINATES_SUB0, SUN_COORDINATES_SUB1
   end interface SUN_COORDINATES

   private

   public :: SUN_COORDINATES

contains

   ! Sun's equation of the center. Page 164
   !  T - Julain centuries of 36525 ephemeris from epoch J2000
   !  M - The mean anomaly of the sun
   real(8) elemental function SUN_CENTER(T, M) result(C)
      use ADC_CONSTANTS, only: DEG2RAD
      implicit none
      real(8), intent(IN) :: T, M
      real(8) :: M_R

      M_R = M*DEG2RAD

      C = (1.914602d0 + T*(-0.004817d0 - 0.000014d0*T))*sin(M_R) + &
          (0.019993d0 - 0.000101d0*T)*sin(2.d0*M_R) + &
          0.000289d0*sin(3.d0*M_R)

   end function SUN_CENTER

   ! The Sun's true longitude. Page 164.
   !  L0 = the geometric mean longitude of the Sun
   !  C  = the Sun's eq1.495978707×1011uation of the center
   !  Dpsi (optional) = nutation
   real(8) elemental function SUN_LON(C, L0, DPsi) result(SLON)
      implicit none
      real(8), intent(IN) :: L0, C
      real(8), optional, intent(IN) :: Dpsi

      SLON = L0 + C ! Sun's true longtiude

      if (present(DPsi)) then
         SLON = SLON - 0.00569d0 + Dpsi
      end if
   end function SUN_LON

   ! True anomaly of the Sun longitude. Page 164.
   !  M  = the mean anomaly of the Sun
   !  C  = the Sun's equation of the center
   real(8) elemental function SUN_NU(M, C) result(nu)
      implicit none
      real(8), intent(IN) :: M, C
      nu = M + C ! its true anomerly
   end function SUN_NU

   ! The Sun's radius. Disance bewteen the centers of the Sub and
   ! the Earth in astronomical units
   !   nu = True anomaly of the Sun's true longitude
   !   eccen = eccentricity of the Earth's orbit
   real(8) elemental function SUN_RADIUS(nu, eccen) result(R)
      use ADC_CONSTANTS, only: DEG2RAD
      implicit none
      real(8), intent(IN) :: nu, eccen
      real(8) :: NUM, DEN, nu_rad

      nu_rad = nu*DEG2RAD

      DEN = 1.d0 + eccen*cos(nu_rad)
      NUM = 1.000001018d0*(1.d0 - eccen*eccen)

      R = NUM/DEN

   end function SUN_RADIUS

   ! Solar's declination. Page 165.
   !  - SLON = the Sun's longtitude
   !  - VarEps = the obliquity of the ecliptic
!  real(8) elemental function SOLAR_DEC(SLON, VarEps) result(DEC)
!    use ADC_CONSTANTS, only: DEG2RAD, RAD2DEG
!    implicit none
!
!    real(8), intent(IN) :: SLON, VarEps
!
!    real(8) :: SLON_RAD, VarEps_RAD
!
!    SLON_RAD = SLON*DEG2RAD
!    VarEps_RAD = VarEps*DEG2RAD
!
!    DEC = sin(SLON_RAD)*sin(VarEps_RAD)
!
!    DEC = asin(DEC)*RAD2DEG
!  end function SOLAR_DEC

   ! Solar's right ascension. Page 165.
   !  - SLON = the Sun's longtitude
   !  - VarEps = the obliquity of the ecliptic
!  real(8) elemental function SOLAR_RA(SLON, VarEps) result(RA)
!    use ADC_CONSTANTS, only: DEG2RAD, RAD2DEG
!    implicit none
!
!    real(8), intent(IN) :: SLON, VarEps
!
!    real(8) :: NUM, DEN
!    real(8) :: SLON_RAD, VarEps_RAD
!
!    SLON_RAD = SLON*DEG2RAD
!    VarEps_RAD = VarEps*DEG2RAD
!
!    NUM = cos(VarEps_RAD)*sin(SLON_RAD)
!    DEN = cos(SLON_RAD)
!
!    RA = mod(atan2(NUM, DEN)*RAD2DEG, 360.d0)
!  end function SOLAR_RA

   !
   ! Chapter 25. Solar coordinates. Page. 163-169
   ! Sun's coordinates.
   ! INPUT:
   !   RA  - Sun's right ascension (Degrees)
   !   DEC - Sun's declination (degrees)
   !   Delta  - Distance between the Earth and the Sun
   ! OUTPUT:
   !   JD - Julain days
   !   ASVAL  (optional) - derived data type storing precomputed
   !                       values for the moon orbit at the time JD
   subroutine SUN_COORDINATES_SUB1(RA, DEC, Delta, JD, ASVAL, NUTATION)
      use mod_astronomic, only: JULIAN_CENTURIES, M_DEG, LP_DEG, OMEGA_DEG, &
                                eccentricity_earth_orbit, L0_DEG, DeltaPsiL, &
                                DeltaVarepsL, varepsilon0_ecliptic, ECLIP2EQ
      implicit none

      real(8), intent(OUT) :: RA, DEC, Delta

      real(8), intent(IN) :: JD
      type(t_astronomic_values), optional, intent(IN) :: ASVAL
      logical, intent(in), optional :: NUTATION

      real(8) :: lambda, beta, T, LP, OG, L0, M, C, eccen
      real(8) :: DPsi, vareps0, DVareps, vareps
      logical :: have_asval
      logical :: use_nutation

      have_asval = .false.
      use_nutation = .true.

      if (present(ASVAL)) then
         if (abs(ASVAL%JD - JD) < 1.0d-9) then
            have_asval = .true.
         end if
      end if

      if (present(NUTATION)) then
         use_nutation = NUTATION
      end if

      if (.not. have_asval) then
         T = JULIAN_CENTURIES(JD)
         M = modulo(M_DEG(T), 360.d0)
         LP = modulo(LP_DEG(T), 360.d0)
         eccen = eccentricity_earth_orbit(T)

         OG = modulo(OMEGA_DEG(T), 360.d0)
         L0 = modulo(L0_DEG(T), 360.d0)
         DPsi = DeltaPsiL(OG, L0, LP)
         DVareps = DeltaVarepsL(OG, L0, LP)

         vareps0 = varepsilon0_ecliptic(T)
         vareps = vareps0 + Dvareps
      else
         T = ASVAL%T
         M = modulo(ASVAL%M, 360.d0)
         L0 = modulo(ASVAL%L0, 360.d0)
         eccen = ASVAL%eccen
         DPsi = ASVAL%DPsi
         vareps = ASVAL%vareps
      end if
      C = SUN_CENTER(T, M)

      if (use_nutation) then
         call SUN_COORDINATES_SUB0(lambda, beta, Delta, L0, M, eccen, C, DPsi)
      else
         call SUN_COORDINATES_SUB0(lambda, beta, Delta, L0, M, eccen, C)
      end if

      ! Apperent right ascension & declination
      call ECLIP2EQ(RA, DEC, lambda, beta, vareps)

   end subroutine SUN_COORDINATES_SUB1

   !
   ! Chapter 25. Solar coordinates. Page. 163-169
   ! Sun's coordinates.
   ! INPUT:
   !   lambda - Sun's longtitude in ecliptic
   !   beta   - Sun's latitude (degrees)
   !   Delta  - Distance between the Earth and the Sun
   ! OUTPUT:
   !   L0 - Mean longitude of the Sun
   !   M  - Mean anomaly of the Sun
   !   C  -  The Sun's equation of the center
   !   eccen - The eccentricity of the Earth's orbit
   !   DPsi (optional) -  The nutation in longitude.
   !                      If present, lambda is the
   !                      aparent lon.
   subroutine SUN_COORDINATES_SUB0(lambda, beta, Delta, L0, M, eccen, C, DPsi)
      implicit none

      real(8), intent(OUT) :: lambda, beta, Delta
      real(8), intent(IN) :: L0, M, C, eccen
      real(8), optional, intent(IN) :: DPsi

      real(8) :: snu

      beta = 0.d0
      if (.not. present(DPsi)) then
         lambda = SUN_LON(C, L0)
      else
         lambda = SUN_LON(C, L0, DPsi)
      end if

      snu = SUN_NU(M, C)
      Delta = SUN_RADIUS(snu, eccen) ! Distance between centers
      ! of the Earth and Sun
      return
   end subroutine SUN_COORDINATES_SUB0

end module mod_sun_coordinates
