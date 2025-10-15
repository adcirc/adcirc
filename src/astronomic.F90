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
! Ref:
!   J. MEEUS, Astronomic algorithms, 2nd Edition, 1991
!
! Implemented by D. Wirasaet, 2024
!
!  NOTE:
!    T - JDE or JD centuries
!
module mod_astronomic

   implicit none

   real(8), parameter :: SEC2DEG = 2.777777777777778d-4, &
                         MIN2DEG = 0.01666666666666d0

   ! Astronomical unit in km
   real(8), parameter :: AUDIST = 1.495978707d+8
   real(8), parameter :: km2AU = 1.d0/AUDIST

   ! Ratio of Sun/Moon mass to Earth mass: Wikipedia !
   real(8), parameter :: MassRatioSunEarth = 332946.0487d0
   real(8), parameter :: MassRatioMoonEarth = 1.d0/81.3005678d0

   ! Keep the value in two saperate numbers: significant and exponent
   real(8), parameter :: EarthRadiuskm(2) = [6371.0088d0, 1.0d0] ! km
   real(8), parameter :: EarthRadiusm(2) = [6.3710088d0, 1.0d+6] ! m
   real(8), parameter :: EarthRadiusAU(2) = [4.258756338033895d0, 1.0d-5] ! in AU

   type t_astronomic_values
      real(8) :: JD ! Julain days.
      real(8) :: T ! Julian centuries from the Epoch J2000.0
      real(8) :: LP ! Moon's mean longtitude, J2000.0 Epoch.
      real(8) :: D ! Mean elongation of the Moon, J2000.0 Epoch.
      real(8) :: M ! Sun's mean anomaly, J2000.0 Epoch.
      real(8) :: MP ! Moon's mean anomaly, J2000 Epoch.
      real(8) :: F ! Moon's argument of latitude, J2000 Epoch.
      real(8) :: OG ! Longtitude of the ascending node of the moon's mean orbit on the ecliptic.
      real(8) :: L0 ! Mean longitude of the Sun, refered to the mean equnoix of the date.
      real(8) :: DPsi ! The nutation in longitude.
      real(8) :: vareps0 ! Mean obliquity of the ecliptic.
      real(8) :: Dvareps ! The nutation in obliquity.
      real(8) :: vareps ! the obliquity of the ecliptic.
      real(8) :: eccen ! the eccentricity of the Earth's orbit

   contains
      procedure, public, pass(self) :: compute_astronomic_values
   end type t_astronomic_values

   interface ECLIP2EQ
      module procedure ECLIP2EQ_S, ECLIP2EQ_V
   end interface ECLIP2EQ

   private

   public ::  t_astronomic_values, JULIAN_CENTURIES, LP_DEG, D_DEG, M_DEG, &
             MP_DEG, F_DEG, OMEGA_DEG, L0_DEG, DeltaPsiL, DeltaVarepsL, &
             varepsilon0_ecliptic, ECLIP2EQ, km2AU, eccentricity_earth_orbit, &
             GMST_DEG, EarthRadiusM, MassRatioSunEarth, MassRatioMoonEarth, &
             EarthRadiusAU, EarthRadiusKM, JulianDay

contains

   subroutine compute_astronomic_values(self, JD)
      implicit none
      class(t_astronomic_values), intent(inout) :: self
      real(8), intent(in) :: JD
      real(8) :: T, OG, L0, LP

      self%JD = JD

      T = JULIAN_CENTURIES(JD)
      self%T = JULIAN_CENTURIES(JD)
      self%D = modulo(D_DEG(T), 360.d0)
      self%M = modulo(M_DEG(T), 360.d0)
      self%MP = modulo(MP_DEG(T), 360.d0)
      self%F = modulo(F_DEG(T), 360.d0)

      OG = modulo(OMEGA_DEG(T), 360.d0)
      L0 = modulo(L0_DEG(T), 360.d0)
      LP = modulo(LP_DEG(T), 360.d0)

      self%LP = LP
      self%L0 = L0
      self%OG = OG
      self%DPsi = DeltaPsiL(OG, L0, LP)
      self%vareps0 = varepsilon0_ecliptic(T)
      self%DVareps = DeltaVarepsL(OG, L0, LP)
      self%vareps = self%vareps0 + self%Dvareps
      self%eccen = eccentricity_earth_orbit(T)

   end subroutine compute_astronomic_values

   ! Julian day, p. 61
   real(8) function JULIANDAY(DD, MM, YYYY, CALENDAR_TYPE) result(JD)
      implicit none

      real(8), intent(IN) :: DD
      integer, intent(in) :: MM, YYYY
      character(LEN=*), intent(in), optional :: CALENDAR_TYPE

      integer :: A, B
      real(8) :: D, M, Y
      character(LEN=9) :: CTYPE
      D = dble(DD)
      M = dble(MM)
      Y = dble(YYYY)

      CTYPE = 'Gregorian'
      if (present(CALENDAR_TYPE)) then
         select case (trim(CALENDAR_TYPE))
         case ('Julian', 'JULIAN', 'julian')
            CTYPE = 'Julian'
         end select
      end if

      if (M <= 2) then
         Y = Y - 1
         M = M + 12
      end if

      A = floor(Y/100.d0)
      B = 2 - A + floor(dble(A)/4.d0)
      select case (trim(CTYPE))
      case ('Julian')
         B = 0
      end select

      JD = dble(floor(365.25d0*(Y + 4716d0))) + dble(floor(30.6001d0*(M + 1d0))) + dble(D) + dble(B) - 1524.50d0

   end function JULIANDAY

   ! - Compute Julian centuries from the Epoch J2000.0 (JDE
   !    2451545)
   function JULIAN_CENTURIES(JDE) result(T)
      implicit none

      real(8) :: T
      real(8), intent(IN) :: JDE

      T = (JDE - 2451545.d0)/36525.d0

   end function JULIAN_CENTURIES

   ! Sidereal time at Greenwich for any JD.                    !
   ! page. 87                                                  !
   !  JD     = Julian days                                     !
   !  DPSi (optional)   = The nutation in longitude (degrees)  !
   !  vareps (optional) = obliquity of the ecliptic.           !
   function GMST_DEG(JD, NUTATION, ASVAL) result(gmst)
      use ADC_CONSTANTS, only: DEG2RAD
      implicit none

      real(8) :: gmst
      real(8), intent(IN) :: JD
      logical, optional, intent(IN) :: NUTATION
      type(t_astronomic_values), intent(in), optional :: ASVAL

      real(8) :: T, JD0, RM
      logical :: include_nutation

      ! If nutation is require !
      real(8) :: LP, L0, OG, DPsi, vareps
      logical :: have_asval

      have_asval = .false.
      include_nutation = .false.

      ! Find JD of the date at UT 0h
      RM = JD - dble(floor(JD))

      JD0 = JD
      if (RM < 0.5d0 - 1.0d-10) then
         JD0 = dble(floor(JD)) - 0.5d0
      else if (RM > 0.5d0 + 1.0d-10) then
         JD0 = dble(floor(JD)) + 0.5d0
      end if
      RM = JD - JD0

      ! Compute T at UT 0h
      T = (JD0 - 2451545.0d0)/36525.d0

      ! \Theta0 EQ. (12.2) page 87
      gmst = 100.46061837d0 + 36000.770053608d0*T &
             + T*T*(0.000387933d0 - T/38710000.d0)

      ! Mean sidereal time at Greenwich for JD. Page 87
      gmst = gmst + 1.00273790935d0*RM*360.d0

      if (present(NUTATION)) include_nutation = NUTATION

      if (include_nutation) then
         if (present(ASVAL)) then
            if (abs(JD - ASVAL%JD) < 1.0d-9) have_asval = .true.
         end if

         if (.not. have_asval) then
            T = JULIAN_CENTURIES(JD)
            LP = modulo(LP_DEG(T), 360.d0)
            L0 = modulo(L0_DEG(T), 360.d0)

            OG = modulo(OMEGA_DEG(T), 360.d0)
            DPsi = DeltaPsiL(OG, L0, LP)
            vareps = varepsilon0_ecliptic(T) + DeltaVarepsL(OG, L0, LP)
         else
            DPsi = ASVAL%DPsi
            vareps = ASVAL%vareps
         end if

         gmst = gmst + DPsi*cos(vareps*DEG2RAD)/15.d0
      end if

!! Eq. 12.4. Yeild slightly differnt results from
!!           the formula above
!!          gmst = 280.46061837D0 + 360.98564736629*(JD -
!!     &       2451545.D0) + T*T*(0.000387933 - T/38710000.D0)
!!
      gmst = modulo(gmst, 360.d0)

   end function GMST_DEG

   ! Moon's mean longtitude, J2000.0 Epoch. page 338
   elemental function LP_DEG(T) result(LP)
      implicit none

      real(8) :: LP
      real(8), intent(IN) :: T ! Julian centuries

      LP = 218.3164477d0 + T*(481267.88123421d0 &
                              + T*(-0.0015786d0 &
                                   + T*(1.d0/538841.d0 &
                                        + T*(-1.d0/65194000.d0))))

   end function LP_DEG

   ! Mean elongation of the Moon, J2000.0 Epoch. page 338
   elemental function D_DEG(T) result(D)
      implicit none

      real(8) :: D
      real(8), intent(IN) :: T

      D = 297.8501921d0 + T*(445267.1114034d0 &
                             + T*(-0.0018819d0 &
                                  + T*(1.d0/545868.d0 &
                                       + T*(-1.d0/113065000.d0))))

   end function D_DEG

   ! Sun's mean anomaly, J2000.0 Epoch, page 338
   elemental function M_DEG(T) result(M)
      implicit none

      real(8) :: M
      real(8), intent(IN) :: T

      M = 357.5291092d0 + T*(35999.0502909d0 &
                             + T*(-0.0001536d0 &
                                  + T*(1.d0/24490000.d0)))
   end function M_DEG

   ! Moon's mean anomaly, J2000 Epoch. Page 338
   elemental function MP_DEG(T) result(MP)
      implicit none

      real(8) :: MP
      real(8), intent(IN) :: T

      MP = 134.9633964d0 + T*(477198.8675055d0 &
                              + T*(0.0087414d0 &
                                   + T*(1.d0/69699.d0 &
                                        + T*(-1.d0/14712000.d0))))
   end function MP_DEG

   ! Moon's argument of latitude, J2000 Epoch. Page 338
   elemental function F_DEG(T) result(F)
      implicit none

      real(8) :: F
      real(8), intent(IN) :: T

      F = 93.2720950d0 + T*(483202.0175233d0 &
                            + T*(-0.0036539d0 &
                                 + T*(-1.d0/3526000d0 &
                                      + T*(1.d0/863310000.d0))))
   end function F_DEG

   ! Longtitude of the ascending node of the moon's mean orbit on
   ! the ecliptic. Page 144
   real(8) elemental function OMEGA_DEG(T) result(OMEGA)
      implicit none
      real(8), intent(IN) :: T
      OMEGA = 125.04452d0 + T*(-1934.136261d0 &
                               + T*(0.0020708d0 &
                                    + T*(1.d0/450000.d0)))
   end function OMEGA_DEG

   ! Mean longitude of the Sun, refered to the mean equnoix of the
   ! date. Page 163
   real(8) elemental function L0_DEG(T) result(L0)
      implicit none
      real(8), intent(IN) :: T
      L0 = 280.46646d0 + T*(36000.76983d0 + T*(0.0003032d0))
   end function L0_DEG

   ! The eccentricity of the Earth's orbit. Page 163
   real(8) elemental function eccentricity_earth_orbit(T) result(eps)
      implicit none
      real(8), intent(IN) :: T
      eps = 0.016708634d0 + T*(-0.000042037d0 &
                               - 0.0000001267d0*T)
   end function eccentricity_earth_orbit

   ! The nutation in longitude. p. 144
   !  -- low accuracy
   ! \Delta \Psi
   real(8) elemental function DeltaPsiL(OMEGA, L, LP) result(DPsi)
      use ADC_CONSTANTS, only: DEG2RAD
      implicit none
      real(8), intent(IN) :: OMEGA, L, LP
      real(8) :: OG_R, L_R, LP_R

      OG_R = OMEGA*DEG2RAD
      L_R = L*DEG2RAD
      LP_R = LP*DEG2RAD

      DPsi = -(17.20d0*sec2deg)*sin(OG_R) &
             - (1.32d0*sec2deg)*sin(2.d0*L_R) &
             - (0.23d0*sec2deg)*sin(2.d0*LP_R) &
             + (0.21d0*sec2deg)*sin(2.d0*OG_R)
   end function DeltaPsiL

   ! The nutation in obliquity. p. 144
   !  -- low accuracy
   ! \Delta \varepsilon
   real(8) elemental function DeltaVarepsL(OMEGA, L, LP) result(DVareps)
      use ADC_CONSTANTS, only: DEG2RAD
      implicit none
      real(8), intent(IN) :: OMEGA, L, LP
      real(8) :: OG_R, L_R, LP_R

      OG_R = OMEGA*DEG2RAD
      L_R = L*DEG2RAD
      LP_R = LP*DEG2RAD

      DVareps = (9.20d0*sec2deg)*cos(OG_R) &
                + (0.57d0*sec2deg)*cos(2.d0*L_R) &
                + (0.10d0*sec2deg)*cos(2.d0*LP_R) &
                - (0.09d0*sec2deg)*cos(2.d0*OG_R)

   end function DeltaVarepsL

   ! Mean obliquity of the ecliptic. Eq. 22.3. Page 147
   real(8) elemental function varepsilon0_ecliptic(T) result(varepsilon0)
      implicit none
      real(8), intent(IN) :: T
      real(8) :: U

      U = 0.01d0*T

      varepsilon0 = 23.d0 + 26.d0*min2deg + 21.448d0*sec2deg

      varepsilon0 = varepsilon0 + U*(-4680.93d0*sec2deg &
                                     + U*(-1.55d0 &
                                          + U*(1999.25d0 &
                                               + U*(-51.38d0 &
                                                    + U*(-249.67d0 &
                                                         + U*(-39.05d00 &
                                                              + U*(7.12d0 &
                                                                   + U*(27.87d0 &
                                                                        + U*(5.79d0 &
                                                                             + U*(2.45d0))))))))))

   end function varepsilon0_ecliptic

   ! Coordinate transformation                            !
   !   Transformation from ecliptical into EQ coordinates !
   subroutine ECLIP2EQ_S(RA, DEC, LAMBDA, BETA, varepsilon)
      use ADC_CONSTANTS, only: deg2rad, rad2deg
      implicit none

      real(8), intent(OUT) :: RA, DEC
      real(8), intent(IN) :: LAMBDA, BETA, varepsilon

      real(8) :: NUM, DEN

      ! Right ascension
      NUM = sin(LAMBDA*deg2rad)*cos(varepsilon*deg2rad) &
            - tan(BETA*deg2rad)*sin(varepsilon*deg2rad)
      DEN = cos(LAMBDA*deg2rad)

      RA = atan2(NUM, DEN)

      ! Declination
      NUM = sin(BETA*deg2rad)*cos(varepsilon*deg2rad) &
            + cos(BETA*deg2rad)*sin(varepsilon*deg2rad)*sin(LAMBDA*deg2rad)

      DEC = asin(NUM)

      ! Convert to degree
      RA = modulo(RA*RAD2DEG, 360.d0)
      DEC = DEC*RAD2DEG

      return
   end subroutine ECLIP2EQ_S

   ! Coordinate transformation                            !
   !   Transformation from ecliptical into EQ coordinates !
   subroutine ECLIP2EQ_V(RA, DEC, LAMBDA, BETA, varepsilon)
      use ADC_CONSTANTS, only: deg2rad, rad2deg
      implicit none

      real(8), dimension(:), intent(OUT) :: RA, DEC
      real(8), dimension(:), intent(IN) :: LAMBDA, BETA, varepsilon

      integer :: len
      real(8), dimension(:), allocatable :: NUM, DEN

      len = size(LAMBDA)
      allocate (NUM(len), DEN(len))

      ! Right ascension
      NUM = sin(LAMBDA*deg2rad)*cos(varepsilon*deg2rad) &
            - tan(BETA*deg2rad)*sin(varepsilon*deg2rad)
      DEN = cos(LAMBDA*deg2rad)

      RA = atan2(NUM, DEN)

      ! Declination
      NUM = sin(BETA*deg2rad)*cos(varepsilon*deg2rad) &
            + cos(BETA*deg2rad)*sin(varepsilon*deg2rad)*sin(LAMBDA*deg2rad)

      DEC = asin(NUM)

      ! Convert to degree
      RA = modulo(RA*RAD2DEG, 360.d0)
      DEC = DEC*RAD2DEG

      deallocate (NUM, DEN)
   end subroutine ECLIP2EQ_V

end module mod_astronomic
