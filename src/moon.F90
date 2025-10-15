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
!
! Postion of the Moon. Chapter 47
! Page. 337 -344
module mod_moon_position

   implicit none

   ! Variables string periodic terms in the the Moon's longtitude,
   ! latitude, and dfistance
   integer, parameter :: NPER = 60 !

   !  TABLE 47. B. page 341. Periodic terms for the latitude of
   !  the Moon (\Sum b). The unit is 1e-6 degree
   real(8), parameter :: LATARG(NPER, 4) = reshape([0.d0, 0.d0, 0.d0, 1.d0, &
                                                    0.d0, 0.d0, 1.d0, 1.d0, &
                                                    0.d0, 0.d0, 1.d0, -1.d0, &
                                                    2.d0, 0.d0, 0.d0, -1.d0, &
                                                    2.d0, 0.d0, -1.d0, 1.d0, &
                                                    2.d0, 0.d0, -1.d0, -1.d0, &
                                                    2.d0, 0.d0, 0.d0, 1.d0, &
                                                    0.d0, 0.d0, 2.d0, 1.d0, &
                                                    2.d0, 0.d0, 1.d0, -1.d0, &
                                                    0.d0, 0.d0, 2.d0, -1.d0, &
                                                    2.d0, -1.d0, 0.d0, -1.d0, &
                                                    2.d0, 0.d0, -2.d0, -1.d0, &
                                                    2.d0, 0.d0, 1.d0, 1.d0, &
                                                    2.d0, 1.d0, 0.d0, -1.d0, &
                                                    2.d0, -1.d0, -1.d0, 1.d0, &
                                                    2.d0, -1.d0, 0.d0, 1.d0, &
                                                    2.d0, -1.d0, -1.d0, -1.d0, &
                                                    0.d0, 1.d0, -1.d0, -1.d0, &
                                                    4.d0, 0.d0, -1.d0, -1.d0, &
                                                    0.d0, 1.d0, 0.d0, 1.d0, &
                                                    0.d0, 0.d0, 0.d0, 3.d0, &
                                                    0.d0, 1.d0, -1.d0, 1.d0, &
                                                    1.d0, 0.d0, 0.d0, 1.d0, &
                                                    0.d0, 1.d0, 1.d0, 1.d0, &
                                                    0.d0, 1.d0, 1.d0, -1.d0, &
                                                    0.d0, 1.d0, 0.d0, -1.d0, &
                                                    1.d0, 0.d0, 0.d0, -1.d0, &
                                                    0.d0, 0.d0, 3.d0, 1.d0, &
                                                    4.d0, 0.d0, 0.d0, -1.d0, &
                                                    4.d0, 0.d0, -1.d0, 1.d0, &
                                                    0.d0, 0.d0, 1.d0, -3.d0, &
                                                    4.d0, 0.d0, -2.d0, 1.d0, &
                                                    2.d0, 0.d0, 0.d0, -3.d0, &
                                                    2.d0, 0.d0, 2.d0, -1.d0, &
                                                    2.d0, -1.d0, 1.d0, -1.d0, &
                                                    2.d0, 0.d0, -2.d0, 1.d0, &
                                                    0.d0, 0.d0, 3.d0, -1.d0, &
                                                    2.d0, 0.d0, 2.d0, 1.d0, &
                                                    2.d0, 0.d0, -3.d0, -1.d0, &
                                                    2.d0, 1.d0, -1.d0, 1.d0, &
                                                    2.d0, 1.d0, 0.d0, 1.d0, &
                                                    4.d0, 0.d0, 0.d0, 1.d0, &
                                                    2.d0, -1.d0, 1.d0, 1.d0, &
                                                    2.d0, -2.d0, 0.d0, -1.d0, &
                                                    0.d0, 0.d0, 1.d0, 3.d0, &
                                                    2.d0, 1.d0, 1.d0, -1.d0, &
                                                    1.d0, 1.d0, 0.d0, -1.d0, &
                                                    1.d0, 1.d0, 0.d0, 1.d0, &
                                                    0.d0, 1.d0, -2.d0, -1.d0, &
                                                    2.d0, 1.d0, -1.d0, -1.d0, &
                                                    1.d0, 0.d0, 1.d0, 1.d0, &
                                                    2.d0, -1.d0, -2.d0, -1.d0, &
                                                    0.d0, 1.d0, 2.d0, 1.d0, &
                                                    4.d0, 0.d0, -2.d0, -1.d0, &
                                                    4.d0, -1.d0, -1.d0, -1.d0, &
                                                    1.d0, 0.d0, 1.d0, -1.d0, &
                                                    4.d0, 0.d0, 1.d0, -1.d0, &
                                                    1.d0, 0.d0, -1.d0, -1.d0, &
                                                    4.d0, -1.d0, 0.d0, -1.d0, &
                                                    2.d0, -2.d0, 0.d0, 1.d0], shape(LATARG), ORDER=[2, 1])

   real(8), parameter :: LAT_SIN_COEF(NPER) = [5128122.d0, &
                                               280602.d0, &
                                               277693.d0, &
                                               173237.d0, &
                                               55413.d0, &
                                               46271.d0, &
                                               32573.d0, &
                                               17198.d0, &
                                               9266.d0, &
                                               8822.d0, &
                                               8216.d0, &
                                               4324.d0, &
                                               4200.d0, &
                                               -3359.d0, &
                                               2463.d0, &
                                               2211.d0, &
                                               2065.d0, &
                                               -1870.d0, &
                                               1828.d0, &
                                               -1794.d0, &
                                               -1749.d0, &
                                               -1565.d0, &
                                               -1491.d0, &
                                               -1475.d0, &
                                               -1410.d0, &
                                               -1344.d0, &
                                               -1335.d0, &
                                               1107.d0, &
                                               1021.d0, &
                                               833.d0, &
                                               777.d0, &
                                               671.d0, &
                                               607.d0, &
                                               596.d0, &
                                               491.d0, &
                                               -451.d0, &
                                               439.d0, &
                                               422.d0, &
                                               421.d0, &
                                               -366.d0, &
                                               -351.d0, &
                                               331.d0, &
                                               315.d0, &
                                               302.d0, &
                                               -283.d0, &
                                               -229.d0, &
                                               223.d0, &
                                               223.d0, &
                                               -220.d0, &
                                               -220.d0, &
                                               -185.d0, &
                                               181.d0, &
                                               -177.d0, &
                                               176.d0, &
                                               166.d0, &
                                               -164.d0, &
                                               132.d0, &
                                               -119.d0, &
                                               115.d0, &
                                               107.d0]

   ! TABLE 47. A.
   !  Periodic terms for the longitude \Sum l and distance \Sum r
   !  of the Moon. The unit is 1e-6 for \Sum l and Kilometer for \Sum r
   real(8), parameter :: LONARG(NPER, 4) = reshape([0.d0, 0.d0, 1.d0, 0.d0, &
                                                    2.d0, 0.d0, -1.d0, 0.d0, &
                                                    2.d0, 0.d0, 0.d0, 0.d0, &
                                                    0.d0, 0.d0, 2.d0, 0.d0, &
                                                    0.d0, 1.d0, 0.d0, 0.d0, &
                                                    0.d0, 0.d0, 0.d0, 2.d0, &
                                                    2.d0, 0.d0, -2.d0, 0.d0, &
                                                    2.d0, -1.d0, -1.d0, 0.d0, &
                                                    2.d0, 0.d0, 1.d0, 0.d0, &
                                                    2.d0, -1.d0, 0.d0, 0.d0, &
                                                    0.d0, 1.d0, -1.d0, 0.d0, &
                                                    1.d0, 0.d0, 0.d0, 0.d0, &
                                                    0.d0, 1.d0, 1.d0, 0.d0, &
                                                    2.d0, 0.d0, 0.d0, -2.d0, &
                                                    0.d0, 0.d0, 1.d0, 2.d0, &
                                                    0.d0, 0.d0, 1.d0, -2.d0, &
                                                    4.d0, 0.d0, -1.d0, 0.d0, &
                                                    0.d0, 0.d0, 3.d0, 0.d0, &
                                                    4.d0, 0.d0, -2.d0, 0.d0, &
                                                    2.d0, 1.d0, -1.d0, 0.d0, &
                                                    2.d0, 1.d0, 0.d0, 0.d0, &
                                                    1.d0, 0.d0, -1.d0, 0.d0, &
                                                    1.d0, 1.d0, 0.d0, 0.d0, &
                                                    2.d0, -1.d0, 1.d0, 0.d0, &
                                                    2.d0, 0.d0, 2.d0, 0.d0, &
                                                    4.d0, 0.d0, 0.d0, 0.d0, &
                                                    2.d0, 0.d0, -3.d0, 0.d0, &
                                                    0.d0, 1.d0, -2.d0, 0.d0, &
                                                    2.d0, 0.d0, -1.d0, 2.d0, &
                                                    2.d0, -1.d0, -2.d0, 0.d0, &
                                                    1.d0, 0.d0, 1.d0, 0.d0, &
                                                    2.d0, -2.d0, 0.d0, 0.d0, &
                                                    0.d0, 1.d0, 2.d0, 0.d0, &
                                                    0.d0, 2.d0, 0.d0, 0.d0, &
                                                    2.d0, -2.d0, -1.d0, 0.d0, &
                                                    2.d0, 0.d0, 1.d0, -2.d0, &
                                                    2.d0, 0.d0, 0.d0, 2.d0, &
                                                    4.d0, -1.d0, -1.d0, 0.d0, &
                                                    0.d0, 0.d0, 2.d0, 2.d0, &
                                                    3.d0, 0.d0, -1.d0, 0.d0, &
                                                    2.d0, 1.d0, 1.d0, 0.d0, &
                                                    4.d0, -1.d0, -2.d0, 0.d0, &
                                                    0.d0, 2.d0, -1.d0, 0.d0, &
                                                    2.d0, 2.d0, -1.d0, 0.d0, &
                                                    2.d0, 1.d0, -2.d0, 0.d0, &
                                                    2.d0, -1.d0, 0.d0, -2.d0, &
                                                    4.d0, 0.d0, 1.d0, 0.d0, &
                                                    0.d0, 0.d0, 4.d0, 0.d0, &
                                                    4.d0, -1.d0, 0.d0, 0.d0, &
                                                    1.d0, 0.d0, -2.d0, 0.d0, &
                                                    2.d0, 1.d0, 0.d0, -2.d0, &
                                                    0.d0, 0.d0, 2.d0, -2.d0, &
                                                    1.d0, 1.d0, 1.d0, 0.d0, &
                                                    3.d0, 0.d0, -2.d0, 0.d0, &
                                                    4.d0, 0.d0, -3.d0, 0.d0, &
                                                    2.d0, -1.d0, 2.d0, 0.d0, &
                                                    0.d0, 2.d0, 1.d0, 0.d0, &
                                                    1.d0, 1.d0, -1.d0, 0.d0, &
                                                    2.d0, 0.d0, 3.d0, 0.d0, &
                                                    2.d0, 0.d0, -1.d0, -2.d0], shape(LONARG), order=[2, 1])

   real(8), parameter :: LON_SIN_COEF(NPER) = [6288774.d0, &
                                               1274027.d0, &
                                               658314.d0, &
                                               213618.d0, &
                                               -185116.d0, &
                                               -114332.d0, &
                                               58793.d0, &
                                               57066.d0, &
                                               53322.d0, &
                                               45758.d0, &
                                               -40923.d0, &
                                               -34720.d0, &
                                               -30383.d0, &
                                               15327.d0, &
                                               -12528.d0, &
                                               10980.d0, &
                                               10675.d0, &
                                               10034.d0, &
                                               8548.d0, &
                                               -7888.d0, &
                                               -6766.d0, &
                                               -5163.d0, &
                                               4987.d0, &
                                               4036.d0, &
                                               3994.d0, &
                                               3861.d0, &
                                               3665.d0, &
                                               -2689.d0, &
                                               -2602.d0, &
                                               2390.d0, &
                                               -2348.d0, &
                                               2236.d0, &
                                               -2120.d0, &
                                               -2069.d0, &
                                               2048.d0, &
                                               -1773.d0, &
                                               -1595.d0, &
                                               1215.d0, &
                                               -1110.d0, &
                                               -892.d0, &
                                               -810.d0, &
                                               759.d0, &
                                               -713.d0, &
                                               -700.d0, &
                                               691.d0, &
                                               596.d0, &
                                               549.d0, &
                                               537.d0, &
                                               520.d0, &
                                               -487.d0, &
                                               -399.d0, &
                                               -381.d0, &
                                               351.d0, &
                                               -340.d0, &
                                               330.d0, &
                                               327.d0, &
                                               -323.d0, &
                                               299.d0, &
                                               294.d0, &
                                               0.d0]

   real(8), parameter :: LON_COS_COEF(NPER) = [-20905355.d0, &
                                               -3699111.d0, &
                                               -2955968.d0, &
                                               -569925.d0, &
                                               48888.d0, &
                                               -3149.d0, &
                                               246158.d0, &
                                               -152138.d0, &
                                               -170733.d0, &
                                               -204586.d0, &
                                               -129620.d0, &
                                               108743.d0, &
                                               104755.d0, &
                                               10321.d0, &
                                               0.d0, &
                                               79661.d0, &
                                               -34782.d0, &
                                               -23210.d0, &
                                               -21636.d0, &
                                               24208.d0, &
                                               30824.d0, &
                                               -8379.d0, &
                                               -16675.d0, &
                                               -12831.d0, &
                                               -10445.d0, &
                                               -11650.d0, &
                                               14403.d0, &
                                               -7003.d0, &
                                               0.d0, &
                                               10056.d0, &
                                               6322.d0, &
                                               -9884.d0, &
                                               5751.d0, &
                                               0.d0, &
                                               -4950.d0, &
                                               4130.d0, &
                                               0.d0, &
                                               -3958.d0, &
                                               0.d0, &
                                               3258.d0, &
                                               2616.d0, &
                                               -1897.d0, &
                                               -2117.d0, &
                                               2354.d0, &
                                               0.d0, &
                                               0.d0, &
                                               -1423.d0, &
                                               -1117.d0, &
                                               -1571.d0, &
                                               -1739.d0, &
                                               0.d0, &
                                               -4421.d0, &
                                               0.d0, &
                                               0.d0, &
                                               0.d0, &
                                               0.d0, &
                                               1165.d0, &
                                               0.d0, &
                                               0.d0, &
                                               8752.d0]

!  REAL(8), private, parameter:: LATARG(NPER,4) = transpose( LATARGTMP ) ;

! This addresses a bug in the nvidia fortran compiler
! See: https://forums.developer.nvidia.com/t/nvfortran-constant-array-initialization-bug/323012
! Instead of of this being a compile time constant, we generate it later
#ifndef __NVCOMPILER
   integer, parameter :: LATMEU(NPER) = int(abs(LATARG(:, 2)))
   integer, parameter :: LONMEU(NPER) = int(abs(LONARG(:, 2)))
#endif

   interface MOON_COORDINATES
      module procedure MOON_COORDINATES_SUB0, MOON_COORDINATES_SUB1
   end interface MOON_COORDINATES

   private

   public :: MOON_COORDINATES

contains

   ! Coefficient muitplying sine and cosine arguments. page 338
   real(8) elemental function E_MUL_COEF(T) result(E)
      implicit none
      real(8), intent(IN) :: T
      E = 1.d0 + T*(-0.002516d0 + T*(-0.0000074d0))
   end function E_MUL_COEF

   ! Page. 338
   real(8) elemental function A1_DEG(T) result(A1)
      implicit none
      real(8), intent(IN) :: T
      A1 = 119.75d0 + 131.849d0*T
   end function A1_DEG

   ! Page. 338
   real(8) elemental function A2_DEG(T) result(A2)
      implicit none
      real(8), intent(IN) :: T
      A2 = 53.09d0 + 479264.290d0*T
   end function A2_DEG

   ! Page. 338
   real(8) elemental function A3_DEG(T) result(A3)
      implicit none
      real(8), intent(IN) :: T
      A3 = 313.45d0 + 481266.484d0*T
   end function A3_DEG

   ! Coefficients for terms contains angle M. See description on
   ! Page 338
   function CAL_EPMUL(MEU, E) result(EPMUL)
      implicit none

      real(8) :: EPMUL(NPER)
      real(8), intent(IN) :: E
      integer, intent(IN) :: MEU(:)

      integer :: I, J

      epmul = 1.d0
      do I = 1, NPER
         do J = 1, MEU(I)
            epmul(I) = epmul(I)*E
         end do
      end do

   end function CAL_EPMUL

   ! Chaper 47.
   ! Output:
   !     RA = Right ascendsion (in Degrees)
   !     DEC = Declination (in Degrees)
   !     Delta = Distance from the Earth to the Moon (in kilometers)
   ! Input:
   !     JD = Julain days
   !     ASVAL (optional): derived data type storing precomputed
   !                       values for the moon orbit at the time JD
   !
   !     ASVAL%LP = Moon's mean longitude
   !     ASVAL%D  = Moon's mean elongation
   !     ASVAL%M  = Sun's mean anomaly
   !     ASVAL%MP = Moon's mean anomaly
   !     ASVAL%F  = Moon's argument of latitude (mean distance of the Moon
   !          from it ascending node)
   !     ASVAL%E  = coefficeint for correcting terms containg M
   !     ASVAL%A1, A2, A3 = coefficeint for additve term accounting
   !                  action of Venus and Jupiter
   ! Intended as a driver subroutine.
   ! it call as a function MOON_COORDINATES_SUB0
   !
   subroutine MOON_COORDINATES_SUB1(RA, DEC, Delta, JD, ASVAL, NUTATION)
      use mod_astronomic, only: t_astronomic_values, JULIAN_CENTURIES, LP_DEG, D_DEG, M_DEG, &
                                MP_DEG, F_DEG, OMEGA_DEG, L0_DEG, DeltaPsiL, DeltaVarepsL, &
                                varepsilon0_ecliptic, ECLIP2EQ
      implicit none

      real(8), intent(OUT) :: RA, DEC, Delta
      real(8), intent(IN) :: JD
      type(t_astronomic_values), optional, intent(IN) :: ASVAL
      logical, intent(in), optional :: NUTATION

      ! local !
      real(8) :: lambda, beta
      real(8) :: T, LP, D, M, MP, F, A1, A2, A3, E
      real(8) :: L0, OG, DPsi, Dvareps, vareps, vareps0
      logical :: have_asval
      logical :: use_nutation
      real(8) :: nutation_mul

      have_asval = .false.
      use_nutation = .false.

      ! Check if the astroval
      if (present(ASVAL)) then
         if (abs(JD - ASVAL%JD) < 1.d-9) have_asval = .true.
      end if
      if (present(NUTATION)) use_nutation = NUTATION

      if (.not. have_asval) then
         T = JULIAN_CENTURIES(JD)
         LP = modulo(LP_DEG(T), 360.d0)
         D = modulo(D_DEG(T), 360.d0)
         M = modulo(M_DEG(T), 360.d0)
         MP = modulo(MP_DEG(T), 360.d0)
         F = modulo(F_DEG(T), 360.d0)

         OG = modulo(OMEGA_DEG(T), 360.d0)
         L0 = modulo(L0_DEG(T), 360.d0)
         DPsi = DeltaPsiL(OG, L0, LP)

         vareps0 = varepsilon0_ecliptic(T)
         DVareps = DeltaVarepsL(OG, L0, LP)
         vareps = vareps0 + Dvareps
      else
         T = ASVAL%T
         LP = modulo(ASVAL%LP, 360.d0)
         D = modulo(ASVAL%D, 360.d0)
         M = modulo(ASVAL%M, 360.d0)
         MP = modulo(ASVAL%MP, 360.d0)
         F = modulo(ASVAL%F, 360.d0)

         DPsi = ASVAL%DPsi
         vareps = ASVAL%vareps
      end if

      A1 = modulo(A1_DEG(T), 360.d0)
      A2 = modulo(A2_DEG(T), 360.d0)
      A3 = modulo(A3_DEG(T), 360.d0)
      E = E_MUL_COEF(T)

      call MOON_COORDINATES_SUB0(lambda, beta, Delta, LP, D, M, MP, F, E, A1, A2, A3)

      ! PRINT*, "-------- RA & DEC --------"
      ! Apperent right ascension & declination
      nutation_mul = 1.d0
      if (.not. use_nutation) nutation_mul = 0.d0

      call ECLIP2EQ(RA, DEC, lambda + nutation_mul*DPsi, beta, vareps)

      ! Geometric right ascension & declination
      ! CALL ECLIP2EQ( RA, DEC, lambda, beta, vareps )
      ! PRINT*, "Geometrical RA  = ", MODULO(RA,360.D0)
      ! PRINT*, "Geometrical DEC = ", DEC

   end subroutine MOON_COORDINATES_SUB1

   ! Chaper 47.
   ! Output:
   !     lambda = Longtitude (in Degrees)
   !     beta = Latitude (in Degrees)
   !     Delta = Distance from the Earth to the Moon (in kilometers)
   ! Input:
   !     LP = Moon's mean longitude
   !     D  = Moon's mean elongation
   !     M  = Sun's mean anomaly
   !     MP = Moon's mean anomaly
   !     F  = Moon's argument of latitude (mean distance of the Moon
   !          from it ascending node)
   !     E  = coefficeint for correcting terms containg M
   !     A1, A2, A3 = coefficeint for additve term accounting
   !                  action of Venus and Jupiter
   ! - Don't account for nutation
   ! - Intended as a lower level subroutine to be called by
   !   higher level functions and subroutine
   subroutine MOON_COORDINATES_SUB0(lambda, beta, Delta, LP, D, M, MP, F, E, A1, A2, A3)
      use ADC_CONSTANTS, only: DEG2RAD
      implicit none

      real(8), intent(OUT) :: lambda, beta, Delta
      real(8), intent(IN) :: LP, D, M, MP, F
      real(8), intent(IN) :: E, A1, A2, A3

      real(8) :: suml, sumr, sumb
      real(8) :: vec(4), argval(NPER), epmul(NPER)

! This addresses a bug in the nvidia fortran compiler
! See: https://forums.developer.nvidia.com/t/nvfortran-constant-array-initialization-bug/323012
! Instead of of this being a compile time constant generated above, we generate it here on the fly
#ifdef __NVCOMPILER
      integer :: LATMEU(NPER)
      integer :: LONMEU(NPER)
      LATMEU(1:NPER) = int(abs(LATARG(1:NPER, 2)))
      LONMEU(1:NPER) = int(abs(LONARG(1:NPER, 2)))
#endif

      ! Convert from degree to radian
      vec = DEG2RAD*[D, M, MP, F]

      ! 1. Compute mean longitude !
      argval = matmul(LONARG, vec)
      epmul = cal_epmul(LONMEU, E)

      ! sine argument
      !  \sum l
      argval = epmul*LON_SIN_COEF*sin(argval)
      suml = sum(argval) + additive_suml()

      ! 2. Compute distance       !
      argval = matmul(LONARG, vec)

      ! cosine argument
      !  \Sum r
      argval = epmul*LON_COS_COEF*cos(argval)
      sumr = sum(argval)

      ! 3. Compute mean latitude !
      argval = matmul(LATARG, vec)
      epmul = cal_epmul(LATMEU, E)

      ! sine argument
      !  \Sum b
      argval = epmul*LAT_SIN_COEF*sin(argval)
      sumb = sum(argval) + additive_sumb()

      ! output
      lambda = LP + suml/1.0d6 ! Longitude in Degree
      beta = sumb/1.d6 ! Lattitude in Degree
      Delta = 385000.56d0 + sumr/1.0d3 ! Distance in km

   contains

      real(8) function additive_suml() result(asuml)
         implicit none
         asuml = 3958.d0*sin(DEG2RAD*A1) &
                 + 1962.d0*sin(DEG2RAD*(LP - F)) &
                 + 318.d0*sin(DEG2RAD*A2)
      end function additive_suml

      real(8) function additive_sumb() result(asumb)
         implicit none
         asumb = -2235.d0*sin(DEG2RAD*LP) &
                 + 382.d0*sin(DEG2RAD*A3) &
                 + 175.d0*sin(DEG2RAD*(A1 - F)) &
                 + 175.d0*sin(DEG2RAD*(A1 + F)) &
                 + 127.d0*sin(DEG2RAD*(LP - MP)) &
                 - 115.d0*sin(DEG2RAD*(LP + MP))

      end function additive_sumb

   end subroutine MOON_COORDINATES_SUB0

end module mod_moon_position
