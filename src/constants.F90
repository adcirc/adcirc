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
!-----------------------------------------------------------------------
!  MODULE ADC_CONSTANTS
!-----------------------------------------------------------------------
!> @author Zachary Cobell, The Water Institute, zcobell@thewaterinstitute.org
!>
!> @copyright Dr. R.A. Luettich and Dr. J.J. Westerink
!>
!> @brief This module contains physical constants and conversion factors used
!> throughout the ADCIRC code
!>
!> In general, all variables within this file should be constants (i.e. they
!> should be declared as a parameter) with exceptions for physical constants
!> that the user can change in the input files
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
module ADC_CONSTANTS

   implicit none

   !---Physical Constants---!

   !...Default gravitational acceleration
   ! Note that the user can update this in the
   ! control file so this isn't made constant
   real(8) :: g = 9.80665

   !...Nominal density of water
   real(8), parameter :: RhoWat0 = 1000.0d0

   !...Sigma T value of reference density
   real(8), parameter :: SigT0 = RHOWAT0 - 1000.0d0

   !...Background Atmospheric pressure in mb
   real(8), parameter :: PRBCKGRND = 1013.0d0
   ! REAL(8), PARAMETER :: PRBCKGRND = 1013.25D0 !...U.S. 1976 Std Atm

   !...Background temperature (kelvin)
   real(8), parameter :: TBCKGRND = 288.15d0

   !...Specific gas constant
   real(8), parameter :: RAir = 287.058d0

   !...Ideal gas law kg/m^3
   !...Note this is the default value but may be modified
   ! by the user at runtime
   real(8) :: RhoAir = 1.293d0
   !...There are other options which may be used but
   ! the default is above
   ! REAL(8), PARAMETER :: RhoAir = 1.15D0 !...Holland
   ! REAL(8), PARAMETER :: RhoAir = 1.1774D0 !... T=79.79, P=1013D0

   real(8), parameter :: PI = 3.141592653589793d0
   real(8), parameter :: TWOPI = PI*2.d0
   real(8), parameter :: HFPI = PI/2.d0
   real(8), parameter :: e = 2.718281828459045d0

   !...Radius of earth (m)
   real(8), parameter :: Rearth = 6378206.4d0

   real(8), parameter :: omega = 7.29212d-5

   !---Time based factors---!
   real(8), parameter :: hour2sec = 3600.0d0
   real(8), parameter :: sec2hour = 1.0d0/hour2sec
   real(8), parameter :: day2hour = 24.0d0
   real(8), parameter :: hour2day = 1.0d0/day2hour
   real(8), parameter :: day2sec = day2hour*hour2sec
   real(8), parameter :: sec2day = 1.0d0/day2sec

   real(8), parameter :: hour2min = 60.d0
   real(8), parameter :: min2hour = 1.d0/hour2min
   real(8), parameter :: min2day = 1.d0*min2hour*hour2day; 
   !---Conversion Factors---!

   !...Wind Reduction Factor
   !  Factor for reducing wind speed from top of planetary
   !  boundary layer to surface
   real(8), parameter :: windReduction = 0.9d0
   ! Other potential factors below
   ! REAL(8), PARAMETER :: windReduction = 1.0D0
   ! REAL(8), PARAMETER :: windReduction = 0.78D0 !...Powell et all 2003

   !...Conversion between 10 minute and 1 minute wind
   real(8), parameter :: one2ten = 0.8928d0 !... Powell et al 1996
   ! REAL(8), PARAMETER :: 1.00D0 !... HWind
   ! REAL(8), PARAMETER :: 0.8787D0 !.. OceanWeather
   real(8), parameter :: ten2one = 1.0d0/one2ten

   !...Conversion between 30 and 1 minute winds
   real(8), parameter :: thirty2one = 1.165d0 !...Luettich
   real(8), parameter :: one2thirty = 1.0d0/thirty2one

   !...Conversion between 30 and 10 minute winds
   real(8), parameter :: thirty2ten = 1.04d0 !...Luettich
   real(8), parameter :: ten2thirty = 1.0d0/thirty2ten

   !...Wave Wind Multiplier. Applied to winds sent to wave models
   real(8) :: waveWindMultiplier = 1.0d0

   !... Wave Stress Gradient Magnitude CAP
   real(8) :: WaveStressGrad_Cap = 1000.0d0

   real(8), parameter :: DEG2RAD = PI/180.0d0
   real(8), parameter :: RAD2DEG = 180.0d0/PI
   real(8), parameter :: MPERDEG = REarth*PI/180.0d0

   !---Length based factors---!
   real(8), parameter :: nm2m = 1852.0d0 ! nautical miles to meters
   real(8), parameter :: m2nm = 1.0d0/nm2m ! meters to nautical miles
   real(8), parameter :: kt2ms = nm2m/3600.0d0 ! knots to m/s
   real(8), parameter :: ms2kt = 1.0d0/kt2ms ! m/s to knots

   !---Pressure---!
   real(8), parameter :: mb2pa = 100.0d0
   real(8), parameter :: pa2mb = 1.0d0/mb2pa

   !---Dispersion---!
   ! These may be modified by the user at runtime
   real(8) :: Bd = 0.23394d0 !...Exponent of Depth
   real(8) :: Ad = 0.0050189d0 !...Coefficient of Depth
   real(8) :: Cs = 1500.0d0 !...Speed of sound in water
   real(8) :: Cs2 !...Speed of sound in water squared
   real(8), parameter :: TwoB = 2d0*0.4779d0
   real(8), parameter :: GM2 = 3.486d0**2.0d0

end module ADC_CONSTANTS
!-----------------------------------------------------------------------
