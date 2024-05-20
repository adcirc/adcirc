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
      MODULE ADC_CONSTANTS
!-----------------------------------------------------------------------

          IMPLICIT NONE

          !---Physical Constants---!
        
          !...Default gravitational acceleration
          ! Note that the user can update this in the
          ! control file so this isn't made constant
          REAL(8) :: g = 9.80665

          !...Nominal density of water
          REAL(8), PARAMETER :: RhoWat0 = 1000.0D0

          !...Sigma T value of reference density
          REAL(8), PARAMETER :: SigT0 = RHOWAT0 - 1000.0D0

          !...Background Atmospheric pressure in mb
          REAL(8), PARAMETER :: PRBCKGRND = 1013.0D0
          ! REAL(8), PARAMETER :: PRBCKGRND = 1013.25D0 !...U.S. 1976 Std Atm

          !...Background temperature (kelvin)
          REAL(8), PARAMETER :: TBCKGRND = 288.15D0

          !...Specific gas constant
          REAL(8), PARAMETER :: RAir = 287.058D0

          !...Ideal gas law kg/m^3
          !...Note this is the default value but may be modified
          ! by the user at runtime
          REAL(8) :: RhoAir = 1.293D0
          !...There are other options which may be used but 
          ! the default is above
          ! REAL(8), PARAMETER :: RhoAir = 1.15D0 !...Holland
          ! REAL(8), PARAMETER :: RhoAir = 1.1774D0 !... T=79.79, P=1013D0

          REAL(8), PARAMETER :: PI = 3.141592653589793D0
          REAL(8), PARAMETER :: TWOPI=PI*2.D0
          REAL(8), PARAMETER :: HFPI=PI/2.D0
          REAL(8), PARAMETER :: e = 2.718281828459045D0

          !...Radius of earth (m)
          REAL(8), PARAMETER :: Rearth = 6378206.4d0

          REAL(8), PARAMETER :: omega = 7.29212d-5

          !---Time based factors---!
          REAL(8), PARAMETER :: hour2sec = 3600.0D0
          REAL(8), PARAMETER :: sec2hour = 1.0D0 / hour2sec
          REAL(8), PARAMETER :: day2hour = 24.0D0
          REAL(8), PARAMETER :: hour2day = 1.0D0 / day2hour
          REAL(8), PARAMETER :: day2sec = day2hour * hour2sec
          REAL(8), PARAMETER :: sec2day = 1.0D0/day2sec 

          !---Conversion Factors---!
          
          !...Wind Reduction Factor
          !  Factor for reducing wind speed from top of planetary
          !  boundary layer to surface
          REAL(8), PARAMETER :: windReduction = 0.9D0
          ! Other potential factors below
          ! REAL(8), PARAMETER :: windReduction = 1.0D0
          ! REAL(8), PARAMETER :: windReduction = 0.78D0 !...Powell et all 2003

          !...Conversion between 10 minute and 1 minute wind 
          REAL(8), PARAMETER :: one2ten = 0.8928D0 !... Powell et al 1996
          ! REAL(8), PARAMETER :: 1.00D0 !... HWind
          ! REAL(8), PARAMETER :: 0.8787D0 !.. OceanWeather
          REAL(8), PARAMETER :: ten2one = 1.0D0/one2ten

          !...Conversion between 30 and 1 minute winds
          REAL(8), PARAMETER :: thirty2one = 1.165D0 !...Luettich
          REAL(8), PARAMETER :: one2thirty = 1.0D0 / thirty2one

          !...Conversion between 30 and 10 minute winds
          REAL(8), PARAMETER :: thirty2ten = 1.04D0 !...Luettich
          REAL(8), PARAMETER :: ten2thirty = 1.0D0 / thirty2ten

          !...Wave Wind Multiplier. Applied to winds sent to wave models
          REAL(8) :: waveWindMultiplier = 1.0D0
          
          !... Wave Stress Gradient Magnitude CAP
          REAL(8) :: WaveStressGrad_Cap = 1000.0D0

          REAL(8), PARAMETER :: DEG2RAD = PI/180.0D0
          REAL(8), PARAMETER :: RAD2DEG = 180.0D0/PI
          REAL(8), PARAMETER :: MPERDEG  = REarth * PI / 180.0d0

          !---Length based factors---!
          REAL(8), PARAMETER :: nm2m = 1852.0D0 ! nautical miles to meters
          REAL(8), PARAMETER :: m2nm = 1.0D0/nm2m ! meters to nautical miles
          REAL(8), PARAMETER :: kt2ms = nm2m / 3600.0D0 ! knots to m/s
          REAL(8), PARAMETER :: ms2kt = 1.0D0 / kt2ms ! m/s to knots

          !---Pressure---!
          REAL(8), PARAMETER :: mb2pa = 100.0D0
          REAL(8), PARAMETER :: pa2mb = 1.0D0 / mb2pa

          !---Dispersion---!
          ! These may be modified by the user at runtime
          REAL(8) :: Bd = 0.23394D0 !...Exponent of Depth
          REAL(8) :: Ad = 0.0050189D0 !...Coefficient of Depth
          REAL(8) :: Cs = 1500.0D0 !...Speed of sound in water
          REAL(8) :: Cs2 !...Speed of sound in water squared
          REAL(8), PARAMETER :: TwoB = 2D0 * 0.4779D0
          REAL(8), PARAMETER :: GM2 = 3.486D0**2.0D0

!-----------------------------------------------------------------------
      END MODULE ADC_CONSTANTS
!-----------------------------------------------------------------------
