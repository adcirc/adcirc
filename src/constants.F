!-----------------------------------------------------------------------
!  MODULE CONSTANTS
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
      MODULE CONSTANTS
!-----------------------------------------------------------------------
          USE SIZES, ONLY: SZ
          IMPLICIT NONE

          !---Physical Constants---!
        
          !...Default gravitational acceleration
          ! Note that the user can update this in the
          ! control file so this isn't made constant
          REAL(SZ) :: g = 9.80665

          !...Nominal density of water
          REAL(SZ), PARAMETER :: RhoWat0 = 1000.0D0

          !...Sigma T value of reference density
          REAL(SZ), PARAMETER :: SigT0 = RHOWAT0 - 1000.0D0

          !...Background Atmospheric pressure in mb
          REAL(SZ), PARAMETER :: PRBCKGRND = 1013.0D0
          ! REAL(SZ), PARAMETER :: PRBCKGRND = 1013.25D0 !...U.S. 1976 Std Atm

          !...Background temperature (kelvin)
          REAL(SZ), PARAMETER :: TBCKGRND = 288.15D0

          !...Specific gas constant
          REAL(SZ), PARAMETER :: RAir = 287.058D0

          !...Ideal gas law kg/m^3
          !...Note this is the default value but may be modified
          ! by the user at runtime
          REAL(SZ) :: RhoAir = 1.293D0
          !...There are other options which may be used but 
          ! the default is above
          ! REAL(SZ), PARAMETER :: RhoAir = 1.15D0 !...Holland
          ! REAL(SZ), PARAMETER :: RhoAir = 1.1774D0 !... T=79.79, P=1013D0

          REAL(SZ), PARAMETER :: PI = 3.141592653589793D0
          REAL(SZ), PARAMETER :: e = 2.718281828459045D0

          !...Radius of earth (m)
          REAL(SZ), PARAMETER :: Rearth = 6378206.4d0

          REAL(SZ), PARAMETER :: omega = 7.29212d-5

          !---Time based factors---!
          REAL(SZ), PARAMETER :: hour2sec = 3600.0D0
          REAL(SZ), PARAMETER :: sec2hour = 1.0D0 / hour2sec
          REAL(SZ), PARAMETER :: day2hour = 24.0D0
          REAL(SZ), PARAMETER :: hour2day = 1.0D0 / day2hour
          REAL(SZ), PARAMETER :: day2sec = day2hour * hour2sec
          REAL(SZ), PARAMETER :: sec2day = 1.0D0/day2sec 

          !---Conversion Factors---!
          
          !...Wind Reduction Factor
          !  Factor for reducing wind speed from top of planetary
          !  boundary layer to surface
          REAL(SZ), PARAMETER :: windReduction = 0.9D0
          ! Other potential factors below
          ! REAL(SZ), PARAMETER :: windReduction = 1.0D0
          ! REAL(SZ), PARAMETER :: windReduction = 0.78D0 !...Powell et all 2003

          !...Conversion between 10 minute and 1 minute wind 
          REAL(SZ), PARAMETER :: one2ten = 0.8928D0 !... Powell et al 1996
          ! REAL(SZ), PARAMETER :: 1.00D0 !... HWind
          ! REAL(SZ), PARAMETER :: 0.8787D0 !.. OceanWeather
          REAL(SZ), PARAMETER :: ten2one = 1.0D0/one2ten

          !...Conversion between 30 and 1 minute winds
          REAL(SZ), PARAMETER :: thirty2one = 1.165D0 !...Luettich
          REAL(SZ), PARAMETER :: one2thirty = 1.0D0 / thirty2one

          !...Conversion between 30 and 10 minute winds
          REAL(SZ), PARAMETER :: thirty2ten = 1.04D0 !...Luettich
          REAL(SZ), PARAMETER :: ten2thirty = 1.0D0 / thirty2ten

          !...Wave Wind Multiplier. Applied to winds sent to wave models
          REAL(SZ) :: waveWindMultiplier = 1.0D0

          REAL(SZ), PARAMETER :: DEG2RAD = PI/180.0D0
          REAL(SZ), PARAMETER :: RAD2DEG = 180.0D0/PI

          !---Length based factors---!
          REAL(SZ), PARAMETER :: nm2m = 1852.0D0 ! nautical miles to meters
          REAL(SZ), PARAMETER :: m2nm = 1.0D0/nm2m ! meters to nautical miles
          REAL(SZ), PARAMETER :: kt2ms = nm2m / 3600.0D0 ! knots to m/s
          REAL(SZ), PARAMETER :: ms2kt = 1.0D0 / kt2ms ! m/s to knots

          !---Pressure---!
          REAL(SZ), PARAMETER :: mb2pa = 100.0D0
          REAL(SZ), PARAMETER :: pa2mb = 1.0D0 / mb2pa

          !---Dispersion---!
          ! These may be modified by the user at runtime
          REAL(SZ) :: Bd = 0.23394D0 !...Exponent of Depth
          REAL(SZ) :: Ad = 0.0050189D0 !...Coefficient of Depth
          REAL(SZ) :: Cs = 1500.0D0 !...Speed of sound in water
          REAL(SZ) :: Cs2 !...Speed of sound in water squared
          REAL(SZ), PARAMETER :: TwoB = 2D0 * 0.4779D0
          REAL(SZ), PARAMETER :: GM2 = 3.486D0**2.0D0

!-----------------------------------------------------------------------
      END MODULE CONSTANTS
!-----------------------------------------------------------------------
