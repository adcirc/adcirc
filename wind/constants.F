      !=================================================================
      !=================================================================
      !=================================================================
      !      =====                                           =====
      !      =====              MODULE constants             =====
      !      =====                                           =====
      !=================================================================
      !=================================================================
      !=================================================================

      !=================================================================
      ! Define and initialize some useful global constants.
      !
      ! Compilation:
      !    Throw C preprocessor macro switches (DEFINES) to select
      !    precision:
      !       f90 -O -c constants.F -DREAL8 -DINT2
      !       xlf90 -O -c constants.F -WF,-DREAL8,-DINT2
      !
      ! Revision history:
      !    Date         Programmer                 Description of change
      !    ----         ----------                 ---------------------
      !    05/23/06     Craig  Mattocks, UNC-CEP   Wrote original code
      !=================================================================

      !=======================================
      ! Define and initialize global constants
      !=======================================
      MODULE constants

         !-------------------------------------------------------------
         ! Import custom precision types for cross-platform portability
         !-------------------------------------------------------------
         USE precision

         !--------------------------------------------
         ! Force explicit declaration of all variables
         !--------------------------------------------
         IMPLICIT NONE

         !--------------------------------------------------------------
         ! Retain values of class instance variables between invocations
         !--------------------------------------------------------------
         SAVE

#ifdef REAL8
         !==================
         ! Double precision:
         !==================

         !----------------------
         ! Commonly used numbers
         !----------------------
         REAL (KIND = dp), PARAMETER :: zero  =  0._dp
         REAL (KIND = dp), PARAMETER :: half  =  .5_dp
         REAL (KIND = dp), PARAMETER :: one   =  1._dp
         REAL (KIND = dp), PARAMETER :: two   =  2._dp
         REAL (KIND = dp), PARAMETER :: three =  3._dp
         REAL (KIND = dp), PARAMETER :: four  =  4._dp
         REAL (KIND = dp), PARAMETER :: five  =  5._dp
         REAL (KIND = dp), PARAMETER :: six   =  6._dp
         REAL (KIND = dp), PARAMETER :: seven =  7._dp
         REAL (KIND = dp), PARAMETER :: eight =  8._dp
         REAL (KIND = dp), PARAMETER :: nine  =  9._dp
         REAL (KIND = dp), PARAMETER :: ten   = 10._dp
         REAL (KIND = dp), PARAMETER :: oneTenth         =     .1_dp
         REAL (KIND = dp), PARAMETER :: oneHundredth     =     .01_dp
         REAL (KIND = dp), PARAMETER :: onePointTwo      =    1.2_dp
         REAL (KIND = dp), PARAMETER :: twoPointFive     =    2.5_dp
         REAL (KIND = dp), PARAMETER :: threePointFive     =  3.5_dp
         REAL (KIND = dp), PARAMETER :: twentyFour       =   24._dp
         REAL (KIND = dp), PARAMETER :: twentyFive       =   25._dp
         REAL (KIND = dp), PARAMETER :: seventyFive      =   75._dp
         REAL (KIND = dp), PARAMETER :: oneHundred       =  100._dp
         REAL (KIND = dp), PARAMETER :: oneEighty        =  180._dp
         REAL (KIND = dp), PARAMETER :: threeSixty       =  360._dp
         REAL (KIND = dp), PARAMETER :: oneThousand      = 1000._dp
         REAL (KIND = dp), PARAMETER :: thirtySixHundred = 3600._dp

         !------------
         ! Error flags
         !------------
         REAL (KIND = dp) , PARAMETER :: error       = -99._dp
         CHARACTER (LEN=8), PARAMETER :: missingChar = "Missing!"
         INTEGER          , PARAMETER :: missingInt  = -99
         REAL (KIND = dp) , PARAMETER :: missingReal = -99._dp

         !---------------------
         ! pi and factors of pi
         !---------------------
!        REAL (KIND = dp), PARAMETER :: pi    = two * DASIN(one)
         REAL (KIND = dp), PARAMETER :: pi    =
     &                    3.141592653589793238462643383279502884197_dp
         REAL (KIND = dp), PARAMETER :: pi2   =
     &                    1.57079632679489661923132169163975144209858_dp
         REAL (KIND = dp), PARAMETER :: twopi =
     &                    6.283185307179586476925286766559005768394_dp 

         !----------------------------------------
         ! Degrees <--> radians conversion factors
         !----------------------------------------
         REAL (KIND = dp), PARAMETER  ::  deg2rad = pi / oneEighty
         REAL (KIND = dp), PARAMETER  ::  rad2deg = oneEighty / pi

         !----------------
         ! Gravity (m/s^2)
         !----------------
         REAL (KIND = dp), PARAMETER :: g = 9.80665_dp

         !-------------------------------
         ! Equatorial radius of earth (m)
         !-------------------------------
         REAL (KIND = dp), PARAMETER :: Rearth = 6378135._dp

         !------------------------------------
         ! Rotation frequency of earth (rad/s)
         !------------------------------------
         REAL (KIND = dp), PARAMETER :: omega = twopi / 86164.2_dp

         !------------------------------------------------
         ! Surface pressure for a standard atmosphere (mb)
         !------------------------------------------------
         REAL (KIND = dp), PARAMETER :: Pzero = 1013.00_dp ! Tropical
!        REAL (KIND = dp), PARAMETER :: Pzero = 1013.25_dp ! U.S. 1976 std. atm.

         !---------------
         ! Density of air
         !---------------
!        REAL (KIND = dp), PARAMETER :: RhoAir = 1.15_dp   ! Holland
!        REAL (KIND = dp), PARAMETER :: RhoAir = 1.168_dp  ! T=25    C, P=1000 mb
         REAL (KIND = dp), PARAMETER :: RhoAir = 1.1774_dp ! T=79.79 F, P=1013 mb

         !----------------------------------------------------------
         ! Factor for reducing wind speed from top of PBL to surface
         !----------------------------------------------------------
!        REAL (KIND = dp), PARAMETER :: windReduction = one      ! no-op
!        REAL (KIND = dp), PARAMETER :: windReduction = .78_dp   ! Powell et al., 2003
         REAL (KIND = dp), PARAMETER :: windReduction = .90_dp   ! Luettich
!        REAL (KIND = dp), PARAMETER :: windReduction = 1.05_dp  ! experimental

         !----------------------------------------------
         ! Wind speed time-averaging conversion factors:
         !    one2ten =  1-minute avg --> 10-minute avg
         !    ten2one = 10-minute avg -->  1-minute avg
         !----------------------------------------------
!        REAL (KIND = dp), PARAMETER :: one2ten = one            ! Hwinds
         REAL (KIND = dp), PARAMETER :: one2ten = .8928_dp       ! Powell et al., 1996
!        REAL (KIND = dp), PARAMETER :: one2ten = .8787_dp       ! Ocean Weather, Inc.
         REAL (KIND = dp), PARAMETER :: ten2one = one/one2ten

         REAL (KIND = dp), PARAMETER :: thirty2one = 1.165_dp    ! Luettich
         REAL (KIND = dp), PARAMETER :: one2thirty = one/thirty2one

         REAL (KIND = dp), PARAMETER :: thirty2ten = 1.04_dp     ! Luettich
         REAL (KIND = dp), PARAMETER :: ten2thirty = one/thirty2ten

         !---------------------------------------------------
         ! Cosine damping distance for translational velocity
         !---------------------------------------------------
         REAL (KIND = dp), PARAMETER :: dampRadii  = five
         REAL (KIND = dp), PARAMETER :: dampLength = ten

         !-------------------------
         ! Nominal density of water
         !-------------------------
         REAL (KIND = dp), PARAMETER :: RhoWat0 = oneThousand
         REAL (KIND = dp), PARAMETER :: RhoWatG = RhoWat0 * g

         !-----------------------------------
         ! Sigma T value of reference density
         !-----------------------------------
         REAL (KIND = dp), PARAMETER :: SigT0 = RhoWat0 - oneThousand

         !------------------------------------
         ! Parameters used in barrier overflow 
         !------------------------------------
         REAL (KIND = dp), PARAMETER :: barmin = .04_dp

         !------------------------
         ! Time conversion factors
         !------------------------
         REAL (KIND = dp), PARAMETER :: hour2sec = thirtySixHundred
         REAL (KIND = dp), PARAMETER :: sec2hour = one / hour2sec
         REAL (KIND = dp), PARAMETER :: day2hour = twentyFour
         REAL (KIND = dp), PARAMETER :: hour2day = one / day2hour
         REAL (KIND = dp), PARAMETER :: day2sec  = day2hour * hour2sec
         REAL (KIND = dp), PARAMETER :: sec2day  = one / day2sec

         !----------------------------------------------
         ! nautical miles <--> meters conversion factors
         !----------------------------------------------
         REAL (KIND = dp), PARAMETER :: nm2m = 1852._dp
         REAL (KIND = dp), PARAMETER :: m2nm = one / nm2m

         !----------------------------------
         ! knots <--> m/s conversion factors
         !----------------------------------
         REAL (KIND = dp), PARAMETER :: kt2ms = nm2m / hour2sec
         REAL (KIND = dp), PARAMETER :: ms2kt = one / kt2ms

         !----------------------------------
         ! mb <--> pascals conversion factor
         !----------------------------------
         REAL (KIND = dp), PARAMETER :: mb2pa = oneHundred
         REAL (KIND = dp), PARAMETER :: pa2mb = oneHundredth

#else
         !==================
         ! Single precision:
         !==================

         !----------------------
         ! Commonly used numbers
         !----------------------
         REAL (KIND = sp), PARAMETER :: zero  =  0._sp
         REAL (KIND = sp), PARAMETER :: half  =  .5_sp
         REAL (KIND = sp), PARAMETER :: one   =  1._sp
         REAL (KIND = sp), PARAMETER :: two   =  2._sp
         REAL (KIND = sp), PARAMETER :: three =  3._sp
         REAL (KIND = sp), PARAMETER :: four  =  4._sp
         REAL (KIND = sp), PARAMETER :: five  =  5._sp
         REAL (KIND = sp), PARAMETER :: six   =  6._sp
         REAL (KIND = sp), PARAMETER :: seven =  7._sp
         REAL (KIND = sp), PARAMETER :: eight =  8._sp
         REAL (KIND = sp), PARAMETER :: nine  =  9._sp
         REAL (KIND = sp), PARAMETER :: ten   = 10._sp
         REAL (KIND = sp), PARAMETER :: oneTenth         =     .1_sp
         REAL (KIND = sp), PARAMETER :: oneHundredth     =     .01_sp
         REAL (KIND = sp), PARAMETER :: onePointTwo      =    1.2_sp
         REAL (KIND = sp), PARAMETER :: twoPointFive     =    2.5_sp
         REAL (KIND = sp), PARAMETER :: threePointFive     =  3.5_sp
         REAL (KIND = sp), PARAMETER :: twentyFour       =   24._sp
         REAL (KIND = sp), PARAMETER :: twentyFive       =   25._sp
         REAL (KIND = sp), PARAMETER :: seventyFive      =   75._sp
         REAL (KIND = sp), PARAMETER :: oneHundred       =  100._sp
         REAL (KIND = sp), PARAMETER :: oneEighty        =  180._sp
         REAL (KIND = sp), PARAMETER :: threeSixty       =  360._sp
         REAL (KIND = sp), PARAMETER :: oneThousand      = 1000._sp
         REAL (KIND = sp), PARAMETER :: thirtySixHundred = 3600._sp

         !------------
         ! Error flags
         !------------
         REAL (KIND = sp) , PARAMETER :: error       = -99._sp
         CHARACTER (LEN=8), PARAMETER :: missingChar = "Missing!"
         INTEGER          , PARAMETER :: missingInt  = -99
         REAL (KIND = sp) , PARAMETER :: missingReal = -99._sp

         !---------------------
         ! pi and factors of pi
         !---------------------
!        REAL (KIND = sp), PARAMETER :: pi    = two *  ASIN(1._sp)
         REAL (KIND = sp), PARAMETER :: pi    =
     &                    3.141592653589793238462643383279502884197_sp
         REAL (KIND = sp), PARAMETER :: pi2   =
     &                    1.57079632679489661923132169163975144209858_sp
         REAL (KIND = sp), PARAMETER :: twopi =
     &                    6.283185307179586476925286766559005768394_sp 

         !----------------------------------------
         ! Degrees <--> radians conversion factors
         !----------------------------------------
         REAL (KIND = sp), PARAMETER :: deg2rad = pi / oneEighty
         REAL (KIND = sp), PARAMETER :: rad2deg = oneEighty / pi

         !----------------
         ! Gravity (m/s^2)
         !----------------
         REAL (KIND = sp), PARAMETER :: g = 9.80665_sp

         !-------------------------------
         ! Equatorial radius of earth (m)
         !-------------------------------
         REAL (KIND = sp), PARAMETER :: Rearth = 6378135._sp

         !------------------------------------
         ! Rotation frequency of earth (rad/s)
         !------------------------------------
         REAL (KIND = sp), PARAMETER :: omega = twopi / 86164.2_sp

         !------------------------------------------------
         ! Surface pressure for a standard atmosphere (mb)
         !------------------------------------------------
         REAL (KIND = sp), PARAMETER :: Pzero = 1013.00_sp ! Tropical
!        REAL (KIND = sp), PARAMETER :: Pzero = 1013.25_sp ! U.S. 1976 std. atm.

         !---------------
         ! Density of air
         !---------------
!        REAL (KIND = sp), PARAMETER :: RhoAir = 1.15_sp   ! Holland
!        REAL (KIND = sp), PARAMETER :: RhoAir = 1.168_sp  ! T=25    C, P=1000 mb
         REAL (KIND = sp), PARAMETER :: RhoAir = 1.1774_sp ! T=79.79 F, P=1013 mb

         !----------------------------------------------------------
         ! Factor for reducing wind speed from top of PBL to surface
         !----------------------------------------------------------
!        REAL (KIND = sp), PARAMETER :: windReduction = one      ! no-op
!        REAL (KIND = sp), PARAMETER :: windReduction = .78_sp   ! Powell et al., 2003
         REAL (KIND = sp), PARAMETER :: windReduction = .90_sp   ! Luettich
!        REAL (KIND = sp), PARAMETER :: windReduction = 1.05_sp  ! experimental

         !----------------------------------------------
         ! Wind speed time-averaging conversion factors:
         !    one2ten =  1-minute avg --> 10-minute avg
         !    ten2one = 10-minute avg -->  1-minute avg
         !----------------------------------------------
!        REAL (KIND = sp), PARAMETER :: one2ten = one            ! Hwinds
         REAL (KIND = sp), PARAMETER :: one2ten = .8928_sp       ! Powell et al., 1996
!        REAL (KIND = sp), PARAMETER :: one2ten = .8787_sp       ! Ocean Weather, Inc.
         REAL (KIND = sp), PARAMETER :: ten2one = one/one2ten

         REAL (KIND = sp), PARAMETER :: thirty2one = 1.165_sp    ! Luettich
         REAL (KIND = sp), PARAMETER :: one2thirty = one/thirty2one

         REAL (KIND = sp), PARAMETER :: thirty2ten = 1.04_sp     ! Luettich
         REAL (KIND = sp), PARAMETER :: ten2thirty = one/thirty2ten

         !---------------------------------------------------
         ! Cosine damping distance for translational velocity
         !---------------------------------------------------
         REAL (KIND = sp), PARAMETER :: dampRadii  = five
         REAL (KIND = sp), PARAMETER :: dampLength = ten

         !-------------------------
         ! Nominal density of water
         !-------------------------
         REAL (KIND = sp), PARAMETER :: RhoWat0 = oneThousand
         REAL (KIND = sp), PARAMETER :: RhoWatG = RhoWat0 * g

         !-----------------------------------
         ! Sigma T value of reference density
         !-----------------------------------
         REAL (KIND = sp), PARAMETER :: SigT0 = RhoWat0 - oneThousand

         !------------------------------------
         ! Parameters used in barrier overflow 
         !------------------------------------
         REAL (KIND = sp), PARAMETER :: barmin = .04_sp

         !------------------------
         ! Time conversion factors
         !------------------------
         REAL (KIND = sp), PARAMETER :: hour2sec = thirtySixHundred
         REAL (KIND = sp), PARAMETER :: sec2hour = one / hour2sec
         REAL (KIND = sp), PARAMETER :: day2hour = twentyFour
         REAL (KIND = sp), PARAMETER :: hour2day = one / day2hour
         REAL (KIND = sp), PARAMETER :: day2sec  = day2hour * hour2sec
         REAL (KIND = sp), PARAMETER :: sec2day  = one / day2sec

         !----------------------------------------------
         ! nautical miles <--> meters conversion factors
         !----------------------------------------------
         REAL (KIND = sp), PARAMETER :: nm2m = 1852._sp
         REAL (KIND = sp), PARAMETER :: m2nm = one / nm2m

         !----------------------------------
         ! knots <--> m/s conversion factors
         !----------------------------------
         REAL (KIND = sp), PARAMETER :: kt2ms = nm2m / hour2sec
         REAL (KIND = sp), PARAMETER :: ms2kt = one / kt2ms

         !----------------------------------
         ! mb <--> pascals conversion factor
         !----------------------------------
         REAL (KIND = sp), PARAMETER :: mb2pa = oneHundred
         REAL (KIND = sp), PARAMETER :: pa2mb = oneHundredth

#endif

         !----------------
         ! Date parameters
         !----------------
         CHARACTER (LEN=9), PARAMETER :: months(12) =
     &              (/'January  ','February ','March    ','April    ',
     &                'May      ','June     ','July     ','August   ',
     &                'September','October  ','November ','December '/)

      END MODULE constants
