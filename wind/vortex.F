      !=================================================================
      !=================================================================
      !=================================================================
      !      =====                                           =====
      !      =====               MODULE vortex               =====
      !      =====                                           =====
      !=================================================================
      !=================================================================
      !=================================================================

      !=================================================================
      ! This module allows one to create an asymmetric hurricane vortex
      ! object from National Hurricane Center forecast advisories, then
      ! compute wind and pressure fields suitable for forcing a storm
      ! surge forecast model. The shape and intensity of the vortex can
      ! be continuously adjusted during the forecast period, as it moves
      ! through the model domain.
      !
      ! Revision history:
      !    Date        Programmer                 Description of change
      !    ----        ----------                 ---------------------
      !    05/23/06    Craig  Mattocks, UNC-CEP   Wrote original code
      !    06/30/06    Cristina Forbes, UNC-CEP   Tested in ADCIRC model
      !    08/25/06    Craig  Mattocks, UNC-CEP   Applied wind reduction:
      !                                           top of SFC layer -> SFC
      !                                           in subroutine uvp.
      !    09/12/06    Craig  Mattocks, UNC-CEP   Subtracted translational
      !                                           wind speed from Vmax in
      !                                           nws9get.
      !    05/12/09    Robert Weaver, UNC-IMS     Modified damping of translational 
      !                                           velocity based on ratio of V/Vmax
      !    05/19/09    Cristina Forbes, UNC-IMS   Implemented old bug fix: units conversion 
      !    05/22/09    Cristina Forbes, UNC-IMS   Reverted damping of translational
      !                                           velocity back to the original formulation
      !                                           developed by Craig Mattocks
      !    05/26/09    Cristina Forbes, UNC-IMS   Fixed implementation of Rick Luettich 
      !                                           V/Vmax tapering formulation (as in NWS=8) 
      !                                           for future experimentation - not activated
      !    06/2009      Robert Weaver, UNC-IMS    changed the gradient wind formula to 
      !                                           use Vmax instead of pressure difference
      !                                           to enable easier manipulation of storm 
      !                                           characteristics
      !    07/2009      Robert Weaver, UNC-IMS    using modules from this code, 
      !                                           wrote a preprocessor to compute Rmax and B
      !                                           prior to running adcirc, and write out an
      !                                           input file like the ATCF file with these
      !                                           variables to make it easier to manipulate
      !                                           the storm parameters.  Changed this code
      !                                           to read in that input file and use the values
      !                                           found there 
      !                                           
      !
      !=================================================================
      MODULE vortex

         !-------------------------------------------------------------
         ! Import custom precision types for cross-platform portability
         !-------------------------------------------------------------
         USE precision

         !------------------------
         ! Import global constants
         !------------------------
         USE constants

         !--------------------------------------------
         ! Force explicit declaration of all variables
         !--------------------------------------------
         IMPLICIT NONE

         !--------------------------------------------------------------
         ! Retain values of class instance variables between invocations
         !--------------------------------------------------------------
         SAVE

         !-----------------------
         ! Set access permissions
         !-----------------------
       PRIVATE :: calcShape, calcRmax, calcRmax_givenRmax,
     &              poly,fitPoly, spline,fitSpline,
     &              setStormName       , setAdvisory        ,
     &              setDate            ,
     &              setAmbientPressure , setCentralPressure ,
     &              setCenterLatitude  , setCenterLongitude ,
     &              setMaxWinds        , setShapeParameter  ,
     &              setCoriolisForce   , setRadialWind      ,
     &              setWindRadii       , setWindRadius      ,
     &              setPolyCoefficients, setPolyCoefficient ,
     &              setRadiiMaxWinds   , setRadiusMaxWinds  ,
     &              setQuadrant        , setRegressionCoefficient,
     &              setFitMethod

        PUBLIC :: newVortex, newVortex_rmax_b_IN, newVortex_rmax_b_OUT,
     &              transVortex, moveVortex, printVortex,
     &              coriolis, Ph, Pasym, Rmw, Vh, Vasym, uvp, uvtrans,
     &              latlon2xy, xy2latlon, Vmaxw ,
     &              getStormName        , getAdvisory        ,
     &              printDate           ,
     &              getYear, getMonth   , getDay, getHour    ,
     &              getForecastIncrement,
     &              getAmbientPressure  , getCentralPressure ,
     &              getCenterLatitude   , getCenterLongitude ,
     &              getMaxWinds         , getShapeParameter  ,
     &              getCoriolisForce    , getRadialWind      ,
     &              getWindRadii        , getWindRadius      ,
     &              getPolyOrder        ,
     &              getPolyCoefficients , getPolyCoefficient ,
     &              getAzimuthalAngles  , getAzimuthalAngle  ,
     &              getRadiiMaxWinds    , getRadiusMaxWinds  ,
     &              getQuadrant         , getRegressionCoefficient,
     &              getFitMethod

         !----------------------------
         ! Hurricane vortex parameters
         !----------------------------
         PRIVATE
         PUBLIC  :: POLY_FIT, SPLINE_FIT

         INTEGER, PARAMETER :: POLY_FIT   = 1    ! Solve with polynomial curve fit
         INTEGER, PARAMETER :: SPLINE_FIT = 2    ! Solve with cubic splines

         INTEGER, PARAMETER :: nQuads  = 4       ! Number of quadrants for
                                                 ! which wind radii are
                                                 ! provided

         INTEGER, PARAMETER :: nPoints = nQuads+2! Number of (theta, Rmax)
                                                 ! points for curve fit

         INTEGER, PARAMETER :: order   = 4       ! Order of polynomial for
                                                 ! Rmax vs. theta curve fit

                                                 ! Azimuthal angles (degrees)
         REAL(sz), DIMENSION(nPoints), PARAMETER :: theta = (/
     &                                                          0._dp,
     &                                                         45._dp,
     &                                                        135._dp,
     &                                                        225._dp,
     &                                                        315._dp,
     &                                                        360._dp
     &                                                      /)

         REAL(sz), DIMENSION(nPoints) :: Vmax_wind   ! maximum wind speed (knots)

         REAL(sz), DIMENSION(nPoints) :: Rmax    ! Radius of maximum winds
                                                 ! (nautical miles)

         CHARACTER(LEN=64) :: stormName          ! Name of storm

         INTEGER  :: advisory                    ! NHC advisory or storm number

         INTEGER  :: year                        ! Year  of advisory  (4 digits)
         INTEGER  :: month                       ! Month of advisory  (2 digits)
         INTEGER  :: day                         ! Day   of advisory  (2 digits)
         INTEGER  :: hour                        ! Hour  of advisory  (2 digits)
         INTEGER  :: hourFcst                    ! Forecast increment (2 digits)

         REAL(sz) :: Pn                          ! Ambient surface pressure (mb)
         REAL(sz) :: Pc                          ! Surface pressure at center of
                                                 ! storm (mb)
         REAL(sz) :: cLat                        ! Latitude  of storm center
                                                 ! (degrees north)
         REAL(sz) :: cLon                        ! Longitude of storm center
                                                 ! (degrees east )
         REAL(sz) :: Vmax                        ! Max sustained wind velocity
                                                 ! in storm (knots)
         REAL(sz) :: B                           ! Exponential shape parameter
         REAL(sz) :: corio                       ! Coriolis force (1/s)
         REAL(sz) :: Vr                          ! Velocity @ wind radii (knots)
         REAL(sz), DIMENSION(nQuads) :: radius   ! Wind radii - the distance
                                                 ! winds of velocity Vr extend
                                                 ! outward from center of storm
                                                 ! (nautical miles)
         REAL(sz), DIMENSION(order+1) :: a       ! Coefficients of the Rmax
                                                 ! vs. theta polynomial
         REAL(sz) :: rSquare                     ! Regression coefficient for
                                                 ! polynomial curve fit
         INTEGER  :: quad                        ! Quadrant counter
         INTEGER  :: fitMethod                   ! Curve fit method
         REAL, DIMENSION(nPoints) :: work        ! Work array
         REAL, DIMENSION(nPoints) :: d2ydx2      ! Cubic spline 2nd derivatives

         CONTAINS

            !=================
            ! Vortex functions
            !=================

            !=================================================================
            ! Create a new Vortex object.
            !
            ! On input:
            !    name         Name of storm
            !    adv          NHC advisory or storm number
            !    yyyy         Year  of advisory  (4 digits)
            !    mm           Month of advisory  (2 digits)
            !    dd           Day   of advisory  (2 digits)
            !    hh           Hour  of advisory  (2 digits)
            !    hf           Forecast increment (2 digits)
            !    Pn           Ambient surface pressure (mb)
            !    Pc           Surface pressure at center of storm (mb)
            !    cLat         Latitude  of storm center (degrees north)
            !    cLon         Longitude of storm center (degrees east )
            !    Vmax         Max sustained wind velocity in storm (knots)
            !    Vr           Velocity at wind radii (knots)
            !    r            Wind radii - the distance winds of velocity Vr
            !                 extend outward from center of storm
            !                 (nautical miles)
            !
            ! On output:
            !    A new vortex is created with essential parameters calculated.
            !  Return Rmax in the 4 quadrants and the Holland B parameter
            !=================================================================
            SUBROUTINE newVortex(name, adv, yyyy,mm,dd,hh,hf,
     &                      Pn,Pc, cLat,cLon, Vmax,Vr, r)
               CHARACTER(LEN=*), INTENT(IN) :: name
               INTEGER , INTENT(IN) :: adv
               INTEGER , INTENT(IN) :: yyyy
               INTEGER , INTENT(IN) :: mm
               INTEGER , INTENT(IN) :: dd
               INTEGER , INTENT(IN) :: hh
               INTEGER , INTENT(IN) :: hf
               REAL(sz), INTENT(IN) :: Pn
               REAL(sz), INTENT(IN) :: Pc
               REAL(sz), INTENT(IN) :: cLat
               REAL(sz), INTENT(IN) :: cLon
               REAL(sz), INTENT(IN) :: Vmax
               REAL(sz), INTENT(IN) :: Vr
               REAL(sz), INTENT(IN), DIMENSION(nQuads) :: r

               CALL setStormName(name)
               CALL setAdvisory(adv)
               CALL setDate(yyyy, mm, dd, hh, hf)

               CALL setAmbientPressure(Pn)
               CALL setCentralPressure(Pc)

               CALL setCenterLatitude (cLat)
               CALL setCenterLongitude(cLon)

               CALL setMaxWinds  (Vmax)
               CALL setRadialWind(Vr)
               CALL setWindRadii (r)

               CALL setCoriolisForce(coriolis(cLat))
               CALL setShapeParameter(calcShape(Pn, Pc, Vmax))
!              CALL setFitMethod(POLY_FIT)
               CALL setFitMethod(SPLINE_FIT)
               CALL calcRmax()

 
            END SUBROUTINE newVortex


            !=================================================================
            ! Create a new Vortex object.
            !
            ! On input:
            !    name         Name of storm
            !    adv          NHC advisory or storm number
            !    yyyy         Year  of advisory  (4 digits)
            !    mm           Month of advisory  (2 digits)
            !    dd           Day   of advisory  (2 digits)
            !    hh           Hour  of advisory  (2 digits)
            !    hf           Forecast increment (2 digits)
            !    Pn           Ambient surface pressure (mb)
            !    Pc           Surface pressure at center of storm (mb)
            !    cLat         Latitude  of storm center (degrees north)
            !    cLon         Longitude of storm center (degrees east )
            !    Vmax         Max sustained wind velocity in storm (knots)
            !    Vr           Velocity at wind radii (knots)
            !    rmaxw        precomputed Radius to Max Winds in the 4
            !                 quadrants (nautical Miiles)
            !    HollB        precomputed Holland B shape parameter 
            !
            ! On output:
            !    A new vortex is created with essential parameters calculated.
            !=================================================================
            SUBROUTINE newVortex_rmax_b_IN(name, adv, yyyy,mm,dd,hh,hf,
     &                      Pn,Pc, cLat,cLon, Vmax,rmaxw, HollB)
               CHARACTER(LEN=*), INTENT(IN) :: name
               INTEGER , INTENT(IN) :: adv
               INTEGER , INTENT(IN) :: yyyy
               INTEGER , INTENT(IN) :: mm
               INTEGER , INTENT(IN) :: dd
               INTEGER , INTENT(IN) :: hh
               INTEGER , INTENT(IN) :: hf
               REAL(sz), INTENT(IN) :: Pn
               REAL(sz), INTENT(IN) :: Pc
               REAL(sz), INTENT(IN) :: cLat
               REAL(sz), INTENT(IN) :: cLon
               REAL(sz), INTENT(IN) :: Vmax
!               REAL(sz), INTENT(IN) :: Vr
!               REAL(sz), INTENT(IN), DIMENSION(nQuads) :: r
               REAL(sz), INTENT(IN) :: HollB
               REAL(sz), INTENT(IN),DIMENSION(nQuads) :: rmaxw

               CALL setStormName(name)
               CALL setAdvisory(adv)
               CALL setDate(yyyy, mm, dd, hh, hf)

               CALL setAmbientPressure(Pn)
               CALL setCentralPressure(Pc)

               CALL setCenterLatitude (cLat)
               CALL setCenterLongitude(cLon)

               CALL setMaxWinds  (Vmax)
C               CALL setRadialWind(Vr)
C               CALL setWindRadii (r)

               CALL setCoriolisForce(coriolis(cLat))
               CALL setShapeParameter(HollB)
!              CALL setFitMethod(POLY_FIT)
               CALL setFitMethod(SPLINE_FIT)
               CALL calcRmax_givenRmax(rmaxw)

            END SUBROUTINE newVortex_rmax_b_IN
                                                                
            !=================================================================
            ! Create a new Vortex object.
            !
            ! On input:
            !    name         Name of storm
            !    adv          NHC advisory or storm number
            !    yyyy         Year  of advisory  (4 digits)
            !    mm           Month of advisory  (2 digits)
            !    dd           Day   of advisory  (2 digits)
            !    hh           Hour  of advisory  (2 digits)
            !    hf           Forecast increment (2 digits)
            !    Pn           Ambient surface pressure (mb)
            !    Pc           Surface pressure at center of storm (mb)
            !    cLat         Latitude  of storm center (degrees north)
            !    cLon         Longitude of storm center (degrees east )
            !    Vmax         Max sustained wind velocity in storm (knots)
            !    Vr           Velocity at wind radii (knots)
            !    r            Wind radii - the distance winds of velocity Vr
            !                 extend outward from center of storm
            !                 (nautical miles)
            !
            ! On output:
            !    A new vortex is created with essential parameters calculated.
            !  Return Rmax in the 4 quadrants and the Holland B parameter
            !   rmaxw        Radus to Max winds in the 4 quadrants 
            !   HollB        Holland B shape parameter
            !=================================================================
            SUBROUTINE newVortex_rmax_b_OUT(name, adv, yyyy,mm,dd,hh,hf,
     &                      Pn,Pc, cLat,cLon, Vmax,Vr, r,rmaxw,HollB)
               CHARACTER(LEN=*), INTENT(IN) :: name
               INTEGER  :: i
               INTEGER , INTENT(IN) :: adv
               INTEGER , INTENT(IN) :: yyyy
               INTEGER , INTENT(IN) :: mm
               INTEGER , INTENT(IN) :: dd
               INTEGER , INTENT(IN) :: hh
               INTEGER , INTENT(IN) :: hf
               REAL(sz), INTENT(IN) :: Pn
               REAL(sz), INTENT(IN) :: Pc
               REAL(sz), INTENT(IN) :: cLat
               REAL(sz), INTENT(IN) :: cLon
               REAL(sz), INTENT(IN) :: Vmax
               REAL(sz), INTENT(IN) :: Vr
               REAL(sz), INTENT(IN), DIMENSION(nQuads) :: r
               REAL(sz), INTENT(OUT) :: HollB
               REAL(sz), DIMENSION(nQuads) :: rmaxw

               CALL setStormName(name)
               CALL setAdvisory(adv)
               CALL setDate(yyyy, mm, dd, hh, hf)

               CALL setAmbientPressure(Pn)
               CALL setCentralPressure(Pc)

               CALL setCenterLatitude (cLat)
               CALL setCenterLongitude(cLon)

               CALL setMaxWinds  (Vmax)
               CALL setRadialWind(Vr)
               CALL setWindRadii (r)

               CALL setCoriolisForce(coriolis(cLat))
               CALL setShapeParameter(calcShape(Pn, Pc, Vmax))
!              CALL setFitMethod(POLY_FIT)
               CALL setFitMethod(SPLINE_FIT)
               CALL calcRmax()

              do i=1,nQuads
               rmaxw(i) = getRadiusMaxWinds(i+1)
               enddo
               HollB=getShapeParameter()


            END SUBROUTINE newVortex_rmax_b_OUT
                                                                
            !=================================================================
            ! Move a Vortex object by changing the location of its center and
            ! any of its intrinsic parameters (except its name).
            !
            ! On input:
            !    adv          NHC advisory or storm number
            !    yyyy         Year  of advisory  (4 digits)
            !    mm           Month of advisory  (2 digits)
            !    dd           Day   of advisory  (2 digits)
            !    hh           Hour  of advisory  (2 digits)
            !    hf           Forecast increment (2 digits)
            !    Pn           Ambient surface pressure (mb)
            !    Pc           Surface pressure at center of storm (mb)
            !    cLat         Latitude  of storm center (degrees north)
            !    cLon         Longitude of storm center (degrees east )
            !    Vmax         Max sustained wind velocity in storm (knots)
            !    Vr           Velocity at wind radii (knots)
            !    r            Wind radii - the distance winds of velocity Vr
            !                 extend outward from center of storm
            !                 (nautical miles)
            !
            ! On output:
            !    The vortex is moved to a new (cLat, cLon) location and all of
            !    its intrinsic parameters are recalculated/updated.
            !=================================================================
            SUBROUTINE moveVortex(adv, yyyy,mm,dd,hh,hf,
     &                            Pn,Pc, cLat,cLon, Vmax,Vr, r)
               INTEGER , INTENT(IN) :: adv
               INTEGER , INTENT(IN) :: yyyy
               INTEGER , INTENT(IN) :: mm
               INTEGER , INTENT(IN) :: dd
               INTEGER , INTENT(IN) :: hh
               INTEGER , INTENT(IN) :: hf
               REAL(sz), INTENT(IN) :: Pn
               REAL(sz), INTENT(IN) :: Pc
               REAL(sz), INTENT(IN) :: cLat
               REAL(sz), INTENT(IN) :: cLon
               REAL(sz), INTENT(IN) :: Vmax
               REAL(sz), INTENT(IN) :: Vr
               REAL(sz), INTENT(IN), DIMENSION(nQuads) :: r

               CALL newVortex(getStormName(), adv, yyyy,mm,dd,hh,hf,
     &                        Pn,Pc, cLat,cLon, Vmax,Vr, r)
            END SUBROUTINE moveVortex

            !==================================================================
            ! Move a vortex object by only changing the location of its center.
            ! The intensity of the storm does not change; the variation in the
            ! Coriolis force may slightly alter the shape of the vortex.
            !
            ! On input:
            !    adv          NHC advisory or storm number
            !    yyyy         Year  of advisory  (4 digits)
            !    mm           Month of advisory  (2 digits)
            !    dd           Day   of advisory  (2 digits)
            !    hh           Hour  of advisory  (2 digits)
            !    hf           Forecast increment (2 digits)
            !    cLat         Latitude  of storm center (degrees north)
            !    cLon         Longitude of storm center (degrees east )
            !
            ! On output:
            !    The vortex is translated to a new (cLat, cLon) location.
            !    Its advisory number and date & time are also updated.
            !==================================================================
            SUBROUTINE transVortex(adv, yyyy,mm,dd,hh,hf, cLat,cLon)
               INTEGER , INTENT(IN) :: adv
               INTEGER , INTENT(IN) :: yyyy
               INTEGER , INTENT(IN) :: mm
               INTEGER , INTENT(IN) :: dd
               INTEGER , INTENT(IN) :: hh
               INTEGER , INTENT(IN) :: hf
               REAL(sz), INTENT(IN) :: cLat
               REAL(sz), INTENT(IN) :: cLon

               CALL setAdvisory(adv)
               CALL setDate(yyyy, mm, dd, hh, hf)

               CALL setCenterLatitude (cLat)
               CALL setCenterLongitude(cLon)

               CALL setCoriolisForce(coriolis(cLat))
               CALL calcRmax()
            END SUBROUTINE transVortex

            !====================
            ! Print vortex object
            !====================
            SUBROUTINE printVortex()
               INTEGER :: n

               WRITE(*, '(64("-"))')
               WRITE(*, '(a, ", Advisory/Storm #", i2)')
     &              TRIM(getStormName()), getAdvisory()
               CALL printDate()
               WRITE(*, '(64("-"))')

               WRITE(*, '("   Ambient pressure = ", F7.2, " mb")')
     &              getAmbientPressure()

               WRITE(*, '("   Central pressure = ", F7.2, " mb")')
     &              getCentralPressure()

!               WRITE(*,'("   Eye location = (",F6.2," N, ",F7.2," E)")')
!     &              getCenterLatitude(), getCenterLongitude()

               WRITE(*,'("   Max winds at surface          = ", I3,
     &                   " knots")')
     &              NINT(getMaxWinds() * windReduction)

               WRITE(*,'("   Max winds at top of SFC layer = ", I3,
     &                   " knots")')
     &              NINT(getMaxWinds())

               WRITE(*,' ("   Wind radii for velocity = ",I3," knots:",
     &                    /, (6X,I3," nm"))')
     &              NINT(getRadialWind()), NINT(getWindRadii())

!               WRITE(*, '("   Coriolis force at ",F6.2," N = ", 1PG10.3,
!     &                    " 1/sec")')
!     &              getCenterLatitude(), getCoriolisForce()

               WRITE(*, '("   Exponential shape parameter B = ", F7.3)')
     &              getShapeParameter()

               WRITE(*, '("   Theta   Rmax")')
               WRITE(*, '("   -----   ----")')
               DO n = 1, nPoints
                   WRITE(111, '(4X, I3, 3X, F7.3, 3X, I4, 3X, F7.3)')
     &                NINT( getAzimuthalAngle(n)), getRadiusMaxWinds(n),
     &                 NINT(getMaxWinds()),getShapeParameter()
                   WRITE(*, '(4X, I3, 3X, F7.3, " nm")')
     &                 NINT(getAzimuthalAngle(n)), getRadiusMaxWinds(n)
               END DO

!               IF (fitMethod == POLY_FIT) THEN
!                  WRITE(*, '(3X, "Polynomial curve fit of Rmax vs. ",
!     &                           "azimuthal angle:")')
!                  WRITE(*, '(6X, "Rmax = ", 1PG10.3, /,
!     &                       6X, "     + ", 1PG10.3, " * theta"  , /,
!     &                       6X, "     + ", 1PG10.3, " * theta^2", /,
!     &                       6X, "     + ", 1PG10.3, " * theta^3", /,
!     &                       6X, "     + ", 1PG10.3, " * theta^4")')
!     &                 getPolyCoefficients()
!                  WRITE(*, '(6X, "Regression coefficient R^2 = ",
!     &                            1PG11.4)')
!     &                 getRegressionCoefficient()
!               ELSE
!                  WRITE(*, '(3X, "Spline fit of Rmax vs. ",
!     &                           "azimuthal angle")')
!               END IF

            END SUBROUTINE printVortex

            !===============================================================
            ! Calculate the Holland exponential shape parameter B:
            !
            !    B = (Vmax)^2 * (density of air) * e / (Pn - Pc)
            !
            ! On input:
            !    Pn           Ambient surface pressure (mb)
            !    Pc           Surface pressure at center of storm (mb)
            !    Vmax         Max sustained wind velocity in storm (knots)
            !
            ! On output:
            !    calcShape    Holland exponential shape parameter (no units)
            !===============================================================
            REAL(sz) FUNCTION calcShape(Pn, Pc, Vmax)
               REAL(sz), INTENT(IN) :: Pn
               REAL(sz), INTENT(IN) :: Pc
               REAL(sz), INTENT(IN) :: Vmax
!              PRINT *, "calcShape: Pn       = ", Pn
!              PRINT *, "calcShape: Pc       = ", Pc
!              PRINT *, "calcShape: Vmax     = ", Vmax
!              PRINT *, "calcShape: kt2ms    = ", kt2ms
!              PRINT *, "calcShape: RhoAir   = ", RhoAir
!              PRINT *, "calcShape: one      = ", one
!              PRINT *, "calcShape: EXP(one) = ", EXP(one)
!              PRINT *, "calcShape: mb2pa    = ", mb2pa
               calcShape = (Vmax*kt2ms)*(Vmax*kt2ms) * RhoAir * EXP(one)
     &                   / ((Pn - Pc) * mb2pa)
               calcShape = MAX( MIN(calcShape,twoPointFive), one)
            END FUNCTION calcShape

            !==================================================
            ! Calculate the Coriolis force at a given latitude.
            !
            !    Coriolis = 2 * omega * sin(latitude)
            !
            ! On input:
            !    lat          Latitude (degrees north)
            !
            ! On output:
            !    coriolis     Coriolis force (1/s)
            !==================================================
            REAL(sz) FUNCTION coriolis(lat)
               REAL(sz), INTENT(IN) :: lat
               coriolis = two * omega * SIN(deg2rad*lat)
!              coriolis = zero
            END FUNCTION coriolis

            !=================================================================
            ! Calculate surface pressure using Holland's hurricane wind model.
            !
            !    Ph = Pc + (Pn - Pc) * exp[ -(Rmax/r)**B ]
            !
            ! On input:
            !    Pc           Surface pressure at center of storm (mb)
            !    Pn           Ambient surface pressure (mb)
            !    Rmax         Radius of maximum winds (NM)
            !    r            Distance from center of storm (NM)
            !    B            Exponential shape parameter (no units)
            !
            ! On output:
            !    Ph           Surface pressure (mb)
            !=================================================================
            REAL(sz) FUNCTION Ph(Pc, Pn, Rmax, r, B)
               REAL(sz), INTENT(IN) :: Pc
               REAL(sz), INTENT(IN) :: Pn
               REAL(sz), INTENT(IN) :: Rmax
               REAL(sz), INTENT(IN) :: r
               REAL(sz), INTENT(IN) :: B
               Ph = Pc + (Pn - Pc) * EXP(-(Rmax/r)**B)
            END FUNCTION Ph
            !=================================================================
            ! Calculate surface pressure using an asymmetric hurricane wind
            ! model, where the radius of maximum winds (Rmax) is a polynomial
            ! function of the azimuthal angle (theta).
            !
            !    Pasym = Pc + (Pn - Pc) * exp[ -(Rmax(theta)/r)**B ]
            !
            ! On input:
            !    Pc           Surface pressure at center of storm (mb)
            !    Pn           Ambient surface pressure (mb)
            !    order        Order of polynomial for Rmax vs. theta curve fit
            !    a            Coefficients of the Rmax(theta) polynomial
            !    theta        Azimuthal angle (degrees)
            !    r            Distance from center of storm (NM)
            !    B            Exponential shape parameter (no units)
            !
            ! On output:
            !    Pasym        Surface pressure (mb)
            !=================================================================
            REAL(sz) FUNCTION Pasym(Pc, Pn, order, a, theta, r, B)
               REAL(sz), INTENT(IN) :: Pc
               REAL(sz), INTENT(IN) :: Pn
               INTEGER , INTENT(IN) :: order
               REAL(sz), INTENT(IN), DIMENSION(order+1) :: a
               REAL(sz), INTENT(IN) :: theta
               REAL(sz), INTENT(IN) :: r
               REAL(sz), INTENT(IN) :: B
               Pasym = Ph(Pc, Pn, Rmw(order,a,theta), r, B)
            END FUNCTION Pasym

            !==================================================================
            ! Calculate radial wind speed using Holland's hurricane wind model.
            !
            !    Vh = SQRT( B/RhoAir * (Rmax/r)**B * (Pn-Pc)*exp[-(Rmax/r)**B]
            !              + (r*coriolis/2)**2 ) - r*coriolis/2
            !
            !    Vh = SQRT( Vmax**2 * (Rmax/r)**B * exp[1 -(Rmax/r)**B]
            !              + (r*coriolis/2)**2 ) - r*coriolis/2
            ! On input:
            !    r            Distance from center of storm (NM)
            !    Rmax         Radius of maximum winds (NM)
            !    B            Exponential shape parameter (no units)
            !    Pn           Ambient surface pressure (mb)
            !    Pc           Surface pressure at center of storm (mb)
            !    cLat         Latitude of storm center (degrees north)
            !    Vmax         Maximum wind velocity
            !
            ! On output:
            !    Vh           Radial wind speed (knots)
            !
            ! Note:
            !    Computation is performed in mks, then converted back to knots.
            !==================================================================
!            REAL(sz) FUNCTION Vh(r, Rmax, B, Pn, Pc, cLat)
            REAL(sz) FUNCTION Vh(r, Rmax, B, Pn, Pc, cLat, Vmax)
               REAL(sz), INTENT(IN) :: r
               REAL(sz), INTENT(IN) :: Rmax
               REAL(sz), INTENT(IN) :: Vmax
               REAL(sz), INTENT(IN) :: B
               REAL(sz), INTENT(IN) :: Pn
               REAL(sz), INTENT(IN) :: Pc
               REAL(sz), INTENT(IN) :: cLat

               REAL(sz) :: deltaP    ! Pn - Pc (Pa)
               REAL(sz) :: Brho      ! B / density of air
               REAL(sz) :: RmaxrB    ! (Rmax/r)**B
               REAL(sz) :: rf2       ! r*f/2
               REAL(sz) :: Vmax_ms       ! get Vmax into m/s (it is passed in in knots)

               Vmax_ms=Vmax*kt2ms  ! get Vmax into m/s (it is passed in in knots)
               deltaP = (Pn - Pc) * mb2pa
               Brho = B / RhoAir
               RmaxrB = (Rmax/r)**B
!              rf2 = nm2m * r * coriolis(cLat) / two
               rf2 = nm2m * r * getCoriolisForce() / two
!               Vh = ms2kt * (SIGN(one, cLat) *
!     &                      SQRT( Brho * RmaxrB * deltaP * EXP(-RmaxrB)
!     &                            + rf2*rf2 ) - rf2)
               Vh = ms2kt * (SIGN(one, cLat) *
     &                      SQRT( Vmax_ms*Vmax_ms * RmaxrB  * EXP(1 - RmaxrB)
     &                            + rf2*rf2 ) - rf2)

!              PRINT *, "Vh:  deltaP = ", deltaP
!              PRINT *, "Vh:  Brho   = ", Brho
!              PRINT *, "Vh:  RmaxrB = ", RmaxrB
!              PRINT *, "Vh:  rf2    = ", rf2
!              PRINT *, "Vh:  Vh     = ", Vh
            END FUNCTION Vh


            !=================================================================
            ! Given a polynomial curve fit of Rmax as a function of azimuthal
            ! angle (theta), calculate the radius of maximum winds.
            !
            ! On input:
            !    order        Order of polynomial for Rmax vs. angle curve fit
            !    a            Coefficients of the Rmax(theta) polynomial
            !    angle        Azimuthal angle (degrees)
            !
            ! On output:
            !    Rmw          Radius of maximum winds (meters) from curve fit
            !                      I DO NOT BELIEVE IT IS IN METERS rjw
            !=================================================================
            REAL(sz) FUNCTION Rmw(order, a, angle)
               INTEGER , INTENT(IN) :: order
               REAL(sz), INTENT(IN), DIMENSION(order+1) :: a
               REAL(sz), INTENT(IN) :: angle

               IF (fitMethod == POLY_FIT) THEN
                  Rmw = poly(order+1, a, angle)
               ELSE
                  Rmw = spline(nPoints, REAL(theta),REAL(Rmax), d2ydx2,
     &                         REAL(angle))
               END IF

               IF (Rmw <= zero) Rmw = oneHundredth

            END FUNCTION Rmw

            !=================================================================
            ! Given a polynomial curve fit of Vmax as a function of azimuthal
            ! angle (theta), calculate the Maximum wind velocity.
            !
            ! On input:
            !    order        Order of polynomial for Vmax vs. angle curve fit
            !    a            Coefficients of the Vmax(theta) polynomial
            !    angle        Azimuthal angle (degrees)
            !
            ! On output:
            !    Vmaxw         maximum winds speed (meters/sec) from curve fit
            !=================================================================
            REAL(sz) FUNCTION Vmaxw(order, a, angle)
               INTEGER , INTENT(IN) :: order
               REAL(sz), INTENT(IN), DIMENSION(order+1) :: a
               REAL(sz), INTENT(IN) :: angle

               IF (fitMethod == POLY_FIT) THEN
                  Vmaxw = poly(order+1, a, angle)
               ELSE
            Vmaxw = spline(nPoints, REAL(theta),REAL(Vmax_wind), d2ydx2,
     &                         REAL(angle))
               END IF

               IF (Vmaxw <= zero) Vmaxw = oneHundredth

            END FUNCTION Vmaxw

            !=================================================================
            ! Calculate radial wind speed using an asymmetric hurricane wind
            ! model, where the radius of maximum winds (Rmax) is a polynomial
            ! function of the azimuthal angle (theta).
            !
            !    Vasym = SQRT( B/RhoAir * (Rmax(theta)/r)**B
            !                 * (Pn-Pc)*exp[-(Rmax(theta)/r)**B]
            !                 + (r*coriolis/2)**2 ) - r*coriolis/2
            !
            ! On input:
            !    r            Distance from center of storm (NM)
            !    order        Order of polynomial for Rmax vs. theta curve fit
            !    a            Coefficients of the Rmax(theta) polynomial
            !    theta        Azimuthal angle (degrees)
            !    B            Exponential shape parameter (no units)
            !    Pn           Ambient surface pressure (mb)
            !    Pc           Surface pressure at center of storm (mb)
            !    cLat         Latitude of storm center (degrees north)
            !    Vmax         Maximum velocity (knots)
            !
            ! On output:
            !    Vasym        Radial wind speed (knots)
            !
            ! Note:
            !    Computation is performed in mks, then converted back to knots
            !=================================================================
            REAL(sz) FUNCTION Vasym(r, order, a, theta, B, Pn, Pc, cLat,Vmax)
               REAL(sz), INTENT(IN) :: r
               INTEGER , INTENT(IN) :: order
               REAL(sz), INTENT(IN), DIMENSION(order+1) :: a
               REAL(sz), INTENT(IN) :: theta
               REAL(sz), INTENT(IN) :: B
               REAL(sz), INTENT(IN) :: Pn
               REAL(sz), INTENT(IN) :: Pc
               REAL(sz), INTENT(IN) :: cLat
               REAL(sz), INTENT(IN) :: Vmax

               Vasym = Vh(r, Rmw(order,a,theta), B, Pn, Pc, cLat,Vmax)

            END FUNCTION Vasym

            !========================================================
            ! Rotate a 2D vector (x,y) by an angle.
            !
            ! On input:
            !    x           x component of vector
            !    y           y component of vector
            !    angle       angle to rotate vector (degrees)
            !    whichWay    direction of rotation:
            !                   - = clockwise, + = counter-clockwise
            !
            ! On output:
            !    xr          x component of rotated vector
            !    yr          y component of rotated vector
            !========================================================
            SUBROUTINE rotate(x,y, angle, whichWay, xr,yr)
               REAL(sz), INTENT(IN)  :: x
               REAL(sz), INTENT(IN)  :: y
               REAL(sz), INTENT(IN)  :: angle
               REAL(sz), INTENT(IN)  :: whichWay

               REAL(sz), INTENT(OUT) :: xr
               REAL(sz), INTENT(OUT) :: yr

               REAL(sz) :: A, cosA, sinA

               A = SIGN(one, whichWay) * deg2rad * angle
               cosA = COS(A)
               sinA = SIN(A)

               xr = x * cosA - y * sinA
               yr = x * sinA + y * cosA
            END SUBROUTINE rotate

            !=======================================================
            ! Compute a wind angle to parameterize frictional inflow
            ! across isobars.
            !
            ! On input:
            !    r           distance from center of storm
            !    rmx         radius of maximum winds
            !
            ! On output:
            !    fang        frictional inflow angle (degrees)
            !=======================================================
            REAL(sz) FUNCTION fang(r, rmx)
               REAL(sz), INTENT(IN) :: r
               REAL(sz), INTENT(IN) :: rmx

               IF (zero <= r .AND. r < rmx) THEN
                  fang = ten * r/rmx
               ELSE IF (rmx <= r .AND. r < onePointTwo*rmx) THEN
                  fang = ten + seventyFive * (r/rmx - one)
               ELSE IF (r >= onePointTwo*rmx) THEN
                  fang = twentyFive
               ELSE
                  fang = zero
               END IF
            END FUNCTION fang

            !=================================================================
            ! Calculate (u,v) wind components and surface pressure from an
            ! asymmetric hurricane wind model.
            !
            ! On input:
            !    lat         Latitude  of nodal point (degrees north)
            !    lon         Longitude of nodal point (degrees east )
            !    uTrans      x component of translational velocity (m/s)
            !    vTrans      y component of translational velocity (m/s)
            !
            ! On output:
            !    u           x component of wind velocity at nodal point (m/s)
            !    v           y component of wind velocity at nodal point (m/s)
            !    p           Surface pressure at nodal point (Pa)
            !
            ! Internal parameters:
            !    dampRadii   How far out (# of Rmax radii) to begin damping
            !                out the translational velocity
            !
            ! Note:
            !    Subroutine directly accesses global class instance variables
            !=================================================================
            SUBROUTINE uvp(lat,lon, uTrans,vTrans, u,v, p)! ,rmx,speed)
! RJW
!           SUBROUTINE uvp(lat,lon, uTrans,vTrans, u,v, p,Vmax)
! MCF 5/26/09
!           SUBROUTINE uvp(lat,lon, uTrans,vTrans, u,v, p,Vmaxp)

               REAL(sz), INTENT(IN)  :: lat
               REAL(sz), INTENT(IN)  :: lon
               REAL(sz), INTENT(IN)  :: uTrans
               REAL(sz), INTENT(IN)  :: vTrans

! RJW
!              REAL(sz), INTENT(IN)  :: Vmax
! MCF 5/26/09
!              REAL(sz), INTENT(IN)  :: Vmaxp

               REAL(sz), INTENT(OUT) :: u
               REAL(sz), INTENT(OUT) :: v
               REAL(sz), INTENT(OUT) :: p

               REAL(sz) :: dx
               REAL(sz) :: dy
               REAL(sz) :: dist
               REAL(sz) :: rmx
               REAL(sz) :: angle
               REAL(sz) :: speed
               REAL(sz) :: uf
               REAL(sz) :: vf
               REAL(sz) :: damp
               !------------------------------------------------------
               ! Calculate distance and angle between eye of hurricane
               ! and input nodal point
               !------------------------------------------------------
               dx = deg2rad * Rearth * (lon - cLon) * COS(deg2rad*cLat)
               dy = deg2rad * Rearth * (lat - cLat)
               dist = SQRT(dx*dx + dy*dy)

               !----------------------------------------
               ! Handle special case at eye of hurricane
               !----------------------------------------
               ! in eye velocity is zero not translational velocity
               !----------------------------------------
               IF (dist < one) THEN
                  u = uTrans
                  v = vTrans
! RJW                  
!                 u = zero
!                 v = zero

                  p = Pc * mb2pa
                  RETURN
               END IF

               dist = m2nm * dist
               angle = threeSixty + rad2deg * ATAN2(dx,dy)
               IF (angle > threeSixty) angle = angle - threeSixty

               !-------------------------------------------------------
               ! Compute the damping factor for translational velocity:
               !    damp = 1 when r < dampRadii
               !    damp = 1 --> 0 as r = dampRadii --> dampRadii + 1
               !    damp = 0 when r > dampRadii + 1
               !-------------------------------------------------------
               rmx =  Rmw(order,a,angle)
!              damp = COS(-pi2*MAX(MIN((dist/rmx)-dampRadii,one),zero))
               damp = COS(-pi2*MAX(MIN(((dist/rmx)-dampRadii)
     &                                  /dampLength,one), zero))
!              damp = one

               !---------------------------------------------------
               ! Compute (u,v) wind velocity components from the
               ! asymmetric hurricane vortex.
               !
               ! Note: the vortex winds are valid at the top of the
               ! surface layer, so reduce the winds to the surface.
               ! Also convert the winds from max sustained 1-minute
               ! averages to 10-minute averages for the storm surge
               ! model.
               !---------------------------------------------------
        speed = kt2ms
     &       * Vasym(dist, order,a, angle, B, Pn,Pc, cLat,getMaxWinds())
     &       * windReduction * one2ten
 
! RJW
!              damp = ABS(speed/Vmax)                                   ! similar to the NWS=8 implementation
! MCF 5/26/09
!              damp = ABS(Vasym(dist, order,a, angle, B, Pn,Pc, cLat)/(Vmaxp))

               u = -speed * COS(deg2rad*angle) + uTrans * damp
               v =  speed * SIN(deg2rad*angle) + vTrans * damp

               !---------------------------------------------------------
               ! Alter wind direction by adding a frictional inflow angle
               !---------------------------------------------------------
               CALL rotate(u,v, fang(dist,rmx), cLat, uf,vf)
               u = uf
               v = vf

               !----------------------------------------------------------
               ! Compute surface pressure from asymmetric hurricane vortex
               !----------------------------------------------------------
               p = mb2pa * Pasym(Pc,Pn, order,a, angle, dist, B)

!              PRINT *, "uvp:  Rearth    = ", Rearth, " meters"
!              PRINT *, "uvp:  lat       = ", lat   , " degrees"
!              PRINT *, "uvp:  lon       = ", lon   , " degrees"
!              PRINT *, "uvp:  cLat      = ", cLat  , " degrees"
!              PRINT *, "uvp:  cLon      = ", cLon  , " degrees"
!              PRINT *, "uvp:  dLat      = ", (lat - cLat), " degrees"
!              PRINT *, "uvp:  dLon      = ", (lon - cLon), " degrees"
!              PRINT *, "uvp:  cos(cLat) = ", COS(deg2rad*cLat), " degrees"
!              PRINT *, "uvp:  dx        = ", dx,     " meters"
!              PRINT *, "uvp:  dy        = ", dy,     " meters"
!              PRINT *, "uvp:  dist      = ", dist,   " nm"
!              PRINT *, "uvp:  angle     = ", angle,  " degrees"
!              PRINT *, "uvp:  speed     = ", speed,  " m/s"
!              PRINT *, "uvp:  u         = ", u,      " m/s"
!              PRINT *, "uvp:  v         = ", v,      " m/s"
!              PRINT *, "uvp:  p         = ", p,      " Pascals"
!              PRINT *, "uvp:  uTrans    = ", uTrans,      " Pascals"
!              PRINT *, "uvp:  uTrans    = ", vTrans,      " Pascals"

            END SUBROUTINE uvp

            !============================================================
            ! Calculate the translational velocity of a moving hurricane.
            !
            ! On input:
            !    latOld      Previous latitude  of center (degrees north)
            !    lonOld      Previous longitude of center (degrees east )
            !    latNew      Current  latitude  of center (degrees north)
            !    lonNew      Current  longitude of center (degrees east )
            !    tOld        Previous time (seconds)
            !    tNew        Current  time (seconds)
            !
            ! On output:
            !    uTrans      x component of translational velocity (m/s)
            !    vTrans      y component of translational velocity (m/s)
            !============================================================
            SUBROUTINE uvtrans(latOld,lonOld, latNew,lonNew, tOld,tNew,
     &                         uTrans,vTrans)

               REAL(sz), INTENT(IN)  :: latOld
               REAL(sz), INTENT(IN)  :: lonOld

               REAL(sz), INTENT(IN)  :: latNew
               REAL(sz), INTENT(IN)  :: lonNew

               REAL(sz), INTENT(IN)  :: tOld
               REAL(sz), INTENT(IN)  :: tNew

               REAL(sz), INTENT(OUT) :: uTrans
               REAL(sz), INTENT(OUT) :: vTrans

               REAL(sz) :: dx
               REAL(sz) :: dy
               REAL(sz) :: dt

               dx = deg2rad * Rearth * (lonNew - lonOld)
     &            * COS(deg2rad*half * (latOld + latNew))

               dy = deg2rad * Rearth * (latNew - latOld)

               dt = tNew - tOld

               uTrans = dx/dt
               vTrans = dy/dt

!              PRINT *, "uvtrans:  Rearth    = ", Rearth, " meters"
!              PRINT *, "uvtrans:  latNew    = ", latNew, " degrees"
!              PRINT *, "uvtrans:  latOld    = ", latOld, " degrees"
!              PRINT *, "uvtrans:  dLat      = ", (latNew - latOld), " degrees"
!              PRINT *, "uvtrans:  lonNew    = ", lonNew, " degrees"
!              PRINT *, "uvtrans:  lonOld    = ", lonOld, " degrees"
!              PRINT *, "uvtrans:  dLon      = ", (lonNew - lonOld), " degrees"
!              PRINT *, "uvtrans:  tNew      = ", tNew,   " seconds"
!              PRINT *, "uvtrans:  tOld      = ", tOld,   " seconds"
!              PRINT *, "uvtrans:  dx        = ", dx,     " meters"
!              PRINT *, "uvtrans:  dy        = ", dy,     " meters"
!              PRINT *, "uvtrans:  dt        = ", dt,     " seconds"
!              PRINT *, "uvtrans:  uTrans    = ", uTrans, " m/s"
!              PRINT *, "uvtrans:  vTrans    = ", vTrans, " m/s"

            END SUBROUTINE uvtrans

            !=================================================================
            ! Transform (lat,lon) --> (x,y) coordinates.
            !
            ! On input:
            !    lat        Latitude  (degrees north)
            !    lon        Longitude (degrees east )
            !    lat0       Latitude  where projection is true (degrees north)
            !    lon0       Longitude where projection is true (degrees east )
            !
            ! On output:
            !    x          x (meters)
            !    y          y (meters)
            !=================================================================
            SUBROUTINE latlon2xy(lat,lon, lat0,lon0, x,y)
               REAL(sz), INTENT(IN)  :: lat ,lon
               REAL(sz), INTENT(IN)  :: lat0,lon0
               REAL(sz), INTENT(OUT) :: x,y

               x = deg2rad * Rearth * (lon - lon0) * COS(deg2rad*lat0)
               y = deg2rad * Rearth * lat
            END SUBROUTINE latlon2xy

            !=================================================================
            ! Transform (x,y) --> (lat,lon) coordinates.
            !
            ! On input:
            !    x          x (meters)
            !    y          y (meters)
            !    lat0       Latitude  where projection is true (degrees north)
            !    lon0       Longitude where projection is true (degrees east )
            !
            ! On output:
            !    lat        Latitude  (degrees north)
            !    lon        Longitude (degrees east )
            !=================================================================
            SUBROUTINE xy2latlon(x,y, lat0,lon0, lat,lon)
               REAL(sz), INTENT(IN)  :: x,y
               REAL(sz), INTENT(IN)  :: lat0,lon0
               REAL(sz), INTENT(OUT) :: lat ,lon

               lat = y / (deg2rad * Rearth)
               lon = lon0 + x / (deg2rad * Rearth * COS(deg2rad*lat0))
            END SUBROUTINE xy2latlon

            !=============================================================
            ! Evaluate a polynomial of degree n-1 at the value x.
            !
            !         n
            !    y = sum [ a(i) * x^(i-1) ]
            !        i=1
            !
            ! On input:
            !    n        number of terms (coefficients) in the polynomial
            !    a        array in which the coefficients are stored
            !    x        value at which to evaluate the polynomial
            !
            ! On output:
            !    poly     scalar result of evaluating the polynomial at x
            !=============================================================
            REAL(sz) FUNCTION poly(n, a, x)
               INTEGER , INTENT(IN) :: n
               REAL(sz), INTENT(IN), DIMENSION(n) :: a
               REAL(sz), INTENT(IN) :: x
               INTEGER :: i

               IF (n <= 0) THEN
                  poly = zero
               ELSE
                  poly = a(n)
                  DO i = n-1, 1, -1
                     poly = x * poly + a(i)
                  END DO
               END IF
            END FUNCTION poly

            !=========================================================
            ! A simple subroutine to fit a polynomial in one variable.
            ! Data must be stored in the form of n data pairs (x,y).
            !
            ! Polynomial to be fitted:
            !    y = a(0) + a(1)*x + a(2)*x^2 + ... + a(power)*x^power
            !
            ! On input:
            !    maxPoints    maximum number of (x,y) data points
            !    n            actual  number of (x,y) data points
            !    x,y          values of (x,y) data points
            !    maxPower     maximum   power of polynomial allowed
            !    power        requested power of polynomial
            !
            ! On output:
            !    a            array containing polynomial coefficients
            !=========================================================
            SUBROUTINE fitPoly(maxPoints,n, x,y, maxPower,power, a)

               !-----------------------------------------
               ! Import least-squares curve-fitting class
               !-----------------------------------------
               USE lsq

               !--------------------------------------------
               ! Force explicit declaration of all variables
               !--------------------------------------------
               IMPLICIT NONE

               INTEGER  , INTENT(IN)  :: maxPoints, n
               REAL (sz), INTENT(IN)  :: x(maxPoints), y(maxPoints)
               INTEGER  , INTENT(IN)  :: maxPower, power

               REAL (sz), INTENT(OUT) :: a(0:maxPower)

               REAL (sz) :: xrow(0:maxPower), sterr(0:maxPower),
     &                      covmat(231), var, rsq
               REAL (sz) :: wt = 1.d0
               INTEGER   :: anyErrors, i, j
               LOGICAL   :: fit_const = .TRUE., lindep(0:maxPower)

               !-----------------
               ! Input error trap
               !-----------------
               IF (n > maxPoints) THEN

                  WRITE(*, '(/, a)') 'FATAL ERROR in SUBROUTINE fitPoly'
                  WRITE(*, '(a, a)') 'Number of input data points ',
     &                               'exceeds maximum number allocated:'
                  WRITE(*, '(i6," > ",i6)') n, maxPoints
                  WRITE(*, '(a, i5)') "Increase 'maxPoints' to ", n
                  CALL killJob('Error trap -- array bounds overflow')

               ELSE IF (power > maxPower) THEN

                  WRITE(*, '(/, a)') 'FATAL ERROR in SUBROUTINE fitPoly'
                  WRITE(*, '(a, a)') 'Requested power of polynomial ',
     &                               'is too high:'
                  WRITE(*, '(i6," > ",i6)') power, maxPower
                  WRITE(*, '(a, i5)') "Increase 'maxPower' to ", power
                  CALL killJob('Error trap -- array bounds overflow')

               END IF

               !---------------------------
               ! Least-squares calculations
               !---------------------------
               CALL startup(power, fit_const)

               DO i = 1, n
                  xrow(0) = 1.d0
                  DO j = 1, power
                     xrow(j) = x(i) * xrow(j-1)
                  END DO
                  CALL includ(wt, xrow, y(i))
               END DO

               CALL sing(lindep, anyErrors)

               IF (anyErrors /= 0) THEN
                  DO i = 0, power
                     IF (lindep(i)) WRITE(*, '(a, i3)')
     &                             'Singularity detected for power: ', i
                  END DO
               END IF

               !-----------------------------------------------
               ! Calculate progressive residual sums of squares
               !-----------------------------------------------
               CALL ss()
               var = rss(power+1) / (n - power - 1)

               !------------------------------------------------
               ! Calculate least-squares regression coefficients
               !------------------------------------------------
               CALL regcf(a, power+1, anyErrors)

               !------------------------------------------------
               ! Calculate covariance matrix, and hence standard
               ! errors of coefficients.
               !------------------------------------------------
               CALL cov(power+1, var, covmat, 231, sterr, anyErrors)

               !-------------------------------
               ! Set regression coefficient R^2
               !-------------------------------
               rsq = (rss(1) - rss(power+1)) / rss(1)
               CALL setRegressionCoefficient(rsq)

               !--------------
               ! Print results
               !--------------
!              WRITE(*, '(36("- "))')
!              WRITE(*, '(21x, a)') 'SUBROUTINE fitPoly Diagnostics'
!              WRITE(*, '(21x, a)') '------------------------------'
!              WRITE(*, '(9x, a, a)') 'y = a(0) + a(1)*x + a(2)*x^2 + ... ', &
!     &                               '+ a(power)*x^power'
!              WRITE(*, '(a,a)') 'Least-squares coefficients and standard ', &
!     &                          'errors for polynomial curve fit:'
!              WRITE(*, *) 'Power Coefficient a         Std.error      ',    &
!                          'Resid.sum of sq.'
!              WRITE(*, *) '----- -------------         ---------      ',    &
!                          '----------------'
!              DO i = 0, power
!                 WRITE(*, '(i4, g20.12, "   ", g14.6, "   ", g14.6)')       &
!                      i, a(i), sterr(i), rss(i+1)
!              END DO
!              WRITE(*, '(a, g20.12)') 'Residual standard deviation = ',     &
!                                      SQRT(var)
!              WRITE(*, '(a, g20.12)') 'Regression coefficient R^2  = ',     &
!                                      rsq
!              WRITE(*, '(36("- "))')

            END SUBROUTINE fitPoly

            !===========================================================
            ! Fit a one-dimensional curve to a set of points
            ! (xgrid,ygrid) using cubic splines under tension.
            !
            ! On input:
            !    ndim      Number of grid points that define the curve.
            !    xgrid     A 1-D array containing the ndim x values that
            !              define the curve.
            !    ygrid     A 1-D array containing the ndim y values that
            !              define the curve.
            !    dydx--    Slopes of the given curve at its end points
            !              (optional).
            !    slopes    An integer flag indicating which and how the
            !              slope data should be used:
            !                 slopes     explanation
            !                 ------     -----------
            !                   0        use input slope data at both
            !                            end points.
            !                   1        use input slope data only at
            !                            lowest end point.
            !                   2        use input slope data only at
            !                            highest end point.
            !                   3        slope data is generated
            !                            internally.
            !    work      1-D array used for work space.
            !
            ! On output:
            !    d2ydx2    2nd derivatives at grid points along the 1-D
            !              curve.
            !===========================================================
            SUBROUTINE fitSpline(ndim, xgrid,ygrid, dydxlo,dydxhi,
     &                           slopes, work, d2ydx2)

               !-----------------------------------------
               ! Import cubic splines curve-fitting class
               !-----------------------------------------
               USE fitpack

               !--------------------------------------------
               ! Force explicit declaration of all variables
               !--------------------------------------------
               IMPLICIT NONE

               INTEGER, INTENT(IN)  :: ndim
               REAL   , INTENT(IN)  :: xgrid (ndim)
               REAL   , INTENT(IN)  :: ygrid (ndim)
               REAL   , INTENT(IN)  :: dydxlo
               REAL   , INTENT(IN)  :: dydxhi
               INTEGER, INTENT(IN)  :: slopes
               REAL   , INTENT(IN)  :: work  (ndim)
               REAL   , INTENT(OUT) :: d2ydx2(ndim)

               REAL, PARAMETER :: tension = 1.75

               INTEGER :: error_vortex

               !--------------------------------------------------------
               ! Generate derivatives at grid points along the 1-D curve
               !--------------------------------------------------------
               CALL curv1 (ndim, xgrid,ygrid, dydxlo,dydxhi,slopes,
     &                     d2ydx2, work, tension, error_vortex)

               !------------------
               ! Error diagnostics
               !------------------
               IF (error_vortex .NE. 0) THEN

                  IF (error_vortex .EQ. 1) THEN
                     CALL killJob("SUBROUTINE fitSpline: " //
     &                            "too few end points specified!")
                  ELSE IF (error_vortex .EQ. 2) THEN
                     CALL killJob("SUBROUTINE fitSpline: x values " //
     &                            "not monotonically increasing!")
                  END IF

               END IF
            END SUBROUTINE fitSpline

            !===========================================================
            ! Fit a one-dimensional curve to a set of points
            ! (xgrid,ygrid) using cubic splines under tension.
            !
            ! On input:
            !    ndim      Number of grid points that define the curve.
            !    xgrid     A 1-D array containing the ndim x values that
            !              define the curve.
            !    ygrid     A 1-D array containing the ndim y values that
            !              define the curve.
            !    d2ydx2    2nd derivatives at grid points along the 1-D
            !              curve.
            !    x         x value where data value is to be
            !              interpolated.
            !
            ! On output:
            !    spline    y value interpolated from cubic spline curve
            !===========================================================
            REAL(sz) FUNCTION spline(ndim, xgrid,ygrid, d2ydx2, x)

               !-----------------------------------------
               ! Import cubic splines curve-fitting class
               !-----------------------------------------
               USE fitpack

               !--------------------------------------------
               ! Force explicit declaration of all variables
               !--------------------------------------------
               IMPLICIT NONE

               INTEGER, INTENT(IN)  :: ndim
               REAL   , INTENT(IN)  :: xgrid (ndim)
               REAL   , INTENT(IN)  :: ygrid (ndim)
               REAL   , INTENT(IN)  :: d2ydx2(ndim)
               REAL   , INTENT(IN)  :: x

               REAL, PARAMETER :: tension = 1.75

               !------------------------------------------------
               ! Perform cubic spline interpolation to determine
               ! y value along the curve at data point x.
               !------------------------------------------------
               spline = REAL(curv2(x, ndim,xgrid,ygrid, d2ydx2,tension),
     &                       KIND = sz)
            END FUNCTION spline

            !==============================================================
            ! Use brute-force marching to find a root the interval [x1,x2].
            !
            ! On input:
            !    func        function f(x)=0 for which root is sough
            !    x1          left  side of interval
            !    x2          right side of interval
            !    dx          x increment for march
            !
            ! On output:
            !    a           left  side of interval that brackets the root
            !    b           right side of interval that brackets the root
            !    findRoot    root returned
            !==============================================================
            FUNCTION findRoot(func, x1,x2, dx, a,b)
               USE CONSTANTS, ONLY : error
               !--------------------------------------------------
               ! External function f(x)=0 for which root is sough
               !--------------------------------------------------
               REAL(sz), EXTERNAL :: func

               !------------------------
               ! Search interval [x1,x2]
               !------------------------
               REAL(sz), INTENT(IN) :: x1, x2

               !-------------------
               ! Marching incremen
               !-------------------
               REAL(sz), INTENT(IN) :: dx

               !---------------------------
               ! x values that bracket root
               !---------------------------
               REAL(sz), INTENT(OUT) :: a, b

               !---------------
               ! The root found
               !---------------
               REAL(sz) :: findRoot

               INTEGER , PARAMETER :: itmax = 1000000 ! Max # of iterations
               INTEGER :: iter                        ! iteration counter
               REAL(sz) :: fa,fb                      ! function values f(x)

               !---------------------------------
               ! Initialize left side of interval
               !---------------------------------
               a  = x1
               fa = func(a)

               !-----------------------------------------
               ! March along interval until root is found
               ! or solution diverges.
               !-----------------------------------------
               findRoot = a
               DO iter = 1,itmax
                  b = x1 + iter * dx
                  fb = func(b)

                  !---------------
                  ! Check progress
                  !---------------
!                 PRINT *, "FUNCTION findRoot:  At iteration # ", iter
!                 PRINT *, "   [a,b]     = [", a, ",", b, "]"
!                 PRINT *, "   f(a),f(b) = ", fa, ", ", fb

                  IF ((fa*fb < zero) .OR. (ABS(fb) > ABS(fa))) THEN
!                    PRINT *, "FUNCTION findRoot:"
!                    PRINT *, "   Root found at iteration # ", iter
!                    PRINT *, "   [a,b]     = [", a, ",", b, "]"
!                    PRINT *, "   f(a),f(b) = ", fa, ", ", fb

                     !------------
                     ! Assign root
                     !------------
                     IF (ABS(fb) > ABS(fa)) THEN
                        findRoot = a
!                       PRINT *, "   root = ", a
                     ELSE
                        findRoot = b
!                       PRINT *, "   root = ", b
                     END IF

                     EXIT
                  END IF

                  !-----------------------------------------------
                  ! Move right search interval values to left side
                  ! for next iteration.
                  !-----------------------------------------------
                  a  = b
                  fa = fb
               END DO

               IF (iter >= itmax) THEN
!                 PRINT *, "FUNCTION findRoot: exceeded max # of iterations"
                  findRoot = error
               END IF
            END FUNCTION findRoot

            !=============================================================
            ! Using Brent's method, find the root of a FUNCTION func known
            ! to lie between x1 and x2.
            !
            ! On input:
            !    func     function f(x)=0 for which root is sough
            !    x1       left  side of interval that brackets the root
            !    x2       right side of interval that brackets the root
            !    tol      user-specified tolerance (accuracy) of the root
            !
            ! On output:
            !    zbrent   root returned will be refined until its accuracy
            !             is tol.
            !
            ! Internal parameters:
            !    itmax    maximum allowed # of iterations
            !    eps      small number at machine floating-point precision
            !
            ! Reference:
            !    http://en.wikipedia.org/wiki/Brent's_method
            !=============================================================
            FUNCTION zbrent(func, x1,x2, tol)
               USE CONSTANTS, ONLY : ERROR
               !--------------------------------------------------
               ! External function f(x)=0 for which root is sough
               !--------------------------------------------------
               REAL(sz), EXTERNAL :: func

!              INTERFACE
!                 FUNCTION func(x)
!                    IMPLICIT NONE
!                    REAL(sz), INTENT(IN) :: x
!                    REAL(sz) :: func
!                 END FUNCTION func
!              END INTERFACE

               REAL(sz), INTENT(IN) :: x1,x2,tol
               REAL(sz) :: zbrent

               INTEGER, PARAMETER :: itmax = 100      ! Max # of iterations

               !--------------------------------------------------------------
               ! The routine zbrent works best when eps is exactly the machine
               ! precision. The Fortran 90 intrinsic function EPSILON allows
               ! us to code this in a portable fashion.
               !--------------------------------------------------------------
               REAL(sz), PARAMETER :: eps = EPSILON(x1)  ! a very small number

               INTEGER :: iter            ! iteration counter
               REAL(sz) :: a,b,c, d,e     ! x values in bounding interval
               REAL(sz) :: fa,fb,fc       ! their function values f(x)
               REAL(sz) :: p,q,r,s        ! inverse quadratic factors
               REAL(sz) :: xm             ! interval midpoin
               REAL(sz) :: tol1           ! computed tolerance

               !--------------------------------------------------------------
               ! Initialize Brent's methods with two points a and b, such tha
               ! f(a) and f(b) have opposite signs. The intermediate value
               ! theorem guarantees a solution between a and b.
               !--------------------------------------------------------------
               a = x1
               b = x2

               fa = func(a)
               fb = func(b)

               !--------------------------------------------------
               ! Oops - function values have same signs! Bail out.
               !--------------------------------------------------
               IF ( (fa > zero .AND. fb > zero) .OR.
     &              (fa < zero .AND. fb < zero) ) THEN
!                 PRINT *, "zbrent: fa and fb have same signs:"
!                 PRINT *, "zbrent: fa = ", fa
!                 PRINT *, "zbrent: fb = ", fb
                  zbrent = error
!                 CALL killJob("FUNCTION zbrent: root must be bracketed")
                  RETURN
               END IF

               !----------------------------------
               ! Save previous value of b and f(b)
               !----------------------------------
               c  = b
               fc = fb

               !--------------------------------------------------------------
               ! Iterate until we find an approximation to the root within the
               ! user-specified tolerance.
               !--------------------------------------------------------------
               DO iter = 1,itmax

                  IF ( (fb > zero .AND. fc > zero) .OR.
     &                 (fb < zero .AND. fc < zero) ) THEN
                     !----------------------------------------------
                     ! Rename a, b, c and adjust bounding interval d
                     !----------------------------------------------
                     c = a
                     fc = fa
                     d = b - a
                     e = d
                  END IF

                  IF (ABS(fc) < ABS(fb)) THEN
                     a = b
                     b = c
                     c = a

                     fa = fb
                     fb = fc
                     fc = fa
                  END IF

                  !------------------
                  ! Convergence check
                  !------------------
                  tol1 = two * eps * ABS(b) + half * tol

                  !-----------------------------
                  ! Compute midpoint of interval
                  !-----------------------------
                  xm = half * (c-b)

                  !---------------
                  ! Check progress
                  !---------------
!                 PRINT *, "At iteration # ", iter
!                 PRINT *, "   [a,b] = [", a, ",", b, "]"
!                 PRINT *, "   root = ", b

                  IF (ABS(xm) <= tol1 .OR. fb == zero) THEN
                     !-----------------------
                     ! Success -- we're done!
                     !-----------------------
                     zbrent = b
                     RETURN
                  END IF

                  IF (ABS(e) >= tol1 .AND. ABS(fa) > ABS(fb)) THEN
                     !----------------------------------------
                     ! Attempt inverse quadratic interpolation
                     !----------------------------------------
                     s = fb/fa

                     IF (a == c) THEN
!                  PRINT *, "   a = c:  secant method"
                        p = two * xm * s
                        q = one - s
                     ELSE
!                  PRINT *, "   inverse quadratic method"
                        q = fa/fc
                        r = fb/fc
                        p = s * (two*xm*q*(q-r) - (b-a)*(r-one))
                        q = (q-one) * (r-one) * (s-one)
                     END IF

                     !------------------------
                     ! Check whether in bounds
                     !------------------------
                     IF (p > zero) q = -q

                     p = ABS(p)

                     IF (two*p < MIN(three*xm*q-ABS(tol1*q), ABS(e*q)))
     &               THEN
                        !----------------------
                        ! Accept interpolation.
                        !----------------------
!                       PRINT *, "   using interpolation"
                        e = d
                        d = p/q
                     ELSE
                        !-------------------------------------
                        ! Interpolation failed; use bisection.
                        !-------------------------------------
!                       PRINT *, "   using bisection"
                        d = xm
                        e = d
                     END IF

                  ELSE
                     !---------------------------------------------
                     ! Bounds decreasing too slowly; use bisection.
                     !---------------------------------------------
                     d = xm
                     e = d
                  END IF

                  !--------------------------
                  ! Move last best guess to a
                  !--------------------------
                  a = b
                  fa = fb

                  !------------------------
                  ! Evaluate new trial root
                  !------------------------
                  b = b + MERGE(d, SIGN(tol1,xm), ABS(d)>tol1)
                  fb = func(b)
               END DO

               PRINT *, "FUNCTION zbrent: exceeded max # of iterations"
               zbrent = b
            END FUNCTION zbrent

            !====================================================
            ! External function f(x) = 0 for which a root is
            ! sought using Brent's root-finding method.
            !
            ! On input:
            !    x       iterative values which converge to root
            !
            ! On output:
            !    func    f(x)
            !
            ! Internal parameters:
            !    vortex instance variables via accessor functions
            !====================================================
            FUNCTION func(x)

               REAL(sz), INTENT(IN) :: x    ! x = Rmax
               REAL(sz) :: func

               REAL(sz) :: VeeH

               !-------------------------
               ! func(x = Rmax) = Vh - Vr
               !-------------------------
!              PRINT *, "func:  Quadrant    = ", getQuadrant()
!              PRINT *, "func:  Wind radius = ", getWindRadius(getQuadrant())
!              PRINT *, "func:  x           = ", x
!              PRINT *, "func:  B           = ", getShapeParameter()
!              PRINT *, "func:  Pn          = ", getAmbientPressure()
!              PRINT *, "func:  Pc          = ", getCentralPressure()
!              PRINT *, "func:  cLat        = ", getCenterLatitude()
!              VeeH = Vh (getWindRadius(getQuadrant()),
!                         x,
!    &                    getShapeParameter (),
!    &                    getAmbientPressure(),
!    &                    getCentralPressure(),
!    &                    getCenterLatitude (),
!     &                   getMaxWinds() )
!              PRINT *, "----------------------------------------"
!              PRINT *, "func:  x           = ", x
!              PRINT *, "func   Vh          = ", Veeh
!              PRINT *, "func:  Vr          = ", getRadialWind()
!              PRINT *, "func:  Vmax        = ", getMaxWinds()
!              PRINT *, "func:  func        = ", Veeh - getRadialWind()

               func = ABS(Vh (getWindRadius(getQuadrant()),
     &                        x,
     &                        getShapeParameter (),
     &                        getAmbientPressure(),
     &                        getCentralPressure(),
     &                        getCenterLatitude (),
     &                        getMaxWinds() ) )
     &              - getRadialWind()

            END FUNCTION func

            !===============================================================
            ! Calculate the radius of maximum winds for all storm quadrants.
            !
            ! On input:
            !    none
            !
            ! On output:
            !    Rmax    radius of maximum winds (nm) in all quadrants, plus
            !            2 extra values to tie down circular periodicity
            !===============================================================
            SUBROUTINE calcRmax()
!              REAL(sz), EXTERNAL  :: func
               REAL(sz)            :: root        ! Radius of maximum winds
               REAL(sz), PARAMETER :: innerRadius = one
               REAL(sz), PARAMETER :: outerRadius = oneHundred
               REAL(sz), PARAMETER :: accuracy    = .0001_dp
               REAL(sz), PARAMETER :: zoom        = oneHundredth
               INTEGER , PARAMETER :: itmax = 3
               REAL(sz)            :: r1,r2, r3,r4, dr
               INTEGER             :: n, iter

              !-----------------------------
              ! Loop over quadrants of storm
              !-----------------------------
               DO n = 1, nQuads
                  !---------------------------------------------
                  ! Find Rmax using Brent root-finding algorithm
                  ! where func(Rmax) = Vh-Vr = 0
                  !---------------------------------------------
                  CALL setQuadrant(n)
                  root = -one
!                 PRINT *, "========================================"
!                 PRINT *, "calcRmax:  Quadrant    = ", getQuadrant()
!                 PRINT *, "calcRmax:  Wind radius = ",
!    &                     getWindRadius(getQuadrant())
                  root = zbrent(func, innerRadius,outerRadius, accuracy)

                  !------------------------------------------
                  ! Sophisticated root-finding method failed,
                  ! use ugly, brute-force method.
                  !------------------------------------------
                  IF (root < zero) THEN
!                    PRINT *, "-------------------------------------------"
!                    PRINT *, "calcRmax: WARNING! zbrent method failed, so"
!                    PRINT *, "          using ugly, brute-force method to"
!                    PRINT *, "          find root Rmax in quadrant ", getQuadrant()
                     r1 = innerRadius
                     r2 = outerRadius
                     dr = one
                     DO iter = 1, itmax
                        root = findRoot(func, r1,r2, dr, r3,r4)
!                       PRINT *, "          ", iter,": ", r1,", ",r2,", ",dr,", ",root
                        r1 = r3
                        r2 = r4
                        dr = dr * zoom
                     END DO
!                    PRINT *, "-------------------------------------------"
                  END IF

                  !-------------------------------------------
                  ! All root-finding methods failed, bail out.
                  !-------------------------------------------
                  IF (root < zero) THEN
                     CALL killJob("SUBROUTINE calcRmax: " //
     &                            "Rmax calculation failed!")
                  END IF

                  CALL setRadiusMaxWinds(n+1, root)


               END DO
               !--------------------------------------------------------
               ! Generate 2 additional (theta,Rmax) points for curve-fit
               !--------------------------------------------------------
               root = half *
     &              (getRadiusMaxWinds(2) + getRadiusMaxWinds(nQuads+1))
               CALL setRadiusMaxWinds(         1, root)
               CALL setRadiusMaxWinds(nQuads + 2, root)
c                DO n = 1, nQuads
c          Write(*,*) "Rmax in Quatrandt", n, " = ", getRadiusMaxWinds(n)
c                enddo

               IF (fitMethod == POLY_FIT) THEN
                  !--------------------------------------------------
                  ! Fit a polynomial to the (theta,Rmax) data points.
                  ! Polynomial coefficients are stored in array 'a'.
                  !--------------------------------------------------
                  CALL fitPoly(nPoints,nPoints,theta,Rmax,order,order,a)
               ELSE
                  !----------------------------------------------------
                  ! Fit a cubic spline to the (theta,Rmax) data points.
                  ! 2nd derivatives are stored in array 'd2ydx2'.
                  !----------------------------------------------------
                  CALL fitSpline(nPoints, REAL(theta),REAL(Rmax),
     &                            0.,0., 3, work, d2ydx2)
               END IF

            END SUBROUTINE calcRmax


            !===============================================================
            ! RJW 07 - 2009 
            ! Calculate the coefficients that fit the given 
            ! radius of maximum winds for all storm quadrants.
            !
            ! On input:
            !   Rmax in all 4 quadrants plus 2 extra values to tie down circular periodicity
            !
            ! On output:
            !    Rmax    radius of maximum winds (nm) in all quadrants, plus
            !            2 extra values to tie down circular periodicity
            !===============================================================
            SUBROUTINE calcRmax_givenRmax(rmaxw)
!              REAL(sz), EXTERNAL  :: func
               REAL(sz)            :: root        ! Radius of maximum winds
C               REAL(sz), PARAMETER :: innerRadius = one
C               REAL(sz), PARAMETER :: outerRadius = oneHundred
C               REAL(sz), PARAMETER :: accuracy    = .0001_dp
C               REAL(sz), PARAMETER :: zoom        = oneHundredth
C               INTEGER , PARAMETER :: itmax = 3
C               REAL(sz)            :: r1,r2, r3,r4, dr
               INTEGER             :: n, iter,i
               REAL(sz), INTENT(IN),DIMENSION(nQuads) :: rmaxw

              !-----------------------------
              ! Loop over quadrants of storm
              !-----------------------------
               DO n = 1, nQuads
                  !---------------------------------------------
                  ! Find Rmax using Brent root-finding algorithm
                  ! where func(Rmax) = Vh-Vr = 0
                  !---------------------------------------------
                  CALL setQuadrant(n)

                  CALL setRadiusMaxWinds(n+1, rmaxw(n))

               END DO
               !--------------------------------------------------------
               ! Generate 2 additional (theta,Rmax) points for curve-fit
               !--------------------------------------------------------
               root = half *
     &              (getRadiusMaxWinds(2) + getRadiusMaxWinds(nQuads+1))
               CALL setRadiusMaxWinds(         1, root)
               CALL setRadiusMaxWinds(nQuads + 2, root)
c                DO n = 1, nQuads+2
c             Write(*,*) "Rmax in Quatrandt", n, " = ", getRadiusMaxWinds(n)
c                enddo

               IF (fitMethod == POLY_FIT) THEN
                  !--------------------------------------------------
                  ! Fit a polynomial to the (theta,Rmax) data points.
                  ! Polynomial coefficients are stored in array 'a'.
                  !--------------------------------------------------
                  CALL fitPoly(nPoints,nPoints,theta,Rmax,order,order,a)
               ELSE
                  !----------------------------------------------------
                  ! Fit a cubic spline to the (theta,Rmax) data points.
                  ! 2nd derivatives are stored in array 'd2ydx2'.
                  !----------------------------------------------------
C                  write(16,85) (Rmax(i),i=1,nPoints) !jgfdebug
C 85               format('Rmax=',6(E14.5,','))
                  CALL fitSpline(nPoints, REAL(theta),REAL(Rmax),
     &                            0.,0., 3, work, d2ydx2)
               END IF

            END SUBROUTINE calcRmax_givenRmax

            !==================================================
            ! Print an error message, then terminate execution.
            !==================================================
            SUBROUTINE killJob(message)
             CHARACTER(LEN = *), INTENT(IN) :: message
              WRITE (*,'(a)') message
               STOP "Program terminated by SUBROUTINE killJob"
            END SUBROUTINE killJob

            !===========================
            ! Setter (mutator) functions
            !===========================
            SUBROUTINE setStormName(name)
               CHARACTER(LEN=*), INTENT(IN) :: name
               stormName = name
            END SUBROUTINE setStormName

            SUBROUTINE setAdvisory(adv)
               INTEGER, INTENT(IN) :: adv
               advisory = adv
            END SUBROUTINE setAdvisory

            SUBROUTINE setDate(yyyy, mm, dd, hh, hf)
               INTEGER, INTENT(IN) :: yyyy
               INTEGER, INTENT(IN) :: mm
               INTEGER, INTENT(IN) :: dd
               INTEGER, INTENT(IN) :: hh
               INTEGER, INTENT(IN) :: hf
               year     = yyyy
               month    = mm
               day      = dd
               hour     = hh
               hourFcst = hf
            END SUBROUTINE setDate

            SUBROUTINE setAmbientPressure(p)
               REAL(sz), INTENT(IN) :: p
               Pn = p
            END SUBROUTINE setAmbientPressure

            SUBROUTINE setCentralPressure(p)
               REAL(sz), INTENT(IN) :: p
               Pc = p
            END SUBROUTINE setCentralPressure

            SUBROUTINE setCenterLatitude(lat)
               REAL(sz), INTENT(IN) :: lat
               cLat = lat
            END SUBROUTINE setCenterLatitude

            SUBROUTINE setCenterLongitude(lon)
               REAL(sz), INTENT(IN) :: lon
               cLon = lon
            END SUBROUTINE setCenterLongitude

            SUBROUTINE setMaxWinds(v)
               REAL(sz), INTENT(IN) :: v
               Vmax = v
C              write(*,*) 'Vmax in vortex = ', Vmax
            END SUBROUTINE setMaxWinds

            SUBROUTINE setShapeParameter(s)
               REAL(sz), INTENT(IN) :: s
               B = s
            END SUBROUTINE setShapeParameter

            SUBROUTINE setCoriolisForce(f)
               REAL(sz), INTENT(IN) :: f
               corio = f
            END SUBROUTINE setCoriolisForce

            SUBROUTINE setRadialWind(v)
               REAL(sz), INTENT(IN) :: v
               Vr = v
            END SUBROUTINE setRadialWind

            SUBROUTINE setWindRadii(r)
               REAL(sz), DIMENSION(nQuads), INTENT(IN) :: r
               radius = r
            END SUBROUTINE setWindRadii

            SUBROUTINE setWindRadius(n, r)
               INTEGER , INTENT(IN) :: n
               REAL(sz), INTENT(IN) :: r
               radius(n) = r
            END SUBROUTINE setWindRadius

            SUBROUTINE setPolyCoefficients(c)
               REAL(sz), DIMENSION(order+1), INTENT(IN) :: c
               a = c
            END SUBROUTINE setPolyCoefficients

            SUBROUTINE setPolyCoefficient(n, c)
               INTEGER , INTENT(IN) :: n
               REAL(sz), INTENT(IN) :: c
               a(n) = c
            END SUBROUTINE setPolyCoefficient

            SUBROUTINE setRadiiMaxWinds(rmx)
               REAL(sz), DIMENSION(nPoints), INTENT(IN) :: rmx
               Rmax = rmx
            END SUBROUTINE setRadiiMaxWinds

            SUBROUTINE setRadiusMaxWinds(n, rmx)
               INTEGER , INTENT(IN) :: n
               REAL(sz), INTENT(IN) :: rmx
               Rmax(n) = rmx
            END SUBROUTINE setRadiusMaxWinds

            SUBROUTINE setRegressionCoefficient(rsq)
               REAL(sz), INTENT(IN) :: rsq
               rSquare = rsq
            END SUBROUTINE setRegressionCoefficient

            SUBROUTINE setQuadrant(q)
               INTEGER, INTENT(IN) :: q
               quad = q
            END SUBROUTINE setQuadrant

            SUBROUTINE setFitMethod(fm)
               INTEGER, INTENT(IN) :: fm
               fitMethod = fm
            END SUBROUTINE setFitMethod

            !============================
            ! Getter (accessor) functions
            !============================
            CHARACTER(LEN=64) FUNCTION getStormName()
               getStormName = stormName
            END FUNCTION getStormName

            INTEGER FUNCTION getAdvisory()
               getAdvisory = advisory
            END FUNCTION getAdvisory

            SUBROUTINE printDate()
               IF (hourFcst > 0) THEN
                  IF (hour > 100) THEN
                     WRITE (*, '(i2, "-hour forecast valid at ", i2,
     &                           ":", i2, " UTC ", a, i3, ", ", i4)')
     &                     hourFcst, hour/100, (hour-100*(hour/100)),
     &                     TRIM(months(month)), day, year
                  ELSE
                     WRITE (*, '(i2, "-hour forecast valid at ", i2,
     &                           ":00 UTC ", a, i3, ", ", i4)')
     &                     hourFcst,hour, TRIM(months(month)), day, year
                  END IF
               ELSE
                  IF (hour > 100) THEN
                     WRITE (*, '("NHC advisory valid at ", i2,
     &                           ":", i2, " UTC ", a, i3, ", ", i4)')
     &                     hour/100, (hour-100*(hour/100)),
     &                     TRIM(months(month)), day, year
                  ELSE
                     WRITE (*, '("NHC advisory valid at ", i2,
     &                           ":00 UTC ", a, i3, ", ", i4)')
     &                     hour, TRIM(months(month)), day, year
                  END IF
               END IF
            END SUBROUTINE printDate

            INTEGER FUNCTION getYear()
               getYear = year
            END FUNCTION getYear

            INTEGER FUNCTION getMonth()
               getMonth = month
            END FUNCTION getMonth

            INTEGER FUNCTION getDay()
               getDay = day
            END FUNCTION getDay

            INTEGER FUNCTION getHour()
               getHour = hour
            END FUNCTION getHour

            INTEGER FUNCTION getForecastIncrement()
               getForecastIncrement = hourFcst
            END FUNCTION getForecastIncrement

            REAL(sz) FUNCTION getAmbientPressure()
               getAmbientPressure = Pn
            END FUNCTION getAmbientPressure

            REAL(sz) FUNCTION getCentralPressure()
               getCentralPressure = Pc
            END FUNCTION getCentralPressure

            REAL(sz) FUNCTION getCenterLatitude()
               getCenterLatitude = cLat
            END FUNCTION getCenterLatitude

            REAL(sz) FUNCTION getCenterLongitude()
               getCenterLongitude = cLon
            END FUNCTION getCenterLongitude

            REAL(sz) FUNCTION getMaxWinds()
               getMaxWinds = Vmax
            END FUNCTION getMaxWinds

            REAL(sz) FUNCTION getShapeParameter()
               getShapeParameter = B
            END FUNCTION getShapeParameter

            REAL(sz) FUNCTION getCoriolisForce()
               getCoriolisForce = corio
            END FUNCTION getCoriolisForce

            REAL(sz) FUNCTION getRadialWind()
               getRadialWind = Vr
            END FUNCTION getRadialWind

            FUNCTION getWindRadii() RESULT(radii)
               REAL(sz), DIMENSION(nQuads) :: radii
               radii = radius
            END FUNCTION getWindRadii

            REAL(sz) FUNCTION getWindRadius(n)
               INTEGER , INTENT(IN) :: n
               getWindRadius = radius(n)
            END FUNCTION getWindRadius

            INTEGER FUNCTION getPolyOrder()
               getPolyOrder = order
            END FUNCTION getPolyOrder

            FUNCTION getPolyCoefficients() RESULT(coefficients)
               REAL(sz), DIMENSION(order+1) :: coefficients
               coefficients = a
            END FUNCTION getPolyCoefficients

            REAL(sz) FUNCTION getPolyCoefficient(n)
               INTEGER , INTENT(IN) :: n
               getPolyCoefficient = a(n)
            END FUNCTION getPolyCoefficient

            FUNCTION getAzimuthalAngles() RESULT(angles)
               REAL(sz), DIMENSION(nPoints) :: angles
               angles = theta
            END FUNCTION getAzimuthalAngles

            REAL(sz) FUNCTION getAzimuthalAngle(n)
               INTEGER , INTENT(IN) :: n
               getAzimuthalAngle = theta(n)
            END FUNCTION getAzimuthalAngle

            FUNCTION getRadiiMaxWinds() RESULT(maxWinds)
               REAL(sz), DIMENSION(nPoints) :: maxWinds
               maxWinds = Rmax
            END FUNCTION getRadiiMaxWinds

            REAL(sz) FUNCTION getRadiusMaxWinds(n)
               INTEGER , INTENT(IN) :: n
               getRadiusMaxWinds = Rmax(n)
            END FUNCTION getRadiusMaxWinds

            REAL(sz) FUNCTION getRegressionCoefficient()
               getRegressionCoefficient = rSquare
            END FUNCTION getRegressionCoefficient

            INTEGER FUNCTION getQuadrant()
               getQuadrant = quad
            END FUNCTION getQuadrant

            INTEGER FUNCTION getFitMethod()
               getFitMethod = fitMethod
            END FUNCTION getFitMethod

      END MODULE vortex
