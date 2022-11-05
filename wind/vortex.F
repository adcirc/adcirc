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
      !    05/12/09    Robert Weaver, UNC-IMS
      !                              Modified damping of translational
      !                              velocity based on ratio of V/Vmax
      !    05/19/09    Cristina Forbes, UNC-IMS
      !                              Implemented old bug fix: units conversion
      !    05/22/09    Cristina Forbes, UNC-IMS
      !                              Reverted damping of translational
      !                              velocity back to the original formulation
      !                              developed by Craig Mattocks
      !    05/26/09    Cristina Forbes, UNC-IMS
      !                              Fixed implementation of Rick Luettich
      !                              V/Vmax tapering formulation (as in NWS=8)
      !                              for future experimentation - not activated
      !    06/2009      Robert Weaver, UNC-IMS
      !                              changed the gradient wind formula to
      !                              use Vmax instead of pressure difference
      !                              to enable easier manipulation of storm
      !                              characteristics
      !    07/2009      Robert Weaver, UNC-IMS
      !                          using modules from this code,
      !                          wrote a preprocessor to compute Rmax and B
      !                          prior to running adcirc, and write out an
      !                          input file like the ATCF file with these
      !                          variables to make it easier to manipulate
      !                          the storm parameters.  Changed this code
      !                          to read in that input file and use the values
      !                          found there
      !    02/2013      Jie Gao, UNC-IMS
      !                          wrote subroutines and functions for
      !                          the Generalized Asymmetric Holland Model (GAHM)
      !
      !=================================================================
      MODULE vortex
      USE SIZES, ONLY : sz
      USE GLOBAL, ONLY : sphericalDistance, deg2rad, rad2deg,
     &                   omega, mb2pa, rhoAir, kt2ms, nm2m, m2nm,
     &                   ms2kt, Rearth, windReduction, one2ten
      IMPLICIT NONE
      SAVE

      PUBLIC :: newVortex, calcRmaxes, fitRmaxes,
     &          Rmw, uvp, uvtrans, fang,
     &          latlon2xy, xy2latlon,
     &          getShapeParameter   ,
     &          setShapeParameter   ,
     &          getLatestRmax,
     &          getLatestAngle      ,
     &          getUseQuadrantVr    , getRmaxes,
     &          setUseQuadrantVr    , setRmaxes,
     &          setIsotachWindSpeed ,
     &          getShapeParameters  , getPhiFactors,
     &          setIsotachWindSpeeds, setIsotachRadii,
     &          setVortex, spInterp,interpR,
     &          setUseVmaxesBL,
     &          getVmaxesBL, setVmaxesBL,
     &          fitRmaxes4,calcRmaxesFull,
     &          uvpr,newVortexFull,
     &          Rmaxes4, QuadFlag4, QuadIr4, Bs4, VmBL4 !jgfdebug

      PRIVATE

      INTEGER, PARAMETER :: nQuads  = 4       ! Number of quadrants for
                                                 ! which wind radii are
                                                 ! provided
      INTEGER, PARAMETER :: nPoints = nQuads+2   ! Number of (theta, Rmax)
                                                 ! points for curve fit
      REAL(sz), DIMENSION(nPoints) :: Rmaxes     ! Radius of maximum winds
      REAL(sz), DIMENSION(nPoints,4) :: Rmaxes4  ! (nautical miles)
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
      REAL(sz) :: Phi
      REAL(sz), DIMENSION(nPoints) :: Bs        
      REAL(sz), DIMENSION(nPoints) :: Phis    ! Correction factor to B and Vh
      REAL(sz), DIMENSION(nPoints) :: VmBL
      REAL(sz), DIMENSION(nPoints,4) :: Bs4        
      REAL(sz), DIMENSION(nPoints,4) :: Phis4    ! Correction factor to B and Vh
      REAL(sz), DIMENSION(nPoints,4) :: VmBL4
      INTEGER,  DIMENSION(nPoints,4) :: QuadFlag4  
      REAL(sz), DIMENSION(nPoints,4) :: QuadIr4        
      REAL(sz), DIMENSION(nQuads) :: VrQuadrant
      REAL(sz), DIMENSION(nQuads) :: radius   ! Wind radii - the distance
                                              ! winds of velocity Vr extend
                                              ! outward from center of storm
                                              ! (nautical miles)
      INTEGER  :: quad                        ! Quadrant counter

      REAL(sz) :: latestRmax                  ! most recently calculated
                                              ! value of fitted rmax
      REAL(sz) :: latestAngle                 ! angle of the most recently
                                              ! calculated node w.r.t. the
                                              ! storm location
      LOGICAL :: useQuadrantVr
      LOGICAL :: useVmaxesBL
      CONTAINS

      !=================
      ! Vortex functions
      !=================

      !=================================================================
      ! Create a new Vortex object.
      !
      ! On input:
      !    Pn           Ambient surface pressure (mb)
      !    Pc           Surface pressure at center of storm (mb)
      !    cLat         Latitude  of storm center (degrees north)
      !    cLon         Longitude of storm center (degrees east )
      !    Vmax         Max sustained wind velocity in storm (knots)
      !
      ! On output:
      !    A new vortex is created with essential parameters calculated.
      !=================================================================
      SUBROUTINE newVortex(pinf, p0, lat, lon, vm)
      REAL(sz), INTENT(IN) :: Pinf
      REAL(sz), INTENT(IN) :: P0
      REAL(sz), INTENT(IN) :: lat
      REAL(sz), INTENT(IN) :: lon
      REAL(sz), INTENT(IN) :: vm
      ! set instance variables
      Pn = pinf
      Pc = p0
      cLat = lat
      cLon = lon
      Vmax = vm
      ! evaluate basic physical params
!P.V 11/04/2022 use absolute coriolis to account for Southern Hemisphere
      corio = ABS(2.0d0 * omega * SIN(deg2rad*cLat))
      B = (Vmax*kt2ms)**2 * RhoAir * EXP(1.d0)
     &             / ((Pn - Pc) * mb2pa)
      B = MAX( MIN(B,2.5d0), 1.0d0) ! limit B to range 1.0->2.5
      ! added for compatibility of calcRmaxes to use with simplified nws20
      Bs(1:6)=B
      VmBL(1:6)=Vmax

      END SUBROUTINE newVortex

      ! A new vortex is created for the full gradient wind balance
      SUBROUTINE newVortexFull(pinf, p0, lat, lon, vm)
      REAL(sz), INTENT(IN) :: Pinf
      REAL(sz), INTENT(IN) :: P0
      REAL(sz), INTENT(IN) :: lat
      REAL(sz), INTENT(IN) :: lon
      REAL(sz), INTENT(IN) :: vm
      ! set instance variables
      Pn = pinf
      Pc = p0
      cLat = lat
      cLon = lon
      Vmax = vm
      ! evaluate basic physical params
!P.V 11/04/2022 use absolute coriolis to account for Southern Hemisphere
      corio = ABS(2.0d0 * omega * SIN(deg2rad*cLat))
      B = (Vmax*kt2ms)**2 * RhoAir * EXP(1.d0)
     &             / ((Pn - Pc) * mb2pa)
      Phi=1.d0
      Bs(1:6)=B
      Phis(1:6)=Phi
      VmBL(1:6)=Vmax
      ! Jie 2013.01 
      ! B = MAX( MIN(B,2.5d0), 1.0d0) ! limit B to range 1.0->2.5

      END SUBROUTINE newVortexFull

      !=================================================================
      ! Set basic parameter for a new Vortex object.
      !
      ! On input:
      !    Pn           Ambient surface pressure (mb)
      !    Pc           Surface pressure at center of storm (mb)
      !    cLat         Latitude  of storm center (degrees north)
      !    cLon         Longitude of storm center (degrees east )
      !
      ! On output:
      !    Aim is to define Pn, Pc, and corio
      !=================================================================
      SUBROUTINE setVortex(pinf, p0, lat, lon)
      REAL(sz), INTENT(IN) :: Pinf
      REAL(sz), INTENT(IN) :: P0
      REAL(sz), INTENT(IN) :: lat
      REAL(sz), INTENT(IN) :: lon
      ! set instance variables
      Pn = pinf
      Pc = p0
      cLat = lat
      cLon = lon
      ! evaluate basic physical params
!P.V 11/04/2022 use absolute coriolis to account for Southern Hemisphere
      corio = ABS(2.0d0 * omega * SIN(deg2rad*cLat))      
      END SUBROUTINE setVortex


      !===============================================================
      ! Calculate the radius of maximum winds for all storm quadrants.
      !
      ! On input:
      !    none
      !
      ! On output:
      !    Rmax    radius of maximum winds (nm) in all quadrants, plus
      !            2 extra values to tie down circular periodicity
      ! Jie 2014.07 Modified with quadrant-varying VmBL, which not only
      ! works for nws19 but for the simplified nws20
      !===============================================================
      SUBROUTINE calcRmaxes()
!              REAL(sz), EXTERNAL  :: func
      REAL(sz)            :: root        ! Radius of maximum winds
      REAL(sz), PARAMETER :: innerRadius = 1.d0
      REAL(sz), PARAMETER :: outerRadius = 400.d0
      REAL(sz), PARAMETER :: accuracy    = 0.0001d0
      REAL(sz), PARAMETER :: zoom        = 0.01d0
      INTEGER , PARAMETER :: itmax = 3
      REAL(sz)            :: r1,r2, r3,r4, dr
      INTEGER             :: n, iter, i
      REAL(sz) :: vicinity

      !-----------------------------
      ! Loop over quadrants of storm
      !-----------------------------
      DO n = 1, nQuads
         ! set B and Vmax values for each quadrant
         ! for nws19, B and Vmax are constant
         ! for simplified nws20, B is constant, while Vmax is not
         B = Bs(n+1)
         Vmax = VmBL(n+1)

         quad = n
         root = -1.d0
         r1 = innerRadius
         r2 = outerRadius
         dr = 1.d0
         DO iter = 1, itmax
            root = findRoot(VhWithCori, r1,r2, dr, r3,r4)
            r1 = r3
            r2 = r4
            dr = dr * zoom
         END DO
         ! determine if Rmax is actually in the vicinity of the
         ! isotach radius that we are using to solve for Rmax,
         ! and if so, take another shot at finding the
         ! Rmax using the gradient wind balance that neglects
         ! coriolis (and is appropriate in the vicinity of Rmax)
         vicinity = abs(root-radius(quad))/root
         if ( (root.lt.0.d0).or.(vicinity.le.0.1d0)) then
            r1 = innerRadius
            r2 = outerRadius
            dr = 1.d0
            DO iter = 1, itmax
               root = findRoot(VhNoCori, r1,r2, dr, r3,r4)
               r1 = r3
               r2 = r4
               dr = dr * zoom
            END DO
         endif
         Rmaxes(n+1) = root
      END DO

      END SUBROUTINE calcRmaxes

      !===============================================================
      ! Calculate the radius of maximum winds for all storm quadrants.
      ! Solving the full gradient wind equation without the assumption
      ! of cyclostrohpic balance.
      !
      ! On input:
      !    none
      !
      ! On output:
      !    Rmax    radius of maximum winds (nm) in all quadrants, plus
      !            2 extra values to tie down circular periodicity
      !
      ! Jie 2013.02 added looping procedures to calculate Bs and Phis
      !===============================================================
      SUBROUTINE calcRmaxesFull()
!              REAL(sz), EXTERNAL  :: func
      REAL(sz)            :: root        ! Radius of maximum winds
      REAL(sz), PARAMETER :: innerRadius = 1.d0
      REAL(sz), PARAMETER :: outerRadius = 500.d0
      REAL(sz), PARAMETER :: accuracy    = 0.0001d0
      REAL(sz), PARAMETER :: zoom        = 0.01d0
      INTEGER , PARAMETER :: itmax = 3
      REAL(sz)            :: r1,r2, r3,r4, dr
      INTEGER             :: n, iter, i, norootflg
      REAL(sz) :: vicinity
      REAL(sz) :: B_new, B_new1
      REAL(sz) :: Phi_new
      INTEGER , PARAMETER :: cont = 400 ! Max # of iterations
      INTEGER :: icont,ibcont                        ! iteration counter
 211  FORMAT(a7,x,i2,x,a38)
      !-----------------------------
      ! Loop over quadrants of storm
      !-----------------------------
      DO n = 1, nQuads
         norootflg=0
      ! initialize B and Phi values for each quadrant
         B = Bs(n+1)
         Phi = Phis(n+1)
         Vmax = VmBL(n+1)
      ! Loop the root-solving process to converge B, for in the
      ! new wind formulation, B is a function of Rmax, Vmax, f, and Phi 
      DO icont = 1, cont ! logical expre. is at the end to exit the loop
         norootflg=0
         quad = n
         root = -1.d0
         r1 = innerRadius
         r2 = outerRadius
         dr = 1.d0
         DO iter = 1, itmax
            root = findRoot(VhWithCoriFull, r1,r2, dr, r3,r4)
            r1 = r3
            r2 = r4
            dr = dr * zoom
         END DO
         ! avoid invalid B value when root is not found        
         if (root.lt.0.d0) then
         !   r1 = innerRadius
         !   r2 = outerRadius
         !   dr = 1.d0
         !   DO iter = 1, itmax
         !      root = findRoot(VhNoCori, r1,r2, dr, r3,r4)
         !      r1 = r3
         !      r2 = r4
         !      dr = dr * zoom
         !   END DO
         root=1.0*radius(quad)
         norootflg=1
         endif
         Rmaxes(n+1) = root  
                
         ! Jie 2013.02
         ! determine if B converges, if yes, break loop and assign 
         ! values to Rmaxes, if not, continue the loop to re-calculate 
         ! root and re-evaluate Bs
         Phi_new = 1+Vmax*kt2ms*root*nm2m*corio/
     &       (B*((Vmax*kt2ms)**2+Vmax*kt2ms*root*nm2m*corio))
         B_new = ((Vmax*kt2ms)**2+Vmax*kt2ms*root*nm2m*corio)
     &       *RhoAir*EXP(Phi_new)/(Phi_new*(Pn-Pc)* mb2pa)
         do ibcont = 1, cont
           B_new1 = B_new
           Phi_new = 1+Vmax*kt2ms*root*nm2m*corio/
     &       (B_new*((Vmax*kt2ms)**2+Vmax*kt2ms*root*nm2m*corio))
           B_new = ((Vmax*kt2ms)**2+Vmax*kt2ms*root*nm2m*corio)
     &       *RhoAir*EXP(Phi_new)/(Phi_new*(Pn-Pc)* mb2pa)
           if (abs(B_new-B_new1).le.0.01d0) then
              EXIT
           endif
         end do
!  debug with aswip          
!         if (ibcont.ge.cont) then
!           write(1111,211) 
!     &       "iquad=",n,  
!     &       "B_new did not fully converge, procede"
!         endif
!  end debug with aswip         
         if (abs(B-B_new).le.0.01d0) then
              EXIT
         endif
         ! update B and Phi for next iteration
         ! warning: modifications made here also affect other subroutines
         B=B_new
         Phi=Phi_new
      END DO !icont = 1, cont 
         ! update to the latest values for aswip output
         Bs(n+1)=B_new
         Phis(n+1)=Phi_new
!  debug with aswip         
!         if (icont.ge.cont) then
!           write(1111,211) 
!     &       "iquad=",n,  
!     &       "B did not fully converge, procede"
!         endif
!  end debug with aswip 
         ! determine if Rmax is actually in the vicinity of the
         ! isotach radius that we are using to solve for Rmax,
         ! and if so, take another shot at finding the
         ! Rmax using the gradient wind equation that neglects
         ! coriolis (and is appropriate in the vicinity of Rmax)
         ! Jie 2013.01
         !vicinity = abs(root-radius(quad))/root
         if (norootflg==1) then
            write(*,*)
     &      "iquad=",n,
     &      "No root found, return dist. to Isotach"  
         endif
 
      END DO  !n = 1, nQuads

      END SUBROUTINE calcRmaxesFull
      
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
      !
      ! Jie 2013.02 Modified to use the full gradient wind eq.
      !====================================================
      FUNCTION VhWithCoriFull(testRmax)

      REAL(sz), INTENT(IN) :: testRmax
      REAL(sz) :: VhWithCoriFull
      REAL(sz) :: thisVr ! the radial wind speed we've been given
      REAL(sz) :: Vh
      !-------------------------
      ! func(x = Rmax) = Vh - Vr
      !-------------------------
      if ( getUseQuadrantVr().eqv..true. ) then
         thisVr = VrQuadrant(quad)
      else
         thisVr = Vr
      endif
      ! Jie 2013.02
      Vh = ms2kt * ( SQRT(
     &            ((Vmax*kt2ms)**2.d0+Vmax*kt2ms*testRmax*nm2m*corio) 
     &                   * (testRmax/radius(quad))**B
     &                   * EXP(Phi*(1.0d0-(testRmax/radius(quad))**B))
     &            + (nm2m*radius(quad)*corio/2.d0)**2.d0)
     &            - nm2m*radius(quad)*corio/2.d0 )
      VhWithCoriFull = Vh - thisVr
      END FUNCTION VhWithCoriFull
      !===============================================================
      
      

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
      FUNCTION VhWithCori(testRmax)

      REAL(sz), INTENT(IN) :: testRmax
      REAL(sz) :: VhWithCori
      REAL(sz) :: thisVr ! the radial wind speed we've been given
      REAL(sz) :: Vh
      !-------------------------
      ! func(x = Rmax) = Vh - Vr
      !-------------------------
      if ( getUseQuadrantVr().eqv..true. ) then
         thisVr = VrQuadrant(quad)
      else
         thisVr = Vr
      endif
      Vh = ms2kt * ( SQRT(
     &            (Vmax*kt2ms)**2.d0 * (testRmax/radius(quad))**B
     &                   * EXP(1.0d0-(testRmax/radius(quad))**B)
     &            + (nm2m*radius(quad)*corio/2.d0)**2.d0)
     &            - nm2m*radius(quad)*corio/2.d0 )
      VhWithCori = Vh - thisVr
      END FUNCTION VhWithCori
      !===============================================================


      !===============================================================
      REAL(sz) FUNCTION VhNoCori(testRmax)
      REAL(sz), INTENT(IN) :: testRmax
      REAL(sz) :: thisVr ! the radial wind speed we've been given
      if ( getUseQuadrantVr().eqv..true. ) then
         thisVr = VrQuadrant(quad)
      else
         thisVr = Vr
      endif
      VhNoCori = ABS(ms2kt *
     &           SQRT( (Vmax*kt2ms)**2.0d0 * (testRmax/radius(quad))**B
     &           * EXP(1 - (testRmax/radius(quad))**B)))
     &           - thisVr

      !===============================================================
      END FUNCTION VhNoCori
      !===============================================================

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
      FUNCTION findRoot(func, x1, x2, dx, a,b)
      REAL(sz), EXTERNAL :: func
      REAL(sz), INTENT(IN) :: x1, x2   ! Search interval [x1,x2]
      REAL(sz), INTENT(IN) :: dx       ! Marching increment
      REAL(sz), INTENT(OUT) :: a, b    ! x values that bracket root
      REAL(sz) :: findRoot             ! The root found

      INTEGER , PARAMETER :: itmax = 400 ! Max # of iterations
      INTEGER :: iter                        ! iteration counter
      REAL(sz) :: fa,fb                      ! function values f(x)
      !
      ! Initialize left side of interval
      a  = x1
      fa = func(a)
      !
      ! March along interval until root is found
      ! or solution diverges.
      findRoot = a
      DO iter = 1,itmax
         b = x1 + iter * dx
         fb = func(b)
         ! Check progress
         IF ((fa*fb < 0.0d0) .OR. (ABS(fb) > ABS(fa))) THEN
            ! Assign root
            IF (ABS(fb) > ABS(fa)) THEN
               findRoot = a
            ELSE
               findRoot = b
            END IF
            EXIT
         END IF

         !
         ! Move right search interval values to left side
         ! for next iteration.
         a  = b
         fa = fb
      END DO

      IF (iter >= itmax) THEN
         PRINT *,
     &              "FUNCTION findRoot: exceeded max # of iterations"
         findRoot = -99999.0
      END IF
      END FUNCTION findRoot
      !=================================================================
      ! Calculate (u,v) wind components and surface pressure from an
      ! asymmetric hurricane wind model.
      !
      ! On input:
      !    lat         Latitude  of nodal point (degrees north)
      !    lon         Longitude of nodal point (degrees east )
      !    uTrans      x component of translational velocity (kts)
      !    vTrans      y component of translational velocity (kts)
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
      SUBROUTINE uvp(lat,lon, uTrans,vTrans, u,v, p)

      REAL(sz), INTENT(IN)  :: lat
      REAL(sz), INTENT(IN)  :: lon
      REAL(sz), INTENT(IN)  :: uTrans
      REAL(sz), INTENT(IN)  :: vTrans

      REAL(sz), INTENT(OUT) :: u
      REAL(sz), INTENT(OUT) :: v
      REAL(sz), INTENT(OUT) :: p

      REAL(sz) :: TransSpdX  !NWS8-style translation speed
      REAL(sz) :: TransSpdY  !NWS8-style translation speed

      REAL(sz) :: avgLat
      REAL(sz) :: dx
      REAL(sz) :: dy
      REAL(sz) :: dist
      REAL(sz) :: rmx
      REAL(sz) :: angle
      REAL(sz) :: speed
      REAL(sz) :: maxspeed
      REAL(sz) :: uf
      REAL(sz) :: vf
      REAL(sz) :: damp
      REAL(sz) :: percentCoriolis
      REAL(sz) :: speedAtRmax
      REAL(sz) :: vmaxFactor
      INTEGER :: i

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
      IF (dist < 1.d0) THEN
         u = 0.0d0
         v = 0.0d0
         p = Pc * mb2pa
         RETURN
      END IF

      dist = m2nm * dist

      angle = 360.0d0 + rad2deg * ATAN2(dx,dy)
      IF (angle > 360.d0) angle = angle - 360.d0
      latestAngle = angle
      rmx =  Rmw(angle)
      latestRmax = rmx

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
      percentCoriolis=1.0d0
      speed = SQRT(
     &           (Vmax*kt2ms)**2.d0 * (rmx/dist)**B
     &                  * EXP(1.d0-(rmx/dist)**B)
     &            + (nm2m*dist*percentCoriolis
     &                           *corio/2.d0)**2.d0)
     &            - nm2m*dist*percentCoriolis*corio/2.d0

      ! calculate the wind speed (m/s) at Rmax, using
      ! equation that includes full coriolis
      speedAtRmax = SQRT(
     &           (Vmax*kt2ms)**2.d0 *
     &                  * EXP(0.d0)
     &             + (nm2m*dist*percentCoriolis
     &                           *corio/2.d0)**2.d0)
     &            - nm2m*dist*percentCoriolis*corio/2.d0

      ! calculate a factor to place the velocity profile so that
      ! it hits vmax
      vmaxFactor = Vmax*kt2ms / speedAtRmax

      ! jgf20111007: Calculate NWS8-like translation speed
      TransSpdX = (abs(speed/speedAtRmax))*uTrans*kt2ms
      TransSpdY = (abs(speed/speedAtRmax))*vTrans*kt2ms

      speed = speed * vmaxFactor

      ! now reduce the wind speed to the surface
      speed = speed * windReduction

!P.V 11/04/2022 to account for Southern Hemisphere
      ! Caution with SH/NH
      if(cLat.lt.0.d0) then
      u =  speed * COS(DEG2RAD * angle)
      v = -speed * SIN(DEG2RAD * angle)
      else
      u = -speed * COS(DEG2RAD * angle)
      v =  speed * SIN(DEG2RAD * angle)
      endif
      !
      ! Alter wind direction by adding a frictional inflow angle
      CALL rotate(u,v, fang(dist,rmx), cLat, uf,vf)
      u = uf
      v = vf
      !
      ! jgf20111007: Add in the translation velocity
      u = u + TransSpdX
      v = v + TransSpdY
      !
      ! convert from 1 minute averaged winds to 10 minute averaged
      ! winds for use in ADCIRC
      u = u * one2ten
      v = v * one2ten

      ! Compute surface pressure from asymmetric hurricane vortex
      p = mb2pa * (Pc + (Pn - Pc) * EXP(-(rmx/dist)**B))

      ! cut off the vortex field after 401nm
      !
      ! TODO: 401nm should be replaced with something less
      ! arbitrary ... and find a better way to blend this
      !if ( dist.gt.401.d0 ) then
      !   u = 0.0d0
      !   v = 0.0d0
      !   p = mb2pa * Pn
      !endif

      END SUBROUTINE uvp

      !=================================================================
      ! Calculate (u,v) wind components and surface pressure from an
      ! asymmetric hurricane wind model.
      !
      ! On input:
      !    pinf         hurricane ambient pressure
      !    p0           hurricane central pressure
      !    idist        dist to hurricane center in nautical mile
      !    irmx         RMW
      !    iangle       Azimuth Angle
      !    iB           Holland B parameter
      !    iVm          vortex maximum velocity at upper boundary
      !    iPhi       vortex correction factor
      !    uTrans        x component of translational velocity (kts)
      !    vTrans        y component of translational velocity (kts)
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
      
      SUBROUTINE uvpr(idist,iangle,irmx,irmx_true,iB,iVm,iPhi,
     &                uTrans,vTrans,geof, u,v,p)
      REAL(sz), INTENT(IN)  :: idist
      REAL(sz), INTENT(IN)  :: iangle
      REAL(sz), INTENT(IN)  :: irmx
      REAL(sz), INTENT(IN)  :: irmx_true
      REAL(sz), INTENT(IN)  :: iB
      REAL(sz), INTENT(IN)  :: iVm
      REAL(sz), INTENT(IN)  :: iPhi
      REAL(sz), INTENT(IN)  :: uTrans
      REAL(sz), INTENT(IN)  :: vTrans
      INTEGER , INTENT(IN)  :: geof

      REAL(sz), INTENT(OUT) :: u
      REAL(sz), INTENT(OUT) :: v
      REAL(sz), INTENT(OUT) :: p

      REAL(sz) :: TransSpdX  !NWS8-style translation speed
      REAL(sz) :: TransSpdY  !NWS8-style translation speed
      REAL(sz) :: rmx
      REAL(sz) :: avgLat
      REAL(sz) :: speed
      REAL(sz) :: uf
      REAL(sz) :: vf
      REAL(sz) :: damp
      REAL(sz) :: percentCoriolis
      INTEGER :: i
      
      rmx=irmx
      B=iB
      Vmax=iVm
      Phi=iPhi
      !----------------------------------------
      ! Handle special case at eye of hurricane
      !----------------------------------------
      ! in eye velocity is zero not translational velocity
      !----------------------------------------
      IF (idist < 1.d0) THEN
         u = 0.0d0
         v = 0.0d0
         p = Pc * mb2pa
         RETURN
      END IF

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
      percentCoriolis=1.0d0
      ! Jie 2014.07
      IF (geof == 1) THEN
         speed = SQRT(
     &           ((Vmax*kt2ms)**2.d0+Vmax*kt2ms*rmx*nm2m
     &             *percentCoriolis*corio) 
     &                  * (rmx/idist)**B
     &                  * EXP(Phi*(1.d0-(rmx/idist)**B))
     &            + (nm2m*idist*percentCoriolis
     &                           *corio/2.d0)**2.d0)
     &            - nm2m*idist*percentCoriolis*corio/2.d0
      ELSE 
         speed = SQRT(
     &           (Vmax*kt2ms)**2.d0 * (rmx/idist)**B
     &            * EXP(1.d0-(rmx/idist)**B)
     &            + (nm2m*idist*percentCoriolis
     &            *corio/2.d0)**2.d0)
     &            - nm2m*idist*percentCoriolis*corio/2.d0      
      ENDIF
      ! jgf20111007: Calculate NWS8-like translation speed
      TransSpdX = (abs(speed/(Vmax*kt2ms)))*uTrans*kt2ms
      TransSpdY = (abs(speed/(Vmax*kt2ms)))*vTrans*kt2ms

      ! now reduce the wind speed to the surface
      speed = speed * windReduction

!P.V 11/05/2022 to account for Southern Hemisphere
      ! Caution with SH/NH
      if(cLat.lt.0.d0) then ! SH
        u =  speed * COS(deg2rad*iangle)
        v = -speed * SIN(deg2rad*iangle)
      else
        u = -speed * COS(deg2rad*iangle)
        v =  speed * SIN(deg2rad*iangle)
      endif
      !
      ! Alter wind direction by adding a frictional inflow angle
      CALL rotate(u,v, fang(idist,irmx_true), cLat, uf,vf)
      u = uf
      v = vf
      !
      ! jgf20111007: Add in the translation velocity
      u = u + TransSpdX
      v = v + TransSpdY
      !
      ! convert from 1 minute averaged winds to 10 minute averaged
      ! winds for use in ADCIRC
      u = u * one2ten
      v = v * one2ten

      ! Compute surface pressure from asymmetric hurricane vortex
      IF (geof == 1) THEN
         p = mb2pa * (Pc + (Pn - Pc) * EXP(-Phi*(rmx/idist)**B))
      ELSE
         p = mb2pa * (Pc + (Pn - Pc) * EXP(-(rmx/idist)**B))
      ENDIF
      ! cut off the vortex field after 401nm
      !
      ! TODO: 401nm should be replaced with something less
      ! arbitrary ... and find a better way to blend this
      !if ( dist.gt.401.d0 ) then
      !   u = 0.0d0
      !   v = 0.0d0
      !   p = mb2pa * Pn
      !endif

      END SUBROUTINE uvpr
      
      !=================================================================
      ! Spatial Interpolation function based on angle and r.
      !
      ! On input:
      !    angle        Azimuthal angle (degrees)
      !    r            Distnace to storm Center (nm)
      !
      ! On output:
      !    interpolated value for Rmax/Vmax/B          
      !
      ! INTEGER valid_Isot is used as a marker to indicate how many isotachs
      ! are available in a certain quadrant
      ! select case(valid_Isot)
      ! case(1): 1 situation
      ! case(2): 3 situations
      ! case(3): 4 situations
      ! case(4): 5 situations                    
      !=================================================================
      REAL(sz) FUNCTION spInterp(angle,dist,opt)
      REAL(sz), INTENT(IN) :: angle,dist  
      INTEGER, INTENT(IN)  :: opt
      REAL(sz), DIMENSION(nPoints,4) :: param 
      REAL(sz) :: temp1, temp2  
      REAL(sz) :: delta_angle
      INTEGER  :: iquad

      if (opt == 1) then
         param = Rmaxes4
      else if (opt == 2) then
         param = Bs4
      else if (opt == 3) then
         param = VmBL4    
      end if   
      
      if (angle.le.45.0d0) then
         iquad = 5
         delta_angle = 45.d0 + angle
      else if (angle.le.135.d0) then
         iquad = 2
         delta_angle = angle - 45.d0
      else if (angle.le.225.d0) then
         iquad = 3
         delta_angle = angle - 135.d0
      else if (angle.le.315.d0) then
         iquad = 4
         delta_angle = angle - 225.d0
      else if (angle.gt.315.d0) then
         iquad = 5
         delta_angle = angle - 315.d0
      endif

      ! nearest neighbor weighted interpolation     
         if ( delta_angle.lt.1.d0 ) then
             spInterp=interpR(param, iquad, dist)
         else if (delta_angle.gt.89.d0) then
             spInterp=interpR(param, iquad+1, dist)
         else
             temp1= interpR(param, iquad, dist)
             temp2= interpR(param, iquad+1, dist)
             spInterp = ( temp1/delta_angle**2.d0
     &             + temp2/(90.0-delta_angle)**2.d0 )
     &             / ( 1.d0/delta_angle**2.d0
     &              + 1.d0/(90.d0-delta_angle)**2.d0 )
         end if

      END FUNCTION spInterp
      
      
      REAL(sz)  FUNCTION interpR(val,quad,r)
      REAL(sz), DIMENSION(nPoints,4),INTENT(IN) :: val 
      INTEGER, INTENT(IN) :: quad
      REAL(sz), INTENT(IN) :: r 
      REAL(sz) :: fac
      INTEGER  :: total_Isot
      
      total_Isot=sum(QuadFlag4(quad,:))
         select case(total_Isot)
         case(1)
           interpR = val(quad,maxloc(QuadFlag4(quad,:),1)) 
         case(2)
            if (r > QuadIr4(quad,1)) then
               interpR = val(quad,1)
            else if (r > QuadIr4(quad,2)) then
               fac = (r-QuadIr4(quad,2))/
     &               (QuadIr4(quad,1)-QuadIr4(quad,2))
               interpR = val(quad,1)* fac+
     &         val(quad,2)* (1-fac)
            else
               interpR = val(quad,2)
            end if
         case(3)
            if (r > QuadIr4(quad,1)) then
               interpR = val(quad,1)
            else if (r > QuadIr4(quad,2)) then
               fac = (r-QuadIr4(quad,2))/
     &               (QuadIr4(quad,1)-QuadIr4(quad,2))
               interpR = val(quad,1)* fac+
     &         val(quad,2)* (1-fac)
            else if (r > QuadIr4(quad,3)) then
               fac = (r-QuadIr4(quad,3))/
     &               (QuadIr4(quad,2)-QuadIr4(quad,3))
               interpR = val(quad,2)* fac+
     &         val(quad,3)* (1-fac)
            else
               interpR = val(quad,3)
            end if
         case(4)
            if (r > QuadIr4(quad,1)) then
               interpR = val(quad,1)
            else if (r > QuadIr4(quad,2)) then
               fac = (r-QuadIr4(quad,2))/
     &               (QuadIr4(quad,1)-QuadIr4(quad,2))
               interpR = val(quad,1)* fac+
     &         val(quad,2)* (1-fac)
            else if (r > QuadIr4(quad,3)) then
               fac = (r-QuadIr4(quad,3))/
     &               (QuadIr4(quad,2)-QuadIr4(quad,3))
               interpR = val(quad,2)* fac+
     &         val(quad,3)* (1-fac)
            else if (r > QuadIr4(quad,4)) then
               fac = (r-QuadIr4(quad,4))/
     &               (QuadIr4(quad,3)-QuadIr4(quad,4))
               interpR = val(quad,3)* fac+
     &         val(quad,4)* (1-fac)
            else
               interpR = val(quad,4)
            end if         
         case default
         ! For whatever reason if our algorithm fails, add the following
         ! line to avoid run-time errors
            interpR = val(quad,maxloc(QuadFlag4(quad,:),1)) 
            write(*,*) "ERROR: interpR failed in nws20get."
         end select  
         
         END function interpR

      !=================================================================
      ! calculate the radius of maximum winds.
      !
      ! On input:
      !    angle        Azimuthal angle (degrees)
      !
      ! On output:
      !    Rmw          Radius of maximum winds (meters) from curve fit
      !                      I DO NOT BELIEVE IT IS IN METERS rjw
      !=================================================================
      REAL(sz) FUNCTION Rmw(angle)
      REAL(sz), INTENT(IN) :: angle
      INTEGER :: base_quadrant
      REAL(sz) :: delta_angle
      if (angle.le.45.0d0) then
         base_quadrant = 5
         delta_angle = 45.d0 + angle
      else if (angle.le.135.d0) then
         base_quadrant = 2
         delta_angle = angle - 45.d0
      else if (angle.le.225.d0) then
         base_quadrant = 3
         delta_angle = angle - 135.d0
      else if (angle.le.315.d0) then
         base_quadrant = 4
         delta_angle = angle - 225.d0
      else if (angle.gt.315.d0) then
         base_quadrant = 5
         delta_angle = angle - 315.d0
      endif
      ! nearest neighbor weighted interpolation
      if ( delta_angle.lt.1.d0 ) then
         Rmw = Rmaxes(base_quadrant)   ! avoid div by zero
      else if ( delta_angle.gt.89.d0 ) then
         Rmw = Rmaxes(base_quadrant+1) ! avoid div by zero
      else
         Rmw = ( Rmaxes(base_quadrant)/delta_angle**2.d0
     &             + Rmaxes(base_quadrant+1)/(90.0-delta_angle)**2.d0 )
     &             / ( 1.d0/delta_angle**2.d0
     &              + 1.d0/(90.d0-delta_angle)**2.d0 )
      endif

      ! linearly interpolate
C            Rmw = (delta_angle/90.0d0)
C     &             *(Rmaxes(base_quadrant+1)-Rmaxes(base_quadrant))
C     &            + Rmaxes(base_quadrant)

      END FUNCTION Rmw

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

      A = SIGN(1.d0, whichWay) * deg2rad * angle
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

      IF (0.0d0 <= r .AND. r < rmx) THEN
         fang = 10.0d0 * r/rmx
      ELSE IF (rmx <= r .AND. r < 1.2d0*rmx) THEN
         fang = 10.0d0 + 75.0d0 * (r/rmx - 1.d0)
      ELSE IF (r >= 1.2d0*rmx) THEN
         fang = 25.0d0
      ELSE
         fang = 0.0d0
      END IF
      END FUNCTION fang


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
C
      REAL(sz) :: dx
      REAL(sz) :: dy
      REAL(sz) :: dt
      REAL(sz) :: avgLat

      avgLat = (latOld + latNew) / 2.d0
      dx = sphericalDistance( ((lonNew-lonOld)*deg2rad),
     &         0.0d0, avgLat, avgLat)
      dy = sphericalDistance(0.0d0,((latNew-latOld)*deg2rad),
     &         latOld,latNew)

      dt = tNew - tOld
      uTrans = sign(dx/dt,(lonNew-lonOld))
      vTrans = sign(dy/dt,(latNew-latOld))

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
      SUBROUTINE fitRmaxes()
      REAL(sz)            :: root        ! Radius of maximum winds
      INTEGER             :: n, iter,i

      !
      ! Generate 2 additional (theta,Rmax) points for curve-fit
      Rmaxes(1) = Rmaxes(5)
      Rmaxes(6) = Rmaxes(2)

      END SUBROUTINE fitRmaxes

      SUBROUTINE fitRmaxes4()
      REAL(sz)            :: root        ! Radius of maximum winds
      INTEGER             :: n, iter,i
      !
      ! Generate 2 additional points for curve-fit
       
      QuadFlag4(1,1:4) = QuadFlag4(5,1:4)
      QuadFlag4(6,1:4) = QuadFlag4(2,1:4)
      
      QuadIr4(1,1:4) = QuadIr4(5,1:4)
      QuadIr4(6,1:4) = QuadIr4(2,1:4)
            
      Rmaxes4(1,1:4) = Rmaxes4(5,1:4)
      Rmaxes4(6,1:4) = Rmaxes4(2,1:4)
      
      Bs4(1,1:4) = Bs4(5,1:4)
      Bs4(6,1:4) = Bs4(2,1:4)
      
      Phis4(1,1:4) = Phis4(5,1:4)
      Phis4(6,1:4) = Phis4(2,1:4)  
      
      VmBL4(1,1:4) = VmBL4(5,1:4)
      VmBL4(6,1:4) = VmBL4(2,1:4) 
 
      END SUBROUTINE fitRmaxes4


      SUBROUTINE setUseQuadrantVr(u)
         LOGICAL, INTENT(IN) :: u
         useQuadrantVr = u
      END SUBROUTINE setUseQuadrantVr

      SUBROUTINE setIsotachWindSpeeds(vrq)
         REAL(sz), DIMENSION(4), INTENT(IN) :: vrq
         VrQuadrant(:) = vrq(:)
      END SUBROUTINE setIsotachWindSpeeds

      SUBROUTINE setIsotachWindSpeed(sp)
         REAL(sz), INTENT(IN) :: sp
         Vr = sp
      END SUBROUTINE setIsotachWindSpeed

      SUBROUTINE setIsotachRadii(ir)
         REAL(sz), DIMENSION(4), INTENT(IN) :: ir
         radius(:) = ir(:)
      END SUBROUTINE setIsotachRadii

      REAL(sz) FUNCTION getShapeParameter()
         getShapeParameter = B
      END FUNCTION getShapeParameter

      SUBROUTINE setShapeParameter(param)
         REAL(sz) :: param
         B = param
      END SUBROUTINE setShapeParameter

      SUBROUTINE getRmaxes(rmaxw)
         REAL(sz), DIMENSION(4), INTENT(OUT) :: rmaxw
         INTEGER :: i
         do i=1,4
            rmaxw(i) = Rmaxes(i+1)
         end do
      END SUBROUTINE getRmaxes

      SUBROUTINE setRmaxes(rmaxw)
         REAL(sz), DIMENSION(4), INTENT(IN) :: rmaxw
         INTEGER :: i
         do i=1,4
            Rmaxes(i+1) = rmaxw(i)
         end do
      END SUBROUTINE setRmaxes

      REAL(sz) FUNCTION getLatestRmax()
         getLatestRmax = latestRmax
      END FUNCTION getLatestRmax

      REAL(sz) FUNCTION getLatestAngle()
         getLatestAngle = latestAngle
      END FUNCTION getLatestAngle

      LOGICAL FUNCTION getUseQuadrantVr()
         getUseQuadrantVr = useQuadrantVr
      END FUNCTION getUseQuadrantVr
      
      !Jie 2013.02
      SUBROUTINE setUseVmaxesBL(u)
         LOGICAL, INTENT(IN) :: u
         useVmaxesBL = u
      END SUBROUTINE setUseVmaxesBL
      
      SUBROUTINE getVmaxesBL(vmaxw)
         REAL(sz), DIMENSION(4), INTENT(OUT) :: vmaxw
         INTEGER :: i
         do i=1,4
            vmaxw(i) = VmBL(i+1)
         end do
      END SUBROUTINE getVmaxesBL

      SUBROUTINE setVmaxesBL(vmaxw)
         REAL(sz), DIMENSION(4), INTENT(IN) :: vmaxw
         INTEGER :: i
         do i=1,4
            VmBL(i+1) = vmaxw(i)
         end do
      END SUBROUTINE setVmaxesBL  
      
      FUNCTION getShapeParameters()
         REAL(sz),DIMENSION(4) :: getShapeParameters
         INTEGER :: i
         do i=1,4
            getShapeParameters(i) = Bs(i+1)
         end do
      END FUNCTION getShapeParameters   
      
      FUNCTION getPhiFactors()
         REAL(sz),DIMENSION(4) :: getPhiFactors
         INTEGER :: i
         do i=1,4
            getPhiFactors(i) = Phis(i+1)
         end do
      END FUNCTION getPhiFactors      

      END MODULE vortex
