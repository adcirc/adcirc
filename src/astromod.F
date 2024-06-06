! Ref: 
!   J. MEEUS, Astronomic algorithms, 2nd Edition, 1991
!  
! Implemented by D. Wirasaet, 2024
!
!  NOTE:  
!    T - JDE or JD centuries 
!
      MODULE ASTROFORMOD
        IMPLICIT NONE

        REAL (8), parameter:: DEG2RAD = 0.017453292519943D0
        REAL (8), parameter:: RAD2DEG = 57.295779513082323D0

        REAL (8), parameter:: SEC2DEG = 2.777777777777778D-4,
     &                        MIN2DEG = 0.01666666666666D0   

        ! Astronomical unit in km
        REAL (8), parameter:: AUDIST = 1.495978707e+11

        INTERFACE ECLIP2EQ
            MODULE PROCEDURE ECLIP2EQ_S, ECLIP2EQ_V    
        END INTERFACE ECLIP2EQ

        TYPE ASTROVAL
           REAL (8):: JD    ! Julain days.
           REAL (8):: T     ! Julian centuries from the Epoch J2000.0 
           REAL (8):: LP    ! Moon's mean longtitude, J2000.0 Epoch. 
           REAL (8):: D     ! Mean elongation of the Moon, J2000.0 Epoch. 
           REAL (8):: M     ! Sun's mean anomaly, J2000.0 Epoch.
           REAL (8):: MP    ! Moon's mean anomaly, J2000 Epoch. 
           REAL (8):: F     ! Moon's argument of latitude, J2000 Epoch. 
           REAL (8):: OG    ! Longtitude of the ascending node of the moon's mean orbit on the ecliptic. 
           REAL (8):: L0    ! Mean longitude of the Sun, refered to the mean equnoix of the date. 
           REAL (8):: DPsi  ! The nutation in longitude.
           REAL (8):: vareps0  ! Mean obliquity of the ecliptic.
           REAL (8):: Dvareps  ! The nutation in obliquity. 
           REAL (8):: vareps ! the obliquity of the ecliptic.
           REAL (8):: eccen ! the eccentricity of the Earth's orbit          
        END TYPE ASTROVAL        
                
      CONTAINS

        FUNCTION COMP_ASTROVAL(  JD  ) RESULT ( ASVAL ) 
          IMPLICIT NONE 

          REAL (8):: JD
          TYPE (ASTROVAL):: ASVAL

          REAL (8):: T, OG, L0, LP
          
          ASVAL%JD = JD ; 

          T = JULIAN_CENTURIES( JD ) ;
          ASVAL%T = JULIAN_CENTURIES( JD )
          ASVAL%D  = MODULO( D_DEG(T ), 360.D0 ) 
          ASVAL%M  = MODULO( M_DEG(T), 360.D0 )
          ASVAL%MP = MODULO( MP_DEG(T), 360.D0 )
          ASVAL%F  = MODULO( F_DEG( T), 360.D0 )

          OG = MODULO( OMEGA_DEG( T ), 360.D0 )
          L0 = MODULO( L0_DEG( T ), 360.D0 )  
          LP = MODULO( LP_DEG( T ), 360.D0 )

          ASVAL%LP = LP
          ASVAL%L0 = L0
          ASVAL%OG = OG 
          ASVAL%DPsi = DeltaPsiL( OG, L0, LP )
          ASVAL%vareps0 = varepsilon0_ecliptic( T ) 
          ASVAL%DVareps = DeltaVarepsL( OG, L0, LP ) 
          ASVAL%vareps = ASVAL%vareps0 + ASVAL%Dvareps 
          ASVAL%eccen =  eccentricity_earth_orbit( T ) 

        END FUNCTION COMP_ASTROVAL

        ! Julian day, p. 61      
        FUNCTION JULIANDAY( DD, MM, YYYY, CALENDAR_TYPE ) RESULT (JD)
          IMPLICIT NONE
         
          REAL(8):: JD
          REAL(8), INTENT(IN) :: DD
          INTEGER:: MM, YYYY
          CHARACTER (LEN=*), optional:: CALENDAR_TYPE
          
          INTEGER:: A, B
          REAL(8):: D, M, Y
          CHARACTER (LEN=9):: CTYPE = 'Gregorian'            
          D = DD 
          M = MM
          Y = YYYY

          IF ( PRESENT(CALENDAR_TYPE) ) THEN
            SELECT CASE(trim(CALENDAR_TYPE))
            CASE ('Julian','JULIAN','julian')
              CTYPE = 'Julian' ;       
            END SELECT        
          ENDIF

          IF ( M <= 2 ) THEN
             Y = Y - 1  ;
             M = M + 12 ;
          ENDIF

          A = FLOOR(Y/100.D0) ;
          B = 2 - A + FLOOR(A/4.D0) ;
          SELECT CASE(trim(CTYPE))
          CASE ('Julian')
            B = 0 ;
          END SELECT 

          JD = FLOOR(365.25D0*(Y + 4716)) + 
     &         FLOOR(30.6001D0*(M+1)) + D + B - 1524.50 ;     
               
        END FUNCTION JULIANDAY
      
        ! - Compute Julian centuries from the Epoch J2000.0 (JDE
        !    2451545)
        FUNCTION JULIAN_CENTURIES(  JDE ) RESULT ( T ) 
          IMPLICIT NONE

          REAL (8):: T
          REAL (8), INTENT(IN):: JDE

          T = (JDE - 2451545.D0)/36525.D0 ;

        END FUNCTION JULIAN_CENTURIES

        ! Sidereal time at Greenwich for any JD.                    !
        ! page. 87                                                  !
        !  JD     = Julian days                                     !
        !  DPSi (optional)   = The nutation in longitude (degrees)  !
        !  vareps (optional) = obliquity of the ecliptic.           !
        FUNCTION GMST_DEG( JD, NUTATION, ASVAL ) RESULT( gmst )
          IMPLICIT NONE

          REAL (8):: gmst
          REAL (8), INTENT(IN):: JD
          LOGICAL, optional, INTENT(IN):: NUTATION
          TYPE (ASTROVAL), optional:: ASVAL 

          REAL (8):: T, JD0, RM
          LOGICAL:: include_nutation = .false. 
          
          ! If nutation is require ! 
          REAL (8)::  LP, L0, OG, DPsi, Dvareps, vareps
          LOGICAL:: have_asval = .false. ; 
          ! Find JD of the date at UT 0h
          RM =  JD - floor(JD) ;

          JD0 = JD ; 
          IF ( RM < 0.5D0 - 1.0e-10 ) THEN
             JD0 = floor(JD) - 0.5 ;        
          ELSE IF ( RM > 0.5D0 + 1.0e-10 ) THEN
             JD0 = floor(JD) + 0.5 ;
          ENDIF
          RM = JD - JD0 ;

          ! Compute T at UT 0h
          T = (JD0 - 2451545.0)/36525.D0 ; 

          ! \Theta0 EQ. (12.2) page 87
          gmst = 100.46061837D0 + 36000.770053608*T
     &         + T*T*(0.000387933 - T/38710000.D0 ) ; 

          ! Mean sidereal time at Greenwich for JD. Page 87
          gmst = gmst + 1.00273790935D0*RM*360.D0 ;

          IF ( present( NUTATION ) ) include_nutation = NUTATION ;

          IF ( include_nutation  ) THEN
            IF ( present(ASVAL) ) THEN
                IF ( abs(JD - ASVAL%JD) < 1.0e-9 ) have_asval = .true.   
            END IF

            IF ( .NOT.  have_asval ) THEN 
              T = JULIAN_CENTURIES( JD ) ; 
              LP = MODULO( LP_DEG( T ), 360.D0 )
              L0 = MODULO( L0_DEG( T ), 360.D0 )  
   
              OG = MODULO( OMEGA_DEG( T ), 360.D0 )
              DPsi = DeltaPsiL( OG, L0, LP )
              vareps = varepsilon0_ecliptic( T )
     &              +   DeltaVarepsL( OG, L0, LP ) 
            ELSE
              DPsi = ASVAL%DPsi ;
              vareps = ASVAL%vareps ;       
            END IF         

            gmst = gmst + DPsi*cos( vareps*DEG2RAD  )/15.D0 ;
          END IF  

!! Eq. 12.4. Yeild slightly differnt results from
!!           the formula above          
!!          gmst = 280.46061837D0 + 360.98564736629*(JD -
!!     &       2451545.D0) + T*T*(0.000387933 - T/38710000.D0) ;  
!!         
          gmst = modulo( gmst, 360.D0 );        
          
        END FUNCTION GMST_DEG 

        ! Moon's mean longtitude, J2000.0 Epoch. page 338
        ELEMENTAL FUNCTION LP_DEG( T ) RESULT ( LP )
          IMPLICIT NONE

          REAL (8):: LP
          REAL (8), INTENT(IN):: T  ! Julian centuries

          LP = 218.3164477D0 + T*( 481267.88123421D0
     &                       + T*(-0.0015786D0 
     &                       + T*( 1.D0/538841.D0  
     &                       + T*(-1.D0/65194000.D0)))) ; 

        END FUNCTION LP_DEG         

        ! Mean elongation of the Moon, J2000.0 Epoch. page 338
        ELEMENTAL FUNCTION D_DEG( T )  RESULT ( D )
          IMPLICIT NONE

          REAL (8):: D
          REAL (8), INTENT(IN):: T

          D = 297.8501921D0 + T*( 445267.1114034D0
     &                      + T*(-0.0018819D0      ! ))  ;
     &                      + T*( 1.D0/545868.D0   ! ))) ;
     &                      + T*(-1.D0/113065000.D0  )))) ;

        END FUNCTION D_DEG         

        ! Sun's mean anomaly, J2000.0 Epoch, page 338
        ELEMENTAL FUNCTION M_DEG( T ) RESULT( M )
          IMPLICIT NONE

          REAL (8):: M
          REAL (8), INTENT(IN):: T

          M = 357.5291092D0 + T*( 35999.0502909
     &                      + T*(-0.0001536D0
     &                      + T*( 1.D0/24490000.D0)))  ;      
        END FUNCTION M_DEG

        ! Moon's mean anomaly, J2000 Epoch. Page 338
        ELEMENTAL FUNCTION MP_DEG( T ) RESULT( MP )
          IMPLICIT NONE

          REAL (8):: MP
          REAL (8), INTENT(IN):: T

          MP = 134.9633964D0 + T*( 477198.8675055D0
     &                       + T*( 0.0087414D0
     &                       + T*( 1.D0/69699.D0
     &                       + T*(-1.D0/14712000.D0)))) ;         
        END FUNCTION MP_DEG

        ! Moon's argument of latitude, J2000 Epoch. Page 338
        ELEMENTAL FUNCTION F_DEG( T ) RESULT (F )
          IMPLICIT NONE

          REAL (8):: F
          REAL (8), INTENT(IN):: T

          F = 93.2720950D0 + T*( 483202.0175233D0
     &                     + T*(-0.0036539
     &                     + T*(-1.D0/3526000
     &                     + T*( 1.D0/863310000.D0)))) ;
        END FUNCTION F_DEG

        ! Longtitude of the ascending node of the moon's mean orbit on
        ! the ecliptic. Page 144
        ELEMENTAL FUNCTION OMEGA_DEG( T ) RESULT( OMEGA )
           IMPLICIT NONE

           REAL (8):: OMEGA
           REAL (8), INTENT(IN):: T

           OMEGA = 125.04452D0 + T*( -1934.136261D0 
     &                         + T*(  0.0020708D0 
     &                         + T*(  1.D0/450000.D0 ))) ;       

        END FUNCTION OMEGA_DEG    
  

        ! Mean longitude of the Sun, refered to the mean equnoix of the
        ! date. Page 163
        ELEMENTAL FUNCTION L0_DEG( T ) RESULT( L0 )
           IMPLICIT NONE

           REAL (8):: L0
           REAL (8), INTENT(IN):: T 

           L0 =  280.46646D0 + T*( 36000.76983D0
     &                       + T*( 0.0003032D0   )) ; 

        END FUNCTION L0_DEG

        ! The eccentricity of the Earth's orbit. Page 163        
        ELEMENTAL FUNCTION eccentricity_earth_orbit( T ) RESULT( eps )
           IMPLICIT NONE 

           REAL(8):: eps
           REAL(8), INTENT(IN):: T

           eps = 0.016708634D0 + T*( -0.000042037D0
     &                                -0.0000001267*T ) ;
        END FUNCTION eccentricity_earth_orbit

        ! The nutation in longitude. p. 144
        !  -- low accuracy
        ! \Delta \Psi 
        ELEMENTAL FUNCTION DeltaPsiL( OMEGA, L, LP ) RESULT( DPsi )
          IMPLICIT NONE

          REAL (8):: DPsi
          REAL (8), INTENT(IN):: OMEGA, L, LP        

          REAL (8):: OG_R, L_R, LP_R

          OG_R = OMEGA*DEG2RAD ;
          L_R = L*DEG2RAD  ;
          LP_R = LP*DEG2RAD ;

          DPsi = - (17.20D0*sec2deg)*sin( OG_R ) 
     &           - (1.32D0*sec2deg)*sin(2.D0*L_R) 
     &           - (0.23D0*sec2deg)*sin(2.D0*LP_R)
     &           + (0.21D0*sec2deg)*sin(2.D0*OG_R) ;     

        END FUNCTION DeltaPsiL        

        ! The nutation in obliquity. p. 144
        !  -- low accuracy
        ! \Delta \varepsilon 
        ELEMENTAL FUNCTION DeltaVarepsL( OMEGA, L, LP ) RESULT( DVareps )
          IMPLICIT NONE

          REAL (8):: DVareps
          REAL (8), INTENT(IN):: OMEGA, L, LP        

          REAL (8):: OG_R, L_R, LP_R

          OG_R = OMEGA*DEG2RAD ;
          L_R = L*DEG2RAD  ;
          LP_R = LP*DEG2RAD ;

          DVareps =  (9.20D0*sec2deg)*cos(OG_R) 
     &             + (0.57D0*sec2deg)*cos(2.D0*L_R) 
     &             + (0.10D0*sec2deg)*cos(2.D0*LP_R) 
     &             - (0.09D0*sec2deg)*cos(2.D0*OG_R) ;    

        END FUNCTION DeltaVarepsL       

        ! Mean obliquity of the ecliptic. Eq. 22.3. Page 147
        ELEMENTAL FUNCTION varepsilon0_ecliptic( T ) RESULT( varepsilon0 )
          IMPLICIT NONE

          REAL (8):: varepsilon0
          REAL (8), INTENT(IN):: T

          REAL (8):: U

          U = 0.01D0*T ;

          varepsilon0 = 23.D0 + 26.D0*min2deg + 21.448D0*sec2deg ;

          varepsilon0 = varepsilon0 + U*( -4680.93*sec2deg
     &                        + U*( -1.55D0 
     &                        + U*( 1999.25D0
     &                        + U*( -51.38D0
     &                        + U*( -249.67D0
     &                        + U*( -39.05D00
     &                        + U*(  7.12D0
     &                        + U*(  27.87D0
     &                        + U*(  5.79D0 
     &                        + U*(  2.45D0 )))))))))) ;   

        END FUNCTION varepsilon0_ecliptic

        
        ! Coordinate transformation                            !
        !   Transformation from ecliptical into EQ coordinates !
        SUBROUTINE ECLIP2EQ_S( RA, DEC, LAMBDA, BETA, varepsilon )
          IMPLICIT NONE

          REAL (8), INTENT(OUT):: RA, DEC
          REAL (8), INTENT(IN):: LAMBDA, BETA, varepsilon

          REAL (8):: NUM, DEN
          
          ! Right ascension 
          NUM = sin( LAMBDA*deg2rad )*cos( varepsilon*deg2rad ) 
     &        -  tan( BETA*deg2rad )*sin( varepsilon*deg2rad ) ;
          DEN = cos( LAMBDA*deg2rad ) ;   

          RA = ATAN2( NUM, DEN ) ; 

          ! Declination
          NUM =  sin( BETA*deg2rad )*cos( varepsilon*deg2rad ) 
     &      +  cos( BETA*deg2rad )*sin( varepsilon*deg2rad)*sin(LAMBDA*deg2rad) ;

          DEC = ASIN( NUM ) ;

          ! Convert to degree 
          RA = MODULO(RA*RAD2DEG,360.D0) ; 
          DEC = DEC*RAD2DEG ; 

          RETURN ;      
        END SUBROUTINE ECLIP2EQ_S


        ! Coordinate transformation                            !
        !   Transformation from ecliptical into EQ coordinates !
        SUBROUTINE ECLIP2EQ_V( RA, DEC, LAMBDA, BETA, varepsilon )
          IMPLICIT NONE

          REAL (8), DIMENSION(:), INTENT(OUT):: RA, DEC
          REAL (8), DIMENSION(:), INTENT(IN):: LAMBDA, BETA, varepsilon

          INTEGER:: len
          REAL (8), DIMENSION(:), ALLOCATABLE:: NUM, DEN
          

          len = size( LAMBDA ) ; 
          ALLOCATE( NUM(len), DEN(len) ) ; 

          ! Right ascension 
          NUM = sin( LAMBDA*deg2rad )*cos( varepsilon*deg2rad ) 
     &        -  tan( BETA*deg2rad )*sin( varepsilon*deg2rad ) ;
          DEN = cos( LAMBDA*deg2rad ) ;   

          RA = ATAN2( NUM, DEN ) ; 

          ! Declination
          NUM =  sin( BETA*deg2rad )*cos( varepsilon*deg2rad ) 
     &      +  cos( BETA*deg2rad )*sin( varepsilon*deg2rad)*sin(LAMBDA*deg2rad) ;

          DEC = ASIN( NUM ) ;

          ! Convert to degree 
          RA = MODULO(RA*RAD2DEG,360.D0) ; 
          DEC = DEC*RAD2DEG ; 

          DEALLOCATE( NUM, DEN ) ; 

          RETURN ;      
        END SUBROUTINE ECLIP2EQ_V

      END MODULE ASTROFORMOD

      
      MODULE SUNPOS
        USE ASTROFORMOD

        IMPLICIT NONE

        INTERFACE SUN_COORDINATES
            MODULE PROCEDURE SUN_COORDINATES_SUB0,
     &                       SUN_COORDINATES_SUB1    
        END INTERFACE SUN_COORDINATES

      CONTAINS

         ! Sun's equation of the center. Page 164
         !  T - Julain centuries of 36525 ephemeris from eoch J2000
         !  M - The mean anomaly of the sun 
         ELEMENTAL FUNCTION SUN_CENTER( T, M ) RESULT( C )
           IMPLICIT NONE

           REAL (8):: C
           REAL (8), INTENT(IN):: T, M
        
           REAL (8):: M_R 

           M_R = M*DEG2RAD ; 

           C = ( 1.914602D0 +  T*(-0.004817D0 - 0.000014D0*T))*sin(M_R) +
     &         ( 0.019993D0 - 0.000101D0*T)*sin(2.D0*M_R) +
     &           0.000289D0*sin(3.D0*M_R) ;

         END FUNCTION SUN_CENTER        


         ! The Sun's true longitude. Page 164.
         !  L0 = the geometric mean longitude of the Sun
         !  C  = the Sun's eq1.495978707×1011uation of the center 
         !  Dpsi (optional) = nutation
         ELEMENTAL FUNCTION SUN_LON( C, L0, DPsi ) RESULT( SLON )
           IMPLICIT NONE

           REAL(8):: SLON
           REAL(8), INTENT(IN):: L0, C
           REAL(8), optional, INTENT(IN):: Dpsi

           SLON =  L0 + C  ;  ! Sun's true longtiude  
       
           IF ( present(DPsi) ) THEN
              SLON = SLON - 0.00569 + Dpsi ;    
           END IF
         END FUNCTION SUN_LON

   
         ! True anomaly of the Sun longitude. Page 164.
         !  M  = the mean anomaly of the Sun
         !  C  = the Sun's equation of the center    
         ELEMENTAL FUNCTION SUN_NU( M, C ) RESULT( nu )
           IMPLICIT NONE

           REAL(8):: nu
           REAL(8), INTENT(IN):: M, C

           nu =  M + C   ;  ! its true anomerly
           
         END FUNCTION SUN_NU


         ! The Sun's radius. Disance bewteen the centers of the Sub and
         ! the Earth in astronomical units 
         !   nu = True anomaly of the Sun's true longitude 
         !   eccen = eccentricity of the Earth's orbit 
         ELEMENTAL  FUNCTION SUN_RADIUS(  nu, eccen ) RESULT( R )
           IMPLICIT NONE 

           REAL (8):: R 
           REAL (8), INTENT(IN):: nu, eccen  
           
           REAL (8):: NUM, DEN, nu_rad

           nu_rad = nu*DEG2RAD ;

           DEN = 1.D0 + eccen*cos(nu_rad) ; 
           NUM = 1.000001018D0*(1.D0 - eccen*eccen) ;

           R = NUM/DEN ; 

        END FUNCTION SUN_RADIUS 

        ! Solar's declination. Page 165.
        !  - SLON = the Sun's longtitude
        !  - VarEps = the obliquity of the ecliptic
        ELEMENTAL FUNCTION SOLAR_DEC( SLON, VarEps ) RESULT( DEC )
          IMPLICIT NONE 

          REAL(8)::  DEC
          REAL(8), INTENT(IN):: SLON, VarEps
           
          REAL (8):: SLON_RAD, VarEps_RAD  

          SLON_RAD = SLON*DEG2RAD ;
          VarEps_RAD = VarEps*DEG2RAD ;
          
          DEC = sin( SLON_RAD )*sin( VarEps_RAD ) ;

          DEC = asin( DEC )*RAD2DEG ;
        END FUNCTION SOLAR_DEC


        ! Solar's right ascension. Page 165.
        !  - SLON = the Sun's longtitude
        !  - VarEps = the obliquity of the ecliptic
        ELEMENTAL FUNCTION SOLAR_RA( SLON, VarEps ) RESULT( RA )
          IMPLICIT NONE

          REAL (8):: RA
          REAL (8), INTENT(IN):: SLON, VarEps

          REAL (8):: NUM, DEN
          REAL (8):: SLON_RAD, VarEps_RAD  

          SLON_RAD = SLON*DEG2RAD ;
          VarEps_RAD = VarEps*DEG2RAD ;

          NUM = cos( VarEps_RAD )*sin( SLON_RAD ) ;
          DEN = cos( SLON_RAD ) 

          RA = MOD(ATAN2( NUM, DEN )*RAD2DEG,360.D0)  ;
        END FUNCTION SOLAR_RA


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
        SUBROUTINE SUN_COORDINATES_SUB1( RA, DEC, Delta, 
     &                                   JD, ASVAL, NUTATION )
          IMPLICIT NONE 

          REAL (8), INTENT(OUT):: RA, DEC, Delta 

          REAL (8), INTENT(IN):: JD
          TYPE (ASTROVAL), optional, INTENT(IN):: ASVAL 
          LOGICAL, optional:: NUTATION  

          REAL (8):: lambda, beta, T, LP, OG, L0, M, C, snu, eccen
          REAL (8):: DPsi, vareps0, DVareps, vareps
          LOGICAL:: have_asval = .false., use_nutation = .true. ; 
        

          IF ( present(ASVAL) ) THEN
             IF ( abs(ASVAL%JD - JD) < 1.0e-9 ) THEN
                have_asval = .true.      
             END IF      
          END IF     

          IF ( present(NUTATION) ) THEN
             use_nutation = NUTATION ;
          END IF   

          IF ( .NOT. have_asval ) THEN
            T = JULIAN_CENTURIES( JD )
            M = MODULO( M_DEG(T), 360.D0 )
            LP = MODULO( LP_DEG(T), 360.D0 ) 
            eccen =  eccentricity_earth_orbit( T ) ;

            OG = MODULO( OMEGA_DEG( T ), 360.D0 )
            L0 = MODULO( L0_DEG( T ), 360.D0 )  
            DPsi = DeltaPsiL( OG, L0, LP ) 
            DVareps = DeltaVarepsL( OG, L0, LP ) ;
   
            vareps0 = varepsilon0_ecliptic( T ) ;
            vareps = vareps0 + Dvareps  ;
          ELSE        
            T = ASVAL%T ; 
            M  = MODULO( ASVAL%M, 360.D0 )
            L0 = MODULO( ASVAL%L0, 360.D0 )
            eccen = ASVAL%eccen  
            DPsi = ASVAL%DPsi 
            vareps = ASVAL%vareps 
          ENDIF     
          C = SUN_CENTER( T, M )

          IF ( use_nutation ) THEN
             CALL SUN_COORDINATES_SUB0( lambda, beta, Delta, 
     &                        L0, M, eccen, C, DPsi )
          ELSE
             CALL SUN_COORDINATES_SUB0( lambda, beta, Delta, 
     &                        L0, M, eccen, C )

          END IF

          ! Apperent right ascension & declination
          CALL ECLIP2EQ( RA, DEC, lambda, beta, vareps ) ;

          RETURN ; 
        END SUBROUTINE SUN_COORDINATES_SUB1        


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
        SUBROUTINE SUN_COORDINATES_SUB0( lambda, beta, Delta, 
     &                        L0, M, eccen, C, DPsi )
          IMPLICIT NONE 

          REAL (8), INTENT(OUT):: lambda, beta, Delta 
          REAL (8), INTENT(IN):: L0, M, C, eccen
          REAL (8), optional, INTENT(IN):: DPsi 
       
          REAL (8):: snu

          beta =  0.D0  ; 
          IF ( .NOT. present(DPsi) ) THEN
            lambda = SUN_LON( C, L0 ) ;
          ELSE  
            lambda = SUN_LON( C, L0, DPsi ) ;  
          END IF    

          snu = SUN_NU( M, C ) ;
          Delta = SUN_RADIUS(  snu, eccen ) ; ! Distance between centers 
                                              ! of the Earth and Sun
          RETURN ; 
        END SUBROUTINE SUN_COORDINATES_SUB0        

      END MODULE SUNPOS

      !
      ! Postion of the Moon. Chapter 47 
      ! Page. 337 -344
      MODULE MOONPOS
        USE ASTROFORMOD

        IMPLICIT NONE

        ! Variables string periodic terms in the the Moon's longtitude,
        ! latitude, and dfistance 
        INTEGER, private, parameter:: NPER = 60 !
        REAL (8), private:: LATARG(NPER,4), LAT_SIN_COEF(NPER)  ! for Moon's Latitude
        REAL (8), private:: LONARG(NPER,4), LON_SIN_COEF(NPER), LON_COS_COEF(NPER) ! For Moon's longitud and distance
    
        INTEGER, private:: LONMEU(NPER), LATMEU(NPER) ! Absolue value of factor of M argument in LATARG, and ONARG  

        PRIVATE:: CAL_EPMUL

        INTERFACE MOON_COORDINATES
            MODULE PROCEDURE MOON_COORDINATES_SUB0,
     &                       MOON_COORDINATES_SUB1    
        END INTERFACE MOON_COORDINATES
      CONTAINS

        ! Coefficient muitplying sine and cosine arguments. page 338
        ELEMENTAL FUNCTION E_MUL_COEF( T ) RESULT( E )
           IMPLICIT NONE
 
           REAL (8):: E
           REAL (8), INTENT(IN):: T
 
           E = 1.D0 + T*( -0.002516D0  
     &              + T*( -0.0000074D0 )) ;   
 
        END FUNCTION E_MUL_COEF        
 
        ! Page. 338 
        ELEMENTAL FUNCTION A1_DEG( T ) RESULT (A1)
           IMPLICIT NONE
 
           REAL(8):: A1
           REAL(8), INTENT(IN):: T
 
           A1 = 119.75D0 + 131.849D0*T ; 
              
        END FUNCTION A1_DEG
 
        ! Page. 338
        ELEMENTAL FUNCTION A2_DEG( T ) RESULT (A2)
           IMPLICIT NONE
 
           REAL(8):: A2
           REAL(8), INTENT(IN):: T
 
           A2 = 53.09D0 + 479264.290D0*T ;
              
        END FUNCTION A2_DEG
 
        ! Page. 338
        ELEMENTAL FUNCTION A3_DEG( T ) RESULT (A3)
           IMPLICIT NONE
 
           REAL(8):: A3
           REAL(8), INTENT(IN):: T
 
           A3 = 313.45D0 + 481266.484D0*T ;
              
        END FUNCTION A3_DEG
         
       ! Coefficients for terms contains angle M. See description on
       ! Page 338   
       FUNCTION CAL_EPMUL( MEU, E ) RESULT( EPMUL )
         IMPLICIT NONE

         REAL (8):: EPMUL(NPER)
         REAL (8), INTENT(IN)::  E
         INTEGER, INTENT(IN)::  MEU(:)

         INTEGER:: I, J

         epmul = 1.D0
         DO I = 1, NPER
           DO J = 1, MEU( I )
             epmul(I) = epmul(I)*E 
           END DO   
         END DO

       END FUNCTION CAL_EPMUL  

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
       SUBROUTINE MOON_COORDINATES_SUB1( RA, DEC, Delta, 
     &                                   JD, ASVAL, NUTATION )
         IMPLICIT NONE

         REAL (8), INTENT(OUT):: RA, DEC, Delta
         REAL (8), INTENT(IN):: JD
         TYPE (ASTROVAL), optional, INTENT(IN):: ASVAL 
         LOGICAL, optional:: NUTATION

         ! local !
         REAL (8):: lambda, beta
         REAL (8):: T, LP, D, M, MP, F, A1, A2, A3, E
         REAL (8):: L0, OG, DPsi, Dvareps, vareps, vareps0
         LOGICAL:: have_asval = .false., use_nutation = .false.  
         REAL (8):: nutation_mul

         ! Check if the astroval 
         IF ( present(ASVAL) ) THEN
            IF ( abs(JD - ASVAL%JD) < 1.e-9 ) have_asval = .true. ; 
         ENDIF
         IF ( present(NUTATION) ) use_nutation = NUTATION ;

         IF ( .NOT. have_asval ) THEN
           T  = JULIAN_CENTURIES( JD )
           LP = MODULO( LP_DEG( T ), 360.D0 )
           D  = MODULO( D_DEG(T ), 360.D0 ) 
           M  = MODULO( M_DEG(T), 360.D0 )
           MP = MODULO( MP_DEG(T), 360.D0 )
           F  = MODULO( F_DEG( T), 360.D0 )
  
           OG = MODULO( OMEGA_DEG( T ), 360.D0 )
           L0 = MODULO( L0_DEG( T ), 360.D0 )  
           DPsi = DeltaPsiL( OG, L0, LP ) 
   
           vareps0 = varepsilon0_ecliptic( T ) ;
           DVareps = DeltaVarepsL( OG, L0, LP ) ;
           vareps = vareps0 + Dvareps  ;
         ELSE        
           T = ASVAL%T ; 
           LP = MODULO( ASVAL%LP, 360.D0 )
           D  = MODULO( ASVAL%D, 360.D0 ) 
           M  = MODULO( ASVAL%M, 360.D0 )
           MP = MODULO( ASVAL%MP, 360.D0 )
           F  = MODULO( ASVAL%F, 360.D0 )

           DPsi = ASVAL%DPsi 
           vareps = ASVAL%vareps 
         ENDIF        

         A1 = MODULO( A1_DEG( T ), 360.D0 ) 
         A2 = MODULO( A2_DEG( T ), 360.D0 )
         A3 = MODULO( A3_DEG( T ), 360.D0 ) 
         E  = E_MUL_COEF( T ) ;
 
         CALL MOON_COORDINATES_SUB0( lambda, beta, Delta, 
     &                LP, D, M, MP, F, E, A1, A2, A3 ) 

         ! PRINT*, "-------- RA & DEC --------" 
         ! Apperent right ascension & declination
         nutation_mul = 1.D0 ;
         IF ( .NOT. use_nutation ) nutation_mul = 0.D0 ; 

         CALL ECLIP2EQ( RA, DEC, lambda + nutation_mul*DPsi, beta, vareps ) ;
      
         ! Geometric right ascension & declination
         ! CALL ECLIP2EQ( RA, DEC, lambda, beta, vareps ) ;
         ! PRINT*, "Geometrical RA  = ", MODULO(RA,360.D0)
         ! PRINT*, "Geometrical DEC = ", DEC
         RETURN ; 
       END SUBROUTINE MOON_COORDINATES_SUB1

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
       SUBROUTINE MOON_COORDINATES_SUB0( lambda, beta, Delta, 
     &                       LP, D, M, MP, F, E, A1, A2, A3 ) 
         IMPLICIT NONE

         REAL (8), INTENT(OUT):: lambda, beta, Delta
         REAL (8), INTENT(IN):: LP, D, M, MP, F
         REAL (8), INTENT(IN):: E, A1, A2, A3

         INTEGER:: I, J
         REAL (8):: suml, sumr, sumb
         REAL (8):: vec(4), argval(NPER), epmul(NPER)

         LOGICAL, save:: first = .TRUE. ;

         ! Initialize the periodic argements 
         IF ( first ) THEN 
            CALL INIT_MOON_POS_ARG() ;
         END IF   
         first = .FALSE. ;

         ! Convert from degree to radian 
         vec = DEG2RAD*(/ D, M, MP, F /) ;

         ! 1. Compute mean longitude !
         argval = MATMUL( LONARG, vec ) ;
         epmul = cal_epmul( LONMEU, E ) ; 

         ! sine argument
         !  \sum l
         argval = epmul*LON_SIN_COEF*sin(argval) ;
         suml = sum(argval) + additive_suml() ; 
  
         ! 2. Compute distance       !
         argval = MATMUL( LONARG, vec ) ; 

         ! cosine argument
         !  \Sum r
         argval = epmul*LON_COS_COEF*cos(argval) ;
         sumr = sum( argval ) ; 

         ! 3. Compute mean latitude !
         argval = MATMUL( LATARG, vec ) ;
         epmul = cal_epmul( LATMEU, E ) ;

         ! sine argument
         !  \Sum b
         argval = epmul*LAT_SIN_COEF*sin(argval) ; 
         sumb = sum( argval ) + additive_sumb() ;    

         ! output   
         lambda = LP + suml/1.0D6 ; ! Longitude in Degree
         beta = sumb/1.D6   ;      ! Lattitude in Degree
         Delta = 385000.56D0 + sumr/1.0D3 ; ! Distance in km  

       CONTAINS
         
         FUNCTION additive_suml( ) RESULT( asuml )
           IMPLICIT NONE

           REAL (8):: asuml 
          
           asuml =  3958.D0*sin( DEG2RAD*A1 )
     &            + 1962.D0*sin( DEG2RAD*(LP - F) ) 
     &            +  318.D0*sin( DEG2RAD*A2 ) ;

         END FUNCTION additive_suml        

         FUNCTION additive_sumb( ) RESULT( asumb )
           IMPLICIT NONE 
           
           REAL (8):: asumb 
 
           asumb = -2235.D0*sin( DEG2RAD*LP ) 
     &             + 382.D0*sin( DEG2RAD*A3  )
     &             + 175.D0*sin( DEG2RAD*(A1 - F) ) 
     &             + 175.D0*sin( DEG2RAD*(A1 + F) )
     &             + 127.D0*sin( DEG2RAD*(LP - MP) )
     &             - 115.D0*sin( DEG2RAD*(LP + MP) ) ;

        END FUNCTION additive_sumb         

       END SUBROUTINE MOON_COORDINATES_SUB0

       ! TABLE 47. A and TABLE 47. B. Pages 339-341
       SUBROUTINE INIT_MOON_POS_ARG( )
          IMPLICIT NONE  

          !  TABLE 47. B. page 341. Periodic terms for the latitude of
          !  the Moon (\Sum b). The unit is 1e-6 degree                
          LATARG( 1,:) = (/    0.D0,   0.D0,   0.D0,   1.D0 /) ;  LAT_SIN_COEF( 1) =   5 128 122.D0 ; 
          LATARG( 2,:) = (/    0.D0,   0.D0,   1.D0,   1.D0 /) ;  LAT_SIN_COEF( 2) =    280 602.D0 ; 
          LATARG( 3,:) = (/    0.D0,   0.D0,   1.D0,  -1.D0 /) ;  LAT_SIN_COEF( 3) =    277693.D0 ; 
          LATARG( 4,:) = (/    2.D0,   0.D0,   0.D0,  -1.D0 /) ;  LAT_SIN_COEF( 4) =    173237.D0 ; 
          LATARG( 5,:) = (/    2.D0,   0.D0,  -1.D0,   1.D0 /) ;  LAT_SIN_COEF( 5) =     55413.D0 ; 
          LATARG( 6,:) = (/    2.D0,   0.D0,  -1.D0,  -1.D0 /) ;  LAT_SIN_COEF( 6) =     46271.D0 ; 
          LATARG( 7,:) = (/    2.D0,   0.D0,   0.D0,   1.D0 /) ;  LAT_SIN_COEF( 7) =     32573.D0 ; 
          LATARG( 8,:) = (/    0.D0,   0.D0,   2.D0,   1.D0 /) ;  LAT_SIN_COEF( 8) =     17198.D0 ; 
          LATARG( 9,:) = (/    2.D0,   0.D0,   1.D0,  -1.D0 /) ;  LAT_SIN_COEF( 9) =      9266.D0 ; 
          LATARG(10,:) = (/    0.D0,   0.D0,   2.D0,  -1.D0 /) ;  LAT_SIN_COEF(10) =      8822.D0 ; 
          LATARG(11,:) = (/    2.D0,  -1.D0,   0.D0,  -1.D0 /) ;  LAT_SIN_COEF(11) =      8216.D0 ; 
          LATARG(12,:) = (/    2.D0,   0.D0,  -2.D0,  -1.D0 /) ;  LAT_SIN_COEF(12) =      4324.D0 ; 
          LATARG(13,:) = (/    2.D0,   0.D0,   1.D0,   1.D0 /) ;  LAT_SIN_COEF(13) =      4200.D0 ; 
          LATARG(14,:) = (/    2.D0,   1.D0,   0.D0,  -1.D0 /) ;  LAT_SIN_COEF(14) =     -3359.D0 ; 
          LATARG(15,:) = (/    2.D0,  -1.D0,  -1.D0,   1.D0 /) ;  LAT_SIN_COEF(15) =      2463.D0 ; 
          LATARG(16,:) = (/    2.D0,  -1.D0,   0.D0,   1.D0 /) ;  LAT_SIN_COEF(16) =      2211.D0 ; 
          LATARG(17,:) = (/    2.D0,  -1.D0,  -1.D0,  -1.D0 /) ;  LAT_SIN_COEF(17) =      2065.D0 ; 
          LATARG(18,:) = (/    0.D0,   1.D0,  -1.D0,  -1.D0 /) ;  LAT_SIN_COEF(18) =     -1870.D0 ; 
          LATARG(19,:) = (/    4.D0,   0.D0,  -1.D0,  -1.D0 /) ;  LAT_SIN_COEF(19) =      1828.D0 ; 
          LATARG(20,:) = (/    0.D0,   1.D0,   0.D0,   1.D0 /) ;  LAT_SIN_COEF(20) =     -1794.D0 ; 
          LATARG(21,:) = (/    0.D0,   0.D0,   0.D0,   3.D0 /) ;  LAT_SIN_COEF(21) =     -1749.D0 ; 
          LATARG(22,:) = (/    0.D0,   1.D0,  -1.D0,   1.D0 /) ;  LAT_SIN_COEF(22) =     -1565.D0 ; 
          LATARG(23,:) = (/    1.D0,   0.D0,   0.D0,   1.D0 /) ;  LAT_SIN_COEF(23) =     -1491.D0 ; 
          LATARG(24,:) = (/    0.D0,   1.D0,   1.D0,   1.D0 /) ;  LAT_SIN_COEF(24) =     -1475.D0 ; 
          LATARG(25,:) = (/    0.D0,   1.D0,   1.D0,  -1.D0 /) ;  LAT_SIN_COEF(25) =     -1410.D0 ; 
          LATARG(26,:) = (/    0.D0,   1.D0,   0.D0,  -1.D0 /) ;  LAT_SIN_COEF(26) =     -1344.D0 ; 
          LATARG(27,:) = (/    1.D0,   0.D0,   0.D0,  -1.D0 /) ;  LAT_SIN_COEF(27) =     -1335.D0 ; 
          LATARG(28,:) = (/    0.D0,   0.D0,   3.D0,   1.D0 /) ;  LAT_SIN_COEF(28) =      1107.D0 ; 
          LATARG(29,:) = (/    4.D0,   0.D0,   0.D0,  -1.D0 /) ;  LAT_SIN_COEF(29) =      1021.D0 ; 
          LATARG(30,:) = (/    4.D0,   0.D0,  -1.D0,   1.D0 /) ;  LAT_SIN_COEF(30) =       833.D0 ; 
          LATARG(31,:) = (/    0.D0,   0.D0,   1.D0,  -3.D0 /) ;  LAT_SIN_COEF(31) =       777.D0 ; 
          LATARG(32,:) = (/    4.D0,   0.D0,  -2.D0,   1.D0 /) ;  LAT_SIN_COEF(32) =       671.D0 ; 
          LATARG(33,:) = (/    2.D0,   0.D0,   0.D0,  -3.D0 /) ;  LAT_SIN_COEF(33) =       607.D0 ; 
          LATARG(34,:) = (/    2.D0,   0.D0,   2.D0,  -1.D0 /) ;  LAT_SIN_COEF(34) =       596.D0 ; 
          LATARG(35,:) = (/    2.D0,  -1.D0,   1.D0,  -1.D0 /) ;  LAT_SIN_COEF(35) =       491.D0 ; 
          LATARG(36,:) = (/    2.D0,   0.D0,  -2.D0,   1.D0 /) ;  LAT_SIN_COEF(36) =      -451.D0 ; 
          LATARG(37,:) = (/    0.D0,   0.D0,   3.D0,  -1.D0 /) ;  LAT_SIN_COEF(37) =       439.D0 ; 
          LATARG(38,:) = (/    2.D0,   0.D0,   2.D0,   1.D0 /) ;  LAT_SIN_COEF(38) =       422.D0 ; 
          LATARG(39,:) = (/    2.D0,   0.D0,  -3.D0,  -1.D0 /) ;  LAT_SIN_COEF(39) =       421.D0 ; 
          LATARG(40,:) = (/    2.D0,   1.D0,  -1.D0,   1.D0 /) ;  LAT_SIN_COEF(40) =      -366.D0 ; 
          LATARG(41,:) = (/    2.D0,   1.D0,   0.D0,   1.D0 /) ;  LAT_SIN_COEF(41) =      -351.D0 ; 
          LATARG(42,:) = (/    4.D0,   0.D0,   0.D0,   1.D0 /) ;  LAT_SIN_COEF(42) =       331.D0 ; 
          LATARG(43,:) = (/    2.D0,  -1.D0,   1.D0,   1.D0 /) ;  LAT_SIN_COEF(43) =       315.D0 ; 
          LATARG(44,:) = (/    2.D0,  -2.D0,   0.D0,  -1.D0 /) ;  LAT_SIN_COEF(44) =       302.D0 ; 
          LATARG(45,:) = (/    0.D0,   0.D0,   1.D0,   3.D0 /) ;  LAT_SIN_COEF(45) =      -283.D0 ; 
          LATARG(46,:) = (/    2.D0,   1.D0,   1.D0,  -1.D0 /) ;  LAT_SIN_COEF(46) =      -229.D0 ; 
          LATARG(47,:) = (/    1.D0,   1.D0,   0.D0,  -1.D0 /) ;  LAT_SIN_COEF(47) =       223.D0 ; 
          LATARG(48,:) = (/    1.D0,   1.D0,   0.D0,   1.D0 /) ;  LAT_SIN_COEF(48) =       223.D0 ; 
          LATARG(49,:) = (/    0.D0,   1.D0,  -2.D0,  -1.D0 /) ;  LAT_SIN_COEF(49) =      -220.D0 ; 
          LATARG(50,:) = (/    2.D0,   1.D0,  -1.D0,  -1.D0 /) ;  LAT_SIN_COEF(50) =      -220.D0 ; 
          LATARG(51,:) = (/    1.D0,   0.D0,   1.D0,   1.D0 /) ;  LAT_SIN_COEF(51) =      -185.D0 ; 
          LATARG(52,:) = (/    2.D0,  -1.D0,  -2.D0,  -1.D0 /) ;  LAT_SIN_COEF(52) =       181.D0 ; 
          LATARG(53,:) = (/    0.D0,   1.D0,   2.D0,   1.D0 /) ;  LAT_SIN_COEF(53) =      -177.D0 ; 
          LATARG(54,:) = (/    4.D0,   0.D0,  -2.D0,  -1.D0 /) ;  LAT_SIN_COEF(54) =       176.D0 ; 
          LATARG(55,:) = (/    4.D0,  -1.D0,  -1.D0,  -1.D0 /) ;  LAT_SIN_COEF(55) =       166.D0 ; 
          LATARG(56,:) = (/    1.D0,   0.D0,   1.D0,  -1.D0 /) ;  LAT_SIN_COEF(56) =      -164.D0 ; 
          LATARG(57,:) = (/    4.D0,   0.D0,   1.D0,  -1.D0 /) ;  LAT_SIN_COEF(57) =       132.D0 ; 
          LATARG(58,:) = (/    1.D0,   0.D0,  -1.D0,  -1.D0 /) ;  LAT_SIN_COEF(58) =      -119.D0 ; 
          LATARG(59,:) = (/    4.D0,  -1.D0,   0.D0,  -1.D0 /) ;  LAT_SIN_COEF(59) =       115.D0 ; 
          LATARG(60,:) = (/    2.D0,  -2.D0,   0.D0,   1.D0 /) ;  LAT_SIN_COEF(60) =       107.D0 ;
          LATMEU(:) = INT( ABS(LATARG(:,2)) ) ;

          ! TABLE 47. A. 
          !  Periodic terms for the longitude \Sum l and distance \Sum r
          !  of the Moon. The unit is 1e-6 for \Sum l and Kilometer for \Sum r 
          LONARG( 1,:) = (/  0.D0,  0.D0,  1.D0,  0.D0 /) ; LON_SIN_COEF( 1) = 6288774.D0 ;  LON_COS_COEF( 1) = -20905355.D0 ; 
          LONARG( 2,:) = (/  2.D0,  0.D0, -1.D0,  0.D0 /) ; LON_SIN_COEF( 2) = 1274027.D0 ;  LON_COS_COEF( 2) =  -3699111.D0 ; 
          LONARG( 3,:) = (/  2.D0,  0.D0,  0.D0,  0.D0 /) ; LON_SIN_COEF( 3) =  658314.D0 ;  LON_COS_COEF( 3) =  -2955968.D0 ; 
          LONARG( 4,:) = (/  0.D0,  0.D0,  2.D0,  0.D0 /) ; LON_SIN_COEF( 4) =  213618.D0 ;  LON_COS_COEF( 4) =   -569925.D0 ; 
          LONARG( 5,:) = (/  0.D0,  1.D0,  0.D0,  0.D0 /) ; LON_SIN_COEF( 5) = -185116.D0 ;  LON_COS_COEF( 5) =     48888.D0 ; 
          LONARG( 6,:) = (/  0.D0,  0.D0,  0.D0,  2.D0 /) ; LON_SIN_COEF( 6) = -114332.D0 ;  LON_COS_COEF( 6) =     -3149.D0 ; 
          LONARG( 7,:) = (/  2.D0,  0.D0, -2.D0,  0.D0 /) ; LON_SIN_COEF( 7) =   58793.D0 ;  LON_COS_COEF( 7) =    246158.D0 ; 
          LONARG( 8,:) = (/  2.D0, -1.D0, -1.D0,  0.D0 /) ; LON_SIN_COEF( 8) =   57066.D0 ;  LON_COS_COEF( 8) =   -152138.D0 ; 
          LONARG( 9,:) = (/  2.D0,  0.D0,  1.D0,  0.D0 /) ; LON_SIN_COEF( 9) =   53322.D0 ;  LON_COS_COEF( 9) =   -170733.D0 ; 
          LONARG(10,:) = (/  2.D0, -1.D0,  0.D0,  0.D0 /) ; LON_SIN_COEF(10) =   45758.D0 ;  LON_COS_COEF(10) =   -204586.D0 ; 
          LONARG(11,:) = (/  0.D0,  1.D0, -1.D0,  0.D0 /) ; LON_SIN_COEF(11) =  -40923.D0 ;  LON_COS_COEF(11) =   -129620.D0 ; 
          LONARG(12,:) = (/  1.D0,  0.D0,  0.D0,  0.D0 /) ; LON_SIN_COEF(12) =  -34720.D0 ;  LON_COS_COEF(12) =    108743.D0 ; 
          LONARG(13,:) = (/  0.D0,  1.D0,  1.D0,  0.D0 /) ; LON_SIN_COEF(13) =  -30383.D0 ;  LON_COS_COEF(13) =    104755.D0 ; 
          LONARG(14,:) = (/  2.D0,  0.D0,  0.D0, -2.D0 /) ; LON_SIN_COEF(14) =   15327.D0 ;  LON_COS_COEF(14) =     10321.D0 ; 
          LONARG(15,:) = (/  0.D0,  0.D0,  1.D0,  2.D0 /) ; LON_SIN_COEF(15) =  -12528.D0 ;  LON_COS_COEF(15) =         0.D0 ; 
          LONARG(16,:) = (/  0.D0,  0.D0,  1.D0, -2.D0 /) ; LON_SIN_COEF(16) =   10980.D0 ;  LON_COS_COEF(16) =     79661.D0 ; 
          LONARG(17,:) = (/  4.D0,  0.D0, -1.D0,  0.D0 /) ; LON_SIN_COEF(17) =   10675.D0 ;  LON_COS_COEF(17) =    -34782.D0 ; 
          LONARG(18,:) = (/  0.D0,  0.D0,  3.D0,  0.D0 /) ; LON_SIN_COEF(18) =   10034.D0 ;  LON_COS_COEF(18) =    -23210.D0 ; 
          LONARG(19,:) = (/  4.D0,  0.D0, -2.D0,  0.D0 /) ; LON_SIN_COEF(19) =    8548.D0 ;  LON_COS_COEF(19) =    -21636.D0 ; 
          LONARG(20,:) = (/  2.D0,  1.D0, -1.D0,  0.D0 /) ; LON_SIN_COEF(20) =   -7888.D0 ;  LON_COS_COEF(20) =     24208.D0 ; 
          LONARG(21,:) = (/  2.D0,  1.D0,  0.D0,  0.D0 /) ; LON_SIN_COEF(21) =   -6766.D0 ;  LON_COS_COEF(21) =     30824.D0 ; 
          LONARG(22,:) = (/  1.D0,  0.D0, -1.D0,  0.D0 /) ; LON_SIN_COEF(22) =   -5163.D0 ;  LON_COS_COEF(22) =     -8379.D0 ; 
          LONARG(23,:) = (/  1.D0,  1.D0,  0.D0,  0.D0 /) ; LON_SIN_COEF(23) =    4987.D0 ;  LON_COS_COEF(23) =    -16675.D0 ; 
          LONARG(24,:) = (/  2.D0, -1.D0,  1.D0,  0.D0 /) ; LON_SIN_COEF(24) =    4036.D0 ;  LON_COS_COEF(24) =    -12831.D0 ; 
          LONARG(25,:) = (/  2.D0,  0.D0,  2.D0,  0.D0 /) ; LON_SIN_COEF(25) =    3994.D0 ;  LON_COS_COEF(25) =    -10445.D0 ; 
          LONARG(26,:) = (/  4.D0,  0.D0,  0.D0,  0.D0 /) ; LON_SIN_COEF(26) =    3861.D0 ;  LON_COS_COEF(26) =    -11650.D0 ; 
          LONARG(27,:) = (/  2.D0,  0.D0, -3.D0,  0.D0 /) ; LON_SIN_COEF(27) =    3665.D0 ;  LON_COS_COEF(27) =     14403.D0 ; 
          LONARG(28,:) = (/  0.D0,  1.D0, -2.D0,  0.D0 /) ; LON_SIN_COEF(28) =   -2689.D0 ;  LON_COS_COEF(28) =     -7003.D0 ; 
          LONARG(29,:) = (/  2.D0,  0.D0, -1.D0,  2.D0 /) ; LON_SIN_COEF(29) =   -2602.D0 ;  LON_COS_COEF(29) =         0.D0 ; 
          LONARG(30,:) = (/  2.D0, -1.D0, -2.D0,  0.D0 /) ; LON_SIN_COEF(30) =    2390.D0 ;  LON_COS_COEF(30) =     10056.D0 ; 
          LONARG(31,:) = (/  1.D0,  0.D0,  1.D0,  0.D0 /) ; LON_SIN_COEF(31) =   -2348.D0 ;  LON_COS_COEF(31) =      6322.D0 ; 
          LONARG(32,:) = (/  2.D0, -2.D0,  0.D0,  0.D0 /) ; LON_SIN_COEF(32) =    2236.D0 ;  LON_COS_COEF(32) =     -9884.D0 ; 
          LONARG(33,:) = (/  0.D0,  1.D0,  2.D0,  0.D0 /) ; LON_SIN_COEF(33) =   -2120.D0 ;  LON_COS_COEF(33) =      5751.D0 ; 
          LONARG(34,:) = (/  0.D0,  2.D0,  0.D0,  0.D0 /) ; LON_SIN_COEF(34) =   -2069.D0 ;  LON_COS_COEF(34) =         0.D0 ; 
          LONARG(35,:) = (/  2.D0, -2.D0, -1.D0,  0.D0 /) ; LON_SIN_COEF(35) =    2048.D0 ;  LON_COS_COEF(35) =     -4950.D0 ; 
          LONARG(36,:) = (/  2.D0,  0.D0,  1.D0, -2.D0 /) ; LON_SIN_COEF(36) =   -1773.D0 ;  LON_COS_COEF(36) =      4130.D0 ; 
          LONARG(37,:) = (/  2.D0,  0.D0,  0.D0,  2.D0 /) ; LON_SIN_COEF(37) =   -1595.D0 ;  LON_COS_COEF(37) =         0.D0 ; 
          LONARG(38,:) = (/  4.D0, -1.D0, -1.D0,  0.D0 /) ; LON_SIN_COEF(38) =    1215.D0 ;  LON_COS_COEF(38) =     -3958.D0 ; 
          LONARG(39,:) = (/  0.D0,  0.D0,  2.D0,  2.D0 /) ; LON_SIN_COEF(39) =   -1110.D0 ;  LON_COS_COEF(39) =         0.D0 ; 
          LONARG(40,:) = (/  3.D0,  0.D0, -1.D0,  0.D0 /) ; LON_SIN_COEF(40) =    -892.D0 ;  LON_COS_COEF(40) =      3258.D0 ; 
          LONARG(41,:) = (/  2.D0,  1.D0,  1.D0,  0.D0 /) ; LON_SIN_COEF(41) =    -810.D0 ;  LON_COS_COEF(41) =      2616.D0 ; 
          LONARG(42,:) = (/  4.D0, -1.D0, -2.D0,  0.D0 /) ; LON_SIN_COEF(42) =     759.D0 ;  LON_COS_COEF(42) =     -1897.D0 ; 
          LONARG(43,:) = (/  0.D0,  2.D0, -1.D0,  0.D0 /) ; LON_SIN_COEF(43) =    -713.D0 ;  LON_COS_COEF(43) =     -2117.D0 ; 
          LONARG(44,:) = (/  2.D0,  2.D0, -1.D0,  0.D0 /) ; LON_SIN_COEF(44) =    -700.D0 ;  LON_COS_COEF(44) =      2354.D0 ; 
          LONARG(45,:) = (/  2.D0,  1.D0, -2.D0,  0.D0 /) ; LON_SIN_COEF(45) =     691.D0 ;  LON_COS_COEF(45) =         0.D0 ; 
          LONARG(46,:) = (/  2.D0, -1.D0,  0.D0, -2.D0 /) ; LON_SIN_COEF(46) =     596.D0 ;  LON_COS_COEF(46) =         0.D0 ; 
          LONARG(47,:) = (/  4.D0,  0.D0,  1.D0,  0.D0 /) ; LON_SIN_COEF(47) =     549.D0 ;  LON_COS_COEF(47) =     -1423.D0 ; 
          LONARG(48,:) = (/  0.D0,  0.D0,  4.D0,  0.D0 /) ; LON_SIN_COEF(48) =     537.D0 ;  LON_COS_COEF(48) =     -1117.D0 ; 
          LONARG(49,:) = (/  4.D0, -1.D0,  0.D0,  0.D0 /) ; LON_SIN_COEF(49) =     520.D0 ;  LON_COS_COEF(49) =     -1571.D0 ; 
          LONARG(50,:) = (/  1.D0,  0.D0, -2.D0,  0.D0 /) ; LON_SIN_COEF(50) =    -487.D0 ;  LON_COS_COEF(50) =     -1739.D0 ; 
          LONARG(51,:) = (/  2.D0,  1.D0,  0.D0, -2.D0 /) ; LON_SIN_COEF(51) =    -399.D0 ;  LON_COS_COEF(51) =         0.D0 ; 
          LONARG(52,:) = (/  0.D0,  0.D0,  2.D0, -2.D0 /) ; LON_SIN_COEF(52) =    -381.D0 ;  LON_COS_COEF(52) =     -4421.D0 ; 
          LONARG(53,:) = (/  1.D0,  1.D0,  1.D0,  0.D0 /) ; LON_SIN_COEF(53) =     351.D0 ;  LON_COS_COEF(53) =         0.D0 ; 
          LONARG(54,:) = (/  3.D0,  0.D0, -2.D0,  0.D0 /) ; LON_SIN_COEF(54) =    -340.D0 ;  LON_COS_COEF(54) =         0.D0 ; 
          LONARG(55,:) = (/  4.D0,  0.D0, -3.D0,  0.D0 /) ; LON_SIN_COEF(55) =     330.D0 ;  LON_COS_COEF(55) =         0.D0 ; 
          LONARG(56,:) = (/  2.D0, -1.D0,  2.D0,  0.D0 /) ; LON_SIN_COEF(56) =     327.D0 ;  LON_COS_COEF(56) =         0.D0 ; 
          LONARG(57,:) = (/  0.D0,  2.D0,  1.D0,  0.D0 /) ; LON_SIN_COEF(57) =    -323.D0 ;  LON_COS_COEF(57) =      1165.D0 ; 
          LONARG(58,:) = (/  1.D0,  1.D0, -1.D0,  0.D0 /) ; LON_SIN_COEF(58) =     299.D0 ;  LON_COS_COEF(58) =         0.D0 ; 
          LONARG(59,:) = (/  2.D0,  0.D0,  3.D0,  0.D0 /) ; LON_SIN_COEF(59) =     294.D0 ;  LON_COS_COEF(59) =         0.D0 ; 
          LONARG(60,:) = (/  2.D0,  0.D0, -1.D0, -2.D0 /) ; LON_SIN_COEF(60) =       0.D0 ;  LON_COS_COEF(60) =      8752.D0 ;
          LONMEU(:) = INT(ABS(LONARG(:,2))) ; 

          RETURN    
       END SUBROUTINE INIT_MOON_POS_ARG         
        
      END MODULE MOONPOS
     
      ! Driver subroutine 
      MODULE MOON_SUN_COORS
         USE MOONPOS
         USE SUNPOS
         USE ASTROFORMOD 

         LOGICAL, private:: INCLUDE_NUTATION = .true. 
         TYPE (ASTROVAL), private:: ASTRO_VALUES

      CONTAINS         

        SUBROUTINE SET_NUTATION( NUTATION )
           IMPLICIT NONE

           LOGICAL:: NUTATION

           INCLUDE_NUTATION = NUTATION 

           RETURN ;     
        END SUBROUTINE SET_NUTATION 

        !
        ! Compute the geocentric coodinates of the Moon, Sun
        !
        ! Output:
        !   MOON_POS = Moon coordinates and distance
        !            = (/ RA, DEC, Delta /)
        !   SUN_POS  = Sun coordinates and distance
        !            = (/ RA, DEC, Delta /)
        ! Input:
        !   JD = Julian days
        ! 
        ! NOTE: 
        !   - By default, the nutation is accounted for and tthe output 
        !     coodinates are the apparent coordinates.  
        !      
        !   - To ignore the nutation and get the geometrical
        !     coordinates, call the subroutine 
        !  
        !         CALL SET_NUTATION( .FALSE. )
        ! 
        !     prior to the use of this subroutine. 
        !        
        SUBROUTINE HEAVENLY_OBJS_COORDS( MOON_POS, SUN_POS, JD )
           IMPLICIT NONE

           REAL (8), INTENT(IN):: JD
           REAL (8), INTENT(OUT):: MOON_POS(3), SUN_POS(3)

           ASTRO_VALUES = COMP_ASTROVAL(  JD  ) ;

           CALL MOON_COORDINATES(  MOON_POS(1), MOON_POS(2), 
     &          MOON_POS(3), JD, ASTRO_VALUES, INCLUDE_NUTATION )

           CALL SUN_COORDINATES( SUN_POS(1), SUN_POS(2), 
     &          SUN_POS(3), JD, ASTRO_VALUES, INCLUDE_NUTATION )

           RETURN ;            
        END SUBROUTINE HEAVENLY_OBJS_COORDS  

      END MODULE MOON_SUN_COORS        
