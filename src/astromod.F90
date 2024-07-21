! Ref: 
!   J. MEEUS, Astronomic algorithms, 2nd Edition, 1991
!  
! Implemented by D. Wirasaet, 2024
!
!  NOTE:  
!    T - JDE or JD centuries 
!
      MODULE ASTROFORMOD

        USE ADC_CONSTANTS       

        IMPLICIT NONE

        REAL (8), parameter:: SEC2DEG = 2.777777777777778D-4,&
                              MIN2DEG = 0.01666666666666D0   

        ! Astronomical unit in km
        REAL (8), parameter:: AUDIST = 1.495978707e+8
        REAL (8), parameter:: km2AU = 1.D0/AUDIST 

        ! Ratio of Sun/Moon mass to Earth mass: Wikipedia ! 
        REAL (8), parameter:: MassRatioSunEarth = 332946.0487       
        REAL (8), parameter:: MassRatioMoonEarth = 1.D0/81.3005678D0 

        ! Keep the value in two saperate numbers: significant and exponent
        REAL (8), parameter:: EarthRadiuskm(2) = (/ 6371.0088D0, 1.0D0  /)   ! km
        REAL (8), parameter:: EarthRadiusm(2)  = (/ 6.3710088D0, 1.0D+6 /)   ! m
        REAL (8), parameter:: EarthRadiusAU(2) = (/ 4.258756338033895D0, 1.0D-5 /)  ! in AU   

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
          
          ASVAL%JD = JD  

          T = JULIAN_CENTURIES( JD ) 
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
        REAL(8) FUNCTION JULIANDAY( DD, MM, YYYY, CALENDAR_TYPE ) RESULT (JD)
          IMPLICIT NONE

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
              CTYPE = 'Julian'        
            END SELECT        
          ENDIF

          IF ( M <= 2 ) THEN
             Y = Y - 1  
             M = M + 12 
          ENDIF

          A = FLOOR(Y/100.D0) 
          B = 2 - A + FLOOR(A/4.D0) 
          SELECT CASE(trim(CTYPE))
          CASE ('Julian')
            B = 0 
          END SELECT 

          JD = FLOOR(365.25D0*(Y + 4716)) + FLOOR(30.6001D0*(M+1)) + D + B - 1524.50      
               
        END FUNCTION JULIANDAY
      
        ! - Compute Julian centuries from the Epoch J2000.0 (JDE
        !    2451545)
        FUNCTION JULIAN_CENTURIES(  JDE ) RESULT ( T ) 
          IMPLICIT NONE

          REAL (8):: T
          REAL (8), INTENT(IN):: JDE

          T = (JDE - 2451545.D0)/36525.D0 

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
          LOGICAL:: have_asval = .false.  
          ! Find JD of the date at UT 0h
          RM =  JD - floor(JD) 

          JD0 = JD  
          IF ( RM < 0.5D0 - 1.0e-10 ) THEN
             JD0 = floor(JD) - 0.5         
          ELSE IF ( RM > 0.5D0 + 1.0e-10 ) THEN
             JD0 = floor(JD) + 0.5 
          ENDIF
          RM = JD - JD0 

          ! Compute T at UT 0h
          T = (JD0 - 2451545.0)/36525.D0  

          ! \Theta0 EQ. (12.2) page 87
          gmst = 100.46061837D0 + 36000.770053608*T &
               + T*T*(0.000387933 - T/38710000.D0 )  

          ! Mean sidereal time at Greenwich for JD. Page 87
          gmst = gmst + 1.00273790935D0*RM*360.D0 

          IF ( present( NUTATION ) ) include_nutation = NUTATION 

          IF ( include_nutation  ) THEN
            IF ( present(ASVAL) ) THEN
                IF ( abs(JD - ASVAL%JD) < 1.0e-9 ) have_asval = .true.   
            END IF

            IF ( .NOT.  have_asval ) THEN 
              T = JULIAN_CENTURIES( JD )  
              LP = MODULO( LP_DEG( T ), 360.D0 )
              L0 = MODULO( L0_DEG( T ), 360.D0 )  
   
              OG = MODULO( OMEGA_DEG( T ), 360.D0 )
              DPsi = DeltaPsiL( OG, L0, LP )
              vareps = varepsilon0_ecliptic( T ) + DeltaVarepsL( OG, L0, LP ) 
            ELSE
              DPsi = ASVAL%DPsi 
              vareps = ASVAL%vareps        
            END IF         

            gmst = gmst + DPsi*cos( vareps*DEG2RAD  )/15.D0 
          END IF  

!! Eq. 12.4. Yeild slightly differnt results from
!!           the formula above          
!!          gmst = 280.46061837D0 + 360.98564736629*(JD -
!!     &       2451545.D0) + T*T*(0.000387933 - T/38710000.D0)   
!!         
          gmst = modulo( gmst, 360.D0 )        
          
        END FUNCTION GMST_DEG 

        ! Moon's mean longtitude, J2000.0 Epoch. page 338
        ELEMENTAL FUNCTION LP_DEG( T ) RESULT ( LP )
          IMPLICIT NONE

          REAL (8):: LP
          REAL (8), INTENT(IN):: T  ! Julian centuries

          LP = 218.3164477D0 + T*( 481267.88123421D0 &
                             + T*(-0.0015786D0 &
                             + T*( 1.D0/538841.D0 &
                             + T*(-1.D0/65194000.D0))))

        END FUNCTION LP_DEG         

        ! Mean elongation of the Moon, J2000.0 Epoch. page 338
        ELEMENTAL FUNCTION D_DEG( T )  RESULT ( D )
          IMPLICIT NONE

          REAL (8):: D
          REAL (8), INTENT(IN):: T

          D = 297.8501921D0 + T*( 445267.1114034D0 &
                            + T*(-0.0018819D0      &
                            + T*( 1.D0/545868.D0   &
                            + T*(-1.D0/113065000.D0  )))) 

        END FUNCTION D_DEG         

        ! Sun's mean anomaly, J2000.0 Epoch, page 338
        ELEMENTAL FUNCTION M_DEG( T ) RESULT( M )
          IMPLICIT NONE

          REAL (8):: M
          REAL (8), INTENT(IN):: T

          M = 357.5291092D0 + T*( 35999.0502909 &
                            + T*(-0.0001536D0   &
                            + T*( 1.D0/24490000.D0)))        
        END FUNCTION M_DEG

        ! Moon's mean anomaly, J2000 Epoch. Page 338
        ELEMENTAL FUNCTION MP_DEG( T ) RESULT( MP )
          IMPLICIT NONE

          REAL (8):: MP
          REAL (8), INTENT(IN):: T

          MP = 134.9633964D0 + T*( 477198.8675055D0 &
                             + T*( 0.0087414D0 &
                             + T*( 1.D0/69699.D0 &
                             + T*(-1.D0/14712000.D0))))          
        END FUNCTION MP_DEG

        ! Moon's argument of latitude, J2000 Epoch. Page 338
        ELEMENTAL FUNCTION F_DEG( T ) RESULT (F )
          IMPLICIT NONE

          REAL (8):: F
          REAL (8), INTENT(IN):: T

          F = 93.2720950D0 + T*( 483202.0175233D0 &
                           + T*(-0.0036539 &
                           + T*(-1.D0/3526000 &
                           + T*( 1.D0/863310000.D0)))) 
        END FUNCTION F_DEG

        ! Longtitude of the ascending node of the moon's mean orbit on
        ! the ecliptic. Page 144
        REAL(8) ELEMENTAL FUNCTION OMEGA_DEG( T ) RESULT( OMEGA )
           IMPLICIT NONE
           REAL (8), INTENT(IN):: T
           OMEGA = 125.04452D0 + T*( -1934.136261D0 &
                               + T*(  0.0020708D0  &
                               + T*(  1.D0/450000.D0 )))        
        END FUNCTION OMEGA_DEG    
  
        ! Mean longitude of the Sun, refered to the mean equnoix of the
        ! date. Page 163
        REAL(8) ELEMENTAL FUNCTION L0_DEG( T ) RESULT( L0 )
           IMPLICIT NONE
           REAL (8), INTENT(IN) :: T 
           L0 =  280.46646D0 + T*( 36000.76983D0 + T*( 0.0003032D0   ))  
        END FUNCTION L0_DEG

        ! The eccentricity of the Earth's orbit. Page 163        
        REAL(8) ELEMENTAL FUNCTION eccentricity_earth_orbit( T ) RESULT( eps )
           IMPLICIT NONE 
           REAL(8), INTENT(IN):: T
           eps = 0.016708634D0 + T*( -0.000042037D0 &
                                     -0.0000001267*T ) 
        END FUNCTION eccentricity_earth_orbit

        ! The nutation in longitude. p. 144
        !  -- low accuracy
        ! \Delta \Psi 
        REAL(8) ELEMENTAL FUNCTION DeltaPsiL( OMEGA, L, LP ) RESULT( DPsi )
          IMPLICIT NONE
          REAL (8), INTENT(IN) :: OMEGA, L, LP        
          REAL (8) :: OG_R, L_R, LP_R

          OG_R = OMEGA*DEG2RAD 
          L_R = L*DEG2RAD  
          LP_R = LP*DEG2RAD 

          DPsi = - (17.20D0*sec2deg)*sin( OG_R ) &
                 - (1.32D0*sec2deg)*sin(2.D0*L_R) &
                 - (0.23D0*sec2deg)*sin(2.D0*LP_R) &
                 + (0.21D0*sec2deg)*sin(2.D0*OG_R)      
        END FUNCTION DeltaPsiL        

        ! The nutation in obliquity. p. 144
        !  -- low accuracy
        ! \Delta \varepsilon 
        REAL(8) ELEMENTAL FUNCTION DeltaVarepsL( OMEGA, L, LP ) RESULT( DVareps )
          IMPLICIT NONE
          REAL (8), INTENT(IN):: OMEGA, L, LP        
          REAL (8) :: OG_R, L_R, LP_R

          OG_R = OMEGA*DEG2RAD 
          L_R = L*DEG2RAD  
          LP_R = LP*DEG2RAD 

          DVareps =  (9.20D0*sec2deg)*cos(OG_R) &
                   + (0.57D0*sec2deg)*cos(2.D0*L_R) &
                   + (0.10D0*sec2deg)*cos(2.D0*LP_R) &
                   - (0.09D0*sec2deg)*cos(2.D0*OG_R)     

        END FUNCTION DeltaVarepsL       

        ! Mean obliquity of the ecliptic. Eq. 22.3. Page 147
        REAL(8) ELEMENTAL FUNCTION varepsilon0_ecliptic( T ) RESULT( varepsilon0 )
          IMPLICIT NONE
          REAL (8), INTENT(IN):: T
          REAL (8) :: U

          U = 0.01D0*T 

          varepsilon0 = 23.D0 + 26.D0*min2deg + 21.448D0*sec2deg 

          varepsilon0 = varepsilon0 + U*( -4680.93*sec2deg &
                              + U*( -1.55D0 &
                              + U*( 1999.25D0 &
                              + U*( -51.38D0 &
                              + U*( -249.67D0 &
                              + U*( -39.05D00 &
                              + U*(  7.12D0 &
                              + U*(  27.87D0 &
                              + U*(  5.79D0 &
                              + U*(  2.45D0 ))))))))))    

        END FUNCTION varepsilon0_ecliptic

        
        ! Coordinate transformation                            !
        !   Transformation from ecliptical into EQ coordinates !
        SUBROUTINE ECLIP2EQ_S( RA, DEC, LAMBDA, BETA, varepsilon )
          IMPLICIT NONE

          REAL (8), INTENT(OUT):: RA, DEC
          REAL (8), INTENT(IN):: LAMBDA, BETA, varepsilon

          REAL (8):: NUM, DEN
          
          ! Right ascension 
          NUM = sin( LAMBDA*deg2rad )*cos( varepsilon*deg2rad ) &
              -  tan( BETA*deg2rad )*sin( varepsilon*deg2rad ) 
          DEN = cos( LAMBDA*deg2rad )    

          RA = ATAN2( NUM, DEN )  

          ! Declination
          NUM =  sin( BETA*deg2rad )*cos( varepsilon*deg2rad ) &
            +  cos( BETA*deg2rad )*sin( varepsilon*deg2rad)*sin(LAMBDA*deg2rad) 

          DEC = ASIN( NUM ) 

          ! Convert to degree 
          RA = MODULO(RA*RAD2DEG,360.D0)  
          DEC = DEC*RAD2DEG  

          RETURN       
        END SUBROUTINE ECLIP2EQ_S


        ! Coordinate transformation                            !
        !   Transformation from ecliptical into EQ coordinates !
        SUBROUTINE ECLIP2EQ_V( RA, DEC, LAMBDA, BETA, varepsilon )
          IMPLICIT NONE

          REAL (8), DIMENSION(:), INTENT(OUT):: RA, DEC
          REAL (8), DIMENSION(:), INTENT(IN):: LAMBDA, BETA, varepsilon

          INTEGER:: len
          REAL (8), DIMENSION(:), ALLOCATABLE:: NUM, DEN
          

          len = size( LAMBDA )  
          ALLOCATE( NUM(len), DEN(len) )  

          ! Right ascension 
          NUM = sin( LAMBDA*deg2rad )*cos( varepsilon*deg2rad ) &
              -  tan( BETA*deg2rad )*sin( varepsilon*deg2rad ) 
          DEN = cos( LAMBDA*deg2rad )    

          RA = ATAN2( NUM, DEN )  

          ! Declination
          NUM =  sin( BETA*deg2rad )*cos( varepsilon*deg2rad ) &
            +  cos( BETA*deg2rad )*sin( varepsilon*deg2rad)*sin(LAMBDA*deg2rad) 

          DEC = ASIN( NUM ) 

          ! Convert to degree 
          RA = MODULO(RA*RAD2DEG,360.D0)  
          DEC = DEC*RAD2DEG  

          DEALLOCATE( NUM, DEN )  

          RETURN       
        END SUBROUTINE ECLIP2EQ_V

      END MODULE ASTROFORMOD

      
      MODULE SUNPOS
        USE ASTROFORMOD

        IMPLICIT NONE

        INTERFACE SUN_COORDINATES
            MODULE PROCEDURE SUN_COORDINATES_SUB0, SUN_COORDINATES_SUB1    
        END INTERFACE SUN_COORDINATES

      CONTAINS

         ! Sun's equation of the center. Page 164
         !  T - Julain centuries of 36525 ephemeris from eoch J2000
         !  M - The mean anomaly of the sun 
         REAL(8) ELEMENTAL FUNCTION SUN_CENTER( T, M ) RESULT( C )
           IMPLICIT NONE
           REAL (8), INTENT(IN):: T, M
           REAL (8):: M_R 

           M_R = M*DEG2RAD  

           C = ( 1.914602D0 +  T*(-0.004817D0 - 0.000014D0*T))*sin(M_R) + &
               ( 0.019993D0 - 0.000101D0*T)*sin(2.D0*M_R) + &
                 0.000289D0*sin(3.D0*M_R) 

         END FUNCTION SUN_CENTER        


         ! The Sun's true longitude. Page 164.
         !  L0 = the geometric mean longitude of the Sun
         !  C  = the Sun's eq1.495978707×1011uation of the center 
         !  Dpsi (optional) = nutation
         REAL(8) ELEMENTAL FUNCTION SUN_LON( C, L0, DPsi ) RESULT( SLON )
           IMPLICIT NONE
           REAL(8), INTENT(IN) :: L0, C
           REAL(8), optional, INTENT(IN):: Dpsi

           SLON =  L0 + C    ! Sun's true longtiude  
       
           IF ( present(DPsi) ) THEN
              SLON = SLON - 0.00569 + Dpsi     
           END IF
         END FUNCTION SUN_LON

   
         ! True anomaly of the Sun longitude. Page 164.
         !  M  = the mean anomaly of the Sun
         !  C  = the Sun's equation of the center    
         REAL(8) ELEMENTAL FUNCTION SUN_NU( M, C ) RESULT( nu )
           IMPLICIT NONE
           REAL(8), INTENT(IN):: M, C
           nu =  M + C     ! its true anomerly
         END FUNCTION SUN_NU


         ! The Sun's radius. Disance bewteen the centers of the Sub and
         ! the Earth in astronomical units 
         !   nu = True anomaly of the Sun's true longitude 
         !   eccen = eccentricity of the Earth's orbit 
         REAL(8) ELEMENTAL FUNCTION SUN_RADIUS(  nu, eccen ) RESULT( R )
           IMPLICIT NONE 
           REAL (8), INTENT(IN):: nu, eccen  
           REAL (8):: NUM, DEN, nu_rad

           nu_rad = nu*DEG2RAD 

           DEN = 1.D0 + eccen*cos(nu_rad)  
           NUM = 1.000001018D0*(1.D0 - eccen*eccen) 

           R = NUM/DEN  

        END FUNCTION SUN_RADIUS 

        ! Solar's declination. Page 165.
        !  - SLON = the Sun's longtitude
        !  - VarEps = the obliquity of the ecliptic
        REAL(8) ELEMENTAL FUNCTION SOLAR_DEC( SLON, VarEps ) RESULT( DEC )
          IMPLICIT NONE 

          REAL(8), INTENT(IN):: SLON, VarEps
           
          REAL (8):: SLON_RAD, VarEps_RAD  

          SLON_RAD = SLON*DEG2RAD 
          VarEps_RAD = VarEps*DEG2RAD 
          
          DEC = sin( SLON_RAD )*sin( VarEps_RAD ) 

          DEC = asin( DEC )*RAD2DEG 
        END FUNCTION SOLAR_DEC


        ! Solar's right ascension. Page 165.
        !  - SLON = the Sun's longtitude
        !  - VarEps = the obliquity of the ecliptic
        REAL(8) ELEMENTAL FUNCTION SOLAR_RA( SLON, VarEps ) RESULT( RA )
          IMPLICIT NONE

          REAL (8), INTENT(IN):: SLON, VarEps

          REAL (8):: NUM, DEN
          REAL (8):: SLON_RAD, VarEps_RAD  

          SLON_RAD = SLON*DEG2RAD 
          VarEps_RAD = VarEps*DEG2RAD 

          NUM = cos( VarEps_RAD )*sin( SLON_RAD ) 
          DEN = cos( SLON_RAD ) 

          RA = MOD(ATAN2( NUM, DEN )*RAD2DEG,360.D0)  
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
        SUBROUTINE SUN_COORDINATES_SUB1( RA, DEC, Delta, JD, ASVAL, NUTATION )
          IMPLICIT NONE 

          REAL (8), INTENT(OUT):: RA, DEC, Delta 

          REAL (8), INTENT(IN):: JD
          TYPE (ASTROVAL), optional, INTENT(IN):: ASVAL 
          LOGICAL, optional:: NUTATION  

          REAL (8):: lambda, beta, T, LP, OG, L0, M, C, snu, eccen
          REAL (8):: DPsi, vareps0, DVareps, vareps
          LOGICAL:: have_asval = .false., use_nutation = .true.  
        

          IF ( present(ASVAL) ) THEN
             IF ( abs(ASVAL%JD - JD) < 1.0e-9 ) THEN
                have_asval = .true.      
             END IF      
          END IF     

          IF ( present(NUTATION) ) THEN
             use_nutation = NUTATION 
          END IF   

          IF ( .NOT. have_asval ) THEN
            T = JULIAN_CENTURIES( JD )
            M = MODULO( M_DEG(T), 360.D0 )
            LP = MODULO( LP_DEG(T), 360.D0 ) 
            eccen =  eccentricity_earth_orbit( T ) 

            OG = MODULO( OMEGA_DEG( T ), 360.D0 )
            L0 = MODULO( L0_DEG( T ), 360.D0 )  
            DPsi = DeltaPsiL( OG, L0, LP ) 
            DVareps = DeltaVarepsL( OG, L0, LP ) 
   
            vareps0 = varepsilon0_ecliptic( T ) 
            vareps = vareps0 + Dvareps  
          ELSE        
            T = ASVAL%T  
            M  = MODULO( ASVAL%M, 360.D0 )
            L0 = MODULO( ASVAL%L0, 360.D0 )
            eccen = ASVAL%eccen  
            DPsi = ASVAL%DPsi 
            vareps = ASVAL%vareps 
          ENDIF     
          C = SUN_CENTER( T, M )

          IF ( use_nutation ) THEN
             CALL SUN_COORDINATES_SUB0( lambda, beta, Delta, L0, M, eccen, C, DPsi )
          ELSE
             CALL SUN_COORDINATES_SUB0( lambda, beta, Delta, L0, M, eccen, C )
          END IF

          ! Apperent right ascension & declination
          CALL ECLIP2EQ( RA, DEC, lambda, beta, vareps ) 

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
        SUBROUTINE SUN_COORDINATES_SUB0( lambda, beta, Delta, L0, M, eccen, C, DPsi )
          IMPLICIT NONE 

          REAL (8), INTENT(OUT):: lambda, beta, Delta 
          REAL (8), INTENT(IN):: L0, M, C, eccen
          REAL (8), optional, INTENT(IN):: DPsi 
       
          REAL (8):: snu

          beta =  0.D0   
          IF ( .NOT. present(DPsi) ) THEN
            lambda = SUN_LON( C, L0 ) 
          ELSE  
            lambda = SUN_LON( C, L0, DPsi )   
          END IF    

          snu = SUN_NU( M, C ) 
          Delta = SUN_RADIUS(  snu, eccen )  ! Distance between centers 
                                              ! of the Earth and Sun
          RETURN  
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

        !  TABLE 47. B. page 341. Periodic terms for the latitude of
        !  the Moon (\Sum b). The unit is 1e-6 degree
        real(8), parameter, private :: LATARG(NPER, 4) = reshape((/0.d0, 0.d0, 0.d0, 1.d0, &
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
                                                                   2.d0, -2.d0, 0.d0, 1.d0/), (/NPER, 4/))

        real(8), parameter, private :: LAT_SIN_COEF(NPER) = (/5128122.d0, &
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
                                                              107.d0/)

        ! TABLE 47. A.
        !  Periodic terms for the longitude \Sum l and distance \Sum r
        !  of the Moon. The unit is 1e-6 for \Sum l and Kilometer for \Sum r
        real(8), parameter, private :: LONARG(NPER, 4) = reshape((/0.d0, 0.d0, 1.d0, 0.d0, &
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
                                                                   2.d0, 0.d0, -1.d0, -2.d0/), (/NPER, 4/))

        real(8), parameter, private :: LON_SIN_COEF(NPER) = (/6288774.d0, &
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
                                                              0.d0/)

        real(8), parameter, private :: LON_COS_COEF(NPER) = (/-20905355.d0, &
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
                                                              8752.d0/)
          integer, private, parameter :: LATMEU(NPER) = INT(ABS(LATARG(:,2))) 
          integer, private, parameter :: LONMEU(NPER) = INT(ABS(LONARG(:,2)))  

        PRIVATE:: CAL_EPMUL

        INTERFACE MOON_COORDINATES
            MODULE PROCEDURE MOON_COORDINATES_SUB0, MOON_COORDINATES_SUB1    
        END INTERFACE MOON_COORDINATES
      CONTAINS

        ! Coefficient muitplying sine and cosine arguments. page 338
        REAL(8) ELEMENTAL FUNCTION E_MUL_COEF( T ) RESULT( E )
           IMPLICIT NONE
           REAL (8), INTENT(IN) :: T
           E = 1.D0 + T*( -0.002516D0 + T*( -0.0000074D0 ))    
        END FUNCTION E_MUL_COEF        
 
        ! Page. 338 
        REAL(8) ELEMENTAL FUNCTION A1_DEG( T ) RESULT (A1)
           IMPLICIT NONE
           REAL(8), INTENT(IN):: T
           A1 = 119.75D0 + 131.849D0*T  
        END FUNCTION A1_DEG
 
        ! Page. 338
        REAL(8) ELEMENTAL FUNCTION A2_DEG( T ) RESULT (A2)
           IMPLICIT NONE
           REAL(8), INTENT(IN):: T
           A2 = 53.09D0 + 479264.290D0*T 
        END FUNCTION A2_DEG
 
        ! Page. 338
        REAL(8) ELEMENTAL FUNCTION A3_DEG( T ) RESULT (A3)
           IMPLICIT NONE
           REAL(8), INTENT(IN):: T
           A3 = 313.45D0 + 481266.484D0*T 
        END FUNCTION A3_DEG
         
       ! Coefficients for terms contains angle M. See description on
       ! Page 338   
       FUNCTION CAL_EPMUL( MEU, E ) RESULT( EPMUL )
         IMPLICIT NONE

         REAL (8):: EPMUL(NPER)
         REAL (8), INTENT(IN) ::  E
         INTEGER, INTENT(IN)  ::  MEU(:)

         INTEGER :: I, J

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
       SUBROUTINE MOON_COORDINATES_SUB1( RA, DEC, Delta, JD, ASVAL, NUTATION )
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
            IF ( abs(JD - ASVAL%JD) < 1.e-9 ) have_asval = .true.  
         ENDIF
         IF ( present(NUTATION) ) use_nutation = NUTATION 

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
   
           vareps0 = varepsilon0_ecliptic( T ) 
           DVareps = DeltaVarepsL( OG, L0, LP ) 
           vareps = vareps0 + Dvareps  
         ELSE        
           T = ASVAL%T  
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
         E  = E_MUL_COEF( T ) 
 
         CALL MOON_COORDINATES_SUB0( lambda, beta, Delta, LP, D, M, MP, F, E, A1, A2, A3 ) 

         ! PRINT*, "-------- RA & DEC --------" 
         ! Apperent right ascension & declination
         nutation_mul = 1.D0 
         IF ( .NOT. use_nutation ) nutation_mul = 0.D0  

         CALL ECLIP2EQ( RA, DEC, lambda + nutation_mul*DPsi, beta, vareps ) 
      
         ! Geometric right ascension & declination
         ! CALL ECLIP2EQ( RA, DEC, lambda, beta, vareps ) 
         ! PRINT*, "Geometrical RA  = ", MODULO(RA,360.D0)
         ! PRINT*, "Geometrical DEC = ", DEC
         RETURN  
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
       SUBROUTINE MOON_COORDINATES_SUB0( lambda, beta, Delta, LP, D, M, MP, F, E, A1, A2, A3 ) 
         IMPLICIT NONE

         REAL (8), INTENT(OUT):: lambda, beta, Delta
         REAL (8), INTENT(IN):: LP, D, M, MP, F
         REAL (8), INTENT(IN):: E, A1, A2, A3

         INTEGER:: I, J
         REAL (8):: suml, sumr, sumb
         REAL (8):: vec(4), argval(NPER), epmul(NPER)

         ! Convert from degree to radian 
         vec = DEG2RAD*(/ D, M, MP, F /) 

         ! 1. Compute mean longitude !
         argval = MATMUL( LONARG, vec ) 
         epmul = cal_epmul( LONMEU, E )  

         ! sine argument
         !  \sum l
         argval = epmul*LON_SIN_COEF*sin(argval) 
         suml = sum(argval) + additive_suml()  
  
         ! 2. Compute distance       !
         argval = MATMUL( LONARG, vec )  

         ! cosine argument
         !  \Sum r
         argval = epmul*LON_COS_COEF*cos(argval) 
         sumr = sum( argval )  

         ! 3. Compute mean latitude !
         argval = MATMUL( LATARG, vec ) 
         epmul = cal_epmul( LATMEU, E ) 

         ! sine argument
         !  \Sum b
         argval = epmul*LAT_SIN_COEF*sin(argval)  
         sumb = sum( argval ) + additive_sumb()     

         ! output   
         lambda = LP + suml/1.0D6  ! Longitude in Degree
         beta = sumb/1.D6         ! Lattitude in Degree
         Delta = 385000.56D0 + sumr/1.0D3  ! Distance in km  

       CONTAINS
         
         REAL(8) FUNCTION additive_suml( ) RESULT( asuml )
           IMPLICIT NONE          
           asuml =  3958.D0*sin( DEG2RAD*A1 ) &
                  + 1962.D0*sin( DEG2RAD*(LP - F) ) &
                  +  318.D0*sin( DEG2RAD*A2 ) 
         END FUNCTION additive_suml        

         REAL(8) FUNCTION additive_sumb( ) RESULT( asumb )
           IMPLICIT NONE  
           asumb = -2235.D0*sin( DEG2RAD*LP ) &
                   + 382.D0*sin( DEG2RAD*A3  ) &
                   + 175.D0*sin( DEG2RAD*(A1 - F) ) &
                   + 175.D0*sin( DEG2RAD*(A1 + F) ) &
                   + 127.D0*sin( DEG2RAD*(LP - MP) ) &
                   - 115.D0*sin( DEG2RAD*(LP + MP) ) 

        END FUNCTION additive_sumb         

       END SUBROUTINE MOON_COORDINATES_SUB0

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

           RETURN      
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
        SUBROUTINE HEAVENLY_OBJS_COORDS_JM( MOON_POS, SUN_POS, JD, IERR )
           IMPLICIT NONE

           REAL (8), INTENT(IN):: JD
           REAL (8), INTENT(OUT):: MOON_POS(3), SUN_POS(3)
           INTEGER:: IERR

           ASTRO_VALUES = COMP_ASTROVAL(  JD  ) 

           CALL MOON_COORDINATES(  MOON_POS(1), MOON_POS(2), MOON_POS(3), JD, ASTRO_VALUES, INCLUDE_NUTATION )
           CALL SUN_COORDINATES( SUN_POS(1), SUN_POS(2), SUN_POS(3), JD, ASTRO_VALUES, INCLUDE_NUTATION )

           IERR = 0 

           RETURN             
        END SUBROUTINE HEAVENLY_OBJS_COORDS_JM  

        FUNCTION GMST_DEG_FN( JDE ) RESULT( GMST )
           IMPLICIT NONE 

           REAL (8):: gmst
           REAL (8), INTENT(IN):: JDE
             
           gmst = GMST_DEG( JDE, INCLUDE_NUTATION, ASTRO_VALUES )  
  
        END FUNCTION GMST_DEG_FN        

      END MODULE MOON_SUN_COORS        
