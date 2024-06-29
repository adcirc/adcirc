!
! Fortran module for computing the full luna-solar equilibrium tides of
! different orders of approximation, 
! 
! Ref:
!   J,-M. Hervouet, Free surface flows: modelling with the finite
!      element method. John Wiley & Sons, Ltd, 2007
!          
!   Z. Kowalik and J. L. Luick, Modern theory and practice of tide
!      analysis and tidal power, Austides consulting, 2019 
!      https://austides.com/downloads/
!   
! DW/AC, 2024
!      
      MODULE TIPFUNCMOD

#ifndef TIPSTANDALONE
       use global, only: DEBUG, INFO, ECHO, WARNING, ERROR,
     &   scratchMessage, logMessage, screenMessage, allMessage,
     &   setMessageSource, openFileForRead, unsetMessageSource,
     &   scratchFormat

        USE ADC_CONSTANTS, only: hour2day, min2day, sec2day, 
     &       DEG2RAD, RAD2DEG 
        USE SIZES, only: MYPROC
#endif

#ifdef DATETIME
        use datetime_module, only: strptime
#endif

        USE ASTROFORMOD
        USE MOON_SUN_COORS, only: HEAVENLY_OBJS_COORDS_JM, 
     &                                 SET_NUTATION, GMST_DEG_FN
        USE ephemri_module, only: HEAVENLY_OBJS_COORDS_FROM_TABLE


        INTERFACE HEAVENLY_OBJS_COORDS
            MODULE PROCEDURE  HEAVENLY_OBJS_COORDS_JM, 
     &             HEAVENLY_OBJS_COORDS_FROM_TABLE            
        END INTERFACE HEAVENLY_OBJS_COORDS

        !
        !  UseFullTIPFormula = T/F (default F)
        !  TIPOrder = 2 (P2), 3(p3), >= 4(exact)  (default 2)
        !  TIPStartDate = 'BaseDate'/'YYYY-MM-DD HH:mm:SS'
        !  MoonSunPostionMethod = 'JM'/'External' (defaul 'JM') 
        !           'JM' = Jean Meeus's approach 
        !           'External' = User supplied Moon's, Sun's RA, DEC and
        !                        distances 
        !  MoonSunCoordFile = filename of an external file storting the
        !                        coordinates of the Moon and the Sun         
        !  IncludeNutation = T/F, Include/exclude the nutation in JM formula
        !  k2value = Love number
        !  h2value = Love number
        !
        LOGICAL:: UseFullTIPFormula = .FALSE. 
        INTEGER:: TIPOrder = 2
        CHARACTER (LEN=24):: TIPStartDate = 'Basedate' 
        CHARACTER (LEN=24):: MoonSunPositionComputeMethod = 'JM' 
        CHARACTER (LEN=80):: MoonSunCoordFile = 'none'
        LOGICAL:: UniformResMoonSunTimeData = .FALSE.  
        LOGICAL:: IncludeNutation = .TRUE.  
        REAL (8):: k2value = 0.302D0 
        REAL (8):: h2value = 0.609D0    
 
        REAL (8), private:: JDE_BEG,  JDE_CURRENT 
        REAL (8), private, dimension(:), allocatable:: AH, cosZA, workarr, workarr1

        REAL (8), private, parameter:: massratio(2) = (/
     &                       MassRatioMoonEarth, MassRatioSunEarth  /) ; 

        REAL (8), private, parameter:: significant_radiusearth(2) = (/
     &                       EarthRadiuskm(1), EarthRadiusAu(1) /) ;
        REAL (8), private, parameter:: exponent_radiusearth(2) = (/
     &                       EarthRadiuskm(2), EarthRadiusAu(2) /) ; 

#ifdef TIPSTANDALONE
        REAL(8), PARAMETER :: hour2sec = 3600.0D0
        REAL(8), PARAMETER :: sec2hour = 1.0D0 / hour2sec
        REAL(8), PARAMETER :: day2hour = 24.0D0
        REAL(8), PARAMETER :: hour2day = 1.0D0 / day2hour
        REAL(8), PARAMETER :: day2sec = day2hour * hour2sec
        REAL(8), PARAMETER :: sec2day = 1.0D0/day2sec

        REAL (8), parameter:: hour2min = 60.D0 
        REAL (8), parameter:: min2hour = 1.D0/hour2min
        REAL (8), parameter:: min2day = 1.D0*min2hour*hour2day ;
#endif        

        PRIVATE:: AINTPOWER, INIT_TIP

       CONTAINS
         
        !
        ! 
        SUBROUTINE COMP_FULL_TIP( TIPVAL, TimeLoc, SLAM, SFEA )  
           IMPLICIT NONE

           REAL (8), INTENT(INOUT):: TIPVAL(:) 
           REAL (8), INTENT(IN):: TimeLoc      ! Time in seconds since TIPStartDate !
           REAL (8), INTENT(IN), DIMENSION(:):: SLAM, SFEA ! lon-lat in radian

           ! local ! 
           INTEGER:: I, IOBJ, iexp, NPP

           REAL (8):: JDELoc, tocgmst
           REAL (8):: MoonSunCoor(3,2), RA, DEC, DELTA
           REAL (8):: radius_div_Delta, KP, C0
           INTEGER:: IERR

           LOGICAL, SAVE:: first = .true.

#ifndef TIPSTANDALONE           
          call setMessageSource("comp_full_tip")
#if defined(ALL_TRACE)
          call allMessage(DEBUG,"Enter.")
#endif
#endif

           NPP = SIZE( TIPVAL ) ; 
           IF ( first ) THEN   
              CALL INIT_TIP( NPP ) ;
              first = .false. ;
           END IF     

           ! Julian day
           JDELoc = JDE_BEG + TimeLoc*sec2day ;            

           ! Compute the postions of the Moon and Sun !
           IF ( TRIM(MoonSunPositionComputeMethod) == 'JM' ) THEN
              ! use algorithms in Jean Meeus's book !
              CALL HEAVENLY_OBJS_COORDS( MoonSunCoor(:,1), 
     &                     MoonSunCoor(:,2), JDELoc, IERR ) ; 
           ELSE
              ! interpolate from an external look up table !
              CALL HEAVENLY_OBJS_COORDS( MoonSunCoor, JDELoc, 
     &          MoonSunCoordFile, IERR, UniformResMoonSunTimeData ) ;
           END IF   
        
           ! Convert RA, DEC from deg --> rad
           MoonSunCoor(1:2,:) = MoonSunCoor(1:2,:)*DEG2RAD ;

           ! Sidereal time at Greenwich ! 
           tocgmst = GMST_DEG_FN( JDELoc )  ;
          
           TIPVAL = 0.D0 ;
           DO IOBJ = 1, 2 
              ! Hour angle !
              AH = DEG2RAD*tocgmst + SLAM - MoonSunCoor(1,IOBJ) ;
              ! AH = modulo( tocgmst + SLAM*RAD2DEG, 360.D0 )*DEG2RAD -  MoonSunCoor(1,IOBJ) ;
              RA = MoonSunCoor(1,IOBJ) ; ! right ascendsion
              DEC = MoonSunCoor(2,IOBJ) ; ! declination
              DELTA = MoonSunCoor(3,IOBJ) ; ! distance (km for the Moon,
                                            ! AU for the sun

              ! Ratio of the radius of the earth and the Moon/Sun                              
              radius_div_delta = significant_radiusearth(IOBJ)/DELTA
              radius_div_delta = radius_div_delta*exponent_radiusearth(IOBJ) 

              !c cosZA = cos( Z )            
              cosZA = sin(SFEA)*sin(DEC) ;
              cosZA = cosZA + cos(SFEA)*cos(DEC)*cos(AH) ;   
              
              SELECT CASE( TipOrder )
              CASE (2,3)
                ! P2  term
                workarr = 3.D0*cosZA*cosZA ; 
                workarr = workarr - 1.D0   ;
                workarr = 0.5D0*workarr ;
 
                KP = AINTPOWER( radius_div_delta, 3 ) ;
                KP = KP*EarthRadiusm(1)*EarthRadiusm(2) ; 
                KP = KP*massratio(IOBJ) ;

                TIPVAL = TIPVAL + KP*workarr ; 

                ! if include P3 term 
                IF ( TipOrder == 3 ) THEN
                   workarr = 5.D0*cosZA*cosZA*cosZA ;
                   workarr = workarr - 3.D0*cosZA ; 
                   workarr = 0.5D0*workarr ;

                   KP = AINTPOWER( radius_div_delta, 4 ) ; 
                   KP = KP*EarthRadiusm(1)*EarthRadiusm(2) ;
                   KP = KP*massratio(IOBJ) ; 
                   ! Offset values is zero for this term c!

                   TIPVAL = TIPVAL + KP*workarr ; 
                END IF   
              CASE (22)    
                ! P2 - Only dirunal + semi-dirunali
                ! Long period
                ! workarr1 = sin(SFEA) ;
                ! workarr = 1.D0 - 3.D0*workarr1*workarr1 ; 
                ! workarr = workarr*(1.D0 - 3.D0*sin(DEC)*sin(DEC))
                ! workarr = (1.D0/4.D0)*workarr ; 
                ! Dirunal
                ! workarr = workarr + (3.D0/4.D0)*sin(2.D0*SFEA)*sin(2.D0*DEC)*cos(AH) ;

                ! Dirunal
                workarr = (3.D0/4.D0)*sin(2.D0*SFEA)*sin(2.D0*DEC)*cos(AH) ;

                ! Semi-Diurnal
                workarr1 = cos( SFEA ) ; 
                workarr = workarr + 
     &                  (3.D0/4.D0)*workarr1*workarr1*cos(DEC)*cos(DEC)*cos(2.D0*AH) ;

 
                KP = AINTPOWER( radius_div_delta, 3 ) ;
                KP = KP*EarthRadiusm(1)*EarthRadiusm(2) ; 
                KP = KP*massratio(IOBJ) ; 

                TIPVAL = TIPVAL + KP*workarr ; 
              CASE DEFAULT
                ! Exact
                workarr = 1.D0 + AINTPOWER( radius_div_delta, 2 ) ;
                workarr = workarr - 2.D0*radius_div_delta*cosZA  ;

                workarr = sqrt( workarr ) ; 
                workarr = 1.D0/workarr ; 

                workarr = workarr - radius_div_delta*cosZA ; 

                ! Offset value, use Proudman's approach, see Kowalik,
                ! page 14
                C0 = 1.D0 + AINTPOWER(radius_div_delta, 2) ;
                C0 =  SQRT(C0 + 2.D0*radius_div_delta) - 
     &                SQRT(C0 - 2.D0*radius_div_delta)  ;
                C0 = -0.5D0*C0/radius_div_delta ;          

                ! NOTE: 
                ! The constant of integration C0 equals to
                !
                ! C0 = -1/2*(1/e)*(sqrt(1 + e**2  + 2*e) - sqrt(1 + e**2 - 2*e))
                ! e =  a/r = radius_div_delta
                ! 
                ! For radius_div_delta (a/r) << 1, it can be shown
                ! through the use of the Taylor series that
                !
                !    sqrt(1 + (a/r)^2 + 2*(a/r)) \sim  1 + (a/r) + High order terms 
                !    sqrt(1 + (a/r)^2 + 2*(a/r)) \sim  1 - (a/r) + High order terms
                !
                ! As a result, a good approximate value of the offset is 
                ! 
                !    C0 \sim -1 
                !             
                workarr = workarr + C0  ;

                KP = AINTPOWER( radius_div_delta, 1 ) ;
                KP = KP*EarthRadiusm(1)*EarthRadiusm(2) ;
                KP = KP*massratio(IOBJ) ; 

                TIPVAL = TIPVAL + KP*workarr ;  
              END SELECT 
           END DO

           TIPVAL = (1.D0 + k2value - h2value)*TIPVAL ;

#ifndef TIPSTANDALONE
#if defined(ALL_TRACE)
           call allMessage(DEBUG,"Return.")
#endif
           call unsetMessageSource()
#endif

           RETURN ;
         CONTAINS

           SUBROUTINE check_err( IERR  )
             IMPLICIT NONE
             
             INTEGER:: IERR 

             IF ( IERR > 1 ) THEN
#ifndef TIPSTANDALONE      
               SELECT CASE( IERR )
               CASE (1)      
                 CALL  screenMessage( ERROR, "Error in the calculation"//
     &             " the postion of the Moon and Sun. ADCIRC is not "//
     &              " built with the -DADCNETCDF flag." ) ;
                CASE (2)
                  CALL  screenMessage( ERROR, "Error in the calculation"//
     &              " the postion of the Moon and Sun. Date not within"//
     &               "  database." ) ;
               END SELECT
#else
               SELECT CASE( IERR )
               CASE (1)
                 WRITE(*,*) "Erorr in the calculation of"//
     &             " the position of the Moon and Sun. ADCIRC is not "//
     &             " built with the -DADCNETCDF flag."  ;
               CASE (2)   
                 WRITE(*,*) "Error in the calculation of"//
     &             " the postion of the Moon and Sun. Date not within"//
     &              "  database."  ;
               END SELECT 
             STOP ;             
#endif            
             END IF
           END SUBROUTINE check_err

        END SUBROUTINE COMP_FULL_TIP 
        !        

        !
        SUBROUTINE SET_TIP_BASEDATE( base_date )
           IMPLICIT NONE

           CHARACTER (LEN=*):: base_date
 
           CALL StringUpper( TIPStartDate ) ; 

           IF ( trim(TIPStartDate) == 'BASEDATE' ) THEN
              TIPStartDate =trim( base_date)    
           END IF

           RETURN ;     
        END SUBROUTINE SET_TIP_BASEDATE

        !
        SUBROUTINE StringUpper( string ) 
           IMPLICIT NONE     
 
           CHARACTER (LEN=*):: string
 
           INTEGER:: I, asciinum
 
           DO I = 1, len(string)
              asciinum = iachar(string(I:I)) ; 
 
              SELECT case(asciinum)
              CASE (97:122)
                string(I:I) = char( asciinum - 32 ) ; 
              END SELECT       
           END DO
 
        END SUBROUTINE StringUpper

        ! A^{N}
        ELEMENTAL FUNCTION AINTPOWER( A, N ) RESULT ( AP )
           IMPLICIT NONE 

           REAL (8):: AP
           REAL (8), INTENT(IN):: A
           INTEGER, INTENT(IN):: N 

           INTEGER:: I

           AP = 1.D0 ; 
           DO I = 1, N
            AP = AP*A ; 
           END DO

        END FUNCTION AINTPOWER

        !
        SUBROUTINE INIT_TIP( NP ) 
          IMPLICIT NONE

          INTEGER, INTENT(IN):: NP

          REAL (8):: DDD
          INTEGER:: YYYY, MM, DD, HH, MMM, SS

          INTEGER:: ii, cpos(2), ivec(3)
          CHARACTER (LEN=80)::  tmparr
          CHARACTER:: delimter(2) = (/ '-', ':' /) ; 

          ! Extract date from base_date
          tmparr = adjustl(TIPStartDate) ;
 
          CALL extractvalues( YYYY, MM, DD, delimter(1) ) ;  
          CALL extractvalues(  HH, MMM, SS, delimter(2) ) ; 

          DDD = DBLE(DD) + DBLE(HH)*hour2day + 
     &       DBLE(MMM)*min2day + DBLE(SS)*sec2day ;
           
          ! Get  correspond Julian days !
          JDE_BEG = JULIANDAY( DDD, MM, YYYY )  ;
          JDE_CURRENT = JDE_BEG ; 

          CALL StringUpper( MoonSunPositionComputeMethod ) ; 
          IF ( TRIM(ADJUSTL(MoonSunPositionComputeMethod)) == 'JM' ) THEN 
            CALL SET_NUTATION(  IncludeNutation  ) ;
          END IF

          CALL ALLOCATEWORKARR( AH, NP ) ; 
          CALL ALLOCATEWORKARR( cosZA, NP )
          CALL ALLOCATEWORKARR( workarr, NP ) ; 

          IF ( TIPOrder > 3 ) THEN
            CALL ALLOCATEWORKARR( workarr1, NP ) ; 
          END IF

          RETURN 
        CONTAINS

           SUBROUTINE ALLOCATEWORKARR( arr, N ) 
              IMPLICIT NONE

              INTEGER:: N
              REAL (8), ALLOCATABLE:: arr(:)

              IF ( ALLOCATED(arr) ) THEN
                DEALLOCATE(arr) ;
              END IF
              ALLOCATE( arr(1:N) ) ; 
              arr = 0.D0 ; 

              RETURN ;
           END SUBROUTINE ALLOCATEWORKARR  

           SUBROUTINE extractvalues( arg1, arg2, arg3, delimc )
             IMPLICIT NONE

             CHARACTER:: delimc
             INTEGER:: arg1, arg2, arg3 

             INTEGER:: ii,  cpos(3), ivec(3) 

             cpos(1) = 1 ;
   
             DO ii = 1, 3
               cpos(2) = index( tmparr, delimc ) ;
              
               if ( cpos(2) == 0 )  cpos(2) = index( tmparr, ' ' ) ; 
  
               if ( cpos(2) > 0  ) then
                  read( tmparr(cpos(1):cpos(2)-1), * ) ivec(ii)  ;
                  
                  tmparr = adjustl(tmparr(cpos(2)+1:) ) ;
               else
                  read( tmparr(cpos(1):), * ) ivec(ii)  ;
               end if
             END DO    
             
             arg1 = ivec(1) ;
             arg2 = ivec(2) ;
             arg3 = ivec(3) ;

             RETURN 
           END SUBROUTINE extractvalues        
           
        END SUBROUTINE INIT_TIP 
        !

        !
        SUBROUTINE INIT_TIP_FTCE( L_N, NP, SFEA ) 
          IMPLICIT NONE

          INTEGER:: NP
          REAL (8):: L_N(:,:), SFEA(:)

          INTEGER:: I 

          L_N = 0.D0 ;

          DO I = 1,NP
C           LONG-TERM SPECIES
            L_N(0,I) = 1.5d0*COS(SFEA(I))*COS(SFEA(I)) - 1d0
C           DIURNAL SPECIES
            L_N(1,I) = SIN(2d0*SFEA(I))
C           SEMIDIURNAL SPECIES
            L_N(2,I) = COS(SFEA(I))*COS(SFEA(I))
          END DO
          RETURN      
        END SUBROUTINE INIT_TIP_FTCE

      END MODULE TIPFUNCMOD        

