      MODULE TIPMOD
         USE MESH, ONLY : NE, SLAM0, SFEA0, 
     &      NP, SLAM, SFEA

        REAL (8), ALLOCATABLE, DIMENSION(:,:):: L_N

        !  UseFullTIPFormula = T/F (default F)
        !  TIPOrder = 2 (P2), 3(p3), >= 4(exact)  (default 2)
        !  TIPStartDate = 'BaseDate'/'YYYY-MM-DD HH:SS'
        !  MoonSunPostionMethod = 'JM'/'External' (defaul 'JM') 
        !           'JM' = Jean Meeus's approach 
        !           'External' = User supplied Moon's, Sun's RA, DEC and
        !                        distances 
        !  MoonSunCoordFile = filename of an external file storting the
        !                        coordinates of the Moon and the Sun         
        !  IncludeNutation = T/F, Include/exclude the nutation in JM formula
        !  k2value = Love number
        !  h2value = Love number
        LOGICAL:: UseFullTIPFormula = .FALSE. 
        INTEGER:: TIPOrder = 2 
        CHARACTER (LEN=24):: TIPStartDate = 'Basedate' 
        CHARACTER (LEN=24):: MoonSunPositionComputeMethod = 'JM' 
        CHARACTER (LEN=80):: MoonSunCoordFile = 'none' 
        LOGICAL:: IncludeNutation = .TRUE. 
        REAL (8):: MSWTINC = 3
        REAL (8):: k2value = 0.302D0 
        REAL (8):: h2value = 0.609D0    
  
      CONTAINS
         
        SUBROUTINE INIT_TIP_FTCE( ) 
          IMPLICIT NONE

          IF ( ALLOCATED(L_N) ) DEALLOCATE(L_N) ;
          ALLOCATE( L_N(0:2,NP) ) ; L_N = 0.D0 ;

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

      END MODULE TIPMOD        

