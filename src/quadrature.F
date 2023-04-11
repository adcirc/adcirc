      MODULE QUADRATURETRI
      USE SIZES
 
      IMPLICIT NONE

      INTEGER, private:: pcub, ncubs
      REAL(8), private, pointer:: cub2d(:,:)
      REAL(8), private, allocatable:: brscub2d(:,:)
      REAL(8), private, allocatable:: workvec(:)

      TYPE MATARR
          INTEGER :: mdim
          REAL(8), POINTER :: ARRVAL(:,:)
      END TYPE MATARR
    
      TYPE VECARR
         INTEGER :: mdim
         REAL(8), POINTER :: VECVAL(:)
      END TYPE VECARR

      CONTAINS

!C=================== High-level subroutines ====================
      SUBROUTINE AllocMatArr( msarr, iidim )
        IMPLICIT NONE
       
        INTEGER      :: iidim
        TYPE(MATARR) :: msarr
 
        msarr%mdim = iidim
        IF ( associated(msarr%arrval) ) THEN
          deallocate(msarr%arrval)
          nullify(msarr%arrval)
        END IF
        
        ALLOCATE(msarr%arrval(iidim,iidim))

        RETURN        
      END SUBROUTINE AllocMatArr

      SUBROUTINE AllocVecArr( msvec, iidim )
        IMPLICIT NONE
       
        INTEGER      :: iidim
        TYPE(VECARR) :: msvec
 
        msvec%mdim = iidim
        IF ( associated(msvec%vecval) ) THEN
           deallocate(msvec%vecval) 
           nullify(msvec%vecval)
        END IF
        
        ALLOCATE(msvec%vecval(iidim))

        RETURN
      END SUBROUTINE AllocVecArr

      SUBROUTINE GetDefaultCub2D( p )
        IMPLICIT NONE
        
        INTEGER:: p

        IF ( associated(cub2d) ) THEN
          deallocate(cub2d) ;
        END IF
        IF ( allocated(brscub2d) ) THEN
          deallocate(brscub2d) ;
        END IF
        IF ( allocated(workvec) ) THEN
          deallocate(workvec) ;
        END IF

        pcub = p ;
        CALL GetQuadrature2DTri( cub2d, pcub ) ;
        
        ncubs = ubound( cub2d, 1 ) ; 
        ALLOCATE( brscub2d(ncubs,3) ) ;

        CALL rsbarycentriccoords( brscub2d, 
     &                 cub2d(:,1), cub2d(:,2), ncubs ) ;

        allocate( workvec(ncubs) ) ; 

        RETURN ;      
      END SUBROUTINE GetDefaultCub2D

      SUBROUTINE FetchDefaultNcub2d( nn )
        IMPLICIT NONE
        
        INTEGER:: nn

        nn = ncubs ;

        RETURN ;
      END SUBROUTINE FetchDefaultNcub2d

      SUBROUTINE InterpDefaultCub2D( fxy, fv )
        IMPLICIT NONE

        REAL(8), dimension(:):: fxy, fv

        CALL barycentricinterp( fxy, fv, brscub2d, ncubs ) ;

        RETURN ;  
      END SUBROUTINE InterpDefaultCub2D
     
      !c
      !c Compute elemental:
      !c
      !c S = (1/(A/2)) \int f_{h} phi_{i} phi_{j}
      !c
      SUBROUTINE CompElmMsfh( S, fv )
         IMPLICIT NONE

         INTEGER:: II, JJ
         REAL(8):: S(3,3), fv(:)

         workvec = 0.D0 ;
         CALL InterpDefaultCub2D( workvec, fv ) ;
         
         DO II = 1, 3
           DO JJ = 1, 3
              S(II,JJ) = sum(workvec*cub2d(:,3)*
     &                brscub2d(:,II)*brscub2d(:,JJ)) ;
           END DO
         END DO
               
         RETURN ;
      END SUBROUTINE CompElmMsfh
C============================================================

      !
      ! Linear interpolation on triangle
      !
      ! given \lambda | (r,s)  return  f(xy)
      SUBROUTINE barycentricinterp( fxy, fv, brs, nrs )
        IMPLICIT NONE

        INTEGER:: nrs
        REAL(8), dimension(:):: fxy, fv
        REAL(8):: brs(:,:)


        fxy(1:nrs) = fv(1)*brs(1:nrs,1) + 
     &          fv(2)*brs(1:nrs,2) + fv(3)*brs(1:nrs,3) ;

        RETURN ;     
      END SUBROUTINE barycentricinterp

      !c RS to lambda
      SUBROUTINE rsbarycentriccoords( brs, r, s, nrs )
        IMPLICIT NONE
         
        INTEGER:: nrs
        REAL(8):: brs(:,:), r(:), s(:)

        INTEGER:: I

        brs(1:nrs,1) = -0.5D0*(r(1:nrs) + s(1:nrs)) ;
        brs(1:nrs,2) =  0.5D0*(r(1:nrs) + 1.D0) ;
        brs(1:nrs,3) =  0.5D0*(s(1:nrs) + 1.D0) ;
      
        RETURN ;
      END SUBROUTINE rsbarycentriccoords

      !c
      !c ref element:   -1 < x,y < 1, x + y = 0 
      SUBROUTINE GetQuadrature2DTri( cub2d, pcub )

      IMPLICIT NONE

      INTEGER:: pcub
      REAL(8), pointer:: cub2d(:,:)


      INTEGER:: p1D, sk, ii, jj, ncub
      REAL(8), allocatable:: x1D(:,:), x2D(:,:)

      SELECT CASE( pcub ) 
      case( 1 )
         allocate(cub2d(1,3)) ; 
         cub2d( 1, :) = (/ -3.333333333333330d-01,  -3.333333333333330d-01, 2.000000000000000d+00 /) ; 
      case( 2 )
         allocate(cub2d(3,3)) ; 
         cub2d( 1, :) = (/ -6.666666666666670d-01,  -6.666666666666670d-01, 6.666666666666670d-01 /) ; 
         cub2d( 2, :) = (/ 3.333333333333330d-01,  -6.666666666666670d-01, 6.666666666666670d-01 /) ; 
         cub2d( 3, :) = (/ -6.666666666666670d-01,  3.333333333333330d-01, 6.666666666666670d-01 /) ; 
      case( 3 )
         allocate(cub2d(6,3)) ; 
         cub2d( 1, :) = (/ -8.168475729804580d-01,  -8.168475729804580d-01, 2.199034873106440d-01 /) ; 
         cub2d( 2, :) = (/ 6.336951459609170d-01,  -8.168475729804590d-01, 2.199034873106440d-01 /) ; 
         cub2d( 3, :) = (/ -8.168475729804590d-01,  6.336951459609170d-01, 2.199034873106440d-01 /) ; 
         cub2d( 4, :) = (/ -1.081030181680700d-01,  -1.081030181680700d-01, 4.467631793560230d-01 /) ; 
         cub2d( 5, :) = (/ -7.837939636638600d-01,  -1.081030181680700d-01, 4.467631793560230d-01 /) ; 
         cub2d( 6, :) = (/ -1.081030181680700d-01,  -7.837939636638600d-01, 4.467631793560230d-01 /) ; 
      case( 4 )
         allocate(cub2d(6,3)) ; 
         cub2d( 1, :) = (/ -8.168475729804580d-01,  -8.168475729804580d-01, 2.199034873106440d-01 /) ; 
         cub2d( 2, :) = (/ 6.336951459609170d-01,  -8.168475729804590d-01, 2.199034873106440d-01 /) ; 
         cub2d( 3, :) = (/ -8.168475729804590d-01,  6.336951459609170d-01, 2.199034873106440d-01 /) ; 
         cub2d( 4, :) = (/ -1.081030181680700d-01,  -1.081030181680700d-01, 4.467631793560230d-01 /) ; 
         cub2d( 5, :) = (/ -7.837939636638600d-01,  -1.081030181680700d-01, 4.467631793560230d-01 /) ; 
         cub2d( 6, :) = (/ -1.081030181680700d-01,  -7.837939636638600d-01, 4.467631793560230d-01 /) ; 
      case( 5 )
         allocate(cub2d(7,3)) ; 
         cub2d( 1, :) = (/ -3.333333333333330d-01,  -3.333333333333330d-01, 4.500000000000000d-01 /) ; 
         cub2d( 2, :) = (/ -5.971587178977000d-02,  -5.971587178977000d-02, 2.647883055770120d-01 /) ; 
         cub2d( 3, :) = (/ -8.805682564204600d-01,  -5.971587178977000d-02, 2.647883055770120d-01 /) ; 
         cub2d( 4, :) = (/ -5.971587178977000d-02,  -8.805682564204600d-01, 2.647883055770120d-01 /) ; 
         cub2d( 5, :) = (/ -7.974269853530870d-01,  -7.974269853530870d-01, 2.518783610896540d-01 /) ; 
         cub2d( 6, :) = (/ 5.948539707061750d-01,  -7.974269853530870d-01, 2.518783610896540d-01 /) ; 
         cub2d( 7, :) = (/ -7.974269853530870d-01,  5.948539707061750d-01, 2.518783610896540d-01 /) ; 
      case( 6 )
         allocate(cub2d(12,3)) ; 
         cub2d( 1, :) = (/ -5.014265096581790d-01,  -5.014265096581790d-01, 2.335725514527590d-01 /) ; 
         cub2d( 2, :) = (/ 2.853019316358000d-03,  -5.014265096581790d-01, 2.335725514527590d-01 /) ; 
         cub2d( 3, :) = (/ -5.014265096581790d-01,  2.853019316358000d-03, 2.335725514527590d-01 /) ; 
         cub2d( 4, :) = (/ -8.738219710169960d-01,  -8.738219710169960d-01, 1.016898127404140d-01 /) ; 
         cub2d( 5, :) = (/ 7.476439420339910d-01,  -8.738219710169960d-01, 1.016898127404140d-01 /) ; 
         cub2d( 6, :) = (/ -8.738219710169960d-01,  7.476439420339910d-01, 1.016898127404140d-01 /) ; 
         cub2d( 7, :) = (/ -3.792950979324310d-01,  -8.937099003103660d-01, 1.657021512367470d-01 /) ; 
         cub2d( 8, :) = (/ -8.937099003103660d-01,  -3.792950979324310d-01, 1.657021512367470d-01 /) ; 
         cub2d( 9, :) = (/ 2.730049982427970d-01,  -8.937099003103660d-01, 1.657021512367470d-01 /) ; 
         cub2d( 10, :) = (/ -8.937099003103660d-01,  2.730049982427970d-01, 1.657021512367470d-01 /) ; 
         cub2d( 11, :) = (/ 2.730049982427970d-01,  -3.792950979324310d-01, 1.657021512367470d-01 /) ; 
         cub2d( 12, :) = (/ -3.792950979324310d-01,  2.730049982427970d-01, 1.657021512367470d-01 /) ; 
      CASE DEFAULT
         PRINT*, "Error: Quadrature for pcub = ", pcub, " has not yet been implemented" ;        
      END SELECT

      END SUBROUTINE GetQuadrature2DTri

      END MODULE QUADRATURETRI
