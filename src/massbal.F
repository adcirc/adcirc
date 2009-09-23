

      SUBROUTINE COMPUTE_RESELEM

      USE GLOBAL, ONLY:
     &  NE,
     &  NM,
     &  RESELEM,
     &  PI,DEG2RAD,
     &  X,Y,SLAM0,SFEA0,ICS

C    &  SLAM,SFEA

      USE SIZES, ONLY:
     &  MYPROC, SZ

      IMPLICIT NONE

      INTRINSIC             :: REAL

      CHARACTER             :: JunkC

      INTEGER, ALLOCATABLE  :: Connectivity(:)
      INTEGER               :: GlobalNode1
      INTEGER               :: GlobalNode2
      INTEGER               :: GlobalNode3
      INTEGER               :: I
      INTEGER               :: J
      INTEGER               :: JunkI
      INTEGER               :: NumNodesGlobal

      REAL(SZ), ALLOCATABLE :: GlobalX(:),GL_X(:)
      REAL(SZ), ALLOCATABLE :: GlobalY(:),GL_Y(:)
      REAL(SZ)              :: Tol = 0.00001d0      

      OPEN(14,FILE='fort.14')
      READ(14,'(A)') JunkC
      READ(14,*) JunkI, NumNodesGlobal
      ALLOCATE(GlobalX(1:NumNodesGlobal),GL_X(1:NumNodesGlobal))
      ALLOCATE(GlobalY(1:NumNodesGlobal),GL_Y(1:NumNodesGlobal))
      DO I=1,NumNodesGlobal
         READ(14,*) JunkI, GlobalX(I), GlobalY(I)
   
C.RJW    need to convert global node lat and lon into CPP projected units for comparison
C...     with local node location x and y values if ICS =2 
        IF (ICS .eq. 1 ) then

           GL_X(I) = GlobalX(I)
           GL_Y(I) = GlobalY(I)
          
        ELSE  !(ICS=2)

            GlobalX(I)=DEG2RAD*GlobalX(I)
            GlobalY(I)=DEG2RAD*GlobalY(I)

         CALL CPP(GL_X(I),GL_Y(I),GlobalX(I),GlobalY(I),SLAM0,SFEA0)

        ENDIF

      ENDDO
      CLOSE(14)
       

      OPEN(80,FILE='partmesh.txt')
      ALLOCATE(Connectivity(1:NumNodesGlobal))
      DO I=1,NumNodesGlobal
         READ(80,*) Connectivity(I)
         Connectivity(I) = Connectivity(I) - 1  ! processors start with 000
      ENDDO
      CLOSE(80)
        
      ALLOCATE(RESELEM(1:NE))
      DO I=1,NE
C...RJW initialize variable
         GlobalNode1=-999
         GlobalNode2=-999
         GlobalNode3=-999    

         DO J=1,NumNodesGlobal

            IF((ABS(GL_X(J)-X(NM(I,1))).LT.Tol).AND.
     &         (ABS(GL_Y(J)-Y(NM(I,1))).LT.Tol))THEN
                GlobalNode1 = J
            ENDIF
            IF((ABS(GL_X(J)-X(NM(I,2))).LT.Tol).AND.
     &         (ABS(GL_Y(J)-Y(NM(I,2))).LT.Tol))THEN
                GlobalNode2 = J
            ENDIF
            IF((ABS(GL_X(J)-X(NM(I,3))).LT.Tol).AND.
     &         (ABS(GL_Y(J)-Y(NM(I,3))).LT.Tol))THEN
                GlobalNode3 = J
            ENDIF

         ENDDO

         IF((Connectivity(GlobalNode1).EQ.MYPROC).AND.
     &      (Connectivity(GlobalNode2).EQ.MYPROC).AND.
     &      (Connectivity(GlobalNode3).EQ.MYPROC))THEN
             RESELEM(I) = 1
         ELSEIF(
     &      (Connectivity(GlobalNode1).EQ.MYPROC).AND.
     &      (Connectivity(GlobalNode2).EQ.MYPROC).AND.
     &      (Connectivity(GlobalNode3).NE.MYPROC))THEN
             RESELEM(I) = 2
         ELSEIF(
     &      (Connectivity(GlobalNode1).EQ.MYPROC).AND.
     &      (Connectivity(GlobalNode2).NE.MYPROC).AND.
     &      (Connectivity(GlobalNode3).EQ.MYPROC))THEN
             RESELEM(I) = 2
         ELSEIF(
     &      (Connectivity(GlobalNode1).NE.MYPROC).AND.
     &      (Connectivity(GlobalNode2).EQ.MYPROC).AND.
     &      (Connectivity(GlobalNode3).EQ.MYPROC))THEN
             RESELEM(I) = 2
         ELSEIF(
     &      (Connectivity(GlobalNode1).EQ.MYPROC).AND.
     &      (Connectivity(GlobalNode2).EQ.
     &       Connectivity(GlobalNode3)))THEN
             RESELEM(I) = 2
         ELSEIF(
     &      (Connectivity(GlobalNode2).EQ.MYPROC).AND.
     &      (Connectivity(GlobalNode3).EQ.
     &       Connectivity(GlobalNode1)))THEN
             RESELEM(I) = 2
         ELSEIF(
     &      (Connectivity(GlobalNode3).EQ.MYPROC).AND.
     &      (Connectivity(GlobalNode1).EQ.
     &       Connectivity(GlobalNode2)))THEN
             RESELEM(I) = 2
         ELSEIF(
     &      (Connectivity(GlobalNode1).EQ.MYPROC).AND.
     &      (Connectivity(GlobalNode2).NE.
     &       Connectivity(GlobalNode3)))THEN
             RESELEM(I) = 3
         ELSEIF(
     &      (Connectivity(GlobalNode2).EQ.MYPROC).AND.
     &      (Connectivity(GlobalNode3).NE.
     &       Connectivity(GlobalNode1)))THEN
             RESELEM(I) = 3
         ELSEIF(
     &      (Connectivity(GlobalNode3).EQ.MYPROC).AND.
     &      (Connectivity(GlobalNode1).NE.
     &       Connectivity(GlobalNode2)))THEN
             RESELEM(I) = 3
         ELSE
             RESELEM(I) = 0
Casey 090304: Perhaps this error message is not helpful.
C            PRINT*, "Yikes!"
         ENDIF
C           PRINT*, "RESELEM = ", RESELEM(I)
      ENDDO

      DEALLOCATE(Connectivity)
      DEALLOCATE(GlobalX,GL_X)
      DEALLOCATE(GlobalY,GL_Y)

      RETURN
      END SUBROUTINE

      SUBROUTINE MASSBAL3D(IT)

      USE GLOBAL, ONLY:      
     &  AREAS,                              ! The area of each element.
     &  DP,                                 ! The bathymetry at each node.
     &  DT,                                 ! The time step.
     &  ETA2,                               ! The most recent elevation at each node.
     &  IFNLFA,                             ! The flag for the nonlinear, finite amplitude term.
     &  RESELEM,
     &  NE,                                 ! Total number of horizontal elements. 
     &  NM,                                 ! Node numbers associated with each element.
     &  X,Y                                 ! Horizontal coordinates of node points.

      USE GLOBAL_3DVS, ONLY: 
     &  A,B,                                ! Extents of sigma coordinates.
     &  MASSCUM,                            ! Cumulative mass balance error.
     &  NFEN,                               ! Number of vertical layers.
     &  Q,                                  ! Horizontal velocities.
     &  SIGMA,                              ! Elevation of each vertical node in the sigma world. 
     &  WZ                                  ! Vertical velocities.

#ifdef CMPI
      USE MESSENGER
#endif

      IMPLICIT NONE

#ifdef CMPI
#ifndef HAVE_MPI_MOD
      INCLUDE 'mpif.h'
#endif
#endif

      INTRINSIC             :: ABS            ! Returns the absolute value of a number. 
      INTRINSIC             :: AIMAG          ! Returns the imaginary part of a complex number.
      INTRINSIC             :: MOD
      INTRINSIC             :: REAL           ! Returns the real part of a complex number.

      INTEGER               :: FACE           ! Face counter.
      INTEGER               :: HE             ! Horizontal element counter.
      INTEGER               :: HN1,HN2,HN3    ! Local horizontal nodes for the element of interest.
      INTEGER               :: IT             ! Current time step.
      INTEGER               :: VE             ! Vertical element counter.
      INTEGER               :: VN1,VN2        ! Local vertical nodes for the element of interest.

      REAL(8)               :: AVGETA 
      REAL(8), DIMENSION(3) :: AVGVEL         ! Average velocity across any face.
      REAL(8)               :: FACEAREA       ! Area of the face.
      REAL(8)               :: FLUX           ! A local net flux. 
      REAL(8)               :: FLUXCUM
      REAL(8)               :: FLUXMIN
      REAL(8)               :: FLUXTOT        ! The global net flux. 
      REAL(8), DIMENSION(3) :: NORM           ! Unit outward normal vector, see function below. 
      REAL(8), EXTERNAL     :: TRIANGLEAREA3D ! See function at end of file.
      REAL(8)               :: ZVN1,ZVN2      ! Z-coordinate of the vertical nodes.
      REAL(8)               :: ZVN3,ZVN4      ! "

      FLUXCUM = 0.D0
      FLUXMIN = 0.D0
      FLUXTOT = 0.D0
      IF( IT.EQ.1 )THEN
          MASSCUM = 0.D0
      ENDIF

Casey: Begin loop over horizontal elements.

      DO HE=1,NE

Casey: Each horizontal element is actually a stack of (NFEN-1) vertical elements,
C      and we have to compute mass balance in each one of them.  Begin loop over
C      vertical elements.

        DO VE=1,NFEN-1

            FLUX = 0.D0

Casey: Each vertical element has five faces: top, bottom, and three sides.
C      Let's compute the flux across each of the sides first.  Loop over
C      the three vertical faces of each element. 

            DO FACE=1,3
 
                VN1 = VE
                VN2 = VE + 1

                HN1 = FACE
                HN2 = FACE + 1
                IF( HN2 .EQ. 4 )THEN
                    HN2 = 1
                ENDIF

Casey: Compute the average velocity on this face.

                AVGVEL(1) = ( 0.25D0 * (  REAL(Q(NM(HE,HN1),VN1))   
     &                                 +  REAL(Q(NM(HE,HN2),VN1))   
     &                                 +  REAL(Q(NM(HE,HN1),VN2))   
     &                                 +  REAL(Q(NM(HE,HN2),VN2)) ) ) 
                AVGVEL(2) = ( 0.25D0 * ( AIMAG(Q(NM(HE,HN1),VN1))   
     &                                 + AIMAG(Q(NM(HE,HN2),VN1))   
     &                                 + AIMAG(Q(NM(HE,HN1),VN2))   
     &                                 + AIMAG(Q(NM(HE,HN2),VN2)) ) )
                AVGVEL(3) = ( 0.25D0 * (      WZ(NM(HE,HN1),VN1)    
     &                                 +      WZ(NM(HE,HN2),VN1)    
     &                                 +      WZ(NM(HE,HN1),VN2)    
     &                                 +      WZ(NM(HE,HN2),VN2)  ) )

Casey: Compute the locations of the vertical nodes in the z-direction.

                ZVN1 = DP(NM(HE,HN1)) + ( ( B - SIGMA(VN1) )/( A - B ) )
     &               * ( DP(NM(HE,HN1)) + IFNLFA * ETA2(NM(HE,HN1)) )
                ZVN2 = DP(NM(HE,HN2)) + ( ( B - SIGMA(VN1) )/( A - B ) )
     &               * ( DP(NM(HE,HN2)) + IFNLFA * ETA2(NM(HE,HN2)) )
                ZVN3 = DP(NM(HE,HN2)) + ( ( B - SIGMA(VN2) )/( A - B ) )
     &               * ( DP(NM(HE,HN2)) + IFNLFA * ETA2(NM(HE,HN2)) )
                ZVN4 = DP(NM(HE,HN1)) + ( ( B - SIGMA(VN2) )/( A - B ) )
     &               * ( DP(NM(HE,HN1)) + IFNLFA * ETA2(NM(HE,HN1)) )

Casey: Compute the outward normal vector on this face.

                CALL CROSSPRODUCT(    X(NM(HE,HN2)) - X(NM(HE,HN1))  ,  
     &                                Y(NM(HE,HN2)) - Y(NM(HE,HN1))  ,  
     &                             -( ZVN2          - ZVN1          ),  
     &                                X(NM(HE,HN1)) - X(NM(HE,HN1))  ,  
     &                                Y(NM(HE,HN1)) - Y(NM(HE,HN1))  ,  
     &                             -( ZVN4          - ZVN1          ),
     &                                NORM                           )  

Casey: Compute the area of this face.

                FACEAREA = TRIANGLEAREA3D(
     &                        X(NM(HE,HN1)),Y(NM(HE,HN1)),ZVN1,
     &                        X(NM(HE,HN2)),Y(NM(HE,HN2)),ZVN2,
     &                        X(NM(HE,HN2)),Y(NM(HE,HN2)),ZVN3)
     &                   + TRIANGLEAREA3D(
     &                        X(NM(HE,HN1)),Y(NM(HE,HN1)),ZVN1,
     &                        X(NM(HE,HN2)),Y(NM(HE,HN2)),ZVN3,
     &                        X(NM(HE,HN1)),Y(NM(HE,HN1)),ZVN4)

Casey: Compute the flux.

                FLUX = FLUX + FACEAREA *  
     &                      ( AVGVEL(1) * NORM(1)   
     &                      + AVGVEL(2) * NORM(2)   
     &                      + AVGVEL(3) * NORM(3) )

#ifdef CMPI
     &                      / REAL(RESELEM(HE))
#endif

Casey: End loop over vertical faces.

            END DO

Casey: Now let's loop over the top and bottom faces. 

            DO FACE=1,2
 
                IF( FACE.EQ.1 )THEN
                    VN1 = VE 
                ELSEIF( FACE.EQ.2 )THEN
                    VN1 = VE + 1
                ENDIF

                HN1 = 1 
                HN2 = 2 
                HN3 = 3

Casey: Compute the average velocity on this face.

                AVGVEL(1) = ( ( 1.D0 / 3.D0 ) *            
     &                        (  REAL(Q(NM(HE,HN1),VN1))   
     &                        +  REAL(Q(NM(HE,HN2),VN1))   
     &                        +  REAL(Q(NM(HE,HN3),VN1)) ) ) 
                AVGVEL(2) = ( ( 1.D0 / 3.D0 ) *            
     &                        ( AIMAG(Q(NM(HE,HN1),VN1))   
     &                        + AIMAG(Q(NM(HE,HN2),VN1))   
     &                        + AIMAG(Q(NM(HE,HN3),VN1)) ) )
                AVGVEL(3) = ( ( 1.D0 / 3.D0 ) *            
     &                        (      WZ(NM(HE,HN1),VN1)   
     &                        +      WZ(NM(HE,HN2),VN1)    
     &                        +      WZ(NM(HE,HN3),VN1)  ) )

Casey: Compute the locations of the vertical nodes in the z-direction.

                ZVN1 = DP(NM(HE,HN1)) + ( ( B - SIGMA(VN1) )/( A - B ) )
     &               * ( DP(NM(HE,HN1)) + IFNLFA * ETA2(NM(HE,HN1)) )
                ZVN2 = DP(NM(HE,HN2)) + ( ( B - SIGMA(VN1) )/( A - B ) )
     &               * ( DP(NM(HE,HN2)) + IFNLFA * ETA2(NM(HE,HN2)) )
                ZVN3 = DP(NM(HE,HN3)) + ( ( B - SIGMA(VN1) )/( A - B ) )
     &               * ( DP(NM(HE,HN3)) + IFNLFA * ETA2(NM(HE,HN3)) )

Casey: Compute the outward normal vector on this face.  If the current face
C      is the bottom of the elemental prism, then the normal vector must
C      point down; if the current face is the top, then the normal vector
C      must point up.

                IF( FACE.EQ.1 )THEN
                    CALL CROSSPRODUCT(  X(NM(HE,HN1))-X(NM(HE,HN2)) ,   
     &                                  Y(NM(HE,HN1))-Y(NM(HE,HN2)) ,   
     &                                -(ZVN1         -ZVN2         ),   
     &                                  X(NM(HE,HN3))-X(NM(HE,HN2)) ,   
     &                                  Y(NM(HE,HN3))-Y(NM(HE,HN2)) ,   
     &                                -(ZVN3         -ZVN2         ),
     &                                  NORM                         )  
                ELSEIF( FACE.EQ.2 )THEN
                    CALL CROSSPRODUCT(  X(NM(HE,HN3))-X(NM(HE,HN2)) ,   
     &                                  Y(NM(HE,HN3))-Y(NM(HE,HN2)) ,   
     &                                -(ZVN3         -ZVN2         ),   
     &                                  X(NM(HE,HN1))-X(NM(HE,HN2)) ,   
     &                                  Y(NM(HE,HN1))-Y(NM(HE,HN2)) ,   
     &                                -(ZVN1         -ZVN2         ),
     &                                  NORM                        )  
                ENDIF

Casey: Compute the area of this face.  Here, we get to use the formula
C      for the area of a triangle in three dimensions.

                FACEAREA = TRIANGLEAREA3D(
     &                     X(NM(HE,HN1)),Y(NM(HE,HN1)),ZVN1,
     &                     X(NM(HE,HN2)),Y(NM(HE,HN2)),ZVN2,
     &                     X(NM(HE,HN3)),Y(NM(HE,HN3)),ZVN3)

Casey: Compute the flux.

                FLUX = FLUX + FACEAREA *  
     &                      ( AVGVEL(1) * NORM(1)   
     &                      + AVGVEL(2) * NORM(2)   
     &                      + AVGVEL(3) * NORM(3) )

#ifdef CMPI
     &                      / REAL(RESELEM(HE))
#endif

Casey: End loop over the top and bottom faces.

            END DO

Casey 070716: Normalize the flux on each element to the deviation from the
C             still water volume in that element.

            AVGETA = (1.D0/3.D0)*(ETA2(NM(HE,1))+
     &                   ETA2(NM(HE,2))+ETA2(NM(HE,3)))
            IF(AVGETA.NE.0.D0)THEN
                FLUX = FLUX / ( AREAS(HE) * AVGETA )
            ELSE
                FLUX = FLUX / ( AREAS(HE) * 0.01D0 )
            ENDIF

            FLUXTOT = FLUXTOT + ABS( FLUX ) 

Casey: End loop over vertical elements.

        END DO

Casey: End loop over horizontal elements.

      END DO

#ifdef CMPI

      FLUXTOT = ABS(FLUXTOT*DT)

      CALL MPI_REDUCE(FLUXTOT, FLUXCUM, 1, MPI_DOUBLE_PRECISION, 
     &                MPI_SUM, 0, COMM_COMP, IERR)
      CALL MPI_REDUCE(FLUXTOT, FLUXMIN, 1, MPI_DOUBLE_PRECISION, 
     &                MPI_MIN, 0, COMM_COMP, IERR)

      IF(MYPROC.EQ.0)THEN

          MASSCUM = MASSCUM + FLUXCUM

          IF(MOD(REAL(IT),100.0).EQ.0)THEN

              PRINT *, "IT: ",IT," FLUXMIN (M^3): ",FLUXMIN,
     &                    " FLUXCUM (M^3): ",FLUXCUM,
     &                    " MASSCUM (M^3/DT): ",(MASSCUM/IT)  

          ENDIF

      ENDIF

#else

      FLUXCUM = ABS(FLUXTOT*DT)

      MASSCUM = MASSCUM + FLUXCUM

      IF(MOD(REAL(IT),100.0).EQ.0)THEN

          PRINT *, "IT: ",IT," FLUXMIN (M^3): ",FLUXMIN,
     &                " FLUXCUM (M^3): ",FLUXCUM,
     &                " MASSCUM (M^3/DT): ",(MASSCUM/IT)  

      ENDIF

#endif

      RETURN
      END SUBROUTINE MASSBAL3D




Casey: The following function calculates the area of a three-dimensional triangle
C      given the coordinates of the three vertices.

      FUNCTION TRIANGLEAREA3D(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3)

      IMPLICIT NONE

      INTRINSIC            :: ABS, SQRT

      REAL(8)              :: TRIANGLEAREA3D
      REAL(8), INTENT(IN)  :: X1, X2, X3
      REAL(8), INTENT(IN)  :: Y1, Y2, Y3
      REAL(8), INTENT(IN)  :: Z1, Z2, Z3

Casey: We get to use the formula for the area of a triangle in three dimensions.

      TRIANGLEAREA3D =  ABS( 0.5D0 * SQRT(                         
     &                    ( ( Y2 * Z3                     
     &                      - Y3 * Z2 )          
     &                    - ( Y1 * Z3                     
     &                      - Y3 * Z1 )          
     &                    + ( Y1 * Z2                   
     &                      - Y2 * Z1 ) )**2.D0   
     &                  + ( ( Z2 * X3           
     &                      - Z3 * X2 )          
     &                    - ( Z1 * X3            
     &                      - Z3 * X1 )          
     &                    + ( Z1 * X2            
     &                      - Z2 * X1 ) )**2.D0   
     &                  + ( ( X2 * Y3            
     &                      - X3 * Y2 )          
     &                    - ( X1 * Y3            
     &                      - X3 * Y1 )          
     &                    + ( X1 * Y2            
     &                      - X2 * Y1 ) )**2.D0 ) )  

      RETURN
      END FUNCTION TRIANGLEAREA3D



