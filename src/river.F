

!-------------------------------------------------------------------------------!
!
! ADCIRC - The ADvanced CIRCulation model
! Copyright (C) 1994-2023 R.A. Luettich, Jr., J.J. Westerink
! 
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU Lesser General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!-------------------------------------------------------------------------------!
C******************************************************************************
C PADCIRC VERSION 56.00 xx/xx/2025   AMAN TEJASWI                             *
C                                                                             *
C                                                                             *
C This module handles normal flux boundary for rivers. Initial conditions    
C(first snap) are called in cold start subroutines and further values are 
C handled in the timestep subroutine.                                         *
C                                                                             *
C******************************************************************************

C-----------------------------------------------------------------------
C                M O D U L E     R E A D _ R I V E R
C-----------------------------------------------------------------------
C
C     CORRECTS RIVER FLUX BOUNDARY CONDITION ACCOUNTING FOR PROJECTION
C                            IN ROTATED CS
C
C.....Qx Qy IMPLEMENTATION TO CONSIDER PROJECTION IN ROTATED CS WHILE NORMAL FLOW
C...... CBOUNDARY IS USED. THE FORMULATION IS LIKE THIS.
C      QN = Qx.Nx + Qy.Ny
C      0 = TAUx.Nx + TAUy.Ny
C-----------------------------------------------------------------------

       MODULE READ_RIVER

      
       USE GLOBAL, ONLY :QN1,QN2,QN0,QX1_R,QY1_R, QX2_R, QY2_R, QN_R, QN2_R,
     &     ScreenUnit, setMessageSource, allMessage, NFFR,NX_R,
     &     NY_R,TAUX_R,TAUY_R, QNIN1, QNIN2  
       USE BOUNDARIES, ONLY : NETA, NFLUXF, NOPE, NVEL, LBCODEI, NPEBC,
     &     CSII, SIII, NVELL, NBV
       USE SIZES, ONLY : MNPROC, GLOBALDIR, MNVEL
       USE MESH, ONLY : SFac, SFMX, SFMY, SFCT,
     &     SFCX, SFCY, YCSFAC, TANPHI  
       IMPLICIT NONE

       INTEGER I,J,K
       REAL(8),allocatable,save :: FORCENODES(:)
       REAL(8),allocatable :: Q(:), Q2(:)     
       PUBLIC 
C-------------------END OF DATA DECLARATION-----C
       CONTAINS

       SUBROUTINE CONVERT_QN(QNIN1,QNIN2)
       
       IMPLICIT NONE
       REAL(8), intent(in) :: QNIN1(NVEL), QNIN2(NVEL)

       

       DO I = 1,NVEL
           Q(:) = 0.d0
           Q2(:) = 0.d0
           QX1_R(:) = 0.D0
           QY1_R(:) = 0.D0
           QX2_R(:) = 0.D0
           QY2_R(:) = 0.D0

       END DO

       DO J = 1,NVEL
         IF ((LBCODEI(J).EQ.22).OR.(LBCODEI(J).EQ.32)) THEN
           Q(J) = QNIN1(J)     
           Q2(J) = QNIN2(J)
!           print *, "QNIN1", QNIN1(J)
!           print *, "QNIN2", QNIN2(J) 
           NX_R(J) = CSII(J)/SFMX(NBV(J))
           NY_R(J) = SIII(J)
           TAUX_R(J) = -1.0d0*NY_R(J)
           TAUY_R(J) = NX_R(J)

 
           
           QY1_R(J) = Q(J)/(((-1.0d0*TAUY_R(J))/TAUX_R(J))*(NX_R(J)) + NY_R(J)) 
           QX1_R(J) = -TAUY_R(J)*QY1_R(J)/TAUX_R(J)    


           QY2_R(J) = Q2(J)/(((-1.0d0*TAUY_R(J))/TAUX_R(J))*(NX_R(J)) + NY_R(J)) 
           QX2_R(J) = -TAUY_R(J)*QY2_R(J)/TAUX_R(J)    
          

       ENDIF  
       END DO
  
         

        END SUBROUTINE CONVERT_QN



C********************************************************
C             FORMULATE QFORCE
C********************************************************
       
       SUBROUTINE FORMULATE_QFORCE(QN_R,QN2_R)

       USE MESH, ONLY : SFac, SFMX, SFMY, SFCT,
     &    SFCX, SFCY, YCSFAC, TANPHI   
      
       USE BOUNDARIES, ONLY: NBOU, NVELL, NVEL, NBV, LBCODEI
       IMPLICIT NONE 
       REAL(8), intent(out) :: QN_R(NVEL), QN2_R(NVEL)
      
       INTEGER P

        DO P = 1,NVEL
         IF ((LBCODEI(P).EQ.22).OR.(LBCODEI(P).EQ.32)) THEN 

          QN_R(P) = SFCX(FORCENODES(P))*QX1_R(P)*CSII(P) +
     &        SFCY(FORCENODES(P))*QY1_R(P)*SIII(P)*YCSFAC(FORCENODES(P))
          
          QN2_R(P) = SFCX(FORCENODES(P))*QX2_R(P)*CSII(P) +
     &        SFCY(FORCENODES(P))*QY2_R(P)*SIII(P)*YCSFAC(FORCENODES(P))
         

        
         ENDIF
       END DO
       
       RETURN
 
       END SUBROUTINE FORMULATE_QFORCE       
      
      
C*********************************************************
C              INITIALIZE RIVER
C*********************************************************      
       
       SUBROUTINE INIT_RIVER()
       
       USE MESH, ONLY : SFac, SFMX, SFMY, SFCT,
     &    SFCX, SFCY, YCSFAC, TANPHI   
      
       USE BOUNDARIES, ONLY: NBOU, NVELL, NVEL, NBV

       IMPLICIT NONE 
       ALLOCATE(Q(MNVEL), Q2(MNVEL),FORCENODES(MNVEL))
       DO J = 1, NVEL
        IF ((LBCODEI(J).EQ.22).OR.(LBCODEI(J).EQ.32)) THEN 

              FORCENODES(J) = NBV(J)
             
        ENDIF
       END DO
      
       END SUBROUTINE INIT_RIVER


      END MODULE READ_RIVER
        
        
      













