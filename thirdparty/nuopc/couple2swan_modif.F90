!>
!! @mainpage ADCIRC NUOPC Cap
!! @author Saeed Moghimi (moghimis@gmail.com)
!! @date 15/1/17 Original documentation
!------------------------------------------------------
!LOG-----------------
!
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      MODULE Couple2Swan_modif

#define ONLY_COMP_FORCES
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      USE mod_logging, ONLY : setMessageSource, unsetMessageSource, &
                allMessage, DEBUG, ECHO, INFO, WARNING, ERROR

      IMPLICIT NONE

!asey 090302: These arrays contain the radiation stresses.
!asey 090820: Be explicit about the size of these REAL variables.
      REAL(8) ,ALLOCATABLE :: ADCIRC_SXX(:,:)
      REAL(8) ,ALLOCATABLE :: ADCIRC_SXY(:,:)
      REAL(8) ,ALLOCATABLE :: ADCIRC_SYY(:,:)

!asey 090302: The interpolation weight controls which value
!             is taken when information is passed between ADCIRC
!             and SWAN.  If InterpoWeight = 0, then the value
!             is taken from the beginning of the coupling interval.
!             If InterpoWeight = 1, then the value is taken from
!             the end of the coupling interval.
      REAL(8) :: InterpoWeight

!-----------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E
!       C O M P U T E  W A V E  D R I V E N  F O R C E S
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      SUBROUTINE ComputeWaveDrivenForces

      USE GLOBAL, ONLY: NODECODE, NOFF, RSNX2, RSNY2
      USE MESH, ONLY : NE, NM, NP, X, Y, AREAS, NODELE, NEITABELE
      USE BOUNDARIES, ONLY : NBDV, NBOU, NBVV, NOPE, NVDLL, NVELL

      IMPLICIT NONE

      INTEGER :: I
      INTEGER :: IE
      INTEGER :: IP
      INTEGER :: K
      INTEGER :: Node1
      INTEGER :: Node2
      INTEGER :: Node3
      INTEGER :: NUMFOUND

      LOGICAL :: Marcel = .FALSE.

      REAL(8),ALLOCATABLE :: DSXXDX(:)
      REAL(8),ALLOCATABLE :: DSXYDY(:)
      REAL(8),ALLOCATABLE :: DSXYDX(:)
      REAL(8),ALLOCATABLE :: DSYYDY(:)

      REAL(8)             :: NCELE

      REAL(8),ALLOCATABLE :: TEMP_SXX(:)
      REAL(8),ALLOCATABLE :: TEMP_SXY(:)
      REAL(8),ALLOCATABLE :: TEMP_SYY(:)

      REAL(8) :: TOTALAREA

      call setMessageSource("ComputeWaveDrivenForces")
#if defined(COUPLE2SWAN_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

!... Check whether radiation stresses have already been computed.
!... If not, then apply forces of zero.
      IF( .FALSE. )THEN

         IF(.NOT.ALLOCATED(RSNX2))THEN
            ALLOCATE(RSNX2(1:NP))
            DO IP=1,NP
               RSNX2(IP) = 0.D0
            ENDDO
         ENDIF
         IF(.NOT.ALLOCATED(RSNY2))THEN
            ALLOCATE(RSNY2(1:NP))
            DO IP=1,NP
               RSNY2(IP) = 0.D0
            ENDDO
         ENDIF

!... If so, then continue to compute wave-driven forces.
      ELSE

!... Allocate arrays for radiation stresses.
         IF(.NOT.ALLOCATED(TEMP_SXX)) ALLOCATE(TEMP_SXX(1:NP))
         IF(.NOT.ALLOCATED(TEMP_SXY)) ALLOCATE(TEMP_SXY(1:NP))
         IF(.NOT.ALLOCATED(TEMP_SYY)) ALLOCATE(TEMP_SYY(1:NP))

!... Loop over all nodes and interpolate the radiation stress for this time step.
         DO IP=1,NP
           TEMP_SXX(IP) = (1.0 - InterpoWeight) * DBLE(ADCIRC_SXX(IP,1)) &
                       + InterpoWeight * DBLE(ADCIRC_SXX(IP,2))
           TEMP_SXY(IP) = (1.0 - InterpoWeight) * DBLE(ADCIRC_SXY(IP,1)) &
                       + InterpoWeight * DBLE(ADCIRC_SXY(IP,2))
           TEMP_SYY(IP) = (1.0 - InterpoWeight) * DBLE(ADCIRC_SYY(IP,1)) &
                       + InterpoWeight * DBLE(ADCIRC_SYY(IP,2))
         ENDDO

!... Allocate arrays for radiation stress gradients.
         IF(.NOT.ALLOCATED(DSXXDX)) ALLOCATE(DSXXDX(1:NE))
         IF(.NOT.ALLOCATED(DSXYDY)) ALLOCATE(DSXYDY(1:NE))
         IF(.NOT.ALLOCATED(DSXYDX)) ALLOCATE(DSXYDX(1:NE))
         IF(.NOT.ALLOCATED(DSYYDY)) ALLOCATE(DSYYDY(1:NE))

!... Loop over all elements and compute the derivatives of Sxx, Sxy and Syy.
!... These derivatives are constant on an element.  Note that the AREAS array
!... actually contains twice the area of each element.
         DO IE=1,NE

!Casey 090707: When using the serial adcswan on Zas, I received memory errors
!... when these calls were nested into the logic below.  Break them out and
!... use these variables to save on the number of calls to memory.
           Node1 = NM(IE,1)
           Node2 = NM(IE,2)
           Node3 = NM(IE,3)

           DSXXDX(IE) = (1.D0/AREAS(IE)) *                         &
                      ( TEMP_SXX(Node1) * (Y(Node2) - Y(Node3))    &
                      + TEMP_SXX(Node2) * (Y(Node3) - Y(Node1))    &
                      + TEMP_SXX(Node3) * (Y(Node1) - Y(Node2)) )

           DSXYDY(IE) = (1.D0/AREAS(IE)) *                         &
                      ( TEMP_SXY(Node1) * (X(Node3) - X(Node2))    &
                      + TEMP_SXY(Node2) * (X(Node1) - X(Node3))    &
                      + TEMP_SXY(Node3) * (X(Node2) - X(Node1)) )

           DSXYDX(IE) = (1.D0/AREAS(IE)) *                         &
                      ( TEMP_SXY(Node1) * (Y(Node2) - Y(Node3))    &
                      + TEMP_SXY(Node2) * (Y(Node3) - Y(Node1))    &
                      + TEMP_SXY(Node3) * (Y(Node1) - Y(Node2)) )

           DSYYDY(IE) = (1.D0/AREAS(IE)) *                         &
                      ( TEMP_SYY(Node1) * (X(Node3) - X(Node2))    &
                      + TEMP_SYY(Node2) * (X(Node1) - X(Node3))    &
                      + TEMP_SYY(Node3) * (X(Node2) - X(Node1)) )

         ENDDO

!... Allocate arrays for wave-driven forces.
         IF(.NOT.ALLOCATED(RSNX2)) ALLOCATE(RSNX2(1:NP))
         IF(.NOT.ALLOCATED(RSNY2)) ALLOCATE(RSNY2(1:NP))

!... Loop over all nodes and compute the wave-driven forces:
!...
!...       Fx = - DSxx/Dx - DSxy/Dy
!...
!...       Fy = - DSxy/Dx - DSyy/Dy
!...
!... We project the element-based radiation stress gradients onto the nodes
!... by taking a weighted average of the gradients in the elements connected
!... to a node.
         outer: DO IP=1,NP

           RSNX2(IP) = 0.D0
           RSNY2(IP) = 0.D0

           TOTALAREA = 0.D0

           IE = 0
           NUMFOUND = 0

           inner: DO

             IE = IE + 1

             IF(NEITABELE(IP,IE).EQ.0)THEN

               CONTINUE

             ELSE

!... Try Marcel's method of zero-ing out the forces at nodes connected to dry nodes/elements.

               NCELE = NODECODE(NM(NEITABELE(IP,IE),1)) &
                     * NODECODE(NM(NEITABELE(IP,IE),2)) &
                     * NODECODE(NM(NEITABELE(IP,IE),3)) &
                     * NOFF(       NEITABELE(IP,IE)   )

               IF(Marcel.AND.(NCELE.EQ.0))THEN

                 RSNX2(IP) = 0.0
                 RSNY2(IP) = 0.0
                 CYCLE outer

               ELSE

                 NUMFOUND = NUMFOUND + 1

                 RSNX2(IP) = RSNX2(IP) + 0.5*AREAS(NEITABELE(IP,IE)) &
                           * ( - DSXXDX(NEITABELE(IP,IE))            &
                               - DSXYDY(NEITABELE(IP,IE)) )
                 RSNY2(IP) = RSNY2(IP) + 0.5*AREAS(NEITABELE(IP,IE)) &
                           * ( - DSXYDX(NEITABELE(IP,IE))            &
                               - DSYYDY(NEITABELE(IP,IE)) )

                 TOTALAREA = TOTALAREA + 0.5*AREAS(NEITABELE(IP,IE))

               ENDIF

             ENDIF

             IF(NUMFOUND.EQ.NODELE(IP))THEN

               EXIT inner

             ENDIF

           ENDDO inner

           RSNX2(IP) = RSNX2(IP) / TOTALAREA
           RSNY2(IP) = RSNY2(IP) / TOTALAREA

         ENDDO outer

!... Try Marcel's method of zero-ing the forces at the boundary nodes.
         IF(Marcel)THEN
           DO K=1,NOPE
             DO I=1,NVDLL(K)
               RSNX2(NBDV(K,I)) = 0.0
               RSNY2(NBDV(K,I)) = 0.0
             ENDDO
           ENDDO
           DO K=1,NBOU
             DO I=1,NVELL(K)
               RSNX2(NBVV(K,I)) = 0.0
               RSNY2(NBVV(K,I)) = 0.0
             ENDDO
           ENDDO
         ENDIF

!... Deallocate the radiation stress gradients.
         IF(ALLOCATED(DSXXDX)) DEALLOCATE(DSXXDX)
         IF(ALLOCATED(DSXYDY)) DEALLOCATE(DSXYDY)
         IF(ALLOCATED(DSXYDX)) DEALLOCATE(DSXYDX)
         IF(ALLOCATED(DSYYDY)) DEALLOCATE(DSYYDY)
         IF(ALLOCATED(TEMP_SXX)) DEALLOCATE(TEMP_SXX)
         IF(ALLOCATED(TEMP_SXY)) DEALLOCATE(TEMP_SXY)
         IF(ALLOCATED(TEMP_SYY)) DEALLOCATE(TEMP_SYY)

      ENDIF

#if defined(COUPLE2SWAN_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN
!-----------------------------------------------------------------------
      END SUBROUTINE ComputeWaveDrivenForces
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      END MODULE Couple2Swan_modif
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
