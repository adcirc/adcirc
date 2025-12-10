!-------------------------------------------------------------------------------!
!
! ADCIRC - The ADvanced CIRCulation model
! Copyright (C) 1994-2025 R.A. Luettich, Jr., J.J. Westerink
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

!-----------------------------------------------------------------------
!  MODULE KDTreeTools
!-----------------------------------------------------------------------
!> @copyright Dr. R.A. Luettich and Dr. J.J. Westerink
!>
!> @author Chris Massey, USACE-ERDC-CHL, Chris.Massey@usace.army.mil
!>            Wrote original code that was in read_global.F
!>         Chris Szpilka, University of Oklahome, cmszpilka@ou.edu
!>            Moved to separate module and added a logical variable to
!>            track whether a search tree has been created.
!>
!> @brief Module to handle setup and searching using KDTREE2 algorithm.
!>
!> @details
!>        tcm v49.48.01 Replaced CoordToEle with KDTREE2 algorithm
!>        cms 56.02  Moved kdtsearch and setup routines out of
!>            read_global.F and created a separate module so that they
!>            can be used for hydrology also. Otherwise, read_global
!>            and hydroprep were codependent.
!>
!> @dependencies
!>        This module uses pre_global, kdtree2_module
!-----------------------------------------------------------------------
      MODULE KDTreeTools

      implicit none

      ! cms 11/25/25 Added for check of KDtree status - still need to
      ! think about adding IF statements around the calls to kdtsearch
      ! will numStations = 0 cause trouble?
      logical :: kdtSetup = .false.

      public :: setup_kdt_search, kdtsearch, kdtSetup
     

      CONTAINS

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!                             P U B L I C
!----------------------------------------------------------------------
!----------------------------------------------------------------------


!----------------------------------------------------------------------
!> @brief Creates lists of element centers and radii for kdtree search.
!----------------------------------------------------------------------
      subroutine setup_kdt_search
!----------------------------------------------------------------------

      use pre_global, only : SLAM,SFEA,NNEG,RMAX,BCXY,NELG

      implicit none
      integer :: i
      integer :: ielm(3)
      real(8) :: xelm(3),yelm(3),shplocal(7)

      do i=1,NELG
         ielm(:) = NNEG((/1,2,3/),i)  !element's node numbers
         xelm(:) = slam(ielm(:))      !element's vertex x-values
         yelm(:) = sfea(ielm(:))      !element's vertex y-values
         shplocal(1)= yelm(2)-yelm(3)
         shplocal(2)= yelm(3)-yelm(1)
         shplocal(3)= yelm(1)-yelm(2)
         shplocal(4)= xelm(3)-xelm(2)
         shplocal(5)= xelm(1)-xelm(3)
         shplocal(6)= xelm(2)-xelm(1)
         shplocal(7)= shplocal(1)*shplocal(5) - shplocal(2)*shplocal(4)

!        compute the radius of the circle that circumscribes
!        the triangle then scale it by 50% larger to allow for
!        a buffer later on
         rmax(i)=1.50D0*( sqrt(shplocal(6)*shplocal(6)+                 &
     &                         shplocal(3)*shplocal(3))*                &
     &                    sqrt(shplocal(4)*shplocal(4)+                 &
     &                         shplocal(1)*shplocal(1))*                &
     &                    sqrt(shplocal(5)*shplocal(5)+                 &
     &                         shplocal(2)*shplocal(2))/                &
     &                     (2.d0*shplocal(7)) )
!        Compute the barycenter of the element
          bcxy(1,i) = sum(xelm)/3.d0
          bcxy(2,i) = sum(yelm)/3.d0
       enddo !i

      kdtSetup = .true.

      return
!----------------------------------------------------------------------
      end subroutine setup_kdt_search
!----------------------------------------------------------------------


!----------------------------------------------------------------------
!> @brief Uses KDTREE2 algorithm to find element a point lies in.
!----------------------------------------------------------------------
      subroutine kdtsearch(XCoords,YCoords,FDEle,NumofStations,Desc)
!----------------------------------------------------------------------

      use PRE_Global
      use kdtree2_module

      implicit NONE
      INTEGER, intent(in) :: NumOfStations                     ! total
      REAL(8), intent(in), dimension(NumOfStations) :: XCoords ! cartesian
      REAL(8), intent(in), dimension(NumOfStations) :: YCoords ! cartesian
      INTEGER, intent(out), dimension(NumOfStations):: FDEle   ! FullDomain
      CHARACTER(len=30), intent(in) :: Desc                    ! description

      INTEGER :: I,itc,iek
      INTEGER :: ielm(3)
      real(8) :: Xsta,Ysta,dist
      real(8) :: x1,x2,x3,y1,y2,y3,A1,A2,A3,AA,AREASK,AE
      real(8) :: elmmin(2),xelm(3),yelm(3)
      TYPE(KDTREE2), POINTER :: TREE
      TYPE(KDTREE2_RESULT), ALLOCATABLE :: KDRESULTS(:)
      LOGICAL ElementFound  ! .true. when a corresponding element is found
      INTEGER ClosestElement  ! element with closest match
      REAL(8), PARAMETER :: Tolerance = 1.0d-5     ! area difference for match

      ElementFound = .false.

!    Be sure the maximum search depth is not larger than
!    the number of elements being kept
      if (nelg.lt.srchdp) srchdp = nelg

!     Create the search tree
      tree => kdtree2_create(bcxy,rearrange=.true.,sort=.true.)

!     allocate space for the search results from the tree
      allocate(KDRESULTS(srchdp))


      do I = 1, NumOfStations
         Xsta = XCoords(I)
         Ysta = YCoords(I)
         ElementFound = .false.

        ! Find the srchdp nearest elements to this point
         call kdtree2_n_nearest(tp=tree,qv=(/Xsta,Ysta/),               &
     &                  nn=srchdp,results=KDRESULTS)
         itc = 1
         ClosestElement = KDRESULTS(itc)%idx
!       Check to see if the point lies with rmax of any of these elements
         elmmin = minval(sqrt(KDRESULTS(1:srchdp)%dis)                  &
     &                              - rmax(KDRESULTS(1:srchdp)%idx) )
         if(elmmin(1).le.0.0D0) then  ! Point lies within search radius of an element
!           loop through the elements in the search list
            do while ((ElementFound.eqv..false.).and.(itc.le.srchdp))
               iek = KDRESULTS(itc)%idx  !Current search element number
!              Get the distance from this point to the barycenter of the
!              current element
               dist = sqrt(KDRESULTS(itc)%dis)
!              If the distance is less than or equal to rmax (rmax=1.5*element radius)
!              Then the point is near the element and might be in it
!              Proceed with the weights test
               if(dist-rmax(iek).le.0.0d0) then
                  !get the shape function for this element
                  ielm(:) = NNEG((/1,2,3/),iek)  !element's node numbers
                  xelm(:) = slam(ielm(:))      !element's vertex x-values
                  yelm(:) = sfea(ielm(:))      !element's vertex y-values
                  X1=xelm(1)
                  X2=xelm(2)
                  X3=xelm(3)
                  Y1=yelm(1)
                  Y2=yelm(2)
                  Y3=yelm(3)
                  A1=(Xsta-X3)*(Y2-Y3)+(X2-X3)*(Y3-Ysta)
                  A2=(Xsta-X1)*(Y3-Y1)-(Ysta-Y1)*(X3-X1)
                  A3=(Ysta-Y1)*(X2-X1)-(Xsta-X1)*(Y2-Y1)
                  AA=ABS(A1)+ABS(A2)+ABS(A3)
                  AREASK=X2*Y3+X1*Y2+X3*Y1-Y1*X2-Y2*X3-Y3*X1
                  AE=ABS(AA-AREASK)/AREASK
               IF (AE.LT.Tolerance) THEN
                  ElementFound = .true.
                  ClosestElement = iek
                  FDEle(I) = ClosestElement
                  else !not in this element keep looking
                    itc = itc + 1
                  endif !End area ratio test
               else !
!                point is too far away from the barycenter of the
!                element to possibly be in the element, so move to
!                the next element
                 itc = itc + 1
               endif !end Radius test
            enddo !end the while loop
         endif
         if(.not.ElementFound) then !We did not find the element
            write(*,1234) Desc, I
            print *, "Please check the coordinates."
            IF (NFOVER.EQ.1) THEN
               print *, "The program will estimate nearest element."
               print *, "WARNING. Distance to nearest element is ",     &
     &                   sqrt(KDRESULTS(1)%dis)
               print *, " "
               FDEle(I) = ClosestElement
            ELSE
               print *, "ERROR. Distance to nearest element is ",       &
     &                   sqrt(KDRESULTS(1)%dis)
               CALL EXIT(1)
            ENDIF
         ENDIF
      enddo !loop over station points

!     Deallocate the tree
      call kdtree2_destroy(tp=tree)

 1234 format(A30,1x,I4,1x,'does not lie in the grid.')

      return

!----------------------------------------------------------------------
      end subroutine kdtsearch
!----------------------------------------------------------------------


      end module KDTreeTools
