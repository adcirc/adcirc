module mod_wind_drag

   implicit none

   real(8) :: WindDragLimit = 0.0025d0 ! maximum wind drag value

   ! jgf50.32: If DragLawString is not set in the fort.15 input file,
   ! it will remain 'default'; this will cause 'Garratt' to be used
   ! in FUNCTION WindDrag, and 'IceCube' to be used in
   ! FUNCTION WindIceDrag.
   integer, parameter :: DRAGLAW_GARRATT = 99990
   integer, parameter :: DRAGLAW_POWELL = 99991
   integer, parameter :: DRAGLAW_SWELL = 99992
   integer, parameter :: DRAGLAW_GFDL = 99993

!   Casey 110518: Added for Mark Powell's sectorial wind drag.
   real(8)  :: EyeLat(3) = 0.d0
   real(8)  :: EyeLon(3) = 0.d0
   logical  :: FoundEye

#ifndef POWELL
   character(len=10) :: DragLawString = 'default' ! controls wind drag formulation
   integer :: dragLawType = DRAGLAW_GARRATT
#else
#warning The flag -DPOWELL is deprecated and will be removed in a future release. \
   Please use the&metcontrol name list functionality to avoid issues

   character(len=10) :: DragLawString = 'Powell'
   integer :: dragLawType = DRAGLAW_POWELL
#endif

   private

   public :: WindDrag, WindIceDrag, checkWindDragType, FoundEye, EyeLat, EyeLon, DragLawType
   public :: DRAGLAW_GARRATT, DRAGLAW_POWELL, DRAGLAW_SWELL, DRAGLAW_GFDL, WindDragLimit, DragLawString

contains

   !----------------------------------------------------------------------
   ! Sets an integer for the drag law to be used later in checks
   ! This is simpler than checking a string
   !----------------------------------------------------------------------
   subroutine checkWindDragType()
      use GLOBAL, only: toLowercase, allMessage, logMessage, INFO, WARNING
      implicit none
      character(len(dragLawString)) :: dragLawLowercase
      draglawlowercase = toLowercase(dragLawString)
      select case (trim(draglawlowercase))
      case ("default", "garratt")
         call logMessage(INFO, "Using Garratt wind drag formulation.")
         draglawtype = DRAGLAW_GARRATT
      case ("powell")
         call logMessage(INFO, "Using Powell sector-based wind drag formulation.")
         draglawtype = DRAGLAW_POWELL
      case ("swell")
         call logMessage(INFO, "Using SWAN-based swell wind drag formulation.")
         draglawtype = DRAGLAW_SWELL
      case ("gfdl")
         call logMessage(INFO, "Using GFDL wind drag formulation.")
         draglawtype = DRAGLAW_GFDL
      case DEFAULT
         call allMessage(WARNING, "Wind drag law '"//trim(DragLawString)//"' not recognized. Using Garratt formulation.")
         draglawtype = DRAGLAW_GARRATT
      end select
   end subroutine checkWindDragType
   !----------------------------------------------------------------------

   !  ---------------------------------------------------------------
   !> @brief Computes wind drag using a function with the signiture:
   !> result = function(windSpeed,nodeNumber)
   !> @param[in] WindSpeed wind magnitude
   !> @param[in] NodeNumber adcirc mesh node index
   !> @result wind drag coefficient
   !  ----------------------------------------------------------------
   real(8) function WindDrag(WindSpeed, NodeNumber)
      use GLOBAL, only: windDragOut
      implicit none
      real(8), intent(in) :: WindSpeed
      integer, intent(in) :: NodeNumber

      select case (dragLawType)
      case (DRAGLAW_GARRATT)
         WindDrag = CappedGarrattWindDrag(WindSpeed)
      case (DRAGLAW_GFDL)
         WindDrag = GFDLWindDrag(WindSpeed)
      case (DRAGLAW_POWELL)
         WindDrag = PowellWindDrag(WindSpeed, NodeNumber)
#ifdef CSWAN
      case (DRAGLAW_SWELL)
         WindDrag = SwanWindDrag(NodeNumber)
#endif
      case default
         WindDrag = CappedGarrattWindDrag(WindSpeed)
      end select

      windDragOut(nodeNumber) = windDrag

   end function WindDrag
   ! ----------------------------------------------------------------

#ifdef CSWAN
!  ----------------------------------------------------------------
!> @brief Use the wind drag coefficients calculated inside SWAN
!> for wind drag inside of adcirc
!> @param[in] WindSpeed wind magnitude
!> @param[in] NodeNumber adcirc mesh node index
!> @result wind drag coefficient
!  ----------------------------------------------------------------
   real(8) function SwanWindDrag(NodeNumber) result(Drag)
      use GLOBAL, only: Swan_WDragCo
      implicit none
      integer, intent(in)  :: NodeNumber
      Drag = SWAN_WDragCo(NodeNumber)
   end function SwanWindDrag
!     ----------------------------------------------------------------
#endif

   !  ----------------------------------------------------------------
   !> @brief Function applying the Garratt wind drag formulation.
   !>
   !> Wind drag is calculated using Garratt, 1977 by:
   !> \f[
   !>    d = \frac{\left( 0.75 + 0.067 \left| \vec{w} \right| \right)}{1000}
   !> \f]
   !>
   !> @param[in] WindSpeed wind magnitude
   !> @param[in] NodeNumber adcirc mesh node index
   !> @result wind drag coefficient
   ! ----------------------------------------------------------------
   real(8) pure function GarrattWindDrag(WindSpeed) result(drag)
      implicit none
      real(8), intent(in)  :: WindSpeed
      drag = 0.001d0*(0.75d0 + 0.067d0*WindSpeed)
   end function GarrattWindDrag

   ! ----------------------------------------------------------------
   !> @brief Applies the Garratt Wind Drag capped at WindDragLimit.
   !> By default the cap applied is 0.0035, however this can be adjusted
   !> using the metControl namelist parameter WindDragLimit
   !> @param[in] WindSpeed wind magnitude
   !> @param[in] NodeNumber adcirc mesh node index
   !> @result wind drag coefficient
   ! ----------------------------------------------------------------
   pure real(8) function CappedGarrattWindDrag(WindSpeed) result(WindDrag)
      implicit none
      real(8), intent(in)  :: WindSpeed
      WindDrag = min(GarrattWindDrag(WindSpeed), WindDragLimit)
   end function CappedGarrattWindDrag
   ! ----------------------------------------------------------------

   ! ----------------------------------------------------------------
   !> @brief Applies Mark Powell's sector-based wind drag. This
   !> applies different wind drag values based upon the sector of the
   !> storm a node lies within. The is done to account for the direction
   !> of wave travel and the coincident direction of the wind vector
   !> @param[in] WindSpeed wind magnitude
   !> @param[in] NodeNumber adcirc mesh node index
   !> @result wind drag coefficient
   ! ----------------------------------------------------------------
   function PowellWindDrag(WindSpeed, NodeNumber) result(PowellDrag)
      use GLOBAL, only: DEG2RAD
      use MESH, only: SFEA, SLAM
      implicit none
      intrinsic :: ATAN2
      intrinsic :: MOD
      real(8), intent(IN) :: WindSpeed
      integer, intent(IN) :: NodeNumber
      integer :: IS
      integer :: IW1, IW2
      logical :: FoundSector
      real(8) :: Dir1
      real(8) :: Dir2
      real(8) :: Drag(3)
      real(8) :: NodeDirection
      real(8) :: NodeLat
      real(8) :: NodeLon
      real(8) :: StormDirection
      real(8) :: Weight(3)
      real(8) :: PowellDrag
      real(8) :: rawGarrattWindDrag
      real(8), parameter :: eps = epsilon(1d0)

      !.. Check whether we have previously found the eye.
      if ((abs(EyeLat(1)) < eps) .or. (abs(EyeLon(1)) < eps) .or. &
          (abs(EyeLat(2)) < eps) .or. (abs(EyeLon(2)) < eps) .or. &
          (abs(EyeLat(3)) < eps) .or. (abs(EyeLon(3)) < eps)) then
         FoundEye = .false.
      end if

      !.. If the eye has not been found then apply Garratt.
      if (.not. FoundEye) then
         PowellDrag = CappedGarrattWindDrag(WindSpeed)
         FoundSector = .false.
         Weight = 0.d0
      else
         !.. We have found the eye!
         !.. Compute the storm's direction.  We use the convention of northward
         !.. being zero degrees, and then cycling clockwise.
         Dir1 = atan2(EyeLat(2) - EyeLat(1), EyeLon(2) - EyeLon(1))
         Dir2 = atan2(EyeLat(3) - EyeLat(2), EyeLon(3) - EyeLon(2))
         StormDirection = atan2(sin(Dir1) + sin(Dir2), cos(Dir1) + cos(Dir2))
         StormDirection = StormDirection/DEG2RAD
         StormDirection = 90.d0 - StormDirection
         if (StormDirection < 0.d0) then
            StormDirection = StormDirection + 360.d0
         end if
         !.. Compute the direction of node relative to the eye.
         NodeLon = SLAM(NodeNumber)/DEG2RAD
         NodeLat = SFEA(NodeNumber)/DEG2RAD
         NodeDirection = atan2(NodeLat - EyeLat(3), NodeLon - EyeLon(3))
         NodeDirection = NodeDirection/DEG2RAD
         NodeDirection = 90.d0 - NodeDirection
         if (NodeDirection < 0.d0) then
            NodeDirection = NodeDirection + 360.d0
         end if
         !.. Compute the weights for the sector in which the node is located.
         !..     0- 40 : Transition from left (3) to right (1).
         !..    40-130 : Right (1).
         !..   130-170 : Transition from right (1) to rear (2).
         !..   170-220 : Rear (2).
         !..   220-260 : Transition from rear (2) to left (3).
         !..   260-  0 : Left (3).
         FoundSector = .false.
         do IS = 1, 6
            if (IS == 1) then
               Dir1 = mod(StormDirection + 0.d0, 360.d0)
               Dir2 = mod(StormDirection + 40.d0, 360.d0)
            elseif (IS == 2) then
               Dir1 = mod(StormDirection + 40.d0, 360.d0)
               Dir2 = mod(StormDirection + 130.d0, 360.d0)
            elseif (IS == 3) then
               Dir1 = mod(StormDirection + 130.d0, 360.d0)
               Dir2 = mod(StormDirection + 170.d0, 360.d0)
            elseif (IS == 4) then
               Dir1 = mod(StormDirection + 170.d0, 360.d0)
               Dir2 = mod(StormDirection + 220.d0, 360.d0)
            elseif (IS == 5) then
               Dir1 = mod(StormDirection + 220.d0, 360.d0)
               Dir2 = mod(StormDirection + 260.d0, 360.d0)
            elseif (IS == 6) then
               Dir1 = mod(StormDirection + 260.d0, 360.d0)
               Dir2 = mod(StormDirection + 0.d0, 360.d0)
            end if

            if (Dir1 > Dir2) then
               if ((Dir1 <= NodeDirection) .and. (NodeDirection <= 360.d0)) then
                  Dir2 = Dir2 + 360.d0
               elseif ((0.d0 <= NodeDirection) .and. (NodeDirection <= Dir2)) then
                  Dir1 = Dir1 - 360.d0
               end if
            end if

            if ((Dir1 <= NodeDirection) .and. (NodeDirection <= Dir2)) then
               if (mod(IS, 2) == 0) then
                  FoundSector = .true.
                  Weight = 0.d0
                  if (is == 2) then
                     Weight(1) = 1.d0
                  elseif (is == 4) then
                     Weight(2) = 1.d0
                  elseif (is == 6) then
                     Weight(3) = 1.d0
                  end if
               else
                  FoundSector = .true.
                  Weight = 0.d0
                  if (IS == 1) then
                     IW1 = 3
                     IW2 = 1
                  elseif (IS == 3) then
                     IW1 = 1
                     IW2 = 2
                  elseif (IS == 5) then
                     IW1 = 2
                     IW2 = 3
                  end if
                  Weight(IW1) = 1.d0 + (0.d0 - 1.d0)/(Dir2 - Dir1)*(NodeDirection - Dir1)
                  Weight(IW2) = 0.d0 + (1.d0 - 0.d0)/(Dir2 - Dir1)*(NodeDirection - Dir1)
               end if
            end if
         end do

         !.. Apply garratt wind drag in case of emergency
         if (.not. FoundSector) then
            PowellDrag = CappedGarrattWindDrag(WindSpeed)
            Weight = 0.d0
         else
            rawGarrattWindDrag = GarrattWindDrag(WindSpeed)
            !.. Determine the wind drag for sector 1 (right).
            Drag(1) = rawGarrattWindDrag
            if (WindSpeed <= 35.0d0) then
               if (Drag(1) > 0.0020d0) then
                  Drag(1) = 0.0020d0
               end if
            elseif (WindSpeed <= 45.0d0) then
               Drag(1) = 0.0020d0 + (0.0030d0 - 0.0020d0)/(45.0d0 - 35.0d0)*(WindSpeed - 35.0d0)
            else
               Drag(1) = 0.0030d0
            end if

            !.. Determine the wind drag for sector 2 (rear).
            Drag(2) = rawGarrattWindDrag
            if (WindSpeed <= 35.0d0) then
               if (Drag(2) > 0.0020d0) then
                  Drag(2) = 0.0020d0
               end if
            elseif (WindSpeed <= 45.0d0) then
               Drag(2) = 0.0020d0 + (0.0010d0 - 0.0020d0)/(45.0d0 - 35.0d0)*(WindSpeed - 35.0d0)
            else
               Drag(2) = 0.0010d0
            end if

            !.. Determine the wind drag for sector 3 (left).
            Drag(3) = rawGarrattWindDrag
            if (Drag(3) > 0.0018d0) then
               if (WindSpeed <= 25.0d0) then
                  Drag(3) = 0.0018d0
               elseif (WindSpeed <= 30.0d0) then
                  Drag(3) = 0.0018d0 + (0.0045d0 - 0.0018d0)/(30.0d0 - 25.0d0)*(WindSpeed - 25.0d0)
               elseif (WindSpeed <= 45.0d0) then
                  Drag(3) = 0.0045d0 + (0.0010d0 - 0.0045d0)/(45.0d0 - 30.0d0)*(WindSpeed - 30.0d0)
               else
                  Drag(3) = 0.0010d0
               end if
            end if

            !.. Apply a weighted average.
            PowellDrag = Weight(1)*Drag(1) + Weight(2)*Drag(2) + Weight(3)*Drag(3)
         end if
      end if
   end function PowellWindDrag
   ! ----------------------------------------------------------------

   ! ----------------------------------------------------------------
   !> @brief Function applying the GFDL wind drag formulation.
   !> @param[in] WindSpeed wind magnitude
   !> @result wind drag coefficient
   ! ----------------------------------------------------------------
   real(8) pure function GFDLWindDrag(WindSpeed) result(Drag)
      implicit none
      real(8), intent(in)  :: WindSpeed

      real(8), parameter :: bs0 = -8.3672761723972770d-12
      real(8), parameter :: bs1 = 1.7398510865876079d-09
      real(8), parameter :: bs2 = -1.3318965783633590d-07
      real(8), parameter :: bs3 = 4.5070552944387270d-06
      real(8), parameter :: bs4 = -6.5086768819069140d-05
      real(8), parameter :: bs5 = 0.00044745137674732834d0
      real(8), parameter :: bs6 = -0.0010745704660847233d0

      real(8), parameter :: cf0 = 2.1151080765239772d-13
      real(8), parameter :: cf1 = -3.2260663894433345d-11
      real(8), parameter :: cf2 = -3.3297059587519610d-10
      real(8), parameter :: cf3 = 1.7648562021709124d-07
      real(8), parameter :: cf4 = 7.1076368256941820d-06
      real(8), parameter :: cf5 = -0.0013914681964973246d0
      real(8), parameter :: cf6 = 0.0406766967657759d0

      if (WindSpeed <= 5.0d0) then
         drag = (0.0185d0/9.8d0*(7.59d-4*WindSpeed**2d0 + 2.46d-2*WindSpeed)**2d0)
      elseif (WindSpeed > 5.0d0 .and. WindSpeed < 10.0d0) then
         drag = .00000235d0*(WindSpeed**2d0 - 25.0d0) + 3.805129199617346d-05
      elseif (WindSpeed >= 10.0d0 .and. WindSpeed < 60.0d0) then
         drag = bs6 + bs5*WindSpeed + bs4*WindSpeed**2 + bs3*WindSpeed**3d0 + bs2*WindSpeed**4d0 &
                + bs1*WindSpeed**5d0 + bs0*WindSpeed**6d0
      else
         drag = cf6 + cf5*WindSpeed + cf4*WindSpeed**2d0 + cf3*WindSpeed**3d0 + cf2*WindSpeed**4d0 &
                + cf1*WindSpeed**5d0 + cf0*WindSpeed**6d0
      end if

   end function GFDLWindDrag

   ! ----------------------------------------------------------------
   !  F U N C T I O N   W I N D  I C E  D R A G
   ! ----------------------------------------------------------------
   !
   ! tcm v49.64.01 -- Added a function to modify wind drag coefficents
   !      based on ice effects from a percentage of ice coverage
   !      value (ice concentration limits 0.0 to 100.0)
   !      If the percentage of ice is less than 1.0% then no
   !      modifications are done to the incoming wind drag
   !      coefficient values.
   !
   ! ----------------------------------------------------------------
   real(8) function WindIceDrag(WindDrag, PercentIce)
#ifdef CMPI
      use messenger, only: msg_fini
#endif
      implicit none
      real(8), intent(in) :: PercentIce
      real(8), intent(in) :: WindDrag
      real(8) :: pic, IceDrag, Cform
      real(8), parameter :: Cskin = 1.5d-3
      real(8), parameter :: Cform_max = 3.67d-3
      real(8), parameter :: beta = 0.6d0

      select case (trim(DragLawString))

      case ("RaysIce")
         !........If there is virtually no ice, then
         !........just use the non-ice drag values
         if (PercentIce < 1.d0) then
            WindIceDrag = WindDrag
         else
            pic = PercentIce*0.01d0
            IceDrag = (0.125d0 + 0.5d0*pic*(1.0d0 - pic))*0.01d0
            WindIceDrag = max(IceDrag, WindDrag)
         end if

      case ("IceCube")
         !........This formula is a cubic function such that
         !....... When pic = 0 the drag coefficient is equal to
         !....... the minimum value for Garratt (.000075), when pic = 50
         !....... it has a value of 0.0025 and a zero gradient, and
         !....... when pic = 100 it has a value of 0.00125.
         !....... This represents a smoother transition between Garratt
         !....... and Ray's Ice for low ice values and low winds

         if (PercentIce < 0.d0) then
            WindIceDrag = WindDrag
         else
            pic = PercentIce*0.01d0
            IceDrag = (0.075d0 + 0.75d0*pic - 0.9d0*pic*pic &
                       + 0.2d0*pic*pic*pic)*0.01d0
            WindIceDrag = max(IceDrag, WindDrag)
         end if

      case ("Lupkes", "default", "Powell")
         !........This formula is a function that explicitly
         !....... combines air-sea drag (Garratt) with sea-ice skin drag and
         !....... sea-ice form drag components with area fraction weights
         !....... see: Lüpkes, et al. (2012), JGR: doi:10.1029/2012JD017630. &
         !....... Joyce et al. (2019), OM: doi:10.1016/j.ocemod.2019.101421
         if (PercentIce < 0.d0) then
            WindIceDrag = WindDrag
         else
            pic = min(PercentIce*0.01d0, 1d0)
            Cform = Cform_max*(1d0 - pic)**beta
            WindIceDrag = WindDrag*(1d0 - pic) + (Cskin + Cform)*pic
         end if

      case DEFAULT
         write (16, *) 'ERROR: Ice drag law not recognized:'
         write (16, '(A10)') DragLawString
         write (16, *) 'Execution will now be terminated.'
#ifdef CMPI
         call msg_fini()
#endif
         call exit(1)

      end select

      if (WindIceDrag > WindDragLimit) WindIceDrag = WindDragLimit

   end function WindIceDrag
! ----------------------------------------------------------------

end module mod_wind_drag
