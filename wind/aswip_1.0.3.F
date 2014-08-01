C-----------------------------------------------------------------------
C
C     M O D U L E    A S Y M M E T R I C  V O R T E X  P R O C E S S I N G
C
C-----------------------------------------------------------------------
      module asymmetricVortexProcessing
         use sizes, only : sz
         use global, only : np, ne, slam, sfea, wvnx2, wvny2, prn2,
     &                   nm, nhy, g, rhowat0, ics, slam0, sfea0,
     &                   deg2rad, rad2deg, ms2kt, ten2one,
     &                   nws, windReduction, sphericalDistance,
     &                   DEBUG, WARNING, ERROR, INFO, screenunit,
     &                   setMessageSource, allMessage, openFileForRead,
     &                   screenMessage, initLogging, scratchMessage
         ! Import asymmetric hurricane vortex class
         USE vortex
         USE wind, only : writeFullCircleRmaxes, vortexLon,
     &                 vortexLat, writeRadialVandP,
     &                 writeSpatialVandP, nws19get, windDrag, nws20get,
     &                 writeSpatialuvp, writeNodaluvp,
     &                 geostrophic_switch, geofactor
         integer :: num_entry, num_cycles
         integer ,dimension(:),allocatable :: advr
         integer,dimension(:),allocatable :: iFcstInc ! hours between forecasts
         integer,dimension(:),allocatable :: iyear, imth, iday, ihr
         character(len=4),dimension(:),allocatable :: castType! hindcast/nowcast or forecast?
         integer,dimension(:),allocatable :: ilat,ilon
         real(sz),dimension(:),allocatable :: rlat,rlon
         integer,dimension(:),allocatable :: ispd,icpress
         integer , dimension(:,:), allocatable :: irad ! working isotach radii
         integer ,dimension(:),allocatable :: ivr, dir, speed
         real(sz) , dimension(:,:), allocatable :: rmaxw
         integer ,dimension(:),allocatable :: ipn
         CHARACTER(1) :: ns ! N or S in latitude column
         CHARACTER(1) :: ew ! E or W in longitude column
         integer  :: ient,icyc,itpc,method,j,iquadrot,nquadrot,approach
         real(sz),dimension(:),allocatable::cycletime
         integer,dimension(:),allocatable::atcfrmw
         real(sz) :: vmax !,speedmax,rmax,rmaxmax
         real(sz) :: vr,speedmax2
         real(sz) :: vmax_pseudo
         real(sz), dimension(4) :: r
         real(sz) ,dimension(:),allocatable :: uTrans, vTrans
         real(sz) ,dimension(:),allocatable :: uTrans2, vTrans2
         REAL(sz),DIMENSION(:),ALLOCATABLE :: HollB
         REAL(sz),DIMENSION(:,:),ALLOCATABLE :: HollBs
         REAL(sz),DIMENSION(:,:),ALLOCATABLE :: PhiFactors
         REAL(sz),DIMENSION(:,:),ALLOCATABLE :: VmaxesBL
         CHARACTER (LEN = 10),DIMENSION(:),ALLOCATABLE :: stormname
         INTEGER, DIMENSION(0:5) :: lookup_radii ! periodic interpolation
         real(sz) :: pn
         real(sz) :: pc
         real(sz) :: clat
         real(sz) :: clon
         real(sz) :: dx,dy,dr
         real(sz) :: rhowatg
         real(sz) :: theta
         real(sz), dimension(4) :: gamma ! factor applied to the StormMotion
         real(sz), dimension(4) :: quadRotateAngle, quadRotateAngle_new
         real(sz), dimension(4) :: rmaxwHighIso
         real(sz), dimension(4) :: epsilonAngle
         integer, dimension(4) :: irr
         logical, dimension(4) :: vioflag
         integer , dimension(:,:), allocatable :: ir   ! given isotach radii
         integer , dimension(:,:), allocatable :: quadflag
         integer ,dimension(:),allocatable :: cycle_num
         integer ,dimension(:),allocatable :: isotachs_per_cycle
         real(sz) ,dimension(:),allocatable :: h_speed, h_dir
         real(sz) :: vmaxbl ! max sustained wind at top of atm. b.l.
         real(sz) :: vrbl   ! wind vel at top of atm. b.l. at wind radius
         real(sz), dimension(4) :: quadrantvr
         real(sz), dimension(4) :: quadrantangles
         real(sz), dimension(4) :: quadrantvecangles
         real(sz), dimension(4) :: VmwBL
         integer,  dimension(4) :: VmwBLflag
         real(sz), allocatable :: cyclesecs(:) ! time in seconds that each cycle
                                               ! corresponds to
         real(sz), allocatable :: cyclelons(:) ! longitudes for each cycle
         real(sz), allocatable :: cyclelats(:) ! latitudes for each cycle
         real(sz), allocatable :: cycledirs(:) ! directions for each cycle
         LOGICAL :: writeNWS19,writeNWS20 ! .false. if we don't want to write it
         REAL(sz) :: rmax_percent  ! desired percentage change in Rmax relative to the
                           ! value calculated from isotach radii data (0 to 100)
         REAL(sz) :: forecastRmax ! constant Rmax value for the duration of the
                               ! forecast (nautical miles)
         LOGICAL :: constantRmax ! .true. if forecastRmax was specified as
                              ! something .gt. zero
                              ! ... takes precedence over percent
                              ! change in Rmax
         LOGICAL :: percentSpecified ! .true. if percent was specified and is
                                  ! something other than 100 percent
         REAL(sz) :: vmax_percent  ! desired percentage change in Vmax relative to the
                                ! given value (original Vmax still used to
                                ! calculate Rmax
         INTEGER :: radii_sum ! record radius values for filling in missing vals
         INTEGER :: num_nonzero ! number of nonzero isotach radii
         INTEGER :: first_entry    ! first entry in the cycle
         INTEGER :: last_entry     ! last entry in the cycle
         LOGICAL :: higher_isotach ! true if a valid Rmax exists at higher isotach
         REAL(sz) :: stormMotion  ! portion of Vmax attributable to storm motio
         REAL(sz) :: stormMotionU ! U portion of Vmax attributable to storm motion
         REAL(sz) :: stormMotionV ! V portion of Vmax attributable to storm motion
         REAL(sz) :: U_Vr
         REAL(sz) :: V_Vr
         REAL(sz) :: temp  

      !-------------------------------------------------------------
      contains
      !-------------------------------------------------------------


      !-------------------------------------------------------------
C     ! S U B R O U T I N E
      !   P R O C E S S   A S Y M M E T R I C   T R A C K   F I L E
      !-------------------------------------------------------------
      subroutine processAsymmetricTrackFile(infile)
      implicit none
      character(1024), intent(in) :: infile

      integer :: errorIO = 0
      integer :: i,j,k

      ! jgf49.30 Direct Rmax specification takes precedence over percentage
      ! variation of Rmax.

      IF ((constantRmax.eqv..true.).and.
     &           (percentSpecified.eqv..true.)) THEN
         write(screenunit,*)
     &     "WARNING: Rmax variation was specified as a percent (-P ",
     &     rmax_percent,") as well as a constant (-R ",forecastRmax,")."
         write(screenunit,*)
     &   "WARNING: The percentage variation for Rmax will not be used."
      ENDIF

C
C     jgf49.30 Check to see if the fort.22 file is there and if not, emit error
C     message and terminate.
      call openFileForRead(22,infile,errorIO)
      if (errorIO.ne.0) then
         call screenMessage(ERROR,"aswip: Could not open file.")
         stop
      endif

      num_entry=0
      DO
         READ(22,*,END=991)
         num_entry=num_entry+1
      ENDDO
991   continue
      rewind(22)

      allocate(advr(num_entry),iyear(num_entry),imth(num_entry))
      allocate(iday(num_entry),ihr(num_entry),castType(num_entry))
      allocate(iFcstInc(num_entry),ilat(num_entry),ilon(num_entry))
      allocate(ispd(num_entry))
      allocate(rlat(num_entry),rlon(num_entry))
      allocate(irad(num_entry,4))
      allocate(icpress(num_entry), ivr(num_entry),ir(num_entry,4))
      allocate(rmaxw(num_entry,4), ipn(num_entry))
      allocate(dir(num_entry),speed(num_entry))
      allocate(cycle_num(num_entry))
      allocate(isotachs_per_cycle(num_entry),CycleTime(num_entry))
      allocate(quadflag(num_entry,4),stormname(num_entry))
      allocate(uTrans(num_entry+1), vTrans(num_entry+1))
      allocate(HollB(num_entry),atcfRMW(num_entry))
      allocate(HollBs(num_entry,4)) !varying B
      allocate(PhiFactors(num_entry,4))!correction factor for HollBs and Vh
      allocate(VmaxesBL(num_entry,4)) ! spatially varying Vmax
      allocate(uTrans2(num_entry), vTrans2(num_entry)) !,HollB(i),atcfRMW(i) )
      allocate(h_speed(num_entry), h_dir(num_entry))

      num_cycles=1
      isotachs_per_cycle(num_cycles)=1

      DO ient=1,num_entry

         READ(22,22) advr(ient), iyear(ient),imth(ient),iday(ient),
     &      ihr(ient), castType(ient),iFcstInc(ient), ilat(ient),ns,
     &      ilon(ient),ew, ispd(ient), icpress(ient), ivr(ient),
     &      (ir(ient,i),i=1,4), ipn(ient), atcfRMW(ient),
     &      dir(ient),speed(ient), stormname(ient)

22       FORMAT(3x, i3, 2x, i4, 3i2, 6x, a4, 2(2x,i3),a1, 1x, i5,
     &         a1, 2x, i3,2x, i4, 6x, i3, 7x, 4(i4,2x), i4, 8x, i3,
     &          27x, 2(i3,2x), a10)


         IF(ew.EQ.'W')THEN
            rlon(ient)=-0.1d0*real(ilon(ient))
         ELSE
            rlon(ient)=0.1d0*real(ilon(ient))
         ENDIF
         IF(ns.EQ.'S')THEN
            rlat(ient)=-0.1d0*real(ilat(ient))
         ELSE
            rlat(ient)=0.1d0*real(ilat(ient))
         ENDIF
C
         ! check to see if this is a new cycle; if so, set the cycle time
         if (ient.ne.1) then
            if ( iFcstInc(ient) .ne. iFcstInc(ient-1)) then
               num_cycles = num_cycles + 1
               isotachs_per_cycle(num_cycles) = 1
            else
               ! increment the number of isotachs for this cycle if this
               ! entry belongs to the same cycle as the last
               isotachs_per_cycle(num_cycles) =
     &             isotachs_per_cycle(num_cycles) + 1
            endif
         endif
         cycle_num(ient) = num_cycles
         cycleTime(ient) = iFcstInc(ient) * 3600.d0

      ENDDO
      CLOSE(22)

      ! calculate overland translation velocity in m/s and knots,
      ! set background pressure
      do ient=1, num_entry

         ! Set ipn to background pressure of 1013
         ipn(ient)=1013
         ! check/set central pressure
         if (icpress(ient).eq.0) then
            if ( ient.eq.1 ) then
               call screenMessage(ERROR,
     &           "aswip: Central pressure set to zero on first line.")
               stop
            else
               call screenMessage(WARNING,
     &        "aswip: Central pressure persisted from previous value.")
               icpress(ient) = icpress(ient-1)
            endif
         endif

         if ( CycleTime(ient).eq.CycleTime(1) ) then
            ! we will come back to the first cycle at the end
            cycle
         endif
         if ( CycleTime(ient).eq.CycleTime(ient-1) ) then
            uTrans(ient)=uTrans(ient-1)
            vTrans(ient)=vTrans(ient-1)
            h_speed(ient)=h_speed(ient-1)
         else
            ! approximate u and v translation velocities
            call uvTrans(rlat(ient-1),rlon(ient-1),
     &                    rlat(ient),rlon(ient),
     &                    cycleTime(ient-1),cycleTime(ient),
     &                    uTrans(ient),vTrans(ient))
            ! get translation speed
            h_speed(ient) = sphericalDistance(
     &         (deg2rad*(rlon(ient)-rlon(ient-1))),
     &         (deg2rad*(rlat(ient)-rlat(ient-1))),
     &         rlat(ient-1),rlat(ient))
            h_speed(ient) = ms2kt *
     &          (h_speed(ient) / (cycleTime(ient)-cycleTime(ient-1)))
         endif

      enddo
      ! now set the translation velocity in the first cycle equal
      ! to the translation velocity in the 2nd cycle, for lack of any
      ! better information
      do ient=2, num_entry
         if ( CycleTime(ient).ne.CycleTime(1) ) then
            uTrans(1:(ient-1)) = uTrans(ient)
            vTrans(1:(ient-1)) = vTrans(ient)
            h_speed(1:(ient-1))= h_speed(ient)
            exit
         endif
      enddo
      !
      ! convert Utrans and Vtrans to speed and direction
      ! direction is in compass coordinates 0 == North
      ! increasing clockwise
      DO ient=1,num_entry
         IF (h_speed(ient) .lt. 1.0d0 ) then
            !jgf50.29: The vortex module can't handle speed and direction
            ! being zero; it will return NaNs as a result. Persist the
            ! direction from the previous cycle, and make the storm translation
            ! speed small but nonzero.
            h_speed(ient) = 1.0d0
            if (ient.gt.1) then
               h_dir(ient) = h_dir(ient-1)
            else
               h_dir(ient)=0.0
            endif
         else
            ! calculate angle in compass coordinates
            h_dir(ient) = rad2deg * ATAN2(uTrans(ient),vTrans(ient))
            if (h_dir(ient) .lt. 0.d0) then
               h_dir(ient) = h_dir(ient) + 360.d0
            endif
         ENDIF
      ENDDO

      !-----------------------------------------
      !  Now using the calculated translational velocities
      !  call the vortex module and compute the Rmax's
      !  to be used in the new input file
      !----------------------------------------------
      irad(:,:) = ir(:,:)
      select case(abs(nws))
      case(19)
      DO ient=1,num_entry
         !
         ! transform variables from integers
         ! to real numbers for hurricane vortex calcualtions.
         Vmax =  real( ispd(ient)  )
         Pn   =  real( ipn(ient)   )
         Pc   =  real( icpress(ient)  )
         cLat = rlat(ient)
         cLon = rlon(ient)
         !
         !  need to get some logic incase Vr is  zero
         !  if so we will also be setting ir(:) to Rmax
         !  ... this happens when storms are at the "invest" stage
         if ( ivr(ient) == 0 ) then
            Vr = Vmax
         else
            Vr = real(ivr(ient))
         endif

         ! jgf50.36: fill missing values based on number of nonzero radii.
         lookup_radii(0)=ir(ient,4)
         lookup_radii(5)=ir(ient,1)
         radii_sum = 0
         num_nonzero = 0
         do i=1,4
            lookup_radii(i) = ir(ient,i)
            radii_sum = radii_sum + ir(ient,i)
            if (ir(ient,i).gt.0) then
               num_nonzero = num_nonzero + 1
               quadflag(ient,i) = 1   ! use the Rmax resulting from this
            else
               quadflag(ient,i) = 0   ! don't use Rmax resulting from this
            endif
         end do
         ! fill missing values based on how many are missing
         select case(num_nonzero)
         case(0) ! no isotachs reported, use overall Rmax; set isotach to Vmax
            quadflag(ient,:) = 1
            if (atcfRMW(ient).ne.0) then
               irad(ient,:) = atcfRMW(ient)
            else
               irad(ient,:) = 40 ! need a nonzero value for Rmax calcs,
                                 ! this val will be thrown away later
            endif
            Vr=Vmax
         case(1) ! set missing radii equal to half the nonzero radius
            where (ir(ient,:).eq.0) irad(ient,:) = 0.5 * radii_sum
         case(2) ! set missing to half the avg of the 2 radii that are given
            where (ir(ient,:).eq.0) irad(ient,:) = 0.5 * radii_sum * 0.5
         case(3) ! set missing radius to half the average of the
                 ! radii on either side
            do i=1,4
               if (ir(ient,i).eq.0) then
                  irad(ient,i) = 0.5*
     &               (lookup_radii(i+1)+lookup_radii(i-1))
               endif
            enddo
         case(4)
            ! use all these radii as-is
         case default
            ! the following error message should be unreachable
            write(16,*) "ERROR: aswip: Number of nonzero radii on ",
     &      "line ",ient," not in range 0 to 4."
         end select

         DO i=1,4
            r(i) = real(irad(ient,i))
         END DO
         HollB(ient)=1.d0
         !-------------------------------------------------------
         ! Create a new asymmetric hurricane vortex.
         !
         ! Note: Subtract translational speed from Vmax, then
         ! scale (Vmax - Vt) and Vr up to the top of the surface,
         ! where the cylcostrophic wind balance is valid.
         !-------------------------------------------------------
         stormMotion = 1.5d0 * h_speed(ient)**0.63d0
         stormMotionU = sin(h_dir(ient)*deg2rad)*stormMotion
         stormMotionV = cos(h_dir(ient)*deg2rad)*stormMotion

         VmaxBL = ( Vmax - stormMotion ) / windReduction

         do i=1,4
            ! quadrant angles are in the radial direction, we need
            ! the tangential direction, b/c that is the direction of Vr
            U_Vr = Vr * cos(quadrantAngles(i) + (deg2rad*90.0d0))
            V_Vr = Vr * sin(quadrantAngles(i) + (deg2rad*90.0d0))
            quadrantVr(i) = sqrt(
     &                   (U_Vr - stormMotionU)**2.0d0
     &                  +(V_Vr - stormMotionV)**2.0d0
     &                  )/windReduction
         end do

         ! check to see if any of the isotach wind speeds end up being
         ! greater than Vmax, and if so, don't apply the translation speed
         do i=1,4
            if ( quadrantVr(i).gt.VmaxBL ) then
               quadrantVr(:) = Vr/windReduction
               VmaxBL = Vmax/windReduction
               exit
            endif
         enddo
         call setUseQuadrantVr(.true.)
         call newVortex(Pn,Pc,cLat,cLon,VmaxBL)
         call setIsotachWindSpeeds(quadrantVr)
         call setIsotachRadii(r)
         call calcRmaxes()
         call getRmaxes(rmaxw(ient,:))
         HollB(ient) = getShapeParameter()

C        reset rmax to zero if there was a zero radius to the isotach for all
C        isotachs EXCEPT the 34 kt isotach.  in that case leave the radius that
C        has been substituted.

C        jgf49.30: The isotach wind speed can sometimes be zero in cases
C        where all radii are zero (this has been observed in the BEST
C        track file for IGOR2010). Including this possibility in the if
C        statement, so that we can avoid setting the quadrant Rmax to zero
C        if ivr was zero.
         DO i=1,4
            if ( (ivr(ient).ne.34).AND.(ivr(ient).ne.0).AND.
     &           (ir(ient,i).eq.0) ) then
               rmaxw(ient,i) = 0.0
            endif
         END DO

      ENDDO ! num_entry (main do loop)
      case(20)
      DO icyc = 1, num_cycles
             last_entry = sum(isotachs_per_cycle(1:icyc))
       DO k = 1, isotachs_per_cycle(icyc)
!            last_entry = sum(isotachs_per_cycle(1:icyc))
!            first_entry = last_entry-(isotachs_per_cycle(icyc)-1)
          ient= last_entry+1-k
          write(*,*) "Start    processing ient = ", ient
         !
         ! transform variables from integers
         ! to real numbers for hurricane vortex calcualtions.
         Vmax =  real( ispd(ient)  )
         Pn   =  real( ipn(ient)   )
         Pc   =  real( icpress(ient)  )
         cLat = rlat(ient)
         cLon = rlon(ient)
         !
         !  need to get some logic incase Vr is  zero
         !  if so we will also be setting ir(:) to Rmax
         !  ... this happens when storms are at the "invest" stage
         if ( ivr(ient) == 0 ) then
            Vr = Vmax
         else
            Vr = real(ivr(ient))
         endif

         ! jgf50.36: fill missing values based on number of nonzero radii.
         lookup_radii(0)=ir(ient,4)
         lookup_radii(5)=ir(ient,1)
         radii_sum = 0
         num_nonzero = 0
         do i=1,4
            lookup_radii(i) = ir(ient,i)
            radii_sum = radii_sum + ir(ient,i)
            if (ir(ient,i).gt.0) then
               num_nonzero = num_nonzero + 1
               quadflag(ient,i) = 1   ! use the Rmax resulting from this
            else
               quadflag(ient,i) = 0   ! don't use Rmax resulting from this
            endif
         end do
         ! fill missing values based on how many are missing
         select case(num_nonzero)
         case(0) ! no isotachs reported, use overall Rmax; set isotach to Vmax
            quadflag(ient,:) = 1
            if (atcfRMW(ient).ne.0) then
               irad(ient,:) = atcfRMW(ient)
            else
               irad(ient,:) = 40 ! need a nonzero value for Rmax calcs,
                                 ! this val will be thrown away later
            endif
            Vr=Vmax
         case(1) ! set missing radii equal to half the nonzero radius
            where (ir(ient,:).eq.0) irad(ient,:) = 0.5 * radii_sum
         case(2) ! set missing to half the avg of the 2 radii that are given
            where (ir(ient,:).eq.0) irad(ient,:) = 0.5 * radii_sum * 0.5
         case(3) ! set missing radius to half the average of the
                 ! radii on either side
                 ! Jie: looks like the average of the radii on either side
            do i=1,4
               if (ir(ient,i).eq.0) then
                  irad(ient,i) = 0.5*
     &               (lookup_radii(i+1)+lookup_radii(i-1))
               endif
            enddo
         case(4)
            ! use all these radii as-is
         case default
            ! the following error message should be unreachable
            write(16,*) "ERROR: aswip: Number of nonzero radii on ",
     &      "line ",ient," not in range 0 to 4."
         end select

         DO i=1,4
            r(i) = real(irad(ient,i))
         END DO
         HollB(ient)=1.d0
         HollBs(ient,1:4)=1.d0
         PhiFactors(ient,1:4)=1.d0
         irr=ir(ient,:)
         !-------------------------------------------------------
         ! Create a new asymmetric hurricane vortex.
         !
         ! Note: Subtract translational speed from Vmax, then
         ! scale (Vmax - Vt) and Vr up to the top of the surface,
         ! where the cylcostrophic wind balance is valid.
         !-------------------------------------------------------
         ! Jie 2013.02
117   FORMAT(8(f6.3,x))

         call setUseVmaxesBL(.true.)

         stormMotion = 1.5d0 * h_speed(ient)**0.63d0
         stormMotionU = sin(h_dir(ient)*deg2rad)*stormMotion
         stormMotionV = cos(h_dir(ient)*deg2rad)*stormMotion
         VmaxBL = ( Vmax - stormMotion ) / windReduction
         
         select case(approach)
         case(1) !Normal approach: assume Vr and quadrantVr vectors 
                 !are both tangential to azimuth 
                 
          do i = 1, 4
            ! quadrant angles are in the radial direction, we need
            ! the tangential direction, b/c that is the direction of Vr

            U_Vr = Vr * cos(quadrantAngles(i) + (deg2rad*90.0d0))
            V_Vr = Vr * sin(quadrantAngles(i) + (deg2rad*90.0d0))  
                      
            ! Jie 2013.01 
            ! eliminate the translational speed based on vortex wind speed
            ! gamma = |quadrantVr| / |VmaxBL|  
            gamma(i) = ((2.0d0*U_Vr*stormMotionU+2.0d0*V_Vr*stormMotionV)
     &          -sqrt((2.0d0*U_Vr*stormMotionU+2.0d0*V_Vr*stormMotionV
     &          )**2.0d0-4.0d0*(stormMotion**2.0d0-VmaxBL**2.0d0
     &          *windReduction**2.0d0)*Vr**2.0d0))/(2.0d0*
     &         (stormMotion**2.0d0-VmaxBL**2.0d0*windReduction**2.0d0))
            gamma(i) = max(min(gamma(i),1.0d0),0.0d0)
            
            quadrantVr(i) = sqrt(
     &                   (U_Vr - gamma(i)*stormMotionU)**2.0d0
     &                  +(V_Vr - gamma(i)*stormMotionV)**2.0d0
     &                  )/windReduction
          end do
         
         ! If violation occurs at any quadrant (quadrantVr(i)>VmaxBL),
         ! re-calculate quadrantVr at those violated quadrants

          do i = 1, 4
          if ( quadrantVr(i).gt.VmaxBL ) then
         ! Jie: replace Vmax with Vr when violation occurs (including 
         ! situations when isotach is not reported at thta quadrant: 
         ! especially at the investment stage or for the highest isotachs
         ! that is not always available.
           U_Vr = Vr * cos(quadrantAngles(i) + (deg2rad*90.0d0))
           V_Vr = Vr * sin(quadrantAngles(i) + (deg2rad*90.0d0))
           if (ir(ient,i).gt.0) then
              Vmax_pseudo = Vr
              VmwBL(i) = sqrt(
     &             (Vmax_pseudo*cos(quadrantAngles(i)+(deg2rad*90.0d0))-
     &             stormMotionU)**2.0d0 +
     &             (Vmax_pseudo*sin(quadrantAngles(i)+(deg2rad*90.0d0))-
     &             stormMotionV)**2.0d0
     &             ) / windReduction
            
            gamma(i) = ((2.0d0*U_Vr*stormMotionU+2.0d0*V_Vr*stormMotionV)
     &          -sqrt((2.0d0*U_Vr*stormMotionU+2.0d0*V_Vr*stormMotionV
     &          )**2.0d0-4.0d0*(stormMotion**2.0d0-VmwBL(i)**2.0d0
     &          *windReduction**2.0d0)*Vr**2.0d0))/(2.0d0*
     &         (stormMotion**2.0d0-VmwBL(i)**2.0d0*windReduction**2.0d0))
            !gamma(i) = max(min(gamma(i),1.0d0),0.0d0)
            
            quadrantVr(i) = sqrt(
     &                   (U_Vr - gamma(i)*stormMotionU)**2.0d0
     &                  +(V_Vr - gamma(i)*stormMotionV)**2.0d0
     &                  )/windReduction
            else
              VmwBL(i) = VmaxBL  
            ! gamma = |quadrantVr| / |VmaxBL|  
            gamma(i) = ((2.0d0*U_Vr*stormMotionU+2.0d0*V_Vr*stormMotionV)
     &          -sqrt((2.0d0*U_Vr*stormMotionU+2.0d0*V_Vr*stormMotionV
     &          )**2.0d0-4.0d0*(stormMotion**2.0d0-VmaxBL**2.0d0
     &          *windReduction**2.0d0)*Vr**2.0d0))/(2.0d0*
     &         (stormMotion**2.0d0-VmaxBL**2.0d0*windReduction**2.0d0))
            gamma(i) = max(min(gamma(i),1.0d0),0.0d0)
            
            quadrantVr(i) = (Vr - gamma(i)*stormMotion)/windReduction !scalar cal.
            endif 
     
           else
              VmwBL(i) = VmaxBL
           endif
           end do 
             
         !VmaxesBL(ient,1:4)=VmwBL(1:4)  
         
         call setUseQuadrantVr(.true.)
         call newVortexFull(Pn,Pc,cLat,cLon,VmaxBL)
         HollB(ient) = getShapeParameter()
         call setIsotachWindSpeeds(quadrantVr)
         call setIsotachRadii(r)
         call setVmaxesBL(VmwBL)
         if (geostrophic_switch .eqv. .true.) then
             call calcRmaxesFull()
         else
             call calcRmaxes()
         endif    
         call getRmaxes(rmaxw(ient,:)) 
     
        case(2) !An updated approach: assume quadrantVr has an 
              ! additional inward angnel quadRotateAngle, and Vr angle is not known
              ! calculate quadRotateAngle for the highest isotach, and then 
              ! use it for other lower isotachs of the same cycle_num
        VmwBLflag = 0
        
        if (k == 1) then   
                 nquadrot = 300
                 quadRotateAngle(:) = 25.0d0 ! initial guess of inward rotation angle (degree)
        else
          do i = 1,4
            quadRotateAngle(i) = fang(r(i),rmaxwHighIso(i))
          enddo
            nquadrot = 1
        endif
        ! Add loop to converge inward rotation angle
        DO iquadrot=1,nquadrot
        vioflag = .false.
          do i = 1, 4
           quadrantvecangles(i) = quadrantAngles(i)+
     &             (90.0d0+quadrotateAngle(i))*deg2rad 
          enddo   
            ! radial direction -> tangential direction -> 
            ! add inward direction -> the direction of quadrantVr
            
          do i = 1, 4
            if ((iquadrot == 1).or.(VmwBLflag(i) == 0)) then
               epsilonAngle(i)= 360.0d0 + rad2deg * 
     &         ATAN2(VmaxBL*sin(quadrantvecangles(i))+ stormMotionV,
     &         VmaxBL*cos(quadrantvecangles(i))+ stormMotionU)
     
               if (epsilonAngle(i) > 360.d0) then 
                  epsilonAngle(i) = epsilonAngle(i) - 360.d0
               endif
                 
               U_Vr = Vr * cos(epsilonAngle(i)/rad2deg)
               V_Vr = Vr * sin(epsilonAngle(i)/rad2deg)  
                      
              ! Jie 2013.01 
              ! eliminate the translational speed based on vortex wind speed
              ! gamma = |quadrantVr| / |VmaxBL|  
              gamma(i) = ((2.0d0*U_Vr*stormMotionU+2.0d0*V_Vr*stormMotionV)
     &          -sqrt((2.0d0*U_Vr*stormMotionU+2.0d0*V_Vr*stormMotionV
     &          )**2.0d0-4.0d0*(stormMotion**2.0d0-VmaxBL**2.0d0
     &          *windReduction**2.0d0)*Vr**2.0d0))/(2.0d0*
     &         (stormMotion**2.0d0-VmaxBL**2.0d0*windReduction**2.0d0))
              gamma(i) = max(min(gamma(i),1.0d0),0.0d0)
            
              quadrantVr(i) = sqrt(
     &                   (U_Vr - gamma(i)*stormMotionU)**2.0d0
     &                  +(V_Vr - gamma(i)*stormMotionV)**2.0d0
     &                  )/windReduction
            endif
          end do
 
         ! write(2222,112) VmaxBL, quadrantVr  
         ! If violation occurs at any quadrant (quadrantVr(i)>VmaxBL),
         ! re-calculate quadrantVr at those violated quadrants
           
          do i = 1, 4
          if (( quadrantVr(i).gt.VmaxBL ).or.(VmwBLflag(i) == 1)) then
         ! Jie: replace Vmax with Vr when violation occurs (including 
         ! situations when isotach is not reported at that quadrant)
            if (iquadrot == 1) VmwBLflag(i) = 1 ! assign violation flags
            if (ir(ient,i).gt.0) then
                     quadrantVr(i) = (-2.d0*(
     &                        stormMotionU*cos(quadrantvecangles(i))  +
     &                        stormMotionV*sin(quadrantvecangles(i))) + 
     &               sqrt(4.d0*(stormMotionU*cos(quadrantvecangles(i))+ 
     &                    stormMotionV*sin(quadrantvecangles(i)))**2.d0
     &                       - 4.d0*(stormMotion**2.d0-Vr**2.d0)))
     &                      /2.d0
     
               epsilonAngle(i)= 360.0d0 + rad2deg * 
     &         ATAN2(quadrantVr(i)*sin(quadrantvecangles(i))+ stormMotionV,
     &         quadrantVr(i)*cos(quadrantvecangles(i))+ stormMotionU)    
     
               if (epsilonAngle(i) > 360.d0) then 
                   epsilonAngle(i) = epsilonAngle(i) - 360.d0
               endif 
               quadrantVr(i)= quadrantVr(i)/windReduction
               VmwBL(i) = quadrantVr(i)
             else
               VmwBL(i) = VmaxBL 
               U_Vr = Vr * sin(h_dir(ient)*deg2rad)
               V_Vr = Vr * sin(h_dir(ient)*deg2rad) 
               ! gamma = |quadrantVr| / |VmaxBL|  
               gamma(i) = ((2.0d0*U_Vr*stormMotionU+2.0d0*V_Vr*stormMotionV)
     &            -sqrt((2.0d0*U_Vr*stormMotionU+2.0d0*V_Vr*stormMotionV
     &            )**2.0d0-4.0d0*(stormMotion**2.0d0-VmaxBL**2.0d0
     &            *windReduction**2.0d0)*Vr**2.0d0))/(2.0d0*
     &           (stormMotion**2.0d0-VmaxBL**2.0d0*windReduction**2.0d0))
               gamma(i) = max(min(gamma(i),1.0d0),0.0d0)
            
               quadrantVr(i) = (Vr - gamma(i)*stormMotion)/windReduction !scalar cal.
             endif 
     
           else
              VmwBL(i) = VmaxBL
           endif
           end do       

         
           call setUseQuadrantVr(.true.)
           call newVortexFull(Pn,Pc,cLat,cLon,VmaxBL)
           HollB(ient) = getShapeParameter()
           call setIsotachWindSpeeds(quadrantVr)
           call setIsotachRadii(r)
           call setVmaxesBL(VmwBL)
           if (geostrophic_switch .eqv. .true.) then
               call calcRmaxesFull()
           else
               call calcRmaxes()
           endif
           call getRmaxes(rmaxw(ient,:))
           ! add deterministic statement to exit the loop when conditions met
          if (k==1) rmaxwHighIso(:) = rmaxw(ient,:)
          do i = 1, 4
            quadRotateAngle_new(i)=fang(r(i),rmaxwHighIso(i))
            if (abs(quadRotateAngle_new(i)- quadRotateAngle(i))>0.2d0) then
                vioflag(i) = .true.
            endif
          enddo

          if ((count(vioflag) >= 1)
     &     .and.(iquadrot<nquadrot)) then
             quadRotateAngle(:) = quadRotateAngle_new(:)
          else              
             exit
          endif   
          
          where (.not.vioflag) irr = 0
          if ((sum(irr(:))== 0).and.(iquadrot==2)) exit
            
        ENDDO  ! iquadrot = 1,nquadrot 
                        
          if ((iquadrot.ge.nquadrot).and.(k==1)
     &       .and.(sum(irr(:))/=0)) then
!          write(*,*) "quadRotateAngle not fully converge, ient=", ient
           write(*,*) "Converge issue at ient = ", ient, " iquadrot = ", iquadrot 
           write(*,*) vioflag, irr
           write(*,117) quadRotateAngle(:) ,quadRotateAngle_new(:)
          else
          write(*,*) "Finished processing ient = ", ient, "iquadrot = ", iquadrot
          end if

       case default
          write(*,*) "Wrong approach #, must be 1 or 2"
       end select  
       
         call getVmaxesBL(VmaxesBL(ient,:))
         HollBs(ient,1:4) = getShapeParameters()
         PhiFactors(ient,1:4)=getPhiFactors()

C        reset rmax to zero if there was a zero radius to the isotach for all
C        isotachs EXCEPT the 34 kt isotach.  in that case leave the radius that
C        has been substituted.

C        jgf49.30: The isotach wind speed can sometimes be zero in cases
C        where all radii are zero (this has been observed in the BEST
C        track file for IGOR2010). Including this possibility in the if
C        statement, so that we can avoid setting the quadrant Rmax to zero
C        if ivr was zero.
         DO i=1,4
            if ( (ivr(ient).ne.34).AND.(ivr(ient).ne.0).AND.
     &           (ir(ient,i).eq.0) ) then
               rmaxw(ient,i) = 0.0
            endif
         END DO
         
      ENDDO ! isotachs_per_cycle
      ENDDO ! num_entry (main do loop)
      case default
        write(*,*) "ERROR: nws=", nws
        write(*,*) "Processing only available for nws=19,20"
      end select
                  
C     !-------------------------------------
      ! Now indicate which isotach quadrant radius
      ! that the user desires ADCIRC to read in
      ! for the final calculation of RMX in the
      ! Asymmetric Holland wind calculations
      !
      ! 34... - 0 0 0 0 ...
      ! 50... - 0 0 1 1 ...
      ! 64... - 1 1 0 0 ...
      !
      ! would indicate -
      ! use NO radii from the 34 kt isotach
      ! use the 3 & 4 radii form the 50 kt isotach
      ! use the 1 & 2 radii form the 64 kt isotach

      ! users can then modify the input file
      ! to indicate which set of radii to use
      ! for each cycle
      !
      !  Loop through each cycle and choose
      !  the isotach radii to use
      !
      !  method 1
      !  use the 34kt isotach only (like original NWS=9)
      !
      !  method 2
      !  use the fancy way of taking the highest
      !  isotach Rmax that exists
      !
      !  method 3
      !  use preferably the 50kt isotach Rmax in each quadrant, 
      !  if not available, use the 34kt one
      ! 
      !  method 4
      !  use all available isothch for each cycle, 
      !  linearly weighted-combination will be performed in 
      !  nws20get module
      !  
      !------------------------------------
      select case(method)
      case(1) ! just use the Rmaxes from the 34kt isotach
         do ient=1,num_entry
            if ((ivr(ient).eq.34).or.(ivr(ient).eq.0)) then
               quadflag(ient,:)=1
            else
               quadflag(ient,:)=0
            endif
         enddo
      case(2) ! use the Rmax from the highest isotach in each quadrant
         do icyc=1,num_cycles
            last_entry = sum(isotachs_per_cycle(1:icyc))
            first_entry = last_entry-(isotachs_per_cycle(icyc)-1)
            if (isotachs_per_cycle(icyc).eq.1) then
               quadflag(last_entry,:)=1
            else
               ! loop over quadrants
               do i=1,4
                  num_nonzero =
     &               count(quadflag(first_entry:last_entry,i).ne.0)
                  select case(num_nonzero)
                     case(0,1) ! none, or only 34kt isotach has a radius value
                        quadflag(first_entry,i) = 1
                     case(2) ! the 34kt and 50kt isotachs have radius value
                        quadflag(first_entry,i) = 0
                     case(3) ! the 34kt, 50kt, and 64kt isotachs have values
                        quadflag(first_entry:first_entry+1,i) = 0
                     case default ! zero isotachs have been flagged
                        write(*,*)
     &                  "ERROR: aswip: ",num_nonzero," isotachs were ",
     &                      "nonzero."
                  end select
               enddo
            endif
         enddo
      case(3) ! use preferably the Rmax from the 50kt isotach in each quadrant
         do icyc=1,num_cycles
            last_entry = sum(isotachs_per_cycle(1:icyc))
            first_entry = last_entry-(isotachs_per_cycle(icyc)-1)
            if (isotachs_per_cycle(icyc).eq.1) then
               quadflag(last_entry,:)=1
            else
               ! loop over quadrants
               do i=1,4
                  num_nonzero =
     &               count(quadflag(first_entry:last_entry,i).ne.0)
                  select case(num_nonzero)
                     case(0,1) ! none, or only 34kt isotach has a radius value
                        quadflag(first_entry,i) = 1
                     case(2) ! the 34kt and 50kt isotachs have radius value
                        quadflag(first_entry,i) = 0
                     case(3) ! the 34kt, 50kt, and 64kt isotachs have values
                        quadflag(first_entry,i) = 0
                        quadflag(last_entry,i) = 0
                     case default ! zero isotachs have been flagged
                        write(*,*)
     &                  "ERROR: aswip: ",num_nonzero," isotachs were ",
     &                      "nonzero."
                  end select
               enddo
            endif
         enddo
      case(4) ! use all available Rmaxes from multiple reported isotachs
         do icyc=1,num_cycles
            last_entry = sum(isotachs_per_cycle(1:icyc))
            first_entry = last_entry-(isotachs_per_cycle(icyc)-1)
            ! since quadflag is previously assigned 1 where (ir(ient,:)>0)
            ! here we only have to deal with situations when only 0 or 34
            ! isotach is reported and with missing ir values
            if (isotachs_per_cycle(icyc).eq.1) then 
                quadflag(last_entry,:)=1
            else
               ! loop over quadrants
               do i=1,4
                  num_nonzero =
     &               count(quadflag(first_entry:last_entry,i).ne.0)
                  select case(num_nonzero)
                     case(0,1) ! none, or only 34kt isotach has a radius value
                        quadflag(first_entry,i) = 1
                     case(2,3) ! the 34kt, 50kt, and/or 64kt isotachs have values
                     
                     case default
                        write(*,*)
     &                  "ERROR: aswip: ",num_nonzero," isotachs were ",
     &                      "nonzero." 
                     end select     
                 enddo       
            endif                
         enddo   
      case default
         write(*,*) "ERROR: method=",method,
     &      " is not valid for setting rmax in quadrants."
         write(*,*) "ERROR: Execution terminated."
      end select

      ! persist last good 34kt Rmax values if all radii are missing
      do icyc=1,num_cycles
         if (isotachs_per_cycle(icyc).eq.1) then
            ient = sum(isotachs_per_cycle(1:icyc))
            if ((all(ir(ient,:).eq.0)).and.(atcfRMW(ient).eq.0)) then
               if ((icyc-1).ge.1) then
                  rmaxw(ient,:) =
     &                     rmaxw(ient-isotachs_per_cycle(icyc-1),:)
               else
                  rmaxw(ient,:) = 25 ! default value when all else fails
               endif
            endif
         endif
      enddo
      !-------------------------------------------------------------
      end subroutine processAsymmetricTrackFile
      !-------------------------------------------------------------


      !-------------------------------------------------------------
      ! S U B R O U T I N E
      !     W R I T E  A S Y M M E T R I C  I N P U T  F I L E
      !-------------------------------------------------------------
      subroutine writeAsymmetricInputFile(nwswrite)

      implicit none
      integer :: i
      integer, intent(in) :: nwswrite

      if (nwswrite == 19) then
        OPEN(222,file="NWS_19_fort.22",status="replace")
      else if (nwswrite == 20) then
         if (geostrophic_switch .eqv. .true.) then
             OPEN(222,file="NWS_20_fort.22",status="replace")
         else
             OPEN(222,file="NWS_20_cyclostr_fort.22",status="replace")
         endif
         OPEN(6666,file="PhiFactors.6666",status="replace")    
116      FORMAT(4(f5.2,x))
      end if
      
      DO ient=1,num_entry
         !-----------------------------------
         ! NHC advisory best-track i/o format
         !-----------------------------------
22          FORMAT(3x, i3, 2x, i4, 3i2, 6x, a4, 2(2x,i3),a1, 1x, i5,
     &         a1, 2x, i3,2x, i4, 6x, i3, 7x, 4(i4,2x), i4, 8x, i3,
     &          27x, 2(i3,2x), a10)
         !-----------------------------------
         ! NHC advisory best-track i/o format
         ! WITH nws=19 ADDED DATA
         !-----------------------------------
26          FORMAT("AL,",1x,i2.2,",",i5,3i2.2,",",3x,",",1x,a4,",",
     &         i4,",",i4,a1,",",i5,a1,",",i4,",",i5,",",3x,",",i4,",",
     &       " NEQ,",5(i5,","),5x,",",i4,", ",5(4x,","),i3,", ",i3,",",
     &         a10," ",1x,",",i4,",",i5,",",1x, 4(i1,",",1x),2x,
     &         4(f6.1,",",1x), 2x, f8.4)!,2(2x,f10.4))
     
         !-----------------------------------
         ! NHC advisory best-track i/o format
         ! WITH nws=20 ADDED DATA
         !-----------------------------------
27          FORMAT("AL,",1x,i2.2,",",i5,3i2.2,",",3x,",",1x,a4,",",
     &         i4,",",i4,a1,",",i5,a1,",",i4,",",i5,",",3x,",",i4,",",
     &       " NEQ,",5(i5,","),5x,",",i4,", ",5(4x,","),i3,", ",i3,",",
     &         a10," ",1x,",",i4,",",i5,",",1x, 4(i1,",",1x),2x,
     &         4(f6.1,",",1x),1x, 8(f8.4,",",1x), f8.4)     

         ! jgf49.30: Modify the Rmax as specified on the command line.
         IF (constantRmax.eqv..true.) THEN
            DO i=1,4
               rmaxw(ient,i) = forecastRmax
            ENDDO
         ELSEIF (percentSpecified.eqv..true.) THEN
            DO i=1,4
               rmaxw(ient,i) = rmaxw(ient,i)
     &                              *(1.d0+rmax_percent/100.d0)
            ENDDO
         ENDIF
         ! modify the Vmax as specified on the command line
         if ( vmax_percent.ne.0.d0 ) then
            ispd(ient) = nint(real(ispd(ient))
     &                              *(1.d0+vmax_percent/100.d0))
         endif
C
         select case(abs(nws))
         case(19)     
           WRITE(222,26) advr(ient), iyear(ient),imth(ient),iday(ient),
     &         ihr(ient), castType(ient),iFcstInc(ient), ilat(ient),ns,
     &         ilon(ient),ew, ispd(ient), icpress(ient), ivr(ient),
     &         (ir(ient,i),i=1,4), ipn(ient), atcfRMW(ient),
     &         nint(h_dir(ient)),nint(h_speed(ient)), stormname(ient),
     &         cycle_num(ient),isotachs_per_cycle(cycle_num(ient)),
     &        (quadflag(ient,i),i=1,4),(rmaxw(ient,i),i=1,4),HollB(ient)
         case(20)
           WRITE(222,27) advr(ient), iyear(ient),imth(ient),iday(ient),
     &         ihr(ient), castType(ient),iFcstInc(ient), ilat(ient),ns,
     &         ilon(ient),ew, ispd(ient), icpress(ient), ivr(ient),
     &         (ir(ient,i),i=1,4), ipn(ient), atcfRMW(ient),
     &         nint(h_dir(ient)),nint(h_speed(ient)), stormname(ient),
     &         cycle_num(ient),isotachs_per_cycle(cycle_num(ient)),
     &        (quadflag(ient,i),i=1,4),(rmaxw(ient,i),i=1,4),
     &         HollB(ient),(HollBs(ient,i),i=1,4),
     &         (VmaxesBL(ient,i),i=1,4)

           WRITE(6666,116)  ! Jie 2013.02
     &        (PhiFactors(ient,i),i=1,4)
         case default
         end select
      ENDDO
      close(222)
      if (abs(nws)==20) close(6666)
      !-------------------------------------------------------------
       end subroutine writeAsymmetricInputFile
      !-------------------------------------------------------------

C-----------------------------------------------------------------------
      end module asymmetricVortexProcessing
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C
C     M O D U L E    D I A G N O S T I C   O U T P U T
C
C-----------------------------------------------------------------------
      module diagnosticOutput

      use sizes, only : sz, globaldir, GBLINPUTDIR
      use global, only : slam, sfea, wvnx2, wvny2, prn2, np, ics, nws,
     &   deg2rad, rad2deg, ten2one, g, rhowat0, openFileForRead,
     &   mnp, ihot, wtime1, wtime2, wtiminc,
     &   m2nm, sphericaldistance  ! Jie 09/2013
      use wind, only : hWindData_t, hWindFiles, nws15init, nws19get,
     &   writeFullCircleRmaxes, writeRadialVandP, writeSpatialVandP,
     &   nws15get, eyeLatR, eyeLonR, nws20get,writeSpatialuvp,
     &   writeNodaluvp, moving_grid, geostrophic_switch, geofactor
      use vortex, only : xy2latlon, latlon2xy, fang
      use asymmetricVortexProcessing, only :
     &  cycleTime, rlat, rlon, quadrantAngles, cycle_num, h_dir,
     &  cycleSecs, cycleLats, cycleLons, cycleDirs, num_cycles,
     &  num_entry, icyc, ient
      use owiwind, only : nws12init, nws12get ! Jie 09/2013 added OWI
      CHARACTER(2048) :: meshFileName ! name of file containing node locations
      LOGICAL :: nodesFromMesh ! .true. if we are supposed to load a fort.14
      LOGICAL :: nodesFromFile ! .true. if we are supposed to load a list of nodes     
C
C     Variables for recording domain extents
      REAL(sz) :: xmin ! smallest x coordinate in mesh
      REAL(sz) :: ymin ! smallest y coordinate in mesh
      REAL(sz) :: xmax ! largest x coordinate in mesh
      REAL(sz) :: ymax ! largest y coordinate in mesh
      INTEGER :: num_outputs ! number of writes to output file
      INTEGER :: iout        ! output counter
      REAL(sz), ALLOCATABLE :: outputSecs(:) ! time in seconds for output
      REAL(sz), ALLOCATABLE :: outputLons(:) ! longitudes at each output time
      REAL(sz), ALLOCATABLE :: outputLats(:) ! latitudes at each output time
      REAL(sz), ALLOCATABLE :: outputDirs(:) ! storm direction at output times
      REAL(sz) :: interp   ! time interpolation factor for time varying params
      REAL(sz) :: trigDir
      REAL(sz) :: azimuth  ! angle of node w.r.t. vortex (radians)
      REAL(sz) :: vtest    ! wind speed along Rmax curve
      REAL(sz) :: vmax_azimuth ! azimuthal angle at which vmax occurs
      INTEGER, PARAMETER :: radial_extent = 400 ! (nm) max dist. to write V,P
      INTEGER, PARAMETER :: quadrant_radial_resolution = 2 ! data points per nm
      INTEGER, PARAMETER :: spatial_radial_increment = 1
      INTEGER, PARAMETER :: azimuthal_increment = 5
      INTEGER :: num_rvalues ! dimension in radial direction for spatial output
      INTEGER :: num_avalues ! dimension in the azimuthal direction for spatial
      INTEGER :: ry, t  ! loop counters in radial and theta directions
      CHARACTER(80) :: rmax_file_name ! name of full circle rmax output
      CHARACTER(80) :: radialvp_file_name ! name of radial V and P output
      CHARACTER(80) :: node_file_name ! name of file containing node locations
      CHARACTER(80) :: station_output_name ! name of output at given stations
      CHARACTER(80) :: spatial_file_name ! file containing spatial output
      CHARACTER(80) :: scratch  ! throw-away line from file
      REAL(sz) :: output_increment    ! time increment for output (seconds)
      LOGICAL :: defaultOutputIncrement = .true. ! equal to cycle time incr.
      INTEGER :: frame ! output frame number to render (0 for all frames)
      INTEGER :: cell_list_size ! size of elements list for vtk output
      LOGICAL :: meshFileReadComplete ! only need to read fort.14 once
      CHARACTER(1024) :: elementTableFileName !needed to post process spatialvp
      LOGICAL :: elementTableProvided = .false. ! changes output format!
      REAL(sz), ALLOCATABLE :: nodex(:) ! cpp projected x coord (meters)
      REAL(sz), ALLOCATABLE :: nodey(:) ! cpp projected x coord (meters)
      INTEGER, ALLOCATABLE :: outer_boundary_nodes(:) ! node numbers
      INTEGER, ALLOCATABLE :: inner_boundary_nodes(:) ! node numbers
      REAL(sz), DIMENSION(:), ALLOCATABLE :: pressure
      REAL(sz), DIMENSION(:), ALLOCATABLE :: radial_speed
      REAL(sz), DIMENSION(:), ALLOCATABLE :: tempspd

      REAL(sz) :: totalWindStress ! magnitude, summed over all nodes
      REAL(sz) :: windMag ! wind speed, m/s
      REAL(sz) :: windDragCoef ! wind drag coefficient
      REAL(sz) :: stressX  ! wind stress in x-dir
      REAL(sz) :: stressY  ! wind stress in y-dir
      REAL(8) :: timeOffset ! like a hotstart time
      REAL(sz) :: owi_skip ! total number of owi snaps to be skipped
      REAL(sz) :: owi_output_factor !frequency to output owi snaps
      real(sz) :: dist     ! (m) distance from storm center to mesh node
      real(sz) :: angle    ! (nm) angle to location of mesh vmax (on subdomain in parallel)
      real(sz) :: maxSpeed ! (m/s) max wind speed on mesh (subdomain in parallel)
      integer  :: maxNode  ! mesh node on which maxSpeed occurs
      !-------------------------------------------------------------
      contains
      !-------------------------------------------------------------

      !-------------------------------------------------------------
      ! S U B R O U T I N E
      !    G E N E R A T E   D I A G N O S T I C   O U T P U T
      !-------------------------------------------------------------
      subroutine initializeDiagnosticOutput()
      use wind, only : numFiles
      implicit none
      integer :: i
      integer :: errorio

!      nodesFromMesh=.false.
!      meshFileReadComplete=.false.
      frame=0
      ics = 2
C
      select case(nws)
      case(8)
         ! open the file we just wrote, so NWS8 can access it
         call openFileForRead(22,'fort.22',errorIO)
         if (errorIO.ne.0) stop
      case(12) ! Jie 2013/09
         ! For ASWIP program, read fort.22* for additional info
         ! Obtain storm center location beforehand, so we can 
         ! form a mesh around it each time
         call openFileForRead(22,'fort.22',errorIO)
         if (errorIO.ne.0) stop
         READ(22,*) !# of files
         READ(22,*) owi_skip ! skip or insert snaps 
         READ(22,*) ! wind multiplier
         ! The above 3 parameters will be treated by nws12init
         READ(22,*) num_entry
         READ(22,*) wtiminc ! in seconds
         READ(22,*) owi_output_factor ! frequency to output owi snaps (default=1)
         close(22)
         num_cycles = num_entry + owi_skip     
         allocate(cycleTime(num_cycles))
         DO i = 1, num_cycles
             cycleTime(i) = wtiminc*(i-1)
         ENDDO
      case(15)
         num_cycles = numFiles
         num_entry = numFiles
         allocate(cycleTime(numFiles))
         ! we can't actually fill the rlat and rlon arrays with
         ! correct values all at once for NWS15, b/c the storm
         ! center locations are in the individual HWind files,
         ! which are read on the fly during the computation,
         ! rather than all at once at the beginning
         cycleTime(:) = hWindFiles(:)%cycleTime ! note: hotstart accounted for nws15
      case(19)
         ! open the file we just wrote, so NWS19 can access it
         call openFileForRead(22,'NWS_19_fort.22',errorIO)
         if (errorIO.ne.0) stop
      case(20)
         ! open the file we just wrote, so NWS20 can access it
         if (geostrophic_switch .eqv. .true.) then
             call openFileForRead(22,'NWS_20_fort.22',errorIO)
         else
             call openFileForRead(22,'NWS_20_cyclostr_fort.22',errorIO)
         endif
         if (errorIO.ne.0) stop
      case default
         write(*,*) "ERROR: aswip: -n ",nws,
     &        " is not an nws value supported by aswip."
      end select
      ! place the time in seconds to which each cycle corresponds,
      ! along with the storm position for each cycle, in new arrays
      ! so that we can conveniently loop through them
      if ((abs(nws).eq.12).and.(owi_output_factor.gt.1)) then
         defaultOutputIncrement = .false. 
      end if  
      if ( defaultOutputIncrement.eqv..true. ) then
         num_outputs = num_cycles
      else
!         num_outputs =
!     &      int((cycleTime(num_entry)-timeOffset)/output_increment)
        if (abs(nws).eq.12) then
          num_outputs = 
     &     int(num_cycles/owi_output_factor)+1  
        else    
          ! add +1 to total number of outputs, so we won't miss the last frame Jie 09/2013
          num_outputs =
     &    int((cycleTime(num_entry)-timeOffset)/output_increment)+1
        endif
      endif
      allocate(outputSecs(num_outputs),cycleSecs(num_cycles))
      allocate(outputLats(num_outputs),cycleLats(num_cycles))
      allocate(outputLons(num_outputs),cycleLons(num_cycles))
      allocate(outputDirs(num_outputs),cycleDirs(num_cycles))
      ! for NWS19 and NWS20, need to go through the array and pick out the
      ! values that correspond to each cycle
      if ((abs(nws).eq.19).or.(abs(nws).eq.20)) then
         cycleSecs(1) = cycleTime(1)
         cycleLons(1) = rlon(1)
         cycleLats(1) = rlat(1)
         cycleDirs(1) = h_dir(1)
         icyc=2
         do ient=2,num_entry
            if ( cycle_num(ient).ne.cycle_num(ient-1) ) then
               cycleSecs(icyc) = CycleTime(ient)
               cycleLons(icyc) = rlon(ient)
               cycleLats(icyc) = rlat(ient)
               cycleDirs(icyc) = h_dir(ient)
               icyc = icyc + 1
            end if
         end do
      endif
      
      ! Jie: for NWS12, need to assign output time according to OWI increments and
      ! owi_output_factor
      if (abs(nws).eq.12) then
         cycleSecs(:) = cycleTime(:) ! in seconds obtained from fort.22*
         DO icyc = 1, num_cycles
             call nws12get(wvnx2, wvny2, prn2, np, rhowat0, g)
             cycleLats(icyc) = EyeLatR(3)
             cycleLons(icyc) = EyeLonR(3)
             cycleDirs(icyc) = 0 ! not activated yet
         ENDDO
         ! write storm center locations to a file
400      FORMAT(f10.2,x,2(f8.4,x)) 
         OPEN(399,FILE='OWI_snap_storm_centers.txt',STATUS='REPLACE')
         do icyc = 1,num_cycles
            write(399,400) cycleSecs(icyc), cycleLats(icyc), 
     &         cycleLons(icyc)
         end do
         close(399)
      endif 
      
      ! Jie: for NWS15, need to assign cycleSecs according to cycleTime;
      ! outputLons and outputLats will be assigned after calling
      ! nws15get to update the position of the storm's center 
      if (abs(nws).eq.15) then
         cycleSecs(:) = cycleTime(:) ! in seconds calcu. by nws15init
      endif
      
      !
      ! if we are analyzing just the cycle times
      if ( defaultOutputIncrement.eqv..true. ) then
         outputSecs(:) = cycleSecs(:)
         if ((abs(nws).eq.19).or.(abs(nws).eq.20).or.(abs(nws).eq.12)) then
            outputLons(:) = cycleLons(:)
            outputLats(:) = cycleLats(:)
            outputDirs(:) = cycleDirs(:)
         endif
      else
         !
         ! we are analyzing at a constant time increment (e.g., hourly)
         ! that does not necessarily correspond with the cycle times
!         outputSecs(1) = CycleTime(1) + timeOffset ! hotstarttime was double-accounted for nws15
!         output_increment not activated for nws12, since a moving grid is used 
         if ((abs(nws).eq.19).or.(abs(nws).eq.20)) then
            outputSecs(1) = cycleSecs(1) + timeOffset
         elseif  ((abs(nws).eq.15).or.(abs(nws).eq.12)) then 
            outputSecs(1) = cycleSecs(1)
         endif
         outputLons(1) = cycleLons(1)
         outputLats(1) = cycleLats(1)
         outputDirs(1) = cycleDirs(1)
 
         if (abs(nws).eq.12) then
            do iout=2,num_outputs
             icyc = (iout - 1) * owi_output_factor + 1 
             outputSecs(iout) = cycleSecs(icyc)
             outputLons(iout) = cycleLons(icyc)
             outputLats(iout) = cycleLats(icyc)
             outputDirs(iout) = cycleDirs(icyc)
            enddo
         else
          icyc = 2
          do iout=2,num_outputs
            outputSecs(iout) = outputSecs(iout-1) + output_increment
            if ((abs(nws).eq.19).or.(abs(nws).eq.20)) then
               if ( outputSecs(iout).gt.cycleSecs(icyc) ) then
                  icyc = icyc + 1
               endif
               interp = (outputSecs(iout) - cycleSecs(icyc-1))
     &                  /(cycleSecs(icyc)-cycleSecs(icyc-1))
               outputLons(iout) = interp *
     &                   (cycleLons(icyc)-cycleLons(icyc-1))
     &                   + cycleLons(icyc-1)
               outputLats(iout) = interp *
     &                   (cycleLats(icyc)-cycleLats(icyc-1))
     &                   + cycleLats(icyc-1)
               outputDirs(iout) = interp *
     &                   (cycleDirs(icyc)-cycleDirs(icyc-1))
     &                   + cycleDirs(icyc-1)
            endif
          end do
         endif ! (abs(nws).eq.12)
      end if ! ( defaultOutputIncrement.eqv..true. )
       !-------------------------------------------------------------
      end subroutine initializeDiagnosticOutput
      !-------------------------------------------------------------


      !-------------------------------------------------------------
      ! S U B R O U T I N E   F U L L   C I R C L E   R M A X E S
      !-------------------------------------------------------------
      subroutine fullCircleRmaxes()
      use wind, only : draglawstring
      implicit none
      integer :: i, j, npp
      real(sz) :: vmax
      real(sz) :: dr
      real(sz) :: rhowatg
      
      
      RhoWatG = RhoWat0 * g
      moving_grid=.true.
      if (abs(nws).eq.12) then
        np = 360*360 ! cannot call nws12get repreatedly
        npp=360
      else
        np = 360
      endif  
      call allocPseudoMeshArrays()
      
      if (abs(nws).eq.12) then  ! Jie 2013/09
         DragLawString='Powell' ! Activate Powell Drag Law to obtain storm centers
         GBLINPUTDIR="." ! nws12init will open GBLINPUTDIR/fort.22
         slam(:) = -79.0d0*deg2rad ! initialize to avoid spurious error messages
         sfea(:) = 34.0d0*deg2rad
         ! initialize arrays
!         wvnx2 = 0.d0, wvny2 = 0.d0, prn2 = 0.d0
         call nws12init(wvnx2, wvny2, prn2, np, RhoWat0, g) !call this to skip unwanted snaps        
      endif
      if (abs(nws).eq.15) then
         mnp = np
         globaldir="." ! nws15init will open globaldir/fort.22 
         call nws15init(timeOffset)
      endif
      call initializeDiagnosticOutput()
      
      write(*,*) ' Initialization Finished Successfully'
      ! for nws12, need to rewind wind file and start over again
      if (abs(nws).eq.12) then  ! Jie 2013/09
         close(221)
         close(222)
         close(223)
         close(224)
         
         slam(:) = -79.0d0*deg2rad ! initialize to avoid spurious error messages
         sfea(:) = 34.0d0*deg2rad
         ! initialize arrays
         ! wvnx2 = 0.d0, wvny2 = 0.d0, prn2 = 0.d0
         call nws12init(wvnx2, wvny2, prn2, np, RhoWat0, g) !call this to skip unwanted snaps        
      endif
      write(*,*) ' Begin to process diagnostic outputs for NWS = ', nws
      ! write storm center locations to a file
402   FORMAT(f10.2,x,2(f8.4,x)) 
      OPEN(401,FILE='Output_storm_centers.txt',STATUS='REPLACE')
      open(1044, file='TC_centerloc.txt',status='replace')
            
      do iout=1,num_outputs
         write(*,*) 'iout = ', iout
         ! if we only want one output frame, and this is not it,
         ! just skip to the next one
         if ( (frame.ne.0).and.(iout.ne.frame) ) then
            cycle
         endif
         if ( abs(nws).eq.15 ) then
            ! call this just to get the updated position of the
            ! storm's center so we can radiate out from that location
            slam(:) = -79.0d0*deg2rad ! initialize to avoid spurious error messages
            sfea(:) = 34.0d0*deg2rad
            call nws15get(wvnx2, wvny2, prn2, outputSecs(iout))
            outputLats(iout) = EyeLatR(3)
            outputLons(iout) = EyeLonR(3)
            write(401,402) OutputSecs(iout), outputLats(iout), 
     &         outputLons(iout)
         endif

         ! open the file where the full circle Rmaxes will be written
         write(rmax_file_name,'(I4.4,".d")') iout
         rmax_file_name = "full_circle_rmaxes_"//
     &                           rmax_file_name//".d"
         OPEN(444,FILE=trim(rmax_file_name),STATUS='REPLACE')
         write(node_file_name,'(I4.4,".d")') iout
         node_file_name = "full_circle_latlon_"//
     &                           node_file_name//".d"
         ! create a set of nodes in a circle around the center of the
         ! storm so that we can get an Rmax in all directions
         do i=1,np
            azimuth = real(i-1) * deg2rad
            slam(i) = outputLons(iout)*deg2rad
     &               + (1.0d0*deg2rad)*cos(azimuth)
            sfea(i) = outputLats(iout)*deg2rad
     &               + (1.0d0*deg2rad)*sin(azimuth)
         enddo
         !
         ! write circular lat/lon values to a file for plotting
         OPEN(333,FILE=trim(node_file_name),STATUS='REPLACE')
         do i=1,np
            write(333,*) slam(i)*rad2deg, sfea(i)*rad2deg
         end do
         close(333)
         !
443      FORMAT(5(f10.4,x)) 
         select case(abs(nws))
         case(12)
         ! Calculate and write the Rmax data to the file here Jie 11/2013
         ! Prepare a circular mesh first
             do i=1,npp
               azimuth = real(i-1) * deg2rad ! 1 deg resolution 0->359 deg
               do j=1,npp ! 1 nm radial resolution out to 360 nm
                  dr = real(j) * 1852.0d0 ! convert nmiles to meters
                  call xy2latlon(dr*cos(azimuth),
     &                 dr*sin(azimuth), outputLats(iout),
     &                 outputLons(iout), 
     &                 sfea(j+(i-1)*npp), slam(j+(i-1)*npp))
               end do               
            end do  
               
               slam(:) = slam(:) * deg2rad
               sfea(:) = sfea(:) * deg2rad + outputLats(iout) * deg2rad

               if (iout.eq.1) then
                  call nws12get(wvnx2, wvny2, prn2, np, RhoWat0, g) 
               else  
                  do j=1, owi_output_factor
                     call nws12get(wvnx2, wvny2, prn2, np, RhoWat0, g) 
                  end do
               endif
               ! if this run is for wind analysis
               ! write out Rmax (nautical miles) and azimuthal angle (degrees),
               ! adapted from nws15get   Jie 11/2013
               if (writeFullCircleRmaxes.eqv..true.) then
                  
                  do i=1,npp
                     azimuth = real(i-1) * deg2rad ! 1 deg resolution 0->359 deg
                     j=1+(i-1)*npp
                     tempspd(:)= sqrt(wvnx2(j:j+npp-1)**2 +
     &                                wvny2(j:j+npp-1)**2)
                     maxSpeed = maxval(tempspd)
                     maxNode = maxloc(tempspd,1)
                     maxNode = maxNode + j-1
                     
                     dist = m2nm * sphericalDistance(slam(maxNode)-
     &                      deg2rad*outputLons(iout),
     &                      sfea(maxNode)-deg2rad*outputLats(iout),
     &                      outputLats(iout), rad2deg*sfea(maxNode))
                  
                     angle = 360.0d0 + rad2deg *
     &               ATAN2((slam(maxNode)-deg2rad*outputLons(iout)),
     &              (sfea(maxNode)-deg2rad*outputLats(iout)))
         
                    if ( angle > 360.d0) angle = angle - 360.d0
                       write(444,443) dist, angle , real(maxNode-j+1),
     &                 azimuth*rad2deg, (ten2one * maxSpeed /0.514444d0)
                  enddo
               endif         
         case(15)
           ! nws15get actually writes the Rmax data to the file
            do i=1,np
               azimuth = real(i-1) * deg2rad ! 1 deg resolution 0->3569 deg
               do j=1,np ! 1 nm radial resolution out to 360 nm
                  dr = real(j) * 1852.0d0 ! convert nmiles to meters
                  call xy2latlon(dr*cos(azimuth),
     &                 dr*sin(azimuth), outputLats(iout),
     &                 outputLons(iout), sfea(j), slam(j))
               end do
               slam(:) = slam(:) * deg2rad
               sfea(:) = sfea(:) * deg2rad + outputLats(iout) * deg2rad
               call nws15get(wvnx2,wvny2,prn2,outputsecs(iout))
            end do
         case(19)
            ! nws19get actually writes the Rmax data to the file
            call nws19get(slam,sfea,wvnx2,wvny2,prn2,np,
     &             outputsecs(iout), ics)    
         case(20)
            ! nws20get actually writes the Rmax data to the file
            call nws20get(slam,sfea,wvnx2,wvny2,prn2,np,
     &             outputsecs(iout), ics)
         case default
            write(*,*) "ERROR: aswip: -n ",nws,
     &           " is not an nws value supported by aswip."
         end select
         CLOSE(444)
         ! search around the storm along Rmax to find the vmax
         ! as well as the azimuth where vmax occurs
         vmax = 0.d0
         do i=1,np
            vtest = sqrt(wvnx2(i)**2 + wvny2(i)**2)
     &                       /0.514444d0
            if ( vtest.gt.vmax ) then
               vmax = vtest
               vmax_azimuth = real(i-1)
            endif
         end do
      end do
      !-------------------------------------------------------------
      end subroutine fullCircleRmaxes
      !-------------------------------------------------------------


    !-------------------------------------------------------------
      ! S U B R O U T I N E    R A D I A L   V   A N D   P
      !-------------------------------------------------------------
      subroutine radialVandP()
      use wind, only : angleNow, draglawstring
      implicit none
      integer :: i, j
      real(sz) :: dr
      real(sz) :: rhowatg
      !
      RhoWatG = RhoWat0 * g
      moving_grid=.true.
      ! 4 quadrants + max
      np = radial_extent * quadrant_radial_resolution * 5
      call allocPseudoMeshArrays()
      if (abs(nws).eq.12) then  ! Jie 2013/09
         DragLawString='Powell' ! Activate Powell Drag Law to obtain storm centers
         GBLINPUTDIR="." ! nws12init will open GBLINPUTDIR/fort.22
         slam(:) = -79.0d0*deg2rad ! initialize to avoid spurious error messages
         sfea(:) = 34.0d0*deg2rad
         ! initialize arrays
!         wvnx2 = 0.d0, wvny2 = 0.d0, prn2 = 0.d0
         call nws12init(wvnx2, wvny2, prn2, np, RhoWat0, g) !call this to skip unwanted snaps        
      endif
      if (abs(nws).eq.15) then
         mnp = np
         globaldir="." ! nws15init will open globaldir/fort.22
         call nws15init(timeOffset)
      endif
      call initializeDiagnosticOutput()
      write(*,*) ' Initialization Finished Successfully'
      ! for nws12, need to rewind wind file and start over again
      if (abs(nws).eq.12) then  ! Jie 2013/09
         close(221)
         close(222)
         close(223)
         close(224)
         
         slam(:) = -79.0d0*deg2rad ! initialize to avoid spurious error messages
         sfea(:) = 34.0d0*deg2rad
         ! initialize arrays
         ! wvnx2 = 0.d0, wvny2 = 0.d0, prn2 = 0.d0
         call nws12init(wvnx2, wvny2, prn2, np, RhoWat0, g) !call this to skip unwanted snaps        
      endif
      
      write(*,*) ' Begin to process diagnostic outputs for NWS = ', nws
      ! write storm center locations to a file
402      FORMAT(f10.2,x,2(f8.4,x)) 
         OPEN(401,FILE='Output_storm_centers.txt',STATUS='REPLACE')
      do iout=1,num_outputs
         write(*,*) 'iout = ', iout
         ! if we only want one output frame, and this is not it,
         ! just skip to the next one
         if ( (frame.ne.0).and.(iout.ne.frame) ) then
            cycle
         endif
         
         if ( abs(nws).eq.15 ) then
            ! call this just to get the updated position of the
            ! storm's center so we can form a mesh at that location
            slam(:) = -79.0d0*deg2rad ! initialize to avoid spurious error messages
            sfea(:) = 34.0d0*deg2rad
            call nws15get(wvnx2, wvny2, prn2, outputSecs(iout))
            outputLats(iout) = EyeLatR(3)
            outputLons(iout) = EyeLonR(3)
         endif
            write(401,402) OutputSecs(iout), outputLats(iout), 
     &         outputLons(iout)
         ! create four strings of nodes in straight lines out from
         ! the center of the storm, one in the center of each quadrant
         dr = 0.0d0
         do i=1,np-5,5
            call xy2latlon(dr*cos(45.d0*deg2rad),
     &                 dr*sin(45.d0*deg2rad), outputLats(iout),
     &                 outputLons(iout), sfea(i), slam(i)) ! NEQ
            call xy2latlon(dr*cos(315.d0*deg2rad),
     &                 dr*sin(315.0d0*deg2rad), outputLats(iout),
     &                 outputLons(iout), sfea(i+1), slam(i+1)) ! SEQ
            call xy2latlon(dr*cos(225.d0*deg2rad),
     &                 dr*sin(225.d0*deg2rad), outputLats(iout),
     &                 outputLons(iout), sfea(i+2),slam(i+2))  ! SWQ
            call xy2latlon(dr*cos(135.d0*deg2rad),
     &                 dr*sin(135.d0*deg2rad),outputLats(iout),
     &                 outputLons(iout), sfea(i+3), slam(i+3)) ! NWQ
            select case(abs(nws))
            case(12) !Not activated yet 09/2013 Jie
               ! compute the translation direction as trigonometric
               ! angle (as opposed to a compass angle)
               trigDir = 360.d0 - outputDirs(iout) + 90.d0
               ! the Vmax should occur in a direction that is
               ! 90 degrees to CCW from this angle
               trigDir = trigDir - 90.d0
            case(15)
               ! this is the angle to the Vmax
               trigDir = angleNow
            case(19)
               ! compute the translation direction as trigonometric
               ! angle (as opposed to a compass angle)
               trigDir = 360.d0 - outputDirs(iout) + 90.d0
               ! the Vmax should occur in a direction that is
               ! 90 degrees to CCW from this angle
               trigDir = trigDir - 90.d0
            case(20)
               ! compute the translation direction as trigonometric
               ! angle (as opposed to a compass angle)
               trigDir = 360.d0 - outputDirs(iout) + 90.d0
               ! the Vmax should occur in a direction that is
               ! 90 degrees to CCW from this angle
               trigDir = trigDir - 90.d0
            end select
            call xy2latlon(
     &                 dr*cos(trigDir*deg2rad),
     &                 dr*sin(trigDir*deg2rad),
     &                 outputLats(iout),outputLons(iout),
     &                 sfea(i+4), slam(i+4)) ! MAX
            dr = dr + 1852.0d0/real(quadrant_radial_resolution)
         end do
         slam(:) = slam(:) * deg2rad
         sfea(:) = sfea(:) * deg2rad + outputLats(iout) * deg2rad
         write(node_file_name,'(I4.4,".d")') iout
         node_file_name = "radialvp_latlon_"//
     &                           node_file_name//".d"
         ! write lat/lon values to a file for plotting
         OPEN(333,FILE=trim(node_file_name),STATUS='REPLACE')
         do i=1,np
            write(333,*) slam(i)*rad2deg, sfea(i)*rad2deg
         end do
         close(333)
         select case(abs(nws))
         case(12)
            ! output each frame only, no jumping, no time-interpolation
            if (iout.eq.1) then
              call nws12get(wvnx2, wvny2, prn2, np, RhoWat0, g) 
            else  
              do j=1, owi_output_factor
                 call nws12get(wvnx2, wvny2, prn2, np, RhoWat0, g) 
              end do
             endif  
         case(15)
            call nws15get(wvnx2, wvny2, prn2, outputSecs(iout))
         case(19)
            call nws19get(slam,sfea,wvnx2,wvny2,prn2,np,
     &             outputsecs(iout), ics)
         case(20)
            ! nws20get actually writes the Rmax data to the file
            call nws20get(slam,sfea,wvnx2,wvny2,prn2,np,
     &             outputsecs(iout), ics)
         case default
            write(*,*) "ERROR: aswip: -n ",nws,
     &           " is not an nws value supported by aswip."
         end select
         ! open the file where the radial V and P will be written
         write(radialvp_file_name,'(I4.4,".d")') iout
         radialvp_file_name = "radialvp_"//
     &                           radialvp_file_name//".d"
         OPEN(555,FILE=trim(radialvp_file_name),STATUS='REPLACE')
         ! convert from m/s velocity to speed in kt and convert
         ! from 10minute averaged winds (needed by ADCIRC) back
         ! to 1minute averaged winds (which match the data from NHC)
         do i=1,np
            radial_speed(i)
     &               = ten2one * sqrt(wvnx2(i)**2 + wvny2(i)**2)
     &                       /0.514444d0
            ! convert from mH2O to mbar
            pressure(i) = prn2(i) * RhoWatG / 100.d0
         enddo
         dr = 0.0d0
         do i=1,np-5,5
            write(555,*) dr,(radial_speed(i+j),j=0,4),
     &                         (pressure(i+j),j=0,4)
            dr = dr + 1.0d0/real(quadrant_radial_resolution)
         end do
         CLOSE(555)
      end do
      close(401)
      !-------------------------------------------------------------
      end subroutine radialVandP
      !-------------------------------------------------------------


      !-------------------------------------------------------------
      ! S U B R O U T I N E   S P A T I A L   V   A N D   P
      !-------------------------------------------------------------
      subroutine spatialVandP()
      use global, only : nm, ne, nhy, slam0, sfea0
      implicit none
      integer :: i,j
      integer :: errorio
      real(sz) :: dr
      real(sz) :: dx
      real(sz) :: dy
      real(sz) :: rhowatg

      RhoWatG = RhoWat0 * g
      if (nodesFromMesh.eqv..false.) then
         np = radial_extent/spatial_radial_increment
     &              * (360/azimuthal_increment)
      else
         call openFileForRead(14,trim(meshFileName),errorIO)
         if (errorIO.ne.0) stop
         read(14,*) scratch  ! throw away AGRID
         read(14,*) ne, np
      endif

      call allocPseudoMeshArrays()
      if (abs(nws).eq.15) then
         mnp = np
         globaldir="."
         call nws15init(timeOffset)
      endif
      call initializeDiagnosticOutput()

      do iout=1,num_outputs
         ! if we only want one output frame, and this is not it,
         ! just skip to the next one
         if ( (frame.ne.0).and.(iout.ne.frame) ) then
            cycle
         endif
         if ( abs(nws).eq.15 ) then
            ! call this just to get the updated position of the
            ! storm's center so we can form a mesh at that location
            slam(:) = -79.0d0*deg2rad ! initialize to avoid spurious error messages
            sfea(:) = 34.0d0*deg2rad
            call nws15get(wvnx2, wvny2, prn2, outputSecs(iout))
            outputLats(iout) = EyeLatR(3)
            outputLons(iout) = EyeLonR(3)
         endif
         if (nodesFromMesh.eqv..false.) then
            ! create a set of nodes in circles around the center of the
            ! storm so that we can get visualize V and P in all directions
            dr = real(spatial_radial_increment)
            num_rvalues = radial_extent/spatial_radial_increment
            num_avalues = 360/azimuthal_increment
            i=1
            ! if we weren't provided with an element table, we
            ! have to write out the nodes to a .nod formatted file
            ! so we an element table can be constructed with an
            ! external program, for use in visualization of
            ! wind velocities (can't do particle tracing without
            ! elements to interpolate over)
            if ( elementTableProvided.eqv..false.) then
               allocate(nodex(np),nodey(np))
               open(777,file="spatial.nod",status="replace")
               write(777,'(i12,A)') np, "  ! NC, number of nodes"
            endif
            ! create outer boundary
            do t=1, num_avalues
               azimuth = real(t*azimuthal_increment) * deg2rad
               dx = real(radial_extent)*cos(azimuth) * 1852.0d0 ! nm to m
               dy = real(radial_extent)*sin(azimuth) * 1852.0d0 ! nm to m
               call xy2latlon(dx, dy,
     &                     outputLats(iout), outputLons(iout),
     &                     sfea(i), slam(i))
               if ( elementTableProvided.eqv..false. ) then
                  nodex(i) = dx
                  nodey(i) = dy
               endif
               i=i+1
            enddo
            ! create inner boundary (must be counter clockwise)
            do t=num_avalues,1,-1
               azimuth = real(t*azimuthal_increment) * deg2rad
               dx = real(spatial_radial_increment)*cos(azimuth)
     &                   * 1852.0d0 ! nm to m
               dy = real(spatial_radial_increment)*sin(azimuth)
     &                   * 1852.0d0 ! nm to m
               call xy2latlon(dx, dy,
     &                     outputLats(iout), outputLons(iout),
     &                     sfea(i), slam(i))
               if ( elementTableProvided.eqv..false. ) then
                  nodex(i) = dx
                  nodey(i) = dy
               endif
               i=i+1
            enddo
            ! create interior (i.e., non-boundary) nodes
            do ry=2, num_rvalues-1
               do t=1, num_avalues
                  azimuth = real(t*azimuthal_increment) * deg2rad
                  dx =
     &                real(ry*spatial_radial_increment)*cos(azimuth)
     &                     * 1852.0d0 ! nm to m
                  dy =
     &                    real(ry*spatial_radial_increment)
     &                     * sin(azimuth) * 1852.0d0 ! nm to m
                  call xy2latlon(dx, dy,
     &                     outputLats(iout), outputLons(iout),
     &                     sfea(i), slam(i))
                  if ( elementTableProvided.eqv..false. ) then
                     nodex(i) = dx
                     nodey(i) = dy
                  endif
                  i=i+1
               enddo
            enddo
            if ( i.ne.(np+1) ) then
               write(*,*) 'ERROR: aswip: i=',i,' while np=',np
            endif
            slam0 = outputLons(iout)
            sfea0 = outputLats(iout)
            slam(:) = slam(:) * deg2rad
            sfea(:) = sfea(:) * deg2rad
     &                 + outputLats(iout) * deg2rad
            ! now write the .nod file if needed
            if (elementTableProvided.eqv..false.) then
               write(777,'(A)') "2   ! NB, number of boundaries"
               ! start with the outer boundary
               write(777,'(I3,A)') num_avalues,
     &                  "   ! num nodes on outer boundry"
               do i=1,num_avalues
                  write(777,'(2(E15.8,2X),A)')
     &                     nodex(i), nodey(i)," 0.0"
               enddo
               ! now the inner boundary
               write(777,'(I3,A)') num_avalues,
     &                   "      ! num nodes on inner boundry"
               do i=num_avalues+1,num_avalues*2
                  write(777,'(2(E15.8,2X),A)')
     &                     nodex(i), nodey(i), " 0.0"
               enddo
               ! then write the list of internal (nonboundary) nodes
               write(777,'(I12)') (np-2*num_avalues)
               do i=(2*num_avalues+1),np
                  write(777,'(2(E15.8,2X),A)')
     &                     nodex(i), nodey(i)," 0.0"
               enddo
               close(777)
            else
               ! we have been provided with an element table
               if ( meshFileReadComplete.eqv..false.) then
                  ! read in the element table (file only has ele table)
                  call openFileForRead(14,elementTableFileName,
     &                      errorIO)
                  if (errorIO.ne.0) stop
                  ! first count the number of elements in the file
                  ne=0
                  do
                     read(14,'(I8)',end=9295) i
                     ne=ne+1
                  enddo
9295                   continue
                  rewind(14)
                  ! now allocate memory and read in the element table
                  allocate(nm(ne,3)) !TODO: this effectively hardcodes nhy=3
                  cell_list_size=0
                  nhy=3
                  do i=1,ne
                     read(14,'(3I8)') nm(i,1),nm(i,2),nm(i,3)
                     cell_list_size = cell_list_size + 1 + nhy
                  enddo
                  close(14)
                  meshFileReadComplete=.true.
               endif
            endif
         else
            ! we are reading the nodal locations from an ADCIRC mesh
            ! file (fort.14)
            if (meshFileReadComplete.eqv..false.) then
               do i=1,np
                  read(14,*) j, dx, dy
                  slam(j) = dx
                  sfea(j) = dy
               enddo
               ! read the element table, just so we can write it out in
               ! vtk format to create streamlines
               allocate(nm(ne,3)) ! TODO: this effectively hardcodes nhy=3
               cell_list_size=0
               do i=1,ne
                  read(14,*) j,nhy,nm(j,1),nm(j,2),nm(j,3)
                  cell_list_size = cell_list_size + 1 + nhy
               enddo
               slam0=outputLons(iout)
               sfea0=outputLats(iout)
               meshFileReadComplete=.true.
               slam(:) = slam(:) * deg2rad
               sfea(:) = sfea(:) * deg2rad
            endif
         endif

         ! open the file where spatial coordinates for V and P
         ! will be written
         write(spatial_file_name,'(I4.4,".d")') iout
         spatial_file_name = "spatial_data_"//
     &                           spatial_file_name//".d"
         open(555,FILE=trim(spatial_file_name),STATUS='REPLACE')
         write(555,'(A)') '# vtk DataFile Version 3.0'
         write(555,'(A)') 'Vortex V and P data'
         write(555,'(A)') 'ASCII'
         !write(555,'(A)') 'DATASET STRUCTURED_GRID'
         write(555,'(A)') 'DATASET UNSTRUCTURED_GRID'
         !write(555,'(A,I3,I3,I3)') 'DIMENSIONS 360 400 1'
         write(555,'(A,I6,A)') 'POINTS ',np,' float'
         ! save coordinates as Cartesian (m) for use in VTK
         xmin=huge(1.d0)
         ymin=huge(1.d0)
         xmax=tiny(1.d0)
         ymax=tiny(1.d0)
         do i=1,np
            call latlon2xy(sfea(i)*rad2deg-sfea0,
     &              slam(i)*rad2deg, sfea0, slam0, dx, dy)
            write(555,'(E15.8,2X,E15.8,2X,E15.8)')
     &              dx, dy, 0.0d0
            if ( dx.gt.xmax ) xmax=dx
            if ( dx.lt.xmin ) xmin=dx
            if ( dy.gt.ymax ) ymax=dy
            if ( dy.lt.ymin ) ymin=dy
         enddo
         !write(*,*) "DEBUG: xmin=",xmin," xmax=",xmax
         !write(*,*) "DEBUG: ymin=",ymin," ymax=",ymax
         if ((nodesFromMesh.eqv..false.).and.
     &             (elementTableProvided.eqv..false.)) then
            ! we write this out as an unstructured point set for vtk
            write(555,'(A,I6,I7)') 'CELLS ',np,np*2
            do i=1,np
               write(555,'(A,I6)') '1 ',i
            enddo
            write(555,'(A,I6)') 'CELL_TYPES ',np
            do i=1,np
               write(555,'(A)') '1'
            enddo
         else
            ! write out cells (triangular elements), provided by
            ! an ADCIRC mesh file or an element table file
            write(555,'(A,I6,2X,I7)') 'CELLS ',ne,cell_list_size
            do i=1,ne
               ! vtk zero-indexes the node numbers, so we have
               ! to subtract 1 from the adcirc node numbers, b/c
               ! they are indexed starting at 1
               write(555,'(A,3(I6,2X))') '3 ',
     &                  nm(i,1)-1, nm(i,2)-1, nm(i,3)-1
            enddo
            write(555,'(A,I6)') 'CELL_TYPES ',ne
            do i=1,ne
               write(555,'(A)') '5'
            enddo
         endif
         select case(abs(nws))
         case(15)
            call nws15get(wvnx2,wvny2,prn2,outputSecs(iout))
         case(19)
             call nws19get(slam,sfea,wvnx2,wvny2,prn2,np,
     &             outputsecs(iout), ics)
         case default
            write(*,*) "ERROR: aswip: -n ",nws,
     &           " is not an nws value supported by aswip."
         end select
         write(555,'(A,I6)') 'POINT_DATA ',np
         write(555,'(A)') 'SCALARS WindSpeed float 1'
         write(555,'(A)') 'LOOKUP_TABLE default'
         do i=1,np
            write(555,'(E15.8)')
     &              (ten2one * sqrt(wvnx2(i)**2 + wvny2(i)**2)
     &                       /0.514444d0)
         enddo
         write(555,'(A)') 'SCALARS Pressure float 1'
         write(555,'(A)') 'LOOKUP_TABLE default'
         do i=1,np
            write(555,'(E15.8)') prn2(i) * RhoWatG / 100.d0
         enddo
         write(555,'(A)') 'VECTORS WindVelocity float'
         do i=1,np
            write(555,'(E15.8,E15.8,E15.8)')
     &              (ten2one*wvnx2(i)/0.514444d0),
     &              (ten2one*wvny2(i)/0.514444d0), 0.0d0
         enddo
         close(555)
      end do
      !-------------------------------------------------------------
      end subroutine spatialVandP
      !-------------------------------------------------------------

   !-------------------------------------------------------------
      ! S U B R O U T I N E   S P A T I A L   V   A N D   P
      ! 
      ! A simple version of spatialVandP to output 
      ! (id, dx, dy, u, v, p) to plot in MATLAB
      ! Added NWS=12 option       Jie 10/2013
      !-------------------------------------------------------------
      subroutine spatialuvp()
      use global, only : nm, ne, nhy, slam0, sfea0
      use wind, only : draglawstring
      implicit none
      integer :: i,j
      integer :: errorio
      real(sz) :: dr
      real(sz) :: dx
      real(sz) :: dy
      real(sz) :: rhowatg

      RhoWatG = RhoWat0 * g
      moving_grid=.true.
      if (nodesFromMesh.eqv..false.) then
         np = radial_extent/spatial_radial_increment
     &              * (360/azimuthal_increment)
      else
         call openFileForRead(14,trim(meshFileName),errorIO)
         if (errorIO.ne.0) stop
         read(14,*) scratch  ! throw away AGRID
         read(14,*) ne, np
      endif

      call allocPseudoMeshArrays()
      if (abs(nws).eq.12) then  ! Jie 2013/09
         DragLawString='Powell' ! Activate Powell Drag Law to obtain storm centers
         GBLINPUTDIR="." ! nws12init will open GBLINPUTDIR/fort.22
         slam(:) = -79.0d0*deg2rad ! initialize to avoid spurious error messages
         sfea(:) = 34.0d0*deg2rad
         ! initialize arrays
!         wvnx2 = 0.d0, wvny2 = 0.d0, prn2 = 0.d0
         call nws12init(wvnx2, wvny2, prn2, np, RhoWat0, g) !call this to skip unwanted snaps        
      endif
      if (abs(nws).eq.15) then
         mnp = np
         globaldir="."
         call nws15init(timeOffset)
      endif
      call initializeDiagnosticOutput()
            write(*,*) ' Initialization Finished Successfully'
      ! for nws12, need to rewind wind file and start over again
      if (abs(nws).eq.12) then  ! Jie 2013/09
         close(221)
         close(222)
         close(223)
         close(224)
         
         slam(:) = -79.0d0*deg2rad ! initialize to avoid spurious error messages
         sfea(:) = 34.0d0*deg2rad
         ! initialize arrays
         ! wvnx2 = 0.d0, wvny2 = 0.d0, prn2 = 0.d0
         call nws12init(wvnx2, wvny2, prn2, np, RhoWat0, g) !call this to skip unwanted snaps        
      endif
      
      write(*,*) ' Begin to process diagnostic outputs for NWS = ', nws
      ! write storm center locations to a file
402   FORMAT(f10.2,x,2(f8.4,x)) 
      OPEN(401,FILE='Output_storm_centers.txt',STATUS='REPLACE')
      open(1044, file='TC_centerloc.txt',status='replace')
      do iout=1,num_outputs
         write(*,*) 'iout = ', iout
         ! if we only want one output frame, and this is not it,
         ! just skip to the next one
         if ( (frame.ne.0).and.(iout.ne.frame) ) then
            cycle
         endif
         if ( abs(nws).eq.15 ) then
            ! call this just to get the updated position of the
            ! storm's center so we can form a mesh at that location
            slam(:) = -79.0d0*deg2rad ! initialize to avoid spurious error messages
            sfea(:) = 34.0d0*deg2rad
            call nws15get(wvnx2, wvny2, prn2, outputSecs(iout))
            outputLats(iout) = EyeLatR(3)
            outputLons(iout) = EyeLonR(3)
         endif
         write(401,402) OutputSecs(iout), outputLats(iout), 
     &         outputLons(iout)
         if (nodesFromMesh.eqv..false.) then
            ! create a set of nodes in circles around the center of the
            ! storm so that we can get visualize V and P in all directions
            dr = real(spatial_radial_increment)
            num_rvalues = radial_extent/spatial_radial_increment
            num_avalues = 360/azimuthal_increment
            i=1
            ! if we weren't provided with an element table, we
            ! have to write out the nodes to a .nod formatted file
            ! so an element table can be constructed with an
            ! external program, for use in visualization of
            ! wind velocities (can't do particle tracing without
            ! elements to interpolate over)
            if ( elementTableProvided.eqv..false.) then
               allocate(nodex(np),nodey(np))
               open(777,file="spatial.nod",status="replace")
               write(777,'(i12,A)') np, "  ! NC, number of nodes"
            endif
            ! create outer boundary
            do t=1, num_avalues
               azimuth = real(t*azimuthal_increment) * deg2rad
               dx = real(radial_extent)*cos(azimuth) * 1852.0d0 ! nm to m
               dy = real(radial_extent)*sin(azimuth) * 1852.0d0 ! nm to m
               call xy2latlon(dx, dy,
     &                     outputLats(iout), outputLons(iout),
     &                     sfea(i), slam(i))
               if ( elementTableProvided.eqv..false. ) then
                  nodex(i) = dx
                  nodey(i) = dy
               endif
               i=i+1
            enddo
            ! create inner boundary (must be counter clockwise)
            do t=num_avalues,1,-1
               azimuth = real(t*azimuthal_increment) * deg2rad
               dx = real(spatial_radial_increment)*cos(azimuth)
     &                   * 1852.0d0 ! nm to m
               dy = real(spatial_radial_increment)*sin(azimuth)
     &                   * 1852.0d0 ! nm to m
               call xy2latlon(dx, dy,
     &                     outputLats(iout), outputLons(iout),
     &                     sfea(i), slam(i))
               if ( elementTableProvided.eqv..false. ) then
                  nodex(i) = dx
                  nodey(i) = dy
               endif
               i=i+1
            enddo
            ! create interior (i.e., non-boundary) nodes
            do ry=2, num_rvalues-1
               do t=1, num_avalues
                  azimuth = real(t*azimuthal_increment) * deg2rad
                  dx =
     &                real(ry*spatial_radial_increment)*cos(azimuth)
     &                     * 1852.0d0 ! nm to m
                  dy =
     &                    real(ry*spatial_radial_increment)
     &                     * sin(azimuth) * 1852.0d0 ! nm to m
                  call xy2latlon(dx, dy,
     &                     outputLats(iout), outputLons(iout),
     &                     sfea(i), slam(i))
                  if ( elementTableProvided.eqv..false. ) then
                     nodex(i) = dx
                     nodey(i) = dy
                  endif
                  i=i+1
               enddo
            enddo
            if ( i.ne.(np+1) ) then
               write(*,*) 'ERROR: aswip: i=',i,' while np=',np
            endif
            slam0 = outputLons(iout)
            sfea0 = outputLats(iout)
            slam(:) = slam(:) * deg2rad
            sfea(:) = sfea(:) * deg2rad
     &                 + outputLats(iout) * deg2rad
            ! write(*,*) "converting (dr,angle) to (lon,lat)"
            ! now write the .nod file if needed
            if (elementTableProvided.eqv..false.) then
               write(777,'(A)') "2   ! NB, number of boundaries"
               ! start with the outer boundary
               write(777,'(I3,A)') num_avalues,
     &                  "   ! num nodes on outer boundry"
               do i=1,num_avalues
                  write(777,'(2(E15.8,2X),A)')
     &                     nodex(i), nodey(i)," 0.0"
               enddo
               ! now the inner boundary
               write(777,'(I3,A)') num_avalues,
     &                   "      ! num nodes on inner boundry"
               do i=num_avalues+1,num_avalues*2
                  write(777,'(2(E15.8,2X),A)')
     &                     nodex(i), nodey(i), " 0.0"
               enddo
               ! then write the list of internal (nonboundary) nodes
               write(777,'(I12)') (np-2*num_avalues)
               do i=(2*num_avalues+1),np
                  write(777,'(2(E15.8,2X),A)')
     &                     nodex(i), nodey(i)," 0.0"
               enddo
               close(777)
            else  !(elementTableProvided.eqv..true.)
               ! we have been provided with an element table
               if ( meshFileReadComplete.eqv..false.) then
                  ! read in the element table (file only has ele table)
                  call openFileForRead(14,elementTableFileName,
     &                      errorIO)
                  if (errorIO.ne.0) stop
                  ! first count the number of elements in the file
                  ne=0
                  do
                     read(14,'(I8)',end=9295) i
                     ne=ne+1
                  enddo
9295                   continue
                  rewind(14)
                  ! now allocate memory and read in the element table
                  allocate(nm(ne,3)) !TODO: this effectively hardcodes nhy=3
                  cell_list_size=0
                  nhy=3
                  do i=1,ne
                     read(14,'(3I8)') nm(i,1),nm(i,2),nm(i,3)
                     cell_list_size = cell_list_size + 1 + nhy
                  enddo
                  close(14)
                  meshFileReadComplete=.true.
               endif
            endif
         else   !(nodesFromMesh.eqv..true.)
            ! we are reading the nodal locations from an ADCIRC mesh
            ! file (fort.14)
            if (meshFileReadComplete.eqv..false.) then
               do i=1,np
                  read(14,*) j, dx, dy
                  slam(j) = dx
                  sfea(j) = dy
               enddo
               ! read the element table, just so we can write it out in
               ! vtk format to create streamlines
               allocate(nm(ne,3)) ! TODO: this effectively hardcodes nhy=3
               cell_list_size=0
               do i=1,ne
                  read(14,*) j,nhy,nm(j,1),nm(j,2),nm(j,3)
                  cell_list_size = cell_list_size + 1 + nhy
               enddo
               slam0=outputLons(iout)
               sfea0=outputLats(iout)
               meshFileReadComplete=.true.
               slam(:) = slam(:) * deg2rad
               sfea(:) = sfea(:) * deg2rad
            endif
         endif

         ! open the file where nodal coordinates, and spatial uvp
         ! will be written
         write(spatial_file_name,'(I4.4,".d")') iout
         spatial_file_name = "spatial_data_"//
     &                           spatial_file_name//".d"
         open(555,FILE=trim(spatial_file_name),STATUS='REPLACE')
!         write(555,'(A)') '# vtk DataFile Version 3.0'
!         write(555,'(A)') 'Vortex V and P data'
!         write(555,'(A)') 'ASCII'
!         !write(555,'(A)') 'DATASET STRUCTURED_GRID'
!         write(555,'(A)') 'DATASET UNSTRUCTURED_GRID'
!         !write(555,'(A,I3,I3,I3)') 'DIMENSIONS 360 400 1'
!         write(555,'(A,I6,A)') 'POINTS ',np,' float'
         ! save coordinates as Cartesian (m) for use in VTK
         xmin=huge(1.d0)
         ymin=huge(1.d0)
         xmax=tiny(1.d0)
         ymax=tiny(1.d0)
         do i=1,np
            call latlon2xy(sfea(i)*rad2deg-sfea0,
     &              slam(i)*rad2deg, sfea0, slam0, dx, dy)
!            write(555,'(E15.8,2X,E15.8,2X,E15.8)')
!     &              dx, dy, 0.0d0
            if ( dx.gt.xmax ) xmax=dx
            if ( dx.lt.xmin ) xmin=dx
            if ( dy.gt.ymax ) ymax=dy
            if ( dy.lt.ymin ) ymin=dy
         enddo       
         !write(*,*) "DEBUG: ymin=",ymin," ymax=",ymax
         if ((nodesFromMesh.eqv..false.).and.
     &             (elementTableProvided.eqv..false.)) then
            ! we write this out as an unstructured point set for vtk
            write(555,'(A,I6,I7)') 'CELLS ',np,np*2
            do i=1,np
               write(555,'(A,I6)') '1 ',i
            enddo
            write(555,'(A,I6)') 'CELL_TYPES ',np
            do i=1,np
               write(555,'(A)') '1'
            enddo
         else
            ! write out cells (triangular elements), provided by
            ! an ADCIRC mesh file or an element table file
!            write(555,'(A,I6,2X,I7)') 'CELLS ',ne,cell_list_size
!            do i=1,ne
!               ! vtk zero-indexes the node numbers, so we have
!               ! to subtract 1 from the adcirc node numbers, b/c
!               ! they are indexed starting at 1
!               write(555,'(A,3(I6,2X))') '3 ',
!     &                  nm(i,1)-1, nm(i,2)-1, nm(i,3)-1
!            enddo
!           write(555,'(A,I6)') 'CELL_TYPES ',ne
!            do i=1,ne
!               write(555,'(A)') '5'
!            enddo
         endif
         select case(abs(nws))
         case(12)
            ! output each frame only, no jumping, no time-interpolation
            if (iout.eq.1) then
              call nws12get(wvnx2, wvny2, prn2, np, RhoWat0, g) 
            else  
              do j=1, owi_output_factor
                 call nws12get(wvnx2, wvny2, prn2, np, RhoWat0, g) 
              end do
            endif           
         case(15)
            call nws15get(wvnx2,wvny2,prn2,outputSecs(iout))
         case(19)
            call nws19get(slam,sfea,wvnx2,wvny2,prn2,np,
     &             outputsecs(iout), ics)
         case(20)
            ! nws20get actually writes the Rmax data to the file
            call nws20get(slam,sfea,wvnx2,wvny2,prn2,np,
     &             outputsecs(iout), ics)
         case default
            write(*,*) "ERROR: aswip: -n ",nws,
     &           " is not an nws value supported by aswip."
         end select
!         write(555,'(A,I6)') 'POINT_DATA ',np
!         write(555,'(A)') 'SCALARS WindSpeed float 1'
!         write(555,'(A)') 'LOOKUP_TABLE default'
!         do i=1,np
!            write(555,'(E15.8)')
!     &              (ten2one * sqrt(wvnx2(i)**2 + wvny2(i)**2)
!     &                       /0.514444d0)
!         enddo
!         write(555,'(A)') 'SCALARS Pressure float 1'
!        write(555,'(A)') 'LOOKUP_TABLE default'
!         do i=1,np
!            write(555,'(E15.8)') prn2(i) * RhoWatG / 100.d0
!         enddo
!         write(555,'(A)') 'VECTORS WindVelocity float'
!         do i=1,np
!            write(555,'(E15.8,E15.8,E15.8)')
!     &              (ten2one*wvnx2(i)/0.514444d0),
!     &              (ten2one*wvny2(i)/0.514444d0), 0.0d0
!         enddo
            ! Jie: output hurricane center locations for use in MATLAB script
104         FORMAT(f6.2,x,2(f6.2,x))            
            write(1044,104) outputsecs(iout)/3600.d0, 
     &       outputLons(iout), outputLats(iout)
!445         FORMAT(a7, x, 2(f6.2,x),a14,x,f6.2)            
!            write(*,445) "TC at [",  slam0, sfea0, 
!     &                   "] at fcsthour=", outputsecs(iout)/3600.d0

!446         FORMAT(a13,x,4(f6.2,x),a1)            
!            write(*,446) "Grid Range: [",minval(slam(:))*rad2deg,
!     &                                   maxval(slam(:))*rad2deg,
!     &                                   minval(sfea(:))*rad2deg,
!     &                                   maxval(sfea(:))*rad2deg, "]"
         do i=1,np
            call latlon2xy(sfea(i)*rad2deg-sfea0,
     &              slam(i)*rad2deg, sfea0, slam0, dx, dy)    
            write(555,'(E15.6,2X,E15.6,2X,E15.6,2X,E15.6,2X,E15.6)')
     &              dx, dy, ! output slam & sfea will make the image distorted
     &              (ten2one*wvnx2(i)/0.514444d0),
     &              (ten2one*wvny2(i)/0.514444d0), 
     &               prn2(i) * RhoWatG / 100.d0
         enddo
         close(555)
      end do
      close(401)
      close(1044)
      !-------------------------------------------------------------
      end subroutine spatialuvp
      
      
!-------------------------------------------------------------
      ! S U B R O U T I N E   N O D A L   V   A N D   P
      ! 
      ! A simple routine to read in nodal locations and output
      ! timeseries uvp for data comparison at meteo stations
      ! Added NWS=12 option       Jie 10/2013
      !-------------------------------------------------------------
      subroutine nodaluvp()
      use global, only : nm, ne, nhy, slam0, sfea0
      use wind, only : draglawstring
      implicit none
      integer :: i,j
      integer :: errorio
      real(sz) :: dr
      real(sz) :: dx
      real(sz) :: dy
      real(sz) :: rhowatg

      RhoWatG = RhoWat0 * g
      if (nodesFromFile.eqv..true.) then
         call openFileForRead(14,trim(meshFileName),errorIO)
         if (errorIO.ne.0) stop
         read(14,*) np  
         write(*,*) "Total station# = ", np
         call allocPseudoMeshArrays()
C      elseif (nodesFromMesh.eqv..true.) then
C         call openFileForRead(14,trim(meshFileName),errorIO)
C         if (errorIO.ne.0) stop
C         read(14,*) scratch  ! throw away AGRID
C         read(14,*) ne, np
              do i=1,np 
                  read(14,*) j, dx, dy
                  slam(i) = dx
                  sfea(i) = dy
                 write(*,*) "lon=", dx, "  lat=", dy
              enddo
              meshFileReadComplete=.true.
         ! open the file where nodal coordinates, and spatial uvp
         ! will be written
         station_output_name = "station_timeseries.d"
         open(555,FILE=trim(station_output_name),STATUS='REPLACE')
         
      elseif (nodesFromMesh.eqv..true.) then
         call openFileForRead(14,trim(meshFileName),errorIO)
         if (errorIO.ne.0) stop
         read(14,*) scratch  ! throw away AGRID
         read(14,*) ne, np  
         write(*,*) "Total number of nodes = ", np
         call allocPseudoMeshArrays()
              do i=1,np 
                  read(14,*) j, dx, dy
                  slam(i) = dx
                  sfea(i) = dy
C                  write(*,*) "lon=", dx, "  lat=", dy
              enddo
              write(*,*) "Read ADCIRC grid successfully"
              meshFileReadComplete=.true.
        write(*,*) "DEBUG: slam(np)=",slam(np)," sfea(np)=",sfea(np)
      endif
      if (abs(nws).eq.12) then  ! Jie 2013/09
         DragLawString='Powell' ! Activate Powell Drag Law to obtain storm centers
         GBLINPUTDIR="." ! nws12init will open GBLINPUTDIR/fort.22
         slam(:) = -79.0d0*deg2rad ! initialize to avoid spurious error messages
         sfea(:) = 34.0d0*deg2rad
         ! initialize arrays
!         wvnx2 = 0.d0, wvny2 = 0.d0, prn2 = 0.d0
         call nws12init(wvnx2, wvny2, prn2, np, RhoWat0, g) !call this to skip unwanted snaps        
      endif
      if (abs(nws).eq.15) then
         mnp = np
         globaldir="."
         call nws15init(timeOffset)
      endif

      call initializeDiagnosticOutput()
      write(*,*) ' Initialization Finished Successfully'
      ! for nws12, need to rewind wind file and start over again
      if (abs(nws).eq.12) then  ! Jie 2013/09
         close(221)
         close(222)
         close(223)
         close(224)
         
         slam(:) = -79.0d0*deg2rad ! initialize to avoid spurious error messages
         sfea(:) = 34.0d0*deg2rad
         ! initialize arrays
         ! wvnx2 = 0.d0, wvny2 = 0.d0, prn2 = 0.d0
         call nws12init(wvnx2, wvny2, prn2, np, RhoWat0, g) !call this to skip unwanted snaps        
      endif
      
      write(*,*) ' Begin to process diagnostic outputs for NWS = ', nws
      open(1044, file='TC_centerloc.txt',status='replace')
      
      write(*,*) "num_outputs = ", num_outputs
      do iout=1,num_outputs
         write(*,*) "iout = ", iout
         ! if we only want one output frame, and this is not it,
         ! just skip to the next one
         if ( (frame.ne.0).and.(iout.ne.frame) ) then
            cycle
         endif
         if ( abs(nws).eq.15 ) then
            ! call this just to get the updated position of the
            ! storm's center so we can form a mesh at that location
            slam(:) = -79.0d0*deg2rad ! initialize to avoid spurious error messages
            sfea(:) = 34.0d0*deg2rad
            call nws15get(wvnx2, wvny2, prn2, outputSecs(iout))
            outputLats(iout) = EyeLatR(3)
            outputLons(iout) = EyeLonR(3)
         endif
         
            ! we are reading the station locations from a text file
         if (meshFileReadComplete.eqv..true.) then
            
!               ! read the element table, just so we can write it out in
!               ! vtk format to create streamlines
!               allocate(nm(ne,3)) ! TODO: this effectively hardcodes nhy=3
!               cell_list_size=0
!               do i=1,ne
!                  read(14,*) j,nhy,nm(j,1),nm(j,2),nm(j,3)
!                  cell_list_size = cell_list_size + 1 + nhy
!               enddo
               slam0=outputLons(iout)
               sfea0=outputLats(iout)
               slam(:) = slam(:) * deg2rad
               sfea(:) = sfea(:) * deg2rad
          endif
            
         xmin=huge(1.d0)
         ymin=huge(1.d0)
         xmax=tiny(1.d0)
         ymax=tiny(1.d0)
         do i=1,np
            call latlon2xy(sfea(i)*rad2deg-sfea0,
     &              slam(i)*rad2deg, sfea0, slam0, dx, dy)
!            write(555,'(E15.8,2X,E15.8,2X,E15.8)')
!     &              dx, dy, 0.0d0
            if ( dx.gt.xmax ) xmax=dx
            if ( dx.lt.xmin ) xmin=dx
            if ( dy.gt.ymax ) ymax=dy
            if ( dy.lt.ymin ) ymin=dy
         enddo    
         write(*,*) slam0,sfea0
!         write(*,*) "DEBUG: ymin=",ymin," ymax=",ymax
         select case(abs(nws))
          case(12)
            ! output each frame only, no jumping, no time-interpolation
            if (iout.eq.1) then
              call nws12get(wvnx2, wvny2, prn2, np, RhoWat0, g) 
            else  
              do j=1, owi_output_factor
                 call nws12get(wvnx2, wvny2, prn2, np, RhoWat0, g) 
              end do
             endif 
         case(15)
            call nws15get(wvnx2,wvny2,prn2,outputSecs(iout))
         case(19)
            call nws19get(slam,sfea,wvnx2,wvny2,prn2,np,
     &             outputsecs(iout), ics)
         case(20)
            write(*,*) "outputsecs=", outputsecs(iout)
            ! nws20get actually writes the Rmax data to the file
            call nws20get(slam,sfea,wvnx2,wvny2,prn2,np,
     &             outputsecs(iout), ics) 
!            write(*,*) "Debug: process only 1000 nodes"
!            call nws20get(slam(1:1000),sfea(1:1000),wvnx2(1:1000),
!     &       wvny2(1:1000),prn2(1:1000),1000,
!     &             outputsecs(iout), ics)    
            write(*,*) "nws=", nws

         case default
            write(*,*) "ERROR: aswip: -n ",nws,
     &           " is not an nws value supported by aswip."
         end select

            ! Jie: output hurricane center locations for use in MATLAB script
105         FORMAT(f6.2,x,2(f6.2,x)) 
!            write(*,105) outputsecs(iout)/3600.d0, 
!     &       outputLons(iout), outputLats(iout)           
!            write(1044,105) outputsecs(iout)/3600.d0, 
!     &       outputLons(iout), outputLats(iout)
!445         FORMAT(a7, x, 2(f6.2,x),a14,x,f6.2)            
!            write(*,445) "TC at [",  slam0, sfea0, 
!     &                   "] at fcsthour=", outputsecs(iout)/3600.d0

!446         FORMAT(a13,x,4(f6.2,x),a1)            
!            write(*,446) "Grid Range: [",minval(slam(:))*rad2deg,
!     &                                   maxval(slam(:))*rad2deg,
!     &                                   minval(sfea(:))*rad2deg,
!     &                                   maxval(sfea(:))*rad2deg, "]"
         if (nodesFromFile.eqv..true.) then
         
                     write(555,'(f12.2,2X,150(E15.6,2X))')
     &               (outputSecs(iout)/3600.d0),
     &              ((ten2one*wvnx2(i)/0.514444d0),i=1,np),
     &              ((ten2one*wvny2(i)/0.514444d0),i=1,np), 
     &              ((prn2(i) * RhoWatG / 100.d0),i=1,np)

          elseif (nodesFromMesh.eqv..true.) then
             write(node_file_name,'(I4.4,".d")') iout
             node_file_name = "mesh_timeseries_"//
     &                           node_file_name//".d"
             ! write lat/lon values to a file for plotting
             OPEN(555,FILE=trim(node_file_name),STATUS='REPLACE')
             do i=1,np
                write(555,'(3(E15.6,2X))') 
     &                       ten2one*wvnx2(i)/0.514444d0,
     &                       ten2one*wvny2(i)/0.514444d0,    
     &                       prn2(i) * RhoWatG / 100.d0       
             end do
             close(555)  
             write(*,*) "Finish writing output for iout = ", iout
          else
           do i=1,np
             call latlon2xy(sfea(i)*rad2deg-sfea0,
     &              slam(i)*rad2deg, sfea0, slam0, dx, dy)    
             write(555,'(E15.6,2X,E15.6,2X,E15.6,2X,E15.6,2X,E15.6)')
     &              dx, dy, ! output slam & sfea will make the image distorted
     &              (ten2one*wvnx2(i)/0.514444d0),
     &              (ten2one*wvny2(i)/0.514444d0), 
     &               prn2(i) * RhoWatG / 100.d0
           enddo

          endif
     
      end do
      if (nodesFromFile.eqv..true.) close(555)
      if (nodesFromFile.eqv..true.) close(1044)
      !-------------------------------------------------------------
      end subroutine nodaluvp
      

      !-------------------------------------------------------------
      ! S U B R O U T I N E
      !     A L L O C   P S E U D O   M E S H   A R R A Y S
      !-------------------------------------------------------------
      subroutine allocPseudoMeshArrays()
      implicit none
      allocate(slam(np),sfea(np))
      allocate(wvnx2(np),wvny2(np),prn2(np))
      allocate(radial_speed(np),pressure(np))
      allocate(tempspd(360)) ! to use with nws12
      !-------------------------------------------------------------
      end subroutine allocPseudoMeshArrays
      !-------------------------------------------------------------


C-----------------------------------------------------------------------
      end module diagnosticOutput
C-----------------------------------------------------------------------


      !=================================================================
      !=================================================================
      !=================================================================
      !      =====                                           =====
      !      =====              PROGRAM ASWIP                =====
      !      =====                                           =====
      !=================================================================
      !=================================================================
      !=================================================================
      program aswip
      use asymmetricVortexProcessing
      use diagnosticOutput
      implicit none
C
C     jgf49.30: New variables to support command line options
      INTEGER :: IARGC   ! function to return number of command line arguments
      INTEGER :: ARGCOUNT ! number of command line arguments
      CHARACTER(2048) :: CMDLINEARG
      CHARACTER(2048) :: infile ! name of input file, default is "fort.22"

      CHARACTER(2) :: CMDLINEOPTION ! the option flag itself
      integer :: i

      call initLogging()
      call setMessageSource("aswip")
#if defined(GLOBAL_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
C     jgf49.30 Initialize to reasonable defaults.

      screenunit=6
      forecastRmax = 25.d0 ! default to a reasonable number, just in case
      constantRmax=.false. ! forecastRmax normally calculated from isotach data
      percentSpecified=.false.
      method=4             ! calc rMaxes based on largest available isotach
                           ! radius for each quadrant, a.k.a. the fancy method
      approach=2           ! set both Vr and quadrantVr tangential to azimuth
      nws=20
      writeNWS19=.true.
      writeNWS20=.true.

      infile="fort.22"
      vmax_percent=0.d0
      timeOffset = 0.0d0
      moving_grid=.false.
      nodesFromMesh=.false.
      nodesFromFile=.false.
      !
      ! jgf49.30: Process command line options, if any.
      ARGCOUNT = IARGC()
      IF (ARGCOUNT.ne.0) THEN
         I=0
         DO WHILE (I.lt.ARGCOUNT)
            I=I+1
            CALL GETARG(I,CMDLINEARG)
            CMDLINEOPTION = CMDLINEARG(1:2)
            SELECT CASE(CMDLINEOPTION)
            CASE("-A","-a") ! Rmw and Azimuth vals req. (360 deg/1deg incr.)
               write(screenunit,*) "INFO: Processing '",CMDLINEOPTION,
     &            "'."
               writeFullCircleRmaxes=.true.
            CASE("-V","-v") ! radial V and P values requested
               write(screenunit,*) "INFO: Processing '",CMDLINEOPTION,
     &            "'."
               writeRadialVandP=.true.
            CASE("-S","-s") ! spatial V and P unstructured mesh vals requested
               write(screenunit,*) "INFO: Processing '",CMDLINEOPTION,
     &            "'."
               writeSpatialVandP=.true.
            CASE("-B","-b") ! uvp at given node location xy
               write(screenunit,*) "INFO: Processing '",CMDLINEOPTION,
     &            "'."
               writeNodaluvp=.true.     
            CASE("-C","-c") ! spatial xy&uvp,unstructured mesh vals requested
               write(screenunit,*) "INFO: Processing '",CMDLINEOPTION,
     &            "'."
               writeSpatialuvp=.true.               
            CASE("-D","-d") ! don't write NWS19 fort.22
               write(screenunit,*) "INFO: Processing '",CMDLINEOPTION,
     &            "'."
               writeNWS19=.false.
            CASE("-E","-e")  ! providing an element table for spatial output
               I=I+1
               CALL GETARG(I,CMDLINEARG)
               write(screenunit,*) "INFO: Processing '",CMDLINEOPTION,
     &           "'."
               elementTableFileName = trim(CMDLINEARG)
               elementTableProvided=.true.
            CASE("-F","-f")  ! selecting a particular frame to render
               I=I+1
               CALL GETARG(I,CMDLINEARG)
               write(screenunit,*) "INFO: Processing '",CMDLINEOPTION,
     &            " ",trim(CMDLINEARG),"'."
               READ(CMDLINEARG,*) frame
            CASE("-G","-g")  ! selecting geostrophic or cyclostrophic
               I=I+1
               CALL GETARG(I,CMDLINEARG)
               write(screenunit,*) "INFO: Processing '",CMDLINEOPTION,
     &            " ",trim(CMDLINEARG),"'."
               READ(CMDLINEARG,*) geofactor
               IF ( geofactor.NE.1 ) THEN
                  geostrophic_switch = .false.
                  geofactor = 0
               ENDIF
            CASE("-I","-i")  ! input fort.14 file containing mesh locations
               I=I+1
               CALL GETARG(I,CMDLINEARG)
               write(screenunit,*) "INFO: Processing '",CMDLINEOPTION,
     &            " ",trim(CMDLINEARG),"'."
               meshFileName = trim(CMDLINEARG)
               nodesFromMesh = .true.
            CASE("-L","-l")  ! input text file containing station locations
               I=I+1
               CALL GETARG(I,CMDLINEARG)
               write(screenunit,*) "INFO: Processing '",CMDLINEOPTION,
     &            " ",trim(CMDLINEARG),"'."
               meshFileName = trim(CMDLINEARG)
               nodesFromFile = .true.               
            CASE("-M","-m")  ! method of selecting quadrant Rmax values
               I=I+1
               CALL GETARG(I,CMDLINEARG)
               write(screenunit,*) "INFO: Processing '",CMDLINEOPTION,
     &            " ",trim(CMDLINEARG),"'."
               READ(CMDLINEARG,*) method
            CASE("-Z","-z")  ! approach of calculating quadrantVr
               I=I+1
               CALL GETARG(I,CMDLINEARG)
               write(screenunit,*) "INFO: Processing '",CMDLINEOPTION,
     &            " ",trim(CMDLINEARG),"'."
               READ(CMDLINEARG,*) approach   
            CASE("-O","-o")  ! time increment for writing output (seconds)
               defaultOutputIncrement = .false.
               I=I+1
               CALL GETARG(I,CMDLINEARG)
               write(screenunit,*) "INFO: Processing '",CMDLINEOPTION,
     &            " ",trim(CMDLINEARG),"'."
               READ(CMDLINEARG,*) output_increment
            CASE("-N","-n")  ! nws parameter (12, 15 or 19, or 20)
               I=I+1
               CALL GETARG(I,CMDLINEARG)
               write(screenunit,*) "INFO: Processing '",CMDLINEOPTION,
     &            " ",trim(CMDLINEARG),"'."
               READ(CMDLINEARG,*) nws
            CASE("-P","-p")  ! percent change of rmax relative to calculated val
               I=I+1
               CALL GETARG(I,CMDLINEARG)
               write(screenunit,*) "INFO: Processing '",CMDLINEOPTION,
     &            " ",trim(CMDLINEARG),"'."
               READ(CMDLINEARG,*) rmax_percent
               IF ( rmax_percent.ne.100 ) THEN
                  percentSpecified = .true.
               ENDIF
            CASE("-X","-x")  ! percent change of vmax relative to given val
               I=I+1
               CALL GETARG(I,CMDLINEARG)
               write(screenunit,*) "INFO: Processing '",CMDLINEOPTION,
     &            " ",trim(CMDLINEARG),"'."
               READ(CMDLINEARG,*) vmax_percent
            CASE("-T","-t")  ! time offset (like hotstart time)
               I=I+1
               CALL GETARG(I,CMDLINEARG)
               write(screenunit,*) "INFO: Processing '",CMDLINEOPTION,
     &            " ",trim(CMDLINEARG),"'."
               READ(CMDLINEARG,*) timeOffset
               ihot = 1
            CASE("-R","-r")  ! constant forecastRmax
               I=I+1
               CALL GETARG(I,CMDLINEARG)
               READ(CMDLINEARG,*) forecastRmax
               write(screenunit,*) "INFO: Processing '",CMDLINEOPTION,
     &            " ",trim(CMDLINEARG),"'."
               IF ( forecastRmax.gt.0 ) THEN
                  constantRmax=.true.
                  write(screenunit,*)
     &               "INFO: All quadrant Rmaxes will be set to ",
     &               forecastRmax," throughout the forecast."
               ENDIF
            CASE("-W","-w") ! the name of the track file to read
                            ! ("fort.22" is default)
               I=I+1
               CALL GETARG(I,CMDLINEARG)
               write(scratchMessage,*) "Processing '",
     &            trim(CMDLINEOPTION)," ",trim(CMDLINEARG),"'."
               call screenMessage(INFO,scratchMessage)
               infile = trim(CMDLINEARG)
            CASE DEFAULT
               write(screenunit,*) "WARNING: The command line option '",
     &             CMDLINEARG(1:2),"' is not valid and will be ignored."
            END SELECT
         END DO
      ENDIF
C
C     initialize azimuth values in quadrants
      azimuth = 45.d0
      do i=1,4
         quadrantAngles(i) = deg2rad*azimuth
         azimuth = azimuth - 90.d0
      end do
C
      if ((abs(nws).eq.19).or.(abs(nws).eq.20)) then
         call processAsymmetricTrackFile(infile)
         if (( writeNWS19.eqv..true. ).and.(nws == 19) )then
            call writeAsymmetricInputFile(nws)
         else if (( writeNWS20.eqv..true. ).and.(nws == 20) )then
            call writeAsymmetricInputFile(nws)
         endif
      endif

      !
      ! If additional diagnostic output was requested, generate it
      !
      if (writeFullCircleRmaxes.eqv..true.) then
         call fullCircleRmaxes()
      endif
      if (writeRadialVandP.eqv..true.) then
         call radialVandP()
      endif
      if (writeSpatialVandP.eqv..true.) then
         call spatialVandP()
      endif
      if (writeSpatialuvp.eqv..true.) then
         call spatialuvp()
      endif
      if (writeNodaluvp.eqv..true.) then
         call nodaluvp()
      endif
      
C----------------------------------------------------------------------
      end program aswip
C----------------------------------------------------------------------
