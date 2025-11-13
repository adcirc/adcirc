!-----------------------------------------------------------------------
!  MODULE HYDROLOGY
!-----------------------------------------------------------------------
!> @author Chris Szpilka, University of Oklahoma, cmszpilka@ou.edu
!> 
!> @copyright Dr. R.A. Luettich and Dr. J.J. Westerink
!>
!> @brief This module handles the I/O and computations related to
!> adding interior hydrology (point) source terms to ADCIRC. These terms
!> allow for interior contributions from hydrologic sources that would
!> not be captured at the boundary.
!>
!> @version These source terms were originally added in v53 but not
!>  merged into the main trunk until v5#.
!>
!>  NOTE: Hydrology source terms are currently added only at "wet" ADCIRC
!>  nodes since ADCIRC cannot route over dry land - read notes below.
!>  ADCIRC will only apply the source discharge at nodes that are wet and
!>  will "ignore" the source terms on nodes that are dry at the current
!>  time step.
!>  
!>  Interior hydrology is modeled by adding point sources at specified
!>  locations within the ADCIRC domain and can be used to add hydrologic
!>  input at any interior point, however they are only applied when
!>  the node is wet. This source term should be used for all interior
!>  rivers (that do not cross a boundary) and any major tributaries that
!>  feed into them. Very small coastal features may be added as well, but
!>  it is advised not to add small tributaries and feeding streams in
!>  upland regions (above 10m topographic contour) since these features
!>  are typically not resolved within the mesh and would likely remain
!>  dry throughout the simulation. Please read the documentation for
!>  best practice when incorporating hydrology within ADCIRC domains.
!>  FULL URL TO DOCUMENTATION ON GITHUB OR ADCIRC.ORG
!>
!>  Point source terms are expected in units of discharge [L^3/T] and
!>  can be provided from any source (USGS gages, NWM, or any hydrology
!>  model). The units must be consistent with the units used for all
!>  other input in the fort.15 file. In most cases, discharges should
!>  be provided in units of [cubic meters per second].
!>
!>-----------------------------------------------------------------------
!>  USAGE:
!>  interiorHydrologyControl namelist at bottom of fort.15 activates
!>  this feature starting at V##.
!>
!>  useInteriorHydrology (logical): default .false.
!>  interiorHydrologyType (integer): default 0 (see options below)
!>  interiorHydrologyTimeInc(real): default 3600.d0 (hourly)
!>     time increment between input - controls reading and interpolation
!>  useInteriorHydrologyRamp (logical): default .false.
!>     time of ramp is set at 2 hours, user can select whether or not to
!>     apply - tests indicate that model typically remains stable even
!>     with large discharges, but option is available if user notices
!>     instability near inland hydrology input locations
!>
!>  Additional hydrology options can be included by adding more
!>  options in the fort.15 namelist and this module. Initially, two
!>  options are available:
!>
!>  interiorHydrologyType = 1:
!>  Explicit coupling at each interiorHydrologyTimeInc
!>     This option can be used for real-world simulations but was
!>     primarily coded for use with small test cases. The user
!>     must specify the ADCIRC node number(s) to apply forcing at
!>     as well as the actual discharges at each location.
!>
!>     Requires two files as input:
!>        fort.420 static location file read once: ADCIRC node#
!>        fort.430 discharge values (m3/s) at each point in location
!>            file for each interiorHydrologyTimeInc - discharges must
!>            be given in same order as the locations in fort.420
!>
!>  interiorHydrologyType = 2:
!>  Static coupling at apriori specified geographic (lon,lat) locations
!>     Static: means that each river/tributary discharge is always
!>             applied at the same location within ADCIRC
!>     Requires two files as input:
!>        fort.428 static location file read once: lon, lat
!>        fort.430 discharge values (m3/s) at each point in location file
!>            for each interiorHydrologyTimeInc - discharges must be
!>            given in same order as the locations in fort.428
!>
!>      Locations are given in fort.428 and read in once during readinput
!>      Only one location may be given for each hydrologic feature; it is
!>      recommended that this be near the main stem where ADCIRC nodes
!>      will be wet. Note that more than one feature may share the same
!>      location to allow for tributaries. The discharges for each
!>      feature are provided in fort.430 for each increment and must be
!>      given in the same order as the fort.428 location file. Discharges
!>      must be accumulated as they are read in since multiple sources
!>      could contribute to the same ADCIRC node at junctions.
!>
!>  interiorHydrologyType = 3: COMING SOON
!>  Dynamic coupling at apriori specified geographic (lon,lat) locations
!>     Dynamic: means that every river/tributary discharge is allowed
!>        to move and be applied at upstream locations as the mesh wets
!>     Requires two files as input:
!>        fort.429 dynamic location file read once: lon, lat, network
!>        fort.430 discharge values (m3/s) at each point in location file
!>            for each interiorHydrologyTimeInc - discharges must be
!>            given in same order as the locations in fort.429
!>     Locations are given in fort.429 and read in once during readinput
!>     Multiple locations may be given for each hydrologic feature.
!>     In addition to the location, the connectivity (or network) is
!>     given in the fort.429 file. See documentation for more details.
!>     It is recommended that the most downstream location be near the
!>     main stem where ADCIRC nodes will be wet (similar to option 2).
!>     ADCIRC will dynamically change to the most upstream location for
!>     each feature during runtime, so each discharge will be read in and
!>     stored separately, in contrast to the static coupling mode. More
!>     than one feature may still share the same location to allow for
!>     tributaries, but discharges for each feature will be stored in an
!>     array and applied in a separate subroutine ___________.
!>     The discharges for each feature are provided in fort.430 for each
!>     increment and must be given in the same order as fort.429 file.
!>
!>  I/O unit numbers are set as described above until these are automated
!>  in the larger ADCIRC scheme.
!>
!>----------------------------------------------------------------------- 
      MODULE Hydrology

      use global, only : screenUnit, setMessageSource,                  &
     &    unsetMessageSource, openFileForRead, scratchMessage,          &
     &    logMessage, allMessage, screenMessage, DEBUG, ECHO, INFO,     &
     &    WARNING, ERROR
      use sizes, only : inputdir, localdir, mnp

      implicit none

      ! Hydrology (point source) timing variables
      real(8) :: hydroTime1 ! time at previous hydrology input
      real(8) :: hydroTime2 ! time at new hydrology input
      real(8) :: hydroRamp=0.08333333  ! optional ramp: 2 hrs

      ! Hydrology arrays for input datasets (sized by numHydro)
      real(8), allocatable :: hydro1(:) ! previous hydrology input
      real(8), allocatable :: hydro2(:) ! new hydrology input
      real(8), allocatable :: hydroX(:) ! Current interpolated
      real(8), allocatable :: hydroOld(:) ! Previous interpolated

      ! Hydrology scalar variables
      integer :: numHydro ! number of hydrologic sources

      ! Hydrology arrays for status & location (sized by numHydro)
      integer, allocatable :: isActive(:)   ! flag for active sources
      integer, allocatable :: hydroNearNode(:)   ! nearest adcirc node

      ! Namelist input from within read_input.F 
      logical :: useInteriorHydrology = .false. ! default none
      integer :: interiorHydrologyType = 0 ! default none
      real(8) :: interiorHydrologyTimeInc = 3600 ! input time (sec)
      logical :: useInteriorHydrologyRamp = .false. ! default none
      namelist /interiorHydrologyControl/ useInteriorHydrology,         &
     &          interiorHydrologyType, InteriorHydrologyTimeInc,        &
     &          useInteriorHydrologyRamp
 
      ! Source terms used in gwce.F (sized by np)
      real(8),allocatable :: hydroCurr(:) ! hydro sources at current ts
      real(8),allocatable :: hydroODT(:) ! hydro sources time derivative

      ! Private variables
      private :: hydro1, hydro2, hydroX, hydroOld
      private :: hydroTime1, hydroTime2
      private :: hydroRamp, numHydro, hydroNearNode, isActive
      ! Private subroutines
      private :: readInteriorHydrologyLocs, kdtSearchHydrology
      private :: hydrologyTerminate, checkActiveHydrology

      ! Public variables
      public :: interiorHydrologyControl
      public :: useInteriorHydrology, interiorHydrologyType
      public :: interiorHydrologyTimeInc, useInteriorHydrologyRamp
      public :: hydroCurr, hydroODT
      ! Public subroutines
      public :: interiorHydrologyCheck, allocInteriorHydrology
      public :: getInteriorHydrology, coldstartInteriorHydrology
      public :: hotstartInteriorHydrology


!>---------------------end of data declarations--------------------------


      CONTAINS


!-----------------------------------------------------------------------
!> Called from read_input.F (when useInteriorHydrology is .true.) to 
!> check validity of hydrology point source type and log messages.
!>
!> The value of interiorHydrologyType is set during readinput for
!> the interiorHydrologyControl namelist if present
!>
!> CMS: edit when change for dynamic unit numbers and filenames
!-----------------------------------------------------------------------
      subroutine interiorHydrologyCheck()
!-----------------------------------------------------------------------

      implicit none
  
      integer :: lun1, lun2
      character(10) :: reftime

      call setMessageSource("interiorHydrologyCheck")
#if defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter")
#endif

      ! Check value of useInteriorHydrology before messages are
      ! logged, otherwise errors may arise when not being used
      ! but namelist was found in read_input.
      if (.not. useInteriorHydrology) interiorHydrologyType=0
      !if (useInteriorHydrology == .false.) interiorHydrologyType=0
      ! Log input messages related to implementation
      select case (interiorHydrologyType)
      case(0) ! no point sources
         write(16,3268) interiorHydrologyType
      case(-1) ! discharges at specified nodes (hotstart time)
         lun1=420
         lun2=430
         reftime='hotstart'
         write(16,3269) interiorHydrologyType, lun1,                    &
     &        lun2, interiorHydrologyTimeInc, lun2, trim(reftime)
      case(1) ! discharges at specified nodes (coldstart time)
         lun1=420
         lun2=430
         reftime='coldstart'
         write(16,3269) interiorHydrologyType, lun1,                    &
     &        lun2, interiorHydrologyTimeInc, lun2, trim(reftime)
      case(-2) ! static (two files) - lon/lat and Q (HS time)
         lun1=428
         lun2=430
         reftime='hotstart'
         write(16,3270) interiorHydrologyType, lun1,                    &
     &        lun2, interiorHydrologyTimeInc, lun2, trim(reftime)
      case(2) ! static (two files) - lon/lat and Q (CS time)
         lun1=428
         lun2=430
         reftime='coldstart'
         write(16,3270) interiorHydrologyType, lun1,                    &
     &        lun2, interiorHydrologyTimeInc, lun2, trim(reftime)
      case default
         write(screenUnit,3271) interiorHydrologyType
         write(16,3271) interiorHydrologyType
         call hydrologyTerminate()
      end select
  
      if (useInteriorHydrology) then
         if(useInteriorHydrologyRamp) then
            write(16,3272)
         else
            write(16,3273)
         end if
      end if


 3268   format(/,5X, 'interiorHydrologyType = ',I3,                     &
     &   /,9X,'Interior hydrology forcing will not be used in',         &
     &   ' the GWC equation') 
 3269   format(/5X,'interiorHydrologyType = ',I3,                       &
     &    /,9X,'Interior hydrology forcing will be used in the GWCE',   &
     &    /,9X,'at specified ADCIRC nodes. Nodal locations will be',    &
     &    /,9X,'read once from ',I3,' and discharge values in units',   &
     &    /,9X,'of (m3/s) will be read from unit ',I3,' at each ',      &
     &    /,9X,F9.1,' interiorHydrologyTimeInc (seconds).',    &
     &    /,9X,'Interpolation in time is done to sync the source term', &
     &    /,9X,'with the model time step. The unit ',I3,' file begins', &
     &    /,9X,'at the time of the ',A,'.')
 3270   format(/5X,'interiorHydrologyType = ',I3                        &
     &    /,9X,'Interior hydrology forcing will be used in the GWCE',   &
     &    /,9X,'at provided lon/lat locations. Geographic locations',   &
     &    /,9X,'will be read once from ',I3,' and discharge values',    &
     &    /,9X,'in units of (m3/s) will be read from unit ',I3,' at',   &
     &    /,9X,' each ',F9.1,' interiorHydrologyTimeInc (seconds).',    &
     &    /,9X,'Interpolation in time is done to sync the source term', &
     &    /,9X,'with the model time step. The unit ',I3,' file begins', &
     &    /,9X,'at the time of the ',A,'.')
 3271   format(/,5X, 'interiorHydrologyType = ',I3,                     &
     &    /,9X,'Your selection (a UNIT 15 input parameter) is not an',  &
     &    /,9X,'allowable value. Execution will be terminated.')
 3272   format(//5X,'useInteriorHydrologyRamp = .true.',                &
     &    /,9X,'Interior hydrology forcing will gradually be ramped',   &
     &    /,9X,'to full strength over the first 2 hours of simulation.')
 3273   format(//5X,'useInteriorHydrologyRamp = .false.',               &
     &    /,9X,'Interior hydrology forcing will be at full strength',   &
     &    /,9X,'from the start of the simulation.')
 
#if defined(ALL_TRACE)
      call allMessage(DEBUG,"Return")
#endif
      call unsetMessageSource()
      return
  
!-----------------------------------------------------------------------
      end subroutine interiorHydrologyCheck
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> Called from read_input.F (when useInteriorHydrology is .true.) to 
!> allocate all arrays necessary for hydrology point source terms.
!> Also calls readInteriorHydrologyLocs() to read location files for
!> interior hydrology options, since this only needs to be done once.
!> Discharge values are processed in ___InteriorHydrologyForcing
!> called from cstart/hstart/timestep.  Only allocation of the arrays
!> and location of interior hydrology source terms are processed here.
!-----------------------------------------------------------------------
      subroutine allocInteriorHydrology()
!-----------------------------------------------------------------------

      implicit none
  
      integer :: ioerr   !IO error for interior hydrology
      integer :: psLun   !unit number for input file
     
      character(8) :: hFile  !default filename for input

      call setMessageSource("allocInteriorHydrology")
#if defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter")
#endif
  
      ! First set of allocate statements are utilized by
      ! all the options for point source terms and do not require
      ! reading hydrology input
      allocate ( hydroCurr(mnp), hydroODT(mnp) )
      hydroCurr(1:mnp) = 0.D0
      hydroODT(1:mnp) = 0.D0
  
      write(16,3010)

      ! Set file name for reading input based on Type
      select case (abs(interiorHydrologyType))
       case(1)
          psLun = 420
          hFile = 'fort.420'
       case(2)
          psLun = 428
          hFile = 'fort.428'
      end select

      ! Now arrays that require reading the hydrology input:
      ! May need to add more arrays with new types.
      call openFileForRead(psLun,trim(inputdir)//'/'//                  &
     &        trim(hFile),ioerr)
      if (ioerr.GT.0) call hydrologyTerminate()
      read(psLun,*) numHydro
      allocate ( hydro1(numHydro), hydro2(numHydro) )
      allocate ( hydroX(numHydro), hydroOld(numHydro) )
      allocate ( hydroNearNode(numHydro) )
      allocate ( isActive(numHydro) )
      hydro1(1:numHydro) = 0.D0
      hydro2(1:numHydro) = 0.D0
      hydroX(1:numHydro) = 0.D0
      hydroOld(1:numHydro) = 0.D0
      hydroNearNode(1:numHydro) = 0
      isActive(1:numHydro) = 0

      ! Subroutine call to read in actual locations
      call readInteriorHydrologyLocs(psLun)
      close(psLun)

! DEBUG: 10/03/2025 Add printing to check hydrology status
!#ifdef CMPI
!      open(888,file=trim(localdir)//'/'//'fort.888', action='write',    &
!     &     status='replace')
!      write(888,*) 'Debugging hydrology sources:'
!      write(888,'(2A)') 'time nearNode isActive divByTotalArea',        &
!     &                  ' hydroX hydroCurr'
!#else
!      open(888,file='fort.888', action='write', status='replace')
!      write(888,*) 'Debugging hydrology sources:'
!      write(888,'(2A)') 'time nearNode isActive divByTotalArea',        &
!     &                  ' hydroX hydroCurr'
!#endif

 3010 format(/," Interior Hydrology Sources")
 
#if defined(ALL_TRACE)
      call allMessage(DEBUG,"Return")
#endif
      call unsetMessageSource()

      RETURN
!-----------------------------------------------------------------------
      end subroutine allocInteriorHydrology
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> Called from allocInteriorHydrology() to read in locations where
!> interior hydrology sources will be applied. Opens and reads the
!> location file. Geographic locations require additional processing
!> using kdtsearchHydrology()
!> CMS: update documentation below
!> What is read in depends upon the value of ps_lun:   
!>  428 -  standard lon, lat input (two columns) 
!>  429 -  National Water Model (NWM) input (five columns) 
!>         N designates endpoint for this feature (point source) 
!>     and T designates endpoint for the next downstream feature
!>  429 -  NOT ACTIVATED YET - code removed
!-----------------------------------------------------------------------     
      subroutine readInteriorHydrologyLocs(ps_lun) 
!-----------------------------------------------------------------------
   
      use global, only : deg2rad, rad2deg, ifsprots
      use mesh, only : ics, slam0, sfea0, drvspcoorsrots, cylindermap

      implicit none

      integer, intent(in) :: ps_lun  !file specifier
  
      integer :: I
      real(8) :: xcoordN, ycoordN, xcoordT, ycoordT
      real(8) :: lonN, latN
      real(8) :: lonrN, latrN   !for transform

      call setMessageSource("readInteriorHydrologyLocs") 
#if defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter")
#endif

      select case (ps_lun)
        case(420) ! adcirc node provided
           write(16,3013)
           do I=1, numHydro
              read(ps_lun,*) hydroNearNode(I)
              write(16,1880) I,hydroNearNode(I)
           end do

        case(428) ! simple lon/lat format - find nearest node
           if(ICS.EQ.1) then
             write(16,3014)
           else if(ICS.EQ.2) then 
             write(16,3015)
           end if
           do I=1,numHydro
              if (ICS.EQ.1) then       
                 read(ps_lun,*) xcoordN, ycoordN
              else
                 read(ps_lun,*) lonN, latN
                 lonN = lonN*deg2rad
                 latN = latN*deg2rad
                 if ( ifsprots .EQ. 1 ) then            
                    call drvspcoorsrots(lonrN, latrN, lonN, latN)
                 else
                    latrN = latN
                    lonrN = lonN
                 end if
                 call cylindermap(xcoordN, ycoordN,                     &
     &                lonrN, latrN, slam0, sfea0, ics)
              end if 
              call kdtSearchHydrology(xcoordN,ycoordN,hydroNearNode(I))
              if (ICS.EQ.1) then
                 write(16,1881) I,hydroNearNode(I),                     &
     &                xcoordN, ycoordN
              else
                 write(16,1882) I,hydroNearNode(I),                     &
     &                lonN*rad2deg,latN*rad2deg
              end if
           end do
      end select
      write(16,3016)

3013  format(/,"TYPE   FEATURE      NODE",/)
3014  format(/,"TYPE   FEATURE      NODE         X             Y",/)
3015  format(/,"TYPE   FEATURE      NODE      LON(deg)     LAT(deg)",/)
3016  format(/,"Warning: if NODE=0, hydrology source is not in domain", &
     &        /,"and will not be activated.",/)
1880  format("Actual",1X,I9,4X,I9)
1881  format("Near",2X,I9,4X,I9,2(4X,F14.2))
1882  format("Near",2X,I9,4X,I9,2(4X,F9.4))

#if defined(ALL_TRACE)
      call allMessage(DEBUG,"Return")
#endif
      call unsetMessageSource()
  
!-----------------------------------------------------------------------
      end subroutine readInteriorHydrologyLocs
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> Called from timestep.F to perform I/O, time interpolation and any
!> other procedures necessary to prepare interior hydrology data at a
!> particular point in time. Add other options by expanding the
!> interiorHydrologyType case options. All discharge values should be
!> provided in units of (m^3/s) unless otherwise stated in each case. 
!-----------------------------------------------------------------------
      subroutine getInteriorHydrology(timeloc)
!-----------------------------------------------------------------------

      use global, only : dt
      use mesh, only : np, divbyTotalArea0
 
      implicit none
 
      real(8), intent(in) :: timeloc ! current simulation time

      integer :: i ! iterator for loop
      integer :: indXH ! node index for hydro source
      real(8) :: hydroRatio ! temporal interpolation weight
      real(8) :: rampPS ! multiplier to time-delay source terms

      call setMessageSource("getInteriorHydrology")
#if defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter")
#endif

      ! Update interior hydrology(PS) ramp
      ! CMS: Decide if we want to keep this feature
      if (useInteriorHydrologyRamp) then
         rampPS = tanh((2.d0*timeloc/86400.d0)/hydroRamp)
      else
         rampPS = 1.d0
      endif 

      ! Update active array for interior hydrology
      call checkActiveHydrology()

      ! Check if need to read more data
      if(timeloc.GT.hydroTime2) then
         hydroTime1 = hydroTime2
         hydroTime2 = hydroTime2 + interiorHydrologyTimeInc
         hydro1(1:numHydro) = hydro2(1:numHydro)
         hydro2(1:numHydro) = 0.d0 ! initialize at zero each input cycle
         ! All interiorHydrologyType options currently read from a
         ! fort.430 file; if this changes, will use:
         ! select case (ABS(interiorHydrologyType))
         do I=1,numHydro
            read(430,*) hydro2(I)
         end do
      endif

      ! v56+ UPDATE: move "dry" determination here instead of gwce.F and
      ! also build the global array using loop through numHydro instead
      ! of np - changed significantly from original implementation
      hydroRatio = (timeloc-hydroTime1)/interiorHydrologyTimeInc
      hydroCurr(1:np) = 0.d0  ! reset since accumulating
      hydroODT(1:np) = 0.d0  ! reset since accumulating
      do I=1,numHydro
         indXH=hydroNearNode(I)  ! mesh node id
         ! individual source contribution (numHydro arrays)
         hydroOld(I) = hydroX(I) ! save current value
         hydroX(I)=rampPS*(hydro1(I)+hydroRatio*(hydro2(I)-hydro1(I)))* &
     &       3.d0*divbyTotalArea0(indXH)
         ! accumulate at shared junctions (np arrays) and apply active
         hydroCurr(indXH) = hydroCurr(indXH) + hydroX(I)*isActive(I)
         hydroODT(indXH) = hydroODT(indXH) + isActive(I)*               &
     &                    (hydroX(I)-hydroOld(I))/dt
! DEBUG: add write statements to track variables - OK since no accum
! fort.888 opened from allocateInteriorHydrology() above
!         write(888,fmt=8880) timeloc, indXH, isActive(I),               &
!     &        divbyTotalArea0(indXH), hydroX(I), hydroCurr(indXH)
      end do
!8880  format(F10.2,1X,I9,1X,I2,2X,F8.6,1X,F8.6,1X,F8.6)

#if defined(ALL_TRACE)
      call allMessage(DEBUG,"Return")
#endif
      call unsetMessageSource()
  
!-----------------------------------------------------------------------
      end subroutine getInteriorHydrology
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> Called from cstart.F to setup interior hydrology I/O for a coldstart
!-----------------------------------------------------------------------
      subroutine coldstartInteriorHydrology()
!-----------------------------------------------------------------------

      use global, only : statim
      use sizes, only : inputdir
      use mesh, only : divByTotalArea0

      implicit none
      
      integer :: i ! node loop counter
      integer :: indXH ! local nearNode index
      integer :: ioerr ! IO error code
      real(8) :: rampPS ! multiplier to time-delay PS terms

      call setMessageSource("coldstartInteriorHydrology")
#if defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter")
#endif

      ! Interior hydrology ramp is independent of all other ramps and has
      ! a default value of 2 hrs (0.08333 days) in case the interior
      ! streamflows are very high at the beginning - testing indicated
      ! that it may not be needed.
      if(useInteriorHydrologyRamp) then
        rampPS = 0.d0
      else
        rampPS = 1.d0
      endif

      ! Arrays are allocated and initialized during readInput but
      ! need to update isActive array after nodecode population
      call checkActiveHydrology()

      ! For all interiorHydrologyType options, discharges are read from
      ! a fort.430 file and locations are processed during readInput in
      ! other subroutines. Only need to read the first two snaps and
      ! initialize the current hydrologic sources and time variables.
      call openFileForRead(430,trim(inputdir)//'/'//               &
     &          'fort.430',ioerr)
      if (ioerr.GT.0) call hydrologyTerminate()
      !Read first set of values
      do I=1,numHydro
         read(430,*) hydro1(I)
      end do
      !Read second set of values
      do I=1,numHydro
         read(430,*) hydro2(I)
      end do
      !Update time dependent arrays
      hydroTime1 = statim*86400.d0
      hydroTime2 = hydroTime1 + interiorHydrologyTimeInc
      do I=1,numHydro
         indXH = hydroNearNode(I)
         ! Discharges must be accumulated since multiple sources could
         ! contribute to the same ADCIRC node at junctions. Must also
         ! check for active status before applying in gwce module.
         ! hydroCurr is initialized to zero during allocation
         hydroCurr(indXH) = hydroCurr(indXH) + rampPS*hydro1(I)*&
     &       3.d0*divbyTotalArea0(indXH)*isActive(I)
      end do 
 
#if defined(ALL_TRACE)
      call allMessage(DEBUG,"Return")
#endif
      call unsetMessageSource()

      end subroutine coldstartInteriorHydrology


!------------------------------------------------------------------------
!> Setup interior hydrology i/o for a hotstart. Called from hstart.F
!>
!> + file starts from coldstart time - read until hotstart time
!>                                     and then read next snap
!> - file starts from hotstart time - read first snap in file
!------------------------------------------------------------------------
      subroutine hotstartInteriorHydrology(timeloc)
!------------------------------------------------------------------------

      use global, only : statim, dt, iths
      use sizes, only : inputdir
      use mesh, only : np, divByTotalArea0

      implicit none
  
      real(8), intent(in) :: timeloc
      integer :: i ! node loop counter
      integer :: indXH ! temp global node ID for hydro
      integer :: it ! temporal loop counter
      integer :: ioerr ! I/O error code
      real(8) :: timeit ! current time
      real(8) :: rampPS ! multiplier to time-delay PS terms

      call setMessageSource("hotstartInteriorHydrology")
#if defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter")
#endif

      ! Initialize the interior hydrology ramp and update isActive()
      ! before reading in and processing input. Note that initial call
      ! to totalAreaCalc() is done within the main hstart module.
      if (useInteriorHydrologyRamp) then
         rampPS = tanh((2.d0*timeloc/86400.d0)/hydroRamp)
      else
         rampPS = 1.d0
      endif 
      call checkActiveHydrology()

      ! Read initial two sets into memory (same for all intHydroTypes)
      call openFileForRead(430,trim(inputdir)//'/'//              &
     &     'fort.430',ioerr)
      if (ioerr.GT.0) call hydrologyTerminate()
      !Read first set of values
      do I=1,numHydro
          read(430,*) hydro1(I)
      end do
      !Read second set of values			
      do I=1,numHydro
          read(430,*) hydro2(I)
      end do

      ! For all interiorHydrologyType options, discharges are read from
      ! a fort.430 file and locations are processed during readInput in
      ! other subroutines. Only need to test if file starts at coldstart
      ! time so can skip snaps and initialize the current hydrologic
      ! sources and time variables.
      if (interiorHydrologyType .gt. 0) then  !coldstart: skip
         !Read until get to current time
         hydroTime1 = statim*86400.d0
         hydroTime2 = hydroTime1 + interiorHydrologyTimeInc
         do IT=1,iths
             timeit=it*dt + statim*86400.d0
             if (timeit.GT.hydroTime2) then 
                 hydroTime1=hydroTime2
                 hydroTime2=hydroTime2+interiorHydrologyTimeInc
                 hydro1(1:numHydro)=hydro2(1:numHydro)
                 hydro2(1:numHydro)=0.D0  ! intialize for each new read
                 do I=1,numHydro
                    read(430,*) hydro2(I)
                 end do 
             end if
         end do
      else  !hotstart: no skip, but set time
         hydroTime1 = timeloc
         hydroTime2 = hydroTime1 + interiorHydrologyTimeInc
      end if

      ! Now update time dependent arrays: need to consider that
      ! multiple discharges may be applied at the same location,
      ! must add not just assign. Also use isActive to make sure that
      ! source shoudl be contributing within the gwce module.
      ! 09/19/25: Do I need to update hydroODT in a hotstart? Have not
      ! done in previous versions - is this a potential error? Read
      ! notes from original implementation
      hydroCurr(1:np) = 0.d0  ! reset to zero before accumulating
      do I=1,numHydro
         indXH = hydroNearNode(I)
         hydroCurr(indXH) = hydroCurr(indXH) + rampPS*hydro1(I)*&
     &       3.d0*divbyTotalArea0(indXH)*isActive(I)
      end do 

#if defined(ALL_TRACE)
      call allMessage(DEBUG,"Return")
#endif
      call unsetMessageSource()

      end subroutine hotstartInteriorHydrology


!------------------------------------------------------------------------
!> Check for active hydrology locations: called at end of wetdry
!> Original implementation was to multiply each nodal hydroCurr value by
!> the nodecode to handle dry nodes. However, that zeroed out all sources
!> that terminate at the same nearNode. The use of isActive during the
!> accumulation in getInteriorHydrology() removes this and lets each
!> source be handled separately. This is only strictly necessary for
!> dynamic locations, but is more flexible for all Types. 
!> Note: totalArea and divByTotalArea are updated at wet/dry step4
!------------------------------------------------------------------------
      subroutine checkActiveHydrology()
!------------------------------------------------------------------------

      use global, only : nodecode

      implicit none

      integer :: I

       ! totalArea and 1/totalArea are updated in wetdry.F
       ! Part4 with a call to totalAreaCalc() during timestepping and
       ! initialized from main cstart and hstart modules to begin.

       ! Add other tests for other interiorHydrologyType values
       ! and call this subroutine after the nodecode update occurs
       ! in wetdry.F

       ! May need an update similar to below for all "active" arrays
       !#ifdef CMPI
       !         CALL UPDATEI(NNODECODE,IDUMY,1)
       !#endif

      select case(abs(interiorHydrologyType))
        case (1,2) ! Static locations - just use nodecode
          do I = 1, numHydro
             isActive(I) = nodecode(hydroNearNode(I))
          end do
        !case (3) ! dynamic locations, need to find most upstream
        ! active node and deactivate the downstream sources
        ! Below lines are from the very first v53 PS code - not useable
        ! as is
        ! isActive(:)=1 ! resets all values to on or active
        ! DO I=1,numPS
        !    IF (DivByTotalArea(nearNode(I)).NE.0.d0) THEN 
             ! node is active and will contribute, turn
             ! off "to" node if exists
        !       IF (toNode(I).NE.0) THEN 
        !          isActive(toNode(I))=0
        !       END IF 
        !    END IF 
        ! END DO 

      end select

!------------------------------------------------------------------------
      end subroutine checkActiveHydrology
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!                        P R I V A T E
!------------------------------------------------------------------------
!------------------------------------------------------------------------

      !cms: Added this subroutine to eliminate dependence on adcirc_mod
      ! module, which should be built last to simplify the build system.
      ! Added additional deallocate statements to clean up arrays
      subroutine hydrologyTerminate(NO_MPI_FINALIZE)

      use global, only : INFO, setMessageSource, allMessage,            &
     &    unsetMessageSource
#ifdef cmpi
      use messenger, only : subdomainFatalError, msg_fini
#endif
      implicit none

      logical, optional :: no_mpi_finalize

      call setMessageSource("hydrologyTerminate")
      call allMessage(INFO,"ADCIRC terminating")

      if (allocated(hydro1)) deallocate(hydro1)
      if (allocated(hydro2)) deallocate(hydro2)
      if (allocated(hydroX)) deallocate(hydroX)
      if (allocated(hydroOld)) deallocate(hydroOld)
      if (allocated(hydroCurr)) deallocate(hydroCurr)
      if (allocated(hydroODT)) deallocate(hydroODT)
      if (allocated(hydroNearNode)) deallocate(hydroNearNode)
      if (allocated(isActive)) deallocate(isActive)
      close (888)

#ifdef cmpi
      subdomainFatalError = .true.
      if (present(NO_MPI_FINALIZE)) then
        call msg_fini(NO_MPI_FINALIZE)
      else
        call msg_fini()
      endif
#endif
      stop 1 

      call unsetMessageSource()

      end subroutine hydrologyTerminate


!------------------------------------------------------------------------
!> Uses the KDTREE2 algorithm for finding which ADCIRC node is closest
!> for a given interior hydrology source. It was necessary to modify
!> the main kdtsearch() routine since if the hydrology source is not
!> within an element, it cannot be added - unlike a recording station
!> a source term must be located within the domain. 
!> 
!> Hydrology sources that were not found during preprocessing are 
!> indicated by lon=lat=-999.0 which should return a "not found" result
!> from kdtree. 
!> 
!> Returns the nearest node instead of the found element. 
!>
!------------------------------------------------------------------------
       subroutine kdtSearchHydrology(InputXCoordinate,InputYCoordinate, &
     &         OutputNode)
!------------------------------------------------------------------------
       use sizes, only : myproc
       use global, only : srchdp, tree, kdresults
       use mesh, only : nm, x, y, rmax, dp
       use kdtree2_module
 
       implicit none

       real(8), intent(in) :: InputXCoordinate    ! cartesian/proj
       real(8), intent(in) :: InputYCoordinate    ! cartesian/proj
       integer, intent(out) :: OutputNode

       integer :: Element
       integer :: ielm(3), itc, iek
       real(8) :: X1, X2, X3, X4, Y1, Y2, Y3, Y4
       real(8) :: Xsta, Ysta
       real(8) :: A1, A2, A3, AE, AREASK, AA
       real(8) :: SD1, SD2, SD3, SD12, SD13, SD23
       real(8) :: DP1, DP2, DP3
       real(8) :: ELMMIN(2), XELM(3), YELM(3), DPELM(3), DIST
       logical :: ElementFound ! .true. = corresponding element is found
       logical :: KeepSearch ! .true. = minimum distance not found yet

       real(8), parameter :: Tolerance = 1.0d-5 ! area diff for match
 
       call setMessageSource("kdtSearchPS")
#if defined(READ_INPUT_TRACE) || defined(ALL_TRACE)
       call allMessage(DEBUG,"Enter.")
#endif
       ElementFound = .false.
       KeepSearch = .false.

       Xsta = InputXCoordinate
       Ysta = InputYCoordinate

       call kdtree2_n_nearest(tp=tree,qv=(/Xsta,Ysta/),                 &
     &                   nn=srchdp,results=kdresults)

       ! Check to see if the points lies within rmax of any of these 
       ! elements.
       itc = 1 
       elmmin = minval(sqrt(kdresults(1:srchdp)%dis)                    &
     &       - rmax(kdresults(1:srchdp)%idx) )

       if (elmmin(1).LE.0.0d0) then 
        ! Point lies within search radius of an element - loop 
        ! through the elements in the search list
         do while ((ElementFound.EQV..false.).AND.(itc.LE.srchdp))
           iek = kdresults(itc)%idx ! Current search element number
           ! Get the distance from this point to the center of the 
           ! current element
           dist = sqrt(kdresults(itc)%dis)
           ! If the distance is less than or equal to rmax
           ! (rmax = 1.5*element radius) then the point is near 
           ! the element and might be in the element
           ! Proceed with the weights test
           if (dist-rmax(iek).LE.0.0d0) then 
             ! Get the shape function for this element
             ielm(:) = NM(iek,(/1,2,3/)) ! element's node number
             xelm(:) = X(ielm(:))  ! element's vertex x-values
             yelm(:) = Y(ielm(:))  ! element's vertex y-values
             X1 = xelm(1)
             X2 = xelm(2)
             X3 = xelm(3)
             Y1 = yelm(1)
             Y2 = yelm(2)
             Y3 = yelm(3)
             A1 = (Xsta-X3)*(Y2-Y3) + (X2-X3)*(Y3-YSta)
             A2 = (Xsta-X1)*(Y3-Y1) - (YSta-Y1)*(X3-X1)
             A3 = (Ysta-Y1)*(X2-X1) - (XSta-X1)*(Y2-Y1)
             AA = ABS(A1) + ABS(A2) + ABS(A3)
             AREASK=X2*Y3+X1*Y2+X3*Y1-Y1*X2-Y2*X3-Y3*X1 
             AE = ABS(AA-AREASK)/AREASK
             if (AE.LT.Tolerance) then 
               ! Do not calculate distances and find nearest 
               ! node until you've found the point within the 
               ! domain. 
               ElementFound = .true.
               dpelm(:) = dp(ielm(:))
               dp1 = dpelm(1)
               dp2 = dpelm(2)
               dp3 = dpelm(3)
               SD1 = (X1-Xsta)**2 + (Y1-Ysta)**2
               SD2 = (X2-Xsta)**2 + (Y2-Ysta)**2
               SD3 = (X3-Xsta)**2 + (Y3-Ysta)**2
               ! First check for any zero lengths - point is ADCIRC node
               if (SD1 .EQ. 0.d0) then
                  OutputNode = ielm(1) 
               else if (SD2 .EQ. 0.d0) then
                  OutputNode = ielm(2) 
               else if (SD3 .EQ. 0.d0) then
                  OutputNode = ielm(3) 
               else
                  KeepSearch = .true.
                  SD12 = ABS(SD1-SD2)/SD1
                  SD13 = ABS(SD1-SD3)/SD1
                  SD23 = ABS(SD2-SD3)/SD2
               end if
               if (KeepSearch) then
                  ! Check each tolerance to see if any lengths are equal
                  if (SD12.LT.Tolerance) then 
                     if (SD13.LT.Tolerance) then 
                     ! PS is equidistant from all nodes, check
                     ! bathymetry - Bathymetry is positive and 
                     ! topography is negative 
                       if ((dp1.GE.dp2).AND.(dp1.GE.dp3)) then 
                         OutputNode = ielm(1) 
                       else if ((dp2.GE.dp1).AND.(dp2.GE.dp3)) then 
                         OutputNode = ielm(2)
                       else
                         OutputNode = ielm(3)
                       end if 
                     else
                     ! PS is equidistant from nodes 1 and 2
                       if (dp1.GE.dp2) then 
                         OutputNode = ielm(1)
                       else
                         OutputNode = ielm(2)
                       end if 
                     end if 
                  else if (SD13.LT.Tolerance) then 
                  ! PS is equidistant from nodes 1 and 3
                     if (dp1.GE.dp3) then 
                       OutputNode = ielm(1)
                     else
                       OutputNode = ielm(3)
                     end if 
                  else if (SD23.LT.Tolerance) then 
                  ! PS is equidistant from nodes 2 and 3
                     if (dp2.GE.dp3) then 
                       OutputNode = ielm(2) 
                     else 
                       OutputNode = ielm(3) 
                     end if 
                  ! There are no equal distances, find the shortest
                  else if ((SD1.LT.SD2).AND.(SD1.LT.SD3)) then 
                     OutputNode = ielm(1) 
                  else if ((SD2.LT.SD1).AND.(SD2.LT.SD3)) then 
                     OutputNode = ielm(2) 
                  else 
                     OutputNode = ielm(3)
                  end if ! End of search for the nearest node
               end if
             else ! not in this element - keep looking     
               itc = itc + 1
             end if ! End area ratio test
           else 
             ! Point is too far away from the barycenter of the 
             ! element to possibly be in the element, so move to 
             ! the next element. 
             itc = itc + 1
           end if !end rmax test
         end DO ! end the while loop
       end if ! end within search radius test

       if (.not. ElementFound) then 
          OutputNode = 0
          if (XSta.GT.0 .AND. YSta.GT.0) then
             ! Negative coordinates indicate that we intentionally
             ! deactivated this feature during preprocessing.
             write(16,9893)
          end if
       end if 

 9893  format(/,1X,'!!!!!  WARNING - NONFATAL INPUT ERROR  !!!!!',/,    &
     &      1X,'Point Source does NOT lie within any element',          &
     &      1X,'and will be deactivated.',/,                            &
     &      1X,'Verify input coordinates for this source.')

#if defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif 
      call unsetMessageSource()
!-----------------------------------------------------------------------
      end subroutine kdtSearchHydrology
!-----------------------------------------------------------------------

      end MODULE hydrology
