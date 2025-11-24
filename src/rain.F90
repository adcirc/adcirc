!-----------------------------------------------------------------------
!  MODULE RAIN
!-----------------------------------------------------------------------
!> @file rain.F90
!>
!> @author Chris Szpilka <cmszpilka@ou.edu>
!> @author Kendra Dresback
!> 
!> @copyright Dr. R.A. Luettich and Dr. J.J. Westerink
!>
!> @brief This module handles the I/O and computations related to
!> adding precipitation source terms to ADCIRC. 
!>
!>  Originally added in v51 but not merged into the main trunk until v56.
!>  Precipitation source terms are only added to "wet" nodes since
!>  ADCIRC cannot route rainfall over dry land. Two accumulation levels
!>  are tracked globally: total rainfall that has fallen (totalRain) and
!>  total rainfall that has been applied to the domain (totalRainApplied).
!>  Rainfall is only applied if the element is active at the current time
!>  step and is accumulated by multiplying by the timestep. Active
!>  status can be chosen by the user using activeRainType in the
!>  namelist.
!>
!>  Precipitation is expected in units of rainfall rate [L/T], where
!>  units must be consistent with the input in the fort.15 control
!>  file. In most cases, rainfall rate should be provided in [mm/hr]
!>  unless otherwise noted in the rainType description.
!>----------------------------------------------------------------------
!>  USAGE:
!>  rainfallControl namelist at bottom of the fort.15 activates this
!>  feature starting at v56.
!>
!>  useRain (logical): default .false.
!>  rainType (integer): default 0 (see options below)
!>  rainTimeInc (float): default 3600.d0 (hourly)
!>     time increment between input - controls reading and interpolation
!>  activeRainType (integer): default 1
!>     1) currently wet nodes are considered active (nodecode)
!>     2) ND implementation to only apply over water - nodes are only
!>        active if they are currently wet and bathymetrically
!>        defined as water (DP > 0)
!>        NOTE: does not allow rainfall to be applied on any nodes that
!>        become wet during the simulation or rivers that are given
!>        initial water elevations but are above MSL.
!>
!>  Additional precipitation source types can be included by adding more
!>  options to the rainType case constructs in this module. Currently,
!>  two types are available:
!>
!>  rainType = 1:
!>  global rainfall over entire domain at rainTimeInc intervals
!>     This option should only be used for small serial test cases
!>
!>     Units for rainfall rate: [m/s]
!>     Requires one file as input:
!>     fort.270 (must contain enough snaps to cover each rainTimeInc)
!>
!>  rainType = 12:
!>  OWI format file structure
!>     Values are provided in global/regional scale domains
!>     similar to OWI winds with rainfall rates given in [mm/hr]
!>     Conversion to units of [m/s] occurs within owi_rain using a
!>       predetermined multiplier. Use of an optional user defined
!>       multiplier for other input units could be provided in the
!>       fort.27 file but is not supported yet in the code.
!>
!>     Requires two to three files as input:
!>     fort.27 main rainfall control file:
!>            follows format of fort.22 for OWI winds
!>       numSets        <-- Number of sets of rain fields (1 or 2)
!>                          1: global file only (fort.271)
!>                          2: global and regional (fort.271 & fort.272)
!>       numBlankSnaps  <-- Number of snaps to add/skip at beginning
!>       rainMultiplier <-- Rain intensity multiplier (if need to
!>                          convert units)
!>                          NOT used yet 
!>            Default assumes that files are in [mm/hr] which are
!>            converted to [m/s] in owi_rain.F
!>
!>     fort.271 global scale domain rainfall
!>     fort.272 regional scale domain rainfall
!>
!>  I/O unit numbers are set as described above until/unless these are
!>  automated in the larger ADCIRC scheme.
!>
!> @todo after enough rainfall has accumulated (Hmin?), turn on
!> the dry nodes and apply the accumulated rainfall. Would only want
!> to do this on flat elements - not steep regions - since there is no
!> internal routing. Precipitation must be added as a rate, so will need
!> to convert the accumulation somehow. Also, add code to support the
!> use of an optional multiplier for every rainType.
!-----------------------------------------------------------------------
      MODULE RAIN
!-----------------------------------------------------------------------

        use global, only : ScreenUnit, setMessageSource, allMessage,&
     &    screenMessage, logMessage, unsetMessageSource,screenunit,& 
     &    ECHO, INFO, DEBUG, openFileForRead, nodecode
        use sizes, only : localdir, inputdir, mnp

        implicit none

        ! Local precipitation timing variables
        real(8) :: rainTime1 ! time associated with previous rain data
        real(8) :: rainTime2 ! time associated with new rain data

        ! Precipitation arrays for input datasets
        real(8), allocatable :: rain1(:) ! previous rainfall dataset
        real(8), allocatable :: rain2(:) ! new rainfall dataset

        ! Precipitation array for tracking active status
        real(8), allocatable :: activeRain(:)

        ! namelist input from read_input.f
        logical :: useRain=.false.  ! default - do not use
        integer :: rainType=0  ! default none
        real(8) :: rainTimeInc=3600.d0 ! rain input interval (sec)
        integer :: activeRainType=1  ! default current wet/dry state
        namelist /rainfallControl/ useRain, rainType, rainTimeInc,      &
     &            activeRainType
  
        ! Source terms used in gwce.F
        real(8),allocatable :: rainCurr(:) ! rainfall at current ts
        real(8),allocatable :: rainODT(:) ! rainfall time derivative
        ! Rainfall used in write_output.F
        real(8),target,allocatable :: totalRain(:)  !accumulated rain
        real(8),target,allocatable :: totalRainApplied(:) !applied rain
        real(8),target,allocatable :: rain00(:) ! station output array
        
        ! Public/private variables
        private :: rainTime1, rainTime2, rain1, rain2
        public :: useRain, rainType, rainTimeInc, activeRainType
        public :: rainfallControl, activeRain
        public :: rainCurr, rainODT, totalRain, totalRainApplied
        public :: rain00

        ! Public/private subroutines
        private :: checkActiveRain, rainTerminate
        public :: rainTypeCheck, allocRain, getRainForcing
        public :: coldstartRainForcing, hotstartRainForcing
!---------------------end of data declarations--------------------------

      CONTAINS

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!                         P U B L I C
!-----------------------------------------------------------------------
!----------------------------------------------------------------------


!----------------------------------------------------------------------
!  Subroutine rainTypeCheck
!-----------------------------------------------------------------------
!> @brief Checks validity of user entered value for rainType
!>
!> Called from read_input (when useRain is .true.) to check
!> validity of precipitation source type and log messages.
!> Moved all message logging for rainfall from read_input.
!>
!> The value of rainType is set during readinput when the 
!> rainfallControl namelist if present, default value is 0
!-----------------------------------------------------------------------
      subroutine rainTypeCheck()
!-----------------------------------------------------------------------

        implicit none
  
        character(10) :: reftime

        call setMessageSource("rainTypeCheck")
#if defined(ALL_TRACE)
        call allMessage(DEBUG,"Enter.")
#endif

        ! Check value of useRain before messages are logged
        ! otherwise errors may arise when not being used
        ! but namelist was found in read_input.
        if (.not. useRain) rainType=0
        ! Echo input messages related to implementation
        select case (rainType)
        case(0)
           write(16,3298) rainType
           activeRainType=0  ! to remove message below
        case(1) ! global rainfall (coldstart time)
           reftime='coldstart'
           write(16,3299) rainType, trim(reftime)
        case(-1) ! global rainfall (hotstart time)
           reftime='hotstart'
           write(16,3299) rainType, trim(reftime)
        case(12) ! OWI format (coldstart time)
           reftime='coldstart'
           write(16,3302) rainType, trim(reftime)
        case(-12) ! OWI format (hotstart time)
           reftime='hotstart'
           write(16,3302) rainType, trim(reftime)
        case default
           write(screenunit,9816) rainType
           write(16,9816) rainType
           call rainTerminate()
        end select

        ! Check value of activeRainType and log messages
        select case (activeRainType)
        case (0)
           !no message just reset to default value
           activeRainType = 1
        case (1)
           write(16,3303) activeRainType
        case (2)
           write(16,3304) activeRainType
        case default
           write(16,9817) activeRainType
           activeRainType = 1
        end select


 3298   format(/,5X, 'rainType = ',I3,&
     &   /,9X,'Precipitation forcing will not be used in',&
     &   ' the GWC equation') 

 3299   format(/5X,'rainType = ',I3,&
     &    /,9X,'Precipitation forcing will be used in the GWCE.',&
     &    /,9X,'Rainfall values are read from a global input file',&
     &    /,9X,'in units of (m/s). Interpolation in time is done',&
     &    /,9X,'to sync the rainfall datasets with the model time',&
     &    /,9X,'step. The input file begins at the time of the',&
     &    /,9X,A,'.')
 
 3302   format(/5X,'rainType = ',I3,&
     &    /,9X,'Precipitation forcing will be used in the GWCE.',&
     &    /,9X,'Rainfall values are read from global and/or regional',&
     &    /,9X,'OWI format files in units of (mm/hr). Interpolation',&
     &    /,9X,'in time is done to sync the rainfall datasets with',&
     &    /,9X,'the model time step. The global and/or regional',&
     &    /,9X,'files begin at the time of the ',A,'.')

 9816   format(/,5X,'rainType = ',I3,&
     &   /,9X,'Your selection (a UNIT 15 input parameter) is not an',&
     &   /,9X,'allowable value. Execution will be terminated.')

 3303   format(/,5X, 'activeRainType = ',I3,&
     &   /,9X,'Precipitation forcing will be added to all wet nodes',&
     &   /,9X,'at each time step.')

 3304   format(/,5X, 'activeRainType = ',I3,&
     &   /,9X,'Precipitation forcing will be added to wet nodes that',&
     &   /,9X,'have positive bathymetry (water) at each time step.')

 9817   format(/,5X,'activeRainType = ',I,&
     &   /,9X,'Your selection (a UNIT 15 input parameter) is not an',&
     &   /,9X,'allowable value. Reverting to the default value [1]: ',&
     &   /,9X,'precipitation will be added to all wet nodes.')
 
#if defined(ALL_TRACE)
        call allMessage(DEBUG,"Return.")
#endif
        call unsetMessageSource()
        return
  
!-----------------------------------------------------------------------
      end subroutine rainTypeCheck
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> Called from read_input.F (when useRain is .true.) to allocate
!> all arrays necessary for precipitation source terms.
!-----------------------------------------------------------------------
      subroutine allocRain()
!-----------------------------------------------------------------------
        use sizes, only : mnp, mnstae

        implicit none
  
        call setMessageSource("allocRain")
#if defined(ALL_TRACE)
        call allMessage(DEBUG,"Enter")
#endif

        allocate ( rain1(mnp), rain2(mnp) )
        allocate ( rainODT(mnp), rainCurr(mnp) )
        allocate ( totalRain(mnp) )
        allocate ( totalRainApplied(mnp) )
        allocate ( activeRain(mnp) )
        rain1(1:mnp) = 0.d0
        rain2(1:mnp) = 0.d0
        rainCurr(1:mnp) = 0.d0
        rainODT(1:mnp) = 0.d0
        totalRain(1:mnp) = 0.d0
        totalRainApplied(1:mnp) = 0.d0
        activeRain(1:mnp) = 0.d0

        ! check if elevation station output requested - rain station
        ! output tied to elevation station output
        if (mnstae .ge. 1) then
           allocate ( rain00(mnstae) )
           rain00(1:mnstae) = 0.d0
        end if

#if defined(ALL_TRACE)
        call allMessage(DEBUG,"Return")
#endif
        call unsetMessageSource()
        return
!-----------------------------------------------------------------------
      end subroutine allocRain
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> @brief
!>  Subroutine to update rain arrays used in gwce
!>
!>  Called from timestep.F (when useRain is .true.) to perform i/o,
!>  time interpolation, and any other procedures necessary to prepare
!>  precipitation data at a particular point in time. Add other options
!>  by expanding the rainType case options. All values returned in (m/s)
!>  and expected units included in comments for each case. Files opened
!>  from cstart or hstart module.
!-----------------------------------------------------------------------
      subroutine getRainForcing(timeloc)
!-----------------------------------------------------------------------

        use global, only : dt, nodecode
        use mesh, only : np
        use owi_rain, only : rain12get

        implicit none

        real(8), intent(in) :: timeloc

        ! Temporary scalar variables for computations within loops
        integer :: nhg ! node number temp
        integer :: i ! node loop counter
        real(8) :: rainRatio ! linear interpolation for current time
        real(8) :: rainOld !rainfall at previous ts for rainODT comp
        real(8) :: rainX ! current value of rainfall after interpolation
        real(8) :: rainAcc ! rain accumulated during timestep

        call setMessageSource("getRainForcing")
#if defined(ALL_TRACE)
        call allMessage(DEBUG,"Enter")
#endif
 
        ! Check the active status of all nodes
        call checkActiveRain()

        ! Check if need to read more data
        if(timeloc.GT.rainTime2) then
            rainTime1 = rainTime2
            rainTime2 = rainTime2 + rainTimeInc
            rain1(1:np) = rain2(1:np)
            select case (abs(rainType))
                case(1)   !global domain coverage for ideal test cases
                    read(270,*) (nhg,rain2(i),i=1,np)

                case(12)  !OWI format precip forcing
                    call rain12get(rain2,np)
            end select
        endif !end check if need more rain

        ! Compute terms used in GWCE and accumulate rainfall
        rainRatio = (timeloc-rainTime1)/rainTimeInc
        do i=1,np
            rainX = rain1(i) + rainRatio*(rain2(i)-rain1(i))
            rainOld = rainCurr(i)
            rainCurr(i) = rainX
            rainODT(i) = (rainX-rainOld)/dt
            rainAcc = rainX*dt
            totalRain(i) = totalrain(i) + rainAcc
            totalRainApplied(i)=totalRainApplied(i)+                    &
                rainAcc*activeRain(i)
        end do

#if defined(ALL_TRACE)
      call allMessage(DEBUG,"Return")
#endif
      call unsetMessageSource()
      return

!-----------------------------------------------------------------------
      end subroutine getRainForcing
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> Called from cstart.F to setup precipitation I/O for a coldstart
!-----------------------------------------------------------------------
      subroutine coldstartRainForcing()
!-----------------------------------------------------------------------

        use global, only : statim, dt, nodecode
        use sizes, only : inputdir
        use mesh, only : np
        use owi_rain, only : rain12init, rain12get

        implicit none
  
        integer :: nhg    ! temporary node number
        integer :: i      ! node loop counter
        integer :: ioerr  ! I/O error flag

        real(8) :: rainAcc ! rain accumulated during timestep

        call setMessageSource("coldStartRainForcing")
#if defined(ALL_TRACE)
        call allMessage(DEBUG,"Enter")
#endif

        ! Check active status of grid points
        call checkActiveRain()

        ! Read first two fields of values
        select case (ABS(rainType))
            case(1)   !global domain coverage for ideal test cases
               call openFileForRead(270,trim(inputdir)//'/'//&
     &              'fort.270',ioerr)
               if (ioerr.GT.0) call rainTerminate()
               read(270,*) (nhg,rain1(I),i=1,np)
               read(270,*) (nhg,rain2(I),i=1,np)
               rainTime1 = statim*86400.D0
               rainTime2 = rainTime1 + rainTimeInc
            case(12)  !OWI format precip
               call rain12init(rain1,np)
               call rain12get(rain1,np)
               call rain12get(rain2,np)
               rainTime1 = statim*86400.D0
               rainTime2 = rainTime1 + rainTimeInc
        end select

        ! Assign temporal values
        !   active status applied in gwce to reduce memory allocation
        !   rainODT assumed to be zero at coldstart
        do i=1,np
           rainCurr(i) = rain1(i)
           rainAcc = rain1(i)*dt
           totalRain(i) = totalRain(i) + rainAcc
           totalRainApplied(i)=totalRainApplied(i)+rainAcc*activeRain(i)
        end do

#if defined(ALL_TRACE)
        call allMessage(DEBUG,"Return")
#endif
        call unsetMessageSource()
        return

!-----------------------------------------------------------------------
      end subroutine coldStartRainForcing
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> Called from hstart.F to setup precipitation I/O for a hotstart
!>
!>     + file starts from coldstart time - read until hotstart time
!>                                         and then read next snap
!>     - file starts from hotstart time - read first snap in file
!>
!> Precipitation expected in other modules with units (m/s) so convert
!> accordingly during I/O as more options are added to rainType.
!-----------------------------------------------------------------------
      subroutine hotstartRainForcing(timeloc)
!-----------------------------------------------------------------------

      use global, only : dt, statim, iths, nodecode
      use sizes, only : inputdir
      use mesh, only : np
      use owi_rain, only : rain12init, rain12get

      implicit none

      real(8), intent(in) :: timeloc

      integer :: ioerr  ! I/O error code
      integer :: nhg    ! temporary node number
      integer :: i      ! node loop counter
      integer :: it     ! hstart snap counter
      real(8) :: timeit ! hstart time check
      real(8) :: rainAcc ! rain accumulated during timestep

      call setMessageSource("hotStartRainForcing")
#if defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter")
#endif

      ! Check active status of grid points
      call checkActiveRain()

      select case (rainType)
         case(1)    !global domain coverage (filestart=cold)
            ! use for ideal test cases
            rainTime1 = statim*86400.d0
            rainTime2 = rainTime1 + rainTimeInc
            call openFileForRead(270,trim(inputdir)//'/'//&
     &           'fort.270',ioerr)
            if (ioerr.GT.0) call rainTerminate()
            ! Read first two precipitation fields
            read(270,*) (nhg,rain1(i),i=1,np)
            read(270,*) (nhg,rain2(i),i=1,np)
            ! Read until get to current time
            do it=1,iths
               timeit=it*dt + statim*86400.d0
               if (timeit.GT.rainTime2) then 
                  rainTime1 = rainTime2
                  rainTime2 = rainTime2 + rainTimeInc
                  rain1(1:np) = rain2(1:np)
                  read(270,*) (nhg,rain2(i),i=1,np) 
               end if 
            end do 

         case(-1)   !global domain coverage (filestart=hot)
            ! use for ideal test cases
            rainTime1 = timeloc
            rainTime2 = rainTime1 + rainTimeInc
            call openFileForRead(270,trim(inputdir)//'/'//&
     &           'fort.270',ioerr)
            if (ioerr.GT.0) call rainTerminate()
            ! Read first two precipitation fields
            read(270,*) (nhg,rain1(i),i=1,np)
            read(270,*) (nhg,rain2(i),i=1,np)

         case(12)   !OWI format precipitation (filestart=cold)
            rainTime1 = statim*86400.d0
            rainTime2 = rainTime1 + rainTimeInc
            call rain12init(rain1,np)
            ! Read first two precipitation fields
            call rain12get(rain1,np)
            call rain12get(rain2,np)
            ! Read until get to current time
            do it=1,iths
               timeit=it*dt + statim*86400.d0
               if (timeit.GT.rainTime2) then 
                  rainTime1 = rainTime2
                  rainTime2 = rainTime2 + rainTimeInc
                  rain1(1:np) = rain2(1:np)
                  call rain12get(rain2,np)
               end if 
            end do 

         case(-12)  !OWI format precipitation (filestart=hot)
            rainTime1 = timeloc
            rainTime2 = rainTime1 + rainTimeInc
            call rain12init(rain1,np)
            call rain12get(rain1,np)
            call rain12get(rain2,np)
      end select

      ! Assign temporal values
      !   active status applied in gwce to reduce memory allocation
      !   rainODT assumed to be zero at hotstart - old values not
      !   saved in the hotstart files
      do i=1,np
         rainCurr(i) = rain1(i)
         rainAcc = rain1(i)*dt
         totalRain(i) = totalRain(i) + rainAcc
         totalRainApplied(i)=totalRainApplied(i)+rainAcc*activeRain(i)
      end do
  
#if defined(ALL_TRACE)
      call allMessage(DEBUG,"Return")
#endif
      call unsetMessageSource()
      return

!----------------------------------------------------------------------
      end subroutine hotstartRainForcing
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!                             P R I V A T E
!----------------------------------------------------------------------
!----------------------------------------------------------------------

!------------------------------------------------------------------------
!> Check for active rainfall locations: called near beginning of
!> cold/hotstartRainForcing() and getRain() 
!> Unlike interior hydrology, the active status is applied directly in
!> the gwce module to reduce the number of mnp arrays that are
!> allocated. NOTE: rainfall that falls on grid in adcirc is only
!> applied on active nodes to prevent instabilities, since ADCIRC
!> cannot route water over dry land.
!>
!> Currently two options:
!> activeRainType=1 (default)
!>    Original implementation: multiply each nodal rainCurr and rainODT 
!>    value by the nodecode to handle dry nodes.
!> activeRainType=2  rain on "water" defined by bathymetry
!>    Additional implementation from ND group to add test for
!>    bathymetric value and only add rainfall for nodes that are 
!>    active and defined as water: BP > 0.d0
!>
!> @todo: monitor accumulation vs H0 and add rain when there is enough
!> probably only on "flat" portions of the domain.
!------------------------------------------------------------------------
      subroutine checkActiveRain()
!------------------------------------------------------------------------

      use global, only : nodecode
      use mesh, only : np, dp

      implicit none

      integer :: i

        call setMessageSource("checkActiveRain")
#if defined(ALL_TRACE)
        call allMessage(DEBUG,"Enter.")
#endif


      select case(activeRainType)
        case (1) ! Currently wet - just use nodecode
           activeRain(1:np) = 1.d0*nodecode(1:np)
        case (2) ! Currently wet and "water only" - ND implementation
                 ! merge(trueV, falseV, test)
                 ! NOTE: will never rain on any nodes that are above
                 ! MSL even if they become wet during the simulation
          do i=1,np
             activeRain(i) = nodecode(i)*                               &
     &            merge( 1.0d0, 0.d0, dp(i) > 0.d0 ) ;
          end do
      end select  

#if defined(ALL_TRACE)
      call allMessage(DEBUG,"Return")
#endif
      call unsetMessageSource()
      return

!------------------------------------------------------------------------
      end subroutine checkActiveRain
!------------------------------------------------------------------------


!> cms: Added this subroutine to eliminate dependence on 
!  adcirc_mod, which should be built last.
!------------------------------------------------------------------------
      subroutine rainTerminate(NO_MPI_FINALIZE)
!------------------------------------------------------------------------

#ifdef cmpi
      use messenger, only : subdomainFatalError, msg_fini
#endif
      implicit none

      logical, optional :: no_mpi_finalize

      call setMessageSource("rainTerminate")
      call allMessage(INFO,"ADCIRC Terminating")

      if(allocated(rain1)) deallocate(rain1)
      if(allocated(rain2)) deallocate(rain2)
      if(allocated(rainCurr)) deallocate(rainCurr)
      if(allocated(rainODT)) deallocate(rainODT)
      if(allocated(totalRain)) deallocate(totalRain)
      if(allocated(totalRainApplied)) deallocate(totalRainApplied)
      if(allocated(rain00)) deallocate(rain00)
      if(allocated(activeRain)) deallocate(activeRain)
      

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
      return

      end subroutine rainTerminate


      end module rain
