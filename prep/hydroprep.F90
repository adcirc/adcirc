!-----------------------------------------------------------------------
!  MODULE HYDROPREP
!-----------------------------------------------------------------------
!> @author Chris Szpilka, University of Oklahoma, cmszpilka@ou.edu
!> 
!> @copyright Dr. R.A. Luettich and Dr. J.J. Westerink
!>
!> @brief This module handles the I/O for adcprep related to adding
!> interior hydrology (point) source terms to ADCIRC. These terms
!> allow for interior contributions from hydrologic sources that
!> would not be captured at the boundary. See src/hydrology.F90 for
!> more details about these source terms in general, as well as
!> input requirements and expected files.
!>
!> @dependencies
!> Uses global, memory_usage, mesh, adc_constants, pre_global,
!>      KDTreeTools
!>-----------------------------------------------------------------------
!>  UPDATES:
!>  Fall 2025: changed from original implementation for type=1 since you
!>  cannot interpolate in time if numHydro keeps changing. Also changed
!>  unit number to beginning of hydro range [420-430] for just type=1;
!>  any new options will fill in gap from 421-427.
!>
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
!>     Requires two files as input:
!>        fort.420 static location file read once: adcirc node
!>        fort.430 discharge values (m3/s) at each point in location file
!>            for each interiorHydrologyTimeInc - discharges must be
!>            given in same order as the locations in fort.420
!>
!>  interiorHydrologyType = 2:
!>     Requires two files as input:
!>        fort.428 static location file read once: lon, lat
!>        fort.430 discharge values (m3/s) at each point in location file
!>            for each interiorHydrologyTimeInc - discharges must be
!>            given in same order as the locations in fort.428
!>
!>  interiorHydrologyType = 3: COMING SOON
!>     Requires two files as input:
!>        fort.429 dynamic location file read once: lon, lat, network
!>        fort.430 discharge values (m3/s) at each point in location file
!>            for each interiorHydrologyTimeInc - discharges must be
!>            given in same order as the locations in fort.429
!>
!>  I/O unit numbers are set as described above until these are automated
!>  in the larger ADCIRC scheme.
!>
!>----------------------------------------------------------------------- 
      MODULE hydroPrep

      use global, only : screenUnit, setMessageSource,                  &
     &    unsetMessageSource, openFileForRead, scratchMessage,          &
     &    logMessage, allMessage, screenMessage, DEBUG, ECHO, INFO,     &
     &    WARNING, ERROR
      use memory_usage

      implicit none

      ! Module variables used in multiple subroutines ("global")
      integer :: numHydro ! # of hydro sources
      integer, allocatable :: numHydroProc(:)  ! # of hydro on PROC
      integer, allocatable :: imap_hyd_LG(:,:) ! list of hydro on PROC

      ! Namelist default values
      logical :: found_hydroControl_nml = .false.
      logical :: useInteriorHydrology = .false. ! default none
      integer :: interiorHydrologyType = 0 ! default none
      real(8) :: interiorHydrologyTimeInc = 3600 ! input time (sec)
      logical :: useInteriorHydrologyRamp = .false. ! default none
      namelist /interiorHydrologyControl/ useInteriorHydrology,         &
     &          interiorHydrologyType, InteriorHydrologyTimeInc,        &
     &          useInteriorHydrologyRamp

      ! --              Public variables
      public :: useInteriorHydrology, interiorHydrologyType
      public :: interiorHydrologyTimeInc, useInteriorHydrologyRamp
      public :: found_hydroControl_nml, interiorHydrologyControl
      ! --              Public subroutines
      public :: prepInteriorHydrology

      ! --              Private variables
      private :: numHydro, numHydroProc, imap_hyd_LG
      ! --              Private subroutines
      private :: prep420, prepHydroLoc, OpenPrepHydro

!>---------------------end of data declarations--------------------------


      CONTAINS


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                             P U B L I C  
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
!> Called from adcprep.F (when useInteriorHydrology is .true.) to 
!> do basic error checking and then create local (PE) interior hydrology
!> input files using the domain decomposition of the ADCIRC grid
!> created by the routine DECOMP.
!> Calls other subroutines for the processing.
!>
!> CMS: edit when change for dynamic unit numbers and filenames
!-----------------------------------------------------------------------
      subroutine prepInteriorHydrology()
!-----------------------------------------------------------------------

      use pre_global, only : nproc

      implicit none

      integer :: I, K, iproc, indX  ! indices
      integer :: sdu(nproc) ! subdomain unit numbers
      real(8), allocatable :: QVal(:)
      logical :: success
  
      call setMessageSource("prepInteriorHydrology")
#if defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter")
#endif
  
      !-------------------------------------------------
      ! Read in and localize the intHydro location files
      !     varies by interiorHydrologyType
      !-------------------------------------------------
      select case (abs(interiorHydrologyType))
      case(1)
         call prep420()
      case(2)
         call prepHydroLoc()
         !update routines so that can call this way
         !call prepHydroLoc(428) and use for multiple types
      case default
         write(screenUnit,3271) interiorHydrologyType
         write(16,3271) interiorHydrologyType
         stop   !terminate 
      end select
 
      !----------------------------------------------------------
      ! Open and localize the global interiorHydrology discharges
      !   valid for all interiorHydrologyTypes (fort.430)
      !----------------------------------------------------------
      call OpenPrepHydro(430, 'point source locations  ',               &
     &     1, nproc, sdu, success)
      if (.not.success) then
         write(*,*) 'ERROR: Unit 430 file not preprocessed.'
         write(*,*) 'Terminating adcprep.'
         stop  !cms: stop since run will fail w/o files
      endif

      ! Allocate space and read the discharge values - not
      ! tracking memory usage in this subroutine
      allocate (QVal(numHydro))
      do  ! Read to end of file and write each snap
        ! Read in values for each snap
        do I=1,numHydro
          read(430,*,end=9999) QVal(I)
        end do 
       ! Write to correct processors
        do iproc=1,nproc
          do K=1,numHydroProc(iproc)
             indX = abs(imap_hyd_LG(K,iproc))
             write(sdu(iproc),*) QVal(indX)
          enddo
        end do
      end do 
      !Close all local and global discharge files
 9999 close(430)
      do iproc=1, nproc
         close (sdu(iproc))
      enddo

      if(allocated(QVal)) deallocate (QVal)

 3271   format(/,2X, 'useInteriorHydrology = T',  & 
     &    /,2X, 'interiorHydrologyType = ',I3,   &
     &    /,2X,'Your selection (a unit 15 parameter) is not valid.',  &
     &    /,2X,'No interior hydrology input was processed.',  &
     &    /,2X,'ERROR: adcprep terminated.')

#if defined(ALL_TRACE)
      call allMessage(DEBUG,"Return")
#endif
      call unsetMessageSource()
      return

!-----------------------------------------------------------------------
      end subroutine prepInteriorHydrology
!-----------------------------------------------------------------------



!------------------------------------------------------------------------
!------------------------------------------------------------------------
!                        P R I V A T E
!------------------------------------------------------------------------
!------------------------------------------------------------------------


!-----------------------------------------------------------------------
!> Called from prepInteriorHydrology inside adcprep.F
!> Creates local (PE) fort.420 interior hydrology input files using
!> the domain decomposition of the ADCIRC grid created by the
!> routine DECOMP.
!>
!> Author: originally written in v53 by kmd at OU, modularized by cms
!>
!> CMS: edit when change for dynamic unit numbers and filenames
!>      092325: removed dynamic numHydro and changed to fixed with
!>      two input files (fort.420/430)
!-----------------------------------------------------------------------
      subroutine prep420()
!-----------------------------------------------------------------------

      use pre_global, only : nproc, itotproc, imap_nod_gl2
      use memory_usage

      implicit none

      integer(8) :: nbytes = 0
      integer :: I,J,Iproc,Iproc2
      integer :: sdu(nproc) ! subdomain unit numbers
      integer :: indX   ! temp global node ID
      integer :: local  ! temp local node ID
      integer :: maxNumHydro = -99  ! maximum number of sources
      character(6) :: PE
      logical :: success    ! .true. if files opened without errors

      integer, allocatable :: locID(:,:)  ! local node IDs


      call setMessageSource("prep420")
#if defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter")
#endif

      call OpenPrepHydro(420, 'point source information  ',             &
     &     1, nproc, sdu, success)
      if (.not.Success) then
         !cms: change to stop since run will fail w/o files
         write(*,*) 'ERROR: Unit 420 file not preprocessed.'
         write(*,*) 'Terminating adcprep.'
         !return
         stop
      endif

      ! Read numHydro for allocation
      read(420,*) numHydro
      ! Allocate arrays
      allocate ( locID(numHydro,nproc) )
      allocate ( imap_hyd_LG(numHydro,nproc) ) ! module variable
      allocate ( numHydroProc(nproc) )  ! module variable
      nbytes = (4+4*2*nproc)*NumHydro ! 3 integer: 2 dims
      call memory_alloc(nbytes)

      !------------------------------------------------------
      ! Read global file for ADCIRC nodes and determine local
      ! nodes on each PROC; allow for ghost nodes:
      !   1) each source may reside on multiple PROC
      !   2) each PROC store locID for fort.420
      !   3) and LGindX for writing QVal
      !------------------------------------------------------
      print *, " "
      print *, "Interior Hydrology Locations (mesh node)"
      print *, "DOMAIN    numHydro    Global      Local"
      write(*,93) "Global",numHydro

      PE(1:6) = 'PE0000'
      numHydroProc(1:nproc) = 0 ! reset counter
      do I=1, numHydro 
         read(420,*) indX  ! full domain node number
         do J=1, itotproc(indX)    ! loop subdomains to find node
            iproc2 = imap_nod_gl2(2*(J-1)+1,indX) ! find next subdomain
            do iproc=1, nproc
               if (iproc.EQ.iproc2) then ! global node maps to this sd
                  numHydroProc(iproc)=numHydroProc(iproc)+1
                  local=imap_nod_gl2(2*(j-1)+2,indX) !local node
                  imap_hyd_LG(numHydroProc(iproc),iproc) = I
                  locID(I,iproc)=local
                  call Iwrite(PE,3,6,iproc-1)
                  write(*,94) PE, I, indX, local
               end if 
            end do
         end do
      enddo ! loop through numHydro
      close (420)

 93   format(1X,A6,3X,I8)
 94   format(1X,A6,3(3X,I8))

      ! --------------------------------
      ! Write local fort.420 on all PROC
      ! --------------------------------
      do iproc=1,nproc
         ! number of hydroLocs
         write(sdu(iproc),*) numHydroProc(iproc) 
         ! local node number
         do I=1,numHydroProc(iproc)
            indX = imap_hyd_LG(I,iproc)
            write(sdu(iproc),*) locID(indX,iproc)
         end do
         close (sdu(iproc))
      end do

      ! Will read in fort.430 and process when return, so only
      ! deallocate local subroutine arrays - not module arrays
      if (allocated(locID)) then
         deallocate (locID)
         nbytes = 4*nproc*numHydro
         call memory_dealloc(nbytes)
         call memory_status()
      endif

      return

!-----------------------------------------------------------------------
      end subroutine prep420
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> Called from prepInteriorHydrology within adcprep.F 
!> Creates local (PE) interior hydrology input files that have
!> locations specified by position instead of node number. It
!> uses the domain decomposition of the ADCIRC grid created by
!> the routine DECOMP.
!> 
!> Current options include interiorHydrologyType (2 or 3) with files
!> fort.428/429/430. Choice set by passing the iounit (428/429).
!>
!> CMS: edit when change for dynamic unit numbers and filenames
!>      can probably optimize and use some routines from hydrology.F
!>      also if have time to clean up a bit and combine 428/429
!-----------------------------------------------------------------------
      subroutine prepHydroLoc()
!-----------------------------------------------------------------------

      !CMS: add required variables for elements instead of nodes
      use mesh, only : cylindermap
      use adc_constants, only : deg2rad
      use pre_global, only : nproc, nelp, imap_el_LG, ics, sl0, sf0
      use KDTreeTools, only : kdtsearch, setup_kdt_search, kdtSetup
      use memory_usage

      implicit none

      integer :: nbytes = 0
      integer :: I,J,K,KK  ! loop indices
      integer :: iproc,iproc2 !PROC indices
      integer :: sdu(nproc) ! subdomain unit numbers
      integer indX
      logical success    ! .true. if files opened without errors
      real(8) :: latN,lonN  ! for coordinate transformation
      character(6) :: PE

      real(8),allocatable, TARGET ::  xps(:),yps(:),slps(:),sfps(:)
      integer,allocatable :: globID(:)  !global element index


      !----------------------------------------------------------
      !Open and read the global interiorHydrology input locations
      !----------------------------------------------------------
      call OpenPrepHydro(428, 'point source locations  ',               &
     &     1, nproc, sdu, success)
      if (.not.Success) then
         !cms: change to stop since run will fail w/o files
         write(*,*) 'ERROR: Unit 428 file not preprocessed.'
         write(*,*) 'Terminating adcprep.'
         !return
         stop
      endif
 
      ! Read fixed number of sources and allocate
      read(428,*) numHydro
      if(ics.eq.1) then
         allocate ( xps(numHydro),yps(numHydro) )
         nbytes = 8*2*numHydro  ! 2 real array
      else
         allocate ( xps(numHydro),yps(numHydro) )
         allocate ( slps(numHydro),sfps(numHydro) )
         nbytes = 8*4*numHydro  ! 4 real arrays
      end if
      allocate ( globID(numHydro) ) 
      allocate ( numHydroProc(nproc) )
      numHydroProc(1:nproc) = 0 ! set initial counter value 
      allocate ( imap_hyd_LG(numHydro,nproc) ) 

      !Now update memory usage
      nbytes = nbytes + 4*(numHydro+nproc+nproc*numHydro) ! integers
      call memory_alloc(nbytes)

      ! Read in locations and transform if necessary
      do I=1,numHydro
         if(ics.EQ.1) then
            read(428,*) xps(I), yps(I)
         else
            read(428,*) slps(I), sfps(I)
            lonN=deg2rad*slps(I)
            latN=deg2rad*sfps(I)
            ! no rotation considered in adcprep
            call cylindermap(xps(I),yps(I),lonN,latN,sl0,sf0,ics)
         end if
      end do 
      close(428) 

      ! Determine which element each source is in
      if (.not. kdtSetup) call setup_kdt_search()
      call kdtsearch(xps,yps,globID,numHydro,                           &
     &          'Interior hydrology Locations ')
            
      ! Determine which PROC each source is on
      print *, " "
      print *, "Interior Hydrology Locations (element index)"
      print *, "DOMAIN    numHydro    Global      Local"
      write(*,93) "Global",numHydro
      do iproc=1,nproc
        ! Write out to log while processing
        PE(1:6) = 'PE0000'
        call Iwrite(PE,3,6,iproc-1)
        do K = 1,numHydro
           do J=1,NELP(iproc)
              indX = abs(imap_el_LG(J,iproc)) ! global elem indx
              if (indX.EQ.globID(K)) then ! found
                 numHydroProc(iproc) = numHydroProc(iproc) + 1
                 imap_hyd_LG(numHydroProc(iproc),iproc) = K
                 write(*,94) PE, K, globID(K), J
                 exit                  
              endif
           enddo  !end loop thru num elements on PROC
        enddo  ! end loop thru numHydro locations
      end do ! end loop thru nproc
      print *, " "
 93   format(1X,A6,3X,I8)
 94   format(1X,A6,3(3X,I8))

      ! Write each PROC fort.428 file and then close
      do iproc=1,nproc
         ! number of hydroLocs
         write(sdu(iproc),*) numHydroProc(iproc) 
         ! XY or lon/lat coordinates
         do K=1,numHydroProc(iproc)
            indX = imap_hyd_LG(K,iproc)
            !indX = abs(imap_staps_LG(K,iproc))
            if(ics.EQ.1) then
               write(sdu(iproc),*) xps(indX), yps(indX)
            else
               write(sdu(iproc),*) slps(indX), sfps(indX)
            end if
         enddo
         close (sdu(iproc))
      end do 

      !----------------------------------------------------------
      ! Cleanup local arrays before exit - NOT module arrays
      !----------------------------------------------------------
      if(ics.EQ.1) then
         if (allocated(xps)) then
            deallocate (xps,yps)
            nbytes = 8*2*numHydro  ! 2 real array
         end if
      else
         if (allocated(slps)) then
            deallocate (xps,yps,slps,sfps)
            nbytes = 8*4*numHydro  ! 4 real arrays
         end if
      end if
      deallocate (globID)
      nbytes = nbytes + 4*numHydro  ! 1 integer array
      call memory_dealloc(nbytes)
      call memory_status()

      return
!-----------------------------------------------------------------------
      end subroutine prepHydroLoc
!-----------------------------------------------------------------------


!> Called from multiple prep files within hydroprep to open global
!> and local files. Copied from prep.F OpenPrepFiles to remove circular
!> dependency when compiling. Modified somewhat to suit needs of 
!> hydrology and add more descriptive messages.
!-----------------------------------------------------------------------
      subroutine OpenPrepHydro(UnitNumber, Description,  &
     &     startProc, endProc, SDU, Success)
!-----------------------------------------------------------------------
      use pre_global
      implicit none

      integer, intent(in) :: UnitNumber     ! i/o unit number to open
      character(*), intent(in) :: Description ! description of file
      integer, intent(in) :: startProc        ! subdomains to start with
      integer, intent(in) :: endProc          ! subdomain to end on
      integer, intent(out), dimension(nproc) :: sdu ! local unit numbers
      logical, intent(out):: Success     ! .true. if files opened w/o errors
      logical Found               !.true. if the full domain file exists
      character(len=80) FileName   ! name of full domain file
      character(len=80) DefaultName! default name of full domain file
      integer :: dnlen  !length of defaultname (7 or 8)
      integer ErrorIO             ! zero if file opened successfully
      integer iproc               ! subdomain index
      character(len=20) sdFileName     ! subdomain file name
      integer:: unitwidth
      character (len=20):: unitnumberstr

      Found = .false.
      Success = .false.
      ErrorIO = 1
      DefaultName(:) = ' '  !initialize to all blanks
      DefaultName(1:5) = 'fort.'

      unitwidth = ceiling(log(dble(UnitNumber))/log(10.0)) ;        
      dnlen = 5 + unitwidth ;   
      write(unitnumberstr,'(I0)') UnitNumber ;
      DefaultName = trim(adjustl(DefaultName))// &
     &              trim(adjustl(unitNumberStr))

      ! Search for global file and give error message if not found
      ! CMS: removed user opportunity to specify non-default filename
      if (Use_Default) then
         FileName = trim(DefaultName(1:dnlen))
      end if

      ! Determine if full domain file exists
      inquire(file=FileName,exist=found)

      ! If it does exist, open it
      if ( found ) then
         write(*,1011) FileName ! found
         open(unit=UnitNumber, file=FileName, iostat=ErrorIO)
         Success = .true.
         if ( ErrorIO .GT. 0 ) then
            write(*,'(2A)') "ERROR: Full domain file exists but",  &
     &       "  cannot be opened."
            Success = .false.
         end if
      else
      ! Provide message to user if does not exist
         write(*,1010) FileName !not found
         Success = .false.
      end if
      ! return if cannot open global file
      if (.not.Success) return

      ! Open each of the subdomain files
      do iproc = startProc, endProc
         sdFileName(:) = ' '
         sdu(iproc) = 505 + (iproc-1)
         sdFileName(1:7) = 'PE0000/'
         sdFileName(8:) = trim(DefaultName)     
#ifdef ADCSWAN
         sdFileName = 'PE0000/'//FileName
#endif
         CALL IWRITE(sdFileName, 3, 6, iproc-1)
         open (unit=sdu(iproc), file=trim(sdFileName), iostat=ErrorIO)
         Success = .true.
         if ( ErrorIO .GT. 0 ) then
            write(*,*) "ERROR: Subdomain file cannot be opened."
            Success = .false.
            return ! failed to open at least one subdomain file
         end if
      end do

 1010 format('File ',A8,/,' was not found!',/)
 1011 format('File ',A8,/,' was found!  Opening & Processing file.',/)
      return
!------------------------------------------------------------------------
      end subroutine OpenPrepHydro
!------------------------------------------------------------------------


      end MODULE hydroPrep
