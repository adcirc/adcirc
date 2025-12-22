#include "logging_macros.h"

!-----------------------------------------------------------------------
!  MODULE owi_rain
!-----------------------------------------------------------------------
!> @file owi_rain.F
!>
!> @author Kendra Dresback, Univ. of Oklahoma, dresback@ou.edu
!> @author Chris Szpilka, Univ. of Oklahoma, cmszpilka@ou.edu
!>
!> @copyright Dr. R.A. Luettich and Dr. J.J. Westerink
!>
!> @brief This modules handles I/O and interpolation related to
!> Oceanweather format input files for precipitation.
!>
!> @details
!> This file was modified from owi_ice.F by Kendra Dresback
!> (12/27/2013), and owi_ice was modified from owiwind by TCM in 2010.
!> Updates to commenting style, conversion to real(8), and data
!> declarations (public/private) were made by Chris Szpilka in 2025.
!> However, the module is still in need of being updated to match more
!> recent updates to owiwind.
!>
!> Rainfall values are read in from fort.271 and fort.272, which must
!> contain intensity values given in units of [mm/hr]. This module
!> converts these values to units of [m/s] with hard-coded conversions.
!>
!> Expected input files:
!>    fort.27:  control file (same format as wind fort.22 file)
!>        numSets         <-- Number of sets of rain fields
!>        numBlankSnaps   <-- Number of snaps to add/skip at beginning
!>        rainMultiplier  <-- Rain intensity converion NOT USED YET
!>    fort.271: basin precipitation values [mm/hr]
!>    fort.272: regional precipitation values [mm/hr] (optional)
!>
!> Precipitation fields should be listed from coarse to fine in the
!> control file, since the priority of field application is reversed
!> relative to the order given in the control file, with the last
!> field having the highest priority.
!>
!> Precipitation fields are interpolated onto each ADCIRC node using a
!> bilinear interpolation scheme. All nodes located outside of the
!> fields receive a default value of 0.0 mm/hr intensity.
!>
!> @todo Need to think how to merge all owi files into a single module
!> since basic process is the same whether wind speed, pressure, ice or
!> rainfall is being processed. The owiwind module has been updated
!> extensively and the owi_rain and owi_ice files work but are out of
!> date. Will have to decide if separate modules for 1 vs 2 variables
!> or just call a single-variable module twice more than once with the
!> correct file unit and variable information. Note that rain and ice
!> require an additional averaging to account for negative values over
!> land. Getting errors when use formatted write to scratchMessage, so
!> kept the screenUnit for formatted messages with variables.
!>
!-----------------------------------------------------------------------
      module owi_rain
!-----------------------------------------------------------------------

      use sizes, only : myproc
      use mod_logging, only : t_log_scope, init_log_scope,
     &    screenUnit, INFO, allMessage, logMessage
!      use mod_io, only : openFileForRead
      use mod_terminate, only : terminate, ADCIRC_EXIT_FAILURE

      implicit none

      ! Values read from control file: used in init and get
      integer :: numSets, numBlankSnaps

      ! Counters for tracking snaps: initialized in init
      ! and incremented in get
      integer :: cntSnaps, isnapR, isnapB

      ! Flags for when to update grid interpolation factors
      ! initialized in init and updated in get
      integer :: updateB, updateR

      ! Partial datetime from header for regional data only: to track
      ! when need to process regional data (do not read minutes)
      integer :: date1R, date2R

      ! Basin/Regional header values: used to determine when
      ! an interpolation update is required
      ! compare to current in get, update in interp
      integer :: iLatR,iLongR,iCYMDHR,iMinR
      real(8)  :: dxR,dyR,swlatR,swlongR
      integer :: iLatB,iLongB,iCYMDHB,iMinB
      real(8)  :: dxB,dyB,swlatB,swlongB

      ! Current OWI header values: used for B and R
      ! read in get, used to update B/R in interp
      integer,save :: iLatc,iLongc,iCYMDHc,iMinc
      real(8)  :: dxc,dyc,swlatc,swlongc

      ! Interpolation: allocated in init, updated in interp, used in get
      integer,allocatable :: swpointsR(:,:)
      integer,allocatable :: swpointsB(:,:)
      real(8),allocatable :: wR(:,:)
      real(8),allocatable :: wB(:,:)

      ! Arrays for OWI data fields and coordinates (may not need to be
      ! in main module??) used in interpRain and getRain
      real(8), dimension(:,:), allocatable :: rainB,rainR
      real(8), dimension(:), allocatable :: latR,longR,latB,longB

      PUBLIC rain12init, rain12get
      PRIVATE interpRainBasin, interpRainRegion
      !all variables are private
      PRIVATE numSets,numBlankSnaps,cntSnaps,isnapR,isnapB
      PRIVATE swpointsR,swpointsB,wR,wB
      PRIVATE updateB,updateR,date1R,date2R
      PRIVATE iLatc,iLongc,iCYMDHc,iMinc,dxc,dyc,swlatc,swlongc
      PRIVATE iLatR,iLongR,iCYMDHR,iMinR,dxR,dyR,swlatR,swlongR 
      PRIVATE iLatB,iLongB,iCYMDHB,iMinB,dxB,dyB,swlatB,swlongB
      PRIVATE latR,longR,latB,longB,rainB,rainR


      CONTAINS

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! P U B L I C   F U N C T I O N S
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> @brief Initializes reading data from Oceanweather format precip
!fields
!> @param[out] rainval  precipitation rate
!> @param[in]  np       number of nodes in adcirc mesh for processor
!-----------------------------------------------------------------------
      SUBROUTINE rain12init(rainval,NP)
!-----------------------------------------------------------------------

      use sizes, only : gblinputdir

      implicit none

      integer, intent(in) :: np
      real(8), intent(out) :: rainval(:)
   
      integer :: i
      integer :: date1rain,date2rain
      integer :: numSkipSnaps
      character(80) :: owiheader

      LOG_SCOPE_TRACED("rain12init", RAIN_TRACING)

      allocate(swpointsB(NP,2),wB(NP,4))
      allocate(swpointsR(NP,2),wR(NP,4))
      
      ! Read meta info -- todo: edit to use openFileForRead
      open(27,file=TRIM(GBLINPUTDIR)//'/'//'fort.27',status='old',
     & action='read')

      ! Read the number of sets for rain files:
      ! If numSets = 1 then ADCIRC requires UNIT 271
      ! If numSets = 2 then ADCIRC requires UNIT 271 and 272.
      read(27,*,err=99999) numSets

      if(numSets.NE.1.and.numSets.NE.2) then
         call terminate(exit_code=ADCIRC_EXIT_FAILURE,
     &      message='Number of sets in rain control (fort.27) must '
     &       //'be either 1 or 2.')
      endif

      ! Read the number of blank snaps to be inserted before start
      read(27,*,err=99999) numBlankSnaps

      if(numBlankSnaps.LT.0) then
        numSkipSnaps = ABS(numBlankSnaps)
        numBlankSnaps = 0
      else
        numSkipSnaps = 0 
      endif

      ! If want to add USER controlled unit multiplier - add read here
      ! before closing the control file
      close(27)

      ! Read basin rain file header -- todo: edit to use openFileForRead
      open(271,file=trim(GBLINPUTDIR)//'/'//'fort.271', status='old')

      ! Read begining/ending dates of rain file
      owiheader(:) = ' '  !set owiheader to blanks before read
      read(271,fmt='(a80)',err=99999) owiheader
      read(owiheader(52:61),'(I10)') date1rain
      read(owiheader(69:78),'(I10)') date2rain

      if (myproc == 0)then
         write(screenunit,*)'Basin rainfall: startDate = ',date1rain
         write(screenunit,*)'Basin rainfall: endDate = ',date2rain
      endif

      ! Check if region scale data exist
      IF(numSets.eq.1) GOTO 100

      ! Read region rain file header -- todo: edit to use
      ! openFileForRead
      open(272,file=TRIM(GBLINPUTDIR)//'/'//'fort.272',status='old')

      ! Read begining/ending dates of rain file
      owiheader(:) = ' '  !set owiheader to blanks before read
      read(272,fmt='(a80)',err=99999) owiheader
      read(owiheader(52:61),'(I10)') date1rain
      read(owiheader(69:78),'(I10)') date2rain

      date1R = date1rain
      date2R = date2rain

      if (myproc == 0)then
         write(screenunit,*)'Region rainfall: startDate = ',date1rain
         write(screenunit,*)'Region rainfall: endDate = ',date2rain
      endif

 100  continue


      ! Initialize flags
      isnapB = 0
      isnapR = 0
      updateB = 1
      updateR = 1
      cntSnaps = 0

      ! Skip snaps if necessary
      do i = 1,numSkipSnaps
        call rain12get(rainval,NP)
      enddo

      return

99999 continue
      call terminate(exit_code=ADCIRC_EXIT_FAILURE,
     &   message='OWI_rain read error during initialization.')

!-----------------------------------------------------------------------
      END SUBROUTINE rain12init
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> @brief Reads rainfall intensity from Oceanweather format file and
!>        interpolates it to the ADCIRC mesh
!> @param[out] rainval  rainfall intensity (assumed mm/hr)
!> @param[in]  np       number of nodes in adcirc mesh for processor
!-----------------------------------------------------------------------
      SUBROUTINE rain12get(rainval,NP)
!-----------------------------------------------------------------------

      implicit none

      integer, intent(in) :: np
      real(8), intent(out) :: rainval(:)

      integer :: I,J,XI,YI,K,ibar
      real(8) :: CV(4),CVbar
      real(8) :: rain1  ! scalar placeholder during interpolation
      real(8) :: w,w1,w2,w3,w4  !local weights
      logical :: regionExists
      character(1024) :: scratchMessage=" "  ! wont work with formatting
 
      LOG_SCOPE_TRACED("rain12get", RAIN_TRACING)


      ! Formats for OWI reading
 11   format(t6,i4,t16,i4,t23,f6.0,t32,f6.0,
     &     t44,f8.0,t58,f8.0,t69,i10,i2)
 22   format(8f10.0)



      ! BEGIN BASIN-SCALE DATA

      ! Increment counter (cntSnaps initialized to zero in rain12init)
      cntSnaps = cntSnaps+1

      ! Put a blank snap for the first 'numBlankSnaps' snaps
      if(cntSnaps.LE.numBlankSnaps) then
        do i=1,NP
          rainval(I)=0.d0
        enddo
        if(myproc.eq.0) write(screenUnit,12) cntSnaps
      endif

      ! Increment counter
      isnapB = isnapB+1

      ! Read grid specifications/date in rain file
      read (271,11,end=10000,err=99999) 
     &      iLatc,iLongc,dxc,dyc,swlatc,swlongc,iCYMDHc,iMinc

      ! Check if header info has changed from the previous snapshot
      if(isnapB.gt.1) then
        if(iLatc.ne.iLatB.or.iLongc.ne.iLongB.or.dxc.ne.dxB.or.
     &       dyc.ne.dyB.or.swlatc.ne.swlatB.or.
     &       swlongc.ne.swlongB) then
          updateB = 1
        else
          updateB = 0
        endif
      endif
      
      iCYMDHB = iCYMDHc
      iMinB = iMinc

      ! Update coordinate mapping coefficients and reallocate arrays
      if(updateB.eq.1) then
        call interpRainBasin(np)
      endif
      
      ! Read rain intensities
      read(271,22,end=9999,err=99999) ((rainB(i,j),i=1,iLongB),j=1,iLatB)


      ! BEGIN REGION-SCALE DATA
      regionExists = .FALSE.

      IF(numSets.EQ.1) GOTO 100
      if(iCYMDHB.lt.date1R) goto 100
      if(iCYMDHB.eq.date2R.and.iMinR.ne.0) goto 100
      if(iCYMDHB.gt.date2R) goto 100

      regionExists = .TRUE.

      ! Increment counter
      isnapR = isnapR+1

      ! Read grid specifications/date in rain file
      read (272,11,end=10000,err=99999) 
     &     iLatc,iLongc,dxc,dyc,swlatc,swlongc,iCYMDHc,iMinc

      ! Check if header info has changed from the previous snapshot
      if(isnapR.gt.1) then
        if(iLatc.ne.iLatR.or.iLongc.ne.iLongR.or.dxc.ne.dxR.or. 
     &       dyc.ne.dyR.or.swlatc.ne.swlatR.or.
     &       swlongc.ne.swlongR) then
          updateR = 1
        else
          updateR = 0
        endif
      endif

      iCYMDHR = iCYMDHc
      iMinR = iMinc

      if(iCYMDHB.ne.iCYMDHR.or.iMinB.ne.iMinR) then
         call terminate(exit_code=ADCIRC_EXIT_FAILURE,
     &     message='Rainfall snapshots in fort.271 and fort.272 '
     &       //'are not syncronized.')
      endif

      ! Update coordinate mapping coefficients and reallocate arrays
      if(updateR.eq.1) then
        call interpRainRegion(np)
      endif

      ! Read rainfall intensities
      read(272,22,end=9999,err=99999) ((rainR(i,j),i=1,iLongR),j=1,iLatR)


 100  CONTINUE

      ! Interpolate onto ADCIRC grid and write to logfile
      ! cms: removed ramp      
      if(regionExists) then
         if (myproc == 0 ) write(screenUnit,13) iCYMDHB, iMinB
      else
         if (myproc == 0 ) write(screenUnit,14) iCYMDHB, iMinB
      endif

      do i=1,NP

        rain1=-9999.9d0

        ! BASIN
        if (swpointsB(i,1).gt.0) then
          xi = swpointsB(i,1)
          yi = swpointsB(i,2)

          w1=wB(i,1)
          w2=wB(i,2)
          w3=wB(i,3)
          w4=wB(i,4)
          ! If a rainfall intensity is less than zero (land flag), then
          ! set to the average non-zero value before interpolating
          CV(1) = rainB(xi,yi)
          CV(2) = rainB(xi+1,yi)
          CV(3) = rainB(xi+1,yi+1)
          CV(4) = rainB(xi,yi+1)
          CVbar = 0.d0
          ibar = 0
          do k=1,4
             if(CV(k) .ge. 0.d0) then
                CVbar = CVbar + CV(k)
                ibar = ibar + 1
             endif
          enddo
          if (ibar.gt.0) CVbar = CVbar/real(ibar,8)  ! is real needed?
          do k=1,4
             if (CV(k).lt.0.d0) CV(k) = CVbar
          enddo

          ! perform the interpolation
          rain1=w1*CV(1)+w2*CV(2)+w3*CV(3)+w4*CV(4)
        endif

        ! REGION ---------------------------------------------------------
        ! rain1 will be overwritten if region data exist.
        if ((regionExists).and.(swpointsR(i,1).gt.0)) then
          xi = swpointsR(i,1)
          yi = swpointsR(i,2)

          w1=wR(i,1)
          w2=wR(i,2)
          w3=wR(i,3)
          w4=wR(i,4)
        ! If a rainfall value is less than zero (land flag)
        ! then set the value to the average non-zero value before
        ! interpolating
          CV(1) = rainR(xi,yi)
          CV(2) = rainR(xi+1,yi)
          CV(3) = rainR(xi+1,yi+1)
          CV(4) = rainR(xi,yi+1)
          CVbar = 0.d0
          ibar = 0

          do k=1,4
             if(CV(k) .ge. 0.d0) then
                CVbar = CVbar + CV(k)
                ibar = ibar + 1
             endif
          enddo
          if (ibar.gt.0) CVbar = CVbar/real(ibar,8)   !is real necessary
          do k=1,4
             if (CV(k).lt.0.d0) CV(k) = CVbar
          enddo

          ! perform the interpolation
          rain1=w1*CV(1)+w2*CV(2)+w3*CV(3)+w4*CV(4)
        endif

        ! COPY TO ARRAYS (with unit conversion)

        if(rain1.eq.-9999.9d0) then
          rainval(I)=0.d0
        else
          ! Convert the rainfall from mm/hr to m/s
          ! More flexibility if add rainMultiplier in fort.27
          rainval(i) = rain1*0.00000027778d0
        end if

      enddo   ! loop through NP

      return


      ! First I/O error type: fatal
9999  continue

      call terminate(exit_code=ADCIRC_EXIT_FAILURE,
     &   message='Unexpected EOF in unit 271 or 272.')

      ! Second I/O error type: non-fatal
10000 continue

      write(scratchMessage,'(" Either unit 271 or 272 ran out of data.",
     &    /," Execution will continue with rainfall set to zero.")')
      call allMessage(INFO,scratchMessage)

      do i=1,NP
        rainval(I)=0.d0
      enddo

      return

      
99999 continue

      call terminate(exit_code=ADCIRC_EXIT_FAILURE,
     &   message='OWI_rain read error during timestepping.')


 12   format(' Inserting a blank rainfall snap, count=',i4)
 13   format(' Processing basin and regional scale rainfall data ',
     & 'for time stamp ',i12,' ',i2)
 14   format(' Processing basin-scale rainfall data for time stamp ',
     &       i12,' ',i2)

!-----------------------------------------------------------------------
      END SUBROUTINE rain12get
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! P R I V A T E     F U N C T I O N S
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> @brief Generates interpolation coefficients for mapping from a
!>        basin-scale Oceanweather formatted file to the ADCIRC mesh
!> @param[in]  np       number of nodes in adcirc mesh for processor
!-----------------------------------------------------------------------
      SUBROUTINE interpRainBasin(np)
!-----------------------------------------------------------------------

      use ADC_CONSTANTS,only : rad2deg
      use MESH, only : SLAM, SFEA

      implicit none

      integer, intent(in) :: np

      integer :: I,J,K,XI,YI
      real(8) :: adcLat,adcLong
      character(1024) :: scratchMessage = ""

      real(8) :: w,w1,w2,w3,w4  !local weights

      LOG_SCOPE_TRACED("interpRainBasin", RAIN_TRACING)

      write(scratchMessage,'("Updating basin-scale rainfall mapping.")')
      call logMessage(INFO,scratchMessage)

      ! update to current values
      iLatB = iLatc
      iLongB = iLongc
      dxB = dxc
      dyB = dyc
      swlatB = swlatc
      swlongB = swlongc
        
      ! Allocate and create matrices
      if(allocated(longB)) deallocate(longB)
      if(allocated(latB))  deallocate(latB)
      allocate(longB(iLongB),latB(iLatB))
      if(allocated(rainB))  deallocate(rainB)
      allocate(rainB(iLongB,iLatB))
        
      ! Generate long&lat on each grid point
      do i=1,iLatB
        latB(i) = swlatB+(i-1)*dyB
      enddo
      do i=1,iLongB
        longB(i) = swlongB+(i-1)*dxB
      enddo
        
      ! Generate interpolation coefficients (south west point and weights)
      do i=1,NP
        adcLat  = rad2deg*SFEA(i)
        adcLong = rad2deg*SLAM(i)

        if (adcLong>=longB(1).and.adcLong<longB(iLongB) .and. 
     &       adcLat>=latB(1).and.adcLat<latB(iLatB)) then
          do j=1,iLongB-1
            if (adcLong>=longB(j) .and.  
     &           adcLong<longB(j+1)) then
              xi=j
              goto 200
            endif
          enddo
 200      continue
            
          do k=1,iLatB-1
            if (adcLat>=latB(k) .and. 
     &           adcLat<latB(k+1)) then
              yi=k
              goto 300
            endif
          enddo
 300      continue
            
          swpointsB(i,1) = xi
          swpointsB(i,2) = yi
            
          w=(longB(xi+1)-longB(xi))*(latB(yi+1)-latB(yi))
          w1=(longB(xi+1)-adcLong)*(latB(yi+1)-adcLat)
          w2=(adcLong-longB(xi))*(latB(yi+1)-adcLat)
          w3=(adcLong-longB(xi))*(adcLat-latB(yi))
          w4=(longB(xi+1)-adcLong)*(adcLat-latB(yi))
            
          wB(i,1)=w1/w
          wB(i,2)=w2/w
          wB(i,3)=w3/w
          wB(i,4)=w4/w

        else
          swpointsB(i,1) = 0
          swpointsB(i,2) = 0
            
        endif
      enddo

      return
!-----------------------------------------------------------------------
      END SUBROUTINE interpRainBasin
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> @brief Generates interpolation coefficients for mapping from a
!>        region-scale Oceanweather formatted file to the ADCIRC mesh
!> @param[in]  np       number of nodes in adcirc mesh for processor
!-----------------------------------------------------------------------
      SUBROUTINE interpRainRegion(np)
!-----------------------------------------------------------------------
      use ADC_CONSTANTS,only : rad2deg
      use MESH, only : SLAM, SFEA
      implicit none

      integer, intent(in) :: np
      integer :: I,J,K,XI,YI
      real(8) :: adcLat,adcLong
      character(1024) :: scratchMessage = ""

      real(8) :: w,w1,w2,w3,w4  !local weights - put in each subroutine

      LOG_SCOPE_TRACED("interpRainRegion", RAIN_TRACING)

      write(scratchMessage,'("Updating region rainfall mapping.")')
      call logMessage(INFO,scratchMessage)

      ! Update to current values
      iLatR = iLatc
      iLongR = iLongc
      dxR = dxc
      dyR = dyc
      swlatR = swlatc
      swlongR = swlongc
 
      ! Allocate and create matrices
      if(allocated(longR)) deallocate(longR)
      if(allocated(latR))  deallocate(latR)
      allocate(longR(iLongR),latR(iLatR))
      if(allocated(rainR))  deallocate(rainR)
      allocate(rainR(iLongR,iLatR))

      ! Generate long&lat on each grid point
      do i=1,iLatR
        latR(i) = swlatR+(i-1)*dyR
      enddo
      do i=1,iLongR
        longR(i) = swlongR+(i-1)*dxR
      enddo

      ! Generate interpolation info (south west point and weights)
      do i=1,NP
        adcLat  = rad2deg*SFEA(i)
        adcLong = rad2deg*SLAM(i)

        if (adcLong>=longR(1).and.adcLong<longR(iLongR) .and. 
     &       adcLat>=latR(1).and.adcLat<latR(iLatR)) then
          do j=1,iLongR-1
            if (adcLong>=longR(j).and.adcLong<longR(j+1)) then
              xi=j
              goto 201
            endif
          enddo
 201      continue

          do k=1,iLatR-1
            if (adcLat>=latR(k).and.adcLat<latR(k+1)) then
              yi=k
              goto 301
            endif
          enddo
 301      continue

          swpointsR(i,1) = xi
          swpointsR(i,2) = yi

          w=(longR(xi+1)-longR(xi))*(latR(yi+1)-latR(yi))
          w1=(longR(xi+1)-adcLong)*(latR(yi+1)-adcLat)
          w2=(adcLong-longR(xi))*(latR(yi+1)-adcLat)
          w3=(adcLong-longR(xi))*(adcLat-latR(yi))
          w4=(longR(xi+1)-adcLong)*(adcLat-latR(yi))

          wR(i,1)=w1/w
          wR(i,2)=w2/w
          wR(i,3)=w3/w
          wR(i,4)=w4/w
        else
          swpointsR(i,1) = 0
          swpointsR(i,2) = 0
        endif
      enddo

      return
!-----------------------------------------------------------------------
      END SUBROUTINE interpRainRegion
!-----------------------------------------------------------------------


      END MODULE OWI_RAIN

