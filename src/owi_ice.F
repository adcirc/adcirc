C******************************************************************************
C  MODULE OWI_ICE
C    Written by Chris Massey, USACE-ERDC-CHL 12/15/2010
C    Modified from owiwind.F Written by s.b. 08/17/2006
C
C    The values read in from the fort.225/fort.227 contain ice percent
C    concentration values that range from 0.0 to 100.0 percent.  A value
C    of -1.0 is used to signal that either the cell is over land, or no ice
C    data was available.
C
C****************************************************************************** 
C 
      MODULE OWI_ICE
      use CkLunMod, only : CK_LUN
C
      USE SIZES,ONLY : SZ,MYPROC
      USE GLOBAL,ONLY : NSCREEN, ScreenUnit

! DW AMPI
!#ifdef CMPI
!      USE MESSENGER,ONLY : MSG_FINI
!#endif 
!
        
      IMPLICIT NONE
      character*100 :: title
!$omp threadprivate(title)

!      tcm v48.4645 Changed the way the header line is read
!      character :: part1*12,part2*7,part3*6
      character :: owiheader*80
!$omp threadprivate(owiheader)

      real(SZ), dimension(:,:), allocatable :: CiceR,CiceB
!$omp threadprivate(CiceR,CiceB)
      real(SZ), dimension(:), allocatable :: latR,longR,latB,longB
!$omp threadprivate(latR,longR,latB,longB)
      real(SZ) :: Along, Alat
!$omp threadprivate(Along,Alat)
      real(SZ) :: ramp,rampfrac
!$omp threadprivate(ramp,rampfrac)
      real(SZ) :: Cice1
!$omp threadprivate(Cice1)
      real(SZ) :: Cice_env
!$omp threadprivate(Cice_env)
      integer(8):: date1R,date2R,date1B,date2B
!$omp threadprivate(date1R,date2R,date1B,date2B)
      integer(8):: date1ice,date2ice
!$omp threadprivate(date1ice,date2ice)
      integer :: iLatR,iLongR,iCYMDHR,iMinR
!$omp threadprivate(iLatR,iLongR,iCYMDHR,iMinR)
      real(SZ)  :: dxR,dyR,swlatR,swlongR
!$omp threadprivate(dxR,dyR,swlatR,swlongR)
      integer :: iLatB,iLongB,iCYMDHB,iMinB
!$omp threadprivate(iLatB,iLongB,iCYMDHB,iMinB)
      real(SZ)  :: dxB,dyB,swlatB,swlongB
!$omp threadprivate(dxB,dyB,swlatB,swlongB)
      integer,save :: iLatc,iLongc,iCYMDHc,iMinc
!$omp threadprivate(iLatc,iLongc,iCYMDHc,iMinc)
      real(SZ)  :: dxc,dyc,swlatc,swlongc
!$omp threadprivate(dxc,dyc,swlatc,swlongc)
      integer :: isnapR,updateR
!$omp threadprivate(isnapR,updateR)
      integer :: isnapB,updateB
!$omp threadprivate(isnapB,updateB)

      logical :: regionExists
!$omp threadprivate(regionExists)

      integer,allocatable :: swpointsR(:,:)
!$omp threadprivate(swpointsR)
      integer,allocatable :: swpointsB(:,:)
!$omp threadprivate(swpointsB)

      real(SZ) :: w,w1,w2,w3,w4
!$omp threadprivate(w,w1,w2,w3,w4)
      real(SZ),allocatable :: wR(:,:)
!$omp threadprivate(wR)
      real(SZ),allocatable :: wB(:,:)
!$omp threadprivate(wB)

      CHARACTER FNAME*60
!$omp threadprivate(FNAME)

      integer :: numSets,numBlankSnaps,cntSnaps,numSkipSnaps
!$omp threadprivate(numSets,numBlankSnaps,cntSnaps,numSkipSnaps)

      PUBLIC

C---------------------end of data declarations--------------------------------C


      CONTAINS


C***********************************************************************
C   SOBROUTINE NCICE1_INIT
C***********************************************************************

      SUBROUTINE NCICE1_INIT(Ciceval,NP)

      USE SIZES, ONLY : SZ,MYPROC, GBLINPUTDIR

      IMPLICIT NONE
      INTEGER   NP,I
      REAL(SZ)  Ciceval(*)

      allocate(swpointsB(NP,2),wB(NP,4))
      allocate(swpointsR(NP,2),wR(NP,4))
      
! Read meta info -------------------------------------------------

      OPEN(25+CK_LUN,FILE=TRIM(GBLINPUTDIR)//'/'//'fort.25',STATUS='OLD',
     & ACTION='READ')

      ! Read the number of sets of .ice files.
      ! If numSets = 1 then ADCIRC requires UNIT 225
      ! If numSets = 2 then ADCIRC requires UNIT 225 and 227.
      !   UNIT 225 and 227 are ice concentration fields

      read(25+CK_LUN,*,err=99999) numSets

      if(numSets.NE.1.AND.numSets.NE.2) then
        if (myproc == 0) write(screenunit+CK_LUN,1004)
        write(16+CK_LUN,1004)\

! DW
!#ifdef CMPI
!        call msg_fini()
!#endif
!
        stop
      endif

      ! Read the number of blank snaps to be inserted before OWI ice fields start
      read(25+CK_LUN,*,err=99999) numBlankSnaps

      ! If numBlankSnaps < 0, ABS(numBlankSnaps) snaps in OWI ice files (UNIT225 and 227) will be skipped.
      if(numBlankSnaps.LT.0) then
        numSkipSnaps = ABS(numBlankSnaps)
        numBlankSnaps = 0

      ! v48.4628 TCM 10/28/2009 -- Added else to initialize numSkipSnaps to be 0 
      else
        numSkipSnaps = 0 
      endif

      close(25+CK_LUN)

! Read basin ice file header ------------------------------------------------

      OPEN(225+CK_LUN,FILE=TRIM(GBLINPUTDIR)//'/'//'fort.225', STATUS='OLD')

      ! Read begining/ending dates of ice file
!      tcm v48.4645 changed the way the header is read
!      read(225,11,err=99999,ADVANCE="NO")part1,part2,part3
! 11   FORMAT(a,1x,a,1x,a)
!      read(225,*)date1ice,date2ice
      owiheader(:) = ' '  !set owiheader to blanks before read
      read(225+CK_LUN,fmt='(a80)',err=99999) owiheader
      read(owiheader(56:65),'(I10)') date1ice
      read(owiheader(71:80),'(I10)') date2ice

      if (myproc == 0)then
         write(screenunit+CK_LUN,*)'date1 in basin ice file = ',date1ice
         write(screenunit+CK_LUN,*)'date2 in basin ice file = ',date2ice
      endif

      date1B = date1ice
      date2B = date2ice

! Check if region scale data exist

      IF(numSets.eq.1) GOTO 100

! Read region ice file header -----------------------------------------------

      OPEN(227+CK_LUN,FILE=TRIM(GBLINPUTDIR)//'/'//'fort.227',STATUS='OLD')

      ! Read begining/ending dates of ice file
!      tcm v48.4645 changed the way the header is read
!      read(227,'(a,1x,a,1x,a)',err=99999,ADVANCE="NO")part1,part2,part3
!      read(227,*)date1ice,date2ice
      owiheader(:) = ' '  !set owiheader to blanks before read
      read(227+CK_LUN,fmt='(a80)',err=99999) owiheader
      read(owiheader(56:65),'(I10)') date1ice
      read(owiheader(71:80),'(I10)') date2ice

      if (myproc == 0)then
         write(screenunit+CK_LUN,*)'date1 in region ice file = ',date1ice
         write(screenunit+CK_LUN,*)'date2 in region ice file = ',date2ice
      endif

      date1R = date1ice
      date2R = date2ice

 100  CONTINUE


! Initialize flags ----------------------------------------------------------
      isnapB = 0
      isnapR = 0
      updateB = 1
      updateR = 1
      cntSnaps = 0

! Skip snaps if necessary ---------------------------------------------------
      do i = 1,numSkipSnaps
        CALL NCICE1_GET(Ciceval,NP)
      enddo

      RETURN


 10   format(a51,i12,a5,i12)
 1004 FORMAT(//,1X,' NUMBER OF SETS WAS SPECIFIED'//
     &       'INCORRECTLY IN UNIT22.'/
     &       ' IT MUST BE ETEHR OF 1 or 2'/
     &       ' EXECUTION WILL BE TERMINATED.'//)

99999 CONTINUE

! DW
!#ifdef CMPI
!      call msg_fini()
!#endif
      STOP 'OWI_ICE READ ERROR (1)'

      END SUBROUTINE

C***********************************************************************
C   SOBROUTINE NCICE1_GET
C***********************************************************************

      SUBROUTINE NCICE1_GET(Ciceval,NP)
      USE SIZES,ONLY : MYPROC,MNPROC

      IMPLICIT NONE

      INTEGER   NP,I,J,XI,YI,K,ibar
      REAL(SZ)  Ciceval(*)
      REAL(SZ)  CV(4),CVbar
      CHARACTER*80 PBLJAGF


! Read basin data  ---------------------------------------------------------

      ! Increment counter (cntSnaps initialized to zero in ncice1_init)
      cntSnaps = cntSnaps+1

      ! Put a blank snap for the first 'numBlankSnaps' snaps
      if(cntSnaps.LE.numBlankSnaps) then
        do i=1,NP
          Ciceval(I)=0.d0
        enddo
        IF(NSCREEN.GE.1) THEN
          if(MYPROC.EQ.0) then
            WRITE(screenunit+CK_LUN,16) cntSnaps
          endif
          WRITE(16+CK_LUN,16) cntSnaps    !TCMv48.4629 (Changed format number from 15 to 16)
        ENDIF
        return
      endif


      ! Increment counter
      isnapB = isnapB+1

      ! Read grid specifications/date in ice file
      read (225+CK_LUN,11,end=10000,err=99999) 
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
      
      ! Update coordinate mapping coefficients if necessary
      if(updateB.eq.1) then
        call ncice1_interp_basin(np)
      endif
      
      ! Read ice concentrations
      read(225+CK_LUN,22,end=9999,err=99999) ((CiceB(i,j),i=1,iLongB),j=1,iLatB)

! Read region data  --------------------------------------------------------

      regionExists = .FALSE.

      IF(numSets.EQ.1) GOTO 100
      if(iCYMDHB.lt.date1R) goto 100
      if(iCYMDHB.eq.date2R.and.iMinR.ne.0) goto 100
      if(iCYMDHB.gt.date2R) goto 100

      regionExists = .TRUE.

      ! Increment counter
      isnapR = isnapR+1

      ! Read grid specifications/date in ice file
      read (227+CK_LUN,11,end=10000,err=99999) 
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
        if (myproc == 0) then
          WRITE(screenunit+CK_LUN,*) 'SNAPSHOTS NOT SYNCRONIZED'
          WRITE(screenunit+CK_LUN,*) '  iCYMDHB=',iCYMDHB, ' iMinB=',iMinB
          WRITE(screenunit+CK_LUN,*) '  iCYMDHR=',iCYMDHR, ' iMinR=',iMinR
          WRITE(screenunit+CK_LUN,*) 'EXECUTION WILL BE TERMINATED'
        endif
        WRITE(16+CK_LUN,*) 'SNAPSHOTS NOT SYNCRONIZED'
        WRITE(16+CK_LUN,*) '  iCYMDHB=',iCYMDHB, ' iMinB=',iMinB
        WRITE(16+CK_LUN,*) '  iCYMDHR=',iCYMDHR, ' iMinR=',iMinR
        WRITE(16+CK_LUN,*) 'EXECUTION WILL BE TERMINATED'

! DW
!#ifdef CMPI
!        call msg_fini()
!#endif
!
        STOP
      endif

      ! Update coordinate mapping coefficients if necessary
      if(updateR.eq.1) then
        call ncice1_interp_region(np)
      endif

      ! Read ice concentrations
      read(223+CK_LUN,22,end=9999,err=99999) ((CiceR(i,j),i=1,iLongR),j=1,iLatR)


 100  CONTINUE

! Interpolate onto ADCIRC grid and write to file -------------------------
      
      rampfrac = isnapB-1
c      if (rampfrac<36) then
c        ramp = tanh(18d0*rampfrac/36d0)
c      end if
      ramp = 1.0

      IF(NSCREEN.GE.1) THEN
        if(regionExists.EQV..TRUE.) then
          if(MYPROC.EQ.0) then
            WRITE(screenunit+CK_LUN,15) iCYMDHB,iMinB
          endif
          WRITE(16+CK_LUN,15) iCYMDHB,iMinB
        else
          if(MYPROC.EQ.0) then
            WRITE(screenunit+CK_LUN,14) iCYMDHB,iMinB
          endif
          WRITE(16+CK_LUN,14) iCYMDHB,iMinB
        endif
      ENDIF

      do i=1,NP

        Cice1=-9999.9D0

        ! BASIN ---------------------------------------------------------

        if (swpointsB(i,1).gt.0) then
          xi = swpointsB(i,1)
          yi = swpointsB(i,2)

          w1=wB(i,1)
          w2=wB(i,2)
          w3=wB(i,3)
          w4=wB(i,4)
!         If an ice concentration value is less than zero (land flag)
!         then set the value to the average non-zero value before
!         interpolating
          CV(1) = CiceB(xi,yi)
          CV(2) = CiceB(xi+1,yi)
          CV(3) = CiceB(xi+1,yi+1)
          CV(4) = CiceB(xi,yi+1)
          CVbar = 0.d0
          ibar = 0
          do k=1,4
             if(CV(k) .ge. 0.d0) then
                CVbar = CVbar + CV(k)
                ibar = ibar + 1
             endif
          enddo
          if (ibar.gt.0) CVbar = CVbar/real(ibar,sz)
          do k=1,4
             if (CV(k).lt.0.d0) CV(k) = CVbar
          enddo

!         perform the interpolation
          Cice1=w1*CV(1)+w2*CV(2)+w3*CV(3)+w4*CV(4)

        endif

        ! REGION ---------------------------------------------------------
        ! Cice1 will be overwritten if region data exist.

        if ((regionExists).and.(swpointsR(i,1).gt.0)) then
          xi = swpointsR(i,1)
          yi = swpointsR(i,2)

          w1=wR(i,1)
          w2=wR(i,2)
          w3=wR(i,3)
          w4=wR(i,4)
!         If an ice concentration value is less than zero (land flag)
!         then set the value to the average non-zero value before
!         interpolating
          CV(1) = CiceR(xi,yi)
          CV(2) = CiceR(xi+1,yi)
          CV(3) = CiceR(xi+1,yi+1)
          CV(4) = CiceR(xi,yi+1)
          CVbar = 0.d0
          ibar = 0
          do k=1,4
             if(CV(k) .ge. 0.d0) then
                CVbar = CVbar + CV(k)
                ibar = ibar + 1
             endif
          enddo
          if (ibar.gt.0) CVbar = CVbar/real(ibar,sz)
          do k=1,4
             if (CV(k).lt.0.d0) CV(k) = CVbar
          enddo

!         perform the interpolation
          Cice1=w1*CV(1)+w2*CV(2)+w3*CV(3)+w4*CV(4)

        endif

        ! COPY TO ARRAYS ----------------------------------------------------------

        if(Cice1.eq.-9999.9D0) then
          Ciceval(I)=0.d0
        else
          if (rampfrac<36) then
            Cice1=Cice_env-(Cice_env-Cice1)*ramp
          endif
          Ciceval(i) = Cice1
        end if

      enddo

      RETURN

9999  continue

      if (myproc == 0) then
      WRITE(screenunit+CK_LUN,*)''
      WRITE(screenunit+CK_LUN,*)' !!!!!FATAL ERROR!!!!!'
      WRITE(screenunit+CK_LUN,*)'EITHER UNIT 225 or 227',
     & ' COULD NOT BE READ'
      WRITE(screenunit+CK_LUN,*)' EXECUTION WILL BE TERMINATED'
      WRITE(screenunit+CK_LUN,*)''
      endif

      WRITE(16+CK_LUN,*) ''
      WRITE(16+CK_LUN,*) ' !!!!!FATAL ERROR!!!!!'
      WRITE(16+CK_LUN,*) ' EITHER UNIT 225 or 227 COULD NOT BE READ'
      WRITE(16+CK_LUN,*) ' EXECUTION WILL BE TERMINATED'
      WRITE(16+CK_LUN,*) ''
!  DW
!#ifdef CMPI
!      call msg_fini()
!#endif
!
      STOP

10000 continue

      IF(MYPROC.EQ.0) THEN
        WRITE(screenunit+CK_LUN,*) ''
        WRITE(screenunit+CK_LUN,*) ' !!! WARNING !!!'
        WRITE(screenunit+CK_LUN,*) ' EITHER UNIT 225 or 227 RAN OUT'
        WRITE(screenunit+CK_LUN,*) ' EXECUTION WILL CONTINUE'
        WRITE(screenunit+CK_LUN,*) ''
      ENDIF

      WRITE(16+CK_LUN,*) ''
      WRITE(16+CK_LUN,*) ' !!! WARNING !!!'
      WRITE(16+CK_LUN,*) ' EITHER UNIT 225 or 227 RAN OUT'
      WRITE(16+CK_LUN,*) ' EXECUTION WILL CONTINUE'
      WRITE(16+CK_LUN,*) ''

      do i=1,NP
        Ciceval(I)=0.d0
      enddo

      RETURN


 11   format(t6,i4,t16,i4,t23,f6.0,t32,f6.0,
     &     t44,f8.0,t58,f8.0,t69,i10,i2)
 12   format(1X,'SNAPSHOT HEADER IN ICE FILES DO NOT MATCH')
 13   format(1X,'EXECUTION WILL BE TERMINATED')
 14   format(/,1X,'PROCESSING BASIN-SCALE ICE DATA',i12,' ',i2)
 15   format(/,1X,'PROCESSING BASIN&REGION-SCALE DATA',i12,' ',i2)
 16   format(/,1X,'INSERTING A BLANK ICE SNAP, COUNT=',i4)
 
 22   format(8f10.0)
      
99999 CONTINUE

! DW
!#ifdef CMPI
!      call msg_fini()
!#endif
      STOP 'OWI_ICE READ ERROR (2)'

      END SUBROUTINE

C***********************************************************************
C   SUBROUTINE NCICE1_INTERP_BASIN
C
C   This generates and saves interpolation coefficients for mapping
C   from a basin-scale OWI to a ADCIRC grid.
C
C***********************************************************************

      SUBROUTINE NCICE1_INTERP_BASIN(NP)
      USE GLOBAL,ONLY : RAD2DEG
      USE MESH, ONLY : SLAM, SFEA
      IMPLICIT NONE

      INTEGER   NP,I,J,K,XI,YI
      REAL(SZ)  adcLat,adcLong

      WRITE(16+CK_LUN,*) ''
      WRITE(16+CK_LUN,*) 'BASIN-SCALE ICE MAPPING UPDATED'
      WRITE(16+CK_LUN,*) ''
        
      iLatB = iLatc
      iLongB = iLongc
      dxB = dxc
      dyB = dyc
      swlatB = swlatc
      swlongB = swlongc
        
      ! Allocate and create matrices
      if(allocated(CiceB))    deallocate(CiceB)
      if(allocated(longB)) deallocate(longB)
      if(allocated(latB))  deallocate(latB)

      allocate(CiceB(iLongB,iLatB))
      allocate(longB(iLongB),latB(iLatB))
        
      ! Generate long&lat on each grid point
      do i=1,iLatB
        latB(i) = swlatB+(i-1)*dyB
      enddo
      do i=1,iLongB
        longB(i) = swlongB+(i-1)*dxB
      enddo
        
      ! Generate interpolation coefficients (south west point and weights)
      do i=1,NP
        adcLat  = RAD2DEG*SFEA(i)
        adcLong = RAD2DEG*SLAM(i)

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
      END SUBROUTINE


C***********************************************************************
C   SUBROUTINE NCICE1_INTERP_REGION
C
C   This generates and saves interpolation coefficients for mapping
C   from a region-scale OWI to a ADCIRC grid.
C
C***********************************************************************

      SUBROUTINE NCICE1_INTERP_REGION(NP)
      USE GLOBAL, ONLY : RAD2DEG
      USE MESH, ONLY : SLAM,SFEA
      IMPLICIT NONE

      INTEGER   NP,I,J,K,XI,YI
      REAL(SZ)  adcLat,adcLong

      WRITE(16+CK_LUN,*) ''
      WRITE(16+CK_LUN,*) 'REGION-SCALE ICE MAPPING UPDATED'
      WRITE(16+CK_LUN,*) ''
        
      iLatR = iLatc
      iLongR = iLongc
      dxR = dxc
      dyR = dyc
      swlatR = swlatc
      swlongR = swlongc
 
      ! Allocate and create matrices
      if(allocated(CiceR))    deallocate(CiceR)
      if(allocated(longR)) deallocate(longR)
      if(allocated(latR))  deallocate(latR)

      allocate(CiceR(iLongR,iLatR))
      allocate(longR(iLongR),latR(iLatR))

      ! Generate long&lat on each grid point
      do i=1,iLatR
        latR(i) = swlatR+(i-1)*dyR
      enddo
      do i=1,iLongR
        longR(i) = swlongR+(i-1)*dxR
      enddo

      ! Generate interpolation info (south west point and weights)
      do i=1,NP
        adcLat  = RAD2DEG*SFEA(i)
        adcLong = RAD2DEG*SLAM(i)

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
      END SUBROUTINE


      END MODULE OWI_ICE

