c***********************************************************************************
c  PROGRAM  owi22
c
c  owi22 will produce ADCIRC wind file from basin-scale/regional wind/pressure field
c  provided by OWI.  It will require a grid file, two wind field files and two
c  pressure field files. The output file will contain wind velocity in knots and 
c  srface atmospheric pressure in hPa.    sb 8/15/2006
c
      program owi22
      
! Declare variables ----------------------------------------------------

      implicit none
      character(100) :: title
      integer :: nnAtl,nnR,nnB,ne,numAtl,i
      integer :: j,k,xi,yi
      real*8, dimension(:,:), allocatable :: uR,vR,pR,uB,vB,pB,Atlonglat
      real*8, dimension(:), allocatable :: latR,longR,latB,longB
      integer, dimension(:), allocatable :: nAtl
      real*8 :: Along, Alat
      real*8 :: ramp,rampfrac
      real*8 :: uu,vv,PP
      real*8 :: Penv
      integer :: date1R,date2R,date1B,date2B
      integer :: date1w,date2w
      integer :: date1p,date2p
      integer :: iLatR,iLongR,iCYMDHR,iMinR
      real*8  :: dxR,dyR,swlatR,swlongR
      integer :: iLatB,iLongB,iCYMDHB,iMinB
      real*8  :: dxB,dyB,swlatB,swlongB
      integer :: iLatw,iLongw,iCYMDHw,iMinw
      real*8  :: dxw,dyw,swlatw,swlongw
      integer :: iLatp,iLongp,iCYMDHp,iMinp
      real*8  :: dxp,dyp,swlatp,swlongp
      integer :: isnapR,updateR
      integer :: isnapB,updateB

      integer :: regionExists

      integer,allocatable :: swpointsR(:,:)
      integer,allocatable :: swpointsB(:,:)

      real*8 :: w,w1,w2,w3,w4
      real*8,allocatable :: wR(:,:)
      real*8,allocatable :: wB(:,:)

      real*8,allocatable :: umax(:),pmin(:)

      CHARACTER FNAME*60
      logical found
      logical SkipRegion
      integer idummy

! Set a value for Penv

      Penv=1013d0         !mb

! Open Atlantic grid file ----------------------------------------------
      FOUND = .FALSE.
      DO WHILE(.NOT.FOUND)
         WRITE(*,*) 'Enter the name of the ADCIRC UNIT 14 file:'
         READ(*,'(A60)') FNAME
         INQUIRE(FILE=FNAME,EXIST=FOUND)
         IF(.NOT.FOUND) WRITE(*,1010) FNAME
      ENDDO
      WRITE(*,1011) FNAME
1010  FORMAT(' File ',A60,/,' WAS NOT FOUND!  Try again',/)
1011  FORMAT(' File ',A60,/,' WAS FOUND!  Opening & Processing file',/)

      open(AMPI_LUN(1),file=fname)
      call AMPI_LUN_REGISTER(AMPI_LUN(1),file=fname)
      read(AMPI_LUN(1),*) title
      read(AMPI_LUN(1),*) ne, nnAtl
      allocate(Atlonglat(nnAtl,2))
      allocate(nAtl(nnAtl))
      allocate(umax(nnAtl),pmin(nnAtl))
      allocate(swpointsB(nnAtl,2),wB(nnAtl,4))
      allocate(swpointsR(nnAtl,2),wR(nnAtl,4))
      do i=1,nnAtl
         read(AMPI_LUN(1),*) numAtl, Along, Alat
         Atlonglat(i,1)=Along
         Atlonglat(i,2)=Alat
         nAtl(i)=numAtl

         umax(i)=0
         pmin(i)=Penv
      end do


! Create ADCIRC output file --------------------------------------------
      
  41  WRITE(*,*) 'Enter the name of the output unit 22 file:'
      READ(*,'(A60)') FNAME

      open(AMPI_LUN(2),file=fname)
      call AMPI_LUN_REGISTER(AMPI_LUN(2),file=fname)
1000  format(i8,3e13.5)

! Ask if there is region-scale files  ----------------------------------

      DO WHILE(.TRUE.)
         WRITE(*,*) 'Do you have region-scale files?  (YES:1,NO:0)'
         READ(*,*) idummy

         IF(idummy.eq.1) THEN
            SkipRegion = .FALSE.
            EXIT
         ELSE IF(idummy.eq.0) THEN
            SkipRegion = .TRUE.
            EXIT
         ENDIF
      ENDDO


! Oceanweather basin win file -------------------------------------------------

      FOUND = .FALSE.
      DO WHILE(.NOT.FOUND)
         WRITE(*,*) 'Enter the name of the OWI Basin .win file:'
         READ(*,'(A60)') FNAME
         INQUIRE(FILE=FNAME,EXIST=FOUND)
         IF(.NOT.FOUND) WRITE(*,1010) FNAME
      ENDDO
      WRITE(*,1011) FNAME

      open(AMPI_LUN(31),file=fname)
      call AMPI_LUN_REGISTER(AMPI_LUN(31),file=fname)

      ! Read in begining/ending dates of win file
      read(AMPI_LUN(31),10) title,date1w,title,date2w

10    format(a55,i10,a5,i10)


! Oceanweather basin pre file ------------------------------------------------

      FOUND = .FALSE.
      DO WHILE(.NOT.FOUND)
         WRITE(*,*) 'Enter the name of the OWI Basin .pre file:'
         READ(*,'(A60)') FNAME
         INQUIRE(FILE=FNAME,EXIST=FOUND)
         IF(.NOT.FOUND) WRITE(*,1010) FNAME
      ENDDO
      WRITE(*,1011) FNAME

      open(AMPI_LUN(41),file=fname)
      call AMPI_LUN_REGISTER(AMPI_LUN(41),file=fname)

      ! Read in begining/ending dates of win file
      read(AMPI_LUN(41),10) title,date1p,title,date2p

      if(date1w.ne.date1p.or.date2w.ne.date2p) then
         print *, 'HEADER INFO IN WIN AND PRE FILES DO NOT MATCH'
         print *, 'EXECUTION WILL BE TERMINATED'
         stop
      endif

      date1B = date1w
      date2B = date2w


! Oceanweather region win file -----------------------------------------------

      IF(.NOT.SkipRegion) THEN

         FOUND = .FALSE.
         DO WHILE(.NOT.FOUND)
            WRITE(*,*) 'Enter the name of the OWI Region .win file:'
            READ(*,'(A60)') FNAME
            INQUIRE(FILE=FNAME,EXIST=FOUND)
            IF(.NOT.FOUND) WRITE(*,1010) FNAME
         ENDDO
         WRITE(*,1011) FNAME

         open(AMPI_LUN(32),file=fname)
         call AMPI_LUN_REGISTER(AMPI_LUN(32),file=fname)

         ! Read in begining/ending dates of .win file
         read(AMPI_LUN(32),10) title,date1w,title,date2w


! Oceanweather region pre file -----------------------------------------------

         FOUND = .FALSE.
         DO WHILE(.NOT.FOUND)
            WRITE(*,*) 'Enter the name of the OWI Region .pre file:'
            READ(*,'(A60)') FNAME
            INQUIRE(FILE=FNAME,EXIST=FOUND)
            IF(.NOT.FOUND) WRITE(*,1010) FNAME
         ENDDO
         WRITE(*,1011) FNAME

         open(AMPI_LUN(42),file=fname)
         call AMPI_LUN_REGISTER(AMPI_LUN(42),file=fname)

         ! Read in begining/ending dates of .pre file
         read(AMPI_LUN(42),10) title,date1p,title,date2p

         if(date1w.ne.date1p.or.date2w.ne.date2p) then
            print *, 'HEADER INFO IN WIN AND PRE FILES DO NOT MATCH'
            print *, 'EXECUTION WILL BE TERMINATED'
            stop
         endif

         date1R = date1w
         date2R = date2w

      ENDIF ! SkipRegion


! Open first file ------------------------------------------------------------
      isnapB = 0
      isnapR = 0
      updateB = 1
      updateR = 1

      do while(.true.)
!Process Basin Info  ---------------------------------------------------------

         ! Increment counter         
         isnapB = isnapB+1

11       format(t6,i4,t16,i4,t23,f6.0,t32,f6.0,
     &        t44,f8.0,t58,f8.0,t69,i10,i2)

         ! Read Grid Specifications/Date
         read (AMPI_LUN(31),11,end=9999) 
     &       iLatw,iLongw,dxw,dyw,swlatw,swlongw,iCYMDHw,iMinw

         ! Read Grid Specifications/Date
         read (AMPI_LUN(41),11,end=9999) 
     &       iLatp,iLongp,dxp,dyp,swlatp,swlongp,iCYMDHp,iMinp

         if(iLatw.ne.iLatp.or.iLongw.ne.iLongp.or.dxw.ne.dxp.or. 
     &         dyw.ne.dyp.or.swlatw.ne.swlatp.or.swlongw.ne.swlongp.or. 
     &         iCYMDHw.ne.iCYMDHp.or.iMinw.ne.iMinp) then
            print *, 'SNAPSHOT HEADER IN WIN AND PRE FILES DO NOT MATCH'
            print *, 'EXECUTION WILL BE TERMINATED'
            stop
         endif

         ! Check if header info has changed from the previous snapshot
         if(isnapB.gt.1) then
            if(iLatw.ne.iLatB.or.iLongw.ne.iLongB.or.dxw.ne.dxB.or.
     &            dyw.ne.dyB.or.swlatw.ne.swlatB.or.
     &            swlongw.ne.swlongB) then
               updateB = 1
            else
               updateB = 0
            endif
         endif


         iCYMDHB = iCYMDHw
         iMinB = iMinw

         if(updateB.eq.1) then
            print *, 'BASIN  COORDINATE UPDATED'

            iLatB = iLatw
            iLongB = iLongw
            dxB = dxw
            dyB = dyw
            swlatB = swlatw
            swlongB = swlongw
 
            ! Number of grid points
            nnB=iLatB*iLongB


            ! Allocate and create matrices
            if(isnapB.ne.1) then
               deallocate(uB,vB,pB,longB,latB)
            endif

            allocate(uB(iLongB,iLatB),vB(iLongB,iLatB),pB(iLongB,iLatB))
            allocate(longB(iLongB),latB(iLatB))

            ! Generate long&lat on each grid point
            do i=1,iLatB
               latB(i) = swlatB+(i-1)*dyB
            enddo
            do i=1,iLongB
               longB(i) = swlongB+(i-1)*dxB
            enddo

            ! Generate interpolation info (south west point and weights)
            do i=1,nnAtl
               if (Atlonglat(i,1)>=longB(1) .and. 
     &                Atlonglat(i,1)<longB(iLongB) .and. 
     &                Atlonglat(i,2)>=latB(1) .and. 
     &                Atlonglat(i,2)<latB(iLatB)) then
                  do j=1,iLongB-1
                     if (Atlonglat(i,1)>=longB(j) .and.  
     &                      Atlonglat(i,1)<longB(j+1)) then
                        xi=j
                        exit
                     endif
                  enddo
                  
                  do k=1,iLatB-1
                     if (Atlonglat(i,2)>=latB(k) .and. 
     &                      Atlonglat(i,2)<latB(k+1)) then
                        yi=k
                        exit
                     endif
                  enddo

                  swpointsB(i,1) = xi
                  swpointsB(i,2) = yi

                  w=(longB(xi+1)-longB(xi))*(latB(yi+1)-latB(yi))
                  w1=(longB(xi+1)-Atlonglat(i,1))*
     &                 (latB(yi+1)-Atlonglat(i,2))
                  w2=(Atlonglat(i,1)-longB(xi))*
     &                 (latB(yi+1)-Atlonglat(i,2))
                  w3=(Atlonglat(i,1)-longB(xi))*
     &                 (Atlonglat(i,2)-latB(yi))
                  w4=(longB(xi+1)-Atlonglat(i,1))*
     &                 (Atlonglat(i,2)-latB(yi))

                  wB(i,1)=w1/w
                  wB(i,2)=w2/w
                  wB(i,3)=w3/w
                  wB(i,4)=w4/w
               else
                  swpointsB(i,1) = 0
                  swpointsB(i,2) = 0

               endif
            enddo
         endif

      
         ! Read u/v components of the wind
12       format(8f10.0)
         read(AMPI_LUN(31),12) ((uB(i,j),i=1,iLongB),j=1,iLatB)
         read(AMPI_LUN(31),12) ((vB(i,j),i=1,iLongB),j=1,iLatB)

         uB(:,:) = 1.943844492d0*uB(:,:)
         vB(:,:) = 1.943844492d0*vB(:,:)

         ! Read pressure
         read(AMPI_LUN(41),12) ((pB(i,j),i=1,iLongB),j=1,iLatB)



!Process Region Info  --------------------------------------------------------
         
         regionExists = 1

!         print *, iCYMDHB, date1R
!         print *, iCYMDHB, date2R
         if(SkipRegion) regionExists = 0
         if(iCYMDHB.lt.date1R) regionExists = 0
         if(iCYMDHB.eq.date2R.and.iMinR.ne.0) regionExists = 0
         if(iCYMDHB.gt.date2R) regionExists = 0

         if (regionExists .eq. 1) then
            ! Increment counter         
            isnapR = isnapR+1

            ! Read Grid Specifications/Date
            read (AMPI_LUN(32),11,end=9999) 
     &       iLatw,iLongw,dxw,dyw,swlatw,swlongw,iCYMDHw,iMinw

            ! Read Grid Specifications/Date
            read (AMPI_LUN(42),11,end=9999) 
     &       iLatp,iLongp,dxp,dyp,swlatp,swlongp,iCYMDHp,iMinp

            if(iLatw.ne.iLatp.or.iLongw.ne.iLongp.or.dxw.ne.dxp.or.
     &         dyw.ne.dyp.or.swlatw.ne.swlatp.or.swlongw.ne.swlongp.or.
     &         iCYMDHw.ne.iCYMDHp.or.iMinw.ne.iMinp) then
               print *, 'SNAPSHOT HEADER IN WIN AND PRE FILES DO NOT MATCH'
               print *, 'EXECUTION WILL BE TERMINATED'
               stop

            endif

            ! Check if header info has changed from the previous snapshot
            if(isnapR.gt.1) then
               if(iLatw.ne.iLatR.or.iLongw.ne.iLongR.or.dxw.ne.dxR.or. 
     &            dyw.ne.dyR.or.swlatw.ne.swlatR.or.
     &            swlongw.ne.swlongR) then
                  updateR = 1
               else
                  updateR = 0
               endif

            endif


            iCYMDHR = iCYMDHw
            iMinR = iMinw

            if(iCYMDHB.ne.iCYMDHR.or.iMinB.ne.iMinR) then
               print *, 'SNAPSHOTS NOT SYNCRONIZED'
               print *, '  iCYMDHB=',iCYMDHB, ' iMinB=',iMinB
               print *, '  iCYMDHR=',iCYMDHR, ' iMinR=',iMinR
               print *, 'EXECUTION WILL BE TERMINATED'
               stop
            
            endif

            if(updateR.eq.1) then
               print *, 'REGION COORDINATE UPDATED'

               iLatR = iLatw
               iLongR = iLongw
               dxR = dxw
               dyR = dyw
               swlatR = swlatw
               swlongR = swlongw
 
               ! Number of grid points
               nnR=iLatR*iLongR


               ! Allocate and create matrices
               if(isnapR.ne.1) then
                  deallocate(uR,vR,pR,longR,latR,swpointsR,wR)
               endif

               allocate(uR(iLongR,iLatR),vR(iLongR,iLatR),pR(iLongR,iLatR))
               allocate(longR(iLongR),latR(iLatR))

               ! Generate long&lat on each grid point
               do i=1,iLatR
                  latR(i) = swlatR+(i-1)*dyR
               enddo
               do i=1,iLongR
                  longR(i) = swlongR+(i-1)*dxR
               enddo

               ! Generate interpolation info (south west point and weights)
               do i=1,nnAtl
                  if (Atlonglat(i,1)>=longR(1) .and. 
     &                Atlonglat(i,1)<longR(iLongR) .and. 
     &                Atlonglat(i,2)>=latR(1) .and. 
     &                Atlonglat(i,2)<latR(iLatR)) then
                     do j=1,iLongR-1
                        if (Atlonglat(i,1)>=longR(j) .and. 
     &                      Atlonglat(i,1)<longR(j+1)) then
                           xi=j
                           exit
                        endif
                     enddo
                  
                     do k=1,iLatR-1
                        if (Atlonglat(i,2)>=latR(k) .and. 
     &                      Atlonglat(i,2)<latR(k+1)) then
                           yi=k
                           exit
                        endif
                     enddo

                     swpointsR(i,1) = xi
                     swpointsR(i,2) = yi

                     w=(longR(xi+1)-longR(xi))*
     &                 (latR(yi+1)-latR(yi))
                     w1=(longR(xi+1)-Atlonglat(i,1))*
     &                 (latR(yi+1)-Atlonglat(i,2))
                     w2=(Atlonglat(i,1)-longR(xi))*
     &                 (latR(yi+1)-Atlonglat(i,2))
                     w3=(Atlonglat(i,1)-longR(xi))*
     &                 (Atlonglat(i,2)-latR(yi))
                     w4=(longR(xi+1)-Atlonglat(i,1))*
     &                 (Atlonglat(i,2)-latR(yi))

                     wR(i,1)=w1/w
                     wR(i,2)=w2/w
                     wR(i,3)=w3/w
                     wR(i,4)=w4/w
                  else
                     swpointsR(i,1) = 0
                     swpointsR(i,2) = 0

                  endif
               enddo
            endif

      
            ! Read u/v components of the wind
            read(AMPI_LUN(32),12) ((uR(i,j),i=1,iLongR),j=1,iLatR) 
            read(AMPI_LUN(32),12) ((vR(i,j),i=1,iLongR),j=1,iLatR)

            uR(:,:) = 1.943844492d0*uR(:,:)
            vR(:,:) = 1.943844492d0*vR(:,:)

            ! Read pressure
            read(AMPI_LUN(42),12) ((pR(i,j),i=1,iLongR),j=1,iLatR)

         endif ! regionExists .eq. 1
         

!Interpolate onto ADCIRC grid and write to file -------------------------
      

         rampfrac = isnapB-1
         if (rampfrac<36) then
            ramp = tanh(18d0*rampfrac/36d0)
         end if
         ramp = 1.0

         write(AMPI_LUN(2),*) '# ', iCYMDHB, iMinB

         if(regionExists.eq.0) then
30          format('PROCESSING BASIN DATA',i12,' ',i2,' RAMPING=',F5.2)
            print 30, iCYMDHB,iMinB,ramp
         else
33          format('PROCESSING BASIN&REGION DATA',i12,' ',
     &          i2,' RAMPING=',F5.2)
            print 33, iCYMDHB,iMinB,ramp
         endif

         do i=1,nnAtl

            uu=-9999.9D0

            ! BASIN ---------------------------------------------------------

            if (swpointsB(i,1).gt.0) then
               xi = swpointsB(i,1)
               yi = swpointsB(i,2)

               w1=wB(i,1)
               w2=wB(i,2)
               w3=wB(i,3)
               w4=wB(i,4)

               uu=w1*uB(xi,yi)+w2*uB(xi+1,yi)+w3*
     &              uB(xi+1,yi+1)+w4*uB(xi,yi+1)
               vv=w1*vB(xi,yi)+w2*vB(xi+1,yi)+w3*
     &              vB(xi+1,yi+1)+w4*vB(xi,yi+1)
               PP=w1*pB(xi,yi)+w2*pB(xi+1,yi)+w3*
     &              pB(xi+1,yi+1)+w4*pB(xi,yi+1)

            endif

            ! REGION ---------------------------------------------------------
            ! uu, vv and PP will be overwritten if region data exist.

            if (regionExists.eq.1.and.swpointsR(i,1).gt.0) then
               xi = swpointsR(i,1)
               yi = swpointsR(i,2)

               w1=wR(i,1)
               w2=wR(i,2)
               w3=wR(i,3)
               w4=wR(i,4)

               uu=w1*uR(xi,yi)+w2*uR(xi+1,yi)+
     &              w3*uR(xi+1,yi+1)+w4*uR(xi,yi+1)
               vv=w1*vR(xi,yi)+w2*vR(xi+1,yi)+
     &              w3*vR(xi+1,yi+1)+w4*vR(xi,yi+1)
               PP=w1*pR(xi,yi)+w2*pR(xi+1,yi)+
     &              w3*pR(xi+1,yi+1)+w4*pR(xi,yi+1)

            endif

            ! OUTPUT ----------------------------------------------------------

            if(uu.ne.-9999.9D0) then
               if (rampfrac<36) then
                  uu=uu*ramp
                  vv=vv*ramp
                  PP=Penv-(Penv-PP)*ramp
               endif

               ! 08/14/2006 sb
               ! Conversion from OWI 30-min avg winds to 10-min avg winds
               ! 1.1076 = 1.04 * 1.065
               uu = uu * 1.1076
               vv = vv * 1.1076

               write(AMPI_LUN(2),1000) nAtl(i),uu,vv,PP
            
               if(umax(i).lt.sqrt(uu*uu+vv*vv)) then
                  umax(i)=sqrt(uu*uu+vv*vv)
               endif

               if(pmin(i).gt.PP) then
                  pmin(i)=PP
               endif

            end if

         enddo

      enddo

9999  continue

      print *, isnapB-1, 'BASIN  SNAPSHOTS WERE READ'
      print *, isnapR-1, 'REGION SNAPSHOTS WERE READ'
      print *, 'DONE'

      ! Close files
      close(AMPI_LUN(1))
      close(AMPI_LUN(2))
      close(AMPI_LUN(31))
      close(AMPI_LUN(41))
      IF(.not.SkipRegion) THEN
        close(AMPI_LUN(32))
        close(AMPI_LUN(42))
      ENDIF

      ! Write peak pressure velocity
      print *, 'Writing peak pressure velocity in peak73.63 ....'
      open(AMPI_LUN(10),file="peak73.63")
      call AMPI_LUN_REGISTER(AMPI_LUN(10),file="peak73.63")
      write(AMPI_LUN(10),1110) title
      WRITE(AMPI_LUN(10),3645) 1,nnAtl,1000,1000,1
      WRITE(AMPI_LUN(10),2120) 1000,1000
      do i=1,nnAtl
         write(AMPI_LUN(10),3330) nAtl(i),pmin(i)
      end do
      close(AMPI_LUN(10))

      ! Write peak wind velocity
      print *, 'Writing peak wind velocity in peak74.63 ....'
      open(AMPI_LUN(10),file="peak74.63")
      call AMPI_LUN_REGISTER(AMPI_LUN(10),file="peak74.63")
      write(AMPI_LUN(10),1110) title
      WRITE(AMPI_LUN(10),3645) 1,nnAtl,1000,1000,1
      WRITE(AMPI_LUN(10),2120) 1000,1000
      do i=1,nnAtl
         write(AMPI_LUN(10),3330) nAtl(i),umax(i)
      end do
      close(AMPI_LUN(10))
1110  format(a80) 
3645  FORMAT(1X,I10,1X,I10,1X,E15.7,1X,I5,1X,I5)
2120  FORMAT(2X,E20.10,5X,I10)
3330  format(2x,i8,2x,e15.8)

      ! Deallocate memory
      deallocate(Atlonglat)
      deallocate(nAtl)
      deallocate(umax,pmin)
      deallocate(uB,vB,pB,longB,latB,swpointsB,wB)
      if(.not.SkipRegion) THEN
        deallocate(uR,vR,pR,longR,latR,swpointsR,wR)
      endif

      end program owi22
