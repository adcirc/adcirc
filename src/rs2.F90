!******************************************************************************
!  MODULE RS2
!    Written by s.b. 09/xx/2006
!------------------------------------------------------------------------------
! This module reads in UNIT 23 file and external files specified in UNIT 23.
! This module was made to have ADCIRC read in a single global UNIT 23 file
! and sets of global STWAVE output files (.sim and .rad files) from each
! processor.
!------------------------------------------------------------------------------
! sb  12/11/2006  RS2GET now takes maximum values if stwave domains overlap.
!                 (It used to take the value of the domain of the first
!                 appearance in a fort.23 file.)
! vjp 10/21/2006  modified logic for type2
!******************************************************************************

    MODULE RS2
    USE SIZES

    USE GLOBAL,ONLY : NSCREEN, ScreenUnit
#ifdef CMPI
    use messenger, only : msg_fini
#endif
    IMPLICIT NONE

! ariables for input file
    integer :: ns,type
    character(len=100),allocatable,dimension(:) :: simfiles ! names of .sim files
    character(len=100),allocatable,dimension(:) :: radfiles ! names of .rad files
    integer,allocatable,dimension(:)            :: fullplanes ! 0: half plane, 1: full plane
    integer,allocatable,dimension(:)            :: nbs ! number of blank snaps
    integer,allocatable,dimension(:)            :: rfids ! id of .rad files
    integer,allocatable,dimension(:) :: nis,njs
    integer,save,allocatable,dimension(:) :: endoffile

    integer,allocatable,dimension(:,:,:) :: rsij
    real(sz),allocatable,dimension(:,:,:) :: rscnf

    real(sz),allocatable,dimension(:) :: dxincs,dyincs
    real(sz),allocatable,dimension(:) :: x0s,y0s,azimuths

    integer :: numSets,numBlankSnaps,cntSnaps,numSkipSnaps
    real(SZ) :: windMultiplier


    PUBLIC

!---------------------end of data declarations--------------------------------C


    CONTAINS


!***********************************************************************
!   SOBROUTINE RS2INIT
!***********************************************************************

    subroutine rs2init(rsnx,rsny,np)
    use sizes
    implicit none
    integer,intent(in) :: np
    real(sz),intent(out) :: rsnx(:),rsny(:)

    if(myproc == 0) then
        write(screenunit,*)
        write(screenunit,*) 'INIT RADIATION STRESS ARRAYS (NRS=2)...'
    endif
    write(16,*)
    write(16,*) 'INITIALIZING RADIATION STRESS ARRAYS (NRS=2)...'

! ead meta info


    OPEN(23,FILE=TRIM(GBLINPUTDIR)//'/'//'fort.23',STATUS='OLD')

! type of fort.23 file, 1: include wave radiation
!                       2: does not include wave radiation
    read(23,*,err=99999) type

    if(type == 1) then
        stop 'RS2INIT: SORRY! TYPE-1 RS FILE NOT YET SUPPORTED.'
    !        call rs2init_type1(rsnx,rsny,np)
    else if(type == 2) then
        call rs2init_type2(rsnx,rsny,np)
    else
        stop 'RS2 TYPE IS NEITHER 1 NOR 2.  TERMINATED.'
    endif

    return

    99999 CONTINUE

#ifdef CMPI
    call msg_fini()
#endif
    STOP 'RADIATION STRESS READ ERROR (RS2-0)'

    END SUBROUTINE


!$$$C***********************************************************************
!$$$C   SOBROUTINE RS2INIT_TYPE1
!$$$C***********************************************************************
!$$$
!$$$      subroutine rs2init_type1(rsnx,rsny,np)
!$$$      use sizes
!$$$      use global,only : slam,sfea,rad2deg
!$$$      implicit none
!$$$      integer,intent(in) :: np
!$$$      real(sz),intent(out) :: rsnx(:),rsny(:)
!$$$      integer ::  c,i,j,k,n,p,s,ni,nj
!$$$      real(sz) :: lon,lat
!$$$      real(sz) :: dxinc,dyinc
!$$$      real(sz) :: stwlon1,stwlat1,stwlon2,stwlat2
!$$$      real(sz) :: stwlon3,stwlat3,stwlon4,stwlat4
!$$$      real(sz) :: x1,y1,x2,y2,x3,y3
!$$$      real(sz) :: subarea1,subarea2,subarea3,subarea4,totalarea
!$$$      real(sz),dimension(:,:,:), allocatable :: stwlonlat
!$$$      real(sz) :: stwlonmin,stwlonmax,stwlatmin,stwlatmax
!$$$      integer,parameter :: nSearchBins = 50
!$$$      integer,allocatable :: nCellsInSearchBins(:)
!$$$      integer,allocatable :: searchBins(:,:,:)
!$$$      real(sz),allocatable :: searchBinPartitions(:)
!$$$      real(sz) :: p1,p2,cmax,cmin,tol
!$$$
!$$$
!$$$      stop 'RS2INIT_TYPE1: SORRY! TYPE-1 RS FILE NOT YET SUPPORTED.'
!$$$
!$$$      read(23,*,err=99999) ns    !number of stwave grid files
!$$$
!$$$      if(ns.le.0) then
!$$$        write(screenunit,1004)
!$$$        write(16,1004)
!$$$#ifdef CMPI
!$$$        call msg_fini()
!$$$#endif
!$$$        stop
!$$$      endif
!$$$
!$$$      allocate(fullplanes(ns),rfids(ns))
!$$$      allocate(nis(ns),njs(ns),dxincs(ns),dyincs(ns))
!$$$      allocate(x0s(ns),y0s(ns),azimuths(ns))
!$$$      allocate(endoffile(ns))
!$$$
!$$$      allocate(rsij(4,np),rscnf(3,np))
!$$$
!$$$      do n=1,np
!$$$        rsij(1,n) = 0
!$$$      enddo
!$$$
!$$$      !compute mapping coefficients
!$$$      do s=1,ns
!$$$        read(23,*,err=99999) ! skip a line
!$$$        read(23,*,err=99999) fullplanes(s)
!$$$
!$$$        if(fullplanes(s).eq.1) then
!$$$          read(23,*,err=99999) nis(s), njs(s), dxincs(s), dyincs(s)
!$$$        else
!$$$          read(23,*,err=99999) nis(s), njs(s), dxincs(s)
!$$$          dyincs(s) = dxincs(s)
!$$$        endif
!$$$
!$$$        ni = nis(s)
!$$$        nj = njs(s)
!$$$
!$$$        allocate(stwlonlat(2,ni,nj))
!$$$        do j=1,nj
!$$$          read(23,*) ((stwlonlat(k,i,j), k=1,2), i=1,ni)
!$$$        enddo
!$$$
!$$$        stwlonmin =  1d10
!$$$        stwlonmax = -1d10
!$$$        stwlatmin =  1d10
!$$$        stwlatmax = -1d10
!$$$        do j=1,nj
!$$$          do i=1,ni
!$$$            if(stwlonlat(1,i,j).eq.0.d0) cycle
!$$$            if(stwlonmin.gt.stwlonlat(1,i,j)) stwlonmin=stwlonlat(1,i,j)
!$$$            if(stwlonmax.lt.stwlonlat(1,i,j)) stwlonmax=stwlonlat(1,i,j)
!$$$            if(stwlatmin.gt.stwlonlat(2,i,j)) stwlatmin=stwlonlat(2,i,j)
!$$$            if(stwlatmax.lt.stwlonlat(2,i,j)) stwlatmax=stwlonlat(2,i,j)
!$$$          enddo
!$$$        enddo
!$$$
!$$$        !prepare search bins
!$$$        allocate(nCellsInSearchBins(nSearchBins))
!$$$        allocate(searchBins(2,(ni-1)*(nj-1),nSearchBins))
!$$$        allocate(searchBinPartitions(nSearchBins+1))
!$$$        searchBinPartitions(1) = stwlonmin
!$$$        if(myproc.eq.0) then
!$$$          write(screenunit,*) 'PREPARING SEARCH BINS FOR GRID ',s,'...'
!$$$        endif
!$$$        tol = abs(stwlonmax-stwlonmin)/real(nSearchBins)*0.01d0
!$$$        do p=1,nSearchBins
!$$$          searchBinPartitions(p+1) =
!$$$     &         stwlonmin + (stwlonmax-stwlonmin)*
!$$$     &         real(p)/real(nSearchBins)
!$$$
!$$$          p1 = searchBinPartitions(p)
!$$$          p2 = searchBinPartitions(p+1)
!$$$
!$$$          p1 = p1 - tol
!$$$          p2 = p2 + tol
!$$$
!$$$          nCellsInSearchBins(p) = 0
!$$$
!$$$          do j=1,nj-1
!$$$            do i=1,ni-1
!$$$              stwlon1 = stwlonlat(1,i,j)
!$$$              stwlon2 = stwlonlat(1,i+1,j)
!$$$              stwlon3 = stwlonlat(1,i+1,j+1)
!$$$              stwlon4 = stwlonlat(1,i,j+1)
!$$$
!$$$              cmax = stwlon1
!$$$              if(cmax.lt.stwlon2) cmax = stwlon2
!$$$              if(cmax.lt.stwlon3) cmax = stwlon3
!$$$              if(cmax.lt.stwlon4) cmax = stwlon4
!$$$
!$$$              cmin = stwlon1
!$$$              if(cmin.gt.stwlon2) cmin = stwlon2
!$$$              if(cmin.gt.stwlon3) cmin = stwlon3
!$$$              if(cmin.gt.stwlon4) cmin = stwlon4
!$$$
!$$$              if(cmax.ge.p1.and.cmin.le.p2) then
!$$$                 nCellsInSearchBins(p) = nCellsInSearchBins(p) + 1
!$$$                 searchBins(1,nCellsInSearchBins(p),p) = i
!$$$                 searchBins(2,nCellsInSearchBins(p),p) = j
!$$$              endif
!$$$            enddo
!$$$          enddo
!$$$        enddo
!$$$        if(myproc.eq.0) then
!$$$          write(screenunit,*) ' SEARCH BINS ARE READY'
!$$$        endif
!$$$
!$$$        do n=1,np
!$$$
!$$$          if(rsij(1,n).ne.0) cycle
!$$$
!$$$          lat = RAD2DEG*SFEA(n)
!$$$          lon = RAD2DEG*SLAM(n)
!$$$
!$$$          if(lon.lt.stwlonmin.or.lon.gt.stwlonmax.or.
!$$$     &       lat.lt.stwlatmin.or.lat.gt.stwlatmax) cycle
!$$$
!$$$
!$$$          do p=1,nSearchBins
!$$$            p1 = searchBinPartitions(p)
!$$$            p2 = searchBinPartitions(p+1)
!$$$
!$$$            if(lon.ge.p1.and.lon.le.p2) exit
!$$$          enddo
!$$$
!$$$          if(p.gt.nSearchBins) cycle
!$$$
!$$$          cellloop: do c=1,nCellsInSearchBins(p)
!$$$            i = searchBins(1,c,p)
!$$$            j = searchBins(2,c,p)
!$$$            stwlon1 = stwlonlat(1,i,j)
!$$$            stwlat1 = stwlonlat(2,i,j)
!$$$            stwlon2 = stwlonlat(1,i+1,j)
!$$$            stwlat2 = stwlonlat(2,i+1,j)
!$$$            stwlon3 = stwlonlat(1,i+1,j+1)
!$$$            stwlat3 = stwlonlat(2,i+1,j+1)
!$$$            stwlon4 = stwlonlat(1,i,j+1)
!$$$            stwlat4 = stwlonlat(2,i,j+1)
!$$$
!$$$            if(stwlon1.eq.0.d0.or.stwlon2.eq.0.d0.or.
!$$$     &           stwlon3.eq.0.d0.or.stwlon4.eq.0.d0) then
!$$$              cycle
!$$$            endif
!$$$
!$$$            !triangle 1 (nodes 1, 2 and 3)
!$$$            x1 = lon
!$$$            y1 = lat
!$$$            x2 = stwlon2
!$$$            y2 = stwlat2
!$$$            x3 = stwlon3
!$$$            y3 = stwlat3
!$$$            subarea1 = abs((x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1))
!$$$
!$$$            x1 = lon
!$$$            y1 = lat
!$$$            x2 = stwlon3
!$$$            y2 = stwlat3
!$$$            x3 = stwlon1
!$$$            y3 = stwlat1
!$$$            subarea2 = abs((x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1))
!$$$
!$$$            x1 = lon
!$$$            y1 = lat
!$$$            x2 = stwlon1
!$$$            y2 = stwlat1
!$$$            x3 = stwlon2
!$$$            y3 = stwlat2
!$$$            subarea3 = abs((x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1))
!$$$
!$$$            x1 = stwlon1
!$$$            y1 = stwlat1
!$$$            x2 = stwlon2
!$$$            y2 = stwlat2
!$$$            x3 = stwlon3
!$$$            y3 = stwlat3
!$$$            totalarea = abs((x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1))
!$$$
!$$$            if((subarea1+subarea2+subarea3).le.
!$$$     &           (totalarea*1.00001d0)) then
!$$$              rsij(1,n) = s
!$$$              rsij(2,n) = 1
!$$$              rsij(3,n) = i
!$$$              rsij(4,n) = j
!$$$              rscnf(1,n)=subarea1/totalarea
!$$$              rscnf(2,n)=subarea2/totalarea
!$$$              rscnf(3,n)=subarea3/totalarea
!$$$              exit cellloop
!$$$            endif
!$$$
!$$$            !triangle 2 (nodes 3, 4 and 1)
!$$$            x1 = lon
!$$$            y1 = lat
!$$$            x2 = stwlon4
!$$$            y2 = stwlat4
!$$$            x3 = stwlon1
!$$$            y3 = stwlat1
!$$$            subarea1 = abs((x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1))
!$$$
!$$$            x1 = lon
!$$$            y1 = lat
!$$$            x2 = stwlon1
!$$$            y2 = stwlat1
!$$$            x3 = stwlon3
!$$$            y3 = stwlat3
!$$$            subarea2 = abs((x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1))
!$$$
!$$$            x1 = lon
!$$$            y1 = lat
!$$$            x2 = stwlon3
!$$$            y2 = stwlat3
!$$$            x3 = stwlon4
!$$$            y3 = stwlat4
!$$$            subarea3 = abs((x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1))
!$$$
!$$$            x1 = stwlon3
!$$$            y1 = stwlat3
!$$$            x2 = stwlon4
!$$$            y2 = stwlat4
!$$$            x3 = stwlon1
!$$$            y3 = stwlat1
!$$$            totalarea = abs((x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1))
!$$$
!$$$            if((subarea1+subarea2+subarea3).le.
!$$$     &           (totalarea*1.00001d0)) then
!$$$              rsij(1,n) = s
!$$$              rsij(2,n) = 2
!$$$              rsij(3,n) = i
!$$$              rsij(4,n) = j
!$$$              rscnf(1,n)=subarea1/totalarea
!$$$              rscnf(2,n)=subarea2/totalarea
!$$$              rscnf(3,n)=subarea3/totalarea
!$$$              exit cellloop
!$$$            endif
!$$$          enddo cellloop
!$$$        enddo
!$$$
!$$$        deallocate(stwlonlat)
!$$$        deallocate(nCellsInSearchBins)
!$$$        deallocate(searchBins)
!$$$        deallocate(searchBinPartitions)
!$$$
!$$$      enddo
!$$$
!$$$      !initialize some variables
!$$$      cntSnaps = 0
!$$$      endoffile(:) = 0
!$$$
!$$$      if(myproc.eq.0) then
!$$$        write(screenunit,*)
!$$$        write(screenunit,*) 'INIT OF RADIATION STRESS ARRAYS IS DONE.'
!$$$      endif
!$$$      write(16,*)
!$$$      write(16,*) 'INIT OF RADIATION STRESS ARRAYS IS DONE.'
!$$$
!$$$      return
!$$$
!$$$ 1004 FORMAT(//,1X,' NUMBER OF SETS WAS SPECIFIED'//
!$$$     &       'INCORRECTLY IN UNIT23.'/
!$$$     &       ' IT MUST BE GREATER THAN 0.'//)
!$$$
!$$$99999 CONTINUE
!$$$
!$$$#ifdef CMPI
!$$$      call msg_fini()
!$$$#endif
!$$$      STOP 'RADIATION STRESS READ ERROR (RS2-1)'
!$$$
!$$$      END SUBROUTINE


!***********************************************************************
!   SOBROUTINE RS2INIT_TYPE2
!***********************************************************************

    subroutine rs2init_type2(rsnx,rsny,np)
    use sizes
    use global,only : rad2deg
    use mesh, only : slam, sfea
    implicit none
    integer,intent(in) :: np
    real(sz),intent(out) :: rsnx(:),rsny(:)
    integer ::  c,i,j,k,n,p,s,ni,nj
    real(sz) :: lon,lat
    real(sz) :: dxinc,dyinc
    real(sz) :: stwlon1,stwlat1,stwlon2,stwlat2
    real(sz) :: stwlon3,stwlat3,stwlon4,stwlat4
    real(sz) :: x1,y1,x2,y2,x3,y3
    real(sz) :: subarea1,subarea2,subarea3,subarea4,totalarea
    real(sz),dimension(:,:,:), allocatable :: stwlonlat
    real(sz) :: stwlonmin,stwlonmax,stwlatmin,stwlatmax
    integer,parameter :: nSearchBins = 50
    integer,allocatable :: nCellsInSearchBins(:)
    integer,allocatable :: searchBins(:,:,:)
    real(sz),allocatable :: searchBinPartitions(:)
    real(sz) :: p1,p2,cmax,cmin,tol


    read(23,*,err=99999) ns    !number of stwave grid files

    if(ns <= 0) then
        write(screenunit,1004)
        write(16,1004)
#ifdef CMPI
        call msg_fini()
#endif
        stop
    endif

    allocate(fullplanes(ns),rfids(ns))
    allocate(nis(ns),njs(ns),dxincs(ns),dyincs(ns))
    allocate(x0s(ns),y0s(ns),azimuths(ns))
    allocate(endoffile(ns))
    allocate(simfiles(ns),radfiles(ns))
    allocate(nbs(ns))

    allocate(rsij(3,np,ns),rscnf(3,np,ns))

    rsij(1,:,:) = 0
             

! ompute mapping coefficients
    do s=1,ns
        read(23,*,err=99999) ! skip a line

        read(23,'(A)',err=99999) simfiles(s)
        simfiles(s) = adjustl(simfiles(s))
        read(23,'(A)',err=99999) radfiles(s)
        radfiles(s) = adjustl(radfiles(s))
        read(23,*,err=99999) nbs(s)
        if (myproc == 0) then
            write(screenunit,*) " NUMBER BLANK RS SNAPS FOR REGION ",s, &
            " = ", nbs(s)
        endif
        read(23,*,err=99999) x0s(s),y0s(s),azimuths(s)

        read(23,*,err=99999) fullplanes(s)

        if(fullplanes(s) == 1) then
            read(23,*,err=99999) nis(s), njs(s), dxincs(s), dyincs(s)
        else
            read(23,*,err=99999) nis(s), njs(s), dxincs(s)
            dyincs(s) = dxincs(s)
        endif
                
        ni = nis(s)
        nj = njs(s)

        allocate(stwlonlat(2,ni,nj))
        do j=1,nj
            read(23,*) ((stwlonlat(k,i,j), k=1,2), i=1,ni)
        enddo
                
        stwlonmin =  1d10
        stwlonmax = -1d10
        stwlatmin =  1d10
        stwlatmax = -1d10
        do j=1,nj
            do i=1,ni
                if(stwlonlat(1,i,j) == 0.d0) cycle
                if(stwlonmin > stwlonlat(1,i,j)) stwlonmin=stwlonlat(1,i,j)
                if(stwlonmax < stwlonlat(1,i,j)) stwlonmax=stwlonlat(1,i,j)
                if(stwlatmin > stwlonlat(2,i,j)) stwlatmin=stwlonlat(2,i,j)
                if(stwlatmax < stwlonlat(2,i,j)) stwlatmax=stwlonlat(2,i,j)
            enddo
        enddo

    ! repare search bins
        allocate(nCellsInSearchBins(nSearchBins))
        allocate(searchBins(2,(ni-1)*(nj-1),nSearchBins))
        allocate(searchBinPartitions(nSearchBins+1))
        searchBinPartitions(1) = stwlonmin
        if(myproc == 0) then
            write(screenunit,*) 'PREPARING SEARCH BINS FOR GRID ',s,'...'
        endif
        tol = abs(stwlonmax-stwlonmin)/real(nSearchBins)*0.01d0
        do p=1,nSearchBins
            searchBinPartitions(p+1) = &
            stwlonmin + (stwlonmax-stwlonmin)* &
            real(p)/real(nSearchBins)

            p1 = searchBinPartitions(p)
            p2 = searchBinPartitions(p+1)

            p1 = p1 - tol
            p2 = p2 + tol

            nCellsInSearchBins(p) = 0

            do j=1,nj-1
                do i=1,ni-1
                    stwlon1 = stwlonlat(1,i,j)
                    stwlon2 = stwlonlat(1,i+1,j)
                    stwlon3 = stwlonlat(1,i+1,j+1)
                    stwlon4 = stwlonlat(1,i,j+1)

                    cmax = stwlon1
                    if(cmax < stwlon2) cmax = stwlon2
                    if(cmax < stwlon3) cmax = stwlon3
                    if(cmax < stwlon4) cmax = stwlon4

                    cmin = stwlon1
                    if(cmin > stwlon2) cmin = stwlon2
                    if(cmin > stwlon3) cmin = stwlon3
                    if(cmin > stwlon4) cmin = stwlon4

                    if(cmax >= p1 .AND. cmin <= p2) then
                        nCellsInSearchBins(p) = nCellsInSearchBins(p) + 1
                        searchBins(1,nCellsInSearchBins(p),p) = i
                        searchBins(2,nCellsInSearchBins(p),p) = j
                    endif
                enddo
            enddo
        enddo
        if(myproc == 0) then
            write(screenunit,*) ' SEARCH BINS ARE READY'
        endif

        do n=1,np

            lat = RAD2DEG*SFEA(n)
            lon = RAD2DEG*SLAM(n)

            if(lon < stwlonmin .OR. lon > stwlonmax .OR. &
            lat < stwlatmin .OR. lat > stwlatmax) cycle
                        

            do p=1,nSearchBins
                p1 = searchBinPartitions(p)
                p2 = searchBinPartitions(p+1)
                            
                if(lon >= p1 .AND. lon <= p2) exit
            enddo

            if(p > nSearchBins) cycle

            cellloop: do c=1,nCellsInSearchBins(p)
                i = searchBins(1,c,p)
                j = searchBins(2,c,p)
                stwlon1 = stwlonlat(1,i,j)
                stwlat1 = stwlonlat(2,i,j)
                stwlon2 = stwlonlat(1,i+1,j)
                stwlat2 = stwlonlat(2,i+1,j)
                stwlon3 = stwlonlat(1,i+1,j+1)
                stwlat3 = stwlonlat(2,i+1,j+1)
                stwlon4 = stwlonlat(1,i,j+1)
                stwlat4 = stwlonlat(2,i,j+1)

                if(stwlon1 == 0.d0 .OR. stwlon2 == 0.d0 .OR. &
                stwlon3 == 0.d0 .OR. stwlon4 == 0.d0) then
                    cycle
                endif

            ! riangle 1 (nodes 1, 2 and 3)
                x1 = lon
                y1 = lat
                x2 = stwlon2
                y2 = stwlat2
                x3 = stwlon3
                y3 = stwlat3
                subarea1 = abs((x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1))
                            
                x1 = lon
                y1 = lat
                x2 = stwlon3
                y2 = stwlat3
                x3 = stwlon1
                y3 = stwlat1
                subarea2 = abs((x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1))
                            
                x1 = lon
                y1 = lat
                x2 = stwlon1
                y2 = stwlat1
                x3 = stwlon2
                y3 = stwlat2
                subarea3 = abs((x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1))
                            
                x1 = stwlon1
                y1 = stwlat1
                x2 = stwlon2
                y2 = stwlat2
                x3 = stwlon3
                y3 = stwlat3
                totalarea = abs((x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1))
                            
                if((subarea1+subarea2+subarea3) <= &
                (totalarea*1.00001d0)) then
                    rsij(1,n,s) = 1
                    rsij(2,n,s) = i
                    rsij(3,n,s) = j
                    rscnf(1,n,s)=subarea1/totalarea
                    rscnf(2,n,s)=subarea2/totalarea
                    rscnf(3,n,s)=subarea3/totalarea
                    exit cellloop
                endif

            ! riangle 2 (nodes 3, 4 and 1)
                x1 = lon
                y1 = lat
                x2 = stwlon4
                y2 = stwlat4
                x3 = stwlon1
                y3 = stwlat1
                subarea1 = abs((x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1))
                            
                x1 = lon
                y1 = lat
                x2 = stwlon1
                y2 = stwlat1
                x3 = stwlon3
                y3 = stwlat3
                subarea2 = abs((x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1))
                            
                x1 = lon
                y1 = lat
                x2 = stwlon3
                y2 = stwlat3
                x3 = stwlon4
                y3 = stwlat4
                subarea3 = abs((x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1))
                            
                x1 = stwlon3
                y1 = stwlat3
                x2 = stwlon4
                y2 = stwlat4
                x3 = stwlon1
                y3 = stwlat1
                totalarea = abs((x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1))
                            
                if((subarea1+subarea2+subarea3) <= &
                (totalarea*1.00001d0)) then
                    rsij(1,n,s) = 2
                    rsij(2,n,s) = i
                    rsij(3,n,s) = j
                    rscnf(1,n,s)=subarea1/totalarea
                    rscnf(2,n,s)=subarea2/totalarea
                    rscnf(3,n,s)=subarea3/totalarea
                    exit cellloop
                endif
            enddo cellloop
        enddo

        deallocate(stwlonlat)
        deallocate(nCellsInSearchBins)
        deallocate(searchBins)
        deallocate(searchBinPartitions)

    enddo

! pen .rad files if necessary
    do s=1,ns
        rfids(s) = 230 + s
        open(rfids(s),file=TRIM(GBLINPUTDIR)//'/'//trim(radfiles(s)), &
        status='old')
        if(fullplanes(s) == 1) then
            read (rfids(s),*) ni, nj, dxinc, dyinc
        else
            read (rfids(s),*) ni, nj, dxinc
            dyinc = dxinc
        endif
        if((ni /= nis(s)) .OR. (nj /= njs(s)) .OR. &
        (dxinc /= dxincs(s)) .OR. (dyinc /= dyincs(s))) then
            if(myproc == 0) then
                write(screenunit,*) ''
                write(screenunit,*) 'FORT.23 INCONSISTENT WITH .RAD FILES'
                write(screenunit,*) 'PROGRAM WILL BE TERMINATED'
                write(screenunit,*) ''
            endif
            write(16,*) ''
            write(16,*) 'FORT.23 SEEMS INCONSISTENT WITH .RAD FILES'
            write(16,*) 'PROGRAM WILL BE TERMINATED'
            write(16,*) ''
#ifdef CMPI
            call msg_fini()
#endif
            stop 'RS2INIT-1'
        endif
    enddo

! kip snaps if nbs(s) < 0
    do s=1,ns
        do
        if(nbs(s) >= 0) exit
        nbs(s) = nbs(s) + 1
        call rs2get(RSNX,RSNY,NP)
    enddo
enddo

! nitialize some variables
    cntSnaps = 0
    endoffile(:) = 0

    if(myproc == 0) then
        write(screenunit,*)
        write(screenunit,*) 'INIT OF RADIATION STRESS ARRAYS IS DONE.'
    endif
    write(16,*)
    write(16,*) 'INIT OF RADIATION STRESS ARRAYS IS DONE.'

    return

    1004 FORMAT(//,1X,' NUMBER OF SETS WAS SPECIFIED'// &
    'INCORRECTLY IN UNIT23.'/ &
    ' IT MUST BE GREATER THAN 0.'//)

    99999 CONTINUE

#ifdef CMPI
    call msg_fini()
#endif
    STOP 'RADIATION STRESS READ ERROR (RS2-1)'

    END SUBROUTINE


!***********************************************************************
!   SUBROUTINE RS2GET
!***********************************************************************

    subroutine rs2get(rsnx,rsny,np)
    implicit none

    integer,intent(in) :: np
    real(sz),intent(out) :: rsnx(:),rsny(:)

    if(type == 1) then
        stop 'RS2GET: SORRY! TYPE-1 RS FILE NOT YET SUPPORTED.'
    !        call rs2get_type1(rsnx,rsny,np)
    else if(type == 2) then
        call rs2get_type2(rsnx,rsny,np)
    else
        stop 'RS2 TYPE IS NEITHER 1 NOR 2.  TERMINATED.'
    endif

    return

    end subroutine

!$$$C***********************************************************************
!$$$C   SUBROUTINE RS2GET_TYPE1
!$$$C***********************************************************************
!$$$
!$$$      subroutine rs2get_type1(rsnx,rsny,np)
!$$$      use sizes
!$$$#ifdef CMPI
!$$$      use messenger, only: msg_fini, msg_ibcast, msg_cbcast, msg_rbcast
!$$$#ifdef HAVE_MPI_MOD
!$$$      use mpi
!$$$      implicit none
!$$$#else
!$$$      implicit none
!$$$      include 'mpif.h'
!$$$#endif
!$$$#else
!$$$      implicit none
!$$$#endif
!$$$      integer,intent(in) :: np
!$$$      real(sz),intent(out) :: rsnx(:),rsny(:)
!$$$      integer :: i,j,k,s,n,ni,nj,fid,ios(1)
!$$$      integer :: isActive(1)
!$$$      character(100) dummyc, spectrum_id(1)
!$$$      real(sz), dimension(:,:,:), allocatable :: rsxy
!$$$      real(sz) :: rsx,rsy,rsx1,rsx2,rsx3,rsy1,rsy2,rsy3,cnf1,cnf2,cnf3
!$$$      real(sz) :: angle_st,pi
!$$$
!$$$      stop 'RS2GET_TYPE1: SORRY! TYPE-1 RS FILE NOT YET SUPPORTED.'
!$$$
!$$$      if(myproc.eq.0) then
!$$$        write(screenunit,*)
!$$$        write(screenunit,*) 'READING IN RADIATION STRESS'
!$$$      endif
!$$$      write(16,*)
!$$$      write(16,*) 'READING IN RADIATION STRESS'
!$$$
!$$$      pi = 4.d0*atan(1.d0)
!$$$
!$$$      !increment counter
!$$$      cntSnaps = cntSnaps+1
!$$$
!$$$      !initialize
!$$$      do i=1,np
!$$$        rsnx(I)=0.d0
!$$$        rsny(I)=0.d0
!$$$      enddo
!$$$
!$$$
!$$$      fid = 23
!$$$#ifdef CMPI
!$$$      if(myproc.eq.0) read(fid,'(A)',iostat=ios(1)) dummyc
!$$$      call msg_ibcast(ios,1)
!$$$      call msg_cbcast(dummyc,100)
!$$$#else
!$$$      read(fid,'(A)',iostat=ios(1)) dummyc
!$$$#endif
!$$$      write(16,*) 'READING RADIATION STRESS SNAP.  '//
!$$$     &              'SNAP ID = ', trim(dummyc)
!$$$      !if fort.23 had reached EOF or got any read error,
!$$$      ! stop reading fort.23 and use insert blank snaps from now on.
!$$$      if(ios(1) < 0) then
!$$$        if(myproc.eq.0) then
!$$$          write(screenunit,*)
!$$$          write(screenunit,*) 'UNIT23 FILE RUN OUT.'
!$$$          write(screenunit,*) 'BLANK SNAP INSERTED FROM NOW ON.'
!$$$        endif
!$$$        write(16,*)
!$$$        write(16,*) 'UNIT23 FILE RUN OUT.'
!$$$        write(16,*) 'BLANK SNAP WILL BE INSERTED FROM NOW ON.'
!$$$        do s=1,ns
!$$$          endoffile(s) = 1
!$$$        enddo
!$$$        return
!$$$      endif
!$$$
!$$$!-----------
!$$$! S-Loop
!$$$!-----------
!$$$
!$$$      do s=1,ns   ! put a blank snap for the first nbs(s) snaps
!$$$        if(cntsnaps.le.nbs(s)) then
!$$$          if(abs(nscreen) >= 1) then
!$$$            if(myproc.eq.0) then
!$$$              write(screenunit,*) 'INSERTING A BLANK RS SNAP, REGION = ',s,
!$$$     &                   ', COUNT = ',cntsnaps
!$$$            endif
!$$$            write(16,*) 'INSERTING A BLANK RS SNAP, REGION = ',s,
!$$$     &                  ', COUNT = ',cntsnaps
!$$$          endif
!$$$          cycle
!$$$        endif
!$$$
!$$$        ni = nis(s)
!$$$        nj = njs(s)
!$$$
!$$$        allocate(rsxy(2,nis(s),njs(s)))
!$$$
!$$$        fid = 23
!$$$        if (myproc.eq.0) read(fid,*,iostat=ios(1)) ! skip a line
!$$$#ifdef CMPI
!$$$        call msg_ibcast(ios,1)
!$$$#endif
!$$$        if (ios(1) > 0) goto 9999
!$$$
!$$$        if(myproc.eq.0) read(fid,*,iostat=ios(1)) isActive(1) ! active flag
!$$$
!$$$#ifdef CMPI
!$$$        call msg_ibcast(ios,1)
!$$$        call msg_ibcast(isActive,1)
!$$$#endif
!$$$
!$$$        if(ios(1) > 0)       goto 9999
!$$$        if(isActive(1).eq.0) goto 100
!$$$
!$$$        if(myproc.eq.0) then
!$$$          do j = nj, 1, -1
!$$$            read (fid,'(5E15.7)',iostat=ios(1))
!$$$     &           ((rsxy(k,i,j),k=1,2),i=1,ni)
!$$$          enddo
!$$$        endif
!$$$
!$$$#ifdef CMPI
!$$$        call msg_ibcast(ios,1)
!$$$        call msg_rbcast(rsxy,2*ni*nj)
!$$$#endif
!$$$        if(ios(1) > 0)       goto 9999
!$$$
!$$$C--------------------------------------------------------------------------
!$$$C   Interpolate RS Snap onto local ADCIRC mesh
!$$$C--------------------------------------------------------------------------
!$$$
!$$$        do n=1,np
!$$$          if(rsij(1,n).ne.s) cycle !node n does not reside in STWAVE grid s
!$$$
!$$$          i = rsij(3,n)
!$$$          j = rsij(4,n)
!$$$
!$$$          if(rsij(2,n).eq.1) then !node n resides in triangle 1 in cell (i,j).
!$$$            rsx1 = rsxy(1,i,j)
!$$$            rsy1 = rsxy(2,i,j)
!$$$            rsx2 = rsxy(1,i+1,j)
!$$$            rsy2 = rsxy(2,i+1,j)
!$$$            rsx3 = rsxy(1,i+1,j+1)
!$$$            rsy3 = rsxy(2,i+1,j+1)
!$$$          else                    !node n resides in triangle 2 in cell (i,j).
!$$$            rsx1 = rsxy(1,i+1,j+1)
!$$$            rsy1 = rsxy(2,i+1,j+1)
!$$$            rsx2 = rsxy(1,i,j+1)
!$$$            rsy2 = rsxy(2,i,j+1)
!$$$            rsx3 = rsxy(1,i,j)
!$$$            rsy3 = rsxy(2,i,j)
!$$$          endif
!$$$
!$$$          !interpolation
!$$$          rsx= rscnf(1,n)*rsx1 + rscnf(2,n)*rsx2 + rscnf(3,n)*rsx3
!$$$          rsy= rscnf(1,n)*rsy1 + rscnf(2,n)*rsy2 + rscnf(3,n)*rsy3
!$$$
!$$$          !need to rotate
!$$$          angle_st = (azimuths(s)*pi)/180.d0
!$$$          rsnx(n) = cos(angle_st)*rsx-sin(angle_st)*rsy
!$$$          rsny(n) = sin(angle_st)*rsx+cos(angle_st)*rsy
!$$$        enddo
!$$$
!$$$ 100    continue
!$$$        deallocate(rsxy)
!$$$
!$$$      enddo       ! end s-loop
!$$$
!$$$
!$$$      if(myproc.eq.0) then
!$$$        write(screenunit,*) 'RADIATION STRESS IS LOADED'
!$$$        write(screenunit,*)
!$$$      endif
!$$$      write(16,*) 'RADIATION STRESS IS LOADED'
!$$$      write(16,*)
!$$$
!$$$      return
!$$$
!$$$ 9999 continue
!$$$      if(myproc.eq.0) then
!$$$        write(screenunit,*) trim(radfiles(s)),' IS NOT COMPLETE.'
!$$$        write(screenunit,*) 'PROGRAM WILL BE TERMINATED.'
!$$$      endif
!$$$      write(16,*) trim(radfiles(s)),' IS NOT COMPLETE.'
!$$$      write(16,*) 'PROGRAM WILL BE TERMINATED.'
!$$$      close(16)
!$$$
!$$$#ifdef CMPI
!$$$      call msg_fini()
!$$$#endif
!$$$      stop
!$$$
!$$$      end subroutine

!***********************************************************************
!   SUBROUTINE RS2GET_TYPE2
!***********************************************************************

    subroutine rs2get_type2(rsnx,rsny,np)
    use sizes
#ifdef CMPI
! Casey 090327: Implement Seizo's changes for buffering the radiation stress gradients.
    use messenger, only: msg_fini, msg_ibcast, msg_cbcast, msg_rbcastd
#ifdef HAVE_MPI_MOD
    use mpi
    implicit none
#else
    implicit none
    include 'mpif.h'
#endif
#else
    implicit none
#endif

    integer,intent(in) :: np
    real(sz),intent(out) :: rsnx(:),rsny(:)
    integer :: i,j,k,s,n,ni,nj,fid,ios(1)
    integer :: isActive(1)
    character(100) dummyc, spectrum_id(1)
    real(sz), dimension(:,:,:), allocatable :: rsxy
    real(sz) :: rsx,rsy,rsx1,rsx2,rsx3,rsy1,rsy2,rsy3,cnf1,cnf2,cnf3
    real(sz) :: angle_st,pi

    if(myproc == 0) then
        write(screenunit,*)
        write(screenunit,*) 'READING IN RADIATION STRESS'
    endif
    write(16,*)
    write(16,*) 'READING IN RADIATION STRESS'

    pi = 4.d0*atan(1.d0)

! ncrement counter
    cntSnaps = cntSnaps+1

! nitialize
    do i=1,np
        rsnx(I)=0.d0
        rsny(I)=0.d0
    enddo


!-----------
! S-Loop
!-----------

    do s=1,ns   ! put a blank snap for the first nbs(s) snaps
        if(cntsnaps <= nbs(s)) then
            if(abs(nscreen) >= 1) then
                if(myproc == 0) then
                    write(screenunit,*) 'INSERTING A BLANK RS SNAP, REGION = ',s, &
                    ', COUNT = ',cntsnaps
                endif
                write(16,*) 'INSERTING A BLANK RS SNAP, REGION = ',s, &
                ', COUNT = ',cntsnaps
            endif
            cycle
        endif

    ! if previous EOF write blank RS SNAP
        if(endoffile(s) == 1) then
            if (myproc == 0) then
                write(screenunit,*) trim(radfiles(s)), &
                ' PREVIOUSLY REACHED EOF.'
                write(screenunit,'(A25,I6,A9,I6)') &
                " BLANK RS SNAP, REGION = ",s, " COUNT = ", cntsnaps
            endif
            cycle
        endif

        ni = nis(s)
        nj = njs(s)

        allocate(rsxy(2,nis(s),njs(s)))

        fid = rfids(s)

        if (endoffile(s) == 0) then
            if (myproc == 0) then
                read (fid,'(A)',iostat=ios(1)) spectrum_id(1) ! spectrum identifier
            endif
#ifdef CMPI
            call msg_ibcast(ios,1)
#endif

            if (ios(1) < 0) then
                endoffile(s) = 1
                if (myproc == 0) then
                    write(screenunit,*) trim(radfiles(s)),' REACHED EOF.'
                    write(screenunit,'(A25,I6,A9,I6)') &
                    " BLANK RS SNAP, REGION = ",s, " COUNT = ", cntsnaps
                endif
                goto 100
            else
                if (myproc == 0) then
                    write(screenunit,'(A)') " spectrum_id = "//trim(spectrum_id(1))
                    write(screenunit,'(A19,I6,A9,I6)') &
                    " RS SNAP, REGION = ",s, " COUNT = ", cntsnaps
                endif
            endif
        endif
                 
#ifdef CMPI
        if(myproc == 0) then
            ios(1) = 0
            do j = nj, 1, -1
                read (fid,*,err=511) &
                ((rsxy(k,i,j),k=1,2),i=1,ni)
            enddo
            go to 512
            511 ios(1) = 1
            512 continue
        endif
        call msg_ibcast(ios,1)
        if(ios(1) == 0) then
        ! Casey 090327: Implement Seizo's changes for buffering the radiation stress gradients.
            call msg_rbcastd(rsxy,2,ni,nj)
        else
            goto 9999      !  read error because of missing snap data
        endif
#else
        do j = nj, 1, -1
            read (fid,*,err=9999) &
            ((rsxy(k,i,j),k=1,2),i=1,ni)
        enddo
#endif

    !--------------------------------------------------------------------------
    !   Interpolate RS Snap onto local ADCIRC mesh
    !--------------------------------------------------------------------------

        do n=1,np
            if(rsij(1,n,s) == 0) cycle !node n does not reside in STWAVE grid s
                        
            i = rsij(2,n,s)
            j = rsij(3,n,s)

            if(rsij(1,n,s) == 1) then !node n resides in triangle 1 in cell (i,j).
                rsx1 = rsxy(1,i,j)
                rsy1 = rsxy(2,i,j)
                rsx2 = rsxy(1,i+1,j)
                rsy2 = rsxy(2,i+1,j)
                rsx3 = rsxy(1,i+1,j+1)
                rsy3 = rsxy(2,i+1,j+1)
            else                    !node n resides in triangle 2 in cell (i,j).
                rsx1 = rsxy(1,i+1,j+1)
                rsy1 = rsxy(2,i+1,j+1)
                rsx2 = rsxy(1,i,j+1)
                rsy2 = rsxy(2,i,j+1)
                rsx3 = rsxy(1,i,j)
                rsy3 = rsxy(2,i,j)
            endif

        ! nterpolation
            rsx= rscnf(1,n,s)*rsx1 + rscnf(2,n,s)*rsx2 + rscnf(3,n,s)*rsx3
            rsy= rscnf(1,n,s)*rsy1 + rscnf(2,n,s)*rsy2 + rscnf(3,n,s)*rsy3

            if((rsnx(n)*rsnx(n)+rsny(n)*rsny(n)) < (rsx*rsx+rsy*rsy))then
            ! eed to rotate
                angle_st = (azimuths(s)*pi)/180.d0
                rsnx(n) = cos(angle_st)*rsx-sin(angle_st)*rsy
                rsny(n) = sin(angle_st)*rsx+cos(angle_st)*rsy
            endif
        enddo

        100 continue
        deallocate(rsxy)

    enddo       ! end s-loop

    if(myproc == 0) then
        write(screenunit,*) 'RADIATION STRESS IS LOADED'
        write(screenunit,*)
    endif
    write(16,*) 'RADIATION STRESS IS LOADED'
    write(16,*)

    return

    9999 continue
    if(myproc == 0) then
        write(screenunit,*) trim(radfiles(s)),' IS NOT COMPLETE.'
        write(screenunit,*) 'PROGRAM WILL BE TERMINATED.'
    endif
    write(16,*) trim(radfiles(s)),' IS NOT COMPLETE.'
    write(16,*) 'PROGRAM WILL BE TERMINATED.'
    close(16)

#ifdef CMPI
    call msg_fini()
#endif
    stop

    end subroutine

    END MODULE RS2
          
          
