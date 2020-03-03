      program buildstwave23

      implicit none

      integer :: i,j,k,e,s,idum,cnt
!$omp threadprivate(i,j,k,e,s,idum,cnt)
      character(len=100) :: inputfile,outputfile,dummyc,title
!$omp threadprivate(inputfile,outputfile,dummyc,title)
      double precision :: x,y,rdum

      !variables for input file
      integer :: type ! type of output file to produce, 1: include wave radiation
!$omp threadprivate(type)
                      !                                 2: does not include wave radiation
      integer :: ns ! number of stwave grids
!$omp threadprivate(ns)
      character(len=100),allocatable,dimension(:) :: simfiles ! names of .sim files
!$omp threadprivate(simfiles)
      character(len=100),allocatable,dimension(:) :: radfiles ! names of .rad files
!$omp threadprivate(radfiles)
      character(len=100),allocatable,dimension(:) :: spgfiles ! names of grid files in state plane coordinates
!$omp threadprivate(spgfiles)
      integer,allocatable,dimension(:)            :: fullplanes ! 0: half plane, 1: full plane
!$omp threadprivate(fullplanes)
      integer,allocatable,dimension(:)            :: nbs ! number of blank snaps
!$omp threadprivate(nbs)
      integer,allocatable,dimension(:)            :: rfids ! id of .rad files
!$omp threadprivate(rfids)

      !variables for a grid in Geophysical coordiantes
      integer :: np,ne
!$omp threadprivate(np,ne)
      character(len=100) :: geoGridFile
!$omp threadprivate(geoGridFile)
      double precision,allocatable,dimension(:,:) :: xy
      integer,allocatable,dimension(:,:) :: nm
!$omp threadprivate(nm)

      !variables for a grid in state plane coordiantes
      integer :: spnp,spne
!$omp threadprivate(spnp,spne)
      double precision,allocatable,dimension(:,:) :: spxy
      integer,allocatable,dimension(:,:) :: spnm
!$omp threadprivate(spnm)

      !variables for .sim files
      double precision :: x0,y0,azimuth
      double precision,allocatable,dimension(:) :: x0s,y0s,azimuths

      !variables for stwave grid
      integer :: ni,nj
!$omp threadprivate(ni,nj)
      integer,allocatable,dimension(:) :: nis,njs
!$omp threadprivate(nis,njs)
      double precision :: dxinc,dyinc
      double precision,allocatable,dimension(:) :: dxincs,dyincs
      double precision,allocatable,dimension(:,:,:) :: stwxy,stwlonlat

      !variables for searching
      integer :: numElemInside
!$omp threadprivate(numElemInside)
      integer,allocatable,dimension(:) :: elemInside
!$omp threadprivate(elemInside)
      integer,allocatable,dimension(:,:) :: elemInsideBins
!$omp threadprivate(elemInsideBins)
      integer,allocatable,dimension(:) :: numElemInsideBins
!$omp threadprivate(numElemInsideBins)
      double precision,allocatable,dimension(:) :: elemInsideBinPartition
      integer :: elem,n1,n2,n3
!$omp threadprivate(elem,n1,n2,n3)
      double precision :: w1,w2,w3,x1,x2,x3,y1,y2,y3
      integer,parameter :: numBins = 600
      double precision :: stwminx,stwmaxx
      real :: tol,exmax,exmin,p1,p2
!$omp threadprivate(tol,exmax,exmin,p1,p2)
      integer :: i1,i2,j1,j2
!$omp threadprivate(i1,i2,j1,j2)

      !variables for snaps
      integer,allocatable,dimension(:) :: endoffile
!$omp threadprivate(endoffile)
      integer :: ios,snapCnt,fid,spectrum_id
!$omp threadprivate(ios,snapCnt,fid,spectrum_id)
      double precision,allocatable,dimension(:,:,:) :: wxyrs


      write(*,*) ""
      write(*,'(A)',ADVANCE='NO') "enter name of input file: "
      read(*,*) inputfile

      write(*,*) ""
      write(*,*) ""
      write(*,'(A)',ADVANCE='NO') "enter name of output file: "
      read(*,*) outputfile
      write(*,*) ""
      write(*,*) ""

      open(1+CK_LUN,file=inputfile,status="old",action="read")
      open(2+CK_LUN,file=outputfile,action="write")

      read(1+CK_LUN,*) type
      if(type.ne.1.and.type.ne.2) then
         print *, 'type must be 1 or 2.'
         print *, 'program will be terminated.'
      endif
      read(1+CK_LUN,*) geoGridFile
      read(1+CK_LUN,*) ns

      write(2+CK_LUN,'(i1,a)') type, 
     &' ! 1: rad vecs included; 2: rad vecs not included'
      write(2+CK_LUN,'(i1,a)') ns, ' ! number of stwave grids'


      allocate(simfiles(ns),radfiles(ns),spgfiles(ns))
      allocate(fullplanes(ns),nbs(ns),rfids(ns))
      allocate(nis(ns),njs(ns),dxincs(ns),dyincs(ns))
      allocate(x0s(ns),y0s(ns),azimuths(ns))
      allocate(endoffile(ns))

      do i=1,ns
         read(1+CK_LUN,*)  ! skip a line
         read(1+CK_LUN,'(A100)') simfiles(i)
         read(1+CK_LUN,'(A100)') radfiles(i)
         read(1+CK_LUN,'(A100)') spgfiles(i)
         read(1+CK_LUN,*) fullplanes(i)
         read(1+CK_LUN,*) nbs(i)
         if(type.eq.2) then
            read(1+CK_LUN, *) x0s(i), y0s(i), azimuths(i)
            if(fullplanes(i).eq.1) then
               read (1+CK_LUN,*) nis(i), njs(i), dxincs(i), dyincs(i)
            else
               read (1+CK_LUN,*) nis(i), njs(i), dxincs(i)
               dyincs(i) = dxincs(i)
            endif
         endif
      enddo

      close(1+CK_LUN)

      !read adcirc grid in geophysical coordinates
      print *, 'reading adcirc grid :', geoGridFile
      open(14+CK_LUN,file=geoGridFile,status='old',action='read')
      read(14+CK_LUN,'(A80)') title
      read(14+CK_LUN,*) ne, np

      print *, ' grid id = ', trim(title)
      print *, ' np = ',np,', ne = ',ne
    
      allocate(xy(2,np),nm(3,ne))

      do i=1,np
         read(14+CK_LUN,*) idum, xy(1,i), xy(2,i), rdum
      enddo

      do i=1,ne
         read(14+CK_LUN,*) idum, idum, nm(1,i), nm(2,i), nm(3,i)
      enddo

      close(14+CK_LUN)

      allocate(elemInside(ne))

      !process each stwave grid header info
      do s=1,ns
         write(*,*) ""
         write(*,*) ""
         print *, 'reading & producing header info of stwave grid ',s

         if(type.eq.1) then
            !simfile
            open(1+CK_LUN,file=simfiles(s),status="old",action="read")
            read (1+CK_LUN, *) dummyc, x0, y0, azimuth

            print *, ' stwave grid info'
            print *, ' ID     = ', dummyc
            print *, ' x0     = ', x0
            print *, ' y0     = ', y0
            print *, ' azimuth= ',azimuth

            close(1+CK_LUN)

            !radfile
            rfids(s) = s+20
            open(rfids(s+CK_LUN),file=radfiles(s),status="old",action="read")
            if(fullplanes(s).eq.1) then
               read (rfids(s+CK_LUN),*) ni, nj, dxinc, dyinc
            else
               read (rfids(s+CK_LUN),*) ni, nj, dxinc
               dyinc = dxinc
            endif

            nis(s) = ni
            njs(s) = nj
         else
            x0 = x0s(s)
            y0 = y0s(s)
            azimuth = azimuths(s)
            ni = nis(s)
            nj = njs(s)
            dxinc = dxincs(s)
            dyinc = dyincs(s)
         endif

         print *, ''
         print *, ' ni = ',ni,', nj = ',nj,', dxinc = ',real(dxinc),
     &     ', dyinc = ',real(dyinc)

         allocate(stwxy(2,ni,nj),stwlonlat(2,ni,nj))

         print *, ""
         print *, 'making stwave grid'

         do j=1,nj
            do i=1,ni
               stwxy(1,i,j)=dxinc*(real(i)-0.5)
               stwxy(2,i,j)=dyinc*(real(j)-0.5)
            end do
         end do

         !read adcirc grid in state plane coordinates
         print *, ''
         print *, 'reading adcirc grid :', trim(spgfiles(s))

         open(14+CK_LUN,file=spgfiles(s),status='old',action='read')
         read(14+CK_LUN,'(A80)') title
         read(14+CK_LUN,*) spne, spnp

         print *, ' grid id = ', trim(title)
         print *, ' np = ',np,', ne = ',ne
    
         allocate(spxy(2,np),spnm(3,ne))
         
         do i=1,np
            read(14+CK_LUN,*) idum, spxy(1,i), spxy(2,i), rdum
         enddo

         do i=1,ne
            read(14+CK_LUN,*) idum, idum, spnm(1,i), spnm(2,i), spnm(3,i)
         enddo

         close(14+CK_LUN)
    

         if(spnp.ne.np) then
            print *, ''
            print *, ' number of nodes in state plane grid does not '//
     &        'match lat-lon grid.'
            print *, ' program will be terminated.'
            stop
         endif

         if(spne.ne.ne) then
            print *, ''
            print *, ' number of elements in state plane grid does '//
     &        'not match lat-lon grid.'
            print *, ' program will be terminated.'
            stop
         endif

         do j=1,ne
            do i=1,3
               if(nm(i,j).ne.spnm(i,j)) then
                  print *, ''
                  print *, ' element-node map in state plane grid '//
     &              'does not match lat-lon grid.'
                  print *, ' program will be terminated.'
                  stop
               endif
            enddo
         enddo

         print *, ''
         print *, 'rotating state plane grid'
         call rotate_adcirc_mesh(x0,y0,azimuth,spnp,spxy)

         print *, ''
         print *, 'initilizing search info'
         call search_init(ni,nj,stwxy,np,ne,spxy,nm,numElemInside,
     &     elemInside)
         print *, '  ', numElemInside,' elements found'

         allocate(elemInsideBins(numBins,numElemInside))
         allocate(numElemInsideBins(numBins))
         allocate(elemInsideBinPartition(numBins+1))

         stwminx = stwxy(1,1,1)
         stwmaxx = stwxy(1,ni,nj)

         elemInsideBinPartition(1) = stwminx
         do i=1,numBins
            elemInsideBinPartition(i+1) = stwminx +
     &       (stwmaxx-stwminx)*real(i)/real(numBins)
            numElemInsideBins(i) = 0
         enddo

         print *, ''
         print *, 'preparing searching bins'
         tol = abs(elemInsideBinPartition(1)-
     &     elemInsideBinPartition(0))*0.01d0
         do k=1,numBins
            p1 = elemInsideBinPartition(k)
            p2 = elemInsideBinPartition(k+1)

            do j=1,numElemInside
               e = elemInside(j)
               exmax = spxy(1,nm(1,e))
               if(exmax.lt.spxy(1,nm(2,e))) exmax = spxy(1,nm(2,e))
               if(exmax.lt.spxy(1,nm(3,e))) exmax = spxy(1,nm(3,e))

               exmin = spxy(1,nm(1,e))
               if(exmin.gt.spxy(1,nm(2,e))) exmin = spxy(1,nm(2,e))
               if(exmin.gt.spxy(1,nm(3,e))) exmin = spxy(1,nm(3,e))

               if(exmax.ge.(p1-tol).and.exmin.le.(p2+tol)) then
                  numElemInsideBins(k) = numElemInsideBins(k) + 1
                  elemInsideBins(k,numElemInsideBins(k)) = e
               endif
            enddo
            print *, ' ',numElemInsideBins(k),
     &        ' elements were found in bin ',k
         enddo

         print *, ''
         print *, 'searching'
         cnt = 0
         do j=1,nj
            do i=1,ni
               x = stwxy(1,i,j)
               y = stwxy(2,i,j)
               call search(x,y,np,ne,spxy,nm,numElemInside,numBins,
     &           elemInsideBinPartition,numElemInsideBins,
     &           elemInsideBins,elem,w1,w2,w3)

               if(elem.eq.-1) then
                  cnt = cnt + 1
                  stwlonlat(1,i,j) = 0
                  stwlonlat(2,i,j) = 0
               else
                  n1 = nm(1,elem)
                  n2 = nm(2,elem)
                  n3 = nm(3,elem)
                  x1 = xy(1,n1)
                  x2 = xy(1,n2)
                  x3 = xy(1,n3)
                  y1 = xy(2,n1)
                  y2 = xy(2,n2)
                  y3 = xy(2,n3)
                  stwlonlat(1,i,j) = w1*x1+w2*x2+w3*x3
                  stwlonlat(2,i,j) = w1*y1+w2*y2+w3*y3
               endif
            enddo
         enddo

         print *, ''
         print *, 'filling blank coordinates'
         cnt = 0
         do j=1,nj
            do i=1,ni
               x = stwlonlat(1,i,j)
               y = stwlonlat(2,i,j)

               if(x.ne.0.d0.or.y.ne.0.d0) cycle

               do i1=i-1,1,-1
                  if(stwlonlat(1,i1,j).ne.0.d0) exit
               enddo
               do i2=i+1,ni,+1
                  if(stwlonlat(1,i2,j).ne.0.d0) exit
               enddo
               do j1=j-1,1,-1
                  if(stwlonlat(1,i,j1).ne.0.d0) exit
               enddo
               do j2=j+1,nj,+1
                  if(stwlonlat(1,i,j2).ne.0.d0) exit
               enddo
               if(i1.eq.0.and.i2.eq.ni+1) exit
               if(j1.eq.0.and.j2.eq.ni+1) exit

               !interpolation of x
               if(i1.eq.0) then
                  if(i2.ge.ni) then
                     x = 0.d0
                  else
                     do i1=i2+1,ni,+1
                        if(stwlonlat(1,i1,j).ne.0.d0) exit
                     enddo
                     if(i1.eq.ni+1) exit
                     x = stwlonlat(1,i2,j) -
     &                 (stwlonlat(1,i1,j) - stwlonlat(1,i2,j))*
     &                 real(i2-i)/real(i1-i2)
                  endif
                  if(x.gt.0) print *, '1 x,y=',x,y,i,j
               else if(i2.eq.(ni+1)) then
                  if(i1.le.0) then
                     x = 0.d0
                  else
                     do i2=i1-1,1,-1
                        if(stwlonlat(1,i2,j).ne.0.d0) exit
                     enddo
                     if(i1.eq.0) exit
                     x = stwlonlat(1,i1,j) +
     &                 (stwlonlat(1,i1,j) - stwlonlat(1,i2,j))*
     &                 real(i-i1)/real(i1-i2)
                  endif
                  if(x.gt.0) then
                     write(0+CK_LUN,*) '2 x,y=',x,y,i,j,stwlonlat(1,i1,j),
     &                 (stwlonlat(1,i1+1,j) - stwlonlat(1,i1,j))*(i-i1),
     &                 (stwlonlat(1,i1+1,j) - stwlonlat(1,i1,j)), (i-i1)
                     stop
                  endif
               else if(i1.ne.0.and.i2.ne.(ni+1)) then
                  w1 = real(i2-i)/(i2-i1)
                  w2 = real(i-i1)/(i2-i1)

                  if(w1.lt.0.d0) stop 'w1.lt.0.d0 (x)'
                  if(w2.lt.0.d0) stop 'w2.lt.0.d0 (x)'

                  x = w1*stwlonlat(1,i1,j) + w2*stwlonlat(1,i2,j)
                  if(x.gt.0) print *, '3 x,y=',x,y,i,j
               else
                  if(i1.eq.0) exit
               endif

               !interpolation of y
               if(j1.eq.0) then
                  if(j2.ge.nj) then
                     y = 0.d0
                  else
                     do j1=j2+1,nj,+1
                        if(stwlonlat(2,i,j1).ne.0.d0) exit
                     enddo
                     if(j1.eq.nj+1) exit
                     y = stwlonlat(2,i,j2) -
     &                 (stwlonlat(2,i,j1) - stwlonlat(2,i,j2))*
     &                 real(j2-j)/real(j1-j2)
                  endif
               else if(j2.eq.(nj+1)) then
                  if(j1.le.0) then
                     y = 0.d0
                  else
                     do j2=j1-1,1,-1
                        if(stwlonlat(2,i,j2).ne.0.d0) exit
                     enddo
                     if(j1.eq.0) exit
                     y = stwlonlat(2,i,j1) +
     &                 (stwlonlat(2,i,j1) - stwlonlat(2,i,j2))*
     &                 real(j-j1)/real(j1-j2)
                  endif
               else if(j1.ne.0.and.j2.ne.(nj+1)) then
                  w1 = real(j2-j)/(j2-j1)
                  w2 = real(j-j1)/(j2-j1)

                  if(w1.lt.0.d0) stop 'w1.lt.0.d0 (y)'
                  if(w2.lt.0.d0) stop 'w2.lt.0.d0 (y)'

                  y = w1*stwlonlat(2,i,j1) + w2*stwlonlat(2,i,j2)
               else
                  if(i1.eq.0) exit
               endif

               if(x.eq.0.or.y.eq.0) then
                  x = 0.d0
                  y = 0.d0
               endif

               stwlonlat(1,i,j) = x
               stwlonlat(2,i,j) = y

            enddo
         enddo

         print *, ''
         print *, 'writing header'

         write(2+CK_LUN,'(a,1x,i2)') '# grid ', s
         if(type.eq.2) then
            write(2+CK_LUN,'(a)') trim(simfiles(s))
            write(2+CK_LUN,'(a)') trim(radfiles(s))
            write(2+CK_LUN,'(i2,1x,a)') nbs(s), 
     &      ' ! number of blank snaps to be inserted at the beginning'
            write(2+CK_LUN,'(3(f14.4,1x),a)') real(x0s(s)), real(y0s(s)), 
     &      real(azimuths(s)),'  !x0 y0 azimuths'
         endif

         write(2+CK_LUN,'(i2,1x,a)') 
     &     fullplanes(s), ' ! 0: half plane 1: full plane'

         if(fullplanes(s).eq.1) then
            write(2+CK_LUN,100) nis(s), njs(s), real(dxincs(s)), 
     &                   real(dyincs(s)),'  !ni nj dxinc dyinc'
         else
            write(2+CK_LUN,101) nis(s), njs(s), real(dxincs(s)),' !ni nj dxinc'
         endif

 100     FORMAT(i6,i6,1x,f10.4,1x,f10.4,1x,a)
 101     FORMAT(i6,i6,1x,f10.4,1x,a)

         do j=1,nj
            write(2+CK_LUN,'(5E15.7)') ((stwlonlat(k,i,j), k=1,2), i=1,ni)
         enddo

         deallocate(stwxy)
         deallocate(stwlonlat)
         deallocate(spxy)
         deallocate(spnm)
         deallocate(elemInsideBins)
         deallocate(numElemInsideBins)
         deallocate(elemInsideBinPartition)
      enddo

      !type 2 format does not include actual wave radiation values
      if(type.eq.2) goto 9999

      !read and write snaps
      print *, ''
      print *, ' writing out snaps'

      snapCnt = 0
      endoffile(:) = 0
      outer: do
         snapCnt = snapCnt + 1
         print *, ''
         print *, ' snap ',snapCnt
         write(2+CK_LUN,*) ' # snap ', snapCnt
         do s=1,ns
            write(2+CK_LUN,*) ' ## grid ', s
            if(nbs(s).gt.0) then
               nbs(s) = nbs(s) - 1
               write(2+CK_LUN,*) 0,'  !this snap is 0: inactive, 1: active'
               cycle
            endif

            fid = rfids(s)

            read (fid+CK_LUN,*,iostat=ios) spectrum_id ! spectrum identifier

            if(ios.ne.0) then
               if(endoffile(s).eq.0) then
                  print *, '  stwave grid ',s,' is done.'
               endif
               endoffile(s) = 1
               write(2+CK_LUN,*) 0,'  !this snap is 0: inactive, 1: active'
               cycle
            endif
               
            write(2+CK_LUN,*) 1,'  !this snap is 0: inactive, 1: active'

            print *, '  spectrum id = ', spectrum_id, '(grid ',s,')'
            
            ni = nis(s)
            nj = njs(s)

            allocate(wxyrs(2,ni,nj))

            do j = nj, 1, -1
                read (fid+CK_LUN,*,iostat=ios) ((wxyrs(k,i,j),k=1,2),i=1,ni)
                if(ios.ne.0) then
                   print *, ' .rad file is not complete.'
                   print *, ' program will be terminated.'
                   stop
                endif
                write(2+CK_LUN,'(5E15.7)',iostat=ios) ((wxyrs(k,i,j),k=1,2),i=1,ni)
                if(ios.ne.0) then
                   print *, ' error deteced in writing output file.'
                   print *, ' program will be terminated.'
                   stop
                endif
            enddo

            deallocate(wxyrs)
         enddo

         do s=1,ns
            if(endoffile(s).eq.0) cycle outer
         enddo

         exit

      enddo outer

 9999 continue

      close(2+CK_LUN)

      print *, ''
      print *, 'done.'

      stop
      end program

      subroutine search(x,y,np,ne,spxy,nm,numElemInside,numBins,
     &  elemInsideBinPartition, numElemInsideBins,elemInsideBins,
     &  elem,w1,w2,w3)

      implicit none

      integer,intent(in) :: np,ne
      integer,intent(in) :: nm(3,ne)
      double precision,intent(in) :: x,y
      double precision,intent(in) :: spxy(2,np)
      integer,intent(in) :: numBins,numElemInside
      double precision,intent(in) :: elemInsideBinPartition(numBins)
      integer,intent(in) :: numElemInsideBins(numBins)
      integer,intent(in) :: elemInsideBins(numBins,numElemInside)

      integer,intent(out) :: elem
      double precision,intent(out) :: w1,w2,w3

      integer :: i,k,e,n1,n2,n3
      double precision :: x1,x2,x3,y1,y2,y3,w

      integer :: min_rat_elem
      double precision :: min_rat

      elem = -1

      min_rat = 1d10

      do k=1,numBins
         if((x.ge.elemInsideBinPartition(k)).and.
     &    (x.le.elemInsideBinPartition(k+1))) exit
      enddo

      do i=1,numElemInsideBins(k)
         e = elemInsideBins(k,i)
         n1 = nm(1,e)
         n2 = nm(2,e)
         n3 = nm(3,e)

         x1 = x
         x2 = spxy(1,n2)
         x3 = spxy(1,n3)
         y1 = y
         y2 = spxy(2,n2)
         y3 = spxy(2,n3)
         w1 = abs((x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1))
         
         x1 = spxy(1,n1)
         x2 = x
         x3 = spxy(1,n3)
         y1 = spxy(2,n1)
         y2 = y
         y3 = spxy(2,n3)
         w2 = abs((x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1))
         
         x1 = spxy(1,n1)
         x2 = spxy(1,n2)
         x3 = x
         y1 = spxy(2,n1)
         y2 = spxy(2,n2)
         y3 = y
         w3 = abs((x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1))

         x1 = spxy(1,n1)
         x2 = spxy(1,n2)
         x3 = spxy(1,n3)
         y1 = spxy(2,n1)
         y2 = spxy(2,n2)
         y3 = spxy(2,n3)
         w  = abs((x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1))

         if(min_rat.gt.(w1+w2+w3)/w) then
            min_rat = (w1+w2+w3)/w
            min_rat_elem = e
         endif

         if((w1+w2+w3).le.(w*1.00001d0)) then
            elem = e
            w1 = w1/w
            w2 = w2/w
            w3 = w3/w
            exit
         endif

      enddo

      end subroutine

      subroutine search_init(ni,nj,stwxy,np,ne,spxy,nm,numElemInside,elemInside)

      implicit none

      integer,intent(in) :: ni,nj,np,ne
      integer,intent(in) :: nm(3,ne)
      double precision,intent(in) :: stwxy(2,ni,nj)
      double precision,intent(in) :: spxy(2,np)
      integer,intent(out) :: numElemInside
      integer,intent(out) :: elemInside(ne)

      integer i,j,n
      double precision :: stwminx,stwminy,stwmaxx,stwmaxy
      double precision :: x,y

      stwminx = stwxy(1,1,1)
      stwminy = stwxy(2,1,1)
      stwmaxx = stwxy(1,ni,nj)
      stwmaxy = stwxy(2,ni,nj)

      print *, '  stwave grid min x = ',real(stwminx),' max '// 
     &' x = ',real(stwmaxx)
      print *, '  stwave grid min y = ',real(stwminy),' max '// 
     &' y = ',real(stwmaxy)

      numElemInside = 0
      do i=1,ne
         do j=1,3
            n = nm(j,i)
            x = spxy(1,n)
            y = spxy(2,n)

            if((x.ge.stwminx).and.(x.le.stwmaxx).and.
     &        (y.ge.stwminy).and.(y.le.stwmaxy)) then
               numElemInside = numElemInside + 1
               elemInside(numElemInside) = i
               exit
            endif
         enddo
      enddo
      end subroutine

      subroutine rotate_adcirc_mesh(x0,y0,azimuth,np,xy)

      implicit none

      double precision,intent(in) :: x0,y0,azimuth
      integer,intent(in) :: np
      double precision,intent(in out) :: xy(2,np)

      integer :: i
      double precision :: x,y,angle,twopi,pi

      twopi = 8.0d0*atan(1.0d0)
      pi = 0.5d0*twopi
      angle = (azimuth*pi)/180.0d0

      do i=1,np
         x = xy(1,i)
         y = xy(2,i)

         xy(1,i)=cos(angle)*(x-x0) + sin(angle)*(y-y0)
         xy(2,i)=-sin(angle)*(x-x0 ) + cos(angle)*(y-y0)
      end do

      end subroutine
