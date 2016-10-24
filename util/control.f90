!-------------------------------------------------------------------------
! control.f90
!-------------------------------------------------------------------------
! Author: Jason Fleming (jason.fleming@seahorsecoastal.com) 
!
! Provide routines for reading, storing, and writing adcirc fort.15 data. 
!
!-------------------------------------------------------------------------
!
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!  M O D U L E   C O N T R O L
!---------------------------------------------------------------------
!---------------------------------------------------------------------
module control
use adcmesh
implicit none
#ifdef REAL4
integer, parameter :: sz = 4
#else
integer, parameter :: sz = 8
#endif
character(len=80) :: rundes
character(len=80) :: runid
integer :: nfover
real :: warnElev
integer :: iWarnElevDump
real :: warnElevDumpLimit
integer :: nabout
integer :: nscreen
integer :: ihot
integer :: ics
integer :: im
integer :: iden
integer :: nolibf
integer :: nolifa
integer :: nolica
integer :: nolicat
integer :: nwp
character(len=1024), allocatable :: attrName(:)
integer :: ncor
integer :: ntip
integer :: nws
integer :: nrs
integer :: ncice
integer :: nramp
real :: g
real :: tau0
real :: tau0FullDomainMin
real :: tau0FullDomainMax
real(8) :: dtdp
real :: statim
real :: reftim
real :: rstiminc
real :: wtiminc
integer :: irefyr
integer :: irefmo
integer :: irefday
integer :: irefhr
integer :: irefmin
integer :: irefsec
real :: refsec
integer :: nwlat
integer :: nwlon
real :: wlatmax
real :: wlonmin
real :: wlatinc
real :: wloninc
real :: cice_timinc
real :: pureVortex
real :: pureBackground
integer :: stormNumber
real :: bladj
real :: rnday
real :: dramp
real :: drampextflux
real :: fluxsettlingtime
real :: drampintflux
real :: drampelev
real :: dramptip
real :: drampmete
real :: drampwrad
real :: dunrampmete
real :: a00
real :: b00
real :: c00
real :: h0
integer :: nodedrymin
integer :: nodewetmin
real :: velmin
real :: tau
real :: cf
real :: hbreak
real :: ftheta
real :: fgamma
real :: eslm
real :: eslc
real :: cori
integer :: ntif
character(len=5), allocatable :: tipotag(:)
real, allocatable :: tpk(:)
real, allocatable :: amigt(:)
real, allocatable :: etrf(:)
real, allocatable :: fft(:)
real, allocatable :: facet(:)
integer :: nbfr
character(len=5), allocatable :: bountag(:)
real, allocatable :: amig(:)
real, allocatable :: ff(:)
real, allocatable :: face(:)
character(len=10), allocatable :: nbfr_alpha(:)
real, allocatable :: emo(:,:)
real, allocatable :: efa(:,:)
real :: anginn
integer :: nffr
character(len=5), allocatable :: fbountag(:)
real, allocatable :: famig(:)
real, allocatable :: fff(:)
real, allocatable :: fface(:)
character(len=10), allocatable :: nffr_alpha(:)
real, allocatable :: qnam(:,:)
real, allocatable :: qnph(:,:)
real, allocatable :: enam(:,:)
real, allocatable :: enph(:,:)
integer :: noute
real :: toutse
real :: toutfe
integer :: nspoole
integer :: nstae
real(sz), allocatable :: xel(:)
real(sz), allocatable :: yel(:)
integer :: noutv
real :: toutsv
real :: toutfv
integer :: nspoolv
integer :: nstav
real(sz), allocatable :: xev(:)
real(sz), allocatable :: yev(:)
integer :: noutc
real :: toutsc
real :: toutfc
integer :: nspoolc
integer :: nstac
real(sz), allocatable :: xec(:)
real(sz), allocatable :: yec(:)
integer :: noutm
real :: toutsm
real :: toutfm
integer :: nspoolm
integer :: nstam
real(sz), allocatable :: xem(:)
real(sz), allocatable :: yem(:)
integer :: noutge
real :: toutsge
real :: toutfge
integer :: nspoolge
integer :: noutgv
real :: toutsgv
real :: toutfgv
integer :: nspoolgv
integer :: noutgc
real :: toutsgc
real :: toutfgc
integer :: nspoolgc
integer :: noutgw
real :: toutsgw
real :: toutfgw
integer :: nspoolgw
integer :: nfreq
character(len=10), allocatable :: namefr(:)
real, allocatable :: hafreq(:)
real, allocatable :: haff(:)
real, allocatable :: haface(:)
real :: thas
real :: thaf
integer :: nhainc
real :: fmv
integer :: nhase
integer :: nhasv
integer :: nhage
integer :: nhagv
integer :: nhstar
integer :: nhsinc
integer :: ititer
integer :: isldia
real :: convcr
integer :: itmax
character(len=80) :: title 
character(len=80) :: institution
character(len=80) :: source
character(len=80) :: history
character(len=80) :: references
character(len=80) :: comments
character(len=80) :: host
character(len=80) :: convention
character(len=80) :: contact
character(len=80) :: base_date
!
logical :: readMetaData
integer, parameter :: metadataRequired(9)     &     
    = (/3, 5, 6, 367, 368, 567, 568, 667, 668/)
integer :: outputSpecifiers(9)
!
integer :: lineNum  ! line number being read from fort.15
integer :: echoLine ! line number being echoed to screen
!
!-----------------------------------
contains !   C O N T A I N S
!-----------------------------------
!
!---------------------------------------------------------------------
! S U B R O U T I N E     R E A D   C O N T R O L  F I L E
!---------------------------------------------------------------------
subroutine readControlFile(controlFileName, verbose)
use asgsio, only : openFileForRead
implicit none
character(len=1024), intent(in) :: controlFileName
logical, intent(in) :: verbose
character(len=80) :: wtimincLine
integer :: numFields ! number of fields read from wtiminc line
integer :: ios     ! i/o status from read 15; used in error messages
integer :: i
integer :: j
integer :: k
!
! initialization
lineNum = 1
!
CALL openFileForRead(15,trim(controlFileName))
!
read(15,fmt='(a80)',err=10,end=20,iostat=ios) rundes
linenum = linenum + 1
read(15,fmt='(a80)',err=10,end=20,iostat=ios) runid
linenum = linenum + 1
read(15,*,err=10,end=20,iostat=ios) nfover
linenum = linenum + 1
read(15,*,err=10,end=20,iostat=ios) nabout
linenum = linenum + 1
read(15,*,err=10,end=20,iostat=ios) nscreen
linenum = linenum + 1
read(15,*,err=10,end=20,iostat=ios) ihot
linenum = linenum + 1
read(15,*,err=10,end=20,iostat=ios) ics
linenum = linenum + 1
read(15,*,err=10,end=20,iostat=ios) im
linenum = linenum + 1
if ((im.eq.21).or.(im.eq.31)) then
   read(15,*,err=10,end=20,iostat=ios) iden
   linenum = linenum + 1
endif
read(15,*,err=10,end=20,iostat=ios) nolibf
linenum = linenum + 1
read(15,*,err=10,end=20,iostat=ios) nolifa
linenum = linenum + 1
read(15,*,err=10,end=20,iostat=ios) nolica
linenum = linenum + 1
read(15,*,err=10,end=20,iostat=ios) nolicat
linenum = linenum + 1
read(15,*,err=10,end=20,iostat=ios) nwp
linenum = linenum + 1
if (nwp.ne.0) then
   allocate(attrName(nwp))
   do i=1, nwp
      read(15,*,err=10,end=20,iostat=ios) attrName(i)
      linenum = linenum + 1
   end do
endif
read(15,*,err=10,end=20,iostat=ios) ncor
linenum = linenum + 1
read(15,*,err=10,end=20,iostat=ios) ntip
linenum = linenum + 1
read(15,fmt=*,err=10,end=20,iostat=ios) nws
linenum = linenum + 1
read(15,*,err=10,end=20,iostat=ios) nramp
linenum = linenum + 1
read(15,*,err=10,end=20,iostat=ios) g
linenum = linenum + 1
read(15,*,err=10,end=20,iostat=ios) tau0
linenum = linenum + 1
if ((tau0.le.-5.0).and.(tau0.ge.-5.99)) then
   read(15,*,err=10,end=20,iostat=ios) tau0fulldomainmin, tau0fulldomainmax
   linenum = linenum + 1
endif
read(15,*,err=10,end=20,iostat=ios) dtdp
linenum = linenum + 1
read(15,*,err=10,end=20,iostat=ios) statim
linenum = linenum + 1
read(15,*,err=10,end=20,iostat=ios) reftim
linenum = linenum + 1
!
! parse the digits of nws 
ncice = 0
nrs = 0
if (abs(nws).ge.1000) then
   ! ice is 1000s place
   ncice = int(abs(nws)/1000)
   nws = int(abs(nws)-ncice*1000)*int(nws/abs(nws))
endif
if (abs(nws).ge.100) then
   ! wave coupling is 100s place
   nrs=int(abs(nws/100))
   nws=int(abs(nws)-nrs*100)*int(nws/abs(nws))
endif
!
! now use the ncice, nrs, and nws values to parse the wtiminc line
if ((nws.ne.0).or.(nrs.ne.0).or.(ncice.ne.0)) then
   read(15,'(a80)',err=10,end=20,iostat=ios) wtimincLine
endif
numFields = 0
select case(abs(nws))
case(0,1) 
   ! do nothing
case(2,4,5,7,12,15,16)   
   read(wtimincLine,*,err=10,end=20,iostat=ios) wtiminc
   numFields = 1
case(3)
   read(wtimincLine,*,err=10,end=20,iostat=ios) &
      irefyr,irefmo,irefday,irefhr,irefmin,refsec
   read(wtimincLine,*,err=10,end=20,iostat=ios) &
      nwlat,nwlon,wlatmax,wlonmin,wlatinc,wloninc,wtiminc
   numfields = 7
case(6)
   read(wtimincLine,*,err=10,end=20,iostat=ios) &
      nwlat,nwlon,wlatmax,wlonmin,wlatinc,wloninc,wtiminc
   numFields = 7
case(8,19)
   read(wtimincLine,*,err=10,end=20,iostat=ios) &
      irefyr,irefmo,irefday,irefhr,stormnumber,bladj
   numFields = 6
case(10)
   nwlat=190
   nwlon=384
   read(wtimincLine,*,err=10,end=20,iostat=ios) wtiminc
   numFields = 1
case(11)
   nwlat=271
   nwlon=181
   wtiminc=10800.
case(29)
   read(wtimincLine,*,err=10,end=20,iostat=ios) &
      irefyr,irefmo,irefday,irefhr,stormnumber,bladj,purevortex,purebackground,wtiminc
   numFields = 9
case default
   write(6,'("ERROR: The value of nws is ",i6," but this is not a valid value.")') nws
   stop
end select
if ((nrs.ne.0).or.(ncice.ne.0)) then
   call parseWaveAndIceTimeIncrements(wtimincLine,numFields)
endif
linenum = linenum + 1
read(15,*,err=10,end=20,iostat=ios) rnday
linenum = linenum + 1
select case(nramp)
case(0,1)
   read(15,*,err=10,end=20,iostat=ios) dramp
case(2) 
   read(15,*,err=10,end=20,iostat=ios) dramp,drampextflux,fluxsettlingtime
case(3) 
   read(15,*,err=10,end=20,iostat=ios) dramp,drampextflux,fluxsettlingtime,drampintflux
case(4) 
   read(15,*,err=10,end=20,iostat=ios) dramp,drampextflux,fluxsettlingtime,drampintflux, &
              drampelev
case(5) 
   read(15,*,err=10,end=20,iostat=ios) dramp,drampextflux,fluxsettlingtime,drampintflux, &
             drampelev,dramptip
case(6) 
   read(15,*,err=10,end=20,iostat=ios) dramp,drampextflux,fluxsettlingtime,drampintflux, &
             drampelev,dramptip,drampmete
case(7) 
   read(15,*,err=10,end=20,iostat=ios) dramp,drampextflux,fluxsettlingtime,drampintflux, &
             drampelev,dramptip,drampmete,drampwrad
case(8) 
   read(15,*,err=10,end=20,iostat=ios) dramp,drampextflux,fluxsettlingtime,drampintflux, &
             drampelev,dramptip,drampmete,drampwrad,dunrampmete
case default
   write(6,'("ERROR: The value of nramp is ",i2," but this is not a valid value.")') nramp
   stop
end select
linenum = linenum + 1
read(15,*,err=10,end=20,iostat=ios) a00,b00,c00
linenum = linenum + 1
if (nolifa.eq.2) then
   read(15,*,err=10,end=20,iostat=ios) h0,nodedrymin,nodewetmin,velmin
else
   read(15,*,err=10,end=20,iostat=ios) h0
endif
linenum = linenum + 1
read(15,*,err=10,end=20,iostat=ios) slam0,sfea0
linenum = linenum + 1
select case(nolibf)
case(0)
   read(15,*,err=10,end=20,iostat=ios) tau
case(1)
   read(15,*,err=10,end=20,iostat=ios) cf
case(2)
   read(15,*,err=10,end=20,iostat=ios) cf,hbreak,ftheta,fgamma
case default
   write(6,'("error: the value of nolibf is ",i2," but this is not a valid value.")') &
      nolibf
   stop
end select
linenum = linenum + 1
if (im.eq.10) then
   read(15,*,err=10,end=20,iostat=ios) eslm,eslc
   linenum = linenum + 1
else
  read(15,*,err=10,end=20,iostat=ios) eslm
  linenum = linenum + 1
endif
read(15,*,err=10,end=20,iostat=ios) cori
linenum = linenum + 1
read(15,*,err=10,end=20,iostat=ios) ntif
linenum = linenum + 1
if (ntif.ne.0) then
   allocate(tipotag(ntif))
   allocate(tpk(ntif))
   allocate(amigt(ntif))
   allocate(etrf(ntif))
   allocate(fft(ntif))
   allocate(facet(ntif))
   do i=1,ntif
      read(15,'(a5)',err=10,end=20,iostat=ios)  tipotag(i)
      linenum = linenum + 1
      read(15,*,err=10,end=20,iostat=ios)  tpk(i),amigt(i),etrf(i),fft(i),facet(i)
      linenum = linenum + 1
   end do
endif
read(15,*,err=10,end=20,iostat=ios) nbfr
linenum = linenum + 1
if (nbfr.ne.0) then
   allocate(bountag(nbfr))
   allocate(amig(nbfr))
   allocate(ff(nbfr)) 
   allocate(face(nbfr))
   do i=1, nbfr
      read(15,'(a5)',err=10,end=20,iostat=ios) bountag(i)
      linenum = linenum + 1
      read(15,*,err=10,end=20,iostat=ios) amig(i),ff(i),face(i)
      linenum = linenum + 1
   end do
   allocate(nbfr_alpha(nbfr))
   allocate(emo(nbfr,neta))
   allocate(efa(nbfr,neta))
   do i=1, nbfr
      read(15,'(a10)',err=10,end=20,iostat=ios) nbfr_alpha(i)
      linenum = linenum + 1
      do j=1,neta
         read(15,*,err=10,end=20,iostat=ios) emo(i,j),efa(i,j)
         linenum = linenum + 1
      end do
   end do
endif
read(15,*,err=10,end=20,iostat=ios) anginn
linenum = linenum + 1
if (nfluxf.ne.0) then
   read(15,*,err=10,end=20,iostat=ios) nffr
   linenum = linenum + 1
endif
if ((nfluxf.ne.0).and.(nffr.gt.0)) then
   allocate(fbountag(nffr))
   allocate(famig(nffr))
   allocate(fff(nffr))
   allocate(fface(nffr))
   allocate(nffr_alpha(nffr))
   allocate(qnam(nffr,nvel))
   allocate(qnph(nffr,nvel))
   do i=1,nffr
      read(15,'(a5)',err=10,end=20,iostat=ios) fbountag(i)
      linenum = linenum + 1
      read(15,*,err=10,end=20,iostat=ios) famig(i),fff(i),fface(i)
      linenum = linenum + 1
   end do 
   do i=1,nffr
      read(15,'(a10)',err=10,end=20,iostat=ios) nffr_alpha(i)
      linenum = linenum + 1
      do j=1,nvel
         select case(lbcodei(j))
         case(2,12,22,52)
            read(15,*,err=10,end=20,iostat=ios) qnam(i,j), qnph(i,j)
            linenum = linenum + 1
         case default
            cycle
         end select
      end do
   end do
endif
read(15,*,err=10,end=20,iostat=ios) noute,toutse,toutfe,nspoole
linenum = linenum + 1
read(15,*,err=10,end=20,iostat=ios) nstae
linenum = linenum + 1
call readStations(nstae, xel, yel)
read(15,*,err=10,end=20,iostat=ios) noutv,toutsv,toutfv,nspoolv
linenum = linenum + 1
read(15,*,err=10,end=20,iostat=ios) nstav
linenum = linenum + 1
call readStations(nstav, xev, yev)
if (im.eq.10) then
   read(15,*,err=10,end=20,iostat=ios) noutc,toutsc,toutfc,nspoolc
   linenum = linenum + 1
   read(15,*,err=10,end=20,iostat=ios) nstac
   linenum = linenum + 1
   call readStations(nstac, xec, yec)
endif
if (nws.ne.0) then
   read(15,*,err=10,end=20,iostat=ios) noutm,toutsm,toutfm,nspoolm
   linenum = linenum + 1
   read(15,*,err=10,end=20,iostat=ios) nstam
   linenum = linenum + 1
   call readStations(nstam, xem, yem)
endif
read(15,*,err=10,end=20,iostat=ios) noutge,toutsge,toutfge,nspoolge
linenum = linenum + 1
read(15,*,err=10,end=20,iostat=ios) noutgv,toutsgv,toutfgv,nspoolgv
linenum = linenum + 1
if (im.eq.10) then
   read(15,*,err=10,end=20,iostat=ios) noutgc,toutsgc,toutfgc,nspoolgc
   linenum = linenum + 1
endif
if (nws.ne.0) then
   read(15,*,err=10,end=20,iostat=ios) noutgw,toutsgw,toutfgw,nspoolgw
   linenum = linenum + 1
endif
read(15,*,err=10,end=20,iostat=ios) nfreq
linenum = linenum + 1
if (nfreq.ne.0) then
   allocate(namefr(nfreq))
   allocate(hafreq(nfreq))
   allocate(haff(nfreq))
   allocate(haface(nfreq))
   do i=1,nfreq
      read(15,'(a10)',err=10,end=20,iostat=ios) namefr(i)
      linenum = linenum + 1
      read(15,*,err=10,end=20,iostat=ios) hafreq(i),haff(i),haface(i)
      linenum = linenum + 1
   end do
endif
read(15,*,err=10,end=20,iostat=ios) thas,thaf,nhainc,fmv
linenum = linenum + 1
read(15,*,err=10,end=20,iostat=ios) nhase,nhasv,nhage,nhagv
linenum = linenum + 1
read(15,*,err=10,end=20,iostat=ios) nhstar,nhsinc
linenum = linenum + 1
read(15,*,err=10,end=20,iostat=ios) ititer,isldia,convcr,itmax
linenum = linenum + 1
!
! read metadata if netcdf or xdmf format was specified in the output 
! identifiers
outputSpecifiers(1) = noute
outputSpecifiers(2) = noutv
outputSpecifiers(3) = noutc
outputSpecifiers(4) = noutm
outputSpecifiers(5) = noutge
outputSpecifiers(6) = noutgv
outputSpecifiers(7) = noutgc
outputSpecifiers(8) = noutgw
outputSpecifiers(9) = nhstar
! check to see if we need to read in metadata
readMetaData = .false.
outputSpec: do i=1,9 
   metadataReq: do j=1,9 
      if (outputSpecifiers(i).eq.metadataRequired(j)) then
         readMetaData = .true.
         exit outputSpec 
      endif
   end do metadataReq
end do outputSpec
! now read metadata if any output format indicated it
if (readMetaData.eqv..true.) then 
   read(15,'(a80)',err=10,end=20,iostat=ios) title
   linenum = linenum + 1
   read(15,'(a80)',err=10,end=20,iostat=ios) institution
   linenum = linenum + 1
   read(15,'(a80)',err=10,end=20,iostat=ios) source 
   linenum = linenum + 1
   read(15,'(a80)',err=10,end=20,iostat=ios) history
   linenum = linenum + 1
   read(15,'(a80)',err=10,end=20,iostat=ios) references
   linenum = linenum + 1
   read(15,'(a80)',err=10,end=20,iostat=ios) comments
   linenum = linenum + 1
   read(15,'(a80)',err=10,end=20,iostat=ios) host
   linenum = linenum + 1
   read(15,'(a80)',err=10,end=20,iostat=ios) convention
   linenum = linenum + 1
   read(15,'(a80)',err=10,end=20,iostat=ios) contact
   linenum = linenum + 1
   read(15,'(a80)',err=10,end=20,iostat=ios) base_date
endif 
!
if (verbose.eqv..true.) then
   call echoControlFile(6)
endif
close(15)
return
10 write(6,'("ERROR: Reading line ",I6," gave the following error code: ",I6,".")') lineNum, ios
call echoControlFile(6)
close(15)
stop
20 write(6,'("ERROR: Reached premature end of file on line ",I6,".")') lineNum
call echoControlFile(6)
close(15)
stop
!---------------------------------------------------------------------
end subroutine readControlFile
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
! S U B R O U T I N E    R E A D   S T A T I O N S
!---------------------------------------------------------------------
subroutine readStations(nsta, stax, stay)
implicit none
integer, intent(in) :: nsta
real(sz), allocatable, intent(out) :: stax(:)
real(sz), allocatable, intent(out) :: stay(:)
integer :: i
integer :: ios
if (nsta.ne.0) then
   allocate(stax(nsta))
   allocate(stay(nsta))
   do i=1,nsta
      read(15,*,err=10,end=20,iostat=ios) stax(i), stay(i)
      linenum = linenum + 1
   end do
endif
return
10 write(6,'("ERROR: Reading line ",I3," gave the following error code: ",I3,".")') lineNum, ios
call echoControlFile(6)
stop
20 write(6,'("ERROR: Reached premature end of file on line ",I3,".")') lineNum
call echoControlFile(6)
stop
!---------------------------------------------------------------------
end subroutine readStations
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!                      S U B R O U T I N E   
!   P A R S E  W A V E  A N D  I C E  T I M E  I N C R E M E N T S
!---------------------------------------------------------------------
! jgf: This subroutine is needed because of the tricky nature of the 
! WTIMINC line in the fort.15 file. The number of data fields on this
! line is dependent on the value of NWS, but the exact formatting of 
! these data fields (width of fields, number of decimal places, etc)
! is unknown. These uncertainties make this line hard to parse.
!---------------------------------------------------------------------
subroutine parseWaveAndIceTimeIncrements(wtimincLine,numFields)
implicit none
character(len=80), intent(in) :: wtimincLine ! line of data from fort.15
integer, intent(inout) :: numFields ! number of fields that have already been read 
!
character(len=80) :: rstimincField    ! character representation of rstiminc val
character(len=80) :: cicetimeincField ! character representation of cicetimeinc val
logical :: foundField       ! .true. if a new data field was found
integer :: fieldStarts(100) ! array of starting positions of data fields
integer :: fieldEnds(100)   ! array of ending positions of data fields
integer :: i  ! index counter for the wtiminc character array
integer :: j  ! index into the array of starting positions for each field
integer :: k  ! index into the array of ending positions for each field
integer :: ios
integer :: istart
integer :: iend
!
rstimincField(:) = ' '
cicetimeincField(:) = ' ' 
foundField = .false.
j=1
k=1
do i=1,len(wtimincLine)
   if (foundField.eqv..false.) then
      ! look for start of a new field
      if (wtimincLine(i:i).ne.' ') then
         ! found the start of a new field
         fieldStarts(j) = i
         j = j + 1 
         foundField = .true.
      endif
   else
      ! we have found a field, look for the end of it
      if (wtimincLine(i:i).eq.' ') then
         ! found the end of this field
         fieldEnds(k) = i
         k = k + 1
         foundField = .false.
      endif
   endif
end do
if (nrs.gt.0) then
   istart = fieldStarts(numFields + 1)
   iend = fieldEnds(numFields + 1)
   rstimincField = wtimincLine(istart:iend)
   read(rstimincField,*,err=10,end=20,iostat=ios) rstiminc
   numFields = numFields + 1 
endif
if (ncice.gt.0) then
   istart = fieldStarts(numFields + 1)
   iend = fieldEnds(numFields + 1)
   cicetimeincField = wtimincLine(istart:iend)
   read(cicetimeincField,*,err=10,end=20,iostat=ios) cice_timinc
   numFields = numFields + 1 
endif
return
10 write(6,'("ERROR: Reading line ",I3," gave the following error code: ",I3,".")') lineNum, ios
call echoControlFile(6)
stop
20 write(6,'("ERROR: Reached premature end of file on line ",I3,".")') lineNum
call echoControlFile(6)
stop
!---------------------------------------------------------------------
end subroutine parseWaveAndIceTimeIncrements
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
! S U B R O U T I N E    E C H O    C O N T R O L   F I L E 
!---------------------------------------------------------------------
subroutine echoControlFile(echoUnit)
implicit none
integer, intent(in) :: echoUnit ! fortran i/o unit number for output
character(len=80) :: wtimincLine
character(len=80) :: rstimincLine
character(len=80) :: cicetimincLine
character(len=80) :: wtimincComment
integer :: controlnws ! nws as encoded with ncice and nrs parameters
integer :: i 
integer :: j
!
! initialization
echoLine = 1
!
write(6,*) 'INFO: Echoing control file after reading ',lineNum,' lines.'
write(echoUnit,fmt='(A,10x,"! rundes")') trim(rundes)
echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return 
write(echoUnit,fmt='(A,10x,"! runid")') trim(runid)
echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
write(echoUnit,fmt='(i0,10x,"! nfover")') nfover
echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
write(echoUnit,fmt='(i0,10x,"! nabout")') nabout
echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
write(echoUnit,fmt='(i0, 10x, "! nscreen")') nscreen
echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
write(echoUnit,fmt='(i0,10x,"! ihot")') ihot
echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
write(echoUnit,fmt='(i0, 10x, "! ics")') ics
echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
write(echoUnit,fmt='(i0, 10x, "! im")') im
echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
if ((im.eq.21).or.(im.eq.31)) then
   write(echoUnit,fmt='(i0, 10x, "! iden")') iden
   echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
endif
write(echoUnit,fmt='(i0, 10x, "! nolibf")') nolibf
echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
write(echoUnit,fmt='(i0, 10x, "! nolifa")') nolifa
echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
write(echoUnit,fmt='(i0, 10x, "! nolica")') nolica
echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
write(echoUnit,fmt='(i0, 10x, "! nolicat")') nolicat
echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
write(echoUnit,fmt='(i0, 10x, "! nwp")') nwp
echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
if (nwp.ne.0) then
   do i=1, nwp
      write(echoUnit,fmt='(a, 10x, "! attrname")') trim(attrName(i))
      echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
   end do
endif
write(echoUnit,fmt='(i0, 10x, "! ncor")') ncor
echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
write(echoUnit,fmt='(i0, 10x, "! ntip")') ntip
echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
controlnws = 1000 * ncice + 100 * nrs + abs(nws)
if (nws.ne.0) then
   controlnws = controlnws * abs(nws)/nws ! to get the sign right
endif
write(echoUnit,fmt='(i0, 10x, "! nws")') controlnws
echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
write(echoUnit,fmt='(i0, 10x, "! nramp")') nramp
echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
write(echoUnit,fmt='(f10.5, 10x, "! g")') g
echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
write(echoUnit,fmt='(f6.3, 10x, "! tau0")') tau0
echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
if ((tau0.le.-5.0).and.(tau0.ge.-5.99)) then
   write(echoUnit,fmt='(f6.3, f6.3, 10x, "! tau0fulldomainmin tau0fulldomainmax")') tau0fulldomainmin, tau0fulldomainmax
   echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
endif
write(echoUnit,fmt='(f15.7, 10x, "1 dtdp")') dtdp
echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
write(echoUnit,fmt='(f15.7, 10x, "! statim")') statim
echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
write(echoUnit,fmt='(f15.7, 10x, "! reftim")') reftim
echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
wtimincLine = ''
wtimincComment = ''
select case(abs(nws))
case(0,1) 
   ! do nothing
case(2,4,5,7,10,11,12,15,16)   
   write(wtimincLine,*) wtiminc 
   wtimincComment = ' ! wtiminc'
case(3)   
   write(wtimincLine,*) irefyr,irefmo,irefday,irefhr,irefmin,refsec 
   wtimincLine = trim(wtimincLine) // ' ! irefyr,irefmo,irefday,irefhr,irefmin,refsec'
   write(echoUnit,'(a)') trim(wtimincLine)
   echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
   write(wtimincLine,*) nwlat,nwlon,wlatmax,wlonmin,wlatinc,wloninc,wtiminc
   wtimincComment = ' ! nwlat,nwlon,wlatmax,wlonmin,wlatinc,wloninc,wtiminc'
case(6)
   write(wtimincLine,*) nwlat,nwlon,wlatmax,wlonmin,wlatinc,wloninc,wtiminc
   wtimincComment = ' ! nwlat,nwlon,wlatmax,wlonmin,wlatinc,wloninc,wtiminc'
case(8,19)
   write(wtimincLine,*) irefyr,irefmo,irefday,irefhr,stormnumber,bladj
   wtimincComment = ' ! irefyr,irefmo,irefday,irefhr,stormnumber,bladj'
case(29)
   write(wtimincLine,*) irefyr,irefmo,irefday,irefhr,stormnumber,bladj,wtiminc, purevortex, purebackground
   wtimincComment = ' ! irefyr,irefmo,irefday,irefhr,stormnumber,bladj,wtiminc, purevortex, purebackground'
case default
   write(echoUnit,'("ERROR: The value of nws is ",i6," but this is not a valid value.")') nws
   stop
end select
if (nrs.gt.0) then
    write(rstimincLine,*) rstiminc
    wtimincLine = trim(wtimincLine) // ' ' // trim(rstimincLine)
    wtimincComment = trim(wtimincComment) // ',rstiminc'
endif
if (ncice.gt.0) then
    write(cicetimincLine,*) cice_timinc
    wtimincLine = trim(wtimincLine) // ' ' // trim(cicetimincLine)
    wtimincComment = trim(wtimincComment) // ',cice_timinc'
endif
if ((nws.ne.0).or.(nrs.ne.0).or.(ncice.ne.0)) then
   write(echoUnit,'(a)') trim(wtimincLine) // trim(wtimincComment)
   echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
endif
write(echoUnit,fmt='(e20.10, 10x, "! rnday")') rnday
echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
!
select case(nramp)
case(0,1)
   write(echoUnit,fmt='(9f6.3,10x,"! dramp")') dramp
case(2) 
   write(echoUnit,fmt='(9f6.3,10x,"! dramp,drampextflux,fluxsettlingtime")') &
      dramp,drampextflux,fluxsettlingtime
case(3) 
   write(echoUnit,fmt='(9f6.3,10x, &
   "! dramp,drampextflux,fluxsettlingtime,drampintflux")')  &
      dramp,drampextflux,fluxsettlingtime,drampintflux
case(4) 
   write(echoUnit,fmt='(9f6.3,10x, &
   "! dramp,drampextflux,fluxsettlingtime,drampintflux, drampelev")') &
      dramp,drampextflux,fluxsettlingtime,drampintflux, drampelev
case(5) 
   write(echoUnit,fmt='(9f6.3,10x, &
   "! dramp,drampextflux,fluxsettlingtime,drampintflux,drampelev,dramptip")') &
      dramp,drampextflux,fluxsettlingtime,drampintflux,drampelev,dramptip
case(6) 
   write(echoUnit,fmt='(9f6.3,10x, &
   "! dramp,drampextflux,fluxsettlingtime,drampintflux,drampelev,dramptip,drampmete")') &
       dramp,drampextflux,fluxsettlingtime,drampintflux,drampelev,dramptip,drampmete
case(7) 
   write(echoUnit,fmt='(9f6.3,10x, &
   "! dramp,drampextflux,fluxsettlingtime,drampintflux, &
             drampelev,dramptip,drampmete,drampwrad")') &
      dramp,drampextflux,fluxsettlingtime,drampintflux, &
             drampelev,dramptip,drampmete,drampwrad
case(8) 
   write(echoUnit,fmt='(9f6.3,10x, &
   "! dramp,drampextflux,fluxsettlingtime,drampintflux, &
             drampelev,dramptip,drampmete,drampwrad,drampunmete")') &
         dramp,drampextflux,fluxsettlingtime,drampintflux, &
             drampelev,dramptip,drampmete,drampwrad,dunrampmete
case default
   write(echoUnit,'("ERROR: The value of nramp is ",i2," but this is not a valid value.")') nramp
   stop
end select
echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
write(echoUnit,fmt='(3f6.3, 10x, "! a00, b00, c00")') a00,b00,c00
echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
if (nolifa.eq.2) then
   write(echoUnit,fmt='(f15.7, i2, i2, f6.3, 10x, "! h0, int, int, velmin")') &
      h0,nodedrymin,nodewetmin,velmin
else
   write(echoUnit,fmt='(f15.7,10x,"! h0")') h0
endif
echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
write(echoUnit,fmt='(2f15.7,10x,"! slam0, sfea0")') slam0,sfea0
echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
select case(nolibf)
case(0)
   write(echoUnit,fmt='(f15.7,10x,"! tau")') tau
case(1)
   write(echoUnit,fmt='(f15.7,10x,"! cf")') cf
case(2)
   write(echoUnit,fmt='(4f15.7,10x,"! cf, hbreak, ftheta, fgamma")')  &
      cf,hbreak,ftheta,fgamma
case default
   write(echoUnit,'("ERROR: the value of nolibf is ",i2," but this is not a valid value.")') &
      nolibf
   stop
end select
echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
if (im.eq.10) then
   write(echoUnit,fmt='(2f15.7,10x,"! eslm, eslc")') eslm,eslc
   echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
else
  write(echoUnit,fmt='(f15.7,10x,"! eslm")') eslm
  echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
endif
write(echoUnit,fmt='(f15.7,10x,"! cori")') cori
echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
write(echoUnit,fmt='(i0,10x, "! ntif")') ntif
echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
do i=1,ntif
   write(echoUnit,'(a5,10x,"! tipotag")')  trim(tipotag(i))
   echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
   write(echoUnit,fmt='(5f15.7,10x,"! tpk(i),amigt(i),etrf(i),fft(i),facet(i)")') &
      tpk(i),amigt(i),etrf(i),fft(i),facet(i)
   echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
end do
write(echoUnit,fmt='(i0,10x,"! nbfr")') nbfr
echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
do i=1, nbfr
   write(echoUnit,'(a5,10x,"! bountag")') trim(bountag(i))
   echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
   write(echoUnit,fmt='(3f15.7,10x,"! amig(i),ff(i),face(i)")') amig(i),ff(i),face(i)
   echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
end do
do i=1, nbfr
   write(echoUnit,'(a10,10x,"! alpha(i)")') nbfr_alpha(i)
   echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
   do j=1,neta
      write(echoUnit,fmt='(2f15.7,10x,"! emo(i,j), efa(i,j)")') emo(i,j),efa(i,j)
      echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
   end do
end do
write(echoUnit,fmt='(f15.7,10x,"! anginn")') anginn
echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
if ((nfluxf.ne.0).and.(nffr.ne.0)) then
   write(echoUnit,fmt='(i0,10x,"! nffr")') nffr
   echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
   if (nffr.gt.0) then
      do i=1,nffr
         write(echoUnit,'(a5,10x,"! fbountag(i)")') fbountag(i)
         echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
         write(echoUnit,fmt='(3f6.3,10x,"! famig(i),fff(i),fface(i)")') &
             famig(i),fff(i),fface(i)
         echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
      end do 
      do i=1,nffr
         write(echoUnit,'(a10,10x,"! nffr_alpha(i)")') nffr_alpha(i)
         echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
         do j=1,nvel
            select case(lbcodei(j))
            case(2,12,22,52)
               write(echoUnit,fmt='(2f6.3,10x,"! qnam(i,j), qnph(i,j)")') qnam(i,j), qnph(i,j)
               echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
            case default
               cycle
            end select
         end do
      end do
   endif
endif
write(echoUnit,fmt='(i0,1x,2f6.3,1x,i0,10x,"! noute,toutse,toutfe,nspoole")') noute,toutse,toutfe,nspoole
echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
write(echoUnit,fmt='(i0,10x,"! nstae")') nstae
echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
call echoStations(nstae, xel, yel, echoUnit)
write(echoUnit,fmt='(i0,1x,2f15.7,1x,i0,10x,"! noutv,toutsv,toutfv,nspoolv")') noutv,toutsv,toutfv,nspoolv
echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
write(echoUnit,fmt='(i0,10x,"! nstav")') nstav
echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
call echoStations(nstav, xev, yev, echoUnit)
if (im.eq.10) then
   write(echoUnit,fmt='(i0,1x,2f15.7,1x,i0,10x,"! noutc,toutsc,toutfc,nspoolc")') noutc,toutsc,toutfc,nspoolc
   echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
   write(echoUnit,fmt='(i0,10x,"! nstac")') nstac
   echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
   call echoStations(nstac, xec, yec, echoUnit)
endif
if (nws.ne.0) then
   write(echoUnit,fmt='(i0,1x,2f15.7,1x,i0,10x,"! noutm,toutsm,toutfm,nspoolm")') noutm,toutsm,toutfm,nspoolm
   echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
   write(echoUnit,fmt='(i0,10x,"! nstam")') nstam
   echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
   call echoStations(nstam, xem, yem, echoUnit)
endif
write(echoUnit,fmt='(i0,1x,2f15.7,1x,i0,10x,"! noutge,toutsge,toutfge,nspoolge")') noutge,toutsge,toutfge,nspoolge
echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
write(echoUnit,fmt='(i0,1x,2f15.7,1x,i0,10x,"! noutgv,toutsgv,toutfgv,nspoolgv")') noutgv,toutsgv,toutfgv,nspoolgv
echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
if (im.eq.10) then
   write(echoUnit,fmt='(i0,1x,2f15.7,1x,i0,10x,"! noutgc,toutsgc,toutfgc,nspoolgc")') noutgc,toutsgc,toutfgc,nspoolgc
   echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
endif
if (nws.ne.0) then
   write(echoUnit,fmt='(i0,1x,2f15.7,1x,i0,10x,"! noutgw,toutsgw,toutfgw,nspoolgw")') noutgw,toutsgw,toutfgw,nspoolgw
   echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
endif
write(echoUnit,fmt='(i0,10x,"! nfreq")') nfreq
echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
if (nfreq.ne.0) then
   do i=1,nfreq
      write(echoUnit,'(a10,10x,"! namefr(i)")') namefr(i)
      echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
      write(echoUnit,fmt='(3f15.7,10x,"! hafreq(i),haff(i),haface(i)")') hafreq(i),haff(i),haface(i)
      echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
   end do
endif
write(echoUnit,fmt='(2f15.7,1x,i0,1x,f15.7,10x,"! thas,thaf,nhainc,fmv")') thas,thaf,nhainc,fmv
echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
write(echoUnit,fmt='(i0,1x,i0,1x,i0,1x,i0,10x,"! nhase,nhasv,nhage,nhagv")') nhase,nhasv,nhage,nhagv
echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
write(echoUnit,fmt='(i0,1x,i0,10x,"! nhstar,nhsinc")') nhstar,nhsinc
echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
write(echoUnit,fmt='(i0,1x,i0,f15.7,1x,i0,10x,"! ititer,isldia,convcr,itmax")') ititer,isldia,convcr,itmax
echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
!
! read metadata if netcdf or xdmf format was specified in the output 
! identifiers
outputSpecifiers(1) = noute
outputSpecifiers(2) = noutv
outputSpecifiers(3) = noutc
outputSpecifiers(4) = noutm
outputSpecifiers(5) = noutge
outputSpecifiers(6) = noutgv
outputSpecifiers(7) = noutgc
outputSpecifiers(8) = noutgw
outputSpecifiers(9) = nhstar
! now write metadata if any output format indicated it
if (readMetaData.eqv..true.) then 
   write(echoUnit,fmt='(a)') trim(title)
   echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
   write(echoUnit,fmt='(a)') trim(institution)
   echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
   write(echoUnit,fmt='(a)') trim(source)
   echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
   write(echoUnit,fmt='(a)') trim(history)
   echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
   write(echoUnit,fmt='(a)') trim(references)
   echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
   write(echoUnit,fmt='(a)') trim(comments)
   echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
   write(echoUnit,fmt='(a)') trim(host)
   echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
   write(echoUnit,fmt='(a)') trim(convention)
   echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
   write(echoUnit,fmt='(a)') trim(contact)
   echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
   write(echoUnit,fmt='(a)') trim(base_date)
endif 
!---------------------------------------------------------------------
end subroutine echoControlFile
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
! S U B R O U T I N E    E C H O   S T A T I O N S
!---------------------------------------------------------------------
subroutine echoStations(nsta, stax, stay, echoUnit)
implicit none
integer, intent(in) :: nsta
real(sz), intent(in) :: stax(:)
real(sz), intent(in) :: stay(:)
integer, intent(in) :: echoUnit ! i/o unit for echoing data
integer :: i
if (nsta.ne.0) then
   do i=1,nsta
      write(echoUnit,fmt='(2f15.7,10x,"! lon lat")') stax(i), stay(i)
      echoLine = echoLine + 1 ; if (echoLine.gt.lineNum) return
   end do
endif
!---------------------------------------------------------------------
end subroutine echoStations
!---------------------------------------------------------------------

!---------------------------------------------------------------------
! S U B R O U T I N E    W R I T E    C O N T R O L   X D M F
!---------------------------------------------------------------------
subroutine writeControlXDMF(xdmfFortranObj)
implicit none
include 'Xdmf.f'
integer*8, intent(in) :: xdmfFortranObj ! XDMF object to receive data
character(len=80) :: wtimincLine
character(len=80) :: rstimincLine
character(len=80) :: cicetimincLine
character(len=80) :: info
character(len=80) :: varname
character(len=80) :: stationType
character(len=80) :: wtimincComment
character(len=80) :: drampComment
integer :: informationID
integer :: controlnws ! nws as encoded with ncice and nrs parameters
integer :: i 
integer :: j
!
write(6,'(a)') 'INFO: Writing control parameters to XDMF.'
informationID = XdmfAddInformation(xdmfFortranObj, 'rundes'//char(0), &
   trim(adjustl(rundes))//char(0))
informationID = XdmfAddInformation(xdmfFortranObj, 'runid'//char(0), &
   trim(adjustl(runid))//char(0))
write(info,fmt='(i0)') nfover
informationID = XdmfAddInformation(xdmfFortranObj, 'nfover'//char(0), & 
   trim(adjustl(info))//char(0))

write(info,fmt='(i0)') nabout
informationID = XdmfAddInformation(xdmfFortranObj, 'nabout'//char(0), & 
   trim(adjustl(info))//char(0))
write(info,fmt='(i0)') nscreen
informationID = XdmfAddInformation(xdmfFortranObj, 'nscreen'//char(0), & 
   trim(adjustl(info))//char(0))
write(info,fmt='(i0)') ihot
informationID = XdmfAddInformation(xdmfFortranObj, 'ihot'//char(0), & 
   trim(adjustl(info))//char(0))
write(info,fmt='(i0)') ics
informationID = XdmfAddInformation(xdmfFortranObj, 'ics'//char(0), & 
   trim(adjustl(info))//char(0))
write(info,fmt='(i0)') im
informationID = XdmfAddInformation(xdmfFortranObj, 'im'//char(0), & 
   trim(adjustl(info))//char(0))
if ((im.eq.21).or.(im.eq.31)) then
   write(info,fmt='(i0)') iden
   informationID = XdmfAddInformation(xdmfFortranObj, 'iden'//char(0), & 
      trim(adjustl(info))//char(0))
endif
write(info,fmt='(i0)') nolibf
informationID = XdmfAddInformation(xdmfFortranObj, 'nolibf'//char(0), & 
   trim(adjustl(info))//char(0))
write(info,fmt='(i0)') nolifa
informationID = XdmfAddInformation(xdmfFortranObj, 'nolifa'//char(0), & 
   trim(adjustl(info))//char(0))
write(info,fmt='(i0)') nolica
informationID = XdmfAddInformation(xdmfFortranObj, 'nolica'//char(0), & 
   trim(adjustl(info))//char(0))
write(info,fmt='(i0)') nolicat
informationID = XdmfAddInformation(xdmfFortranObj, 'nolicat'//char(0), & 
   trim(adjustl(info))//char(0))
write(info,fmt='(i0)') nwp
informationID = XdmfAddInformation(xdmfFortranObj, 'nwp'//char(0), & 
   trim(adjustl(info))//char(0))

if (nwp.ne.0) then
   do i=1, nwp
      write(info,'("attribute_name(",i0,")")') i
      informationID = XdmfAddInformation(xdmfFortranObj, & 
         trim(adjustl(info))//char(0), trim(adjustl(attrName(i)))//char(0))
   end do
endif
write(info,fmt='(i0)') ncor
informationID = XdmfAddInformation(xdmfFortranObj, 'ncor'//char(0), & 
   trim(adjustl(info))//char(0))
write(info,fmt='(i0)') ntip
informationID = XdmfAddInformation(xdmfFortranObj, 'ntip'//char(0), & 
   trim(adjustl(info))//char(0))
controlnws = 1000 * ncice + 100 * nrs + abs(nws)
if (nws.ne.0) then
   controlnws = controlnws * abs(nws)/nws ! to get the sign right
endif
write(info,fmt='(i0)') controlnws
informationID = XdmfAddInformation(xdmfFortranObj, 'nws'//char(0), & 
   trim(adjustl(info))//char(0))
write(info,fmt='(i0)') nramp
informationID = XdmfAddInformation(xdmfFortranObj, 'nramp'//char(0), & 
   trim(adjustl(info))//char(0))
write(info,fmt='(f10.5)') g
informationID = XdmfAddInformation(xdmfFortranObj, 'g'//char(0), & 
   trim(adjustl(info))//char(0))
write(info,fmt='(f6.3)') tau0
informationID = XdmfAddInformation(xdmfFortranObj, 'tau0'//char(0), & 
   trim(adjustl(info))//char(0))
if ((tau0.le.-5.0).and.(tau0.ge.-5.99)) then
   write(info,fmt='(f6.3,f6.3)') tau0fulldomainmin,tau0fulldomainmax
   informationID = XdmfAddInformation(xdmfFortranObj, &
      'tau0fulldomainmin,tau0fulldomainmax'//char(0), & 
      trim(adjustl(info))//char(0))
endif
write(info,fmt='(f15.7)') dtdp
informationID = XdmfAddInformation(xdmfFortranObj,'dtdp'//char(0), & 
   trim(adjustl(info))//char(0))
write(info,fmt='(f15.7)') statim
informationID = XdmfAddInformation(xdmfFortranObj,'statim'//char(0), & 
   trim(adjustl(info))//char(0))
write(info,fmt='(f15.7)') reftim
informationID = XdmfAddInformation(xdmfFortranObj,'reftim'//char(0), & 
   trim(adjustl(info))//char(0))

wtimincLine = ''
wtimincComment = ''
select case(abs(nws))
case(0,1) 
   ! do nothing
case(2,4,5,7,10,11,12,15,16)   
   wtimincComment = 'wtiminc'
case(3)   
   wtimincComment = 'irefyr,irefmo,irefday,irefhr,irefmin,refsec'
   write(wtimincLine,*) irefyr,irefmo,irefday,irefhr,irefmin,refsec 
   informationID = XdmfAddInformation(xdmfFortranObj,trim(wtimincComment)//char(0), & 
      trim(adjustl(wtimincLine))//char(0))
   wtimincComment = 'nwlat,nwlon,wlatmax,wlonmin,wlatinc,wloninc,wtiminc'
   write(wtimincLine,*) nwlat,nwlon,wlatmax,wlonmin,wlatinc,wloninc,wtiminc
case(6)
   wtimincComment = 'nwlat,nwlon,wlatmax,wlonmin,wlatinc,wloninc,wtiminc'
   write(wtimincLine,*) nwlat,nwlon,wlatmax,wlonmin,wlatinc,wloninc,wtiminc
case(8,19)
   wtimincComment = 'irefyr,irefmo,irefday,irefhr,stormnumber,bladj'
   write(wtimincLine,*) irefyr,irefmo,irefday,irefhr,stormnumber,bladj
case(29)
   wtimincComment = 'irefyr,irefmo,irefday,irefhr,stormnumber,bladj,wtiminc,purevortex,purebackground'
   write(wtimincLine,*) irefyr,irefmo,irefday,irefhr,stormnumber,bladj,wtiminc,purevortex, purebackground
case default
   write(info,'("ERROR: The value of nws is ",i6," but this is not a valid value.")') nws
   stop
end select
if (nrs.gt.0) then
    wtimincComment = trim(wtimincComment) // ',rstiminc'
    write(rstimincLine,*) rstiminc
    wtimincLine = trim(wtimincLine) // ' ' // trim(rstimincLine)
endif
if (ncice.gt.0) then
    wtimincComment = trim(wtimincComment) // ',cice_timinc'
    write(cicetimincLine,*) cice_timinc
    wtimincLine = trim(wtimincLine) // ' ' // trim(cicetimincLine)
endif
if ((nws.ne.0).or.(nrs.ne.0).or.(ncice.ne.0)) then
   informationID = XdmfAddInformation(xdmfFortranObj,trim(wtimincComment)//char(0), & 
      trim(adjustl(wtimincLine))//char(0))
endif


write(info,fmt='(e20.10)') rnday
   informationID = XdmfAddInformation(xdmfFortranObj,'rnday'//char(0), & 
      trim(adjustl(info))//char(0))
!
select case(nramp)
case(0,1)
   write(info,fmt='(9f6.3)') dramp  
   drampComment = 'dramp'
case(2) 
   write(info,fmt='(9f6.3)') &
      dramp,drampextflux,fluxsettlingtime
   drampComment = 'dramp,drampextflux,fluxsettlingtime'
case(3) 
   write(info,fmt='(9f6.3)')  &
      dramp,drampextflux,fluxsettlingtime,drampintflux
   drampComment = 'dramp,drampextflux,fluxsettlingtime,drampintflux'
case(4) 
   write(info,fmt='(9f6.3)') &
      dramp,drampextflux,fluxsettlingtime,drampintflux, drampelev
   drampComment = 'dramp,drampextflux,fluxsettlingtime,drampintflux, drampelev'
case(5) 
   write(info,fmt='(9f6.3)') &
      dramp,drampextflux,fluxsettlingtime,drampintflux,drampelev,dramptip
   drampComment = 'dramp,drampextflux,fluxsettlingtime,drampintflux,drampelev,dramptip'
case(6) 
   write(info,fmt='(9f6.3)') &
       dramp,drampextflux,fluxsettlingtime,drampintflux,drampelev,dramptip,drampmete
   drampComment = 'dramp,drampextflux,fluxsettlingtime,drampintflux,drampelev,dramptip,drampmete'
case(7) 
   write(info,fmt='(9f6.3)') &
      dramp,drampextflux,fluxsettlingtime,drampintflux, &
             drampelev,dramptip,drampmete,drampwrad
   drampComment = 'dramp,drampextflux,fluxsettlingtime,drampintflux,' // &
             'drampelev,dramptip,drampmete,drampwrad'
case(8) 
   write(info,fmt='(9f6.3)') &
         dramp,drampextflux,fluxsettlingtime,drampintflux, &
             drampelev,dramptip,drampmete,drampwrad,dunrampmete
   drampcomment = 'dramp,drampextflux,fluxsettlingtime,drampintflux,' // &
             'drampelev,dramptip,drampmete,drampwrad,dunrampmete'
case default
   write(info,'("ERROR: The value of nramp is ",i2," but this is not a valid value.")') nramp
   stop
end select
informationID = XdmfAddInformation(xdmfFortranObj,trim(drampComment)//char(0), & 
   trim(adjustl(info))//char(0))
!
write(info,fmt='(3f6.3)') a00,b00,c00
informationID = XdmfAddInformation(xdmfFortranObj,'a00,b00,c00'//char(0), & 
   trim(adjustl(info))//char(0))
if (nolifa.eq.2) then
   write(info,fmt='(f15.7, i2, i2, f6.3)') &
      h0,nodedrymin,nodewetmin,velmin
   informationID = XdmfAddInformation(xdmfFortranObj, &
      'h0,nodedrymin,nodewetmin,velmin'//char(0), & 
      trim(adjustl(info))//char(0))
else
   write(info,fmt='(f15.7)') h0
   informationID = XdmfAddInformation(xdmfFortranObj,'h0'//char(0), & 
      trim(adjustl(info))//char(0))
endif
write(info,fmt='(2f15.7)') slam0,sfea0
   informationID = XdmfAddInformation(xdmfFortranObj,'slam0,sfea0'//char(0), & 
      trim(adjustl(info))//char(0))
select case(nolibf)
case(0)
   write(info,fmt='(f15.7)') tau
   informationID = XdmfAddInformation(xdmfFortranObj,'tau'//char(0), & 
      trim(adjustl(info))//char(0))
case(1)
   write(info,fmt='(f15.7)') cf
   informationID = XdmfAddInformation(xdmfFortranObj,'cf'//char(0), & 
      trim(adjustl(info))//char(0))
case(2)
   write(info,fmt='(4f15.7)') cf,hbreak,ftheta,fgamma
   informationID = XdmfAddInformation(xdmfFortranObj,'cf,hbreak,ftheta,fgamma'//char(0), & 
      trim(adjustl(info))//char(0))
case default
   write(info,'("ERROR: the value of nolibf is ",i2," but this is not a valid value.")') &
      nolibf
   stop
end select
if (im.eq.10) then
   write(info,fmt='(2f15.7)') eslm,eslc
   informationID = XdmfAddInformation(xdmfFortranObj,'eslm,eslc'//char(0), & 
      trim(adjustl(info))//char(0))
else
  write(info,fmt='(f15.7)') eslm
  informationID = XdmfAddInformation(xdmfFortranObj,'eslm'//char(0), & 
      trim(adjustl(info))//char(0))
endif
write(info,fmt='(f15.7)') cori
informationID = XdmfAddInformation(xdmfFortranObj,'cori'//char(0), & 
   trim(adjustl(info))//char(0))
write(info,fmt='(i0)') ntif
informationID = XdmfAddInformation(xdmfFortranObj,'ntif'//char(0), & 
   trim(adjustl(info))//char(0))
do i=1,ntif
   write(info,'(a)')  trim(tipotag(i))
   write(varname,'("tipotag(",i0,")")') i
   informationID = XdmfAddInformation(xdmfFortranObj, & 
      trim(adjustl(varname))//char(0), trim(adjustl(info))//char(0))     
   write(info,fmt='(5f15.7)') tpk(i),amigt(i),etrf(i),fft(i),facet(i)
   write(varname,'("tpk(",i0,"),amigt(",i0,"),etrf(",i0,"),fft(",i0,"),facet(",i0,")")' ) &
       i, i, i, i, i
   informationID = XdmfAddInformation(xdmfFortranObj, & 
      trim(adjustl(varname))//char(0), trim(adjustl(info))//char(0))
end do

write(info,fmt='(i0)') nbfr
informationID = XdmfAddInformation(xdmfFortranObj,'nbfr'//char(0), & 
   trim(adjustl(info))//char(0))
do i=1, nbfr
   write(info,'(a)') trim(bountag(i))
   write(varname,'("boundtag(",i0,")")') i
   informationID = XdmfAddInformation(xdmfFortranObj, & 
      trim(adjustl(varname))//char(0), trim(adjustl(info))//char(0))   
   write(info,fmt='(3f15.7)') amig(i),ff(i),face(i)
   write(varname,'("amig(",i0,"),ff(",i0,"),face(",i0,")")' ) &
       i, i, i
   informationID = XdmfAddInformation(xdmfFortranObj, & 
      trim(adjustl(varname))//char(0), trim(adjustl(info))//char(0))   
end do
do i=1, nbfr
   write(info,'(a)') trim(nbfr_alpha(i))
   write(varname,'("nbfr_alpha(",i0,")")') i
   informationID = XdmfAddInformation(xdmfFortranObj, & 
      trim(adjustl(varname))//char(0), trim(adjustl(info))//char(0))   
   do j=1,neta
      write(varname,'("emo(",i0,",",i0,"),efa(",i0,",",i0,")")' ) &
       i, j, i, j
      write(info,fmt='(2f15.7)') emo(i,j),efa(i,j)
         informationID = XdmfAddInformation(xdmfFortranObj, & 
            trim(adjustl(varname))//char(0), trim(adjustl(info))//char(0))
   end do
end do
write(info,fmt='(f15.7)') anginn
informationID = XdmfAddInformation(xdmfFortranObj, & 
   'anginn'//char(0), trim(adjustl(info))//char(0))
if ((nfluxf.ne.0).and.(nffr.ne.0)) then
   write(info,fmt='(i0)') nffr
   informationID = XdmfAddInformation(xdmfFortranObj, & 
      'nffr'//char(0), trim(adjustl(info))//char(0))
   if (nffr.gt.0) then
      do i=1,nffr
         write(info,'(a)') trim(fbountag(i))
         write(varname,'("fboundtag(",i0,")")') i
         informationID = XdmfAddInformation(xdmfFortranObj, & 
            trim(adjustl(varname))//char(0), trim(adjustl(info))//char(0))
         write(info,fmt='(3f6.3)') &
            famig(i),fff(i),fface(i)
         write(varname,'("famig(",i0,"),fff(",i0,"),fface(",i0,")")' ) &
            i, i, i
         informationID = XdmfAddInformation(xdmfFortranObj, & 
            trim(adjustl(varname))//char(0), trim(adjustl(info))//char(0))
      end do 
      do i=1,nffr
         write(info,'(a)') trim(nffr_alpha(i))
         write(varname,'("nffr_alpha(",i0,")")') i
         informationID = XdmfAddInformation(xdmfFortranObj, & 
            trim(adjustl(varname))//char(0), trim(adjustl(info))//char(0))
         do j=1,nvel
            select case(lbcodei(j))
            case(2,12,22,52)
               write(info,fmt='(2f6.3)') qnam(i,j), qnph(i,j)
               write(varname,'("qnam(",i0,",",i0,"),qnph(",i0,",",i0,")")' ) &
                  i, j, i, j
               informationID = XdmfAddInformation(xdmfFortranObj, & 
                  trim(adjustl(varname))//char(0), trim(adjustl(info))//char(0))
            case default
               cycle
            end select
         end do
      end do
   endif
endif

write(info,fmt='(i0,1x,2f6.3,1x,i0)') noute,toutse,toutfe,nspoole
informationID = XdmfAddInformation(xdmfFortranObj,'noute,toutse,toutfe,nspoole'//char(0), & 
   trim(adjustl(info))//char(0))
write(info,fmt='(i0)') nstae
informationID = XdmfAddInformation(xdmfFortranObj,'nstae'//char(0), & 
   trim(adjustl(info))//char(0))
stationType = 'sea_surface_height_above_geoid'

call writeStationsXDMF(nstae, xel, yel, stationType, xdmfFortranObj)

write(info,fmt='(i0,1x,2f15.7,1x,i0)') noutv,toutsv,toutfv,nspoolv
informationID = XdmfAddInformation(xdmfFortranObj,'noutv,toutsv,toutfv,nspoolv'//char(0), & 
   trim(adjustl(info))//char(0))
write(info,fmt='(i0)') nstav
informationID = XdmfAddInformation(xdmfFortranObj,'nstav'//char(0), & 
   trim(adjustl(info))//char(0))
stationType = 'water_column_vertically_averaged_velocity'
call writeStationsXDMF(nstav, xev, yev, stationType, xdmfFortranObj)
if (im.eq.10) then
   write(info,fmt='(i0,1x,2f15.7,1x,i0)') noutc,toutsc,toutfc,nspoolc
   informationID = XdmfAddInformation(xdmfFortranObj,'noutc,toutsc,toutfc,nspoolc'//char(0), trim(adjustl(info))//char(0))
   write(info,fmt='(i0)') nstac
   informationID = XdmfAddInformation(xdmfFortranObj,'nstac'//char(0), & 
      trim(adjustl(info))//char(0))
   stationType = 'water_column_vertically_averaged_concentration'   
   call writeStationsXDMF(nstav, xev, yev, stationType, xdmfFortranObj)
endif
if (nws.ne.0) then
   write(info,fmt='(i0,1x,2f15.7,1x,i0)') noutm,toutsm,toutfm,nspoolm
   informationID = XdmfAddInformation(xdmfFortranObj,'noutm,toutsm,toutfm,nspoolm'//char(0), trim(adjustl(info))//char(0))   
   write(info,fmt='(i0)') nstam
   informationID = XdmfAddInformation(xdmfFortranObj,'nstam'//char(0), & 
      trim(adjustl(info))//char(0))   
   stationType = 'atmospheric_pressure_and_wind_velocity_at_sea_level'   
   call writeStationsXDMF(nstam, xem, yem, stationType, xdmfFortranObj)
endif
write(info,fmt='(i0,1x,2f15.7,1x,i0)') noutge,toutsge,toutfge,nspoolge
informationID = XdmfAddInformation(xdmfFortranObj,'noutge,toutsge,toutfge,nspoolge'//char(0), trim(adjustl(info))//char(0))
write(info,fmt='(i0,1x,2f15.7,1x,i0)') noutgv,toutsgv,toutfgv,nspoolgv
informationID = XdmfAddInformation(xdmfFortranObj,'noutgv,toutsgv,toutfgv,nspoolgv'//char(0), trim(adjustl(info))//char(0))
if (im.eq.10) then
   write(info,fmt='(i0,1x,2f15.7,1x,i0)') noutgc,toutsgc,toutfgc,nspoolgc
   informationID = XdmfAddInformation(xdmfFortranObj,'noutgc,toutsgc,toutfgc,nspoolgc'//char(0), trim(adjustl(info))//char(0))
endif
if (nws.ne.0) then
   write(info,fmt='(i0,1x,2f15.7,1x,i0)') noutgw,toutsgw,toutfgw,nspoolgw
   informationID = XdmfAddInformation(xdmfFortranObj,'noutgw,toutsgw,toutfgw,nspoolgw'//   char(0), trim(adjustl(info))//char(0))   
endif

write(info,fmt='(i0)') nfreq
informationID = XdmfAddInformation(xdmfFortranObj,'nfreq'//char(0), & 
   trim(adjustl(info))//char(0))
if (nfreq.ne.0) then
   do i=1,nfreq
      write(info,'(a10)') namefr(i)
      write(varname,'("namefr(",i0,")")') i
         informationID = XdmfAddInformation(xdmfFortranObj, & 
            trim(adjustl(varname))//char(0), trim(adjustl(info))//char(0))      
      write(info,fmt='(3f15.7)') hafreq(i),haff(i),haface(i)
      write(varname,'("hafreq(",i0,"),haff(",i0,"),haface(",i0,")")' ) &
         i, i, i
      informationID = XdmfAddInformation(xdmfFortranObj, & 
         trim(adjustl(varname))//char(0), trim(adjustl(info))//char(0))
   end do
endif
write(info,fmt='(2f15.7,1x,i0,1x,f15.7)') thas,thaf,nhainc,fmv
informationID = XdmfAddInformation(xdmfFortranObj,'thas,thaf,nhainc,fmv' &
   //char(0), trim(adjustl(info))//char(0))
write(info,fmt='(i0,1x,i0,1x,i0,1x,i0)') nhase,nhasv,nhage,nhagv
informationID = XdmfAddInformation(xdmfFortranObj,'nhase,nhasv,nhage,nhagv' &
   //char(0), trim(adjustl(info))//char(0))
write(info,fmt='(i0,1x,i0)') nhstar,nhsinc
informationID = XdmfAddInformation(xdmfFortranObj,'nhstar,nhsinc' &
   //char(0), trim(adjustl(info))//char(0))
write(info,fmt='(i0,1x,i0,f15.7,1x,i0)') ititer,isldia,convcr,itmax
informationID = XdmfAddInformation(xdmfFortranObj,'ititer,isldia,convcr,itmax' &
   //char(0), trim(adjustl(info))//char(0))
!
! read metadata if netcdf or xdmf format was specified in the output 
! identifiers
outputSpecifiers(1) = noute
outputSpecifiers(2) = noutv
outputSpecifiers(3) = noutc
outputSpecifiers(4) = noutm
outputSpecifiers(5) = noutge
outputSpecifiers(6) = noutgv
outputSpecifiers(7) = noutgc
outputSpecifiers(8) = noutgw
outputSpecifiers(9) = nhstar
! now write metadata if any output format indicated it
if (readMetaData.eqv..true.) then 
   informationID = XdmfAddInformation(xdmfFortranObj,'title'//char(0), &
      trim(adjustl(title))//char(0))   
   informationID = XdmfAddInformation(xdmfFortranObj,'institution'//char(0), &
      trim(adjustl(institution))//char(0))
   informationID = XdmfAddInformation(xdmfFortranObj,'source'//char(0), &
      trim(adjustl(source))//char(0))   
   informationID = XdmfAddInformation(xdmfFortranObj,'history'//char(0), &
      trim(adjustl(history))//char(0))   
   informationID = XdmfAddInformation(xdmfFortranObj,'references'//char(0), &
      trim(adjustl(references))//char(0))   
   informationID = XdmfAddInformation(xdmfFortranObj,'comments'//char(0), &
      trim(adjustl(comments))//char(0))   
   informationID = XdmfAddInformation(xdmfFortranObj,'host'//char(0), &
      trim(adjustl(host))//char(0))   
   informationID = XdmfAddInformation(xdmfFortranObj,'convention'//char(0), &
      trim(adjustl(convention))//char(0))   
   informationID = XdmfAddInformation(xdmfFortranObj,'contact'//char(0), &
      trim(adjustl(contact))//char(0))   
   informationID = XdmfAddInformation(xdmfFortranObj,'base_date'//char(0), &
      trim(adjustl(base_date))//char(0))
endif 
!---------------------------------------------------------------------
end subroutine writeControlXDMF
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
! S U B R O U T I N E    W R I T E   S T A T I O N S   X D M F
!---------------------------------------------------------------------
subroutine writeStationsXDMF(nsta, stax, stay, stationType, xdmfFortranObj)
implicit none
include 'Xdmf.f'
integer, intent(in) :: nsta
real(sz), intent(in) :: stax(:)
real(sz), intent(in) :: stay(:)
character(len=80), intent(in) :: stationType
integer*8, intent(in) :: xdmfFortranObj ! XDMF object to receive data
integer :: informationID ! information index, used as ID to insert array
!
! add the station locations as an array related to the station information object
informationID = XdmfAddInformation(xdmfFortranObj,'station longitude'//char(0), & 
   trim(adjustl(stationType))//char(0))
! arguments are: fortran obj, information index, values, numValues, arrayType
call XdmfAddInformationArray(xdmfFortranObj, informationID, stax, nsta, &
   XDMF_ARRAY_TYPE_FLOAT64)
informationID = XdmfAddInformation(xdmfFortranObj,'station latitude'//char(0), & 
   trim(adjustl(stationType))//char(0))
call XdmfAddInformationArray(xdmfFortranObj, informationID, stay, nsta, &
   XDMF_ARRAY_TYPE_FLOAT64)
!---------------------------------------------------------------------
end subroutine writeStationsXDMF
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
end module control
!---------------------------------------------------------------------
