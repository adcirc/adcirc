module swn_outnc
!
!   --|-----------------------------------------------------------|--
!     | BMT ARGOSS                                                |
!     | Voorsterweg 28, 8316 PT Vollenhove                        |
!     | http://www.bmtargoss.com                                  |
!     |                                                           |
!     |                                                           |
!     | Programmer: A.Th.C. Hulst                                 |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2016  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!   Authors
!
!    0.01: A.Th.C.Hulst
!
!   Updates
!
!    0.01, Feb 2012: New Module
!
!   Purpose
!
!   Write SWAN output into netcdf files
!
!   Method
!
!   MODULE construct
!
!   Modules used
!
!   NETCDF          : Standard netcdf library
!   AGIONCMD        : Set of utilities for writing (metocean) data to netcdf
!   NCTABLEMD       : description table for metocean parameters
    use agioncmd
    use M_PARALL
    use NETCDF
    use nctablemd, only: nctable, nctable_record, get_nctable_record
    use OUTP_DATA, only: ORQDAT, MAX_OUTP_REQ, LCOMPGRD, OUTP_FILES, NREOQ
    use OCPCOMM2,  only: PROJNR, PROJID, VERTXT
    use OCPCOMM4
    use SWCOMM1,   only: CHTIME, OVEXCV
    use SWCOMM2,   only: XOFFS, YOFFS, OPTG, EXCFLD
    use SWCOMM3,   only: NSTATM, ALPC, MDC, MSC, MCGRD, MXC, MYC, DNORTH, &
                         PI2, ICUR, BNAUT
    use SWCOMM4
    use TIMECOMM,  only: TFINC, TIMCO, TINIC, DT
!
    implicit none
!
!   Module parameters
!
!      oqi(1)       :   saved alias for nref
!      oqi(2)       :   index of request
!      oqi(3)       :   number of output variables
!      oqi(4)       :   idla of output request
!      oqr(1)       :   tbegin
!      oqr(2)       :   delta
!
!   Module variables
    character*40                :: STNAMES(71,2) = ''
    logical                     :: stnames_initialized = .false., &
                                   skip_range_error = .true.
    type(recordaxe_type), save  :: recordaxe(MAX_OUTP_REQ)
    integer,save                :: ncoffset(MAX_OUTP_REQ) = 0

!PUN    logical                     :: PUNSWAN = .true.
    logical                     :: PUNSWAN = .false.
    type spcaux_type
        real, dimension(:), allocatable          :: hs, depth, ux, uy, wndx, wndy, &
                                                    xc, yc, xp, yp, f, theta
        integer, dimension(:), allocatable       :: ips
        real, dimension(:,:), allocatable        :: E
        logical                                  :: is2d = .false., relative = .false.
        integer                                  :: mip = 0, ndir = 0, nfreq = 0
    end type spcaux_type

!
!   Source text
!

contains
    subroutine swn_outnc_spec(RTYPE, OQI, OQR, MIP, VOQR, VOQ, AC2, &
                              SPCSIG, SPCDIR, DEP2, KGRPNT, CROSS, IONOD)
!
!
!   --|-----------------------------------------------------------|--
!     | BMT ARGOSS                                                |
!     | Voorsterweg 28, 8316 PT Vollenhove                        |
!     | http://www.bmtargoss.com                                  |
!     |                                                           |
!     |                                                           |
!     | Programmer: A.Th.C. Hulst                                 |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2016  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!
!  1. Updates
!
!
!  2. Purpose
!
!     Write action density spectrum to netcdf file when serial OR
!     binary file when PARLL or PUNSWAN
!
!  3. Method
!
!     ---
!
!  4. Argument variables
!
! IO..name....type...size.......description
! i   RTYPE   char   *          type of output request: 'SPEC' for 2-D spectral
!                                                       'SPE1' for 1-D freq. spectrum
! i   IONOD   int    MIP        array indicating in which subdomain output points
!                               are located
! i   SPCDIR  real   MDC        (*,1); spectral directions (radians)
!                               (*,2); cosine of spectral directions
!                               (*,3); sine of spectral directions
!                               (*,4); cosine^2 of spectral directions
!                               (*,5); cosine*sine of spectral directions
!                               (*,6); sine^2 of spectral directions
! i   SPCSIG real   MSC         Relative frequencies in computational domain in sigma-space
! i   VOQR   int    *           Adminstrative list tracks which column in VOQ has which
!                               variable
! i   CROSS  bool   4,MIP
! i   KGRPNT int    MXC,MYC
! i   MIP    int    1           number of output points
! i   VOQ    real   *           Array with "integrated" parameters
! i   AC2    real   MDC,MSC,    2D spectra
!                   MCGRD
! i   DEP2   real   MCGRD
!
      real,                 intent(   in) :: SPCDIR(MDC,6), SPCSIG(MSC)
      logical,              intent(   in) :: CROSS(1:4,1:MIP)
      character (len=*),    intent(   in) :: RTYPE
      integer,              intent(   in) :: VOQR(*), &
                                             KGRPNT(MXC,MYC), IONOD(*), &
                                             MIP
      integer,              intent(inout) :: OQI(4)
      real,                 intent(   in) :: OQR(2)
      real,                 intent(   in) :: VOQ(MIP,*), AC2(MDC,MSC,MCGRD), &
                                             DEP2(MCGRD)
!
!  5. Local variables
!
      integer                             :: irq
      integer, save                       :: IENT=0
      type(spcaux_type)                   :: spcaux
!  8. Subroutines used
!
!     swn_outnc_spcaux
!     swn_outnc_openspecfile
!     swn_outnc_appendspc
!     nccheck
!     nf90_close
!     swn_outnc_deallocate_spcaux
!
      logical STPNOW
!
!  9. Subroutines calling
!
!     SWOUTP (SWAN/OUTP)
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!
! 12. Structure.
!
!    a. compute spectra and auxilary data (writes binary file if parll or punswan)
!    if serial
!      b. open / create netcdf file
!      c. write spectra to netcdf
!      d. close ncfile
!
! 13. Source text
!
      if (LTRACE) call STRACE (IENT,'swn_outnc_spec')
      irq = OQI(2)
      call swn_outnc_spcaux(RTYPE, OQI, MIP, VOQR, VOQ, AC2, &
                            SPCSIG, SPCDIR, DEP2, KGRPNT, CROSS, IONOD, &
                            spcaux)
      if (STPNOW()) return

!
! b. open / create netcdf file (generic output file)
!
      if ( .not. (PARLL .or. PUNSWAN) ) then

        ! opening / creating the netcdf file. Sets module variable recordaxe(irq)
        call swn_outnc_openspecfile(OUTP_FILES(irq), spcaux, OQI, OQR)
        if (STPNOW()) return

        call swn_outnc_appendspc(OQI, spcaux)
        if (STPNOW()) return

        call nccheck(nf90_close(oqi(1) + ncoffset(irq) ))
        oqi(1) = 0
      end if

      call swn_outnc_deallocate_spcaux( spcaux )

    end subroutine swn_outnc_spec

    subroutine swn_outnc_colspc ( RTYPE, OQI, OQR, MIP, KGRPGL )
      use SwanGriddata, only: xcugrdgl, ycugrdgl, nverts
!PUN      USE OCPCOMM2, ONLY: LENFNM, DIRCH2
!PUN      USE SIZES, ONLY: GLOBALDIR, LOCALDIR
!
!
!   --|-----------------------------------------------------------|--
!     | BMT ARGOSS                                                |
!     | Voorsterweg 28, 8316 PT Vollenhove                        |
!     | http://www.bmtargoss.com                                  |
!     |                                                           |
!     |                                                           |
!     | Programmer: A.Th.C. Hulst                                 |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2016  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!
!  1. Updates
!
!
!  2. Purpose
!
!     Collect spectra from binary files and store in netcdf file
!
!  3. Method
!
!     ---
!
!  4. Argument variables
      character (len=*),    intent(   in) :: RTYPE * 4
      integer,              intent(inout) :: OQI(4)
      real,                 intent(   in) :: OQR(2)
      integer,              intent(   in) :: MIP, KGRPGL(MXCGL,MYCGL)
!
! IO..name....type...size.......description
! i   RTYPE   char   *          type of output request: 'SPEC' for 2-D spectral
!                                                       'SPE1' for 1-D freq. spectrum
! i   XP, YP  real   MIP        coordinates of all output points
! i   MIP    int    1           number of output points
!
!  5. Local variables
!
      logical                             :: lopen, read_header, file_exists, &
                                             EQREAL, STPNOW
      integer                             :: ip, ierr, otype, xpctime2, &
                                             i, tmip, binnr, irq, &
                                             tip, ips(MIP), iproc, nref, ilpos
      integer, save                       :: IENT=0
      character*80                        :: binfile, ncfile, basefile
      character*256                       :: errmsg
      type(spcaux_type)                   :: spcaux, lspcaux
      real                                :: xref, yref

!
!  9. Subroutines calling
!
!     SWCOLOUT
!
! 10. Error messages
!
!
! 11. Remarks
!
!
! 12. Structure
!
!
! 13. Source text
!
        if (LTRACE) call STRACE (IENT,'swn_outnc_colspc')
        ierr    = 0
!
! a. Prepare spectra/auxilary data struct
!
        irq        = OQI(2)
        spcaux%mip = MIP
        allocate(spcaux%hs(MIP), spcaux%depth(MIP), spcaux%ux(MIP), spcaux%uy(MIP), &
                 spcaux%wndx(MIP), spcaux%wndy(MIP), spcaux%xc(MIP), spcaux%yc(MIP))
        if ( LCOMPGRD .and. PUNSWAN) then
            allocate(spcaux%xp(nverts), spcaux%yp(nverts))
        else
            allocate(spcaux%xp(MIP), spcaux%yp(MIP))
        end if

        spcaux%hs      = OVEXCV(10)
        spcaux%depth   = OVEXCV(4)
        spcaux%ux      = OVEXCV(5)
        spcaux%uy      = OVEXCV(5)
        spcaux%wndx    = OVEXCV(26)
        spcaux%wndy    = OVEXCV(26)

        if (RTYPE(3:3) .eq. 'R') spcaux%relative = .true.
        if (RTYPE(4:4) .eq. 'C') then
            allocate(spcaux%E(MSC * MDC, MIP))
            spcaux%is2d = .true.
        else
            allocate(spcaux%E(3 * MSC, MIP))
            spcaux%is2d = .false.
        end if
        spcaux%E     = 0
        spcaux%nfreq = MSC
        spcaux%ndir  = MDC

        allocate( spcaux%f(MSC), spcaux%theta(MDC))
!
! open and read binary files
!
        binfile = OUTP_FILES(irq)

        ! nref is unitnr set while writing the binary file.
        nref = HIOPEN + IRQ
        basefile = outp_files(irq)
        ilpos  = index(basefile, '-0')
        if (ilpos.ne.0) basefile = basefile(1:ilpos-1)
!PUN        ilpos   = len(trim(LOCALDIR))+2
!PUN        basefile = basefile(ilpos:LENFNM)
!PUN        ilpos   = index(LOCALDIR, '0')-1

        do iproc = 1, NPROC
            binnr = HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+iproc

!PUN              write(binfile, '(A, I4.4, A, A)') LOCALDIR(1:ilpos), iproc - 1, DIRCH2, trim(basefile)
              write(binfile, '(A,"-",I3.3)') trim(basefile), iproc
            inquire(unit=binnr, opened=lopen)

            if ( .not. lopen) then
                open(file=binfile, unit=binnr, form='unformatted', IOSTAT=ierr)
                if (ierr /= 0) then
                    write(errmsg, '("Opening ",A, " on ",I3, " failed with IOSTAT ",I3)') trim(binfile), binnr, ierr
                    call MSGERR(4, errmsg)
                    return
                end if
                read_header = .true.
            else
                read_header = .false.
            end if
            inquire(unit=binnr, opened=lopen)
            if ( .not. lopen ) then
                write(errmsg,'("unit ",I3," is NOT open")'), binnr
                call MSGERR(4, errmsg)
                return
            end if

            ! In the final collocation the spatial, spectral and directional grids are unknown.
            ! They are exchanged via the binary file of the master process.
            if ( read_header .and. iproc == 1 ) then
                read(binnr, IOSTAT=ierr) spcaux%xp
                if ( ierr /= 0 ) then
                  write(errmsg,'("Reading grid data from ", A, " returned IOSTAT ", I3.3)') trim(binfile), ierr
                  call MSGERR(4, errmsg)
                  return
                end if
                read(binnr) spcaux%yp
                if ( PUNSWAN .and. LCOMPGRD ) then
                    deallocate(spcaux%xp, spcaux%yp)
                    allocate(spcaux%xp(MIP), spcaux%yp(MIP))
                    spcaux%xp = xcugrdgl
                    spcaux%yp = ycugrdgl
                end if
                read(binnr) spcaux%f
                read(binnr) spcaux%theta
            end if

            read(binnr, IOSTAT=ierr) xpctime2, tmip
            if ( ierr /= 0 ) then
                write(errmsg,'("Reading time and tmip from ", A, "(",I3,") , time ", A, ", returned IOSTAT ", I3.3)') &
                     trim(binfile), binnr, CHTIME, ierr
                call MSGERR(4, errmsg)
                return
            end if
            if ( tmip == 0 ) cycle

            read(binnr) ips(1:tmip)
            do tip = 1, tmip
              ip = ips(tip)
              read(binnr, iostat = ierr) xref, yref, spcaux%E(:, ip), spcaux%depth(ip), &
                          spcaux%ux(ip), spcaux%uy(ip), spcaux%hs(ip), spcaux%wndx(ip), spcaux%wndy(ip)

              if ( ierr /= 0 ) then
                write(errmsg,'("Reading IP ",I3.3, ", TIP ", I3.3, " in ", A, " returned IOSTAT ", I3.3)') &
                  ip, tip, binfile, ierr
                call MSGERR(4, errmsg)
                return
              end if
            end do ! points
        end do !nproc
        !
        ! MPI with a cold start shows some interesting initial values
        ! remove them
        !
        where(spcaux%E < epsilon(1.)) spcaux%E = 0.
!
! b. netcdf IO
!

        ncfile = outp_files(irq)
        ilpos  = index(ncfile, '-0')
        if (ilpos.ne.0) ncfile = ncfile(1:ilpos-1)
!PUN         ilpos   = len(trim(LOCALDIR))+2
!PUN         ncfile = ncfile(ilpos:LENFNM)
!PUN         ncfile = trim(GLOBALDIR)//DIRCH2//trim(ncfile)

        call swn_outnc_openspecfile(ncfile, spcaux, OQI, OQR)
        if (STPNOW()) return

        call swn_outnc_appendspc(OQI, spcaux, xpctime2)
        if (STPNOW()) return

        call swn_outnc_deallocate_spcaux( spcaux )

        if ( NSTATM == 0 .or. TIMCO >= TFINC ) then
            call close_ncfile(OQI(1) + ncoffset(irq) )
            OQI(1) = 0

            ! The intermediate binary files are closed and deleted in
            !    swanparll.ftn -> SWCOLOUT

        else
            call nccheck( NF90_SYNC( OQI(1) + ncoffset(irq)  ))
        end if
    end subroutine swn_outnc_colspc

    subroutine swn_outnc_spcaux(RTYPE, OQI, MIP, VOQR, VOQ, AC2, &
                              SPCSIG, SPCDIR, DEP2, KGRPNT, CROSS, IONOD, &
                              lspcaux)
      USE OCPCOMM2
      USE TIMECOMM, only: TINIC, TFINC, TIMCO
!PUN      USE SIZES, ONLY: MYPROC
!
!
!   --|-----------------------------------------------------------|--
!     | BMT ARGOSS                                                |
!     | Voorsterweg 28, 8316 PT Vollenhove                        |
!     | http://www.bmtargoss.com                                  |
!     |                                                           |
!     |                                                           |
!     | Programmer: A.Th.C. Hulst                                 |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2016  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!
!  1. Updates
!
!
!  2. Purpose
!
!     Write action density spectrum to netcdf file
!
!  3. Method
!
!     ---
!
!  4. Argument variables
!
! IO..name....type...size.......description
! i   RTYPE   char   *          type of output request: 'SPEC' for 2-D spectral
!                                                       'SPE1' for 1-D freq. spectrum
! i   IONOD   int    MIP        array indicating in which subdomain output points
!                               are located
! i   SPCDIR  real   MDC        (*,1); spectral directions (radians)
!                               (*,2); cosine of spectral directions
!                               (*,3); sine of spectral directions
!                               (*,4); cosine^2 of spectral directions
!                               (*,5); cosine*sine of spectral directions
!                               (*,6); sine^2 of spectral directions
! i   SPCSIG real   MSC         Relative frequencies in computational domain in sigma-space
! i   VOQR   int    *           Adminstrative list tracks which column in VOQ has which
!                               variable
! i   CROSS  bool   4,MIP
! i   KGRPNT int    MXC,MYC
! i   MIP    int    1           number of output points
! i   VOQ    real   *           Array with "integrated" parameters
! i   AC2    real   MDC,MSC,    2D spectra
!                   MCGRD
! i   DEP2   real   MCGRD
!
      real,                            intent(   in) :: SPCDIR(MDC,6), SPCSIG(MSC)
      logical,                         intent(   in) :: CROSS(1:4,1:MIP)
      character (len=*),               intent(   in) :: RTYPE
      integer,                         intent(   in) :: OQI(4), VOQR(*), &
                                                        KGRPNT(MXC,MYC), IONOD(*), &
                                                        MIP
      real,                            intent(   in) :: VOQ(MIP,*), AC2(MDC,MSC,MCGRD), &
                                                        DEP2(MCGRD)
      type(spcaux_type),               intent(  out) :: lspcaux
!
!  5. Local variables
!
      logical                                        :: EQREAL, lopen, do_open_files, write_header
      integer                                        :: ip, ierr, otype, xpctmp(2), xpctime, &
                                                        pnr, ri, i, tmip, irq, ips(MIP), iproc, &
                                                        binnr, xi, yi, npnts, kgrpnt1d(MXC*MYC)
      integer, save                                  :: IENT=0
      character*80                                   :: outfile
      character*256                                  :: errmsg
!  8. Subroutines used
!
!     SWCMSP
!
      logical STPNOW
!
!  9. Subroutines calling
!
!     SWOUTP (SWAN/OUTP)
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     - PUNSWAN support was added later on. In punswan, PARLL is always false but the
!       IAMMASTER is set. For most output purposes, PARLL should be considered
!       true. Hence a logical PUNSWAN is set at the top if this module. In the structure
!       below PARLL means if PARLL .or. PUNSWAN.
!
! 12. Structure.
!    a. decide wether relative and 2D spectra are required
!    b. collect depth, current, xc, yc, f, theta
!    c  if PARLL open unformatted binary
!    d. set output type (1D or 2D, relative or absolute) and allocate Energy array
!    e. if PARLL, write time and number of output points in this subgrid to
!       unformatted binary file
!    for each point...
!      f.1. compute relative / absolute 1D or 2D energy
!      f.2. if PARLL spectra and auxillary data to unformatted binary file
!
! 13. Source text
!
      iproc = INODE
!PUN      iproc = MYPROC
      if (LTRACE) call STRACE (IENT,'swn_outnc_spcaux')
      lspcaux%mip = MIP
      allocate(lspcaux%hs(MIP), lspcaux%depth(MIP), lspcaux%ux(MIP), lspcaux%uy(MIP), &
               lspcaux%wndx(MIP), lspcaux%wndy(MIP), lspcaux%xc(MIP), lspcaux%yc(MIP), &
               lspcaux%xp(MIP), lspcaux%yp(MIP))
!
! a. decide wether relative and 2D spectra are required
!
      irq      = OQI(2)
      ierr     = 0

      if (RTYPE(3:3) .eq. 'R') lspcaux%relative = .true.
      if (RTYPE(4:4) .eq. 'C') lspcaux%is2d     = .true.
!
! b. collect auxillary data / set default (task e. set hs and wind values)
!
      lspcaux%depth = VOQ(:,VOQR(4))
      lspcaux%hs    = 0
      lspcaux%wndx  = 0
      lspcaux%wndy  = 0
      if (ICUR > 0) then
          lspcaux%ux = VOQ(:,VOQR(5))
          lspcaux%uy = VOQ(:,VOQR(5)+1)
      else
          lspcaux%ux = 0.
          lspcaux%uy = 0.
      end if

      lspcaux%xc = OVEXCV(1)
      lspcaux%yc = OVEXCV(2)
      lspcaux%xp = OVEXCV(1)
      lspcaux%yp = OVEXCV(2)

      do ip = 1, MIP
        if (OPTG.ne.5) then
          if (.not.EQREAL(VOQ(ip,VOQR(24)),OVEXCV(1) )) then
            lspcaux%xc(ip) = VOQ(ip,VOQR(24))
            lspcaux%yc(ip) = VOQ(ip,VOQR(25))
          end if
        else
          if (.not.EQREAL(VOQ(ip,1),OVEXCV(1))) then
            lspcaux%xc(ip) = VOQ(ip,1) - XOFFS
            lspcaux%yc(ip) = VOQ(ip,2) - YOFFS
          end if
        end if
        lspcaux%xp(ip) = VOQ(ip,1)
        lspcaux%yp(ip) = VOQ(ip,2)
      end do
      lspcaux%nfreq = MSC
      lspcaux%ndir  = MDC

      allocate( lspcaux%f(MSC), lspcaux%theta(MDC))
      lspcaux%f     = SPCSIG/PI2
      if ( BNAUT ) then
        ! meteorological convention
        lspcaux%theta = PI + DNORTH*PI/180 - SPCDIR(:,1)
      else
        lspcaux%theta = SPCDIR(:,1)
      end if
!
! c if PARLL or PUNSWAN open unformatted binary if not already opened
!
      if ( PARLL .or. PUNSWAN) then
          if ( oqi(1) == 0 ) then
              call FOR(oqi(1), OUTP_FILES(irq), 'UU', ierr)
              if ( ierr > 0 ) then
                  write(errmsg,'("File ", A, " open returned IOSTAT ", I3)'), trim(OUTP_FILES(irq)), ierr
                  call MSGERR(4, errmsg)
                  return
              end if
              write_header = .true.
          else
              write_header = .false.
          end if

      end if

! d. set output type (1D or 2D, relative or absolute) and allocate Energy array
      if ( lspcaux%is2d ) then
          allocate(lspcaux%E(MSC * MDC, MIP))
          if ( lspcaux%relative ) then
              otype = 2
          else
              otype = -2
          end if
      else
          allocate(lspcaux%E(3 * MSC, MIP))
          if ( lspcaux%relative ) then
              otype = 1
          else
              otype = -1
          end if
      end if
      lspcaux%E = 0

! e. write time and number of output points in this subgrid to unformatted binary file
      if ( PARLL .or. PUNSWAN) then
        tmip = 0
        ! if part of this subgrid and an active point
        do ip = 1, MIP
            if ( IONOD(ip).eq.IPROC ) then
                tmip      = tmip + 1
                ips(tmip) = ip
            end if
        end do
        if ( NSTATM > 0 ) then
            READ (chtime, '(I8,1X,I6)') (xpctmp(i), i=1,2)
            xpctime = seconds_since_epoch(xpctmp(1), xpctmp(2))
        else
            xpctime = 0
        end if

        ! In the final collocation the spatial, spectral and directional grids are unknown.
        ! They are exchanged via the binary file of the master on the first timestep.
        !
        ! note the statement below and its reading counterpart imply that the master is
        ! always IPROC 1
        if ( write_header .and. IAMMASTER ) then
            write(OQI(1)) lspcaux%xp
            write(OQI(1)) lspcaux%yp
            write(OQI(1)) lspcaux%f
            write(OQI(1)) lspcaux%theta
        end if

        write(OQI(1), IOSTAT=ierr) xpctime, tmip
        if (ierr /= 0) then
          write(errmsg, '(A, "-> IOSTAT=", I3, " time=", I10, " TMIP=", I4, " OQI(1)=", I3)') &
                  trim(outfile), ierr, xpctime, tmip, oqi(1)
          call MSGERR(4, errmsg)
          return
        end if
        if ( tmip /= 0 ) write(OQI(1)) ips

      end if

! for each point, compute spectra and collect remaining auxilary
      do ip = 1, MIP
        if ( (.not. PARLL .and. .not. PUNSWAN) .or. &
             IONOD(ip).eq.INODE ) then
!PUN             IONOD(ip).eq.MYPROC ) then
          ! if depth > 0 and not a dummy value
          if ( lspcaux%depth(ip) > epsilon(1.) .and. abs(lspcaux%depth(ip) - OVEXCV(4)) > epsilon(1.) ) then
            call SWCMSP (otype, lspcaux%xc(ip), lspcaux%yc(ip), AC2, lspcaux%E(:, ip), SPCSIG, &
                         lspcaux%depth(ip) , dep2, lspcaux%ux(ip), lspcaux%uy(ip),        &
                         SPCDIR(1,2), SPCDIR(1,3), 1., KGRPNT,  &
                         CROSS(1,IP), ierr)

            ! when ierr > 0 then energy == 0. That's fine.

            lspcaux%hs(ip)   = VOQ(ip,VOQR(10))
            lspcaux%wndx(ip) = VOQ(ip,VOQR(26))
            lspcaux%wndy(ip) = VOQ(ip,VOQR(26)+1)
          end if

          if ( PARLL .or. PUNSWAN ) then
            ! write to binary
            write(OQI(1)) lspcaux%xp(ip), lspcaux%yp(ip), lspcaux%E(:, ip), lspcaux%depth(ip), &
                          lspcaux%ux(ip), lspcaux%uy(ip), lspcaux%hs(ip), lspcaux%wndx(ip), &
                          lspcaux%wndy(ip)
          end if
        end if
      end do

      if ( PARLL .or. PUNSWAN ) then
        ! in parallel mode, close binary file when required
        if (NSTATM == 0 .or. TIMCO >= TFINC) then
          inquire(unit=OQI(1), opened=lopen)
          if ( lopen ) then
            close(unit=OQI(1))
          else
            write(errmsg, '("File ", A, " is already closed. This is a code error")') trim(OUTP_FILES(irq))
            call MSGERR(4, errmsg)
            return
          end if
        end if
      end if

    end subroutine swn_outnc_spcaux

    subroutine swn_outnc_spcaux_on_wetnodes( lkgrpnt, lspcaux, spcaux )
        ! lkgrpnt is KGRPNT when called from shared memory mode or KGRPGL
        ! from colspc in the case of distributed memory. In both cases, it covers
        ! the full grid
        integer,                         intent(   in) :: lkgrpnt(MXCGL,MYCGL)
        type(spcaux_type),               intent(   in) :: lspcaux
        type(spcaux_type),               intent(  out) :: spcaux

        integer                                        :: ip, indx, npnts, ix, iy
        integer, save                                  :: IENT=0
        character*256                                  :: errmsg
        if (LTRACE) call STRACE(IENT,'swn_outnc_spcaux_on_wetnodes')

        npnts = MCGRDGL - 1
        allocate(spcaux%hs(npnts), spcaux%depth(npnts), spcaux%ux(npnts), spcaux%uy(npnts), &
                spcaux%wndx(npnts), spcaux%wndy(npnts), spcaux%xp(npnts), spcaux%yp(npnts), &
                spcaux%xc(npnts), spcaux%yc(npnts), spcaux%ips(npnts), &
                spcaux%E(size(lspcaux%E,1), npnts) )

        do ix = 1, MXCGL
          do iy = 1, MYCGL
            ip = (iy-1)*MXCGL + ix
            if ( lkgrpnt(ix,iy) /= 1 ) then
                indx               = lkgrpnt(ix,iy) - 1
                spcaux%ips(indx)   = ip
                spcaux%hs(indx)    = lspcaux%hs(ip)
                spcaux%depth(indx) = lspcaux%depth(ip)
                spcaux%ux(indx)    = lspcaux%ux(ip)
                spcaux%uy(indx)    = lspcaux%uy(ip)
                spcaux%wndx(indx)  = lspcaux%wndx(ip)
                spcaux%wndy(indx)  = lspcaux%wndy(ip)
                spcaux%xc(indx)    = lspcaux%xc(ip)
                spcaux%yc(indx)    = lspcaux%yc(ip)
                spcaux%xp(indx)    = lspcaux%xp(ip)
                spcaux%yp(indx)    = lspcaux%yp(ip)
                spcaux%E(:,indx)   = lspcaux%E(:,ip)
            end if
          end do
        end do

        if ( indx /= npnts ) then
            write(errmsg, '("Inconsistent number of wet nodes:", i5, i5)') indx, npnts
            call MSGERR(4, errmsg)
            return
        end if

        spcaux%is2d     = lspcaux%is2d
        spcaux%relative = lspcaux%relative

        spcaux%nfreq = lspcaux%nfreq
        spcaux%ndir  = lspcaux%ndir

        allocate( spcaux%f(spcaux%nfreq), spcaux%theta(spcaux%ndir))
        spcaux%f     = lspcaux%f
        spcaux%theta = lspcaux%theta
        spcaux%mip   = npnts

    end subroutine swn_outnc_spcaux_on_wetnodes

    subroutine swn_outnc_deallocate_spcaux( lspcaux )
        type(spcaux_type), intent(inout) :: lspcaux
        integer, save                    :: IENT=0

        if (LTRACE) call STRACE (IENT,'swn_outnc_deallocate_spcaux')

        ! yesyes, derived type fields are deallocated when out of scope
        ! but sometimes variables stay in scope
        if (allocated(lspcaux%hs))    deallocate(lspcaux%hs)
        if (allocated(lspcaux%depth)) deallocate(lspcaux%depth)
        if (allocated(lspcaux%ux))    deallocate(lspcaux%ux)
        if (allocated(lspcaux%uy))    deallocate(lspcaux%uy)
        if (allocated(lspcaux%wndx))  deallocate(lspcaux%wndx)
        if (allocated(lspcaux%wndy))  deallocate(lspcaux%wndy)
        if (allocated(lspcaux%xc))    deallocate(lspcaux%xc)
        if (allocated(lspcaux%yc))    deallocate(lspcaux%yc)
        if (allocated(lspcaux%xp))    deallocate(lspcaux%xp)
        if (allocated(lspcaux%yp))    deallocate(lspcaux%yp)
        if (allocated(lspcaux%f))     deallocate(lspcaux%f)
        if (allocated(lspcaux%theta)) deallocate(lspcaux%theta)
        if (allocated(lspcaux%E))     deallocate(lspcaux%E)
        if (allocated(lspcaux%ips))   deallocate(lspcaux%ips)

        lspcaux%mip   = 0
        lspcaux%ndir  = 0
        lspcaux%nfreq = 0
    end subroutine swn_outnc_deallocate_spcaux

    subroutine swn_outnc_appendspc(oqi, spcaux, xpctime2)
      USE OCPCOMM2
!
!
!   --|-----------------------------------------------------------|--
!     | BMT ARGOSS                                                |
!     | Voorsterweg 28, 8316 PT Vollenhove                        |
!     | http://www.bmtargoss.com                                  |
!     |                                                           |
!     |                                                           |
!     | Programmer: A.Th.C. Hulst                                 |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2016  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!
!  1. Updates
!
!
!  2. Purpose
!
!     Append spectra and auxilary variables to netcdf file.
!
!  3. Method
!
!     ---
!
!  4. Argument variables
!
! IO..name....type...size.......description
!
      integer,                      intent(   in) :: oqi(4)
      type(spcaux_type), target,    intent(inout) :: spcaux
      integer,           optional,  intent(   in) :: xpctime2
!
!  5. Local variables
!
      integer                             :: ierr, xpctmp(2), xpctime, pnr, ri, irq, i, &
                                             ncid
      character*256                       :: errmsg
      integer, save                       :: IENT=0
      logical                             :: spc_as_map, wetnode_list, noaux
      type(spcaux_type), target           :: lspcaux
      type(spcaux_type), pointer          :: pspcaux

!  8. Subroutines used
!
!     SWCMSP
!
      logical STPNOW
!
!  9. Subroutines calling
!
!     SWOUTP (SWAN/OUTP)
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!     Auxillary data is retrieved from the VOQ array. Please modify SWORDC when
!     VOQ need to be extended
!
! 12. Structure
!    a. Determine the index in the record axe of ncfile
!    b. add spectra to ncfile
!    c. add auxilary variables (depth, current and Hs) to ncfile
!
! 13. Source text
!
      if (LTRACE) call STRACE (IENT,'swn_outnc_appendspc')
      irq  = oqi(2)
      ncid = oqi(1) + ncoffset(irq)

      call swn_outnc_spcflags(oqi(4), spc_as_map, wetnode_list, spcaux, noaux=noaux)

!
! Determine the index in the record axe
!
      if ( recordaxe(irq)%nstatm ) then
          READ (chtime, '(I8,1X,I6)') (xpctmp(i), i=1,2)
          xpctime = seconds_since_epoch(xpctmp(1), xpctmp(2))

          ! If provided, check the time against a reference time
          if ( present(xpctime2) ) then
            if ( xpctime2 /= xpctime ) then
                write(errmsg, '(I12, " does not equal ", I12)') xpctime2, xpctime
                call MSGERR(4, errmsg)
                return
            end if
          end if
          call record_index(ncid, recordaxe(irq), xpctime, ri)
      else
          if ( string_is_int(PROJNR) ) then
            read( PROJNR, '(i16)' ) pnr
            if ( pnr == 0 ) then
                call MSGERR(3, 'The projnr number may not be zero')
                return
            end if
          else
            pnr = 1
          end if
          call record_index(ncid, recordaxe(irq), pnr, ri )
      end if

      if ( wetnode_list ) then
          ! remove dry points. In memory there's a copy of spcaux with dry points and
          ! lspcaux without them. If this remains problematic than we need to do something clever
          ! with KGRPGL at write time
          call swn_outnc_spcaux_on_wetnodes(KGRPGL, spcaux, lspcaux )
          if (STPNOW()) return

          ! pointer to wetnode list
          pspcaux => lspcaux

          ! deallocate list with dry points
          call swn_outnc_deallocate_spcaux( spcaux )

      else
          pspcaux => spcaux
      end if

!
! add spectra
!
      where(pspcaux%E < epsilon(1.)) pspcaux%E = 0.
      if ( pspcaux%is2d ) then
          ! factor 2*PI to account for transition from rad/s to Hz
          pspcaux%E = pspcaux%E * 2 * pi
          call agnc_add_spcdata_density(ncid, ri, pspcaux%E, spc_as_map)
      else
          call agnc_add_spcdata(ncid, ri, pspcaux%E(1:MSC,:), &
                                          pspcaux%E(MSC+1:2*MSC,:), &
                                          pspcaux%E(2*MSC+1:3*MSC,:), &
                                          spc_as_map)
      end if
!
! add auxilary variables (depth, current and Hs).
!
!     hs is provided so that user can check his integration routines. His / hers
!     hs should match the one computed by SWAN.
      if ( .not.noaux ) then
        if ( spc_as_map ) then
            call agnc_add_mapdata(ncid, 'depth', ri, pspcaux%depth, OVEXCV( 4), skip_range_error)
            call agnc_add_mapdata(ncid, 'xcur',  ri, pspcaux%ux,    OVEXCV( 5), skip_range_error)
            call agnc_add_mapdata(ncid, 'ycur',  ri, pspcaux%uy,    OVEXCV( 5), skip_range_error)
            call agnc_add_mapdata(ncid, 'hs',    ri, pspcaux%hs,    OVEXCV(10), skip_range_error)
            call agnc_add_mapdata(ncid, 'xwnd',  ri, pspcaux%wndx,  OVEXCV(26), skip_range_error)
            call agnc_add_mapdata(ncid, 'ywnd',  ri, pspcaux%wndy,  OVEXCV(26), skip_range_error)
        else
            call agnc_add_pntdata(ncid, 'depth', ri, pspcaux%depth, OVEXCV( 4), skip_range_error)
            call agnc_add_pntdata(ncid, 'xcur',  ri, pspcaux%ux,    OVEXCV( 5), skip_range_error)
            call agnc_add_pntdata(ncid, 'ycur',  ri, pspcaux%uy,    OVEXCV( 5), skip_range_error)
            call agnc_add_pntdata(ncid, 'hs',    ri, pspcaux%hs,    OVEXCV(10), skip_range_error)
            call agnc_add_pntdata(ncid, 'xwnd',  ri, pspcaux%wndx,  OVEXCV(26), skip_range_error)
            call agnc_add_pntdata(ncid, 'ywnd',  ri, pspcaux%wndy,  OVEXCV(26), skip_range_error)
        end if
    end if

      ! garbage collect
      nullify(pspcaux)
      if ( wetnode_list ) then
        call swn_outnc_deallocate_spcaux( lspcaux )
      else
        call swn_outnc_deallocate_spcaux( spcaux )
      end if

    end subroutine swn_outnc_appendspc

    subroutine swn_outnc_openspecfile(ncfile, spcaux, oqi, oqr)
!
! Structure:
!
!    1) input arguments
!    2) local variables
!    if file does not exist
!      3) Define new record axe (run or time)
!      If spectra requested on COMPGRID
!        if not rotated
!          4.a) define mspcgrid with 1D coordinate axes
!        else
!          4.b) define mspcgrid with 2D coordinate axes
!        end
!        4.c) create ncfile
!      else
!        5.a) define spcgrid
!        if WETCGRID
!          5.b) filter pointlist
!        end if
!        5.b) create ncfile
!      end if
!      6.a) add global attributes
!      6.b) create variables
!      6.c) close define mode and reopen write
!      6.d) fill dimension variables and close
!    end if file does not exist
!    7) open file in write mode
!
! 1) input arguments
!
        character*80,                  intent(   in) :: ncfile
        type(spcaux_type), target,     intent(   in) :: spcaux
        integer,                       intent(   in) :: oqi(4)
        real,                          intent(   in) :: oqr(2)
!
! 2) local variables
!
        logical                               :: file_exists, nc_debug, &
                                                 Escale, monthly, spc_as_map, &
                                                 wetnode_list, STPNOW, noaux
        integer                               :: ncid, xpctmp(2), i, irq
        type(spcgrid_type)                    :: spcgrid
        type(mapgrid_type)                    :: mapgrid
        type(pntgrid_type)                    :: pntgrid
        type(spcaux_type), target             :: lspcaux
        type(spcaux_type), pointer            :: pspcaux
        integer, save                         :: IENT=0

        file_exists  = .false.
        nc_debug     = .false.

        if (LTRACE) call STRACE (IENT,'swn_outnc_openspecfile')

        call swn_outnc_spcflags(oqi(4), spc_as_map, wetnode_list, spcaux, Escale, noaux, monthly)
        irq = OQI(2)

        if ( .not. stnames_initialized ) call stnames_init()

        inquire( FILE=ncfile, EXIST=file_exists )

!
! if file does not exists
!
        if ( .not. file_exists) then
!
! 3) Define new record axe (run or time) and spectral grid
!
            allocate(recordaxe(irq)%content(1))
            recordaxe(irq)%ncontent = 1

            if ( NSTATM > 0 ) then
                READ (chtime, '(I8,1X,I6)') (xpctmp(i), i=1,2)
                ! should be seconds since epoch.
                recordaxe(irq)%content(1) = seconds_since_epoch(xpctmp(1), xpctmp(2))
                recordaxe(irq)%delta      = maxval((/oqr(2), DT/))
                recordaxe(irq)%nstatm     = .true.
            else
                ! stationary, record axe is "run"
                if ( string_is_int(PROJNR) ) then
                    read( PROJNR, '(i16)' ) recordaxe(irq)%content(1)
                else
                    recordaxe(irq)%content(1) = 1
                end if
                recordaxe(irq)%delta      = 1
                recordaxe(irq)%nstatm     = .false.
            end if

            allocate (spcgrid%frequency(spcaux%nfreq))
            spcgrid%nfreq     = spcaux%nfreq
            spcgrid%frequency = spcaux%f

            if ( spcaux%is2d ) then
                spcgrid%ndir  = spcaux%ndir
                allocate (spcgrid%direction(spcaux%ndir))
                spcgrid%direction = spcaux%theta
            else
                spcgrid%ndir  = 0
            end if

            spcgrid%relative  = spcaux%relative
!
! If spectra requested on COMPGRID
!
            if ( spc_as_map .and. OPTG /= 5 ) then
                ! COMPGRID is a special case. By specifying this spectra are 4D
                ! written to the netcdf instead of only the wet points
!
! if not rotated
!
                ! XGRDGL, YGRDGL, MXCGL and MYCGL from M_PARALL
                if ( abs(ALPC) < epsilon(1.) .and. OPTG == 1 ) then
!
! 4.a) define mapgrid with 1D coordinate axes
!
                    allocate (mapgrid%longitude(MXCGL,1))
                    allocate (mapgrid%latitude (1,MYCGL))
                    if ( PARLL ) then
                        mapgrid%longitude(:,1) = xgrdgl(:,1) + XOFFS
                        mapgrid%latitude(1,:)  = ygrdgl(1,:) + YOFFS
                    else
                        mapgrid%longitude(:,1) = xgrdgl(:,1)
                        mapgrid%latitude(1,:)  = ygrdgl(1,:)
                    end if
                    mapgrid%mdc  = .false.
                    mapgrid%alpc = 0.
                else
!
! 4.b) define mapgrid with 2D coordinate axes
!
                    allocate (mapgrid%longitude(MXCGL, MYCGL))
                    allocate (mapgrid%latitude (MXCGL, MYCGL))
                    if ( PARLL ) then
                        mapgrid%longitude = XGRDGL + XOFFS
                        mapgrid%latitude  = YGRDGL + YOFFS
                    else
                        mapgrid%longitude = XGRDGL
                        mapgrid%latitude  = YGRDGL
                    end if
                    mapgrid%mdc  = .true.
                    mapgrid%alpc = ALPC
                end if
                ! curvilinear grids may contain the EXCEPTION value of the bottom
                ! Replace with the hardcoded float _FillValue
                where ( abs(mapgrid%longitude - EXCFLD(1)) < epsilon(1.) ) &
                    mapgrid%longitude = NF90_FILL_FLOAT
                where ( abs(mapgrid%latitude - EXCFLD(1)) < epsilon(1.) ) &
                    mapgrid%latitude = NF90_FILL_FLOAT

                if ( KSPHER == 0 ) mapgrid%lunit = 'meter'

                mapgrid%nx = MXCGL
                mapgrid%ny = MYCGL
!
! 4.c) create ncfile
!

                ! Definition mode
                call create_ncfile(ncfile, ncid, recordaxe(irq), spcgrid=spcgrid, mapgrid=mapgrid, Escale=Escale, nautical=BNAUT)

            else
!
! 5.a) define pntgrid, eg, spectra along a point list
!
                if ( wetnode_list ) then
                    ! remove dry points
                    call swn_outnc_spcaux_on_wetnodes(KGRPGL, spcaux, lspcaux )
                    if (STPNOW()) return
                    pspcaux => lspcaux
                else
                    pspcaux => spcaux
                end if
                allocate (pntgrid%longitude(pspcaux%mip))
                allocate (pntgrid%latitude(pspcaux%mip))
                pntgrid%npoints   = pspcaux%mip
                pntgrid%longitude = pspcaux%xp
                pntgrid%latitude  = pspcaux%yp
                if ( wetnode_list ) then
                    allocate(pntgrid%ips(pspcaux%mip))
                    pntgrid%ips     = pspcaux%ips
                    pntgrid%xdimlen = MXCGL
                    pntgrid%ydimlen = MYCGL
                end if
                nullify(pspcaux)
                if ( wetnode_list ) then
                    call swn_outnc_deallocate_spcaux( lspcaux )
                end if

                if ( KSPHER == 0 ) pntgrid%lunit = 'meter'

!
! 5.b) create point list ncfile
!
                ! Definition mode
                call create_ncfile(ncfile, ncid, recordaxe(irq), spcgrid=spcgrid, pntgrid=pntgrid, Escale=Escale, nautical=BNAUT)

            end if
!
! 6.a) add global attributes
!
            call nccheck ( nf90_put_att( ncid, NF90_GLOBAL, 'project', PROJID) )
            call nccheck ( nf90_put_att( ncid, NF90_GLOBAL, 'model',   VERTXT) )
            if ( NSTATM > 0 ) then
                call nccheck ( nf90_put_att( ncid, NF90_GLOBAL, 'run', PROJNR) )
            end if

!
! 6.b) create variables
!
            if ( .not.noaux ) call create_auxvariables(ncid, spc_as_map)


! 6.c) close define mode and reopen write
!
            ! write mode, e.g, fill dimension variables. Note that record dimension
            ! variable is updated by open_ncfile
            call close_ncfile(ncid)
            call open_ncfile( ncfile, "write", ncid, recordaxe(irq), monthly)
!
! 6.d) fill dimension variables and close
!
            call agnc_set_spcgrid (ncid, spcgrid)
            if ( spc_as_map .and. OPTG /= 5 ) then
                call agnc_set_mapgrid (ncid, mapgrid)
            else
                call agnc_set_pntgrid (ncid, pntgrid)
            end if

            ! close the freshly created file so all dimension variable data is written
            call close_ncfile(ncid)
        end if
!
! 7) open file in write mode
!
        call open_ncfile( ncfile, "write", ncid, recordaxe(irq))

        ! if only only record found, set delta to delta specified in command file (or 1 on stationary)
        if ( recordaxe(irq)%ncontent == 1 ) then
            if ( NSTATM > 0 ) then
                recordaxe(irq)%delta = maxval((/oqr(2), DT/))
            else
                recordaxe(irq)%delta = 1
            end if
        end if
        ncoffset(irq) = ncid - oqi(1)

    end subroutine swn_outnc_openspecfile

    subroutine swn_outnc_spcflags(oqi, spc_as_map, wetnode_list, spcaux, Escale, noaux, monthly)
        integer,                intent(   in) :: oqi
        type(spcaux_type),      intent(   in) :: spcaux
        logical,                intent(  out) :: spc_as_map, wetnode_list
        logical, optional,      intent(  out) :: Escale, monthly, noaux

        integer                               :: iocode

        iocode       = oqi
        wetnode_list = .false.
        spc_as_map   = .false.
        noaux        = .false.

        if ( iocode >= 16 ) then
            noaux  = .true.
            iocode = iocode - 16
        end if

        ! MDGRID is the special keyword to store spectra (including dry points) on the
        ! full geographic grid instead of along a pointlist
        ! On the other hand, when COMPGRID is requested on any grid but unstructured,
        ! store the point indices in the computational grid in the netcdf file
        if ( spcaux%mip == MCGRDGL-1 .or. spcaux%mip == (MXCGL * MYCGL) ) then
            if (iocode >= 8 ) then
                spc_as_map = .true.
                iocode     = iocode - 8
            elseif ( iocode >= 4 .and. OPTG /= 5  ) then
                ! request to remove dry points from COMPGRID pointlist, e.g., compress
                wetnode_list = .true.
            end if
        end if

        ! catch conflicting MDGRID + COMPRESS
        iocode       = mod(iocode, 4)

        ! Enery scaling is of creating a new variable that stores a scale per timestep and point
        ! the energy itself is then stored as int16 instead of float
        if ( present( EScale) ) then
            Escale = .false.
            if ( iocode > 1) Escale = .true.
        end if

         ! By setting the keyword MONTH a netcdf file is created that has a time axis
        ! starting at the first of the month.
        if ( present( monthly ) ) then
            monthly = .false.
            if ( mod(iocode,2) == 1 .and. NSTATM == 1) monthly = .true.
        end if

    end subroutine swn_outnc_spcflags

    subroutine swn_outnc_appendblock(myk, mxk, ivtype, nref, irq, data, excv, col)
        ! note that ivtyp is an array and ivtype is ivtyp(JVAR)
        ! col is used for writing vector data.
        integer,                intent(   in) :: irq, ivtype, nref, col, myk, mxk
        real,                   intent(   in) :: data(mxk * myk), excv

        integer                               :: ri, pnr, xpctmp(2), i, &
                                                 xpctime, ilpos, ncid
        integer, save                         :: IENT=0
        if (LTRACE) call STRACE (IENT,'swn_outnc_appendblock')
        if ( ivtype == 40 ) return

        ncid = nref + ncoffset(irq)

        if ( recordaxe(irq)%nstatm ) then
            READ (chtime, '(I8,1X,I6)') (xpctmp(i), i=1,2)
            xpctime = seconds_since_epoch(xpctmp(1), xpctmp(2))
            call record_index(ncid, recordaxe(irq), xpctime, ri)
        else
            if ( string_is_int(PROJNR) ) then
                read( PROJNR, '(i16)' ) pnr
                if ( pnr == 0 ) then
                    call MSGERR(3, 'The projnr number may not be zero')
                    return
                end if
            else
                pnr = 1
            end if
            call record_index(ncid, recordaxe(irq), pnr, ri )
        end if

        if (myk == 1) then
            call agnc_add_pntdata(ncid, STNAMES(ivtype,col), ri, data, excv, skip_range_error)
        else
            call agnc_add_mapdata(ncid, STNAMES(ivtype,col), ri, data, excv, skip_range_error)
        end if

    end subroutine swn_outnc_appendblock

    subroutine swn_outnc_close_on_end(nref, irq)
        integer,           intent(inout) :: nref
        integer, optional, intent(   in) :: irq
        integer                          :: ncid

        integer, save                    :: IENT=0

        ncid = nref
        if ( present(irq) ) ncid = nref + ncoffset(irq)

        ! Decided to sync data to file even in non stationary mode because
        ! otherwise data is buffered during the entire computation,
        ! increasing the chance on file corruption
        if (LTRACE) call STRACE (IENT,'swn_outnc_close_on_end')
        if ( nstatm == 0 .or. TIMCO >= TFINC ) then
            call close_ncfile(ncid)
            nref = 0
        else
            call nccheck( NF90_SYNC(ncid) )
        end if
    end subroutine swn_outnc_close_on_end

    subroutine record_index(ncid, recordaxe, xpctime, ri)
        integer,                              intent(   in) :: ncid, xpctime
        type(recordaxe_type),                 intent(inout) :: recordaxe
        integer,                              intent(  out) :: ri
        integer, allocatable, dimension(:)                  :: tmp_time
        integer, save                                       :: IENT=0

        if (LTRACE) call STRACE (IENT,'record_index')
        ri = timeindex(recordaxe%content, xpctime)

        if ( ri < 0 ) then
            call MSGERR(4, 'agioncmd cannot prepend data, only overwrite '// &
                           'and append. Consider setting option in order '// &
                           'to generate monthly files. You''ll need to rerun '// &
                           'the first period too though')
            call close_ncfile(ncid)
        end if
        if ( ri == 0) then
            ! record needs to be appended
            allocate(tmp_time(recordaxe%ncontent))
            tmp_time = recordaxe%content
            deallocate(recordaxe%content);
            allocate(recordaxe%content(recordaxe%ncontent+1))
            recordaxe%content(1:recordaxe%ncontent) = tmp_time
            recordaxe%content(recordaxe%ncontent+1) = xpctime
            recordaxe%ncontent = recordaxe%ncontent+1
            !
            ! Update netcdf file with new time axe
            !
            call agnc_set_recordaxe (ncid, recordaxe)
            ri = recordaxe%ncontent

        end if

    end subroutine record_index

    subroutine swn_outnc_openblockfile(ncfile, myk, mxk, ovlnam, xgrdgl,   &
                                       ygrdgl, oqi, oqr, ivtyp ,  irq  )
        character*80,           intent(   in) :: ncfile
        integer,                intent(   in) :: myk, mxk, irq
        integer, dimension(:),  intent(   in) :: oqi, ivtyp
        real,    dimension(:),  intent(   in) :: oqr
        character*40,           intent(   in) :: ovlnam(:)
        real,                   intent(   in) :: xgrdgl(mxk, myk), &
                                                 ygrdgl(mxk, myk)

        ! local variables
        logical                               :: file_exists = .false., &
                                                 nc_debug    = .false., &
                                                 monthly     = .false.
        integer                               :: ncid, xpctmp(2), i
        type(mapgrid_type)                    :: mapgrid
        type(pntgrid_type)                    :: pntgrid
        integer, save                         :: IENT=0

        if (LTRACE) call STRACE (IENT,'swn_outnc_openblockfile')

        if ( .not. stnames_initialized ) call stnames_init()

        ! By setting the idla to 5 a netcdf file is created that has a time axis
        ! that starts at the first of the month. Monthly is only valid on non-stationary
        ! run that has time axis. Declaring a monthly file is required when computing
        ! hindcasts in several chunks in time.
        if ( oqi(4) == 5 .and. nstatm == 1 ) monthly = .true.

        inquire( FILE=ncfile, EXIST=file_exists )

        if ( .not. file_exists) then
            allocate(recordaxe(irq)%content(1))

            recordaxe(irq)%ncontent = 1
            if ( nstatm > 0 ) then
                READ (chtime, '(I8,1X,I6)') (xpctmp(i), i=1,2)
                ! should be seconds since epoch. This is probably not it.
                recordaxe(irq)%content(1) = seconds_since_epoch(xpctmp(1), xpctmp(2))
                recordaxe(irq)%delta      = maxval((/oqr(2), DT/))
                recordaxe(irq)%nstatm     = .true.
            else
                ! stationary, record axe is "run"
                if ( string_is_int(PROJNR) ) then
                    read( PROJNR, '(i16)' ) recordaxe(irq)%content(1)
                else
                    recordaxe(irq)%content(1) = 1
                end if
                recordaxe(irq)%delta      = 1
                recordaxe(irq)%nstatm     = .false.
            end if

            if (myk == 1) then
                ! unstructured mesh is a pointlist. TABLE output also uses this section
                pntgrid%npoints = mxk
                allocate (pntgrid%longitude(mxk))
                allocate (pntgrid%latitude(mxk))
                pntgrid%latitude  = ygrdgl(:,1)
                pntgrid%longitude = xgrdgl(:,1)
                if ( KSPHER == 0 ) pntgrid%lunit = 'meter'

                ! Definition mode
                call create_ncfile(ncfile, ncid, recordaxe(irq), pntgrid=pntgrid, nautical=BNAUT)
                call nccheck ( nf90_put_att( ncid, NF90_GLOBAL, 'project', PROJID) )
                call nccheck ( nf90_put_att( ncid, NF90_GLOBAL, 'run', PROJNR) )
                call create_variables(ncid, ivtyp, ovlnam, .false.)

                ! write mode, e.g, fill dimension variables. Note that record dimension
                ! variable is updated by open_ncfile
                call close_ncfile(ncid)
                call open_ncfile( ncfile, "write", ncid, recordaxe(irq), monthly)
                call agnc_set_pntgrid (ncid, pntgrid)
            else
                ! Regular mesh. Depending on the angle of the computational grid,
                ! longitude and latitude can be a vector or a matrix
                ! curvilinear grids always have multi-dimension coordinate fields

                ! Definition mode
                if ( abs(ALPC) < epsilon(1.) .and. OPTG == 1 ) then
                    allocate (mapgrid%longitude(mxk,1))
                    allocate (mapgrid%latitude (1,myk))
                    if ( PARLL ) then
                        mapgrid%longitude(:,1) = xgrdgl(:,1) + XOFFS
                        mapgrid%latitude(1,:)  = ygrdgl(1,:) + YOFFS
                    else
                        mapgrid%longitude(:,1) = xgrdgl(:,1)
                        mapgrid%latitude(1,:)  = ygrdgl(1,:)
                    end if
                    mapgrid%mdc  = .false.
                else
                    allocate (mapgrid%longitude(mxk, myk))
                    allocate (mapgrid%latitude (mxk, myk))
                    if ( PARLL ) then
                        mapgrid%longitude = xgrdgl + XOFFS
                        mapgrid%latitude  = ygrdgl + YOFFS
                    else
                        mapgrid%longitude = xgrdgl
                        mapgrid%latitude  = ygrdgl
                    end if
                    mapgrid%mdc  = .true.
                end if
                ! curvilinear grids may contain the EXCEPTION value of the bottom
                ! Replace with the hardcoded float _FillValue
                where ( abs(mapgrid%longitude - EXCFLD(1)) < epsilon(1.) ) &
                    mapgrid%longitude = NF90_FILL_FLOAT
                where ( abs(mapgrid%latitude - EXCFLD(1)) < epsilon(1.) ) &
                    mapgrid%latitude = NF90_FILL_FLOAT

                mapgrid%nx   = mxk
                mapgrid%ny   = myk
                mapgrid%alpc = ALPC
                if ( KSPHER == 0 ) mapgrid%lunit = 'meter'

                call create_ncfile(ncfile, ncid, recordaxe(irq), mapgrid=mapgrid, nautical=BNAUT)
                call nccheck ( nf90_put_att( ncid, NF90_GLOBAL, 'project', PROJID) )
                call nccheck ( nf90_put_att( ncid, NF90_GLOBAL, 'run', PROJNR) )
                call create_variables(ncid, ivtyp, ovlnam, .true.)

                ! write mode, e.g, fill dimension variables. Note that record dimension
                ! variable is updated by open_ncfile
                call close_ncfile(ncid)
                call open_ncfile( ncfile, "write", ncid, recordaxe(irq), monthly)
                call agnc_set_mapgrid (ncid, mapgrid)
            end if

            ! close the freshly created file so all buffers are emptied
            call close_ncfile(ncid)
        end if

        call open_ncfile( ncfile, "write", ncid, recordaxe(irq))

        ! set the offset to translate between swan's HIOPEN + IRQ and
        ! the netcdf 65xxx range
        ncoffset(irq) = ncid - oqi(1)

        ! if only only record found, set delta to delta specified in command file (or 1 on stationary)
        if ( recordaxe(irq)%ncontent == 1 ) then
            if ( NSTATM > 0 ) then
                recordaxe(irq)%delta = maxval((/oqr(2), DT/))
            else
                recordaxe(irq)%delta = 1
            end if
        end if

    end subroutine swn_outnc_openblockfile

    subroutine create_auxvariables(ncid, spc_as_map)
        integer,                intent( in) :: ncid
        logical,                intent( in) :: spc_as_map
        integer                             :: ivtyp(4)
        character*40                        :: dummy(4)

        ivtyp = (/4, 5, 10, 26/)
        call create_variables(ncid, ivtyp, dummy, spc_as_map)

    end subroutine create_auxvariables

    subroutine create_variables(ncid, ivtyp, ovlnam, spc_as_map)
        integer,                intent( in) :: ncid
        integer,dimension(:),   intent( in) :: ivtyp
        character*40,           intent( in) :: ovlnam(:)
        logical,                intent( in) :: spc_as_map
        integer                             :: i

        character*40                        :: stname
        integer, save                       :: IENT=0
        if (LTRACE) call STRACE (IENT,'create_variables')

        do i=1, size(ivtyp)
            stname = STNAMES(ivtyp(i),1)
            if ( stname /= '' ) then
                if ( spc_as_map ) then
                    call agnc_define_mapvariable( ncid, stname)
                else
                    call agnc_define_pntvariable( ncid, stname)
                end if

                ! if vector, write second variable
                stname = STNAMES(ivtyp(i),2)
                if ( stname /= '' ) then
                    if ( spc_as_map ) then
                        call agnc_define_mapvariable( ncid, stname)
                   else
                        call agnc_define_pntvariable( ncid, stname)
                    end if
                end if
            else
                ! See list at the end of the source file
                if ( ivtyp(i) > 3) call MSGERR(1, "Field '"//trim(ovlnam(ivtyp(i)))//"' not supported (yet).")
            end if
        end do

    end subroutine create_variables

    subroutine stnames_init()
        STNAMES( 4, 1) = 'depth'
        STNAMES( 5, 1) = 'xcur'
        STNAMES( 5, 2) = 'ycur'
        STNAMES( 6, 1) = 'ubot'
        STNAMES(10, 1) = 'hs'
        STNAMES(11, 1) = 'tm01'
        STNAMES(12, 1) = 'tp'
        STNAMES(13, 1) = 'theta0'
        STNAMES(14, 1) = 'thetap'
        STNAMES(16, 1) = 'spread'
        STNAMES(17, 1) = 'L'
        STNAMES(26, 1) = 'xwnd'
        STNAMES(26, 2) = 'ywnd'
        STNAMES(28, 1) = 'rtm01'
        STNAMES(30, 1) = 'dhs'
        STNAMES(31, 1) = 'dtm'
        STNAMES(32, 1) = 'tm02'
        STNAMES(34, 1) = 'urms'
        STNAMES(35, 1) = 'ustar'
        STNAMES(38, 1) = 'cdrag'
        STNAMES(44, 1) = 'hswe'
        STNAMES(47, 1) = 'tmm10'
        STNAMES(48, 1) = 'rtmm10'
        STNAMES(51, 1) = 'ssh'
        STNAMES(52, 1) = 'botl'
        STNAMES(53, 1) = 'tps'
        stnames_initialized = .true.

    end subroutine stnames_init

    function string_is_int(str) result (is_int)
        character(len=*),  intent( in)  :: str
        logical                         :: is_int
        character(len=10)               :: nums = "0123456789"
        is_int = .true.
        if (  verify(str, nums) > 0 ) is_int = .false.
    end function string_is_int

end module swn_outnc

                !  ( 3) [              'distance along output curve']
                !  ( 7) [                       'Energy dissipation']
                !  ( 8) [                  'Fraction breaking waves']
                !  ( 9) [     'Energy leak over spectral boundaries']
                !  (15) [        'direction of the energy transport']
                !  (18) [                           'Wave steepness']
                !  (19) [                    'Wave energy transport']
                !  (20) [       'Wave driven force per unit surface']
                !  (21) [                  'spectral action density']
                !  (22) [                  'spectral energy density']
                !  (23) [                       'auxiliary variable']
                !  (24) [          'X computational grid coordinate']
                !  (25) [          'Y computational grid coordinate']
                !  (27) [              'Bottom friction coefficient']
                !  (29) [ 'energy density integrated over direction']
                !  (33) [         'Frequency spectral width (Kappa)']
                !  (35) [                        'Friction velocity']
                !  (40) [                                'Date-time']
                !  (41) [      'Time in seconds from reference time']
                !  (36) ['Zero velocity thickness of boundary layer']
                !  (37) [                                     '    ']
                !  (39) [                       'Setup due to waves']
                !  (49) [                    'Diffraction parameter']
                !  (50) [                       'Bottom wave period']
                !  (54) [              'Bottom friction dissipation']
                !  (55) [                'Surf breaking dissipation']
                !  (56) [                 'Whitecapping dissipation']
                !  (57) [                   'Vegetation dissipation']
                !  (58) [                               'Peakedness']
                !  (59) [                      'Benjamin-Feir index']
                !  (60) [                        'Energy generation']
                !  (61) [                         'Wind source term']
                !  (62) [                    'Energy redistribution']
                !  (63) [        'Total absolute 4-wave interaction']
                !  (64) [        'Total absolute 3-wave interaction']
                !  (65) [                       'Energy propagation']
                !  (66) [                           'xy-propagation']
                !  (67) [                        'theta-propagation']
                !  (68) [                        'sigma-propagation']
                !  (69) [                         'Radiation stress']
                !  (70) [                            'Plants per m2']
                !  (71) [                         'Peak wave length']
