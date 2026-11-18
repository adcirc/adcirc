!-------------------------------------------------------------------------------!
!
! ADCIRC - The ADvanced CIRCulation model
! Copyright (C) 1994-2025 R.A. Luettich, Jr., J.J. Westerink
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!-------------------------------------------------------------------------------!
!
!                      ADCIRC - SUBDOMAIN MODELING MODULE
!
!    ========================================================================
!    |                                                                      |
!    |   This file contains the subroutines required by Subdomain Modeling, |
!    |   an approach to reduce the total runtime of a series of hurricane   |
!    |   storm surge simulations on a smaller grid, a subdomain grid,       |
!    |   extracted from the original grid.                                  |
!    |                                                                      |
!    |   Written by Alper Altuntas, aaltunt@ncsu.edu                        |
!    |   North Carolina State University,                                   |
!    |   2013                                                               |
!    |                                                                      |
!    ========================================================================

#include "logging_macros.h"

module mod_subdomain
   use mod_sizes, only: globaldir
   use mod_logging, only: DEBUG, ECHO, INFO, WARNING, ERROR, &
                          logMessage, allMessage, t_log_scope, &
                          init_log_scope
   implicit none

   private

   public :: subdomainOn, NOutGS, enforceBN, nobn, obn, &
             readFort015, openFort019H, readFort019, readFort020, &
             readFort021, writeFort065, writeFort066, writeFort067, &
             enforceEcb, enforceEob, enforceEib, enforceUVcb, &
             enforceUVob, enforceWDcb, enforceWDob, checkChange

   logical :: subdomainOn
   integer :: NOUTGS, NSPOOLGS, enforceBN
   integer :: ncbnr, nobnr, nibnr
   integer, allocatable :: cbnr(:), obnr(:), ibnr(:)
#ifdef CMPI
   integer :: localncbnr, localnobnr, localnibnr
   integer, allocatable :: localcbnr(:), localobnr(:), localibnr(:)
#endif
   integer :: sbtiminc, nobn, nibn, ncbn
   real(8), allocatable :: eobn2(:), uobn2(:), vobn2(:), eibn2(:)
   real(8), allocatable :: ecbn1(:), ucbn1(:), vcbn1(:)
   real(8), allocatable :: ecbn2(:), ucbn2(:), vcbn2(:)
   real(8), allocatable :: setEob(:), setUob(:), setEcb(:), setUcb(:)
   real(8), allocatable :: setVob(:), setVcb(:), setEib(:)
   integer, allocatable :: setWDob(:), setWDcb(:)
   integer, allocatable :: wdobn1(:), wdobn2(:), wdcbn1(:), wdcbn2(:)
   integer, allocatable :: cbn(:), obn(:), ibn(:)
   integer :: bchange

contains

   subroutine readFort015()
      !    ========================================================================
      !      This subroutine reads in additional modeling parameters from fort.015
      !      and initiates the Subdomain Modeling
      !
      !       - NOUTGS: Type of the run:
      !     NOUTGS=0 => full run
      !     NOUTGS=1 => subdomain run (old file formatting)
      !     NOUTGS=2 => subdomain run (new file formatting)
      !       - NSPOOLGS: The number of timesteps at which information is
      !                   written to output files fort.06*
      !       - enforceBN: Boundary enforcing flag.
      !     enforceBN = 0 => no forcing (full domain)
      !     enforceBN = 1 => forcing (subdomain, old file formatting)
      !     enforceBN = 2 => forcing (subdomain, new file formatting)
      !       - nobnr: The number of outer boundary nodes of subdomain grids to be
      !                recorded to fort.06* during a full run.
      !       - obnr(i): i.th outer boundary node to be recorded to fort.065 and
      !                  fort.066
      !       - NIBNR: The number of inner boundary nodes of subdomain grids to be
      !                recorded to fort.065 and fort.066 during a full run.
      !       - IBNR(i): i.th inner boundary node to be recorded to fort.065 and
      !                  fort.066
      !    ========================================================================

      use mod_sizes, only: myproc
      use mod_terminate, only: terminate, ADCIRC_EXIT_FAILURE
#ifdef CMPI
      use global, only: nodes_lg
      use mesh, only: np
#endif
      implicit none
      integer :: i
#ifdef CMPI
      integer :: j
#endif
      logical :: fileFound

      LOG_SCOPE_TRACED("readFort015", SUBDOMAIN_TRACING)

      fileFound = .false.
      inquire (file=trim(globaldir)//'/'//'fort.015', Exist=fileFound)
      if (fileFound .eqv. .false.) then
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                        message='Subdomain modeling was activated but the ' &
                        //'fort.015 file was not found in "' &
                        //trim(globaldir)//'".')
      end if

      if (myproc == 0) print *, "Subdomain Active"
      open (1015, file=trim(globaldir)//'/'//'fort.015')
      read (1015, *) NOUTGS
      read (1015, *) NSPOOLGS
      read (1015, *) enforceBN
      select case (noutgs)
      case (0)
         ! subdomain run (no b.c. recording)
      case (1)
         ! open ocean boundary nodes:
         read (1015, *) ncbnr
         allocate (cbnr(ncbnr))
         do i = 1, ncbnr
            read (1015, *) cbnr(i)
         end do
      case (2)
         ! outer boundary nodes:
         read (1015, *) nobnr
         allocate (obnr(nobnr))
         do i = 1, nobnr
            read (1015, *) obnr(i)
         end do
         ! inner boundary nodes:
         read (1015, *) nibnr
         allocate (ibnr(nibnr))
         do i = 1, nibnr
            read (1015, *) ibnr(i)
         end do
      case default
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                        message='The NOUTGS value is invalid.')
      end select
#ifdef CMPI
! localize global record-node arrays:
      select case (noutgs)
      case (0)
! subdomain run (no b.c. recording)
      case (1)
! open ocean boundary nodes:
         localncbnr = 0
         do i = 1, np
            if (any(cbnr == nodes_lg(i))) then
               localncbnr = localncbnr + 1
            end if
         end do
         allocate (localcbnr(localncbnr))
         j = 1
         do i = 1, np
            if (any(cbnr == nodes_lg(i))) then
               localcbnr(j) = i
               j = j + 1
            end if
         end do
      case (2)
! outer boundary nodes:
         localnobnr = 0
         do i = 1, np
            if (any(obnr == nodes_lg(i))) then
               localnobnr = localnobnr + 1
            end if
         end do
         allocate (localobnr(localnobnr))
         j = 1
         do i = 1, np
            if (any(obnr == nodes_lg(i))) then
               localobnr(j) = i
               j = j + 1
            end if
         end do
! inner boundary nodes:
         localnibnr = 0
         do i = 1, np
            if (any(ibnr == nodes_lg(i))) then
               localnibnr = localnibnr + 1
            end if
         end do
         allocate (localibnr(localnibnr))
         j = 1
         do i = 1, np
            if (any(ibnr == nodes_lg(i))) then
               localibnr(j) = i
               j = j + 1
            end if
         end do
      case default
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                        message='The NOUTGS value is invalid.')
      end select
#endif

   end subroutine readFort015

   subroutine writeFort065(it)

      !     | This subroutine writes ETA2, x/y velocity and wet/dry states of
      !     | open ocean boundary nodes of subdomain grids to fort.065 file,
      !     | during a full run. This subroutine is called within timestep.F

#ifdef CMPI
      use mod_sizes, only: myproc
#endif
      use GLOBAL, only: dt, rnday, eta2, uu2, vv2, nodecode
#ifdef CMPI
      use global, only: nodes_lg
#endif
      implicit none
#ifdef CMPI
      character(6) :: procLoc
      integer :: gn
#endif
      integer, intent(in) :: it
      integer :: i, n

      ! open fort.065 at the first timestep
      if (it == 1) then
#ifdef CMPI
         write (procLoc, '(A2,I4.4)') "PE", myproc
         open (1065, file=procLoc//'/'//'fort.065', status='REPLACE')
         write (1065, *) 'Subdomain Modeling CBN output'
         write (1065, *) nspoolgs, localncbnr, int(rnday*86400.d0/(dt*dble(nspoolgs))), &
            " ! NSPOOLGS,lncbnr,tsteps"
#else
         open (1065, file='fort.065', status='REPLACE')
         write (1065, *) 'Subdomain Modeling CBN output'
         write (1065, *) nspoolgs, ncbnr, int(rnday*86400.d0/(dt*dble(nspoolgs))), &
            " ! NSPOOLGS,ncbnr,tsteps"
#endif
         close (1065)
      end if

      if (mod(it, NSPOOLGS) == 0) then
#ifdef CMPI
         write (procLoc, '(A2,I4.4)') "PE", myproc
         open (1065, file=procLoc//'/'//'fort.065', &
               access='SEQUENTIAL', position='APPEND')
         write (1065, *) it, "! Timestep"
         do i = 1, localncbnr
            n = localcbnr(i)
            gn = nodes_lg(n)
            write (1065, *) gn, eta2(n), uu2(n)
            write (1065, *) vv2(n), nodecode(n)
         end do
#else
         open (1065, file='fort.065', access='SEQUENTIAL', &
               position='APPEND')
         write (1065, *) it, "! Timestep"
         do i = 1, ncbnr
            n = cbnr(i)
            write (1065, *) n, eta2(n), uu2(n)
            write (1065, *) vv2(n), nodecode(n)
         end do
#endif
         close (1065)
      end if

   end subroutine writeFort065

   subroutine writeFort066(it)

      !     | This subroutine writes ETAS, x/y velocity and wet/dry states of
      !     | outer boundary nodes of subdomain grids to fort.066 file,
      !     | during a full run. This subroutine is called within timestep.F

#ifdef CMPI
      use mod_sizes, only: myproc
#endif
      use GLOBAL, only: dt, rnday, etas, uu2, vv2, nodecode
#ifdef CMPI
      use global, only: nodes_lg
#endif
      implicit none
#ifdef CMPI
      character(6) :: procLoc
      integer :: gn
#endif
      integer, intent(in) :: it
      integer :: i, n

      ! open fort.066 at the first timestep
      if (it == 1) then
#ifdef CMPI
         write (procLoc, '(A2,I4.4)') "PE", myproc
         open (1066, file=procLoc//'/'//'fort.066', status='REPLACE')
         write (1066, *) 'Subdomain Modeling OBN output'
         write (1066, *) nspoolgs, localnobnr, int(rnday*86400.d0/(dt*dble(nspoolgs))), &
            " ! NSPOOLSG,lnobnr,tsteps"
#else
         open (1066, file='fort.066', status='REPLACE')
         write (1066, *) 'Subdomain Modeling OBN output'
         write (1066, *) nspoolgs, nobnr, int(rnday*86400.d0/(dt*dble(nspoolgs))), &
            " ! NSPOOLGS,nobnr,tsteps"
#endif
         close (1066)
      end if

      if (mod(it, NSPOOLGS) == 0) then
#ifdef CMPI
         write (procLoc, '(A2,I4.4)') "PE", myproc
         open (1066, file=procLoc//'/'//'fort.066', &
               access='SEQUENTIAL', position='APPEND')
         write (1066, *) it, "! Timestep"
         do i = 1, localnobnr
            n = localobnr(i)
            gn = nodes_lg(n)
            write (1066, *) gn, etas(n), uu2(n)
            write (1066, *) vv2(n), nodecode(n)
         end do
#else
         open (1066, file='fort.066', access='SEQUENTIAL', &
               position='APPEND')
         write (1066, *) it, "! Timestep"
         do i = 1, nobnr
            n = obnr(i)
            write (1066, *) n, etas(n), uu2(n)
            write (1066, *) vv2(n), nodecode(n)
         end do
#endif
         close (1066)
      end if

   end subroutine writeFort066

   subroutine writeFort067(it)

      !     | This subroutine writes ETAS of inner boundary nodes of subdomain
      !     | grids to fort.067 file during a full run. This subroutine is called
      !     | within timestep.F

#ifdef CMPI
      use mod_sizes, only: myproc
#endif
      use global, only: dt, rnday, etas
#ifdef CMPI
      use global, only: nodes_lg
#endif
      implicit none
#ifdef CMPI
      character(6) :: procLoc
      integer :: gn
#endif
      integer, intent(in) :: it
      integer :: i, n

      ! open fort.067 at the first timestep
      if (it == 1) then
#ifdef CMPI
         write (procLoc, '(A2,I4.4)') "PE", myproc
         open (1067, file=procLoc//'/'//'fort.067', status='REPLACE')
         write (1067, *) 'Subdomain Modeling IBN output'
         write (1067, *) nspoolgs, localnibnr, int(rnday*86400.d0/(dt*dble(nspoolgs))), &
            " ! NSPOOLGS,lnibnr,tsteps"
#else
         open (1067, file='fort.067', status='REPLACE')
         write (1067, *) 'Subdomain Modeling IBN output'
         write (1067, *) nspoolgs, nibnr, int(rnday*86400.d0/(dt*dble(nspoolgs))), &
            " ! NSPOOLGS,nibnr,tsteps"
#endif
         close (1067)
      end if

      if (mod(it, NSPOOLGS) == 0) then
#ifdef CMPI
         write (procLoc, '(A2,I4.4)') "PE", myproc
         open (1067, file=procLoc//'/'//'fort.067', &
               access='SEQUENTIAL', position='APPEND')
         write (1067, *) it, "! Timestep"
         do i = 1, localnibnr
            n = localibnr(i)
            gn = nodes_lg(n)
            write (1067, *) gn, etas(n)
         end do
#else
         open (1067, file='fort.067', access='SEQUENTIAL', &
               position='APPEND')
         write (1067, *) it, "! Timestep"
         do i = 1, nibnr
            n = ibnr(i)
            write (1067, *) n, etas(n)
         end do
#endif
         close (1067)
      end if

   end subroutine writeFort067

   subroutine openFort019H()
      !     | This subroutine opens fort.019 for a hot-started model.
      !     | fort.019 is read from its first record; the hot-start record skip
      !     | was removed in v51.45.

      use MOD_SIZES, only: myproc
      use BOUNDARIES, only: neta
      use mod_terminate, only: terminate, ADCIRC_EXIT_FAILURE
      implicit none
#ifdef CMPI
      character(6) :: procLoc
#endif
      integer :: i, n
      character(1024) :: scratchMessage
#ifdef CMPI
      write (procLoc, '(A2,I4.4)') "PE", myproc
      open (1019, file=procLoc//'/'//'fort.019')
#else
      open (1019, file=trim(globaldir)//'/'//'fort.019')
#endif
      read (1019, *) ! Header
      read (1019, *) sbtiminc, ncbn
      allocate (cbn(ncbn))
      do i = 1, ncbn
         read (1019, *) cbn(i)
      end do
      if (abs(ncbn - neta) > 1) then
         write (scratchMessage, '(A,I0,A,I0,A,I0)') &
            'ncbn != neta: ncbn=', ncbn, ', neta=', neta, &
            ', myproc=', myproc
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                        message=trim(scratchMessage))
      end if
      allocate (ecbn1(ncbn), ucbn1(ncbn), vcbn1(ncbn), wdcbn1(ncbn))
      allocate (ecbn2(ncbn), ucbn2(ncbn), vcbn2(ncbn), wdcbn2(ncbn))
      allocate (setEcb(ncbn), setUcb(ncbn))
      allocate (setVcb(ncbn), setWDcb(ncbn))
      do i = 1, ncbn
         ecbn1(i) = 0.d0
         ucbn1(i) = 0.d0
         vcbn1(i) = 0.d0
         wdcbn1(i) = 0
      end do
      ! Read the first set of boundary conditions:
      read (1019, *) n
      do i = 1, ncbn
         read (1019, *) n, ecbn2(i), ucbn2(i)
         read (1019, *) vcbn2(i), wdcbn2(i)
         wdcbn1(i) = wdcbn2(i) !ensure that the wet nodes begin wet
      end do

   end subroutine openFort019H

   subroutine openFort019C()
      !     | This subroutine opens fort.019 for a cold-started model

      use MOD_SIZES, only: myproc
      use BOUNDARIES, only: neta
      use mod_terminate, only: terminate, ADCIRC_EXIT_FAILURE
      implicit none
#ifdef CMPI
      character(6) :: procLoc
#endif
      integer :: i, n
      character(1024) :: scratchMessage

#ifdef CMPI
      write (procLoc, '(A2,I4.4)') "PE", myproc
      open (1019, file=procLoc//'/'//'fort.019')
#else
      open (1019, file=trim(globaldir)//'/'//'fort.019')
#endif
      read (1019, *) ! Header
      read (1019, *) sbtiminc, ncbn
      allocate (cbn(ncbn))
      do i = 1, ncbn
         read (1019, *) cbn(i)
      end do
      if (abs(ncbn - neta) > 1) then
         write (scratchMessage, '(A,I0,A,I0,A,I0)') &
            'ncbn != neta: ncbn=', ncbn, ', neta=', neta, &
            ', myproc=', myproc
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                        message=trim(scratchMessage))
      end if
      allocate (ecbn1(ncbn), ucbn1(ncbn), vcbn1(ncbn), wdcbn1(ncbn))
      allocate (ecbn2(ncbn), ucbn2(ncbn), vcbn2(ncbn), wdcbn2(ncbn))
      allocate (setEcb(ncbn), setUcb(ncbn))
      allocate (setVcb(ncbn), setWDcb(ncbn))
      do i = 1, ncbn
         ecbn1(i) = 0.d0
         ucbn1(i) = 0.d0
         vcbn1(i) = 0.d0
         wdcbn1(i) = 0
      end do

      ! Read the first set of boundary conditions:
      read (1019, *) n
      do i = 1, ncbn
         read (1019, *) n, ecbn2(i), ucbn2(i)
         read (1019, *) vcbn2(i), wdcbn2(i)
         wdcbn1(i) = wdcbn2(i) !ensure that the wet nodes begin wet
      end do

   end subroutine openFort019C

   subroutine readFort019(it)
      !     | This subroutine reads in ETA2, x/y velocities and wet/dry states of
      !     | outer boundary nodes of the subdomain grid from fort.019 file during
      !     | a subdomain run. This subroutine is called within timestep.F.

      use GLOBAL, only: ihot
      implicit none
      integer, intent(in) :: it
      double precision :: rateTS
      integer :: i, n

      if (it == 1 .and. ihot == 0) then ! coldstart
         call openFort019C()
      end if

      if (mod(it, sbtiminc) == 0) then
         read (1019, *) n
         !             print *, "Reading fort.019 at timestep",n
         do i = 1, ncbn
            ecbn1(i) = ecbn2(i)
            ucbn1(i) = ucbn2(i)
            vcbn1(i) = vcbn2(i)
            wdcbn1(i) = wdcbn2(i)
            read (1019, *) n, ecbn2(i), ucbn2(i)
            read (1019, *) vcbn2(i), wdcbn2(i)
         end do
      end if

      ! Iteration
      rateTS = dble(mod(it, sbtiminc))/dble(sbtiminc)
      do i = 1, ncbn
         setEcb(i) = ecbn1(i) + (ecbn2(i) - ecbn1(i))*ratets
         setUcb(i) = ucbn1(i) + (ucbn2(i) - ucbn1(i))*ratets
         setVcb(i) = vcbn1(i) + (vcbn2(i) - vcbn1(i))*ratets
         setWDcb(i) = wdcbn1(i)
      end do

   end subroutine readFort019

   subroutine readFort020(it)

      !     | This subroutine reads in ETAS, x/y velocities and wet/dry states of
      !     | outer boundary nodes of the subdomain grid from fort.020 file during
      !     | a subdomain run. This subroutine is called within timestep.F.

#ifdef CMPI
      use mod_sizes, only: myproc
#endif
      implicit none
      integer, intent(in) :: it
#ifdef CMPI
      character(6) :: procLoc
#endif
      integer :: i, n

      if (it == 1) then
#ifdef CMPI
         write (procLoc, '(A2,I4.4)') "PE", myproc
         open (1020, file=procLoc//'/'//'fort.020')
#else
         open (1020, file=trim(globaldir)//'/'//'fort.020')
#endif
         read (1020, *) ! Header
         read (1020, *) sbtiminc, nobn
         allocate (obn(nobn))
         do i = 1, nobn
            read (1020, *) obn(i)
         end do
         allocate (eobn2(nobn), uobn2(nobn), vobn2(nobn))
         allocate (wdobn1(nobn), wdobn2(nobn))
         allocate (setEob(nobn), setUob(nobn))
         allocate (setVob(nobn), setWDob(nobn))
      end if

      !if (mod(it,sbtiminc).eq.0) then
      read (1020, *) n
      do i = 1, nobn
         wdobn1(i) = wdobn2(i)
         read (1020, *) n, eobn2(i), uobn2(i)
         read (1020, *) vobn2(i), wdobn2(i)
         setEob(i) = eobn2(i)
         setUob(i) = uobn2(i)
         setVob(i) = vobn2(i)
         setWDob(i) = wdobn2(i)
      end do
      !endif

   end subroutine readFort020

   subroutine readFort021(it)

      !     | This subroutine reads in ETAS of inner boundary nodes of the
      !     | subdomain grid from fort.021 file during a subdomain run.
      !     | This subroutine is called within timestep.F.

#ifdef CMPI
      use mod_sizes, only: myproc
#endif
      implicit none
      integer, intent(in) :: it
#ifdef CMPI
      character(6) :: procLoc
#endif
      integer :: i, n

      if (it == 1) then
#ifdef CMPI
         write (procLoc, '(A2,I4.4)') "PE", myproc
         open (1021, file=procLoc//'/'//'fort.021')
#else
         open (1021, file=trim(globaldir)//'/'//'fort.021')
#endif
         read (1021, *) ! Header
         read (1021, *) sbtiminc, nibn
         allocate (ibn(nibn))
         do i = 1, nibn
            read (1021, *) ibn(i)
         end do
         allocate (eibn2(nibn), setEib(nibn))
      end if

      !         if (mod(it,sbtiminc).eq.0) then
      read (1021, *) n
      !             print *, "Reading fort.021 at timestep",n
      do i = 1, nibn
         read (1021, *) n, eibn2(i)
         setEib(i) = eibn2(i)
      end do
      !         endif

   end subroutine readFort021

   subroutine enforceEcb()

      !     | This subroutine enforces ETA2 of opean ocean boundary nodes
      !     | using setEcb array. Called within timestep.F

      use GLOBAL, only: eta2
      implicit none
      integer :: i, n

      do i = 1, ncbn
         n = cbn(i)
         ETA2(n) = setEcb(i)
      end do

   end subroutine enforceEcb

   subroutine enforceEob()

      !     | This subroutine enforces ETAS of outer boundary nodes
      !     | using setEob array. Called within timestep.F

      use GLOBAL, only: etas
      implicit none
      integer :: i, n

      do i = 1, nobn
         n = obn(i)
         ETAS(n) = setEob(i)
      end do

   end subroutine enforceEob

   subroutine enforceUVcb()

      !     | This subroutine enforces U and V velocities of open ocean boundary
      !     | nodes using setUcb and setVcb arrays. Called within timestep.F

      use GLOBAL, only: uu2, vv2
      implicit none
      integer :: i, n

      do i = 1, ncbn
         n = cbn(i)
         UU2(n) = setUcb(i)
         VV2(n) = setVcb(i)
      end do

   end subroutine enforceUVcb

   subroutine enforceUVob()

      !     | This subroutine enforces U and V velocities of outer boundary
      !     | nodes using setUob and setVob arrays. Called within timestep.F

      use GLOBAL, only: uu2, vv2
      implicit none
      integer :: i, n

      do i = 1, nobn
         n = obn(i)
         UU2(n) = setUob(i)
         VV2(n) = setVob(i)
      end do

   end subroutine enforceUVob

   subroutine enforceWDcb()

      !     | This subroutine enforces wet/dry flags of outer boundary
      !     | nodes using setWDcb. Called within timestep.F

      use GLOBAL, only: nnodecode
      implicit none
      integer :: i, n

      do i = 1, ncbn
         n = cbn(i)
         NNODECODE(n) = setWDcb(i)
      end do

   end subroutine enforceWDcb

   subroutine enforceWDob()

      !     | This subroutine enforces wet/dry flags of outer boundary
      !     | nodes using setWDob. Called within timestep.F

      use GLOBAL, only: nnodecode
      implicit none
      integer :: i, n

      do i = 1, nobn
         n = obn(i)
         NNODECODE(n) = setWDob(i)
      end do

   end subroutine enforceWDob

   subroutine enforceEib()

      !     | This subroutine enforces ETAS of inner boundary nodes
      !     | using setEib array. Called within timestep.F

      use GLOBAL, only: etas
      implicit none
      integer :: i, n

      do i = 1, nibn
         n = ibn(i)
         ETAS(n) = setEib(i)
      end do

   end subroutine enforceEib

   subroutine checkChange()

      !     | This subroutine checks if w/d status of any o.b.n. has changed.
      !     | If a change has detected, ncchange is set to 1 at the next timestep
      !     | to recalculate COEF of outer-inner boundary nodes

      use global, only: ncchange
#ifdef CMPI
      use global, only: ilump
      use messenger, only: WetDrySum
#endif
      implicit none
      integer :: i

      if (bchange == 1) then
         ncchange = 1
         bchange = 0
      end if
      do i = 1, nobn
         if (wdobn2(i) /= wdobn1(i)) then
            bchange = 1
         end if
      end do
#ifdef CMPI
      if (ILump == 0) then
         call WetDrySum(NCCHANGE)
      end if
#endif

   end subroutine checkChange

end module mod_subdomain
