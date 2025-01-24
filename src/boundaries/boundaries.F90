!-------------------------------------------------------------------------------!
!
! ADCIRC - The ADvanced CIRCulation model
! Copyright (C) 1994-2023 R.A. Luettich, Jr., J.J. Westerink
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
! ------------------------------------------------------------------
!              M O D U L E   B O U N D A R I E S
! ------------------------------------------------------------------
! jgf51.21.10: Created this boundaries module in order to modularize
! boundary definition data and related subroutines. This module also
! introduces new data structures that facilitate the reading and
! writing of boundaries in other formats, such as NetCDF and XDMF.
!
! The variables and subroutines in this module were refactored
! out of the other parts of the code, particularly from the global
! module.
! ------------------------------------------------------------------
module boundaries

   implicit none

   real(8), allocatable :: BNDLEN2O3(:)
   real(8), allocatable :: CSII(:), SIII(:)
   real(8), allocatable :: CSIICN(:), SIIICN(:) ! normal direction at condensed nodes 09/22/2022 SB
   integer, allocatable :: ME2GW(:)
   integer, allocatable :: NBV(:), LBCODEI(:)

   integer, allocatable :: NVDLL(:), NBD(:)
   integer, allocatable :: NBDV(:, :)
   integer, allocatable :: NVELL(:)
   integer, allocatable :: NBVV(:, :)
   integer, allocatable :: IBTYPE(:)
   integer, allocatable :: IBTYPEE(:)
   integer, allocatable :: ibtype_orig(:)
   integer, allocatable :: NEleZNG(:)
   real(8), allocatable :: ZNGIF1(:), ZNGIF2(:), ZNGIF3(:)
   !.....for internal barrier boundaries with flowthrough pipes
   real(8), allocatable :: BARLANHT(:), BARLANCFSP(:)
   real(8), allocatable :: BARINHT(:), BARINCFSB(:), BARINCFSP(:)
   real(8), allocatable :: PIPEHT(:), PIPECOEF(:), PIPEDIAM(:)
   integer, allocatable :: IBCONN(:)
   integer, allocatable :: ISSUBMERGED64(:) ! used to identify submerged ibtype=64 weir  SB
   integer, allocatable :: ISSUBMERGED64P(:)
   ! these arrays will only be needed temporarily, between reading
   ! the ascii mesh file and initializing the boundaries
   real(8), allocatable :: BARLANHTR(:, :), BARLANCFSPR(:, :)
   real(8), allocatable :: BARINHTR(:, :), BARINCFSBR(:, :), BARINCFSPR(:, :)
   real(8), allocatable :: PIPEHTR(:, :), PIPECOEFR(:, :), PIPEDIAMR(:, :)
   integer, allocatable :: IBCONNR(:, :)

   integer NOPE
   ! kmd - added for rivers in baroclinic simulation
   integer, allocatable :: BCRNBVV(:, :)
   integer, allocatable :: BCRNVELL(:)
   integer :: totalbcrivernodes ! nodes for the river boundary conditions for baroclinic
   integer :: bcriversec ! river boundary information
   logical :: BndBCRiver = .false.
   logical :: NPEBC = .false. !WJP 03.272018 Non-periodic elevation boundary condition flag

   integer :: NBOU
   integer :: NETA
   integer :: NBVI, NBVJ

   integer :: NVEL, NVELEXT, NVELME
   integer :: IM, JGW, JKI, JME
   integer :: NFLUXIBP, NPIPE
   integer :: nfluxb, nfluxf, nfluxib, nfluxrbc, nfluxgbc
   integer :: nfluxib64, nfluxib64_gbl !sb flag to indicate that ibtype=64 (vertical element wall) boundary exists
   real(8) :: costset
   real(8) :: anginn

   real(8), allocatable :: BTRAN3(:), BTRAN4(:), BTRAN5(:)
   real(8), allocatable :: BTRAN6(:), BTRAN7(:), BTRAN8(:)
   integer, allocatable :: NTRAN1(:), NTRAN2(:)

   ! elevation boundaries and flux boundaries where
   ! ibtype = 0,1,2,10,11,12,20,21,22,30,52
   type simpleBoundary_t
      integer :: indexNum ! order within the fort.14 file (used to get IBTYPE etc)
      integer :: informationID ! xdmf ID for IBTYPEE or IBTYPE info
      integer :: setID ! xdmf ID for node numbers
      integer, allocatable :: nodes(:) ! node numbers on boundary
      integer, allocatable :: xdmf_nodes(:) ! 0 offset node numbers on boundary
      real(8), allocatable :: bGeom(:) ! coordinates for visualization
   end type simpleBoundary_t
   ! variable holding elevation boundaries
   type(simpleBoundary_t), allocatable :: elevationBoundaries(:)

   ! variable holding flux boundaries with
   ! ibtype = 0, 1, 2, 10, 11, 12, 20, 21, 22, 30, 52, 102, 112, 122, 152
   type(simpleBoundary_t), allocatable :: simpleFluxBoundaries(:)
   integer :: numSimpleFluxBoundaries ! for memory allocation
   integer :: sfCount ! index into the simpleFluxBoundaries array

   ! flux boundaries where ibtype = 3, 13, 23
   type externalFluxBoundary_t
      integer :: indexNum ! order within the fort.14 file
      integer :: informationID ! xdmf ID for IBTYPE info
      integer :: setID ! xdmf ID for node numbers
      integer :: numAttributes = 2
      integer :: attributeIDs(2) ! xdmf IDs for parameters
      integer, allocatable :: nodes(:)
      real(8), allocatable :: barlanht(:)
      real(8), allocatable :: barlancfsp(:)
      real(8), allocatable :: leveeGeom(:)
      integer, allocatable :: xdmf_nodes(:)
   end type externalFluxBoundary_t
   type(externalFluxBoundary_t), allocatable :: externalFluxBoundaries(:)
   integer :: numExternalFluxBoundaries
   integer :: efCount ! index into the externalFluxBoundaries array

   ! flux boundaries where ibtype = 4, 24, 64
   type internalFluxBoundary_t
      integer :: indexNum ! order within the fort.14 file
      integer :: informationID ! xdmf ID for IBTYPE info
      integer :: setID ! xdmf ID for node numbers
      integer :: numAttributes = 4
      integer :: attributeIDs(4) ! xdmf IDs for parameters
      integer, allocatable :: nodes(:)
      integer, allocatable :: ibconn(:)
      real(8), allocatable :: barinht(:)
      real(8), allocatable :: barincfsb(:)
      real(8), allocatable :: barincfsp(:)
      real(8), allocatable :: leveeGeom(:)
      integer, allocatable :: xdmf_nodes(:)
      integer, allocatable :: xdmf_ibconn(:)
      integer, allocatable :: ibTypeAttribute(:) ! used for visualization
   end type internalFluxBoundary_t
   type(internalFluxBoundary_t), allocatable :: internalFluxBoundaries(:)
   integer :: numInternalFluxBoundaries
   integer :: ifCount ! index into the internalFluxBoundaries array

   ! flux boundaries where ibtype = 5, 25
   type internalFluxBoundaryWithPipes_t
      integer :: indexNum ! order within the fort.14 file
      integer :: informationID ! xdmf ID for IBTYPE info
      integer :: setID ! xdmf ID for node numbers
      integer :: numAttributes = 7
      integer :: attributeIDs(7) ! xdmf IDs for parameters
      integer, allocatable :: nodes(:)
      integer, allocatable :: ibconn(:)
      real(8), allocatable :: barinht(:)
      real(8), allocatable :: barincfsb(:)
      real(8), allocatable :: barincfsp(:)
      real(8), allocatable :: pipeht(:)
      real(8), allocatable :: pipecoef(:)
      real(8), allocatable :: pipediam(:)
      real(8), allocatable :: leveeGeom(:)
      integer, allocatable :: xdmf_nodes(:)
      integer, allocatable :: xdmf_ibconn(:)
   end type internalFluxBoundaryWithPipes_t
   type(internalFluxBoundaryWithPipes_t), allocatable :: internalFluxBoundariesWithPipes(:)
   integer :: numInternalFluxBoundariesWithPipes
   integer :: ifwpCount ! index into the internalFluxBoundariesWithPipes array

   integer, parameter :: specifiedFluxBoundaryTypes(5) = (/2, 12, 22, 32, 52/)

   ! all info needed for self describing dataset in XDMF
   type xdmfMetaData_t
      logical :: createdIDs ! .true. if infoIDs have been created
      character(80) :: variable_name
      integer :: variable_name_id
      character(80) :: long_name
      integer :: long_name_id
      character(80) :: standard_name
      integer :: standard_name_id
      character(80) :: coordinates
      integer :: coordinates_id
      character(80) :: units
      integer :: units_id
      character(80) :: positive
      integer :: positive_id
      integer :: ndset ! number of data sets
   end type xdmfMetaData_t

contains

   !------------------------------------------------------------------
   !                   S U B R O U T I N E
   ! A L L O C A T E   E L E V A T I O N   B O U N D A R Y   L E N G T H S
   !------------------------------------------------------------------
   ! Allocate the arrays that hold the number of nodes on each elevation
   ! boundary segment
   !------------------------------------------------------------------
   subroutine allocateElevationBoundaryLengths()
      use sizes, only: mnope
      implicit none
      allocate (nvdll(mnope)) ! number of nodes on each elevation boundary segment
      allocate (ibtypee(mnope)) ! type of each elevation boundary segment

      ! initialize to something troublesome to make it easy to spot issues
      ibtypee = -99999
      nvdll = -99999

   end subroutine allocateElevationBoundaryLengths
   !------------------------------------------------------------------

   !------------------------------------------------------------------
   !                   S U B R O U T I N E
   !      A L L O C A T E  F L U X  B O U N D A R Y   L E N G T H S
   !------------------------------------------------------------------
   ! Allocate the arrays that hold the number of nodes (primary nodes
   ! in the case of paired node boundaries like levees) on each flux
   ! boundary segment
   !------------------------------------------------------------------
   subroutine allocateFluxBoundaryLengths()
      use sizes, only: mnbou
      implicit none
      allocate (nvell(mnbou)) ! number of nodes on each flux boundary segment
      allocate (ibtype_orig(mnbou))
      allocate (ibtype(mnbou))

      ! initialize to something troublesome to make it easy to spot issues
      nvell = -99999
      ibtype_orig = -99999
      ibtype = -99999

   end subroutine allocateFluxBoundaryLengths
   !------------------------------------------------------------------

   !------------------------------------------------------------------
   !                   S U B R O U T I N E
   !      A L L O C A T E   A D C I R C  E L E V A T I O N
   !               B O U N D A R Y  A R R A Y S
   !------------------------------------------------------------------
   ! Allocate space for elevation boundary-related variables
   !------------------------------------------------------------------
   subroutine allocateAdcircElevationBoundaryArrays()
      use sizes, only: mnope, mneta
      implicit none
      allocate (nbdv(mnope, mneta))
      allocate (nbd(mneta))

      ! initialize to something troublesome to make it easy to spot issues
      nbdv = -99999
      nbd = -99999

   end subroutine allocateAdcircElevationBoundaryArrays
   !------------------------------------------------------------------

   !------------------------------------------------------------------
   !                       S U B R O U T I N E
   !              A L L O C A T E  A D C I R C   F L U X
   !                   B O U N D A R Y  A R R A Y S
   !------------------------------------------------------------------
   !     jgf51.21.11 Allocate space for flux boundary-related variables
   !------------------------------------------------------------------
   subroutine allocateAdcircFluxBoundaryArrays()
      use sizes, only: mnbou, mnvel
      implicit none
      allocate (me2gw(mnvel))
      allocate (csii(mnvel), siii(mnvel))
      allocate (nbv(mnvel), lbcodei(mnvel))
      allocate (barlanht(mnvel), barlancfsp(mnvel))
      allocate (barinht(mnvel), barincfsb(mnvel), barincfsp(mnvel))
      allocate (pipeht(mnvel), pipecoef(mnvel), pipediam(mnvel))
      allocate (ibconn(mnvel))
      allocate (nbvv(mnbou, 0:mnvel))
      allocate (issubmerged64(mnvel), issubmerged64p(mnvel))
      allocate (bndlen2o3(mnvel))
      allocate (ntran1(mnvel), ntran2(mnvel))
      allocate (btran3(mnvel), btran4(mnvel), btran5(mnvel))
      allocate (btran6(mnvel), btran7(mnvel), btran8(mnvel))
      allocate (nelezng(mnvel))
      allocate (zngif1(mnvel), zngif2(mnvel), zngif3(mnvel))

      ! kmd - added for rivers in baroclinic simulation
      allocate (bcrnbvv(mnbou, 0:mnvel))
      allocate (bcrnvell(mnbou))
      !
      ! initialize to something troublesome to make it easy to spot issues
      nbv = -99999
      lbcodei = -99999
      barlanht = -99999.d0
      barlancfsp = -99999.d0
      barinht = -99999.d0
      barincfsb = -99999.d0
      barincfsp = -99999.d0
      pipeht = -99999.d0
      pipecoef = -99999.d0
      pipediam = -99999.d0
      ibconn = -99999
      nbvv = -99999
      bndlen2o3 = -99999.0
      bcrnbvv = -99999
      bcrnvell = -99999

   end subroutine allocateAdcircFluxBoundaryArrays
   !------------------------------------------------------------------

   !------------------------------------------------------------------
   !                   S U B R O U T I N E
   !        A L L O C A T E  B O U N D A R Y  A R R A Y S
   !------------------------------------------------------------------
   ! Allocate space for boundary-related variables
   !------------------------------------------------------------------
   subroutine allocateBoundaryArrays()
      use sizes, only: mnope
      use global, only: nabout
      implicit none
      integer :: i

      allocate (elevationBoundaries(mnope))
      do i = 1, nope
         allocate (elevationBoundaries(i)%nodes(nvdll(i)))
         allocate (elevationBoundaries(i)%xdmf_nodes(nvdll(i)))
      end do
      allocate (simpleFluxBoundaries(numSimpleFluxBoundaries))
      allocate (externalFluxBoundaries(numExternalFluxBoundaries))
      allocate (internalFluxBoundaries(numInternalFluxBoundaries))
      allocate (internalFluxBoundariesWithPipes(numInternalFluxBoundariesWithPipes))
      sfCount = 1
      efCount = 1
      ifCount = 1
      ifwpCount = 1
      do i = 1, nbou
         select case (ibtype_orig(i))
         case (0, 1, 2, 10, 11, 12, 20, 21, 22, 30, 32, 52, 94)
            allocate (simpleFluxBoundaries(sfCount)%nodes(nvell(i)))
            allocate (simpleFluxBoundaries(sfCount)%xdmf_nodes(nvell(i)))
            sfCount = sfCount + 1
         case (3, 13, 23)
            allocate (externalFluxBoundaries(efCount)%nodes(nvell(i)))
            allocate (externalFluxBoundaries(efCount)%barlanht(nvell(i)))
            allocate (externalFluxBoundaries(efCount)%barlancfsp(nvell(i)))
            allocate (externalFluxBoundaries(efCount)%xdmf_nodes(nvell(i)))
            efCount = efCount + 1
         case (4, 24, 64)
            allocate (internalFluxBoundaries(ifCount)%nodes(nvell(i)))
            allocate (internalFluxBoundaries(ifCount)%ibconn(nvell(i)))
            allocate (internalFluxBoundaries(ifCount)%barinht(nvell(i)))
            allocate (internalFluxBoundaries(ifCount)%barincfsb(nvell(i)))
            allocate (internalFluxBoundaries(ifCount)%barincfsp(nvell(i)))
            allocate (internalFluxBoundaries(ifCount)%xdmf_nodes(nvell(i)))
            allocate (internalFluxBoundaries(ifCount)%xdmf_ibconn(nvell(i)))
            ifCount = ifCount + 1
         case (5, 25)
            allocate (internalFluxBoundariesWithPipes(ifwpCount)%nodes(nvell(i)))
            allocate (internalFluxBoundariesWithPipes(ifwpCount)%ibconn(nvell(i)))
            allocate (internalFluxBoundariesWithPipes(ifwpCount)%barinht(nvell(i)))
            allocate (internalFluxBoundariesWithPipes(ifwpCount)%barincfsb(nvell(i)))
            allocate (internalFluxBoundariesWithPipes(ifwpCount)%barincfsp(nvell(i)))
            allocate (internalFluxBoundariesWithPipes(ifwpCount)%pipeht(nvell(i)))
            allocate (internalFluxBoundariesWithPipes(ifwpCount)%pipecoef(nvell(i)))
            allocate (internalFluxBoundariesWithPipes(ifwpCount)%pipediam(nvell(i)))
            allocate (internalFluxBoundariesWithPipes(ifwpCount)%xdmf_nodes(nvell(i)))
            allocate ( &
               internalFluxBoundariesWithPipes(ifwpCount)%xdmf_ibconn(nvell(i)))
            ifwpCount = ifwpCount + 1
         case default
            write (6, '(a,i0,a)') "ERROR: The boundary type ", ibtype_orig(i), &
               " was found in the file but is not valid."
            call exit(1)
         end select
      end do
      ! initialize to something troublesome to make it easy to spot issues
      do i = 1, nope
         elevationBoundaries(i)%nodes(:) = -99999
         elevationBoundaries(i)%xdmf_nodes(:) = -99999
      end do
      do i = 1, numSimpleFluxBoundaries
         simpleFluxBoundaries(i)%nodes(:) = -99999
         simpleFluxBoundaries(i)%xdmf_nodes(:) = -99999
      end do
      do i = 1, numExternalFluxBoundaries
         externalFluxBoundaries(i)%nodes(:) = -99999
         externalFluxBoundaries(i)%xdmf_nodes(:) = -99999
         externalFluxBoundaries(i)%barlanht(:) = -99999.d0
         externalFluxBoundaries(i)%barlancfsp(:) = -99999.d0
      end do
      do i = 1, numInternalFluxBoundaries
         internalFluxBoundaries(i)%nodes(:) = -99999
         internalFluxBoundaries(i)%ibconn(:) = -99999
         internalFluxBoundaries(i)%barinht(:) = -99999.d0
         internalFluxBoundaries(i)%barincfsb(:) = -99999.d0
         internalFluxBoundaries(i)%barincfsp(:) = -99999.d0
         internalFluxBoundaries(i)%xdmf_nodes(:) = -99999
         internalFluxBoundaries(i)%xdmf_ibconn(:) = -99999
      end do
      do i = 1, numInternalFluxBoundariesWithPipes
         internalFluxBoundariesWithPipes(i)%nodes(:) = -99999
         internalFluxBoundariesWithPipes(i)%ibconn(:) = -99999
         internalFluxBoundariesWithPipes(i)%barinht(:) = -99999.d0
         internalFluxBoundariesWithPipes(i)%barincfsb(:) = -99999.d0
         internalFluxBoundariesWithPipes(i)%barincfsp(:) = -99999.d0
         internalFluxBoundariesWithPipes(i)%pipeht(:) = -99999.d0
         internalFluxBoundariesWithPipes(i)%pipecoef(:) = -99999.d0
         internalFluxBoundariesWithPipes(i)%pipediam(:) = -99999.d0
         internalFluxBoundariesWithPipes(i)%xdmf_nodes(:) = -99999
         internalFluxBoundariesWithPipes(i)%xdmf_ibconn(:) = -99999
      end do

   end subroutine allocateBoundaryArrays
   !------------------------------------------------------------------

   ! ------------------------------------------------------------------
   !                 S U B R O U T I N E
   ! A L L O C A T E  E L E V A T I O N  B O U N D A R Y  A R R A Y S
   ! ------------------------------------------------------------------
   ! jgf51.21.11 Allocate space for elevation boundary-related variables
   ! ------------------------------------------------------------------
   subroutine allocateElevationBoundaryArrays()
      use sizes, only: mneta, mnope
      implicit none
      allocate (nbdv(mnope, mneta))
      allocate (nvdll(mnope), nbd(mneta), ibtypee(mnope))
   end subroutine allocateElevationBoundaryArrays
   !------------------------------------------------------------------

   !------------------------------------------------------------------
   !                  S U B R O U T I N E
   !    A L L O C A T E  F L U X  B O U N D A R Y  A R R A Y S
   !------------------------------------------------------------------
   !jgf51.21.11 Allocate space for flux boundary-related variables
   !------------------------------------------------------------------
   subroutine allocateFluxBoundaryArrays()
      use sizes, only: mnvel, mnbou
      implicit none
      allocate (nbv(mnvel), lbcodei(mnvel))
      allocate (bndlen2o3(mnvel))
      allocate (me2gw(mnvel))
      allocate (csii(mnvel), siii(mnvel))

      allocate (barlanht(mnvel), barlancfsp(mnvel))
      allocate (barinht(mnvel), barincfsb(mnvel), barincfsp(mnvel))
      allocate (pipeht(mnvel), pipecoef(mnvel), pipediam(mnvel))
      allocate (ibconn(mnvel))
      allocate (issubmerged64(mnvel))

      allocate (nbvv(mnbou, 0:mnvel))
      allocate (nvell(mnbou), ibtype(mnbou), ibtype_orig(mnbou))
      allocate (ntran1(mnvel), ntran2(mnvel))
      allocate (btran3(mnvel), btran4(mnvel), btran5(mnvel))
      allocate (btran6(mnvel), btran7(mnvel), btran8(mnvel))
      allocate (nelezng(mnvel))
      allocate (zngif1(mnvel), zngif2(mnvel), zngif3(mnvel))
      ! kmd - added for rivers in baroclinic simulation
      allocate (bcrnbvv(mnbou, 0:mnvel))
      allocate (bcrnvell(mnbou))

   end subroutine allocateFluxBoundaryArrays
   !------------------------------------------------------------------

   !-----------------------------------------------------------------
   !                   S U B R O U T I N E
   !       G E T   B O U N D A R Y  S I Z E S  F O R   P R E P
   !-----------------------------------------------------------------
   ! jgf51.21.28: Called by adcprep, which uses this module to read
   ! boundary sizes in XDMF format.
   !-----------------------------------------------------------------
   subroutine getBoundarySizesForPrep(p_nbou, p_nvel, p_neta, p_nope, nweir, exist_flux)
      implicit none

      integer, intent(out) :: p_nbou, p_nope, p_nvel, p_neta, &
                              nweir
      integer, intent(out) :: exist_flux
      integer :: j, k

      p_nbou = nbou
      p_nvel = nvel
      p_neta = neta
      p_nope = nope
      !
      ! count the number of nodes on river flux boundaries and allocate
      ! memory for holding these node numbers so that adcprep can split
      ! up the fort.20 files across subdomains
      exist_flux = 0
      do k = 1, nbou
         select case (ibtype(k))
         case (2, 12, 22, 32, 52)
            exist_flux = exist_flux + size(simpleFluxBoundaries(k)%nodes)
         case default
            ! other boundary types don't count
         end select
      end do

      j = 0
      nweir = 0
      do k = 1, nbou
         select case (ibtype(k))
         case (4, 24, 64)
            nweir = nweir + nvell(k)
         case (5, 25)
            print *, 'error: ibtype 5 and 25 not supported by adcprep.'
            call exit(1)
         case default
            ! no action required for other boundary types
         end select
      end do

   end subroutine getBoundarySizesForPrep
   !-----------------------------------------------------------------

   !-----------------------------------------------------------------
   !                     S U B R O U T I N E
   !        G E T   B O U N D A R I E S   F O R   P R E P
   !-----------------------------------------------------------------
   ! jgf51.21.28: Called by adcprep, which uses this module to read
   ! mesh files in XDMF format. Because of the duplication between
   ! the pre_global module and the variables used here, this subroutine
   ! is necessary to extract parameters for adcprep. Someday adcprep
   ! will be integrated with the rest of the code and this sub will
   ! no longer be necessary.
   !-----------------------------------------------------------------
   subroutine getBoundariesForPrep(p_nvdll, p_nvell, p_ibtypee, p_nbdv, lbcode, &
                                   p_ibtype, p_nbvv, p_ibconnr, bar1, bar2, bar3, weir, &
                                   weird, exist_flux, flux14_ary)
      implicit none
      integer, intent(out) :: p_nvdll(:)
      integer, intent(out) :: p_nvell(:)
      integer, intent(out) :: p_ibtypee(:)
      integer, intent(out) :: p_nbdv(:, :), p_nbvv(:, :)
      integer, intent(out) :: p_ibtype(:)
      integer, intent(out) :: lbcode(:)
      integer, intent(out) :: p_ibconnr(:, :)
      real(8), intent(out) :: bar1(:, :)
      real(8), intent(out) :: bar2(:, :)
      real(8), intent(out) :: bar3(:, :)
      integer, intent(out) :: weir(:)
      integer, intent(out) :: weird(:)
      integer, allocatable, intent(out) :: flux14_ary(:)
      integer, intent(inout) :: exist_flux
      integer :: i, j, k, nweir

      if (nope /= 0) then
         p_nvdll = nvdll
         p_nbdv = nbdv
      end if
      if (nbou /= 0) then
         p_nvell = nvell
         P_nbvv = nbvv
      end if
      p_ibtypee = ibtypee
      p_ibtype = ibtype
      p_ibconnr = ibconnr
      lbcode = lbcodei
      allocate (flux14_ary(exist_flux))
      j = 0
      nweir = 0
      exist_flux = 0
      do k = 1, nbou
         select case (ibtype(k))
         case (0, 10, 20, 30, 40, 1, 11, 21, 41)
            do i = 1, nvell(k)
               p_ibconnr(k, i) = 0
            end do
         case (2, 12, 22, 32, 52)
            do i = 1, nvell(k)
               ibconnr(k, i) = 0
               flux14_ary(exist_flux) = nbvv(k, i)
               exist_flux = exist_flux + 1
            end do
            bndbcriver = .false.
         case (3, 13, 23)
            do i = 1, nvell(k)
               bar1(k, i) = barlanhtr(k, i)
               bar2(k, i) = barlancfspr(k, i)
               p_ibconnr(k, :) = 0
            end do
         case (4, 24)
            do i = 1, nvell(k)
               bar1(k, i) = barinhtr(k, i)
               bar2(k, i) = barincfsbr(k, i)
               bar3(k, i) = barincfspr(k, i)
               !--construct list of weir nodes and their duals
               nweir = nweir + 1
               weir(nweir) = nbvv(k, i)
               weird(nweir) = ibconnr(k, i)
            end do
         case (5, 25)
            print *, 'error: ibtype 5 and 25 not supported by adcprep.'
            call exit(1)
         case default
            write (6, '(a,i0,a)') 'ERROR: ibtype ', ibtype(k), &
               ' not supported by adcprep.'
         end select
      end do

   end subroutine getBoundariesForPrep
   !-----------------------------------------------------------------

   !------------------------------------------------------------------
   !                  S U B R O U T I N E
   !    A L L O C A T E  F L U X  B O U N D A R Y  A R R A Y
   !                 T E M P O R A R I E S
   !------------------------------------------------------------------
   !jgf51.21.17 Temporarily allocate space for flux boundary parameters
   !that will not be needed after initialization.
   !------------------------------------------------------------------
   subroutine allocateFluxBoundaryArrayTemporaries()
      use sizes, only: mnvel, mnbou
      implicit none
      allocate (barlanhtr(mnbou, mnvel), barlancfspr(mnbou, mnvel))
      allocate (barinhtr(mnbou, mnvel), barincfsbr(mnbou, mnvel))
      allocate (barincfspr(mnbou, mnvel))
      allocate (pipehtr(mnbou, mnvel), pipecoefr(mnbou, mnvel))
      allocate (pipediamr(mnbou, mnvel))
      allocate (ibconnr(mnbou, mnvel))

   end subroutine allocateFluxBoundaryArrayTemporaries
   !------------------------------------------------------------------

   !------------------------------------------------------------------
   !                  S U B R O U T I N E
   !        F R E E  F L U X  B O U N D A R Y  A R R A Y
   !                 T E M P O R A R I E S
   !------------------------------------------------------------------
   !jgf51.21.17 Free space required for flux boundary parameters after
   !boundaries are initialized.
   !------------------------------------------------------------------
   subroutine freeFluxBoundaryArrayTemporaries()
      implicit none
      deallocate (barlanhtr, barlancfspr)
      deallocate (barinhtr, barincfsbr)
      deallocate (barincfspr)
      deallocate (pipehtr, pipecoefr)
      deallocate (pipediamr)
      deallocate (ibconnr)

   end subroutine freeFluxBoundaryArrayTemporaries
   !------------------------------------------------------------------

end module boundaries
!------------------------------------------------------------------
