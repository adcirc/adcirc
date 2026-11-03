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
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
#include "logging_macros.h"
module Couple2Swan
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

   use WRITE_OUTPUT, only: terminate
   use mod_logging, only: allMessage, DEBUG, ECHO, INFO, &
                          WARNING, ERROR, t_log_scope, init_log_scope
#ifdef CSWAN
   use SWCOMM1, only: NMOVAR
#endif

   implicit none

   private

#ifdef CSWAN
!Cobell 20261103: Entities used by ADCIRC (adcirc.F, driver.F,
!                 read_input.F, timestep.F).
   public :: PADCSWAN_INIT, PADCSWAN_RUN, PADCSWAN_FINAL, &
             ComputeWaveDrivenForces, CouplingInterval, InterpoWeight, &
             SWAN_WX2, SWAN_WY2, WriteSwanHotStart, SwanHotStartUnit

!Cobell 20261103: Entities used by SWAN through the !ADC lines in
!                 thirdparty/swan.
   public :: ComputeRadiationStresses, ComputeSwanWindDrag, SwanOutput, &
             SWAN_ETA2, SWAN_UU2, SWAN_VV2, SWAN_Z0, &
             COUPCUR, COUPWIND, COUPWLV, COUPFRIC

!Cobell 20261103: Interfaces for procedures provided by SWAN (Fortran 77
!                 code, no module interfaces). These mirror the declarations
!                 in the SWAN sources after switch.pl -adcirc processing and
!                 must be kept in sync with them.
   interface
!... swanser.ftn
      subroutine KSCIP1(MMT, SIG, D, K, CG, N, ND)
         implicit none
         integer, intent(in) :: MMT
         real(4), intent(in) :: SIG(MMT), D
         real(4), intent(out) :: K(MMT), CG(MMT), N(MMT), ND(MMT)
      end subroutine KSCIP1

!... ocpmix.ftn
      subroutine DTRETI(TSTRNG, IOPT, TIMESC)
         implicit none
         character(LEN=*), intent(in) :: TSTRNG
         integer, intent(in) :: IOPT
         real(8), intent(out) :: TIMESC
      end subroutine DTRETI

!... swanmain.ftn
      subroutine SWMAIN(ITIME, IT)
         implicit none
         integer, intent(in) :: ITIME, IT
      end subroutine SWMAIN

!... swanparll.ftn
      subroutine SWINITMPI()
         implicit none
      end subroutine SWINITMPI
      subroutine SWEXITMPI()
         implicit none
      end subroutine SWEXITMPI

!... swanout1.ftn
      subroutine SWOEXD(RTYPE, OQPROC, MIP, XC, YC, VOQR, VOQ, COMPDA, &
                        KGRPNT, FORCE, CROSS, IONOD, IRQ)
         use SWCOMM1, only: NMOVAR
         use SWCOMM3, only: MCGRD, MCMVAR, MXC, MYC
         use SwanGriddata, only: nverts
         implicit none
         character(LEN=4), intent(in) :: RTYPE
         logical, intent(in) :: OQPROC(NMOVAR)
         integer, intent(in) :: MIP
         real(4), intent(in) :: XC(MIP), YC(MIP)
         integer, intent(in) :: VOQR(NMOVAR)
!... SWAN declares VOQ(MIP,*); its columns are addressed through VOQR.
         ! allow(C071)
         real(4), intent(inout) :: VOQ(MIP, *)
         real(4), intent(in) :: COMPDA(MCGRD, MCMVAR)
         integer, intent(in) :: KGRPNT(MXC, MYC)
         real(4), intent(in) :: FORCE(nverts, 2)
         logical, intent(in) :: CROSS(4, MIP)
         integer, intent(inout) :: IONOD(MIP)
         integer, intent(in) :: IRQ
      end subroutine SWOEXD
      subroutine SWOEXA(OQPROC, BKC, MIP, XC, YC, VOQR, VOQ, AC2, &
                        ACLOC, SPCSIG, WK, CG, SPCDIR, NE, NED, &
                        KGRPNT, DEPXY, CROSS)
         use SWCOMM1, only: NMOVAR
         use SWCOMM3, only: MCGRD, MDC, MSC, MXC, MYC
         implicit none
         logical, intent(in) :: OQPROC(NMOVAR)
         integer, intent(in) :: BKC
         integer, intent(in) :: MIP
         real(4), intent(in) :: XC(MIP), YC(MIP)
         integer, intent(in) :: VOQR(NMOVAR)
!... SWAN declares VOQ(MIP,*); its columns are addressed through VOQR.
         ! allow(C071)
         real(4), intent(inout) :: VOQ(MIP, *)
         real(4), intent(in) :: AC2(MDC, MSC, MCGRD)
         real(4), intent(inout) :: ACLOC(MDC, MSC)
         real(4), intent(in) :: SPCSIG(MSC)
         real(4), intent(inout) :: WK(MSC), CG(MSC), NE(MSC), NED(MSC)
         real(4), intent(in) :: SPCDIR(MDC, 6)
         integer, intent(in) :: KGRPNT(MXC, MYC)
         real(4), intent(in) :: DEPXY(MCGRD)
         logical, intent(in) :: CROSS(4, MIP)
      end subroutine SWOEXA
   end interface

!Casey 090302: The coupling interval controls how ADCIRC and SWAN
!             take turns during the simulation.  It is set
!             at the end of PADCSWAN_INIT.
   integer :: CouplingInterval
   integer, save :: SwanTimeStep ! Counter for SWAN time steps.

!Casey 090302: These logical variables will be reset to .TRUE.
!             if the water levels in SWAN are coupled to ADCIRC.
   logical :: COUPCUR
   logical :: COUPWIND
   logical :: COUPWLV
!Casey 091216: Added another coupling variable.
   logical :: COUPFRIC

!Casey 090302: These arrays contain the radiation stresses.
!Casey 090820: Be explicit about the size of these REAL variables.
   real(4), allocatable :: ADCIRC_SXX(:, :)
   real(4), allocatable :: ADCIRC_SXY(:, :)
   real(4), allocatable :: ADCIRC_SYY(:, :)

!Casey 090302: The interpolation weight controls which value
!             is taken when information is passed between ADCIRC
!             and SWAN.  If InterpoWeight = 0, then the value
!             is taken from the beginning of the coupling interval.
!             If InterpoWeight = 1, then the value is taken from
!             the end of the coupling interval.
   real(8) :: InterpoWeight

!Casey 090302: The following arrays will contain the ADCIRC values
!             that may be passed to SWAN.
   real(8), allocatable :: SWAN_ETA2(:, :)
   real(8), allocatable :: SWAN_UU2(:, :)
   real(8), allocatable :: SWAN_VV2(:, :)
   real(8), allocatable :: SWAN_WX2(:, :)
   real(8), allocatable :: SWAN_WY2(:, :)
!Casey 091216: Added another array.
   real(8), allocatable :: SWAN_Z0(:, :)

!Casey 090302: These variables control SWAN output.
! jgf51.21.25: Added target attribute to file position counters
! so that they could be pointed to from the OutputDataDescript_t
! for each output type and updated in subroutine writeOutArray in
! the write_output module.
   integer, allocatable, target :: IGSW(:) ! file position counters
   integer              :: NumSwanOutput
   logical              :: Processed26 = .false.
   logical              :: SWAN_OQPROC(NMOVAR)
   integer              :: SWAN_VOQR(NMOVAR)

!Casey 090303: Add variable for maximum SWAN output.
   real(8), allocatable        :: SWAN_MAX(:, :)

!Casey 100205: Add variables for writing of SWAN hot-start files.
   integer              :: SwanHotStartUnit
   logical              :: WriteSwanHotStart = .false.

!Casey 120302: Merging ice changes from Taylor Asher.
!TGA 110912:  Implementing SWAN ice handling.  This variable is set to
!             0 wherever ice reaches/exceeds a threshold value
   integer, allocatable :: ICY(:)
   real(8), allocatable :: SWAN_CICE2(:, :)

!-----------------------------------------------------------------------
contains
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E
!       C O M P U T E  R A D I A T I O N  S T R E S S E S
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   subroutine ComputeRadiationStresses(AC2, SPCDIR, SPCSIG)

!Casey 090616: On Jade, there was a problem with getting water depths
!              from the COMPDA array.  Use ETA2 and DEPTH instead.
      use GLOBAL, only: ETA2
      use M_GENARR, only: DEPTH
      use SwanGridData, only: nverts
      use SWCOMM3, only: DDIR, DEPMIN, FRINTF, GRAV, MDC, MSC

      implicit none

!Casey 090616: These REAL variables must be size(4) to match SWAN.
      real(4), intent(in) :: AC2(MDC, MSC, nverts)
      real(4), intent(in) :: SPCDIR(MDC, 6)
      real(4), intent(in) :: SPCSIG(MSC)

      integer :: I, ID, IERR, IS
      real(4) :: DEPLOC
      real(4) :: WK(MSC), CG(MSC), NE(MSC), NED(MSC)

      LOG_SCOPE_TRACED("ComputeRadiationStresses", COUPLE2SWAN_TRACING)
!... (Re)allocate the arrays.
      if (.not. allocated(ADCIRC_SXX)) then
         allocate (ADCIRC_SXX(1:nverts, 1:2), STAT=IERR)
         do I = 1, nverts
            ADCIRC_SXX(I, 1) = 0.0
            ADCIRC_SXX(I, 2) = 0.0
         end do
      end if
      if (.not. allocated(ADCIRC_SXY)) then
         allocate (ADCIRC_SXY(1:nverts, 1:2))
         do I = 1, nverts
            ADCIRC_SXY(I, 1) = 0.0
            ADCIRC_SXY(I, 2) = 0.0
         end do
      end if
      if (.not. allocated(ADCIRC_SYY)) then
         allocate (ADCIRC_SYY(1:nverts, 1:2))
         do I = 1, nverts
            ADCIRC_SYY(I, 1) = 0.0
            ADCIRC_SYY(I, 2) = 0.0
         end do
      end if
!... Loop over the nodes and compute radiation stresses.
      do I = 1, nverts
!... Transfer the stresses from the new time level to the old time level.
         ADCIRC_SXX(I, 1) = ADCIRC_SXX(I, 2)
         ADCIRC_SXY(I, 1) = ADCIRC_SXY(I, 2)
         ADCIRC_SYY(I, 1) = ADCIRC_SYY(I, 2)
!... Initialize the radiation stresses as zero at all of the nodes.
         ADCIRC_SXX(I, 2) = 0.0
         ADCIRC_SXY(I, 2) = 0.0
         ADCIRC_SYY(I, 2) = 0.0
!... Compute ratio of group and phase velocity.
         DEPLOC = real(real(DEPTH(I), 8) + ETA2(I), 4)
         if (DEPLOC <= DEPMIN) then
            cycle
         end if
         call KSCIP1(MSC, SPCSIG, DEPLOC, WK, CG, NE, NED)
!... Loop over all sigma and theta.
         do ID = 1, MDC
            do IS = 1, MSC
!... Sum contributions to radiation stresses.
               ADCIRC_SXX(I, 2) = ADCIRC_SXX(I, 2) &
                                  + (NE(IS)*SPCDIR(ID, 4) + NE(IS) - 0.5) &
                                  *SPCSIG(IS)*SPCSIG(IS)*AC2(ID, IS, I)
               ADCIRC_SXY(I, 2) = ADCIRC_SXY(I, 2) &
                                  + NE(IS)*SPCDIR(ID, 5) &
                                  *SPCSIG(IS)*SPCSIG(IS)*AC2(ID, IS, I)
               ADCIRC_SYY(I, 2) = ADCIRC_SYY(I, 2) &
                                  + (NE(IS)*SPCDIR(ID, 6) + NE(IS) - 0.5) &
                                  *SPCSIG(IS)*SPCSIG(IS)*AC2(ID, IS, I)
            end do
         end do
!... Multiply summed radiation stresses by the stuff outside of the sums.
!Casey 080602: ADCIRC accepts wave-driven stresses "in units of velocity squared
!    (consistent with the units of gravity).  Stress in these units is obtained
!    by dividing stress in units of force/area by the reference density of water."
!    SO WE MUST DIVIDE BY RHO! (WJP: just ommitted RHO multiple in this
!    line to avoid the division by RHO below)
         ADCIRC_SXX(I, 2) = GRAV*ADCIRC_SXX(I, 2)*DDIR*FRINTF
         ADCIRC_SXY(I, 2) = GRAV*ADCIRC_SXY(I, 2)*DDIR*FRINTF
         ADCIRC_SYY(I, 2) = GRAV*ADCIRC_SYY(I, 2)*DDIR*FRINTF
!... End loop over wet nodes.  The radiation stresses are ready for ADCIRC.
      end do

      return
!-----------------------------------------------------------------------
   end subroutine ComputeRadiationStresses
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E
!       C O M P U T E  W A V E  D R I V E N  F O R C E S
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   subroutine ComputeWaveDrivenForces

      use GLOBAL, only: NODECODE, NOFF, RSNX2, RSNY2, IFSFM, &
                        SWAN_RSNX1, SWAN_RSNY1, SWAN_RSNX2, SWAN_RSNY2
      use MESH, only: NE, NM, NP, AREAS, NODELE, NEITABELE, &
                      FDXE, FDYE, SFacEle, SFMYEle, SFMXEle
      use BOUNDARIES, only: NBDV, NBOU, NBVV, NOPE, NVDLL, NVELL
      implicit none

      integer :: I
      integer :: IE
      integer :: IP
      integer :: K
      integer :: Node1
      integer :: Node2
      integer :: Node3
      integer :: NUMFOUND

!...  To zero out forces on boundary/dry nodes
#ifdef MARCELSWAN
      logical, parameter :: Marcel = .true.
#else
      logical, parameter :: Marcel = .false.
#endif
      real(8), allocatable :: DSXXDX(:)
      real(8), allocatable :: DSXYDY(:)
      real(8), allocatable :: DSXYDX(:)
      real(8), allocatable :: DSYYDY(:)

      integer             :: NCELE

      real(8), allocatable :: TEMP_SXX(:)
      real(8), allocatable :: TEMP_SXY(:)
      real(8), allocatable :: TEMP_SYY(:)

      real(8) :: TOTALAREA
      real(8) :: FDX1, FDX2, FDX3
      real(8) :: FDY1, FDY2, FDY3
      real(8) :: SFacAvg, SFmxAvg, SFmyAvg, sfdxfac, sfdyfac

!     Upper limit on Radiation stresses
!     Default - don't limit (set to high value) but may want to use: UL = 0.1D0
!     (0.1 m^2/s^2 is equivalent to 175 m/s wind speeds!)
      real(8), parameter :: UL = 999d0

      LOG_SCOPE_TRACED("ComputeWaveDrivenForces", COUPLE2SWAN_TRACING)
!... Check whether radiation stresses have already been computed.
!... If not, then apply forces of zero.
      if (.false.) then

         if (.not. allocated(RSNX2)) then
            allocate (RSNX2(1:NP))
            RSNX2(:) = 0.d0
         end if
         if (.not. allocated(RSNY2)) then
            allocate (RSNY2(1:NP))
            RSNY2(:) = 0.d0
         end if

!... If so, then continue to compute wave-driven forces.
      else

!... Allocate arrays for radiation stresses.
         if (.not. allocated(TEMP_SXX)) allocate (TEMP_SXX(1:NP))
         if (.not. allocated(TEMP_SXY)) allocate (TEMP_SXY(1:NP))
         if (.not. allocated(TEMP_SYY)) allocate (TEMP_SYY(1:NP))

!... Loop over all nodes and interpolate the radiation stress for this time step.
         do IP = 1, NP
            TEMP_SXX(IP) = (1.d0 - InterpoWeight)*dble(ADCIRC_SXX(IP, 1)) &
                           + InterpoWeight*dble(ADCIRC_SXX(IP, 2))
            TEMP_SXY(IP) = (1.d0 - InterpoWeight)*dble(ADCIRC_SXY(IP, 1)) &
                           + InterpoWeight*dble(ADCIRC_SXY(IP, 2))
            TEMP_SYY(IP) = (1.d0 - InterpoWeight)*dble(ADCIRC_SYY(IP, 1)) &
                           + InterpoWeight*dble(ADCIRC_SYY(IP, 2))
         end do

!... Allocate arrays for radiation stress gradients.
         if (.not. allocated(DSXXDX)) allocate (DSXXDX(1:NE))
         if (.not. allocated(DSXYDY)) allocate (DSXYDY(1:NE))
         if (.not. allocated(DSXYDX)) allocate (DSXYDX(1:NE))
         if (.not. allocated(DSYYDY)) allocate (DSYYDY(1:NE))

!... Loop over all elements and compute the derivatives of Sxx, Sxy and Syy.
!... These derivatives are constant on an element.  Note that the AREAS array
!... actually contains twice the area of each element.
         do IE = 1, NE

!Casey 090707: When using the serial adcswan on Zas, I received memory errors
!... when these calls were nested into the logic below.  Break them out and
!... use these variables to save on the number of calls to memory.
            Node1 = NM(IE, 1)
            Node2 = NM(IE, 2)
            Node3 = NM(IE, 3)
!     WJP: tke into account correction factors
            SFacAvg = SFacEle(IE)
!..... BEG DW/WJP
            SFmxAvg = SFMXEle(IE); 
            SFmyAvg = SFMYEle(IE); 
            sfdxfac = dble(1 - IFSFM)*SFacAvg + dble(IFSFM)*SFmxAvg; 
            sfdyfac = dble(1 - IFSFM)*1.0d0 + dble(IFSFM)*SFmyAvg; 
            FDX1 = FDXE(1, IE)*sfdxfac; !c FDX1=(Y(NM2)-Y(NM3))*MX !b1*mx
            FDX2 = FDXE(2, IE)*sfdxfac; !c FDX2=(Y(NM3)-Y(NM1))*MX !b2*mx
            FDX3 = FDXE(3, IE)*sfdxfac; !c FDX3=(Y(NM1)-Y(NM2))*MX !b3*mx
            FDY1 = FDYE(1, IE)*sfdyfac; !c FDY1=(X(NM3)-X(NM2))*MY !a1*my
            FDY2 = FDYE(2, IE)*sfdyfac; !c FDY2=(X(NM1)-X(NM3))*MY !a2*my
            FDY3 = FDYE(3, IE)*sfdyfac; !c FDY3=(X(NM2)-X(NM1))*MY !a3*my
!..... END DW/WJP

            DSXXDX(IE) = (1.d0/AREAS(IE))* &
                         (TEMP_SXX(Node1)*FDX1 &
                          + TEMP_SXX(Node2)*FDX2 &
                          + TEMP_SXX(Node3)*FDX3)

            DSXYDY(IE) = (1.d0/AREAS(IE))* &
                         (TEMP_SXY(Node1)*FDY1 &
                          + TEMP_SXY(Node2)*FDY2 &
                          + TEMP_SXY(Node3)*FDY3)

            DSXYDX(IE) = (1.d0/AREAS(IE))* &
                         (TEMP_SXY(Node1)*FDX1 &
                          + TEMP_SXY(Node2)*FDX2 &
                          + TEMP_SXY(Node3)*FDX3)

            DSYYDY(IE) = (1.d0/AREAS(IE))* &
                         (TEMP_SYY(Node1)*FDY1 &
                          + TEMP_SYY(Node2)*FDY2 &
                          + TEMP_SYY(Node3)*FDY3)

         end do

!... Allocate arrays for wave-driven forces.
         if (.not. allocated(RSNX2)) allocate (RSNX2(1:NP))
         if (.not. allocated(RSNY2)) allocate (RSNY2(1:NP))

!... Loop over all nodes and compute the wave-driven forces:
!...
!...       Fx = - DSxx/Dx - DSxy/Dy
!...
!...       Fy = - DSxy/Dx - DSyy/Dy
!...
!... We project the element-based radiation stress gradients onto the nodes
!... by taking a weighted average of the gradients in the elements connected
!... to a node.
         outer: do IP = 1, NP

            SWAN_RSNX1(IP) = SWAN_RSNX2(IP)
            SWAN_RSNY1(IP) = SWAN_RSNY2(IP)
            SWAN_RSNX2(IP) = 0.d0
            SWAN_RSNY2(IP) = 0.d0

            TOTALAREA = 0.d0

            IE = 0
            NUMFOUND = 0

            inner: do

               IE = IE + 1

               if (IE > size(NEITABELE, 2)) exit inner

               if (NEITABELE(IP, IE) == 0) then

                  continue

               else

!... Try Marcel's method of zero-ing out the forces at nodes connected to dry nodes/elements.

                  NCELE = NODECODE(NM(NEITABELE(IP, IE), 1)) &
                          *NODECODE(NM(NEITABELE(IP, IE), 2)) &
                          *NODECODE(NM(NEITABELE(IP, IE), 3)) &
                          *NOFF(NEITABELE(IP, IE))

                  if (Marcel .and. (NCELE == 0)) then

                     RSNX2(IP) = 0.0d0
                     RSNY2(IP) = 0.0d0
                     cycle outer

                  else

                     NUMFOUND = NUMFOUND + 1

                     SWAN_RSNX2(IP) = SWAN_RSNX2(IP) + 0.5d0*AREAS(NEITABELE(IP, IE)) &
                                      *(-DSXXDX(NEITABELE(IP, IE)) &
                                        - DSXYDY(NEITABELE(IP, IE)))
                     SWAN_RSNY2(IP) = SWAN_RSNY2(IP) + 0.5d0*AREAS(NEITABELE(IP, IE)) &
                                      *(-DSXYDX(NEITABELE(IP, IE)) &
                                        - DSYYDY(NEITABELE(IP, IE)))

                     TOTALAREA = TOTALAREA + 0.5d0*AREAS(NEITABELE(IP, IE))

                  end if

               end if

               if (NUMFOUND == NODELE(IP)) then

                  exit inner

               end if

            end do inner

            SWAN_RSNX2(IP) = SWAN_RSNX2(IP)/TOTALAREA
            SWAN_RSNY2(IP) = SWAN_RSNY2(IP)/TOTALAREA

!    WJPs artifical limiter; UL is to a really large by default, user
!    can change UL to be a strict value if they please
            SWAN_RSNX2(IP) = min(UL, max(-UL, SWAN_RSNX2(IP)))
            SWAN_RSNY2(IP) = min(UL, max(-UL, SWAN_RSNY2(IP)))

         end do outer

!... Try Marcel's method of zero-ing the forces at the boundary nodes.
         if (Marcel) then
            do K = 1, NOPE
               do I = 1, NVDLL(K)
                  SWAN_RSNX2(NBDV(K, I)) = 0.0d0
                  SWAN_RSNY2(NBDV(K, I)) = 0.0d0
               end do
            end do
            do K = 1, NBOU
               do I = 1, NVELL(K)
                  SWAN_RSNX2(NBVV(K, I)) = 0.0d0
                  SWAN_RSNY2(NBVV(K, I)) = 0.0d0
               end do
            end do
         end if

!... Deallocate the radiation stress gradients.
         if (allocated(DSXXDX)) deallocate (DSXXDX)
         if (allocated(DSXYDY)) deallocate (DSXYDY)
         if (allocated(DSXYDX)) deallocate (DSXYDX)
         if (allocated(DSYYDY)) deallocate (DSYYDY)
         if (allocated(TEMP_SXX)) deallocate (TEMP_SXX)
         if (allocated(TEMP_SXY)) deallocate (TEMP_SXY)
         if (allocated(TEMP_SYY)) deallocate (TEMP_SYY)

      end if

      return
!-----------------------------------------------------------------------
   end subroutine ComputeWaveDrivenForces
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E
!       C O M P U T E  S W A N  W I N D  D R A G
!-----------------------------------------------------------------------
!     @mattbilskie
!-----------------------------------------------------------------------
   subroutine ComputeSwanWindDrag(CDRAG, NodeNumber)
      use GLOBAL, only: NP_G, SWAN_WDragCo
      implicit none
      real(4), intent(in) :: CDRAG
      integer, intent(in) :: NodeNumber
      if (.not. allocated(SWAN_WDragCo)) allocate (SWAN_WDragCo(NP_G))
      SWAN_WDragCo(NodeNumber) = real(CDRAG, 8)
!-----------------------------------------------------------------------
   end subroutine ComputeSwanWindDrag
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E  M A N N I N G  2  M A D S E N
!-----------------------------------------------------------------------
!Casey 091216: This routine will convert the ADCIRC Manning's n values
!             into roughness lengths that can be used with the Madsen
!             friction formulation inside SWAN.
!-----------------------------------------------------------------------
   subroutine Manning2Madsen
      use ADC_CONSTANTS, only: G
      use MESH, only: NP, DP
      use NodalAttributes, only: ManningsN

      implicit none

      integer  :: IN

      real(8) :: H
      real(8), parameter :: K = 0.4d0
      real(8) :: N
      real(8) :: Z0

      LOG_SCOPE_TRACED("Manning2Madsen", COUPLE2SWAN_TRACING)
      do IN = 1, NP

         H = SWAN_ETA2(IN, 2) + DP(IN)
         N = ManningsN(IN)

!Casey 110518: Enforce a lower limit on the Manning's n seen by SWAN.
         if (N < 0.02d0) N = 0.02d0

         Z0 = (H)*exp(-1.d0*(1.d0 + K*H**(1.d0/6.d0) &
                             /(N*sqrt(G))))*30.0d0

!Casey 091216: If we get a junk number, then use the default value.
         if (Z0 <= 0.d0) then
            Z0 = 0.05d0
         end if

         SWAN_Z0(IN, 2) = Z0

      end do

!-----------------------------------------------------------------------
   end subroutine Manning2Madsen
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E    S W A N  O U T P U T
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   subroutine SwanOutput

      use Couple2Adcirc, only: COMPDA
      use MESH, only: DP, NP
      use GLOBAL, only: NP_G, ETA2, NODECODE, NOUTGW, UU2, VV2, h0, &
                        SWAN_OutputHS, &
                        SWAN_OutputDIR, &
                        SWAN_OutputTM01, &
                        SWAN_OutputTPS, &
                        SWAN_OutputWIND, &
                        SWAN_OutputTM02, &
                        SWAN_OutputTMM10, &
                        SWAN_OutputAgg, &
                        Swan_HSOut, Swan_TPSOut, Swan_TM01Out, Swan_DirOut, &
                        Swan_TMM10Out, Swan_WindXOut, Swan_WindYOut, &
                        Swan_HSMaxOut, Swan_TPSMaxOut, Swan_TM01MaxOut, &
                        Swan_TMM10MaxOut, Swan_WindMaxOut, Swan_DirMaxOut, &
                        Swan_TM02Out, Swan_TM02MaxOut
      use GLOBAL_IO, only: Header73, Header74, OPEN_GBL_FILE
      use M_GENARR, only: AC2, &
                          KGRPNT, &
                          SPCDIR, &
                          SPCSIG
      use SIZES, only: GLOBALDIR
      use SwanGriddata, only: nverts, &
                              xcugrd, &
                              ycugrd
      use SWCOMM1, only: COSCQ, &
                         OVEXCV, &
                         SINCQ
      use SWCOMM3, only: MDC, &
                         MSC

      implicit none

      intrinsic                     :: ALLOCATED
      intrinsic                     :: TRIM

      character(LEN=15), allocatable :: FileName(:)
      character(LEN=10), allocatable :: Names(:)
!Cobell 20261103: SWOEXD only acts on RTYPE == 'BLKV' (VTK block output),
!                 which does not apply here, so pass blanks.
      character(LEN=4), parameter   :: SWAN_RTYPE = '    '
!... SWAN output variable index (OVEXCV/SWAN_VOQR) for each ADCIRC output slot IW:
!    HS, DIR, TM01, TPS, WIND, TM02, TMM10
      integer, parameter            :: SWAN_IVTYPE(7) = [10, 13, 11, 53, 26, 32, 47]

      integer                       :: IO
      integer                       :: IONOD(NP)
      integer                       :: IP
      integer                       :: IS
      integer                       :: IVTYPE
      integer                       :: IW
      integer                       :: SWAN_BKC
      integer                       :: UnitNumber
      integer                       :: UpdateMax(NP)

      logical                       :: CROSS(4, NP)

!Casey 090820: Sapphire doesn't like it if these variables
!              are declared as only REAL.  They must be declared
!              as REAL(4) to interface correctly with the
!              SWAN output subroutines.
      real(4)                       :: ACLOC(MDC, MSC)
      real(4)                       :: DEPXY(NP)
      real(4)                       :: FORCE(NP, 2)
      real(4)                       :: SWAN_CG(MSC)
      real(4)                       :: SWAN_NE(MSC)
      real(4)                       :: SWAN_NED(MSC)
      real(4), allocatable           :: SWAN_VOQ(:, :)
      real(4)                       :: SWAN_WK(MSC)
      real(8), pointer              :: SwanOut(:)
      real(8), pointer              :: SwanOut2(:)
      real(8), pointer              :: SwanMaxOut(:)
      real(4)                       :: XC(NP)
      real(4)                       :: YC(NP)

      LOG_SCOPE_TRACED("SwanOutput", COUPLE2SWAN_TRACING)
!... Initialize the Swan output variables.
      do IO = 1, NMOVAR
         SWAN_OQPROC(IO) = .false.
      end do
      NumSwanOutput = 0
!... For significant wave heights (HS).
      NumSwanOutput = NumSwanOutput + 1
      OVEXCV(10) = -99999.
      SWAN_OQPROC(10) = .true.
      SWAN_VOQR(10) = 7 + NumSwanOutput
!... For mean wave directions (DIR).
      NumSwanOutput = NumSwanOutput + 1
      OVEXCV(13) = -99999.
      SWAN_OQPROC(13) = .true.
      SWAN_VOQR(13) = 7 + NumSwanOutput
!... For mean wave periods (TM01).
      NumSwanOutput = NumSwanOutput + 1
      OVEXCV(11) = -99999.
      SWAN_OQPROC(11) = .true.
      SWAN_VOQR(11) = 7 + NumSwanOutput
!... For peak wave periods (TPS).
      NumSwanOutput = NumSwanOutput + 1
      OVEXCV(53) = -99999.
      SWAN_OQPROC(53) = .true.
      SWAN_VOQR(53) = 7 + NumSwanOutput
!... For wind speeds (WX2 and WY2).
      NumSwanOutput = NumSwanOutput + 1
      OVEXCV(26) = 0.
      SWAN_OQPROC(26) = .true.
      SWAN_VOQR(26) = 7 + NumSwanOutput
!... For mean wave periods (TM02).
      NumSwanOutput = NumSwanOutput + 1
      OVEXCV(32) = -99999.
      SWAN_OQPROC(32) = .true.
      SWAN_VOQR(32) = 7 + NumSwanOutput + 1
!... For mean wave periods (TMM10).
      NumSwanOutput = NumSwanOutput + 1
      OVEXCV(47) = -99999.
      SWAN_OQPROC(47) = .true.
      SWAN_VOQR(47) = 7 + NumSwanOutput + 1
      COSCQ = cos(0.)
      SINCQ = sin(0.)
      if (.not. allocated(IGSW)) allocate (IGSW(1:NumSwanOutput))
      if (.not. allocated(Names)) allocate (Names(1:NumSwanOutput))
      Names(1) = "HS"
      Names(2) = "DIR"
      Names(3) = "TM01"
      Names(4) = "TPS"
      Names(5) = "WIND"
      Names(6) = "TM02"
      Names(7) = "TMM10"
      if (.not. allocated(FileName)) allocate (FileName(1:NumSwanOutput))
      FileName(1) = "swan_"//trim(Names(1))//".63"//"       "
      FileName(2) = "swan_"//trim(Names(2))//".63"//"       "
      FileName(3) = "swan_"//trim(Names(3))//".63"//"       "
      FileName(4) = "swan_"//trim(Names(4))//".63"//"       "
      FileName(5) = "swan_"//trim(Names(5))//".64"//"       "
      FileName(6) = "swan_"//trim(Names(6))//".63"//"       "
      FileName(7) = "swan_"//trim(Names(7))//".63"//"       "

      if (.not. Processed26) then
!Cobell 20120510: At first pass, we set up the output data arrays
!.....Arrange on/off for output into array
         SWAN_OutputAgg(1) = SWAN_OutputHS
         SWAN_OutputAgg(2) = SWAN_OutputDIR
         SWAN_OutputAgg(3) = SWAN_OutputTM01
         SWAN_OutputAgg(4) = SWAN_OutputTPS
         SWAN_OutputAgg(5) = SWAN_OutputWIND
         SWAN_OutputAgg(6) = SWAN_OutputTM02
         SWAN_OutputAgg(7) = SWAN_OutputTMM10

         if (.not. allocated(SWAN_MAX)) allocate (SWAN_MAX(1:NP, 1:NumSwanOutput))
         do IP = 1, nverts
            do IW = 1, NumSwanOutput
               IVTYPE = SWAN_IVTYPE(IW)
               SWAN_MAX(IP, IW) = real(OVEXCV(IVTYPE), 8)
            end do
         end do

!        !...Open the global output files
         do IW = 1, NumSwanOutput
            if (.not. SWAN_OutputAgg(IW)) cycle
            UnitNumber = 300 + IW

            if ((abs(NOUTGW) == 1) .or. (abs(NOUTGW) == 4)) then
               if (IW /= 5) then
                  call OPEN_GBL_FILE(UnitNumber, trim(GLOBALDIR)//'/'//trim(FileName(IW)), &
                                     NP_G, NP, Header73)
               else
                  call OPEN_GBL_FILE(UnitNumber, trim(GLOBALDIR)//'/'//trim(FileName(IW)), &
                                     NP_G, NP, Header74)
               end if
            end if
         end do

         Processed26 = .true.
      end if

!... I don't know what BKC does.  Everything seems to work okay
!... when it is set to two, though.
      SWAN_BKC = 2
!... If Swan needs these values for wave numbers, etc.,
!... then maybe it solves for them inside SWOEXA.
      do IS = 1, MSC
         SWAN_CG(IS) = 0.
         SWAN_NE(IS) = 0.
         SWAN_NED(IS) = 0.
         SWAN_WK(IS) = 0.
      end do
!... Allocate VOQ and set up the first seven entries inside it.
      if (.not. allocated(SWAN_VOQ)) allocate (SWAN_VOQ(1:NP, 1:(7 + NumSwanOutput + 1)))
      do IP = 1, nverts
!... X,Y in problem coordinate system.
         SWAN_OQPROC(1) = .true.
         SWAN_OQPROC(2) = .true.
         SWAN_VOQR(1) = 1
         SWAN_VOQR(2) = 2
         SWAN_VOQ(IP, SWAN_VOQR(1)) = xcugrd(IP)
         SWAN_VOQ(IP, SWAN_VOQR(2)) = ycugrd(IP)
!... X,Y in adjusted coordinate system.
         SWAN_OQPROC(24) = .true.
         SWAN_OQPROC(25) = .true.
         SWAN_VOQR(24) = 3
         SWAN_VOQR(25) = 4
         SWAN_VOQ(IP, SWAN_VOQR(24)) = 0.
         SWAN_VOQ(IP, SWAN_VOQR(25)) = 0.
!... Depth of water.
         SWAN_OQPROC(4) = .true.
         SWAN_VOQR(4) = 5
         SWAN_VOQ(IP, SWAN_VOQR(4)) = real(DP(IP)) + real(ETA2(IP))

! jgf51.46: Fix from Casey to eliminate floating point
! exception in SWAN (swanser.F lines 703 and 726) when depths
! are negative.
         if (SWAN_VOQ(IP, SWAN_VOQR(4)) < real(H0)) then
            SWAN_VOQ(IP, SWAN_VOQR(4)) = real(H0)
         end if

!... Water current velocities.
         SWAN_OQPROC(5) = .true.
         SWAN_VOQR(5) = 6
         SWAN_VOQ(IP, SWAN_VOQR(5)) = real(UU2(IP))
         SWAN_VOQ(IP, SWAN_VOQR(5) + 1) = real(VV2(IP))
!... Assemble other arrays?
         DEPXY(IP) = real(DP(IP)) + real(ETA2(IP))
         XC(IP) = xcugrd(IP)
         YC(IP) = ycugrd(IP)
         FORCE(IP, 1) = 0.
         FORCE(IP, 2) = 0.
         IONOD(IP) = -999
      end do

!... Call the Swan subroutine to interpolate internal Swan quantities.
      call SWOEXD(SWAN_RTYPE, SWAN_OQPROC, nverts, XC, YC, SWAN_VOQR, &
                  SWAN_VOQ, COMPDA, KGRPNT, FORCE, CROSS, IONOD, -999)

!... Call the Swan subroutine to compute the output quantities.
      call SWOEXA(SWAN_OQPROC, SWAN_BKC, nverts, XC, YC, &
                  SWAN_VOQR, SWAN_VOQ, AC2, ACLOC, SPCSIG, &
                  SWAN_WK, SWAN_CG, SPCDIR, SWAN_NE, SWAN_NED, &
                  KGRPNT, DEPXY, CROSS)

      UpdateMax = 0
      do IW = 1, NumSwanOutput
         if (.not. SWAN_OutputAgg(IW)) cycle
         UnitNumber = 300 + IW
!... Assign the output variable to our data structure.
         IVTYPE = SWAN_IVTYPE(IW)
         select case (IW)
         case (1)
            SwanOut => Swan_HSOut
            SwanMaxOut => Swan_HSMaxOut
         case (2)
            SwanOut => Swan_DIROut
            SwanMaxOut => Swan_DIRMaxOut
         case (3)
            SwanOut => Swan_TM01Out
            SwanMaxOut => Swan_TM01MaxOut
         case (4)
            SwanOut => Swan_TPSOut
            SwanMaxOut => Swan_TPSMaxOut
         case (5)
            SwanOut => Swan_WindXOut
            SwanOut2 => Swan_WindYOut
            SwanMaxOut => Swan_WindMaxOut
         case (6)
            SwanOut => Swan_TM02Out
            SwanMaxOut => Swan_TM02MaxOut
         case (7)
            SwanOut => Swan_TMM10Out
            SwanMaxOut => Swan_TMM10MaxOut
         case DEFAULT
            cycle
         end select
         do IP = 1, NP
            SwanOut(IP) = dble(SWAN_VOQ(IP, SWAN_VOQR(IVTYPE)))
            if (IW == 5) then
               SwanOut2(IP) = dble(SWAN_VOQ(IP, SWAN_VOQR(IVTYPE) + 1))
            end if
!           ! jgf51.21.27: Move this here so it can be used to write
!           ! data in various formats below.
!           ! SwanDescript % alternate_value =  DBLE(OVEXCV(IVTYPE))

            if (IW /= 5) then
!Casey 110518: Ensure that dry nodes are written with default values.
!Casey 120522: Correct logic for the files that don't contain significant wave heights.
               if (NODECODE(IP) == 1) then
                  if (IW == 1) then
                     if (SwanOut(IP) > SwanMaxOut(IP)) then
                        UpdateMax(IP) = 1
                        SwanMaxOut(IP) = SwanOut(IP)
                     end if
                  elseif (UpdateMax(IP) == 1) then
                     SwanMaxOut(IP) = SwanOut(IP)
                  end if
               end if
            else
               if ((NODECODE(IP) == 1) .and. &
                   (sqrt(SwanOut(IP)*SwanOut(IP) + SwanOut2(IP)*SwanOut2(IP)) > SwanMaxOut(IP))) then
                  if (UpdateMax(IP) == 1) then
                     SwanMaxOut(IP) = sqrt(SwanOut(IP)*SwanOut(IP) + SwanOut2(IP)*SwanOut2(IP))
                  end if
               end if
            end if
         end do
      end do

!Casey 090620: More deallocations.
      if (allocated(FileName)) deallocate (FileName)
      if (allocated(IGSW)) deallocate (IGSW)
      if (allocated(Names)) deallocate (Names)
      if (allocated(SWAN_VOQ)) deallocate (SWAN_VOQ)

      return
!-----------------------------------------------------------------------
   end subroutine SwanOutput
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E    P A D C S W A N _ I N I T
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   subroutine PADCSWAN_INIT

!     !USE, INTRINSIC :: IEEE_ARITHMETIC !jgfdebug ieee_is_nan()
      use GLOBAL, only: DT, &
                        ETA2, UU2, VV2, &
                        NCICE, &
                        CICE1, &
                        CICE2 &
                        , RunStartDateTime
      use WIND, only: WVNX1, &
                      WVNX2, &
                      WVNY1, &
                      WVNY2
      use MESH, only: NP, DP
!TGA 111003:  Adding ability to read in ice info from fort.26 file
      use SWCOMM3, only: SWAN_WBICETH => WBICETH, &
                         SWAN_IICE => IICE
      use TIMECOMM, only: SWAN_DT => DT &
                                                         , SWAN_TINIC => TINIC

      implicit none

      intrinsic :: ALLOCATED

      integer   :: IN ! Node counter.
!Casey 251212: Adding time control.
      real(8) :: AdcircStartTime

      LOG_SCOPE_TRACED("PADCSWAN_INIT", COUPLE2SWAN_TRACING)
!Casey 090302: Allocate memory for the ADCIRC values that may be
!             passed to SWAN.
      if (.not. allocated(SWAN_ETA2)) allocate (SWAN_ETA2(1:NP, 1:2))
      if (.not. allocated(SWAN_UU2)) allocate (SWAN_UU2(1:NP, 1:2))
      if (.not. allocated(SWAN_VV2)) allocate (SWAN_VV2(1:NP, 1:2))
      if (.not. allocated(SWAN_WX2)) allocate (SWAN_WX2(1:NP, 1:2))
      if (.not. allocated(SWAN_WY2)) allocate (SWAN_WY2(1:NP, 1:2))
!Casey 091216: Added allocation for friction coupling.
      if (.not. allocated(SWAN_Z0)) allocate (SWAN_Z0(1:NP, 1:2))

!Casey 120302: Merging ice changes from Taylor Asher.
!TGA 110912:  Full SWAN ice implementation.  Allocating ice arrays.
      if (NCICE /= 0) then
         if (.not. allocated(ICY)) allocate (ICY(1:NP))
         if (.not. allocated(SWAN_CICE2)) allocate (SWAN_CICE2(1:NP, 1:2))
         ICY = 0
         SWAN_CICE2 = 0.d0
      end if

!Casey 090302: Initialize the arrays.
      do IN = 1, NP
         SWAN_ETA2(IN, 1) = ETA2(IN)
         SWAN_ETA2(IN, 2) = ETA2(IN)
         SWAN_UU2(IN, 1) = UU2(IN) ! sb251003: Initialize velocity arrays
         SWAN_UU2(IN, 2) = UU2(IN) ! with current velocities, not zero.
         SWAN_VV2(IN, 1) = VV2(IN) ! This should improve transition at
         SWAN_VV2(IN, 2) = VV2(IN) ! hotstart.
      end do

!Casey 090302: Need to pass correct first wind snaps to SWAN.
      do IN = 1, NP
         SWAN_WX2(IN, 1) = WVNX1(IN)
         SWAN_WY2(IN, 1) = WVNY1(IN)
         SWAN_WX2(IN, 2) = WVNX2(IN)
         SWAN_WY2(IN, 2) = WVNY2(IN)
      end do

!Casey 091216: Initialize arrays for friction coupling and compute first set of values.
!sb   251003: Remove IF(COUPFRIC) so that SWAN_Z0 is treated consistently with SWAN_ETA2, SWAN_UU2, and SWAN_VV2.
      call Manning2Madsen
      do IN = 1, NP
         SWAN_Z0(IN, 1) = SWAN_Z0(IN, 2)
      end do

!Casey 120302: Merging ice changes from Taylor Asher.
!TGA 110912:  Full SWAN ice implementation.  Applying ice coverage mask
!             wherever ice concentration reaches/exceeds threshold.
!             Mask is applied by making total water depth (water
!             elevation plus bathymetry) equal to zero.  Same as what
!             is done for dry nodes.
!             The array ICY is used so it can be output in the future.
!TGA111003:  Adding ability to read in ice parameters from fort.26 file
!TGA      IF(NCICE.NE.0) THEN
      if ((NCICE /= 0) .and. (SWAN_IICE == 2)) then
         ICY = 0
         do IN = 1, NP
            SWAN_CICE2(IN, 1) = CICE1(IN)
            SWAN_CICE2(IN, 2) = CICE2(IN)
!TGA 111003:  Adding ability to read in ice parameters from fort.26 file
!TGA            IF (SWAN_CICE2(IN,1).GE.70) ICY(IN) = 1
            if (SWAN_CICE2(IN, 1) >= real(SWAN_WBICETH, 8)) ICY(IN) = 1
            if (ICY(IN) == 1) then
!TGA 110915:  Ensure SWAN considers conditions to be too icy for wave
!             generation for the current timestep by setting both
!             values of SWAN_*(IN,:) to the "dry" values.
               SWAN_ETA2(IN, :) = -DP(IN)
               SWAN_UU2(IN, :) = 0.d0
               SWAN_VV2(IN, :) = 0.d0
!TGA110915               SWAN_ETA2(IN,2) = - DP(IN)
!TGA110915               SWAN_UU2 (IN,2) = 0.D0
!TGA110915               SWAN_VV2 (IN,2) = 0.D0
            end if
         end do
      end if
! GML 2020-08-21 enable pass ACIRC ice concentration to SWAN
      if ((NCICE /= 0) .and. (SWAN_IICE == 3)) then
!         ICY=0
         do IN = 1, NP
            SWAN_CICE2(IN, 1) = CICE1(IN)
            SWAN_CICE2(IN, 2) = CICE2(IN)
!            IF (SWAN_CICE2(IN,1).GE.70) ICY(IN) = 1
!            IF(ICY(IN).EQ.1) THEN
!               SWAN_ETA2(IN,:) = - DP(IN)
!               SWAN_UU2 (IN,:) = 0.D0
!               SWAN_VV2 (IN,:) = 0.D0
!            ENDIF
         end do
      end if
!Casey 090302: Allow SWAN to initialize.
      call SWINITMPI
      call SWMAIN(0, 0)

!Casey 090302: Advance the wind speeds now so that they do not
!             get written over during the time step loop.
!sb   251003: Removing this copy of future winds to the arrays for past winds.
!             SWAN_WX/Y2(:,1) is supposed to store the past wind speed, not the future wind speed.
!             SWAN_WX/Y2(:,2) will be updated in timestep.F at each timestep. If this copy is kept,
!             Interpolation of the wind speeds conducted using InterpoWeight in SWAN would
!             be between the future (SWAN_WX/Y2(:,1)) and the current (SWAN_WX/Y2(:,2)) wind speeds,
!             which is not correct.
!      DO IN=1,NP
!         SWAN_WX2(IN,1) = SWAN_WX2(IN,2)
!         SWAN_WY2(IN,1) = SWAN_WY2(IN,2)
!      ENDDO

!Casey 090302: Compute the coupling interval as equivalent
!             to the SWAN time step.
      CouplingInterval = nint(SWAN_DT/DT)
!Casey 251212: Adding time control. Determine how many SWAN time steps
!             to skip before the start of its calculations.
      if (index(RunStartDateTime, "-") > 0) then
         SwanTimeStep = 0
      else
!        ! I will mimic code from thirdparty/swan/swanpre1.ftn
!        ! Convert datetime to real time since reference.
!Nicole 251218: fixed hardcoded date to datetime from namelist
         call DTRETI(RunStartDateTime, 1, AdcircStartTime)
!        ! Initialize the starting time step.
         SwanTimeStep = nint((AdcircStartTime - SWAN_TINIC)/SWAN_DT)
      end if

!-----------------------------------------------------------------------
   end subroutine PADCSWAN_INIT
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E    P A D C S W A N _ R U N
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   subroutine PADCSWAN_RUN(ITIME)

      use GLOBAL, only: DT, &
                        ETA2, &
                        NODECODE, &
                        UU2, &
                        VV2, &
                        STATIM, &
                        NCICE, &
                        CICE1, &
                        CICE2, &
                        CICE_TIME1, &
                        CICE_TIMINC &
                        , Swan_HSOut, &
                        Swan_DIROut, &
                        Swan_TM01Out, &
                        Swan_TPSOut, &
                        Swan_WindXOut, &
                        Swan_WindYOut, &
                        Swan_TM02Out, &
                        Swan_TMM10Out
!TGA 111003:  Adding ability to read in ice info from fort.26 file
      use MESH, only: NP, DP
      use SWCOMM3, only: SWAN_IICE => IICE, &
                         SWAN_WBICETH => WBICETH &
                                         , SWAN_MTC => MTC

      implicit none

      integer, intent(in)      :: ITIME ! Time step from ADCIRC.

      intrinsic    :: real

      integer                  :: IN ! Node counter.
!      INTEGER,SAVE :: SwanTimeStep ! Counter for SWAN time steps.

!Casey 120302: Merging ice changes from Taylor Asher.
!TGA 110912:  Full SWAN ice implementation.  Adding variables.
      real(8)             :: SWAN_PIC
      real(8)             :: ICETIMEFRAC

      LOG_SCOPE_TRACED("PADCSWAN_RUN", COUPLE2SWAN_TRACING)
      do IN = 1, NP
!Casey 090302: Move the 'old' values into the first columns.
         SWAN_ETA2(IN, 1) = SWAN_ETA2(IN, 2)
         SWAN_UU2(IN, 1) = SWAN_UU2(IN, 2)
         SWAN_VV2(IN, 1) = SWAN_VV2(IN, 2)
!Casey 090302: Find the 'new' values for the second columns.
         if (NODECODE(IN) == 1) then
            SWAN_ETA2(IN, 2) = ETA2(IN)
            SWAN_UU2(IN, 2) = UU2(IN)
            SWAN_VV2(IN, 2) = VV2(IN)
         else
            SWAN_ETA2(IN, 2) = -DP(IN)
            SWAN_UU2(IN, 2) = 0.d0
            SWAN_VV2(IN, 2) = 0.d0
         end if
      end do

!Casey 120302: Merging ice changes from Taylor Asher.
!TGA 110912:  Full SWAN ice implementation.  Applying ice coverage mask
!             wherever ice concentration reaches/exceeds threshold.
!             Mask is applied by making total water depth (water
!             elevation plus bathymetry) equal to zero.  Same as what
!             is done for dry nodes.
!             The array ICY is used so it can be output in the future.
!TGA111003:  Adding ability to read in ice parameters from fort.26 file
!TGA      IF(NCICE.NE.0) THEN
      if ((NCICE /= 0) .and. (SWAN_IICE == 2)) then
         ICY = 0
         ICETIMEFRAC = ((STATIM*86400.d0 + dble(ITIME)*DT) - CICE_TIME1)/CICE_TIMINC
         do IN = 1, NP
            SWAN_CICE2(IN, 1) = CICE1(IN)
            SWAN_CICE2(IN, 2) = CICE2(IN)
            SWAN_PIC = ICETIMEFRAC*(SWAN_CICE2(IN, 2) - SWAN_CICE2(IN, 1)) &
                       + SWAN_CICE2(IN, 1)
!TGA 111003:  Adding ability to read in ice parameters from fort.26 file
!TGA            IF (SWAN_PIC.GE.70) ICY(IN) = 1
            if (SWAN_PIC >= real(SWAN_WBICETH, 8)) ICY(IN) = 1
            if (ICY(IN) == 1) then
!TGA 110915:  Ensure SWAN considers conditions to be too icy for wave
!             generation for the current timestep by setting both
!             values of SWAN_*(IN,:) to the "dry" values.
               SWAN_ETA2(IN, :) = -DP(IN)
               SWAN_UU2(IN, :) = 0.d0
               SWAN_VV2(IN, :) = 0.d0
!TGA110915               SWAN_ETA2(IN,1) = - DP(IN)
!TGA110915               SWAN_UU2 (IN,1) = 0.D0
!TGA110915               SWAN_VV2 (IN,1) = 0.D0
!TGA110915               SWAN_ETA2(IN,2) = - DP(IN)
!TGA110915               SWAN_UU2 (IN,2) = 0.D0
!TGA110915               SWAN_VV2 (IN,2) = 0.D0
            end if
         end do
      end if
! GML 2020-08-21 enable pass ACIRC ice concentration to SWAN
      if ((NCICE /= 0) .and. (SWAN_IICE == 3)) then
!         ICY=0
         ICETIMEFRAC = ((STATIM*86400.d0 + dble(ITIME)*DT) - CICE_TIME1)/CICE_TIMINC
         do IN = 1, NP
            SWAN_CICE2(IN, 1) = CICE1(IN)
            SWAN_CICE2(IN, 2) = CICE2(IN)
            SWAN_PIC = ICETIMEFRAC*(SWAN_CICE2(IN, 2) - SWAN_CICE2(IN, 1)) &
                       + SWAN_CICE2(IN, 1)
!            IF (SWAN_PIC.GE.70) ICY(IN) = 1
!            IF(ICY(IN).EQ.1) THEN
!               SWAN_ETA2(IN,:) = - DP(IN)
!               SWAN_UU2 (IN,:) = 0.D0
!               SWAN_VV2 (IN,:) = 0.D0
!            ENDIF
         end do
      end if

!Casey 091216: Move the 'old' values into the first column,
!             and then compute the 'new' values.
!sb   251003: Remove IF(COUPFRIC) so that SWAN_Z0 is always initialized. COUPFRIC is not initialized yet here.
      do IN = 1, NP
         SWAN_Z0(IN, 1) = SWAN_Z0(IN, 2)
      end do
      call Manning2Madsen

!     !jgf50.28: Removed DO ... when SWAN coupling was under development,
!     ! this was a loop so that SWAN could do multiple time steps per coupling
!     ! interval, but it became clear that there was not a need for that and
!     ! that new information should be used at the start of each SWAN step.

!Casey 090729: Average values over the time step.
      InterpoWeight = 0.5d0
      SwanTimeStep = SwanTimeStep + 1
!Casey 251212: Adding time control.
      if ((SwanTimeStep >= 1) .and. (SwanTimeStep <= SWAN_MTC)) then
         call SWMAIN(ITIME, SwanTimeStep)
      end if

!Casey 251212: Adding time control. The call to 'SwanOutput' is made from
!             the SWAN side. If we did not call SWAN, then we need to
!             'compute' or at least initialize the output quantities.
      if ((SwanTimeStep < 1) .or. (SwanTimeStep > SWAN_MTC)) then
         Swan_HSOut(:) = -99999.d0
         Swan_DIROut(:) = -99999.d0
         Swan_TM01Out(:) = -99999.d0
         Swan_TPSOut(:) = -99999.d0
         Swan_WindXOut(:) = 0.d0
         Swan_WindYOut(:) = 0.d0
         Swan_TM02Out(:) = 0.d0
         Swan_TMM10Out(:) = -99999.d0
      end if
!Casey 251212: End changes for output.

!Casey 251212: Adding time control. The call to 'ComputeRadiationStresses'
!             is made from the SWAN side. If we did not call SWAN,
!             then we need to 'compute' the radiation stresses here.
      if ((SwanTimeStep < 1) .or. (SwanTimeStep > SWAN_MTC)) then
!... (Re)allocate the arrays.
         if (.not. allocated(ADCIRC_SXX)) then
            allocate (ADCIRC_SXX(1:NP, 1:2))
            ADCIRC_SXX(:, :) = 0.0
         end if
         if (.not. allocated(ADCIRC_SXY)) then
            allocate (ADCIRC_SXY(1:NP, 1:2))
            ADCIRC_SXY(:, :) = 0.0
         end if
         if (.not. allocated(ADCIRC_SYY)) then
            allocate (ADCIRC_SYY(1:NP, 1:2))
            ADCIRC_SYY(:, :) = 0.0
         end if
         do IN = 1, NP
!... Transfer the stresses from the new time level to the old time level.
            ADCIRC_SXX(IN, 1) = ADCIRC_SXX(IN, 2)
            ADCIRC_SXY(IN, 1) = ADCIRC_SXY(IN, 2)
            ADCIRC_SYY(IN, 1) = ADCIRC_SYY(IN, 2)
!... Initialize the radiation stresses as zero at all of the nodes.
            ADCIRC_SXX(IN, 2) = 0.0
            ADCIRC_SXY(IN, 2) = 0.0
            ADCIRC_SYY(IN, 2) = 0.0
         end do
      end if
!Casey 251212: End changes for radiation stresses.

!Casey 090302: Move the 'old' values into the first columns.
!             We do this for the winds after the SWAN time steps
!             because the 'new' values are written over in TIMESTEP.
!             We need to move them to the 'old' values now,
!             so that we do not lose them during the ADCIRC time stepping.
      do IN = 1, NP
         SWAN_WX2(IN, 1) = SWAN_WX2(IN, 2)
         SWAN_WY2(IN, 1) = SWAN_WY2(IN, 2)
      end do

!-----------------------------------------------------------------------
   end subroutine PADCSWAN_RUN
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E    P A D C S W A N _ F I N A L
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   subroutine PADCSWAN_FINAL
      implicit none

      LOG_SCOPE_TRACED("PADCSWAN_FINAL", COUPLE2SWAN_TRACING)
      call SWEXITMPI()

      if (allocated(SWAN_ETA2)) deallocate (SWAN_ETA2)
      if (allocated(SWAN_UU2)) deallocate (SWAN_UU2)
      if (allocated(SWAN_VV2)) deallocate (SWAN_VV2)
      if (allocated(SWAN_WX2)) deallocate (SWAN_WX2)
      if (allocated(SWAN_WY2)) deallocate (SWAN_WY2)

!Casey 091216: Add deallocation for friction coupling.
      if (allocated(SWAN_Z0)) deallocate (SWAN_Z0)

!Casey 120302: Merging ice changes from Taylor Asher.
!TGA 110912:
      if (allocated(ICY)) deallocate (ICY)
      if (allocated(SWAN_CICE2)) deallocate (SWAN_CICE2)

!-----------------------------------------------------------------------
   end subroutine PADCSWAN_FINAL
!-----------------------------------------------------------------------
#endif

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
end module Couple2Swan
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
