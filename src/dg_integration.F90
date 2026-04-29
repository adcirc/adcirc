
module dg_integration
   use sizes, only: MNE
   use NodalAttributes, only: GeoidOffset, LoadGeoidOffset
   use global, only: noff, nodecode, uu1, vv1, qtime1
   use mesh, only: NM
   use DG, only: ZE, RHS_ZE, NEDSD, NEDEL, ATVD, BTVD, DTVD, NEEDN, QNPH_DG, QNAM_DG, &
                 NFEDN, WDFLG, COSNX, SINNX, XLEN, MAX_BOA_DT, neled, hb, nedno, u_modal, &
                 v_modal, niedn, phi_corner, efa_dg, emo_dg, nfeds, pa, dofh, needs, edgeq, xegp, wegp, &
                 m_inv, phi_edge, phi_area, xfac, yfac, bathed, sfaced, negp, bath, srfac, ncele, nagp, &
                 bath, dbathdx, dbathdy, sfac_elem, nrk, leq, nieds, nleq, prep_DG
   use ADC_CONSTANTS, only: G

   implicit none

   private
   integer, parameter :: sz = 8
   real(sz) :: etratio, rampdg

   public :: DG_HYDRO_TIMESTEP
contains

   subroutine DG_HYDRO_TIMESTEP(IT, timeh)
      !! Set `ETA1 := ETA2`, and
      !! compute `ETA2` at the next timestep using DG formulation

      use SLOPELIMITERS, only: SLOPELIMITER
#ifdef CMPI
      use GLOBAL, only: dumy1, dumy2, &
                        DTDP, STATIM, RampExtFlux, NRAMP, DRampExtFlux, &
                        DRAMP, NFFR, NBFR, FTIMINC, QNIN1, QNIN2, &
                        ESBIN1, ESBIN2, ETA2, ETA1, noff
#else
      use GLOBAL, only: DTDP, STATIM, RampExtFlux, NRAMP, DRampExtFlux, &
                        DRAMP, NFFR, NBFR, FTIMINC, QNIN1, QNIN2, &
                        ESBIN1, ESBIN2, ETA2, ETA1, noff
#endif
      use SIZES, only: MNE
      use BOUNDARIES, only: NVEL, LBCODEI, NFLUXF, NOPE, NETA, NBD
      use GWCE, only: ETIME1, ETIME2, ETIMINC
#ifdef CMPI
      use MESSENGER_ELEM, only: updater_elem_mod
      use messenger, only: updateR, updatei
#endif


      implicit none
      integer, intent(in) :: IT
      !! Current time step number
      real(8), intent(in) :: timeh
      !! Current time in seconds (including reference time)

      integer :: timestepper, NQEDS
      integer :: I, J, K, KK, NBDI, IRK
      real(SZ) :: ARK, BRK, time_a, timedg, qtratio

      if (it == 1) then
         call prep_DG()
      end if

      TIME_A = real(IT,8)*DTDP + STATIM*86400.d0

      eta1 = eta2

      call projectMomentum()
      call positive_depth()
#ifdef CMPI
      call UPDATER(UU1, VV1, DUMY2, 2)
      call UPDATER_elem_mod(ze, ze, ze, 1, 1)
#endif
      call update_ncele()
      WDFLG = noff

#ifndef NDEBUG
      call check_bathy(IT)
      call check_element_depth(IT)
      call check_edge_depth(IT)
#endif

!.....Begin RK time stepper

      timestepper = 1
      do IRK = 1, NRK

         TIMEDG = TIME_A - DTDP

!.....Compute the ramping

         RAMPDG = 1.d0
         RAMPExtFlux = 1.0d0
         if (NRAMP >= 1) then
            if (NRAMP == 1) then
               RAMPDG = tanh((2.d0*((real(IT,8) - 1) + DTVD(IRK))*DTDP/86400.d0)/DRAMP)
               RAMPExtFlux = tanh((2.d0*((real(IT,8) - 1) + DTVD(IRK))*DTDP/86400.d0)/DRAMPExtFlux)
            end if
            if (NRAMP == 2) then
               RAMPDG = tanh((2.d0*((real(IT,8) - 1) + DTVD(IRK))*DTDP/86400.d0)/DRAMP)
               RAMPExtFlux = tanh((2.d0*((real(IT,8) - 1) + DTVD(IRK))*DTDP/86400.d0)/DRAMPExtFlux)
            end if
            if (NRAMP == 3) then
               write (*, *) 'NRAMP = 3 not supported '
               stop
            end if
         end if

!.......For non-periodic flux bcs (fort.20)
         ! Copy this from timestep.F
         if (NFLUXF == 1) then
            if ((NFFR == 0) .or. (NFFR == -1)) then
               QTRATIO = (TIMEDG - QTIME1)/FTIMINC
               NQEDS = 0
               do I = 1, NVEL
                  if ((LBCODEI(I) == 2) .or. (LBCODEI(I) == 12) .or. (LBCODEI(I) == 22)) then
                     NQEDS = NQEDS + 1
                     if (NQEDS <= NFEDS) then
                        QNAM_DG(1, NQEDS, 1) = RAMPDG*(QNIN1(I) + QTRATIO*(QNIN2(I) - QNIN1(I)))
                        QNPH_DG(1, NQEDS, 1) = 0.d0
                        QNAM_DG(1, NQEDS, 2) = RAMPDG*(QNIN1(I + 1) + QTRATIO*(QNIN2(I + 1) - QNIN1(I + 1)))
                        QNPH_DG(1, NQEDS, 2) = 0.d0
                     end if
                  end if
               end do
            end if

         end if ! NFLUXF

!........For non-periodic elevation BCs (fort.19). Will be used in ocean_edge_hydro

         if ((NBFR == 0) .and. (NOPE > 0)) then
            if (TIME_A > ETIME2) then
               ETIME1 = ETIME2
               ETIME2 = ETIME1 + ETIMINC
               do J = 1, NETA
                  ESBIN1(J) = ESBIN2(J)
                  read (19, *) ESBIN2(J)
               end do
            end if
            ETRATIO = (TIMEDG - ETIME1)/ETIMINC
            do I = 1, NETA
               NBDI = NBD(I)
               ETA2(NBDI) = RAMPDG*(ESBIN1(I) + ETRATIO*(ESBIN2(I) - ESBIN1(I)))
            end do
         end if
         if (LoadGeoidOffset) then
            do I = 1, NETA
               ETA2(NBD(I)) = ETA2(NBD(I)) + GeoidOffset(NBD(I))
            end do
         end if

         if (NFEDS > 0) call FLOW_EDGE_HYDRO(irk, timedg)

         if (NEEDS > 0) call OCEAN_EDGE_HYDRO(IRK, timedg, rampdg)

         call INTERNAL_EDGE_HYDRO(IRK)

         call RHS_DG_HYDRO(IRK)

!.......SSP Runge-Kutta Time Scheme

         do I = 1, IRK
            ARK = ATVD(IRK, I)
            BRK = BTVD(IRK, I)*DTDP
            do J = 1, MNE
               do K = 1, DOFH
                  ZE(K, J, IRK + 1) = ZE(K, J, irk + 1) + ARK*ZE(K, J, I) + BRK*RHS_ZE(K, J, I)
               end do
            end do
         end do

#ifdef CMPI
         call updater_elem_mod(ZE, ZE, ZE, IRK + 1, 1)
#endif

         call slopelimiter(IRK)

#ifdef CMPI
         call updater_elem_mod(ZE, ZE, ZE, IRK + 1, 1)
#endif
      end do ! IRK = NRK

!.....RK stage calculations done. Set the new state as the last RK stage
      ! Also store the previous timestep in ze0
      do J = 1, MNE
         do K = 1, DOFH
            ZE(K, J, 1) = ZE(K, J, NRK + 1)
         end do
      end do

!.....Zero out the RK stage arrays
      do KK = 2, NRK + 1
         do J = 1, MNE
            do K = 1, DOFH
               ZE(K, J, KK) = 0.d0
               RHS_ZE(K, J, KK - 1) = 0.d0
            end do
         end do
      end do

      call write_results()

#ifdef CMPI
      call UPDATER(ETA2, DUMY1, DUMY2, 1)
#endif

      call computeOceanPressure(timeh, .false.)

   end subroutine DG_HYDRO_TIMESTEP

   subroutine FLOW_EDGE_HYDRO(irk, timedg)

      use sizes, only: mnffr
      use mesh, only: areas
      use global, only: nffr, fper, famig, fface, fff

      implicit none

      integer, intent(in) :: irk
      real(sz), intent(in) :: timedg

      integer :: L, LED, GED, i, k, jj,   el_in, ncyc
      real(sz) :: q_n_ext, q_t_ext, argj, rff, qnam_gp, qnph_gp, arg
      real(SZ) :: TX, TY
      real(SZ) :: W_IN
      real(sz) :: ze_in, ze_ex, qx_in, qx_ex, qy_in, qy_ex, &
                  hb_in, hb_ex, sfac_in, sfac_ex, nx, ny
      real(sz) :: f_hat

      do 1000 L = 1, NFEDS

!.....Retrieve the global and local edge number

         GED = NFEDN(L)
         LED = NEDSD(1, GED)

!.....Retrieve the element to which the edge belongs

         EL_IN = NEDEL(1, GED)

!.....If the element is dry then skip the edge calculation

         if (WDFLG(EL_IN) == 0) goto 1000

!.....Retrieve the components of the normal vector to the edge

         NX = COSNX(GED)
         NY = SINNX(GED)

!.....Set the components for the tangential vector to the edge

         TX = -NY
         TY = NX

!.....Compute ZE, QX, QY, and HB at each edge Gauss quadrature point

         do I = 1, 3

            ZE_IN = 0.D0
            QX_IN = 0.D0
            QY_IN = 0.D0

            ZE_EX = 0.D0
            QX_EX = 0.D0
            QY_EX = 0.D0

            HB_IN = 0.D0
            SFAC_IN = SFACED(I, LED, EL_IN, pa)

!.....Compute the specified flow boundaries for the exterior state

            Q_N_EXT = 0.D0
            do JJ = 1, MNFFR
               if (NFFR == 0) then
                  ARGJ = 0.d0
                  RFF = RAMPDG
               elseif (abs(FPER(JJ)) <= 1d-15) then
                  NCYC = 0
                  ARGJ = FAMIG(JJ)*(TIMEDG - real(NCYC,8)*FPER(JJ)) + FFACE(JJ)
                  RFF = FFF(JJ)*RAMPDG
               else
                  NCYC = int(TIMEDG/FPER(JJ))
                  ARGJ = FAMIG(JJ)*(TIMEDG - real(NCYC,8)*FPER(JJ)) + FFACE(JJ)
                  RFF = FFF(JJ)*RAMPDG
               end if

               QNAM_GP = 0.5d0*(QNAM_DG(JJ, L, 1) + QNAM_DG(JJ, L, 2)) + &
                         0.5d0*(QNAM_DG(JJ, L, 2) - QNAM_DG(JJ, L, 1))*XEGP(I, PA)
               QNPH_GP = 0.5d0*(QNPH_DG(JJ, L, 1) + QNPH_DG(JJ, L, 2)) + &
                         0.5d0*(QNPH_DG(JJ, L, 2) - QNPH_DG(JJ, L, 1))*XEGP(I, PA)

               ARG = ARGJ - QNPH_GP

               Q_N_EXT = Q_N_EXT + QNAM_GP*RFF*cos(ARG)
               Q_T_EXT = 0.d0

               QX_EX = -(TY*Q_N_EXT - NY*Q_T_EXT)/(NX*TY - NY*TX)
               QY_EX = -(-TX*Q_N_EXT + NX*Q_T_EXT)/(NX*TY - NY*TX)

            end do

!.....Compute the solution at the interior state

            do K = 1, 3
               ZE_IN = ZE_IN + ZE(K, EL_IN, IRK)*PHI_EDGE(K, I, LED, pa)
               HB_IN = HB_IN + HB(K, EL_IN, 1)*PHI_EDGE(K, I, LED, pa)
            end do

            qy_ex = qy_ex/(hb_in + ze_in)
            qx_ex = qx_ex/(hb_in + ze_in)

!.....Set the exterior bed and surface elevation equal to the interior

            ZE_EX = ZE_IN
            HB_EX = HB_IN
            SFAC_EX = SFAC_IN

            f_hat = llf_flux(ZE_IN, ZE_EX, HB_IN, HB_EX, qx_ex, qy_ex, &
                             qx_ex, qy_ex, NX, NY, SFAC_IN, SFAC_EX)

!.....Compute the edge integral

            do K = 1, 3

               W_IN = 2.d0*M_INV(K, pa)/AREAS(EL_IN)*XLEN(GED)* &
                      PHI_EDGE(K, I, LED, pa)*WEGP(I, pa)

               RHS_ZE(K, EL_IN, IRK) = RHS_ZE(K, EL_IN, IRK) - W_IN*F_HAT

            end do

         end do

1000     continue

         end subroutine FLOW_EDGE_HYDRO
!
!     SUBROUTINE INTERNAL_EDGE_HYDRO( )
!
!     This subroutine does the following:
!
!     1.  Calculates the values of the necessary variables at the edge
!     gauss points for INTERNAL edges
!     2.  Calls the appropriate subroutine to compute the flux at
!     these points.
!     3.  Calls the appropriate subroutine to compute the boundary
!     integrals.
!
!     Written by Ethan Kubatko (06-11-2004)
!
!     01-10-2011 - cem - adapted for p_enrichment and multicomponent
!     11-11-2011 - cem - adapted for layered sediment
!
!-----------------------------------------------------------------------
!
!     01-02-2007, sb, Modified for LDG
!     C***********************************************************************

         subroutine INTERNAL_EDGE_HYDRO(IRK)

!.....Use appropriate modules

            use mesh, only: AREAS
            implicit none

            integer, intent(in) :: IRK
      !! Current RK stage

            real(sz) :: ze_ex, hb_ex, sfac_ex, ze_in
            real(sz) :: sfac_in, hb_in,  nx, ny
            integer :: el_in, el_ex, el
            integer :: n1, n2
            real(sz) :: U_EDGE, V_EDGE, f_hat
            integer :: L, LED_IN, LED_EX, GED, GP_IN, GP_EX, k, i
            !REAL(SZ), PARAMETER :: ZERO = 1.D-12
            real(SZ) :: TX, TY, W_IN, W_EX
            real(SZ) :: EDFAC_IN, EDFAC_EX
            real(SZ) :: XLEN_EL_IN, XLEN_EL_EX
            real(SZ) :: MASS_EL_IN, MASS_EL_EX


            do L = 1, NIEDS

!.......Retrieve the global and local edge number

               GED = NIEDN(L)
               LED_IN = NEDSD(1, GED)
               LED_EX = NEDSD(2, GED)

!.......Retrieve the elements which share the edge

               EL_IN = NEDEL(1, GED)
               EL_EX = NEDEL(2, GED)

               EL = EL_EX

!.......If both elements on either side of edge are dry then skip

               wet: if ((WDFLG(EL_IN) == 1) .or. (WDFLG(EL_EX) == 1)) then

!.....Compute the sum of the lengths of three edges

                  XLEN_EL_IN = XLEN(NELED(1, EL_IN))
                  XLEN_EL_IN = XLEN_EL_IN + XLEN(NELED(2, EL_IN))
                  XLEN_EL_IN = XLEN_EL_IN + XLEN(NELED(3, EL_IN))

                  XLEN_EL_EX = XLEN(NELED(1, EL_EX))
                  XLEN_EL_EX = XLEN_EL_EX + XLEN(NELED(2, EL_EX))
                  XLEN_EL_EX = XLEN_EL_EX + XLEN(NELED(3, EL_EX))

!.....Compute the total mass in the elements

                  MASS_EL_IN = (ZE(1, EL_IN, IRK) + HB(1, EL_IN, 1))*AREAS(EL_IN)*0.5d0
                  MASS_EL_EX = (ZE(1, EL_EX, IRK) + HB(1, EL_EX, 1))*AREAS(EL_EX)*0.5d0

!.....Retrieve the components of the normal vector to the edge

                  NX = COSNX(GED)
                  NY = SINNX(GED)

                  N1 = NEDNO(1, GED)
                  N2 = NEDNO(2, GED)

!.....Set the components for the tangential vector to the edge

                  TX = -NY
                  TY = NX

                  EDFAC_IN = XLEN(GED)/AREAS(EL_IN)
                  EDFAC_EX = XLEN(GED)/AREAS(EL_EX)

!.....Compute ZE, QX, QY, and HB at each edge Gauss quadrature point
! namo - for now, use constant velocities across edge, obtained
                  ! by averaging the values at the 2 nodes of the edge

                  do I = 1, 3

                     GP_IN = I
                     GP_EX = NEGP(pa) - I + 1

                     HB_IN = BATHED(GP_IN, LED_IN, EL_IN, pa)
                     SFAC_IN = SFACED(GP_IN, LED_IN, EL_IN, pa)

                     HB_EX = HB_IN
                     SFAC_EX = SFACED(GP_EX, LED_EX, EL_EX, pa)

                     ZE_IN = 0D0
                     ZE_EX = 0D0

                     U_EDGE = 0D0
                     V_EDGE = 0D0
                     do K = 1, 3
                        U_EDGE = U_EDGE + U_modal(K, EL_IN)*PHI_EDGE(K, GP_IN, LED_IN, pa)
                        V_EDGE = V_EDGE + V_modal(K, EL_IN)*PHI_EDGE(K, GP_IN, LED_IN, pa)
                        ZE_IN = ZE_IN + ZE(K, EL_IN, IRK)*PHI_EDGE(K, GP_IN, LED_IN, pa)
                        ZE_EX = ZE_EX + ZE(K, EL_EX, IRK)*PHI_EDGE(K, GP_EX, LED_EX, pa)
                     end do

!DIR$ FORCEINLINE
                     f_hat = llf_flux(ZE_IN, ZE_EX, HB_IN, HB_EX, U_EDGE, V_EDGE, &
                                      U_EDGE, V_EDGE, NX, NY, SFAC_IN, SFAC_EX)

!.....Check if the flux is large enough to dry up the elements
!.....1.01D0 is a safty factor.

                     if ((1.01d0*F_HAT*XLEN_EL_IN*MAX_BOA_DT(IRK) >= MASS_EL_IN) &
                         .or. (1.01d0*F_HAT*XLEN_EL_EX*MAX_BOA_DT(IRK)*(-1.d0) >= &
                               MASS_EL_EX)) then

                        cycle
                     end if

!........Check to make sure mass flux is not coming from a dry element
                     if (abs(f_hat) > 1.d-12) then
!
                        if (wdflg(el_in) == 0) then
! el_in is dry !
                           if (f_hat > 0) then
! flux going from the dry element (in)
! on the wet side (ex): reflect boundary
                              cycle
                           end if

                        elseif (wdflg(el_ex) == 0) then

! el_ex is dry
                           if (f_hat < 0) then
! flux comming from dry size (ex)
! on the wet side (in): reflect boundary
                              cycle
                           end if
                        end if
                     end if
                     do K = 1, 3

                        W_IN = EDFAC_IN*EDGEQ(K, GP_IN, LED_IN, pa)
                        W_EX = EDFAC_EX*EDGEQ(K, GP_EX, LED_EX, pa)

                        RHS_ZE(K, EL_IN, IRK) = RHS_ZE(K, EL_IN, IRK) - W_IN*F_HAT
                        RHS_ZE(K, EL_EX, IRK) = RHS_ZE(K, EL_EX, IRK) + W_EX*F_HAT
                     end do
                  end do
               end if wet
            end do

         end subroutine INTERNAL_EDGE_HYDRO

         subroutine RHS_DG_HYDRO(IRK)
!! This subroutine computes the area integrals for the DG hydro and
!! adds them into the RHS. For each element E, these terms appear as
!! $$(\nabla v, F)_{E} - \langle \hat{F}_n, v \rangle_{\partial E}$$

            use GLOBAL, only: dtdp
            use sizes, only: MNE
            use precipitation, only: elem_rain

            implicit none

            integer, intent(in) :: IRK

            integer :: L, k, i
            real(SZ) :: source_r
            real(SZ) :: DEPTH, F1_NL
            real(SZ) :: DTDPH, SFACQUAD

            real(sz) :: ze_in, hb_in, dhb_x, dhb_y
            real(sz) :: fx_in, fy_in
            real(sz) :: u_quad, v_quad

            DTDPH = 1.d0/DTDP

!$omp simd
            do L = 1, MNE

!.......namo - rain source before wet/dry check

               source_r = elem_rain(L)
               do K = 1, 3
                  do I = 1, 3
                     RHS_ZE(K, L, IRK) = RHS_ZE(K, L, IRK) + SRFAC(K, I, L, pa)*SOURCE_R
                  end do
               end do

!.......If element is dry then skip calculations

               wet: if (WDFLG(L) == 1) then

!.......Compute ZE, QX, QY, and HB at each area Gauss quadrature point

                  do I = 1, 3

                     U_QUAD = 0D0
                     V_QUAD = 0D0

                     ZE_IN = 0D0
                     HB_IN = BATH(I, L, pa)
                     DHB_X = DBATHDX(I, L, pa)
                     DHB_Y = DBATHDY(I, L, pa)

                     SFACQUAD = SFAC_ELEM(I, L, pa)

                     do K = 1, 3
                        U_QUAD = U_QUAD + U_modal(K, L)*PHI_AREA(K, I, pa)
                        V_QUAD = V_QUAD + V_modal(K, L)*PHI_AREA(K, I, pa)
                        ZE_IN = ZE_IN + ZE(K, L, IRK)*PHI_AREA(K, I, pa)

                     end do

                     DEPTH = ZE_IN + HB_IN

!.........Compute continuity fluxes

                     F1_NL = NLEQ + LEQ*HB_IN

!FX_IN = QX_IN*F1_NL*SFACQUAD
                     FX_IN = U_QUAD*HB_IN*SFACQUAD*LEQ + U_QUAD*DEPTH*SFACQUAD*NLEQ
!     for linear swe
                     FY_IN = V_QUAD*HB_IN*LEQ + V_QUAD*DEPTH*NLEQ

                     do K = 2, 3

                        RHS_ZE(K, L, IRK) = RHS_ZE(K, L, IRK) + XFAC(K, I, L, pa)*FX_IN &
                                            + YFAC(K, I, L, pa)*FY_IN
                     end do
                  end do
               end if wet
            end do

            return
         end subroutine RHS_DG_HYDRO

         subroutine ocean_edge_hydro( irk, timedg, rampdg)

! adcirc new stuff
            use boundaries, only: nope
            use mesh, only: areas
            use GLOBAL, only: NBFR, PER, AMIG, FF, H0, IFNLFA, FACE, ETA2

            implicit none

      !! Current RK stage
            integer, intent(in) :: irk
            real(sz), intent(in) :: timedg, rampdg

!.....Declare local variables

            real(sz) :: u_edge, v_edge, nx, ny, ze_in, ze_ex, hb_in, hb_ex
            real(sz) :: f_hat, sfac_ex, sfac_in, w_in
            integer ::  el_in, n1, n2

            integer :: L, LED, GED, i, k, jj
            integer:: IPT
            real(SZ):: ZEFREQ
            real(SZ), dimension(2):: EFA_GPT, EMO_GPT, ARG_GPT, ZE_GPT
            integer :: ncyc
            real(sz) :: argj, rff

            do L = 1, needs

!.....Retrieve the global and local edge number

               GED = NEEDN(L)
               LED = NEDSD(1, GED)

!.....Retrieve the elements which share the edge

               EL_IN = NEDEL(1, GED)

!.....If the element is dry then skip the edge calculation

               wet: if (WDFLG(EL_IN) == 1) then

!.....Retrieve the components of the normal vector to the edge

                  NX = COSNX(GED)
                  NY = SINNX(GED)

!.....Retrieve the nodes of the edge

                  N1 = NEDNO(1, GED)
                  N2 = NEDNO(2, GED)

!.....Compute ZE, QX, QY, and HB at each edge Gauss quadrature point

                  do I = 1, 3
                     ZE_EX = 0.d0

                     HB_IN = BATHED(I, LED, EL_IN, pa)
                     SFAC_IN = SFACED(I, LED, EL_IN, pa)

!.....Compute the specified open ocean elevation

                     do JJ = 1, NBFR

                        if (abs(PER(JJ)) <= 1d-15) then
                           NCYC = 0
                        else
                           NCYC = int(TIMEDG/PER(JJ))
                        end if

!...........Surface Elevation

                        ARGJ = AMIG(JJ)*(TIMEDG - real(NCYC,8)*PER(JJ)) + FACE(JJ)
                        RFF = FF(JJ)*RAMPDG

!..............linearly interpolate from the time-series of the harmonic forcing
                        do IPT = 1, 2
                           EFA_GPT(IPT) = EFA_DG(JJ, L, IPT); 
                           EMO_GPT(IPT) = EMO_DG(JJ, L, IPT); 
                           ARG_GPT(IPT) = ARGJ - EFA_GPT(IPT); 
                           ZE_GPT(IPT) = EMO_GPT(IPT)*RFF*cos(ARG_GPT(IPT)); 
                        end do

                        ZEFREQ = 0.5d0*(ZE_GPT(1) + ZE_GPT(2)) + &
                                 0.5d0*(ZE_GPT(2) - ZE_GPT(1))*XEGP(I, pa); 
                        ZE_EX = ZE_EX + ZEFREQ; 
                     end do

!.....For non-periodic elevation bcs

                     if ((NBFR == 0) .and. (NOPE > 0)) then
                        ZE_EX = 0.5d0*(ETA2(N1) + ETA2(N2)) &
                                + 0.5d0*(ETA2(N2) - ETA2(N1))*XEGP(I, pa)
                     end if

!.....Compute the solution at the interior state

                     ZE_IN = 0D0
                     U_EDGE = 0D0
                     V_EDGE = 0D0
                     do K = 1, 3

                        U_EDGE = U_EDGE + U_modal(K, EL_IN)*PHI_EDGE(K, I, LED, pa)
                        V_EDGE = V_EDGE + V_modal(K, EL_IN)*PHI_EDGE(K, I, LED, pa)
                        ZE_IN = ZE_IN + ZE(K, EL_IN, IRK)*PHI_EDGE(K, I, LED, pa)

                     end do

!.....Set the exterior value of the bathymetry equal to the interior

                     HB_EX = HB_IN
                     SFAC_EX = SFAC_IN

!$$$            IF (LoadGeoidOffset) then
!$$$               ZE_EX = ZE_EX + .5*(GeoidOffset(N1)+GeoidOffset(N2))
!$$$            endif

                     ! Eirik's fix
                     if ((ZE_EX*real(IFNLFA,8) + HB_EX) <= 0.d0) then
                        ZE_EX = abs(HB_EX) + H0
                     end if
!DIR$ FORCEINLINE
                     F_hat = llf_flux(ZE_IN, ZE_EX, HB_IN, HB_EX, U_EDGE, V_EDGE, &
                                      U_EDGE, V_EDGE, NX, NY, SFAC_IN, SFAC_EX)

!.....Compute the edge integral
                     do K = 1, 3

                        W_IN = 2.0d0*M_INV(K, pa)/AREAS(EL_IN)*XLEN(GED)* &
                               PHI_EDGE(K, I, LED, pa)*WEGP(I, pa)

                        RHS_ZE(K, EL_IN, IRK) = RHS_ZE(K, EL_IN, IRK) - W_IN*F_HAT

                     end do
                  end do

               end if wet
            end do

         end subroutine ocean_edge_hydro

!
!     SUBROUTINE:  WRITE_RESULTS
!
!     Taken from the timestep subroutine.  Modified to include
!     interpolation of DG modal degrees of freedom to the nodes.  These
!     multiple nodal values are then averaged to a single value.
!
!     Aug ??, 2005, sb, Modifications for parallel runs
!     Jan 01, 2007, sb, Files are forced to be written if FORCE_WRITE = .TRUE.
!
!***********************************************************************

         subroutine WRITE_RESULTS()

!.....Use appropriate modules

            use GLOBAL, only: etamax, eta2
            use MESH, only: NM, AREAS
            use sizes, only: MNP


            integer ::   kk,  i, n1, n2, n3
            real(sz) ::  ze1, ze2, ze3
            real(sz) :: node_area(MNP), node_ze(MNP)

!.....Transform from modal coordinates to nodal coordinates and average
!.....to single nodal values
            node_area = 0.d0
            node_ze = 0.d0
            do I = 1, MNE
               if (ncele(I) == 1) then
                  N1 = NM(I, 1)
                  N2 = NM(I, 2)
                  N3 = NM(I, 3)

                  node_area(n1) = node_area(n1) + 0.5d0*areas(I)
                  node_area(n2) = node_area(n2) + 0.5d0*areas(I)
                  node_area(n3) = node_area(n3) + 0.5d0*areas(I)

                  ze1 = 0
                  ze2 = 0
                  ze3 = 0
                  do KK = 1, DOFH
                     ZE1 = ZE1 + PHI_CORNER(KK, 1, 1)*ZE(KK, I, 1)
                     ZE2 = ZE2 + PHI_CORNER(KK, 2, 1)*ZE(KK, I, 1)
                     ZE3 = ZE3 + PHI_CORNER(KK, 3, 1)*ZE(KK, I, 1)
                  end do

                  node_ze(N1) = node_ze(n1) + ze1*0.5d0*areas(i)
                  node_ze(N2) = node_ze(n2) + ze2*0.5d0*areas(i)
                  node_ze(N3) = node_ze(n3) + ze3*0.5d0*areas(i)
               end if
            end do

!$omp simd
            do I = 1, MNP
               if (node_area(i) > 0) then
                  eta2(i) = node_ze(i)/node_area(i)
               else
                  eta2(i) = 0.d0
               end if
               etamax(i) = max(etamax(i), eta2(i))
            end do

         end subroutine WRITE_RESULTS

         subroutine projectMomentum()
!! Convert nodal `UU2`,`VV2` into DG modal representation
!!`U_modal`, `V_modal`

            use global, only: UU2, VV2
            use mesh, only: NM
            use sizes, only: MNE
            implicit none

            integer :: N1, N2, N3, J
            real(8) :: u1, u2, u3, v1, v2, v3

            do J = 1, MNE
               N1 = NM(J, 1)
               N2 = NM(J, 2)
               N3 = NM(J, 3)

               u1 = UU2(N1)
               u2 = UU2(N2)
               u3 = UU2(N3)

               v1 = VV2(N1)
               v2 = VV2(N2)
               v3 = VV2(N3)

               U_modal(1, J) = (1.d0/3.d0*(u1 + u2 + u3))
               U_modal(2, J) = (-1.d0/6.d0*(u1 + u2) + 1.d0/3.d0*u3)
               U_modal(3, J) = (-0.5d0*u1 + 0.5d0*u2)

               V_modal(1, J) = (1.d0/3.d0*(v1 + v2 + v3))
               V_modal(2, J) = (-1.d0/6.d0*(v1 + v2) + 1.d0/3.d0*v3)
               V_modal(3, J) = (-0.5d0*v1 + 0.5d0*v2)

            end do

! Rearrange loops

         end subroutine projectMomentum

!     SUBROUTINE COMPUTE OCEAN PRESSURE
!--------------------------------------
!     Compute the elevation at the ocean boundary nodes specifically
!     for use in computing the barotropic pressure in the momentum equations.
!     Note that ETA2 from DG itself is not changed, but the pressure terms
!     requires exact values at the ocean boundary (and not just weak enforcement).
!--------------------------------------
!     Modifies: PETA1, PETA2 - these are ETA1 and ETA2 except that the
!               ocean boundary terms are exact as given
!******************************************************************************
         subroutine computeOceanPressure(timeh, forceFlag)
            use global, only: peta1, peta2, PER, FACE, FF, EMO, EFA, rampElev, &
                              nbfr, AMIG, ETA2, ESBIN1, ESBIN2
            use boundaries, only: NETA, NBD, NOPE
            use mesh, only: NM

            real(sz), intent(in) :: timeh
            logical, intent(in) :: forceFlag

            ! Local vars
            integer :: NCYC, NBDI, J, I, n1, n2, n3, L, LED, GED
            real(sz) :: ARG, ARGJ, RFF

!.....If we have fort.19
            if ((NBFR == 0) .and. (NOPE > 0)) then
               PETA1(:) = PETA2(:)
               PETA2(:) = ETA2(:)

               do I = 1, NETA
                  NBDI = NBD(I)
                  PETA2(NBDI) = rampdg*(ESBIN1(I) + ETRATIO*(ESBIN2(I) - ESBIN1(I)))
               end do

            else

!.....Else if we use periodic elevation BC

               ! Save the current state and set it to ETA2 initially
               Peta1(:) = Peta2(:)
               Peta2(:) = ETA2(:)

               ! Zero the elevation-specified nodes
               do J = 1, NETA
                  PETA2(NBD(J)) = 0.d0
               end do

               !Compute the elevation-specified values at timeh
               do J = 1, NBFR
                  if (abs(PER(J)) <= 1d-15) then
                     NCYC = 0
                  else
                     NCYC = int(timeh/PER(J))
                  end if
                  ARGJ = AMIG(J)*(timeh - real(NCYC,8)*PER(J)) + FACE(J)
                  RFF = FF(J)*RampElev
                  do I = 1, NETA
                     ARG = ARGJ - EFA(J, I)
                     NBDI = NBD(I)
                     PETA2(NBDI) = PEta2(NBDI) + EMO(J, I)*RFF*cos(ARG)
                  end do
               end do

            end if
            if (LoadGeoidOffset) then
               do I = 1, NETA
                  PETA2(NBD(I)) = PETA2(NBD(I)) + GeoidOffset(NBD(I))
               end do
            end if
            ETA2(:) = PETA2(:)

!     if forceFlag, enforce the actual ETA2 and modify the DG basis to
            ! have the same values
            if (forceFlag) then
               ETA2(:) = PETA2(:)

               ! loop through the elements at the ocean boundary
               do L = 1, needs
!.....Retrieve the global and local edge number
                  GED = NEEDN(L)
                  LED = NEDSD(1, GED)

!.....Retrieve the element
                  J = NEDEL(1, GED)

                  N1 = NM(J, 1)
                  N2 = NM(J, 2)
                  N3 = NM(J, 3)

                  ZE(1, J, 1) = 1.d0/3.d0*(ETA2(N1) + ETA2(N2) + ETA2(N3))
                  ZE(2, J, 1) = -1.d0/6.d0*(ETA2(N1) + ETA2(N2)) + 1.d0/3.d0*ETA2(N3)
                  ZE(3, J, 1) = -0.5d0*ETA2(N1) + 0.5d0*ETA2(N2)

               end do
            end if

         end subroutine computeOceanPressure

         subroutine check_element_depth(it)
!! Loop through elements and check if the depth at any AREA
!! quadrature point is negative, in which case stop the program.

            use sizes, only: MNE
            implicit none

            integer, intent(in) :: it
            integer :: l, i, k
            real(sz) :: ze_in, hb_in, depth

            do l = 1, MNE
               if (.true.) then
                  do I = 1, NAGP(pa)

                     ze_in = 0.d0
                     HB_IN = BATH(I, L, pa)

                     do k = 1, DOFH
                        ZE_IN = ZE_IN + ZE(K, L, 1)*PHI_AREA(K, I, pa)
                     end do

                     depth = ze_in + hb_in

                     if (depth <= 0) then
                        print *, 'negative element depth at timestep ', it
                        print *, 'element ', l
                        print *, 'ze = ', ze_in
                        print *, 'hb_in = ', hb_in
                        stop 5
                     end if

                  end do
               end if
            end do

         end subroutine check_element_depth

         subroutine check_edge_depth(it)
!! Loop through internal edges and check if the depth at any EDGE
!! quadrature point is negative, in which case stop the program.

            implicit none

            integer, intent(in), value :: it
            real(sz) :: depth_in, depth_ex
            real(sz) :: ze_ex, hb_ex, sfac_ex, ze_in
            real(sz) :: sfac_in, hb_in
            !real(sz) q_n_ext,q_t_ext,q_n_int,q_t_int
            integer :: el_in, el_ex
! function
            !real(sz) llf_flux_coupling

            integer :: L, LED_IN, LED_EX, GED, GP_IN, GP_EX, k, i
            !REAL(SZ), PARAMETER :: ZERO = 1.D-12

            do L = 1, NIEDS

!.......Retrieve the global and local edge number

               GED = NIEDN(L)
               LED_IN = NEDSD(1, GED)
               LED_EX = NEDSD(2, GED)

!.......Retrieve the elements which share the edge

               EL_IN = NEDEL(1, GED)
               EL_EX = NEDEL(2, GED)

               do I = 1, NEGP(pa)

                  GP_IN = I
                  GP_EX = NEGP(pa) - I + 1

                  HB_IN = BATHED(GP_IN, LED_IN, EL_IN, pa)
                  SFAC_IN = SFACED(GP_IN, LED_IN, EL_IN, pa)

                  HB_EX = BATHED(GP_EX, LED_EX, EL_EX, pa)
                  SFAC_EX = SFACED(GP_EX, LED_EX, EL_EX, pa)

                  ZE_IN = 0D0
                  ZE_EX = 0D0

                  do K = 1, 3
                     ZE_IN = ZE_IN + ZE(K, EL_IN, 1)*PHI_EDGE(K, GP_IN, LED_IN, pa)
                     ZE_EX = ZE_EX + ZE(K, EL_EX, 1)*PHI_EDGE(K, GP_EX, LED_EX, pa)
                  end do

                  depth_in = ze_in + hb_in
                  depth_ex = ze_ex + hb_ex

                  if (depth_in <= 0 .or. depth_ex <= 0) then
                     print *, 'negative edge depth at timestep ', it
                     print *, 'internal edge ', L
                     print *, 'ze = ', ze_in
                     print *, 'hb_in = ', hb_in
                     stop 5
                  end if
               end do
            end do

         end subroutine check_edge_depth

         subroutine positive_depth()
!! Enforce ZE to have positive depth using the algorithm in
!! Shintaro's 2008 paper. There, it is referred to as the operator \(M\Pi_h\).

            use global, only: NOFF, nodecode, uu1, vv1
            use global, only: H0
            use mesh, only: NM, DP
            use sizes, only: MNE

            implicit none

            integer :: j, kk, k, m1, m2, m3, inds(3)
            real(sz) :: zevertex(3), depth(3), ze_hat(3), depth_avg, depth2
            real(sz) :: H1

            H1 = 2.d0*H0

            do j = 1, MNE
               zevertex = 0.d0

               do KK = 1, 3
                  ZEVERTEX(1) = ZEVERTEX(1) + PHI_CORNER(KK, 1, 1)*ZE(kk, j, 1)
                  ZEVERTEX(2) = ZEVERTEX(2) + PHI_CORNER(KK, 2, 1)*ZE(kk, j, 1)
                  ZEVERTEX(3) = ZEVERTEX(3) + PHI_CORNER(KK, 3, 1)*ZE(kk, j, 1)
               end do

               do k = 1, 3
                  depth(k) = zevertex(k) + DP(NM(j, k))
               end do

               depth_avg = sum(depth)/3.d0

               if (all(depth > H1)) then
                  NOFF(j) = 1
                  nodecode(NM(j, :)) = 1
                  cycle ! move on to the next element
               elseif (depth_avg < 0) then
                  ze_hat(:) = H0 - DP(nm(j, :))
                  NOFF(j) = 0
                  nodecode(NM(j, :)) = 0
                  UU1(NM(j, :)) = 0.d0
                  VV1(NM(j, :)) = 0.d0
               elseif (depth_avg <= H1) then
! If mean value is less than H1, then set the whole element to that depth
                  ze_hat(:) = depth_avg - DP(nm(j, :))
                  !if (LoadGeoidOffset) ze_hat = ze_hat + GeoidOffset(NM(j,1))
                  !ze_hat(:) = H0*1.1 - DP(nm(j,:))
                  NOFF(j) = 0
                  nodecode(NM(j, :)) = 0
                  UU1(NM(j, :)) = 0.d0
                  VV1(NM(j, :)) = 0.d0
               else
                  call sort(3, depth, inds)
                  m1 = inds(1)
                  m2 = inds(2)
                  m3 = inds(3)
                  ze_hat(m1) = H1 - DP(nm(j, m1))

                  ze_hat(m2) = max(H1, depth(2) - 0.5d0*(H1 - depth(1))) - DP(nm(j, m2))
                  depth2 = ze_hat(m2) + dp(nm(j, m2))

                  ze_hat(m3) = depth(3) - (H1 - depth(1)) - (depth2 - depth(2)) - DP(nm(j, m3))
                  UU1(NM(j, :)) = 0.d0
                  VV1(NM(j, :)) = 0.d0
               end if

! Reproject vertex values into DG modes
               ZE(1, J, 1) = 1.d0/3.d0*(ze_hat(1) + ze_hat(2) + ze_hat(3))
               ZE(2, J, 1) = -1.d0/6.d0*(ze_hat(1) + ze_hat(2)) + 1.d0/3.d0*ze_hat(3)
               ZE(3, J, 1) = -0.5d0*ze_hat(1) + 0.5d0*ze_hat(2)

            end do

         end subroutine positive_depth

         subroutine check_bathy(IT)
      !! Check if the bathymetry in the DG basis matches the nodal DP
            use mesh, only: DP, NM
            implicit none

            integer, value :: it
            real(sz) :: vertex(3), dps(3)
            integer :: j, kk, i

            do j = 1, MNE
               vertex = 0.d0
               do i = 1, 3
                  dps(i) = DP(NM(J, i))
               end do

               do KK = 1, 3
                  VERTEX(1) = VERTEX(1) + PHI_CORNER(KK, 1, 1)*HB(kk, j, 1)
                  VERTEX(2) = VERTEX(2) + PHI_CORNER(KK, 2, 1)*HB(kk, j, 1)
                  VERTEX(3) = VERTEX(3) + PHI_CORNER(KK, 3, 1)*HB(kk, j, 1)
               end do

               do i = 1, 3
                  if (abs(vertex(i) - dps(i)) > 1d-8) then
                     print *, 'Bathymetry not matching at element, vertex ', J, i
                     print *, 'hb, dp = ', vertex(i), dps(i)
                     print *, 'at time step ', it
                     stop 6
                  end if
               end do
            end do

         end subroutine check_bathy

         subroutine sort(n, a, is)
!! Sort the input array `a` and write the corresponding indices into `is`.
            implicit none
            integer, intent(in) :: n
            real(sz), intent(inout) :: a(n)
            integer, intent(out) :: is(n)

            integer :: i, j
            real(sz) :: x

            do i = 1, n
               is(i) = i
            end do

            do i = 2, n
               x = a(i)
               j = i - 1
               do while (j >= 1)
                  if (a(j) <= x) exit
                  a(j + 1) = a(j)
                  is(j + 1) = is(j)
                  j = j - 1
               end do
               a(j + 1) = x
               is(j + 1) = i
            end do
         end subroutine sort

         subroutine update_ncele()
            implicit none

            integer :: j, n1, n2, n3

            do j = 1, MNE
               n1 = nm(j, 1)
               n2 = nm(j, 2)
               n3 = nm(j, 3)
               ncele(j) = noff(j)*nodecode(n1)*nodecode(n2)*nodecode(n3)
            end do

         end subroutine update_ncele
!
         pure real(sz) function LLF_FLUX(ZE_IN, ZE_EX, HB_IN, HB_EX, U_IN, V_IN, &
                                         U_EX, V_EX, NX, NY, SFAC_IN, SFAC_EX)

            implicit none

            real(sz), intent(in), value :: ZE_IN, ZE_EX, HB_IN, HB_EX
            real(sz), intent(in), value :: U_IN, V_IN, U_EX, V_EX, NX, NY, SFAC_IN, SFAC_EX

!.....Declare local variables
            real(SZ) :: EIGVALS(6), EIGMAX, JUMP, HT_IN, HT_EX, QX_IN, QY_IN, QX_EX, QY_EX
            real(SZ) :: C_EX, C_IN, F1_NL, FY1_IN, FY1_EX, FX1_IN, FX1_EX, F1_AVG
            real(sz) :: Un_in, Un_ex

!.....Compute the jump in the variables.

            JUMP = ZE_EX - ZE_IN

!.....Compute the total height of the water column

            HT_IN = ZE_IN*NLEQ + HB_IN
            HT_EX = ZE_EX*NLEQ + HB_EX

!.....Compute the momentum flux from the velocities
            QX_IN = HT_IN*U_IN*NLEQ + U_IN*LEQ
            QY_IN = HT_IN*V_IN*NLEQ + V_IN*LEQ
            QX_EX = HT_EX*U_EX*NLEQ + U_EX*LEQ
            QY_EX = HT_EX*V_EX*NLEQ + V_EX*LEQ

!.....Compute continuity fluxes at interior state

            F1_NL = NLEQ + LEQ*HB_IN ! this is just 1 if nolica, nolicat = 1
            FX1_IN = QX_IN*F1_NL*SFAC_IN
            FY1_IN = QY_IN*F1_NL

!.....Compute continuity fluxes at exterior state

            F1_NL = NLEQ + LEQ*HB_EX
            FX1_EX = QX_EX*F1_NL*SFAC_EX
            FY1_EX = QY_EX*F1_NL

!.....Compute the average flux function

            F1_AVG = 0.5d0*((FX1_IN + FX1_EX)*NX + (FY1_IN + FY1_EX)*NY)

!.....Evaluate the eigenvalues at the interior and exterior states

            UN_IN = (U_IN*NX + V_IN*NY)*NLEQ
            C_IN = sqrt(G*HT_IN*(NY**2 + (NX*SFAC_IN)**2)) !srb - spherical coordinate correction
            EIGVALS(1) = abs(UN_IN + C_IN)
            EIGVALS(2) = abs(UN_IN)
            EIGVALS(3) = abs(UN_IN - C_IN)

            UN_EX = (U_EX*NX + V_EX*NY)*NLEQ
            C_EX = sqrt(G*HT_EX*(NY**2 + (NX*SFAC_EX)**2)) !srb - spherical coordinate correction
            EIGVALS(4) = abs(UN_EX + C_EX)
            EIGVALS(5) = abs(UN_EX)
            EIGVALS(6) = abs(UN_EX - C_EX)

!.....Find the maximum eigenvalue (in absolute value)

            EIGMAX = max(EIGVALS(1), EIGVALS(2), EIGVALS(3), &
                         EIGVALS(4), EIGVALS(5), eigvals(6))

!.....Compute the Local Lax Friedrichs Fluxes

            llf_flux = F1_AVG - 0.5d0*EIGMAX*(JUMP)


         end function llf_flux

         end module dg_integration
