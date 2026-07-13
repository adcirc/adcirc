
MODULE DG

   USE SIZES, only: mne, mnp, mneta, mnbfr, mnei, mnffr, mnbfr

   implicit none
   private
   public :: ze, U_modal, V_modal, nedel, wdflg, neled, ncele
   public :: rhs_ze, g2root, slopeflag
   public :: phi_corner, el_count, eletab, leq, nleq, nleqg
   public :: atvd, btvd, dtvd, nrk
   public :: needn, qnph_dg, qnam_dg, MAX_BOA_DT, nedno
   public :: sinnx, cosnx, hb, xlen, nedsd, niedn, nagp, xfac, yfac
   public :: dbathdx, bath, srfac, sfac_elem, phi_area, negp, phi_edge
   public :: wegp, bathed, sfaced, xegp, emo_dg, efa_dg
   public :: M_inv, dbathdy, needs, nieds, dofh, pa
   public :: psi2, psi3, psi1, etiminc_dg
   public :: edgeq, nfeds, qtratio
   public :: nfedn

   protected :: nedel, neled
   protected :: g2root, slopeflag
   protected :: phi_corner, el_count, eletab, leq, nleq, nleqg
   protected :: atvd, btvd, dtvd, nrk
   protected :: needn, MAX_BOA_DT, nedno
   protected :: sinnx, cosnx, hb, xlen, nedsd, niedn, nagp, xfac, yfac
   protected :: dbathdx, bath, srfac, sfac_elem, phi_area, negp, phi_edge
   protected :: wegp, bathed, sfaced, xegp, emo_dg, efa_dg
   protected :: M_inv, dbathdy, needs, nieds, dofh, pa
   protected :: psi2, psi3, psi1, etiminc_dg
   protected :: edgeq, nfeds, qtratio
   protected :: nfedn

   public :: prep_DG, nodal_to_modal, nodal_to_quad_points

   integer, parameter, private :: sz = 8

  !! Listing of public variables

   REAL(SZ), ALLOCATABLE :: ZE(:, :, :)
  !! `ze(i,j,k)` = ith mode of water elevation on element j at RK stage k
   REAL(SZ), ALLOCATABLE :: U_modal(:, :)
  !! `u_modal(i,j)` = ith mode of x-velocity on element j
   real(sz), ALLOCATABLE :: V_modal(:, :)
  !! `v_modal(i,j)` = ith mode of y-velocity on element j
   INTEGER, ALLOCATABLE :: WDFLG(:)
  !! Same as `NOFF(j)`: 1 if element j is wet, 0 if dry
   INTEGER, ALLOCATABLE :: NCELE(:)
  !! `NCELE` = `NOFF*NODECODE(n1)*NODECODE(n2)*NODECODE(n3)`

   integer :: dofh
  !! DG polynomial order (degrees of freedom). Currently only 1 is supported.

   integer :: NAGP(8)
  !! `NAGP(i)` = number of area quadrature points for dofh = i
   integer :: NEGP(8)
  !! `NEGP(i)` = number of edge quadrature points for dofh = i

   INTEGER :: DOF, dofl, dofx
   INTEGER :: EL
   INTEGER :: J1, J2, J3 
   INTEGER :: NCHECK(8), NEDGES, NRK !42 hardwires for ph=7
   INTEGER :: NIEDS, NLEDS, NEEDS, NFEDS, NREDS, NEBEDS, NIBEDS
   INTEGER :: NIBSEG, NEBSEG
   INTEGER :: MNED
   INTEGER, TARGET :: MODAL_IC
   INTEGER, TARGET :: SLOPEFLAG
   INTEGER, TARGET :: RK_STAGE, RK_ORDER
   Integer, TARGET :: padapt, pl, ph, px
   INTEGER :: pa
   !

   !.....Declare real variables

   REAL(SZ) :: C13, C16
   REAL(SZ) :: DOT, DPHIDX, DPHIDY
   REAL(SZ) :: EL_ANG
   REAL(SZ) :: FG_L
   REAL(SZ) :: MAG1, MAG2
   REAL(SZ) :: S1, S2, SAV

   REAL(SZ), ALLOCATABLE :: ATVD(:, :), BTVD(:, :), CTVD(:, :)
   REAL(SZ), ALLOCATABLE :: DTVD(:), MAX_BOA_DT(:)
   !     sb...Wetting and drying
   INTEGER, ALLOCATABLE :: WDFLG_TMP(:)
   INTEGER, ALLOCATABLE :: DOFW(:)
   INTEGER, ALLOCATABLE :: EL_UPDATED(:)
   REAL(SZ), ALLOCATABLE :: LEDGE_NVEC(:, :, :)

   ! LEDGE_NVEC(1,1:3,1:NELEM) = whether a node is on the land boundary
   ! LEDGE_NVEC(2,1:3,1:NELEM) = x component of the normal vector
   ! LEDGE_NVEC(3,1:3,1:NELEM) = y component of the normal vector

   !Declare some stuff for function parsing for bed load


   !.....Declare real variable arrays

   REAL(SZ) :: DRPSI(3), DSPSI(3)
   REAL(SZ) :: VEC1(2), VEC2(2)

   !.....Declare allocatable integer arrays

   INTEGER, ALLOCATABLE :: DOFS(:), PCOUNT(:), PDG(:)
   INTEGER, ALLOCATABLE :: NCOUNT(:)
   INTEGER, ALLOCATABLE :: NEDEL(:, :), NEDNO(:, :), NEDSD(:, :)
   INTEGER, ALLOCATABLE :: NIEDN(:), NLEDN(:), NEEDN(:)
   INTEGER, ALLOCATABLE :: NFEDN(:), NREDN(:), NEBEDN(:), NIBEDN(:)
   INTEGER, ALLOCATABLE :: NIBSEGN(:, :)
   INTEGER, ALLOCATABLE :: NEBSEGN(:)
   INTEGER, ALLOCATABLE :: EL_NBORS(:, :)
   INTEGER, ALLOCATABLE :: BACKNODES(:, :)
   INTEGER, ALLOCATABLE :: MARK(:)

   !.....Declare allocatable real arrays

   REAL(SZ), ALLOCATABLE :: BATH(:, :, :), DBATHDX(:, :, :), DBATHDY(:, :, :)
   REAL(SZ), ALLOCATABLE :: SFAC_ELEM(:, :, :)
   REAL(SZ), ALLOCATABLE :: BATHED(:, :, :, :), SFACED(:, :, :, :)
   REAL(SZ), ALLOCATABLE :: COSNX(:), SINNX(:)
   REAL(SZ), ALLOCATABLE :: DP_NODE(:, :, :)
   REAL(SZ), ALLOCATABLE :: DP_VOL(:, :)
   REAL(SZ), ALLOCATABLE :: DRPHI(:, :, :), DSPHI(:, :, :)
   REAL(SZ), ALLOCATABLE :: DRDX(:), DSDX(:), DRDY(:), DSDY(:)
   REAL(SZ), ALLOCATABLE :: EFA_DG(:, :, :), EMO_DG(:, :, :)
   REAL(SZ), ALLOCATABLE :: UFA_DG(:, :, :), UMO_DG(:, :, :)
   REAL(SZ), ALLOCATABLE :: VFA_DG(:, :, :), VMO_DG(:, :, :)
   REAL(SZ), ALLOCATABLE :: XLEN(:)
   REAL(SZ), ALLOCATABLE :: HB(:, :, :)
   REAL(SZ), ALLOCATABLE :: MANN(:, :)
   REAL(SZ), ALLOCATABLE :: M_INV(:, :)
   REAL(SZ), ALLOCATABLE :: PHI_AREA(:, :, :), PHI_EDGE(:, :, :, :)
   REAL(SZ), ALLOCATABLE :: PHI_CENTER(:, :), PHI_CORNER(:, :, :)
   REAL(SZ), ALLOCATABLE :: PHI_CHECK(:, :, :)
   REAL(SZ), ALLOCATABLE :: PHI_MID(:, :, :)
   REAL(SZ), ALLOCATABLE :: PHI_INTEGRATED(:, :)
   REAL(SZ), ALLOCATABLE :: PSI_CHECK(:, :)
   REAL(SZ), ALLOCATABLE :: PSI1(:, :), PSI2(:, :), PSI3(:, :)
   REAL(SZ), ALLOCATABLE :: Q_HAT(:)
   REAL(SZ), ALLOCATABLE :: QIB(:)
   REAL(SZ), ALLOCATABLE :: CORI_EL(:), FRIC_EL(:)
   REAL(SZ), ALLOCATABLE :: ZE_MAX(:), ZE_MIN(:), DPE_MIN(:)
   REAL(SZ), ALLOCATABLE :: WATER_DEPTH_OLD(:, :), WATER_DEPTH(:, :)
   REAL(SZ), ALLOCATABLE :: ADVECTQX(:), ADVECTQY(:)
   REAL(SZ), ALLOCATABLE :: SOURCEQX(:), SOURCEQY(:)
   REAL(SZ), ALLOCATABLE :: QNAM_DG(:, :, :), QNPH_DG(:, :, :)
   REAL(SZ), ALLOCATABLE :: RHS_ZE(:, :, :)
   REAL(SZ), ALLOCATABLE :: XAGP(:, :), YAGP(:, :), WAGP(:, :)
   REAL(SZ), ALLOCATABLE :: XEGP(:, :),  WEGP(:, :)
   REAL(SZ), ALLOCATABLE :: SL3(:, :)
   REAL(SZ), ALLOCATABLE :: XBC(:), YBC(:)
   REAL(SZ), ALLOCATABLE :: XFAC(:, :, :, :), YFAC(:, :, :, :), &
                            SRFAC(:, :, :, :)
   REAL(SZ), ALLOCATABLE :: EDGEQ(:, :, :, :)
   REAL(SZ), ALLOCATABLE :: PHI(:), DPHIDZ1(:), DPHIDZ2(:)
   REAL(SZ), ALLOCATABLE :: PHI_STAE(:, :), PHI_STAV(:, :)

   !.....These (below) are defined in prep_slopelim.F

   Integer, Allocatable :: fact(:), focal_neigh(:, :), focal_up(:), bi(:), &
                           bj(:)
   Real(SZ), Allocatable :: XBCb(:), YBCb(:), xi1(:, :), xi2(:, :)
   Real(SZ), Allocatable :: xtransform(:, :), ytransform(:, :)
   Real(SZ), Allocatable :: xi1BCb(:), xi2BCb(:), xi1vert(:, :)
   Real(SZ), Allocatable :: xi2vert(:, :), xtransformv(:, :), &
                            ytransformv(:, :)
   Real(SZ), Allocatable :: Area_integral(:, :, :)
   Real(SZ), Allocatable :: f(:, :, :, :), g0(:, :, :, :), varsigma0(:, :, :, :)
   Real(SZ), Allocatable :: fv(:, :, :, :), g0v(:, :, :, :), &
                            varsigma0v(:, :, :, :)
   Real(SZ), Allocatable :: var2sigmag(:, :, :), var2sigmav(:, :, :)
   Real(SZ), Allocatable :: Nmatrix(:, :, :, :), NmatrixInv(:, :, :, :)
   Real(SZ), Allocatable :: deltx(:), delty(:), pmatrix(:, :, :)

   !.....These (below) are defined in slopelimiter.F

   Real(SZ), Allocatable :: ZEmin(:, :), ZEmax(:, :), QXmin(:, :)
   Real(SZ), Allocatable :: QXmax(:, :), QYmin(:, :), QYmax(:, :)

   ! namo - for ADCIRC -----------------------------------------------


   real(sz) :: G2ROOT

   REAL(SZ) :: etiminc_dg


   ! init in prep_DG
   INTEGER, ALLOCATABLE ::   pdg_el(:)

   ! initialized in prep_DG.F
   REAL(SZ), ALLOCATABLE ::    UMO(:, :), UFA(:, :), VMO(:, :), VFA(:, :)

   ! initialized in prep_DG.F
   INTEGER, ALLOCATABLE :: EL_COUNT(:)


   ! initialized in hstart.F
   REAL(SZ), ALLOCATABLE ::   DP0(:)

   ! initialized in write_results.F (temp variable?)
   real(sz), allocatable :: DPe(:)

   ! init in prep_DG.F
   real(sz) :: qtratio

   ! initialized in prep_DG.F
   INTEGER, ALLOCATABLE :: NNOEL(:, :)

   ! initialized in read_input.F
   ! same as ibtype_orig in ADCIRC
   !INTEGER,ALLOCATABLE ::    SEGTYPE(:)

   ! initialized in prep_DG.F
   REAL(SZ), ALLOCATABLE :: ANGTAB(:, :), CENTAB(:, :)
   integer, allocatable :: ELETAB(:, :)
   !INTEGER,ALLOCATABLE ::    NEITAB(:,:), neigh_elem(:,:)

   ! initialized in read_input.F
   ! this is probably the same as ndel?
   !integer, allocatable :: NEIGH_ELEM(:,:)

   ! init in prep_DG.F
   real(sz) :: h0l, h0h

   ! init in prep_DG.F
   real(sz), allocatable :: ydub(:, :, :)

   ! init in prep_DG.F
   real(sz) :: habsmin
   real(8) :: x1, x2, x3, y1, y2, y3

   ! init in read_input.F / prep_DG.F
   integer :: nstartdry

   ! init in create_edge_data.F
   INTEGER, ALLOCATABLE ::    EDFLG(:, :)

   ! init in create_edge_data.F
   INTEGER, ALLOCATABLE ::    NELED(:, :)

   ! init in read_input.F
   ! NOTE: nndel = nneighele and ndel = neitabele
   !$$$      integer, allocatable :: nndel(:) ! number of elements associated with a node
   !$$$      integer, allocatable :: ndel(:,:) ! NDEL(I,J) = number of Jth element associated with node I
   !$$$      integer MNNDEL            ! max number of elements associated with a node

   ! init in create_edge_data.F
   INTEGER, ALLOCATABLE ::    NOT_AN_EDGE(:), weir_buddy_node(:, :)
   integer :: jnmm
   INTEGER, ALLOCATABLE ::    ONE_OR_TWO(:)

   ! more variables not in prep_DG but somewhere else in dgswem
   ! init in adcirc.F
   REAL(SZ) :: NLEQ, LEQ, NLEQG

   real(sz), allocatable :: dg_ang(:), dp_dg(:)

   !**********************END OF DATA DECLARATIONS ***********************

CONTAINS

   !namo - allocate DG variables not in ADCIRC global.F

   subroutine alloc_adcirc()
      ALLOCATE (UMO(MNBFR, MNETA), UFA(MNBFR, MNETA))
      ALLOCATE (VMO(MNBFR, MNETA), VFA(MNBFR, MNETA))
      !ALLOCATE ( EL_COUNT(MNP) )
      ALLOCATE (DP0(MNP), DPe(MNE))
      !ALLOCATE ( NNOEL(MNP,mnei),CENTAB(MNP,mnei+1) )
      !ALLOCATE (ELETAB(MNP,mnei+1),ANGTAB(MNP,mnei+1))
      ALLOCATE (YDUB(36, MNE, 8))
      ALLOCATE (EDFLG(3, MNE))
      ALLOCATE (NELED(3, MNE))
      ALLOCATE (NOT_AN_EDGE(MNP))
      ALLOCATE (WEIR_BUDDY_NODE(MNP, 2))
      allocate (one_or_two(mnp))
      allocate (PDG_EL(MNE))

   end subroutine alloc_adcirc

   subroutine ALLOC_NNOEL1()
      allocate (el_count(mnp))
   end subroutine ALLOC_NNOEL1

   SUBROUTINE ALLOC_NNOEL2(maxel)
     integer, intent(in) :: maxel
      ALLOCATE (DP_DG(MAXEL), DG_ANG(MAXEL))
      ALLOCATE (NNOEL(MNP, MAXEL), CENTAB(MNP, MAXEL + 1))
      ALLOCATE (ELETAB(MNP, MAXEL + 1), ANGTAB(MNP, MAXEL + 1))
   end subroutine ALLOC_NNOEL2

   !.....Set edge array sizes


   SUBROUTINE ALLOC_EDGES1()
      ALLOCATE (NEDNO(2, MNED), NEDEL(2, MNED), NEDSD(2, MNED))
      ALLOCATE (NIBSEGN(2, MNED))
      ALLOCATE (NEBSEGN(MNED))
      ALLOCATE (NIEDN(MNED), NLEDN(MNED), NEEDN(MNED))
      ALLOCATE (NFEDN(MNED), NREDN(MNED), NIBEDN(MNED), NEBEDN(MNED))
      ALLOCATE (NCOUNT(MNED))
      ALLOCATE (COSNX(MNED), SINNX(MNED), XLEN(MNED))
      ALLOCATE (Q_HAT(MNED))
      RETURN
   end subroutine ALLOC_EDGES1

   !.....Set DG SWE array sizes

   SUBROUTINE ALLOC_DG1()
      ALLOCATE (EFA_DG(MNBFR, NEEDS + 2, 2), EMO_DG(MNBFR, NEEDS + 2, 2))
      ALLOCATE (UFA_DG(MNBFR, NEEDS + 2, 2), UMO_DG(MNBFR, NEEDS + 2, 2))
      ALLOCATE (VFA_DG(MNBFR, NEEDS + 2, 2), VMO_DG(MNBFR, NEEDS + 2, 2))
      RETURN
   end subroutine ALLOC_DG1

   SUBROUTINE ALLOC_DG2()
      ALLOCATE (QNAM_DG(MNFFR, NFEDS, 2), QNPH_DG(MNFFR, NFEDS, 2))
      RETURN
   end subroutine ALLOC_DG2

   SUBROUTINE ALLOC_DG3()
      ALLOCATE (QIB(MNP))
      RETURN
   end subroutine ALLOC_DG3

   SUBROUTINE ALLOC_DG4()
      !sb-20070228 NRK+1-->NRK+2 --- XX(:,:,NRK+2) will be used by slope limiter
      ALLOCATE (HB(DOFH, MNE, NRK + 2))
      ALLOCATE (MANN(DOFH, MNE)) !,arrayfix(DOFH,MNE,NRK+2) )
      ALLOCATE (U_modal(DOFH, MNE), V_modal(DOFH, MNE))
      ALLOCATE (ZE(DOFH, MNE, NRK + 2))
      !sb-20060711 For wet/dry
      ALLOCATE (ZE_MAX(MNE), ZE_MIN(MNE), DPE_MIN(MNE))
      ALLOCATE (WATER_DEPTH(MNE, 3), WATER_DEPTH_OLD(MNE, 3))
      ALLOCATE (ADVECTQX(MNE), ADVECTQY(MNE), &
                SOURCEQX(MNE), SOURCEQY(MNE))
      ALLOCATE (MARK(MNE))
      !em-2012 for sediment

      !--
      !sb-20070101
      !--
      ALLOCATE (RHS_ZE(DOFH, MNE, NRK))
      ALLOCATE (DRDX(MNE), DSDX(MNE), DRDY(MNE), DSDY(MNE))
      ALLOCATE (CORI_EL(MNE), FRIC_EL(MNE))
      ALLOCATE (PHI(DOFH), DPHIDZ1(DOFH), DPHIDZ2(DOFH))
      ALLOCATE (DOFS(MNE), PCOUNT(MNE))
      ALLOCATE (PDG(MNP))
      RETURN
   end subroutine ALLOC_DG4

   !.....Set sizes for arrays used in orthobasis


   !.....Set sizes for arrays for area integrals

   SUBROUTINE ALLOC_AREA_GAUSS()
      ALLOCATE (XAGP(NAGP(ph), ph), YAGP(NAGP(ph), ph), WAGP(NAGP(ph), ph))
      ALLOCATE (PHI_AREA(DOFH, NAGP(ph) + 1, ph))
      ALLOCATE (DSPHI(DOFH, NAGP(ph) + 1, ph), DRPHI(DOFH, NAGP(ph) + 1, ph))
      ALLOCATE (PHI_CORNER(DOFH, 3, ph), PHI_MID(DOFH, 3, ph))
      ALLOCATE (PHI_CENTER(DOFH, DOFH))
      ALLOCATE (PSI1(NAGP(ph), ph), PSI2(NAGP(ph), ph), PSI3(NAGP(ph), ph))
      ALLOCATE (BATH(NAGP(ph), MNE, ph), DBATHDX(NAGP(ph), MNE, ph))
      Allocate (DBATHDY(NAGP(ph), MNE, ph))
      ALLOCATE (SFAC_ELEM(NAGP(ph), MNE, ph))
      ALLOCATE (XFAC(DOFH, NAGP(ph), MNE, ph), YFAC(DOFH, NAGP(ph), MNE, ph))
      ALLOCATE (SRFAC(DOFH, NAGP(ph), MNE, ph))
      RETURN
   end subroutine ALLOC_AREA_GAUSS

   !.....Set sizes for arrays for edge integrals

   SUBROUTINE ALLOC_EDGE_GAUSS()
      ALLOCATE (XEGP(NEGP(ph), ph), WEGP(NEGP(ph), ph))
      ALLOCATE (PHI_EDGE(DOFH, NEGP(ph) + 1, 3, ph))
      ALLOCATE (M_INV(DOFH, ph))
      ALLOCATE (BATHED(NEGP(ph), 3, MNE, ph), SFACED(NEGP(ph), 3, MNE, ph))
      ALLOCATE (EDGEQ(DOFH, NEGP(ph), 3, ph))
      RETURN
   end subroutine ALLOC_EDGE_GAUSS

   !.....Set sizes for the arrays for the slope limiter
   !.....slopelim arrays

   SUBROUTINE ALLOC_SLOPELIM()
      ALLOCATE (XBC(MNE), YBC(MNE))
      ALLOCATE (EL_NBORS(4, MNE))
      ALLOCATE (SL3(3, MNE))

      !.....These are defined in prep_slopelim.F

      Allocate (fact(0:ph), focal_neigh(MNE, 3*MNEI), focal_up(MNE), &
                bi(dofh), bj(dofh))
      Allocate (XBCb(MNE), YBCb(MNE), xi1(MNE, NAGP(ph)), &
                xi2(MNE, NAGP(ph)))
      Allocate (xtransform(MNE, NAGP(ph)), ytransform(MNE, NAGP(ph)))
      Allocate (xi1BCb(MNE), xi2BCb(MNE), xi1vert(MNE, 3))
      Allocate (xi2vert(MNE, 3), xtransformv(MNE, 3), ytransformv(MNE, 3))
      allocate (Area_integral(MNE, 0:ph, 0:ph))
      Allocate (f(MNE, NAGP(ph), 0:ph, 0:ph), g0(MNE, NAGP(ph), 0:ph, 0:ph))
      Allocate (fv(MNE, 3, 0:ph, 0:ph), g0v(MNE, 3, 0:ph, 0:ph))
      Allocate (varsigma0(MNE, NAGP(ph), 0:ph, 0:ph))
      Allocate (varsigma0v(MNE, 3, 0:ph, 0:ph))
      Allocate (pmatrix(MNE, dofh, dofh), var2sigmag(MNE, NAGP(ph), dofh))
      Allocate (Nmatrix(MNE, dofh, dofh, dofh), &
                NmatrixInv(MNE, dofh, dofh, dofh))
      Allocate (deltx(MNE), delty(MNE), var2sigmav(MNE, 3, dofh))

      !.....These (below) are defined in slopelimiter.F (slopelimiter4)

      Allocate (ZEmin(MNP, dofh), ZEmax(MNP, dofh), QXmin(MNP, dofh))
      Allocate (QXmax(MNP, dofh), QYmin(MNP, dofh), QYmax(MNP, dofh))

   end subroutine ALLOC_SLOPELIM

   !sb...Set sizes for arrays for wetting and drying
   SUBROUTINE ALLOC_DG_WETDRY()
      ALLOCATE (WDFLG(MNE), DOFW(MNE), ncele(MNE))
      ALLOCATE (EL_UPDATED(MNE))
      ALLOCATE (WDFLG_TMP(MNE))
      ALLOCATE (LEDGE_NVEC(3, 3, MNE))
      ALLOCATE (DP_VOL(MNE, ph))
      ALLOCATE (PHI_INTEGRATED(DOFH, ph))
      ALLOCATE (PHI_CHECK(DOFH, NCHECK(ph), ph))
      ALLOCATE (DP_NODE(NCHECK(ph), MNE, ph))
      ALLOCATE (PSI_CHECK(3, 12*3))
      RETURN
   end subroutine ALLOC_DG_WETDRY

   SUBROUTINE ALLOC_STAE(L)
      integer, intent(in) :: L
      ALLOCATE (PHI_STAE(DOFH, L))
      RETURN
   end subroutine ALLOC_STAE

   SUBROUTINE ALLOC_STAV(L)
      integer, intent(in) :: L
      ALLOCATE (PHI_STAV(DOFH, L))
      RETURN
   end subroutine ALLOC_STAV

   SUBROUTINE PREP_DG()
    !! Allocate and initialize DG variables for solving the
    !! continuity equation

      use adc_constants, only: G, rad2deg
      use sizes, only: myproc, mnffr
      USE GLOBAL, only: ftiminc, eta2, efa, emo, noff, &
                        qnin1, qnam, qnph, qnin2, qtime1, h0, ifwind, nbfr, nffr, &
                        nstae, nstav, corif, xel, xev, yel, yev, nne, nnv,  peta1, peta2, &
                        IM, nolica, nolicat, nolifa, ihot, statim
      USE wetdry, only: computeWettingAndDrying
      USE NodalAttributes, ONLY: STARTDRY, FRIC, GeoidOffset, &
                                 LoadGeoidOffset, LoadManningsN, ManningsN
#ifdef CMPI
      use MESSENGER_ELEM, only: msg_table_elem, message_start_elem
#endif
      use BOUNDARIES, only: NOPE, NVDLL, nvell, nbou, nvel, ibtype_orig, lbcodei
      use mesh, only: NE, NM, neitab, nneigh, ics,  sfea0,  &
                     X, Y, areas, DP, SFAC

      IMPLICIT NONE

      ! dummy
      integer :: maxel

      !.....Declare local variables
      logical :: wetflag
      real(sz) :: col
      real(sz) :: qtratio_dg

      integer :: n1, n2, n3
      INTEGER :: II, l, P_0, DOF_0, j, k, kk, jj, i, chi,  Q, M, P, SZ2, w, III
      real(sz) :: ireal, jreal, kreal, jjreal
      CHARACTER(LEN=8) :: REGION
      REAL(SZ) :: AREA,  DP_MIN, XP, YP, timedg
      REAL(SZ) :: XI, YI, ZE1, ZE2, ZE3
      REAL(SZ), Allocatable :: BARY(:), VERT(:, :), BASIS(:), DBASIS(:, :)
      REAL(SZ), Allocatable :: PTS(:, :), WTS(:), PT(:)
      real(sz) :: R
      integer :: ELEM, ADDGP, NEDGS
      integer :: phh, DIM, NQEDS
      Real(SZ), allocatable :: XBCbt(:), YBCbt(:), radial(:), XB(:), YB(:), &
                               l2e(:)
      Real(SZ), allocatable :: hbo(:, :, :), &
                               ydubo(:, :)
      Real(SZ), allocatable :: YELEM(:), YED(:), HB1(:, :, :, :), zeo(:, :, :)

      !     ....................................................................

      !.....Allocate variables from DG not already in ADCIRC
      call alloc_adcirc()

      !.....Define variables from read fort.dg routine (will be added later)
      TIMEDG = statim
      DIM = 2
      RK_STAGE = 1
      RK_ORDER = 1
      NRK = RK_STAGE
      PADAPT = 0
      PL = 1
      PH = 1
      SLOPEFLAG = 6
      G2ROOT = SQRT(G/2.d0)

      !.....Set nonlinear flags
      if (nolica == 0 .or. nolicat == 0) then
         NLEQ = 0
         LEQ = 1
      else
         NLEQ = 1
         LEQ = 0
      end if
      NLEQG = NLEQ*G
      FG_L = LEQ*G

      IFWIND = 1
      IF (IM == 1) IFWIND = 0
      !     ....................................................................

      if (IM == 0) then
         DIM = 2
      elseif (IM == 1) then
         DIM = 3
      elseif (IM == 2) then
         DIM = 3
      end if
      Allocate (XBCbt(MNE), YBCbt(MNE), radial(MNE), XB(MNE), YB(MNE), &
                l2e(MNE))
      Allocate ( hbo(36, MNE, 1), &
                ydubo(36, mne))
      Allocate (YELEM(ph), YED(ph), hb1(36, mne, 1, ph), zeo(36, mne, 1), &
                BARY(DIM))
      Allocate (PT(DIM))

      C13 = 1.D0/3.D0
      C16 = 1.D0/6.D0
      R = 6378206.4d0


      !.....Obtain RK time scheme parameters

      CALL RK_TIME()

      !.....Compute the degrees of freedom per element

      DOF = (pl + 1)*(pl + 2)/2
      dofx = (px + 1)*(px + 2)/2 ! dofx for variable functions f=f(x)
      P_0 = pl
      DOF_0 = (pl + 1)*(pl + 2)/2 ! dof at lowest order when p!=0
      dofh = (ph + 1)*(ph + 2)/2

      !.....Allocate some DG stuff

      IF (PADAPT == 1) THEN

         dofh = (ph + 1)*(ph + 2)/2
         dofl = (pl + 1)*(pl + 2)/2
         pa = pl

      elseif (padapt == 0) then
         dofh = dofh
         dofl = DOF_0
         pa = pl

      end if


      !.....Compute the number of gauss points needed for the edge integrals

      CALL ALLOC_DG4() !moved here 6.28.10, for p_adapt because
      !of messenger_elem
      dofs(:) = dofl
      PDG_EL(:) = pl
      PDG(:) = pl
      PCOUNT(:) = 0
      pa = pl

      do chi = pl, ph
         NEGP(chi) = chi + 1
      end do

      IF (pl == 0) THEN

         PDG_EL(:) = 1
         PDG(:) = 1
         DOF = 3
         pl = 1
         dofl = 3
         P_0 = 0
         DOF_0 = 1
         NEGP(pl) = 2

      END IF

#ifdef CMPI
      CALL MSG_TABLE_ELEM() ! Read Message-Passing Tables
#endif

      !.....Create the edge based data

      IF (MYPROC == 0) THEN
         PRINT *, 'CREATING EDGE DATA...'
         PRINT *, ''
      END IF
      CALL CREATE_EDGE_DATA()
      IF (MYPROC == 0) THEN
         print *, 'CREATING EDGE DATA DONE'
         print *, ''
      END IF

#ifdef CMPI
      CALL MESSAGE_START_ELEM() ! Startup persistent message passing
#endif

      !.....Re-arrange elevation specified boundary segment data for DG

      IF (NEEDS > 0) THEN
         CALL ALLOC_DG1()
         II = 1
         JJ = 1
         DO I = 1, NBFR
            DO J = 1, NOPE
               DO K = 1, NVDLL(J) - 1
                  EMO_DG(I, II, 1) = EMO(I, JJ)
                  EMO_DG(I, II, 2) = EMO(I, JJ + 1)
                  EFA_DG(I, II, 1) = EFA(I, JJ)
                  EFA_DG(I, II, 2) = EFA(I, JJ + 1)
                  UMO_DG(I, II, 1) = UMO(I, JJ)
                  UMO_DG(I, II, 2) = UMO(I, JJ + 1)
                  UFA_DG(I, II, 1) = UFA(I, JJ)
                  UFA_DG(I, II, 2) = UFA(I, JJ + 1)
                  VMO_DG(I, II, 1) = VMO(I, JJ)
                  VMO_DG(I, II, 2) = VMO(I, JJ + 1)
                  VFA_DG(I, II, 1) = VFA(I, JJ)
                  VFA_DG(I, II, 2) = VFA(I, JJ + 1)
                  II = II + 1
                  JJ = JJ + 1
               END DO
               JJ = JJ + 1
            END DO
            II = 1
            JJ = 1
         END DO
      END IF

      !.....Re-arrange non-zero flow specified boundary segment data for DG

      IF (NFEDS > 0) THEN
         CALL ALLOC_DG2()
         II = 1
         JJ = 1
         DO I = 1, MNFFR
            IF (NFFR == 0) THEN
               QTRATIO_DG = (TIMEDG - QTIME1)/FTIMINC
               NQEDS = 0
               DO J = 1, NVEL
                  IF ((LBCODEI(J) == 2) .OR. (LBCODEI(J) == 12) &
                      .OR. (LBCODEI(J) == 22)) THEN
                     NQEDS = NQEDS + 1
                     IF (NQEDS <= NFEDS) THEN
                        QNAM_DG(1, NQEDS, 1) = (QNIN1(J) + &
                                                QTRATIO_DG*(QNIN2(J) - QNIN1(J)))
                        QNPH_DG(1, NQEDS, 1) = 0.D0
                        QNAM_DG(1, NQEDS, 2) = (QNIN1(J + 1) + &
                                                QTRATIO_DG*(QNIN2(J + 1) - QNIN1(J + 1)))
                        QNPH_DG(1, NQEDS, 2) = 0.D0
                     END IF
                  END IF
               END DO
            ELSE
               DO J = 1, NBOU
                  IF ((ibtype_orig(J) == 2) .OR. (ibtype_orig(J) == 12) &
                      .OR. (ibtype_orig(J) == 22)) THEN
                     DO K = 1, NVELL(J) - 1
                        QNAM_DG(I, II, 1) = QNAM(I, JJ)
                        QNAM_DG(I, II, 2) = QNAM(I, JJ + 1)
                        QNPH_DG(I, II, 1) = QNPH(I, JJ)
                        QNPH_DG(I, II, 2) = QNPH(I, JJ + 1)
                        II = II + 1
                        JJ = JJ + 1
                     END DO
                     JJ = JJ + 1
                  END IF
               END DO
               II = 1
               JJ = 1
            END IF
         END DO
      END IF

      !.....If there are internal barriers allocate some stuff

      IF (NIBEDS /= 0) CALL ALLOC_DG3()

      !.....Allocate the array for node to element table

      CALL ALLOC_NNOEL1()

      !.....Determine the number of elements connected at each node

      EL_COUNT = 0
      MAXEL = 1
      DO K = 1, 3
         DO J = 1, MNE
            N1 = NM(J, K)
            EL_COUNT(N1) = EL_COUNT(N1) + 1
         END DO
      END DO
      MAXEL = MAXVAL(EL_COUNT)

      !.....Allocate the array for the node to element table

      CALL ALLOC_NNOEL2(maxel)

      !.....Construct node to element table

      EL_COUNT = 0
      DO K = 1, 3
         DO J = 1, MNE
            N1 = NM(J, K)
            NNOEL(N1, 1 + EL_COUNT(N1)) = J
            EL_COUNT(N1) = EL_COUNT(N1) + 1
         END DO
      END DO

      !.....Construct node to element angle table

      DO I = 1, MNP
         !ETAMAX(I) = -99999
         KK = 1
         ELETAB(I, 1) = I
         S1 = SFAC(I)
         J1 = NEITAB(I, 1)
         DO K = 1, NNEIGH(I) - 1
            J2 = NEITAB(I, 1 + K)
            IF (K < (NNEIGH(I) - 1)) THEN
               J3 = NEITAB(I, 2 + K)
            ELSE
               J3 = NEITAB(I, 2)
            END IF
            DO J = 1, EL_COUNT(I)
               EL = NNOEL(I, J)
               N1 = NM(EL, 1)
               N2 = NM(EL, 2)
               N3 = NM(EL, 3)
               IF ((J1 == N1) .OR. (J1 == N2) .OR. (J1 == N3)) THEN
                  IF ((J2 == N1) .OR. (J2 == N2) .OR. (J2 == N3)) THEN
                     IF ((J3 == N1) .OR. (J3 == N2) .OR. (J3 == N3)) THEN
                        ELETAB(I, 1 + KK) = EL
                        S2 = SFAC(J2)
                        SAV = (S1 + S2)/2.D0
                        VEC1(1) = X(J1) - X(J2)
                        VEC1(2) = SAV*(Y(J1) - Y(J2))
                        VEC2(1) = X(J1) - X(J3)
                        VEC2(2) = SAV*(Y(J1) - Y(J3))
                        MAG1 = SQRT(VEC1(1)**2 + VEC1(2)**2)
                        MAG2 = SQRT(VEC2(1)**2 + VEC2(2)**2)
                        DOT = DOT_PRODUCT(VEC1, VEC2)
                        EL_ANG = ACOS(DOT/(MAG1*MAG2))
                        ANGTAB(I, KK + 1) = RAD2DEG*EL_ANG
                        KK = KK + 1
                        exit
                     END IF
                  END IF
               END IF
            END DO
         END DO
      END DO

      !     C.....Allocate some DG stuff

      !     CALL ALLOC_DG4()

      !.....Initialize the DG arrays

      ZE = 0.D0
      zeo = 0.D0
      hbo = 0.D0
      hb = 0.D0
      MARK = 0
      RHS_ZE = 0.D0
      !.....If using modal initial conditions transform the bathymetry from
      !.....nodal coordinates to modal dof

      ! namo - set eta to constant to test projection to modal
      ! ETA2 = 0.1
      ! hotstart = .false.

      DO J = 1, MNE
         N1 = NM(J, 1)
         N2 = NM(J, 2)
         N3 = NM(J, 3)
         hbo(1, J, 1) = 1.D0/3.D0*(DP(N1) + DP(N2) + DP(N3))
         hbo(2, J, 1) = -1.D0/6.D0*(DP(N1) + DP(N2)) + 1.D0/3.D0*DP(N3)
         hbo(3, J, 1) = -0.5D0*DP(N1) + 0.5D0*DP(N2)

         ydubo(1, J) = 1.D0/3.D0*(Y(N1) + Y(N2) + Y(N3))
         ydubo(2, J) = -1.D0/6.D0*(Y(N1) + Y(N2)) &
                       + 1.D0/3.D0*Y(N3)
         ydubo(3, J) = -0.5D0*Y(N1) + 0.5D0*Y(N2)

         ! namo - if hotstart from adcirc

         if (IHOT /= 0) then

            ZE(1, J, 1) = 1.D0/3.D0*(ETA2(N1) + ETA2(N2) + ETA2(N3))
            ZE(2, J, 1) = -1.D0/6.D0*(ETA2(N1) + ETA2(N2)) + 1.D0/3.D0*ETA2(N3)
            ZE(3, J, 1) = -0.5D0*ETA2(N1) + 0.5D0*ETA2(N2)

            do I = 1, 3
               col = ze(i, j, 1) + hbo(i, j, 1)
            end do
         end if
      END DO

      if (LoadManningsN) then
         DO J = 1, NE
            N1 = NM(J, 1)
            N2 = NM(J, 2)
            N3 = NM(J, 3)
            MANN(1, J) = 1.D0/3.D0*(ManningsN(N1) &
                                    + ManningsN(N2) + ManningsN(N3))
            MANN(2, J) = -1.D0/6.D0*(ManningsN(N1) &
                                     + ManningsN(N2)) + 1.D0/3.D0*ManningsN(N3)
            MANN(3, J) = -0.5D0*ManningsN(N1) + 0.5D0*ManningsN(N2)
         END DO
      end if

      IF (MODAL_IC == 0) THEN
         !     this assumes a cold start
         if (LoadGeoidOffset) then
            DO J = 1, NE
               N1 = NM(J, 1)
               N2 = NM(J, 2)
               N3 = NM(J, 3)
               zeo(1, J, 1) = 1.d0/3.d0*(GeoidOffset(N1) + GeoidOffset(N2) + &
                                         GeoidOffset(N3))
               IF (dof_0 /= 1) THEN
                  zeo(2, J, 1) = -1.d0/6.d0*(GeoidOffset(N1) + GeoidOffset(N2)) &
                                 + 1.d0/3.d0*GeoidOffset(N3)
                  zeo(3, J, 1) = -.5d0*GeoidOffset(N1) + .5d0*GeoidOffset(N2)
               END IF
            END DO
         end if
      END IF

      !--

      !As part of initializing the system, let us determine the partials
      !of the sediment discharge equation, fed in by fort.dg


      IF (MYPROC == 0) THEN
         print *, 'PREP FOR WET/DRY BEGINS...'
      END IF

      !.....1. Set initial surface elevation above the bed elevation
      !.....if wetting-and-drying is enabled and the initial water depth is
      !.....not specified by fort.12
      !.....2. Set wet-and-dry elemental flags
      !.....3. Set the DOF at dry elements = 1
      do j = 1, ph
         jj = 2*j
         NEGP(j) = CEILING(real((jj + 3),8)/2.0d0)
      end do

      NCHECK(1) = 3
      if (ph > 1) then
         do chi = 2, ph
            NCHECK(chi) = NCHECK(1) + 3*negp(chi)
         end do
      end if

      CALL ALLOC_DG_WETDRY()
      PHI_CHECK = 0.D0
      PSI_CHECK = 0.D0
      H0L = H0
      H0H = H0*1.0d0
      HABSMIN = H0*1.0d0

      !.....Retrieve the normals to the edges

      CALL CALC_NORMAL()
      !***********************************************************************
      ! dmw 03-14-2017 until next line of *'s
      !.....Master element type for quadrature
      !.....(1 = triangle, 2 = quad)
      ELEM = 1 ! For future release using quads. For now only tris work

      !.....Retrieve the area integral gauss quadrature points

      phh = -1 ! dummy initial value
      do j = 1, ph
         if (j == 1) then
            phh = 2*ph
            if (ELEM == 1 .and. DIM == 2) then
               !               Allocate( character(8) :: REGION )
               REGION = 'TRIANGLE'
               if (phh == 0) then
                  NAGP(ph) = 1
               elseif (phh <= 2) then
                  NAGP(ph) = 3
               elseif (phh <= 4) then
                  NAGP(ph) = 6
               elseif (phh <= 6) then
                  NAGP(ph) = 12
               elseif (phh <= 8) then
                  NAGP(ph) = 16
               elseif (phh <= 10) then
                  NAGP(ph) = 25
               elseif (phh <= 12) then
                  NAGP(ph) = 33
               elseif (phh <= 14) then
                  NAGP(ph) = 42
               elseif (phh <= 16) then
                  NAGP(ph) = 55
               elseif (phh <= 18) then
                  NAGP(ph) = 72
               elseif (phh <= 20) then
                  NAGP(ph) = 88
               end if
            elseif (ELEM == 2 .and. DIM == 2) then
               !               Allocate( character(6) :: REGION )
               REGION = 'SQUARE  '
               if (phh == 0) then
                  NAGP(ph) = 1
               elseif (phh <= 2) then
                  NAGP(ph) = 3
               elseif (phh <= 4) then
                  NAGP(ph) = 6
               elseif (phh <= 6) then
                  NAGP(ph) = 10
               elseif (phh <= 8) then
                  NAGP(ph) = 16
               elseif (phh <= 10) then
                  NAGP(ph) = 22
               elseif (phh <= 12) then
                  NAGP(ph) = 31
               elseif (phh <= 14) then
                  NAGP(ph) = 44
               elseif (phh <= 16) then
                  NAGP(ph) = 56
               elseif (phh <= 18) then
                  NAGP(ph) = 68
               elseif (phh <= 20) then
                  NAGP(ph) = 81
               elseif (phh <= 22) then
                  NAGP(ph) = 100
               end if
            end if
            !...........Allocate XaGP, YaGP, and WaGP

            call ALLOC_AREA_GAUSS()
         end if
         jj = 2*j
         sz2 = -1 ! dummy initial value
         if (ELEM == 1 .and. DIM == 2) then
            if (jj == 0) then
               SZ2 = 1
            elseif (jj <= 2) then
               SZ2 = 3
            elseif (jj <= 4) then
               SZ2 = 6
            elseif (jj <= 6) then
               SZ2 = 12
            elseif (jj <= 8) then
               SZ2 = 16
            elseif (jj <= 10) then
               SZ2 = 25
            elseif (jj <= 12) then
               SZ2 = 33
            elseif (jj <= 14) then
               SZ2 = 42
            elseif (jj <= 16) then
               SZ2 = 55
            elseif (jj <= 18) then
               SZ2 = 72
            elseif (jj <= 20) then
               SZ2 = 88
            end if
         elseif (ELEM == 2 .and. DIM == 2) then
            if (jj == 0) then
               SZ2 = 1
            elseif (jj <= 2) then
               SZ2 = 3
            elseif (jj <= 4) then
               SZ2 = 6
            elseif (jj <= 6) then
               SZ2 = 10
            elseif (jj <= 8) then
               SZ2 = 16
            elseif (jj <= 10) then
               SZ2 = 22
            elseif (jj <= 12) then
               SZ2 = 31
            elseif (jj <= 14) then
               SZ2 = 44
            elseif (jj <= 16) then
               SZ2 = 56
            elseif (jj <= 18) then
               SZ2 = 68
            elseif (jj <= 20) then
               SZ2 = 81
            elseif (jj <= 22) then
               SZ2 = 100
            end if
         end if

         Allocate (PTS(SZ2, DIM), WTS(SZ2))

         call QUADRATURE(jj, REGION, PTS, WTS, NAGP(j), SZ2, 2)
         XAGP(1:NAGP(j), j) = PTS(1:NAGP(j), 1)
         YAGP(1:NAGP(j), j) = PTS(1:NAGP(j), 2)

         WAGP(1:NAGP(j), j) = WTS(1:NAGP(j))
         Deallocate (PTS, WTS)
      end do
      !.....Retrieve the edge integral gauss quadrature points

      do j = 1, ph

         if (j == 1) then

            NEGP(ph) = CEILING(real((phh + 3),8)/2.0d0)
            CALL ALLOC_EDGE_GAUSS()

         end if
         jj = 2*j
         SZ2 = CEILING(real((jj + 3),8)/2.0d0)
         NEGP(j) = SZ2

         Allocate (PTS(SZ2, 1), WTS(SZ2))
         CALL quad_rules_general(NEGP(j), .true., 0.0d0, 0.0d0, PTS, WTS)

         XEGP(1:NEGP(j), j) = PTS(1:NEGP(j), 1)
         WEGP(1:NEGP(j), j) = WTS(1:NEGP(j))
         Deallocate (PTS, WTS)
      end do

      !.....Evaluate the orthogonal basis and its derivatives at the area
      !.....gauss quadrature points

      !.....Determine vertices and element barycenter by elements type

      if (ELEM == 1 .and. DIM == 2) then
         Allocate (VERT(3, DIM))
         BARY(:) = [-1.D0/3.D0, -1.D0/3.D0]
         VERT(1, :) = [-1.D0, -1.D0]
         VERT(2, :) = [1.D0, -1.D0]
         VERT(3, :) = [-1.D0, 1.D0]
         ADDGP = 3
      elseif (ELEM == 2 .and. DIM == 2) then
         Allocate (VERT(4, DIM))
         BARY(:) = [0.D0, 0.D0]
         VERT(1, :) = [-1.D0, -1.D0]
         VERT(2, :) = [1.D0, -1.D0]
         VERT(3, :) = [1.D0, 1.D0]
         VERT(4, :) = [-1.D0, 1.D0]
         ADDGP = 4
      end if

      !.....Loop over orders of polynomials required
      do w = 1, ph
         P = w
         L = w
         SZ2 = (P + 2)*(P + 1)**(DIM - 1)/2
         Allocate (BASIS(SZ2), DBASIS(SZ2, DIM))

         !.......Calculate orthogonal basis functions and derivatives at
         !.......gauss area integral points, element barycenter and vertices

         !....Loop over area gauss points, then element barycenter, then vertices
         do Q = 1, NAGP(L) + 1 + ADDGP
            if (Q <= NAGP(L)) then
               PT(1) = XAGP(Q, L)
               PT(2) = YAGP(Q, L)
            elseif (Q == NAGP(L) + 1) then
               PT(1) = BARY(1)
               PT(2) = BARY(2)
            else
               PT(1) = VERT(Q - NAGP(L) - 1, 1)
               PT(2) = VERT(Q - NAGP(L) - 1, 2)
            end if

            !...........Evaluate basis functions for specified element at current
            !...........point, order basis functions hierarchically
            call ORTHOGONAL_BASIS(ELEM, PT, P, DIM, BASIS, DBASIS)
            if (Q <= NAGP(L)) then
               PHI_AREA(1:SZ2, Q, P) = BASIS
               DRPHI(1:SZ2, Q, P) = DBASIS(:, 1)
               DSPHI(1:SZ2, Q, P) = DBASIS(:, 2)
            elseif (Q == NAGP(L) + 1) then
               PHI_CENTER(:, P) = BASIS
            else
               PHI_CORNER(:, Q - NAGP(L) - 1, P) = BASIS
            end if
         end do
         Deallocate (BASIS, DBASIS)
      end do
      Deallocate (BARY, VERT)

      !.....Evaluate the orthogonal basis at the edge gauss quadrature points

      if (ELEM == 1 .and. DIM == 2) then
         NEDGS = 3
      elseif (ELEM == 2 .and. DIM == 2) then
         NEDGS = 4
      elseif (ELEM == 1 .and. DIM == 3) then
         NEDGS = 9
      elseif (ELEM == 2 .and. DIM == 3) then
         NEDGS = 12
      else
         print *, '  ****** ERROR!!! INVALID CHOICE OF ELEMENT/DIMENSION &
  &            PARAMETER!  *******  '
         print *, '  Execution terminated in subroutine prep_DG  '
         STOP
      end if

      do w = 1, ph
         P = w
         L = w
         SZ2 = (P + 2)*(P + 1)**(DIM - 1)/2
         Allocate (BASIS(SZ2), DBASIS(SZ2, DIM))

         !........Calculate entries for inverse mass stiffness matrix

         if (ELEM == 1 .and. DIM == 2) then
            !........For triangular elements
            M = 1
            do J = 0, P
               do I = 0, J
                  ireal = real(i,8)
                  JJ = J - I
                  jjreal = real(jj,8)
                  M_INV(M, P) = ((2.D0*Ireal + 1.D0)*(2.D0*JJreal + 2.D0*Ireal + 2.D0)/4.D0)
                  M = M + 1
               end do
            end do
         elseif (ELEM == 2 .and. DIM == 2) then
            !........For rectangular elements
            M = 1
            do J = 0, P
               do I = 0, J
                  ireal = real(i,8)
                  JJ = J - I
                  jjreal = real(jj, 8)
                  M_INV(M, P) = ((2.D0*Ireal + 1.D0)*(2.D0*jjreal + 1.D0)/4.D0)
                  M = M + 1
               end do
            end do
         elseif (ELEM == 1 .and. DIM == 3) then
            !........For triangular prism elements
            M = 1
            do k = 0, P
               kreal = real(k,8)
               do J = 0, P
                  jreal = real(j,8)
                  do I = 0, J
                     ireal = real(i,8)
                     JJ = J - I
                     jjreal = real(jj, 8)
                     M_INV(M, P) = ((2.D0*Ireal + 1.D0)*(2.D0*JJreal + 2.D0*Ireal + 2.D0)* &
                                    (2.D0*kreal + 1.D0)/8.D0)
                     M = M + 1
                  end do
               end do
            end do
         elseif (ELEM == 2 .and. DIM == 3) then
            !........For rectangular hexahedron elements
            M = 1
            do k = 0, P
               kreal = real(k,8)
               do J = 0, P
                  jreal = real(j,8)
                  do I = 0, J
                     ireal = real(i,8)
                     JJ = J - I
                     jjreal = real(jj,8)
                     M_INV(M, P) = ((2.D0*Ireal + 1.D0)*(2.D0*JJreal + 1.D0)* &
                                    (2.D0*kreal + 1.D0)/8.D0)
                     M = M + 1
                  end do
               end do
            end do
         end if

         !........Loop over edges
         do II = 1, NEDGS
            !..........Loop over edge gauss points
            do Q = 1, NEGP(L) + 1
               !..............Get edge gauss points (and midpoint) for triangle edges
               if (ELEM == 1 .and. DIM == 2) then
                  if (II == 1) then
                     if (Q < (NEGP(L) + 1)) then
                        PT(1) = -XEGP(Q, L)
                        PT(2) = XEGP(Q, L)
                     else
                        PT(1) = 0.D0
                        PT(2) = 0.D0
                     end if
                  elseif (II == 2) then
                     if (Q < (NEGP(L) + 1)) then
                        PT(1) = -1.D0
                        PT(2) = -XEGP(Q, L)
                     else
                        PT(1) = -1.D0
                        PT(2) = 0.D0
                     end if
                  elseif (II == 3) then
                     if (Q < (NEGP(L) + 1)) then
                        PT(1) = XEGP(Q, L)
                        PT(2) = -1.D0
                     else
                        PT(1) = 0.D0
                        PT(2) = -1.D0
                     end if
                  end if
                  !..............Get edge gauss points (and midpoint) for rectangle edges
               elseif (ELEM == 2 .and. DIM == 2) then
                  if (II == 1) then
                     if (Q < (NEGP(L) + 1)) then
                        PT(1) = -XEGP(Q, L)
                        PT(2) = 1.D0
                     else
                        PT(1) = 0.D0
                        PT(2) = 1.D0
                     end if
                  elseif (II == 2) then
                     if (Q < (NEGP(L) + 1)) then
                        PT(1) = -1.D0
                        PT(2) = -XEGP(Q, L)
                     else
                        PT(1) = -1.D0
                        PT(2) = 0.D0
                     end if
                  elseif (II == 3) then
                     if (Q < (NEGP(L) + 1)) then
                        PT(1) = XEGP(Q, L)
                        PT(2) = -1.D0
                     else
                        PT(1) = 0.D0
                        PT(2) = -1.D0
                     end if
                  elseif (II == 4) then
                     if (Q < (NEGP(L) + 1)) then
                        PT(1) = 1.D0
                        PT(2) = XEGP(Q, L)
                     else
                        PT(1) = 1.D0
                        PT(2) = 0.D0
                     end if
                  end if
               end if

               !..............Evaluate basis functions for specified element at current
               !..............point, order basis functions hierarchically
               call ORTHOGONAL_BASIS(ELEM, PT, P, DIM, BASIS, DBASIS)
               if (Q <= NEGP(L)) then
                  PHI_EDGE(:, Q, II, P) = BASIS
               else
                  PHI_MID(:, II, P) = BASIS
               end if
            end do

            if (ELEM == 1 .and. DIM == 2) then
               XP = PT(1)
               YP = PT(2)
               III = 3
               do Q = 1, NEGP(L) + 1
                  if (P > 1) then
                     if (Q < (NEGP(L) + 1)) then
                        PSI_CHECK(1, III) = -1.D0/2.D0*(XP + YP)
                        PSI_CHECK(2, III) = 1.D0/2.D0*(XP + 1.D0)
                        PSI_CHECK(3, III) = 1.D0/2.D0*(YP + 1.D0)
                        III = III + 1
                     end if
                  end if
               end do
            elseif (ELEM == 2 .and. DIM == 2) then
               XP = PT(1)
               YP = PT(2)
               III = 4
               do Q = 1, NEGP(L) + 1
                  if (P > 1) then
                     if (Q < (NEGP(L) + 1)) then
                        PSI_CHECK(1, III) = 1.D0/2.D0*(XP - 1.D0)
                        PSI_CHECK(2, III) = 1.D0/2.D0*(YP - 1.D0)
                        PSI_CHECK(3, III) = 1.D0/2.D0*(XP + 1.D0)
                        PSI_CHECK(4, III) = 1.D0/2.D0*(YP + 1.D0)
                        III = III + 1
                     end if
                  end if
               end do
            end if
         end do
         Deallocate (BASIS, DBASIS)
      end do
      !***********************************************************************

      !.....Do the L2-projection of the initial conditions

      !hb = 0.D0
      !qx = 0.D0
      !qy = 0.D0
      !ze = 0.D0
      !ydub = 0.d0
      !hb1 = 0.D0
      !
      hb(1:dofh, :, 1) = hbo(1:dofh, :, 1)


      if (IHOT == 0) then
         ze(1:dofh, :, 1) = zeo(1:dofh, :, 1)
      end if

      ydub(1:dofh, :, 1) = ydubo(1:dofh, :)
      do chi = pl, ph
         hb1(1:dofh, :, 1, chi) = hbo(1:dofh, :, 1)
      end do
      !....if layers are on, distribute them evenly across the total bed load
      !....other approaches are clearly available, this is a simple choice
      !....adapt for higher order initial data


      DO J = 1, NE
         DPE_MIN(J) = MIN(DP(NM(J, 1)), DP(NM(J, 2)), DP(NM(J, 3)))
      END DO

      !.....Compute the values of the nodal basis functions at the
      !.....area gauss quadrature points, at every p level chi

      do chi = 1, ph
         do I = 1, NAGP(chi)
            PSI1(I, chi) = -1.D0/2.D0*(XAGP(I, chi) + YAGP(I, chi))
            PSI2(I, chi) = 1.D0/2.D0*(XAGP(I, chi) + 1.D0)
            PSI3(I, chi) = 1.D0/2.D0*(YAGP(I, chi) + 1.D0)
         end do
      end do

      !.....Store the derivatives of the (linear) nodal basis functions

      DRPSI(1) = -1.D0/2.D0
      DRPSI(2) = 1.D0/2.D0
      DRPSI(3) = 0.D0
      DSPSI(1) = -1.D0/2.D0
      DSPSI(2) = 0.D0
      DSPSI(3) = 1.D0/2.D0

      !.....Pre-compute the derivatives of the coordinate transformation for
      !.....each element

      DO J = 1, MNE

         !.....Retrieve the global node numbers for the element

         N1 = NM(J, 1)
         N2 = NM(J, 2)
         N3 = NM(J, 3)
         x1 = x(n1)
         y1 = y(n1)
         x2 = x(n2)
         y2 = y(n2)
         x3 = x(n3)
         y3 = y(n3)
         AREA = (X1 - X3)*(Y2 - Y3) + (X3 - X2)*(Y1 - Y3)
         area = area*.5d0

         !.....Compute the derivatives of the coordinate transformation

         DRDX(J) = 1.D0/AREA*(Y(N3) - Y(N1))
         DSDX(J) = 1.D0/AREA*(Y(N1) - Y(N2))

         DRDY(J) = 1.D0/AREA*(X(N1) - X(N3))
         DSDY(J) = 1.D0/AREA*(X(N2) - X(N1))

         !.......Compute elemental Coriolis and friction terms

         CORI_EL(J) = (CORIF(N1) + CORIF(N2) + CORIF(N3))/3.D0
         FRIC_EL(J) = (FRIC(N1) + FRIC(N2) + FRIC(N3))/3.D0

         !.......Pre-compute the bathymetry and the gradient of the bathymetry at
         !.......the quadrature points and compute volume of water

         DP_VOL(J, :) = 0.D0
         SFAC_ELEM(:, J, :) = 0.D0

         chi = 1
         area_quad: do I = 1, NAGP(chi) ! Area quadrature points

            BATH(I, J, chi) = 0.D0
            DBATHDX(I, J, chi) = 0.D0
            DBATHDY(I, J, chi) = 0.D0
            YELEM(chi) = 0.d0

            DO K = 1, (chi + 1)*(chi + 2)/2

               DPHIDX = DRPHI(K, I, chi)*DRDX(J) + DSPHI(K, I, chi)*DSDX(J)
               DPHIDY = DRPHI(K, I, chi)*DRDY(J) + DSPHI(K, I, chi)*DSDY(J)
               XFAC(K, I, J, chi) = M_INV(K, chi)*WAGP(I, chi)*DPHIDX
               YFAC(K, I, J, chi) = M_INV(K, chi)*WAGP(I, chi)*DPHIDY
               SRFAC(K, I, J, chi) = M_INV(K, chi)*WAGP(I, chi)* &
                                     PHI_AREA(K, I, chi)
               BATH(I, J, chi) = BATH(I, J, chi) + HB(K, J, 1)* &
                                 PHI_AREA(K, I, chi)
               YELEM(chi) = YELEM(chi) + YDUB(K, J, chi)* &
                            PHI_AREA(K, I, chi)

               IF (ICS == 1) THEN
                  SFAC_ELEM(I, J, chi) = 1.0D0
               ELSE
                  SFAC_ELEM(I, J, chi) = COS(SFEA0)/COS(YELEM(chi)/R)
               END IF

               DBATHDX(I, J, chi) = DBATHDX(I, J, chi) + HB(K, J, 1)*DPHIDX
               DBATHDY(I, J, chi) = DBATHDY(I, J, chi) + HB(K, J, 1)*DPHIDY
               DP_VOL(J, chi) = DP_VOL(J, chi) + WAGP(I, chi)*HB(K, J, 1)* &
                                PHI_AREA(K, I, chi)

            END DO

         END DO area_quad

         do chi = 1, ph
            DP_VOL(J, chi) = 0.25D0*AREAS(J)*DP_VOL(J, chi)
         end do

         do chi = 1, ph

            if (chi >= 1) then

               DO L = 1, 3
                  do I = 1, NEGP(chi) ! Edge quadrature points

                     BATHED(I, L, J, chi) = 0.D0
                     yed(chi) = 0.d0

                     DO K = 1, (chi + 1)*(chi + 2)/2

                        BATHED(I, L, J, chi) = BATHED(I, L, J, chi) + HB(K, J, 1)* &
                                               PHI_EDGE(K, I, L, chi)

                        YED(chi) = YED(chi) + YDUB(K, J, chi)* &
                                   PHI_EDGE(K, I, L, chi)

                        IF (ICS == 1) THEN
                           SFACED(I, L, J, chi) = 1.0d0
                        ELSE
                           SFACED(I, L, J, chi) = COS(SFEA0)/COS(YED(chi)/R)
                        END IF

                        EDGEQ(K, I, L, chi) = 2.0d0*M_INV(K, chi)* &
                                              PHI_EDGE(K, I, L, chi)*WEGP(I, chi)

                     END DO
                  end do
               END DO

            else

               DO L = 1, 3
                  do I = 1, NEGP(chi) ! Edge quadrature points

                     BATHED(I, L, J, chi) = 0.D0
                     yed(chi) = 0.d0

                     DO K = 1, (chi + 1)*(chi + 2)/2

                        BATHED(I, L, J, chi) = BATHED(I, L, J, chi) + HB(K, J, 1)* &
                                               PHI_EDGE(K, I, L, chi)

                        YED(chi) = YED(chi) + YDUB(K, J, chi)* &
                                   PHI_EDGE(K, I, L, chi)

                        IF (ICS == 1) THEN
                           SFACED(I, L, J, chi) = 1.0d0
                        ELSE
                           SFACED(I, L, J, chi) = COS(SFEA0)/COS(YED(chi)/R)
                        END IF

                        EDGEQ(K, I, L, chi) = 2.0d0*M_INV(K, chi)* &
                                              PHI_EDGE(K, I, L, chi)*WEGP(I, chi)

                     END DO
                  end do
               END DO

            end if

         end do

         !........Store bathymetry at triangular vertices and edge gauss points
         !........for wet-dry

         do chi = 1, ph

            DO I = 1, 3
               DP_NODE(I, J, chi) = DP(NM(J, I))
            END DO

            IF (NCHECK(chi) > 3) THEN
               II = 4
               DO L = 1, 3
                  DO I = 1, NEGP(chi)
                     DP_NODE(II, J, chi) = BATHED(I, L, J, chi)
                     II = II + 1
                  END DO
               END DO
            END IF
         end do
      END DO

      DO I = 1, 3
         IF (I == 1) THEN
            XI = -1.D0
            YI = -1.D0
         ELSEIF (I == 2) THEN
            XI = 1.D0
            YI = -1.D0
         ELSE
            XI = -1.D0
            YI = 1.D0
         END IF
         PSI_CHECK(1, I) = -1.D0/2.D0*(XI + YI)
         PSI_CHECK(2, I) = 1.D0/2.D0*(XI + 1.D0)
         PSI_CHECK(3, I) = 1.D0/2.D0*(YI + 1.D0)
         do chi = 1, ph
            DO K = 1, (chi + 1)*(chi + 2)/2
               PHI_CHECK(K, I, chi) = PHI_CORNER(K, I, chi)
            END DO
         end do
      END DO

      do chi = 1, ph
         IF (NCHECK(chi) > 3) THEN
            II = 4
            DO L = 1, 3
               DO I = 1, NEGP(chi)
                  DO K = 1, (chi + 1)*(chi + 2)/2

                     PHI_CHECK(K, II, chi) = PHI_EDGE(K, I, L, chi)

                  END DO
                  II = II + 1
               END DO
            END DO
         END IF
      end do
      !.....Integrate the basis functions

      PHI_INTEGRATED = 0.D0
      do chi = 1, ph
         DO I = 1, NAGP(chi)
            DO K = 1, (chi + 1)*(chi + 2)/2
               PHI_INTEGRATED(K, chi) = PHI_INTEGRATED(K, chi) + WAGP(I, chi)* &
                                        PHI_AREA(K, I, chi)
            end do
         END DO
      END DO
      !.....Wetting and drying is not turned on

      ! namo - disable this part on wet/dry

      wetflag = .true.

      NSTARTDRY = 0

      if (wetflag) then
         IF (NOLIFA == 0 .OR. NOLIFA == 1) THEN
            DO J = 1, MNE
               WDFLG(J) = 1
               !DOFS(J) = 3
            END DO

            !.....Wetting and drying is turned on but there are no dry nodes below
            !.....geoid

         ELSEIF (NOLIFA == 2 .AND. NSTARTDRY == 0) THEN

            DO J = 1, MNE

               !.........Check to see if there are initially any dry nodes

               !ZE1 = 0.D0
               !ZE2 = 0.D0
               !ZE3 = 0.D0

               ZE1 = ze(1, J, 1)
               ZE2 = ze(1, J, 1)
               ZE3 = ze(1, J, 1)
               IF (DP(NM(J, 1)) < H0) ZE1 = max(ze1, H0 - DP(NM(J, 1)))
               IF (DP(NM(J, 2)) < H0) ZE2 = max(ze2, H0 - DP(NM(J, 2)))
               IF (DP(NM(J, 3)) < H0) ZE3 = max(ze3, H0 - DP(NM(J, 3)))

               !.........If so set initial surface elevation values

               IF (abs((ZE1 + ZE2 + ZE3)/3.D0 - ze(1, J, 1)) >= 1d-12) THEN
                  IF (p_0 == 0) THEN
                     DP_MIN = MIN(DP(NM(J, 1)), DP(NM(J, 2)), DP(NM(J, 3)))
                     ze(1, J, 1) = max(ze(1, j, 1), H0 - DP_MIN)
                  ELSE
                     IF (ze(1, J, 1) > (ZE1 + ZE2 + ZE3)/3.d0) THEN
                        IF (dof_0 /= 1) THEN
                           ze(2, J, 1) = 0.d0
                           ze(3, J, 1) = 0.d0
                           ze(4:dofh, J, 1) = 0.D0 ! forced again for
                           ! transparency
                        END IF
                     ELSE
                        ze(1, J, 1) = (ZE1 + ZE2 + ZE3)/3.D0
                        IF (DOF_0 /= 1) THEN
                           ze(2, J, 1) = -1.D0/6.D0*(ZE1 + ZE2) + 1.D0/3.D0*ZE3
                           ze(3, J, 1) = -0.5D0*ZE1 + 0.5D0*ZE2
                           ze(4:dofh, J, 1) = 0.D0 ! forced again for
                           ! transparency
                        END IF
                     END IF
                  END IF

                  !............and set wet/dry flag (0 = dry, 1 = wet)

                  WDFLG(J) = 0
                  DOFS(J) = dof_0
               ELSE
                  WDFLG(J) = 1
                  DOFS(J) = dofs(J)
               END IF

            END DO

            !.....If there are dry nodes below geoid

         ELSEIF (NOLIFA == 2 .AND. NSTARTDRY == 1) THEN

            !.......Loop over elements

            DO J = 1, NE

               !........Retrieve global node numbers for element

               N1 = NM(J, 1)
               N2 = NM(J, 2)
               N3 = NM(J, 3)

               !........Check to see if nodes are initially dry

               ZE1 = 0.d0
               ZE2 = 0.d0
               ZE3 = 0.d0
               IF (int(STARTDRY(N1)) == 1) ZE1 = H0 - DP(N1)
               IF (DP(N1) < H0) ZE1 = H0 - DP(N1)
               IF (int(STARTDRY(N2)) == 1) ZE2 = H0 - DP(N2)
               IF (DP(N2) < H0) ZE2 = H0 - DP(N2)
               IF (int(STARTDRY(N3)) == 1) ZE3 = H0 - DP(N3)
               IF (DP(N3) < H0) ZE3 = H0 - DP(N3)

               IF (MODAL_IC == 3) THEN
                  IF (int(STARTDRY(N1)) == -88888) then
                     ZE1 = H0 - DP(N1)
                  else
                     ZE1 = STARTDRY(N1)
                  end if
                  IF (int(STARTDRY(N2)) == -88888) then
                     ZE2 = H0 - DP(N2)
                  else
                     ZE2 = STARTDRY(N2)
                  end if
                  IF (int(STARTDRY(N3)) == -88888) then
                     ZE3 = H0 - DP(N3)
                  else
                     ZE3 = STARTDRY(N3)
                  end if
               END IF

               !.........If so set initial surface elevation values

               IF (abs(ZE1 + ZE2 + ZE3) >= 0.d0) THEN
                  IF (P_0 == 0) THEN
                     DP_MIN = MIN(DP(NM(J, 1)), DP(NM(J, 2)), DP(NM(J, 3)))
                     ze(1, J, 1) = max(ze(1, j, 1), H0 - DP_MIN)
                  ELSE
                     IF (ze(1, J, 1) > (ZE1 + ZE2 + ZE3)/3.d0) THEN
                        IF (DOF_0 /= 1) THEN
                           ze(2, J, 1) = 0.d0
                           ze(3, J, 1) = 0.d0
                           ze(4:dofh, J, 1) = 0.D0 ! forced again for
                           ! transparency
                        END IF
                     ELSE
                        ze(1, J, 1) = (ZE1 + ZE2 + ZE3)/3.D0
                        IF (DOF_0 /= 1) THEN
                           ze(2, J, 1) = -1.D0/6.D0*(ZE1 + ZE2) + 1.D0/3.D0*ZE3
                           ze(3, J, 1) = -0.5D0*ZE1 + 0.5D0*ZE2
                           ze(4:dofh, J, 1) = 0.D0 ! forced again for
                           ! transparency
                        END IF
                     END IF
                  END IF

                  !...........Set wet/dry flag (0 = dry, 1 = wet)

                  WDFLG(J) = 0
                  DOFS(J) = DOF_0
               ELSE
                  WDFLG(J) = 1
                  DOFS(J) = DOFS(J)
               END IF
            END DO
         END IF

         IF (MYPROC == 0) THEN
            print *, 'DONE'
            print *, ''
         END IF

      end if


      !.....Set p back to original value if p = 0

      IF (P_0 /= pl) THEN
         PDG_EL(:) = 0
         PDG(:) = 0
         DOF = 1
         DOFL = 1
         DOFS(:) = 1
         DOF_0 = 1
         pl = 0
         !p = 0
         pa = 0
      END IF

      !.....Compute basis functions at stations

      IF (NSTAE > 0) THEN ! Elevation stations
         CALL ALLOC_STAE(NSTAE)
         DO I = 1, NSTAE
            CALL STA_BASIS(ELEM, DIM, XEL(I), YEL(I), NNE(I), &
                           PHI_STAE(:, I))
         END DO
      END IF

      IF (NSTAV > 0) THEN ! Velocity Stations
         CALL ALLOC_STAV(NSTAV)
         DO I = 1, NSTAV
            CALL STA_BASIS(ELEM, DIM, XEV(I), YEV(I), NNV(I), &
                           PHI_STAV(:, I))
         END DO
      END IF

      !.....Prep the slopelimiter

      IF (SLOPEFLAG /= 0) THEN
         IF (MYPROC == 0) THEN
            print *, 'Slope limiting prep begins, "kshanti"'
            print *, 'Using slope limiter ', SLOPEFLAG
            print *, 'H0 = ', H0
         END IF
         CALL ALLOC_SLOPELIM()
         CALL PREP_SLOPELIM()
         IF (MYPROC == 0) THEN
            print *, 'Finished'
         END IF
      END IF
      !--
      ! 2/1/24 - Hot starting Only works when everything is wet
      ! ETA2, UU2, and VV2 have been read from a hotstart file
      ! after the call to HOTSTART() in adcirc.F
      IF (IHOT /= 0) THEN
         print *, 'Hotstarting mode'
         DO J = 1, NE
            N1 = NM(J, 1)
            N2 = NM(J, 2)
            N3 = NM(J, 3)
            ZE(1, J, 1) = 1.D0/3.D0*(ETA2(N1) + ETA2(N2) + ETA2(N3))
            ZE(2, J, 1) = -1.D0/6.D0*(ETA2(N1) + ETA2(N2)) + 1.D0/3.D0*ETA2(N3)
            ZE(3, J, 1) = -0.5D0*ETA2(N1) + 0.5D0*ETA2(N2)
         END DO
      END IF

      NOFF = WDFLG
      peta2 = 0.D0
      peta1 = 0.D0

      if (LoadGeoidOffset) then
         DO J = 1, NE
            if (WDFLG(j) == 0) then
               N1 = NM(J, 1)
               N2 = NM(J, 2)
               N3 = NM(J, 3)
               ze(1, J, 1) = ze(1, J, 1) + 1.d0/3.d0*(GeoidOffset(N1) + GeoidOffset(N2) + &
                                                      GeoidOffset(N3))
               ze(2, J, 1) = ze(2, J, 1) + (-1.d0/6.d0*(GeoidOffset(N1) + GeoidOffset(N2)) &
                                            + 1.d0/3.d0*GeoidOffset(N3))
               ze(3, J, 1) = ze(3, J, 1) + (-.5d0*GeoidOffset(N1) + .5d0*GeoidOffset(N2))
            end if
         END DO
      end if
   END SUBROUTINE PREP_DG

   SUBROUTINE CALC_NORMAL()

      !.....Use appropriate modules

      use mesh, only: X, Y, sfac

      IMPLICIT NONE

      !.....Declare local variables

      INTEGER :: i, n1, n2

      !.....Loop over the edges

      DO I = 1, NEDGES

         !.....Retrieve the node numbers for the given edge

         N1 = NEDNO(1, I)
         N2 = NEDNO(2, I)

         !.....Compute an average SFAC to adjust normal for CPP coordinates

         SAV = (SFAC(N1) + SFAC(N2))/2.d0

         !.....Compute the length of the given egde

         XLEN(I) = SQRT((Y(N2) - Y(N1))**2.D0 &
                        + (X(N2) - X(N1))**2.D0)

         !.....Compute the components of the normal vector

         !COSNX(I) = SAV*(Y(N2) - Y(N1))/XLEN(I)
         COSNX(I) = (Y(N2) - Y(N1))/XLEN(I)
         SINNX(I) = -(X(N2) - X(N1))/XLEN(I)
      END DO

      RETURN
   end subroutine CALC_NORMAL

   !
   !     SUBROUTINE  CREATE_EDGE_DATA.F
   !
   !     This program takes the original ADCIRC data sets given in files
   !     'fort.14' and 'fort.15' and generates edge based data structures
   !
   !     Written by: Srinivas Chippada
   !
   !     Turned into a subroutine by Clint Dawson, May 2002
   !
   !     Modifications for different boundary types and notational changes
   !     made by Ethan Kubatko, March 2005
   !
   !     Modifications for parallel runs
   !     made by Shintaro Bunya, Aug 2005
   !
   !     Bug fix
   !     made by Shintaro Bunya, Aug 26, 2005
   !
   !     Modifications to skip an edge on a barriar overlapping on
   !     the external boundary.  Boundary edge creation in a group
   !     used to be aborted when such an edge is found.
   !     made by Shintaro Bunya, Sep  1, 2005
   !
   !     Speed up in generating edge pairs
   !     made by Shintaro Bunya, Feb 26, 2006
   !
   !***********************************************************************

   SUBROUTINE CREATE_EDGE_DATA()

      !.....Use appropriate modules

      !USE GLOBAL
      use sizes, only: myproc
      use mesh, only: NM, NE, neitabele, nneighele
      use boundaries, only: nvel, ibconn, nvell, nbou, nope, nvdll, nbdv, &
                            nbv, nbvv, lbcodei, ibtype_orig
      IMPLICIT NONE

      !.....Declare local variables

      INTEGER :: IED, JED, JJED, IEL, IEL1, IEL2, JJ1, JJ2
      !sb--
      INTEGER :: JEL, JJEL, ED_ID,   I, J, K
      !--
      INTEGER :: I1, I2,  II, LED(2, 3), n1, n2, n3
      INTEGER :: NERR, NN1, NN2
      INTEGER :: IERROR

      ! namo - for adcirc
      !nndel = nodele

      !.....Compute maximum number of edges

      MNED = 3*MNE

      !.....Allocate the edge data arrays

      CALL ALLOC_EDGES1()

      !.....Generate the edge connectivity

      DO J = 1, MNE
         EDFLG(1, J) = 0
         EDFLG(2, J) = 0
         EDFLG(3, J) = 0
         NELED(1, J) = 0
         NELED(2, J) = 0
         NELED(3, J) = 0
      END DO

      NEDNO(:, :) = 0
      NEDEL(:, :) = 0

      NEDGES = 0

      IF (MYPROC == 0) THEN
         PRINT *, 'CREATING EDGE PAIRS...'
      END IF

      DO IEL = 1, NE
         N1 = NM(IEL, 1)
         N2 = NM(IEL, 2)
         N3 = NM(IEL, 3)
         LED(1, 1) = N2
         LED(2, 1) = N3
         LED(1, 2) = N3
         LED(2, 2) = N1
         LED(1, 3) = N1
         LED(2, 3) = N2

         DO IED = 1, 3
            IF (EDFLG(IED, IEL) == 1) cycle

            I1 = LED(1, IED)
            I2 = LED(2, IED)

            NEDGES = NEDGES + 1
            ED_ID = NEDGES
            NELED(IED, IEL) = ED_ID
            NEDNO(1, ED_ID) = I1
            NEDNO(2, ED_ID) = I2
            NEDEL(1, ED_ID) = IEL
            EDFLG(IED, IEL) = 1

            DO JJEL = 1, nneighele(I1)
               JEL = neitabele(I1, JJEL)
               IF (JEL == IEL) cycle
               DO JED = 1, 3
                  J1 = NM(JEL, MOD(JED + 0, 3) + 1)
                  J2 = NM(JEL, MOD(JED + 1, 3) + 1)
                  IF (((J1 == I1) .AND. (J2 == I2)) .OR. &
                      ((J1 == I2) .AND. (J2 == I1))) THEN

                     IF (EDFLG(JED, JEL) == 1) THEN
                        PRINT *, 'POSSIBLE DUPLICATE ELEMENT'
                        PRINT *, 'MYPROC=', MYPROC
                        PRINT *, 'EL1=', JEL, ' EL2=', IEL, ', J1=', J1, ', J2=', J2
                        PRINT *, '  (CREATE_EDGE_DATA.F)'
                        PRINT *, 'EXECUTION WILL BE TERMINATED'
                        PRINT *, '! CHECK THE GRID CAREFULLY !'
                        PRINT *, '--------------------------------------'
                     ELSE

                        NELED(JED, JEL) = ED_ID
                        NEDEL(2, ED_ID) = JEL
                        EDFLG(JED, JEL) = 1
                        exit
                     END IF

                  END IF
               END DO
            end do
         end do
      end do

      DEALLOCATE (EDFLG)

      IF (MYPROC == 0) THEN
         PRINT *, 'DONE'
         PRINT *, ''
      END IF

      !.....An index to keep track of the edges

      DO I = 1, NEDGES
         NCOUNT(I) = -1
      END DO

      !.....Zero out internal edge counter and array

      NIEDS = 0
      NIEDN = 0

      !.....Zero out land (no-normal flow) edge counter and array

      NLEDS = 0
      NLEDN = 0

      !.....Zero out elevation specified open edge counter and array

      NEEDS = 0
      NEEDN = 0

      !.....Zero out flow specified open edge counter and array

      NFEDS = 0
      NFEDN = 0

      !.....Zero out the radiation edge counter and array

      NREDS = 0
      NREDN = 0

      !.....Zero out internal and external barrier counters and arrays

      NIBEDS = 0
      NEBEDS = 0

      NIBSEG = 0
      NEBSEG = 0

      NIBEDN = 0
      NEBEDN = 0
      NIBSEGN = 0
      NEBSEGN = 0

      !.....Find node pairs that are not edges of elements (eg. an internal
      !.....barrier against a land boundary)

      NOT_AN_EDGE = 0
      WEIR_BUDDY_NODE = 0
      DO II = 1, NVEL
         IF ((LBCODEI(II) == 4) .OR. (LBCODEI(II) == 24) .OR. &
             (LBCODEI(II) == 5) .OR. (LBCODEI(II) == 25)) THEN
            J1 = NBV(II) ! GLOBAL NODE NUMBER ON BACK SIDE OF BARRIER
            J2 = IBCONN(II) ! GLOBAL NODE NUMBER ON FRONT SIDE OF BARRIER
            DO K = 1, NBOU
               IF ((IBTYPE_ORIG(K) == 4) .OR. (IBTYPE_ORIG(K) == 24)) THEN
                  IF (WEIR_BUDDY_NODE(J1, 1) == 0) THEN
                     WEIR_BUDDY_NODE(J1, 1) = J2
                  ELSE
                     WEIR_BUDDY_NODE(J1, 2) = J2
                  END IF
               END IF
               IF ((IBTYPE_ORIG(K) == 0) .OR. (IBTYPE_ORIG(K) == 3) .OR. &
                   (IBTYPE_ORIG(K) == 2) .OR. (IBTYPE_ORIG(K) == 12) .OR. &
                   (IBTYPE_ORIG(K) == 13) .OR. (IBTYPE_ORIG(K) == 20) .OR. &
                   (IBTYPE_ORIG(K) == 22) .OR. (IBTYPE_ORIG(K) == 23)) THEN
                  DO IED = 1, NVELL(K) - 1
                     N1 = NBVV(K, IED)
                     N2 = NBVV(K, IED + 1)
                     IF ((N1 == J1) .OR. (N1 == J2)) THEN
                        IF ((N2 == J1) .OR. (N2 == J2)) THEN
                           NOT_AN_EDGE(N1) = 1
                           NOT_AN_EDGE(N2) = 1
                        END IF
                     END IF
                  END DO
               ELSEIF ((IBTYPE_ORIG(K) == 1) .OR. (IBTYPE_ORIG(K) == 11) .OR. &
                       (IBTYPE_ORIG(K) == 21)) THEN
                  DO IED = 1, NVELL(K) - 1
                     N1 = NBVV(K, IED)
                     N2 = NBVV(K, IED + 1)
                     IF ((N1 == J1) .OR. (N1 == J2)) THEN
                        IF ((N2 == J1) .OR. (N2 == J2)) THEN
                           NOT_AN_EDGE(N1) = 1
                           NOT_AN_EDGE(N2) = 1
                        END IF
                     END IF
                  END DO
               END IF
            END DO
         END IF
      end do

      !.....Find the interior edges

      DO I = 1, NEDGES
         IEL1 = NEDEL(1, I)
         IEL2 = NEDEL(2, I)
         IF ((IEL1 /= 0) .AND. (IEL2 /= 0)) THEN
            NIEDS = NIEDS + 1
            NIEDN(NIEDS) = I
            NCOUNT(I) = 0
         END IF
      END DO

      !.....Find elevation specified boundary edges

      DO K = 1, NOPE
         DO IED = 1, NVDLL(K) - 1
            N1 = NBDV(K, IED)
            N2 = NBDV(K, IED + 1)
            IERROR = 0
            DO JED = 1, NEDGES
               J1 = NEDNO(1, JED)
               J2 = NEDNO(2, JED)
               IF ((J1 == N1) .OR. (J1 == N2)) THEN
                  IF ((J2 == N1) .OR. (J2 == N2)) THEN
                     NEEDS = NEEDS + 1
                     NEEDN(NEEDS) = JED
                     NCOUNT(JED) = 4
                     IERROR = 1
                  END IF
               END IF
            END DO
            !sb-PDG1
#ifdef CMPI
            IF (IERROR == 0) THEN
               WRITE (16, *) &
                  'ERROR IN PROCESSING OPEN OCEAN BOUNDARY CONDITIONS'
            END IF
#else
            IF (IERROR == 0) THEN
               STOP 'ERROR IN PROCESSING OPEN OCEAN BOUNDARY CONDITIONS'
            END IF
#endif
            !--
         END DO
      end do

      !.....Find flux specified boundary edges

      JNMM = 0
      ONE_OR_TWO = 1
      DO K = 1, NBOU
         !sb-
         DO IED = 1, NVELL(K) - 1
            !--
            N1 = NBVV(K, IED)
            N2 = NBVV(K, IED + 1)

            IERROR = 0
            DO JED = 1, NEDGES
               J1 = NEDNO(1, JED)
               J2 = NEDNO(2, JED)
               IF ((NOT_AN_EDGE(N1) == 1) .AND. (NOT_AN_EDGE(N2) == 1)) THEN
                  !              PRINT*,'NODES ',N1,' AND ',N2,' MAKE UP A BOUNDARY',
                  !     &                  'SEGMENT THAT IS NOT AN EDGE TO AN ELEMENT; ',
                  !sb-
                  !     &                  'THIS EDGE IS SKIPPED'
                  cycle
                  !--
               END IF
               IF ((N1 == J1) .OR. (N1 == J2)) THEN
                  IF ((N2 == J1) .OR. (N2 == J2)) THEN

                     IERROR = 1
                     NCOUNT(JED) = 1

                     !.....Determine the different boundary types

                     IF ((IBTYPE_ORIG(K) == 0) .OR. (IBTYPE_ORIG(K) == 10) .OR. &
                         (IBTYPE_ORIG(K) == 20)) THEN
                        NLEDS = NLEDS + 1
                        NLEDN(NLEDS) = JED
                     END IF

                     IF ((IBTYPE_ORIG(K) == 1) .OR. (IBTYPE_ORIG(K) == 11) .OR. &
                         (IBTYPE_ORIG(K) == 21)) THEN
                        NLEDS = NLEDS + 1
                        NLEDN(NLEDS) = JED
                     END IF

                     IF ((IBTYPE_ORIG(K) == 2) .OR. (IBTYPE_ORIG(K) == 12) .OR. &
                         (IBTYPE_ORIG(K) == 22)) THEN
                        NFEDS = NFEDS + 1
                        NFEDN(NFEDS) = JED
                     END IF

                     IF ((IBTYPE_ORIG(K) == 3) .OR. (IBTYPE_ORIG(K) == 13) .OR. &
                         (IBTYPE_ORIG(K) == 23)) THEN
                        NEBSEG = NEBSEG + 1
                        NEBEDS = NEBEDS + 1
                        NEBEDN(NEBEDS) = JED
                        NEBSEGN(NEBSEG) = JED
                     END IF

                     IF ((IBTYPE_ORIG(K) == 4) .OR. (IBTYPE_ORIG(K) == 14) .OR. &
                         (IBTYPE_ORIG(K) == 24)) THEN
                        NIBSEG = NIBSEG + 1
                        NIBEDS = NIBEDS + 1
                        NIBEDN(NIBEDS) = JED
                        NIBSEGN(1, NIBSEG) = JED

                        !.....Find the edge on the opposite side of the barrier

                        NN1 = BACKNODES(1, NIBSEG)
                        NN2 = BACKNODES(2, NIBSEG)
                        !                  if (myproc.eq.24) then
                        !                     write(200+myproc,*) jed,nn1,nn2
                        !                  endif
                        DO JJED = 1, NEDGES
                           JJ1 = NEDNO(1, JJED)
                           JJ2 = NEDNO(2, JJED)
                           IF ((NN1 == JJ1) .OR. (NN1 == JJ2)) THEN
                              IF ((NN2 == JJ1) .OR. (NN2 == JJ2)) THEN
                                 NIBEDS = NIBEDS + 1
                                 NIBEDN(NIBEDS) = JJED
                                 NIBSEGN(2, NIBSEG) = JJED
                                 NCOUNT(JJED) = 1
                                 ONE_OR_TWO(N1) = 2
                                 ONE_OR_TWO(N2) = 2
                              END IF
                           END IF
                        END DO
                        !                  if (nibsegn(2,nibseg).eq.0) then
                        !                     write(200+myproc,*) 'error in create_edge_data'
                        !                     write(200+myproc,*) myproc,nibseg,jed,nn1,nn2
                        !                  endif
                     END IF

                     IF ((IBTYPE_ORIG(K) == 30)) THEN
                        NREDS = NREDS + 1
                        NREDN(NREDS) = JED
                     END IF

                  END IF
               END IF
               !sb-
            END DO
            !--
            !sb-PDG1
#ifdef CMPI
            IF (IERROR == 0) THEN
               WRITE (16, *) 'NODE PAIR (', N1, ',', N2, ') IS NOT AN EDGE.'
               WRITE (16, *) 'ERROR IN PROCESSING LAND SEGMENT'
               WRITE (16, *) ''
            END IF
#else
            IF (IERROR == 0) then
               WRITE (*, *) 'NODE PAIR (', N1, ',', N2, ') IS NOT AN EDGE.'
               WRITE (*, *) ''
               WRITE (*, *) 'ERROR IN PROCESSING LAND SEGMENT'
               WRITE (*, *) ''
               !             STOP 'ERROR IN PROCESSING LAND SEGMENT'
            end if
#endif
            !--
         end do
      end do

      !.....Check the order of the nodes assigned to an edge - important in
      !.....the calculation of the unit normal

      DO I = 1, NEDGES
         N1 = NEDNO(1, I)
         N2 = NEDNO(2, I)
         IEL = NEDEL(1, I)
         IF ((N1 == NM(IEL, 2)) .AND. (N2 == NM(IEL, 1))) THEN
            WRITE (6, *) 'THE ORDER OF NODES ASSIGNED IS WRONG FOR EDGE', I
            WRITE (16, *) 'THE ORDER OF NODES ASSIGNED IS WRONG FOR EDGE', I
            STOP
         END IF
         IF ((N1 == NM(IEL, 3)) .AND. (N2 == NM(IEL, 2))) THEN
            WRITE (6, *) 'THE ORDER OF NODES ASSIGNED IS WRONG FOR EDGE', I
            WRITE (16, *) 'THE ORDER OF NODES ASSIGNED IS WRONG FOR EDGE', I
            STOP
         END IF
         IF ((N1 == NM(IEL, 1)) .AND. (N2 == NM(IEL, 3))) THEN
            WRITE (6, *) 'THE ORDER OF NODES ASSIGNED IS WRONG FOR EDGE', I
            WRITE (16, *) 'THE ORDER OF NODES ASSIGNED IS WRONG FOR EDGE', I
            STOP
         END IF
      END DO

      !.....Check for missing edges

      !sb-PDG1 modified
      NERR = NEDGES - (NIEDS + NLEDS + NEEDS + NFEDS + NIBEDS + NEBEDS + &
                       NREDS)

      IF (MYPROC == 0) THEN
         WRITE (6, *) '  '
         WRITE (6, *) 'TOTAL NO. OF EDGES = ', NEDGES
         WRITE (6, *) '  '
         WRITE (6, *) 'NO. OF INTERNAL (NON-BOUNDARY) EDGES = ', NIEDS
         WRITE (6, *) 'NO. OF NO-NORMAL FLOW EDGES = ', NLEDS
         WRITE (6, *) 'NO. OF NON-ZERO NORMAL FLOW EDGES = ', NFEDS
         WRITE (6, *) 'NO. OF ELEVATION SPECIFIED EDGES = ', NEEDS
         WRITE (6, *) 'NO. OF EXTERNAL BARRIER EDGES = ', NEBEDS
         WRITE (6, *) 'NO. OF INTERNAL BARRIER EDGES = ', NIBEDS
         WRITE (6, *) 'NO. OF RADIATION EDGES = ', NREDS
         WRITE (6, *) &
            '-----------------------------------------------------'
         WRITE (6, *) 'NO. OF MISSING EDGES = ', NERR
         WRITE (6, *) ''
      END IF
      WRITE (16, *) '  '
      WRITE (16, *) 'TOTAL NO. OF EDGES = ', NEDGES
      WRITE (16, *) '  '
      WRITE (16, *) 'NO. OF INTERNAL (NON-BOUNDARY) EDGES = ', NIEDS
      WRITE (16, *) 'NO. OF NO-NORMAL FLOW EDGES = ', NLEDS
      WRITE (16, *) 'NO. OF NON-ZERO NORMAL FLOW EDGES = ', NFEDS
      WRITE (16, *) 'NO. OF ELEVATION SPECIFIED EDGES = ', NEEDS
      WRITE (16, *) 'NO. OF EXTERNAL BARRIER EDGES = ', NEBEDS
      WRITE (16, *) 'NO. OF INTERNAL BARRIER EDGES = ', NIBEDS
      WRITE (16, *) 'NO. OF RADIATION EDGES = ', NREDS
      WRITE (16, *) &
         '-----------------------------------------------------'
      WRITE (16, *) 'NO. OF MISSING EDGES = ', NERR
      WRITE (16, *) ''

      DO I = 1, NEDGES
         IF (NCOUNT(I) < 0) THEN
            N1 = NEDNO(1, I)
            N2 = NEDNO(2, I)
            !sb-PDG1
#ifdef CMPI
            !          WRITE(6,*)' '
            !          WRITE(6,*)'EDGE ',I,' IS MADE UP OF NODES ',N1,' AND ',N2
            !          WRITE(6,*)'EDGE ',I, ' IS NEITHER AN INTERNAL NOR A BOUNDARY',
            !     &       'MAKE SURE IF THIS IS DUE TO THE DOMAIN DECOMPOSITION'
#else
            WRITE (16, *) ' '
            WRITE (16, *) 'EDGE', I, ',MADE UP OF NODES', N1, 'AND', N2, ', IS NOT'
            WRITE (16, *) 'AN INTERNAL(NON-BOUNDARY) EDGE OR A BOUNDARY EDGE'
            WRITE (16, *) 'ASSUMING EDGE', I, 'IS A NO-NORMAL FLOW EDGE !!!'
            !          STOP
#endif
            !--
            NLEDS = NLEDS + 1
            NLEDN(NLEDS) = I
         END IF
      END DO

      !.....Add internal barrier edges to land edge table for wet-dry
      !.....post-processing

      DO I = 1, NIBEDS
         NLEDN(NLEDS + I) = NIBEDN(I)
      END DO

      !.....Construct global edge to local edge (1,2, or 3) table

      NEDSD(:, :) = 0
      DO I = 1, MNED
         DO K = 1, 2
            IF (NEDEL(K, I) /= 0) THEN
               N1 = NELED(1, NEDEL(K, I))
               N2 = NELED(2, NEDEL(K, I))
               N3 = NELED(3, NEDEL(K, I))
               IF (N1 == I) NEDSD(K, I) = 1
               IF (N2 == I) NEDSD(K, I) = 2
               IF (N3 == I) NEDSD(K, I) = 3
            END IF
         END DO
      END DO

      RETURN
   END SUBROUTINE CREATE_EDGE_DATA

   !
   !     Subroutine QUADRATURE
   !
   !     Written by Ethan Kubatko (02-05-2010)
   !
   !     This subroutine returns the weights and points of a D degree
   !     quadrature rule for a given region.  Regions may be one of the
   !     following:
   !
   !       EDGE:       -1 <= X <= 1
   !       TRIANGLE:   -1 <= X =< 1, -1 <= Y < 1,  X + Y <= 0
   !       SQUARE:     -1 <= X =< 1, -1 <= Y < 1,
   !       TRIPRISM:   -1 <= X <= 1, -1 <= Y <= 1, -1 <= Z <= 1, X + Y <= 0
   !
   !     Edge rules are standard Gauss-Legendre rules.  The remaining rules
   !     are obtained from a number of different sources listed below. The
   !     references for particular rules are noted below.
   !
   !-----------------------------------------------------------------------
   !
   !     Input:
   !     ------
   !       D:       Polynomial degree of rule
   !       REGION:  EDGE, TRIANGLE, SQUARE, or TRIPRISM (prism with a tri-
   !                angular base)
   !
   !     Output:
   !     -------
   !       PTS: Array of quadrature points of size QUAD_PTS by DIM where
   !            QUAD_PTS is the number of points for a given rule and DIM
   !            (=1 ,2, or 3) is the dimension of the region.
   !       WTS: Array of quadrature weights of size QUAD_PTS by 1.
   !
   !-----------------------------------------------------------------------

   SUBROUTINE QUADRATURE(D, REGION, PTS, WTS, N, SZ2, SZ3)
      !USE SIZES

      IMPLICIT NONE

      !.....Declare subroutine input and output

      INTEGER, intent(in) :: D, SZ2, SZ3
      CHARACTER(LEN=8), intent(in) :: REGION
      INTEGER, intent(out) :: N
      REAL(SZ), intent(inout) :: PTS(SZ2, SZ3)
      REAL(SZ), intent(inout) :: WTS(SZ2)

      INTEGER :: I, J, L
      REAL(SZ) :: AREA
      REAL(SZ) :: V1(2), V2(2), V3(2)
      REAL(SZ), ALLOCATABLE :: A(:), B(:), C(:),  W(:)
      integer, allocatable :: M(:)

      !-----------------------------------------------------------------------
      !                        EDGE QUADRATURE RULES
      !-----------------------------------------------------------------------
      !
      !     Notes: (1) These are standard Gauss-Legendre rules.
      !            (2) An n point Gauss-Legendre quadrature rule will
      !                integrate up to a 2n-1 degree polynomial exactly.
      !
      !-----------------------------------------------------------------------

      IF (REGION == 'EDGE    ') THEN
         N = CEILING(real((D + 1),8)/2d0)
         L = CEILING(real(N,8)/2d0)
         !-----------------------------------------------------------------------
         IF (D <= 1) THEN
            !-----------------------------------------------------------------------
            PTS(1, 1) = 0.000000000000000D0 ! Point
            WTS(1) = 1.00000000000000D0 ! Weight
            !-----------------------------------------------------------------------
         ELSEIF (D <= 3) THEN
            !-----------------------------------------------------------------------
            PTS(1, 1) = -0.57735026918963D0 ! Points
            WTS(1) = 1.00000000000000D0 ! Weights
            !-----------------------------------------------------------------------
         ELSEIF (D <= 5) THEN
            !-----------------------------------------------------------------------
            PTS(1, 1) = -0.77459666924148D0 ! Points
            PTS(2, 1) = 0.00000000000000D0
            WTS(1) = 0.55555555555556D0 ! Weights
            WTS(2) = 0.88888888888888D0
            !-----------------------------------------------------------------------
         ELSEIF (D <= 7) THEN
            !-----------------------------------------------------------------------
            PTS(1, 1) = -0.86113631159405D0 ! Points
            PTS(2, 1) = -0.33998104358486D0
            WTS(1) = 0.34785484513745D0 ! Weights
            WTS(2) = 0.65214515486255D0
            !-----------------------------------------------------------------------
         ELSEIF (D <= 9) THEN
            !-----------------------------------------------------------------------
            PTS(1, 1) = -0.90617984593866D0 ! Points
            PTS(2, 1) = -0.53846931010568D0
            PTS(3, 1) = 0.00000000000000D0
            WTS(1) = 0.23692688505619D0 ! Weights
            WTS(2) = 0.47862867049937D0
            WTS(3) = 0.56888888888889D0
            !-----------------------------------------------------------------------
         ELSEIF (D <= 11) THEN
            !-----------------------------------------------------------------------
            PTS(1, 1) = -0.93246951420315D0 ! Points
            PTS(2, 1) = -0.66120938646626D0
            PTS(3, 1) = -0.23861918608320D0
            WTS(1) = 0.17132449237917D0 ! Weights
            WTS(2) = 0.36076157304814D0
            WTS(3) = 0.46791393457269D0
            !-----------------------------------------------------------------------
         ELSEIF (D <= 13) THEN
            !-----------------------------------------------------------------------
            PTS(1, 1) = -0.94910791234276D0 ! Points
            PTS(2, 1) = -0.74153118559939D0
            PTS(3, 1) = -0.40584515137740D0
            PTS(4, 1) = 0.D0
            WTS(1) = 0.12948496616887D0 ! Weights
            WTS(2) = 0.27970539148928D0
            WTS(3) = 0.38183005050512D0
            WTS(4) = 0.41795918367347D0
            !-----------------------------------------------------------------------
         ELSEIF (D <= 15) THEN
            !-----------------------------------------------------------------------
            PTS(1, 1) = -0.96028985649754D0 ! Points
            PTS(2, 1) = -0.79666647741363D0
            PTS(3, 1) = -0.52553240991633D0
            PTS(4, 1) = -0.18343464249565D0
            WTS(1) = 0.10122853629038D0 ! Weights
            WTS(2) = 0.22238103445337D0
            WTS(3) = 0.31370664587789D0
            WTS(4) = 0.36268378337836D0
            !-----------------------------------------------------------------------
         ELSEIF (D <= 17) THEN
            !-----------------------------------------------------------------------
            PTS(1, 1) = -0.96816023950763D0 ! Points
            PTS(2, 1) = -0.83603110732664D0
            PTS(3, 1) = -0.61337143270059D0
            PTS(4, 1) = -0.32425342340381D0
            PTS(5, 1) = 0.00000000000000D0
            WTS(1) = 0.08127438836163D0 ! Weights
            WTS(2) = 0.18064816069483D0
            WTS(3) = 0.26061069640294D0
            WTS(4) = 0.31234707704000D0
            WTS(5) = 0.33023935500126D0
            !-----------------------------------------------------------------------
         ELSEIF (D <= 19) THEN
            !-----------------------------------------------------------------------
            PTS(1, 1) = -0.97390652851717D0 ! Points
            PTS(2, 1) = -0.86506336668898D0
            PTS(3, 1) = -0.67940956829902D0
            PTS(4, 1) = -0.43339539412925D0
            PTS(5, 1) = -0.14887433898163D0
            WTS(1) = 0.06667134430869D0 ! Weights
            WTS(2) = 0.14945134915058D0
            WTS(3) = 0.21908636251598D0
            WTS(4) = 0.26926671931000D0
            WTS(5) = 0.29552422471475D0
            !-----------------------------------------------------------------------
         ELSEIF (D <= 21) THEN
            !-----------------------------------------------------------------------
            PTS(1, 1) = -0.97822865814606D0 ! Points
            PTS(2, 1) = -0.88706259976810D0
            PTS(3, 1) = -0.73015200557405D0
            PTS(4, 1) = -0.51909612920681D0
            PTS(5, 1) = -0.26954315595234D0
            PTS(6, 1) = 0.00000000000000D0
            WTS(1) = 0.05566856711627D0 ! Weights
            WTS(2) = 0.12558036946485D0
            WTS(3) = 0.18629021092774D0
            WTS(4) = 0.23319376459199D0
            WTS(5) = 0.26280454451025D0
            WTS(6) = 0.27292508677790D0
            !-----------------------------------------------------------------------
         ELSEIF (D <= 23) THEN
            !-----------------------------------------------------------------------
            PTS(1, 1) = -0.98156063424672D0 ! Points
            PTS(2, 1) = -0.90411725637047D0
            PTS(3, 1) = -0.76990267419430D0
            PTS(4, 1) = -0.58731795428662D0
            PTS(5, 1) = -0.36783149899818D0
            PTS(6, 1) = -0.12523340851147D0
            WTS(1) = 0.04717533638677D0 ! Weights
            WTS(2) = 0.10693932599520D0
            WTS(3) = 0.16007832854334D0
            WTS(4) = 0.20316742672308D0
            WTS(5) = 0.23349253653835D0
            WTS(6) = 0.24914704581340D0
            !-----------------------------------------------------------------------
         ELSEIF (D <= 25) THEN
            !-----------------------------------------------------------------------
            PTS(1, 1) = -0.98418305471859D0 ! Points
            PTS(2, 1) = -0.91759839922298D0
            PTS(3, 1) = -0.80157809073331D0
            PTS(4, 1) = -0.64234933944034D0
            PTS(5, 1) = -0.44849275103645D0
            PTS(6, 1) = -0.23045831595513D0
            PTS(7, 1) = 0.00000000000000D0
            WTS(1) = 0.04048400476532D0 ! Weights
            WTS(2) = 0.09212149983773D0
            WTS(3) = 0.13887351021979D0
            WTS(4) = 0.17814598076195D0
            WTS(5) = 0.20781604753689D0
            WTS(6) = 0.22628318026290D0
            WTS(7) = 0.23255155323087D0
            !-----------------------------------------------------------------------
         ELSE
            !-----------------------------------------------------------------------
            call quad_rules_general(N, .false., 0.D0, 0.D0, pts, wts)

         END IF
         !-----------------------------------------------------------------------
         J = MOD(N, 2)
         DO I = L, N
            PTS(I + 1, 1) = -PTS(L - J, 1)
            WTS(I + 1) = WTS(L - J)
            J = J + 1
         END DO
         !-----------------------------------------------------------------------
         !                        TRIANGLE QUADRATURE RULES
         !-----------------------------------------------------------------------
         !
         !     All rules have positive weights with all points located inside the
         !     triangle (so-called PI rules). Rules marked by * are optimal.
         !
         !     References:

         !     [1] D.A. Dunavant, "High Degree Efficient Symmetrical Gaussian
         !         Quadrature Rules for the Triangle", International Journal for
         !         Numerical Methods in Engineering, 21, 1129--1148, 1985.
         !
         !     [2] L. Zhang, T. Cui, and Hui Liu, "A Set of Symmetric Quadrature
         !         Rules on Triangles and Tetrahedra", Journal of Computational
         !         Mathematics, 27, 89--96, 2009.
         !
         !     [3] P. Hillion, "Numerical Integration on a triangle",
         !         International Journal for Numerical Methods in Engineering,
         !         11, 797--815, 1977.
         !
         !     [4] K. Gatermann, "The construction of symmetric cubature formulas
         !         for the square and the triangle, Computing, 40, 229--240,
         !         1988.
         !
         !-----------------------------------------------------------------------

      ELSEIF (REGION == 'TRIANGLE') THEN

         !-----------------------------------------------------------------------
         IF (D <= 1) THEN ! Ref [1]*
            !-----------------------------------------------------------------------
            L = 1
            ALLOCATE (A(L), B(L), C(L), M(L), W(L))
            M(1:L) = [1] ! Multip.
            A(1:L) = [0.333333333333333D0] ! Points
            W(1:L) = [1.00000000000000D0] ! Weights
            !-----------------------------------------------------------------------
         ELSEIF (D <= 2) THEN ! Ref [1]*
            !-----------------------------------------------------------------------
            L = 1
            ALLOCATE (A(L), B(L), C(L), M(L), W(L))
            M(1:L) = [3] ! Multip.
            A(1:L) = [0.166666666666666D0] ! Points
            W(1:L) = [0.333333333333333D0] ! Weights
            !-----------------------------------------------------------------------
            !        ELSEIF (D.LE.3) THEN                                  ! Ref [3]*
            !-----------------------------------------------------------------------
            !-----------------------------------------------------------------------
         ELSEIF (D <= 4) THEN ! Ref [2]*
            !-----------------------------------------------------------------------
            L = 2
            ALLOCATE (A(L), B(L), C(L), M(L), W(L))
            M(1:L) = [3, 3] ! Multip.
            A(1:L) = [0.445948490915964D0, &
                      0.0915762135097707D0]
            W(1:L) = [0.223381589678011D0, &
                      0.109951743655321D0]
            !-----------------------------------------------------------------------
         ELSEIF (D <= 5) THEN ! Ref [1]*
            !-----------------------------------------------------------------------
            L = 3
            ALLOCATE (A(L), B(L), C(L), M(L), W(L))
            M(1:L) = [1, 3, 3] ! Multip.
            A(1:L) = [0.333333333333333D0, &
                      0.470142064105115D0, &
                      0.101286507323456D0]
            W(1:L) = [0.225000000000000D0, &
                      0.132394152788506D0, &
                      0.125939180544827D0]
            !-----------------------------------------------------------------------
         ELSEIF (D <= 6) THEN ! Ref [1]
            !-----------------------------------------------------------------------
            L = 3
            ALLOCATE (A(L), B(L), C(L), M(L), W(L))
            M(1:L) = [3, 3, 6] ! Multip.
            A(1:L) = [0.0630890144915022D0, &
                      0.249286745170910D0, &
                      0.0531450498448169D0]
            B(3:3) = [0.310352451033784D0]
            W(1:L) = [0.0508449063702068D0, &
                      0.116786275726379D0, &
                      0.0828510756183735D0]
            !-----------------------------------------------------------------------
            !        ELSEIF (D.LE.7) THEN                                  ! Ref [4]*
            !-----------------------------------------------------------------------
            !          L = 4
            !          ALLOCATE( A(L), B(L), C(L), M(L), W(L) )
            !          M(1:L) = (/ 3,3,3,3 /)
            !          A(1:L) = (/ /)
            !          B(1:L) = (/ /)
            !          W(1:L) = (/ /)
            !-----------------------------------------------------------------------
         ELSEIF (D <= 8) THEN ! Ref [1]
            !-----------------------------------------------------------------------
            L = 5
            ALLOCATE (A(L), B(L), C(L), M(L), W(L))
            M(1:L) = [1, 3, 3, 3, 6] ! Multip.
            A(1:L) = [0.333333333333333D0, &
                      0.170569307751760D0, &
                      0.0505472283170309D0, &
                      0.459292588292723D0, &
                      0.263112829634638D0]
            B(5:5) = [0.00839477740995760D0]
            W(1:L) = [0.144315607677787D0, &
                      0.103217370534718D0, &
                      0.0324584976231980D0, &
                      0.0950916342672846D0, &
                      0.0272303141744349D0]
            !-----------------------------------------------------------------------
         ELSEIF (D <= 9) THEN ! Ref [1]
            !-----------------------------------------------------------------------
            L = 6
            ALLOCATE (A(L), B(L), C(L), M(L), W(L))
            M(1:L) = [1, 3, 3, 3, 3, 6] ! Multip.
            A(1:L) = [0.333333333333333D0, &
                      0.489682519198737D0, &
                      0.0447295133944527D0, &
                      0.437089591492936D0, &
                      0.188203535619032D0, &
                      0.741198598784498D0]
            B(6:6) = [0.221962989160765D0]
            W(1:L) = [0.0971357962827988D0, &
                      0.0313347002271390D0, &
                      0.0255776756586980D0, &
                      0.0778275410047742D0, &
                      0.0796477389272102D0, &
                      0.0432835393772893D0]
            !-----------------------------------------------------------------------
         ELSEIF (D <= 10) THEN ! Ref [2]
            !-----------------------------------------------------------------------
            L = 7
            ALLOCATE (A(L), B(L), C(L), M(L), W(L))
            M(1:L) = [1, 3, 3, 3, 3, 6, 6] ! Multip.
            A(1:L) = [0.333333333333333D0, &
                      0.427273178846775D0, &
                      0.183099222448675D0, &
                      0.490434019701130D0, &
                      0.0125724455515805D0, &
                      0.654268667920066D0, &
                      0.122804577068559D0]
            B(6:7) = [0.308046001685247D0, &
                      0.0333718337393047D0]
            W(1:L) = [0.0809374287976228D0, &
                      0.0772985880029631D0, &
                      0.0784576386123717D0, &
                      0.0174691679959294D0, &
                      0.00429237418483282D0, &
                      0.0374688582104676D0, &
                      0.0269493525918799D0]
            !-----------------------------------------------------------------------
         ELSEIF (D <= 11) THEN ! Ref [2]
            !-----------------------------------------------------------------------
            L = 8
            ALLOCATE (A(L), B(L), C(L), M(L), W(L))
            M(1:L) = [1, 3, 3, 3, 3, 3, 6, 6] ! Multip.
            A(1:L) = [0.333333333333333D0, &
                      0.0309383552454307D0, &
                      0.436498181134128D0, &
                      0.498984763702593D0, &
                      0.214688197958594D0, &
                      0.113683104042113D0, &
                      0.825618766164862D0, &
                      0.640472310134865D0]
            B(7:8) = [0.159742304591850D0, &
                      0.311783715709599D0]
            W(1:L) = [0.0811779602968671D0, &
                      0.0123240435069094D0, &
                      0.0628280097444101D0, &
                      0.0122203790493645D0, &
                      0.0677013489528115D0, &
                      0.0402196936288516D0, &
                      0.0147622727177161D0, &
                      0.0407279964582990D0]
            !-----------------------------------------------------------------------
         ELSEIF (D <= 12) THEN ! Ref [2]
            !-----------------------------------------------------------------------
            L = 8
            ALLOCATE (A(L), B(L), C(L), M(L), W(L))
            M(1:L) = [3, 3, 3, 3, 3, 6, 6, 6] ! Multip.
            A(1:L) = [0.0213173504532103D0, &
                      0.271210385012115D0, &
                      0.127576145541585D0, &
                      0.439724392294460D0, &
                      0.488217389773804D0, &
                      0.695836086787803D0, &
                      0.858014033544072D0, &
                      0.608943235779787D0]
            B(6:8) = [0.281325580989939D0, &
                      0.116251915907597D0, &
                      0.275713269685514D0]
            W(1:L) = [0.00616626105155901D0, &
                      0.0628582242178851D0, &
                      0.0347961129307089D0, &
                      0.0436925445380384D0, &
                      0.0257310664404553D0, &
                      0.0223567732023034D0, &
                      0.0173162311086588D0, &
                      0.0403715577663809D0]
            !-----------------------------------------------------------------------
         ELSEIF (D <= 13) THEN ! Ref [2]
            !-----------------------------------------------------------------------
            L = 9
            ALLOCATE (A(L), B(L), C(L), M(L), W(L))
            M(1:L) = [1, 3, 3, 3, 3, 6, 6, 6, 6] ! Multip.
            A(1:L) = [0.333333333333333D0, &
                      0.426941414259800D0, &
                      0.221372286291832D0, &
                      0.0215096811088431D0, &
                      0.489076946452539D0, &
                      0.623545995553675D0, &
                      0.864707770295442D0, &
                      0.748507115899952D0, &
                      0.722357793124187D0]
            B(6:9) = [0.308441760892117D0, &
                      0.110922042803463D0, &
                      0.163597401067850D0, &
                      0.272515817773429D0]
            W(1:L) = [0.0679600365868316D0, &
                      0.0556019675304533D0, &
                      0.0582784851191999D0, &
                      0.00605233710353917D0, &
                      0.0239944019288947D0, &
                      0.0346412761408483D0, &
                      0.0149654011051656D0, &
                      0.0241790398115938D0, &
                      0.00959068100354326D0]
            !-----------------------------------------------------------------------
         ELSEIF (D <= 14) THEN ! Ref [1]
            !-----------------------------------------------------------------------
            L = 10
            ALLOCATE (A(L), B(L), C(L), M(L), W(L))
            M(1:L) = [3, 3, 3, 3, 3, 3, 6, 6, 6, 6] ! Multipl.
            A(1:L) = [0.488963910362179D0, &
                      0.417644719340454D0, &
                      0.273477528308839D0, &
                      0.177205532412543D0, &
                      0.0617998830908730D0, &
                      0.0193909612487010D0, &
                      0.172266687821356D0, &
                      0.336861459796345D0, &
                      0.298372882136258D0, &
                      0.118974497696957D0]
            B(7:L) = [0.770608554774996D0, &
                      0.570222290846683D0, &
                      0.686980167808088D0, &
                      0.879757171370171D0]
            W(1:L) = [0.0218835813694290D0, &
                      0.0327883535441250D0, &
                      0.0517741045072920D0, &
                      0.0421625887369930D0, &
                      0.0144336996697770D0, &
                      0.00492340360240000D0, &
                      0.0246657532125640D0, &
                      0.0385715107870610D0, &
                      0.0144363081135340D0, &
                      0.00501022883850100D0]
            !-----------------------------------------------------------------------
         ELSEIF (D <= 15) THEN ! Ref [2]
            !-----------------------------------------------------------------------
            L = 13
            ALLOCATE (A(L), B(L), C(L), M(L), W(L))
            M(1:L) = [1, 3, 3, 3, 3, 3, 3, 3, 6, 6, 6, 6, 6] ! Multip.
            A(1:L) = [0.333333333333333D0, &
                      0.227332218819142D0, &
                      0.497162577431887D0, &
                      0.478849735348954D0, &
                      0.404986039098271D0, &
                      0.0159312166717444D0, &
                      0.165583262426081D0, &
                      0.0731336047192287D0, &
                      0.665260733072213D0, &
                      0.712521987242545D0, &
                      0.559648362235393D0, &
                      0.810476597619076D0, &
                      0.916075644031731D0]
            B(9:L) = [0.316352839344947D0, &
                      0.0934607511499175D0, &
                      0.344229017582193D0, &
                      0.171047248314257D0, &
                      0.0730559964791864D0]
            W(1:L) = [0.0440387108784342D0, &
                      0.0461847871820269D0, &
                      0.00649890661733271D0, &
                      0.0179936142526584D0, &
                      0.0417731050391413D0, &
                      0.00305954760911646D0, &
                      0.00201243505255864D0, &
                      0.0167756109305091D0, &
                      0.0154607491897142D0, &
                      0.0284998903395474D0, &
                      0.0320943504834895D0, &
                      0.0115085816368707D0, &
                      0.00461430652896710D0]
            !-----------------------------------------------------------------------
         ELSEIF (D <= 16) THEN ! Ref [2]
            !-----------------------------------------------------------------------
            L = 13
            ALLOCATE (A(L), B(L), C(L), M(L), W(L))
            M(1:L) = [1, 3, 3, 3, 3, 3, 3, 6, 6, 6, 6, 6, 6] ! Multip.
            A(1:L) = [0.333333333333333D0, &
                      0.0817949831313738D0, &
                      0.165300601969779D0, &
                      0.468592105349461D0, &
                      0.0144388134454166D0, &
                      0.241784285391783D0, &
                      0.495310342987769D0, &
                      0.650513402661352D0, &
                      0.604011281495997D0, &
                      0.802168257574741D0, &
                      0.756505606442828D0, &
                      0.465938438714118D0, &
                      0.906394843992041D0]
            B(8:L) = [0.331399744537089D0, &
                      0.303247162749942D0, &
                      0.188028059521237D0, &
                      0.183504668522296D0, &
                      0.359645948797504D0, &
                      0.0771943712957554D0]
            W(1:L) = [0.0480221886803770D0, &
                      0.0147091003068019D0, &
                      0.0295445865493192D0, &
                      0.0261250173510883D0, &
                      0.00278038735239000D0, &
                      0.0318217730005366D0, &
                      0.00864583434950965D0, &
                      0.0143003329044953D0, &
                      0.0278497772036008D0, &
                      0.00704167340663609D0, &
                      0.0178998382599337D0, &
                      0.0274582003843497D0, &
                      0.00729979693943176D0]
            !-----------------------------------------------------------------------
         ELSEIF (D <= 17) THEN ! Ref [2]
            !-----------------------------------------------------------------------
            L = 14
            ALLOCATE (A(L), B(L), C(L), M(L), W(L))
            M(1:L) = [1, 3, 3, 3, 3, 3, 3, 6, 6, 6, 6, 6, 6, 6] ! Multip.
            A(1:L) = [0.333333333333333D0, &
                      0.0956985088627109D0, &
                      0.170138639678775D0, &
                      0.418020685867954D0, &
                      0.496581480506624D0, &
                      0.0416621148288076D0, &
                      0.467932905729423D0, &
                      0.969531198903722D0, &
                      0.759724387538624D0, &
                      0.295499316968301D0, &
                      0.625606382157697D0, &
                      0.872174447233184D0, &
                      0.747512319440006D0, &
                      0.598868790883238D0]
            B(8:L) = [0.0289250916202182D0, &
                      0.234441755263568D0, &
                      0.495911246660753D0, &
                      0.353417694541497D0, &
                      0.112728641814219D0, &
                      0.199070278797857D0, &
                      0.303585183071326D0]
            W(1:L) = [0.0447568714443446D0, &
                      0.0173668850267477D0, &
                      0.0305993480761035D0, &
                      0.0285877085785997D0, &
                      0.00664743192975369D0, &
                      0.00747618940201851D0, &
                      0.0250499865038387D0, &
                      0.00147981089211964D0, &
                      0.00512113624674810D0, &
                      0.0273173593695928D0, &
                      0.0140057286759092D0, &
                      0.00780927569745836D0, &
                      0.0181657284597916D0, &
                      0.0274443739924583D0]
            !-----------------------------------------------------------------------
         ELSEIF (D <= 18) THEN ! Ref [2]
            !-----------------------------------------------------------------------
            L = 14
            ALLOCATE (A(L), B(L), C(L), M(L), W(L))
            M(1:L) = [3, 3, 3, 3, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6] ! Multipl.
            A(1:L) = [0.0732708864643828D0, &
                      0.00391774898322823D0, &
                      0.467597318988711D0, &
                      0.417916210967411D0, &
                      0.165381693360289D0, &
                      0.287500894405783D0, &
                      0.125889314319824D0, &
                      0.0632219159465026D0, &
                      0.0789102274540205D0, &
                      0.0380580535067857D0, &
                      0.0142903521304540D0, &
                      0.0129672723432531D0, &
                      0.00764859482084089D0, &
                      0.0127104605722554D0]
            B(5:L) = [0.563696705660870D0, &
                      0.286042326139204D0, &
                      0.696043218642461D0, &
                      0.760545551887682D0, &
                      0.592019631271758D0, &
                      0.683681259635999D0, &
                      0.851704037137055D0, &
                      0.574732492888149D0, &
                      0.735510440830729D0, &
                      0.939345087643731D0]
            W(1:L) = [0.0139778616452860D0, &
                      0.000554906979213213D0, &
                      0.0210268138197046D0, &
                      0.0340182121799276D0, &
                      0.0279101658047749D0, &
                      0.0182146861271508D0, &
                      0.0142670236581097D0, &
                      0.0142371230906750D0, &
                      0.0192575838546747D0, &
                      0.00970513228438064D0, &
                      0.00762978813433212D0, &
                      0.0106187391363503D0, &
                      0.00571066980327583D0, &
                      0.00432685746087641D0]
            !-----------------------------------------------------------------------
         ELSEIF (D <= 19) THEN ! Ref [2]
            !-----------------------------------------------------------------------
            L = 17
            ALLOCATE (A(L), B(L), C(L), M(L), W(L))
            M(1:L) = [1, 3, 3, 3, 3, 3, 3, 3, 3, 6, 6, 6, 6, 6, 6, 6, 6] ! Multip.
            A(1:L) = [0.333333333333333D0, &
                      0.489609987073006D0, &
                      0.454536892697892D0, &
                      0.401416680649431D0, &
                      0.255551654403097D0, &
                      0.177077942152129D0, &
                      0.110061053227951D0, &
                      0.0555286242518396D0, &
                      0.0126218637772286D0, &
                      0.600633794794645D0, &
                      0.134466754530779D0, &
                      0.720987025817365D0, &
                      0.594527068955870D0, &
                      0.839331473680838D0, &
                      0.223861424097915D0, &
                      0.822931324069856D0, &
                      0.924344252620784D0]
            B(10:L) = [0.395754787356942D0, &
                       0.557603261588783D0, &
                       0.264566948406520D0, &
                       0.358539352205950D0, &
                       0.157807405968594D0, &
                       0.701087978926173D0, &
                       0.142421601113383D0, &
                       0.0654946280829377D0]
            W(1:L) = [0.0329063313889186D0, &
                      0.0103307318912720D0, &
                      0.0223872472630163D0, &
                      0.0302661258694680D0, &
                      0.0304909678021977D0, &
                      0.0241592127416409D0, &
                      0.0160508035868008D0, &
                      0.00808458026178406D0, &
                      0.00207936202748478D0, &
                      0.00388487690498138D0, &
                      0.0255741606120219D0, &
                      0.00888090357333805D0, &
                      0.0161245467617313D0, &
                      0.00249194181749067D0, &
                      0.0182428401189505D0, &
                      0.0102585637361985D0, &
                      0.00379992885530191D0]
            !-----------------------------------------------------------------------
         ELSEIF (D <= 20) THEN ! Ref [2]
            !-----------------------------------------------------------------------
            L = 18
            ALLOCATE (A(L), B(L), C(L), M(L), W(L))
            M(1:L) = [1, 3, 3, 3, 3, 3, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6] ! Multip.
            A(1:L) = [0.333333333333333D0, &
                      0.215874305932991D0, &
                      0.0753767665297472D0, &
                      0.0103008281372217D0, &
                      0.493602211298700D0, &
                      0.461550938106925D0, &
                      0.328621406424236D0, &
                      0.260480361786568D0, &
                      0.137074235846455D0, &
                      0.146726945872299D0, &
                      0.0269989777425532D0, &
                      0.0618717859336170D0, &
                      0.0477243674276219D0, &
                      0.120600515186364D0, &
                      0.00269714779670978D0, &
                      0.00301563327794236D0, &
                      0.0299053757884570D0, &
                      0.00675665422246098D0]
            B(7:L) = [0.429340570258210D0, &
                      0.101577534280969D0, &
                      0.710065973001130D0, &
                      0.498545477678414D0, &
                      0.0491867226725820D0, &
                      0.779660146540569D0, &
                      0.370491539149547D0, &
                      0.863346948754752D0, &
                      0.0561949381877455D0, &
                      0.208675006748421D0, &
                      0.721151240912034D0, &
                      0.640055441940541D0]
            W(1:L) = [0.0125376079944966D0, &
                      0.0274718698764242D0, &
                      0.00976527227705142D0, &
                      0.00139841953539182D0, &
                      0.00929210262518518D0, &
                      0.0165778760323669D0, &
                      0.0206677623486650D0, &
                      0.0208222355211545D0, &
                      0.00956863841984906D0, &
                      0.0244527709689724D0, &
                      0.00315573063063053D0, &
                      0.0121367963653212D0, &
                      0.0149664801438864D0, &
                      0.00632759332177773D0, &
                      0.00134256031206369D0, &
                      0.00277607691634755D0, &
                      0.0107398444741849D0, &
                      0.00536780573818745D0]
            !-----------------------------------------------------------------------
         ELSEIF (D <= 21) THEN ! Ref [2]
            !-----------------------------------------------------------------------
            L = 19
            ALLOCATE (A(L), B(L), C(L), M(L), W(L))
            M(1:L) = [1, 3, 3, 3, 3, 3, 3, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6] ! Multip.
            A(1:L) = [0.333333333333333D0, &
                      0.200935277065085D0, &
                      0.437659165961927D0, &
                      0.00343395649059617D0, &
                      0.0466434847753067D0, &
                      0.386422251763071D0, &
                      0.0954354711085309D0, &
                      0.955513803350456D0, &
                      0.886638813428868D0, &
                      0.784262845880434D0, &
                      0.882923955050200D0, &
                      0.668991964441077D0, &
                      0.552072121035560D0, &
                      0.797592965596568D0, &
                      0.677514715119771D0, &
                      0.542997415589091D0, &
                      0.705459905569968D0, &
                      0.574800573066508D0, &
                      0.471778808504614D0]
            B(8:L) = [0.0357186278731633D0, &
                      0.108143224915646D0, &
                      0.207464449599876D0, &
                      0.0856847087203169D0, &
                      0.321494003014288D0, &
                      0.437942218793341D0, &
                      0.161916453063577D0, &
                      0.274504767401994D0, &
                      0.405335998075006D0, &
                      0.187737680656435D0, &
                      0.305696834766055D0, &
                      0.312144466870890D0]
            W(1:L) = [0.0275622569528764D0, &
                      0.0220602154134885D0, &
                      0.0234600159386714D0, &
                      0.000326889595047190D0, &
                      0.00326531946293996D0, &
                      0.0117564629154127D0, &
                      0.0117807684199115D0, &
                      0.00226881081880114D0, &
                      0.00259601096443632D0, &
                      0.00463452978587186D0, &
                      0.00479433605454885D0, &
                      0.00571247883672361D0, &
                      0.00586582760432212D0, &
                      0.00941376305909158D0, &
                      0.0134149437966564D0, &
                      0.0157169180920832D0, &
                      0.0168636830144369D0, &
                      0.0213900270853200D0, &
                      0.0230767921894926D0]
            !-----------------------------------------------------------------------
         ELSE
            !-----------------------------------------------------------------------
            PRINT *, '  ******** ERROR!!! D must be 1 <= D <= 21 ********  '
            PRINT *, '  ********      for REGION = TRIANGLE      ********  '
            PRINT *, '    Execution terminated in subroutine quadrature    '
            STOP
         END IF

         !.......Vertices and area for the master triangular element

         V1 = [-1.D0, -1.D0]
         V2 = [1.D0, -1.D0]
         V3 = [-1.D0, 1.D0]
         AREA = 2.D0

         !.......Count up number of quadrature points

         N = 0
         DO I = 1, L
            N = N + M(I)
         END DO
         !      ALLOCATE ( PTS(N,2), WTS(N) )
         !.......Transform quadrature points to master element coordinates

         J = 1
         DO I = 1, L
            IF (M(I) == 1) THEN
               B(I) = A(I)
               C(I) = B(I)
               PTS(J, 1) = V1(1)*A(I) + V2(1)*B(I) + V3(1)*C(I)
               PTS(J, 2) = V1(2)*A(I) + V2(2)*B(I) + V3(2)*C(I)
               WTS(J) = AREA*W(I)
               J = J + 1
            ELSEIF (M(I) == 3) THEN
               B(I) = A(I)
               C(I) = 1.D0 - (A(I) + B(I))
               PTS(J, 1) = V1(1)*A(I) + V2(1)*B(I) + V3(1)*C(I)
               PTS(J, 2) = V1(2)*A(I) + V2(2)*B(I) + V3(2)*C(I)
               WTS(J) = AREA*W(I)
               PTS(J + 1, 1) = V1(1)*B(I) + V2(1)*C(I) + V3(1)*A(I)
               PTS(J + 1, 2) = V1(2)*B(I) + V2(2)*C(I) + V3(2)*A(I)
               WTS(J + 1) = AREA*W(I)
               PTS(J + 2, 1) = V1(1)*C(I) + V2(1)*A(I) + V3(1)*B(I)
               PTS(J + 2, 2) = V1(2)*C(I) + V2(2)*A(I) + V3(2)*B(I)
               WTS(J + 2) = AREA*W(I)
               J = J + 3
            ELSEIF (M(I) == 6) THEN
               C(I) = 1.D0 - (A(I) + B(I))
               PTS(J, 1) = V1(1)*A(I) + V2(1)*B(I) + V3(1)*C(I)
               PTS(J, 2) = V1(2)*A(I) + V2(2)*B(I) + V3(2)*C(I)
               WTS(J) = AREA*W(I)
               PTS(J + 1, 1) = V1(1)*B(I) + V2(1)*C(I) + V3(1)*A(I)
               PTS(J + 1, 2) = V1(2)*B(I) + V2(2)*C(I) + V3(2)*A(I)
               WTS(J + 1) = AREA*W(I)
               PTS(J + 2, 1) = V1(1)*C(I) + V2(1)*A(I) + V3(1)*B(I)
               PTS(J + 2, 2) = V1(2)*C(I) + V2(2)*A(I) + V3(2)*B(I)
               WTS(J + 2) = AREA*W(I)
               PTS(J + 3, 1) = V1(1)*A(I) + V2(1)*C(I) + V3(1)*B(I)
               PTS(J + 3, 2) = V1(2)*A(I) + V2(2)*C(I) + V3(2)*B(I)
               WTS(J + 3) = AREA*W(I)
               PTS(J + 4, 1) = V1(1)*B(I) + V2(1)*A(I) + V3(1)*C(I)
               PTS(J + 4, 2) = V1(2)*B(I) + V2(2)*A(I) + V3(2)*C(I)
               WTS(J + 4) = AREA*W(I)
               PTS(J + 5, 1) = V1(1)*C(I) + V2(1)*B(I) + V3(1)*A(I)
               PTS(J + 5, 2) = V1(2)*C(I) + V2(2)*B(I) + V3(2)*A(I)
               WTS(J + 5) = AREA*W(I)
               J = J + 6
            END IF
         END DO

         !-----------------------------------------------------------------------
         !                        SQUARE QUADRATURE RULES
         !-----------------------------------------------------------------------
         !
         !     All rules have positive weights (with the exception of rule 23,
         !     which has 1 negative weight) with all points located inside the
         !     square (so-called PI rules). Rules marked by * are optimal.  The
         !     rest are the best PI rules currently known.
         !
         !     References:

         !     [1] A.H. Stroud, "Approximate calculation of multiple integrals",
         !         Prentice-Hall, Englewood Cliffs, N.J., 1971.
         !
         !     [2] J.W. Wissman and T. Becker, "Partially Symmetric Cubature
         !         Formulas for Even Degrees of Exactness", SIAM Journal on
         !         Numerical Analysis, 23, 676--685, 1986.
         !
         !     [3] Construction of Cubature Formulas of Degree Seven and Nine
         !         Using Symmetric Planar Regions, Using Orthogonal Polynmials",
         !         SIAM Journal on Numerical Analysis, 14, 492--508, 1977.
         !
         !     [4] I.P. Omelyan and V.B. Solovyan, "Improved cubature formulae of
         !         high degrees of exactness for the square", Journal of Compu-
         !         tational and Applied Mathematics, 188, 190--204, 2006.
         !
         !     [5] H.M. Moller, "Minimum-Point Cubature Formula", Numerische
         !         Mathematik, 25, 185--200, 1976. (in German)
         !
         !-----------------------------------------------------------------------

      ELSEIF (REGION == 'SQUARE  ') THEN

         !-----------------------------------------------------------------------
         IF (D <= 1) THEN ! Ref [1]*
            !-----------------------------------------------------------------------
            L = 1
            ALLOCATE (A(L), B(L), M(L), W(L))
            M(1:L) = [1] ! Multiplicity
            A(1:L) = [0.000000000000000D0] ! Point
            B(1:L) = [0.000000000000000D0]
            W(1:L) = [4.00000000000000D0] ! Weight
            !-----------------------------------------------------------------------
         ELSEIF (D <= 2) THEN ! Ref [1]
            !-----------------------------------------------------------------------
            L = 3
            ALLOCATE (A(L), B(L), M(L), W(L))
            M(1:L) = [1, 1, 1]
            A(1:L) = [0.81649658092773D0, &
                      -0.40824829046386D0, &
                      -0.40824829046386D0]
            B(1:L) = [0.00000000000000D0, &
                      0.70710678118655D0, &
                      -0.70710678118655D0]
            W(1:L) = [1.33333333333333D0, &
                      1.33333333333333D0, &
                      1.33333333333333D0]
            !-----------------------------------------------------------------------
         ELSEIF (D <= 3) THEN ! Ref [1]*
            !-----------------------------------------------------------------------
            L = 1
            ALLOCATE (A(L), B(L), M(L), W(L))
            M(1:L) = [4] ! Multiplicity
            A(1:L) = [0.577350269189626D0] ! Points
            B(1:L) = [0.577350269189626D0]
            W(1:L) = [1.00000000000000D0] ! Weights
            !-----------------------------------------------------------------------
         ELSEIF (D <= 4) THEN ! Ref [2]*
            !-----------------------------------------------------------------------
            L = 6
            ALLOCATE (A(L), B(L), M(L), W(L))
            M(1:L) = [1, 1, 1, 1, 1, 1] ! Multipl.
            A(1:L) = [0.000000000000000D0, &
                      0.000000000000000D0, &
                      0.851914653304601D0, &
                      0.630912788976754D0, &
                      -0.851914653304601D0, &
                      -0.630912788976754D0]
            B(1:L) = [0.000000000000000D0, &
                      0.966091783079296D0, &
                      0.455603727836193D0, &
                      -0.731629951573135D0, &
                      0.455603727836193D0, &
                      -0.731629951573135D0]
            W(1:L) = [1.14285714285714D0, &
                      0.439560439560440D0, &
                      0.566072207007532D0, &
                      0.642719001783677D0, &
                      0.566072207007532D0, &
                      0.642719001783677D0]
            !-----------------------------------------------------------------------
         ELSEIF (D <= 5) THEN ! Ref [1]*
            !-----------------------------------------------------------------------
            L = 7
            ALLOCATE (A(L), B(L), M(L), W(L))
            M(1:L) = [1, 1, 1, 1, 1, 1, 1]
            A(1:L) = [0.00000000000000D0, &
                      0.00000000000000D0, &
                      0.00000000000000D0, &
                      0.77459666924148D0, &
                      0.77459666924148D0, &
                      -0.77459666924148D0, &
                      -0.77459666924148D0]
            B(1:L) = [0.00000000000000D0, &
                      0.96609178307930D0, &
                      -0.96609178307930D0, &
                      0.57735026918963D0, &
                      -0.57735026918963D0, &
                      0.57735026918963D0, &
                      -0.57735026918963D0]
            W(1:L) = [1.14285714285714D0, &
                      0.31746031746032D0, &
                      0.31746031746032D0, &
                      0.55555555555556D0, &
                      0.55555555555556D0, &
                      0.55555555555556D0, &
                      0.55555555555556D0]
            !-----------------------------------------------------------------------
         ELSEIF (D <= 6) THEN ! Ref [2]*
            !-----------------------------------------------------------------------
            L = 10
            ALLOCATE (A(L), B(L), M(L), W(L))
            M(1:L) = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
            A(1:L) = [0.000000000000000D0, &
                      0.000000000000000D0, &
                      0.888764014654765D0, &
                      0.604857639464685D0, &
                      0.955447506641064D0, &
                      0.565459993438754D0, &
                      -0.888764014654765D0, &
                      -0.604857639464685D0, &
                      -0.955447506641064D0, &
                      -0.565459993438754D0]
            B(1:L) = [0.836405633697626D0, &
                      -0.357460165391307D0, &
                      0.872101531193131D0, &
                      0.305985162155427D0, &
                      -0.410270899466658D0, &
                      -0.872869311156879D0, &
                      0.872101531193131D0, &
                      0.305985162155427D0, &
                      -0.410270899466658D0, &
                      -0.872869311156879D0]
            W(1:L) = [0.455343245714174D0, &
                      0.827395973202966D0, &
                      0.144000884599645D0, &
                      0.668259104262665D0, &
                      0.225474004890679D0, &
                      0.320896396788441D0, &
                      0.144000884599645D0, &
                      0.668259104262665D0, &
                      0.225474004890679D0, &
                      0.320896396788441D0]
            !-----------------------------------------------------------------------
         ELSEIF (D <= 7) THEN ! Ref [3]*
            !-----------------------------------------------------------------------
            L = 4
            ALLOCATE (A(L), B(L), M(L), W(L))
            M(1:L) = [2, 2, 2, 2] ! Multiplicity
            A(1:L) = [0.917117822312770D0, &
                      0.611268766465328D0, &
                      0.529422802042655D0, &
                      0.00000000000000000000D0]
            B(1:L) = [0.547931206828092D0, &
                      0.938843256658858D0, &
                      0.00000000000000000000D0, &
                      0.627041373780395D0]
            W(1:L) = [0.213057211620949D0, &
                      0.174009488946895D0, &
                      0.635853883443279D0, &
                      0.590012715421030D0]
            !-----------------------------------------------------------------------
         ELSEIF (D <= 8) THEN ! Ref [2]
            !-----------------------------------------------------------------------
            L = 16
            ALLOCATE (A(L), B(L), M(L), W(L))
            M(1:L) = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
            A(1:L) = [0.000000000000000D0, &
                      0.000000000000000D0, &
                      0.000000000000000D0, &
                      0.000000000000000D0, &
                      0.906662316102171D0, &
                      0.560348127669843D0, &
                      0.518695856299547D0, &
                      0.839267525316398D0, &
                      0.553832724273449D0, &
                      0.342287447682940D0, &
                      -0.906662316102171D0, &
                      -0.560348127669843D0, &
                      -0.518695856299547D0, &
                      -0.839267525316398D0, &
                      -0.553832724273449D0, &
                      -0.342287447682940D0]
            B(1:L) = [0.000000000000000D0, &
                      0.953321175521807D0, &
                      -0.882458098900126D0, &
                      0.582334188120499D0, &
                      0.294592444333741D0, &
                      -0.771253032094644D0, &
                      0.713923598834010D0, &
                      -0.272694549383947D0, &
                      0.179951160534772D0, &
                      -0.471118254595022D0, &
                      0.294592444333741D0, &
                      -0.771253032094644D0, &
                      0.713923598834010D0, &
                      -0.272694549383947D0, &
                      0.179951160534772D0, &
                      -0.471118254595022D0]
            W(1:L) = [0.334521439965580D0, &
                      0.087911387704639D0, &
                      0.167628039106469D0, &
                      0.305874815913735D0, &
                      0.087911387704639D0, &
                      0.087911387704639D0, &
                      0.167628039106469D0, &
                      0.167628039106469D0, &
                      0.305874815913735D0, &
                      0.305874815913735D0, &
                      0.087911387704639D0, &
                      0.087911387704639D0, &
                      0.167628039106469D0, &
                      0.167628039106469D0, &
                      0.305874815913735D0, &
                      0.305874815913735D0]
            !-----------------------------------------------------------------------
         ELSEIF (D <= 9) THEN ! Ref [5]*
            !-----------------------------------------------------------------------
            L = 5
            ALLOCATE (A(L), B(L), M(L), W(L))
            M(1:L) = [1, 4, 4, 4, 4]
            A(1:L) = [0.00000000000000000000d0, &
                      0.968849966361977D0, &
                      0.750277099978900D0, &
                      0.523735820214429D0, &
                      0.0762083281926171D0]
            B(1:L) = [0.00000000000000000000d0, &
                      0.630680119731668D0, &
                      0.927961645959569D0, &
                      0.453339821135647D0, &
                      0.852615729333662D0]
            W(1:L) = [0.526748971193415D0, &
                      0.0888793781701987D0, &
                      0.112099602129596D0, &
                      0.398282439262070D0, &
                      0.269051337639780D0]
            !-----------------------------------------------------------------------
         ELSEIF (D <= 10) THEN ! Ref []
            !-----------------------------------------------------------------------
            L = 22
            ALLOCATE (A(L), B(L), M(L), W(L))
            M(1:L) = 1
            A(1:L) = [4.73248988492765D-01, &
                      -3.50726726089189D-01, &
                      -4.71139214907016D-01, &
                      3.21100231203865D-02, &
                      1.07332278651087D-01, &
                      8.10374922601918D-01, &
                      -7.78825415983185D-01, &
                      -6.47635484262675D-01, &
                      3.92474875396096D-01, &
                      7.60506550713973D-01, &
                      -8.11715106016487D-01, &
                      -1.14295173642237D-01, &
                      7.52965632479960D-01, &
                      -6.24024379589847D-01, &
                      -2.06501346198872D-01, &
                      5.08913190429606D-01, &
                      -9.81711926404796D-01, &
                      -9.40618557199211D-01, &
                      -9.22548168257411D-01, &
                      9.53801922342551D-01, &
                      9.66342083687358D-01, &
                      9.57749591600075D-01]
            B(1:L) = [1.65578525100383D-01, &
                      1.84471720621219D-01, &
                      -6.66647330598211D-01, &
                      -3.18793575936407D-01, &
                      6.18866191392992D-01, &
                      6.11596783034925D-01, &
                      -2.10527389148215D-01, &
                      6.47494698175254D-01, &
                      -7.63111493924383D-01, &
                      -3.66313916780679D-01, &
                      -9.24668424290535D-01, &
                      -9.49219131408870D-01, &
                      -9.70718373967774D-01, &
                      9.85383311931460D-01, &
                      9.11958871035734D-01, &
                      9.21529075578982D-01, &
                      -6.25866193532396D-01, &
                      3.18845359683929D-01, &
                      8.79234804399032D-01, &
                      -7.55126920614355D-01, &
                      1.04312325566363D-01, &
                      9.26210500125838D-01]
            W(1:L) = [3.71717649308961D-01, &
                      3.88114474024408D-01, &
                      2.83958422182789D-01, &
                      4.08241977261545D-01, &
                      3.14988782212311D-01, &
                      2.00183206202775D-01, &
                      2.65871134771260D-01, &
                      2.30570045533700D-01, &
                      2.58279394103428D-01, &
                      2.60423976819168D-01, &
                      8.86862022169754D-02, &
                      1.18617672074659D-01, &
                      5.91101951503511D-02, &
                      3.33871292470730D-02, &
                      1.34609773861980D-01, &
                      1.33377311922401D-01, &
                      6.32692727611109D-02, &
                      1.19841585323912D-01, &
                      6.25379411875521D-02, &
                      7.18058987605166D-02, &
                      9.79404294841319D-02, &
                      3.44675255889838D-02]
            !-----------------------------------------------------------------------
         ELSEIF (D <= 11) THEN ! Ref [5]*
            !-----------------------------------------------------------------------
            L = 6
            ALLOCATE (A(L), B(L), M(L), W(L))
            M(1:L) = [4, 4, 4, 4, 4, 4]
            A(1:L) = [0.982639223540855D+00, &
                      0.825775835902963D+00, &
                      0.188586138718641D+00, &
                      0.812520548304813D+00, &
                      0.525320250364547D+00, &
                      0.416580719120223D-01]
            B(1:L) = [0.698076104549567D+00, &
                      0.939486382816736D+00, &
                      0.953539528201532D+00, &
                      0.315623432915254D+00, &
                      0.712001913075336D+00, &
                      0.424847248848669D+00]
            W(1:L) = [0.480207633507238D-01, &
                      0.660713291645505D-01, &
                      0.973867773586681D-01, &
                      0.211736349998948D+00, &
                      0.225626061728863D+00, &
                      0.351158718398245D+00]
            !-----------------------------------------------------------------------
         ELSEIF (D <= 12) THEN ! Ref [ ]
            !-----------------------------------------------------------------------
            L = 31
            ALLOCATE (A(L), B(L), M(L), W(L))
            M(1:L) = 1
            A(1:L) = [-3.81311194591487D-01, &
                      6.45731911094539D-02, &
                      -1.08343213291947D-02, &
                      -3.18827030208519D-01, &
                      4.91240514400928D-01, &
                      -2.09067994476707D-01, &
                      1.89483371892721D-01, &
                      -6.97400941145688D-01, &
                      -5.84025135895536D-01, &
                      4.47013772379773D-01, &
                      3.48679122083780D-01, &
                      -8.25626979146408D-01, &
                      8.28397946692489D-02, &
                      -4.03519746475846D-01, &
                      6.54010420564062D-01, &
                      -9.02643772854603D-01, &
                      -8.52953975519270D-01, &
                      8.22167862959016D-01, &
                      7.45105557704970D-01, &
                      -9.41339863623151D-01, &
                      -5.44611873876121D-01, &
                      7.44073229669089D-01, &
                      -7.22548854814926D-01, &
                      6.75244065693527D-01, &
                      9.38076769751955D-01, &
                      -9.91714087647107D-01, &
                      -9.73911493587488D-01, &
                      -9.22312952371136D-01, &
                      9.55577080195839D-01, &
                      9.82165912950434D-01, &
                      9.25572123351298D-01]
            B(1:L) = [-1.53980007841999D-01, &
                      -4.97948715686572D-01, &
                      1.22244789716577D-01, &
                      -7.32948988938657D-01, &
                      -7.84092576887272D-01, &
                      6.56110514972985D-01, &
                      8.99118055367382D-01, &
                      -4.22481270338075D-01, &
                      3.31006453787163D-01, &
                      -2.26596649736778D-01, &
                      4.38604049551637D-01, &
                      -7.83392232755809D-01, &
                      -9.42307986950857D-01, &
                      9.64828810963431D-01, &
                      7.38570980286750D-01, &
                      -1.58818522113793D-02, &
                      4.01632516022314D-01, &
                      -5.26607486922562D-01, &
                      1.09885498152123D-01, &
                      -9.54141481429644D-01, &
                      -9.55300363791224D-01, &
                      -9.77448354606736D-01, &
                      7.96622306496305D-01, &
                      9.96261104859096D-01, &
                      8.87292440266582D-01, &
                      -5.28539454774543D-01, &
                      5.99553405482386D-01, &
                      9.50825923437282D-01, &
                      -8.49033534185472D-01, &
                      -1.63697113637744D-01, &
                      4.74258804343308D-01]
            W(1:L) = [2.35987012929336D-01, &
                      2.41108136008478D-01, &
                      2.70482242276395D-01, &
                      1.86481668808841D-01, &
                      1.72419583938898D-01, &
                      2.41094420928215D-01, &
                      1.41253914532472D-01, &
                      1.82920254549737D-01, &
                      2.22229975818643D-01, &
                      2.49045377857700D-01, &
                      2.51123622294624D-01, &
                      9.28113599988029D-02, &
                      9.82383371730872D-02, &
                      6.85838177438892D-02, &
                      1.55069384444496D-01, &
                      1.31823208813660D-01, &
                      6.00685183552832D-02, &
                      1.57660185197732D-01, &
                      1.99308606026434D-01, &
                      2.17191093790923D-02, &
                      6.45960987588242D-02, &
                      4.05845199091051D-02, &
                      1.26033410079390D-01, &
                      3.00069124767560D-02, &
                      4.49747533292613D-02, &
                      3.87212445772882D-02, &
                      4.77180569223790D-02, &
                      3.00561896911893D-02, &
                      4.44992828465288D-02, &
                      5.65660932238902D-02, &
                      9.68147011095613D-02]
            !-----------------------------------------------------------------------
         ELSEIF (D <= 13) THEN ! Ref [5]
            !-----------------------------------------------------------------------
            L = 9
            ALLOCATE (A(L), B(L), M(L), W(L))
            M(1:L) = [1, 4, 4, 4, 4, 4, 4, 4, 4]
            A(1:L) = [0.00000000000000000000d+00, &
                      0.778809711554419D+00, &
                      0.957297699786307D+00, &
                      0.138183459862465D+00, &
                      0.941327225872925D+00, &
                      0.475808625218275D+00, &
                      0.755805356572081D+00, &
                      0.696250078491749D+00, &
                      0.342716556040406D+00]
            B(1:L) = [0.00000000000000000000d+00, &
                      0.983486682439872D+00, &
                      0.859556005641638D+00, &
                      0.958925170287534D+00, &
                      0.390736216129461D+00, &
                      0.850076673699748D+00, &
                      0.647821637187010D+00, &
                      0.707415089964449D-01, &
                      0.409304561694038D+00]
            W(1:L) = [0.300382115431225D+00, &
                      0.299918388644991D-01, &
                      0.381744213170836D-01, &
                      0.604249238177499D-01, &
                      0.774927385331053D-01, &
                      0.118844667300595D+00, &
                      0.129763550370002D+00, &
                      0.213341581457189D+00, &
                      0.256870749481967D+00]
            !-----------------------------------------------------------------------
         ELSEIF (D <= 15) THEN ! Ref [4]
            !-----------------------------------------------------------------------
            L = 11
            ALLOCATE (A(L), B(L), M(L), W(L))
            M(1:L) = [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4] ! Multipl.
            A(1:L) = [0.987984566507718D0, &
                      0.908159496006570D0, &
                      0.679283658334533D0, &
                      0.509113734117583D0, &
                      0.976753324669101D0, &
                      0.756199367191492D0, &
                      0.897785693286338D0, &
                      0.205993070742521D0, &
                      0.451443125112991D0, &
                      0.666838245383608D0, &
                      0.742957047557658D-1]
            B(1:L) = [0.771268212238755D0, &
                      0.957031834346906D0, &
                      0.882601975930872D0, &
                      0.971203129741836D0, &
                      0.835598626087816D-1, &
                      0.756199367191492D0, &
                      0.466762659237964D0, &
                      0.840794484540785D0, &
                      0.562456862332199D0, &
                      0.190466302435717D0, &
                      0.323977022497530D0]
            W(1:L) = [0.208814702044975D-1, &
                      0.255459015744972D-1, &
                      0.312038666249333D-1, &
                      0.380107615950748D-1, &
                      0.414490618524261D-1, &
                      0.793204070040833D-1, &
                      0.889012657587515D-1, &
                      0.120169821582060D0, &
                      0.168820434106397D0, &
                      0.169871624973361D0, &
                      0.215825384723915D0]
            !-----------------------------------------------------------------------
         ELSEIF (D <= 17) THEN ! Ref [4]
            !-----------------------------------------------------------------------
            L = 14
            ALLOCATE (A(L), B(L), M(L), W(L))
            M(1:L) = [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4] ! Multipl.
            A(1:L) = [0.833959144220507D0, &
                      0.967012407603778D0, &
                      0.986514410860330D0, &
                      0.170350608089954D-1, &
                      0.425238965224530D0, &
                      0.881941090892153D0, &
                      0.699789077190586D0, &
                      0.908584129588383D0, &
                      0.739337592052920D0, &
                      0.212461558378854D0, &
                      0.156441720958463D0, &
                      0.503985638194279D0, &
                      0.582438950744672D0, &
                      0.289703230655412D0]
            B(1:L) = [0.996901349982582D0, &
                      0.921015650153696D0, &
                      0.562517802446672D0, &
                      0.982049142568430D0, &
                      0.957964532368651D0, &
                      0.741176885697325D0, &
                      0.887280242932557D0, &
                      0.258431221518207D0, &
                      0.472955832572976D0, &
                      0.00000000000000000000000000000000D0, &
                      0.812592125239123D0, &
                      0.682010932977925D0, &
                      0.118465445606478D0, &
                      0.408769859537943D0]
            W(1:L) = [0.106934839869745D-1, &
                      0.167716229893254D-1, &
                      0.215208348031730D-1, &
                      0.248932015326650D-1, &
                      0.424632584720309D-1, &
                      0.537112650376450D-1, &
                      0.545794796933824D-1, &
                      0.673756536224613D-1, &
                      0.980252828851022D-1, &
                      0.983256515846666D-1, &
                      0.105768983196657D0, &
                      0.105934535745754D0, &
                      0.149904054849169D0, &
                      0.150032691600992D0]
            !-----------------------------------------------------------------------
         ELSEIF (D <= 19) THEN ! Ref [4]
            !-----------------------------------------------------------------------
            L = 17
            ALLOCATE (A(L), B(L), M(L), W(L))
            M(1:L) = [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4] ! Multipl.
            A(1:L) = [0.937909440601746D0, &
                      0.597879665191571D0, &
                      0.968657817477988D0, &
                      0.988717132764473D0, &
                      0.983985341323383D0, &
                      0.365857759345555D0, &
                      0.783921491760966D0, &
                      0.869371448989578D0, &
                      0.918893777775738D0, &
                      0.140783968044568D0, &
                      0.563796668154465D0, &
                      0.762002742930703D0, &
                      0.697453883421912D0, &
                      0.523103923394044D0, &
                      0.305668369039291D0, &
                      0.448986426282880D0, &
                      0.982534187598351D-1]
            B(1:L) = [0.999985460728529D0, &
                      0.987328479417815D0, &
                      0.898374951635325D0, &
                      0.622596343895302D0, &
                      0.811094634878249D-1, &
                      0.964624004228919D0, &
                      0.937564445443781D0, &
                      0.743808660345971D0, &
                      0.384108758227378D0, &
                      0.879879573079810D0, &
                      0.818792106076360D0, &
                      0.157499201227322D0, &
                      0.546575274601777D0, &
                      0.00000000000000000000000000000000D0, &
                      0.650176702706879D0, &
                      0.342404653806802D0, &
                      0.235096216291155D0]
            W(1:L) = [0.421571893124572D-2, &
                      0.992376010147412D-2, &
                      0.150786788795815D-1, &
                      0.151214968648229D-1, &
                      0.238216210473395D-1, &
                      0.287464372521896D-1, &
                      0.313487155038614D-1, &
                      0.497628526667171D-1, &
                      0.555347751596040D-1, &
                      0.650137104321739D-1, &
                      0.738190689007317D-1, &
                      0.903858259681506D-1, &
                      0.918818333364250D-1, &
                      0.949671426388560D-1, &
                      0.104981747241028D0, &
                      0.118713368509280D0, &
                      0.126683246566517D0]
            !-----------------------------------------------------------------------
         ELSEIF (D <= 21) THEN ! Ref [4]
            !-----------------------------------------------------------------------
            L = 21
            ALLOCATE (A(L), B(L), M(L), W(L))
            M(1:L) = [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 1]
            A(1:L) = [0.997428443180714D0, &
                      0.900922057228577D0, &
                      0.981375816611523D0, &
                      0.475623988460619D0, &
                      0.683962178755243D-1, &
                      0.964451711702909D0, &
                      0.967895340677605D0, &
                      0.712046261839021D0, &
                      0.866865035635049D0, &
                      0.293888957654902D0, &
                      0.882452897205096D0, &
                      0.883564655497776D0, &
                      0.537143303287965D0, &
                      0.732977372786888D0, &
                      0.743408158979943D0, &
                      0.143523798672578D0, &
                      0.577163748064870D0, &
                      0.380106531055197D0, &
                      0.493803047507045D0, &
                      0.233006942769649D0, &
                      0.00000000000000000000000000000000D0]
            B(1:L) = [0.523493335403422D0, &
                      0.992136246111987D0, &
                      0.935777665005192D0, &
                      0.988065039814063D0, &
                      0.982891421213467D0, &
                      0.730795641927920D0, &
                      0.284169437040225D0, &
                      0.950017375773952D0, &
                      0.857545149044262D0, &
                      0.914783753744201D0, &
                      0.783105713193479D-1, &
                      0.506500278197455D0, &
                      0.821116623215359D0, &
                      0.678114741579903D0, &
                      0.258912124059179D0, &
                      0.737142844448782D0, &
                      0.426051378322523D0, &
                      0.591909801960054D0, &
                      0.469684383899153D-1, &
                      0.278540549928700D0, &
                      0.00000000000000000000000000000000D0]
            W(1:L) = [0.592452899102747D-2, &
                      0.631818799935309D-2, &
                      0.776543567718858D-2, &
                      0.133132365246493D-1, &
                      0.180280236320003D-1, &
                      0.207592110848114D-1, &
                      0.243455946869629D-1, &
                      0.272918298131574D-1, &
                      0.313137244893206D-1, &
                      0.430286844636254D-1, &
                      0.474351074977605D-1, &
                      0.495352327910999D-1, &
                      0.590866097375227D-1, &
                      0.615988987555731D-1, &
                      0.686628827371055D-1, &
                      0.775110257604497D-1, &
                      0.776559008613981D-1, &
                      0.841651591432533D-1, &
                      0.116643957423565D0, &
                      0.125870779667014D0, &
                      0.134983953052639D0]
            !-----------------------------------------------------------------------
         ELSEIF (D <= 23) THEN ! Ref [4]
            !-----------------------------------------------------------------------
            L = 25
            ALLOCATE (A(L), B(L), M(L), W(L))
            M(1:L) = [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]
            A(1:L) = [0.224757764352695D0, &
                      0.874207930876896D0, &
                      0.976022672613640D0, &
                      0.989010257580283D0, &
                      0.993269209240316D0, &
                      0.564995923166212D0, &
                      0.124799332346678D0, &
                      0.902404399320343D0, &
                      0.966034160999610D0, &
                      0.743274033812573D0, &
                      0.943971582343373D0, &
                      0.364735933238215D0, &
                      0.900253092874432D0, &
                      0.480346658395002D-1, &
                      0.549191162143093D0, &
                      0.864875376908999D0, &
                      0.787696213093480D0, &
                      0.576710284209456D0, &
                      0.200134082627620D0, &
                      0.746027287830323D0, &
                      0.641772804586726D0, &
                      0.402364113097523D0, &
                      0.469939831200515D0, &
                      0.224917494381230D0, &
                      0.228317533864552D0]
            B(1:L) = [0.295629261234400D0, &
                      0.991644841702323D0, &
                      0.960386443899881D0, &
                      0.794761706954390D0, &
                      0.452699068266705D0, &
                      0.989634439635466D0, &
                      0.988916345115731D0, &
                      0.866817471658503D0, &
                      0.221288379835874D0, &
                      0.937392797891401D0, &
                      0.627343738112143D0, &
                      0.936540535531906D0, &
                      0.250709903414279D-1, &
                      0.653046876919900D0, &
                      0.00000000000000000000000000000000D0, &
                      0.420217007311403D0, &
                      0.714861371796371D0, &
                      0.826496431062503D0, &
                      0.808242673492922D0, &
                      0.205819811656463D0, &
                      0.515166376877061D0, &
                      0.648907818198540D0, &
                      0.277147991514298D0, &
                      0.403755035172687D0, &
                      0.561167909823553D-1]
            W(1:L) = [-0.224991445901807D-1, &
                      0.567577972797097D-2, &
                      0.612945416323857D-2, &
                      0.773073997169215D-2, &
                      0.798510306890593D-2, &
                      0.107477249090711D-1, &
                      0.140737229019146D-1, &
                      0.229501666562158D-1, &
                      0.231045036685527D-1, &
                      0.248865834606209D-1, &
                      0.254956087107813D-1, &
                      0.351487914573924D-1, &
                      0.399676124712756D-1, &
                      0.422109061992908D-1, &
                      0.456293089176021D-1, &
                      0.471849315439838D-1, &
                      0.480043605117840D-1, &
                      0.514071322819658D-1, &
                      0.542104950316212D-1, &
                      0.639802264708508D-1, &
                      0.719998547139675D-1, &
                      0.741303684244850D-1, &
                      0.870755978570734D-1, &
                      0.100800083238107D0, &
                      0.111970088231815D0]
            !-----------------------------------------------------------------------
         ELSE
            !-----------------------------------------------------------------------
            PRINT *, '  ******** ERROR!!! D must be 1 <= D <= 23 ********  '
            PRINT *, '  ********       for REGION = SQUARE       ********  '
            PRINT *, '    Execution terminated in subroutine quadrature    '
            STOP
         END IF

         !.......Count up number of quadrature points

         N = 0
         DO I = 1, L
            N = N + M(I)
         END DO

         !        ALLOCATE ( PTS(N,2), WTS(N) )
         !.......Transform quadrature points to master element coordinates

         J = 1
         DO I = 1, L
            IF (M(I) == 1) THEN
               PTS(J, 1) = A(I)
               PTS(J, 2) = B(I)
               WTS(J) = W(I)
               J = J + M(I)
            ELSEIF (M(I) == 2) THEN
               PTS(J, 1) = A(I)
               PTS(J, 2) = B(I)
               WTS(J) = W(I)
               PTS(J + 1, 1) = -A(I)
               PTS(J + 1, 2) = -B(I)
               WTS(J + 1) = W(I)
               J = J + M(I)
            ELSEIF (M(I) == 4) THEN
               PTS(J, 1) = A(I)
               PTS(J, 2) = B(I)
               WTS(J) = W(I)
               PTS(J + 1, 1) = -B(I)
               PTS(J + 1, 2) = A(I)
               WTS(J + 1) = W(I)
               PTS(J + 2, 1) = -A(I)
               PTS(J + 2, 2) = -B(I)
               WTS(J + 2) = W(I)
               PTS(J + 3, 1) = B(I)
               PTS(J + 3, 2) = -A(I)
               WTS(J + 3) = W(I)
               J = J + M(I)
            END IF
         END DO

         !-----------------------------------------------------------------------
         !                        TRIANGULAR PRISM QUADRATURE RULES
         !-----------------------------------------------------------------------
         !
         !     All rules have positive weights with all points located inside the
         !     triangle (so-called PI rules). Rules marked by * are optimal.  The
         !     rest are the best PI rules currently known.
         !
         !     [1] E.J. Kubatko, "Efficient cubature rules for triangular prisms"
         !         in preparation.
         !
         !-----------------------------------------------------------------------

      ELSEIF (REGION == 'TRIPRISM') THEN

         !-----------------------------------------------------------------------
         IF (D <= 1) THEN ! Ref [1]*
            !-----------------------------------------------------------------------
            L = 1
            ALLOCATE (A(L), B(L), C(L), M(L), W(L))
            M(1:L) = [1]
            A(1:L) = [-0.333333333333333D0]
            B(1:L) = [-0.333333333333333D0]
            C(1:L) = [0.0000000000000000000D0]
            W(1:L) = [4.00000000000000D0]
            !-----------------------------------------------------------------------
         ELSEIF (D <= 2) THEN ! Ref [1]*
            !-----------------------------------------------------------------------
            L = 4
            ALLOCATE (A(L), B(L), C(L), M(L), W(L))
            M(1:L) = [1, 1, 1, 1]
            A(1:L) = [0.24401693585629D0, &
                      -0.91068360252296D0, &
                      -0.33333333333333D0, &
                      -0.33333333333333D0]
            B(1:L) = [-0.33333333333333D0, &
                      0.24401693585629D0, &
                      -0.91068360252296D0, &
                      -0.91068360252296D0]
            C(1:L) = [0.00000000000000D0, &
                      0.00000000000000D0, &
                      1.00000000000000D0, &
                      -1.00000000000000D0]
            W(1:L) = [1.33333333333333D0, &
                      1.33333333333333D0, &
                      0.66666666666667D0, &
                      0.66666666666667D0]
            !-----------------------------------------------------------------------
         ELSEIF (D <= 3) THEN ! Ref [1]
            !-----------------------------------------------------------------------
            !-----------------------------------------------------------------------
         ELSE
            !-----------------------------------------------------------------------
            PRINT *, '  ******** ERROR!!! D must be 1 <= D <= 4  ********  '
            PRINT *, '  ********       for REGION = SQUARE       ********  '
            PRINT *, '    Execution terminated in subroutine quadrature    '
            STOP
         END IF

         !.......Count up number of quadrature points

         N = 0
         DO I = 1, L
            N = N + M(I)
         END DO
         !        ALLOCATE ( PTS(N,3), WTS(N) )

         !.......Transform quadrature points to master element coordinates

         J = 1
         DO I = 1, N
            IF (M(I) == 1) THEN
               PTS(J, 1) = A(I)
               PTS(J, 2) = B(I)
               PTS(J, 3) = C(I)
               WTS(J) = W(I)
               J = J + M(I)
            END IF
         END DO
         !-----------------------------------------------------------------------
      ELSE
         !-----------------------------------------------------------------------
         PRINT *, '  ********  ERROR!!! REGION must be EDGE,  ********  '
         PRINT *, '  ********  TRIANGLE, SQUARE OR TRIPRISM   ********  '
         PRINT *, '    Execution terminated in subroutine quadrature    '
         STOP
      END IF
      RETURN
   END SUBROUTINE QUADRATURE

   !
   !     Subroutine quad_rules_general
   !
   !     Written by Dylan Wood (03-10-2017)
   !
   !     This subroutine returns the points of a n point quadrature rule
   !     for the 1-dimensional domain -1 to 1, via Golub-Welsch algorithm.
   !     Legendre weights are determined from algebraic formulas.
   !     Jacobi weights are computed from Gaussian elimination.
   !
   !     Quadrature type/weight function is Jacobi and determined from
   !     input values of alpha and beta.
   !     A logical input indicates whether Lobatto rules are desired.
   !
   !-----------------------------------------------------------------------
   !
   !     Input:
   !     ------
   !       n:       Number of quadrature points/weights desired.
   !       lobatyn: Logicial value indicating whether a Lobatto type rule
   !                 is desired, i.e. rule including the endpoints -1 and 1
   !                 as quadrature points.
   !       alph:    Value for alpha in Jacobi weight,
   !                  (1 - x)**alpha*(1 + x)**beta.
   !       bet:     Value for beta in Jacobi weight.
   !
   !       If alph = bet = 0 and lobatyn = .false.,
   !          then a Gauss-Legendre type rule is generated.
   !       If either alph or bet does not equal 0, and lobatyn = .false.,
   !          then a Gauss-Jacobi type rule is generated.
   !       If lobatyn = .true.,
   !          then alph and bet MUST BOTH equal 0 for a quadrature rule to
   !          be generated by this subroutine. In this case the rule is of
   !          type Gauss-Legendre-Lobatto.
   !       Gauss-Jacobi-Lobatto rules are not currently generated by this
   !        subroutine.
   !
   !     Output:
   !     -------
   !       pts: Array of quadrature points of size n by 1
   !       wts: Array of quadrature weights of size n by 1.
   !
   !-----------------------------------------------------------------------
   subroutine quad_rules_general(n, lobatyn, alph, bet, pts, wts)

      implicit none

      integer, intent(inout) :: n
      integer :: nin
      logical, intent(in) :: lobatyn
      integer :: j, m, quadtype, nn, p, ierr, i
      real(8), intent(in) :: alph, bet
      real(8), dimension(n), intent(out) :: pts, wts
      real(8), dimension(n) :: aj, bj
      real(8) :: Jac(n, n), augJ(n + 1, n), eig(n), v(n)
      real(8), dimension(0:n + 1) :: polyo
      real(8) :: eps, jrl, dum, ajp1, bjp1, x, prl, nv, lambda
      real(8) :: summ, al, be, nrl, B0

      !.....Determine whether Lobatto rules are needed. If so,
      !.....prepare elements for modified Jacobi matrix, and
      !.....places these elements in the necessary locations.
      if (lobatyn) then
         al = alph
         be = bet
         if (.not. (abs(alph) <= 1d-15 .and. abs(bet) <= 1d-15)) then
            PRINT *, ' ERROR ENCOUNTERED IN quad_rules_general.f! '
            PRINT *, ' INVALID CHOICE(S) FOR ALPHA/BETA FOR QUADRATURE! '
            PRINT *, ' IF LOBATTO RULES ARE NEEDED, CURRENTLY ONLY '
            PRINT *, ' LEGENDRE WEIGHTS ARE SUPPORTED FOR THESE RULES '
            PRINT *, ' NOW SETTING ALPHA = BETA = 0. '
            al = 0.0d0
            be = 0.0d0
            PRINT *, ' CONTINUING EXECUTION. '
         end if

         n = n - 1
         nn = n - 1
         nrl = real(nn, kind=8)
         Jac(:, :) = 0.0d0
         ajp1 = (al - be)/(2.0d0*nrl + al + be + 2.0d0)
         bjp1 = 4.0d0*(nrl + al + 1.0d0)*(nrl + be + 1.0d0) &
                *(nrl + al + be + 1.0d0) &
                /((2.0d0*nrl + al + be + 1.0d0) &
                  *(2.0d0*nrl + al + be + 2.0d0) &
                  *(2.0d0*nrl + al + be + 2.0d0))
         bjp1 = sqrt(bjp1)
         Jac(n + 1, n + 1) = ajp1
         Jac(n, n + 1) = bjp1
         Jac(n + 1, n) = bjp1
         nn = n + 1

      else
         al = alph
         be = bet
         Jac(:, :) = 0.0d0
         nn = n
      end if

      !.....Determine whether Legendre or Jacobi weights are desired and
      !.....calculate entries of Jacobi matrix accordingly.
      if (abs(al) <= 1d-15 .and. abs(be) <= 1d-15) then
         quadtype = 1
         aj(1:n) = 0.0d0
         do j = 1, n - 1
            jrl = real(j, kind=8)
            bj(j) = jrl*jrl/(4.0d0*jrl*jrl - 1.0d0)
         end do
      elseif (al > -1.0d0 .and. be > -1.0d0) then
         quadtype = 2
         do j = 1, n
            jrl = real(j, kind=8)
            dum = 2.0d0*jrl + al + be
            aj(j) = (be - al)*(be + al)/(dum*(dum - 2.0d0))
         end do
         do j = 1, n - 1
            jrl = real(j, kind=8)
            dum = 2.0d0*jrl + al + be
            bj(j) = 4.0d0*jrl*(al + jrl) &
                    *(be + jrl)*(al + be + jrl) &
                    /(dum*dum*(dum + 1.0d0)*(dum - 1.0d0))
         end do
      else
         PRINT *, ' ERROR ENCOUNTERED IN quad_rules_general.f! '
         PRINT *, ' INVALID CHOICE(S) FOR ALPHA/BETA FOR QUADRATURE! '
         PRINT *, ' ALPHA AND BETA MUST BOTH BE GREATER THAN -1! '
         PRINT *, ' EXECUTION TERMINATED IN SUBROUTINE &
  &            quad_rules_general. '
         STOP
      end if
      bj(n) = 0.0d0
      bj(1:n) = sqrt(bj(1:n))

      !.....Construct Jacobi matrix
      do j = 1, n
         Jac(j, j) = aj(j)
      end do
      do j = 2, n
         Jac(j, j - 1) = bj(j - 1)
         Jac(j - 1, j) = bj(j - 1)
      end do

      !.....Determine machine precision
      eps = epsilon(aj(1))

      !.....Find eigenvalues of Jacobi matrix via QR iteration
      call QRdecomp(Jac, nn, eps, eig)

      !.....If rule is symmetric, duplicate values across the symmetry.
      if (abs(al - be) > 1d-14) then
         if (mod(nn, 2) == 1) then
            m = ceiling(real(nn, kind=8)/2.d0)
            pts(1:m) = eig(1:m)
            pts(m + 1:nn) = -pts(m - 1:1:-1)
         else
            m = nn/2
            pts(1:m) = eig(1:m)
            pts(m + 1:nn) = -pts(m:1:-1)
         end if
      else
         pts(:) = eig(:)
      end if

      !.....If a point is very close to 0, then set it equal to 0.
      do j = 1, nn
         if (abs(pts(j)) < 1.0d-12) pts(j) = 0.0d0
      end do
      !.....If Lobatto endpoints are close to 1 in absolute value, then
      !.....set endpoints equal to 1.
      if (lobatyn) then
         if (abs(pts(1) - (-1.0d0)) < 1.0d-12) pts(1) = -1.0d0
         if (abs(pts(nn) - 1.0d0) < 1.0d-12) pts(nn) = 1.0d0
      end if

      !.....Find weights from algebraic formulas
      do m = 1, nn
         x = pts(m)
         if (quadtype == 1 .and. .not. lobatyn) then
            !...........Legendre weights
            polyo(0) = 1.0d0
            polyo(1) = x
            !...........Compute values of Legendre polys. at current point from
            !...........recurrence relation.
            do p = 1, n
               prl = real(p, kind=8)
               polyo(p + 1) = ((2.0d0*prl + 1.0d0)*x*polyo(p) - prl &
                               *polyo(p - 1))/(prl + 1.0d0)
            end do
            !...........Compute weights
            nrl = real(n, kind=8)
            wts(m) = 2.0d0*(1.0d0 - x*x)/((nrl + 1.0d0)*(nrl + 1.0d0) &
                                          *polyo(n + 1)*polyo(n + 1))
         elseif (quadtype == 2 .and. .not. lobatyn) then
            !...........Jacobi weights, computed from eigenvectors, which are
            !...........computed by Guassian elimination.
            lambda = pts(m)

            do j = 1, nn
               do i = 1, nn
                  augJ(i, j) = Jac(i, j)
               end do
               augJ(j, j) = augJ(j, j) - lambda
               augJ(nn + 1, j) = 0.0d0
            end do

            call Gauss_Elim(augJ, nn, ierr)

            v = augJ(nn + 1, :)

            summ = 0.0d0
            do j = 1, nn
               summ = summ + v(j)*v(j)
            end do
            nv = sqrt(summ)

            if (abs(nv) < eps) nv = eps
            v = v/nv

            B0 = 2.0d0**(alph + bet + 1.0d0) &
                 *gammaf(alph + 1.0d0)*gammaf(bet + 1.0d0) &
                 /gammaf(alph + bet + 2.0d0)
            wts(m) = B0*v(1)*v(1)

         elseif (lobatyn) then
            if (quadtype == 1) then
               !..............Lobatto rule weights with Lengendre weight function.
               nin = nn - 2
               nrl = real(nin, kind=8)
               if (m == 1) then
                  !..............Compute weights at endpoints.
                  wts(m) = 2.0d0/((nrl + 1.0d0)*(nrl + 2.0d0))
               elseif (m == nn) then
                  wts(m) = wts(1)
               else
                  !..............Compute weights at interior quad points.
                  polyo(0) = 1.0d0
                  polyo(1) = x
                  !.................Compute values of Legendre polys. at current point
                  !.................from recurrence relation.
                  do p = 1, nin
                     prl = real(p, kind=8)
                     polyo(p + 1) = ((2.0d0*prl + 1.0d0)*x*polyo(p) - prl &
                                     *polyo(p - 1))/(prl + 1.0d0)
                  end do
                  !.................Compute weights
                  wts(m) = 2.0d0/((nrl + 1.0d0)*(nrl + 2.d0) &
                                  *polyo(nin + 1)*polyo(nin + 1))
               end if
            elseif (quadtype == 2) then
               !..............This is where code for weights of Guass-Jacobi-Lobatto
               !..............type rules would go. The author has no practical need
               !..............for these rules, so they are currently excluded.
            end if
         end if

      end do
      if (lobatyn) then
         n = n + 1
      end if

   end subroutine quad_rules_general

   subroutine QRdecomp(A, n, eps, eig)
      !.....QR decomposition with shifting
      implicit none
      integer, intent(in)::n
      integer :: i, k
      real(8), dimension(n, n), intent(in):: A
      real(8), dimension(n, n) :: A2, oldA2, Q, QT, R, EYE
      real(8) :: res, sig, normA
      real(8), intent(in):: eps
      real(8), dimension(n) :: vec1, vec2
      real(8), dimension(n), intent(out)::eig

      EYE(:, :) = 0.0d0
      do i = 1, n
         EYE(i, i) = 1.0d0
      end do

      oldA2(:, :) = 1.0d0
      A2(:, :) = A(:, :)

      do while (DoIterate(A2, oldA2, eps, n))
         sig = 0.0d0
         do i = 1, n
            normA = 0.0d0
            do k = 1, n
               normA = normA + A2(i, k)*A2(i, k)
            end do
            normA = sqrt(normA)
            sig = max(normA, sig)
         end do

         oldA2(:, :) = A2(:, :)
         A2(:, :) = A2(:, :) - sig*EYE(:, :)

         Q(:, :) = A2(:, :)
         do i = 1, n
            do k = 1, i - 1
               vec1 = Q(:, i)
               vec2 = Q(:, k)
               call dotprod(vec1, vec2, n, res)
               Q(:, i) = Q(:, i) - res*Q(:, k)
            end do
            call dotprod(Q(:, i), Q(:, i), n, res)
            if (abs(res) < eps) res = eps
            Q(:, i) = Q(:, i)/sqrt(res)
         end do
         QT = transpose(Q)
         R = matmul(QT, A2)
         A2 = matmul(R, Q) + sig*EYE(:, :)
      end do

      do i = 1, n
         eig(i) = A2(i, i)
      end do
   end subroutine QRdecomp

   function gammaf(n)
      ! Compute the gamma function for an input
      implicit none
      real(8), intent(in) :: n
      integer :: i, m
      real(8) :: ans, gammaf

      ! Subtract 1 from input number (should be an integer as a real)
      m = int(n - 1.0d0)
      ans = 1.0d0
      ! Calculate the factorial of n-1
      do i = 1, m
         ans = ans*dble(i)
      end do
      gammaf = ans
      return

   end function gammaf

   subroutine dotprod(vec1, vec2, n, prod)
      ! Compute the dot product of two vector inputs
      implicit none
      integer, intent(in) ::n
      integer :: i
      real(8), intent(in), dimension(n)::vec1, vec2
      real(8), intent(out) ::prod
      prod = 0.0d0
      do i = 1, n
         prod = prod + vec1(i)*vec2(i)
      end do

   end subroutine dotprod

   pure function DoIterate(A1, A2, eps, n)
      ! Evaluate whether or not to iterate QR decomp.
      implicit none
      integer, intent(in) ::n
      integer :: i, j
      real(8), intent(in), dimension(n, n)::A1, A2
      real(8), intent(in) ::eps
      logical::DoIterate

      DoIterate = .false.
      do i = 1, n
         do j = 1, n
            if (abs(abs(A1(i, j)) - abs(A2(i, j))) > eps) then
               DoIterate = .true.
               exit
            end if
         end do
      end do
      return

   end function DoIterate

   subroutine Gauss_Elim(M, n, ierr)
      ! Perform Gaussian elimination on an "augmented matrix"
      implicit none
      integer, intent(in)::n
      integer, intent(out)::ierr
      real(8), dimension(n + 1, n), intent(inout)::M
      real(8)::tempM, val, eps
      integer::i, j, k, c, a, b, temp
      logical::broken

      eps = 1.0d-8
      broken = .false.
      ierr = 0
      do i = 1, n
         temp = i
         val = abs(M(i, i))
         do k = i + 1, n
            if (val < abs(M(i, k))) then
               val = abs(M(i, k))
               temp = k
            end if
         end do
         if (temp /= i) then
            do k = i, n + 1
               tempM = M(k, i)
               M(k, i) = M(k, temp)
               M(k, temp) = tempM
            end do
         end if
         val = M(i, i)
         if (abs(val) < eps) then
            broken = .true.
            exit
         end if
         do k = i, n + 1
            M(k, i) = M(k, i)/val
         end do
         do c = i + 1, n
            val = M(i, c)
            do k = i, n + 1
               M(k, c) = M(k, c) - val*M(k, i)
            end do
         end do
      end do

      if (broken) then
         ierr = -i
         do a = i, n + 1
            do b = i, n
               if (abs(M(a, b)) > eps) then
                  ierr = a
                  return
               end if
               M(a, b) = 0.0d0
            end do
         end do
         do b = i, n
            M(b, b) = 1.0d0
            M(n + 1, b) = real(n - b + 1, kind=8)
         end do
      end if

      do i = n, 1, -1
         do j = i - 1, 1, -1
            M(n + 1, j) = M(n + 1, j) - M(i, j)*M(n + 1, i)
            M(i, j) = M(i, j) - M(i, j)*M(i, i)
         end do
      end do

   end subroutine Gauss_Elim

   !
   !     Subroutine ORTHOGONAL_BASIS
   !
   !     Written by Ethan Kubatko (11-12-2009)
   !
   !     This subroutine evaluates the orthogonal basis functions and its
   !     derivatives at a given point for a 2D triangle, a 3D
   !     triangular prism, a 2D quadrilateral or a 3D
   !     hexahedral.  The recurrence relations used in the 2D triangle
   !     are from [1].
   !
   !-----------------------------------------------------------------------
   !
   !     Input:
   !     ------
   !       ELEM:    Element type (1 for Triangle 2 for Quad)
   !       PT:      Coordinate point (of dimension 2 or 3)
   !       D:       Polynomial degree of basis
   !       DIM:     Dimension (2 or 3)
   !
   !     Output:
   !     -------
   !       BASIS:  Array of length (D+2)*(D+1)**(DIM-1)/2 containing the
   !               values of the basis functions evaluated at PT.
   !       DBASIS: Array of size (D+2)*(D+1)**(DIM-1)/2 by DIM containing
   !               the values of the gradients of the basis functions
   !               evaluated at PT.
   !
   !-----------------------------------------------------------------------
   !
   !     Reference: [1] Robert C. Kirby, "Singularity-free evaluation of
   !                    collapsed-coordinate orthogonal polynomials", ACM
   !                    Transactions on Mathematical Software, Volume 37
   !                    Issue 1, January 2010, Article No. 5, Pages 5:1--
   !                    5:16.
   !
   !                    http://doi.acm.org/10.1145/1644001.1644006
   !
   !-----------------------------------------------------------------------
   !
   !     Updates:
   !
   !     Quadrilateral and hexahedral element bases added
   !       - Ashley Maggi, 07-13-2010
   !
   !-----------------------------------------------------------------------

   SUBROUTINE ORTHOGONAL_BASIS(ELEM, PT, D, DIM, BASIS, DBASIS)

      IMPLICIT NONE

      !.....Declare subroutine input and output

      INTEGER, INTENT(IN)  :: D, DIM, ELEM
      REAL(8), INTENT(IN)  :: PT(DIM)
      REAL(8), INTENT(OUT) :: BASIS((D + 2)*(D + 1)**(DIM - 1)/2)
      REAL(8), INTENT(OUT) :: DBASIS((D + 2)*(D + 1)**(DIM - 1)/2, DIM)

      !.....Explicitly declare remaining variables

      INTEGER :: P, Q, R
      REAL(8) :: AQ, BQ, CQ, X, Y, Z, preal, qreal
      REAL(8), DIMENSION(0:D) :: LEGENDRE
      REAL(8), DIMENSION(0:D) :: LEGENDRE_X, LEGENDRE_Y
      REAL(8), DIMENSION(0:D) :: DLEGENDRE
      REAL(8), DIMENSION(0:D) :: DLEGENDRE_X, DLEGENDRE_Y
      REAL(8), DIMENSION(0:D, 0:D) :: PHI
      REAL(8), DIMENSION(0:D, 0:D, 2) :: DPHI
      REAL(8), DIMENSION((D + 2)*(D + 1)/2) :: PHI_2D
      REAL(8), DIMENSION((D + 2)*(D + 1)/2, 2) :: DPHI_2D

      !.....Check to make sure input is correct

      IF (D < 0) THEN
         PRINT *, '  ****** ERROR!!! D must be a positve integer *******  '
         PRINT *, '  Execution terminated in subroutine orthogonal_basis  '
         STOP
      ELSEIF ((DIM > 3) .OR. (DIM < 2)) THEN
         PRINT *, '  *********** ERROR!!! DIM must be 2 or 3 ***********  '
         PRINT *, '  Execution terminated in subroutine orthogonal_basis  '
         STOP
      ELSEIF ((ELEM /= 1) .AND. (ELEM /= 2)) THEN
         PRINT *, '  **** ERROR!!! ELEM must be TRIA = 1 or QUAD = 2 ****  '
         PRINT *, '  Execution terminated in subroutine orthogonal_basis  '
         STOP
      END IF

      !.....Assign X, Y, and Z coordinates

      X = PT(1)
      Y = PT(2)
      IF (DIM == 3) THEN
         Z = PT(3)
      ELSE
         Z = 0.D0
      END IF

      !.....Case D = 0 for 2D or 3D

      IF (D == 0) THEN
         BASIS = 1.D0
         DBASIS = 0.D0
         RETURN
      END IF

      !-----------------------------------------------------------------------
      !
      !     2D Triangular basis
      !
      !-----------------------------------------------------------------------
      IF (ELEM == 1) THEN

         !.......Check to make sure point PT is contained within the reference
         !.......element.  Issue warning if not but continue execution.
#if 0
         IF ((ABS(X) <= 1.D0) .AND. (ABS(Y) <= 1.D0) .AND. &
             ((X + Y) <= 0.D0) .AND. (ABS(Z) <= 1.D0)) THEN
            CONTINUE
         ELSE
            PRINT *, ' ***** WARNING!! in subroutine orthogonal_basis ***** '
            PRINT *, '       Point is outside of reference element!!        '
            PRINT *, '               Execution will continue                '
         END IF
#endif

         !.......First basis function and derivatives

         PHI(0, 0) = 1.D0
         DPHI(0, 0, 1) = 0.D0
         DPHI(0, 0, 2) = 0.D0

         !.......Second basis function and derivatives

         PHI(1, 0) = (1.D0 + 2.D0*X + Y)/2.D0
         DPHI(1, 0, 1) = 1.D0
         DPHI(1, 0, 2) = 1.D0/2.D0

         !.......Recurrence relations

         DO P = 1, D - 1
            preal = real(P,8)
            PHI(P + 1, 0) = (2.D0*Preal + 1.D0)/(Preal + 1.D0)*PHI(1, 0)*PHI(P, 0) &
                            - Preal/(Preal + 1.D0)*((1.D0 - Y)/2.D0)**2.D0*PHI(P - 1, 0)
            DPHI(P + 1, 0, 1) = (2.D0*Preal + 1.D0)/(Preal + 1.D0)*(DPHI(1, 0, 1)*PHI(P, 0) &
                                                            + PHI(1, 0)*DPHI(P, 0, 1)) &
                                - Preal/(Preal + 1.D0)*((1.D0 - Y)/2.D0)**2.D0*DPHI(P - 1, 0, 1)
            DPHI(P + 1, 0, 2) = (2.D0*Preal + 1.D0)/(Preal + 1.D0)*(DPHI(1, 0, 2)*PHI(P, 0) &
                                                            + PHI(1, 0)*DPHI(P, 0, 2)) &
                                - Preal/(Preal + 1.D0)*((Y - 1)/2 &
                                                *PHI(P - 1, 0) &
                                                + ((1.D0 - Y)/2.D0)**2.D0*DPHI(P - 1, 0, 2))
         END DO
         DO P = 0, D - 1
            preal = real(P,8)
            PHI(P, 1) = PHI(P, 0)*(1.D0 + 2.D0*Preal + (3.D0 + 2.D0*Preal)*Y)/2.D0
            DPHI(P, 1, 1) = DPHI(P, 0, 1)*(1.D0 + 2.D0*Preal + (3.D0 + 2.D0*Preal)*Y)/2.D0
            DPHI(P, 1, 2) = DPHI(P, 0, 2)*(1.D0 + 2.D0*Preal + (3.D0 + 2.D0*Preal)*Y)/2.D0 &
                            + PHI(P, 0)*(3.D0 + 2.D0*Preal)/2.D0
         END DO

         DO Q = 1, D - 1
            qreal = real(q,8)
            DO P = 0, D - Q - 1
               preal = real(p,8)
               AQ = (2.D0*Qreal + 2.D0 + 2.D0*Preal)*(2.D0*Qreal + 3.D0 + 2.D0*Preal)/(2.D0*Qreal + 2.D0) &
                    /(Qreal + 2.D0 + 2.D0*Preal)
               BQ = (2.D0*Preal + 1.D0)**2.D0*(2.D0*Qreal + 2.D0 + 2.D0*Preal)/ &
                    (2.D0*Qreal + 2.D0)/(2.D0*Qreal + 2.D0*Preal + 1.D0)/(Qreal + 2.D0 + 2.D0*Preal)
               CQ = (Qreal + 2.D0*Preal + 1.D0)*Qreal*(2.D0*Qreal + 3.D0 + 2.D0*Preal)/ &
                    (Qreal + 1.D0)/(Qreal + 2.D0 + 2.D0*Preal)/(2.D0*Qreal + 2.D0*Preal + 1.D0)
               PHI(P, Q + 1) = (AQ*Y + BQ)*PHI(P, Q) - CQ*PHI(P, Q - 1)
               DPHI(P, Q + 1, 1) = (AQ*Y + BQ)*DPHI(P, Q, 1) - CQ*DPHI(P, Q - 1, 1)
               DPHI(P, Q + 1, 2) = (AQ*Y + BQ)*DPHI(P, Q, 2) + AQ*PHI(P, Q) &
                                   - CQ*DPHI(P, Q - 1, 2)
            END DO
         END DO

         !-----------------------------------------------------------------------
         !
         !     2D Quadrilateral basis
         !
         !-----------------------------------------------------------------------

      ELSEIF (ELEM == 2) THEN

         !.......Check to make sure point PT is contained within the reference
         !.......element.  Issue warning if not but continue execution.

         IF ((ABS(X) <= 1.D0) .AND. (ABS(Y) <= 1.D0) .AND. &
             (ABS(Z) <= 1.D0)) THEN
            CONTINUE
         ELSE
            PRINT *, '  ***** WARNING in suroutine orthogonal_basis *****  '
            PRINT *, '       Point is outside of reference element!!       '
            PRINT *, '               Execution will continue               '
         END IF

         !.......Legendre Polynomials in X

         LEGENDRE_X(0) = 1.D0
         DLEGENDRE_X(0) = 0.D0
         LEGENDRE_X(1) = X
         DLEGENDRE_X(1) = 1.D0
         DO P = 1, D - 1
            preal = real(P, 8)
            LEGENDRE_X(P + 1) = ((2.D0*Preal + 1.D0)*X*LEGENDRE_X(P) - Preal* &
                                 LEGENDRE_X(P - 1))/(Preal + 1.D0)
            DLEGENDRE_X(P + 1) = ((2.D0*Preal + 1.D0)*(X*DLEGENDRE_X(P) + &
                                                   LEGENDRE_X(P)) - Preal*DLEGENDRE_X(P - 1))/(Preal + 1.D0)
         END DO

         !.......Legendre Polynomials in Y

         LEGENDRE_Y(0) = 1.D0
         DLEGENDRE_Y(0) = 0.D0
         LEGENDRE_Y(1) = Y
         DLEGENDRE_Y(1) = 1.D0
         DO P = 1, D - 1
            preal = real(p, 8)
            LEGENDRE_Y(P + 1) = ((2.D0*Preal + 1.D0)*Y*LEGENDRE_Y(P) - Preal* &
                                 LEGENDRE_Y(P - 1))/(Preal + 1.D0)
            DLEGENDRE_Y(P + 1) = ((2.D0*Preal + 1.D0)*(Y*DLEGENDRE_Y(P) + &
                                                   LEGENDRE_Y(P)) - Preal*DLEGENDRE_Y(P - 1))/(Preal + 1.D0)
         END DO

         !.......Forming the basis functions
         DO P = 0, D
            DO Q = 0, D - P
               PHI(Q, P) = LEGENDRE_X(Q)*LEGENDRE_Y(P)
               !.......DPHI with respect to X
               DPHI(Q, P, 1) = DLEGENDRE_X(Q)*LEGENDRE_Y(P)
               !.......DPHI with respect to Y
               DPHI(Q, P, 2) = LEGENDRE_X(Q)*DLEGENDRE_Y(P)
            END DO
         END DO
      END IF

      !.....Re-order basis functions hierarchically

      R = 1
      DO Q = 0, D
         DO P = 0, Q
            PHI_2D(R) = PHI(P, Q - P)
            DPHI_2D(R, 1) = DPHI(P, Q - P, 1)
            DPHI_2D(R, 2) = DPHI(P, Q - P, 2)
            IF (ABS(PHI_2D(R)) < 1.0d-15) PHI_2D(R) = 0.D0
            IF (ABS(DPHI_2D(R, 1)) < 1.0d-15) DPHI_2D(R, 1) = 0.D0
            IF (ABS(DPHI_2D(R, 2)) < 1.0d-15) DPHI_2D(R, 2) = 0.D0
            R = R + 1
         END DO
      END DO

      IF (DIM == 2) THEN
         BASIS = PHI_2D
         DBASIS = DPHI_2D
      ELSEIF (DIM == 3) THEN

         !-----------------------------------------------------------------------
         !
         !     3D Bases
         !
         !-----------------------------------------------------------------------

         LEGENDRE(0) = 1.D0
         DLEGENDRE(0) = 0.D0
         LEGENDRE(1) = Z
         DLEGENDRE(1) = 1.D0
         DO P = 1, D - 1
            preal = real(p, 8)
            LEGENDRE(p+ 1) = ((2.D0*preal + 1.D0)*Z*LEGENDRE(p) - preal*LEGENDRE(p-1)) &
                              /(preal + 1.D0)
            DLEGENDRE(p+ 1) = ((2.D0*preal + 1.D0)* &
                                (Z*DLEGENDRE(p) + LEGENDRE(p)) - preal*DLEGENDRE(p-1))/(preal + 1.D0)
         END DO
         R = 1
         DO Q = 0, D
            DO P = 1, (D + 2)*(D + 1)/2
               BASIS(R) = PHI_2D(P)*LEGENDRE(Q)
               DBASIS(R, 1) = DPHI_2D(P, 1)*LEGENDRE(Q)
               DBASIS(R, 2) = DPHI_2D(P, 2)*LEGENDRE(Q)
               DBASIS(R, 3) = PHI_2D(P)*DLEGENDRE(Q)
               IF (ABS(BASIS(R)) < 1.0d-15) BASIS(R) = 0.D0
               IF (ABS(DBASIS(R, 1)) < 1.0d-15) DBASIS(R, 1) = 0.D0
               IF (ABS(DBASIS(R, 2)) < 1.0d-15) DBASIS(R, 2) = 0.D0
               IF (ABS(DBASIS(R, 3)) < 1.0d-15) DBASIS(R, 3) = 0.D0
               R = R + 1
            END DO
         END DO
      END IF

      RETURN
   END SUBROUTINE ORTHOGONAL_BASIS

   !
   !     SUBROUTINE PREP_SLOPELIM()
   !
   !     This subroutine does preparatory stuff for the slope limiter
   !
   !     Written by Ethan Kubatko (07-13-2005)
   !     01-10-2011 - cem - adapted for slew of new limiters
   !     NOTE: This is a very expensive step and can be saved in a file
   !     for multiple runs; needs to be optimized still
   !
   !***********************************************************************

   SUBROUTINE PREP_SLOPELIM()

      !.....Use appropriate modules

      !namo - change neigh_elem to neightabele for now, not sure
      use mesh, only: NM, X, Y

      IMPLICIT NONE

      !.....Declare local variables

      INTEGER :: GED, i, j, n1, n2, n3

      Real(SZ), Allocatable :: tempmat(:, :), tempInv(:, :), tempag(:, :)
      Real(SZ), Allocatable :: AreaV_integral(:, :, :, :), A(:, :)
      Real(SZ), Allocatable :: Full_M_inv(:, :), temp_p(:, :), temp_t(:, :)
      Real(SZ), Allocatable :: Taylor_mass(:, :, :)

      Allocate (tempmat(dofh, dofh), tempInv(dofh, dofh), tempag(dofh, dofh))
      Allocate (AreaV_integral(MNE, 0:ph, 0:ph, 3), A(dofh, dofh))
      Allocate (Full_M_inv(dofh, dofh), temp_p(dofh, dofh))
      Allocate (Taylor_mass(MNE, dofh, dofh), temp_t(dofh, dofh))

      !.....Initialize SL3

      SL3 = 0

      !.....Loop over the elements

      DO J = 1, MNE

         !.....Retrieve the nodal coordinates of the given element

         N1 = NM(J, 1)
         N2 = NM(J, 2)
         N3 = NM(J, 3)

         !.....Compute the barycenter coordinates of the element and store

         XBC(J) = 1.D0/3.D0*(X(N1) + X(N2) + X(N3))
         YBC(J) = 1.D0/3.D0*(Y(N1) + Y(N2) + Y(N3))

         X1 = XBC(J)
         Y1 = YBC(J)

         !.....Loop over the edges to fnd the neighboring elements

         DO I = 1, 3

            !.....Retrieve the global edge number of the element for the given edge

            GED = NELED(I, J)

            !.....Retrieve the neighboring element number and store into an array

            EL_NBORS(I, J) = NEDEL(1, GED)

            IF (EL_NBORS(I, J) == J) EL_NBORS(I, J) = NEDEL(2, GED)

            !.....If the element has an edge that is on boundary go to next element
            !     sb-2007/07/27 commented out

            !     IF (EL_NBORS(I,J).EQ.0) GOTO 111

         END DO

         !.....Set the 4th neighbroing element equal to the first

         EL_NBORS(4, J) = EL_NBORS(1, J)

         !.....Now loop over the three edges of the element again

         DO I = 1, 3

            IF (EL_NBORS(I, J) == 0 .OR. EL_NBORS(I + 1, J) == 0) cycle

            !.....Compute the barycenter coordinates of two neighboring elements

            N1 = NM(EL_NBORS(I, J), 1)
            N2 = NM(EL_NBORS(I, J), 2)
            N3 = NM(EL_NBORS(I, J), 3)

            X2 = 1.D0/3.D0*(X(N1) + X(N2) + X(N3))
            Y2 = 1.D0/3.D0*(Y(N1) + Y(N2) + Y(N3))

            N1 = NM(EL_NBORS(I + 1, J), 1)
            N2 = NM(EL_NBORS(I + 1, J), 2)
            N3 = NM(EL_NBORS(I + 1, J), 3)

            X3 = 1.D0/3.D0*(X(N1) + X(N2) + X(N3))
            Y3 = 1.D0/3.D0*(Y(N1) + Y(N2) + Y(N3))

            !.....Compute the time independent planar constant

            SL3(I, J) = X1*(Y2 - Y3) + X2*(Y3 - Y1) + X3*(Y1 - Y2)

            IF (SL3(I, J) <= 0 .AND. SLOPEFLAG /= 0) then
               WRITE (16, *) 'WARNING. SL3(', I, ',', J, ') =', SL3(I, J), ' <= 0.', &
                  '    ELEMENT ', J, &
                  ' WILL NOT BE CONSIDERED IN SLOPE LIMITING.'
               WRITE (16, *)
               WRITE (*, *) 'WARNING. SL3(', I, ',', J, ') =', SL3(I, J), ' <= 0.', &
                  '    ELEMENT ', J, &
                  ' WILL NOT BE CONSIDERED IN SLOPE LIMITING.'
               WRITE (*, *)
               SL3(I, J) = 0.D0
            END IF

         END DO

      END DO

   END subroutine prep_slopelim

   !
   !     SUBROUTINE STA_BASIS()
   !
   !     This subroutine computes the basis functions at a given station
   !
   !     Written by Ethan Kubatko
   !
   !***********************************************************************

   SUBROUTINE STA_BASIS(ELEM, DIM, XSTA, YSTA, ELSTA, PHI_STA)

      !$$$      USE GLOBAL
      !$$$      USE SIZES
      use mesh, only: X, Y, areas, NM

      IMPLICIT NONE

      integer, intent(in) :: elem, elsta
      real(sz), intent(in) :: xsta, ysta
      INTEGER :: SZ2
      INTEGER, intent(in) :: DIM
      REAL(SZ), intent(out) :: PHI_STA(DOF)
      REAL(SZ) ::  AREA
      REAL(SZ) :: Z1, Z2
      REAL(8) :: PT(DIM)
      REAL(8), Allocatable  :: BASIS(:), DBASIS(:, :)
      INTEGER :: i

      !.....Retrieve element coordinates and area

      X1 = X(NM(ELSTA, 1))
      X2 = X(NM(ELSTA, 2))
      X3 = X(NM(ELSTA, 3))

      Y1 = Y(NM(ELSTA, 1))
      Y2 = Y(NM(ELSTA, 2))
      Y3 = Y(NM(ELSTA, 3))

      AREA = 0.5D0*AREAS(ELSTA)

      !.....Transform to local Z element coordinates

      Z1 = -1.D0/(2.D0*AREA)*(XSTA*(2.D0*Y1 - 2.D0*Y3) + &
                              YSTA*(2.D0*X3 - 2.D0*X1) + &
                              X1*Y2 + X1*Y3 - X2*Y1 + X2*Y3 - X3*Y1 - X3*Y2)

      Z2 = 1.D0/(2.D0*AREA)*(XSTA*(2.D0*Y1 - 2.D0*Y2) + &
                             YSTA*(2.D0*X2 - 2.D0*X1) + &
                             X1*Y2 + X1*Y3 - X2*Y1 - X2*Y3 - X3*Y1 + X3*Y2)

      !.....Compute the basis functions at that point and store
      PT(1) = Z1
      PT(2) = Z2

      SZ2 = (PDG_EL(ELSTA) + 2)*(PDG_EL(ELSTA) + 1)**(DIM - 1)/2

      Allocate (BASIS(SZ2), DBASIS(SZ2, DIM))

      do i = 1, DIM
         if (abs(real(NINT(PT(i)), 8)) - PT(i) <= 1.0d-12) then
            PT(i) = REAL(NINT(PT(i)), 8)
         end if
      end do

      CALL ORTHOGONAL_BASIS(ELEM, PT, PDG_EL(ELSTA), DIM, BASIS, DBASIS)

      PHI_STA(:) = BASIS

      RETURN
   END SUBROUTINE STA_BASIS

   SUBROUTINE ALLOC_RK()

      ALLOCATE (ATVD(NRK, NRK), BTVD(NRK, NRK), CTVD(NRK, NRK))
      ALLOCATE (DTVD(NRK), MAX_BOA_DT(NRK))

      RETURN
   END SUBROUTINE ALLOC_RK

   SUBROUTINE RK_TIME()
      use global, only: dt

      IMPLICIT NONE

      INTEGER :: L, i, j, k, IRK
      REAL(SZ) :: ARK, BRK, CASUM, MAX_BOA

      !print *, 'Running RK_TIME()'
      !.....Allocate the time stepping arrays

      !NRK = RK_STAGE
      CALL ALLOC_RK()

      !.....The forward Euler method

      IF ((RK_STAGE == 1) .AND. (RK_ORDER == 1)) THEN

         ATVD(:, :) = 0.D0
         BTVD(:, :) = 0.D0
         CTVD(:, :) = 0.D0
         DTVD(:) = 0.D0

         ATVD(1, 1) = 1.D0
         BTVD(1, 1) = 1.D0

         !.....SSP(s,2) schemes
         !.....THESE ARE DG-OPTIMIZED METHODS!
         !..... Kubatko, Ethan J., Benjamin A. Yeager, and David I. Ketcheson.
         !....."Optimal strong-stability-preserving Runge–Kutta time
         !..... discretizations for discontinuous Galerkin methods."
         !..... Journal of Scientific Computing 60.2 (2014): 313-344.
      ELSEIF (RK_ORDER == 2) THEN
         ATVD(:, :) = 0.D0
         BTVD(:, :) = 0.D0
         CTVD(:, :) = 0.D0
         DTVD(:) = 0.D0
         IF (RK_STAGE == 3) THEN

            ATVD(1, 1) = 1.00000000000000D0
            ATVD(2, 1) = 0.087353119859156D0
            ATVD(3, 1) = 0.344956917166841D0
            ATVD(2, 2) = 0.912646880140844D0
            ATVD(3, 3) = 0.655043082833159D0

            BTVD(1, 1) = 0.528005024856522D0
            BTVD(3, 1) = 0.022826837460491D0
            BTVD(2, 2) = 0.481882138633993D0
            BTVD(3, 3) = 0.345866039233415D0

         ELSEIF (RK_STAGE == 4) THEN

            ATVD(1, 1) = 1.00000000000000D0
            ATVD(2, 1) = 0.394806441339829D0
            ATVD(3, 1) = 0.002797307087390D0
            ATVD(4, 1) = 0.252860909354373D0
            ATVD(2, 2) = 0.605193558660171D0
            ATVD(3, 3) = 0.997202692912610D0
            ATVD(4, 4) = 0.747139090645627D0

            BTVD(1, 1) = 0.406584463657504D0
            BTVD(3, 1) = 0.013637216641451D0
            BTVD(4, 1) = 0.016453567333598D0
            BTVD(2, 2) = 0.246062298456822D0
            BTVD(3, 3) = 0.405447122055692D0
            BTVD(4, 4) = 0.303775146447707D0

         ELSEIF (RK_STAGE == 5) THEN

            ATVD(1, 1) = 1.00000000000000D0
            ATVD(2, 1) = 0.235593265061659D0
            ATVD(3, 1) = 0.174017972351526D0
            ATVD(4, 1) = 0.235264368870758D0
            ATVD(5, 1) = 0.141720372339803D0
            ATVD(2, 2) = 0.764406734938341D0
            ATVD(4, 2) = 0.000058643383967D0
            ATVD(5, 2) = 0.095374613155521D0
            ATVD(3, 3) = 0.825982027648475D0
            ATVD(5, 3) = 0.000311763705780D0
            ATVD(4, 4) = 0.764676987745275D0
            ATVD(5, 5) = 0.762593250798895D0

            BTVD(1, 1) = 0.324840618151514D0
            BTVD(3, 1) = 0.108822380501601D0
            BTVD(4, 1) = 0.054392262422093D0
            BTVD(5, 1) = 0.000000180291569D0
            BTVD(2, 2) = 0.248310356296551D0
            BTVD(4, 2) = 0.000019049753098D0
            BTVD(5, 2) = 0.030981548293401D0
            BTVD(3, 3) = 0.268312512443371D0
            BTVD(5, 3) = 0.000101273514903D0
            BTVD(4, 4) = 0.248398145385413D0
            BTVD(5, 5) = 0.247721262987686D0
         ELSE

            ATVD(:, :) = 0.D0
            BTVD(:, :) = 0.D0
            CTVD(:, :) = 0.D0
            DTVD(:) = 0.D0

            DO I = 1, NRK
               DO J = 0, NRK - 1

                  IF ((J == (I - 1)) .AND. (I < NRK)) THEN
                     ATVD(I, J + 1) = 1.D0
                     BTVD(I, J + 1) = 1.D0/(real(NRK,8) - 1)
                  ELSEIF ((J == 0) .AND. (I == NRK)) THEN
                     ATVD(I, J + 1) = 1.D0/real(NRK,8)
                  ELSEIF ((J == (NRK - 1)) .AND. (I == NRK)) THEN
                     ATVD(I, J + 1) = (real(NRK,8) - 1.D0)/real(NRK,8)
                     BTVD(I, J + 1) = 1.D0/real(NRK, 8)
                  END IF

               END DO
            END DO
         END IF
         !.....SSP(3,3) scheme

      ELSEIF ((RK_STAGE == 3) .AND. (RK_ORDER == 3)) THEN

         ATVD(:, :) = 0.D0
         BTVD(:, :) = 0.D0
         CTVD(:, :) = 0.D0
         DTVD(:) = 0.D0

         ATVD(1, 1) = 1.D0
         ATVD(2, 1) = 3.D0/4.D0
         ATVD(2, 2) = 1.D0/4.D0
         ATVD(3, 1) = 1.D0/3.D0
         ATVD(3, 3) = 2.D0/3.D0

         BTVD(1, 1) = 1.D0
         BTVD(2, 2) = 1.D0/4.D0
         BTVD(3, 3) = 2.D0/3.D0

         !.....SSP(4,3) scheme

      ELSEIF ((RK_STAGE == 4) .AND. (RK_ORDER == 3)) THEN

         ATVD(:, :) = 0.D0
         BTVD(:, :) = 0.D0
         CTVD(:, :) = 0.D0
         DTVD(:) = 0.D0

         ATVD(1, 1) = 1.D0
         ATVD(2, 2) = 1.D0
         ATVD(3, 1) = 2.D0/3.D0
         ATVD(3, 3) = 1.D0/3.D0
         ATVD(4, 4) = 1.D0

         BTVD(1, 1) = 1.D0/2.D0
         BTVD(2, 2) = 1.D0/2.D0
         BTVD(3, 3) = 1.D0/6.D0
         BTVD(4, 4) = 1.D0/2.D0

         !.....SSP(5,3) scheme

      ELSEIF ((RK_STAGE == 5) .AND. (RK_ORDER == 3)) THEN

         ATVD(:, :) = 0.D0
         BTVD(:, :) = 0.D0
         CTVD(:, :) = 0.D0
         DTVD(:) = 0.D0

         ATVD(1, 1) = 1.D0
         ATVD(2, 2) = 1.D0
         ATVD(3, 1) = 0.355909775063327D0
         ATVD(3, 3) = 0.644090224936674D0
         ATVD(4, 1) = 0.367933791638137D0
         ATVD(4, 4) = 0.632066208361863D0
         ATVD(5, 3) = 0.237593836598569D0
         ATVD(5, 5) = 0.762406163401431D0

         BTVD(1, 1) = 0.377268915331368D0
         BTVD(2, 2) = 0.377268915331368D0
         BTVD(3, 3) = 0.242995220537396D0
         BTVD(4, 4) = 0.238458932846290D0
         BTVD(5, 5) = 0.287632146308408D0

         !.....SSP(6,3) scheme

      ELSEIF ((RK_STAGE == 6) .AND. (RK_ORDER == 3)) THEN

         ATVD(:, :) = 0.D0
         BTVD(:, :) = 0.D0
         CTVD(:, :) = 0.D0
         DTVD(:) = 0.D0

         ATVD(1, 1) = 1.D0
         ATVD(2, 2) = 1.D0
         ATVD(3, 3) = 1.D0
         ATVD(4, 1) = 0.476769811285196D0
         ATVD(4, 2) = 0.098511733286064D0
         ATVD(4, 4) = 0.424718455428740D0
         ATVD(5, 5) = 1.D0
         ATVD(6, 3) = 0.155221702560091D0
         ATVD(6, 6) = 0.844778297439909D0

         BTVD(1, 1) = 0.284220721334261D0
         BTVD(2, 2) = 0.284220721334261D0
         BTVD(3, 3) = 0.284220721334261D0
         BTVD(4, 4) = 0.120713785765930D0
         BTVD(5, 5) = 0.284220721334261D0
         BTVD(6, 6) = 0.240103497065900D0

         !.....SSP(7,3) scheme

      ELSEIF ((RK_STAGE == 7) .AND. (RK_ORDER == 3)) THEN

         ATVD(:, :) = 0.D0
         BTVD(:, :) = 0.D0
         CTVD(:, :) = 0.D0
         DTVD(:) = 0.D0

         ATVD(1, 1) = 1.D0
         ATVD(2, 2) = 1.D0
         ATVD(3, 3) = 1.D0
         ATVD(4, 1) = 0.184962588071072D0
         ATVD(4, 4) = 0.815037411928928D0
         ATVD(5, 1) = 0.180718656570380D0
         ATVD(5, 2) = 0.314831034403793D0
         ATVD(5, 5) = 0.504450309025826D0
         ATVD(6, 6) = 1.D0
         ATVD(7, 4) = 0.120199000000000D0
         ATVD(7, 7) = 0.879801000000000D0

         BTVD(1, 1) = 0.233213863663009D0
         BTVD(2, 2) = 0.233213863663009D0
         BTVD(3, 3) = 0.233213863663009D0
         BTVD(4, 4) = 0.190078023865845D0
         BTVD(5, 5) = 0.117644805593912D0
         BTVD(6, 6) = 0.233213863663009D0
         BTVD(7, 7) = 0.205181790464579D0

         !.....SSP(8,3) scheme

      ELSEIF ((RK_STAGE == 8) .AND. (RK_ORDER == 3)) THEN

         ATVD(:, :) = 0.D0
         BTVD(:, :) = 0.D0
         CTVD(:, :) = 0.D0
         DTVD(:) = 0.D0

         ATVD(1, 1) = 1.D0
         ATVD(2, 2) = 1.D0
         ATVD(3, 3) = 1.D0
         ATVD(4, 4) = 1.D0
         ATVD(5, 1) = 0.421366967085359D0
         ATVD(5, 2) = 0.005949401107575D0
         ATVD(5, 5) = 0.572683631807067D0
         ATVD(6, 2) = 0.004254010666365D0
         ATVD(6, 6) = 0.995745989333635D0
         ATVD(7, 3) = 0.104380143093325D0
         ATVD(7, 4) = 0.243265240906726D0
         ATVD(7, 7) = 0.652354615999950D0
         ATVD(8, 8) = 1.D0

         BTVD(1, 1) = 0.195804015330143D0
         BTVD(2, 2) = 0.195804015330143D0
         BTVD(3, 3) = 0.195804015330143D0
         BTVD(4, 4) = 0.195804015330143D0
         BTVD(5, 5) = 0.112133754621673D0
         BTVD(6, 6) = 0.194971062960412D0
         BTVD(7, 7) = 0.127733653231944D0
         BTVD(8, 8) = 0.195804015330143D0

         !.....SSP(5,4) scheme

      ELSEIF ((RK_STAGE == 5) .AND. (RK_ORDER == 4)) THEN

         ATVD(:, :) = 0.D0
         BTVD(:, :) = 0.D0
         CTVD(:, :) = 0.D0
         DTVD(:) = 0.D0

         ATVD(1, 1) = 1.D0
         ATVD(2, 1) = 0.44437049406734D0
         ATVD(2, 2) = 0.55562950593266D0
         ATVD(3, 1) = 0.62010185138540D0
         ATVD(3, 3) = 0.37989814861460D0
         ATVD(4, 1) = 0.17807995410773D0
         ATVD(4, 4) = 0.82192004589227D0
         ATVD(5, 1) = 0.00683325884039D0
         ATVD(5, 3) = 0.51723167208978D0
         ATVD(5, 4) = 0.12759831133288D0
         ATVD(5, 5) = 0.34833675773694D0

         BTVD(1, 1) = 0.39175222700392D0
         BTVD(2, 2) = 0.36841059262959D0
         BTVD(3, 3) = 0.25189177424738D0
         BTVD(4, 4) = 0.54497475021237D0
         BTVD(5, 4) = 0.08460416338212D0
         BTVD(5, 5) = 0.22600748319395D0

         !.....SSP(6,4) scheme

      ELSEIF ((RK_STAGE == 6) .AND. (RK_ORDER == 4)) THEN

         ATVD(:, :) = 0.D0
         BTVD(:, :) = 0.D0
         CTVD(:, :) = 0.D0
         DTVD(:) = 0.D0

         ATVD(1, 1) = 1.00000000000000D0
         ATVD(2, 1) = 0.30948026455053D0
         ATVD(2, 2) = 0.69051973544947D0
         ATVD(3, 1) = 0.54205244285557D0
         ATVD(3, 3) = 0.45794755714443D0
         ATVD(4, 1) = 0.35984960863377D0
         ATVD(4, 4) = 0.64015039136623D0
         ATVD(5, 5) = 1.00000000000000D0
         ATVD(6, 1) = 0.05776282890116D0
         ATVD(6, 3) = 0.44216432622405D0
         ATVD(6, 5) = 0.10115567086469D0
         ATVD(6, 6) = 0.39891717401009D0

         BTVD(1, 1) = 0.39270746575722D0
         BTVD(2, 2) = 0.30154043149172D0
         BTVD(3, 3) = 0.19997937335132D0
         BTVD(4, 4) = 0.27954483459696D0
         BTVD(5, 5) = 0.43668618869443D0
         BTVD(6, 3) = 0.09150931531680D0
         BTVD(6, 5) = 0.04417328437472D0
         BTVD(6, 6) = 0.14911300530736D0

         !.....SSP(7,4) scheme

      ELSEIF ((RK_STAGE == 7) .AND. (RK_ORDER == 4)) THEN

         ATVD(:, :) = 0.D0
         BTVD(:, :) = 0.D0
         CTVD(:, :) = 0.D0
         DTVD(:) = 0.D0

         ATVD(1, 1) = 1.D0
         ATVD(2, 1) = 0.20161507213829D0
         ATVD(2, 2) = 0.79838492786171D0
         ATVD(3, 1) = 0.19469598207921D0
         ATVD(3, 3) = 0.80530401792079D0
         ATVD(4, 1) = 0.58143386885601D0
         ATVD(4, 4) = 0.41856613114399D0
         ATVD(5, 1) = 0.01934367892154D0
         ATVD(5, 5) = 0.98065632107846D0
         ATVD(6, 6) = 1.D0
         ATVD(7, 1) = 0.06006304558847D0
         ATVD(7, 3) = 0.30152730794242D0
         ATVD(7, 4) = 0.10518998496676D0
         ATVD(7, 5) = 0.01483791154585D0
         ATVD(7, 7) = 0.51838174995650D0

         BTVD(1, 1) = 0.30111872706068D0
         BTVD(2, 2) = 0.24040865318216D0
         BTVD(3, 3) = 0.24249212077315D0
         BTVD(4, 4) = 0.12603810060080D0
         BTVD(5, 5) = 0.29529398308716D0
         BTVD(6, 6) = 0.30111872706068D0
         BTVD(7, 3) = 0.09079551914158D0
         BTVD(7, 4) = 0.02888359354880D0
         BTVD(7, 7) = 0.15609445267839D0

         !.....SSP(8,4) scheme

      ELSEIF ((RK_STAGE == 8) .AND. (RK_ORDER == 4)) THEN

         ATVD(:, :) = 0.D0
         BTVD(:, :) = 0.D0
         CTVD(:, :) = 0.D0
         DTVD(:) = 0.D0

         ATVD(1, 1) = 1.D0
         ATVD(2, 1) = 0.10645325745007D0
         ATVD(2, 2) = 0.89354674254993D0
         ATVD(3, 3) = 1.D0
         ATVD(4, 1) = 0.57175518477257D0
         ATVD(4, 4) = 0.42824481522743D0
         ATVD(5, 1) = 0.19161667219044D0
         ATVD(5, 5) = 0.80838332780956D0
         ATVD(6, 6) = 1.D0
         ATVD(7, 7) = 1.D0
         ATVD(8, 1) = 0.02580435327923D0
         ATVD(8, 3) = 0.03629901341774D0
         ATVD(8, 4) = 0.31859181340256D0
         ATVD(8, 5) = 0.05186768980103D0
         ATVD(8, 6) = 0.03944076217320D0
         ATVD(8, 7) = 0.00511633747411D0
         ATVD(8, 8) = 0.52288003045213D0

         BTVD(1, 1) = 0.24120020561311D0
         BTVD(2, 2) = 0.21552365802797D0
         BTVD(3, 3) = 0.24120020561311D0
         BTVD(4, 4) = 0.10329273748560D0
         BTVD(5, 5) = 0.19498222488188D0
         BTVD(6, 6) = 0.24120020561311D0
         BTVD(7, 7) = 0.24120020561311D0
         BTVD(8, 3) = 0.00875532949991D0
         BTVD(8, 4) = 0.06195575835101D0
         BTVD(8, 6) = 0.00951311994571D0
         BTVD(8, 8) = 0.12611877085604D0

      END IF

      !.....Compute the time dependent parameters

      DO K = 0, RK_STAGE - 1
         DO I = 1, RK_STAGE
            CASUM = 0.D0
            DO L = K + 1, I - 1
               CASUM = CASUM + CTVD(L, K + 1)*ATVD(I, L + 1)
            END DO
            CTVD(I, K + 1) = BTVD(I, K + 1) + CASUM
         END DO
      END DO

      DO K = 1, RK_STAGE - 1
         DTVD(K + 1) = 0.D0
         DO L = 0, K - 1
            DTVD(K + 1) = DTVD(K + 1) + CTVD(K, L + 1)
         END DO
      END DO

      !.....Compute the maximum beta over alpha ratio at each stage

      DO IRK = 1, RK_STAGE
         MAX_BOA = 0.D0
         DO I = 1, IRK
            ARK = ATVD(IRK, I)
            BRK = BTVD(IRK, I)
            IF (abs(ARK) >= 1d-15) THEN
               IF (MAX_BOA < BRK/ARK) MAX_BOA = BRK/ARK
            END IF
         END DO
         MAX_BOA_DT(IRK) = MAX_BOA*DT
      END DO

      !-----------------------------------------------------------
      !.... Compute the Runge-Kutta Chebyshev (RKC) version

      RETURN
   END SUBROUTINE RK_TIME

   subroutine nodal_to_modal(A, B)
    !! Compute DG modal representation of nodal A and store result in B

      use mesh, only: NM

      implicit none
      real(sz), intent(in) :: A(:)
      real(sz), intent(out) :: B(:, :)

      integer :: j, n1, n2, n3

      DO J = 1, MNE
         N1 = NM(J, 1)
         N2 = NM(J, 2)
         N3 = NM(J, 3)

         B(1, J) = 1.D0/3.D0*(A(N1) + A(N2) + A(N3))
         B(2, J) = -1.D0/6.D0*(A(N1) + A(N2)) + 1.D0/3.D0*A(N3)
         B(3, J) = -0.5D0*A(N1) + 0.5D0*A(N2)
      END DO
   end subroutine nodal_to_modal

   subroutine modal_to_area_quad(A, B)
    !! Compute the area quadrature values of a given modal representation

      implicit none
      real(sz), intent(in) :: A(:, :)
      real(sz), intent(out) :: B(:, :, :)

      integer :: i, j, k, chi

      DO J = 1, MNE
         do chi = 1, ph
            do I = 1, NAGP(chi) ! Area quadrature points

               B(I, J, chi) = 0.D0

               DO K = 1, (chi + 1)*(chi + 2)/2
                  B(I, J, chi) = B(I, J, chi) + A(K, J)*PHI_AREA(K, I, chi)
               END DO
            END DO
         end do
      END DO
   end subroutine modal_to_area_quad

   subroutine modal_to_edge_quad(A, B)
    !! Compute the edge quadrature values of a given modal representation

      implicit none
      real(sz), intent(in) :: A(:, :)
      real(sz), intent(out) :: B(:, :, :, :)

      integer :: i, j, k, L, chi

      do J = 1, MNE
         do chi = 1, ph
            DO L = 1, 3
               do I = 1, NEGP(chi) ! Edge quadrature points

                  B(I, L, J, chi) = 0.D0

                  DO K = 1, (chi + 1)*(chi + 2)/2
                     B(I, L, J, chi) = B(I, L, J, chi) + A(K, J)*PHI_EDGE(K, I, L, chi)
                  END DO
               end do
            END DO
         END DO
      end do
   end subroutine modal_to_edge_quad

   subroutine nodal_to_quad_points(A, B, C, D)
    !! Given a nodal array A, compute its modal representation into B,
    !! its area quad points into C, and its edge quadrature values into D
      implicit none
      real(sz), intent(in) :: A(:)
      real(sz), intent(out) :: B(:, :)
      real(sz), intent(out) :: C(:, :, :)
      real(sz), intent(out) :: D(:, :, :, :)

      call nodal_to_modal(A, B)
      call modal_to_area_quad(B, C)
      call modal_to_edge_quad(B, D)
   end subroutine nodal_to_quad_points

END MODULE DG
