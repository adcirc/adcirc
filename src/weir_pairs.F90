! ADCIRC: internal weir pairs at the new time level
!-----------------------------------------------------------------------
!  M O D U L E   W E I R _ P A I R S
!-----------------------------------------------------------------------
!  Internal weir pairs (IBTYPE 4 and 24) at the new time level.
!
!  An internal weir is two strings of boundary nodes facing each other
!  across a crest; each facing couple is a pair, side a on the first
!  string and side b on the second. A pair moves water from a to b at a
!  rate q per unit width, and each side's normal flux is QN = -q on a
!  and +q on b (positive into the domain).
!
!  The pair's flow term in the GWCE is f = c+ q(n+1) - c- q(n-1), with
!  c+- = 1/(2 dt) +- tau0/2 (tau0 the mean of the two nodes), through
!  one length per pair, l = (l_a + l_b)/2, where each side's l_j is half
!  of each boundary edge next to it on its own string (the chord across
!  the gap is left out). What a loses, b gains. With the flux taken at
!  the new time level, the flux and the levels it moves are at the same
!  time level, so the weir cannot overcorrect from one step to the next.
!
!  Water crosses only between wet nodes. A side may spill only when its
!  node and a neighbour along the weir are wet. Wet first, flow next: a
!  pair with a dry side passes nothing; a wet source above the crest and
!  above the dry node's level, at least H0 above the lowest level that
!  covers every corner of one of the dry node's elements (its reach),
!  marks the dry node for wetting (NIBNODECODE), and the pair flows from
!  the next step. A receiving node of a flowing pair is marked as well,
!  as before, which keeps its elements wet.
!
!  Each side's own tau0: the pair is solved with the mean of the two
!  nodes' tau0, which keeps the consistent matrix symmetric, and each
!  side adds (tau0_j - mean) q(n) explicitly, as the rest of its GWCE row
!  carries its own tau0; pairs with unequal tau0 then conserve water.
!
!  Lumped GWCE (WEIR_PAIRS_SOLVE_LUMPED): after the diagonal division,
!  each pair solves its flow term f against the levels without the weir,
!  zeta_a* and zeta_b*:
!     G(f) = f + c- q(n-1) - c+ q(zeta_a* - r_a f, zeta_b* + r_b f) = 0,
!  r_j = k_j l / D_j, with k_j the number of pairs at node j and D_j its
!  diagonal. G increases with f, so the root is unique; it lies between
!  the f that stops the flow and the f that brings the two sides level.
!  It is found by false position with the Illinois weighting in
!  x = sqrt(distance from the levelling end), where the law is smooth.
!  Each node then takes every pair's term through its whole diagonal.
!
!  Consistent GWCE (WEIR_PAIRS_LINEARIZE / WEIR_PAIRS_FINISH): the solve
!  is for the increments, so each pair's flow at n+1 is written as a line
!  in the two increments, q = q0 + s (dz_a - dz_b), the law's chord
!  across the span the last step moved the level difference (at least
!  1 micrometre), about the levels that last increment extrapolates to.
!  The part c+ l s (dz_a - dz_b) goes into the matrix as the symmetric
!  positive semidefinite block c+ l s [1 -1; -1 1] in the rows and
!  columns of a and b, the rest, l (c+ q0 - c- q(n-1)), into the load.
!  One solve per step; the flux is read back at the solved increments.
!  The two nodes share no element, so the matrix gets a column for the
!  partner (WEIR_PAIRS_GRAPH); the block is added to the assembled matrix
!  before the solve and removed after it.
!-----------------------------------------------------------------------
module WEIR_PAIRS
   use ADC_CONSTANTS, only: G
   implicit none
   private

   !...Solve the IBTYPE 4/24 pairs at the new time level (default) or, with
   !   &boundaryControl IMPLICIT_WEIR = .false., take their flux from the
   !   previous step's levels (WEIR_FLUX)
   logical, public :: IMPLICIT_WEIR = .true.
   !...The spilling side reads its level over the crest (default) or, with
   !   &boundaryControl WEIR_TOTAL_HEAD = .true., its total head: the level
   !   plus the velocity head of the node's velocity toward the weir, from the
   !   momentum equations at time level n (here and in WEIR_FLUX)
   logical, public :: WEIR_TOTAL_HEAD = .false.
   !...consistent GWCE: ITPACK's JCG with the blocks of the connected pairs (WEIR_BLOCKS)
   !   in place of the diagonal alone (default) or, with &boundaryControl
   !   WEIR_PAIR_BLOCK_JACOBI = .false., the diagonal alone
   logical, public :: WEIR_PAIR_BLOCK_JACOBI = .true.
   !...set by the GWCE for the solve that uses them
   logical, public :: WEIR_BLOCKS_ON = .false.

   integer :: NPAIR = -1 !...-1: not set up yet
   logical, public :: ANY_PAIRS = .false. !...some subdomain has a pair
   integer, allocatable :: PA(:), PB(:) !...boundary entries of sides a, b
   integer, allocatable :: NA_NEI(:, :), NB_NEI(:, :) !...nodes next along the weir (0: none)
   real(8), allocatable :: PLEN(:) !...the pair's length l [m]
   integer, allocatable :: KAT(:) !...number of pairs at each node
   real(8), allocatable :: REACH(:) !...lowest level covering one of a node's elements
   real(8), allocatable :: VHA(:), VHB(:) !...velocity head toward the weir on sides a, b [m]
   !...boundary entries whose forcing the pairs supply (skip in QFORCE)
   logical, allocatable, public :: WEIR_PAIR_ENTRY(:)

   !...each boundary entry's crest, kept current by the weir code (time-varying weirs)
   real(8), allocatable, public :: WEIR_PAIR_CREST(:)

   !...consistent GWCE: the matrix's neighbour table with the weir partners
   integer, public :: MNEIW = 0
   integer, allocatable, public :: NEITABW(:, :), NNEIGHW(:)
   integer, allocatable :: SLOTA(:), SLOTB(:) !...column of b in a's row, of a in b's row
   real(8), allocatable :: LSLOPE(:), LBASE(:) !...the pair's line
   real(8), allocatable :: LBLOCK(:) !...its coefficient in the matrix
   integer, allocatable :: LRA(:), LRB(:) !...the rows it is in
   real(8), allocatable :: SAVED(:, :) !...matrix entries before the blocks
   logical, allocatable :: LACTIVE(:)

   !...the pairs' blocks (WEIR_BLOCKS): components of the flowing pairs with a resident
   !   node; component c has the nodes BNODE(BPTR(c)+1:BPTR(c+1)) and the matrices C and
   !   C^-1/2 by rows in BC, BMAT(BMPTR(c)+1:BMPTR(c+1))
   integer, parameter :: MAXC = 4 !...largest component whose rows go to ghost nodes
   integer :: NCOMP = 0, BMODE = 0
   integer, allocatable :: BNODE(:), BPTR(:), BMPTR(:)
   real(8), allocatable :: BC(:), BMAT(:)
   integer, allocatable :: BLIST(:), BROOT(:), BPOS(:) !...work: nodes of the flowing pairs, their roots
   real(8), allocatable :: BWT(:), BWS(:), BWY(:) !...work of WEIR_BLOCKS_BT, _DT (at the blocks' nodes only)
#ifdef CMPI
   !...ghost nodes with a row: columns BGIDX(BGPTR(k)+1:BGPTR(k+1)), entries BGVAL
   integer :: NGROW = 0
   integer, allocatable :: BGNODE(:), BGPTR(:), BGIDX(:)
   real(8), allocatable :: BGVAL(:)
   integer, allocatable :: BGHOST(:) !...the subdomain's ghost nodes
   integer, allocatable :: BSFULL(:), BSLOC(:) !...full-domain node numbers, sorted, and their nodes
   integer :: NBSENT = 0
   integer, allocatable :: BSENT(:) !...nodes whose rows are in BRVAL, BRCOL
   real(8), allocatable :: BRVAL(:, :), BRCOL(:, :) !...rows sent: entries, full-domain column numbers
#endif

   public :: WEIR_PAIRS_SETUP, WEIR_PAIRS_SYNC_MARKS, WEIR_PAIRS_SOLVE_LUMPED, WEIR_OVERFLOW, &
             WEIR_PAIRS_GRAPH, WEIR_PAIRS_LINEARIZE, WEIR_PAIRS_FINISH, &
             WEIR_BLOCKS_SETUP, WEIR_BLOCKS_HALF, WEIR_BLOCKS_FULL, WEIR_BLOCKS_LOCAL, &
             WEIR_BLOCKS_BT, WEIR_BLOCKS_DT
#ifdef CMPI
   public :: WEIR_BLOCKS_EXCHANGE
#endif

contains

   !-----------------------------------------------------------------------
   !  WEIR_PAIRS_SETUP: build the pair table from the IBTYPE 4/24 boundaries
   !-----------------------------------------------------------------------
   subroutine WEIR_PAIRS_SETUP()
      use BOUNDARIES, only: NBOU, NVELL, NVEL, NBV, LBCODEI, BndLen2O3, BARINHT
      use MESH, only: NP
#ifdef CMPI
      use MESSENGER, only: msg_imax
#endif
      implicit none
      integer :: I0, K, J, NV, P, IA, IB
      real(8), allocatable :: elen(:)

      if (NPAIR >= 0) return
      allocate (WEIR_PAIR_ENTRY(NVEL), elen(NVEL), KAT(NP))
      if (.not. allocated(WEIR_PAIR_CREST)) then
         allocate (WEIR_PAIR_CREST(NVEL))
         if (NVEL > 0) WEIR_PAIR_CREST = BARINHT(1:NVEL)
      end if
      WEIR_PAIR_ENTRY = .false.
      elen = 0.d0
      KAT = 0
      NPAIR = 0
      if (.not. IMPLICIT_WEIR) return
      !...count the pairs and set each entry's length
      I0 = 0
      do K = 1, NBOU
         NV = NVELL(K)
         select case (LBCODEI(I0 + 1))
         case (4, 24)
            do J = 2, NV
               elen(I0 + J - 1) = elen(I0 + J - 1) + 0.75d0*BndLen2O3(I0 + J - 1)
               elen(I0 + J) = elen(I0 + J) + 0.75d0*BndLen2O3(I0 + J - 1)
               elen(I0 + NV + J - 1) = elen(I0 + NV + J - 1) + 0.75d0*BndLen2O3(I0 + NV + J - 1)
               elen(I0 + NV + J) = elen(I0 + NV + J) + 0.75d0*BndLen2O3(I0 + NV + J - 1)
            end do
            do J = 1, NV
               if (NBV(I0 + J) > 0 .and. NBV(I0 + 2*NV + 1 - J) > 0) NPAIR = NPAIR + 1
            end do
            I0 = I0 + 2*NV
         case (5, 25)
            I0 = I0 + 2*NV
         case default
            I0 = I0 + NV
         end select
      end do
      allocate (PA(NPAIR), PB(NPAIR), PLEN(NPAIR), NA_NEI(2, NPAIR), NB_NEI(2, NPAIR))
      NA_NEI = 0
      NB_NEI = 0
      P = 0
      I0 = 0
      do K = 1, NBOU
         NV = NVELL(K)
         select case (LBCODEI(I0 + 1))
         case (4, 24)
            do J = 1, NV
               IA = I0 + J
               !...the back face is stored in reverse order (mesh.F)
               IB = I0 + 2*NV + 1 - J
               if (NBV(IA) <= 0 .or. NBV(IB) <= 0) cycle
               P = P + 1
               PA(P) = IA
               PB(P) = IB
               PLEN(P) = 0.5d0*(elen(IA) + elen(IB))
               if (J > 1) then
                  NA_NEI(1, P) = NBV(IA - 1)
                  NB_NEI(1, P) = NBV(IB + 1)
               end if
               if (J < NV) then
                  NA_NEI(2, P) = NBV(IA + 1)
                  NB_NEI(2, P) = NBV(IB - 1)
               end if
               KAT(NBV(IA)) = KAT(NBV(IA)) + 1
               KAT(NBV(IB)) = KAT(NBV(IB)) + 1
               WEIR_PAIR_ENTRY(IA) = .true.
               WEIR_PAIR_ENTRY(IB) = .true.
            end do
            I0 = I0 + 2*NV
         case (5, 25)
            I0 = I0 + 2*NV
         case default
            I0 = I0 + NV
         end select
      end do
      deallocate (elen)
      call SET_REACH()
      ANY_PAIRS = NPAIR > 0
#ifdef CMPI
      ANY_PAIRS = msg_imax(merge(1, 0, ANY_PAIRS)) > 0
#endif
   end subroutine WEIR_PAIRS_SETUP

   !-----------------------------------------------------------------------
   !  WEIR_PAIRS_SYNC_MARKS: a ghost node takes its owner's wetting mark.
   !  The owner holds all of the node's elements, its reach and the pairs
   !  along its weir; another subdomain can see only some of them.
   !  Called on every subdomain after the GWCE solve.
   !-----------------------------------------------------------------------
   subroutine WEIR_PAIRS_SYNC_MARKS()
#ifdef CMPI
      use GLOBAL, only: NIBNODECODE
      use MESSENGER, only: UPDATEI
      implicit none
      integer :: idum(1)
      if (ANY_PAIRS) call UPDATEI(NIBNODECODE, idum, 1)
#endif
   end subroutine WEIR_PAIRS_SYNC_MARKS

   !-----------------------------------------------------------------------
   !  SET_REACH: for each node, the lowest water level that covers every
   !  corner of one of its elements (a dry node has nowhere to put water
   !  below it)
   !-----------------------------------------------------------------------
   subroutine SET_REACH()
      use MESH, only: NP, NE, NM, DP
      implicit none
      integer :: IE, K
      real(8) :: top
      allocate (REACH(NP))
      REACH = huge(1.d0)
      do IE = 1, NE
         top = max(-DP(NM(IE, 1)), -DP(NM(IE, 2)), -DP(NM(IE, 3)))
         do K = 1, 3
            REACH(NM(IE, K)) = min(REACH(NM(IE, K)), top)
         end do
      end do
   end subroutine SET_REACH

   !-----------------------------------------------------------------------
   !  MARK_PAIR: mark the receiving node of a flowing pair, and the dry node
   !  of a pair whose wet source can wet one of its elements, for wetting
   !  (NIBNODECODE); la, lb are the two sides' levels, q the flow a -> b
   !-----------------------------------------------------------------------
   subroutine MARK_PAIR(P, la, lb, q, NODECODE)
      use BOUNDARIES, only: NBV
      use GLOBAL, only: NIBNODECODE, H0
      implicit none
      integer, intent(IN) :: P
      real(8), intent(IN) :: la, lb, q
      integer, intent(IN) :: NODECODE(:)
      integer :: NA, NB
      real(8) :: crest
      NA = NBV(PA(P))
      NB = NBV(PB(P))
      if (q > 0.d0) NIBNODECODE(NB) = 1
      if (q < 0.d0) NIBNODECODE(NA) = 1
      crest = WEIR_PAIR_CREST(PA(P))
      if (NODECODE(NB) == 0 .and. la > max(crest, lb) .and. la >= REACH(NB) + H0) then
         if (CAN_SOURCE(NA, NA_NEI(:, P), NODECODE)) NIBNODECODE(NB) = 1
      end if
      if (NODECODE(NA) == 0 .and. lb > max(crest, la) .and. lb >= REACH(NA) + H0) then
         if (CAN_SOURCE(NB, NB_NEI(:, P), NODECODE)) NIBNODECODE(NA) = 1
      end if
   end subroutine MARK_PAIR

   !-----------------------------------------------------------------------
   !  WEIR_OVERFLOW: overflow per unit width over a crest [m^2/s] for the
   !  heads on the higher (upper) and lower side, measured from the crest.
   !  Free flow while the lower head is at most 2/3 of the upper; above
   !  that the lesser of the free and submerged laws, which keeps the flow
   !  rising with the upper head and never rising with the lower one.
   !-----------------------------------------------------------------------
   real(8) pure function WEIR_OVERFLOW(upper, lower, csub, cfree) result(q)
      implicit none
      real(8), intent(IN) :: upper, lower, csub, cfree
      real(8) :: hm
      q = 0.d0
      if (upper <= 0.d0) return
      hm = 2.d0*upper/3.d0
      q = cfree*hm*sqrt(G*hm)
      if (lower > hm) q = min(q, csub*lower*sqrt(2.d0*G*(upper - lower)))
   end function WEIR_OVERFLOW

   !-----------------------------------------------------------------------
   !  PAIR_FLOW: flow per unit width from side a to side b at levels ea, eb;
   !  va, vb: each side's velocity head toward the weir. The spilling side
   !  reads its total head, its level over the crest plus its velocity head;
   !  the lower side keeps its level, since the velocity that leaves over the
   !  crest is lost behind it. The velocity heads are fixed during a step, so
   !  the flow still never falls as the spilling side's level rises.
   !-----------------------------------------------------------------------
   real(8) pure function PAIR_FLOW(ea, eb, crest, csub, cfree, ramp, asrc, bsrc, va, vb) result(q)
      implicit none
      real(8), intent(IN) :: ea, eb, crest, csub, cfree, ramp, va, vb
      logical, intent(IN) :: asrc, bsrc
      q = 0.d0
      if (ea - crest >= eb - crest) then
         if (asrc .and. ea > crest) q = ramp*WEIR_OVERFLOW(ea - crest + va, eb - crest, csub, cfree)
      else
         if (bsrc .and. eb > crest) q = -ramp*WEIR_OVERFLOW(eb - crest + vb, ea - crest, csub, cfree)
      end if
   end function PAIR_FLOW

   !-----------------------------------------------------------------------
   !  SET_VELOCITY_HEADS: each side's velocity head toward the weir, from
   !  the node's velocity at time level n along the side's outward normal
   !  (zero for a velocity away from the weir, or without WEIR_TOTAL_HEAD)
   !-----------------------------------------------------------------------
   subroutine SET_VELOCITY_HEADS()
      use GLOBAL, only: UU1, VV1
      use BOUNDARIES, only: NBV, CSII, SIII
      implicit none
      integer :: P
      real(8) :: vn
      if (.not. allocated(VHA)) allocate (VHA(NPAIR), VHB(NPAIR))
      VHA = 0.d0
      VHB = 0.d0
      if (.not. WEIR_TOTAL_HEAD) return
      do P = 1, NPAIR
         vn = max(UU1(NBV(PA(P)))*CSII(PA(P)) + VV1(NBV(PA(P)))*SIII(PA(P)), 0.d0)
         VHA(P) = vn*vn/(2.d0*G)
         vn = max(UU1(NBV(PB(P)))*CSII(PB(P)) + VV1(NBV(PB(P)))*SIII(PB(P)), 0.d0)
         VHB(P) = vn*vn/(2.d0*G)
      end do
   end subroutine SET_VELOCITY_HEADS

   !-----------------------------------------------------------------------
   !  CAN_SOURCE: a side may spill when its node and a neighbour along the
   !  weir are wet
   !-----------------------------------------------------------------------
   logical function CAN_SOURCE(node, nei, nodecode)
      implicit none
      integer, intent(IN) :: node, nei(2)
      integer, intent(IN) :: nodecode(:)
      CAN_SOURCE = .false.
      if (nodecode(node) == 0) return
      if (nei(1) > 0) then
         if (nodecode(nei(1)) /= 0) CAN_SOURCE = .true.
      end if
      if (nei(2) > 0) then
         if (nodecode(nei(2)) /= 0) CAN_SOURCE = .true.
      end if
   end function CAN_SOURCE

   !-----------------------------------------------------------------------
   !  WEIR_PAIRS_SOLVE_LUMPED: solve every pair at the new time level after
   !  the lumped GWCE division and correct the increments ETAS of its two
   !  nodes; the solved flux goes to QN2 (positive into the domain).
   !  D is the lumped diagonal, ETA1 the level at n, QN0 the flux at n-1.
   !-----------------------------------------------------------------------
   subroutine WEIR_PAIRS_SOLVE_LUMPED(ETAS, D, ETA1, TAU0V, DT, QN0, QN1, QN2, NODECODE, RAMP)
      use BOUNDARIES, only: NBV, BARINCFSB, BARINCFSP
      use MESH, only: NP
      implicit none
      real(8), intent(INOUT) :: ETAS(:)
      real(8), intent(IN) :: D(:), ETA1(:), TAU0V(:), QN0(:), QN1(:)
      real(8), intent(IN) :: DT, RAMP
      real(8), intent(INOUT) :: QN2(:)
      integer, intent(IN) :: NODECODE(:)
      real(8), allocatable :: change(:)
      integer :: P, NA, NB
      real(8) :: la, lb, ra, rb, cp, cm, qold, carried, crest, csub, cfree, f, q, tau, qn, xa, xb
      logical :: asrc, bsrc

      if (NPAIR <= 0) return
      call SET_VELOCITY_HEADS()
      allocate (change(NP))
      change = 0.d0
      do P = 1, NPAIR
         NA = NBV(PA(P))
         NB = NBV(PB(P))
         QN2(PA(P)) = 0.d0
         QN2(PB(P)) = 0.d0
         if (NODECODE(NA) == 0 .or. NODECODE(NB) == 0) then
            !...wet first, flow next
            la = ETA1(NA)
            lb = ETA1(NB)
            if (NODECODE(NA) /= 0) la = la + ETAS(NA)
            if (NODECODE(NB) /= 0) lb = lb + ETAS(NB)
            call MARK_PAIR(P, la, lb, 0.d0, NODECODE)
            cycle
         end if
         if (D(NA) <= 0.d0 .or. D(NB) <= 0.d0) cycle
         !...each side's own tau0: (tau0_j - mean) q(n), explicit
         tau = 0.5d0*(TAU0V(NA) + TAU0V(NB))
         qn = 0.5d0*(QN1(PB(P)) - QN1(PA(P)))
         xa = -PLEN(P)*(TAU0V(NA) - tau)*qn/D(NA)
         xb = PLEN(P)*(TAU0V(NB) - tau)*qn/D(NB)
         change(NA) = change(NA) + xa
         change(NB) = change(NB) + xb
         la = ETA1(NA) + ETAS(NA) + xa
         lb = ETA1(NB) + ETAS(NB) + xb
         !...a node on k pairs answers each pair k times as strongly
         ra = dble(KAT(NA))*PLEN(P)/D(NA)
         rb = dble(KAT(NB))*PLEN(P)/D(NB)
         cp = 0.5d0/DT + 0.5d0*tau
         cm = 0.5d0/DT - 0.5d0*tau
         qold = 0.5d0*(QN0(PB(P)) - QN0(PA(P))) ! a -> b at n-1
         carried = cm*qold
         crest = WEIR_PAIR_CREST(PA(P))
         csub = BARINCFSB(PA(P))
         cfree = BARINCFSP(PA(P))
         asrc = CAN_SOURCE(NA, NA_NEI(:, P), NODECODE)
         bsrc = CAN_SOURCE(NB, NB_NEI(:, P), NODECODE)
         f = solve_pair()
         q = (f + carried)/cp
         QN2(PA(P)) = -q
         QN2(PB(P)) = q
         change(NA) = change(NA) - PLEN(P)*f/D(NA)
         change(NB) = change(NB) + PLEN(P)*f/D(NB)
         call MARK_PAIR(P, la - ra*f, lb + rb*f, q, NODECODE)
      end do
      do NA = 1, NP
         if (abs(change(NA)) > 0.d0) ETAS(NA) = ETAS(NA) + change(NA)
      end do
      deallocate (change)

   contains

      real(8) function resid(term)
         real(8), intent(IN) :: term
         resid = term + carried - cp*PAIR_FLOW(la - ra*term, lb + rb*term, crest, &
                                               csub, cfree, RAMP, asrc, bsrc, VHA(P), VHB(P))
      end function resid

      !...the root of resid, as in Cocoa's solve_weir_pair
      real(8) function solve_pair() result(term)
         real(8) :: rn, rl, levelling, span, kept, rkept, latest, rlatest, nxt, rnxt
         integer :: it
         term = -carried
         rn = resid(term)
         if (abs(rn) <= 0.d0) return
         levelling = (la - lb)/(ra + rb)
         rl = resid(levelling)
         if ((rn > 0.d0) .eqv. (rl > 0.d0)) then
            if (abs(rl) < abs(rn)) term = levelling
            return
         end if
         !...x = 0 is the levelling end, x = sqrt|span| the no-flow end
         span = levelling + carried
         kept = sqrt(abs(span))
         rkept = rn
         latest = 0.d0
         rlatest = rl
         do it = 1, 64
            if (abs(rlatest) <= 0.d0 .or. abs(latest - kept) <= 1.d-14*(abs(latest) + abs(kept))) exit
            nxt = (kept*rlatest - latest*rkept)/(rlatest - rkept)
            rnxt = resid(levelling - sign(nxt*nxt, span))
            if ((rnxt > 0.d0) .eqv. (rlatest > 0.d0)) then
               rkept = 0.5d0*rkept
            else
               kept = latest
               rkept = rlatest
            end if
            latest = nxt
            rlatest = rnxt
         end do
         term = levelling - sign(latest*latest, span)
      end function solve_pair

   end subroutine WEIR_PAIRS_SOLVE_LUMPED

   !-----------------------------------------------------------------------
   !  WEIR_PAIRS_GRAPH: the consistent GWCE matrix's neighbour table, the
   !  mesh table plus each weir node's partner(s) (columns appended after
   !  the mesh's, so the stored element positions stay valid). Unused
   !  slots point to the node itself and carry a zero coefficient.
   !  Call after WEIR_PAIRS_SETUP.
   !-----------------------------------------------------------------------
   subroutine WEIR_PAIRS_GRAPH(NEITAB, NNEIGH, MNEI)
      use BOUNDARIES, only: NBV
      use MESH, only: NP
      implicit none
      integer, intent(IN) :: NEITAB(:, :), NNEIGH(:), MNEI
      integer, allocatable :: extra(:, :), nextra(:)
      integer :: P, I, J, NA, NB, mx
      logical :: found

      allocate (nextra(NP), extra(NP, max(1, maxval(KAT))))
      nextra = 0
      do P = 1, NPAIR
         NA = NBV(PA(P))
         NB = NBV(PB(P))
         call add(NA, NB)
         call add(NB, NA)
      end do
      mx = 0
      do I = 1, NP
         mx = max(mx, NNEIGH(I) + nextra(I))
      end do
      MNEIW = max(MNEI, mx)
      allocate (NEITABW(size(NEITAB, 1), MNEIW), NNEIGHW(size(NEITAB, 1)))
      do I = 1, size(NEITAB, 1)
         NEITABW(I, :) = I
      end do
      NNEIGHW = 0
      do I = 1, NP
         NEITABW(I, 1:NNEIGH(I)) = NEITAB(I, 1:NNEIGH(I))
         do J = 1, nextra(I)
            NEITABW(I, NNEIGH(I) + J) = extra(I, J)
         end do
         NNEIGHW(I) = NNEIGH(I) + nextra(I)
      end do
      allocate (SLOTA(NPAIR), SLOTB(NPAIR), LSLOPE(NPAIR), LBASE(NPAIR), SAVED(4, NPAIR), &
                LACTIVE(NPAIR), LBLOCK(NPAIR), LRA(NPAIR), LRB(NPAIR))
      do P = 1, NPAIR
         NA = NBV(PA(P))
         NB = NBV(PB(P))
         SLOTA(P) = slot(NA, NB)
         SLOTB(P) = slot(NB, NA)
      end do
      LACTIVE = .false.
      deallocate (nextra, extra)

   contains

      subroutine add(i, j)
         integer, intent(IN) :: i, j
         integer :: k
         do k = 1, NNEIGH(i)
            if (NEITAB(i, k) == j) return
         end do
         found = .false.
         do k = 1, nextra(i)
            if (extra(i, k) == j) found = .true.
         end do
         if (found) return
         nextra(i) = nextra(i) + 1
         extra(i, nextra(i)) = j
      end subroutine add

      integer function slot(i, j)
         integer, intent(IN) :: i, j
         integer :: k
         slot = 0
         do k = 2, NNEIGHW(i)
            if (NEITABW(i, k) == j) then
               slot = k
               return
            end if
         end do
      end function slot

   end subroutine WEIR_PAIRS_GRAPH

   !-----------------------------------------------------------------------
   !  WEIR_PAIRS_LINEARIZE: before the consistent GWCE solve, add each wet
   !  pair's block to the matrix COEF (saving the entries it changes) and
   !  its constant part to the load GWCE_LV. ETA1 is the level at n and
   !  ETAS the last step's increment. Pairs at an elevation-specified
   !  boundary node are left out (their rows are fixed).
   !-----------------------------------------------------------------------
   subroutine WEIR_PAIRS_LINEARIZE(COEF, GWCE_LV, ETA1, ETAS, TAU0V, DT, QN0, QN1, NODECODE, RAMP)
      use BOUNDARIES, only: NBV, BARINCFSB, BARINCFSP, NETA, NBD
      use MESH, only: NP
      implicit none
      real(8), intent(INOUT) :: COEF(:, :), GWCE_LV(:)
      real(8), intent(IN) :: ETA1(:), ETAS(:), TAU0V(:), QN0(:), QN1(:), DT, RAMP
      integer, intent(IN) :: NODECODE(:)
      logical, allocatable, save :: ISOBC(:)
      integer :: P, NA, NB
      real(8) :: ea, eb, span, up, dn, q0, s, cp, cm, tau, qold, m, rhs, crest, csub, cfree, qn
      logical :: asrc, bsrc

      if (NPAIR <= 0) return
      call SET_VELOCITY_HEADS()
      if (.not. allocated(ISOBC)) then
         allocate (ISOBC(NP))
         ISOBC = .false.
         do P = 1, NETA
            if (NBD(P) > 0) ISOBC(NBD(P)) = .true.
         end do
      end if
      do P = 1, NPAIR
         LACTIVE(P) = .false.
         LSLOPE(P) = 0.d0
         LBASE(P) = 0.d0
         NA = NBV(PA(P))
         NB = NBV(PB(P))
         if (NODECODE(NA) == 0 .or. NODECODE(NB) == 0) cycle
         if (ISOBC(NA) .or. ISOBC(NB)) cycle
         if (SLOTA(P) == 0 .or. SLOTB(P) == 0) cycle
         crest = WEIR_PAIR_CREST(PA(P))
         csub = BARINCFSB(PA(P))
         cfree = BARINCFSP(PA(P))
         asrc = CAN_SOURCE(NA, NA_NEI(:, P), NODECODE)
         bsrc = CAN_SOURCE(NB, NB_NEI(:, P), NODECODE)
         !...the levels the last increment extrapolates to, and the span the
         !   last step moved their difference
         ea = ETA1(NA) + ETAS(NA)
         eb = ETA1(NB) + ETAS(NB)
         span = max(abs(ETAS(NA) - ETAS(NB)), 1.d-6)
         tau = 0.5d0*(TAU0V(NA) + TAU0V(NB))
         cp = 0.5d0/DT + 0.5d0*tau
         cm = 0.5d0/DT - 0.5d0*tau
         if (ea > crest .or. eb > crest) then
            up = PAIR_FLOW(ea + 0.5d0*span, eb - 0.5d0*span, crest, csub, cfree, RAMP, asrc, bsrc, &
                           VHA(P), VHB(P))
            dn = PAIR_FLOW(ea - 0.5d0*span, eb + 0.5d0*span, crest, csub, cfree, RAMP, asrc, bsrc, &
                           VHA(P), VHB(P))
            s = (up - dn)/(2.d0*span)
            !...the line's flow at the difference at n
            q0 = 0.5d0*(up + dn) + s*((ETA1(NA) - ETA1(NB)) - (ea - eb))
         else
            s = 0.d0
            q0 = 0.d0
         end if
         qold = 0.5d0*(QN0(PB(P)) - QN0(PA(P)))
         m = cp*PLEN(P)*s
         rhs = PLEN(P)*(cp*q0 - cm*qold)

         SAVED(1, P) = COEF(NA, 1)
         SAVED(2, P) = COEF(NA, SLOTA(P))
         SAVED(3, P) = COEF(NB, 1)
         SAVED(4, P) = COEF(NB, SLOTB(P))
         COEF(NA, 1) = COEF(NA, 1) + m
         COEF(NA, SLOTA(P)) = COEF(NA, SLOTA(P)) - m
         COEF(NB, 1) = COEF(NB, 1) + m
         COEF(NB, SLOTB(P)) = COEF(NB, SLOTB(P)) - m
         GWCE_LV(NA) = GWCE_LV(NA) - rhs
         GWCE_LV(NB) = GWCE_LV(NB) + rhs
         !...each side's own tau0: (tau0_j - mean) q(n), explicit
         qn = 0.5d0*(QN1(PB(P)) - QN1(PA(P)))
         GWCE_LV(NA) = GWCE_LV(NA) - PLEN(P)*(TAU0V(NA) - tau)*qn
         GWCE_LV(NB) = GWCE_LV(NB) + PLEN(P)*(TAU0V(NB) - tau)*qn
         LSLOPE(P) = s
         LBASE(P) = q0
         LBLOCK(P) = m
         LRA(P) = NA
         LRB(P) = NB
         LACTIVE(P) = .true.
      end do
   end subroutine WEIR_PAIRS_LINEARIZE

   !-----------------------------------------------------------------------
   !  WEIR_PAIRS_FINISH: after the consistent GWCE solve, remove the blocks
   !  from COEF (in reverse order, restoring the saved entries) and write
   !  each pair's flux at the solved increments ETAS to QN2.
   !-----------------------------------------------------------------------
   subroutine WEIR_PAIRS_FINISH(COEF, ETAS, ETA1, QN2, NODECODE)
      use BOUNDARIES, only: NBV
      implicit none
      real(8), intent(INOUT) :: COEF(:, :), QN2(:)
      real(8), intent(IN) :: ETAS(:), ETA1(:)
      integer, intent(IN) :: NODECODE(:)
      real(8) :: la, lb
      integer :: P, NA, NB
      real(8) :: q

      if (NPAIR <= 0) return
      do P = NPAIR, 1, -1
         NA = NBV(PA(P))
         NB = NBV(PB(P))
         QN2(PA(P)) = 0.d0
         QN2(PB(P)) = 0.d0
         la = ETA1(NA)
         lb = ETA1(NB)
         if (NODECODE(NA) /= 0) la = la + ETAS(NA)
         if (NODECODE(NB) /= 0) lb = lb + ETAS(NB)
         if (.not. LACTIVE(P)) then
            call MARK_PAIR(P, la, lb, 0.d0, NODECODE)
            cycle
         end if
         COEF(NB, SLOTB(P)) = SAVED(4, P)
         COEF(NB, 1) = SAVED(3, P)
         COEF(NA, SLOTA(P)) = SAVED(2, P)
         COEF(NA, 1) = SAVED(1, P)
         q = LBASE(P) + LSLOPE(P)*(ETAS(NA) - ETAS(NB))
         QN2(PA(P)) = -q
         QN2(PB(P)) = q
         call MARK_PAIR(P, la, lb, q, NODECODE)
      end do
   end subroutine WEIR_PAIRS_FINISH

   !-----------------------------------------------------------------------
   !  The pairs' blocks for ITPACK's JCG (&boundaryControl
   !  WEIR_PAIR_BLOCK_JACOBI). JCG scales the matrix to a unit diagonal. In
   !  the scaled matrix a flowing pair ties its two nodes by one negative
   !  term; the nodes joined by flowing pairs are a component (2 nodes, or
   !  more where pairs share a node), and C is the identity with each
   !  component's unit diagonal and these terms. C is symmetric and positive
   !  definite. JCG then solves C^-1/2 A C^-1/2 instead of A, which takes
   !  the pairs' stiffness out of the iterations:
   !     WEIR_BLOCKS_SETUP  the components of this solve, and their C^-1/2
   !     WEIR_BLOCKS_HALF   C^-1/2 x
   !     WEIR_BLOCKS_FULL   C x
   !  WEIR_BLOCKS_ON is set by the GWCE for its solve, on every subdomain or
   !  none. Components with a resident node are kept; a ghost node takes its
   !  owner's row of C^-1/2 (BMODE 0), so that HALF is right at ghost nodes
   !  without an exchange of its result. If a row cannot be used (a component
   !  of more than MAXC nodes, a node not held) the result is exchanged
   !  (BMODE 1); if the subdomains see a component differently, the components
   !  across subdomains are left out as well (BMODE 2).
   !-----------------------------------------------------------------------
   subroutine WEIR_BLOCKS_SETUP(N, WK)
#ifdef CMPI
      use MESSENGER, only: RESNODE, msg_imax, UPDATER_SUBSET
      use GLOBAL, only: NODES_LG
#endif
      use GLOBAL, only: NPERSEG, PERBC_IDN_MAP
      implicit none
      integer, intent(IN) :: N
      real(8), intent(IN) :: WK(:) !...1/sqrt(diagonal) at every node (ghost: its owner's)
      integer :: P, K, c, ib, jb, kb, n1, nlist, nmat
      real(8) :: t

      if (.not. allocated(BROOT)) call INIT()
      !...the nodes of the flowing pairs, joined into components
      nlist = 0
      do P = 1, NPAIR
         if (.not. PAIR_USED(P)) cycle
         call ADD_NODE(LRA(P))
         call ADD_NODE(LRB(P))
         ib = BFIND(LRA(P))
         jb = BFIND(LRB(P))
         if (ib /= jb) BROOT(max(ib, jb)) = min(ib, jb)
      end do
      !...number the components that have a resident node (BPOS of the root)
      NCOMP = 0
      do K = 1, nlist
         if (.not. OWNED(BLIST(K))) cycle
         c = BFIND(BLIST(K))
         if (BPOS(c) == 0) then
            NCOMP = NCOMP + 1
            BPOS(c) = NCOMP
         end if
      end do
      !...their nodes: BPTR(c+1) counts, then is the next free place, then the end
      BPTR(1:NCOMP + 1) = 0
      do K = 1, nlist
         c = BPOS(BFIND(BLIST(K)))
         if (c > 0) BPTR(c + 1) = BPTR(c + 1) + 1
      end do
      BMPTR(1) = 0
      do c = 1, NCOMP
         BMPTR(c + 1) = BMPTR(c) + BPTR(c + 1)**2
         BPTR(c + 1) = BPTR(c + 1) + BPTR(c)
      end do
      do c = NCOMP, 1, -1
         BPTR(c + 1) = BPTR(c)
      end do
      do K = 1, nlist
         c = BPOS(BFIND(BLIST(K)))
         if (c == 0) cycle
         BPTR(c + 1) = BPTR(c + 1) + 1
         BNODE(BPTR(c + 1)) = BLIST(K)
      end do
      !...C of each component, by rows
      nmat = BMPTR(NCOMP + 1)
      if (allocated(BC)) then
         if (size(BC) < nmat) deallocate (BC, BMAT)
      end if
      if (.not. allocated(BC)) allocate (BC(nmat), BMAT(nmat))
      BC(1:nmat) = 0.d0
      do c = 1, NCOMP
         n1 = BPTR(c + 1) - BPTR(c)
         do ib = 1, n1
            BC(BMPTR(c) + (ib - 1)*n1 + ib) = 1.d0
         end do
      end do
      do P = 1, NPAIR
         if (.not. PAIR_USED(P)) cycle
         c = BPOS(BFIND(LRA(P)))
         if (c == 0) cycle
         n1 = BPTR(c + 1) - BPTR(c)
         ib = 0
         jb = 0
         do kb = 1, n1
            if (BNODE(BPTR(c) + kb) == LRA(P)) ib = kb
            if (BNODE(BPTR(c) + kb) == LRB(P)) jb = kb
         end do
         t = LBLOCK(P)*WK(LRA(P))*WK(LRB(P))
         kb = BMPTR(c) + (ib - 1)*n1 + jb
         BC(kb) = BC(kb) - t
         kb = BMPTR(c) + (jb - 1)*n1 + ib
         BC(kb) = BC(kb) - t
      end do
      !...C^-1/2 of each component
      BMAT(1:nmat) = BC(1:nmat)
      do c = 1, NCOMP
         call INVERSE_ROOT(BMAT(BMPTR(c) + 1:BMPTR(c + 1)), BPTR(c + 1) - BPTR(c))
      end do
      !...clear the work arrays for the next solve
      do K = 1, nlist
         BROOT(BLIST(K)) = 0
         BPOS(BLIST(K)) = 0
      end do
      BMODE = 0
#ifdef CMPI
      call GHOST_ROWS()
#endif

   contains

      !...arrays of the run
      subroutine INIT()
         integer :: np2
#ifdef CMPI
         integer :: i1, ng
#endif
         np2 = max(NPAIR, 0)
         allocate (BROOT(N), BPOS(N))
         BROOT = 0
         BPOS = 0
         allocate (BLIST(2*np2), BNODE(2*np2), BPTR(2*np2 + 1), BMPTR(2*np2 + 1))
         allocate (BWT(N), BWS(N), BWY(N))
#ifdef CMPI
         ng = 0
         do i1 = 1, N
            if (.not. RESNODE(i1)) ng = ng + 1
         end do
         allocate (BGHOST(ng), BGNODE(ng), BGPTR(ng + 1), BGIDX(MAXC*ng), BGVAL(MAXC*ng))
         ng = 0
         do i1 = 1, N
            if (RESNODE(i1)) cycle
            ng = ng + 1
            BGHOST(ng) = i1
         end do
         allocate (BSENT(2*np2), BRVAL(N, MAXC), BRCOL(N, MAXC), BSFULL(N), BSLOC(N))
         BRVAL = 0.d0
         BRCOL = 0.d0
         do i1 = 1, N
            BSLOC(i1) = i1
            BSFULL(i1) = abs(NODES_LG(i1))
         end do
         call SORT_FULL(N)
         call PAIR_LISTS()
#endif
      end subroutine INIT

#ifdef CMPI
      !...the nodes of the pairs, the only ones where the values the blocks
      !   exchange (UPDATER_SUBSET) differ from the ones the ghost nodes have;
      !   the ghost nodes whose rows GHOST_ROWS takes are then these nodes only
      subroutine PAIR_LISTS()
         use MESSENGER, only: MSG_SUBSET_SETUP
         use BOUNDARIES, only: NBV
         real(8), allocatable :: mk(:)
         integer :: p1, i1, ng
         allocate (mk(N))
         mk = 0.d0
         do p1 = 1, NPAIR
            mk(NBV(PA(p1))) = 1.d0
            mk(NBV(PB(p1))) = 1.d0
         end do
         call MSG_SUBSET_SETUP(mk, 2*MAXC)
         if (NPERSEG > 0) return
         ng = 0
         do i1 = 1, size(BGHOST)
            if (mk(BGHOST(i1)) < 0.5d0) cycle
            ng = ng + 1
            BGHOST(ng) = BGHOST(i1)
         end do
         BGHOST = BGHOST(1:ng)
      end subroutine PAIR_LISTS
#endif

      logical function PAIR_USED(P1)
         integer, intent(IN) :: P1
         PAIR_USED = .false.
         if (.not. LACTIVE(P1)) return
         PAIR_USED = LBLOCK(P1) > 0.d0
         !...a periodic slave node's row (IBTYPE 94) is cut off from the rest, so
         !   its pair needs no block
         if (NPERSEG > 0 .and. PAIR_USED) then
            if (PERBC_IDN_MAP(LRA(P1)) /= LRA(P1) .or. PERBC_IDN_MAP(LRB(P1)) /= LRB(P1)) PAIR_USED = .false.
         end if
      end function PAIR_USED

      subroutine ADD_NODE(J1)
         integer, intent(IN) :: J1
         if (BROOT(J1) > 0) return
         BROOT(J1) = J1
         nlist = nlist + 1
         BLIST(nlist) = J1
      end subroutine ADD_NODE

      recursive integer function BFIND(J1) result(R1)
         integer, intent(IN) :: J1
         R1 = J1
         if (BROOT(J1) /= J1) then
            R1 = BFIND(BROOT(J1))
            BROOT(J1) = R1
         end if
      end function BFIND

      logical function OWNED(I1)
         integer, intent(IN) :: I1
#ifdef CMPI
         OWNED = RESNODE(I1)
#else
         OWNED = I1 > 0
#endif
      end function OWNED

      !...A^-1/2 of a symmetric positive definite matrix, by Jacobi rotations
      subroutine INVERSE_ROOT(A1, n2)
         integer, intent(IN) :: n2
         real(8), intent(INOUT) :: A1(n2, n2)
         real(8) :: V(n2, n2), W1(n2), t, cs, sn, th, x1, x2, off
         integer :: i1, j1, k1, it
         V = 0.d0
         do i1 = 1, n2
            V(i1, i1) = 1.d0
         end do
         do it = 1, 60
            off = 0.d0
            do i1 = 1, n2 - 1
               do j1 = i1 + 1, n2
                  off = max(off, abs(A1(i1, j1)))
               end do
            end do
            if (off <= 1.d-15) exit
            do i1 = 1, n2 - 1
               do j1 = i1 + 1, n2
                  if (abs(A1(i1, j1)) <= 1.d-16) cycle
                  th = 0.5d0*(A1(j1, j1) - A1(i1, i1))/A1(i1, j1)
                  t = sign(1.d0, th)/(abs(th) + sqrt(1.d0 + th*th))
                  cs = 1.d0/sqrt(1.d0 + t*t)
                  sn = t*cs
                  do k1 = 1, n2
                     x1 = A1(k1, i1)
                     x2 = A1(k1, j1)
                     A1(k1, i1) = cs*x1 - sn*x2
                     A1(k1, j1) = sn*x1 + cs*x2
                  end do
                  do k1 = 1, n2
                     x1 = A1(i1, k1)
                     x2 = A1(j1, k1)
                     A1(i1, k1) = cs*x1 - sn*x2
                     A1(j1, k1) = sn*x1 + cs*x2
                  end do
                  do k1 = 1, n2
                     x1 = V(k1, i1)
                     x2 = V(k1, j1)
                     V(k1, i1) = cs*x1 - sn*x2
                     V(k1, j1) = sn*x1 + cs*x2
                  end do
               end do
            end do
         end do
         do i1 = 1, n2
            W1(i1) = 1.d0/sqrt(A1(i1, i1))
         end do
         do i1 = 1, n2
            do j1 = 1, n2
               A1(i1, j1) = 0.d0
               do k1 = 1, n2
                  A1(i1, j1) = A1(i1, j1) + V(i1, k1)*W1(k1)*V(j1, k1)
               end do
            end do
         end do
      end subroutine INVERSE_ROOT

#ifdef CMPI
      !...send the resident members' rows, take the ghost nodes' rows, and check
      !   that the subdomains agree on the components (BMODE)
      subroutine GHOST_ROWS()
         integer :: i1, c1, j1, k1, m1, n2, l, nrow
         logical :: held
         do k1 = 1, NBSENT
            BRVAL(BSENT(k1), :) = 0.d0
            BRCOL(BSENT(k1), :) = 0.d0
         end do
         NBSENT = 0
         do c1 = 1, NCOMP
            n2 = BPTR(c1 + 1) - BPTR(c1)
            if (n2 > MAXC) BMODE = 1
            if (n2 > MAXC) cycle
            do k1 = 1, n2
               i1 = BNODE(BPTR(c1) + k1)
               if (.not. OWNED(i1)) cycle
               NBSENT = NBSENT + 1
               BSENT(NBSENT) = i1
               do j1 = 1, n2
                  BRVAL(i1, j1) = BMAT(BMPTR(c1) + (k1 - 1)*n2 + j1)
                  BRCOL(i1, j1) = dble(abs(NODES_LG(BNODE(BPTR(c1) + j1))))
               end do
            end do
         end do
         call UPDATER_SUBSET(N, BRVAL, MAXC, BRCOL, MAXC)
         NGROW = 0
         BGPTR(1) = 0
         do k1 = 1, size(BGHOST)
            i1 = BGHOST(k1)
            nrow = BGPTR(NGROW + 1)
            held = .true.
            do j1 = 1, MAXC
               if (BRCOL(i1, j1) < 0.5d0) exit
               l = LOCAL_NODE(nint(BRCOL(i1, j1)))
               if (l == 0) then
                  held = .false.
                  exit
               end if
               nrow = nrow + 1
               BGIDX(nrow) = l
               BGVAL(nrow) = BRVAL(i1, j1)
            end do
            if (.not. held) BMODE = max(BMODE, 1)
            if (.not. held .or. nrow == BGPTR(NGROW + 1)) cycle
            NGROW = NGROW + 1
            BGNODE(NGROW) = i1
            BGPTR(NGROW + 1) = nrow
         end do
         !...the ghost members of the local components must have the owner's rows
         do c1 = 1, NCOMP
            n2 = BPTR(c1 + 1) - BPTR(c1)
            do k1 = 1, n2
               i1 = BNODE(BPTR(c1) + k1)
               if (OWNED(i1)) cycle
               if (n2 > MAXC) then
                  BMODE = 2
                  cycle
               end if
               do j1 = 1, n2
                  l = 0
                  do m1 = 1, MAXC
                     if (nint(BRCOL(i1, m1)) == abs(NODES_LG(BNODE(BPTR(c1) + j1)))) l = m1
                  end do
                  if (l == 0) then
                     BMODE = 2
                  else if (abs(BRVAL(i1, l) - BMAT(BMPTR(c1) + (k1 - 1)*n2 + j1)) > 1.d-8*abs(BRVAL(i1, l))) then
                     BMODE = 2
                  end if
               end do
               if (n2 < MAXC) then
                  if (BRCOL(i1, n2 + 1) > 0.5d0) BMODE = 2
               end if
            end do
         end do
         BMODE = msg_imax(BMODE)
         if (BMODE < 2) return
         !...leave the components across subdomains out
         do c1 = 1, NCOMP
            n2 = BPTR(c1 + 1) - BPTR(c1)
            held = .true.
            do k1 = 1, n2
               if (.not. OWNED(BNODE(BPTR(c1) + k1))) held = .false.
            end do
            if (held) cycle
            BMAT(BMPTR(c1) + 1:BMPTR(c1 + 1)) = 0.d0
            BC(BMPTR(c1) + 1:BMPTR(c1 + 1)) = 0.d0
            do k1 = 1, n2
               BMAT(BMPTR(c1) + (k1 - 1)*n2 + k1) = 1.d0
               BC(BMPTR(c1) + (k1 - 1)*n2 + k1) = 1.d0
            end do
         end do
      end subroutine GHOST_ROWS

      !...subdomain node of a full-domain node, 0 if not held
      integer function LOCAL_NODE(IFULL)
         integer, intent(IN) :: IFULL
         integer :: lo, hi, mid
         LOCAL_NODE = 0
         lo = 1
         hi = N
         do while (lo <= hi)
            mid = (lo + hi)/2
            if (BSFULL(mid) == IFULL) then
               LOCAL_NODE = BSLOC(mid)
               return
            else if (BSFULL(mid) < IFULL) then
               lo = mid + 1
            else
               hi = mid - 1
            end if
         end do
      end function LOCAL_NODE

      subroutine SORT_FULL(N2)
         integer, intent(IN) :: N2
         integer :: k1, m1, t1, t2
         do k1 = N2/2, 1, -1
            call SIFT(k1, N2)
         end do
         do m1 = N2, 2, -1
            t1 = BSFULL(1)
            t2 = BSLOC(1)
            BSFULL(1) = BSFULL(m1)
            BSLOC(1) = BSLOC(m1)
            BSFULL(m1) = t1
            BSLOC(m1) = t2
            call SIFT(1, m1 - 1)
         end do
      end subroutine SORT_FULL

      subroutine SIFT(ilo, ihi)
         integer, intent(IN) :: ilo, ihi
         integer :: a1, b1, t1, t2
         t1 = BSFULL(ilo)
         t2 = BSLOC(ilo)
         a1 = ilo
         b1 = 2*a1
         do while (b1 <= ihi)
            if (b1 < ihi) then
               if (BSFULL(b1 + 1) > BSFULL(b1)) b1 = b1 + 1
            end if
            if (BSFULL(b1) <= t1) exit
            BSFULL(a1) = BSFULL(b1)
            BSLOC(a1) = BSLOC(b1)
            a1 = b1
            b1 = 2*a1
         end do
         BSFULL(a1) = t1
         BSLOC(a1) = t2
      end subroutine SIFT
#endif

   end subroutine WEIR_BLOCKS_SETUP

#ifdef CMPI
   !-----------------------------------------------------------------------
   !  WEIR_BLOCKS_EXCHANGE: ghost nodes take their owners' values, as JCG's
   !  exchange (with periodic boundaries, a ghost slave node its primary's)
   !-----------------------------------------------------------------------
   subroutine WEIR_BLOCKS_EXCHANGE(V1, V2, V3, NMSG)
      use MESSENGER, only: UPDATER, UPDATER_W_PERBC
      use GLOBAL, only: NPERSEG
      implicit none
      integer, intent(IN) :: NMSG
      real(8), intent(INOUT) :: V1(:), V2(:), V3(:)
      if (NPERSEG > 0) then
         call UPDATER_W_PERBC(V1, V2, V3, NMSG)
      else
         call UPDATER(V1, V2, V3, NMSG)
      end if
   end subroutine WEIR_BLOCKS_EXCHANGE
#endif

   !-----------------------------------------------------------------------
   !  WEIR_BLOCKS_HALF: Y = C^-1/2 X; X valid at ghost nodes too
   !-----------------------------------------------------------------------
   subroutine WEIR_BLOCKS_HALF(N, X, Y)
#ifdef CMPI
      use MESSENGER, only: RESNODE, UPDATER_SUBSET
#endif
      implicit none
      integer, intent(IN) :: N
      real(8), intent(IN) :: X(:)
      real(8), intent(INOUT) :: Y(:)
      integer :: c, ii, jj, n1
      real(8) :: t
      real(8) :: dum(1)
      dum = 0.d0
      Y(1:N) = X(1:N)
      do c = 1, NCOMP
         n1 = BPTR(c + 1) - BPTR(c)
         do ii = 1, n1
#ifdef CMPI
            if (.not. RESNODE(BNODE(BPTR(c) + ii))) cycle
#endif
            t = 0.d0
            do jj = 1, n1
               t = t + BMAT(BMPTR(c) + (ii - 1)*n1 + jj)*X(BNODE(BPTR(c) + jj))
            end do
            Y(BNODE(BPTR(c) + ii)) = t
         end do
      end do
#ifdef CMPI
      if (BMODE > 0) then
         call UPDATER_SUBSET(N, Y, 1, dum, 0)
      else
         do c = 1, NGROW
            t = 0.d0
            do jj = BGPTR(c) + 1, BGPTR(c + 1)
               t = t + BGVAL(jj)*X(BGIDX(jj))
            end do
            Y(BGNODE(c)) = t
         end do
      end if
#endif
   end subroutine WEIR_BLOCKS_HALF

   !-----------------------------------------------------------------------
   !  WEIR_BLOCKS_LOCAL: C^-1/2 is right at ghost nodes without an exchange
   !  (BMODE 0), so that JCG's iteration can use WEIR_BLOCKS_BT and _DT
   !-----------------------------------------------------------------------
   logical function WEIR_BLOCKS_LOCAL()
      WEIR_BLOCKS_LOCAL = BMODE == 0
   end function WEIR_BLOCKS_LOCAL

   !-----------------------------------------------------------------------
   !  WEIR_BLOCKS_BT, WEIR_BLOCKS_DT: JCG's DTWD = D - C^-1/2 A C^-1/2 D,
   !  with work only at the blocks' nodes, where C^-1/2 is not the identity:
   !     call WEIR_BLOCKS_BT(D)           D = C^-1/2 D (saving D)
   !     DTWD = D - A D                   JCG's product, then its exchange
   !     call WEIR_BLOCKS_DT(D, DTWD)     D back, and DTWD at the blocks' nodes
   !  Elsewhere C^-1/2 A C^-1/2 D and A C^-1/2 D are the same. BMODE 0 only.
   !-----------------------------------------------------------------------
   subroutine WEIR_BLOCKS_BT(D)
#ifdef CMPI
      use MESSENGER, only: RESNODE
#endif
      implicit none
      real(8), intent(INOUT) :: D(:)
      integer :: c, ii, i1
      do c = 1, NCOMP
         do ii = 1, BPTR(c + 1) - BPTR(c)
            i1 = BNODE(BPTR(c) + ii)
#ifdef CMPI
            if (.not. RESNODE(i1)) cycle
#endif
            BWT(i1) = BLOCK_ROW(c, ii, D)
         end do
      end do
#ifdef CMPI
      do c = 1, NGROW
         BWT(BGNODE(c)) = GHOST_ROW(c, D)
      end do
#endif
      do c = 1, NCOMP
         do ii = 1, BPTR(c + 1) - BPTR(c)
            i1 = BNODE(BPTR(c) + ii)
#ifdef CMPI
            if (.not. RESNODE(i1)) cycle
#endif
            BWS(i1) = D(i1)
            D(i1) = BWT(i1)
         end do
      end do
#ifdef CMPI
      do c = 1, NGROW
         i1 = BGNODE(c)
         BWS(i1) = D(i1)
         D(i1) = BWT(i1)
      end do
#endif
   end subroutine WEIR_BLOCKS_BT

   subroutine WEIR_BLOCKS_DT(D, DTWD)
#ifdef CMPI
      use MESSENGER, only: RESNODE
#endif
      implicit none
      real(8), intent(INOUT) :: D(:), DTWD(:)
      integer :: c, ii, i1
      !...A C^-1/2 D at the rows' columns (all of them nodes of the blocks)
      do c = 1, NCOMP
         do ii = BPTR(c) + 1, BPTR(c + 1)
            BWY(BNODE(ii)) = D(BNODE(ii)) - DTWD(BNODE(ii))
         end do
      end do
#ifdef CMPI
      do c = 1, NGROW
         i1 = BGNODE(c)
         BWY(i1) = D(i1) - DTWD(i1)
         do ii = BGPTR(c) + 1, BGPTR(c + 1)
            BWY(BGIDX(ii)) = D(BGIDX(ii)) - DTWD(BGIDX(ii))
         end do
      end do
#endif
      do c = 1, NCOMP
         do ii = 1, BPTR(c + 1) - BPTR(c)
            i1 = BNODE(BPTR(c) + ii)
#ifdef CMPI
            if (.not. RESNODE(i1)) cycle
#endif
            D(i1) = BWS(i1)
            DTWD(i1) = D(i1) - BLOCK_ROW(c, ii, BWY)
         end do
      end do
#ifdef CMPI
      do c = 1, NGROW
         i1 = BGNODE(c)
         D(i1) = BWS(i1)
         DTWD(i1) = D(i1) - GHOST_ROW(c, BWY)
      end do
#endif
   end subroutine WEIR_BLOCKS_DT

   !...row ii of component c of C^-1/2, times X
   real(8) function BLOCK_ROW(C, II, X)
      integer, intent(IN) :: C, II
      real(8), intent(IN) :: X(:)
      integer :: jj, n1
      n1 = BPTR(C + 1) - BPTR(C)
      BLOCK_ROW = 0.d0
      do jj = 1, n1
         BLOCK_ROW = BLOCK_ROW + BMAT(BMPTR(C) + (II - 1)*n1 + jj)*X(BNODE(BPTR(C) + jj))
      end do
   end function BLOCK_ROW

#ifdef CMPI
   !...ghost row c of C^-1/2, times X
   real(8) function GHOST_ROW(C, X)
      integer, intent(IN) :: C
      real(8), intent(IN) :: X(:)
      integer :: jj
      GHOST_ROW = 0.d0
      do jj = BGPTR(C) + 1, BGPTR(C + 1)
         GHOST_ROW = GHOST_ROW + BGVAL(jj)*X(BGIDX(jj))
      end do
   end function GHOST_ROW
#endif

   !-----------------------------------------------------------------------
   !  WEIR_BLOCKS_FULL: Y = C X at the resident nodes
   !-----------------------------------------------------------------------
   subroutine WEIR_BLOCKS_FULL(N, X, Y)
      implicit none
      integer, intent(IN) :: N
      real(8), intent(IN) :: X(:)
      real(8), intent(INOUT) :: Y(:)
      integer :: c, ii, jj, n1
      real(8) :: t
      Y(1:N) = X(1:N)
      do c = 1, NCOMP
         n1 = BPTR(c + 1) - BPTR(c)
         do ii = 1, n1
            t = 0.d0
            do jj = 1, n1
               t = t + BC(BMPTR(c) + (ii - 1)*n1 + jj)*X(BNODE(BPTR(c) + jj))
            end do
            Y(BNODE(BPTR(c) + ii)) = t
         end do
      end do
   end subroutine WEIR_BLOCKS_FULL

end module WEIR_PAIRS
