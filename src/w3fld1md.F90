!/ ------------------------------------------------------------------- /
      Module W3FLD1MD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III      NOAA/NCEP/NOPP |
!/                  |           B. G. Reichl            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         22-Mar-2021 |
!/                  +-----------------------------------+
!/
!/    01-Jul-2013 : Origination.                        ( version 4.xx )
!/    18-Mar-2015 : Clean-up/prepare for distribution   ( version 5.12 )
!/    15-Jan-2016 : Updates                             ( version 5.12 )
!/                                                      ( B. G. Reichl )
!/    27-Jul-2016 : Added Charnock output  (J.Meixner)  ( version 5.12 )
!/    22-Jun-2018 : updated SIG2WN subroutine (X.Chen)  ( version 6.06 )
!/                  modified the range of wind profile computation;
!/                  corrected direction of the shortest waves
!/    22-Mar-2021 : Consider DAIR a variable            ( version 7.xx )
!/    xx-xxx-2025 : Edited by A.Papandreou for computing variables 
!/                  necessary for applying the URI tail to the pi term, 
!/                  instead of wind stress. Part of the code was taken
!/                  from Jie's code.
!/                  (Clean-up by S.Bunya, 02-Mar-2026)
!/
!/
!/    Copyright 2009 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS.
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     This Module contains routines to compute the wind stress vector
!     from the wave spectrum, the wind vector, and the lower atmospheric
!     stability (the form included here is for neutral conditions, but
!     the structure needed to include stability is included in comments).
!     The stress calculated via this subroutine is
!     intended for coupling to serve as the boundary condition
!     between the ocean and atmosphere, and (for now)
!     and has no impact on the wave spectrum calculated.
!     The calculation in w3fld1 is based on the method
!     presented in Reichl, Hara, and Ginis (2014), "Sea State Dependence
!     of the Wind Stress under Hurricane Winds."
!
!  2. Variables and types :
!
!     Not applicable.
!
!  3. Subroutines and functions :
!
!      Name       Type  Scope    Description
!     ----------------------------------------------------------------
!      W3FLD1     Subr. Public   Reichl et al. 2014 stress calculation
!      INFLD1     Subr. Public   Corresponding initialization routine.
!      APPENDTAIL Subr. Public   Modification of tail for calculation
!      SIG2WN     Subr. Public   Depth-dependent dispersion relation
!      WND2Z0M    Subr. Public   Wind to roughness length
!      W3FLD1MD_INIT Subr. Public   Initialization of tail parameters
!     ----------------------------------------------------------------
!
!  4. Subroutines and functions used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!     ----------------------------------------------------------------
!
!  5. Remarks :
!
!  6. Switches :
!
!     !/S  Enable subroutine tracing.
!     !/
!
!  7. Source code :
!/
!/ ------------------------------------------------------------------- /
!/
!
      PUBLIC

      INTEGER :: Tail_Choice ! Tail_Choice: Chose the method to 
                             ! determine the level of the tail
      REAL    :: Tail_Level  !if Tail_Choice=0, tail is constant
      REAL    :: Tail_transition_ratio1 ! freq/fpi where tail
                                        ! adjustment begins
      REAL    :: Tail_transition_ratio2 ! freq/fpi where tail
                                        ! adjustment ends
!/
      PRIVATE :: Tail_Choice, Tail_Level, &
                 Tail_transition_ratio1, Tail_transition_ratio2
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3FLD1MD_INIT(CHOICE, LEVEL, TRANS1, TRANS2)
!/
!/    02-Mar-2026 : Added initialization of parameters (S.Bunya)
!/
!   1. Purpose :
!     Initialize the parameters for the tail calculation
!   2. Method :
!     See source code.
!   3. Parameters :
!     ----------------------------------------------------------------
!       CHOICE  Integer, optional  I   Choice of method to determine 
!         the level of the tail. Default is 10 m wind level.
!         For the URI tail you need to set this parameter to 1. 
!         This will make the tail dependent on the wind speed magnitude, 
!         according to the lines at the end of the file. 
!         If you set it to 0, the tail will be constant for any wind speed.
!         Choices are:
!         0: Constant tail level
!         1: Wind speed magnitude dependent tail level
!       LEVEL   Real, optional     I   Level of the tail if Tail_Choice=0.
!         Default is 0.1. This value is tentative. This value should be 
!         set properly when CHOICE=0. Not used when CHOICE=1.
!       TRANS1  Real, optional     I   Frequency/fpi where tail adjustment
!         begins. Default is 1.25, adopted from Jie's code.
!       TRANS2  Real, optional     I   Frequency/fpi where tail adjustment
!         ends. Default is 3.0, adopted from Jie's code.
!     ----------------------------------------------------------------
!/ ------------------------------------------------------------------- /
      IMPLICIT NONE
      INTEGER, INTENT(IN), OPTIONAL :: CHOICE
      REAL, INTENT(IN), OPTIONAL :: LEVEL
      REAL, INTENT(IN), OPTIONAL :: TRANS1
      REAL, INTENT(IN), OPTIONAL :: TRANS2

      ! Initialize the parameters for the tail calculation
      if (present(CHOICE)) then
        Tail_Choice = CHOICE
      else
        Tail_Choice = 1
      endif
      if (present(LEVEL)) then
        Tail_Level = LEVEL
      else
        Tail_Level = 0.1
      endif
      if (present(TRANS1)) then
        Tail_transition_ratio1 = TRANS1
      else
        Tail_transition_ratio1 = 1.25
      endif
      if (present(TRANS2)) then
        Tail_transition_ratio2 = TRANS2
      else
        Tail_transition_ratio2 = 3.0
      endif

      RETURN
      
      END SUBROUTINE W3FLD1MD_INIT
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3FLD1( ASPC, FPI, WNDX, WNDY,               &
                    DEPTH, UST, USTD, TAUX, TAUY)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III      NOAA/NCEP/NOPP |
!/                  |           B. G. Reichl            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         22-Mar-2021 |
!/                  +-----------------------------------+
!/
!/    01-Jul-2013 : Origination.                        ( version 4.xx )
!/    18-Mar-2015 : Prepare for submission              ( version 5.12 )
!/    22-Mar-2021 : Consider DAIR a variable            ( version 7.xx )
!/
!  1. Purpose :
!
!     Diagnostic stress vector calculation from wave spectrum, lower
!     atmosphere stability, and wind vector (at some given height).
!     The height of wind vector is assumed to be within the constant
!     stress layer.  These parameterizations are meant to be performed
!     at wind speeds > 10 m/s, and may not converge for extremely young
!     seas (i.e. starting from flat sea conditions).
!
!  2. Method :
!     See Reichl et al. (2014).
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       ASPC    Real   I   1-D Wave action spectrum.
!       FPI     Real   I   Peak input frequency.
!       WNDX    Real   I   X-dir wind (assumed referenced to current)
!       WNDY    Real   I   Y-dir wind (assumed referenced to current)
!       ZWND    Real   I   Wind height.
!       DEPTH   Real   I   Water depth.
!       RIB     REAL   I   Bulk Richardson in lower atmosphere
!                          (for determining stability in ABL to get
!                          10 m neutral wind)
!       DAIR    REAL   I   Air density
!       TAUNUX  Real   0   X-dir viscous stress (guessed from prev.)
!       TAUNUY  Real   0   Y-dir viscous stress (guessed from prev.)
!       UST     Real   O   Friction velocity.
!       USTD    Real   O   Direction of friction velocity.
!       Z0      Real   O   Surface roughness length
!       CHARN   Real   O,optional    Charnock parameter
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3ASIM    Subr. W3ASIMMD Air-sea interface module.
!      W3EXPO    Subr.   N/A    Point output post-processor.
!      GXEXPO    Subr.   N/A    GrADS point output post-processor.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!       None.
!
!  7. Remarks :
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/S  Enable subroutine tracing.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!      USE CONSTANTS, ONLY: GRAV, DWAT, TPI, PI, KAPPA
!      USE W3GDATMD, ONLY: NK, NTH, NSPEC, SIG, DTH, XFR, TH
!      USE W3ODATMD, ONLY: NDSE
!      USE W3SERVMD, ONLY: EXTCDE
      USE SWCOMM3, ONLY       : GRAV, RHO, PWIND, PI2, PI, &
                                MSC, MDC, DDIR !, MMCGR
      USE M_GENARR, ONLY      : SPCDIR, SPCSIG !relative frequency in sigma-space!!!


!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL, INTENT(IN)        :: ASPC(MDC*MSC), WNDX, WNDY,   &
                                  DEPTH, FPI
      REAL, INTENT(OUT)       :: UST, USTD
!      REAL, INTENT(OUT), OPTIONAL :: CHARN
      REAL, INTENT(INOUT)     :: TAUX, TAUY
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      ! jie 08/2015 For coupling with SWAN
      REAL                    :: KAPPA
      REAL                    :: DWAT
      REAL                    :: DAIR
      REAL                    :: TPI
      INTEGER                 :: NK
      INTEGER                 :: NTH
      INTEGER                 :: NSPEC
      REAL                    :: DTH
      REAL                    :: XFR
      REAL, ALLOCATABLE, DIMENSION(:) :: SIG
      REAL, ALLOCATABLE, DIMENSION(:) :: TH      
      REAL, PARAMETER         ::  NU=0.105/10000.0
      REAL, PARAMETER         ::  DELTA=0.03
      ! Commonly used parameters
      REAL                    ::  wnd_in_mag, wnd_in_dir
      !For Calculating Tail
      REAL                    ::  KMAX, KTAILA, KTAILB, KTAILC
      REAL                    ::  SAT, z01, z02, u10
      LOGICAL                 ::  ITERFLAG
      INTEGER                 ::  COUNT
      !For Iterations
      REAL                    ::  DTX, DTY, iter_thresh, &
                                  USTSM, Z0SM, Z1
      !For stress calculation
      REAL                    ::  WAGE, CBETA, BP, CD,       &
                                  USTRB, ANGDIF, USTAR, ZNU, &
                                  TAUT, BETAG, TAUDIR, &
                                  TAUDIRB
      !For wind profile calculation
      REAL                    ::  UPROFV, VPROFV
      !For wind profile iteration
      REAL                    ::  WND_1X, WND_1Y, &
                                  WND_2X, WND_2Y, &
                                  WND_3X, WND_3Y, &
                                  DIFU10XX, DIFU10YX, DIFU10XY, DIFU10YY, &
                                  FD_A, FD_B, FD_C, FD_D, &
                                  DWNDX, DWNDY, &
                                  APAR, CH,UITV, VITV,USTL,&
                                  CK
      !For adding stability to wind profile
      REAL                    ::  WND_TOP, ANG_TOP, WND_PA, WND_PE,   &
                                  WND_PEx, WND_PEy, WND_PAx, WND_PAy, &
                                  CDM
      INTEGER                 ::  NKT, K, T, Z2, ITER, ZI, ZII, &
                                  I, CTR, ITERATION, KA1, KA2, &
                                  KA3, KB
      ! For defining extended spectrum with appended tail.
      REAL, ALLOCATABLE, DIMENSION(:)   :: WN, DWN, CP,SIG2
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SPC2
      REAL, ALLOCATABLE, DIMENSION(:)   :: TLTN, TLTE, TAUD, &
                                           TLTND, &
                                           TLTED, ZOFK, UPROF, VPROF, &
                                           FTILDE, UP1, VP1, UP, VP, &
                                           TLTNA, TLTEA
      LOGICAL                 :: FSFL1,FSFL2, CRIT1, CRIT2
      LOGICAL                 :: IT_FLAG1, IT_FLAG2
      LOGICAL, SAVE           :: FIRST = .TRUE.
!/
!/ ------------------------------------------------------------------- /
!/
!
! 0.  Initializations ------------------------------------------------ *
!
!     **********************************************************
!     ***    The initialization routine should include all   ***
!     *** initialization, including reading data from files. ***
!     **********************************************************
!
! jie initialize necessary varialbes and arrays
       KAPPA = PWIND(15)
       DWAT = PWIND(17) !RHO
       DAIR = PWIND(16)
       TPI = PI2
       NK = MSC
       NTH = MDC
       NSPEC = MDC*MSC !MMCGR
       DTH = DDIR

! jie  allocate necessary arrays
       allocate(SIG(NK))
       SIG = SPCSIG !relative frequency in sigma-space(radian)!!!
       XFR = SIG(NK)/SIG(NK-1) !should be constant
       allocate(TH(NTH))
       DO T = 1, NTH
        TH(T) = SPCDIR(T,1)
       ENDDO
!      IF ( FIRST ) THEN
!          CALL INFLD
!          FIRST  = .FALSE.
!       END IF
       wnd_in_mag = sqrt( wndx**2 + wndy**2 )
       wnd_in_dir = atan2(wndy, wndx)
       !----------------------------------------------------------+
       ! Assume wind input is neutral 10 m wind.  If wind input   +
       ! is not 10 m, tail level will need to be calculated based +
       ! on esimation of 10 m wind.                               +
       !----------------------------------------------------------+
       u10 = wnd_in_mag
       ! - Get tail level
       if (Tail_Choice.eq.0) then
          SAT=Tail_Level
       elseif (Tail_Choice.eq.1) then
          CALL WND2SAT(U10,SAT)
       endif
!
! 1.  Attach Tail ---------------------------------------------------- *
!
! If the depth remains constant, the allocation could be limited to the
!   first time step.  Since this code is designed for coupled
!   implementation where the water level can change, I keep it the
!   allocation on every time step.  When computational efficiency is
!   important, this process may be rethought.
!
      ! i. Find maximum wavenumber of input spectrum
      call sig2wn(sig(nk),depth,kmax)
      NKT = NK
      ! ii. Find additional wavenumber bins to extended to cm scale waves
      DO WHILE ( KMAX .LT. 366.0 )
        NKT = NKT + 1
        KMAX = ( KMAX * XFR**2 )
      ENDDO!K<366
      ! iii. Allocate new "extended" spectrum
      ALLOCATE( WN(NKT), DWN(NKT), CP(NKT), SIG2(NKT),SPC2(NKT,NTH), &
                TLTN(NKT), TLTE(NKT), TAUD(NKT), &
                TLTND(NKT), TLTED(NKT), ZOFK(NKT), UPROF(NKT+1),&
                VPROF(NKT+1), FTILDE(NKT), UP1(NKT+1),VP1(NKT+1), &
                UP(NKT+1), VP(NKT+1), TLTNA(NKT),TLTEA(NKT))
!
! 1a. Build Discrete Wavenumbers for defining extended spectrum on---- *
!
      !i. Copy existing sig to extended sig2, calculate phase speed.
      DO K = 1, NK !existing spectrum
         call sig2wn(sig(k),depth,wn(k))
         CP(K) = ( SIG(K) / WN(K) )
         sig2(k) = sig(k)
      ENDDO!K
      !ii. Calculate extended sig2 and phase speed.
      DO K = ( NK + 1 ), ( NKT) !extension
         sig2(k) = sig2(k-1) *XFR
         call sig2wn(sig2(k),depth,wn(k))
         CP(K) = SIG2(K) / WN(K)
      ENDDO!K
      !iii. Calculate dk's for integrations.
      DO K = 1, NKT-1
        DWN(K) = WN(K+1) - WN(K)
      ENDDO
      DWN(NKT) = WN(NKT)*XFR**2 - WN(NKT)
!
! 1b. Attach initial tail--------------------------------------------- *
!
      !i. Convert action spectrum to variance spectrum
      !   SPC(k,theta) = A(k,theta) * sig(k)
      ! This could be redone for computational efficiency
      I=0
      DO K=1, NK
        DO T=1, NTH
          I = I + 1
          SPC2(K,T) = ASPC(I)  * SIG(K)
        ENDDO!T
      ENDDO!K
      !ii. Extend k^-3 tail to extended spectrum
      DO K=NK+1, NKT
        DO T=1, NTH
          SPC2(K,T)=SPC2(NK,T)*WN(NK)**3.0/WN(K)**(3.0)
        ENDDO!T
      ENDDO!K
!
! 1c. Calculate transitions for new (constant saturation ) tail ------ *
!
      !i. Find wavenumber for beginning spc level transition to tail
      call sig2wn (FPI*TPI*tail_transition_ratio1,depth,ktaila )
      !ii. Find wavenumber for ending spc level transition to tail
      call sig2wn (FPI*TPI*tail_transition_ratio2,depth,ktailb )
      !iii. Find wavenumber for ending spc direction transition to tail
      KTAILC= KTAILB * 2.0
      !iv. Find corresponding indices of wavenumber transitions
      KA1 = 2     ! Do not modify 1st wavenumber bin
      DO WHILE ( ( KTAILA .GE. WN(KA1) ) .AND. (KA1 .LT. NKT-6) )
        KA1 = KA1 + 1
      ENDDO
      KA2 = KA1+2
      DO WHILE ( ( KTAILB .GE. WN(KA2) ) .AND. (KA2 .LT. NKT-4) )
        KA2 = KA2 + 1
      ENDDO
      KA3 = KA2+2
      DO WHILE ( ( KTAILC .GE. WN(KA3)) .AND. (KA3 .LT. NKT-2) )
        KA3 = KA3 + 1
      ENDDO
      !v. Call subroutine to perform actually tail truncation
      ! only if there is some energy in spectrum
      CALL APPENDTAIL(SPC2, WN, NKT, KA1, KA2, KA3,&
           wnd_in_dir, SAT)
      ! Spectrum is now set for stress integration
      !
!
! 2.  Prepare for iterative calculation of wave-form stress----------- *
!   Edited by A.Papandreou in 2025 for computing variables necessary for
!   applying the URI tail to the pi term, instead of wind stress.  
!
      Z2 = NKT+1
      DO ZI=1, NKT
        Z2 = Z2 - 1
        TLTED(zi)=0.0
        TLTND(zi)=0.0
        DO T=1, NTH
          TLTND(zi) = tltnd(zi) + SIN(TH(T)) * SIG2(Z2) * WN(Z2) &
                      * 2 * SPC2(Z2,T) * dth
          TLTED(zi) = tlted(zi) + COS(TH(T)) * SIG2(Z2) * WN(Z2) &
                      * 2 * SPC2(Z2,T) * dth              
        IF (isnan(TLTND(zi))) THEN
        TLTND(zi) = 0.0
        ENDIF
        IF (isnan(TLTED(zi))) THEN
        TLTED(zi) = 0.0
        ENDIF
        ENDDO
        IF (ZI .EQ. 1) THEN
          TLTNA(ZI)=TLTND(ZI)*DWN(Z2)*0.5
          TLTEA(ZI)=TLTED(ZI)*DWN(Z2)*0.5
        ELSE
          TLTNA(ZI)=(TLTND(ZI)+TLTND(ZI-1))*0.5*DWN(Z2)
          TLTEA(ZI)=(TLTED(ZI)+TLTED(ZI-1))*0.5*DWN(Z2)
        ENDIF
        IF (ZI.GT.1) then
            TLTN(ZI)=TLTN(ZI-1)+TLTNA(ZI)
            TLTE(ZI)=TLTE(ZI-1)+TLTEA(ZI)
        else
            TLTN(ZI)=TLTNA(ZI)
            TLTE(ZI)=TLTEA(ZI)
        endif
      ENDDO
      TAUY=TLTN(NKT) !by NKT full stress is entirely
      TAUX=TLTE(NKT) !from turbulent stress
      USTAR=SQRT(TAUY**2.0+TAUX**2.0)
      TAUDIR=atan2(TAUY, TAUX)
!|---------------------------------------------------------------------|
      USTD = ATAN2(TAUY,TAUX)
      UST = SQRT( TAUX**2 + TAUY**2 )
      !Z0 = 0
      !CHARN = TAUX
      DEALLOCATE(WN, DWN, CP,SIG2, SPC2, TLTN, TLTE, TAUD, &
                 TLTND, TLTED, ZOFK, UPROF, &
                 VPROF, FTILDE, UP1, VP1, UP, VP, TLTNA, TLTEA)
!/ End of W3FLD1 ----------------------------------------------------- /
!/
      RETURN
!
      END SUBROUTINE W3FLD1
!/ ------------------------------------------------------------------- /
!      SUBROUTINE INFLD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           B. G. Reichl            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         15-Jan-2016 |
!/                  +-----------------------------------+
!/
!/    15-Jan-2016 : Origination.                        ( version 5.12 )
!/
!  1. Purpose :
!
!     Initialization for w3fld1 (also used by w3fld2)
!
!  2. Method :
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3FLDX    Subr. W3FLDXMD Corresponding source term.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!       None.
!
!  7. Remarks :
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/S  Enable subroutine tracing.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!      USE W3ODATMD, ONLY: NDSE
!      USE W3GDATMD, ONLY: TAIL_ID, TAIL_LEV, TAIL_TRAN1, TAIL_TRAN2
!      USE W3SERVMD, ONLY: EXTCDE
!/
!      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
!/
!/ ------------------------------------------------------------------- /
!/
!
! 1.  .... ----------------------------------------------------------- *
!
!      Tail_Choice=Tail_ID
!      Tail_Level=TAIL_Lev
!      Tail_transition_ratio1 = TAIL_TRAN1
!      Tail_transition_ratio2 = TAIL_TRAN2
 
!
!      RETURN
!
! Formats
!
 
!/
!/ End of INFLD1 ----------------------------------------------------- /
!/
!      END SUBROUTINE INFLD
!/
!/ ------------------------------------------------------------------- /
      SUBROUTINE APPENDTAIL(INSPC, WN2, NKT, KA1, KA2, KA3, WNDDIR,SAT)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           B. G. Reichl            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         15-Jan-2016 |
!/                  +-----------------------------------+
!/
!/    15-Jan-2016 : Origination.                        ( version 5.12 )
!/
!  1. Purpose :
!
!     Set tail for stress calculation.
!
!  2. Method :
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3FLD1    Subr. W3FLD1MD Corresponding source term.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!       None.
!
!  7. Remarks :
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/S  Enable subroutine tracing.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!      USE CONSTANTS, ONLY: TPI, PI
!      USE W3GDATMD, ONLY: NTH, TH, DTH
!      USE W3ODATMD, ONLY: NDSE
!      USE W3SERVMD, ONLY: EXTCDE
      USE SWCOMM3, ONLY       : PI2, PI, &
                                MDC, DDIR
      USE M_GENARR, ONLY      : SPCDIR
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN) :: NKT, KA1, KA2, KA3
      REAL, INTENT(IN)    :: WN2(NKT), WNDDIR,SAT
      REAL, INTENT(INOUT)   :: INSPC(NKT,MDC)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      ! jie 08/2015 For coupling with SWAN
      REAL                    :: TPI
      INTEGER                 :: NTH
      REAL                    :: DTH
      REAL, ALLOCATABLE, DIMENSION(:) :: TH

      REAL                :: BT(NKT), IC, ANGLE2, ANG(NKT),&
                             NORMSPC(MDC), AVG, ANGDIF, M, MAXANG, &
                             MAXAN, MINAN
      INTEGER             :: MAI, I, K, T
      REAL, ALLOCATABLE, DIMENSION(:)  :: ANGLE1
!/
!/ ------------------------------------------------------------------- /
!/
!
! 1.  .... ----------------------------------------------------------- *
!
      !|###############################################################|
      !|##1. Get the level of the saturation spectrum in transition
      !|##   region A
      !|###############################################################|
      !-------------------------------------------
      ! 1a, get saturation level at KA1 (1.25xFPI)
      !-------------------------------------------
! jie initialize necessary varialbes and arrays
       TPI = PI2
       NTH = MDC
       DTH = DDIR
! jie  allocate necessary arrays
       allocate(TH(NTH))
       DO T = 1, NTH
        TH(T) = SPCDIR(T,1)
       ENDDO

      BT(KA1) = 0
      ANG = 0.0
      DO T=1, NTH
        BT(KA1)=BT(KA1)+INSPC(KA1,T)*WN2(KA1)**3.0*DTH
      ENDDO
      !-----------------------------------------------
      ! 1b, Set saturation level at KA2 (3xFPI) to SAT
      !-----------------------------------------------
      BT(KA2) = SAT
      !-------------------------------------------------------------
      ! 1c, Find slope of saturation spectrum in transition region A
      !-------------------------------------------------------------
      M = ( BT(KA2) - BT(KA1) ) / ( WN2(KA2) - WN2(KA1) )
      !----------------------------------------------------------------
      ! 1d, Find intercept of saturation spectrum in transition region
      !     A
      !----------------------------------------------------------------
      IC = BT(KA1) - M * WN2(KA1)
      !------------------------------------------------------
      ! 1e, Calculate saturation level for all wavenumbers in
      !     transition region A
      !------------------------------------------------------
      DO K=KA1,KA2
        BT(K)=M*WN2(K)+IC
      ENDDO
      !|###############################################################|
      !|##2. Determine the directionality at each wavenumber in
      !|##   transition region B
      !|###############################################################|
      !-----------------------------------------------
      ! 2a, Find angle of spectral peak at KA2 (3xFPI)
      !-----------------------------------------------
      MAXANG = 0.0
      DO T=1, NTH
        IF (INSPC(KA2,T) .GT. MAXANG) THEN
          MAXANG=INSPC(KA2,T)
        ENDIF
      ENDDO
      !-------------------------------
      ! 2b, Check if peak spans 2 bins
      !-------------------------------
      !MAI = total number of angles of peak (if it spans more than 1)
      MAI = 0
      DO T=1, NTH
        IF (MAXANG .EQ. INSPC(KA2,T)) THEN
          MAI = MAI+1
        ENDIF
      ENDDO
      !ANGLE1 = angles that correspond to peak (array)
      MAI = MAX(1,MAI)
      ALLOCATE(ANGLE1(MAI))
      !-----------------------------------------------------
      ! 2c, If peak spans 2 or more bins it must be averaged
      !-----------------------------------------------------
      K=1
      DO T=1, NTH
        IF (MAXANG .EQ. INSPC(KA2,T)) THEN
          ANGLE1(K) = TH(T)
          K=K+1
        ENDIF
      ENDDO
      DO K=1, MAI
        DO WHILE (ANGLE1(K) .LT. 0.0)
          ANGLE1(K) = ANGLE1(K) + TPI
        ENDDO
        DO WHILE (ANGLE1(K) .GE. TPI)
          ANGLE1(K) = ANGLE1(K) - TPI
        ENDDO
      ENDDO
      IF (MAI .GT. 1) THEN
        MAXAN = ANGLE1(1)
        MINAN = ANGLE1(1)
        DO I=2, MAI
          IF (MAXAN .LT. ANGLE1(I) )THEN
            MAXAN = ANGLE1(I)
          ENDIF
          IF (MINAN .GT. ANGLE1(I) )THEN
            MINAN = ANGLE1(I)
          ENDIF
        ENDDO
      !------------------------------------------------------
      !  Need to distinguish if mean cross the origin (0/2pi)
      !------------------------------------------------------
        IF (MAXAN-MINAN .GT. PI) THEN
          DO I=1, MAI
            IF (MAXAN - ANGLE1(I) .GT. PI) THEN
              ANGLE1(I) = ANGLE1(I) + TPI
            ENDIF
          ENDDO
          ANGLE2=SUM(ANGLE1)/MAX(REAL(MAI),1.)
        ELSE
           ANGLE2=SUM(ANGLE1)/MAX(REAL(MAI),1.)
        ENDIF
      ELSE
        ANGLE2=ANGLE1(1)
      ENDIF
      DO WHILE (ANGLE2 .LT. 0.0)
        ANGLE2 = ANGLE2 + TPI
      ENDDO
      DO WHILE (ANGLE2 .GE. TPI)
        ANGLE2 = ANGLE2 - TPI
      ENDDO
      !
      !---------------------------------------------------
      ! This deals with angles that are less than 90
      !---------------------------------------------------
      if (cos(angle2-wnddir) .ge. 0.) then  !Less than 90
        m=asin(sin(wnddir-angle2))/(wn2(ka3)-wn2(ka2))
        ic=angle2
        do k=ka2, ka3
          ang(k)=ic +m*(wn2(k)-wn2(ka2))
        enddo
      else
      !----------------------------------------------------
      !  This deals with angels that turn clockwise
      !----------------------------------------------------
        if (sin(wnddir-angle2).GE.0) then
         m=acos(cos(wnddir-angle2))/(wn2(ka3)-wn2(ka2))
          ic=angle2
          do k=ka2, ka3
            ang(k)=ic+m*(wn2(k)-wn2(ka2))
          enddo
        else
      !-----------------------------------------------------
      !  This deals with angels that cross counter-clockwise
      !-----------------------------------------------------
          m=acos(cos(wnddir-angle2))/(wn2(ka3)-wn2(ka2))
          ic=angle2
          do k=ka2, ka3
            ang(k)=ic-m*(wn2(k)-wn2(ka2))
          enddo
        endif
      endif
      !----------------------------------------------
      ! Region A, Saturation level decreased linearly
      ! while direction is maintained
      !----------------------------------------------
      DO K=KA1, KA2-1
         AVG=SUM(INSPC(K,:))/MAX(REAL(NTH),1.)
         DO T=1,NTH
            INSPC(K,T)=BT(K)*INSPC(K,T)/TPI/(WN2(K)**3.0)/AVG
         ENDDO
      ENDDO
      !-----------------------------------------------------------
      ! Region B, Saturation level left flat while spectrum turned
      ! to direction of wind
      !-----------------------------------------------------------
      DO K = KA2, KA3
        DO T=1, NTH
          angdif=th(t)-ang(k)
          IF (COS(ANGDIF) .GT. 0.0) THEN
            NORMSPC(T) = COS(ANGDIF)**2.0
          ELSE
            NORMSPC(T)=0.0
          ENDIF
       ENDDO
        AVG=SUM(NORMSPC)/MAX(REAL(NTH),1.)
        DO T=1, NTH
          INSPC(K,T) = SAT * NORMSPC(T)/TPI/(WN2(K)**3.0)/AVG
        ENDDO
      ENDDO
      DO T=1, NTH
        angdif=th(t)-wnddir
        IF (COS(ANGDIF) .GT. 0.0) THEN
          NORMSPC(T) = COS(ANGDIF)**2.0
        ELSE
          NORMSPC(T) = 0.0
        ENDIF
      ENDDO
      AVG=SUM(NORMSPC)/MAX(REAL(NTH),1.)!1./4.
      DO K=KA3+1, NKT
        DO T=1, NTH
         INSPC(K,T)=NORMSPC(T)*(SAT)/TPI/(WN2(K)**3.0)/AVG
        ENDDO
      ENDDO
      DEALLOCATE(ANGLE1)
!
! Formats
!
!/
!/ End of APPENDTAIL ----------------------------------------------------- /
!/
      RETURN
!
      END SUBROUTINE APPENDTAIL
!/ ------------------------------------------------------------------- /
!/
!/ ------------------------------------------------------------------- /
      SUBROUTINE SIG2WN(SIG,DEPTH,WN)
!/ ------------------------------------------------------------------- /
!Author: Brandon Reichl (GSO/URI)
!Origination  : 2013
!Update       : March - 18 - 2015
!             : June -22 -2018  (XYC)
!             : February -27 -2026 (sb) Hardened the Newton's method 
!                 iteration by introducing a few parameters and 
!                 checking for NaN values
!Puropse      : Convert from angular frequency to wavenumber
!               using full gravity wave dispersion relation
!               if tanh(kh)<0.99, otherwise uses deep-water
!               approximation.
!NOTE: May be a better version internal to WW3 that can replace this.
!      Improved by using newton's method for iteration.(2018)
!/ ------------------------------------------------------------------- /
!/
!        use constants, only: GRAV
        USE SWCOMM3, ONLY       : GRAV
!/
        implicit none
!/
        REAL,INTENT(IN)    :: SIG,DEPTH
        REAL,INTENT(OUT)   :: WN
!/
        real    :: wn1,wn2 !,sig1,sig2,dsigdk
        real    :: fk, fk_slp, denom
        integer :: iter
        real, parameter :: rel_tol = 1.0e-4
        real, parameter :: abs_tol = 1.0e-10
        real, parameter :: slope_eps = 1.0e-12
        integer, parameter :: max_iter = 100
!/ ------------------------------------------------------------------- /
        wn1=sig**2/GRAV
        wn2=wn1
!/ Updated code with Newton's method by XYC:
      if (wn1 .gt. 0.0 .and. depth .gt. 0.0 .and. tanh(wn1*depth) .LT. 0.99) then
         do iter=1,max_iter
           fk=GRAV*wn1*tanh(wn1*depth) - sig**2
           fk_slp = GRAV*tanh(wn1*depth) + GRAV*wn1*depth/(cosh(wn1*depth))**2

           if (abs(fk_slp) .lt. slope_eps .or. fk_slp .ne. fk_slp) exit

           wn2=wn1 - fk/fk_slp
           if (wn2 .ne. wn2 .or. wn2 .le. 0.0) exit

           denom = max(abs(wn1), abs_tol)
           if (abs(wn2-wn1)/denom .LT. rel_tol ) then
              wn1 = wn2
              exit
           endif

           wn1=wn2
         enddo
      endif
      WN=max(wn1,0.0)
!/ END of update
!/
!/ Previous code by BR:
!/ ------------------------------------------------------------------- /
!        wn1=sig**2/GRAV
!        wn2=wn1+0.00001
!        SWITCH=.true.
!/ ------------------------------------------------------------------- /
!        if (tanh(wn1*depth).LT.0.99) then
!           do i=1,5
!              if (SWITCH) then
!                 sig1=sqrt(GRAV*wn1*tanh(wn1*depth))
!                 sig2=sqrt(GRAV*wn2*tanh(wn2*depth))
!                 if (sig1.lt.sig*.99999.or.sig1.gt.sig*1.00001) then
!                    dsigdk=(sig2-sig1)/(wn2-wn1)
!                    WN1=WN1+(SIG2-SIG1)/dsigdk
!                    wn2=wn1+wn1*0.00001
!                 else
!                    SWITCH = .FALSE.
!                 endif
!              endif
!           enddo
!        endif
!/
!        WN=WN1
!/
        RETURN
      END SUBROUTINE SIG2WN
!/ ------------------------------------------------------------------- /
!/
!/ ------------------------------------------------------------------- /
      SUBROUTINE WND2Z0M( W10M , ZNOTM )
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           B. G. Reichl            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         04-Aug-2016 |
!/                  +-----------------------------------+
!/
!/    09-Apr-2014 : Last Update.                        ( version 5.12 )
!/    15-Aug-2016 : Updated for 2016 HWRF z0            ( J. Meixner   )
!/
!  1. Purpose :
!
!     Get bulk momentum z0 from 10-m wind.
!          Bulk stress corresponds to 2015 GFDL Hurricane model
!          Not published yet, but contact Brandon Reichl or
!          Isaac Ginis (Univ. of Rhode Island) for further info
!
!  2. Method :
!          This has now been updated for 2016 HWRF z0 using routines
!          from HWRF  znot_m_v1, Biju Thomas, 02/07/2014
!           and       znot_wind10m Weiguo Wang, 02/24/2016
!
!  3. Parameters :
!       Name  Unit  Type      Description
!     ----------------------------------------------------------------
!       W10M   m/s  input    10 m neutral wind [m/s]
!       ZNOTM  m    output   Roughness scale for momentum
!     ----------------------------------------------------------------
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3FLD1    Subr. W3FLD1MD Corresponding source term.
!      W3FLD2    Subr. W3FLD2MD Corresponding source term.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!       None.
!
!  7. Remarks :
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/S  Enable subroutine tracing.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
      REAL, INTENT(IN) :: W10M
      REAL, INTENT(OUT):: ZNOTM
 
      !Parameters from znot_m_v1
      REAL, PARAMETER :: bs0 = -8.367276172397277e-12
      REAL, PARAMETER :: bs1 = 1.7398510865876079e-09
      REAL, PARAMETER :: bs2 = -1.331896578363359e-07
      REAL, PARAMETER :: bs3 = 4.507055294438727e-06
      REAL, PARAMETER :: bs4 = -6.508676881906914e-05
      REAL, PARAMETER :: bs5 = 0.00044745137674732834
      REAL, PARAMETER :: bs6 = -0.0010745704660847233
 
      REAL, PARAMETER :: cf0 = 2.1151080765239772e-13
      REAL, PARAMETER :: cf1 = -3.2260663894433345e-11
      REAL, PARAMETER :: cf2 = -3.329705958751961e-10
      REAL, PARAMETER :: cf3 = 1.7648562021709124e-07
      REAL, PARAMETER :: cf4 = 7.107636825694182e-06
      REAL, PARAMETER :: cf5 = -0.0013914681964973246
      REAL, PARAMETER :: cf6 = 0.0406766967657759
 
      !Variables from znot_wind10m
      REAL            :: Z10, U10,AAA,TMP
 
 
      !Values as set in znot_wind10m
      Z10=10.0
      U10=W10M
      if (U10 > 85.0) U10=85.0
 
      !Calculation of z0 as in znot_m_v1
      IF ( U10 .LE. 5.0 ) THEN
        ZNOTM = (0.0185 / 9.8*(7.59e-4*U10**2+2.46e-2*U10)**2)
      ELSEIF (U10 .GT. 5.0 .AND. U10 .LT. 10.0) THEN
        ZNOTM =.00000235*(U10**2 - 25 ) + 3.805129199617346e-05
      ELSEIF ( U10 .GE. 10.0  .AND. U10 .LT. 60.0) THEN
        ZNOTM = bs6 + bs5*U10 + bs4*U10**2 + bs3*U10**3 + bs2*U10**4 +&
              bs1*U10**5 + bs0*U10**6
      ELSE
        ZNOTM = cf6 + cf5*U10 + cf4*U10**2 + cf3*U10**3 + cf2*U10**4 +&
              cf1*U10**5 + cf0*U10**6
      END IF
 
      !Modifications as in znot_wind10m for icoef_sf=4
 
      !for wind<20, cd similar to icoef=2 at 10m, then reduced
      TMP=0.4*0.4/(ALOG(10.0/ZNOTM))**2   ! cd at zlev
      AAA=0.75
      IF (U10 < 20) THEN
        AAA=0.99
      ELSEIF(U10 < 45.0) THEN
        AAA=0.99+(U10-20)*(0.75-0.99)/(45.0-20.0)
      END IF
      ZNOTM=Z10/EXP( SQRT(0.4*0.4/(TMP*AAA)) )
 
      END SUBROUTINE WND2Z0M
!/ ------------------------------------------------------------------- /
      SUBROUTINE WND2SAT(WND10,SAT)
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           B. G. Reichl            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         04-Aug-2016 |
!/                  +-----------------------------------+
!/
!/    15-Jan-2016 : Origination.                        ( version 5.12 )
!/    04-Aug-2016 : Updated for 2016 HWRF CD/U10 curve  ( J. Meixner   )
!/
!  1. Purpose :
!
!     Gives level of saturation spectrum to produce
!         equivalent Cd as in wnd2z0m (for neutral 10m wind)
!         tuned for method of Reichl et al. 2014
!
!  2. Method :
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
      !Input: WND10 - 10 m neutral wind [m/s]
      !Output: SAT  - Level 1-d saturation spectrum in tail [non-dim]
!     ----------------------------------------------------------------
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3FLD1    Subr. W3FLD1MD Corresponding source term.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!       None.
!
!  7. Remarks :
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/S  Enable subroutine tracing.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
        IMPLICIT NONE
!/
        REAL, INTENT(IN) :: WND10
        REAL, INTENT(OUT) :: SAT
!
!/ Old HWRF 2015 and ST2
!        SAT =0.000000000001237 * WND10**6 +&
!             -0.000000000364155 * WND10**5 +&
!             0.000000037435015 * WND10**4 +&
!             -0.000001424719473 * WND10**3 +&
!             0.000000471570975 * WND10**2 +&
!             0.000778467452178 * WND10**1 +&
!             0.002962335390055
!
!     SAT values based on
!     HWRF 2016 CD curve, created with  fetch limited cases ST4 physics
!        IF (WND10<20.0) THEN
!          SAT = -0.000018541921682*WND10**2  &
!                +0.000751515452434*WND10     &
!                +0.002466529381004
!        ELSEIF (WND10<45) THEN
!          SAT = -0.000000009060349*WND10**4  &
!                +0.000001276678367*WND10**3  &
!                -0.000068274393789*WND10**2  &
!                +0.001418180888868*WND10     &
!                +0.000262277682984
!        ELSE
!          SAT = -0.000155976275073*WND10     &
!                +0.012027763023184
!        ENDIF
!
!        SAT = min(max(SAT,0.002),0.014)
      IF (WND10<54.0) THEN
         SAT = -0.0000000000933045564 * WND10**4  &
               +0.000000368296046  * WND10**3  &
               -0.0000418658707 * WND10**2  &
               +0.00130595281  * WND10     &
               -0.000527661260
      ELSE
        SAT = 0.0000000742701634 * WND10**3  &
              -0.0000100486405 * WND10**2  &
              +0.000393347926 * WND10     &
              +0.00143171024
      ENDIF
      END SUBROUTINE WND2SAT
!
!/ End of module W3FLD1MD -------------------------------------------- /
!/
      END MODULE W3FLD1MD
