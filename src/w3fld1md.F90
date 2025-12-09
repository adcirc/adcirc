!/ ------------------------------------------------------------------- /
      Module W3FLD1MD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III      NOAA/NCEP/NOPP |
!/                  |           B. G. Reichl            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         18-Mar-2015 |
!/                  +-----------------------------------+
!/
!/    01-Jul-2013 : Origination.                        ( version 4.xx )
!/    18-Mar-2015 : Clean-up/prepare for distribution   ( version 5.xx )
!/                                                      ( B. G. Reichl )
!/
!/    Copyright 2009 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS.
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     This Module computes the wind stress vector from the wave
!     spectrum, the wind vector, and the lower atmospheric
!     stability.  The stress calculated via this subroutine is
!     intended for coupling to serve as the boundary condition
!     between the ocean and atmosphere, and (for now)
!     and has no impact on the wave spectrum calculated.
!     The calculation in w3fld1 is based on the method
!     presented in Reichl et al. (2014), "Sea State Dependence
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
!      APPENDTAIL Subr. Public  Modification of tail for calculation
!      SIG2WN     Subr. Public   Depth-dependent dispersion relation
!      WND2Z0M    Subr. Public   Wind to roughness length
!      MFLUX      Subr. Public   MO stability correction
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
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3FLD1( ASPC, FPI, WNDX,WNDY, ZWND, &
                         DEPTH, RIB, UST, USTD, Z0,TAUX,TAUY)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III      NOAA/NCEP/NOPP |
!/                  |           B. G. Reichl            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         18-Mar-2015 |
!/                  +-----------------------------------+
!/
!/    01-Jul-2013 : Origination.                        ( version 4.xx )
!/    18-Mar-2015 : Prepare for submission              ( version 5.xx )
!/
!  1. Purpose :
!
!     Diagnostic stress vector calculation from wave spectrum, lower
!     atmosphere stability, and wind vector (at some given height).
!     The height of wind vector is assumed to be within the constant
!     stress layer.
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
!       TAUX    Real   0   X-dir total stress (guessed from prev.)
!       TAUY    Real   0   Y-dir total stress (guessed from prev.)
!       UST     Real   O   Friction velocity.
!       USTD    Real   O   Direction of friction velocity.
!       Z0      Real   O   Surface roughness length
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
!     jie 08/2015 modify USE dependencies
!      USE CONSTANTS, ONLY: GRAV, DWAT, DAIR, TPI, PI, KAPPA
!      USE W3GDATMD, ONLY: NK, NTH, NSPEC, SIG, DTH, XFR, TH
      
!      USE GLOBAL, ONLY        : KAPPA !rhoAir,RHOWAT0
      USE SWCOMM3, ONLY       : GRAV, RHO, PWIND, PI2, PI, &
                                MSC, MDC, DDIR !, MMCGR
      USE M_GENARR, ONLY      : SPCDIR, SPCSIG !relative frequency in sigma-space!!!
      
!/      USE W3ODATMD, ONLY: NDSE
!/      USE W3SERVMD, ONLY: EXTCDE
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
!/    jie modified dependencies
      REAL, INTENT(IN)        :: ASPC(MDC*MSC), FPI, WNDX, WNDY,  &
                                 ZWND, DEPTH, RIB
      REAL, INTENT(OUT)       :: UST, USTD, Z0
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
      !parameters
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
                                  TAUT, TAUNUX, TAUNUY, BETAG
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
      LOGICAL                 :: FSFL1,FSFL2
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
!/      IF ( FIRST ) THEN
!/          CALL INFLD1
!/          FIRST  = .FALSE.
!/       END IF
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
       !debug only 
       !write(16,*) 'KAPPA=',KAPPA,'DWAT=',DWAT,'DAIR=',DAIR,'NSPEC=',NSPEC,'DTH=',DTH,'XFR=',XFR   
       !write(16,*) 'sig(1)=',SIG(1),'sig(NK)=',SIG(NK)!,'TH(1)=',TH(1),'TH(NTH)=',TH(NTH)
       !write(16,*) 'FPI = ',FPI, "FPI*2pi = ", FPI*TPI
                
       wnd_in_mag = sqrt( wndx**2 + wndy**2 )
       wnd_in_dir = atan2(wndy, wndx)
       !Get guess for 10m (use bulk neutral)
       IF (abs(zwnd-10.).gt.1.) THEN
          IterFLAG=.true.
          COUNT = 1 !COUNT is now counting iteration over z0
          do while(IterFLAG)
             u10=wnd_in_mag*log(10./z01)/log(zwnd/z01)
             CALL wnd2z0m(u10,z02)
             if ( (abs(z01/z02-1.).GT.0.001) .AND. &
                  (COUNT.LT.10))THEN
                z01 = z02
             else
                IterFLAG = .false.
             endif
             COUNT = COUNT + 1
          enddo
       ELSE
          u10 = wnd_in_mag
       ENDIF
       CALL WND2SAT(U10,SAT)
!
! 1.  Attach Tail ---------------------------------------------------- *
!
      ! i. Find maximum wavenumber of input spectrum
      call sig2wn(sig(nk),depth,kmax)
      NKT = NK
      ! ii. Find additional wavenumber bins to extended to cm scale waves
      DO WHILE ( KMAX .LT. 366.0 )
        NKT = NKT + 1
        KMAX = ( KMAX * XFR**2 )
      ENDDO !K<366
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
      I=0
      DO K=1, NK
        DO T=1, NTH
          I = I + 1
          SPC2(K,T) = ASPC(I) * SIG(K)
        ENDDO !T
      ENDDO !K
      !ii. Extend k^-3 tail to extended spectrum
      DO K=NK+1, NKT
        DO T=1, NTH
          SPC2(K,T)=SPC2(NK,T)*WN(NK)**3.0/WN(K)**(3.0)
        ENDDO !T
      ENDDO !K
!
! 1c. Calculate transitions for new (constant saturation ) tail ------ *
!
      !i. Find wavenumber for beginning spc level transition to tail
      call sig2wn (FPI*TPI*1.25,depth,ktaila )
      !ii. Find wavenumber for ending spc level transition to tail
      call sig2wn ( FPI*TPI*3.0,depth,ktailb )
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
      CALL APPENDTAIL(SPC2, WN, NKT, KA1, KA2, KA3,&
                      wnd_in_dir, SAT)
      ! Spectrum is now set for stress integration
      ! Jie DEBUG only to be continued
      !write(16,*) "U10=",U10, "SAT=", SAT, "depth=", depth 
      !write(16,*) "NK=", NK, "NKT=", NKT
      !write(16,*) "KA1=", KA1,"KA2=", KA2,"KA3=", KA3
      !write(16,*) "KTAILA=",KTAILA,"KTAILB=",KTAILB
      !write(16,*) "SPC2(KA1,1:NTH):", SPC2(KA1,1:NTH)
      !write(16,*) "SPC2(KA2-1,1:NTH):", SPC2(KA2-1,1:NTH)

!
! 2.  Prepare for iterative calculation of wave-form stress----------- *
!
      DTX = 0.00005
      DTY = 0.00005
      iter_thresh = 0.001
!
! 2a. Calculate initial guess for viscous stress from smooth-law------ *
! (Would be preferable to use prev. step)
!
      Z0SM = 0.001  !Guess
      IT_FLAG1 = .true.
      ITERATION = 0
      DO WHILE( IT_FLAG1 )
        ITERATION = ITERATION + 1
        Z1 = Z0SM
        USTSM = KAPPA * wnd_in_mag / ( LOG( ZWND / Z1 ) )
        Z0SM = 0.132 * NU / USTSM
        IF ( (ABS( Z0SM - Z1 ) .LT. 10.0**(-6)) .OR.&
             ( ITERATION .GT. 5 )) THEN
          IT_FLAG1 = .false.
        ENDIF
      ENDDO
      ITERATION = 1
      ! Guessed values of viscous stress
      TAUNUX = USTSM**2 * DAIR * wndx / wnd_in_mag
      TAUNUY = USTSM**2 * DAIR * wndy / wnd_in_mag
!
! 3.  Enter iterative calculation of wave form/skin stress----------  *
!
      IT_FLAG1 = .true.
      DO WHILE (IT_FLAG1)
        DO ITER=1, 3 !3 loops for TAUNU iteration
          Z2 = NKT
          ! First : TAUNUX + DX
          IF (ITER .EQ. 1) THEN
            TAUNUX = TAUNUX + DTX
          ! Second : TAUNUY + DY
          ELSEIF (ITER .EQ. 2) THEN
            TAUNUX = TAUNUX - DTX
            TAUNUY = TAUNUY + DTY
          ! Third : unmodified
          ELSEIF (ITER .EQ. 3) THEN
            TAUNUY = TAUNUY - DTY
          ENDIF
          ! Near surface turbulent stress = taunu
          TLTN(1) = TAUNUY
          TLTE(1) = TAUNUX
!|---------------------------------------------------------------------|
!|-----Calculate first guess at growth rate and local turbulent stress-|
!|-----for integration as a function of wavedirection------------------|
!|---------------------------------------------------------------------|
          DO ZI = 2, NKT
            USTL=0.0
            TLTND(zi)=0.0
            TLTED(zi)=0.0
            Z2 = Z2 - 1
            ! Use value of prev. wavenumber/height
            TAUD(ZI) = ATAN2( TLTN(ZI-1), TLTE(ZI-1))
            USTL = SQRT( SQRT( TLTN(ZI-1)**2 + TLTE(ZI-1)**2 )/ DAIR )
            DO T = 1, NTH
              ANGDIF=TAUD(ZI)-TH(T) !stress/wave angle
              IF ( COS( ANGDIF ) .GE. 0.0 ) THEN !Waves aligned
                WAGE = CP(Z2) / (USTL)
                ! First, waves much slower than wind.
                IF ( WAGE .LT. 10 ) THEN
                  CBETA = 25.0
                ! Transition from waves slower than wind to faster
                ELSEIF ( ( WAGE .GE. 10.0 ) .AND. &
                          ( WAGE .LE. 25.0 ) ) THEN
                  CBETA = 10.0 + 15.0 * COS( PI * ( WAGE - 10.0 ) &
                          / 15.0 )
                ! Waves faster than wind
                ELSEIF ( WAGE .GT. 25.0 ) THEN
                  CBETA = -5.0
                ENDIF
              ! Waves opposing wind
              ELSE
                CBETA = -25.0
              ENDIF
              !Integrate turbulent stress !jie/backward wn^2??
              TLTND(ZI) =TLTND(ZI)+( SIN( TH(T) ) * COS( ANGDIF )**2)&
                          * CBETA * SPC2(Z2,T) * &
                          SQRT( TLTE(ZI-1)**2 + TLTN(ZI-1)**2.0 ) &
                          * ( WN(Z2)**2.0 )*DTH
              TLTED(ZI) = TLTED(ZI)+(COS( TH(T) ) * COS( ANGDIF )**2)&
                          * CBETA * SPC2(Z2,T) * &
                          SQRT( TLTE(ZI-1)**2 + TLTN(ZI-1)**2.0 ) &
                          * ( WN(Z2)**2.0 )*DTH
            ENDDO !T
!|---------------------------------------------------------------------|
!|-----Complete the integrations---------------------------------------|
!|---------------------------------------------------------------------|
            IF (ZI .EQ. 2) THEN
              !First turbulent stress bin above taunu
              TLTNA(ZI) = TLTND(ZI) * DWN(Z2) * 0.5
              TLTEA(ZI) = TLTED(ZI) * DWN(Z2) * 0.5
            ELSE
              TLTNA(ZI)=(TLTND(ZI)+TLTND(ZI-1))*0.5*DWN(Z2)
              TLTEA(ZI)=(TLTED(ZI)+TLTED(ZI-1))*0.5*DWN(Z2)
            ENDIF
            TLTN(ZI)=TLTN(ZI-1)+TLTNA(ZI)
            TLTE(ZI)=TLTE(ZI-1)+TLTEA(ZI)
          ENDDO !ZI
          TAUY=TLTN(NKT)
          TAUX=TLTE(NKT)
          ! This is the first guess at the stress.
!|---------------------------------------------------------------------|
!|----Iterate til convergence------------------------------------------|
!|---------------------------------------------------------------------|
          USTRB=SQRT(SQRT(TAUY**2.0+TAUX**2.0)/DAIR)
          IT_FLAG2 = .TRUE.
          CTR=1
          DO WHILE ( (IT_FLAG2) .AND. ( CTR .LT. 10 ) )
           Z2=NKT+1
            DO ZI=1, NKT
              Z2=Z2-1
              USTL=0.0
              TLTED(zi)=0.0
              TLTND(zi)=0.0
              FTILDE(z2)=0.0
              TAUD(ZI) = ATAN2(TLTN(ZI),TLTE(ZI))
              USTL = SQRT(SQRT(TLTN(ZI)**2+TLTE(ZI)**2)/DAIR)
              DO T=1, NTH
                BETAG=0.0
                ANGDIF = TAUD(ZI)-TH(T)
                IF ( COS( ANGDIF ) .GE. 0.0 ) THEN
                  WAGE = CP(Z2)  / (USTL)
                 IF ( WAGE .LT. 10 ) THEN
                    CBETA = 25.0
                 ELSEIF ( ( WAGE .GE. 10.0 ) .AND. &
                          ( WAGE .LE. 25.0 ) ) THEN
                   CBETA = 10.0 + 15.0 * COS( PI * ( WAGE - 10.0 ) &
                          / 15.0 )
                 ELSEIF ( WAGE .GT. 25.0 ) THEN
                   CBETA = -5.0
                  ENDIF
                ELSE
                  CBETA = -25.0
                ENDIF
                BP = SQRT( (COS( TH(T) ) * COS( ANGDIF )**2.0)**2.0 &
                     + (SIN( TH(T) ) * COS( ANGDIF )**2.0)**2.0 )
                BETAG=BP*CBETA*SQRT(TLTE(ZI)**2.0+TLTN(ZI)**2.0) &
                           /(DWAT)*SIG2(Z2)/CP(Z2)**2
                FTILDE(Z2) = FTILDE(Z2) + BETAG * DWAT * GRAV &
                           * SPC2(Z2,T) * DTH
                TLTND(zi) =tltnd(zi)+ (SIN( TH(T) ) * COS( ANGDIF )**2.0)&
                            * CBETA * SPC2(Z2,T) * SQRT( &
                            TLTE(ZI)**2.0 + TLTN(ZI)**2.0 ) * &
                            ( WN(Z2)**2.0 )*dth
                TLTED(zi) = tlted(zi)+(COS( TH(T) ) * COS( ANGDIF )**2.0)&
                            * CBETA * SPC2(Z2,T) * SQRT( &
                            TLTE(ZI)**2.0 + TLTN(ZI)**2.0 ) * &
                            ( WN(Z2)**2.0 )*dth
              ENDDO !T
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
                 TLTN(ZI)=TAUNUY+TLTNA(ZI)
                 TLTE(ZI)=TAUNUX+TLTEA(ZI)
              endif
            ENDDO !ZI
            TAUY=TLTN(NKT) !by NKT full stress is entirely
            TAUX=TLTE(NKT) !from turbulent stress
            TAUT=SQRT(TAUY**2.0+TAUX**2.0)
            USTAR=SQRT(SQRT(TAUY**2.0+TAUX**2.0)/DAIR)
            IF ((ABS(USTAR-USTRB)*100.0)/((USTAR+USTRB)*0.5) .GT. 0.1) THEN
               USTRB=USTAR
              CTR=CTR+1
            ELSE
              IT_FLAG2 = .FALSE.
            ENDIF
          ENDDO
          KB=1
          DO WHILE(((TLTN(KB)**2+TLTE(KB)**2)/(TAUX**2+TAUY**2)).LE. &
             .99)
             KB=KB+1
          ENDDO
!|---------------------------------------------------------------------|
!|----Now begin work on wind profile-----------------------------------|
!|---------------------------------------------------------------------|
          DO I=1,NKT
            ZOFK(I)=DELTA/WN(I)
          ENDDO
          ZNU=0.1 * 1.45E-5 / SQRT(SQRT(TAUNUX**2.0+TAUNUY**2.0)/DAIR)
          UPROF(1:NKT)=0.0
          VPROF(1:NKT)=0.0
          UPROFV=0.0
          VPROFV=0.0
          ZI=1
          Z2=NKT
          UP1(ZI) = ( ( ( WN(Z2)**2 / DELTA ) * FTILDE(z2) ) + &
                    ( DAIR / ( ZOFK(Z2) * KAPPA ) ) * ( SQRT( &
                    TLTN(ZI)**2 + TLTE(ZI)**2 ) / DAIR )**(3/2) ) &
                    * ( TLTE(ZI) ) / ( TLTE(ZI) * TAUX &
                    + TLTN(ZI) * TAUY )
          VP1(ZI) = ( ( ( WN(Z2)**2 / DELTA ) * FTILDE(z2) ) + &
                    ( DAIR / ( ZOFK(Z2) * KAPPA ) ) * ( SQRT ( &
                    TLTN(ZI)**2 + TLTE(ZI)**2 ) / DAIR )**(3/2) ) &
                    * ( TLTN(ZI) ) / ( TLTE(ZI) * TAUX &
                    + TLTN(ZI) * TAUY )
          UP(ZI) = UP1(ZI)
          VP(ZI) = VP1(ZI)
          UPROF(ZI) = DAIR / KAPPA * ( SQRT( TAUNUX**2.0 + TAUNUY**2.0 ) &
                      / DAIR )**(1.5) * ( TAUNUX / ( TAUX * &
                      TAUNUX + TAUY * TAUNUY ) ) * LOG( &
                      ZOFK(Z2) / ZNU )
          VPROF(ZI) = DAIR / KAPPA * ( SQRT( TAUNUX**2.0 + TAUNUY**2.0 ) &
                      / DAIR )**(1.5) * ( TAUNUY / ( TAUX * &
                      TAUNUX + TAUY * TAUNUY ) ) * LOG( &
                      ZOFK(Z2) / ZNU )
          DO ZI=2, KB
            Z2 = Z2 - 1
            UP1(ZI) = ( ( ( WN(Z2)**2.0 / DELTA ) * FTILDE(Z2) ) + &
                      ( DAIR / ( ZOFK(Z2) * KAPPA ) ) * ( SQRT( &
                      TLTN(ZI)**2.0 + TLTE(ZI)**2.0 ) / DAIR )**(1.5) ) &
                      * ( TLTE(ZI) ) / ( TLTE(ZI) * TAUX + &
                      TLTN(ZI) * TAUY )
            VP1(ZI) = ( ( ( WN(Z2)**2.0 / DELTA ) * FTILDE(Z2) ) + &
                      ( DAIR / ( ZOFK(Z2) * KAPPA ) ) * ( SQRT( &
                      TLTN(ZI)**2.0 + TLTE(ZI)**2.0 ) / DAIR )**(1.5) ) &
                      * ( TLTN(ZI) ) / ( TLTE(ZI) * TAUX + &
                      TLTN(ZI) * TAUY )
            UP(ZI) = UP1(ZI) * 0.5 + UP1(ZI-1) * 0.5
            VP(ZI) = VP1(ZI) * 0.5 + VP1(ZI-1) * 0.5
            UPROF(ZI) = UPROF(ZI-1) + UP(ZI) * ( ZOFK(Z2) - ZOFK(Z2+1) )
            VPROF(ZI) = VPROF(ZI-1) + VP(ZI) * ( ZOFK(Z2) - ZOFK(Z2+1) )
          ENDDO
!|---------------------------------------------------------------------|
!|----Iteration completion/checks--------------------------------------|
!|---------------------------------------------------------------------|
          !ZI = ( KB + 1 )
 
          UPROF(KB+1) = UPROF(KB) + ( SQRT( SQRT( TAUY**2.0 + &
                      TAUX**2.0 ) / DAIR ) ) / KAPPA * TAUX &
                      / SQRT( TAUY**2.0 +TAUX**2.0 ) * LOG( ZWND &
                      / ZOFK(Z2) )
          VPROF(KB+1) = VPROF(KB) + ( SQRT( SQRT( TAUY**2.0 + &
                      TAUX**2.0 ) / DAIR ) ) / KAPPA * TAUY &
                      / SQRT( TAUY**2.0 +TAUX**2.0 ) * LOG( ZWND &
                      / ZOFK(Z2) )
          !---------------------------
          !  Adding stability effects
          !----------------------------
          !Get Wind at top of wave boundary layer
          WND_TOP=SQRT(UPROF(KB)**2+VPROF(KB)**2)
          ! Get Wind Angle at top of wave boundary layer
          ANG_TOP=ATAN2(VPROF(KB),UPROF(KB))
          ! Stress and direction
          USTD = ATAN2(TAUY,TAUX)
          UST = SQRT( SQRT( TAUX**2 + TAUY**2 ) / DAIR)
          ! Calclate along (PA) and across (PE) wind components
          WND_PA=WND_TOP*COS(ANG_TOP-USTD)
          WND_PE=WND_TOP*SIN(ANG_TOP-USTD)
          !Calculate cartesion across wind
          WND_PEx=WND_PE*cos(ustd+pi/2.)
          WND_PEy=WND_PE*sin(ustd+pi/2.)
          ! Assume neutral inside WBL calculate Z0
          Z0=ZOFK(Z2)*EXP(-WND_PA*kappa/UST)
          ! Using that Z0 calculate
          CALL MFLUX(WND_IN_MAG*COS(WND_IN_DIR-USTD),ZWND,Z0,RIB,CDM)
          WND_PA=UST/SQRT(CDM)
          WND_PAx=WND_PA*cos(ustd)
          WND_PAy=WND_PA*sin(USTd)
          IF (ITER .EQ. 3) THEN
            WND_1X = WND_PAx+WND_PEx
            WND_1Y = WND_PAy+WND_PEy
          ELSEIF (ITER .EQ. 2) THEN
            WND_2X = WND_PAx+WND_PEx
            WND_2Y = WND_PAy+WND_PEy
          ELSEIF (ITER .EQ. 1) THEN
            WND_3X = WND_PAx+WND_PEx
            WND_3Y = WND_PAy+WND_PEy
          ENDIF
        ENDDO !ITER=1,3
        ITERATION = ITERATION + 1
        CALL APPENDTAIL(SPC2, WN, NKT, KA1, KA2, KA3,&
                      atan2(VPROF(KB),UPROF(KB)), SAT)
        DIFU10XX = WND_3X - WND_1X
        DIFU10YX = WND_3Y - WND_1Y
        DIFU10XY = WND_2X - WND_1X
        DIFU10YY = WND_2Y - WND_1Y
        FD_A = DIFU10XX / DTX
        FD_B = DIFU10XY / DTY
        FD_C = DIFU10YX / DTX
        FD_D = DIFU10YY / DTY
        DWNDX = - WNDX + WND_1X
        DWNDY = - WNDY + WND_1Y
        UITV = ABS( DWNDX )
        VITV = ABS( DWNDY )
        CH = SQRT( UITV**2.0 + VITV**2.0 )
        IF (CH .GT. 15.) THEN
          APAR = 0.5 / ( FD_A * FD_D - FD_B * FD_C )
        ELSE
          APAR = 1.0 / ( FD_A * FD_D - FD_B * FD_C )
        ENDIF
        CK=4.
        IF (((VITV/MAX(ABS(WNDY),CK) .GT. iter_thresh) .OR. &
             (UITV/MAX(ABS(WNDX),CK) .GT. iter_thresh)) .AND. &
             (ITERATION .LT. 2)) THEN
          TAUNUX = TAUNUX - APAR * ( FD_D * DWNDX - FD_B * DWNDY )
          TAUNUY = TAUNUY - APAR * ( -FD_C * DWNDX +FD_A * DWNDY )
        ELSEIF (((VITV/MAX(ABS(WNDY),CK) .GT. iter_thresh) .OR. &
             (UITV/MAX(ABS(WNDX),CK) .GT. iter_thresh)) .AND. &
             (ITERATION .LT. 24)) THEN
           iter_thresh = 0.001
           TAUNUX = TAUNUX - APAR * ( FD_D * DWNDX - FD_B * DWNDY )
           TAUNUY = TAUNUY - APAR * ( -FD_C * DWNDX +FD_A * DWNDY )
        ELSEIF (((VITV/MAX(ABS(WNDY),CK) .GT. iter_thresh) .OR. &
             (UITV/MAX(ABS(WNDX),CK) .GT. iter_thresh)) .AND. &
             (ITERATION .LT. 26)) THEN
           iter_thresh = 0.01
           TAUNUX = TAUNUX - APAR * ( FD_D * DWNDX - FD_B * DWNDY )
           TAUNUY = TAUNUY - APAR * ( -FD_C * DWNDX +FD_A * DWNDY )
        ELSEIF (((VITV/MAX(ABS(WNDY),CK) .GT. iter_thresh) .OR. &
             (UITV/MAX(ABS(WNDX),CK) .GT. iter_thresh)) .AND. &
             (ITERATION .LT. 30)) THEN
           iter_thresh = 0.05
           TAUNUX = TAUNUX - APAR * ( FD_D * DWNDX - FD_B * DWNDY )
           TAUNUY = TAUNUY - APAR * ( -FD_C * DWNDX +FD_A * DWNDY )
        ELSEIF (ITERATION .GE. 30) THEN
           write(16,*)'Attn: W3FLD1 not converged.'
           write(16,*)'      Wind (X/Y): ',WNDX,WNDY, 'At water depth=',DEPTH 
           IT_FLAG1 = .FALSE.
           UST=-999
           TAUNUX=0.
           TAUNUY=0.
           !jie 12/2015 in this case return zero stress
           !TAUX = 0.
           !TAUY = 0.
        ELSEIF (((VITV/MAX(ABS(WNDY),CK) .LT. iter_thresh) .AND.&
             (UITV/MAX(ABS(WNDX),CK) .LT. iter_thresh)) .AND. &
             (ITERATION .GE. 2)) THEN
           IT_FLAG1 = .FALSE.
        ENDIF
     ENDDO !IT_FLAG1
!|---------------------------------------------------------------------|
!|----Finish-----------------------------------------------------------|
!|---------------------------------------------------------------------|
      USTD = ATAN2(TAUY,TAUX)
      UST = SQRT( SQRT( TAUX**2 + TAUY**2 ) / DAIR)
      CD = UST**2 / wnd_in_mag**2
      !Jie modified range from [0.0005 0.01] to [0.0005 0.005]
      !FSFL1=.not.((CD .LT. 0.01).AND.(CD .GT. 0.0005))
      FSFL1=.not.((CD .LT. 0.005).AND.(CD .GT. 0.0005))
      FSFL2=.not.(cos(wnd_in_dir-ustd).GT.0.9)
      IF (FSFL1 .or. FSFL2) THEN
         !Fail safe to bulk
         write(16,*)'Attn: W3FLD1 failed (cd out of range or dir-ustd > 0.9), need to use bulk...'
         write(16,*)'      Wind/Cd that failed: ',wnd_in_mag,CD
         UST = wnd_in_mag*SQRT(0.0015)
         USTD = wnd_in_dir
         !jie 12/2015 in this case return zero stress
         TAUX = 0.
         TAUY = 0.
      ENDIF
      DEALLOCATE(WN, DWN, CP,SIG2, SPC2, TLTN, TLTE, TAUD, &
                 TLTND, TLTED, ZOFK, UPROF, &
                 VPROF, FTILDE, UP1, VP1, UP, VP, TLTNA, TLTEA, &
                 SIG, TH) !jie 08/2015
!/ End of W3FLD1 ----------------------------------------------------- /
!/
      RETURN
!
      END SUBROUTINE W3FLD1
!/ ------------------------------------------------------------------- /
!/      SUBROUTINE INFLD1
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         01-Jul-2006 |
!/                  +-----------------------------------+
!/
!/    01-Jul-2006 : Origination.                        ( version 3.09 )
!/
!  1. Purpose :
!
!     Initialization for source term routine.
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
!/      USE W3ODATMD, ONLY: NDSE
!/      USE W3SERVMD, ONLY: EXTCDE
!/
!/      IMPLICIT NONE
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
!/      RETURN
!
! Formats
!
 
!/
!/ End of INFLD1 ----------------------------------------------------- /
!/
!/      END SUBROUTINE INFLD1
!/
!/ ------------------------------------------------------------------- /
      SUBROUTINE APPENDTAIL(INSPC, WN2, NKT, KA1, KA2, KA3, WNDDIR,SAT)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         01-Jul-2006 |
!/                  +-----------------------------------+
!/
!/    01-Jul-2006 : Origination.                        ( version 3.09 )
!/
!  1. Purpose :
!
!     Initialization for source term routine.
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
!jie 08/2015 modified module dependencies
!      USE CONSTANTS, ONLY: TPI, PI
!      USE W3GDATMD, ONLY: NTH, TH, DTH
      USE SWCOMM3, ONLY       : PI2, PI, &
                                MDC, DDIR
      USE M_GENARR, ONLY      : SPCDIR    
            
!/      USE W3ODATMD, ONLY: NDSE
!/      USE W3SERVMD, ONLY: EXTCDE
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/    jie modified dependencies
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
      AVG=1./4.!SUM(NORMSPC)/MAX(REAL(NTH),1.)!1./4.
      DO K=KA3+1, NKT
        DO T=1, NTH
         INSPC(K,T)=NORMSPC(T)*(SAT)/TPI/(WN2(K)**3.0)/AVG
        ENDDO
      ENDDO
      DEALLOCATE(TH)
      DEALLOCATE(ANGLE1)
      !jie DEBUG only
      !write(16,*), 'BT(KA1)=',BT(KA1),'BT(KA2-1)=',BT(KA2-1)
      !write(16,*), 'BT(KA2)=',BT(KA2),'BT(KA3)=',BT(KA3)
      !write(16,*), 'BT(KA3+1)=',BT(KA3+1),'BT(NKT)=',BT(NKT)
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
!Puropse      : Convert from angular frequency to wavenumber
!               using full gravity wave dispersion relation
!               if tanh(kh)<0.99, otherwise uses deep-water
!               approximation.
!/ ------------------------------------------------------------------- /
!/
!       jie change dependencies
!       use constants, only: GRAV
        USE SWCOMM3, ONLY       : GRAV
!/
        implicit none
!/
        REAL,INTENT(IN)    :: SIG,DEPTH
        REAL,INTENT(OUT)   :: WN
!/
        real    :: wn1,wn2,sig1,sig2,dsigdk
        integer :: i
        logical :: SWITCH
!/ ------------------------------------------------------------------- /
        wn1=sig**2/GRAV
        wn2=wn1+0.00001
        SWITCH=.true.
!/ ------------------------------------------------------------------- /
        if (tanh(wn1*depth).LT.0.99) then
           do i=1,5
              if (SWITCH) then
                 sig1=sqrt(GRAV*wn1*tanh(wn1*depth))
                 sig2=sqrt(GRAV*wn2*tanh(wn2*depth))
                 if (sig1.lt.sig*.99999.or.sig1.gt.sig*1.00001) then
                    dsigdk=(sig2-sig1)/(wn2-wn1)
                    WN1=WN1+(SIG-SIG1)/dsigdk
                    wn2=wn1+wn1*0.00001
                 else
                    SWITCH = .FALSE.
                 endif
              endif
           enddo
        endif
!/
        WN=WN1
!/
        RETURN
      END SUBROUTINE SIG2WN
!/ ------------------------------------------------------------------- /
!/
!/ ------------------------------------------------------------------- /
      SUBROUTINE WND2Z0M( U10 , ZNOTM )
!/ ------------------------------------------------------------------- /
!/ Date: Feb - 04 - 2014
!/ Author: Biju Thomas URI-GSO
!/ Updated: Apr - 09 - 2014
!/ Author: Brandon Reichl URI-GSO
!/ Title: WND2Z0M
!/ Purpose: Get bulk momentum z0 from 10-m wind.
!/ Input:  U10   - 10-m wind (magnitude) in m/s
!/ Output: ZNOTM - momentum z0 in m
!/ ------------------------------------------------------------------- /
!/
!/
        IMPLICIT NONE
!/ ------------------------------------------------------------------- /
      REAL, INTENT(IN)   :: U10
      REAL, INTENT(OUT)  :: ZNOTM
!/ ------------------------------------------------------------------- /
      REAL, PARAMETER :: a1=1.044183210405817e-12
      REAL, PARAMETER :: a2=-5.707116220939218e-11
      REAL, PARAMETER :: a3=8.005722172810571e-10
      REAL, PARAMETER :: a4=6.322045801589353e-09
      REAL, PARAMETER :: a5=-2.422002988137712e-07
      REAL, PARAMETER :: a6=2.269200594753249e-06
      REAL, PARAMETER :: a7=-6.029592778169796e-06
      REAL, PARAMETER :: a8=8.882284703541603e-06
      REAL, PARAMETER :: a9=-2.371341185499601e-06
!/
      REAL, PARAMETER :: b1=1.814407011197660e-15
      REAL, PARAMETER :: b2=-1.602907562918788e-13
      REAL, PARAMETER :: b3=-3.351205313520358e-11
      REAL, PARAMETER :: b4= 6.036179295940524e-09
      REAL, PARAMETER :: b5=-3.725481686822030e-07
      REAL, PARAMETER :: b6= 1.059761705898929e-05
      REAL, PARAMETER :: b7=-1.375241830530252e-04
      REAL, PARAMETER :: b8= 8.538858261732818e-04
      REAL, PARAMETER :: b9=-1.936638976963742e-03
!/ ------------------------------------------------------------------- /
      IF ( U10 .LE. 0.4) THEN
         ZNOTM = 4E-7
      ELSEIF( U10 .LE. 9.3 ) THEN
         ZNOTM = A9 + A8*U10 + A7*U10**2 + A6*U10**3 + A5*U10**4 + &
              A4*U10**5 + A3*U10**6 + A2*U10**7 + A1*U10**8
      ELSEIF (U10 .LT. 60.0) THEN
         ZNOTM = B9 + B8*U10 + B7*U10**2 + B6*U10**3 + B5*U10**4 + &
              B4*U10**5 + B3*U10**6 +B2*U10**7 + B1*U10**8
      ELSE
        ZNOTM = 1.3025e-3
      END IF
!/
      RETURN
    END SUBROUTINE WND2Z0M
!/ ------------------------------------------------------------------- /
    SUBROUTINE WND2Z0T(U_IN,Z_IN,ZNOTT)
!/ ------------------------------------------------------------------- /
!/ Date: Apr - 09 - 2014
!/ Author: Biju Thomas
!/ Adapted into WW3 by: Brandon Reichl URI-GSO
!/ Title: WWwnd2z0t
!/ Purpose: Get bulk thermal z0 from Za-meter wind speed.
!/ Input: uu - wind in m/s
!         za - height of wind in m
!/ Output: znott - thermal z0 in m
      IMPLICIT NONE
!/
      REAL, INTENT(IN) :: U_IN,Z_IN
      REAL, INTENT(OUT):: ZNOTT
!/
      REAL :: U10, U35
      REAL :: Z10 = 10.0
      REAL :: Z35 = 35.0
      REAL :: Z0I,Z0F
!/ TR
      REAL, PARAMETER  :: tr0 = 6.451939325286488e-08
      REAL, PARAMETER  :: tr1 = -7.306388137342143e-07
      REAL, PARAMETER  :: tr2 = -1.3709065148333262e-05
      REAL, PARAMETER  :: tr3 = 0.00019109962089098182
!/ TO
      REAL, PARAMETER  :: to0 = 1.4379320027061375e-08
      REAL, PARAMETER  :: to1 = -2.0674525898850674e-07
      REAL, PARAMETER  :: to2 = -6.8950970846611e-06
      REAL, PARAMETER  :: to3 = 0.00012199648268521026
!/ TN
      REAL, PARAMETER  :: tn0 = 1.4023940955902878e-10
      REAL, PARAMETER  :: tn1 = -1.4752557214976321e-08
      REAL, PARAMETER  :: tn2 = 5.90998487691812e-07
      REAL, PARAMETER  :: tn3 = -1.0920804077770066e-05
      REAL, PARAMETER  :: tn4 = 8.898205876940546e-05
      REAL, PARAMETER  :: tn5 = -0.00021123340439418298
!/ TT
      REAL, PARAMETER  :: tt0 = 1.92409564131838e-12
      REAL, PARAMETER  :: tt1 = -5.765467086754962e-10
      REAL, PARAMETER  :: tt2 = 7.276979099726975e-08
      REAL, PARAMETER  :: tt3 = -5.002261599293387e-06
      REAL, PARAMETER  :: tt4 = 0.00020220445539973736
      REAL, PARAMETER  :: tt5 = -0.0048088230565883
      REAL, PARAMETER  :: tt6 = 0.0623468551971189
      REAL, PARAMETER  :: tt7 = -0.34019193746967424
!/ TA
      REAL, PARAMETER  :: ta0 = -1.7787470700719361e-10
      REAL, PARAMETER  :: ta1 = 4.4691736529848764e-08
      REAL, PARAMETER  :: ta2 = -3.0261975348463414e-06
      REAL, PARAMETER  :: ta3 = -0.00011680322286017206
      REAL, PARAMETER  :: ta4 = 0.024449377821884846
      REAL, PARAMETER  :: ta5 = -1.1228628619105638
      REAL, PARAMETER  :: ta6 = 17.358026773905973
!/
      INTEGER :: LOOP
!/ ------------------------------------------------------------------- /
      !covert from uref to u35
      ! Get guess at z0
      CALL WND2Z0M(U_IN,Z0F)
      Z0I=100.0
      !/ Quick and dirty iteration to get approx. u35
      DO LOOP=1,5
         IF(ABS(Z0F-Z0I)/Z0F .GT. 0.01 ) THEN
            U10 = LOG(Z10/Z0F) / log(Z_IN/Z0I) * U_IN
            Z0I=Z0F
            CALL WND2Z0M(U10,Z0F)
         ENDIF
      ENDDO
      U35 = LOG(Z35/Z0F) / log(Z_IN/Z0I) * U_IN
      IF ( U35 .LE. 7.0 ) THEN
         !Not replacing constants (so znott remains smooth)
         ZNOTT = (0.0185 / 9.8*(7.59E-4*U35**2+2.46E-2*U35)**2)
      ELSEIF ( U35  .GE. 7.0 .AND. U35 .LT. 12.5 ) THEN
         ZNOTT =  TR3 + TR2*U35 + TR1*U35**2 + TR0*U35**3
      ELSEIF ( U35  .GE. 12.5 .AND. U35 .LT. 15.0 ) THEN
         ZNOTT =  TO3 + TO2*U35 + TO1*U35**2 + TO0*U35**3
      ELSEIF ( U35 .GE. 15.0  .AND. U35 .LT. 30.0) THEN
         ZNOTT =  TN5 + TN4*U35 + TN3*U35**2 + TN2*U35**3 + TN1*U35**4 + &
              TN0*U35**5
      ELSEIF ( U35 .GE. 30.0  .AND. U35 .LT. 60.0) THEN
         ZNOTT = TT7 + TT6*U35 + TT5*U35**2  + TT4*U35**3 + TT3*U35**4 + &
              TT2*U35**5 + TT1*U35**6 + TT0*U35**7
      ELSE
         ZNOTT =  TA6 + TA5*U35 + TA4*U35**2  + TA3*U35**3 + TA2*U35**4 + &
              TA1*U35**5 + TA0*U35**6
      END IF
      return
    end subroutine wnd2z0t
!/ ------------------------------------------------------------------- /
    SUBROUTINE MFLUX(WND,ZH,Z0,RIB,CD)
!/
!/ DATE: 03/17/2015
!/ Imported by: Brandon Reichl (GSO)
!/ Purpose:
!/   Calculate stability dependent fluxes
!/   Needed when coupled to atmosphere.
!/ Inputs:
!/   WND - Wind speed (m/s)
!/   ZH  - Atm. height of WND (m)
!/   Z0  - Surface roughness length (m)
!/   RIB - Lower atm. bulk Richardson (non-dimensional)
!/ Outputs:
!/   CD  - Drag coefficient for ZH (Tau/(rho_a U_ZH^2))
!/ ------------------------------------------------------------------- /
!/
!     jie change dependencies 
!     USE CONSTANTS, only: GRAV, TPI, KAPPA
      !USE GLOBAL, ONLY        : KAPPA 
      USE SWCOMM3, ONLY       : GRAV, PI2, PWIND
!}
      IMPLICIT NONE
      REAL, INTENT(IN) :: WND,ZH,Z0,RIB
      REAL, INTENT(OUT):: CD
      real, parameter :: c=0.76
      real, parameter :: a=5.0
      REAL :: ZETA2, zetat, zeta0,ZETA
      REAL :: FS,FMZ, FM0, FHZ, FH0,FM,FH
      REAL :: ZETAo, ZETA2o, DZD, DZ,ZT,DT,CDN
      INTEGER :: I, IT
      REAL :: TPI, KAPPA
!/
      TPI = PI2
      KAPPA = PWIND(15)
      ! Use Rb to calculate zeta
      ZETA=kappa*Rib/0.03
      !Estimate stability functions
      ZETA0=ZETA/ZH*Z0
      call wnd2z0t(WND,ZH,ZT)
      ZETAT=ZETA/ZH*ZT/100.
      ZETAo=1.0e+10
      ZETA2o=0.0e+00
      IT=1
      IF (ABS(ZETA).LT.0.001) THEN
         IT=0
         FMZ = 1./kappa * log(ZH)
         FM0 = 1./kappa * log(Z0)
         FHZ = FMZ
         FH0 = 1./kappa * log(Z0)
         FM=FMZ-FM0
         FH=FHZ-FH0
      ENDIF
      DO I=1,30
         ZETA0=ZETA/ZH*Z0
         ZETAT=ZETA/ZH*Z0
         IF (IT.EQ.1) THEN
            IF (ZETA.GT.10.) THEN
               FMZ = 1./kappa*(c*ZETA+17.193+0.5*a-10.0*c)
               FHZ = FMZ
            ELSEIF((ZETA.GT..5).AND.(ZETA.LE.10.)) THEN
               FMZ = 1./kappa*(8.*log(zeta)+4.25*ZETA**(-1) &
                    - 0.5*ZETA**(-2) + a*.5- 1.648 )
               FHZ = FMZ
            ELSEIF((ZETA.GT.0).AND.(ZETA.LE..5))THEN
               FMZ = 1./kappa*( log(zeta) + a*zeta)
               FHZ = FMZ
            ELSEIF(ZETA.LT.0)THEN
               FMZ = 1./kappa* ( log(zh)-log(((1.-16.*zeta)**(1./2.)+1.)&
                    *((1-16*zeta)**(1./4.)+1.)**2)+2.*atan((1.-16.*zeta)**(1./4.)))
               FM0 =1./kappa* ( log(abs(z0))-log(((1.-16.*zeta0)**(1./2.)+1.)&
                    *((1-16*zeta0)**(1./4.)+1.)**2)+2.*atan((1.-16.*zeta0)**(1./4.)))
               FHZ = 1./kappa* (log(abs(zh))-2.*log(1.+(1.-16.*zeta)**(1./2.)) )
               FH0 = 1./kappa* (log(abs(zetat))-2.*log(1.+(1.-16.*zetat)**(1./2.)) )
               FS = (1.-16.*zeta0)**(-1./2.)
            ENDIF
            IF (ZETA0.GT.10.) THEN
               FM0 = 1./kappa*(c*ZETA0+17.193+0.5*a-10.0*c)
               FH0 = 1./kappa*(c*ZETAT+17.193+0.5*a-10.0*c)
               FS = c*ZETA0
            ELSEIF((ZETA0.GT..5).AND.(ZETA0.LE.10.)) THEN
               FM0 = 1./kappa*(8.*log(zeta0)+4.25*ZETA0**(-1) &
                    - 0.5*ZETA0**(-2) + a*.5- 1.648 )
               FH0 = 1./kappa*(8.*log(zetaT)+4.25*ZETAT**(-1.) &
                    - 0.5*ZETAT**(-2.) +a*.5 - 1.648 )
               FS = 8. - 4.25* zeta0**(-1.) +.5*zeta0**(-2)
            ELSEIF((ZETA0.GT.0).AND.(ZETA0.LE..5))THEN
               FM0 = 1./kappa*( log(zeta0) + a*zeta0)
               FH0 =1./kappa*( log(zetat) + a*zetat)
               FS = 1. + a*zeta0
            ENDIF
            FM=FMz-FM0
            FH=FHZ-FH0
            ZETA2=kappa*Rib*FM**2/(FH+2/kappa*FS)
            DZ=(ZETA-ZETA2)/ZETA
            if (.NOT.((abs(DZ).LT.1.0e-4).AND.(I.LT.30))) then
               DZD=(ZETA2-ZETA2o)/(ZETA-ZETAo)
               ZETAo=ZETA
               ZETA=(ZETA2-ZETAo*DZD)/(1.-DZD)
               ZETA2o=ZETA2
            elseif (abs(DZ).LE.1.0e-3) THEN
               it=0
            elseif (abs(DZ).GT.1.0e-3.AND.(I.EQ.30)) THEN
               ! If MFLUX fails use neutral value.
               print*,'MFLUX FAIL, use neutral '
               !print*,'Wind: ',WND,'RIB: ',rib
               FMZ = 1./kappa * log(ZH)
               FM0 = 1./kappa * log(Z0)
               FHZ = FMZ
               FH0 = 1./kappa * log(Z0)
               FM=FMZ-FM0
               FH=FHZ-FH0
            endif
         endif
         CD=1/FM**2
      ENDDO
      FMZ = 1./kappa * log(ZH)
      FM0 = 1./kappa * log(Z0)
      FM=FMZ-FM0
      CDN=1/FM**2
      IF (abs(CD/CDN-1.).gt.0.5) THEN
         CD=CDN
      ENDIF
      RETURN
    END SUBROUTINE MFLUX
    
    
    !jie to be continued / test different SAT values
    SUBROUTINE WND2SAT(WND10,SAT)
        IMPLICIT NONE
        REAL, INTENT(IN) :: WND10
        REAL, INTENT(OUT) :: SAT
 
        SAT =0.000000000001237 * WND10**6 +&
             -0.000000000364155 * WND10**5 +&
             0.000000037435015 * WND10**4 +&
             -0.000001424719473 * WND10**3 +&
             0.000000471570975 * WND10**2 +&
             0.000778467452178 * WND10**1 +&
             0.002962335390055
        SAT = min(max(SAT,0.002),0.014)
      END SUBROUTINE WND2SAT
 
 
!/ End of module W3FLD1MD -------------------------------------------- /
!/
      END MODULE W3FLD1MD
