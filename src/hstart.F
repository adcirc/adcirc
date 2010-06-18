C******************************************************************************
C PADCIRC VERSION 45.12 03/17/2006                                            *
C  last changes in this file VERSION 45.07                                    *
C                                                                             *
C  mod history                                                                *
C  v45.07     - 11/07/05 - jgf- from 45.06 - removed EP from hot start file   *
C  v44.02     - 09/30/03 - rl - from 44.01 - additional changes associated    *
C                                          with removing local variables from *
C                                          global module.                     *
C  v44.01     - 10/06/03 - rl - from 43.03 - changed names of load vectors in *
C                                          2D solution                        *
C  v43.03     - 05/20/03 - rl - from 43.02 - parallel wind stuff (m.brown)    *
C                                          output buffer flush (m.cobb)       *
C                                          3D fixes (k.dresback)              *
C                                          drop MNPROC in fort.15 (t.campbell)*
C                                          various bug fixes in RBCs          *
C                                          ZSURFBUOY/BCPG calc                *
C  v43.00a    - sum  /02 - tc - from 36.01 (3D) & 41.12? (2D), create F90/    *
C                                               parallel unified 2D/3D source *
C  v41.11     - 09/14/01 - rl - from 41.10 - added NWS = -2 capability        *
C  v41.09     - 06/30/01 - jw - from 41.08 - minor mods per vp version 41.05  *
C  v40.02m002 - 12/22 - jjw/vjp - Vic suggested to avoid compiler conflicts   *
C  v40.02m001 - 12/21 - jjw - add cross barrier pipes cjjwm001                *
C                                                                             *
C******************************************************************************
C 
      SUBROUTINE HOTSTART()
C     
C**************************************************************************
C     
C     HOT START PROGRAM SETUP ROUTINE 
C     
C**************************************************************************
C     
      USE GLOBAL
      USE GLOBAL_IO
      USE GLOBAL_3DVS, ONLY: I3DSD, I3DSV, I3DST, I3DGD, I3DGV, I3DGT, qsurf


      USE SIZES
      USE HARM
      USE WIND
      USE OWIWIND,ONLY : NWS12INIT,NWS12GET   ! sb46.28sb01 added 09/xx/2006
      USE RS2,ONLY : RS2INIT,RS2GET           ! sb46.28sb01 added 09/xx/2006
C   Start Addition for netCDF by MCF 5/18/08
#ifdef NETCDF
      USE NETCDFIO, ONLY :
     &     timenc, timenc_len, fileunit, irtype,
     &     xenc, yenc, xvnc, yvnc, xpnc, ypnc, xnc, ync, xwnc, ywnc,
     &     et00nc, UU00nc, VV00nc, RMP00nc, RMU00nc, RMV00nc, nbvvnc,
     &     eta22, uu22, vv22, pr22, WNCX22, WNCY22,
     &     eta11,etadd, nodecodenc, noffnc, uu2, vv2,
     &     nbounc,nvelnc,nvellnc,
     &     nopenc,netanc,nvdllnc,ibtypenc,
     &     writenc61, writenc62, writenc63, writenc64,
     &     writenc71, writenc72, writenc73, writenc74, readnc67
#endif
C   End Addition for netCDF by MCF 5/18/08
      USE NodalAttributes, ONLY : 
     &     nolibf, nwp, tau0, cf, eslm,                       ! for netcdf
     &     LoadDirEffRLen,
     &     ApplyDirectionalWindReduction,
     &     OutputTau0 !jgf47.06
      IMPLICIT NONE

      INTEGER IT
      INTEGER I, J              !local loop counters
      INTEGER NCyc, NA, InputFileFmtVn
C   Start Addition for netCDF by MCF 4/8/08
      INTEGER ncmode
C   End Addition for netCDF by MCF 4/8/08

C   kmd48.33bc - add variable for timestep change 
      INTEGER NTHS  

      REAL(SZ) HollandTime
      REAL(SZ) ArgT, ArgTP, ArgSAlt
      REAL(SZ) H2     
      REAL(SZ) CCSFEA
      REAL(SZ) QTRatio
      REAL(SZ) SAltMul, S2SFEA
      REAL(SZ) TPMul
      REAL(SZ) WTRatio, WindX, WindY, WindMag, WDragCo
      REAL(SZ) WindDragLimit    !jgf46.01 Cap on max wind drag coefficient
      REAL(8) Time
      REAL(SZ) PRDIFF,PRBCKGRND_MH2O     ! tcm v49.16 20100617 added

C     jgf46.08 Fine grained ramp functions
      REAL(SZ) RampExtFlux1   ! Ramp for external flux b.c.s @ITHS-1
      REAL(SZ) RampExtFlux2   ! Ramp for external flux b.c.s @ITHS
      REAL(SZ) RampIntFlux1   ! Ramp for internal flux b.c.s @ITHS-1
      REAL(SZ) RampIntFlux2   ! Ramp for internal flux b.c.s @ITHS
      REAL(SZ) RampTip2    ! Ramp for tidal potential @ITHS
      REAL(SZ) RampMete2   ! Ramp for wind and atmospheric pressure @ITHS 
      REAL(SZ) RampWRad2   ! Ramp for wave radiation stress @ITHS
C     jgf46.14 Check for file existence before attempting to open.
      LOGICAL FileFound    ! .True. if the the file exists
      INTEGER ErrorIO      ! zero if file opened successfully
      INTEGER NP_G_IN, NE_G_IN  ! Global and
      INTEGER NP_A_IN, NE_A_IN  ! active number nodes, elements.
      INTEGER IHOTNC            ! variable to keep -67 and -68 in parallel
C
C     jgf48.4628 Add capability to turn off solution if only met output
C     was requested.
      LOGICAL metOutputOn     ! .true. if met output was requested
      LOGICAL nonMetOutputOff ! .true. no other output was requested
C
C     jgf48.4636 Add capability to read min/max files upon hotstart, 
C     if they are available
      CHARACTER(len=12) minMaxFN ! name of the min or max file to read in
      INTEGER node               ! node number of the min/max data
      CHARACTER(len=80) skipline ! dummy variable for min/max header data
      INTEGER numLinesInMinMaxFile ! count the lines to report to log file
      LOGICAL tooFewMinMaxLines  ! true if couldn't read enough lines from file

C   kmd48.33bc - add variables for heat flux 
      INTEGER :: NOD
      REAL(SZ),ALLOCATABLE :: TMP(:,:)
      REAL(SZ) :: CD, CDQ, QWIND
      CHARACTER(80) :: CDUM80
      REAL(SZ), PARAMETER :: RHOAIR=1.225d0
      REAL(SZ) :: RNDAYHS
      metOutputOn = .FALSE.
      nonMetOutputOff = .FALSE.

C
C     ASSIGN value to NCMODE according to the number of processors for NETCDF - MCF
C
      IF (MNPROC.gt.1) ncmode=2
      IF (MNPROC.eq.1) ncmode=1

C...  
C......READ IN 2DDI HOT START INITIAL CONDITION OVER WHOLE DOMAIN
C......THIS FILE ALWAYS HAS A RECL=8 BECAUSE IT IS ASSUMED THAT THE HARMONIC
C......ANALYSIS IS ALWAYS DONE IN 64 BITS, EVEN ON A WORKSTATION
C...  

C     Determine if the hot start file exists.
      FileFound=.False.
!      IF(myproc.eq.0) PRINT *,"IHOT = ", IHOT, " ncmode = ",ncmode
      IF (IHOT.eq.67.OR.
     & (IHOT.eq.-67.and.ncmode.eq.2) ) THEN
!      IF(myproc.eq.0) PRINT *,"1- IHOT = ", IHOT, " ncmode = ",ncmode
         INQUIRE(FILE=TRIM(INPUTDIR)//'/'//'fort.67',EXIST=FileFound)
      ENDIF
      IF (IHOT.eq.68.OR.
     & ( IHOT.eq.-68.and.ncmode.eq.2) ) THEN
!      IF(myproc.eq.0) PRINT *,"2- IHOT = ", IHOT, " ncmode = ",ncmode
         INQUIRE(FILE=TRIM(INPUTDIR)//'/'//'fort.68',EXIST=FileFound)
      ENDIF

#ifdef NETCDF
      if (ncmode.eq.1.or.MNPROC.EQ.1) THEN
!      IF(myproc.eq.0) PRINT *,"======serial netcdf========",ncmode
!      IF(myproc.eq.0) PRINT *,"3- IHOT = ", IHOT, " ncmode = ",ncmode
      IF(IHOT.EQ.-67)THEN 
!      IF(myproc.eq.0) PRINT *,"4- IHOT = ", IHOT, " ncmode = ",ncmode
!        PRINT *, "OPEN NETCDF HOTSTART FILE: fort.67.nc"
!        INQUIRE(FILE=TRIM(INPUTDIR)//'/'//'fort.67.nc',EXIST=FileFound)
        INQUIRE(FILE='fort.67.nc',EXIST=FileFound)
      ELSEIF(IHOT.EQ.-68) THEN
!      IF(myproc.eq.0) PRINT *,"5- IHOT = ", IHOT, " ncmode = ",ncmode
!        PRINT *, "OPEN NETCDF HOTSTART FILE: fort.68.nc"
!        INQUIRE(FILE=TRIM(INPUTDIR)//'/'//'fort.68.nc',EXIST=FileFound)
        INQUIRE(FILE='fort.68.nc',EXIST=FileFound)
      ENDIF
      ENDIF
#endif

      IF (.not.FileFound) THEN
         WRITE(16,1001)         ! hot start file  
         WRITE(16,1011)         ! was not found.
         WRITE(16,9973)         ! execution terminated
         IF (NScreen.ne.0.and.MyProc.eq.0) THEN
            WRITE(ScreenUnit,1001)       
            WRITE(ScreenUnit,1011)       
            WRITE(16,9973)      ! execution terminated
         ENDIF
         STOP
      ENDIF
C
C     Now open the hot start (unit 67 or 68) file.

      IF (IHOT.EQ.67.OR.IHOT.EQ.68. or.
     &    (ncmode.eq.2.and.
     &    ( IHOT.EQ.-67.or. IHOT.EQ.-68 ) ) )THEN                                 ! start reading binary hot start file

      WRITE(16,240) 
 240  FORMAT(/,9X,'Hot start file was found.',
     &     ' Opening file.') 
!      IF(IHOT.EQ.67) OPEN(67,FILE=TRIM(INPUTDIR)//'/'//'fort.67',
!     &     ACCESS='DIRECT',RECL=8,IOSTAT=ErrorIO)
!      IF(IHOT.EQ.68) OPEN(68,FILE=TRIM(INPUTDIR)//'/'//'fort.68',
!     &     ACCESS='DIRECT',RECL=8,IOSTAT=ErrorIO)


      IF(IHOT.EQ.67.or.(ncmode.eq.2.and.IHOT.EQ.-67) ) THEN 
           OPEN(67,FILE=TRIM(INPUTDIR)//'/'//'fort.67',
     &     ACCESS='DIRECT',RECL=8,IOSTAT=ErrorIO)
      
      ELSEIF(IHOT.EQ.68.or.(ncmode.eq.2.and.IHOT.EQ.-68) ) THEN
           OPEN(68,FILE=TRIM(INPUTDIR)//'/'//'fort.68',
     &     ACCESS='DIRECT',RECL=8,IOSTAT=ErrorIO)
      ENDIF

      IF (ErrorIO.GT.0) THEN 
         WRITE(16,1001)         ! hot start file
         WRITE(16,1005)         ! exists but can't be opened
         WRITE(16,9973)         ! execution terminated
         IF (NScreen.ne.0.and.MyProc.eq.0) THEN
            WRITE(ScreenUnit,1001) 
            WRITE(ScreenUnit,1005)
            WRITE(ScreenUnit,9973) 
         ENDIF
         STOP                   ! We're toast.
      ENDIF
      
      ENDIF
C
 1001 FORMAT('ERROR: The hot start file')
 1011 FORMAT('was not found.') 
 1012 FORMAT('was a nonmatching version') 
 1005 FORMAT('exists but cannot be opened.')
 9973 FORMAT(/,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)

C   Start Addition for netCDF by MCF 5/8/08
C
C   READ NETCDF HOTSTART FILE IF SERIAL VERSION IS BEING USED

#ifdef NETCDF
      
      IF(ncmode.eq.1.and.(IHOT.EQ.-67.or.IHOT.eq.-68)) THEN                                              !serial
!      PRINT *, "hstart: IHOT",IHOT
!      PRINT *, "hstart: ncmode",ncmode
!                PRINT *, "MNP = ",MNP
!                PRINT *, "NP = ",NP
!                PRINT *, "MNE = ",MNE
!                PRINT *, "NE = ",NE
!                PRINT *," NBOU = ",NBOU
!                PRINT *," NOPE = ",NOPE
!                PRINT *," NVEL = ",NVEL
!                PRINT *," NETA = ",NETA
       netanc=neta
       nvelnc=nvel
!                PRINT *," netanc = ",netanc
!                PRINT *," nvelnc = ",nvelnc
       ENDIF         
      IF(ncmode.eq.1) THEN                                              !serial
      
      
      IF (IHOT.EQ.-67.OR.IHOT.EQ.-68)THEN                               !start reading netCDF hot start file
         timenc_len=1
         IF(.NOT. ALLOCATED(timenc))  ALLOCATE(timenc(timenc_len))
!         timenc(timenc_len)=time
         ncmode=1                                                       !serial version

        IF(.NOT.ALLOCATED(eta11))ALLOCATE(eta11(np,timenc_len))
        IF(.NOT.ALLOCATED(eta22))ALLOCATE(eta22(np,timenc_len))
        IF(.NOT.ALLOCATED(etadd))ALLOCATE(etadd(np,timenc_len))
        IF(.NOT.ALLOCATED(uu22))ALLOCATE(uu22(np,timenc_len))
        IF(.NOT.ALLOCATED(vv22))ALLOCATE(vv22(np,timenc_len))
        IF(.NOT.ALLOCATED(nodecodenc))ALLOCATE(nodecodenc(np))
        IF(.NOT.ALLOCATED(noffnc))ALLOCATE (noffnc(ne))
        IF(.NOT. ALLOCATED(xnc)) ALLOCATE(xnc(np))
        IF(.NOT. ALLOCATED(ync)) ALLOCATE(ync(np))
        
        IF(.NOT. ALLOCATED(nbvvnc)) ALLOCATE(nbvvnc(nbou,nvel))
        IF(.NOT. ALLOCATED(NBDV)) ALLOCATE(NBDV(NOPE,NETA))
        IF(.NOT. ALLOCATED(nvellnc)) ALLOCATE(nvellnc(nbou))
        IF(.NOT. ALLOCATED(nvdllnc)) ALLOCATE(nvdllnc(nope))
        IF(.NOT. ALLOCATED(ibtypenc)) ALLOCATE(ibtypenc(nbou))
            IF(.NOT. ALLOCATED(nopenc)) ALLOCATE(nopenc(netanc))
            IF(.NOT. ALLOCATED(nbounc)) ALLOCATE(nbounc(nvelnc))

!     iestp,nscoue,ivstp,nscouv, icstp, 
!     &                  nscouc, ipstp, iwstp, nscoum, igep, nscouge,
!     &                  igvp, nscougv, igcp, nscougc, igpp, igwp,
!     &                  nscougw,

!     &                    nciestp,nccoue,ncivstp,nccouv, icstp,
!     &                    nscouc, ncipstp, nciwstp, nccoum, 
!     &                    ncigep, nccouge, ncigvp, nccougv, 
!     &                    igcp, nscougc, ncigpp, ncigwp,
!     &                    nccougw,


         call readnc67 (nhsinc,ABS(IHOT),timenc_len,np,timenc,
     &                  nt, nhsinc, uu22,vv22,
     &                  eta11, eta22, EtaDd, nodecodenc,noffnc,
     &                  imhs, iths,
     &                  iestp,nscoue,ivstp,nscouv, icstp, 
     &                  nscouc, ipstp, iwstp, nscoum, igep, nscouge,
     &                  igvp, nscougv, igcp, nscougc, igpp, igwp,
     &                  nscougw,
     &                  rundes, runid, agrid, irtype,
     &                  dtdp, ihot, ics, nolibf, nolifa,
     &                  nolica, nolicat, nwp, 
     &                  ncor, ntip, nws, nramp, tau0,
     &                  statim, reftim, rnday, dramp,
     &                  a00, b00, c00, h0,
     &                  slam0, sfea0, cf ,eslm, cori,
     &                  ntif, nbfr, ne, nhy, nm,
     &                  xnc, ync, dp,
     &                  nope, neta, nbdv,nopenc,netanc,
     &                  nvdllnc,
     &                  nbou, nvel, nbvvnc,nbounc,nvelnc,
     &                  nvellnc,ibtypenc,
     &                  ncmode,
     &                  title, institution, source, history,
     &                  references, comments, host, convention, contact,
     &                  base_date)
       time= timenc(timenc_len)

             DO I=1,MNP
                eta1 (I)=eta11(I, 1)
!                PRINT *, "ETA1(I) = ",I, ETA1(I)
             END DO
             DO I=1,MNP
                eta2 (I)=eta22(I, 1)
!                PRINT *, "ETA2(I) = ",I, ETA2(I)
             END DO
             DO I=1,MNP
                EtaDisc (I)=etadd(I, 1)
!                PRINT *, "EtaDisc(I) = ",I, EtaDisc(I)
             END DO
             DO I=1,MNP
                UU2 (I)=UU22(I, 1)
!                PRINT *, "UU2(I) = ",I, UU2(I)
             END DO
             DO I=1,MNP
                VV2 (I)=VV22(I, 1)
!                PRINT *, "VV22(I) = ",I, VV2(I)
            END DO
             DO I=1,MNP
                nnodecode (I)=nodecodenc(I)
!                PRINT *, "nodecodenc(I) = ",I, nnodecode(I)
            END DO
             DO I=1,MNE
                noff (I)=noffnc(I)
!                PRINT *, "noff(I) = ",I, noff(I)
             END DO
             TIME=timenc(timenc_len)
             NP_G_IN = mnp
             NE_G_IN = mne
             NP_A_IN = mnp
             NE_A_IN = mne
!                   PRINT *,"  ----------START--NETCDF after read ------"
!                   PRINT *,"  IMHS = ", IMHS
!                   PRINT *,"  TIME = ", TIME
!                   PRINT *,"  TIME = ", timenc
!                   PRINT *,"  ITHS = ", ITHS
!                   PRINT *,"  ETA11 = ", ETA1(11)
!                   PRINT *,"  ETA22 = ", ETA2(11)
!                   PRINT *,"  ETADD = ", EtaDisc(11)
!                   PRINT *,"  uu22 = ", uu2(11)
!                   PRINT *,"  vv22 = ", vv2(11)
!                   PRINT *,"  nodecode = ", nodecode(11)
!                   PRINT *,"  noff = ", noff(11)
!                   PRINT *,"  IESTP = ",IESTP
!                   PRINT *,"  NSCOUE = ",NSCOUE
!                   PRINT *,"  IVSTP = ",IVSTP
!                   PRINT *,"  NSCOUV = ",NSCOUV
!                   PRINT *,"  ICSTP = ",ICSTP
!                   PRINT *,"  NSCOUC = ",NSCOUC
!                   PRINT *,"  IPSTP = ",IPSTP
!                   PRINT *,"  IWSTP = ",IWSTP
!                   PRINT *,"  NSCOUM = ",NSCOUM
!                   PRINT *,"  IGEP = ",IGEP
!                   PRINT *,"  NSCOUGE = ",NSCOUGE
!                   PRINT *,"  IGVP = ",IGVP
!                   PRINT *,"  NSCOUGV = ",NSCOUGV
!                   PRINT *,"  IGCP = ",IGCP
!                   PRINT *,"  NSCOUGC = ",NSCOUGC
!                   PRINT *,"  IGPP = ",IGPP
!                   PRINT *,"  IGWP = ",IGWP
!                   PRINT *,"  NSCOUGW = ",NSCOUGW
!             PRINT *,"  ---------END-------NETCDF--------------------"
      ENDIF                                                             ! finish reading netCDF hot start file
      
      ENDIF     

C   End Addition for netCDF by MCF 5/8/08
#endif


C
C     Now read the hot start file.

       IF(( IHOT.EQ.-67.or. IHOT.EQ.-68 )) IHOTNC=IHOT
       
      IF (IHOT.EQ.67.OR.IHOT.EQ.68. or.
     &    (ncmode.eq.2.and.
     &    ( IHOT.EQ.-67.or. IHOT.EQ.-68 ) ) )THEN        ! start reading binary hot start file

       IHOT=ABS(IHOT)

      IHOTSTP=1
      READ(IHOT,REC=IHOTSTP) InputFileFmtVn ; IHOTSTP = IHOTSTP + 1

      IF(IHOTNC.GE.0)THEN
      if (.not.CMP_VERSION_NUMBERS(InputFileFmtVn,FileFmtVersion))then
!      PRINT *,"IHOTNC ", IHOTNC
!      PRINT *,"InputFileFmtVn",InputFileFmtVn
!      PRINT *,"FileFmtVersion",FileFmtVersion
        WRITE(16, 1001)
        write(16, 1012)
        write(16, 9973)
        if (NScreen /= 0 .and. myproc == 0) then
          WRITE(ScreenUnit, 1001)
          write(ScreenUnit, 1012)
          write(ScreenUnit, 9973)
        endif
        stop
      endif
      ENDIF
      READ(IHOT,REC=IHOTSTP) IMHS         ; IHOTSTP = IHOTSTP + 1
      READ(IHOT,REC=IHOTSTP) TIME         ; IHOTSTP = IHOTSTP + 1
      READ(IHOT,REC=IHOTSTP) ITHS         ; IHOTSTP = IHOTSTP + 1

C   kmd48.33bc - implemented a varying Dt for 3D baroclinic
      IF((IDEN.GE.1).AND.(C3D).AND.(CBAROCLINIC)) THEN 
         CHOTHS=.TRUE.
      END IF 
      IF (CHOTHS) THEN 
         DTDPHS=TIME/ITHS  ! Set up time step information for hotstart file. 
         IF (DTDPHS.NE.DTDP) THEN
            RNDAYHS=TIME/86400.d0
            NTHS=INT((RNDAY-RNDAYHS)*(86400.d0/DTDP)+0.5d0)
            NT=ITHS+NTHS
            Write(16,*) "The number of time steps has been reset
     &             due to a change in the timestep interval"
            WRITE(16,*) "It has been changed to ", NT
         END IF 
      ELSE 
        DTDPHS=DTDP
      END IF 


      READ(IHOT,REC=IHOTSTP) NP_G_IN      ; IHOTSTP = IHOTSTP + 1
      READ(IHOT,REC=IHOTSTP) NE_G_IN      ; IHOTSTP = IHOTSTP + 1
      READ(IHOT,REC=IHOTSTP) NP_A_IN      ; IHOTSTP = IHOTSTP + 1
      READ(IHOT,REC=IHOTSTP) NE_A_IN      ; IHOTSTP = IHOTSTP + 1
      
!      IF(ncmode.eq.2.and.MyProc.eq.0)PRINT *,"IMHS = ",IMHS
!      IF(ncmode.eq.2.and.MyProc.eq.0)PRINT *,"TIME = ",TIME
!      IF(ncmode.eq.2.and.MyProc.eq.0)PRINT *,"ITHS = ",ITHS
!      IF(ncmode.eq.2.and.MyProc.eq.0)PRINT *,"NP_G_IN = ",NP_G_IN
!      IF(ncmode.eq.2.and.MyProc.eq.0)PRINT *,"NE_G_IN = ",NE_G_IN
!      IF(ncmode.eq.2.and.MyProc.eq.0)PRINT *,"NP_A_IN = ",NP_A_IN
!      IF(ncmode.eq.2.and.MyProc.eq.0)PRINT *,"NE_A_IN = ",NE_A_IN
      
      DO I=1, NP
         READ(IHOT,REC=IHOTSTP) ETA1(I)   ; IHOTSTP = IHOTSTP + 1









      END DO
      DO I=1, NP
         READ(IHOT,REC=IHOTSTP) ETA2(I)   ; IHOTSTP = IHOTSTP + 1









      END DO
C     jgf46.34 Added support for IBTYPE=52.
      DO I=1, NP
         READ(IHOT,REC=IHOTSTP) EtaDisc(I); IHOTSTP = IHOTSTP + 1
      END DO
      DO I=1, NP
         READ(IHOT,REC=IHOTSTP) UU2(I)    ; IHOTSTP = IHOTSTP + 1
      END DO
      DO I=1, NP
         READ(IHOT,REC=IHOTSTP) VV2(I)    ; IHOTSTP = IHOTSTP + 1
      END DO
      IF(IMHS.EQ.10) THEN
        DO I=1, NP
           READ(IHOT,REC=IHOTSTP) CH1(I)  ; IHOTSTP = IHOTSTP + 1
        END DO
      ENDIF
      DO I=1, NP
         READ(IHOT,REC=IHOTSTP) NNODECODE(I)  ; IHOTSTP = IHOTSTP + 1
      END DO
      DO I=1,NE
         READ(IHOT,REC=IHOTSTP) NOFF(I)       ; IHOTSTP = IHOTSTP + 1
      END DO

      ENDIF                                                             ! finish reading binary hot start file
!      if(myproc.eq.0) THEN
!                   PRINT *,"  ETA1(1) = ", ETA1(1)
!                   PRINT *,"  ETA2(1) = ", ETA2(1)
!      ENDIF
      IF (DTDP.NE.DTDPHS) THEN
         DO I=1,NP
            ETA1(I)=ETA2(I) + DTDP*((ETA1(I)-ETA2(I))/DTDPHS)
         END DO 
      END IF 

      DO I=1, NP
         ETAS(I)=ETA2(I)-ETA1(I)
         H2=DP(I)+IFNLFA*ETA2(I)
         QX2(I)=UU2(I)*H2
         QY2(I)=VV2(I)*H2
      END DO
     
C     jgf46.08 Fine grained ramp functions
C     jgf46.21 Split flux into internal and external, added support
C     for IBTYPE=52.
!Kendra: Changed these to DTDPHS because if timestep is changed then
!       the value of DTDP will not match with the ITHS.  
!       
      IF (NRamp.eq.0) THEN
         RampExtFlux1=1.0d0
         RampExtFlux2=1.0d0
         RampIntFlux1=1.0d0
         RampIntFlux2=1.0d0
         RampTip2=1.0d0
         RampMete2=1.0d0
         RampWRad2=1.0d0
      ELSE
         RampExtFlux1=TANH((2.D0*(ITHS-1)*DTDPHS/86400.D0)/DRampExtFlux)
         RampExtFlux2=TANH((2.D0*(ITHS)*DTDPHS/86400.D0)/DRampExtFlux)
         RampIntFlux1=TANH((2.D0*(ITHS-1)*DTDPHS/86400.D0)/DRampIntFlux)
         RampIntFlux2=TANH((2.D0*(ITHS)*DTDPHS/86400.D0)/DRampIntFlux)
         RampTip2=TANH((2.D0*(ITHS)*DTDPHS/86400.D0)/DRampTip)
         RampMete2=TANH((2.D0*(ITHS)*DTDPHS/86400.D0)/DRampMete)
         RampWRad2=TANH((2.D0*(ITHS)*DTDPHS/86400.D0)/DRampWRad)
      ENDIF
#ifdef IBM
      FluxSettlingIT=INT(FluxSettlingTime*86400.d0/DTDPHS,KIND(0.0D0))
#else
      FluxSettlingIT=INT(FluxSettlingTime*86400.d0/DTDPHS)
#endif
      IF(ITHS.LT.(FluxSettlingIT+10)) THEN
         RampIntFlux1=0.0d0
         RampIntFlux2=0.0d0
         RampTip2=0.0d0
         RampMete2=0.0d0
         RampWRad2=0.0d0
      ENDIF
     
Cjgf46.02 comm.out for Katrina  IF(NRAMP.EQ.1) THEN
Cjgf46.02 comm.out for Katrina    RAMP1=TANH((2.D0*(ITHS-1)*DTDP/86400.D0)/DRAMP)
Cjgf46.02 comm.out for Katrina    RAMP2=TANH((2.D0*ITHS*DTDP/86400.D0)/DRAMP)
Cjgf46.02 comm.out for Katrina  ENDIF

C     jgf48.4627 If user only requests meteorological output, then 
C     set a flag that turns off the calculation of the GWCE and
C     momentum equations
      IF (NWS.NE.0) THEN
         IF ((NOUTE.EQ.0) .AND. (NOUTV.EQ.0) .AND. (NOUTC.EQ.0) .AND.
     &       (NOUTGE.EQ.0) .AND. (NOUTGV.EQ.0) .AND. (NOUTGC.EQ.0 ).AND.
     &       (I3DSD.EQ.0) .AND. (I3DSV.EQ.0) .AND. (I3DST.EQ.0) .AND.
     &       (I3DGD.EQ.0) .AND. (I3DGV.EQ.0) .AND. (I3DGT.EQ.0) .AND.
     &       (NHASE.EQ.0) .AND. (NHASV.EQ.0) .AND.
     &       (NHAGE.EQ.0) .AND. (NHAGV.EQ.0) )
     &       THEN
            nonMetOutputOff = .TRUE.
         ENDIF
         IF ( (NOUTM.NE.0) .OR. (NOUTGW.NE.0) ) THEN
            metOutputOn = .TRUE.
         ENDIF
         IF ( metOutputOn .AND. nonMetOutputOff ) THEN
            IF ((NSCREEN.NE.0).AND.MYPROC.EQ.0) THEN
               WRITE(ScreenUnit,*)
     &      'INFO: Only meterological output was requested. ADCIRC'
               WRITE(ScreenUnit,*)
     &      'will not solve the GWCE or momentum equations.'
             ENDIF
             WRITE(16,*)
     &      'INFO: Only meterological output was requested. ADCIRC'
             WRITE(16,*)
     &      'will not solve the GWCE or momentum equations.'
             METONLY = .TRUE. ! set flag
         ENDIF
      ENDIF
C     
C.... SET POSITIONS IN BOUNDARY CONDITION, WIND AND OUTPUT FILES
C     
      WRITE(16,1112)
      WRITE(16,1794)
 1794 FORMAT(//,' INFORMATION ABOUT RE-STARTING THE TIME SERIES',
     &     ' OUTPUT FILES (UNITS 61-64,71-74,81,83),',
     &     /,' WIND/PRESSURE FILE (UNIT 22) AND FLOW BOUNDARY',
     &     ' CONDITION FILE (UNIT 20)',//)

C......INITIALLY, ZERO OUT THE NORMAL FLOW ON ALL BOUNDARIES

      DO I=1,NVEL
         QN2(I)=0.D0
         QN1(I)=0.D0
         QN0(I)=0.D0
      END DO

C.... FIND THE PROPER PLACE IN THE APERIODIC ELEVATION SPECIFIED
C.... BOUNDARY CONDITION FILE IF IT IS REQURIED.

      IF((NOPE.GT.0).AND.(NBFR.EQ.0)) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1112)
         WRITE(16,1112)
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1976)
         WRITE(16,1976)
 1976    FORMAT(/,1X,'LOCATING ELEVATION SPECIFIED INFORMATION IN ',
     &        'UNIT 19',/)
         OPEN(19,FILE=TRIM(INPUTDIR)//'/'//'fort.19')
         READ(19,*) ETIMINC
         ETIME1=STATIM*86400.D0
         ETIME2=ETIME1+ETIMINC
         DO J=1,NETA
            READ(19,*) ESBIN1(J)
         END DO
         DO J=1,NETA
            READ(19,*) ESBIN2(J)
         END DO
         DO IT=1,ITHS-1
            TIMEIT=IT*DTDPHS + STATIM*86400.D0 ! kmd48.33bc - changed DTDP to DTDPHS
            IF(TIMEIT.GT.ETIME2) THEN
               ETIME1=ETIME2
               ETIME2=ETIME1+ETIMINC
               DO J=1,NETA
                  ESBIN1(J)=ESBIN2(J)
                  READ(19,*) ESBIN2(J)
               END DO
            ENDIF
         END DO
         IF(TIME.GT.ETIME2) THEN
            ETIME1=ETIME2
            ETIME2=ETIME1+ETIMINC
            DO J=1,NETA
               ESBIN1(J)=ESBIN2(J)
               READ(19,*) ESBIN2(J)
            END DO
         ENDIF
         ETRATIO=(TIMEIT-ETIME1)/ETIMINC
      ENDIF

C......FIND PROPER PLACE IN THE APERIODIC NORMAL FLOW BOUNDARY CONDITION
C......FILE IF IT IS REQUIRED

      IF((NFLUXF.EQ.1).AND.(NFFR.EQ.0)) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1112)
         WRITE(16,1112)
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1978)
         WRITE(16,1978)
 1978    FORMAT(/,1X,'LOCATING NORMAL FLOW INFORMATION IN UNIT 20',/)
         OPEN(20,FILE=TRIM(INPUTDIR)//'/'//'fort.20')
         READ(20,*) FTIMINC
         QTIME1=STATIM*86400.D0
         QTIME2=QTIME1+FTIMINC
         DO J=1,NVEL
            QNIN1(J)=0.D0
            IF((LBCODEI(J).EQ.2).OR.(LBCODEI(J).EQ.12)
     &           .OR.(LBCODEI(J).EQ.22))
     &           READ(20,*) QNIN1(J)
         END DO
         DO J=1,NVEL
            QNIN2(J)=0.D0
            IF((LBCODEI(J).EQ.2).OR.(LBCODEI(J).EQ.12)
     &           .OR.(LBCODEI(J).EQ.22))
     &           READ(20,*) QNIN2(J)
         END DO
         DO IT=1,ITHS-1
            TIMEIT=IT*DTDPHS + STATIM*86400.D0 ! kmd48.33bc - changed DTDP to DTDPHS
            IF(TIMEIT.GT.QTIME2) THEN
               QTIME1=QTIME2
               QTIME2=QTIME2+FTIMINC
               DO J=1,NVEL
                  IF((LBCODEI(J).EQ.2).OR.(LBCODEI(J).EQ.12)
     &                 .OR.(LBCODEI(J).EQ.22)) THEN
                     QNIN1(J)=QNIN2(J)
                     READ(20,*) QNIN2(J)
                  ENDIF
               END DO
            ENDIF
         END DO
         QTRATIO=(TIMEIT-QTIME1)/FTIMINC
         DO I=1,NVEL
            QN1(I)=RampExtFlux1*(QNIN1(I)+QTRATIO*(QNIN2(I)-QNIN1(I)))
         END DO
         IF(TIME.GT.QTIME2) THEN
            QTIME1=QTIME2
            QTIME2=QTIME1+FTIMINC
            DO J=1,NVEL
               IF((LBCODEI(J).EQ.2).OR.(LBCODEI(J).EQ.12)
     &              .OR.(LBCODEI(J).EQ.22)) THEN
                  QNIN1(J)=QNIN2(J)
                  READ(20,*) QNIN2(J)
               ENDIF
            END DO
         ENDIF
         QTRATIO=(TIME-QTIME1)/FTIMINC
         DO I=1,NVEL
            QN2(I)=RampExtFlux2*(QNIN1(I)+QTRATIO*(QNIN2(I)-QNIN1(I)))
         END DO
      ENDIF

C......RESTART THE PERIODIC NORMAL FLOW BOUNDARY CONDITION

      IF((NFLUXF.EQ.1).AND.(NFFR.GT.0)) THEN
         DO J=1,NFFR
            IF(FPER(J).EQ.0.) THEN
               NCYC=0.
            ELSE
#ifdef IBM
               NCYC=INT(TIME/FPER(J),KIND(0.0d0))
#else
               NCYC=INT(TIME/FPER(J))
#endif
            ENDIF
            ARGJ1=FAMIG(J)*(TIME-DTDPHS-NCYC*FPER(J))+FFACE(J) ! kmd48.33bc - changed from DTDP to DTDPHS
            ARGJ2=FAMIG(J)*(TIME-NCYC*FPER(J))+FFACE(J)
            RFF1=FFF(J)*RampExtFlux1
            RFF2=FFF(J)*RampExtFlux2
            DO I=1,NVEL
               ARG1=ARGJ1-QNPH(J,I)
               ARG2=ARGJ2-QNPH(J,I)
               QN1(I)=QN1(I)+QNAM(J,I)*RFF1*COS(ARG1)
               QN2(I)=QN2(I)+QNAM(J,I)*RFF2*COS(ARG2)
            END DO
         END DO
      ENDIF

C...  
C...  RESTART SUPERCRITICAL OUTWARD NORMAL FLOW OVER SPECIFIED
C.... EXTERNAL BARRIER BOUNDARY NODES
C...  
      IF(NFLUXB.EQ.1) THEN
         DO I=1,NVEL
            IF((LBCODEI(I).EQ.3).OR.(LBCODEI(I).EQ.13)
     &           .OR.(LBCODEI(I).EQ.23)) THEN
               NNBB=NBV(I)
               RBARWL=2.D0*(ETA1(NNBB)-BARLANHT(I))/3.D0
               IF(RBARWL.GT.0.0D0) THEN
                  QN1(I)=-RampIntFlux1
     &                 *BARLANCFSP(I)*RBARWL*(RBARWL*G)**0.5D0
               ENDIF
               RBARWL=2.D0*(ETA2(NNBB)-BARLANHT(I))/3.D0
               IF(RBARWL.GT.0.0D0) THEN
                  QN2(I)=-RampIntFlux2
     &                 *BARLANCFSP(I)*RBARWL*(RBARWL*G)**0.5D0
               ENDIF
            ENDIF
         END DO
      ENDIF

C...  
C...  RESTART INWARD/OUTWARD NORMAL FLOW OVER SPECIFIED
cjj   wm001 - modified/added the following 3 lines
C.... INTERNAL BARRIERS AND FOR INTERNAL BARRIER BOUNDARIES
C.... WITH CROSS BARRIER PIPES
C.... THIS SECTION ONLY RESTARTS THE OVER BARRIER FLOW COMPONENT
C...  
cjj   wm001 - modified following line
      IF(NFLUXIB.EQ.1) THEN
         DO I=1,NVEL
cjj   wm001 - modified following 2 lines
            IF((LBCODEI(I).EQ.4).OR.(LBCODEI(I).EQ.24) 
     &           .OR.(LBCODEI(I).EQ.5).OR.(LBCODEI(I).EQ.25)) THEN
               NNBB1=NBV(I)     ! GLOBAL NODE NUMBER ON THIS SIDE OF BARRIER
               NNBB2=IBCONN(I)  ! GLOBAL NODE NUMBER ON OPPOSITE SIDE OF BARRIER
C.............RESET INFORMATION FOR K-1 TIME LEVEL
               RBARWL1=ETA1(NNBB1)-BARINHT(I)
               RBARWL2=ETA1(NNBB2)-BARINHT(I)
               RBARWL1F=2.0D0*RBARWL1/3.0D0
               RBARWL2F=2.0D0*RBARWL2/3.0D0
               IF((RBARWL1.LT.0.0).AND.(RBARWL2.LT.0.0)) THEN ! WATER LEVEL BELOW BARRIER
                  QN1(I)=0.0D0  ! NO FLOW
                  GOTO 1998
               ENDIF
               IF(RBARWL1.EQ.RBARWL2) THEN ! WATER LEVEL EQUAL ON BOTH SIDES OF BARRIER
                  QN1(I)=0.0D0  ! NO FLOW
                  GOTO 1998
               ENDIF
               IF(RBARWL1.GT.RBARWL2) THEN ! WATER LEVEL GREATER ON THIS SIDE OF BARRIER
                  IF(RBARWL2.GT.RBARWL1F) THEN ! OUTWARD SUBCRITICAL FLOW
                     QN1(I)=-RampIntFlux1*BARINCFSB(I)*RBARWL2*
     &                    (2.d0*G*(RBARWL1-RBARWL2))**0.5D0
                     GOTO 1998
                  ELSE          ! OUTWARD SUPERCRITICAL FLOW
                     QN1(I)=-RampIntFlux1
     &                    *BARINCFSP(I)*RBARWL1F*(RBARWL1F*G)**0.5D0
                     GOTO 1998
                  ENDIF
               ENDIF
               IF(RBARWL2.GT.RBARWL1) THEN ! WATER LEVEL LOWER ON THIS SIDE OF BARRIER
                  IF(RBARWL1.GT.RBARWL2F) THEN ! INWARD SUBCRITICAL FLOW
                     QN1(I)=RampIntFlux1*BARINCFSB(I)*RBARWL1*
     &                    (2.d0*G*(RBARWL2-RBARWL1))**0.5D0
                     GOTO 1998
                  ELSE          ! INWARD SUPERCRITICAL FLOW
                     QN1(I)=RampIntFlux1
     &                    *BARINCFSP(I)*RBARWL2F*(RBARWL2F*G)**0.5D0
                     GOTO 1998
                  ENDIF
               ENDIF
 1998          CONTINUE
C.............RESET INFORMATION FOR K TIME LEVEL
               RBARWL1=ETA2(NNBB1)-BARINHT(I)
               RBARWL2=ETA2(NNBB2)-BARINHT(I)
               RBARWL1F=2.0D0*RBARWL1/3.0D0
               RBARWL2F=2.0D0*RBARWL2/3.0D0
               IF((RBARWL1.LT.0.0).AND.(RBARWL2.LT.0.0)) THEN ! WATER LEVEL BELOW BARRIER
                  QN2(I)=0.0D0  ! NO FLOW
                  GOTO 1999
               ENDIF
               IF(RBARWL1.EQ.RBARWL2) THEN ! WATER LEVEL EQUAL ON BOTH SIDES OF BARRIER
                  QN2(I)=0.0D0  ! NO FLOW
                  GOTO 1999
               ENDIF
               IF(RBARWL1.GT.RBARWL2) THEN ! WATER LEVEL GREATER ON THIS SIDE OF BARRIER
                  IF(RBARWL2.GT.RBARWL1F) THEN ! OUTWARD SUBCRITICAL FLOW
                     QN2(I)=-RampIntFlux2*BARINCFSB(I)*RBARWL2*
     &                    (2.d0*G*(RBARWL1-RBARWL2))**0.5D0
                     GOTO 1999
                  ELSE          ! OUTWARD SUPERCRITICAL FLOW
                     QN2(I)=-RampIntFlux2
     &                    *BARINCFSP(I)*RBARWL1F*(RBARWL1F*G)**0.5D0
                     GOTO 1999
                  ENDIF
               ENDIF
               IF(RBARWL2.GT.RBARWL1) THEN !WATER LEVEL LOWER ON THIS SIDE OF BARRIER
                  IF(RBARWL1.GT.RBARWL2F) THEN ! INWARD SUBCRITICAL FLOW
                     QN2(I)=RampIntFlux2*BARINCFSB(I)*RBARWL1*
     &                    (2.d0*G*(RBARWL2-RBARWL1))**0.5D0
                     GOTO 1999
                  ELSE          ! INWARD SUPERCRITICAL FLOW
                     QN2(I)=RampIntFlux2
     &                    *BARINCFSP(I)*RBARWL2F*(RBARWL2F*G)**0.5D0
                     GOTO 1999
                  ENDIF
               ENDIF
 1999          CONTINUE
            ENDIF
         END DO
      ENDIF

cjj   wm001 - start add
C...  
C...  RESTART INWARD/OUTWARD NORMAL FLOW OVER SPECIFIED
C.... INTERNAL BARRIERS WITH CROSS BARRIER PIPES
C.... THIS SECTION RESTARTS THE PIPE FLOW COMPONENT
C.... NOTE THAT PIPE FLOW COMPONENT IS ADDED INTO BARRIER FLOW COMPONENT
C.... THAT WAS PREVIOUSLY SET
C...  
      IF(NFLUXIBP.EQ.1) THEN
         DO I=1,NVEL
            IF((LBCODEI(I).EQ.5).OR.(LBCODEI(I).EQ.25)) THEN
               NNBB1=NBV(I)     ! GLOBAL NODE NUMBER ON THIS SIDE OF BARRIER
               NNBB2=IBCONN(I)  ! GLOBAL NODE NUMBER ON OPPOSITE SIDE OF BARRIER
C.............RESET INFORMATION FOR K-1 TIME LEVEL
               RBARWL1=ETA1(NNBB1)-PIPEHT(I)
               RBARWL2=ETA1(NNBB2)-PIPEHT(I)
               IF((RBARWL1.LT.0.0).AND.(RBARWL2.LT.0.0)) THEN ! WATER LEVEL BELOW PIPE
                  QN1(I)=QN1(I)+0.0D0 ! NO FLOW
                  GOTO 2002
               ENDIF
               IF(RBARWL1.EQ.RBARWL2) THEN ! WATER LEVEL EQUAL ON BOTH SIDES OF PIPE
                  QN1(I)=QN1(I)+0.0D0 ! NO FLOW
                  GOTO 2002
               ENDIF
               IF(RBARWL1.GT.RBARWL2) THEN ! WATER LEVEL GREATER ON THIS SIDE OF PIPE
                  IF(RBARWL2.LE.0) THEN ! OUTWARD FREE DISCHARGE 
                     QN1(I)=QN1(I)-RampIntFlux1
     &                    *0.25D0*PI*(PIPEDIAM(I))**2
     &                    *(2.D0*G*RBARWL1/(1+PIPECOEF(I)))**0.5D0
                     GOTO 2002
                  ELSE          ! OUTWARD SUBMERGED DISCHARGE
                     QN1(I)=QN1(I)-RampIntFlux1
     &                    *0.25D0*PI*(PIPEDIAM(I))**2
     &                    *(2.D0*G*(RBARWL1-RBARWL2)/PIPECOEF(I))**0.5D0
                     GOTO 2002
                  ENDIF
               ENDIF
               IF(RBARWL2.GT.RBARWL1) THEN ! WATER LEVEL LOWER ON THIS SIDE OF PIPE
                  IF(RBARWL1.LE.0) THEN ! INWARD FREE DISCHARGE 
                     QN1(I)=QN1(I)+RampIntFlux1
     &                    *0.25D0*PI*(PIPEDIAM(I))**2
     &                    *(2.D0*G*RBARWL2/(1+PIPECOEF(I)))**0.5D0
                     GOTO 2002
                  ELSE          ! INWARD SUBMERGED DISCHARGE
                     QN1(I)=QN1(I)+RampIntFlux1
     &                    *0.25D0*PI*(PIPEDIAM(I))**2
     &                    *(2.D0*G*(RBARWL2-RBARWL1)/PIPECOEF(I))**0.5D0
                     GOTO 2002
                  ENDIF
               ENDIF
 2002          CONTINUE
C.............RESET INFORMATION FOR K TIME LEVEL
               RBARWL1=ETA2(NNBB1)-PIPEHT(I)
               RBARWL2=ETA2(NNBB2)-PIPEHT(I)
               IF((RBARWL1.LT.0.0).AND.(RBARWL2.LT.0.0)) THEN ! WATER LEVEL BELOW PIPE
                  QN2(I)=QN2(I)+0.0D0 ! NO FLOW
                  GOTO 2003
               ENDIF
               IF(RBARWL1.EQ.RBARWL2) THEN ! WATER LEVEL EQUAL ON BOTH SIDES OF PIPE
                  QN2(I)=QN2(I)+0.0D0 ! NO FLOW
                  GOTO 2003
               ENDIF
               IF(RBARWL1.GT.RBARWL2) THEN ! WATER LEVEL GREATER ON THIS SIDE OF PIPE
                  IF(RBARWL2.LE.0) THEN ! OUTWARD FREE DISCHARGE
                     QN2(I)=QN2(I)-RampIntFlux2
     &                    *0.25D0*PI*(PIPEDIAM(I))**2
     &                    *(2.D0*G*RBARWL1/(1+PIPECOEF(I)))**0.5D0
                     GOTO 2003
                  ELSE          ! OUTWARD SUBMERGED DISCHARGE
                     QN2(I)=QN2(I)-RampIntFlux2
     &                    *0.25D0*PI*(PIPEDIAM(I))**2
     &                    *(2.D0*G*(RBARWL1-RBARWL2)/PIPECOEF(I))**0.5D0
                     GOTO 2003
                  ENDIF
               ENDIF
               IF(RBARWL2.GT.RBARWL1) THEN !WATER LEVEL LOWER ON THIS SIDE OF PIPE
                  IF(RBARWL1.LE.0) THEN ! INWARD FREE DISCHARGE
                     QN2(I)=QN2(I)+RampIntFlux2
     &                    *0.25D0*PI*(PIPEDIAM(I))**2
     &                    *(2.D0*G*RBARWL2/(1+PIPECOEF(I)))**0.5D0
                     GOTO 2003
                  ELSE          ! INWARD SUBMERGED DISCHARGE
                     QN2(I)=QN2(I)+RampIntFlux2
     &                    *0.25D0*PI*(PIPEDIAM(I))**2
     &                    *(2.D0*G*(RBARWL2-RBARWL1)/PIPECOEF(I))**0.5D0
                     GOTO 2003
                  ENDIF
               ENDIF
 2003          CONTINUE
            ENDIF
         END DO
      ENDIF
cjjwm001 - end add    

C------------------------MET FORCING---------------------------------------
C
C......RESTART WIND AND PRESSURE INFORMATION
C
C.... tcm v49.16 20100617
C.... convert background pressure from millibars to meters of water
      PRBCKGRND_MH2O = 100.0D0*PRBCKGRND/(RHOWAT0*G)

Casey 091217: Lower the wind drag limit from 0.0035 to 0.0020.
      WindDragLimit = 0.0035
C     No wind, radiation stress or atmospheric pressure forcings are used.
      IF(NWS.EQ.1) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         DO J=1,ITHS
            DO I=1,NP
               READ(22,*) NHG,WSX2(I),WSY2(I),PR2(I)
            END DO
         END DO
         DO I=1,NP
            WSX2(I)=RampMete2*WSX2(I)
            WSY2(I)=RampMete2*WSY2(I)
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*PR2(I)
            PRDIFF = RampMete2*(PR2(I)-PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
         END DO
      ENDIF
      
C     Wind stress and atmospheric pressure are read in at all grid nodes
C     at a time interval that does not equal the model time
C     step. Interpolation in time is used to synchronize the wind and
C     pressure information with the model time step.
      IF(NWS.EQ.2) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         WTIME1 = STATIM*86400.D0
         WTIME2 = WTIME1 + WTIMINC
         READ(22,*) (NHG,WVNX1(I),WVNY1(I),PRN1(I),I=1,NP)  !TCM v49.02 Added read for rec. 1
         READ(22,*) (NHG,WVNX2(I),WVNY2(I),PRN2(I),I=1,NP)
         DO IT=1,ITHS
            TIMEIT=IT*DTDPHS + STATIM*86400.D0 ! kmd48.33bc - changed the DTDP to DTDPHS
            IF(TIMEIT.GT.WTIME2) THEN
               WTIME1=WTIME2
               WTIME2=WTIME2+WTIMINC
               DO I=1,NP
                  WVNX1(I)=WVNX2(I)
                  WVNY1(I)=WVNY2(I)
                  PRN1(I)=PRN2(I)
                  READ(22,*) NHG,WVNX2(I),WVNY2(I),PRN2(I)
               END DO
            ENDIF
         END DO
         WTRATIO=(TIME-WTIME1)/WTIMINC
         DO I=1,NP
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WSX2(I) = RampMete2*WINDX
            WSY2(I) = RampMete2*WINDY
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*(PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
            PRDIFF = RampMete2*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                          -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
         END DO
      ENDIF

      IF(NWS.EQ.-2) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         WTIME1 = TIME
         WTIME2 = WTIME1 + WTIMINC
         READ(22,*) (NHG,WVNX1(I),WVNY1(I),PRN1(I),I=1,NP)
         READ(22,*) (NHG,WVNX2(I),WVNY2(I),PRN2(I),I=1,NP)
         DO I=1,NP
            WSX2(I) = RampMete2*WVNX1(I)
            WSY2(I) = RampMete2*WVNY1(I)
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*PRN1(I)
            PRDIFF = RampMete2*(PRN1(I)-PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
         END DO
      ENDIF

C     Wind velocity in US Navy Fleet Numeric format interpolated in
C     space onto the ADCIRC grid and in time to synchronize the wind and
C     pressure information with the model time step. Garratt's formula
C     is used to compute wind stress from the wind velocity.
      IF(NWS.EQ.3) THEN
C....    tcm_v49.04 20091124 -- Changed from local inputdir to global inputdir
C....                           to correspond with having only a global wind file.
!         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         OPEN(22,FILE=TRIM(GBLINPUTDIR)//'/'//'fort.22')
 2223    CALL NWS3GET(X,Y,SLAM,SFEA,WVNX2,WVNY2,IWTIME,IWYR,WTIMED,NP,
     &        NWLON,NWLAT,WLATMAX,WLONMIN,WLATINC,WLONINC,ICS,
     &        NScreen,ScreenUnit)
         IF(IWYR.NE.IREFYR) THEN
            IWTIMEP=IWTIME
            DO I=1,NP
               WVNX1(I)=WVNX2(I)
               WVNY1(I)=WVNY2(I)
            END DO
            GOTO 2223
         ENDIF
         IF(WTIMED.LE.WREFTIM) THEN
            IWTIMEP=IWTIME
            DO I=1,NP
               WVNX1(I)=WVNX2(I)
               WVNY1(I)=WVNY2(I)
            END DO
            GOTO 2223
         ENDIF
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) 
     &        WRITE(ScreenUnit,*)'FOUND WIND DATA AT TIME= ',IWTIMEP
         WRITE(16,*) 'FOUND WIND DATA AT TIME =',IWTIMEP
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) 
     &        WRITE(ScreenUnit,*)'FOUND WIND DATA AT TIME= ',IWTIME
         WRITE(16,*) 'FOUND WIND DATA AT TIME =',IWTIME
         WTIME2=WTIMED-WREFTIM  !CAST INTO MODEL TIME REFERENCE
         WTIME1=WTIME2-WTIMINC
         DO IT=1,ITHS
            TIMEIT=IT*DTDPHS + STATIM*86400.D0 ! kmd48.33bc - changed the DTDP to DTDPHS
            IF(TIMEIT.GT.WTIME2) THEN
               WTIME1=WTIME2
               WTIME2=WTIME2+WTIMINC
               DO I=1,NP
                  WVNX1(I)=WVNX2(I)
                  WVNY1(I)=WVNY2(I)
               END DO
               CALL NWS3GET(X,Y,SLAM,SFEA,WVNX2,WVNY2,IWTIME,
     &              IWYR,WTIMED,NP,NWLON,NWLAT,WLATMAX,WLONMIN,
     &              WLATINC,WLONINC,ICS,NScreen,ScreenUnit)
               IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) 
     &              WRITE(ScreenUnit,*) 
     &              'WIND FILE ADVANCED TO TIME',' = ', IWTIME
               WRITE(16,*) 'WIND FILE ADVANCED TO TIME = ',IWTIME
            ENDIF
         END DO
         WTRATIO=(TIME-WTIME1)/WTIMINC
         DO I=1,NP              !INTERPOLATE IN TIME
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WINDMAG=SQRT(WINDX*WINDX+WINDY*WINDY)
            WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
C     jgf46.01 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY) 
               WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            ENDIF
            WSX2(I)=RampMete2*0.001293d0*WDRAGCO*WINDX*WINDMAG
            WSY2(I)=RampMete2*0.001293d0*WDRAGCO*WINDY*WINDMAG
         END DO
      ENDIF

C     Wind velocity and atmospheric pressure are read in (PBL/JAG
C     format) at selected ADCIRC grid nodes. Interpolation in time is
C     used to synchronize the wind and pressure information with the
C     model time step. Garratt's formula is used to compute wind stress
C     from wind velocity.
      IF(NWS.EQ.4) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         WTIME1 = STATIM*86400.D0
         WTIME2 = WTIME1 + WTIMINC
         CALL NWS4GET(WVNX1,WVNY1,PRN1,NP,RHOWAT0,G)  !TCM v49.02 Added read for rec. 1
         CALL NWS4GET(WVNX2,WVNY2,PRN2,NP,RHOWAT0,G)
         DO IT=1,ITHS
            TIMEIT=IT*DTDPHS + STATIM*86400.D0 ! kmd48.33bc - changed the DTDP to DTDPHS
            IF(TIMEIT.GT.WTIME2) THEN
               WTIME1=WTIME2
               WTIME2=WTIME2+WTIMINC
               DO I=1,NP
                  WVNX1(I)=WVNX2(I)
                  WVNY1(I)=WVNY2(I)
                  PRN1(I)=PRN2(I)
               END DO
               CALL NWS4GET(WVNX2,WVNY2,PRN2,NP,RHOWAT0,G)
            ENDIF
         END DO
         WTRATIO=(TIME-WTIME1)/WTIMINC
         DO I=1,NP
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
C     jgf46.01 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY) 
               WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            ENDIF
            WSX2(I) = RampMete2*0.001293d0*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RampMete2*0.001293d0*WDRAGCO*WINDY*WINDMAG
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*(PRN1(I)+WTRATIO*(PRN2(I)-PRN1(I)))
            PRDIFF = RampMete2*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                           -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
         END DO
      ENDIF

      IF(NWS.EQ.-4) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         WTIME1 = TIME
         WTIME2 = WTIME1 + WTIMINC
         CALL NWS4GET(WVNX1,WVNY1,PRN1,NP,RHOWAT0,G)
         CALL NWS4GET(WVNX2,WVNY2,PRN2,NP,RHOWAT0,G)
         DO I=1,NP
            WINDX = WVNX1(I)
            WINDY = WVNY1(I)
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
C     jgf46.01 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY) 
               WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            ENDIF
            WSX2(I) = RampMete2*0.001293d0*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RampMete2*0.001293d0*WDRAGCO*WINDY*WINDMAG
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*PRN1(I)
            PRDIFF = RampMete2*(PRN1(I)-PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
         END DO
      ENDIF

C     Wind velocity and atmospheric pressure are read in at all grid
C     nodes. Interpolation in time is used to synchronize the wind and
C     pressure information with the model time step. Garratt's formula
C     is used to compute wind stress from wind velocity.
      IF(NWS.EQ.5) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         WTIME1 = STATIM*86400.D0
         WTIME2 = WTIME1 + WTIMINC
         READ(22,*) (NHG,WVNX1(I),WVNY1(I),PRN1(I),I=1,NP)  !TCM v49.02 Added read for rec. 1
         READ(22,*) (NHG,WVNX2(I),WVNY2(I),PRN2(I),I=1,NP)
         DO IT=1,ITHS
            TIMEIT=IT*DTDPHS + STATIM*86400.D0 ! kmd48.33bc - changed the DTDP to DTDPHS
            IF(TIMEIT.GT.WTIME2) THEN
               WTIME1=WTIME2
               WTIME2=WTIME2+WTIMINC
               DO I=1,NP
                  WVNX1(I)=WVNX2(I)
                  WVNY1(I)=WVNY2(I)
                  PRN1(I)=PRN2(I)
                  READ(22,*) NHG,WVNX2(I),WVNY2(I),PRN2(I)
               END DO
            ENDIF
         END DO
         WTRATIO=(TIME-WTIME1)/WTIMINC
         DO I=1,NP
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
C     jgf46.01 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY) 
               WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            ENDIF
            WSX2(I) = RampMete2*0.001293d0*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RampMete2*0.001293d0*WDRAGCO*WINDY*WINDMAG
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*(PRN1(I)+WTRATIO*(PRN2(I)-PRN1(I)))
            PRDIFF = RampMete2*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                           -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF 
C     jjw-42.06j wrote, jgf46.00 added following two lines:
            WVNXOUT(I)=RampMete2*WINDX
            WVNYOUT(I)=RampMete2*WINDY
         END DO
      ENDIF

      IF(NWS.EQ.-5) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         WTIME1 = TIME
         WTIME2 = WTIME1 + WTIMINC
         READ(22,*) (NHG,WVNX1(I),WVNY1(I),PRN1(I),I=1,NP)
         READ(22,*) (NHG,WVNX2(I),WVNY2(I),PRN2(I),I=1,NP)
         WTRATIO=(TIME-WTIME1)/WTIMINC !jgf46.00 added
         DO I=1,NP
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I)) 
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I)) 
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
C     jgf46.01 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY) 
               WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            ENDIF
            WSX2(I) = RampMete2*0.001293d0*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RampMete2*0.001293d0*WDRAGCO*WINDY*WINDMAG
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*PRN1(I)
            PRDIFF = RampMete2*(PRN1(I)-PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            WVNXOUT(I)=RampMete2*WINDX !jgf46.00 added  ! tcm 20100618 v49.16 changed RampMete to RampMete2
            WVNYOUT(I)=RampMete2*WINDY !jgf46.00 added  ! tcm 20100618 v49.16 changed RampMete to RampMete2
         END DO
      ENDIF

C     Wind velocity and atmospheric pressure are read in for a
C     rectangular grid (either in Longitude, Latitude or Cartesian
C     coordinates, consistent with the grid coordinates) and
C     interpolated in space onto the ADCIRC grid and in time to
C     synchronize the wind and pressure information with the model time
C     step. Garratt's formula is used to compute wind stress from the
C     wind velocity.
      IF(NWS.EQ.6) THEN
C....    tcm_v49.04 20091124 -- Changed from local inputdir to global inputdir
C....                           to correspond with having only a global wind file.
!         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         OPEN(22,FILE=TRIM(GBLINPUTDIR)//'/'//'fort.22')
C     The following 3 lines are a hardwire to allow a non standard met
C     file to be read in at time zero in a hot start.  They should be
C     eliminated or commented out for normal operation
c     OPEN(199,FILE=TRIM(INPUTDIR)//'/'//'fort.199')
c     READ(199,*) (NHG,PRN1(I),WVNX1(I),WVNY1(I),I=1,NP)
c     CLOSE(199)
C     The following CALL statement should be uncommented for normal operation
         CALL NWS6GET(X,Y,SLAM,SFEA,WVNX1,WVNY1,PRN1,NP,NWLON,NWLAT,
     &        WLATMAX,WLONMIN,WLATINC,WLONINC,ICS,RHOWAT0,G)
         CALL NWS6GET(X,Y,SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,NWLON,NWLAT,
     &        WLATMAX,WLONMIN,WLATINC,WLONINC,ICS,RHOWAT0,G)
         WTIME1=TIME
         WTIME2=WTIME1+WTIMINC
         WTRATIO=(TIME-WTIME1)/WTIMINC
         DO I=1,NP
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
C     jgf46.00 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY) 
               WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            ENDIF
            WSX2(I) = RampMete2*0.001293d0*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RampMete2*0.001293d0*WDRAGCO*WINDY*WINDMAG
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*(PRN1(I)+WTRATIO*(PRN2(I)-PRN1(I)))
            PRDIFF = RampMete2*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                           -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
         END DO
      ENDIF

C     jgf46.01 New option to read in surface wind stress and atmospheric
C     pressure for a rectangular grid (either in Longitude, Latitude or
C     Cartesian coordinates, consistent with the grid coordinates) and
C     interpolate in space onto the ADCIRC grid. Interpolation in time
C     is used to synchronize the wind and pressure information with the
C     model time step.
C
      IF(NWS.EQ.7) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         WTIME1 = STATIM*86400.D0
         CALL NWS7GET(X,Y,SLAM,SFEA,WVNX1,WVNY1,PRN1,NP,NWLON,NWLAT,
     &        WLATMAX,WLONMIN,WLATINC,WLONINC,ICS,RHOWAT0,G)  !TCMv49.02 Changed (WVNX2,WVNY2,PRN2) to (WVNX1,WVNY1,PRN1)
         CALL NWS7GET(X,Y,SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,NWLON,NWLAT,
     &        WLATMAX,WLONMIN,WLATINC,WLONINC,ICS,RHOWAT0,G)
         WTIME1=TIME
         WTIME2=WTIME1+WTIMINC
         WTRATIO=(TIME-WTIME1)/WTIMINC
         DO I=1,NP
            WindX = WVNX1(I) + WTRatio*(WVNX2(I)-WVNX1(I))
            WindY = WVNY1(I) + WTRatio*(WVNY2(I)-WVNY1(I))
            WSX2(I) = RampMete2*WindX !apply ramp
            WSY2(I) = RampMete2*WindY
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*(PRN1(I)+WTRATIO*(PRN2(I)-PRN1(I)))
            PRDIFF = RampMete2*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                           -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            wvnxout(i)=WSX2(i) !for met recording sta. output
            wvnyout(i)=WSY2(i)             
         END DO
      ENDIF

      IF(NWS.EQ.-7) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         WTIME1 = STATIM*86400.D0
         CALL NWS7GET(X,Y,SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,NWLON,NWLAT,
     &        WLATMAX,WLONMIN,WLATINC,WLONINC,ICS,RHOWAT0,G)
         CALL NWS7GET(X,Y,SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,NWLON,NWLAT,
     &        WLATMAX,WLONMIN,WLATINC,WLONINC,ICS,RHOWAT0,G)
         WTIME1=TIME
         WTIME2=WTIME1+WTIMINC
         WTRATIO=(TIME-WTIME1)/WTIMINC
C     jgfdebug46.01 How do we convert these marine wind stresses into
C     directional land surface stresses?
         DO I=1,NP
            WindX = WVNX1(I) + WTRatio*(WVNX2(I)-WVNX1(I))
            WindY = WVNY1(I) + WTRatio*(WVNY2(I)-WVNY1(I))
            WSX2(I) = RampMete2*WindX !apply ramp
            WSY2(I) = RampMete2*WindY
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*(PRN1(I)+WTRATIO*(PRN2(I)-PRN1(I)))
            PRDIFF = RampMete2*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                           -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            wvnxout(i)=WSX2(i) !for met recording sta. output
            wvnyout(i)=WSY2(i)             
         END DO
      ENDIF

C     jgf46.02 New option to read in hurricane locations and generate 
C     generate hurricane winds from the Holland Wind Model.
      IF(ABS(NWS).EQ.8) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         HollandTime = TIME
         CALL HollandGet(X,Y,SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,
     &        ICS,RHOWAT0,G,HollandTime,NSCREEN,ScreenUnit)
         DO I=1,NP
            WindX = WVNX2(I)
            WindY = WVNY2(I)
            WindMag = SQRT(WindX*WindX+WindY*WindY)
            WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY) 
               WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            ENDIF
            WSX2(I) = RampMete2*0.001293d0*WDragCo*WindX*WindMag
            WSY2(I) = RampMete2*0.001293d0*WDragCo*WindY*WindMag
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*PRN2(I)
            PRDIFF = RampMete2*(PRN2(I)-PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            WVNXOUT(I)=RampMete2*WindX
            WVNYOUT(I)=RampMete2*WindY
         ENDDO
      ENDIF 

!     cf & cm added nws = 9: asymmetric hurricane winds

      IF(ABS(NWS).EQ.9) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
             CALL NWS9GET(SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,TIME, ICS)
!         IF(NWS.EQ.9) THEN
!             CALL NWS9GET(SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,TIME, ICS)
!         ENDIF
!         IF(NWS.EQ.-9) THEN
!             CALL NWS9GET(SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,TIME, ICS)
!         ENDIF
         DO I=1,NP
            WindX = WVNX2(I)
            WindY = WVNY2(I)
            WindMag = SQRT(WindX*WindX+WindY*WindY)
            WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY)
               WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            ENDIF
            WSX2(I) = RampMete2*0.001293d0*WDragCo*WindX*WindMag
            WSY2(I) = RampMete2*0.001293d0*WDragCo*WindY*WindMag
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*PRN2(I)
            PRDIFF = RampMete2*(PRN2(I)-PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            WVNXOUT(I)=RampMete2*WindX
            WVNYOUT(I)=RampMete2*WindY
         ENDDO
      ENDIF

!     rjw added nws = 19: asymmetric hurricane winds v2.0

      IF(ABS(NWS).EQ.19) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
             CALL NWS19GET(SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,TIME, ICS)
!         IF(NWS.EQ.9) THEN
!             CALL NWS9GET(SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,TIME, ICS)
!         ENDIF
!         IF(NWS.EQ.-9) THEN
!             CALL NWS9GET(SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,TIME, ICS)
!         ENDIF
         DO I=1,NP
            WindX = WVNX2(I)
            WindY = WVNY2(I)
            WindMag = SQRT(WindX*WindX+WindY*WindY)
            WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY)
               WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            ENDIF
            WSX2(I) = RampMete2*0.001293d0*WDragCo*WindX*WindMag
            WSY2(I) = RampMete2*0.001293d0*WDragCo*WindY*WindMag
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*PRN2(I)
            PRDIFF = RampMete2*(PRN2(I)-PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            WVNXOUT(I)=RampMete2*WindX
            WVNYOUT(I)=RampMete2*WindY
         ENDDO
      ENDIF
C
C     Wind velocity (10 m) and atmospheric pressure are read in from a
C     sequence of National Weather Service (NWS) Aviation (AVN) model
C     output files. Each AVN file is assumed to contain data on a
C     Gaussian longitude, latitude grid at a single time. Consecutive
C     files in the sequence are separated by N hours in time. Garratt's
C     formula is used to compute wind stress from the wind velocity.
      IF(NWS.EQ.10) THEN
         WTIME1=TIME
         WTIME2=WTIME1+WTIMINC
         NWSGGWI=-1
         CALL NWS10GET(NWSGGWI,SLAM,SFEA,WVNX1,WVNY1,PRN1,NP,RHOWAT0,G,
     &        NWLON,NWLAT,WTIMINC)
         NWSGGWI=0
         CALL NWS10GET(NWSGGWI,SLAM,SFEA,WVNX1,WVNY1,PRN1,NP,RHOWAT0,G,
     &        NWLON,NWLAT,WTIMINC)
         NWSGGWI=1
         CALL NWS10GET(NWSGGWI,SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,RHOWAT0,G,
     &        NWLON,NWLAT,WTIMINC)
         WTRATIO=(TIME-WTIME1)/WTIMINC
         DO I=1,NP
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
C     jgf46.00 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY) 
               WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            ENDIF
            WSX2(I) = RampMete2*0.001293d0*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RampMete2*0.001293d0*WDRAGCO*WINDY*WINDMAG
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*(PRN1(I)+WTRATIO*(PRN2(I)-PRN1(I)))
            PRDIFF = RampMete2*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                           -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
         END DO
      ENDIF

C     Wind velocity (10 m) and atmospheric pressure are read in from a
C     sequence of stripped down National Weather Service (NWS) ETA 29km
C     model output files. Each ETA file is assumed to contain data on an
C     E grid for a single day (8 data sets, one every 3 hours, beginning
C     @ 03:00 and continuing through 24:00 of the given day). The wind
C     data is converted to an east-west, north-south coordinate system
C     inside ADCIRC. Garratt's formula is used to compute wind stress
C     from the wind velocity.
      IF(NWS.EQ.11) THEN
         WTIME1=TIME
         WTIME2=WTIME1+WTIMINC
         NWSEGWI=0
         IDSETFLG=0
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1197)
         WRITE(16,1197)
         CALL NWS11GET(NWSEGWI,IDSETFLG,SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,
     &        RHOWAT0,G)        !JUST COMPUTE INTERPOLATING FACTORS
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1198)
         WRITE(16,1198)
         NWSEGWI=0
         IDSETFLG=8
         CALL NWS11GET(NWSEGWI,IDSETFLG,SLAM,SFEA,WVNX1,WVNY1,PRN1,NP,
     &        RHOWAT0,G)        !NOW COMPUTE HOTSTART WIND FILED
         NWSEGWI=1
         IDSETFLG=1
         CALL NWS11GET(NWSEGWI,IDSETFLG,SLAM,SFEA,WVNX2,WVNY2,PRN2,NP,
     &        RHOWAT0,G)        !NOW COMPUTE NEXT WIND FIELD
         WTRATIO=(TIME-WTIME1)/WTIMINC
         DO I=1,NP
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
C     jgf46.00 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY) 
               WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
            ENDIF
            WSX2(I) = RampMete2*0.001293d0*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RampMete2*0.001293d0*WDRAGCO*WINDY*WINDMAG
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*(PRN1(I)+WTRATIO*(PRN2(I)-PRN1(I)))
            PRDIFF = RampMete2*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                           -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
         END DO
      ENDIF

C.....sb46.28sb01 NWS=12,-12 were added to deal with raw OWI files.  09/xx/2006
      IF(NWS.EQ.12) THEN
          WTIME1 = STATIM*86400.D0
          WTIME2 = WTIME1 + WTIMINC
          CALL NWS12INIT(WVNX1,WVNY1,PRN1,NP,RHOWAT0,G)  ! This just initializes some variables and might set WVNX1,etc... if skipping ahead.
          CALL NWS12GET(WVNX1,WVNY1,PRN1,NP,RHOWAT0,G)  !TCM v49.02 Added read for rec. 1
          CALL NWS12GET(WVNX2,WVNY2,PRN2,NP,RHOWAT0,G)
          DO IT=1,ITHS
             TIMEIT=IT*DTDPHS + STATIM*86400.D0 ! kmd48.33bc - changed the DTDP to DTDPHS
             IF(TIMEIT.GT.WTIME2) THEN
                WTIME1=WTIME2
                WTIME2=WTIME2+WTIMINC
                DO I=1,NP
                   WVNX1(I)=WVNX2(I)
                   WVNY1(I)=WVNY2(I)
                   PRN1(I)=PRN2(I)
                 END DO
                 CALL NWS12GET(WVNX2,WVNY2,PRN2,NP,RHOWAT0,G)
              ENDIF
          END DO
          WTRATIO=(TIME-WTIME1)/WTIMINC
          DO I=1,NP
             WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
             WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
             WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
             WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
             IF (LoadDirEffRLen) THEN
                CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
                WindMag = SQRT(WindX*WindX+WindY*WindY) 
                WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
             ENDIF
             WSX2(I) = RampMete2*0.001293d0*WDRAGCO*WINDX*WINDMAG
             WSY2(I) = RampMete2*0.001293d0*WDRAGCO*WINDY*WINDMAG
C  tcm v49.16 20100617
!            PR2(I) = RampMete2*PRN1(I)
             PRDIFF = RampMete2*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                           -PRBCKGRND_MH2O)
             PR2(I) = PRBCKGRND_MH2O + PRDIFF
             WVNXOUT(I)=RampMete2*WINDX
             WVNYOUT(I)=RampMete2*WINDY
          END DO
       ENDIF

       IF(NWS.EQ.-12) THEN
          WTIME1 = TIME
          WTIME2 = WTIME1 + WTIMINC
          CALL NWS12INIT(WVNX1,WVNY1,PRN1,NP,RHOWAT0,G)
          CALL NWS12GET(WVNX1,WVNY1,PRN1,NP,RHOWAT0,G)
          CALL NWS12GET(WVNX2,WVNY2,PRN2,NP,RHOWAT0,G)
          DO I=1,NP
             WINDX = WVNX1(I)
             WINDY = WVNY1(I)
             WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
             WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
             IF (LoadDirEffRLen) THEN
                CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), H0, G, WindX, WindY)
                WindMag = SQRT(WindX*WindX+WindY*WindY) 
                WDragCo = WindDrag(WindMag, WindDragLimit, "Garratt   ")
             ENDIF
             WSX2(I) = RampMete2*0.001293d0*WDRAGCO*WINDX*WINDMAG
             WSY2(I) = RampMete2*0.001293d0*WDRAGCO*WINDY*WINDMAG
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*PRN1(I)
             PRDIFF = RampMete2*(PRN1(I)-PRBCKGRND_MH2O)
             PR2(I) = PRBCKGRND_MH2O + PRDIFF
             WVNXOUT(I)=RampMete2*WINDX  ! tcm 20100618 v49.16 changed RampMete to RampMete2
             WVNYOUT(I)=RampMete2*WINDY  ! tcm 20100618 v49.16 changed RampMete to RampMete2
          END DO
       ENDIF

C
C--------------------END MET FORCING---------------------------------------

C......RESTART THE WAVE RADIATION STRESS

        IF(NRS.GE.1) THEN ! sb46.28sb03
          RSTIME1 = TIME
          RSTIME2 = RSTIME1 + RSTIMINC
          ! NRS=2 was added. sb46.28sb03 09/xx/2006
          IF(NRS.EQ.1) THEN
            OPEN(23,FILE=TRIM(INPUTDIR)//'/'//'fort.23')
            CALL RSGET(RSNX1,RSNY1,NP)
            CALL RSGET(RSNX2,RSNY2,NP)
          ENDIF
          IF(NRS.EQ.2) THEN
            CALL RS2INIT(RSNX1,RSNY1,NP)
            CALL RS2GET(RSNX1,RSNY1,NP)
            CALL RS2GET(RSNX2,RSNY2,NP)
          ENDIF
          DO I=1,NP
            WSX2(I) = WSX2(I)+RampWRad2*RSNX1(I)
            WSY2(I) = WSY2(I)+RampWRad2*RSNY1(I)
            END DO
          ENDIF

       if (CTIP) then
Cjromo 11-01-00  Initialize TIP2 for HOTSTART
          DO I=1,NP
             TIP2(I)=0.0
          END DO
CTIP  LINES TO USE TIDAL POTENTIAL FORCING
       IF(NTIP.GE.1) THEN
         DO J=1,NTIF
           IF(PERT(J).EQ.0.) THEN
             NCYC=0
             ELSE
#ifdef IBM                
             NCYC=INT(TIME/PERT(J),KIND(0.0d0))
#else
             NCYC=INT(TIME/PERT(J))
#endif
             ENDIF
           ARGT=AMIGT(J)*(TIME-NCYC*PERT(J))+FACET(J)
           TPMUL=RampTip2*ETRF(J)*TPK(J)*FFT(J)
           SALTMUL=RampTip2*FFT(J)
#ifdef IBM
           NA=NINT(0.00014/AMIGT(J),KIND(0.0d0))
#else
           NA=NINT(0.00014/AMIGT(J))
#endif
           IF(NA.EQ.1) THEN                        !SEMI-DIURNAL SPECIES
             DO I=1,NP
               ARGTP=ARGT+2.*SLAM(I)
               ARGSALT=ARGT-SALTPHA(J,I)
               CCSFEA=COS(SFEA(I))
               CCSFEA=CCSFEA*CCSFEA
               TIP2(I)=TIP2(I)+TPMUL*CCSFEA*COS(ARGTP)
     &                 +SALTMUL*SALTAMP(J,I)*COS(ARGSALT)
               END DO
             ENDIF
           IF(NA.EQ.2) THEN
             DO I=1,NP
               ARGTP=ARGT+SLAM(I)
               ARGSALT=ARGT-SALTPHA(J,I)
cjjw/vjpm002 - modified/added the following 5 lines
#ifdef REAL8
               S2SFEA=SIN(2.d0*SFEA(I))
#else
               S2SFEA=SIN(2.e0*SFEA(I))
#endif
               TIP2(I)=TIP2(I)+TPMUL*S2SFEA*COS(ARGTP)
     &                +SALTMUL*SALTAMP(J,I)*COS(ARGSALT)
               END DO
             ENDIF
           END DO
         ENDIF
      endif     !   CTIP


C...
C....SET UP TO RESTART TIMESERIES OUTPUT FILES
C....
C...
      IF(NBYTE.EQ.4) ITEMPSTP=20
      IF(NBYTE.EQ.8) ITEMPSTP=10

C...
C....IF RESTARTING THE ELEVATION STATION OUTPUT FILE, GO TO THE PROPER PLACE
C....IN THE FILE.  OTHERWISE ZERO OUT NSCOUE.
C...

      IF (IHOT.EQ.67.OR.IHOT.EQ.68. OR.
     &    (ncmode.eq.2. and. (IHOT.EQ.-67.OR.IHOT.EQ.-68) ) )THEN

      READ(IHOT,REC=IHOTSTP) IESTP ; IHOTSTP = IHOTSTP + 1
      READ(IHOT,REC=IHOTSTP) NSCOUE ; IHOTSTP = IHOTSTP + 1
      WRITE(16,1040) IESTP,NSCOUE
 1040 FORMAT(//,1X,I6,' LINES OR RECORDS WRITTEN IN ELEVATION ',
     &                'STATION FILE BY THE TIME OF THE HOT START',
     &           /,8X,'SPOOL COUNTER = ',I6)
     
      ENDIF 
      
      IF (ncmode.ne.0. and. (IHOTNC.EQ.-67.OR.IHOTNC.EQ.-68) ) then
         NCIESTP=IESTP
         NCCOUE=NSCOUE
      ENDIF
    
      IF(NOUTE.LT.0) THEN
        IESTP=0
        NSCOUE=0
        IF((NTCYSE.LT.ITHS).AND.(NSPOOLE.GT.0)) THEN
          NTCYSE=NTCYSE+((REAL(ITHS)-REAL(NTCYSE))/
     &           REAL(NSPOOLE))*REAL(NSPOOLE)
          IF(NTCYSE.LT.ITHS) NTCYSE=NTCYSE+NSPOOLE
          IF(NSPOOLE.NE.0) NTRSPE=(NTCYFE-NTCYSE)/NSPOOLE
        ENDIF
        WRITE(16,1041)
 1041   FORMAT(//,' A NEW ELEVATION STATION FILE WILL BE STARTED')
      ENDIF
      IF(NOUTE.EQ.-2) THEN
        OPEN(61,FILE=TRIM(LOCALDIR)//'/'//'fort.61',
     &                      ACCESS='DIRECT',RECL=NBYTE)
        IF(NBYTE.EQ.4) THEN
          DO I=1,8
            WRITE(61,REC=IESTP+I) RDES4(I)
          ENDDO
          IESTP=IESTP+8
          DO I=1,6
            WRITE(61,REC=IESTP+I) RID4(I)
          ENDDO
          IESTP=IESTP+6
          DO I=1,6
            WRITE(61,REC=IESTP+I) AID4(I)
          ENDDO
          IESTP=IESTP+6
        ENDIF
        IF(NBYTE.EQ.8) THEN
          DO I=1,4
            WRITE(61,REC=IESTP+I) RDES8(I)
          ENDDO
          IESTP=IESTP+4
          DO I=1,3
            WRITE(61,REC=IESTP+I) RID8(I)
          ENDDO
          IESTP=IESTP+3
          DO I=1,3
            WRITE(61,REC=IESTP+I) AID8(I)
          ENDDO
          IESTP=IESTP+3
        ENDIF
        WRITE(61,REC=IESTP+1) NTRSPE
        WRITE(61,REC=IESTP+2) NSTAE
        WRITE(61,REC=IESTP+3) DT*NSPOOLE
        WRITE(61,REC=IESTP+4) NSPOOLE
        WRITE(61,REC=IESTP+5) 1
        IESTP=IESTP+5
        CLOSE(61)                    ! DO THIS TO FLUSH THE WRITE BUFFER
      ENDIF
      IF(NOUTE.EQ.-1) THEN
        CALL OPEN_GBL_FILE(61, TRIM(GLOBALDIR)//'/'//'fort.61',
     $    NSTAE_G, NSTAE, HEADER61) 
        IESTP=2
      ENDIF

      IF(NOUTE.EQ.2) THEN
        OPEN(61,FILE=TRIM(LOCALDIR)//'/'//'fort.61',
     &                      ACCESS='DIRECT',RECL=NBYTE)
c       WRITE(61,REC=ITEMPSTP+1) NTRSPE   ! ALLOW ADDITIONAL OUTPUT DATA TO BE WRITTEN
        CLOSE(61)                         ! DO THIS TO FLUSH THE WRITE BUFFER
      ENDIF

C...
C....GO TO THE PROPER PLACE IN THE VELOCITY STATION OUTPUT FILE
C...

      IF (IHOT.EQ.67.OR.IHOT.EQ.68.OR.
     &    (ncmode.eq.2. and. (IHOT.EQ.-67.OR.IHOT.EQ.-68) ) )THEN

      READ(IHOT,REC=IHOTSTP) IVSTP ; IHOTSTP = IHOTSTP + 1
      READ(IHOT,REC=IHOTSTP) NSCOUV ; IHOTSTP = IHOTSTP + 1
      WRITE(16,1042) IVSTP,NSCOUV
 1042 FORMAT(//,1X,I6,' LINES OR RECORDS WRITTEN IN VELOCITY ',
     &                'STATION FILE BY THE TIME OF THE HOT START',
     &           /,8X,'SPOOL COUNTER =',I6)
     
      ENDIF
     
      IF (ncmode.ne.0. and. (IHOTNC.EQ.-67.OR.IHOTNC.EQ.-68) ) then
         NCIVSTP=IVSTP
         NCCOUV=NSCOUV
      ENDIF
     
      IF(NOUTV.LT.0) THEN
        IVSTP=0
        NSCOUV=0
        IF((NTCYSV.LT.ITHS).AND.(NSPOOLV.GT.0)) THEN
          NTCYSV=NTCYSV+((REAL(ITHS)-REAL(NTCYSV))/
     &           REAL(NSPOOLV))*REAL(NSPOOLV)
          IF(NTCYSV.LT.ITHS) NTCYSV=NTCYSV+NSPOOLV
          NTRSPV=(NTCYFV-NTCYSV)/NSPOOLV
        ENDIF
        WRITE(16,1043)
 1043   FORMAT(//,' A NEW VELOCITY STATION FILE WILL BE STARTED')
      ENDIF

      IF(NOUTV.EQ.-2) THEN
        OPEN(62,FILE=TRIM(LOCALDIR)//'/'//'fort.62',
     &                      ACCESS='DIRECT',RECL=NBYTE)
        IF(NBYTE.EQ.4) THEN
          DO I=1,8
            WRITE(62,REC=IVSTP+I) RDES4(I)
          ENDDO
          IVSTP=IVSTP+8
          DO I=1,6
            WRITE(62,REC=IVSTP+I) RID4(I)
          ENDDO
          IVSTP=IVSTP+6
          DO I=1,6
            WRITE(62,REC=IVSTP+I) AID4(I)
          ENDDO
          IVSTP=IVSTP+6
          ENDIF
        IF(NBYTE.EQ.8) THEN
          DO I=1,4
            WRITE(62,REC=IVSTP+I) RDES8(I)
          ENDDO
          IVSTP=IVSTP+4
          DO I=1,3
            WRITE(62,REC=IVSTP+I) RID8(I)
          ENDDO
          IVSTP=IVSTP+3
          DO I=1,3
            WRITE(62,REC=IVSTP+I) AID8(I)
          ENDDO
          IVSTP=IVSTP+3
        ENDIF
        WRITE(62,REC=IVSTP+1) NTRSPV
        WRITE(62,REC=IVSTP+2) NSTAV
        WRITE(62,REC=IVSTP+3) DT*NSPOOLV
        WRITE(62,REC=IVSTP+4) NSPOOLV
        WRITE(62,REC=IVSTP+5) 2
        IVSTP=IVSTP+5
        CLOSE(62)                    ! DO THIS TO FLUSH THE WRITE BUFFER
      ENDIF
      IF(NOUTV.EQ.-1) THEN
        CALL OPEN_GBL_FILE(62, TRIM(GLOBALDIR)//'/'//'fort.62',
     $    NSTAV_G, NSTAV, HEADER62)
        IVSTP=2
      ENDIF
      IF(NOUTV.EQ.2) THEN
        OPEN(62,FILE=TRIM(LOCALDIR)//'/'//'fort.62',
     &                      ACCESS='DIRECT',RECL=NBYTE)
c       WRITE(62,REC=ITEMPSTP+1) NTRSPV   ! ALLOW ADDITIONAL OUTPUT DATA TO BE WRITTEN
        CLOSE(62)                         ! DO THIS TO FLUSH THE WRITE BUFFER
      ENDIF

C...
C....GO TO THE PROPER PLACE IN THE CONCENTRATION STATION OUTPUT FILE
C...

      IF (IHOT.EQ.67.OR.IHOT.EQ.68. OR.
     &    (ncmode.eq.2. and. (IHOT.EQ.-67.OR.IHOT.EQ.-68) ) )THEN

      READ(IHOT,REC=IHOTSTP) ICSTP ; IHOTSTP = IHOTSTP + 1
      READ(IHOT,REC=IHOTSTP) NSCOUC ; IHOTSTP = IHOTSTP + 1
      WRITE(16,1044) ICSTP,NSCOUC
 1044 FORMAT(//,1X,I6,' LINES OR RECORDS WRITTEN IN CONCENTRATION ',
     &                'STATION FILE BY THE TIME OF THE HOT START',
     &           /,8X,'SPOOL COUNTER = ',I6)
     
      ENDIF
      
      IF(NOUTC.LT.0) THEN
        ICSTP=0
        NSCOUC=0
        IF((NTCYSC.LT.ITHS).AND.(NSPOOLC.GT.0)) THEN
          NTCYSC=NTCYSC+((REAL(ITHS)-REAL(NTCYSC))/
     &           REAL(NSPOOLC))*REAL(NSPOOLC)
          IF(NTCYSC.LT.ITHS) NTCYSC=NTCYSC+NSPOOLC
          NTRSPC=(NTCYFC-NTCYSC)/NSPOOLC
        ENDIF
        WRITE(16,1045)
 1045   FORMAT(//,' A NEW CONCENTRATION STATION FILE WILL BE STARTED')
      ENDIF

      IF(NOUTC.EQ.-2) THEN
        OPEN(81,FILE=TRIM(LOCALDIR)//'/'//'fort.81',
     &                      ACCESS='DIRECT',RECL=NBYTE)
        IF(NBYTE.EQ.4) THEN
          DO I=1,8
            WRITE(81,REC=ICSTP+I) RDES4(I)
          ENDDO
          ICSTP=ICSTP+8
          DO I=1,6
            WRITE(81,REC=ICSTP+I) RID4(I)
          ENDDO
          ICSTP=ICSTP+6
          DO I=1,6
            WRITE(81,REC=ICSTP+I) AID4(I)
          ENDDO
          ICSTP=ICSTP+6
        ENDIF
        IF(NBYTE.EQ.8) THEN
          DO I=1,4
            WRITE(81,REC=ICSTP+I) RDES8(I)
          ENDDO
          ICSTP=ICSTP+4
          DO I=1,3
            WRITE(81,REC=ICSTP+I) RID8(I)
          ENDDO
          ICSTP=ICSTP+3
          DO I=1,3
            WRITE(81,REC=ICSTP+I) AID8(I)
          ENDDO
          ICSTP=ICSTP+3
        ENDIF
        WRITE(81,REC=ICSTP+1) NTRSPC
        WRITE(81,REC=ICSTP+2) NSTAC
        WRITE(81,REC=ICSTP+3) DT*NSPOOLC
        WRITE(81,REC=ICSTP+4) NSPOOLC
        WRITE(81,REC=ICSTP+5) 1
        ICSTP=ICSTP+5
        CLOSE(81)                    ! DO THIS TO FLUSH THE WRITE BUFFER
      ENDIF
      IF(NOUTC.EQ.-1) THEN
        CALL OPEN_GBL_FILE(81, TRIM(GLOBALDIR)//'/'//'fort.81',
     $    NSTAC_G, NSTAC, HEADER81) 
        ICSTP=2
        CLOSE(81)
      ENDIF
      IF(NOUTC.EQ.2) THEN
        OPEN(81,FILE=TRIM(LOCALDIR)//'/'//'fort.81',
     &                      ACCESS='DIRECT',RECL=NBYTE)
c       WRITE(81,REC=ITEMPSTP+1) NTRSPC   ! ALLOW ADDITIONAL OUTPUT DATA TO BE WRITTEN
        CLOSE(81)                         ! DO THIS TO FLUSH THE WRITE BUFFER
      ENDIF

C...
C....GO TO THE PROPER PLACE IN THE METEOROLOGICAL STATION OUTPUT FILE
C...

      IF (IHOT.EQ.67.OR.IHOT.EQ.68.OR.
     &    (ncmode.eq.2. and. (IHOT.EQ.-67.OR.IHOT.EQ.-68) ) )THEN

      READ(IHOT,REC=IHOTSTP) IPSTP ; IHOTSTP = IHOTSTP + 1
      READ(IHOT,REC=IHOTSTP) IWSTP ; IHOTSTP = IHOTSTP + 1
      READ(IHOT,REC=IHOTSTP) NSCOUM ; IHOTSTP = IHOTSTP + 1
      WRITE(16,1038) IWSTP,IPSTP,NSCOUM
 1038 FORMAT(//,1X,I6,' LINES OR RECORDS WRITTEN IN THE WIND STATION',
     &                ' FILE BY THE TIME OF THE HOT START',
     &        /,1X,I6,' LINES OR RECORDS WRITTIN IN THE PRES STATION',
     &                ' FILE BY THE TMIE OF THE HOT START',
     &           /,8X,'SPOOL COUNTER = ',I6)
     
      ENDIF
      IF (ncmode.ne.0. and. (IHOTNC.EQ.-67.OR.IHOTNC.EQ.-68) ) then
         NCIWSTP=IWSTP
         NCIPSTP=IPSTP
         NCCOUM=NSCOUM
      ENDIF

     
      IF(NOUTM.LT.0) THEN
        IPSTP=0
        IWSTP=0
        NSCOUM=0
        IF((NTCYSM.LT.ITHS).AND.(NSPOOLM.GT.0)) THEN
          NTCYSM=NTCYSM+((REAL(ITHS)-REAL(NTCYSM))/
     &           REAL(NSPOOLM))*REAL(NSPOOLM)
          IF(NTCYSM.LT.ITHS) NTCYSM=NTCYSM+NSPOOLM
          NTRSPM=(NTCYFM-NTCYSM)/NSPOOLM
        ENDIF
        WRITE(16,1039)
 1039   FORMAT(//,' A NEW METEOROLOGICAL STATION FILE WILL BE STARTED')
      ENDIF

      IF(NOUTM.EQ.-2) THEN
        OPEN(71,FILE=TRIM(LOCALDIR)//'/'//'fort.71',
     &                      ACCESS='DIRECT',RECL=NBYTE)
        OPEN(72,FILE=TRIM(LOCALDIR)//'/'//'fort.72',
     &                      ACCESS='DIRECT',RECL=NBYTE)
        IF(NBYTE.EQ.4) THEN
          DO I=1,8
            WRITE(71,REC=IPSTP+I) RDES4(I)
            WRITE(72,REC=IWSTP+I) RDES4(I)
          ENDDO
          IPSTP=IPSTP+8
          IWSTP=IWSTP+8
          DO I=1,6
            WRITE(71,REC=IPSTP+I) RID4(I)
            WRITE(72,REC=IWSTP+I) RID4(I)
          ENDDO
          IPSTP=IPSTP+6
          IWSTP=IWSTP+6
          DO I=1,6
            WRITE(71,REC=IPSTP+I) AID4(I)
            WRITE(72,REC=IWSTP+I) AID4(I)
          ENDDO
          IPSTP=IPSTP+6
          IWSTP=IWSTP+6
        ENDIF
        IF(NBYTE.EQ.8) THEN
          DO I=1,4
            WRITE(71,REC=IPSTP+I) RDES8(I)
            WRITE(72,REC=IWSTP+I) RDES8(I)
          ENDDO
          IPSTP=IPSTP+4
          IWSTP=IWSTP+4
          DO I=1,3
            WRITE(71,REC=IPSTP+I) RID8(I)
            WRITE(72,REC=IWSTP+I) RID8(I) 
          ENDDO
          IPSTP=IPSTP+3
          IWSTP=IWSTP+3
          DO I=1,3
            WRITE(71,REC=IPSTP+I) AID8(I)
            WRITE(72,REC=IWSTP+I) AID8(I)
          ENDDO
          IPSTP=IPSTP+3
          IWSTP=IWSTP+3
        ENDIF
        WRITE(71,REC=IPSTP+1) NTRSPM
        WRITE(71,REC=IPSTP+2) NSTAM
        WRITE(71,REC=IPSTP+3) DT*NSPOOLM
        WRITE(71,REC=IPSTP+4) NSPOOLM
        WRITE(71,REC=IPSTP+5) 1
        WRITE(72,REC=IWSTP+1) NTRSPM
        WRITE(72,REC=IWSTP+2) NSTAM
        WRITE(72,REC=IWSTP+3) DT*NSPOOLM
        WRITE(72,REC=IWSTP+4) NSPOOLM
        WRITE(72,REC=IWSTP+5) 2
        IPSTP=IPSTP+5
        IWSTP=IWSTP+5
        CLOSE(71)                    ! DO THIS TO FLUSH THE WRITE BUFFER
        CLOSE(72)
      ENDIF
      IF(NOUTM.EQ.-1) THEN
        CALL OPEN_GBL_FILE(71, TRIM(GLOBALDIR)//'/'//'fort.71',
     $    NSTAM_G, NSTAM, HEADER71)
        IPSTP=2
        CALL OPEN_GBL_FILE(72, TRIM(GLOBALDIR)//'/'//'fort.72',
     $    NSTAM_G, NSTAM, HEADER72)
        IWSTP=2
      ENDIF
      IF(NOUTM.EQ.2) THEN
        OPEN(71,FILE=TRIM(LOCALDIR)//'/'//'fort.71',
     &                      ACCESS='DIRECT',RECL=NBYTE)
        OPEN(72,FILE=TRIM(LOCALDIR)//'/'//'fort.72',
     &                      ACCESS='DIRECT',RECL=NBYTE)
c         WRITE(71,REC=ITEMPSTP+1) NTRSPM   ! ALLOW ADDITIONAL OUTPUT DATA TO BE WRITTEN
c         WRITE(72,REC=ITMEPSTP+1) NTRSPM
        CLOSE(71)                         ! DO THIS TO FLUSH THE WRITE BUFFER
        CLOSE(72)
      ENDIF

C...
C....GO TO THE PROPER PLACE IN THE GLOBAL ELEVATION OUTPUT FILE
C...

      IF (IHOT.EQ.67.OR.IHOT.EQ.68. OR.
     &    (ncmode.eq.2. and. (IHOT.EQ.-67.OR.IHOT.EQ.-68) ) )THEN

!      PRINT *,"START READING HOT START FILE - BINARY IGEP,NSCOUGE"


        READ(IHOT,REC=IHOTSTP) IGEP ; IHOTSTP = IHOTSTP + 1
        READ(IHOT,REC=IHOTSTP) NSCOUGE ; IHOTSTP = IHOTSTP + 1
        WRITE(16,1046) IGEP,NSCOUGE
 1046   FORMAT(//,1X,I6,' LINES OR RECORDS WRITTEN IN THE GLOBAL ',
     &                  'ELEVATION FILE BY THE TIME OF THE HOT START',
     &          /,8X,'SPOOL COUNTER =',I6)
     
!         IF (MYPROC.EQ.0) PRINT *,"B-HS: IGEP == ", IGEP
!         IF (MYPROC.EQ.0) PRINT *,"B-HS: NSCOUGE == ", NSCOUGE
!         IF (MYPROC.EQ.0) PRINT *,"B-HS: NCIGEP == ", NCIGEP
!         IF (MYPROC.EQ.0) PRINT *,"B-HS: NCCOUGE == ", NCCOUGE
!         IF (MYPROC.EQ.0) PRINT *,"B-HS: NCMODE == ", ncmode
!         IF (MYPROC.EQ.0) PRINT *,"B-HS: IHOT == ", IHOT

       ENDIF
      IF (ncmode.ne.0. and. (IHOTNC.EQ.-67.OR.IHOTNC.EQ.-68) ) then
         NCIGEP=IGEP
         NCCOUGE=NSCOUGE
!         IF (MYPROC.EQ.0) PRINT *,"N-HS: NCIGEP == ", NCIGEP
!         IF (MYPROC.EQ.0) PRINT *,"N-HS: NCCOUGE == ", NCCOUGE
      ENDIF
      
        IF(NOUTGE.LT.0) THEN
          IGEP=0
          NSCOUGE=0
          IF((NTCYSGE.LT.ITHS).AND.(NSPOOLGE.GT.0)) THEN
            NTCYSGE=NTCYSGE+((REAL(ITHS)-REAL(NTCYSGE))/
     &              REAL(NSPOOLGE))*REAL(NSPOOLGE)
            IF(NTCYSGE.LT.ITHS) NTCYSGE=NTCYSGE+NSPOOLGE
C     kmd48.33bc - changed NDSETSE for use of new DTDPHS
            IF(DTDP.NE.DTDPHS) THEN 
              NDSETSE=((NTCYFGE*DTDP)-(NTCYSGE*DTDPHS))/
     &                (NSPOOLGE*DTDP)
            ELSE
              NDSETSE=(NTCYFGE-NTCYSGE)/NSPOOLGE
            END IF 
            ENDIF
          WRITE(16,1047)
 1047     FORMAT(//,' A NEW GLOBAL ELEVATION FILE WILL BE STARTED')
          ENDIF

        IF(NOUTGE.EQ.-2) THEN
          OPEN(63,FILE=TRIM(LOCALDIR)//'/'//'fort.63',
     &          ACCESS='DIRECT',RECL=NBYTE)
          IF(NBYTE.EQ.4) THEN
            DO I=1,8
              WRITE(63,REC=IGEP+I) RDES4(I)
              ENDDO
            IGEP=IGEP+8
            DO I=1,6
              WRITE(63,REC=IGEP+I) RID4(I)
              ENDDO
            IGEP=IGEP+6
            DO I=1,6
              WRITE(63,REC=IGEP+I) AID4(I)
              ENDDO
            IGEP=IGEP+6
            ENDIF
          IF(NBYTE.EQ.8) THEN
            DO I=1,4
              WRITE(63,REC=IGEP+I) RDES8(I)
              ENDDO
            IGEP=IGEP+4
            DO I=1,3
              WRITE(63,REC=IGEP+I) RID8(I)
              ENDDO
            IGEP=IGEP+3
            DO I=1,3
              WRITE(63,REC=IGEP+I) AID8(I)
              ENDDO
            IGEP=IGEP+3
            ENDIF
          WRITE(63,REC=IGEP+1) NDSETSE
          WRITE(63,REC=IGEP+2) NP
          WRITE(63,REC=IGEP+3) DT*NSPOOLGE
          WRITE(63,REC=IGEP+4) NSPOOLGE
          WRITE(63,REC=IGEP+5) 1
          IGEP=IGEP+5
          CLOSE(63)                    ! DO THIS TO FLUSH THE WRITE BUFFER
        ENDIF

Casey 090310: Added the possibility for NOUTGE = -4.
        IF((NOUTGE.EQ.-1).OR.(NOUTGE.EQ.-4)) THEN
           CALL OPEN_GBL_FILE(63, TRIM(GLOBALDIR)//'/'//'fort.63',
     $         NP_G, NP, HEADER63)
C
C        jgf47.06 Tau0 output is produced on the same schedule as
C        global elevation
           IF (OUTPUTTAU0) THEN
Casey 090310: Corrected the unit number from 10 to 90.
               CALL writeDomainHeader(90,
     &              TRIM(GLOBALDIR)//'/'//'fort.tau0',
     $               NP_G, NP, 'Tau0      ')
           ENDIF

C kmd48.33bc - information for writing out the sponge layer
           IF(OUTPUTSPONGE) THEN 
              CALL writeDomainHeader(92,
     &             TRIM(GLOBALDIR)//'/'//'fort.92',
     $             NP_G, NP, 'sponge      ')
           END IF 

           IGEP=2
        ENDIF
        IF(NOUTGE.EQ.2) THEN
          OPEN(63,FILE=TRIM(LOCALDIR)//'/'//'fort.63',
     &          ACCESS='DIRECT',RECL=NBYTE)
c         WRITE(63,REC=ITEMPSTP+1) NDSETSE  ! ALLOW ADDITIONAL OUTPUT DATA TO BE WRITTEN
          CLOSE(63)                         ! DO THIS TO FLUSH THE WRITE BUFFER
        ENDIF

C...
C....GO TO THE PROPER PLACE IN THE GLOBAL VELOCITY OUTPUT FILE
C...

        IF (IHOT.EQ.67.OR.IHOT.EQ.68. OR.
     &    (ncmode.eq.2. and. (IHOT.EQ.-67.OR.IHOT.EQ.-68) ) )THEN

!        PRINT *,"START READING HOT START FILE - BINARY IGVP,NSCOUGV"



        READ(IHOT,REC=IHOTSTP) IGVP ; IHOTSTP = IHOTSTP + 1
        READ(IHOT,REC=IHOTSTP) NSCOUGV ; IHOTSTP = IHOTSTP + 1
        WRITE(16,1048) IGVP,NSCOUGV
 1048   FORMAT(//,1X,I6,' LINES OR RECORDS WRITTEN IN THE GLOBAL ',
     &                  'VELOCITY FILE BY THE TIME OF THE HOT START',
     &          /,8X,'SPOOL COUNTER =',I6)
     
        ENDIF
      IF (ncmode.ne.0. and. (IHOTNC.EQ.-67.OR.IHOTNC.EQ.-68) ) then
         NCIGVP=IGVP
         NCCOUGV=NSCOUGV
      ENDIF
     
        IF(NOUTGV.LT.0) THEN
          IGVP=0
          NSCOUGV=0
          IF((NTCYSGV.LT.ITHS).AND.(NSPOOLGV.GT.0)) THEN
            NTCYSGV=NTCYSGV+((REAL(ITHS)-REAL(NTCYSGV))/
     &              REAL(NSPOOLGV))*REAL(NSPOOLGV)
            IF(NTCYSGV.LT.ITHS) NTCYSGV=NTCYSGV+NSPOOLGV
C  kmd48.33bc - change the NDSETSV for the change in the DTDP
            IF(DTDP.NE.DTDPHS) THEN 
              NDSETSV=((NTCYFGV*DTDP)-(NTCYSGV*DTDPHS))/
     &                (NSPOOLGV*DTDP)
            ELSE 
              NDSETSV=(NTCYFGV-NTCYSGV)/NSPOOLGV
            END IF 
            ENDIF
          WRITE(16,1049)
 1049     FORMAT(//,' A NEW GLOBAL VELOCITY FILE WILL BE STARTED')
          ENDIF

        IF(NOUTGV.EQ.-2) THEN
          OPEN(64,FILE=TRIM(LOCALDIR)//'/'//'fort.64',
     &          ACCESS='DIRECT',RECL=NBYTE)
          IF(NBYTE.EQ.4) THEN
            DO I=1,8
              WRITE(64,REC=IGVP+I) RDES4(I)
              ENDDO
            IGVP=IGVP+8
            DO I=1,6
              WRITE(64,REC=IGVP+I) RID4(I)
              ENDDO
            IGVP=IGVP+6
            DO I=1,6
              WRITE(64,REC=IGVP+I) AID4(I)
              ENDDO
            IGVP=IGVP+6
            ENDIF
          IF(NBYTE.EQ.8) THEN
            DO I=1,4
              WRITE(64,REC=IGVP+I) RDES8(I)
              ENDDO
            IGVP=IGVP+4
            DO I=1,3
              WRITE(64,REC=IGVP+I) RID8(I)
              ENDDO
            IGVP=IGVP+3
            DO I=1,3
              WRITE(64,REC=IGVP+I) AID8(I)
              ENDDO
            IGVP=IGVP+3
            ENDIF
          WRITE(64,REC=IGVP+1) NDSETSV
          WRITE(64,REC=IGVP+2) NP
          WRITE(64,REC=IGVP+3) DT*NSPOOLGV
          WRITE(64,REC=IGVP+4) NSPOOLGV
          WRITE(64,REC=IGVP+5) 2
          IGVP=IGVP+5
          CLOSE(64)                    ! DO THIS TO FLUSH THE WRITE BUFFER
        ENDIF
        IF((NOUTGV.EQ.-1).OR.(NOUTGV.EQ.-4)) THEN
          CALL OPEN_GBL_FILE(64, TRIM(GLOBALDIR)//'/'//'fort.64',
     $      NP_G, NP, HEADER64)
          IGVP=2
        ENDIF
        IF(NOUTGV.EQ.2) THEN
          OPEN(64,FILE=TRIM(LOCALDIR)//'/'//'fort.64',
     &          ACCESS='DIRECT',RECL=NBYTE)
c         WRITE(64,REC=ITEMPSTP+1) NDSETSV  ! ALLOW ADDITIONAL OUTPUT DATA TO BE WRITTEN
          CLOSE(64)                         ! DO THIS TO FLUSH THE WRITE BUFFER
        ENDIF

C...
C....GO TO THE PROPER PLACE IN THE GLOBAL CONCENTRATION OUTPUT FILE
C...

        IF (IHOT.EQ.67.OR.IHOT.EQ.68)THEN

!        PRINT *,"START READING HOT START FILE - BINARY IGCP,NSCOUGC"

        READ(IHOT,REC=IHOTSTP) IGCP ; IHOTSTP = IHOTSTP + 1
        READ(IHOT,REC=IHOTSTP) NSCOUGC ; IHOTSTP = IHOTSTP + 1
        WRITE(16,1053) IGCP,NSCOUGC
 1053   FORMAT(//,1X,I6,' LINES OR RECORDS WRITTEN IN THE GLOBAL ',
     &                'CONCENTRATION FILE BY THE TIME OF THE HOT START',
     &          /,8X,'SPOOL COUNTER =',I6)
     
        ENDIF
     
        IF(NOUTGC.LT.0) THEN
          IGCP=0
          NSCOUGC=0
          IF((NTCYSGC.LT.ITHS).AND.(NSPOOLGC.GT.0)) THEN
            NTCYSGC=NTCYSGC+((REAL(ITHS)-REAL(NTCYSGC))/
     &              REAL(NSPOOLGC))*REAL(NSPOOLGC)
            IF(NTCYSGC.LT.ITHS) NTCYSGC=NTCYSGC+NSPOOLGC
C  kmd48.33bc - change the NDSETSC for the change in the DTDP
            IF(DTDP.NE.DTDPHS) THEN 
              NDSETSC=((NTCYFGC*DTDP)-(NTCYSGC*DTDPHS))/
     &                (NSPOOLGC*DTDP)
            ELSE 
              NDSETSC=(NTCYFGC-NTCYSGC)/NSPOOLGC
            END IF 
            ENDIF
          WRITE(16,1054)
 1054     FORMAT(//,' A NEW GLOBAL CONCENTRATION FILE WILL BE STARTED')
        ENDIF

        IF(NOUTGC.EQ.-2) THEN
          OPEN(83,FILE=TRIM(LOCALDIR)//'/'//'fort.83',
     &           ACCESS='DIRECT',RECL=NBYTE)
          IF(NBYTE.EQ.4) THEN
            DO I=1,8
              WRITE(83,REC=IGCP+I) RDES4(I)
              ENDDO
            IGCP=IGCP+8
            DO I=1,6
              WRITE(83,REC=IGCP+I) RID4(I)
              ENDDO
            IGCP=IGCP+6
            DO I=1,6
              WRITE(83,REC=IGCP+I) AID4(I)
              ENDDO
            IGCP=IGCP+6
            ENDIF
          IF(NBYTE.EQ.8) THEN
            DO I=1,4
              WRITE(83,REC=IGCP+I) RDES8(I)
              ENDDO
            IGCP=IGCP+4
            DO I=1,3
              WRITE(83,REC=IGCP+I) RID8(I)
              ENDDO
            IGCP=IGCP+3
            DO I=1,3
              WRITE(83,REC=IGCP+I) AID8(I)
              ENDDO
            IGCP=IGCP+3
            ENDIF
          WRITE(83,REC=IGCP+1) NDSETSC
          WRITE(83,REC=IGCP+2) NP
          WRITE(83,REC=IGCP+3) DT*NSPOOLGC
          WRITE(83,REC=IGCP+4) NSPOOLGC
          WRITE(83,REC=IGCP+5) 1
          IGCP=IGCP+5
          CLOSE(83)                    ! DO THIS TO FLUSH THE WRITE BUFFER
        ENDIF
        IF(NOUTGC.EQ.-1) THEN
          OPEN(83,FILE=TRIM(LOCALDIR)//'/'//'fort.83')
          WRITE(83,3220) RUNDES,RUNID,AGRID
          WRITE(83,3645) NDSETSC,NP,DTDP*NSPOOLGC,NSPOOLGC,1
          IGCP=2
          CLOSE(83)
        ENDIF
        IF(NOUTGC.EQ.1) THEN
          OPEN(83,FILE=TRIM(LOCALDIR)//'/'//'fort.83')
          DO I=1,IGCP           !I DON'T KNOW OF A PRACTICAL WAY TO CHANGE NDSETSC
            READ(83,1050)
 1050       FORMAT(1X)
          ENDDO
          ENDFILE(83)
          CLOSE(83)
        ENDIF
        IF(NOUTGC.EQ.2) THEN
          OPEN(83,FILE=TRIM(LOCALDIR)//'/'//'fort.83',
     &          ACCESS='DIRECT',RECL=NBYTE)
c         WRITE(83,REC=ITEMPSTP+1) NDSETSC  ! ALLOW ADDITIONAL OUTPUT DATA TO BE WRITTEN
          CLOSE(83)                         ! DO THIS TO FLUSH THE WRITE BUFFER
        ENDIF

C...
C....GO TO THE PROPER PLACE IN THE GLOBAL METEOROLOGICAL OUTPUT FILES
C...
      IF (IHOT.EQ.67.OR.IHOT.EQ.68. OR.
     &    (ncmode.eq.2. and. (IHOT.EQ.-67.OR.IHOT.EQ.-68) ) )THEN
!        PRINT *,"START READING HOT START FILE- BINARY IGWP,IGPP,NSCOUGW"

        READ(IHOT,REC=IHOTSTP) IGPP ; IHOTSTP = IHOTSTP + 1
        READ(IHOT,REC=IHOTSTP) IGWP ; IHOTSTP = IHOTSTP + 1
        READ(IHOT,REC=IHOTSTP) NSCOUGW ; IHOTSTP = IHOTSTP + 1
        WRITE(16,1055) IGWP,IGPP,NSCOUGW
 1055   FORMAT(//,1X,I6,' LINES OR RECORDS WRITTEN IN THE GLOBAL ',
     &                  'WIND FILE BY THE TIME OF THE HOT START',
     &          /,1X,I6,'LINES OR RECORDS WRITTEN IN THE GLOBAL ',
     &                  'PRESSURE FILE BY THE TIME OF THE HOT START',   
     &          /,8X,'SPOOL COUNTER =',I6)
     
      ENDIF
      IF (ncmode.ne.0. and. (IHOTNC.EQ.-67.OR.IHOTNC.EQ.-68) ) then
         NCIGWP=IGWP
         NCIGPP=IGPP
         NCCOUGW=NSCOUGW
      ENDIF
     
        IF(NOUTGW.LT.0) THEN
          IGPP=0
          IGWP=0
          NSCOUGW=0
          IF((NTCYSGW.LT.ITHS).AND.(NSPOOLGW.GT.0)) THEN
            NTCYSGW=NTCYSGW+((REAL(ITHS)-REAL(NTCYSGW))/
     &              REAL(NSPOOLGW))*REAL(NSPOOLGW)
            IF(NTCYSGW.LT.ITHS) NTCYSGW=NTCYSGW+NSPOOLGW
C  kmd48.33bc - change the NDSETSC for the change in the DTDP
            IF(DTDP.NE.DTDPHS) THEN 
              NDSETSW=((NTCYFGW*DTDP)-(NTCYSGW*DTDPHS))/
     &                (NSPOOLGW*DTDP)
            ELSE 
              NDSETSW=(NTCYFGW-NTCYSGW)/NSPOOLGW
            END IF 
          ENDIF
          WRITE(16,1056)
 1056     FORMAT(//,' NEW GLOBAL WIND & pressure FILEs WILL BE STARTED')
        ENDIF

        IF(NOUTGW.EQ.-2) THEN
          OPEN(73,file=trim(LOCALDIR)//'/'//'fort.73',
     &          access='DIRECT',recl=nbyte)
          OPEN(74,FILE=TRIM(LOCALDIR)//'/'//'fort.74',
     &          ACCESS='DIRECT',RECL=NBYTE)
          IF(NBYTE.EQ.4) THEN
            DO I=1,8
              write(73,rec=igpp+i) rdes4(i)
              WRITE(74,REC=IGWP+I) RDES4(I)
              ENDDO
            igpp=igpp+8
            IGWP=IGWP+8
            DO I=1,6
              write(73,rec=igpp+i) rid4(i)
              WRITE(74,REC=IGWP+I) RID4(I)
              ENDDO
            igpp=igpp+6
            IGWP=IGWP+6
            DO I=1,6
              write(73,rec=igpp+i) aid4(i)
              WRITE(74,REC=IGWP+I) AID4(I)
              ENDDO
            igpp=igpp+6
            IGWP=IGWP+6
            ENDIF
          IF(NBYTE.EQ.8) THEN
            DO I=1,4
              write(73,rec=igpp+i) rdes8(i)
              WRITE(74,REC=IGWP+I) RDES8(I)
              ENDDO
            igpp=igpp+4
            IGWP=IGWP+4
            DO I=1,3
              write(73,rec=igpp+i) rid8(i)
              WRITE(74,REC=IGWP+I) RID8(I)
              ENDDO
            igpp=igpp+3
            IGWP=IGWP+3
            DO I=1,3
              write(73,rec=igpp+i) aid8(i)
              WRITE(74,REC=IGWP+I) AID8(I)
              ENDDO
            igpp=igpp+3
            IGWP=IGWP+3
            ENDIF
          WRITE(73,rec=igpp+1) ndsetsw
          WRITE(73,rec=igpp+2) np
          WRITE(73,rec=igpp+3) dt*nspoolgw
          WRITE(73,rec=igpp+4) nspoolgw
          WRITE(73,rec=igpp+5) 2
          IGPP=IGPP+5
          CLOSE(73)                    ! DO THIS TO FLUSH THE WRITE BUFFER
          WRITE(74,REC=IGWP+1) NDSETSW
          WRITE(74,REC=IGWP+2) NP
          WRITE(74,REC=IGWP+3) DT*NSPOOLGW
          WRITE(74,REC=IGWP+4) NSPOOLGW
          WRITE(74,REC=IGWP+5) 2
          IGWP=IGWP+5
          CLOSE(74)                    ! DO THIS TO FLUSH THE WRITE BUFFER
        ENDIF

Casey 090310: Added the possibility for NOUTGW = -4.
        IF((NOUTGW.EQ.-1).OR.(NOUTGW.EQ.-4)) THEN
          CALL OPEN_GBL_FILE(73, TRIM(GLOBALDIR)//'/'//'fort.73',
     $      NP_G, NP, HEADER73)
          IGPP=2
          CALL OPEN_GBL_FILE(74, TRIM(GLOBALDIR)//'/'//'fort.74',
     $      NP_G, NP, HEADER74)
          IGWP=2
        ENDIF

#ifdef CSWAN
Casey 090310: Added the output of radiation stress gradients.
        IF(NOUTGW.LT.0)THEN
          IGRadS=0
        ENDIF
        IF(NOUTGW.EQ.-2)THEN
          OPEN(164,FILE=TRIM(LOCALDIR)//'/'//'rads.64',
     &          ACCESS='DIRECT',RECL=NBYTE)
          IF(NBYTE.EQ.4) THEN
            DO I=1,8
              WRITE(164,REC=IGRadS+I) RDES4(I)
              ENDDO
            IGRadS=IGRadS+8
            DO I=1,6
              WRITE(164,REC=IGRadS+I) RID4(I)
              ENDDO
            IGRadS=IGRadS+6
            DO I=1,6
              WRITE(164,REC=IGRadS+I) AID4(I)
              ENDDO
            IGRadS=IGRadS+6
            ENDIF
          IF(NBYTE.EQ.8) THEN
            DO I=1,4
              WRITE(164,REC=IGRadS+I) RDES8(I)
              ENDDO
            IGRadS=IGRadS+4
            DO I=1,3
              WRITE(164,REC=IGRadS+I) RID8(I)
              ENDDO
            IGRadS=IGRadS+3
            DO I=1,3
              WRITE(164,REC=IGRadS+I) AID8(I)
              ENDDO
            IGRadS=IGRadS+3
            ENDIF
          WRITE(164,REC=IGRadS+1) NDSETSW
          WRITE(164,REC=IGRadS+2) NP
          WRITE(164,REC=IGRadS+3) DT*NSPOOLGW
          WRITE(164,REC=IGRadS+4) NSPOOLGW
          WRITE(164,REC=IGRadS+5) 2
          IGRadS=IGRadS+5
          CLOSE(164)                    ! DO THIS TO FLUSH THE WRITE BUFFER
        ENDIF
        IF((NOUTGW.EQ.-1).OR.(NOUTGW.EQ.-4)) THEN
          CALL OPEN_GBL_FILE(164, TRIM(GLOBALDIR)//'/'//'rads.64',
     $      NP_G, NP, HEADER74)
          IGRadS=2
        ENDIF
#endif

        IF(NOUTGW.EQ.2) THEN
          open(73,file=TRIM(LOCALDIR)//'/'//'fort.73',
     &          access='DIRECT',recl=nbyte)
c         WRITE(73,REC=ITempStp+1) ndsetsw  ! ALLOW ADDITIONAL OUTPUT DATA TO BE WRITTEN
          close(73)                         ! DO THIS TO FLUSH THE WRITE BUFFER
          OPEN(74,FILE=TRIM(LOCALDIR)//'/'//'fort.74',
     &          ACCESS='DIRECT',RECL=NBYTE)
c         WRITE(74,REC=ITEMPSTP+1) NDSETSW  ! ALLOW ADDITIONAL OUTPUT DATA TO BE WRITTEN
          CLOSE(74)                         ! DO THIS TO FLUSH THE WRITE BUFFER
        ENDIF

!       IF (ncmode.eq.1. or. MYPROC.eq. 0) THEN
!                 PRINT *,"  ----------START--PARAMETERS-after 74------"
!                   PRINT *,"  IMHS = ", IMHS
!                   PRINT *,"  TIME = ", TIME
!                   PRINT *,"  ITHS = ", ITHS
!                   PRINT *,"  ETA1 = ", ETA1(11)
!                   PRINT *,"  ETA2 = ", ETA2(11)
!                   PRINT *,"  ETADD = ", EtaDisc(11)
!                   PRINT *,"  uu2 = ", uu2(11)
!                   PRINT *,"  vv2 = ", vv2(11)
!                   PRINT *,"  nodecode = ", nnodecode(11)
!                   PRINT *,"  noff = ", noff(11)
!                   PRINT *,"  IESTP = ",IESTP
!                   PRINT *,"  NSCOUE = ",NSCOUE
!                   PRINT *,"  IVSTP = ",IVSTP
!                   PRINT *,"  NSCOUV = ",NSCOUV
!                   PRINT *,"  ICSTP = ",ICSTP
!                   PRINT *,"  NSCOUC = ",NSCOUC
!                   PRINT *,"  IPSTP = ",IPSTP
!                   PRINT *,"  IWSTP = ",IWSTP
!                   PRINT *,"  NSCOUM = ",NSCOUM
!                   PRINT *,"  IGEP = ",IGEP
!                   PRINT *,"  NSCOUGE = ",NSCOUGE
!                   PRINT *,"  IGVP = ",IGVP
!                   PRINT *,"  NSCOUGV = ",NSCOUGV
!                   PRINT *,"  IGCP = ",IGCP
!                   PRINT *,"  NSCOUGC = ",NSCOUGC
!                   PRINT *,"  IGPP = ",IGPP
!                   PRINT *,"  IGWP = ",IGWP
!                   PRINT *,"  NSCOUGW = ",NSCOUGW
!                   PRINT *,
!     &             "  ---------END-------PARAMETERS--------------------"
!     
!      ENDIF
C...
C...  HOT START INFORMATION FOR A 3D RUN
C...

C  kmd48.33bc - added time and time step to be passed to the hotstart file
      IF(C3D) THEN
        CALL HOTSTART_3D(TIME,ITHS)
      ENDIF


C...
C......HOT START INFORMATION FOR HARMONIC ANALYSIS
C...
      IF(IHARIND.EQ.1) THEN
        IHABEG=ITHAS+NHAINC
C...
C........IF HARMONIC ANALYSIS HAS NOT BEGUN, COLD START THE HARMONIC ANALYSIS
C...
        IF(ITHS.LT.IHABEG) THEN
          ICHA=0
          CALL HACOLDS(HAFREQ)
          IF(NHASE.EQ.1) CALL HACOLDSES(NSTAE)
          IF(NHASV.EQ.1) CALL HACOLDSVS(NSTAV)
          IF(NHAGE.EQ.1) CALL HACOLDSEG(NP)
          IF(NHAGV.EQ.1) CALL HACOLDSVG(NP)
          IF (CHARMV) THEN
            DO I=1,NP
               ELAV(I)=0.D0
               XVELAV(I)=0.D0
               YVELAV(I)=0.D0
               ELVA(I)=0.D0
               XVELVA(I)=0.D0
               YVELVA(I)=0.D0
            END DO
          ENDIF   !   charmv

        ENDIF

C...
C........IF HARMONIC ANALYSIS HAS ALREADY BEGUN, READ IN HOT START
C........HARMONIC ANALYSIS, MEAN AND SQUARE INFO
C...
!        PRINT *,"START READING HOT START FILE- HARMONIC"

      IF (IHOT.EQ.67.OR.IHOT.EQ.68)THEN
!        PRINT *,"START READING HOT START FILE- BINARY ICHA"

        IF(ITHS.GT.ITHAS) THEN
          READ(IHOT,REC=IHOTSTP) ICHA ; IHOTSTP = IHOTSTP + 1
        ENDIF
        IF(ITHS.GE.IHABEG) THEN
          CALL HAHOTS(NSTAE,NSTAV,NP,NHASE,NHASV,NHAGE,NHAGV,
     &          NSCREEN,ScreenUnit,IHOTSTP,IHOT,MYPROC)
          IF(NHASE.EQ.1) CALL HAHOTSES(NSTAE,IHOTSTP,IHOT)
          IF(NHASV.EQ.1) CALL HAHOTSVS(NSTAV,IHOTSTP,IHOT)
          IF(NHAGE.EQ.1) CALL HAHOTSEG(NP,IHOTSTP,IHOT)
          IF(NHAGV.EQ.1) CALL HAHOTSVG(NP,IHOTSTP,IHOT)
        ENDIF

      ENDIF  

C..Read in Means and Squares
!        PRINT *,"START READING HOT START FILE - MEANS"

      IF (IHOT.EQ.67.OR.IHOT.EQ.68)THEN
!        PRINT *,"START READING HOT START FILE- BINARY ELAV"

        IF(CHARMV) THEN
          IF((FMV.NE.0.).AND.(ITHS.GT.ITMV)) THEN
            READ(IHOT,REC=IHOTSTP) NTSTEPS ; IHOTSTP = IHOTSTP + 1
            IF(NHAGE.EQ.1) THEN
              DO I=1,NP
                READ(IHOT,REC=IHOTSTP) ELAV(I) ; IHOTSTP = IHOTSTP + 1
                READ(IHOT,REC=IHOTSTP) ELVA(I) ; IHOTSTP = IHOTSTP + 1
              ENDDO
            ENDIF
            IF(NHAGV.EQ.1) THEN
              DO I=1,NP
                READ(IHOT,REC=IHOTSTP) XVELAV(I) ; IHOTSTP = IHOTSTP + 1
                READ(IHOT,REC=IHOTSTP) YVELAV(I) ; IHOTSTP = IHOTSTP + 1
                READ(IHOT,REC=IHOTSTP) XVELVA(I) ; IHOTSTP = IHOTSTP + 1
                READ(IHOT,REC=IHOTSTP) YVELVA(I) ; IHOTSTP = IHOTSTP + 1
              ENDDO
            ENDIF
          ENDIF
        ENDIF   !  charmv
      ENDIF    !HARIND

      ENDIF
C
C     jgf48.4635 If the maxele.63 (etc) files are present, we need to
C     load them up, to preserve high water marks across hotstarts. 
C     If they are not present, thats okay too.  
      if ((MNPROC.gt.1).and.(MyProc.eq.0)
     &     .and.(.not.WRITE_LOCAL_FILES)) then
         ALLOCATE(HOT_ETAMAX_G(NP_G))
         HOT_ETAMAX_G(:) = -99999.d0  ! initialize the value
         CALL readMinMax('maxele.63',HOT_ETAMAX_G)
         ALLOCATE(HOT_UMAX_G(NP_G))
         HOT_UMAX_G(:) = 0.d0  ! initialize the value
         CALL readMinMax('maxvel.63',HOT_UMAX_G)
         ALLOCATE(HOT_WVNOUTMAX_G(NP_G))
         HOT_WVNOUTMAX_G(:) = 0.d0  ! initialize the value
         CALL readMinMax('maxwvel.63',HOT_WVNOUTMAX_G)
         ALLOCATE(HOT_PrMin_G(NP_G))
         HOT_PrMin_G(:) = 99999.d0  ! initialize the value
         CALL readMinMax('minpr.63',HOT_PrMin_G)
         ALLOCATE(HOT_RSNMAX_G(NP_G))
         HOT_RSNMAX_G(:) = 0.d0  ! initialize the value
         CALL readMinMax('maxrs.63',HOT_RSNMAX_G)
      endif
      if ( (MNPROC.eq.1).or.(WRITE_LOCAL_FILES) ) then  
         ALLOCATE(HOT_ETAMAX(NP))
         HOT_ETAMAX(:) = -99999.d0    ! initialize the value
         CALL readMinMax('maxele.63',HOT_ETAMAX)
         ALLOCATE(HOT_UMAX(NP))
         HOT_UMAX(:) = 0.d0  ! initialize the value
         CALL readMinMax('maxvel.63',HOT_UMAX)
         ALLOCATE(HOT_WVNOUTMAX(NP))
         HOT_WVNOUTMAX(:) = 0.d0  ! initialize the value
         CALL readMinMax('maxwvel.63',HOT_WVNOUTMAX)
         ALLOCATE(HOT_PrMin(NP))
         HOT_PrMin(:) = 99999.d0  ! initialize the value
         CALL readMinMax('minpr.63',HOT_PrMin)
         ALLOCATE(HOT_RSNMAX(NP))
         HOT_RSNMAX(:) = 0.d0  ! initialize the value
         CALL readMinMax('maxrs.63',HOT_RSNMAX)
      endif
C
      IF(IHOT.EQ.67.OR.IHOT.EQ.68)THEN
!         IF(Myproc.eq.0) PRINT *, "IHOT == ",IHOT
         CLOSE(IHOT)
      ENDIF
C
 1112 FORMAT(/,1X,79('_'))
 1197 FORMAT(/,1X,'THE E29 MET GRID INTERPOLATING FACTORS ARE ',
     &                'BEING COMPUTED ')
 1198     FORMAT(1X,'FINISHED COMPUTING E29 INTERPOLATING FACTORS',/)
 3220 FORMAT(1X,A32,2X,A24,2X,A24)
 3645 FORMAT(1X,I10,1X,I10,1X,E15.7,1X,I5,1X,I5)
C
      RETURN
      END

C******************************************************************************
C   Subroutine to initialize the 3D routines for a hot start including        *
C   reading the 3D portion of the HOTSTART file and finding proper places in  *
C   3D output files.                                                          *
C                                                                             *
C                                                                             *
c******************************************************************************

      SUBROUTINE HOTSTART_3D(TIME,ITHS)

      USE GLOBAL_3DVS

C   kmd48.33bc - variable information needed from the global.F file
      USE GLOBAL, ONLY : RES_BC_FLAG, RBCTIME1, RBCTIME2,
     &                   RBCTIMEINC, SBCTIME1, SBCTIME2,
     &                   SBCTIMEINC, TBCTIME1, TBCTIME2,
     &                   TBCTIMEINC, BCSTATIM, RBCRATIO,
     &                   LNM_BC1, LNM_BC2, TIMEIT, STATIM,
     &                   SBCSTATIM, TBCSTATIM, SBCRATIO,
     &                   TBCRATIO, DTDPHS, LNM_BC,
     &                   BCFLAG_LNM, BCFLAG_TEMP, TTBCTIME1,
     &                   TTBCTIME2, TTBCSTATIM, TTBCTIMEINC,
     &                   q_heat1, q_heat2, HFLUX, TTBCRATIO,
     &                   NBD

      IMPLICIT NONE

      INTEGER :: IRType                  !number of fields in output file
      INTEGER :: NH, N                   !horizontal & vertical loop counters
!      REAL(SZ) :: RealPartOfQ     
      REAL(SZ) :: ImaginaryPartOfQ

C  kmd48.33bc - new variables
      INTEGER :: NumofBCNodes, J, NOD, IT
      INTEGER, INTENT(IN) :: ITHS
      REAL(SZ), INTENT(IN) :: TIME
      INTEGER :: IDenHS
      CHARACTER(80) :: CDUM80
      REAL(SZ),ALLOCATABLE :: RealQ_Tmp(:,:)
      INTEGER :: INDEX
      INTEGER :: NumofNodes, K
      REAL(SZ),ALLOCATABLE :: TMP(:,:)

      REAL(SZ) :: tempp, temp1, temp2

C...
C...  Define format statements used to initialize 3D output files
C...
 499  FORMAT(1X,A32,2X,A24,2X,A24)
C.. RJW bug fix in format 498 (kendra found this)
 498  FORMAT(1X,I10,1X,I10,1X,E15.7,I10,1X,I10,1X,I3)
 497  FORMAT(5X,'UNIT ',I2,' FORMAT WILL BE ASCII')
 496  FORMAT(5X,'UNIT ',I2,' FORMAT WILL BE BINARY')

C  kmd48.33bc - allocate realQ_tmp
      ALLOCATE(RealQ_Tmp(NP,NFEN) )
C...
C...  Read in 3D portion of hotstart file
C...
      READ(IHOT,REC=IHOTSTP) IDenHS
      IHOTSTP=IHOTSTP+1
      READ(IHOT,REC=IHOTSTP) N3DSD           !3D station density output time counter
      IHOTSTP=IHOTSTP+1
      READ(IHOT,REC=IHOTSTP) I3DSDRec
      IHOTSTP=IHOTSTP+1
      READ(IHOT,REC=IHOTSTP) N3DSV           !3D station velocity output time counter
      IHOTSTP=IHOTSTP+1
      READ(IHOT,REC=IHOTSTP) I3DSVRec
      IHOTSTP=IHOTSTP+1
      READ(IHOT,REC=IHOTSTP) N3DST           !3D station turbulence output time counter
      IHOTSTP=IHOTSTP+1
      READ(IHOT,REC=IHOTSTP) I3DSTRec
      IHOTSTP=IHOTSTP+1
      READ(IHOT,REC=IHOTSTP) N3DGD
      IHOTSTP=IHOTSTP+1
      READ(IHOT,REC=IHOTSTP) I3DGDRec
      IHOTSTP=IHOTSTP+1
      READ(IHOT,REC=IHOTSTP) N3DGV
      IHOTSTP=IHOTSTP+1
      READ(IHOT,REC=IHOTSTP) I3DGVRec
      IHOTSTP=IHOTSTP+1
      READ(IHOT,REC=IHOTSTP) N3DGT
      IHOTSTP=IHOTSTP+1
      READ(IHOT,REC=IHOTSTP) I3DGTRec
      IHOTSTP=IHOTSTP+1

C  kmd48.33bc - need to determine if IDEN parameter in hotstart is different than IDEN from
C               fort.15 file. If so, determine which one should be used in the simulation. 
      IF(MYPROC.EQ.0) THEN 
        PRINT *, "IDEN = ", IDEN
        PRINT *, "IDENHS = ", IDENHS
      END IF 
      IF(IDenHS.NE.IDEN) THEN 
        IF ((CBAROCLINIC).AND.(RES_BC_FLAG.GE.1)) THEN 
            ! IDEN value from fort.15 should be used 
        ELSE IF ((CBAROCLINIC).AND.(RES_BC_FLAG.LE.0)) THEN 
            IDEN=IDenHS
        END IF 
      END IF 

C  kmd48.33bc - changed to the format used with the global IO
      DO NH=1,NP
        READ(IHOT,REC=IHOTSTP) DUU(NH)
        IHOTSTP=IHOTSTP+1
      END DO 
      DO NH=1,NP
        READ(IHOT,REC=IHOTSTP) DUV(NH)
        IHOTSTP=IHOTSTP+1
      END DO 
      DO NH=1,NP
        READ(IHOT,REC=IHOTSTP) DVV(NH)
        IHOTSTP=IHOTSTP+1
      END DO 
      DO NH=1,NP
        READ(IHOT,REC=IHOTSTP) UU(NH)
        IHOTSTP=IHOTSTP+1
      END DO 
      DO NH=1,NP
        READ(IHOT,REC=IHOTSTP) VV(NH)
        IHOTSTP=IHOTSTP+1
      END DO 
      DO NH=1,NP
        READ(IHOT,REC=IHOTSTP) BSX(NH)
        IHOTSTP=IHOTSTP+1
         IF((MYPROC.EQ.23).AND.(NH.EQ.15784)) THEN
             PRINT *, "BSX(346561) = ", BSX(15784)
         END IF 
      END DO 
      DO NH=1,NP
        READ(IHOT,REC=IHOTSTP) BSY(NH)
        IHOTSTP=IHOTSTP+1
         IF((MYPROC.EQ.23).AND.(NH.EQ.15784)) THEN
             PRINT *, "BSY(346561) = ", BSY(15784)
         END IF 
      ENDDO

C  kmd48.33bc - changed to format used in global IO, layer output
      DO N=1,NFEN
        DO NH=1,NP
          READ(IHOT,REC=IHOTSTP) RealQ_Tmp(NH,N)
          IHOTSTP=IHOTSTP+1
         IF((MYPROC.EQ.23).AND.(NH.EQ.15785)) THEN
             PRINT *, "RealQ_Tmp(346562,1) = ", RealQ_Tmp(15785,1)
         END IF 
        END DO 
      END DO 
      DO N=1,NFEN
        DO NH=1,NP
          READ(IHOT,REC=IHOTSTP) ImaginaryPartOfQ
          IHOTSTP=IHOTSTP+1
          Q(NH,N) = RealQ_Tmp(NH,N) + I*ImaginaryPartOfQ 
        END DO 
      END DO 
      DO N=1,NFEN
        DO NH=1,NP
          READ(IHOT,REC=IHOTSTP) WZ(NH,N)
          IHOTSTP=IHOTSTP+1
        END DO 
      END DO 
      DO N=1,NFEN
        DO NH=1,NP
          READ(IHOT,REC=IHOTSTP) q20(NH,N)
          IHOTSTP=IHOTSTP+1
        END DO 
      END DO 
      DO N=1,NFEN
        DO NH=1,NP
          READ(IHOT,REC=IHOTSTP) l(NH,N)
          IHOTSTP=IHOTSTP+1
        END DO 
      END DO 

      IF(ABS(IDen).EQ.1) THEN 
        DO N=1,NFEN
          DO NH=1,NP
            READ(IHOT,REC=IHOTSTP) SigT(NH,N)
            IHOTSTP=IHOTSTP+1           
          END DO
        END DO 
      ENDIF
      IF(ABS(IDen).EQ.2) THEN 
        DO N=1,NFEN
          DO NH=1,NP
            READ(IHOT,REC=IHOTSTP) Sal(NH,N)
            IHOTSTP=IHOTSTP+1           
          END DO
        END DO 
      ENDIF
      IF(ABS(IDen).EQ.3) THEN 
        DO N=1,NFEN
          DO NH=1,NP
            READ(IHOT,REC=IHOTSTP) Temp(NH,N)
            IHOTSTP=IHOTSTP+1
          END DO
        END DO 
      ENDIF
      IF(ABS(IDen).EQ.4) THEN 
        DO N=1,NFEN
          DO NH=1,NP
            READ(IHOT,REC=IHOTSTP) Sal(NH,N)
            IHOTSTP=IHOTSTP+1
         IF((MYPROC.EQ.23).AND.(NH.EQ.15785)) THEN
             PRINT *, "Sal(346562,1) = ", Sal(15785,1)
         END IF 
          END DO 
        END DO 
        DO N=1,NFEN
          DO NH=1,NP
            READ(IHOT,REC=IHOTSTP) Temp(NH,N)
            IHOTSTP=IHOTSTP+1
         IF((MYPROC.EQ.23).AND.(NH.EQ.15785)) THEN
             PRINT *, "Temp(346562,1) = ", Temp(15785,1)
         END IF 
          ENDDO
        ENDDO
      END IF 

      DEALLOCATE(RealQ_Tmp)

C  kmd48.33bc - restart boundary condtions for baroclinic simulations
      IF (MYPROC.EQ.0) THEN 
        PRINT *, "RES_BC_FLAG = ", RES_BC_FLAG
        PRINT *, "BCFLAG_LNM = ", BCFLAG_LNM
      END IF 
      IF((RES_BC_FLAG.GT.0).AND.(CBAROCLINIC)) THEN
         IF (ABS(RES_BC_FLAG).GE.1) THEN
            IF (BCFLAG_LNM.EQ.1) THEN
              OPEN(35,FILE=TRIM(INPUTDIR)//'/'//'fort.35')
              RBCTIME1=BCSTATIM*86400.d0
              RBCTIME2=RBCTIME1+RBCTIMEINC
              READ(35,'(A)') CDUM80
              DO NumofBCNodes=1,NETA
                 READ(35,*) NOD,LNM_BC1(NumofBCNodes)
              END DO
              READ(35,'(A)') CDUM80
              DO NumofBCNodes=1,NETA
                 READ(35,*) NOD,LNM_BC2(NumofBCNodes)
              END DO
              DO IT=1,ITHS
                 TIMEIT=IT*DTDPHS + STATIM*86400.D0
                 IF(TIMEIT.GT.RBCTIME2) THEN
                    RBCTIME1=RBCTIME2
                    RBCTIME2=RBCTIME1+RBCTIMEINC
                    READ(35,'(A)') CDUM80
                    DO NumofBCNodes=1,NETA
                       LNM_BC1(NumofBCNodes)=LNM_BC2(NumofBCNodes)
                       READ(35,*) NOD,LNM_BC2(NumofBCNodes)
                    END DO
                 END IF
              END DO
              IF(TIME.GT.RBCTIME2) THEN
                 RBCTIME1=RBCTIME2
                 RBCTIME2=RBCTIME2+RBCTIMEINC
                 READ(35,'(A)') CDUM80
                 DO NumofBCNodes=1,NETA
                   LNM_BC1(NumofBCNodes)=LNM_BC2(NumofBCNodes)
                   READ(35,*) NOD,LNM_BC2(NumofBCNodes)
                 END DO
              END IF
              RBCRATIO=(TIMEIT-RBCTIME1)/RBCTIMEINC
            ELSE IF (BCFLAG_LNM.EQ.2) THEN
              ! to be added later but currently information
              ! is taken in two ways from the HYCOM results
            ELSE IF (BCFLAG_LNM.EQ.3) THEN
              DO NumofBCNodes=1,NETA
                 INDEX=NBD(NumofBCNodes)
                 LNM_BC1(NumofBCNodes)=ETA2(INDEX)
                 LNM_BC2(NumofBCNodes)=LNM_BC1(NumofBCNodes)
              END DO
            END IF
         END IF
         IF ((ABS(RES_BC_FLAG).EQ.2).OR.(ABS(RES_BC_FLAG).EQ.4)) THEN
            OPEN(36,FILE=TRIM(INPUTDIR)//'/'//'fort.36')
            SBCTIME1=SBCSTATIM*86400.D0
            SBCTIME2=SBCTIME1+SBCTIMEINC
            READ(36,'(A)') CDUM80
            DO NumofBCNodes=1,NETA
               READ(36,*) NOD,(RESSAL1(NumofBCNodes,J),J=1,NFEN)
            END DO
            READ(36,'(A)') CDUM80
            DO NumofBCNodes=1,NETA
               READ(36,*) NOD,(RESSAL2(NumofBCNodes,J),J=1,NFEN)
            END DO
            DO IT=1,ITHS
               TIMEIT=IT*DTDPHS + STATIM*86400.D0
               IF(TIMEIT.GT.SBCTIME2) THEN
                  SBCTIME1=SBCTIME2
                  SBCTIME2=SBCTIME1+SBCTIMEINC
                  READ(36,'(A)') CDUM80
                  DO NumofBCNodes=1,NETA
                     DO J=1,NFEN
                       RESSAL1(NumofBCNodes,J)=RESSAL2(NumofBCNodes,J)
                     END DO
                     READ(36,*) NOD,(RESSAL2(NumofBCNodes,J),J=1,NFEN)
                  END DO
               END IF
            END DO
            IF(TIME.GT.SBCTIME2) THEN
               SBCTIME1=SBCTIME2
               SBCTIME2=SBCTIME2+SBCTIMEINC
               READ(36,'(A)') CDUM80
               DO NumofBCNodes=1,NETA
                 DO J=1,NFEN
                   RESSAL1(NumofBCNodes,J)=RESSAL2(NumofBCNodes,J)
                 END DO
                 READ(36,*) NOD,(RESSAL2(NumofBCNodes,J),J=1,NFEN)
               END DO
            END IF
            SBCRATIO=(TIMEIT-SBCTIME1)/SBCTIMEINC
         END IF
         IF ((ABS(RES_BC_FLAG).EQ.3).OR.
     &            (ABS(RES_BC_FLAG).EQ.4)) THEN
            OPEN(37,FILE=TRIM(INPUTDIR)//'/'//'fort.37')
            TBCTIME1=TBCSTATIM*86400.D0
            TBCTIME2=TBCTIME1+TBCTIMEINC
            READ(37,'(A)') CDUM80
            DO NumofBCNodes=1,NETA
               READ(37,*) NOD,(RESTEMP1(NumofBCNodes,J),J=1,NFEN)
            END DO
            READ(37,'(A)') CDUM80
            DO NumofBCNodes=1,NETA
               READ(37,*) NOD,(RESTEMP2(NumofBCNodes,J),J=1,NFEN)
            END DO
            DO IT=1,ITHS
               TIMEIT=IT*DTDPHS + STATIM*86400.D0
               IF(TIMEIT.GT.TBCTIME2) THEN
                  TBCTIME1=TBCTIME2
                  TBCTIME2=TBCTIME1+TBCTIMEINC
                  READ(37,'(A)') CDUM80
                  DO NumofBCNodes=1,NETA
                     DO J=1,NFEN
                       RESTEMP1(NumofBCNodes,J)=RESTEMP2(NumofBCNodes,J)
                     END DO
                     READ(37,*) NOD,(RESTEMP2(NumofBCNodes,J),J=1,NFEN)
                  END DO
               END IF
            END DO
            IF(TIME.GT.TBCTIME2) THEN
               TBCTIME1=TBCTIME2
               TBCTIME2=TBCTIME2+TBCTIMEINC
               READ(37,'(A)') CDUM80
               DO NumofBCNodes=1,NETA
                 DO J=1,NFEN
                   RESTEMP1(NumofBCNodes,J)=RESTEMP2(NumofBCNodes,J)
                 END DO
                 READ(37,*) NOD,(RESTEMP2(NumofBCNodes,J),J=1,NFEN)
               END DO
            END IF
            TBCRATIO=(TIMEIT-TBCTIME1)/TBCTIMEINC
            IF (BCFLAG_TEMP.EQ.1) THEN ! read in file with multiple values but no calcs
               OPEN(38,FILE=TRIM(INPUTDIR)//'/'//'fort.38')
               TTBCTIME1=TTBCSTATIM*86400.D0
               TTBCTIME2=TTBCTIME1+TTBCTIMEINC
               DO NumofNodes=1,NP
                  READ(38,*) NOD, q_heat1(Numofnodes)
               END DO
               DO NumofNodes=1,NP
                  READ(38,*) NOD, q_heat2(NumofNodes)
               END DO
               DO IT=1,ITHS
                  TIMEIT=IT*DTDPHS + STATIM*86400.D0
                  IF(TIMEIT.GT.TTBCTIME2) THEN
                     TTBCTIME1=TTBCTIME2
                     TTBCTIME2=TTBCTIME1+TTBCTIMEINC
                     DO NumofNodes=1,NH
                       q_heat1(NumofNodes)=q_heat2(NumofNodes)
                       READ(38,*) NOD, q_heat2(NumofNodes)
                     END DO
                  END IF
               END DO
               TTBCRATIO=(TIMEIT-TTBCTIME1)/TTBCTIMEINC
               DO NumofNodes=1,NP
                  HFLUX(NumofNodes) = q_heat1(NumofNodes)+TTBCRATIO*
     &                (q_heat2(NumofNodes)-q_heat1(NumofNodes))
                  qsurf(NumofNodes)=-HFLUX(NumofNodes)
                  qsurfkp1(NumofNodes)=-HFLUX(NumofNodes)
               END DO
            ELSE IF (BCFLAG_TEMP.EQ.2) THEN ! read in from file with 6 components
               ALLOCATE (TMP(NP,6))
               OPEN(38,FILE=TRIM(INPUTDIR)//'/'//'fort.38')
               TTBCTIME1=TTBCSTATIM*86400.D0
               TTBCTIME2=TTBCTIME1+TTBCTIMEINC
               READ(38,*) (NOD,(TMP(K,J),J=1,6),K=1,NP)
               DO NumofNodes=1,NP
                  q_heat1(NumofNodes)=-TMP(NumofNodes,1)-
     &               TMP(NumofNodes,2)+TMP(NumofNodes,3)-
     &               TMP(NumofNodes,5)+TMP(NumofNodes,4)-
     &               TMP(NumofNodes,6)
               END DO
               READ(38,*) (NOD,(TMP(K,J),J=1,6),K=1,NP)
               DO NumofNodes=1,NP
                  q_heat2(NumofNodes)=-TMP(NumofNodes,1)-
     &               TMP(NumofNodes,2)+TMP(NumofNodes,3)-
     &               TMP(NumofNodes,5)+TMP(NumofNodes,4)-
     &               TMP(NumofNodes,6)
               END DO
               DO IT=1,ITHS
                  TIMEIT=IT*DTDPHS + STATIM*86400.D0
                  IF(TIMEIT.GT.TTBCTIME2) THEN
                     TTBCTIME1=TTBCTIME2
                     TTBCTIME2=TTBCTIME1+TTBCTIMEINC
                     DO NumofNodes=1,NP
                        q_heat1(NumofNodes)=q_heat2(NumofNodes)
                     END DO
                     READ(38,*) (NOD,(TMP(K,J),J=1,6),K=1,NP)
                     DO NumofNodes=1,NP
                        q_heat2(NumofNodes)=-TMP(NumofNodes,1)-
     &                     TMP(NumofNodes,2)+TMP(NumofNodes,3)-
     &                     TMP(NumofNodes,5)+TMP(NumofNodes,4)-
     &                     TMP(NumofNodes,6)
                     END DO
                  END IF
               END DO
               TTBCRATIO=(TIMEIT-TTBCTIME1)/TTBCTIMEINC
               DO NumofNodes=1,NP
                  HFLUX(NumofNodes) = q_heat1(NumofNodes)+TTBCRATIO*
     &                (q_heat2(NumofNodes)-q_heat1(NumofNodes))
                  qsurf(NumofNodes)=-HFLUX(NumofNodes)
                  qsurfkp1(NumofNodes)=-HFLUX(NumofNodes)
               END DO
               DEALLOCATE(TMP)
            ELSE IF (BCFLAG_TEMP.EQ.3) THEN ! read in from file with 4 components
               ALLOCATE (TMP(NP,4))
               OPEN(38,FILE=TRIM(INPUTDIR)//'/'//'fort.38')
               TTBCTIME1=TTBCSTATIM*86400.D0
               TTBCTIME2=TTBCTIME1+TTBCTIMEINC
               READ(38,*) (NOD,(TMP(K,J),J=1,4),K=1,NP)
               DO NumofNodes=1,NP
                  q_heat1(NumofNodes)=TMP(NumofNodes,4)+
     &               TMP(NumofNodes,3)-TMP(NumofNodes,1)+
     &               TMP(NumofNodes,2)
               END DO
               READ(38,*) (NOD,(TMP(K,J),J=1,4),K=1,NP)
               DO NumofNodes=1,NP
                  q_heat2(NumofNodes)=TMP(NumofNodes,4)+
     &               TMP(NumofNodes,3)-TMP(NumofNodes,1)+
     &               TMP(NumofNodes,2)
               END DO
               DO IT=1,ITHS
                  TIMEIT=IT*DTDPHS + STATIM*86400.D0
                  IF(TIMEIT.GT.TTBCTIME2) THEN
                     TTBCTIME1=TTBCTIME2
                     TTBCTIME2=TTBCTIME1+TTBCTIMEINC
                     DO NumofNodes=1,NP
                        q_heat1(NumofNodes)=q_heat2(NumofNodes)
                     END DO
                     READ(38,*) (NOD,(TMP(K,J),J=1,4),K=1,NP)
                     DO NumofNodes=1,NP
                        q_heat2(NumofNodes)=TMP(NumofNodes,4)+
     &                     TMP(NumofNodes,3)-TMP(NumofNodes,1)+
     &                     TMP(NumofNodes,2)
                     END DO
                  END IF
               END DO
               TTBCRATIO=(TIMEIT-TTBCTIME1)/TTBCTIMEINC
               DO NumofNodes=1,NP
                  HFLUX(NumofNodes) = q_heat1(NumofNodes)+TTBCRATIO*
     &                (q_heat2(NumofNodes)-q_heat1(NumofNodes))
                  qsurf(NumofNodes)=-HFLUX(NumofNodes)
                  qsurfkp1(NumofNodes)=-HFLUX(NumofNodes)
               END DO
               DEALLOCATE(TMP)
            END IF ! end if for heat flux file
         END IF ! end if for temperature boundary and heat flux
      ELSE IF ((RES_BC_FLAG.LT.0).AND.(CBAROCLINIC)) THEN
         IF (BCFLAG_LNM.EQ.1) THEN 
            OPEN(35,FILE=TRIM(INPUTDIR)//'/'//'fort.35')
            READ(35,'(A)') CDUM80
            DO NumofBCNodes=1,NETA
                READ(35,*) NOD,LNM_BC(NumofBCNodes)
            END DO
         ELSE IF (BCFLAG_LNM.EQ.2) THEN 
            ! add later
         ELSE IF (BCFLAG_LNM.EQ.3) THEN 
            DO NumofBCNodes=1,NETA
               INDEX=NBD(NumofBCNodes)
               LNM_BC(NumofBCNodes)=ETA2(INDEX)
         !IF((MYPROC.EQ.23).AND.(INDEX.EQ.15785)) THEN
         !    PRINT *, "ETA2(346562) = ", ETA2(15785)
         !END IF 
            END DO
         END IF
      END IF 

C.... Initialize the station 3D density output file
!kmd48.33bc - added in the information in front of all 
!             the other output files for resetting the 
!             information. 
      IF (I3DSD.LT.0) THEN 
        N3DSD=0
        I3DSDRec=0 
        IF ((NTO3DSDS.LT.ITHS).AND.(NSpo3DSD.GT.0)) THEN 
           NTO3DSDS=NTO3DSDS+((REAL(ITHS)-REAL(NTO3DSDS))/
     &              REAL(NSpo3DSD))*REAL(NSpo3DSD)
           IF (NTO3DSDS.LT.ITHS) NTO3DSDS=NTO3DSDS+NSpo3DSD
           IF (DTDP.NE.DTDPHS) THEN
              NDSet3DSD=((NTO3DSDF*DTDP)-(NTO3DSDS*DTDPHS))/
     &                  (NSpo3DSD*DTDP)
           ELSE 
              NDSet3DSD=(NTO3DSDF-NTO3DSDS)/NSpo3DSD
           END IF
        END IF 
      END IF 
            

      IF(IDen.EQ.1) IRType=1
      IF((IDen.EQ.2).OR.(IDen.EQ.3)) IRType=2
      IF(IDen.EQ.4) IRType=3

      IF(I3DSD.EQ.-1) THEN             !start a new ASCII file
        WRITE(16,497) 41
        OPEN(41,FILE=TRIM(LOCALDIR)//'/'//'fort.41')
        WRITE(41,499) RUNDes,RunID,AGrid
        WRITE(41,498) NDSet3DSD,NSta3DD,DTDP*NSpo3DSD,NSpo3DSD,NFEN,
     &                                                            IRType
        CLOSE(41)
        I3DSDRec=2
      ENDIF

      IF(I3DSD.EQ.-2) THEN             !start a new BINARY file
        WRITE(16,496) 41
        OPEN(41,FILE=TRIM(LOCALDIR)//'/'//'fort.41',ACCESS='DIRECT',
     &       RECL=NByte)
        I3DSDRec=0
        IF(NByte.EQ.4) THEN
          DO N=1,8
            I3DSDRec=I3DSDRec+1
            WRITE(41,REC=I3DSDRec) RDES4(N)
          ENDDO
          DO N=1,6
            I3DSDREC=I3DSDREC+1
            WRITE(41,REC=I3DSDRec) RID4(N)
          ENDDO
          DO N=1,6
            I3DSDRec=I3DSDRec+1
            WRITE(41,REC=I3DSDRec) AID4(N)
          ENDDO
        ENDIF
        IF(NByte.EQ.8) THEN
          DO N=1,4
            I3DSDRec=I3DSDRec+1
            WRITE(41,REC=I3DSDRec) RDES8(N)
          ENDDO
          DO N=1,3
            I3DSDRec=I3DSDRec+1
            WRITE(41,REC=I3DSDRec) RID8(N)
          ENDDO
          DO N=1,3
            I3DSDRec=I3DSDRec+1
            WRITE(41,REC=I3DSDRec) AID8(N)
          ENDDO
        ENDIF
        WRITE(41,REC=I3DSDRec+1) NDSet3DSD
        WRITE(41,REC=I3DSDRec+2) NSta3DD
        WRITE(41,REC=I3DSDRec+3) DTDP*NSpo3DSD
        WRITE(41,REC=I3DSDRec+4) NSpo3DSD
        WRITE(41,REC=I3DSDRec+5) NFEN   
        WRITE(41,REC=I3DSDRec+6) IRType
        I3DSDRec=I3DSDRec+6
        CLOSE(41)
      ENDIF

      IF(I3DSD.EQ.1) THEN          !find spot in existing ASCII file
        WRITE(16,497) 41
        OPEN(41,FILE=TRIM(INPUTDIR)//'/'//'fort.41')
        DO N=1,I3DSDRec
          READ(41,*) 
        ENDDO
        ENDFILE(41)
        CLOSE(41)
      ENDIF

      IF(I3DSD.EQ.2) THEN          !find spot in existing BINARY file
        WRITE(16,496) 41
        OPEN(41,FILE=TRIM(LOCALDIR)//'/'//'fort.41',ACCESS='DIRECT',
     &       RECL=NByte)
        CLOSE(41)
      ENDIF

C.... Initialize the 3D velocity station output file (Unit 42)
!kmd48.33bc - added in the information in front of all 
!             the other output files for resetting the 
!             information. 
      IF (I3DSV.LT.0) THEN 
        N3DSV=0
        I3DSVRec=0 
        IF ((NTO3DSVS.LT.ITHS).AND.(NSpo3DSV.GT.0)) THEN 
           NTO3DSVS=NTO3DSVS+((REAL(ITHS)-REAL(NTO3DSVS))/
     &              REAL(NSpo3DSV))*REAL(NSpo3DSV)
           IF (NTO3DSVS.LT.ITHS) NTO3DSVS=NTO3DSVS+NSpo3DSV
           IF (DTDP.NE.DTDPHS) THEN
              NDSet3DSV=((NTO3DSVF*DTDP)-(NTO3DSVS*DTDPHS))/
     &                  (NSpo3DSV*DTDP)
           ELSE 
              NDSet3DSV=(NTO3DSVF-NTO3DSVS)/NSpo3DSV
           END IF
        END IF 
      END IF 

      IRType=3

      IF(I3DSV.EQ.-1) THEN             !start a new ASCII file
        WRITE(16,497) 42
        OPEN(42,FILE=TRIM(INPUTDIR)//'/'//'fort.42')
        WRITE(42,499) RUNDES,RUNID,AGRID
        WRITE(42,498) NDSet3DSV,NSta3DV,DTDP*NSpo3DSV,NSpo3DSV,NFEN,
     &                                                            IRType
        CLOSE(42)
        I3DSVRec=2
      ENDIF

      IF(ABS(I3DSV).EQ.-2) THEN        !start a new BINARY file
        WRITE(16,496) 42
        OPEN(42,FILE=TRIM(INPUTDIR)//'/'//'fort.42',ACCESS='DIRECT',
     &       RECL=NBYTE)
        I3DSVRec=0
        IF(NByte.EQ.4) THEN
          DO N=1,8
            I3DSVRec=I3DSVRec+1
            WRITE(42,REC=I3DSVRec) RDES4(N)
          ENDDO
          DO N=1,6
            I3DSVRec=I3DSVRec+1
            WRITE(42,REC=I3DSVRec) RID4(N)
          ENDDO
          DO N=1,6
            I3DSVRec=I3DSVRec+1
            WRITE(42,REC=I3DSVRec) AID4(N)
          ENDDO
        ENDIF
        IF(NByte.EQ.8) THEN
          DO N=1,4
            I3DSVRec=I3DSVRec+1
            WRITE(42,REC=I3DSVRec) RDES8(N)
          ENDDO
          DO N=1,3
            I3DSVRec=I3DSVRec+1
            WRITE(42,REC=I3DSVRec) RID8(N)
          ENDDO
          DO N=1,3
            I3DSVRec=I3DSVRec+1
            WRITE(42,REC=I3DSVRec) AID8(N)
          ENDDO
        ENDIF
        WRITE(42,REC=I3DSVRec+1) NDSet3DSV
        WRITE(42,REC=I3DSVRec+2) NSta3DV
        WRITE(42,REC=I3DSVRec+3) DTDP*NSpo3DSV
        WRITE(42,REC=I3DSVRec+4) NSpo3DSV
        WRITE(42,REC=I3DSVRec+5) NFEN
        WRITE(42,REC=I3DSVRec+6) IRType
        I3DSVRec=I3DSVRec+6
        CLOSE(42)
      ENDIF

      IF(I3DSV.EQ.1) THEN              !find spot in existing ASCII file
        WRITE(16,497) 42
        OPEN(42,FILE=TRIM(INPUTDIR)//'/'//'fort.42')
        DO N=1,I3DSVRec
          READ(42,*) 
        ENDDO
        ENDFILE(42)
        CLOSE(42)
      ENDIF

      IF(I3DSV.EQ.2) THEN              !find spot in existing BINARY file
        WRITE(16,496) 42
        OPEN(42,FILE=TRIM(INPUTDIR)//'/'//'fort.42',ACCESS='DIRECT',
     &       RECL=NByte)
        CLOSE(42)
      ENDIF

C.... Initialize the 3D turbulence station output file (Unit 43)
!kmd48.33bc - added in the information in front of all 
!             the other output files for resetting the 
!             information. 
      IF (I3DST.LT.0) THEN 
        N3DST=0
        I3DSTRec=0 
        IF ((NTO3DSTS.LT.ITHS).AND.(NSpo3DST.GT.0)) THEN 
           NTO3DSTS=NTO3DSTS+((REAL(ITHS)-REAL(NTO3DSTS))/
     &              REAL(NSpo3DST))*REAL(NSpo3DST)
           IF (NTO3DSTS.LT.ITHS) NTO3DSTS=NTO3DSTS+NSpo3DST
           IF (DTDP.NE.DTDPHS) THEN
              NDSet3DST=((NTO3DSTF*DTDP)-(NTO3DSTS*DTDPHS))/
     &                  (NSpo3DST*DTDP)
           ELSE 
              NDSet3DST=(NTO3DSTF-NTO3DSTS)/NSpo3DST
           END IF
        END IF 
      END IF 

      IRType=3

C  kmd48.33bc - changed NDSet3DST to use DTDPHS
      IF (DTDP.NE.DTDPHS) THEN
         NDSet3DST=((NTO3DSTF*DTDP)-(NTO3DSTS*DTDPHS))/
     &              (NSpo3DST*DTDP)
      END IF

      IF(I3DST.EQ.-1) THEN             !start a new ASCII file
        WRITE(16,497) 43
        OPEN(43,FILE=TRIM(LOCALDIR)//'/'//'fort.43')
        WRITE(43,499) RUNDES,RUNID,AGRID
        WRITE(43,498) NDSet3DST,NSta3DT,DTDP*NSpo3DST,NSpo3DST,NFEN,
     &                                                            IRType
        CLOSE(43)
        I3DSTRec=2
      ENDIF

      IF(ABS(I3DST).EQ.2) THEN         !start a new BINARY file
        WRITE(16,496) 43
        OPEN(43,FILE=TRIM(LOCALDIR)//'/'//'fort.43',ACCESS='DIRECT',
     &       RECL=NByte)
        I3DSTRec=0
        IF(NByte.EQ.4) THEN
          DO N=1,8
            I3DSTRec=I3DSTRec+1
            WRITE(43,REC=I3DSTRec) RDES4(N)
          ENDDO
          DO N=1,6
            I3DSTRec=I3DSTRec+1
            WRITE(43,REC=I3DSTRec) RID4(N)
          ENDDO
          DO N=1,6
            I3DSTRec=I3DSTRec+1
            WRITE(43,REC=I3DSTRec) AID4(N)
          ENDDO
        ENDIF
        IF(NByte.EQ.8) THEN
          DO N=1,4
            I3DSTRec=I3DSTRec+1
            WRITE(43,REC=I3DSTRec) RDES8(N)
          ENDDO
          DO N=1,3
            I3DSTRec=I3DSTRec+1
            WRITE(43,REC=I3DSTRec) RID8(N)
          ENDDO
          DO N=1,3
            I3DSTRec=I3DSTRec+1
            WRITE(43,REC=I3DSTRec) AID8(N)
          ENDDO
        ENDIF
        WRITE(43,REC=I3DSTRec+1) NDSet3DST
        WRITE(43,REC=I3DSTRec+2) NSta3DT
        WRITE(43,REC=I3DSTRec+3) DTDP*NSpo3DST
        WRITE(43,REC=I3DSTRec+4) NSpo3DST
        WRITE(43,REC=I3DSTRec+5) NFEN
        WRITE(43,REC=I3DSTRec+6) IRType
        I3DSTRec=I3DSTRec+6
        CLOSE(43)
      ENDIF

      IF(I3DST.EQ.1) THEN          !find spot in existing ASCII file
        WRITE(16,497) 43
        OPEN(43,FILE=TRIM(INPUTDIR)//'/'//'fort.43')
        DO N=1,I3DSTRec
          READ(43,*) 
        ENDDO
        ENDFILE(43)
        CLOSE(43)
      ENDIF

      IF(I3DST.EQ.2) THEN          !find spot in existing BINARY file
        WRITE(16,496) 43
        OPEN(43,FILE=TRIM(INPUTDIR)//'/'//'fort.43',ACCESS='DIRECT',
     &       RECL=NByte)
        CLOSE(43)
      ENDIF

C.... Initialize the global 3D density output file (Unit 44)

C  kmd48.33bc - change NDSet3DGD to use DTDPHS
      IF(I3DGD.LT.0) THEN
         N3DGD=0
         I3DGDRec=0
         IF((NTO3DGDS.LT.ITHS).AND.(NSpo3DGD.GT.0)) THEN
           NTO3DGDS=NTO3DGDS+((REAL(ITHS)-REAL(NTO3DGDS))/
     &               REAL(NSpo3DGD))*REAL(NSpo3DGD)
           IF(NTO3DGDS.LT.ITHS) NTO3DGDS=NTO3DGDS+NSpo3DGD
           IF (DTDP.NE.DTDPHS) THEN
             NDSet3DGD=((NTO3DGDF*DTDP)-(NTO3DGDS*DTDPHS))/
     &                  (NSpo3DGD*DTDP)
           ELSE
             NDSet3DGD=(NTO3DGDF-NTO3DGDS)/NSpo3DGD
           ENDIF
         END IF
      END IF

      IF(IDen.EQ.1) IRType=1
      IF((IDen.EQ.2).OR.(IDen.EQ.3)) IRType=2
      IF(IDen.EQ.4) IRType=3
C     jgf46.27 Replaced IRType with IDen
!kmd48.33bc - added in the information to make a file with all the
!             top temperature boundary condition information in
      IF(I3DGD.EQ.-1) THEN         !start a new ASCII file
        WRITE(16,497) 44
        OPEN(44,FILE=TRIM(LOCALDIR)//'/'//'fort.44')
        WRITE(44,499) RUNDES,RUNID,AGRID
        WRITE(44,498) NDSet3DGD,NP,DTDP*NSpo3DGD,NSpo3DGD,NFEN,IDen
        CLOSE(44)
        IF((IDEN.EQ.3).OR.(IDEN.EQ.4)) THEN
          IF (BCFLAG_TEMP.NE.0) THEN
            OPEN(47,FILE=TRIM(LOCALDIR)//'/'//'fort.47')
            WRITE(47,499) RUNDES,RUNID,AGRID
            WRITE(47,498) NDSet3DGD,NP,DTDP*NSpo3DGD,NSpo3DGD,1
            CLOSE(47)
          END IF
        END IF
        I3DGDRec=2
      ENDIF

      IF(I3DGD.EQ.-2) THEN         !start a new BINARY file
        WRITE(16,496) 44
        OPEN(44,FILE=TRIM(LOCALDIR)//'/'//'fort.44',ACCESS='DIRECT',
     &       RECL=NByte)
        I3DGDRec=0
        IF(NByte.EQ.4) THEN
          DO N=1,8
            I3DGDRec=I3DGDRec+1
            WRITE(44,REC=I3DGDRec) RDES4(N)
          ENDDO
          DO N=1,6
            I3DGDRec=I3DGDRec+1
            WRITE(44,REC=I3DGDRec) RID4(N)
          ENDDO
          DO N=1,6
            I3DGDRec=I3DGDRec+1
            WRITE(44,REC=I3DGDRec) AID4(N)
          ENDDO
        ENDIF
        IF(NByte.EQ.8) THEN
          DO N=1,4
            I3DGDRec=I3DGDRec+1
            WRITE(44,REC=I3DGDRec) RDES8(N)
          ENDDO
          DO N=1,3
            I3DGDRec=I3DGDRec+1
            WRITE(44,REC=I3DGDRec) RID8(N)
          ENDDO
          DO N=1,3
            I3DGDRec=I3DGDRec+1
            WRITE(44,REC=I3DGDRec) AID8(N)
          ENDDO
        ENDIF
        WRITE(44,REC=I3DGDRec+1) NDSet3DGD
        WRITE(44,REC=I3DGDRec+2) NP
        WRITE(44,REC=I3DGDRec+3) DTDP*NSpo3DGD
        WRITE(44,REC=I3DGDRec+4) NSpo3DGD
        WRITE(44,REC=I3DGDRec+5) NFEN
        WRITE(44,REC=I3DGDRec+6) IDEN
        I3DGDRec=I3DGDRec+6
        CLOSE(44)
      ENDIF

      IF(I3DGD.EQ.1) THEN          !find spot in existing ASCII file
        WRITE(16,497) 44
        OPEN(44,FILE=TRIM(INPUTDIR)//'/'//'fort.44')
        DO N=1,I3DGDRec
          READ(44,*) 
        ENDDO
        ENDFILE(44)
        CLOSE(44)
      ENDIF

      IF(I3DGD.EQ.2) THEN          !find spot in existing BINARY file
        WRITE(16,496) 44
        OPEN(44,FILE=TRIM(INPUTDIR)//'/'//'fort.44',ACCESS='DIRECT',
     &       RECL=NByte)
        CLOSE(44)
      ENDIF
      
c.... Initialize the global 3D velocity output file (Unit 45)

C  kmd48.33bc - changed NDSet3DGV to use DTDPHS and added REALs to 
C              calculation of the start time due to timestep change
      IF(I3DGV.LT.0) THEN
        N3DGV=0
        I3DGVRec=0
        IF((NTO3DGVS.LT.ITHS).AND.(NSpo3DGV.GT.0)) THEN
          NTO3DGVS=NTO3DGVS+((REAL(ITHS)-REAL(NTO3DGVS))/
     &             REAL(NSpo3DGV))*REAL(NSpo3DGV)
          IF(NTO3DGVS.LT.ITHS) NTO3DGVS=NTO3DGVS+NSpo3DGV
          IF (DTDP.NE.DTDPHS) THEN
            NDSet3DGV=((NTO3DGVF*DTDP)-(NTO3DGVS*DTDPHS))/
     &                 (NSpo3DGV*DTDP)
          ELSE
            NDSet3DGV=(NTO3DGVF-NTO3DGVS)/NSpo3DGV
          END IF
        ENDIF
      ENDIF

      IF(I3DGV.EQ.-1) THEN         !start a new ASCII file
        WRITE(16,497) 45
        OPEN(45,FILE=TRIM(LOCALDIR)//'/'//'fort.45')
        WRITE(45,499) RUNDES,RUNID,AGRID
        WRITE(45,498) NDSet3DGV,NP,DTDP*NSpo3DGV,NSpo3DGV,NFEN,IRType
        CLOSE(45)
        I3DGVRec=2
      ENDIF

      IF(I3DGV.EQ.-2) THEN        !start a new BINARY file
        WRITE(16,496) 45
        OPEN(45,FILE=TRIM(LOCALDIR)//'/'//'fort.45',ACCESS='DIRECT',
     &       RECL=NBYTE)
        I3DGVRec=0
        IF(NByte.EQ.4) THEN
          DO N=1,8
            I3DGVRec=I3DGVRec+1
            WRITE(45,REC=I3DGVRec) RDES4(N)
          ENDDO
          DO N=1,6
            I3DGVRec=I3DGVRec+1
            WRITE(45,REC=I3DGVRec) RID4(N)
          ENDDO
          DO N=1,6
            I3DGVRec=I3DGVRec+1
            WRITE(45,REC=I3DGVRec) AID4(N)
          ENDDO
        ENDIF
        IF(NByte.EQ.8) THEN
          DO N=1,4
            I3DGVRec=I3DGVRec+1
            WRITE(45,REC=I3DGVRec) RDES8(N)
          ENDDO
          DO N=1,3
            I3DGVRec=I3DGVRec+1
            WRITE(45,REC=I3DGVRec) RID8(N)
          ENDDO
          DO N=1,3
            I3DGVRec=I3DGVRec+1
            WRITE(45,REC=I3DGVRec) AID8(N)
          ENDDO
        ENDIF
        WRITE(45,REC=I3DGVRec+1) NDSet3DGV
        WRITE(45,REC=I3DGVRec+2) NP
        WRITE(45,REC=I3DGVRec+3) DTDP*NSpo3DGV
        WRITE(45,REC=I3DGVRec+4) NSpo3DGV
        WRITE(45,REC=I3DGVRec+5) NFEN
        WRITE(45,REC=I3DGVRec+6) IRType
        I3DGVRec=I3DGVRec+6
        CLOSE(45)
      ENDIF

      IF(I3DGV.EQ.1) THEN          !find spot in existing ASCII file
        WRITE(16,497) 45
        OPEN(45,FILE=TRIM(INPUTDIR)//'/'//'fort.45')
        DO N=1,I3DGVRec
          READ(45,*) 
        ENDDO
        ENDFILE(45)
        CLOSE(45)
      ENDIF

      IF(I3DGV.EQ.2) THEN          !find spot in existing BINARY file
        WRITE(16,496) 45
        OPEN(45,FILE=TRIM(INPUTDIR)//'/'//'fort.45',ACCESS='DIRECT',
     &       RECL=NByte)
        CLOSE(45)
      ENDIF

c.... Initialize the global 3D turbulence output file (Unit 46)

C  kmd48.33bc - changed NDSet3DGT to use DTDPHS
      IF(I3DGT.LT.0) THEN
        N3DGT=0
        I3DGTRec=0
        IF((NTO3DGTS.LT.ITHS).AND.(NSpo3DGT.GT.0)) THEN
          NTO3DGTS=NTO3DGTS+((REAL(ITHS)-REAL(NTO3DGTS))/
     &             REAL(NSpo3DGT))*REAL(NSpo3DGT)
          IF(NTO3DGTS.LT.ITHS) NTO3DGTS=NTO3DGTS+NSpo3DGT
          IF (DTDP.NE.DTDPHS) THEN
            NDSet3DGT=((NTO3DGTF*DTDP)-(NTO3DGTS*DTDPHS))/
     &                 (NSpo3DGT*DTDP)
          ELSE
            NDSet3DGT=(NTO3DGTF-NTO3DGTS)/NSpo3DGT
          END IF
        ENDIF
      ENDIF

      IF(I3DGT.EQ.-1) THEN         !start a new ASCII file
        WRITE(16,497) 46
        OPEN(46,FILE=TRIM(LOCALDIR)//'/'//'fort.46')
        WRITE(46,499) RUNDES,RUNID,AGRID
        WRITE(46,498) NDSet3DGT,NP,DTDP*NSpo3DGT,NSpo3DGT,NFEN,IRType
        CLOSE(46)
        I3DGTRec=2
      ENDIF

      IF(I3DGT.EQ.-2) THEN         !start a new BINARY file
        WRITE(16,496) 46
        OPEN(46,FILE=TRIM(LOCALDIR)//'/'//'fort.46',ACCESS='DIRECT',
     &       RECL=NByte)
        I3DGTRec=0
        IF(NByte.EQ.4) THEN
          DO N=1,8
            I3DGTRec=I3DGTRec+1
            WRITE(46,REC=I3DGTRec) RDES4(N)
          ENDDO
          DO N=1,6
            I3DGTRec=I3DGTRec+1
            WRITE(46,REC=I3DGTRec) RID4(N)
          ENDDO
          DO N=1,6
            I3DGTRec=I3DGTRec+1
            WRITE(46,REC=I3DGTRec) AID4(N)
          ENDDO
        ENDIF
        IF(NByte.EQ.8) THEN
          DO N=1,4
            I3DGTRec=I3DGTRec+1
            WRITE(46,REC=I3DGTRec) RDES8(N)
          ENDDO
          DO N=1,3
            I3DGTRec=I3DGTRec+1
            WRITE(46,REC=I3DGTRec) RID8(N)
          ENDDO
          DO N=1,3
            I3DGTRec=I3DGTRec+1
            WRITE(46,REC=I3DGTRec) AID8(N)
          ENDDO
        ENDIF
        WRITE(46,REC=I3DGTRec+1) NDSet3DGT
        WRITE(46,REC=I3DGTRec+2) NP
        WRITE(46,REC=I3DGTRec+3) DTDP*NSpo3DGT
        WRITE(46,REC=I3DGTRec+4) NSpo3DGT
        WRITE(46,REC=I3DGTRec+5) NFEN
        WRITE(46,REC=I3DGTRec+6) IRType
        I3DGTRec=I3DGTRec+6
        CLOSE(46)
      ENDIF

      IF(I3DGT.EQ.1) THEN          !find spot in existing ASCII file
        WRITE(16,497) 46
        OPEN(46,FILE=TRIM(INPUTDIR)//'/'//'fort.46')
        DO N=1,I3DGTRec
          READ(46,*) 
        ENDDO
        ENDFILE(46)
        CLOSE(46)
      ENDIF

      IF(I3DGT.EQ.2) THEN          !find spot in existing BINARY file
        WRITE(16,496) 46
        OPEN(46,FILE=TRIM(INPUTDIR)//'/'//'fort.46',ACCESS='DIRECT',
     &       RECL=NByte)
        CLOSE(46)
      ENDIF

C.... Set up a few final odds and ends for a 3D run

      CALL VSSTUP ()
        
      RETURN
      END SUBROUTINE HOTSTART_3D
