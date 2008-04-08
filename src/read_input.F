C******************************************************************************
C PADCIRC VERSION 46.00 xx/xx/2006                                            *
C  last changes in this file VERSION 46.00                                    *
C                                                                             *
C This module handles most of the model input.  The primary 2d input is read  *
C in subroutine READ_INPUT.  The primary 3D input is read in subroutine       *
C READ_INPUT_3D.  Initial conditions that are read in for a cold start are    *
C handled in the cold start subroutines.                                      *
C                                                                             *
C******************************************************************************

C-----------------------------------------------------------------------
C     S U B R O U T I N E     R E A D _ I N P U T
C-----------------------------------------------------------------------
C
C     READS INPUT FILES
C
C-----------------------------------------------------------------------
      SUBROUTINE READ_INPUT()
      USE SIZES               !since GLOBAL uses SIZES, is this necessary?
      USE GLOBAL
      USE GLOBAL_IO
      USE GLOBAL_3DVS, ONLY: C3D_PTrans
      USE HARM
      USE WIND
      USE ITPACKV
      USE VERSION
!     jgf46.00
      USE NodalAttributes, ONLY : 
     &     NoLiBF, NWP, Tau0, HBreak, FTheta, FGamma, Tau, CF,
     &     InitNAModule, ReadNodalAttr, InitNodalAttr, ESLM, ESLC,
     &     Tau0FullDomainMin, Tau0FullDomainMax
#ifdef CMPI
      USE MESSENGER, ONLY : msg_fini
#endif
      IMPLICIT NONE
      INTEGER I, J, K, N, ICK, IPRBI                 !local loop counters
      INTEGER NIBP,IBN1,IK,NDISC,NBBN,NVEL2
      INTEGER IDUM80
      INTEGER IMDig1,IMDig2,IMDig3,IMDig4,IMDig5,IMDig6
      REAL(SZ) RampVal
      REAL(SZ) ZNGFicNodeDist,ZNGFicNodeDistTemp
      REAL(SZ) TOUTSC
      REAL(SZ) TOUTSE,TOUTFE
      REAL(SZ) TOUTSV,TOUTFV
      REAL(SZ) TOUTFGC, TOUTFGE, TOUTFGV, TOUTFGW
      REAL(SZ) TOUTFM, TOUTSGC, TOUTSGE, TOUTSGV, TOUTSGW, TOUTSM
      CHARACTER CDUM80

C      MCF ADDITION FOR READING AND WRITING STATION INFO IN NETCDF
C
      REAL, EXTERNAL              :: a2f
      CHARACTER(132) STATLINE
      CHARACTER(50) LVAR(3)
C
C     jgf46.00 Zero out all the variables in the Nodal Attributes
C     Module. 
      CALL InitNAModule() 
C...
C...OPEN STATEMENT FOR UNIT 14 AND 15 INPUT FILES
C...  
      OPEN(14,FILE=TRIM(INPUTDIR)//'/'//'fort.14')
      OPEN(15,FILE=TRIM(INPUTDIR)//'/'//'fort.15')
C...
C...OPEN STATEMENT FOR UNIT 16 OUTPUT FILE
C...
      OPEN(16,FILE=TRIM(LOCALDIR)//'/'//'fort.16')
C
C...GENERAL PURPOSE FORMAT STATEMENTS
C...  
 1112 FORMAT(/,1X,79('_'))
 9972 FORMAT(////,1X,'!!!!!!!!!! INPUT ERROR !!!!!!!!!',/)
 9973 FORMAT(/,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
 9974 FORMAT(/,1X,'!!!!!! EXECUTION WILL CONTINUE !!!!!!!!',//)

C...
C...PRINT OUT HEADER FOR OUTPUT INCLUDING VERSION NUMBER AND COPYRIGHT
C...
      WRITE(16,1112)
      WRITE(16,1112)
      WRITE(16,1114) TRIM(ADC_VERSION)
      WRITE(16,1112)

 1114 FORMAT(//,19X,'PROGRAM ADCIRC   VERSION ',A,
     &     //,5X,'AN ADVANCED CIRCULATION MODEL FOR SHELVES, COASTAL ',
     &     'SEAS AND ESTUARIES',
     &     ///,7X,'-  DEVELOPED BY',
     &     //,10X,'R.A. LUETTICH, JR',
     &     /,12X,'UNIVERSITY OF NORTH CAROLINA AT CHAPEL HILL',
     &     /,12X,'INSTITUTE OF MARINE SCIENCES',
     &     //,10X,'J.J. WESTERINK ',
     &     /,12X,
     &     'DEPARTMENT OF CIVIL ENGINEERING AND GEOLOGICAL SCIENCES',
     &     /,12X,'UNIVERSITY OF NOTRE DAME',
     &     ///,7X,'-  THE ADCIRC SOURCE CODE IS COPYRIGHTED BY',
     &     //,10X,'R.A. LUETTICH, JR. AND J.J. WESTERINK, 1994-2006',
     &     //,7X,
     &     'NO PART OF THIS CODE MAY BE REPRODUCED OR REDISTRIBUTED',
     &     /,10X,'WITHOUT THE WRITTEN PERMISSION OF THE AUTHORS',//)
      
C...  
C...  WRITE OUT HEADER INFORMATION DESCRIBING HOW THE CODE HAS BE SET UP
C...  
      WRITE(16,1210)
 1210 FORMAT(//,1X,'THE ADCIRC SOURCE CODE HAS BEEN CONFIGURED ',
     &             'BY THE PREPROCESSOR AS FOLLOWS:',/)
  
#ifdef CMACHSUN
      WRITE(16,*) '      - CODE SETUP TO RUN ON SUN 4 OR SPARC ',
     &            'COMPUTERS'
#endif

#ifdef REAL4  
      WRITE(16,*) '      - CODE SETUP TO RUN WITH 4 byte REALS'
      WRITE(16,4312)
4312  FORMAT(/,'**** W A R N I N G ****',/
            ' You are running ADCIRC in SINGLE PRECISION! ',/
            ' It is always recommended that you run ADCIRC in ',/
            ' DOUBLE PRECISION ONLY.'
#else
      WRITE(16,*) '      - CODE SETUP TO RUN WITH 8 byte REALS'
#endif

#ifdef CVEC
      WRITE(16,*) '      - CODE OPTIMIZED FOR A VECTOR COMPUTER'
#endif

#ifdef CSCA
      WRITE(16,*) '      - CODE OPTIMIZED FOR A SCALAR COMPUTER'
#endif

      WRITE(16,*) '      - CODE WILL USE JCG ITERATIVE GWCE SOLVER'
      WRITE(16,1112)

C...  
C...  INPUT FROM UNIT 15 AND OUTPUT TO UNIT 16 RUN DESCRIPTION AND RUN
C...  IDENTIFICATION
C...  
      READ(15,'(A80)') RUNDES
      READ(15,'(A80)') RUNID
      WRITE(16,1) RUNDES
  1   FORMAT(//,1X,'RUN DESCRIPTION : ',A80)
      WRITE(16,209) RUNID
209   FORMAT(/,1X,'RUN IDENTIFICATION : ',A80)

	do I=1,20
         J=(I-1)*4+1
         RDES4(I)=RUNDES(J:J+3)
         RID4(I) =RUNID (J:J+3)
      end do
	do I=1,10
         J=(I-1)*8+1
         RDES8(I)=RUNDES(J:J+7)
         RID8(I) =RUNID (J:J+7)
      end do
C
C     jgf46.00 Relocated the reading of the number of nodes and the
C     number of elements from the grid file to this position in the read
C     input subroutine. It is needed in the processing of the nodal
C     attributes file.
      READ(14,'(A80)') AGRID

      do I=1,20
         J=(I-1)*4+1
         AID4(I)=AGRID(J:J+3)
      end do
      do I=1,10
         J=(I-1)*8+1
         AID8(I)=AGRID(J:J+7)
      end do

      READ(14,*) NE,NP
      MNP  = NP
      MNE  = NE
C...
C... READ AND PROCESS NFOVER - NONFATAL ERROR OVERRIDE OPTION
C...
C
C     jgf46.10 Add user-controllable warning, output, and stop criteria
C     for elevations. Initialize default values.
      WarnElev = 20.0         ! default 
      iWarnElevDump = 0       ! init
      WarnElevDump = .False.  ! default
      WarnElevDumpLimit = 50  ! default
      WarnElevDumpCounter = 0 ! init
      ErrorElev = 1000.0      ! default
C
#ifndef DEBUG_WARN_ELEV
      READ(15,*) NFOVER
#else
      READ(15,*) NFOVER, WarnElev, iWarnElevDump, WarnElevDumpLimit, 
     &  ErrorElev
#endif
C
      WRITE(16,1112)
      WRITE(16,1250)
 1250 FORMAT(//,1X,'GENERAL RUN INFORMATION',/)
      IF(NFOVER.EQ.1) THEN
         WRITE(16,1951) NFOVER
 1951    FORMAT(5X,'NFOVER = ',I2,
     &        /,9X,'IF NON-FATAL INPUT ERRORS ARE DETECTED, THEY WILL ',
     &        'BE CORRECTED AND EXECUTION CONTINUED')
      ELSE
         WRITE(16,1952) NFOVER
 1952    FORMAT(/,5X,'NFOVER = ',I3,
     &        /,9X,'NON-FATAL INPUT ERRORS WILL STOP EXECUTION ',/)
      ENDIF
#ifdef DEBUG_WARN_ELEV
      IF (iWarnElevDump .ne. 0) WarnElevDump = .True. 
      WRITE(16,1953) WarnElev,WarnElevDump,WarnElevDumpLimit,ErrorElev
 1953 FORMAT(//,5X,
     &     'A warning will be issued if elevation exceeds WarnElev = ',
     &     e16.8,
     &     /,5X,'A global elevation file (fort.69) will be written if '
     &     /,5X,'WarnElev is exceeded and WarnElevDump is true: ',L2,
     &     /,5X,'Execution will be terminated if ',
     &     '(WarnElevDumpLimit = 'I3,') ', 
     &     /,5X,'global elevation files have been written as warning.'
     &     /,5X,'Execution will be terminated if elevation exceeds' 
     &     ' ErrorElev =',e16.8)
#endif
C...
C...  READ AND PROCESS NABOUT - ABBREVIATED UNIT 16 OUTPUT OPTION
C...
      READ(15,*) NABOUT
      IF(NABOUT.EQ.1) THEN
         WRITE(16,3501) NABOUT
 3501    FORMAT(5X,'NABOUT = ',I2,
     &        /,9X,'ABREVIATED OUTPUT WILL BE PROVIDED TO UNIT 16',/,9X,
     &        'UNIT 14, 21, 22 INPUT DATA WILL NOT BE ECHO PRINTED',/)
      ELSE
         WRITE(16,3502) NABOUT
 3502    FORMAT(/,5X,'NABOUT = ',I3,
     &        /,9X,'DETAILED OUTPUT WILL BE PROVIDED TO UNIT 16',/,9X,
     &        'UNIT 14, 15, 21, 22 INPUT DATA WILL BE ECHO PRINTED',/)
      ENDIF

C...
C...  READ AND PROCESS NSCREEN - SCREEN OUTPUT OPTION
C...
C
C     jgf46.00 Added option to output data to the screen every NSCREEN
C     time steps, rather than on every time step.

C     jgf46.19 Added option to output "screen" data to fort.999 file
C     rather than the screen. This can be superior to a shell redirect
C     because on some platforms, the redirected log file is not
C     available until the run is complete.
C
      READ(15,*) NSCREEN
      ScreenUnit=0
      IF(NSCREEN.GT.0.AND.MYPROC.EQ.0) THEN
         ScreenUnit=6
         WRITE(16,3561) NSCREEN
 3561    FORMAT(5X,'NSCREEN = ',I6,
     &        /,9X,'SCREEN OUTPUT WILL BE PROVIDED TO UNIT 6',
     &        /,9X,'EVERY NSCREEN TIME STEPS.',/)
      ELSEIF(NSCREEN.LT.0.AND.MYPROC.EQ.0) THEN
         ScreenUnit=999
         WRITE(16,3562) NSCREEN
 3562    FORMAT(/,5X,'NSCREEN = ',I6,
     &        /,9X,'SCREEN OUTPUT WILL BE PROVIDED TO adcirc.log',/)
         OPEN(ScreenUnit,FILE=TRIM(GLOBALDIR)//'/'//'adcirc.log')
      ELSE
         WRITE(16,3563) NSCREEN
 3563    FORMAT(/,5X,'NSCREEN = ',I6,
     &        /,9X,'SCREEN OUTPUT WILL NOT BE PROVIDED',/)
      ENDIF

      IF (MYPROC.EQ.0.and.NScreen.ne.0) THEN
         WRITE(ScreenUnit,1112)
         WRITE(ScreenUnit,1114) TRIM(ADC_VERSION)
         WRITE(ScreenUnit,1112)
      ENDIF
C...
C...  READ AND PROCESS IHOT - HOT START OPTION
C...  
      READ(15,*) IHOT
      IF (IHOT < 0) THEN
        hss % lun      = 67
        hss % filename = 'fort.67'
      ELSE IF ((IHOT == 0).OR.(IHOT == 68)) then
        hss % lun      = 67
        hss % filename = 'fort.67'
      ELSE IF (IHOT == 67) then
        hss % lun      = 68
        hss % filename = 'fort.68'
      ENDIF

      IF(IHOT.NE.0) THEN
         WRITE(16,9733) IHOT
 9733    FORMAT(/,5X,'ADCIRC will be hot started using information ',
     &        'on UNIT ',I2)
      ELSE
         WRITE(16,9734)
 9734    FORMAT(/,5X,'ADCIRC will be cold started')
      ENDIF
C...
C...  READ AND PROCESS ICS - CARTESIAN/SPHERICAL COORDINATE OPTION
C...
      READ(15,*) ICS
      IF((ICS.NE.1).AND.(ICS.NE.2)) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
            WRITE(ScreenUnit,9972)
            WRITE(ScreenUnit,*) 'ICS =',ICS
            WRITE(ScreenUnit,9735)
            WRITE(ScreenUnit,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'ICS =',ICS
         WRITE(16,9735)
         WRITE(16,9973)
 9735    FORMAT(/,1X,'Your selection of ICS (a UNIT 15 input ',
     &        'parameter) is not an allowable value')
         STOP
      ENDIF
      IF(ICS.EQ.1) THEN
         WRITE(16,9736) ICS
 9736    FORMAT(/,5X,'ICS = ',I2,
     &        /,9X,'Governing equations are in Cartesian coordinates')
      ELSE
         WRITE(16,9737) ICS
 9737    FORMAT(/,5X,'ICS = ',I2,
     &        /,9X,'Governing equations are in Spherical coordinates',
     &        /,9X,'mapped using a CPP projection')
      ENDIF

C...  
C...  READ AND PROCESS IM - 2D/3D MODEL FORMULATION OPTION
C...  
      READ(15,*) IM
      WRITE(16,*) ' '
      WRITE(16,*) '    IM = ',IM
      WRITE(16,*) ' '
c     - - - - - - - - - - - - - - - - - - - - - - - -
      IF (IM.LT.100) THEN
c     jgf Set defaults for model type (IM). All LOGICAL variables are
c     initialized to .FALSE. when declared in global.F.  To use original
c     GWCE equation, uncomment the following line and comment out the
c     following three lines. 
c         CGWCE_Orig       = .TRUE.  !uncomment for original GWCE
         CGWCE_New        = .TRUE. !comment out for original GWCE
         CGWCE_LS_KGQ     = .TRUE. !comment out for original GWCE
         CGWCE_Advec_NC   = .TRUE. !comment out for original GWCE
c     jgf To use original momentum equations, uncomment the following
c     line and comment out the following two lines.  
c     CME_Orig         = .TRUE. !uncomment for original momentum eqs  
         CME_New_NC       = .TRUE. !comment out for original momentum eqs  
         CME_LS_IBPV      = .TRUE. !comment out for original momentum eqs  
         CME_AreaInt_Corr = .TRUE.  
      ENDIF
c     - - - - - - - - - - - - - - - - - - - - - - - -
      IF (IM.EQ.0) THEN
         C2DDI         = .TRUE.
         ILump=0
         WRITE(16,*) '    ADCIRC is configured for a Barotropic 2DDI ',
     &        'run using: New GWCE and Momentum Eq formulations'
      ELSEIF (IM.EQ.10) THEN
         C2DDI         = .TRUE.
         C2D_PTrans    = .TRUE.
         ILump=0
         WRITE(16,*) '    ADCIRC is configured for a Barotropic 2DDI ',
     &        'run using: New GWCE and Momentum Eq '
         WRITE(16,*) '          formulations + Passive Scalar Transport'
      ELSEIF (IM.EQ.20) THEN
         C2DDI         = .TRUE.
         CBaroclinic   = .TRUE.
         ILump=0
         WRITE(16,*) '    ADCIRC is configured for a Baroclinic 2DDI ',
     &        'run using: New GWCE and Momentum Eq formulations'
      ELSEIF (IM.EQ.30) THEN
         C2DDI         = .TRUE.
         C2D_PTrans    = .TRUE.
         CBaroclinic   = .TRUE.
         ILump=0
         WRITE(16,*) '    ADCIRC is configured for a Baroclinic 2DDI ',
     &        'run using: New GWCE and Momentum Eq '
         WRITE(16,*) '          formulations + Passive Scalar Transport'
      ELSEIF (IM.EQ.1) THEN
         C3D           = .TRUE.
         C3DVS         = .TRUE.
         ILump=0
         WRITE(16,*) '    ADCIRC is configured for a Barotropic 3D ',
     &        'run using: New GWCE and velocity based ',
     &        'Momentum Eqs.'
      ELSEIF (IM.EQ.11) THEN
         C3D           = .TRUE.
         C3DVS         = .TRUE.
         C3D_PTrans    = .TRUE.
         ILump=0
         WRITE(16,*) '    ADCIRC is configured for a Barotropic 3D ',
     &        'run using: New GWCE and velocity based '
         WRITE(16,*) '          Momentum Eqs + Passive Scalar Transport'
      ELSEIF (IM.EQ.21) THEN
         C3D           = .TRUE.
         C3DVS         = .TRUE.
         CBaroclinic   = .TRUE.
         ILump=0
         WRITE(16,*) '    ADCIRC is configured for a Baroclinic 3D ',
     &        'run using: New GWCE and velocity based ',
     &        'Momentum Eqs.'      
      ELSEIF (IM.EQ.31) THEN
         C3D           = .TRUE.
         C3DVS         = .TRUE.
         C3D_PTrans    = .TRUE.
         CBaroclinic   = .TRUE.
         ILump=0
         WRITE(16,*) '    ADCIRC is configured for a Baroclinic 3D ',
     &        'run using: New GWCE and velocity based '
         WRITE(16,*) '          Momentum Eqs + Passive Scalar Transport'
      ELSEIF (IM.EQ.2) THEN
         C3D           = .TRUE.
         C3DDSS        = .TRUE.
         ILump=0
         WRITE(16,*) '    ADCIRC is configured for a 3D run using',
     &        ': New GWCE and stress based Momentum Eqs.'
c     - - - - - - - - - - - - - - - - - - - - - - - -
c     f i n e   g r a i n e d   o p t i o n s ( i m )
c     - - - - - - - - - - - - - - - - - - - - - - - -
      ELSEIF ((IM.GE.111111).AND.(IM.LE.534322)) THEN
         IMDig1 = IM/100000
         IMDig2 = (IM - 100000*IMDig1)/10000
         IMDig3 = (IM - 100000*IMDig1 - 10000*IMDig2)/1000
         IMDig4 = (IM - 100000*IMDig1 - 10000*IMDig2 - 1000*IMDig3)/100
         IMDig5 = (IM - 100000*IMDig1 - 10000*IMDig2 - 1000*IMDig3
     &        -  100*IMDig4)/10
         IMDig6 =  IM - 100000*IMDig1 - 10000*IMDig2 - 1000*IMDig3
     &        - 100*IMDig4 -   10*IMDig5

         C2DDI     = .TRUE.
         CGWCE_New = .TRUE.                       
         WRITE(16,*) '    ADCIRC is configured for a 2DDI run using'
         WRITE(16,*) '    the new GWCE routine and:'
         IF(IMDig1.EQ.1) THEN
            CGWCE_LS_KGQ     = .TRUE. !jgf This is the default.
            WRITE(16,*) '        Kolar-Gray, flux based lateral ',
     &           'stress in GWCE'
         ELSEIF(IMDig1.EQ.2) THEN
            CGWCE_LS_2PartQ  = .TRUE.
            WRITE(16,*) '        2 Part, flux based lateral ',
     &           'stress in GWCE'
         ELSEIF(IMDig1.EQ.3) THEN
            CGWCE_LS_2PartV  = .TRUE.
            WRITE(16,*) '        2 Part, velocity based lateral ',
     &           'stress in GWCE'
         ELSEIF(IMDig1.EQ.4) THEN
            CGWCE_LS_2PartSQ  = .TRUE.
            WRITE(16,*) '        2 Part, flux based lateral ',
     &           'symmetric stress in GWCE'
         ELSEIF(IMDig1.EQ.5) THEN
            CGWCE_LS_2PartSV  = .TRUE.
            WRITE(16,*) '        2 Part, velocity based lateral ',
     &           'symmetric stress in GWCE'
         ENDIF
c     - - - - - - - - - - - - - - - - - - - - - - - -
         IF(IMDig2.EQ.1) THEN
            CGWCE_Advec_NC   = .TRUE. !jgf This is the default.
            WRITE(16,*) '        Non conservative advection in GWCE'
         ELSEIF(IMDig2.EQ.2) THEN
            CGWCE_Advec_C1   = .TRUE.
            WRITE(16,*) '        Conservative form 1 advection in GWCE'
         ELSEIF(IMDig2.EQ.3) THEN
            CGWCE_Advec_C2   = .TRUE.
            WRITE(16,*) '        Conservative form 2 advection in GWCE'
         ENDIF
c     - - - - - - - - - - - - - - - - - - - - - - - -
         IF(IMDig3.EQ.1) THEN
            CME_LS_IBPV      = .TRUE. !jgf This is the default.
            WRITE(16,*) '        Integration by parts, velocity based ',
     &           'lateral stress in Momentum Eqs.'
         ELSEIF(IMDig3.EQ.2) THEN
            CME_LS_IBPQ      = .TRUE. 
            WRITE(16,*) '        Integration by parts, flux based ',
     &           'lateral stress in Momentum Eqs.'
         ELSEIF(IMDig3.EQ.3) THEN
            CME_LS_IBPSV      = .TRUE.
            WRITE(16,*) '        Integration by parts, velocity based ',
     &           'symmetric lateral stress in Momentum Eqs.'
         ELSEIF(IMDig3.EQ.4) THEN
            CME_LS_IBPSQ      = .TRUE. 
            WRITE(16,*) '        Integration by parts, flux based ',
     &           'symmetric lateral stress in Momentum Eqs.'
         ELSEIF(IMDig3.EQ.5) THEN
            CME_LS_2PartV    = .TRUE.
            WRITE(16,*) '        2 Part, velocity based lateral ',
     &           'stress in Momentum Eqs.'
         ELSEIF(IMDig3.EQ.6) THEN
            CME_LS_2PartQ    = .TRUE.
            WRITE(16,*) '        2 Part, flux based lateral ',
     &           'stress in Momentum Eqs.'
         ENDIF
c     - - - - - - - - - - - - - - - - - - - - - - - -
         IF(IMDig4.EQ.1) THEN
            CME_New_NC     = .TRUE. !jgf This is the default.
            WRITE(16,*) '        Non conservative advection in ',
     &           'Momentum Eqs.'
         ELSEIF(IMDig4.EQ.2) THEN
            CME_New_C1       = .TRUE.
            WRITE(16,*) '        Conservative form 1 advection in ',
     &           'Momentum Eqs.'
         ELSEIF(IMDig4.EQ.3) THEN
            CME_New_C2     = .TRUE.
            WRITE(16,*) '        Conservative form 2 advection in ',
     &           'Momentum Eqs.'
         ENDIF
c     - - - - - - - - - - - - - - - - - - - - - - - -
         IF(IMDig5.EQ.1) THEN
            CME_AreaInt_Corr = .TRUE. !jgf This is the default.
            WRITE(16,*) '        Original Area Integration in ',
     &           'Momentum Eqs.'
         ELSEIF(IMDig5.EQ.2) THEN
            CME_AreaInt_Orig = .TRUE. 
            WRITE(16,*) '        Corrected Area Integration in ',
     &           'Momentum Eqs.'
         ENDIF
c     - - - - - - - - - - - - - - - - - - - - - - - -
         IF(IMDig6.EQ.2) THEN
            ILump=1
            CGWCE_Lump = .TRUE.
            WRITE(16,*) '        Lumped GWCE mass matrix'       
         ELSE
            ILump=0
            WRITE(16,*) '        Consistent GWCE mass matrix'       
         ENDIF
c     - - - - - - - - - - - - - - - - - - - - - - - -
      ELSE
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
            WRITE(ScreenUnit,9972)
            WRITE(ScreenUnit,*) 'IM =',IM
            WRITE(ScreenUnit,9721)
            WRITE(ScreenUnit,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'IM =',IM
         WRITE(16,9721)
         WRITE(16,9973)
 9721    FORMAT(/,1X,'Your selection of IM (a UNIT 15 input ',
     &        'parameter) is not an allowable value')
         STOP
      ENDIF


      IDEN=0
      IF (CBaroclinic) READ(15,*) IDEN
      IF (IDEN.EQ.0) THEN
      ELSEIF (IDEN.EQ.1) THEN
         C2D_BTrans    = .TRUE.
         WRITE(16,*) '    2DDI Prognostic Baroclinic ADCIRC run with ',
     &               'SigmaT forcing'
      ELSEIF (IDEN.EQ.-1) THEN
         WRITE(16,*) '    2DDI Diagnostic Baroclinic ADCIRC run with ',
     &               'SigmaT forcing'
      ELSEIF (IDEN.EQ.2) THEN
         C2D_BTrans    = .TRUE.
         WRITE(16,*) '    2DDI Prognostic Baroclinic ADCIRC run with ',
     &               'Salinity forcing'
      ELSEIF (IDEN.EQ.-2) THEN
         WRITE(16,*) '    2DDI Diagnostic Baroclinic ADCIRC run with ',
     &               'Salinity forcing'
      ELSEIF (IDEN.EQ.3) THEN
         C2D_BTrans    = .TRUE.
         WRITE(16,*) '    2DDI Prognostic Baroclinic ADCIRC run with ',
     &               'Temperature forcing'
      ELSEIF (IDEN.EQ.-3) THEN
         WRITE(16,*) '    2DDI Diagnostic Baroclinic ADCIRC run with ',
     &               'Temperature forcing'
      ELSEIF (IDEN.EQ.4) THEN
         C2D_BTrans    = .TRUE.
         WRITE(16,*) '    2DDI Prognostic Baroclinic ADCIRC run with ',
     &               'Salinity and Temperature forcing'
      ELSEIF (IDEN.EQ.-4) THEN
         WRITE(16,*) '    2DDI Diagnostic Baroclinic ADCIRC run with ',
     &               'Salinity and Temperature forcing'
      ENDIF


      WRITE(16,*) ' '

      WRITE(16,*) '     The ADCIRC logical variables are set to:'
      WRITE(16,*) '         C2DDI            = ',C2DDI
      WRITE(16,*) '         C3D              = ',C3D
      WRITE(16,*) '         C3DDSS           = ',C3DDSS
      WRITE(16,*) '         C3DVS            = ',C3DVS
      WRITE(16,*) '         C2D_BTrans       = ',C2D_BTrans
      WRITE(16,*) '         C2D_PTrans       = ',C2D_PTrans
c     WRITE(16,*) '         C3D_BTrans       = ',C3D_BTrans            !haven't yet read 3D input
      WRITE(16,*) '         C3D_PTrans       = ',C3D_PTrans
      WRITE(16,*) '         CBaroclinic      = ',CBaroclinic
      WRITE(16,*) '         CGWCE_Orig       = ',CGWCE_Orig
      WRITE(16,*) '         CGWCE_Lump       = ',CGWCE_Lump
      WRITE(16,*) '         CGWCE_LS_KGQ     = ',CGWCE_LS_KGQ
      WRITE(16,*) '         CGWCE_LS_2PartQ  = ',CGWCE_LS_2PartQ
      WRITE(16,*) '         CGWCE_LS_2PartV  = ',CGWCE_LS_2PartV
      WRITE(16,*) '         CGWCE_LS_2PartSQ = ',CGWCE_LS_2PartSQ
      WRITE(16,*) '         CGWCE_LS_2PartSV = ',CGWCE_LS_2PartSV
      WRITE(16,*) '         CGWCE_Advec_NC   = ',CGWCE_Advec_NC
      WRITE(16,*) '         CGWCE_Advec_C1   = ',CGWCE_Advec_C1
      WRITE(16,*) '         CGWCE_Advec_C2   = ',CGWCE_Advec_C2
      WRITE(16,*) '         CME_Orig         = ',CME_Orig
      WRITE(16,*) '         CME_New_NC       = ',CME_New_NC
      WRITE(16,*) '         CME_New_C1       = ',CME_New_C1
      WRITE(16,*) '         CME_New_C2       = ',CME_New_C2
      WRITE(16,*) '         CME_LS_IBPQ      = ',CME_LS_IBPQ
      WRITE(16,*) '         CME_LS_IBPV      = ',CME_LS_IBPV
      WRITE(16,*) '         CME_LS_IBPSQ     = ',CME_LS_IBPSQ
      WRITE(16,*) '         CME_LS_IBPSV     = ',CME_LS_IBPSV
      WRITE(16,*) '         CME_LS_2PartQ    = ',CME_LS_2PartQ
      WRITE(16,*) '         CME_LS_2PartV    = ',CME_LS_2PartV
      WRITE(16,*) '         CME_AreaInt_Orig = ',CME_AreaInt_Orig
      WRITE(16,*) '         CME_AreaInt_Corr = ',CME_AreaInt_Corr
      WRITE(16,*) '         CTIP             = ',CTIP
      WRITE(16,*) '         CHARMV           = ',CHARMV

C...
C...  READ AND PROCESS NOLIBF - NONLINEAR BOTTOM FRICTION OPTION
C...
      READ(15,*) NOLIBF
      IF((NOLIBF.LT.0).OR.(NOLIBF.GT.2)) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
            WRITE(ScreenUnit,9972)
            WRITE(ScreenUnit,*) 'NOLIBF =',NOLIBF
            WRITE(ScreenUnit,9722)
            WRITE(ScreenUnit,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'NOLIBF =',NOLIBF
         WRITE(16,9722)
         WRITE(16,9973)
 9722    FORMAT(/,1X,'Your selection of NOLIBF (a UNIT 15 input ',
     &        'parameter) is not an allowable value')
         STOP
      ENDIF
      WRITE(16,9845) NOLIBF
 9845 FORMAT(/,5X,'NOLIBF = ',I3)
      IF (NOLIBF.EQ.0) WRITE(16,2050)
 2050 FORMAT(9X,'THE MODEL WILL USE LINEAR BOTTOM FRICTION')
      IF (NOLIBF.EQ.1) WRITE(16,2051)
 2051 FORMAT(9X,'THE MODEL WILL USE NONLINEAR BOTTOM FRICTION')
      IF (NOLIBF.EQ.2) WRITE(16,2052)
 2052 FORMAT(9X,'THE MODEL WILL USE STANDARD QUADRATIC BOTTOM FRICTION',
     &     'IN DEEP WATER ',
     &     /,9X,'AND A FRICTION FACTOR THAT INCREASES AS THE DEPTH ',
     &     'DECREASES IN SHALLOW WATER')
C...
C... READ AND PROCESS NOLIFA - NONLINEAR FINITE AMPLITUDE OPTION
C...
      READ(15,*) NOLIFA
      IF((NOLIFA.LT.0).OR.(NOLIFA.GT.2)) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
            WRITE(ScreenUnit,9972)       ! input error
            WRITE(ScreenUnit,*) 'NOLIFA =',NOLIFA
            WRITE(ScreenUnit,9723)       ! not allowable
         ENDIF
         WRITE(16,9972)         ! input error
         WRITE(16,*) 'NOLIFA =',NOLIFA
         WRITE(16,9723)         ! not allowable 
 9723    FORMAT(/,1X,'Your selection of NOLIFA (a UNIT 15 input ',
     &        'parameter) is not an allowable value')
         IF (NoLiFA.eq.3.and.NFOver.eq.1) THEN
            WRITE(16,8735)
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,8735)
 8735       FORMAT(/,1X,'WARNING: The StartDry file was replaced ',
     &           'by surface_submergence_state in ',
     &           /,1X,'the Nodal Attributes file (unit 13).'
     &           //,1X,'ACTION: NOLIFA will be corrected to 2; the ',
     &           'loading of StartDry data will not ',
     &           /,1X,'be triggered now, although ',
     &           'it may be triggered later in the NWP section.',/)
            NoLiFA = 2
         ELSE
            IF (NoLiFA.eq.3.and.NFOver.eq.0) THEN
               WRITE(16,7624)   ! startdry replaced
               IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,7624)
            ENDIF
 7624       FORMAT(/,1X,'ERROR: NOLIFA=3 formerly triggered the ',
     &           'loading of StartDry data.'
     &           /,1X,'However, the StartDry file was replaced ',
     &           'by surface_submergence_state in ',
     &           /,1X,'the Nodal Attributes file (unit 13). Please ',
     &           'use NWP to load this data.',/)
            WRITE(16,9973)      ! execution will terminate
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9973)
            STOP
         ENDIF
      ENDIF
      WRITE(16,9846) NOLIFA
 9846 FORMAT(/,5X,'NOLIFA = ',I3)
      IF(NOLIFA.EQ.0) WRITE(16,2053)
 2053 FORMAT(9X,'THE MODEL WILL NOT USE FINITE AMPLITUDE TERMS OR ',
     &     'WETTING AND DRYING')
      IF(NOLIFA.EQ.1) WRITE(16,2054)
 2054 FORMAT(9X,'THE MODEL WILL USE FINITE AMPLITUDE TERMS BUT NO ',
     &     'WETTING AND DRYING')
      IF(NOLIFA.EQ.2) WRITE(16,2049)
 2049 FORMAT(9X,'THE MODEL WILL USE FINITE AMPLITUDE TERMS AND ',
     &     'WETTING AND DRYING')
C...
C...  READ AND PROCESS NOLICA - ADVECTIVE TERM SPATIAL GRADIENT
C...
      READ(15,*) NOLICA
      IF((NOLICA.LT.0).OR.(NOLICA.GT.1)) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
            WRITE(ScreenUnit,9972)
            WRITE(ScreenUnit,*) 'NOLICA =',NOLICA
            WRITE(ScreenUnit,9724)
            WRITE(ScreenUnit,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'NOLICA =',NOLICA
         WRITE(16,9724)
         WRITE(16,9973)
 9724    FORMAT(/,1X,'Your selection of NOLICA (a UNIT 15 input ',
     &        'parameter) is not an allowable value')
         STOP
      ENDIF
      WRITE(16,9847) NOLICA
 9847 FORMAT(/,5X,'NOLICA = ',I3)
      IF(NOLICA.EQ.0) WRITE(16,2055)
 2055 FORMAT(9X,'THE MODEL WILL NOT USE SPATIAL DERIVATIVE ',
     &     'COMPONENTS OF THE ADVECTIVE TERMS')
      IF(NOLICA.EQ.1) WRITE(16,2056)
 2056 FORMAT(9X,'THE MODEL WILL USE SPATIAL DERIVATIVE ',
     &     'COMPONENTS OF THE ADVECTIVE TERMS')

C...
C...  READ AND PROCESS NOLICAT - GWCE ADVECTIVE TERM TIME DERIVATIVE
C...
      READ(15,*) NOLICAT
      IF((NOLICAT.LT.0).OR.(NOLICAT.GT.1)) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
            WRITE(ScreenUnit,9972)
            WRITE(ScreenUnit,*) 'NOLICAT =',NOLICAT
            WRITE(ScreenUnit,9725)
            WRITE(ScreenUnit,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'NOLICAT =',NOLICAT
         WRITE(16,9725)
         WRITE(16,9973)
 9725    FORMAT(/,1X,'Your selection of NOLICAT (a UNIT 15 input ',
     &        'parameter) is not an allowable value')
         STOP
      ENDIF

      IF((NOLIFA.GE.1).AND.(NOLICAT.EQ.0)) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
            WRITE(ScreenUnit,9972)
            WRITE(ScreenUnit,*) 'NOLICAT =',NOLICAT
            WRITE(ScreenUnit,9726)
            IF(NFOVER.EQ.1) THEN
               WRITE(ScreenUnit,9974)
            ELSE
               WRITE(ScreenUnit,9973)
            ENDIF
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'NOLICAT =',NOLICAT
         WRITE(16,9726)
         WRITE(16,9974)
 9726    FORMAT(/,1X,'Your selection of NOLICAT (a UNIT 15 input ',
     &        'parameter) is inconsistent with your ',
     &        /,1X,'selection of NOLIFA and may lead to mass ',
     &        'balance problems')
         IF(NFOVER.EQ.1) THEN
            if (myproc == 0) WRITE(ScreenUnit,9974)
         ELSE
            if (myproc == 0) WRITE(ScreenUnit,9973)
#ifdef CMPI
            call msg_fini()
#endif
            STOP
         ENDIF
      ENDIF

      IF((NOLIFA.EQ.0).AND.(NOLICAT.EQ.1)) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
            WRITE(ScreenUnit,9972)
            WRITE(ScreenUnit,*) 'NOLICAT =',NOLICAT
            WRITE(ScreenUnit,9726)
            IF(NFOVER.EQ.1) THEN
               WRITE(ScreenUnit,9974)
            ELSE
               WRITE(ScreenUnit,9973)
            ENDIF
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'NOLICAT =',NOLICAT
         WRITE(16,9726)
         WRITE(16,9974)
         IF(NFOVER.EQ.1) THEN
            if (myproc == 0) WRITE(ScreenUnit,9974)
         ELSE
            if (myproc == 0) WRITE(ScreenUnit,9973)
#ifdef CMPI
            call msg_fini()
#endif
            STOP
         ENDIF
      ENDIF

      IF(NOLICA.NE.NOLICAT) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
            WRITE(ScreenUnit,9972)
            WRITE(ScreenUnit,*) 'NOLICAT =',NOLICAT
            WRITE(ScreenUnit,9727)
            IF(NFOVER.EQ.1) THEN
               WRITE(ScreenUnit,9974)
            ELSE
               WRITE(ScreenUnit,9973)
            ENDIF
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'NOLICAT =',NOLICAT
         WRITE(16,9727)
         WRITE(16,9974)
 9727    FORMAT(/,1X,'Your selection of NOLICAT (a UNIT 15 input ',
     &        'parameter) is inconsistent with your ',
     &        /,1X,'selection of NOLICA and may lead to mass ',
     &        'balance problems')
         IF(NFOVER.EQ.1) THEN
            if (myproc == 0) WRITE(ScreenUnit,9974)
         ELSE
            if (myproc == 0) WRITE(ScreenUnit,9973)
#ifdef CMPI
            call msg_fini()
#endif
            STOP
         ENDIF
      ENDIF

      WRITE(16,9848) NOLICAT
 9848 FORMAT(/,5X,'NOLICAT = ',I3)
      IF(NOLICAT.EQ.0) WRITE(16,2057)
 2057 FORMAT(9X,'THE MODEL WILL NOT USE TIME DERIVATIVE COMPONENTS ',
     &     /,9X,'OF THE ADVECTIVE TERMS IN THE GWCE')
      IF(NOLICAT.EQ.1) WRITE(16,2058)
 2058 FORMAT(9X,'THE MODEL WILL USE TIME DERIVATIVE COMPONENTS ',
     &     /,9X,'OF THE ADVECTIVE TERMS IN THE GWCE')
C
C     READ AND PROCESS NWP jgf46.00 Read in nodal attributes such as
C     tau0, bottom friction, directional wind speed reduction factor,
C     startdry, etc. The full initialization and error checking of these
C     data must wait until the grid has been read in from unit 14.
      READ(15,*) NWP
      CALL ReadNodalAttr(NScreen, ScreenUnit, MyProc, NAbOut)
C...
C...  READ AND PROCESS NCOR - SPATIALLY VARYING CORIOLIS PARAMETER
C...
      READ(15,*) NCOR
      IF((NCOR.NE.0).AND.(NCOR.NE.1)) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
            WRITE(ScreenUnit,9972)
            WRITE(ScreenUnit,*) 'NCOR =',NCOR
            WRITE(ScreenUnit,9729)
            WRITE(ScreenUnit,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'NCOR =',NCOR
         WRITE(16,9729)
         WRITE(16,9973)
 9729    FORMAT(/,1X,'Your selection of NCOR (a UNIT 15 input ',
     &        'parameter) is not an allowable value')
         STOP
      ENDIF
      IF((ICS.EQ.1).AND.(NCOR.EQ.1)) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
            WRITE(ScreenUnit,9972)
            WRITE(ScreenUnit,*) 'NCOR =',NCOR
            WRITE(ScreenUnit,9730)
            WRITE(ScreenUnit,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'NCOR =',NCOR
         WRITE(16,9730)
         WRITE(16,9973)
 9730    FORMAT(/,1X,'Your selection of NCOR (a UNIT 15 input ',
     &        'parameter) is inconsistent with your ',
     &        /,1X,'selection of coordinate systems.  Spatially ',
     &        'variable Coriolis should be used only with ',
     &        /,1X,'Spherical coordinates')
#ifdef CMPI
         call msg_fini()
#endif
         STOP
      ENDIF
      IF(NCOR.EQ.0) THEN
         WRITE(16,233) NCOR
 233     FORMAT(/,5X,'NCOR = ',I2,
     &    /,9X,'A CONSTANT VALUE OF THE CORIOLIS PARAMETER WILL BE ',
     &    /,9X,'USED THROUGHOUT THE DOMAIN')
      ELSE
         WRITE(16,234) NCOR
 234     FORMAT(/,5X,'NCOR = ',I2,
     &    /,9X,'SPATIALLY VARYING CORIOLIS VALUES WILL BE COMPUTED ',
     &    'FROM INPUT LATITUDES')
      ENDIF

C...
C...  READ AND PROCESS NTIP - TIDAL POTENTIAL FORCING
C...
      READ(15,*) NTIP
      IF((NTIP.LT.0).OR.(NTIP.GT.2)) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
            WRITE(ScreenUnit,9972)
            WRITE(ScreenUnit,*) 'NTIP =',NTIP
            WRITE(ScreenUnit,9710)
            WRITE(ScreenUnit,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'NTIP =',NTIP
         WRITE(16,9710)
         WRITE(16,9973)
 9710    FORMAT(/,1X,'Your selection of NTIP (a UNIT 15 input ',
     &        'parameter) is not an allowable value')
#ifdef CMPI
         call msg_fini()
#endif
         STOP
      ENDIF
      IF((ICS.EQ.1).AND.(NTIP.GE.1)) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
            WRITE(ScreenUnit,9972)
            WRITE(ScreenUnit,*) 'NTIP =',NTIP
            WRITE(ScreenUnit,9711)
            WRITE(ScreenUnit,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'NTIP =',NTIP
         WRITE(16,9711)
         WRITE(16,9973)
 9711    FORMAT(/,1X,'Your selection of NTIP (a UNIT 15 input ',
     &        'parameter) is inconsistent with your ',
     &        /,1X,'selection of coordinate systems.  Tidal',
     &        'potential forcing should be used only with ',
     &        /,1X,'Spherical coordinates')
         STOP
      ENDIF
      IF (NTIP.NE.0) CTIP = .TRUE.
      IF(NTIP.EQ.0) THEN
         WRITE(16,235) NTIP
 235     FORMAT(/,5X,'NTIP = ',I2,/,9X,
     &        'TIDAL POTENTIAL FORCING IS NOT USED IN THE COMPUTATION')
      ENDIF
      IF(NTIP.GE.1) THEN
         WRITE(16,236) NTIP
 236     FORMAT(/,5X,'NTIP = ',I2,
     &    /,9X,'TIDAL POTENTIAL FORCING IS USED IN THE COMPUTATION ',
     &   'BASED ON INPUT LONGITUDES/LATITUDES')
      ENDIF
      IF(NTIP.EQ.2) THEN
         WRITE(16,239)
 239     FORMAT(9X,'SELF ATTRACTION/LOAD TIDE FORCING IS ALSO USED ',
     &        'IN THE COMPUTATION')
      ENDIF
C...  
C...  READ AND PROCESS NWS - WIND AND PRESSURE FORCING & WAVE RADIATION
C...  STRESS FORCING
C...  
C     jgf46.00 Added NWS=7 (direct surface stress).
C     jgf46.03 Added NWS=8 (Holland wind model)
C     jgf46.16 Merged: 
C     cf & cm  Added NWS=9 (asymmetric hurricane winds)
C     sb46.28sb01 Added NWS=12 (OWI format) 09/xx/2006
C     sb46.28sb03 Added NWS=2xx for STWAVE output direct read 09/xx/2006
      READ(15,*) NWS
      IF((NWS.NE.0)    .AND.(NWS.NE.1)       .AND.(ABS(NWS).NE.2)  .AND.
     &   (NWS.NE.3)    .AND.(ABS(NWS).NE.4)  .AND.(ABS(NWS).NE.5)  .AND.
     &   (NWS.NE.6)    .AND.(NWS.NE.7)       .AND.(NWS.NE.8)       .AND.
     &   (NWS.NE.9)    .AND.(NWS.NE.10)      .AND.(NWS.NE.11)      .AND.
     &   (ABS(NWS).NE.12)                                          .AND.
     &   (NWS.NE.100)  .AND.(NWS.NE.101)     .AND.(NWS.NE.102)     .AND.
     &   (NWS.NE.103)  .AND.(ABS(NWS).NE.104).AND.(ABS(NWS).NE.105).AND.
     &   (NWS.NE.106)  .AND.(NWS.NE.110)     .AND.(NWS.NE.111)     .AND.
     &   (ABS(NWS).NE.112)                                         .AND.
     &   (NWS.NE.200)  .AND.(NWS.NE.201)     .AND.(NWS.NE.202)     .AND.
     &   (NWS.NE.203)  .AND.(ABS(NWS).NE.204).AND.(ABS(NWS).NE.205).AND.
     &   (NWS.NE.206)  .AND.(NWS.NE.210)     .AND.(NWS.NE.211)     .AND.
     &   (ABS(NWS).NE.212)) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
            WRITE(ScreenUnit,9972)
            WRITE(ScreenUnit,*) 'NWS =',NWS
            WRITE(ScreenUnit,9712)
            WRITE(ScreenUnit,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'NWS =',NWS
         WRITE(16,9712)
         WRITE(16,9973)
 9712    FORMAT(/,1X,'Your selection of NWS (a UNIT 15 input ',
     &        'parameter) is not an allowable value')
#ifdef CMPI
         call msg_fini()
#endif
         STOP
      ENDIF

C.... SET WAVE RADIATION STRESS FLAG AND ADJUST NWS ACCORDINGLY

      NRS=0
      IF(ABS(NWS/100).EQ.1) THEN ! sb46.28sb03
         NRS=1
         NWS=(ABS(NWS)-100)*(NWS/ABS(NWS))
      ENDIF
C     sb46.28sb03 Added NWS=2xx for STWAVE output direct read 09/xx/2006
      IF(ABS(NWS/100).EQ.2) THEN
         NRS=2
         NWS=(ABS(NWS)-200)*(NWS/ABS(NWS))
      ENDIF
      IF(NWS.EQ.0) THEN
         WRITE(16,237) NWS
 237     FORMAT(/,5X,'NWS = ',I2,/,9X,
     &        'WIND STRESS OR SURFACE PRESSURE ARE NOT USED TO FORCE',
     &        'THE COMPUTATION')
      ENDIF
      IF(NWS.EQ.1) THEN
         WRITE(16,238) NWS
 238     FORMAT(/,5X,'NWS = ',I2,
     &        /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',
     &        /,9X,' THE COMPUTATION',
     &        /,9X,'VALUES ARE READ AT ADCIRC GRID NODES FROM UNIT 22',
     &        /,9X,' EVERY MODEL TIME STEP')
      ENDIF
      IF(NWS.EQ.2) THEN
         WRITE(16,2381) NWS
 2381    FORMAT(/,5X,'NWS = ',I2,
     &    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',
     &    /,9X,' THE COMPUTATION',
     &    /,9X,'VALUES ARE READ AT ADCIRC GRID NODES FROM UNIT 22',
     &    /,9X,'THE UNIT 22 FILE BEGINS AT TIME=STATIM.',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.')
      ENDIF
      IF(NWS.EQ.-2) THEN
         WRITE(16,2380) NWS
 2380    FORMAT(/,5X,'NWS = ',I2,
     &    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',
     &    /,9X,' THE COMPUTATION',
     &    /,9X,'VALUES ARE READ AT ADCIRC GRID NODES FROM UNIT 22',
     &    /,9X,'THE UNIT 22 FILE BEGINS AT THE TIME OF THE HOT START.',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.')
      ENDIF
      IF(NWS.EQ.3) THEN
         WRITE(16,2382) NWS
 2382    FORMAT(/,5X,'NWS = ',I2,
     &    /,9X,'WIND STRESS ONLY IS USED TO FORCE THE COMPUTATION.',
     &    /,9X,'WIND SPEEDS AND DIRECTIONS ARE READ FROM A FLEET ',
     &    /,9X,'NUMERIC FORMAT FILE AT UNIT 22 AND INTERPOLATED TO',
     &    /,9X,'THE ADCIRC GRID. ',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.',
     &    /,9X,'WIND SPEEDS ARE CONVERTED TO STRESS USING THE GARRET ',
     &     'DRAG LAW.')
      ENDIF
      IF(NWS.EQ.4) THEN
         WRITE(16,2383) NWS
 2383    FORMAT(/,5X,'NWS = ',I2,
     &    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',
     &    /,9X,' THE COMPUTATION',
     &    /,9X,'WIND VELOCITY AND PRESSURE VALUES ARE READ AT SELECTED',
     &    /,9X,'ADCIRC GRID NODES FROM A PBL FILE AT UNIT 22.',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.',
     &    /,9X,'THE UNIT 22 FILE BEGINS AT TIME=STATIM.',
     &    /,9X,'WIND SPEEDS ARE CONVERTED TO STRESS USING THE GARRET ',
     &    'DRAG LAW.')
      ENDIF
      IF(NWS.EQ.-4) THEN
         WRITE(16,2388) NWS
 2388    FORMAT(/,5X,'NWS = ',I2,
     &    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',
     &    /,9X,' THE COMPUTATION',
     &    /,9X,'WIND VELOCITY AND PRESSURE VALUES ARE READ AT SELECTED',
     &    /,9X,'ADCIRC GRID NODES FROM A PBL FILE AT UNIT 22.',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.',
     &    /,9X,'THE UNIT 22 FILE BEGINS AT THE TIME OF THE HOT START.',
     &    /,9X,'WIND SPEEDS ARE CONVERTED TO STRESS USING THE GARRET ',
     &    'DRAG LAW.')
      ENDIF
      IF(NWS.EQ.5) THEN
         WRITE(16,2384) NWS
 2384    FORMAT(/,5X,'NWS = ',I2,
     &    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',
     &    /,9X,' THE COMPUTATION',
     &    /,9X,'WIND VELOCITY AND PRESSURE VALUES ARE READ AT ADCIRC ',
     &    /,9X,'GRID NODES FROM UNIT 22',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.',
     &    /,9X,'THE UNIT 22 FILE BEGINS AT TIME=STATIM.',
     &    /,9X,'WIND SPEEDS ARE CONVERTED TO STRESS USING THE GARRET ',
     &    'DRAG LAW.')
      ENDIF
      IF(NWS.EQ.-5) THEN
         WRITE(16,2389) NWS
 2389    FORMAT(/,5X,'NWS = ',I2,
     &    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',
     &    /,9X,' THE COMPUTATION',
     &    /,9X,'WIND VELOCITY AND PRESSURE VALUES ARE READ AT ADCIRC ',
     &    /,9X,'GRID NODES FROM UNIT 22',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.',
     &    /,9X,'THE UNIT 22 FILE BEGINS AT THE TIME OF THE HOT START.',
     &    /,9X,'WIND SPEEDS ARE CONVERTED TO STRESS USING THE GARRET ',
     &    'DRAG LAW.')
      ENDIF
      IF(NWS.EQ.6) THEN
         WRITE(16,2385) NWS
 2385    FORMAT(/,5X,'NWS = ',I2,
     &    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',
     &    /,9X,' THE COMPUTATION',
     &    /,9X,'WIND VELOCITY AND PRESSURE VALUES ARE READ FROM A ',
     &    /,9X,'REGULARLY SPACED GRID FROM UNIT 22',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP AND IN SPACE TO BRING THE ',
     &    /,9X,'MET DATA FROM A REGULAR GRID TO THE ADCIRC GRID.'
     &    /,9X,'WIND SPEEDS ARE CONVERTED TO STRESS USING THE GARRET ',
     &    'DRAG LAW.')
      ENDIF
C     jgf46.00 Added NWS=7 (direct surface stress).
      IF(NWS.EQ.7) THEN
         WRITE(16,1234) NWS
 1234    FORMAT(/,5X,'NWS = ',I2,
     &    /,9X,'SURFACE STRESS AND PRESSURE ARE USED TO FORCE',
     &    /,9X,' THE COMPUTATION',
     &    /,9X,'VALUES ARE READ AT ADCIRC GRID NODES FROM UNIT 22',
     &    /,9X,'THE UNIT 22 FILE BEGINS AT TIME=STATIM.',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE STRESS DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.')
      ENDIF
      IF(NWS.EQ.-7) THEN
         WRITE(16,1235) NWS
 1235    FORMAT(/,5X,'NWS = ',I2,
     &    /,9X,'SURFACE STRESS AND SURFACE PRESSURE ARE USED TO FORCE',
     &    /,9X,' THE COMPUTATION',
     &    /,9X,'VALUES ARE READ AT ADCIRC GRID NODES FROM UNIT 22',
     &    /,9X,'THE UNIT 22 FILE BEGINS AT THE TIME OF THE HOT START.',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE STRESS DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.')
      ENDIF
C     jgf46.03 Added NWS=8 (Holland wind model).
      IF(NWS.EQ.8) THEN
         WRITE(16,1237) NWS
 1237    FORMAT(/,5X,'NWS = ',I2,
     &    /,9X,'HURRICANE PARAMETERS AND THE HOLLAND WIND MODEL',
     &    /,9X,'  ARE USED TO FORCE THE COMPUTATION',
     &    /,9X,'VALUES ARE READ AT FOR THE STORM FROM UNIT 22',
     &    /,9X,'THE UNIT 22 FILE BEGINS AT TIME=STATIM.',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE STORM DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.')
      ENDIF
      IF(NWS.EQ.-8) THEN
         WRITE(16,1238) NWS
 1238    FORMAT(/,5X,'NWS = ',I2,
     &    /,9X,'HURRICANE PARAMETERS AND THE HOLLAND WIND MODEL',
     &    /,9X,'  ARE USED TO FORCE THE COMPUTATION',
     &    /,9X,'VALUES ARE READ AT FOR THE STORM FROM UNIT 22',
     &    /,9X,'THE UNIT 22 FILE BEGINS AT THE TIME OF THE HOT START.',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE STORM DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.')
      ENDIF
C     jgf46.16 Merged:
C     cf & cm added nws = 9: asymmetric hurricane winds
      IF(NWS.EQ.9) THEN                                           
         WRITE(16,2399) NWS
 2399    FORMAT(/,5X,'NWS = ',I2,
     &    /,9X,'HURRICANE PARAMETERS AND THE ASYMMETRIC WIND MODEL',
     &    /,9X,'  ARE USED TO FORCE THE COMPUTATION',
     &    /,9X,'VALUES ARE READ AT FOR THE STORM FROM UNIT 22',
     &    /,9X,'THE UNIT 22 FILE BEGINS AT TIME=STATIM.',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE STORM DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.')
      ENDIF
      IF(NWS.EQ.-9) THEN
         WRITE(16,2398) NWS
 2398    FORMAT(/,5X,'NWS = ',I2,
     &    /,9X,'HURRICANE PARAMETERS AND THE ASYMMETRIC WIND MODEL',
     &    /,9X,'  ARE USED TO FORCE THE COMPUTATION',
     &    /,9X,'VALUES ARE READ AT FOR THE STORM FROM UNIT 22',
     &    /,9X,'THE UNIT 22 FILE BEGINS AT THE TIME OF THE HOT START',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE STORM DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.')
      ENDIF
      IF(NWS.EQ.10) THEN
         WRITE(16,2386) NWS
 2386    FORMAT(/,5X,'NWS = ',I2,
     &    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',
     &    /,9X,' THE COMPUTATION',
     &    /,9X,'WIND VELOCITY AND PRESSURE VALUES ARE READ EVERY N',
     &    /,9X,' HOURS FROM A DIFFERENT FILE AT UNITS 200, 200+N,',
     &   ' 200+2N, ETC.',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP AND IN SPACE TO BRING THE ',
     &    /,9X,'MET DATA FROM A GAUSSIAN GRID TO THE ADCRIC GRID.',
     &    /,9X,'WIND SPEEDS ARE CONVERTED TO STRESS USING THE GARRET ',
     &    'DRAG LAW.')
      ENDIF
      IF(NWS.EQ.11) THEN
         WRITE(16,2387) NWS
 2387    FORMAT(/,5X,'NWS = ',I2,
     &    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',
     &    /,9X,' THE COMPUTATION',
     &    /,9X,'WIND VELOCITY AND PRESSURE VALUES ARE READ EVERY 3 ',
     &    /,9X,'HOURS FROM ETA-29 FILES AT UNITS 200, 201, 202, ETC.',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP AND IN SPACE TO BRING THE ',
     &    /,9X,'WIND DATA FROM THE 29 KM E GRID TO THE ADCRIC GRID.',
     &    /,9X,'WIND SPEEDS ARE CONVERTED TO STRESS USING THE GARRET ',
     &   'DRAG LAW.')
      ENDIF
C.....sb46_28sb01 added for NWS=-12,12 09/xx/2006
      IF(NWS.EQ.12) THEN
         WRITE(16,12384) NWS
12384    FORMAT(/,5X,'NWS = ',I2,
     &    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',
     &    /,9X,' THE COMPUTATION',
     &    /,9X,'WIND VELOCITY AND PRESSURE VALUES ARE READ FROM RAW ',
     &    /,9X,'OWI DATA FILES (UNIT 221-224).',
     &    /,9X,'META DATA IS READ FROM UNIT 220.',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.',
     &    /,9X,'THE UNIT 22 FILE BEGINS AT TIME=STATIM.',
     &    /,9X,'WIND SPEEDS ARE CONVERTED TO STRESS USING THE GARRET ',
     &         'DRAG LAW.')
      ENDIF
      IF(NWS.EQ.-12) THEN
         WRITE(16,12389) NWS
12389    FORMAT(/,5X,'NWS = ',I2,
     &    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',
     &    /,9X,' THE COMPUTATION',
     &    /,9X,'WIND VELOCITY AND PRESSURE VALUES ARE READ FROM RAW ',
     &    /,9X,'OWI DATA FILES (UNIT 221-224).',
     &    /,9X,'META DATA IS READ FROM UNIT 220.',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.',
     &    /,9X,'THE UNIT 22 FILE BEGINS AT THE TIME OF THE HOT START.',
     &    /,9X,'WIND SPEEDS ARE CONVERTED TO STRESS USING THE GARRET ',
     &         'DRAG LAW.')
      ENDIF
      IF(NRS.EQ.0) THEN
         WRITE(16,2390) NRS
 2390    FORMAT(/,5X,'NRS = ',I2,
     &        /,9X,'WAVE RADIATION STRESS IS NOT USED TO FORCE THE ',
     &        'COMPUTATION')
      ENDIF
      IF(NRS.EQ.1) THEN
         WRITE(16,2391) NRS
 2391    FORMAT(/,5X,'NRS = ',I2,
     &    /,9X,'WAVE RADIATION STRESS IS USED TO FORCE THE COMPUTATION',
     &    /,9X,'STRESSES ARE READ AT SELECTED ADCIRC GRID NODES FROM A',
     &    /,9X,'PBL TYPE FILE AT UNIT 23.  INTERPOLATION IN TIME IS ',
     &    /,9X,'DONE TO SYNC THE STRESS DATA WITH THE MODEL TIME STEP.',
     &    /,9X,'FOR A COLD START, THE UNIT 23 FILE BEGINS AT THE TIME ',
     &    /,9X,'OF THE COLD START.  FOR A HOT START, THE UNIT 23 FILE ',
     &    /,9X,'BEGINS AT THE TIME OF THE HOT START.')
      ENDIF
C     sb46.28sb03 Added NWS=2xx for STWAVE output direct read 09/xx/2006
      IF(NRS.EQ.2) THEN
         WRITE(16,2392) NRS
 2392    FORMAT(/,5X,'NRS = ',I2,
     &    /,9X,'WAVE RADIATION STRESS IS USED TO FORCE THE COMPUTATION',
     &    /,9X,'STRESSES ARE READ AT SELECTED ADCIRC GRID NODES FROM A',
     &    /,9X,'UNIT 23 FILE GENERATED BY UTIL/BUILDSTWAVE23. ',
     &    /,9X,'NO DECOMPOSITION IS NEEDED FOR THIS UNIT 23 FILE EVEN',
     &    /,9X,'IN A PARALLEL EXCECUTION.  INTERPOLATION IN TIME IS ',
     &    /,9X,'DONE TO SYNC THE STRESS DATA WITH THE MODEL TIME STEP.',
     &    /,9X,'FOR A COLD START, THE UNIT 23 FILE BEGINS AT THE TIME ',
     &    /,9X,'OF THE COLD START.  FOR A HOT START, THE UNIT 23 FILE ',
     &    /,9X,'BEGINS AT THE TIME OF THE HOT START.')
      ENDIF
      
C...  
C...  READ AND PROCESS NRAMP - WHETHER A RAMP FUNCTION WILL BE USED
C...  
C     jgf46.08 Change to the number of ramp functions that will be used.
      READ(15,*) NRAMP
      IF((NRAMP.NE.0).AND.(NRAMP.gt.7)) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
            WRITE(ScreenUnit,9972)
            WRITE(ScreenUnit,*) 'NRAMP =',NRAMP
            WRITE(ScreenUnit,9713)
            WRITE(ScreenUnit,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'NRAMP =',NRAMP
         WRITE(16,9713)
         WRITE(16,9973)
 9713    FORMAT(/,1X,'Your selection of NRAMP (a UNIT 15 input ',
     &        'parameter) is not an allowable value')
#ifdef CMPI
         call msg_fini()
#endif
         STOP
      ENDIF
      IF(NRAMP.EQ.0) THEN
         WRITE(16,240) NRAMP
 240     FORMAT(/,5X,'NRAMP = ',I2,
     &        /,9X,'NO RAMP FUNCTION IS USED IN THE COMPUTATION')
      ELSE
         WRITE(16,241) NRAMP
 241     FORMAT(/,5X,'NRAMP = ',I2,
     &    /,9X,'HYPERBOLIC TANGENT RAMP(S) WILL BE APPLIED TO THE ',
     &    'FORCING FUNCTIONS')
      ENDIF
C...
C...  PROCESS G - GRAVITY
C...
      READ(15,*) G
      IF((ICS.EQ.2).AND.(abs(G-9.81).gt.0.01)) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
            WRITE(ScreenUnit,9972)
            WRITE(ScreenUnit,*) 'G =',G
            WRITE(ScreenUnit,9714)
            WRITE(ScreenUnit,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'G =',G
         WRITE(16,9714)
         WRITE(16,9973)
 9714    FORMAT(/,1X,'Your specification of the gravitational ',
     &        'constant, G, (a UNIT 15 input) is not ',
     &        /,1X,'consistant with the use of spherical coordinates.',
     &        '  G must be in units of m/s^2')
#ifdef CMPI
         call msg_fini()
#endif
         STOP
      ENDIF
      WRITE(16,5) G
    5 FORMAT(///,5X,'GRAVITATIONAL CONSTANT G =',F10.5,/)

C...
C...  READ AND PROCESS TAU0 - WEIGHTING COEFFICIENT IN THE GWCE
C
C...  jgf45.12 Added three-tier tau0 scheme.  
C
C     jgf46.00 Added user specified spatially varying tau0 scheme using
C     the Nodal Attributes File (unit 13) and removed three tier tau0
C     scheme. Since Tau0 is a nodal attribute, it is initialized later
C     in this subroutine using a CALL to the InitNodalAttr subroutine in
C     the NodalAttributes module.
      READ(15,*) Tau0
C
C     jgf47.11 Added a line to the fort.15 file to read in the min
C     and max tau0 values, if the user has chosen to use the time
C     varying tau0.
      IF ( (Tau0.eq.-3.d0).or.(Tau0.eq.-4.d0) ) THEN
        READ(15,*) Tau0FullDomainMin, Tau0FullDomainMax
      ENDIF
C...  
C...  INPUT FROM UNIT 15 AND OUTPUT TO UNIT 16 TIME INTEGRATION
C...  INFORMATION INCLUDING DT,STATIM,REFTIM,AND RNDAY
C...  
      WRITE(16,1112)
      WRITE(16,245)
  245 FORMAT(//,1X,'TIME INTEGRATION INFORMATION',//)

C...
C...  READ AND PROCESS DT - MODEL TIME STEP
C...
ckmd	 Changed the time step to allow for negative values in 
ckmd 	 order to turn on the predictor-corrector algorithm. 
      READ(15,*) DTDP
      IF(DTDP.LT.0.d0) THEN
         CPRECOR    = .TRUE.
         CGWCE_Orig = .FALSE.
         CGWCE_New  = .FALSE. !jgf Turn off the default.
         CME_New_NC = .FALSE. !jgf Turn off the default.
         CME_New_C1 = .FALSE. 
         CME_New_C2 = .FALSE. 
         DT=-DTDP
         DTDP=DT
         WRITE(16,*) ' ADCIRC is configured for a 2DDI run using'
         WRITE(16,*) ' the predictor-corrector algorithm and'
         WRITE(16,*) ' the ADCIRC logical variable is set to:    '
         WRITE(16,*) '         CPRECOR           = ',CPRECOR
      ELSE IF(DTDP.GT.0.d0) THEN
         DT=DTDP
         WRITE(16,*) ' ADCIRC is configured for a 2DDI run '
         WRITE(16,*) ' without the predictor-corrector algorithm and'
         WRITE(16,*) ' the ADCIRC logical variable is set to:  '   
         WRITE(16,*) '         CPRECOR           = ',CPRECOR
      END IF
      WRITE(16,9) DTDP
    9 FORMAT(5X,'TIME STEP =',F12.6,5X,'SECONDS',/)

C...
C...  READ AND PROCESS STATIM - SIMULATION STARTING TIME
C...  
      READ(15,*) STATIM
      WRITE(16,1113) STATIM
 1113 FORMAT(5X,'STARTING TIME FOR SIMULATION = ',F14.6,' DAYS',/)

C...
C...  READ AND PROCESS REFTIM - Harmonic REFERNCE TIME
C...
      READ(15,*) REFTIM
      WRITE(16,1115) REFTIM
 1115 FORMAT(5X,'Harmonic REFERENCE TIME = ',F14.6,' DAYS',/)

C...
C...  Read in and process additional timing information for wind.
C...
      IF((NWS.EQ.0).AND.(NRS.GE.1)) READ(15,*) RSTIMINC ! sb46.28sb03
      IF((NWS.EQ.1).AND.(NRS.GE.1)) READ(15,*) RSTIMINC ! sb46.28sb03
      IF(ABS(NWS).EQ.2) THEN
         IF(NRS.EQ.0) READ(15,*) WTIMINC
         IF(NRS.GE.1) READ(15,*) WTIMINC,RSTIMINC ! sb46.28sb03
      ENDIF
      IF(NWS.EQ.3) THEN
         READ(15,*) IREFYR,IREFMO,IREFDAY,IREFHR,IREFMIN,REFSEC
         WRITE(16,1116) IREFMO,IREFDAY,IREFYR,IREFHR,IREFMIN,REFSEC
 1116    FORMAT(5X,'WIND REFERENCE TIME FOR SIMULATION = ',
     &        I2,'/',I2,'/',I2,'  ',I2,':',I2,':',f7.4,/)
         CALL TIMECONV(IREFYR,IREFMO,IREFDAY,IREFHR,IREFMIN,REFSEC,
     &        WREFTIM, MyProc, NScreen, ScreenUnit)
         IF(NRS.EQ.0) READ(15,*) NWLAT,NWLON,WLATMAX,WLONMIN,WLATINC,
     &        WLONINC,WTIMINC
         IF(NRS.GE.1) READ(15,*) NWLAT,NWLON,WLATMAX,WLONMIN,WLATINC,  ! sb46.28sb03
     &        WLONINC,WTIMINC,RSTIMINC
      ENDIF
      IF(ABS(NWS).EQ.4) THEN
         IF(NRS.EQ.0) READ(15,*) WTIMINC
         IF(NRS.GE.1) READ(15,*) WTIMINC,RSTIMINC ! sb46.28sb03
      ENDIF
      IF(ABS(NWS).EQ.5) THEN
         IF(NRS.EQ.0) READ(15,*) WTIMINC
         IF(NRS.GE.1) READ(15,*) WTIMINC,RSTIMINC ! sb46.28sb03
      ENDIF
      IF(NWS.EQ.6) THEN
         IF(NRS.EQ.0) READ(15,*) NWLAT,NWLON,WLATMAX,WLONMIN,WLATINC,
     &        WLONINC,WTIMINC
         IF(NRS.GE.1) READ(15,*) NWLAT,NWLON,WLATMAX,WLONMIN,WLATINC, ! sb46.28sb03
     &        WLONINC,WTIMINC,RSTIMINC
      ENDIF
C     jgf46.00 Added NWS=7 (direct surface stress).
      IF(ABS(NWS).EQ.7) THEN
         IF(NRS.EQ.0) READ(15,*) WTIMINC
         IF(NRS.GE.1) READ(15,*) WTIMINC,RSTIMINC ! sb46.28sb03
      ENDIF
C     jgf46.05 Added NWS=8 (Holland Wind Model).
C     jgf46.28 Changed WTIMINC to StormNumber for activating 
C     wind multiplier to final wind speeds from Holland model.
      IF(ABS(NWS).EQ.8) THEN
         IF(NRS.EQ.0) THEN
            READ(15,*) IREFYR,IREFMO,IREFDAY,IREFHR,StormNumber,BLAdj
         ELSEIF (NRS.EQ.1) THEN
            READ(15,*) IREFYR,IREFMO,IREFDAY,IREFHR,StormNumber,BLAdj,
     &           RSTIMINC
            WRITE(16,6111) IREFMO,IREFDAY,IREFYR,IREFHR
 6111       FORMAT(5X,'WIND REFERENCE TIME FOR SIMULATION = ',
     &           I2,'/',I2,'/',I2,'  ',I2,'H',/)
         ENDIF
         CALL TIMECONV(IREFYR,IREFMO,IREFDAY,IREFHR,0,0.0d0,
     &        WindRefTime, MyProc, NScreen, ScreenUnit)
      ENDIF
C     jgf46.16 Merged:
C     cf & cm added nws = 9: asymmetric hurricane winds
      IF(ABS(NWS).EQ.9) THEN
         IF(NRS.EQ.0) READ(15,*) WTIMINC
         IF(NRS.GE.1) READ(15,*) WTIMINC,RSTIMINC ! sb46.28sb03
         WRITE(16,1117) WTIMINC
      ENDIF
      IF(NWS.EQ.10) THEN
         NWLAT=190
         NWLON=384
         IF(NRS.EQ.0) READ(15,*) WTIMINC
         IF(NRS.GE.1) READ(15,*) WTIMINC,RSTIMINC ! sb46.28sb03
      ENDIF
      IF(NWS.EQ.11) THEN
         NWLAT=271
         NWLON=181
         WTIMINC=10800.
         IF(NRS.GE.1) READ(15,*) RSTIMINC ! sb46.28sb03
         !READ(15,*) NWLAT,NWLON,WTIMINC
      ENDIF
C     sb46.28sb01 Added NWS=12 (OWI format)      
      IF(ABS(NWS).EQ.12) THEN
         IF(NRS.EQ.0) READ(15,*) WTIMINC
         IF(NRS.GE.1) READ(15,*) WTIMINC,RSTIMINC ! sb46.28sb03
      ENDIF

      IF(NWS.NE.0) WRITE(16,1117) WTIMINC
 1117 FORMAT(5X,'WIND TIME INCREMENT (SEC) = ',F10.2,/)
      IF(NRS.NE.0) WRITE(16,1118) RSTIMINC
 1118 FORMAT(5X,'RADIATION STRESS TIME INCREMENT (SEC) = ',F10.2,/)

C...
C...  READ AND PROCESS RNDAY - SIMULATION DURATION IN DAYS
C...
      READ(15,*) RNDAY
      WRITE(16,10) RNDAY
 10   FORMAT(5X,'TOTAL LENGTH OF NUMERICAL SIMULATION =',F12.4,
     &       5X,'DAYS',/)

C...
C...  COMPUTE TOTAL NUMBER OF TIME STEPS NT
C...
#ifdef IBM
      NT=INT(RNDAY*(86400.D0/DTDP)+0.5d0,KIND(0.0d0))
#else
      NT=INT(RNDAY*(86400.D0/DTDP)+0.5d0)
#endif
      WRITE(16,1920) NT
 1920 FORMAT(5X,'NUMBER OF TIME STEPS  =',I8,/)
C...  
C...  READ AND PROCESS EFFECTIVE LENGTH OF THE HYPERBOLIC TANGENT RAMP(S)
C...  IN DAYS
C...  
C     jgf46.08 Add fine-grained ramp functions.
C     jgf46.21 Add FluxSettlingTime for IBTYPE=52 to accomodate
C     MS river during Katrina, split ramps for flux b.c.s into internal
C     and external.
      FluxSettlingTime = 0.0d0
      DRamp = 1.0d0
      SELECT CASE(NRamp)
C     ---------      
      CASE(0,1)! Either no ramp, or same ramp for all forcings
C     ---------      
         READ(15,*) DRamp
         DRampIntFlux=DRamp
         DRampExtFlux=DRamp
         DRampElev=DRamp
         DRampTip=DRamp
         DRampMete=DRamp
         DRampWRad=DRamp
C     -------
      CASE(2) ! Ramp for external flux boundary conditions.
C     -------
         READ(15,*) DRamp,DRampExtFlux,FluxSettlingTime
         DRampIntFlux=DRamp
         DRampElev=DRamp
         DRampTip=DRamp
         DRampMete=DRamp
         DRampWRad=DRamp
C     -------
      CASE(3) ! Ramp for internal flux boundary conditions.
C     -------
         READ(15,*) DRamp,DRampExtFlux,FluxSettlingTime,DRampIntFlux
         DRampElev=DRamp
         DRampTip=DRamp
         DRampMete=DRamp
         DRampWRad=DRamp
C     -------
      CASE(4) ! Ramp for surface elevation specified boundary conditions.
C     -------
         READ(15,*) DRamp,DRampExtFlux,FluxSettlingTime,DRampIntFlux,
     &        DRampElev
         DRampTip=DRamp
         DRampMete=DRamp
         DRampWRad=DRamp
C     -------
      CASE(5) ! Ramp for tidal potential 
C     -------
         READ(15,*) DRamp,DRampExtFlux,FluxSettlingTime,DRampIntFlux,
     &        DRampElev,DRampTip
         DRampMete=DRamp
         DRampWRad=DRamp
C     -------
      CASE(6) ! Ramp for wind and atmospheric pressure  
C     -------
         READ(15,*) DRamp,DRampExtFlux,FluxSettlingTime,DRampIntFlux,
     &        DRampElev,DRampTip,DRampMete
         DRampWRad=DRamp
C     -------
      CASE(7) ! Ramp for wave radiation stress
C     -------
         READ(15,*) DRamp,DRampExtFlux,FluxSettlingTime,DRampIntFlux,
     &        DRampElev,DRampTip,DRampMete,DRampWRad
C     ------------
      CASE DEFAULT ! fall-through 
C     ------------
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
            WRITE(ScreenUnit,9972)
            WRITE(ScreenUnit,*) 'NRAMP =',NRAMP
            WRITE(ScreenUnit,9713)
            WRITE(ScreenUnit,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'NRAMP =',NRAMP
         WRITE(16,9713)
         WRITE(16,9973)
#ifdef CMPI
         call msg_fini()
#endif
         STOP
C
      END SELECT
C     
      IF(NRAMP.NE.0) THEN
         WRITE(16,8763) DRAMP
 8763    FORMAT(/,5X,'VALUE FOR DRAMP USED IN RAMP EVALUATION =',F12.4,
     &        5X,'DAYS',/)

         WRITE(16,5841)
 5841    FORMAT(11X,' DAYS OF SIMULATION',2X,' TIME  ',6X,'  RAMP VALUE'
     &                                                               ,/)
         Day=0.0d0
 999     RampVal=TANH(Day*2.d0/DRAMP)
         WRITE(16,5845) Day,Day+StaTim,RampVal
 5845    FORMAT(15X,F8.2,6X,F8.2,2X,F15.7)
         DAY=DAY+0.5d0
         IF(Day.LT.DRAMP*1.25) GOTO 999
      ENDIF
C...
C...  READ GWCE TIME WEIGHTING FACTORS
C...
      READ(15,*) A00,B00,C00
      WRITE(16,14)
14    FORMAT(//,5X,'TIME WEIGHTING FACTORS IN THE WAVE EQUATION :'/)
      WRITE(16,15) A00,B00,C00
15    FORMAT(9X,'AT TIME LEVEL K+1 : ',F8.5,
     &  /,9X,'AT TIME LEVEL K   : ',F8.5,
     &  /,9X,'AT TIME LEVEL K-1 : ',F8.5,/)
C...
C...  READ MINIMUM DEPTH OR WET/DRY PARAMETERS FROM UNIT 15
C...
      IF(NOLIFA.NE.2) THEN
         READ(15,*) H0
         WRITE(16,16) H0
 16      FORMAT(//,5X,'THE BATHYMETRIC DEPTH AT ALL NODES WILL BE ',
     &               'INCREASED TO H0= ',F12.4,' IF NECESSARY'/)
      ENDIF
      IF(NOLIFA.EQ.2) THEN
         READ(15,*) H0,NODEDRYMIN,NODEWETMIN,VELMIN
         WRITE(16,17) H0,NODEWETMIN,VELMIN,NODEDRYMIN
 17      FORMAT(//,5X,'DRYING WILL OCCUR WHEN THE WATER DEPTH < H0',
     &          /,5X,'H0 = ',F10.6,
     &          /,5X,'AND NODEREP > NODEWETMIN = ',I6,' TIME STEPS',
     &          /,5X,'NODEREP = NUMBER OF TIME STEPS SINCE A NODE ',
     &               'CHANGED STATE (EITHER WETTED OR DRIED)',
     &         //,5X,'WETTING WILL OCCUR WHEN THERE IS A FAVORABLE ',
     &               'PRESSURE GRADIENT THAT',
     &          /,5X,'WOULD DRIVE A STEADY VELOCITY TOWARDS A DRY NODE',
     &          /,5X,'THAT IS GREATER THAN VELMIN = ',F10.5,
     &          /,5X,'AND NODEREP > NODEDRYMIN = ',I6,' TIME STEPS',/)
      ENDIF
C
C     jgf46.00 Read longitude and latitude on which the CPP coordinate
C     projection is centered (in degrees) if ICS = 2. (The reading of
C     the top of the grid file, including NE and NP, was moved nearer to
C     the beginning of this subroutine.)
C
      READ(15,*) SLAM0,SFEA0
      SLAM0=SLAM0*DEG2RAD
      SFEA0=SFEA0*DEG2RAD
      WRITE(16,1112)
      WRITE(16,246)
 246  FORMAT(//,1X,'GRID INFORMATION',//)

C ALLOCATE ARRAYS Dimensioned by MNP and MNE
      call alloc_main1a()
C...
C...  IF ICS=1 INPUT NODAL COORDINATES AND BATHYMETRY FROM UNIT 14
C.... IF EITHER NTIP=1 OR NCOR=1, COMPUTE THE INVERSE CPP PROJECTION
C...  IF ICS=2 INPUT NODAL COORDINATES AND BATHYMETRY FROM UNIT 14
C.... AND COMPUTE CPP PROJECTED COORDINATES
C...
      IF(ICS.EQ.1) THEN
         DO I=1,NP
            READ(14,*) JKI,X(JKI),Y(JKI),DP(JKI)
            IF(JKI.NE.I) THEN
               IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,99801)
               WRITE(16,99801)
99801          FORMAT(////,1X,'!!!!!!!!!!  WARNING - NONFATAL ',
     &              'INPUT ERROR  !!!!!!!!!',
     &              //,1X,'YOUR NODE NUMBERING IS NOT SEQUENTIAL ',
     &              'CHECK YOUR UNIT 14 INPUT FILE CAREFULLY',//)
            ENDIF
            IF((NTIP.GE.1).OR.(NCOR.EQ.1)) THEN
               CALL INVCP(X(JKI),Y(JKI),SLAM(JKI),SFEA(JKI),SLAM0,SFEA0)
            ENDIF
         END DO
      ENDIF
      IF(ICS.EQ.2) THEN
         DO I=1,NP
            READ(14,*) JKI,SLAM(JKI),SFEA(JKI),DP(JKI)
            IF(JKI.NE.I) THEN
               IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,99801)
               WRITE(16,99801)
            ENDIF
            SLAM(JKI)=DEG2RAD*SLAM(JKI)
            SFEA(JKI)=DEG2RAD*SFEA(JKI)
            CALL CPP(X(JKI),Y(JKI),SLAM(JKI),SFEA(JKI),SLAM0,SFEA0)
         END DO
      ENDIF
C...
C...  IF ICS=1 SET THE SFAC VECTOR EQUAL TO UNITY
C...  IF ICS=2 COMPUTE THE SFAC VECTOR IN ORDER TO ADJUST EQUATIONS TO CPP
C...  COORDINATES
C...
      IF(ICS.EQ.1) THEN
         DO I=1,NP
            SFAC(I)=1.0d0
         END DO
      ENDIF
      IF(ICS.EQ.2) THEN
         DO I=1,NP
            SFAC(I)=COS(SFEA0)/COS(SFEA(I))
         END DO
      ENDIF

C...  
C...  IF WETTING AND DRYING WILL NOT BE USED, MAKE SURE ALL BATHYMETRIC
C...  DEPTHS ARE > OR = TO H0.
C...  
      IF((NOLIFA.EQ.0).OR.(NOLIFA.EQ.1)) THEN
         DO I=1,NP
            IF(DP(I).LT.H0) DP(I)=H0
         END DO
      ENDIF

C...
C...  Read and process the global connectivity table from unit 14
C...
      DO I=1,NE
         READ(14,*) JKI,NHY,NM(JKI,1),NM(JKI,2),NM(JKI,3)
         IF(JKI.NE.I) THEN
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,99802)
            WRITE(16,99802)
99802       FORMAT(////,1X,'!!!!!!!!!!  WARNING - NONFATAL ',
     &                    'INPUT ERROR  !!!!!!!!!',
     &              //,1X,'YOUR ELEMENT NUMBERING IS NOT SEQUENTIAL ',
     &               /,1X,'CHECK YOUR UNIT 14 INPUT FILE CAREFULLY',//)
         ENDIF

C     check that sufficient accuracy is provided by the code to handle
C     the input grid

         X1=X(NM(JKI,1))
         X2=X(NM(JKI,2))
         X3=X(NM(JKI,3))
         Y1=Y(NM(JKI,1))
         Y2=Y(NM(JKI,2))
         Y3=Y(NM(JKI,3))
         AVGXY=(ABS(X1)+ABS(X2)+ABS(X3)+ABS(Y1)+ABS(Y2)+ABS(Y3))/6.d0
         DIF1R=AVGXY/(((X2-X1)**2+(Y2-Y1)**2)**0.5d0)
         DIF2R=AVGXY/(((X3-X2)**2+(Y3-Y2)**2)**0.5d0)
         DIF3R=AVGXY/(((X3-X1)**2+(Y3-Y1)**2)**0.5d0)
         DIF1R=LOG10(DIF1R)
         DIF2R=LOG10(DIF2R)
         DIF3R=LOG10(DIF3R)
         IF((DIF1R.GT.NPREC).OR.(DIF2R.GT.NPREC).OR.
     &                          (DIF3R.GT.NPREC))THEN
           IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9898) JKI
           WRITE(16,9898) JKI
 9898      FORMAT(////,1X,'!!!!!!!!!!  WARNING  !!!!!!!!!',
     &              //,1X,'IF THE GRID COORDINATES HAVE 32 BITS ',
     &                    '(APPROX 7 DIGITS) OF PRECISION',
     &               /,1X,'A ROBUST MODEL SOLUTION CAN NOT BE ',
     &                    'GUARANTEED AT ELEMENT NO. ',I10,
     &              //,1X,'MORE PRECISION MUST BE USED IN THE GRID',//)
           ENDIF

C     compute and store 2 x element areas
         AREAS(JKI)=(X1-X3)*(Y2-Y3)+(X3-X2)*(Y1-Y3)

C     check to insure that the nodal ordering in the connectivity table is 
C     counter clockwise around the elements
         IF(AREAS(JKI).LT.0.0) THEN
           IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9899) JKI
           WRITE(16,9899) JKI
 9899      FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',
     &              //,1X,'IN THE CONNECTIVITY TABLE, THE NODES AROUND',
     &                    ' ELEMENT ',I6,
     &               /,1X,'MUST BE SPECIFIED IN COUNTERCLOCKWISE ORDER',
     &                    ' - CHECK INPUT ',
     &              //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED ',
     &                    '!!!!!!',//)
           STOP
           ENDIF
        END DO
C...
C...  Compute neighbor tables.  
C...

C     estimate the number of neighbor nodes around any node = number of
C     elements containing that node.  This is correct for non-boundary
C     nodes and one too small for boundary nodes

      DO I=1, NP
        NNeigh(I) = 0
        ENDDO
      
      DO I=1,NE
        NNeigh(NM(I,1))=NNeigh(NM(I,1))+1
        NNeigh(NM(I,2))=NNeigh(NM(I,2))+1
        NNeigh(NM(I,3))=NNeigh(NM(I,3))+1
        ENDDO

C     determine the maximum NNeigh

      MNei=0
      DO I=1,NP
        IF(NNeigh(I).GT.MNei) MNei=NNeigh(I)
        ENDDO

C     estimate the maximum array space needed for the neighbor table by 
C     increasing this number by 2, to provide array space for the node itself
C     and in case the maximum number of nodes occurs at a boundary node 

      MNei = MNei+2

C     allocate space for neighbor tables

      call alloc_main1b()

C     compute the neighbor table and redo NNeigh array

      IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1196)
      WRITE(16,1196)
 1196 FORMAT(/,1X,'THE NEIGHBOR TABLE IS BEING COMPUTED ',/)
 
      CALL NEIGHB(NE,NP,NM,NNeigh,NeiTab,NeiTabEle,NeiMin,NeiMax,
     &     X,Y,NScreen)

      IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) 
     &     WRITE(ScreenUnit,1195) NEIMIN,NEIMAX,NEIMAX
      WRITE(16,1195) NEIMIN,NEIMAX,NEIMAX
 1195 FORMAT(1X,'THE NEIGHBOR TABLE IS COMPLETED ',
     &   /,5X,'THE MINIMUM NUMBER OF NEIGHBORS FOR ANY NODE = ',I3,
     &   /,5X,'1+THE MAXIMUM NUMBER OF NEIGHBORS FOR ANY NODE = ',I3,
     &   /,5X,'THE PARAMETER MNEI CAN BE SET AS SMALL AS ',I3,/)

C...
C...IF A 2DDI BAROCLINIC RUN, READ IN INITIAL CONDITION DENSITY FIELDS
C...
      IF((C2DDI).AND.(CBaroclinic)) THEN
         OPEN(11,FILE=TRIM(INPUTDIR)//'/'//'fort.11')
         READ(11,*)
         READ(11,*)
         READ(11,*) NP2
         IF(NP2.NE.NP) THEN
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9943)
            WRITE(16,9943)
 9943       FORMAT(////,' !!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',
     &               //,' THE NUMBER OF NODES IN THE BAROCLINIC',
     &                  ' INITIAL CONDITION FILE (UNIT 11) ',
     &                /,' MUST EQUAL THE NUMBER OF NODES (NP) IN ',
     &                  'THE ADCIRC GRID FILE (UNIT 14)'
     &               //,' !!!!! EXECUTION WILL NOW BE TERMINATED !!!!!')
            STOP
         ENDIF
 
         IF    (ABS(IDEN).EQ.1) THEN
            DO I=1,NP
               READ(11,*) JKI,DASigT(JKI)
            END DO
         ELSEIF(ABS(IDEN).EQ.2) THEN
            DO I=1,NP
               READ(11,*) JKI,DASalt(JKI)
            END DO
c           CALL CALC_SIGMAT_2D()      !need to activate
         ELSEIF(ABS(IDEN).EQ.3) THEN
            DO I=1,NP
               READ(11,*) JKI,DATemp(JKI)
c           CALL CALC_SIGMAT_2D()      !need to activate
            END DO
         ELSEIF(ABS(IDEN).EQ.4) THEN
            DO I=1,NP
               READ(11,*) JKI,DATemp(JKI),DASalt(JKI)
            END DO
c           CALL CALC_SIGMAT_2D()      !need to activate
         ENDIF
         CLOSE(11)
      ENDIF

C...
C...READ IN 2DDI PASSIVE SCALAR TRANSPORT INITIAL CONDITIONS
C...
      IF((C2DDI).AND.(C2D_PTrans)) THEN
         OPEN(10,FILE=TRIM(INPUTDIR)//'/'//'fort.10')
         READ(10,*)
         READ(10,*)
         READ(10,*) NP2
         IF(NP2.NE.NP) THEN
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9944)
            WRITE(16,9943)
 9944       FORMAT(////,' !!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',
     &               //,' THE NUMBER OF NODES IN THE SCALAR ',
     &                  ' INITIAL CONDITION FILE (UNIT 10) ',
     &                /,' MUST EQUAL THE NUMBER OF NODES (NP) IN ',
     &                  'THE ADCIRC GRID FILE (UNIT 14)'
     &               //,' !!!!! EXECUTION WILL NOW BE TERMINATED !!!!!')
            STOP
         ENDIF

         DO I=1,NP
            READ(10,*) JKI,DAConc(JKI)
         END DO
      ENDIF
C...
C...READ INFORMATION CONCERNING BOTTOM FRICTION COEFFICIENT
C
C     jgf46.00 If some type of spatially varying bottom friction is
C     specified in the NWP section, these inputs are ignored, and the
C     friction coefficients that are read in from the nodal attributes
C     file will take precedence.
C
C     jgf47.04 If ManningsN is loaded from the nodal attributes (fort.13)
C     file, the value of BFCdLLimit is set to CF (see nodal attributes
C     module).

      WRITE(16,1112)

      IF(NOLIBF.EQ.0) READ(15,*) TAU
      IF(NOLIBF.EQ.1) READ(15,*) CF
      IF(NOLIBF.EQ.2) READ(15,*) CF,HBREAK,FTHETA,FGAMMA
C     
      WRITE(16,2045)
 2045 FORMAT(//,' BOTTOM FRICTION INFORMATION',//)
      IF(NOLIBF.EQ.0) THEN
         WRITE(16,106) TAU
 106     FORMAT(5X,'LINEAR BOTTOM FRICTION TAU =',F12.8,5X,'1/sec'/)
         IF(TAU.NE.TAU0) THEN   !CHECK TAU VALUE AGAINST TAU0
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9951)
            WRITE(16,9951)
 9951       FORMAT(////,1X,'!!!!!!!!!!  WARNING - NONFATAL ',
     &           'INPUT ERROR  !!!!!!!!!',
     &           //,1X,'TYPICALLY YOUR INPUT VALUE FOR ',
     &           'TAU0 SHOULD BE SET EQUAL TO TAU')
         ENDIF
      ENDIF
      IF(NOLIBF.EQ.1) THEN
         WRITE(16,8) CF
 8       FORMAT(5X,'NONLINEAR FRICTION COEFFICIENT CF =',F12.8,/)
      ENDIF
      IF(NOLIBF.EQ.2) THEN
         WRITE(16,101) CF,HBREAK,FTHETA,FGAMMA
 101     FORMAT(5X,'HYBRID FRICTION RELATIONSHIP PARAMTERS, CFMIN =',
     &        F12.8,'  HBREAK = ',F8.2,
     &        /,5X,'FTHETA = ',F8.2,'  FGAMMA = ',F10.4,//)
      ENDIF
C     
C     jgf46.00 Bottom friction coefficients are initialized along with other
C     nodal attributes in the InitNodalAttr subroutine of the
C     NodalAttributes module.
C...
C...OUTPUT TO UNIT 16 GRID INFORMATION INCLUDING AGRID,NE,NP
C....H0 AND NODAL COORDINATES AND BATHYMETRY
C
C     jgf46.00 Modified this output routine so that it does not print a
C     STARTDRY column, whether the STARTDRY array has been loaded from a
C     file or not.
C
      WRITE(16,2039) AGRID
 2039 FORMAT(/,5X,'GRID IDENTIFICATION : ',A80,/)
      WRITE(16,3) NP
    3 FORMAT(5X,'TOTAL NUMBER OF NODES =',I6,/)
      WRITE(16,4) NE
    4 FORMAT(5X,'TOTAL NUMBER OF ELEMENTS =',I6,/)
      IF(ICS.EQ.2) WRITE(16,13) SLAM0*RAD2DEG,SFEA0*RAD2DEG
 13   FORMAT(5X,'LONGITUDE ABOUT WHICH CPP PROJECTION IS CENTERED',
     &     '  SLAM0 = ',F9.4,' DEGREES',
     &     /,5X,'LATITUDE  ABOUT WHICH CPP PROJECTION IS CENTERED',
     &     '  SFEA0 = ',F9.4,' DEGREES',/)
      IF(NABOUT.NE.1) THEN
         WRITE(16,24)
 24      FORMAT(/,1X,'NODAL COORDINATES AND BATHYMETRY :')
         IF(ICS.EQ.1) THEN
            IF((NTIP.EQ.0).AND.(NCOR.EQ.0)) THEN
               WRITE(16,25)
 25            FORMAT(/,10X,'NODE NO.',10X,'X',20X,'Y',15X,'DP',/)
               DO I=1,NP
                  WRITE (16,2008) I,X(I),Y(I),DP(I)
 2008             FORMAT(5X,I6,2(2X,F20.2),2X,F12.2)
               END DO
            ELSE
               WRITE(16,9195)
 9195          FORMAT(/,1X,'   NODE ',7X,'X',14X,'Y',9X,
     &              'LAMBDA(DEG)',6X,'FEA(DEG)',9X,'DP',/)
               DO I=1,NP
                  WRITE (16,9197) I,X(I),Y(I),SLAM(I)*RAD2DEG,
     &                 SFEA(I)*RAD2DEG,DP(I)
 9197             FORMAT(1X,I6,2(1X,F14.1),1X,2(1X,E15.7),1X,F8.2)
               END DO
            ENDIF
         ELSE
            WRITE(16,9225)
 9225       FORMAT(/,1X,'   NODE ',2X,'LAMBDA(DEG)',5X,'FEA(DEG)',11X,
     &           'XCP',14X,'YCP',11X,'DP',/)
            DO I=1,NP
               WRITE (16,9228) I,SLAM(I)*RAD2DEG,SFEA(I)*RAD2DEG,
     &              X(I),Y(I),DP(I)
 9228          FORMAT(1X,I6,2(1X,F14.8),2(1X,F15.1),1X,F10.2)
            END DO
         ENDIF
      ELSE
         WRITE(16,3511)
 3511    FORMAT(/,5X,'NODAL COORDINATES AND BATHYMETRY',
     &        ' INFORMATION IS AVAILABLE IN THE',
     &        /,6X,'UNIT 14 INPUT FILE')
      ENDIF
C...
C...OUTPUT TO UNIT 16 THE GLOBAL CONNECTIVITY TABLE (NODE NUMBERS FOR ELEMENTS)
C...
      IF(NABOUT.NE.1) THEN
         WRITE(16,26)
 26      FORMAT(//,5X,'GLOBAL NODE NUMBERS FOR EACH ELEMENT :')
         WRITE(16,27)
 27      FORMAT(/,9X,'ELEMENT',8X,'N1',9X,'N2',10X,'N3',/)
         DO I=1,NE
            WRITE(16,2009) I,NM(I,1),NM(I,2),NM(I,3)
 2009       FORMAT(8X,4(I7,4X))
         END DO
      ELSE
         WRITE(16,3512)
 3512    FORMAT(/,5X,'THE GLOBAL CONNECTIVITY TABLE',
     &        ' INFORMATION IS AVAILABLE IN THE',
     &        /,6X,'UNIT 14 INPUT FILE')
      ENDIF
C...
C...READ IN AND WRITE OUT EDDY VISCOSITY/DIFFUSIVITY COEFFICIENTS
C...
C     jgf46.18 Made EVM and EVC nodal attributes. Their values are
C     initialized in the call to InitNodalAttr. 
      IF (IM.EQ.10) THEN
        READ(15,*) ESLM,ESLC
        WRITE(16,111) ESLM,ESLC
 111    FORMAT(5X,'EVM, EDDY VISCOSITY COEFFICIENT =',E15.8,/,
     &         5X,'EVC, EDDY DIFFUSIVITY COEFFICIENT =',E15.8,//)
      ELSE
        READ(15,*) ESLM
        IF(ESLM.LT.0.) THEN
          CSmag_Eh= .TRUE.
          ESLM=ABS(ESLM)
          WRITE(16,1111) ESLM
 1111     FORMAT(5X,'Smagorinski lateral stress coefficient with ',
     &              'constant =',E15.8,//)           
          IF((CGWCE_Orig).OR.(CGWCE_LS_KGQ).OR.(CME_Orig)) THEN
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
              WRITE(screenunit,1112)
              WRITE(screenunit,9972)
              WRITE(screenunit,9973) 
              WRITE(screenunit,1110)
              ENDIF
            WRITE(16,1112)
            WRITE(16,9972)
            WRITE(16,9973) 
            WRITE(16,1110)
1110        FORMAT(5X,'The Smagorinski lateral stress coefficient is ',
     &                'not compatible with the original ADCIRC GWCE',//
     &             5X,'or Momentum Eqn formulations or with the ',
     &                'Kolar & Gray GWCE lateral stress formulation',//)
#ifdef CMPI
            call msg_fini()
#endif
            STOP
          ENDIF
        ELSE
          WRITE(16,11) ESLM
 11       FORMAT(5X,'Constant lateral stress coefficient =',E15.8,//)          
        ENDIF
      ENDIF
C
C     jgf46.00 Initialize nodal attributes, now that grid has been read
C     in from unit 14 file.
      CALL InitNodalAttr(DP, NP, G, NScreen, ScreenUnit, MyProc, NAbOut)
C...
C...  READ CORIOLIS INFORMATION AND COMPUTE THE CORIOLIS VECTOR
C...  OUTPUT RESULTING CORIOLIS INFORMATION
C...
      WRITE(16,1112)
      WRITE(16,2090)
 2090 FORMAT(//,1X,'CORIOLIS INFORMATION ',//)

      READ(15,*) CORI
      IF(NCOR.EQ.0) THEN
         DO I=1,NP
            CORIF(I)=CORI
         END DO
      ENDIF
      IF(NCOR.EQ.1) THEN
         DO I=1,NP
            CORIF(I)=2.0d0*7.29212d-5*SIN(SFEA(I))
         END DO
      ENDIF

      IF(NCOR.EQ.0) THEN
         WRITE(16,12) CORI
 12      FORMAT(5X,'CONSTANT CORIOLIS COEFFICIENT =',E15.8,5X,'1/SEC',/)
      ENDIF
      IF(NCOR.EQ.1) THEN
         WRITE(16,3604)
 3604    FORMAT(/,5X,'LATITUDES ARE USED TO COMPUTE VARIABLE CORIOLIS',
     &        /,7X,'AND ARE BASED ON INPUT NODAL COORDINATES',/)
         IF(NABOUT.NE.1) THEN
            WRITE(16,2092)
 2092       FORMAT(/,10X,' NODE ',5X,'NODAL CORIOLIS CORIF',/)
            DO I=1,NP
               WRITE(16,2096) I,CORIF(I)
 2096          FORMAT(7X,I6,10X,E15.9)
            END DO
         ENDIF
      ENDIF

C...
C...  READ AND PROCESS INFORMATION ABOUT THE TIDAL POTENTIAL CONSTITUENTS
C...
      READ(15,*) NTIF
      mntif = ntif
      if (ntif .eq. 0) mntif = 1

C...  allocate tidal potential arrays
      call alloc_main4a()
C...  READ TIDAL POTENTIAL AMPLITUDE, FREQUENCIES, NODAL FACTORS,
C...  EQUILIBRIUM ARGUMENTS AND ALPHANUMERIC LABEL
C....
      DO I=1,NTIF
         READ(15,'(A5)')  TIPOTAG(I)
         READ(15,*)  TPK(I),AMIGT(I),ETRF(I),FFT(I),FACET(I)
         IF(AMIGT(I).EQ.0.) THEN
            PERT(I)=0.
         ELSE
            PERT(I)=2.D0*PI/AMIGT(I)
         ENDIF
      END DO

C...  LINES TO USE EARTH LOAD/SELF-ATTRACTION PART OF TIDAL POTENTIAL FORCING

      CALL ALLOC_MAIN4b()
      IF(NTIP.EQ.2) THEN
         OPEN(24,FILE='fort.24')
         DO I=1,NTIF
            READ(24,9930)
 9930       FORMAT(///)
            DO J=1,NP
               READ(24,*) JJ,SALTAMP(I,JJ),SALTPHA(I,JJ)
               SALTPHA(I,JJ)=SALTPHA(I,JJ)*DEG2RAD
            END DO
         END DO
      ELSE
         DO I=1,NTIF
            DO J=1,NP
               SALTAMP(I,J)=0.d0
               SALTPHA(I,J)=0.d0
            END DO
         END DO
         CLOSE(24)
      ENDIF

C...
C...  OUTPUT TO UNIT 16 INFORMATION ABOUT TIDAL POTENTIAL FORCING
C...  OUTPUT WILL VARY DEPENDING ON VALUES OF NTIP,NTIF AND NCOR
C...
      WRITE(16,1112)
      WRITE(16,2102)
 2102 FORMAT(//,1X,'TIDAL POTENTIAL FORCING INFORMATION ',//)
      WRITE(16,22) NTIF
 22   FORMAT(/,1X,'TIDAL POTENTIAL IS FORCED FOR ',I5,
     &     ' CONSTITUENT(S) ')
      IF(NTIF.GT.0) WRITE(16,23)
 23   FORMAT(/,1X,'AMPLITUDE',4X,'FREQUENCY',5X,
     &     '    ETRF      ','NODAL FACTOR',2X,
     &     'EQU.ARG(DEG)',1X,'CONSTITUENT',/)
      DO I=1,NTIF
         WRITE(16,2107) TPK(I),AMIGT(I),ETRF(I),FFT(I),FACET(I),
     &        TIPOTAG(I)
 2107    FORMAT(1X,F9.7,1X,F15.12,2X,F10.7,5X,F10.7,1X,F10.3,7X,A5)
      END DO
C...
C...  CONVERT FACET(I) VALUES FROM DEGREES TO RADIANS
C...
      DO I=1,NTIF
         FACET(I)=FACET(I)*DEG2RAD
      END DO
C...
C...  CHECK CONSISTENCY OF INPUT PARAMETERS NTIF AND NTIP
C...

      IF(((NTIP.EQ.0).AND.(NTIF.NE.0)).OR.((NTIP.NE.0).AND.
     &     (NTIF.EQ.0))) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9961)
         WRITE(16,9961)
 9961    FORMAT(////,1X,'!!!!!!!!!!  WARNING - NONFATAL ',
     &        'INPUT ERROR  !!!!!!!!!',
     &        //,1X,'YOUR SELECTION OF NTIF AND NTIP (UNIT 15 INPUT ',
     &        'PARAMETERS) IS INCONSISTENT',
     &        /,1X,'PLEASE CHECK THESE VALUES')
         IF(NFOVER.EQ.1) THEN
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9987)
            WRITE(16,9987)
 9987       FORMAT(/,1X,'PROGRAM WILL OVERRIDE THE SPECIFIED ',
     &           'INPUT AND NEGLECT TIDAL POTENTIAL TERMS',
     &           /,1X,' AND/OR RESET NTIP = 0',
     &           //,1X,'!!!!!! EXECUTION WILL CONTINUE !!!!!!',//)
            NTIP=0
         ELSE
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9973)
            WRITE(16,9973)
            STOP
         ENDIF
         GOTO 1893
      ENDIF
C...  
C...  PRINT OUT LAT/LON VALUES TO BE USED IN COMPUTING TIDAL POTENTIAL
C...  IF NOT ALREADY DONE SO IN CORIOLIS SECTION AND TIDAL POTENTIAL IS
C...  ACTIVATED WITH NTIP=1
C...  
        IF(NTIP.GE.1) THEN
           IF(ICS.EQ.1) THEN
              WRITE(16,3605)
 3605         FORMAT(/,5X,'LONGITUDES AND LATITUDES ARE USED TO',
     &             ' COMPUTE THE TIDAL POTENTIAL FUNCTION',
     &             /,7X,'AND ARE BASED ON AN INVERSE CPP PROJECTION ',
     &             'OF THE INPUT COORDINATES',/)
           ELSE
              WRITE(16,2109)
 2109         FORMAT(/,5X,'LONGITUDES AND LATITUDES ARE USED TO',
     &             ' COMPUTE THE TIDAL POTENTIAL FUNCTION',
     &             /,7X,'AND ARE BASED ON INPUT NODAL COORDINATES ',/)
           ENDIF
        ENDIF
C...  
C...  INPUT FROM UNIT 15 THE TIDAL FORCING FREQUENCIES ON THE ELEVATION
C...  SPECIFIED BOUNDARIES: INCLUDING NBFR, FREQUENCIES, NODAL FACTORS,
C...  EQUILIBRIUM ARGUMENTS AND AN ELEVATION BOUNDARY CONDITION
C...  ALPHANUMERIC DESCRIPTOR
C...  
 1893   READ(15,*) NBFR
        MNBFR = NBFR
        IF (NBFR.EQ.0) MNBFR = 1

C     - Allocate arrays dimensioned by MNBFR
        call alloc_main5()

        WRITE(16,1112)
        WRITE(16,2106)
 2106   FORMAT(//,1X,'ELEVATION SPECIFIED BOUNDARY FORCING INFORMATION '
     &       ,//)
        WRITE(16,20) NBFR
 20     FORMAT(/,5X,'NUMBER OF PERIODIC, ELEVATION SPECIFIED ',
     &       'CONSTITUENTS =',I5)
        IF(NBFR.GE.1) WRITE(16,21)
 21     FORMAT(/,7X,'CONSTITUENT #',4X,'FREQUENCY',4X,'NODAL FACTOR',
     &       3X,'EQU.ARG (DEG)',2X,'CONSTITUENT',/)
        DO I=1,NBFR
           READ(15,'(A5)') BOUNTAG(I)
           READ(15,*) AMIG(I),FF(I),FACE(I)
           WRITE(16,1850) I,AMIG(I),FF(I),FACE(I),BOUNTAG(I)
 1850      FORMAT(12X,I2,6X,F16.12,2X,F10.7,2X,F10.3,10X,A5)
           FACE(I)=FACE(I)*DEG2RAD
           IF(AMIG(I).EQ.0.) THEN
              PER(I)=0.
           ELSE
              PER(I)=2.D0*PI/AMIG(I)
           ENDIF
        END DO
C...  
C...  INPUT ELEVATION BOUNDARY FORCING NODE NUMBER INFORMATION FROM
C...  UNIT 14 AND OUTPUT TO UNIT 16
C...  
C...  INPUT THE TOTAL NUMBER OF ELEVATION BOUNDARY SEGMENTS
C...  
      READ(14,*) NOPE

      WRITE(16,1852) NOPE
1852  FORMAT(///,5X,'TOTAL NUMBER OF ELEVATION BOUNDARY FORCING',
     &         ' SEGMENTS ',' = ',I5)
C...
C...INPUT THE TOTAL NUMBER OF ELEVATION BOUNDARY NODES
C...
      READ(14,*) NETA
      WRITE(16,1854) NETA
 1854 FORMAT(/,5X,'TOTAL NUMBER OF ELEVATION SPECIFIED BOUNDARY NODES ='
     &        ,I6)

C allocate arrays dimensioned by NOPE and NETA
      MNOPE = NOPE
      IF (NOPE.EQ.0) MNOPE = 1
      MNETA = NETA
      IF (NETA.EQ.0) MNETA = 1
      call alloc_main2()     
C...
C...  INPUT THE NODE NUMBERS ON EACH ELEVATION BOUNDARY FORCING SEGMENT
C...  NOTE: In the future, if a boundary type is included for each
C...  elevation segment it would be possible to distinguish between
C...  periodic and aperiodic B.C.s This would require another NBD type
C...  array.
C...
      JNMM=0
      DO K=1,NOPE
        READ(14,*) NVDLL(K)
        WRITE(16,281) K,NVDLL(K)
 281    FORMAT(//,5X,'TOTAL NUMBER OF NODES ON ELEVATION SPECIFIED ',
     &               'BOUNDARY SEGMENT ',2X,I2,2X,'=',1X,I5,/)
        DO I=1,NVDLL(K)
          READ(14,*) NBDV(K,I)
          WRITE(16,1855) NBDV(K,I)
 1855     FORMAT(7X,I7)
          NBD(JNMM+I)=NBDV(K,I)
          ENDDO
        JNMM=JNMM+NVDLL(K)
        ENDDO

C...
C...  CHECK TO MAKE SURE THAT JNMM EQUALS NETA
C...
       IF(NETA.NE.JNMM) THEN
          IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9945)
          WRITE(16,9945)
 9945     FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL INPUT ERROR ',
     &     '!!!!!!!!!',
     &    //,1X,'THE INPUT PARAMETER NETA FROM UNIT 14 DOES NOT MATCH ',
     &    'THE TOTAL NUMBER OF BOUNDARY NODES',
     &    /,1X,' FROM ALL THE SPECIFIED SEGMENTS COMBINED')

          IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9973)
          WRITE(16,9973)
          STOP
       ENDIF
      
C...
C...  SET UP TO READ IN TIME SERIES ELEVATION SPECIFIED BOUNDARY
C...  CONDITIONS IF APPROPRIATE
C...
        IF((NBFR.EQ.0).AND.(NOPE.GT.0)) THEN
           WRITE(16,1871)
 1871      FORMAT(/,5X,'TIME SERIES ELEVATION SPECIFIED VALUES ',
     &          'WILL BE READ FROM UNIT 19',
     &          /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE ',
     &          /,9X,'ELEVATION DATA WITH THE MODEL TIME STEP.')
        ENDIF

C...  INPUT FORCING CONDITIONS ON PERIODIC ELEVATION SPECIFIED
C...  BOUNDARIES FOR EACH OF THE ELEVATION FORCING FREQUENCIES FROM UNIT
C...  15 AND OUTPUT TO UNIT 16
C...  
        DO I=1,NBFR
           WRITE(16,29) I,BOUNTAG(I)
 29        FORMAT(////,5X,'ELEVATION BOUNDARY TIDAL FORCING FOR',
     &          ' CONSTITUENT NUMBER',I4,1X,'DESIGNATED : ',A5)
           READ(15,'(A10)') ALPHA
           WRITE(16,31) ALPHA
 31        FORMAT(9X,'VERIFICATION OF CONSTITUENT : ',A10,/)
           WRITE(16,30)
 30        FORMAT(14X,'NODE',11X,'AMPL.',9X,'PHASE(DEG)',/)
           DO J=1,NETA
              READ(15,*) EMO(I,J),EFA(I,J)
              WRITE(16,1870) NBD(J),EMO(I,J),EFA(I,J)
 1870         FORMAT(10X,I8,4X,F14.5,4X,F12.3)
              EFA(I,J)=EFA(I,J)*DEG2RAD
           END DO
        END DO

C.....READ THE MINIMUM INNER ANGLE FOR WHICH VELOCITY AT FLOW BOUNDARY NODES
C.....WILL BE ZEROED IN THE TANGENTIAL DIRECTIONS WHEN NORMAL FLOW IS AN
C.....ESSENTIAL B.C.

      READ(15,*) ANGINN
      WRITE(16,1112)
      WRITE(16,7654) ANGINN
7654  FORMAT(//,5X,'ANGINN = ',F8.2,' DEGREES',
     &  /,5X,'ALL FLOW BOUNDARY NODES WITH NORMAL FLOW AS AN ',
     &       'ESSENTIAL B.C. AND ',
     &  /,9X,'INNER ANGLES LESS THAN ANGINN WILL HAVE BOTH NORMAL ',
     &  /,9X,'AND TANGENTIAL VELOCITY COMPONENTS ZEROED',/)
      COSTSET=COS(ANGINN*DEG2RAD)

C... 
C... INPUT FLOW BOUNDARY INFORMATION FROM UNIT 14 AND OUTPUT TO UNIT 16
C...  

C...  INTERIOR NODES, LBCODE=-1, COS=0, SIN=1
C...  BOUNDARY NODES, LBCODE=LBCODEI=IBTYPE,
C...  COS & SIN DETERMINED FROM NORMAL DIRECTION IN ALL CASES, ALTHOUGH THIS
C...  INFORMATION IS ONLY USED WHEN NORMAL FLOW IS AN ESSENTIAL B.C. AND
C...  FREE TANGENTIAL SLIP IS ALLOWED.

C...  INPUT THE TOTAL NUMBER OF FLOW BOUNDARY SEGMENTS

      WRITE(16,1112)
      WRITE(16,1878)
1878  FORMAT(//,1X,'FLOW BOUNDARY INFORMATION ',/)
      READ(14,*) NBOU

      WRITE(16,1879) NBOU
 1879 FORMAT(//,5X,'THE TOTAL NUMBER OF FLOW BOUNDARY SEGMENTS = ',I5)

C.....INPUT THE TOTAL NUMBER OF FLOW BOUNDARY NODES

      READ(14,*) NVEL
      WRITE(16,1881) NVEL
1881  FORMAT(/,5X,'THE TOTAL NUMBER OF FLOW BOUNDARY NODES = ',I5)

      MNBOU = NBOU
      IF (NBOU.EQ.0) MNBOU = 1
      MNVEL = NVEL*2    !Cvjp  -  11/28/99 -  upper bound guess for MNVEL

C...  Allocate space for nonperiodic zero and nonzero normal flow
C...  boundary arrays including barriers

      call alloc_main3()

C.....INPUT THE NUMBER OF NODES IN THE NEXT FLOW BOUNDARY SEGMENT
C.....AND THE BOUNDARY TYPE

      JGW=0
      JME=0
      NFLUXF=0
      NFLUXB=0
      NFLUXIB=0
      NFLUXIBP=0
      NFLUXRBC=0 
      NFLUXGBC=0 
      NVELEXT=0

      DO K=1,NBOU
         READ(14,*) NVELL(K),IBTYPE
         IBTYPENC(K)=IBTYPE
C     check that IBTYPE parameter has been set properly
C     jgf46.21 Added IBTYPE=52.

         IF((IBTYPE.NE.0).AND.(IBTYPE.NE.10).AND.(IBTYPE.NE.20)
     &        .AND.(IBTYPE.NE.1).AND.(IBTYPE.NE.11).AND.(IBTYPE.NE.21)
     &        .AND.(IBTYPE.NE.2).AND.(IBTYPE.NE.12).AND.(IBTYPE.NE.22)
     &        .AND.(IBTYPE.NE.52)
     &        .AND.(IBTYPE.NE.3).AND.(IBTYPE.NE.13).AND.(IBTYPE.NE.23)
     &        .AND.(IBTYPE.NE.4).AND.(IBTYPE.NE.24)
     &        .AND.(IBTYPE.NE.5).AND.(IBTYPE.NE.25)
     &        .AND.(IBType.NE.30).AND.(IBType.NE.32)
     &        .AND.(IBType.NE.40).AND.(IBType.NE.41)               )THEN
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9985) K
            WRITE(16,9985) K
 9985       FORMAT(////,1X,'!!!!!!!!!! WARNING - FATAL ERROR !!!!!!!!!',
     &           //,1X,'THE FLOW BOUNDARY TYPE PARAMETER IBTYPE ',
     &           'HAS NOT BEEN SET CORRECTLY FOR ',
     &           /,1X,'FLOW BOUNDARY SEGMENT NO. ',I8,
     &           /,1X,'USER MUST CORRECT UNIT 14 INPUT FILE',
     &           //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED ',
     &           '!!!!!!',//)
            STOP
         ENDIF

C     write out flow boundary information to unit 16

         IF((IBTYPE.EQ.4).OR.(IBTYPE.EQ.24)) THEN
            WRITE(16,28) K,NVELL(K),K,2*NVELL(K)
 28         FORMAT(///,5X,'TOTAL NUMBER OF PAIRS FOR FLOW BOUNDARY',
     &           ' SEGMENT',2X,I2,2X,'=',2X,I5,
     &           /,5X,'TOTAL NUMBER OF NODES FOR FLOW BOUNDARY',
     &           ' SEGMENT',2X,I2,2X,'=',2X,I5)
         ELSE
            WRITE(16,128) K,NVELL(K)
 128        FORMAT(///,5X,'TOTAL NUMBER OF NODES FOR FLOW BOUNDARY',
     &           ' SEGMENT',2X,I2,2X,'=',2X,I5)
         ENDIF

C     continue processing flow boundary information

         IF(IBTYPE.EQ.0) THEN
            WRITE(16,2340)
 2340       FORMAT(5X,'THIS SEGMENT IS AN EXTERNAL BOUNDARY WITH:',/,
     &           7X,'NO NORMAL FLOW AS AN ESSENTIAL B.C.',/,
     &           7X,'AND FREE TANGENTIAL SLIP',/)
         ENDIF
         IF(IBTYPE.EQ.1) THEN
            WRITE(16,2341)
 2341       FORMAT(5X,'THIS SEGMENT IS AN INTERNAL BOUNDARY WITH:',/,
     &           7X,'NO NORMAL FLOW AS AN ESSENTIAL B.C.',/,
     &           7X,'AND FREE TANGENTIAL SLIP',/)
         ENDIF
         IF(IBTYPE.EQ.2) THEN
            NFLUXF=1
            WRITE(16,2342)
 2342       FORMAT(5X,'THIS SEGMENT IS AN EXTERNAL BOUNDARY WITH:',/,
     &           7X,'SPECIFIED NORMAL FLOW AS AN ESSENTIAL B.C.',/,
     &           7X,'AND FREE TANGENTIAL SLIP',/)
         ENDIF
         IF(IBTYPE.EQ.3) THEN
           NFLUXB=1
           WRITE(16,2344)
 2344      FORMAT(5X,'THIS SEGMENT IS AN EXTERNAL BOUNDARY WITH:',/,
     &            7X,'A BARRIER WHICH ALLOWS FREE SURFACE',
     &               ' SUPERCRITICAL OUTFLOW',/,
     &            7X,'FROM THE DOMAIN ONCE THE BARRIER HAS BEEN',
     &               ' OVERTOPPED',/,
     &            7X,'AND FREE TANGENTIAL SLIP',/)
           ENDIF
         IF(IBTYPE.EQ.4) THEN
           NFLUXIB=1
           WRITE(16,2345)
 2345      FORMAT(5X,'THIS SEGMENT IS AN INTERNAL BARRIER BOUNDARY:',/,
     &            7X,'WITH CROSS BARRIER FLOW TREATED AS AN ESSENTIAL ',
     &               ' NORMAL FLOW BOUNDARY CONDITION',/,
     &            7X,'WHICH LEAVES/ENTERS THE DOMAIN ON ONE SIDE OF ',
     &               ' THE BARRIER AND ENTERS/LEAVES THE DOMAIN ON THE',
     &          /,7X,'CORRESPONDING OPPOSITE SIDE OF THE BARRIER ',/,
     &            7X,'FLOW RATE AND DIRECTION ARE BASED ON BARRIER ',
     &               ' HEIGHT, SURFACE WATER ELEVATION',/,
     &            7X,'ON BOTH SIDES OF THE BARRIER, BARRIER COEFFICIEN',
     &               'T AND THE APPROPRIATE BARRIER FLOW FORMULA',/,
     &            7X,'FREE TANGENTIAL SLIP IS ALLOWED',/)
           ENDIF
         IF(IBTYPE.EQ.5) THEN
           NFLUXIB=1
           NFLUXIBP=1
           WRITE(16,2347)
 2347      FORMAT(5X,'THIS SEGMENT IS AN INTERNAL BARRIER BOUNDARY:',/,
     &            7X,'WITH ADDITIONAL CROSS BARRIER PIPES ',
     &               'LOCATED UNDER THE CROWN ',/,
     &            7X,'CROSS BARRIER FLOW IS TREATED AS AN ESSENTIAL',
     &               ' NORMAL FLOW BOUNDARY CONDITION',/,
     &            7X,'WHICH LEAVES/ENTERS THE DOMAIN ON ONE SIDE OF ',
     &               ' THE BARRIER AND ENTERS/LEAVES THE DOMAIN ON THE',
     &          /,7X,'CORRESPONDING OPPOSITE SIDE OF THE BARRIER ',/,
     &            7X,'FLOW RATE AND DIRECTION ARE BASED ON BARRIER ',
     &               ' HEIGHT, SURFACE WATER ELEVATION',/,
     &            7X,'ON BOTH SIDES OF THE BARRIER, BARRIER COEFFICIEN',
     &               'T AND THE APPROPRIATE BARRIER FLOW FORMULA',/,
     &            7X,'IN ADDITION CROSS BARRIER PIPE FLOW RATE AND ',
     &               ' DIRECTION ARE BASED ON PIPE CROWN HEIGHT, ',/,
     &            7X,'SURFACE WATER ELEVATION ON BOTH SIDES OF THE B',
     &               'ARRIER, PIPE FRICTION COEFFICIENT, PIPE DIAMETER',
     &          /,7X,' AND THE APPROPRIATE PIPE FLOW FORMULA',/,
     &            7X,'FREE TANGENTIAL SLIP IS ALLOWED',/)
           ENDIF
         IF(IBTYPE.EQ.10) THEN
           WRITE(16,2350)
 2350      FORMAT(5X,'THIS SEGMENT IS AN EXTERNAL BOUNDARY WITH:',/,
     &            7X,'NO NORMAL FLOW AS AN ESSENTIAL B.C.',/,
     &            7X,'AND NO TANGENTIAL SLIP',/)
           ENDIF
         IF(IBTYPE.EQ.11) THEN
           WRITE(16,2351)
 2351      FORMAT(5X,'THIS SEGMENT IS AN INTERNAL BOUNDARY WITH:',/,
     &            7X,'NO NORMAL FLOW AS AN ESSENTIAL B.C.',/,
     &            7X,'AND NO TANGENTIAL SLIP',/)
           ENDIF
         IF(IBTYPE.EQ.12) THEN
           NFLUXF=1
           WRITE(16,2352)
 2352      FORMAT(5X,'THIS SEGMENT IS AN EXTERNAL BOUNDARY WITH:',/,
     &            7X,'SPECIFIED NORMAL FLOW AS AN ESSENTIAL B.C.',/,
     &            7X,'AND NO TANGENTIAL SLIP',/)
           ENDIF
         IF(IBTYPE.EQ.13) THEN
           NFLUXB=1
           WRITE(16,2353)
 2353      FORMAT(5X,'THIS SEGMENT IS AN EXTERNAL BOUNDARY WITH:',/,
     &            7X,'A BARRIER WHICH ALLOWS FREE SURFACE',
     &               ' SUPERCRITICAL OUTFLOW',/,
     &            7X,'FROM THE DOMAIN ONCE THE BARRIER HAS BEEN',
     &               ' OVERTOPPED',/,
     &            7X,'AND NO TANGENTIAL SLIP',/)
           ENDIF
         IF(IBTYPE.EQ.20) THEN
           WRITE(16,2354)
 2354      FORMAT(5X,'THIS SEGMENT IS AN EXTERNAL BOUNDARY WITH:',/,
     &            7X,'NO NORMAL FLOW AS A NATURAL B.C.',/,
     &            7X,'AND FREE TANGENTIAL SLIP',/)
           ENDIF
         IF(IBTYPE.EQ.21) THEN
           WRITE(16,2355)
 2355      FORMAT(5X,'THIS SEGMENT IS AN INTERNAL BOUNDARY WITH:',/,
     &            7X,'NO NORMAL FLOW AS A NATURAL B.C.',/,
     &            7X,'AND FREE TANGENTIAL SLIP',/)
           ENDIF
         IF(IBTYPE.EQ.22) THEN
           NFLUXF=1
           WRITE(16,2356)
 2356      FORMAT(5X,'THIS SEGMENT IS A EXTERNAL BOUNDARY WITH:',/,
     &            7X,'SPECIFIED NORMAL FLOW AS A NATURAL B.C.',/,
     &            7X,'AND FREE TANGENTIAL SLIP',/)
           ENDIF
         IF(IBTYPE.EQ.23) THEN
           NFLUXB=1
           WRITE(16,2357)
 2357      FORMAT(5X,'THIS SEGMENT IS AN EXTERNAL BOUNDARY WITH:',/,
     &            7X,'A BARRIER WHICH ALLOWS FREE SURFACE',
     &               ' SUPERCRITICAL OUTFLOW',/,
     &            7X,'FROM THE DOMAIN ONCE THE BARRIER HAS BEEN',
     &               ' OVERTOPPED',/,
     &            7X,' IMPLEMENTED AS A NATURAL BOUNDARY CONDITION',
     &            7X,'FREE TANGENTIAL SLIP IS ALSO ALLOWED',/)
           ENDIF
         IF(IBTYPE.EQ.24) THEN
           NFLUXIB=1
           WRITE(16,2358)
 2358      FORMAT(5X,'THIS SEGMENT IS AN INTERNAL BARRIER BOUNDARY:',/,
     &            7X,'WITH CROSS BARRIER FLOW TREATED AS A NATURAL ',
     &               ' NORMAL FLOW BOUNDARY CONDITION',/,
     &            7X,'WHICH LEAVES/ENTERS THE DOMAIN ON ONE SIDE OF ',
     &               ' THE BARRIER AND ENTERS/LEAVES THE DOMAIN ON THE',
     &          /,7X,'CORRESPONDING OPPOSITE SIDE OF THE BARRIER ',/,
     &            7X,'FLOW RATE AND DIRECTION ARE BASED ON BARRIER ',
     &               ' HEIGHT, SURFACE WATER ELEVATION',/,
     &            7X,'ON BOTH SIDES OF THE BARRIER, BARRIER COEFFICIEN',
     &               'T AND THE APPROPRIATE BARRIER FLOW FORMULA',/,
     &            7X,'FREE TANGENTIAL SLIP IS ALLOWED',/)
           ENDIF
         IF(IBTYPE.EQ.25) THEN
           NFLUXIB=1
           NFLUXIBP=1
           WRITE(16,2359)
 2359      FORMAT(5X,'THIS SEGMENT IS AN INTERNAL BARRIER BOUNDARY:',/,
     &            7X,'WITH ADDITIONAL CROSS BARRIER PIPES ',
     &               'LOCATED UNDER THE CROWN ',/,
     &            7X,'CROSS BARRIER FLOW IS TREATED AS A NATURAL',
     &               ' NORMAL FLOW BOUNDARY CONDITION',/,
     &            7X,'WHICH LEAVES/ENTERS THE DOMAIN ON ONE SIDE OF ',
     &               ' THE BARRIER AND ENTERS/LEAVES THE DOMAIN ON THE',
     &          /,7X,'CORRESPONDING OPPOSITE SIDE OF THE BARRIER ',/,
     &            7X,'FLOW RATE AND DIRECTION ARE BASED ON BARRIER ',
     &               ' HEIGHT, SURFACE WATER ELEVATION',/,
     &            7X,'ON BOTH SIDES OF THE BARRIER, BARRIER COEFFICIEN',
     &               'T AND THE APPROPRIATE BARRIER FLOW FORMULA',/,
     &            7X,'IN ADDITION CROSS BARRIER PIPE FLOW RATE AND ',
     &               ' DIRECTION ARE BASED ON PIPE CROWN HEIGHT, ',/,
     &            7X,'SURFACE WATER ELEVATION ON BOTH SIDES OF THE B',
     &               'ARRIER, PIPE FRICTION COEFFICIENT, PIPE DIAMETER',
     &          /,7X,'AND THE APPROPRIATE PIPE FLOW FORMULA',/,
     &            7X,'FREE TANGENTIAL SLIP IS ALLOWED',/)
           ENDIF
         IF(IBType.EQ.30) THEN
           NFluxRBC=1
           WRITE(16,2360) 
 2360      FORMAT(5X,'This segment is an outward radiating boundary:',/,
     &            7X,'the GWCE is forced with normal flux.  A wave ',
     &               'radiation condition',
     &            7X,'is used to relate the time derivative of the ',
     &               'normal flux to',/,
     &            7X,'the time derivative of the elevation.  The ',
     &               'momentum equations',/,
     &            7X,'are used to compute the velocity field the same ',
     &               'as for a nonboundary node.',/)
           ENDIF
         IF(IBType.EQ.32) THEN
           NFluxRBC=1
           NFluxF=1
           WRITE(16,2356)
 2362      FORMAT(5X,'This segment is a combined specified normal ',
     &               'flux and outward radiating boundary.',/,
     &            7X,'The GWCE is forced with the total normal flux ',
     &               'computed by adding the specified',/,
     &            7X,'normal flux and the flux associated with the ',
     &               'outward radiating wave.',/,
     &            7X,'The latter is determine from a Sommerfeld type ',
     &               'condition, flux=celerity*wave elevation',/,
     &            7X,'The momentum equations are used to compute the',/,
     &            7X,'velocity field the same as for a nonboundary ',
     &               'node',/)
           ENDIF
         IF(IBType.EQ.40) THEN
           NFluxGBC=1
           WRITE(16,2370)
 2370      FORMAT(5X,'This segment is a zero normal velocity gradient ', 
     &               'boundary:',/,
     &            7X,'the GWCE is forced with normal flux',/,
     &            7X,'the momentum eqs are sacrificed in favor of ',  
     &               'setting the velocity at a',/,
     &            7X,'boundary node equal to the value at a fictitious',
     &               ' point inside the domain.',/,
     &            7X,'The fictitious point lies on the inward directed',
     &               ' normal to the boundary',/,
     &            7X,'a distance equal to the distance from the ',
     &               'boundary node to its farthest',/,
     &            7X,'neighbor. This should ensure that the fictitious',
     &               ' point does not',/,
     &            7X,'fall into an element that contains the boundary ',
     &               'node.',/,
     &            7X,'The velocity at the fictitious point is ',
     &               'determined by interpolation.',/)
           ENDIF
         IF(IBType.EQ.41) THEN
           NFluxGBC=1
           WRITE(16,2371)
 2371      FORMAT(5X,'This segment is a zero normal velocity gradient ', 
     &               'boundary:',/,
     &            7X,'the GWCE is forced with normal flux',/,
     &            7X,'the momentum eqs are sacrificed in favor of eqs ',
     &               'that set the velocity',/,
     &            7X,'gradient normal to the boundary equal to zero ',
     &               'in the Galerkin sense.',/)
        ENDIF
        IF(IBTYPE.EQ.52) THEN
           NFLUXF=1
           IBTYPE52=.TRUE.
           WRITE(16,6532)
 6532      FORMAT(5X,'This segment uses a specified periodic normal ',
     &          'flux as a natural b.c. combined with outward ',
     &          7X,'radiating boundary. The GWCE is forced with the ,'
     &          'total normal flux computed by adding the',/,
     &          7X,'specified normal flux and the flux associated ',
     &          'with the outward radiating wave.',/,
     &          7X,'The latter is determined from ',
     &          'flux=celerity*wave elevation',/,
     &          7X,'The momentum equations are used to compute the',/,
     &          7X,'velocity field the same as for a nonboundary ',
     &          'node.',/)
        ENDIF

C...  For NON BARRIER type boundaries, read the node numbers for the Kth
C...  boundary segment

         IF((IBTYPE.NE. 3).AND.(IBTYPE.NE.13).AND.(IBTYPE.NE.23).AND.
     &      (IBTYPE.NE. 4).AND.(IBTYPE.NE.24).AND.
     &      (IBTYPE.NE. 5).AND.(IBTYPE.NE.25)) THEN
           DO I=1,NVELL(K)
             READ(14,*) NBVV(K,I)
             END DO
           NPRBI=1
           NPIPE=0
           ENDIF

C...  For EXTERNAL BARRIER type boundaries, read the node number, the
C...  barrier elevation above the Geoid and the coefficient of free
C...  surface super-critical flow for the nodes in the Kth boundary
C...  segment.

         IF((IBTYPE.EQ.3).OR.(IBTYPE.EQ.13).OR.(IBTYPE.EQ.23)) THEN
           DO I=1,NVELL(K)
             READ(14,*) NBVV(K,I),BARLANHTR(I),BARLANCFSPR(I)
             END DO
           NPRBI=1
           NPIPE=0
           ENDIF

C...  For INTERNAL BARRIER type boundaries, read the node number, the
C...  paired node number, the barrier elevation above the Geoid and the
C...  coefficients of supercritical and subcritical flow for the nodes
C...  in the Kth boundary segment.

         IF((IBTYPE.EQ.4).OR.(IBTYPE.EQ.24)) THEN
            DO I=1,NVELL(K)
               READ(14,*) NBVV(K,I),IBCONNR(I),BARINHTR(I),BARINCFSBR(I)
     &              ,BARINCFSPR(I)
            END DO
            NPRBI=2
            NPIPE=0
         ENDIF

C...  For INTERNAL BARRIER WITH CROSS BARRIER PIPES type boundaries,
C...  read the node number, the paired node number, the barrier
C...  elevation above the Geoid, the coefficients of supercritical and
C...  subcritical flow, the cross barrier pipe height, the cross barrier
C...  pipe coefficient and the cross barrier pipe diameter for the nodes
C...  in the Kth boundary segment.

         IF((IBTYPE.EQ.5).OR.(IBTYPE.EQ.25)) THEN
            DO I=1,NVELL(K)
         READ(14,*) NBVV(K,I),IBCONNR(I),BARINHTR(I),BARINCFSBR(I),
     &             BARINCFSPR(I),PIPEHTR(I),PIPECOEFR(I),
     &              PIPEDIAMR(I)
            END DO
            NPRBI=2
            NPIPE=1
         ENDIF
         
C...  PROCESS INFORMATION FOR VARIOUS TYPES OF FLOW BOUNDARY SEGMENTS

         DO IPRBI=1,NPRBI
            
C...  LOAD PAIRED NODES INTO PRIMARY PROCESSING VECTORS AND RESET
C...  CONNECTING NODES FOR BACK FACE THUS BACK/CONNECTING NODES ARE
C...  BEING LOADED AS PRIMARY NODES AND FRONT NODES ARE RELOADED AS
C...  CONNECTING NODES NOTE THAT THE CLOCKWISE ORIENTATION OF ISLAND
C...  TYPE BOUNDARIES IS BEING MAINTAINED WHEN BACK NODES ARE RELOADED
C...  AS PRIMARY NODES ADDITIONAL INTERNAL BARRIER BOUNDARY INFORMATION
C...  IS ALSO RESET

            IF(IPRBI.EQ.2) THEN
               DO I=1,NVELL(K)
                  NTRAN1(I)=NBVV(K,I)
                  NTRAN2(I)=IBCONNR(I)
                  BTRAN3(I)=BARINHTR(I)
                  BTRAN4(I)=BARINCFSBR(I)
                  BTRAN5(I)=BARINCFSPR(I)
                  IF(NPIPE.EQ.1) THEN
                     BTRAN6(I)=PIPEHTR(I)
                     BTRAN7(I)=PIPECOEFR(I)
                     BTRAN8(I)=PIPEDIAMR(I)
                  ENDIF
               END DO
               DO I=1,NVELL(K)
                  NBVV(K,I)=NTRAN2(NVELL(K)+1-I)
                  IBCONNR(I)=NTRAN1(NVELL(K)+1-I)
                  BARINHTR(I)=BTRAN3(NVELL(K)+1-I)
                  BARINCFSBR(I)=BTRAN4(NVELL(K)+1-I)
                  BARINCFSPR(I)=BTRAN5(NVELL(K)+1-I)
                  IF(NPIPE.EQ.1) THEN
                     PIPEHTR(I)=BTRAN6(NVELL(K)+1-I)
                     PIPECOEFR(I)=BTRAN7(NVELL(K)+1-I)
                     PIPEDIAMR(I)=BTRAN8(NVELL(K)+1-I)
                  ENDIF
               END DO
            ENDIF

C...  WRITE OUT ADDITIONAL HEADER FOR INTERNAL BARRIER BOUNDARIES

            IF((IBTYPE.EQ.4).OR.(IBTYPE.EQ.24)) THEN
               IF(IPRBI.EQ.1) THEN
                  WRITE(16,1842)
 1842          FORMAT(/,5X,'FRONT FACE OF INTERNAL BARRIER BOUNDARY',/)
               ELSE
                  WRITE(16,1843)
 1843          FORMAT(/,5X,'BACK FACE OF INTERNAL BARRIER BOUNDARY',/)
               ENDIF
            ENDIF

C...  WRITE OUT ADDITIONAL HEADER FOR INTERNAL BARRIER BOUNDARIES WITH
C...  CROSS BARRIER PIPES

            IF((IBTYPE.EQ.5).OR.(IBTYPE.EQ.25)) THEN
               IF(IPRBI.EQ.1) THEN
                  WRITE(16,1844)
 1844             FORMAT(/,5X,'FRONT FACE OF INTERNAL BARRIER BOUNDARY',
     &                 ' WITH CROSS BARRIER PIPES',/)
               ELSE
                  WRITE(16,1845)
 1845             FORMAT(/,5X,'BACK FACE OF INTERNAL BARRIER BOUNDARY',
     &                 ' WITH CROSS BARRIER PIPES',/)
               ENDIF
            ENDIF

C...  WRITE OUT GENERAL HEADER FOR BOUNDARY INFORMATION

            WRITE(16,1841)
 1841 FORMAT('    JGW    JME    ME2GW   NODE #  BNDRY CODE   INNER',
     &      ' ANGLE',7X,'COS',13X,'SIN',9X,'0.667*BNDRY LEN',/)

C...  COMPLETE THE BOUNDARY ARRAY FOR THE Kth FLOW BOUNDARY SEGMENT

            NBVV(K,0)=NBVV(K,1) !UNCLOSED EXTERNAL
            IF((IBTYPE.EQ.1).OR.(IBTYPE.EQ.11).OR.(IBTYPE.EQ.21)) THEN
               IF(NBVV(K,NVELL(K)).NE.NBVV(K,1)) THEN !CLOSE AN UNCLOSED INTERNAL
                  NVELL(K)=NVELL(K)+1
                  NBVV(K,NVELL(K))=NBVV(K,1)
               ENDIF
            ENDIF
            IF(NBVV(K,NVELL(K)).EQ.NBVV(K,1)) THEN !CLOSED EXTERNAL OR INTERNAL
               NBVV(K,0)=NBVV(K,NVELL(K)-1)
            ENDIF
            NBVV(K,NVELL(K)+1)=NBVV(K,NVELL(K))

C...  PUT BOUNDARY INFORMATION INTO 2 TYPES OF ARRAYS, ONE FOR THE GWCE
C...  B.C. AND ONE FOR THE MOMENTUM EQUATION B.C.
c
C...  THE GWCE ARRAYS INCLUDE EVERY NODE IN THE UNIT 14 FILE, I.E.,
C...  NODES ARE REPEATED WHERE SPECIFIED NORMAL FLOW AND NO NORMAL FLOW
C...  BOUNDARIES MEET AND AT THE BEGINNING AND END OF CLOSED EXTERNAL
C...  BOUNDARIES AND ISLANDS.
c
C...  THE MOMENTUM EQUATION ARRAYS ARE KEYED TO THE GWCE ARRAYS VIA THE
C...  ARRAY ME2GW WHICH INDICATES THE LOCATION IN THE GWCE ARRAYS THAT
C...  THE APPROPRIATE M.E. VALUE LIES.
c
C...  THE M.E. ARRAYS DO NOT REPEAT NODES THAT ARE DUPLICATED IN THE
C...  UNIT 14 FILE, I.E., WHEN SPECIFIED NORMAL FLOW AND NO NORMAL FLOW
C...  BOUNDARIES MEET, THE SPECIFIED NORMAL FLOW BOUNDARY CONDITION
C...  TAKES PRECEDENT.  ALSO THE BEGINNING AND ENDING NODES OF CLOSED
C...  EXTERNAL AND ISLAND BOUNDARIES ARE NOT REPEATED.

            DO I=1,NVELL(K)

C...  SET UP THE GWCE BOUNDARY ARRAYS WHICH CONSIST OF
C...  BOUNDARY NODE NUMBERS
C...  BOUNDARY CODES
C...  0.66667*LENGTH OF EACH BOUNDARY SEGMENT.  NOTE, THE LENGTH OF THE
C...  LAST BOUNDARY SEGMENT ON EACH BOUNDARY SHOULD BE ZERO
C     jgf46.21 Added IBTYPE=52.
               JGW=JGW+1
               IF((IBTYPE.EQ.0).OR.(IBTYPE.EQ.10).OR.(IBTYPE.EQ.20).OR.
     &             (IBTYPE.EQ.2).OR.(IBTYPE.EQ.12).OR.(IBTYPE.EQ.22).OR.
     &             (IBTYPE.EQ.52).OR.
     &             (IBTYPE.EQ.3).OR.(IBTYPE.EQ.13).OR.(IBTYPE.EQ.23).OR.
     &             (IBType.EQ.30).OR.(IBType.EQ.32).OR.
     &             (IBType.EQ.40).OR.(IBType.EQ.41)) THEN
                  NVELEXT=NVELEXT+1
               ENDIF
               NBV(JGW)=NBVV(K,I)
               NBVI=NBVV(K,I)
               NBVJ=NBVV(K,I+1)
               DELX=X(NBVJ)-X(NBVI)
               DELY=Y(NBVJ)-Y(NBVI)
               BNDLEN2O3(JGW)=2.D0*(SQRT(DELX*DELX+DELY*DELY))/3.D0

C...  COMPUTE THE INCLUDED ANGLE AND TEST TO DETERMINE WHETHER TO ZERO
C...  TANGENTIAL VELOCITIES

C...  NOTE:.IMPLEMENTATION FOR ICS=2 REQUIRES COMPUTING ALL COORDINATES
C...  IN A LOCALIZED SYSTEM (I.E. THE TRANSFORMATION IS CENTERED AT
C...  X0,Y0)

               IF(ICS.EQ.1) THEN
                  XL0=X(NBVV(K,I))
                  XL1=X(NBVV(K,I-1))
                  XL2=X(NBVV(K,I+1))
                  YL0=Y(NBVV(K,I))
                  YL1=Y(NBVV(K,I-1))
                  YL2=Y(NBVV(K,I+1))
               ELSE
                  CALL CPP(XL0,YL0,SLAM(NBVV(K,I)),SFEA(NBVV(K,I)),
     &                 SLAM(NBVV(K,I)),SFEA(NBVV(K,I)))
                  CALL CPP(XL1,YL1,SLAM(NBVV(K,I-1)),SFEA(NBVV(K,I-1)),
     &                 SLAM(NBVV(K,I)),SFEA(NBVV(K,I)))
                  CALL CPP(XL2,YL2,SLAM(NBVV(K,I+1)),SFEA(NBVV(K,I+1)),
     &                 SLAM(NBVV(K,I)),SFEA(NBVV(K,I)))
               ENDIF

C...  NOTE: INTERIOR ANGLE AT ENDS OF BOUNDARIES MUST BE EQUAL, EITHER:
C...  A FICTITIOUSLY LARGE VALUE IF THE BOUNDARY IS NOT CLOSED OR A TRUE
C...  VALUE IF THE BOUNDARY IS CLOSED

               THETA=0.d0
               IF((I.EQ.1).AND.(NBVV(K,I).EQ.NBVV(K,I-1))) THEN
                  THETA1=-9999999.d0
                  THETA=THETA1
                  COSTHETA1=COSTSET
                  COSTHETA=COSTHETA1
                  CROSS1=0.d0
                  CROSS=CROSS1
               ENDIF
               IF(I.EQ.NVELL(K)) THEN
                  THETA=THETA1
                  COSTHETA=COSTHETA1
                  CROSS=CROSS1
               ENDIF
               IF(THETA.EQ.0.) THEN
                  VL1X=XL1-XL0
                  VL1Y=YL1-YL0
                  VL2X=XL2-XL0
                  VL2Y=YL2-YL0
                  DOTVEC=VL1X*VL2X+VL1Y*VL2Y
               VECNORM=(SQRT(VL1X**2+VL1Y**2))*(SQRT(VL2X**2+VL2Y**2))
                  COSTHETA=DOTVEC/VECNORM
                  IF(COSTHETA.GT.1.0d0) COSTHETA=1.0d0
                  IF(COSTHETA.LT.-1.0d0) COSTHETA=-1.0d0
                  THETA=RAD2DEG*ACOS(COSTHETA)
                  CROSS=-VL1X*VL2Y+VL2X*VL1Y
                  IF(CROSS.LT.0) THETA=360.d0-THETA
                  IF(I.EQ.1) THEN
                     THETA1=THETA
                     COSTHETA1=COSTHETA
                     CROSS1=CROSS
                  ENDIF
               ENDIF

C...  CHECK WHETHER ANGLE IS LESS THAN MINIMUM ANGLE, IF SO CHANGE THE
C...  BOUNDARY CODE TO ZERO TANGENTIAL VELOCITIES

               LBCODEI(JGW)=IBTYPE
               IF((COSTHETA.GT.COSTSET).AND.(CROSS.GT.0.0)) THEN
                  IF(IBTYPE.EQ.0) LBCODEI(JGW)=10
                  IF(IBTYPE.EQ.1) LBCODEI(JGW)=11
                  IF(IBTYPE.EQ.2) LBCODEI(JGW)=12
                  IF(IBTYPE.EQ.3) LBCODEI(JGW)=13
                  IF((IBTYPE.GE.0).AND.(IBTYPE.LE.3)) THEN
                     WRITE(16,1856) NBVV(K,I),THETA
 1856                FORMAT(2X,I7,4X,'THE INNER ANGLE = ',F8.2,1X,
     &                    'TANGENTIAL SLIP WILL BE ZEROED')
                  ENDIF
               ENDIF

C...  COMPUTE COS AND SIN OF OUTWARD NORMAL REGARDLESS OF BOUNDARY TYPE

               X1=X(NBVV(K,I-1))
               X2=X(NBVV(K,I+1))
               Y1=Y(NBVV(K,I-1))
               Y2=Y(NBVV(K,I+1))
               XL=SQRT((X1-X2)**2+(Y1-Y2)**2)
               CSII(JGW)=SFAC(NBVV(K,I))*(Y2-Y1)/XL
               SIII(JGW)=(X1-X2)/XL

C...  Compute the location of and the element containing the fictitious
C...  node used for a zero normal velocity gradient boundary condition
C...  (type 40)

               IF(IBType.EQ.40) Then
                  ZNGFicNodeDist=0.d0
                  DO N=2,NNeigh(NBVI)
                  ZNGFicNodeDistTemp=SQRT((X(NeiTab(NBVI,N))-X(NBVI))**2
     &                    +(Y(NeiTab(NBVI,N))-Y(NBVI))**2)
                     IF(ZNGFicNodeDist.LT.ZNGFicNodeDistTemp)
     &                    ZNGFicNodeDist=ZNGFicNodeDistTemp
                  ENDDO
                  ZNGFicNodeDist=1.001d0*ZNGFicNodeDist
                  X4=X(NBVI)-ZNGFicNodeDist*CSII(JGW)
                  Y4=Y(NBVI)-ZNGFicNodeDist*SIII(JGW)   
                  AEMIN=1.0E+25
                  DO N=1,NE
                     N1=NM(N,1)
                     N2=NM(N,2)
                     N3=NM(N,3)
                     X1=X(N1)
                     X2=X(N2)
                     X3=X(N3)
                     Y1=Y(N1)
                     Y2=Y(N2)
                     Y3=Y(N3)
                     A1=(X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4)
                     A2=(X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1)
                     A3=(Y4-Y1)*(X2-X1)-(X4-X1)*(Y2-Y1)
                     AA=ABS(A1)+ABS(A2)+ABS(A3)
                     AE=ABS(AA-AREAS(N))/AREAS(N)
                     IF(AE.LT.AEMIN) THEN
                        AEMIN=AE
                        NEleZNG(JGW)=N
                     ENDIF
                  ENDDO
                  
                  IF(AEMIN.GT.1.0E-5) THEN
                     if (NSCREEN.ne.0.and.myproc.eq.0)WRITE(ScreenUnit,9770) NBVI
                     WRITE(16,9770) NBVI
 9770                FORMAT(///,1X,'!!!!!!!!!!  ERROR !!!!!!!!!',//,
     &                1X,'Zero Normal Velocity Gradient Boundary Node ',
     &                   'Number ',I10,
     &                3X,'does not appear to have a fictitious node ',
     &                   'located within the domain.',
     &                3X,'This should be checked.')
                  ENDIF

C...  Pre-compute the information required to interpolate at zero normal
C...  gradient fictitious nodes

                  N1=NM(NEleZNG(JGW),1)
                  N2=NM(NEleZNG(JGW),2)
                  N3=NM(NEleZNG(JGW),3)
                  X1=X(N1)
                  X2=X(N2)
                  X3=X(N3)
                  Y1=Y(N1)
                  Y2=Y(N2)
                  Y3=Y(N3)
                  ZNGIF1(JGW)=((X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4))
     &                 /AREAS(NEleZNG(JGW))
                  ZNGIF2(JGW)=((X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1))
     &                 /AREAS(NEleZNG(JGW))
                  ZNGIF3(JGW)=(-(X4-X1)*(Y2-Y1)+(Y4-Y1)*(X2-X1))
     &                 /AREAS(NEleZNG(JGW))
c             WRITE(2,*) 'For Boundary Node # ',NBVV(K,I)
c             WRITE(2,*) 'The fictitious point is ',X4,Y4
c             WRITE(2,*) 'which is located in element # ',NEleZNG(JGW)
c             WRITE(2,*) 'The interpolating factors are ',ZNGIF1(JGW),
c    &                                           ZNGIF2(JGW),ZNGIF3(JGW)
c             WRITE(2,*) ' '
               ENDIF

C...
C...  SET UP THE MOMENTUM EQUATION BOUNDARY ARRAY WHICH CONSISTS OF A
C...  KEY TO THE GWCE BOUNDARY CONDITION ARRAY
C...
               IF(I.EQ.1) THEN  !DEAL WITH FIRST NODE IN L.B. SEG
                  IF(JGW.EQ.1) THEN !VERY FIRST L.B. SEG
                     JME=JME+1  !M.E. USES IT
                     ME2GW(JME)=JGW
                  ENDIF
                  IF(JGW.NE.1) THEN
                     IF(NBV(JGW).NE.NBV(JGW-1)) THEN !L.B. SEGS DON'T OVERLAP
                        JME=JME+1 !M.E. USES IT
                        ME2GW(JME)=JGW
                     ENDIF
                     IF(NBV(JGW).EQ.NBV(JGW-1)) THEN !L.B. SEGS OVERLAP
                        IF((LBCODEI(JGW).EQ.2) .OR. ! M.E. USES IT ONLY
     &                       (LBCODEI(JGW).EQ.12).OR. ! IF IT IS
     &                       (LBCODEI(JGW).EQ.22).OR. ! SPECIFIED FLOW,
     &                       (LBCODEI(JGW).EQ.3) .OR. ! AN OVERFLOW BARRIER
     &                       (LBCODEI(JGW).EQ.13).OR. ! OR A RADIATION
     &                       (LBCODEI(JGW).EQ.23).OR. ! BOUNDARY
     &                       (LBCODEI(JGW).EQ.30).OR.
     &                       (LBCODEI(JGW).EQ.32).OR. 
     &                       (LBCODEI(JGW).EQ.40).OR. 
     &                       (LBCODEI(JGW).EQ.41).OR.
     &                       (LBCODEI(JGW).EQ.52)) ME2GW(JME)=JGW
                     ENDIF
                  ENDIF
               ENDIF
               IF((I.GT.1).AND.(I.LT.NVELL(K))) THEN !IF NOT FIRST OR
                  JME=JME+1     !LAST NODE
                  ME2GW(JME)=JGW !M.E. USES IT
               ENDIF
               IF(I.EQ.NVELL(K)) THEN !DEAL WITH LAST NODE ON BOUNDARY
                  IF((NBV(JGW).NE.NBVV(K,1)).AND. !IF UNCLOSED BOUNDARY
     &                 (NBV(JGW).NE.NBV(1))) THEN !M.E. USES IT
                     JME=JME+1
                     ME2GW(JME)=JGW
                  ENDIF
                  IF(NBVV(K,I).EQ.NBV(1)) THEN !IF OVERLAPS WITH VERY FIRST
                     IF((LBCODEI(JGW).EQ.2) .OR. ! L.B. NODE
     &                    (LBCODEI(JGW).EQ.12).OR. ! M.E. USES IT ONLY IF IT IS
     &                    (LBCODEI(JGW).EQ.22).OR. ! SPECIFIED FLOW,
     &                    (LBCODEI(JGW).EQ.3) .OR. ! AN OVERFLOW BARRIER OR
     &                    (LBCODEI(JGW).EQ.13).OR. ! A RADIATION
     &                    (LBCODEI(JGW).EQ.23).OR. ! BOUNDARY
     &                    (LBCODEI(JGW).EQ.30).OR. 
     &                    (LBCODEI(JGW).EQ.32).OR.
     &                    (LBCODEI(JGW).EQ.40).OR.
     &                    (LBCODEI(JGW).EQ.41).OR.
     &                    (LBCODEI(JGW).EQ.52)) ME2GW(1)=JGW
                  ENDIF
               ENDIF

C...........LOAD EXTERNAL BARRIER BOUNDARY INFORMATION INTO THE CORRECT VECTORS
               IF((IBTYPE.EQ.3).OR.(IBTYPE.EQ.13)
     &          .OR.(IBTYPE.EQ.23)) THEN
                  BARLANHT(JGW)=BARLANHTR(I)
                  BARLANCFSP(JGW)=BARLANCFSPR(I)
               ENDIF

C...........LOAD INTERNAL BARRIER BOUNDARY INFORMATION INTO THE CORRECT VECTORS
               IF((IBTYPE.EQ.4).OR.(IBTYPE.EQ.24)) THEN
                  IBCONN(JGW)=IBCONNR(I)
                  BARINHT(JGW)=BARINHTR(I)
                  BARINCFSB(JGW)=BARINCFSBR(I)
                  BARINCFSP(JGW)=BARINCFSPR(I)
               ENDIF
               
C...........LOAD INTERNAL BARRIER WITH PIPES BOUNDARY INFORMATION INTO 
C............THE CORRECT VECTORS
               IF((IBTYPE.EQ.5).OR.(IBTYPE.EQ.25)) THEN
                  IBCONN(JGW)=IBCONNR(I)
                  BARINHT(JGW)=BARINHTR(I)
                  BARINCFSB(JGW)=BARINCFSBR(I)
                  BARINCFSP(JGW)=BARINCFSPR(I)
                  PIPEHT(JGW)=PIPEHTR(I)
                  PIPECOEF(JGW)=PIPECOEFR(I)
                  PIPEDIAM(JGW)=PIPEDIAMR(I)
               ENDIF

C...........WRITE OUT BOUNDARY CONDITION ARRAY INFORMATION

               WRITE(16,1857) JGW,JME,ME2GW(JME),NBV(JGW),LBCODEI(JGW),
     &              THETA,CSII(JGW),SIII(JGW),BNDLEN2O3(JGW)
 1857          FORMAT(1X,I6,1X,I6,1X,I6,3X,I6,3X,I4,9X,F8.2,2X,E16.8,1X,
     &              E16.8,2X,E16.8)

C...........CHECK EXTERNAL BARRIER HEIGHTS AGAINST DEPTHS
               IF((IBTYPE.EQ.3).OR.(IBTYPE.EQ.13)
     &          .OR.(IBTYPE.EQ.23)) THEN
                  IF(BARLANHT(JGW).LT.-DP(NBV(JGW))) THEN
                     IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
                        WRITE(ScreenUnit,8367) 
     &                    JGW,NBV(JGW),BARLANHT(JGW),DP(NBV(JGW))
                        ENDIF
                     WRITE(16,8367) JGW,NBV(JGW),BARLANHT(JGW),
     &                DP(NBV(JGW))
 8367            FORMAT(////,1X,'!!!!!!!!!!  FATAL INPUT ERROR   !!!'
     &           ,'!!!!!!!!!',//,
     &           1X,'AT BOUNDARY NODE NO.',I6,' (GLOBAL NODE NO.',
     &           I6, ' AND OF EXTERNAL BARRIER TYPE) ',/,
     &           2X,'THE EXTERNAL BARRIER HEIGHT = ',E12.5,
     &           2X,'IS EXCEEDED BY THE DEPTH SPECIFIED AT ',/,2X
     &           ,'THE ASSOCIATED GLOBAL NODE = ',E12.5,/,2X,
     &           'USER MUST SPECIFY CONSISTENT BARRIER HEIGHTS',
     &           ' AND DEPTHS')
                     IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
                        WRITE(ScreenUnit,9973)
                     ENDIF
                     WRITE(16,9973)
#ifdef CMPI
                     call msg_fini()
#endif
                     STOP
                  ENDIF
               ENDIF
C...........CHECK INTERNAL BARRIER HEIGHTS AGAINST DEPTHS
               IF((IBTYPE.EQ.4).OR.(IBTYPE.EQ.24)) THEN
                  IF(BARINHT(JGW).LT.-DP(NBV(JGW))) THEN
                     IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
                        WRITE(ScreenUnit,8368) 
     &                    JGW,NBV(JGW),BARINHT(JGW),DP(NBV(JGW))
                     ENDIF
                     WRITE(16,8368) JGW,NBV(JGW),BARINHT(JGW),
     &                    DP(NBV(JGW))
 8368            FORMAT(////,1X,'!!!!!!!!!!  FATAL INPUT ERROR   !!!'
     &           ,'!!!!!!!!!',//,
     &           1X,'AT BOUNDARY NODE NO.',I6,' (GLOBAL NODE NO. ',
     &           I6,' AND OF INTERNAL BARRIER TYPE) ',/,
     &           2X,'THE INTERNAL BARRIER HEIGHT = ',E12.5,
     &           2X,'IS EXCEEDED BY THE DEPTH SPECIFIED AT ',/,2X
     &           ,'THE ASSOCIATED GLOBAL NODE = ',E12.5,/,2X,
     &           'USER MUST SPECIFY CONSISTENT BARRIER HEIGHTS',
     &           ' AND DEPTHS')
                     IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
                        WRITE(ScreenUnit,9973)
                     ENDIF
                     WRITE(16,9973)
                     STOP
                  ENDIF
               ENDIF
                        
C...........CHECK INTERNAL BARRIER WITH PIPES BARRIER HEIGHTS AGAINST DEPTHS
               IF((IBTYPE.EQ.5).OR.(IBTYPE.EQ.25)) THEN
                  IF(BARINHT(JGW).LT.-DP(NBV(JGW))) THEN
                     IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
                        WRITE(ScreenUnit,8370) 
     &                    JGW,NBV(JGW),BARINHT(JGW),DP(NBV(JGW))
                     ENDIF
                     WRITE(16,8370) JGW,NBV(JGW),BARINHT(JGW),
     &                    DP(NBV(JGW))
 8370            FORMAT(////,1X,'!!!!!!!!!!  FATAL INPUT ERROR   !!!'
     &           ,'!!!!!!!!!',//,
     &           1X,'AT BOUNDARY NODE NO.',I6,' (GLOBAL NODE NO. ',
     &           I6,' AND OF INTERNAL BARRIER TYPE) ',/,
     &           2X,'THE INTERNAL BARRIER HEIGHT = ',E12.5,
     &           2X,'IS EXCEEDED BY THE DEPTH SPECIFIED AT ',/,2X
     &           ,'THE ASSOCIATED GLOBAL NODE = ',E12.5,/,2X,
     &           'USER MUST SPECIFY CONSISTENT BARRIER HEIGHTS',
     &           ' AND DEPTHS')
                     IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
                        WRITE(ScreenUnit,9973)
                     ENDIF
                     WRITE(16,9973)
#ifdef CMPI
                     call msg_fini()
#endif
                     STOP
                  ENDIF
               ENDIF
C...  CHECK INTERNAL BARRIER WITH PIPES PIPE HEIGHTS AGAINST DEPTHS
               IF((IBTYPE.EQ.5).OR.(IBTYPE.EQ.25)) THEN
                  IF(PIPEHT(JGW).LT.-DP(NBV(JGW))) THEN
                     IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
                        WRITE(ScreenUnit,8372) 
     &                    JGW,NBV(JGW),BARINHT(JGW),DP(NBV(JGW))
                     ENDIF
                     WRITE(16,8372) JGW,NBV(JGW),BARINHT(JGW),
     &                    DP(NBV(JGW))
 8372            FORMAT(////,1X,'!!!!!!!!!!  FATAL INPUT ERROR   !!!'
     &           ,'!!!!!!!!!',//,
     &           1X,'AT BOUNDARY NODE NO.',I6,' (GLOBAL NODE NO. ',
     &           I6,' AND OF INTERNAL BARRIER TYPE) ',/,
     &           2X,'THE BARRIER PIPE HEIGHT = ',E12.5,
     &           2X,'IS EXCEEDED BY THE DEPTH SPECIFIED AT ',/,2X
     &           ,'THE ASSOCIATED GLOBAL NODE = ',E12.5,/,2X,
     &           'USER MUST SPECIFY CONSISTENT PIPE HEIGHTS',
     &           ' AND DEPTHS')
                     IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
                        WRITE(ScreenUnit,9973)
                     ENDIF
                     WRITE(16,9973)
#ifdef CMPI
                     call msg_fini()
#endif
                     STOP
                  ENDIF
               ENDIF

C...  CHECK FOR OVERLAPPING OF AN INTERNAL BARRIER BOUNDARY WITH
C...  ANY EXTERNAL BARRIER BOUNDARY. IF THIS DOES OCCUR, TAKE
C...  APPROPRIATE ACTION
               IF((IBTYPE.EQ.4).OR.(IBTYPE.EQ.24).OR.(IBTYPE.EQ.5)
     &              .OR.(IBTYPE.EQ.25)) THEN
                  DO ICK=1,NVELEXT
C...  CHECK IF OVERLAP EXISTS
                     IF(NBV(ICK).EQ.NBV(JGW)) THEN
C...  CHECK FOR ILLEGAL OVERLAPS
                        IF((LBCODEI(ICK).EQ.2).OR.
     &                   (LBCODEI(ICK).EQ.3).OR.(LBCODEI(ICK).EQ.12).OR.
     &                   (LBCODEI(ICK).EQ.13)) THEN 
                           IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
                              WRITE(ScreenUnit,8567)
     &                             JGW,NBV(JGW),ICK,NBV(ICK)
                           ENDIF
                           WRITE(16,8567) JGW,NBV(JGW),ICK,NBV(ICK)
 8567                      FORMAT(////,1X,
     &                          '!!!!!!!!!!  FATAL INPUT ERROR   !!!'
     &               ,'!!!!!!!!!',//,
     &               1X,'BOUNDARY NODE NO. ',I6,' (GLOBAL NODE NO. ',
     &               I9, 'AND OF INTERNAL BARRIER TYPE) ',/,
     &               2X,'OVERLAPS BOUNDARY NODE NO.',I6,' (GLOBAL NODE'
     &               ,' NO.',I6,' )',/,
     &               2X,'THIS IS AN ILLEGAL TYPE OVERLAP !! - INTERNAL '
     &               ,'BARRIER BOUNDARIES CAN ONLY OVERLAP WITH ',
     &               'NO NORMAL FLOW EXTERNAL BOUNDARIES',/
     &               2X,'(I.E. IBTYPE=0,10,20)')
                           IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
                              WRITE(ScreenUnit,9973)
                           ENDIF
                           WRITE(16,9973)
                           STOP
                        ENDIF
C...  CHECK FOR OVERLAPS WHICH REQUIRE ADJUSTMENTS OF BOUNDARY
C...  CODE ON THE EXTERNAL BOUNDARY
                        IF(((IBTYPE.EQ.4).AND.(LBCODEI(ICK).EQ.0))
     &                       .OR.((IBTYPE.EQ.5)
     &                       .AND.(LBCODEI(ICK).EQ.0))) THEN
                           WRITE(16,8568) JGW,ICK,ICK
 8568                      FORMAT(1X,'DUE TO LEGAL OVERLAPPING OF ',
     &               'BOUNDARY NODE',I7,' (WHICH IS AN ESSENTIAL INTER'
     &               ,'NAL BARRIER BOUNDARY NODE)', /,2X,
     &               'AND BOUNDARY NODE',I7,' (WHICH IS AN ESSENTIAL ',
     &               'EXTERNAL NO NORMAL FLOW WITH SLIP BOUNDARY',
     &               ' NODE),',/,2X,
     &               'THE BOUNDARY TYPE FOR BOUNDARY NODE ',I7,
     &               ' IS BEING RESET TO IBTYPE=20',/,2X,
     &               '(NATURAL NO NORMAL FLOW WITH SLIP BOUNDARY) ')
                           LBCODEI(ICK)=20
                        ENDIF
                        IF(((IBTYPE.EQ.4).AND.(LBCODEI(ICK).EQ.10)) 
     &                       .OR.((IBTYPE.EQ.5)
     &                       .AND.(LBCODEI(ICK).EQ.10))) THEN
                           WRITE(16,8569) JGW,ICK,ICK
 8569                      FORMAT(1X,'DUE TO LEGAL OVERLAPPING OF ',
     &               'BOUNDARY NODE ',I7,' (WHICH IS AN ESSENTIAL INTER'
     &               ,'NAL BARRIER BOUNDARY NODE)', /,2X,
     &               'AND BOUNDARY NODE',I7,' (WHICH IS AN ESSENTIAL ',
     &               'EXTERNAL NO NORMAL FLOW WITH NO SLIP BOUNDARY',
     &               ' NODE),',/,2X,
     &               'THE BOUNDARY TYPE FOR BOUNDARY NODE ',I7,
     &               ' IS BEING RESET TO IBTYPE=20',/,2X,
     &               '(NATURAL NO NORMAL FLOW WITH SLIP BOUNDARY) ')
                           LBCODEI(ICK)=20
                        ENDIF
                        IF(((IBTYPE.EQ.24).AND.(LBCODEI(ICK).EQ.10))  
     &                       .OR.((IBTYPE.EQ.25)
     &                       .AND.(LBCODEI(ICK).EQ.10))) THEN
                           WRITE(16,8570) JGW,ICK,ICK
 8570                      FORMAT(1X,'DUE TO LEGAL OVERLAPPING OF ',
     &               'BOUNDARY NODE',I7,' (WHICH IS A NATURAL INTERNAL'
     &               ,' BARRIER BOUNDARY NODE)', /,2X,
     &              'AND BOUNDARY NODE',I7,' (WHICH IS AN ESSENTIAL ',
     &              'EXTERNAL NO NORMAL FLOW WITH NO SLIP BOUNDARY',
     &              ' NODE),',/,2X,
     &              'THE BOUNDARY TYPE FOR BOUNDARY NODE',I7,
     &                   ' IS BEING RESET TO IBTYPE=0',/,2X,
     &               '(ESSENTIAL NO NORMAL FLOW WITH SLIP BOUNDARY) ')
                           LBCODEI(ICK)=0
                        ENDIF
                     ENDIF
                  END DO
               ENDIF
               
            END DO
         END DO
      END DO

C...  ONCE ALL FLOW BOUNDARY NODES HAVE BEEN PROCESSED, CHECK TO MAKE
C...  SURE THAT JGW LE MNVEL.  NOTE, JME MUST BE < JGW.

      IF(MNVEL.LT.JGW) THEN
        IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9947)
        WRITE(16,9947)
 9947   FORMAT(////,1X,'!!!!!!!!!!  FATAL INPUT ERROR   !!!!!!!!!!!!',
     &    //,1X,'THE DIMENSION PARAMETER MNVEL IS LESS THAN ',
     &          'THE TOTAL NUMBER OF FLOW BOUNDARY NODES',
     &    /,1X,'FROM ALL THE SPECIFIED FLOW SEGMENTS COMBINED',/)
        IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9973)
          WRITE(16,9973)
#ifdef CMPI
          call msg_fini()
#endif
          STOP
        ENDIF

      NVEL=JGW
      NVELME=JME

C.....TRANSFER FLOW BOUNDARY INFORMATION INTO NODAL ARRAYS
 
      DO I=1,NP
         LBArray_Pointer(I)=-1
         CSI(I)=0.
         SII(I)=1.
      END DO

      DO I=1,NVELME
         J=ME2GW(I)
         LBArray_Pointer(NBV(J))=J
         CSI(NBV(J))=CSII(J)
         SII(NBV(J))=SIII(J)
      END DO

C...IF ANY NON ZERO NORMAL FLOW BOUNDARIES WERE SPECIFIED, (NFLUXF=1)
C.....READ FORCING INFORMATION FROM UNIT 15 FILE

      NFFR = 0
      IF(NFLUXF.EQ.1) THEN

C.....INPUT FROM THE NUMBER OF FREQUENCIES PRESENT IN NORMAL FLOW FORCING
C......DATA.  IF THIS = 0, NORMAL FLOW DATA IS READ IN FROM THE FORT.20 FILE.

        READ(15,*) NFFR
        MNFFR = NFFR
        IF (NFFR.EQ.0) MNFFR = 1

C.....Allocate space for periodic normal flow boundary conditions
        call alloc_main6()
C
        DO I=1,NVELME
          J=ME2GW(I)
          QNAM(1,J)=0.
          QNPH(1,J)=0.
          END DO

C.....READ IN AND WRITE OUT INFO ON SPECIFIED NORMAL FLOW BOUNDARIES

        WRITE(16,1112)
        WRITE(16,2200)
 2200   FORMAT(//,1X,'NORMAL FLOW BOUNDARY FORCING INFORMATION ',//)
        IF(NFFR.EQ.0) THEN
          WRITE(16,2201)
 2201     FORMAT(/,5X,'NORMAL FLOW VALUES WILL BE READ FROM UNIT 20 ',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE FLOW DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.')
          ENDIF
        IF(NFFR.NE.0) THEN
          WRITE(16,2202) NFFR
 2202     FORMAT(/,5X,'NUMBER OF PERIODIC NORMAL FLOW CONSTITUENTS =',
     &                                                               I5)
          WRITE(16,2203)
 2203     FORMAT(/,7X,'CONSTITUENT #',4X,'FREQUENCY',4X,'NODAL FACTOR',
     &         3X,'EQU.ARG (DEG)',2X,'CONSTITUENT',/)
          DO I=1,NFFR
            READ(15,'(A5)') FBOUNTAG(I)
            READ(15,*) FAMIG(I),FFF(I),FFACE(I)
            WRITE(16,2204) I,FAMIG(I),FFF(I),FFACE(I),FBOUNTAG(I)
 2204       FORMAT(12X,I2,6X,F16.12,2X,F10.7,2X,F10.3,10X,A5)
            FFACE(I)=FFACE(I)*DEG2RAD
            IF(FAMIG(I).EQ.0.) THEN
              FPER(I)=0.
              ELSE
              FPER(I)=2.D0*PI/FAMIG(I)
              ENDIF
            END DO

C.......INPUT PERIODIC NORMAL FLOW FORCING CONDITIONS ON DESIGNATED FLOW BOUNDARIES
C........FOR EACH OF THE FORCING FREQUENCIES FROM UNIT 15 AND OUTPUT TO UNIT 16

          DO I=1,NFFR
            WRITE(16,2206) I,FBOUNTAG(I)
 2206       FORMAT(////,5X,'PERIODIC NORMAL FLOW CONSTITUENT ',
     &                     'NUMBER',I4,1X,'DESIGNATED : ',A5)
            READ(15,'(A10)') ALPHA
            WRITE(16,31) ALPHA
            WRITE(16,30)
            DO J=1,NVEL
               IF((LBCODEI(J).EQ.2).OR.(LBCODEI(J).EQ.12).OR.
     &              (LBCODEI(J).EQ.22).OR.(LBCODEI(J).EQ.52)) THEN
                  READ(15,*) QNAM(I,J),QNPH(I,J)
                WRITE(16,2205) NBV(J),QNAM(I,J),QNPH(I,J)
 2205           FORMAT(10X,I8,4X,F14.5,4X,F12.3)
                QNPH(I,J)=QNPH(I,J)*DEG2RAD
                ENDIF
              IF(LBCODEI(J).EQ.32) THEN
                READ(15,*) QNAM(I,J),QNPH(I,J),ENAM(I,J),ENPH(I,J)
                WRITE(16,2207) NBV(J),QNAM(I,J),QNPH(I,J),ENAM(I,J),
     &                         ENPH(I,J)
 2207           FORMAT(10X,I8,4X,F14.5,4X,F12.3,4X,F14.5,4X,F12.3)
                QNPH(I,J)=QNPH(I,J)*DEG2RAD
                ENPH(I,J)=ENPH(I,J)*DEG2RAD
                ENDIF
              END DO
            END DO
          ENDIF
        ENDIF

C...IF ANY EXTERNAL BARRIER BOUNDARIES WERE SPECIFIED, (NFLUXB=1)
C.....WRITE OUT EXTERNAL BARRIER BOUNDARY INFORMATION TO UNIT 16 FILE
C.....NOTE THAT THIS INFORMATION WAS READ IN FROM THE UNIT 14 FILE

      IF(NFLUXB.EQ.1) THEN

C.....WRITE OUT INFO ON SPECIFIED EXTERNAL BARRIER BOUNDARIES

        WRITE(16,1112)
        WRITE(16,2220)
 2220   FORMAT(//,1X,'EXTERNAL BARRIER BOUNDARY INFORMATION ',/)

C.......OUTPUT ELEVATION OF EXTERNAL BARRIER NODES ABOVE THE GEOID AND
C........THE COEFFICIENT OF FREE SURFACE SUPERCRITICAL FLOW AT
C........DESIGNATED EXTERNAL BARRIER BOUNDARY NODES TO UNIT 16
        WRITE(16,2224)
 2224   FORMAT(//,9X,'NODE',10X,'BARRIER HEIGHT',
     &                     6X,'SUPER-CRIT. EXTERNAL BAR. COEF.',/)
        DO J=1,NVEL
          IF((LBCODEI(J).EQ.3).OR.(LBCODEI(J).EQ.13)
     &       .OR.(LBCODEI(J).EQ.23)) THEN
            WRITE(16,2225) NBV(J),BARLANHT(J),BARLANCFSP(J)
 2225       FORMAT(5X,I8,6X,F14.5,15X,F12.3)
          ENDIF
        END DO
      ENDIF

C...IF ANY INTERNAL BARRIER BOUNDARIES WERE SPECIFIED, (NFLUXIB=1)
C.....WRITE INTERNAL BARRIER BOUNDARY INFORMATION TO UNIT 16 FILE

      IF(NFLUXIB.EQ.1) THEN

C.....WRITE OUT INFO ON SPECIFIED INTERNAL BARRIER BOUNDARIES

        WRITE(16,1112)
        WRITE(16,2320)
 2320   FORMAT(//,1X,'INTERNAL BARRIER BOUNDARY INFORMATION ',/)

C.......WRITE CONNECTION NODE NUMBER AND ELEVATION OF THE INTERNAL BARRIER
C........NODES ABOVE THE GEOID AND THE COEFFICIENTS OF FREE SURFACE SUPERCRITICAL
C........AND SUBCRITICAL FLOW AT DESIGNATED INTERNAL BARRIER BOUNDARY NODES
C........TO UNIT 16 (NOTE THAT THIS INFORMATION WAS INPUT FROM THE UNIT 14
C........FILE WITH BOUNDARY NODE INFORMATION)
        WRITE(16,2324)
 2324   FORMAT(//,9X,'NODE',6X,'CONNECTED NODE',6X,'BARRIER HEIGHT',
     &                4X,'SUB-CRIT. INT. BAR. COEF.',
     &                4X,'SUPER-CRIT. INT. BAR. COEF.',/)
        DO J=1,NVEL
          IF((LBCODEI(J).EQ.4).OR.(LBCODEI(J).EQ.24)) THEN
            WRITE(16,2325) NBV(J),IBCONN(J),BARINHT(J),
     &                     BARINCFSB(J),BARINCFSP(J)
 2325       FORMAT(5X,I8,7X,I8,6X,F14.5,12X,F12.3,17X,F12.3)
          ENDIF
       END DO
      ENDIF

cjjwm001 - begin add                        
C...IF ANY INTERNAL BARRIER BOUNDARIES WITH CROSS BARRIER PIPES
C.....WERE SPECIFIED, (NFLUXIBP=1)
C.....WRITE INTERNAL BARRIER BOUNDARY INFORMATION WITH CROSS 
C.....BARRIER PIPE INFORMATION TO UNIT 16 FILE

      IF(NFLUXIBP.EQ.1) THEN

C.....WRITE OUT INFO ON SPECIFIED INTERNAL BARRIER BOUNDARIES

        WRITE(16,1112)
        WRITE(16,2326)
 2326   FORMAT(//,1X,'INTERNAL BARRIER BOUNDARY WITH CROSS BARRIER',
     &   ' PIPE INFORMATION ',/)

C.......WRITE CONNECTION NODE NUMBER AND ELEVATION OF THE INTERNAL BARRIER
C........NODES ABOVE THE GEOID AND THE COEFFICIENTS OF FREE SURFACE SUPERCRITICAL
C........AND SUBCRITICAL FLOW AT DESIGNATED INTERNAL BARRIER BOUNDARY NODES
C........IN ADDITION TO CROSS BARRIER PIPE CROWN HEIGHT, CROSS BARRIER PIPE
C........COEFFICIENT AND CROSS BARRIER PIPE DIAMETER TO UNIT 16
C........(NOTE THAT THIS INFORMATION WAS INPUT FROM THE UNIT 14 FILE WITH 
C........BOUNDARY NODE INFORMATION)
        WRITE(16,2327)
 2327   FORMAT(//,7X,'NODE',4X,'CONNECTED NODE',4X,'BARRIER HEIGHT',
     &                4X,'SUB-CRIT INT BAR COEF',
     &                4X,'SUPER-CRIT INT BAR COEF',
     &                4X,'PIPEHT  ',
     &                4X,'PIPECOEF',
     &                4X,'PIPEDIAM',/)
        DO J=1,NVEL
          IF((LBCODEI(J).EQ.5).OR.(LBCODEI(J).EQ.25)) THEN
            WRITE(16,2328) NBV(J),IBCONN(J),BARINHT(J),
     &                     BARINCFSB(J),BARINCFSP(J),
     &                     PIPEHT(J),PIPECOEF(J),PIPEDIAM(J)
 2328       FORMAT(3X,I8,5X,I8,4X,F14.5,8X,F12.3,12X,F12.3,
     &              2X,F10.5,2X,F10.5,2X,F10.5)
          ENDIF
       END DO
      ENDIF
cjjwm001 - end add                        

C...
C...READ IN INFORMATION CONCERNING OUTPUT REQUIREMENTS FROM UNIT 15 AND
C...OUTPUT THIS TO UNIT 16
C...
      WRITE(16,1112)
      WRITE(16,3000)
3000  FORMAT(//,1X,'OUTPUT INFORMATION WILL BE PROVIDED AS'
     &  ,' FOLLOWS :')

C...
C...INPUT INFORMATION FOR ELEVATION RECORDING STATIONS
C...

C....READ IN NOUTE,TOUTSE,TOUTFE,NSPOOLE : IF ABS(NOUTE)>0, INTERPOLATED
C....ELEVATIONS AT ELEVATION STATIONS ARE SPOOLED TO UNIT 61 EVERY NSPOOLE
C....TIME STEPS BETWEEN TIMES TOUTSE AND TOUTFE

      READ(15,*) NOUTE,TOUTSE,TOUTFE,NSPOOLE
      WRITE(16,3001) NOUTE
 3001 FORMAT(///,1X,'ELEVATION RECORDING STATION OUTPUT : ',
     &        //,5X,'NOUTE = ',I2)

C....CHECK INPUT PARAMETER NOUTE

      IF(ABS(NOUTE).GT.3) THEN
        IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,3002)
        WRITE(16,3002)
 3002   FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',
     &           //,1X,'YOUR SELECTION OF THE UNIT 15 INPUT PARAMETER',
     &                 ' NOUTE',
     &            /,1X,'IS NOT AN ALLOWABLE VALUE.  CHECK YOUR INPUT!!')
        IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9973)
        WRITE(16,9973)
#ifdef CMPI
        call msg_fini()
#endif
        STOP
        ENDIF

C....IF STATION ELEVATION OUTPUT WILL NOT BE GENERATED

      IF(NOUTE.EQ.0) THEN
        WRITE(16,3003)
 3003   FORMAT(/,5X,'NO OUTPUT WILL BE SPOOLED AT ELEVATION ',
     &              'RECORDING STATIONS')
        ENDIF

C....IF STATION ELEVATION OUTPUT WILL BE GENERATED

      IF(NOUTE.NE.0) THEN

C......COMPUTE NTCYSE, NTCYFE, WHICH = TOUTSE AND TOUTFE IN TIMESTEPS

#ifdef IBM
        NTCYSE=INT((TOUTSE-STATIM)*(86400.D0/DTDP)+0.5d0,KIND(0.0d0))
        NTCYFE=INT((TOUTFE-STATIM)*(86400.D0/DTDP)+0.5d0,KIND(0.0d0))
#else
        NTCYSE=INT((TOUTSE-STATIM)*(86400.D0/DTDP)+0.5d0)
        NTCYFE=INT((TOUTFE-STATIM)*(86400.D0/DTDP)+0.5d0)
#endif
        IF(NTCYFE.GT.NT) NTCYFE=NT

C......COMPUTE NTRSPE = THE NO. OF DATA SETS TO BE SPOOLED TO UNIT 61

        IF(NSPOOLE.EQ.0) NTRSPE=0
#ifdef IBM
        IF(NSPOOLE.NE.0) NTRSPE=INT((NTCYFE-NTCYSE)/NSPOOLE,KIND(0.0d0))
#else
        IF(NSPOOLE.NE.0) NTRSPE=INT((NTCYFE-NTCYSE)/NSPOOLE)
#endif
C......WRITE TOUTSE,TOUTFE,NTCYSE,NTCYFE,NSPOOLE TO UNIT 16

        WRITE(16,3004) TOUTSE,NTCYSE,TOUTFE,NTCYFE,NSPOOLE
 3004   FORMAT(/,5X,'DATA RECORDS WILL START AFTER TOUTSE =',F8.3,
     &              ' DAY(S) RELATIVE',
     &         /,9X,'TO THE STARTING TIME OR',I9,
     &              ' TIME STEPS INTO THE SIMULATION',
     &        //,5X,'DATA RECORDS WILL STOP AFTER TOUTFE =',F8.3,
     &              ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &           I9,' TIME STEPS INTO THE SIMULATION',
     &        //,5X,'INFORMATION WILL BE SPOOLED TO UNIT 61 EVERY',
     &              ' NSPOOLE =',I8,' TIME STEPS')
        IF(ABS(NOUTE).EQ.1) WRITE(16,3005)
 3005   FORMAT(/,5X,'UNIT 61 FORMAT WILL BE ASCII')
        IF(ABS(NOUTE).EQ.2) WRITE(16,3006)
 3006   FORMAT(/,5X,'UNIT 61 FORMAT WILL BE BINARY')
        ENDIF

C....REGARDLESS OF WHETHER NOUTE=0, READ IN THE NUMBER OF ELEVATION
C....RECORDING STATIONS

      READ(15,*) NSTAE
      WRITE(16,3007) NSTAE
 3007 FORMAT(///,5X,'NUMBER OF INPUT ELEVATION RECORDING STATIONS = ',
     &              I5)

      IF(NSTAE.GT.0) THEN
        IF(ICS.EQ.1) WRITE(16,3008)
 3008   FORMAT(/,7X,'STATION #   ELEMENT',9X,'X',13X,'Y',/)
        IF(ICS.EQ.2) WRITE(16,3009)
 3009   FORMAT(/,5X,'STATION   ELEMENT',3X,'LAMBDA(DEG)',
     &           4X,'FEA(DEG)',10X,'XCP',12X,'YCP',/)
        MNSTAE = NSTAE
      ENDIF
      IF (NSTAE.EQ.0) MNSTAE = 1

C  Allocate arrays dimensioned by MNSTAE
       call alloc_main7()


C....INPUT COORDINATES OF ELEVATION RECORDING STATIONS THEN COMPUTE
C....THE ELEMENT NO. THE STATION LIES IN

C     jgf47.12 From cf, memory allocation for netcdf.
      IF(.NOT. ALLOCATED(STATNAME))ALLOCATE(STATNAME(NSTAE,SNAMLEN))

      DO I=1,NSTAE
        NNE(I)=0
        IF(ICS.EQ.1) THEN

C         jgf47.12 From cf, additions for netcdf.
C    START ADDITION BY CF FOR NETCDF OPTION

             IF(NOUTE.NE.3)THEN
                READ(15,*) XEL(I),YEL(I)
              ELSE
                READ(15,'(A132)') STATLINE
                call parse(STATLINE, LVAR)
                XEL(I)=a2f(LVAR(1))
                YEL(I)=a2f(LVAR(2))
                STATNAME(I,SNAMLEN)=LVAR(3)
                WRITE(6,9911) XEL(I), YEL(I), STATNAME(I,SNAMLEN)
             ENDIF
          ELSE

             IF(NOUTE.NE.3)THEN
                 READ(15,*) SLEL(I),SFEL(I)
             ELSE
                 READ(15,'(A132)') STATLINE
                 call parse(STATLINE, LVAR)
                 SLEL(I)=a2f(LVAR(1))
                 SFEL(I)=a2f(LVAR(2))
                 STATNAME(I,SNAMLEN)=LVAR(3)
             ENDIF

C         END OF ADDITION BY CF FOR NETCDF OPTION

          SLEL(I)=SLEL(I)*DEG2RAD
          SFEL(I)=SFEL(I)*DEG2RAD
          CALL CPP(XEL(I),YEL(I),SLEL(I),SFEL(I),SLAM0,SFEA0)
        ENDIF
 9911            FORMAT(F12.3,2X,F12.3, 6X, A50)
 9111            FORMAT(F8.3,2X,F8.3, 6X, A50)
        CALL CoordinateToElement(XEL(I), YEL(I),
     &       NNE(I), I, 'ELEVATION RECORDING STATION   ')

        IF(ICS.EQ.1) THEN
          WRITE(16,1880) I,NNE(I),XEL(I),YEL(I)
1880      FORMAT(8X,I3,6X,I7,2(2X,F14.2))
          ELSE
          WRITE(16,1883) I,NNE(I),SLEL(I)*RAD2DEG,
     &                   SFEL(I)*RAD2DEG,XEL(I),YEL(I)
1883      FORMAT(6X,I3,4X,I7,2(2X,F13.8),2X,2(1X,F13.2))
          ENDIF

C....PRE-COMPUTE INFORMATION REQUIRED TO INTERPOLATE AT ELEV. RECORDING STATIONS
          CALL ComputeInterpolatingFactors(XEL(I), YEL(I), NNE(I),
     &         STAIE1(I), STAIE2(I), STAIE3(I))

        END DO

C...
C...INPUT INFORMATION FOR VELOCITY RECORDING STATIONS
C...

C....READ IN NOUTV,TOUTSV,TOUTFV,NSPOOLV : IF NOUTV<>0,INTERPOLATED VELOCITIES AT
C....VELOCITY STATIONS ARE SPOOLED TO UNIT 62 EVERY NSPOOLV TIME STEPS BETWEEN
C....TIMES TOUTSV AND TOUTFV; IF ABS(NOUTV)=2, OUTPUT WILL BE BINARY

      READ(15,*) NOUTV,TOUTSV,TOUTFV,NSPOOLV
      WRITE(16,3101) NOUTV
 3101 FORMAT(////,1X,'VELOCITY RECORDING STATION OUTPUT : ',
     &         //,5X,'NOUTV = ',I2)

C....CHECK INPUT PARAMETER NOUTV

      IF(ABS(NOUTV).GT.3) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,3102)
        WRITE(16,3102)
 3102   FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',
     &           //,1X,'YOUR SELECTION OF THE UNIT 15 INPUT PARAMETER',
     &                 ' NOUTV',
     &            /,1X,'IS NOT AN ALLOWABLE VALUE.  CHECK YOUR INPUT!!')
        IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9973)
        WRITE(16,9973)
#ifdef CMPI
        call msg_fini()
#endif
        STOP
        ENDIF

C....IF STATION VELOCITY OUTPUT WILL NOT BE GENERATED

      IF(NOUTV.EQ.0) THEN
        WRITE(16,3103)
 3103   FORMAT(///,5X,'NO OUTPUT WILL BE SPOOLED AT VELOCITY',
     &                ' RECORDING STATIONS')
        ENDIF

C....IF STATION VELOCITY OUTPUT WILL BE GENERATED

      IF(NOUTV.NE.0) THEN

C......COMPUTE NTCYSV, NTCYFV, WHICH = TOUTSV AND TOUTFV IN TIME STEPS
#ifdef IBM
        NTCYSV=INT((TOUTSV-STATIM)*(86400.D0/DTDP) + 0.5d0,KIND(0.0d0))
        NTCYFV=INT((TOUTFV-STATIM)*(86400.D0/DTDP) + 0.5d0,KIND(0.0d0))
#else
        NTCYSV=INT((TOUTSV-STATIM)*(86400.D0/DTDP) + 0.5d0)
        NTCYFV=INT((TOUTFV-STATIM)*(86400.D0/DTDP) + 0.5d0)
#endif
        IF(NTCYFV.GT.NT) NTCYFV=NT

C......CALCULATE NTRSPV = THE NO. OF DATA SETS TO BE SPOOLED TO UNIT 62

        IF(NSPOOLV.EQ.0) NTRSPV=0
#ifdef IBM
        IF(NSPOOLV.NE.0) NTRSPV=INT((NTCYFV-NTCYSV)/NSPOOLV,KIND(0.0d0))
#else
        IF(NSPOOLV.NE.0) NTRSPV=INT((NTCYFV-NTCYSV)/NSPOOLV)
#endif

C......WRITE NOUTV,TOUTSV,TOUTFV,NTCYSV,NTCYFV,NSPOOLV TO UNIT 16

        WRITE(16,3104) TOUTSV,NTCYSV,TOUTFV,NTCYFV,NSPOOLV
 3104   FORMAT(/,5X,'DATA RECORDS WILL START AFTER TOUTSV =',F8.3,
     &              ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &           I9,' TIME STEPS INTO THE SIMULATION',
     &        //,5X,'DATA RECORDS WILL STOP AFTER TOUTFV =',F8.3,
     &              ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &           I9,' TIME STEPS INTO THE SIMULATION',
     &        //,5X,'INFORMATION WILL BE SPOOLED TO UNIT 62 EVERY ',
     &              ' NSPOOLV =',I8,' TIME STEPS')
        IF(ABS(NOUTV).EQ.1) WRITE(16,3105)
 3105   FORMAT(/,5X,'UNIT 62 FORMAT WILL BE ASCII')
        IF(ABS(NOUTV).EQ.2) WRITE(16,3106)
 3106   FORMAT(/,5X,'UNIT 62 FORMAT WILL BE BINARY')
        ENDIF

C....REGARDLESS OF WHETHER NOUTV=0, READ IN THE NUMBER OF VELOCITY
C....RECORDING STATIONS

      READ(15,*) NSTAV
      WRITE(16,3107) NSTAV
 3107 FORMAT(////,5X,'NUMBER OF INPUT VELOCITY RECORDING STATIONS = ',
     &            I5)

      IF(NSTAV.GT.0) THEN
        IF(ICS.EQ.1) WRITE(16,3108)
 3108   FORMAT(/,7X,'STATION #   ELEMENT',9X,'X',13X,'Y',/)
        IF(ICS.EQ.2) WRITE(16,3109)
 3109   FORMAT(/,5X,'STATION   ELEMENT',3X,'LAMBDA(DEG)',
     &             4X,'FEA(DEG)',10X,'XCP',12X,'YCP',/)
        MNSTAV = NSTAV
      ENDIF
      IF (NSTAV.EQ.0) MNSTAV = 1

C      Allocate arrays dimensioned by MNSTAV
       call alloc_main8()

C....INPUT COORDINATES OF VELOCITY RECORDING STATIONS
C....THEN COMPUTE ELEMENT NO. WITHIN WHICH STATION LIES

C     jgf47.12 From cf, memory allocation for netcdf.
      IF(.NOT. ALLOCATED(STATNUMBV))ALLOCATE(STATNUMBV(NSTAV,SNUMLEN))
      IF(.NOT. ALLOCATED(STATNAMEV))ALLOCATE(STATNAMEV(NSTAV,SNAMLEN))
!      WRITE(6,*)"DO LOOP - BEFORE ************"

      DO I=1,NSTAV
        NNV(I)=0
        IF(ICS.EQ.1) THEN

C       jgf47.12 From cf,
CC     START ADDITION BY CF FOR NETCDF OPTION

           IF(NOUTV.NE.3) THEN
              READ(15,*) XEV(I),YEV(I)
           ELSE
              READ(15,'(A132)') STATLINE
              call parse(STATLINE, LVAR)
              XEV(I)=a2f(LVAR(1))
              YEV(I)=a2f(LVAR(2))
              STATNAMEV(I,SNAMLEN)=LVAR(3)
              WRITE(6,9911) XEV(I), YEV(I), STATNAMEV(I,SNAMLEN)
           ENDIF
        ELSE
           IF(NOUTV.NE.3)THEN
              READ(15,*) SLEV(I),SFEV(I)
           ELSE
              WRITE(6,*)"STATION INFO - BEFORE ************"
              READ(15,'(A132)') STATLINE
              call parse(STATLINE, LVAR)
              SLEV(I)=a2f(LVAR(1))
              SFEV(I)=a2f(LVAR(2))
              STATNAMEV(I,SNAMLEN)=LVAR(3)
              WRITE(6,*)"parsed stations"
              WRITE(6,9111) SLEV(I),SFEV(I), STATNAMEV(I,SNAMLEN)
              WRITE(6,*)"STATION INFO - AFTER *************",
     &           SLEV(I),SFEV(I), STATNAMEV(I,SNAMLEN)
         ENDIF

CC    END OF ADDITION BY CF FOR NETCDF OPTION

          SLEV(I)=SLEV(I)*DEG2RAD
          SFEV(I)=SFEV(I)*DEG2RAD
          CALL CPP(XEV(I),YEV(I),SLEV(I),SFEV(I),SLAM0,SFEA0)
          ENDIF

          CALL CoordinateToElement(XEV(I), YEV(I),
     &         NNV(I), I, 'VELOCITY RECORDING STATION    ')

          IF(ICS.EQ.1) THEN
             WRITE(16,1880) I,NNV(I),XEV(I),YEV(I)
          ELSE
             WRITE(16,1883) I,NNV(I),SLEV(I)*RAD2DEG,SFEV(I)*RAD2DEG,
     &            XEV(I),YEV(I)
          ENDIF

C....PRE-COMPUTE INFORMATION REQUIRED TO INTERPOLATE AT VEL. RECORDING STATIONS
          CALL ComputeInterpolatingFactors(XEV(I), YEV(I), NNV(I),
     &         STAIV1(I), STAIV2(I), STAIV3(I))

      END DO

C...
C...  IF TRANSPORT IS INCLUDED IN THE RUN, INPUT INFORMATION FOR CONCENTRATION
C...  RECORDING STATIONS
C...
      NOUTC=0
      IF(IM.EQ.10) THEN

C...  READ IN NOUTC,TOUTSC,TOUTFC,NSPOOLC : IF NOUTC<>0,INTERPOLATED
C...  CONCENTRATIONS ARE SPOOLED TO UNIT 81 EVERY NSPOOLC TIME STEPS
C...  BETWEEN TIMES TOUTSC AND TOUTFC; IF ABS(NOUTC)=2, OUTPUT WILL BE BINARY

         READ(15,*) NOUTC,TOUTSC,TOUTFC,NSPOOLC
         WRITE(16,3201) NOUTC
 3201    FORMAT(///,1X,'CONCENTRATION RECORDING STATION OUTPUT : ',
     &          //,5X,'NOUTC = ',I2)
         
C...  CHECK INPUT PARAMETER NOUTC

         IF(ABS(NOUTC).GT.2) THEN
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,3202)
            WRITE(16,3202)
 3202       FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!',
     &           //,1X,'YOUR SELECTION OF THE UNIT 15 INPUT PARAMETER',
     &           ' NOUTC',
     &           /,1X,'IS NOT AN ALLOWABLE VALUE.  CHECK YOUR INPUT!!')
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9973)
            WRITE(16,9973)
#ifdef CMPI
            call msg_fini()
#endif
            STOP
         ENDIF

C...  IF STATION CONCENTRATION OUTPUT WILL NOT BE GENERATED

         IF(NOUTC.EQ.0) THEN
            WRITE(16,3203)
 3203       FORMAT(/,5X,'NO OUTPUT WILL BE SPOOLED AT CONCENTRATION',
     &           ' RECORDING STATIONS')
         ENDIF
         
C...  IF STATION CONCENTRATION OUTPUT WILL BE GENERATED

         NSTAC = 0
         IF(NOUTC.NE.0) THEN
            
C...  COMPUTE NTCYSC, NTCYFC, WHICH = TOUTSC AND TOUTFC IN TIMESTEPS
#ifdef IBM
            NTCYSC=INT((TOUTSC-STATIM)*(86400.D0/DTDP) + 0.5d0, 
     &           KIND(0.0d0))
            NTCYFC=INT((TOUTFC-STATIM)*(86400.D0/DTDP) + 0.5d0, 
     &           KIND(0.0d0))
#else
            NTCYSC=INT((TOUTSC-STATIM)*(86400.D0/DTDP) + 0.5d0) 
            NTCYFC=INT((TOUTFC-STATIM)*(86400.D0/DTDP) + 0.5d0) 
#endif
            IF(NTCYFC.GT.NT) NTCYFC=NT

C...  COMPUTE NTRSPC = THE NO. OF DATA SETS TO BE SPOOLED TO UNIT 81
            
            IF(NSPOOLC.EQ.0) NTRSPC=0
#ifdef IBM
            IF(NSPOOLC.NE.0) NTRSPC=INT((NTCYFC-NTCYSC)/NSPOOLC,
     &           KIND(0.0d0))
#else
            IF(NSPOOLC.NE.0) NTRSPC=INT((NTCYFC-NTCYSC)/NSPOOLC)
#endif

C...  WRITE TOUTSC,TOUTFC,NTCYSC,NTCYFC,NSPOOLC TO UNIT 16

            WRITE(16,3204) TOUTSC,NTCYSC,TOUTFC,NTCYFC,NSPOOLC
 3204       FORMAT(/,5X,'DATA RECORDS WILL START AFTER TOUTSC =',F8.3,
     &                ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &             I9,' TIME STEPS INTO THE SIMULATION',
     &          //,5X,'DATA RECORDS WILL STOP AFTER TOUTFC =',F8.3,
     &                ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &             I9,' TIME STEPS INTO THE SIMULATION',
     &          //,5X,'INFORMATION WILL BE SPOOLED TO UNIT 81 EVERY',
     &                ' NSPOOLC =',I8,' TIME STEPS')
            IF(ABS(NOUTC).EQ.1) WRITE(16,3205)
 3205       FORMAT(/,5X,'UNIT 81 FORMAT WILL BE ASCII')
            IF(ABS(NOUTC).EQ.2) WRITE(16,3206)
 3206       FORMAT(/,5X,'UNIT 81 FORMAT WILL BE BINARY')
         ENDIF

C...  REGARDLESS OF WHETHER NOUTC=0, READ IN THE NUMBER OF CONCENTRATION
C...  RECORDING STATIONS

         READ(15,*) NSTAC
         WRITE(16,3207) NSTAC
 3207    FORMAT(///,5X,'NUMBER OF INPUT CONCENTRATION RECORDING ',
     &        'STATIONS = ',I5)

         IF(NSTAC.GT.0) THEN
            IF(ICS.EQ.1) WRITE(16,3208)
 3208       FORMAT(/,7X,'STATION #   ELEMENT',9X,'X',13X,'Y',/)
            IF(ICS.EQ.2) WRITE(16,3209)
 3209       FORMAT(/,5X,'STATION   ELEMENT',3X,'LAMBDA(DEG)',
     &           4X,'FEA(DEG)',10X,'XCP',12X,'YCP',/)
            MNSTAC = NSTAC
         ENDIF
        
C  Allocate arrays dimensioned by MNSTAC
         call alloc_main9()

C...  INPUT COORDINATES OF CONCENTRATION RECORDING STATIONS
C...  THEN COMPUTE ELEMENT NO. WITHIN WHICH STATION LIES

         DO I=1,NSTAC
            NNC(I)=0
            IF(ICS.EQ.1) THEN
               READ(15,*) XEC(I),YEC(I)
            ELSE
               READ(15,*) SLEC(I),SFEC(I)
               SLEC(I)=SLEC(I)*DEG2RAD
               SFEC(I)=SFEC(I)*DEG2RAD
               CALL CPP(XEC(I),YEC(I),SLEC(I),SFEC(I),SLAM0,SFEA0)
            ENDIF

            CALL CoordinateToElement(XEC(I), YEC(I),
     &         NNC(I), I, 'CONCENTRATION REC. STATION    ')

            IF(ICS.EQ.1) THEN
               WRITE(16,1880) I,NNC(I),XEC(I),YEC(I)
            ELSE
               WRITE(16,1883) I,NNC(I),SLEC(I)*RAD2DEG,
     &              SFEC(I)*RAD2DEG,XEC(I),YEC(I)
            ENDIF

C...  PRE-COMPUTE INFORMATION REQUIRED TO INTERPOLATE AT CONCENTRATION
C...  RECORDING STATIONS
            CALL ComputeInterpolatingFactors(XEC(I), YEC(I), NNC(I),
     &           STAIC1(I), STAIC2(I), STAIC3(I))
            
         END DO
      ENDIF
      IF (NSTAC.EQ.0) MNSTAC = 1
C...  
C...  IF METEOROLOICAL FORCING IS INCLUDED IN THE RUN, INPUT
C...  INFORMATION FOR MET RECORDING STATIONS - OUTPUT
C...  
      NOUTM=0
      NSTAM = 0
C
      IF(NWS.NE.0) THEN

C...  READ IN NOUTM,TOUTSM,TOUTFM,NSPOOLM : IF NOUTM<>0,INTERPOLATED
C...  MET DATA ARE SPOOLED TO UNITS 71&72 EVERY NSPOOLM TIME STEPS
C...  BETWEEN TIMES TOUTSM AND TOUTFM; IF ABS(NOUTM)=2, OUTPUT WILL BE BINARY

         READ(15,*) NOUTM,TOUTSM,TOUTFM,NSPOOLM
         WRITE(16,3211) NOUTM
 3211    FORMAT(///,1X,'METEOROLOGICAL RECORDING STATION OUTPUT : ',
     &        //,5X,'NOUTM = ',I2)

C...  CHECK INPUT PARAMETER NOUTM

         IF(ABS(NOUTM).GT.3) THEN
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,3212)
            WRITE(16,3202)
 3212       FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!',
     &           //,1X,'YOUR SELECTION OF THE UNIT 15 INPUT PARAMETER',
     &           ' NOUTC',
     &           /,1X,'IS NOT AN ALLOWABLE VALUE.  CHECK YOUR INPUT!!')
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9973)
            WRITE(16,9973)
            STOP
         ENDIF

C...  IF STATION METEOROLOGICAL OUTPUT WILL NOT BE GENERATED
         
         IF(NOUTM.EQ.0) THEN
            WRITE(16,3213)
 3213       FORMAT(/,5X,'NO OUTPUT WILL BE SPOOLED AT METEOROLOGICAL',
     &           ' RECORDING STATIONS')
         ENDIF
         
C...  IF STATION MET OUTPUT WILL BE GENERATED

         IF(NOUTM.NE.0) THEN

C...  COMPUTE NTCYSM, NTCYFM, WHICH = TOUTSM AND TOUTFM IN TIMESTEPS
#ifdef IBM
            NTCYSM=INT((TOUTSM-STATIM)*(86400.D0/DTDP) + 0.5d0,
     &           KIND(0.0d0))
            NTCYFM=INT((TOUTFM-STATIM)*(86400.D0/DTDP) + 0.5d0,
     &           KIND(0.0d0))
#else
            NTCYSM=INT((TOUTSM-STATIM)*(86400.D0/DTDP))
            NTCYFM=INT((TOUTFM-STATIM)*(86400.D0/DTDP))
#endif

            IF(NTCYFM.GT.NT) NTCYFM=NT

C...  COMPUTE NTRSPM = THE NO. OF DATA SETS TO BE SPOOLED TO UNITS 71&72
            
            IF(NSPOOLM.EQ.0) NTRSPM=0
#ifdef IBM
            IF(NSPOOLM.NE.0) NTRSPM=INT((NTCYFM-NTCYSM)/NSPOOLM,
     &           KIND(0.0d0))
#else
            IF(NSPOOLM.NE.0) NTRSPM=INT((NTCYFM-NTCYSM)/NSPOOLM)
#endif

C...  WRITE TOUTSM,TOUTFM,NTCYSM,NTCYFM,NSPOOLM TO UNIT 16
            
            WRITE(16,3214) TOUTSM,NTCYSM,TOUTFM,NTCYFM,NSPOOLM
 3214       FORMAT(/,5X,'DATA RECORDS WILL START AFTER TOUTSM =',F8.3,
     &           ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &           I9,' TIME STEPS INTO THE SIMULATION',
     &           //,5X,'DATA RECORDS WILL STOP AFTER TOUTFM =',F8.3,
     &           ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &           I9,' TIME STEPS INTO THE SIMULATION',
     &           //,5X,'INFORMATION WILL BE SPOOLED TO UNITS 71&72',
     &           ' EVERY NSPOOLM =',I8,' TIME STEPS')
            IF(ABS(NOUTM).EQ.1) WRITE(16,3215)
 3215       FORMAT(/,5X,'UNITS 71&72 FORMAT WILL BE ASCII')
            IF(ABS(NOUTM).EQ.2) WRITE(16,3216)
 3216       FORMAT(/,5X,'UNITS 71&72 FORMAT WILL BE BINARY')
         ENDIF

C...  REGARDLESS OF WHETHER NOUTM=0, READ IN THE NUMBER OF METEOROLOGICAL
C...  RECORDING STATIONS

         READ(15,*) NSTAM
         WRITE(16,3217) NSTAM
 3217    FORMAT(///,5X,'NUMBER OF INPUT METEOROLOGICAL RECORDING ',
     &        'STATIONS = ',I5)
         
         IF(NSTAM.GT.0) THEN
            IF(ICS.EQ.1) WRITE(16,3218)
 3218       FORMAT(/,7X,'STATION #   ELEMENT',9X,'X',13X,'Y',/)
            IF(ICS.EQ.2) WRITE(16,3219)
 3219       FORMAT(/,5X,'STATION   ELEMENT',3X,'LAMBDA(DEG)',
     &           4X,'FEA(DEG)',10X,'XCP',12X,'YCP',/)
            MNSTAM = NSTAM
         ENDIF
         
C  Allocate arrays dimensioned by MNSTAM
         call alloc_main10()

C...  INPUT COORDINATES OF METEOROLOGICAL RECORDING STATIONS
C...  THEN COMPUTE ELEMENT NO. WITHIN WHICH STATION LIES
         
C        jgf47.12 from cf, memory allocation for netcdf.         
         IF(.NOT.ALLOCATED(STATNAMEM))ALLOCATE(STATNAMEM(NSTAM,SNAMLEN))


         DO I=1,NSTAM
            NNM(I)=0
            IF(ICS.EQ.1) THEN

C           jgf47.12 from cf,
C    START ADDITION BY CF FOR NETCDF OPTION

             IF(NOUTM.NE.3)THEN
                READ(15,*) XEM(I),YEM(I)
             ELSE
                READ(15,'(A132)') STATLINE
                call parse(STATLINE, LVAR)
                XEM(I)=a2f(LVAR(1))
                YEM(I)=a2f(LVAR(2))
                STATNAMEM(I,SNAMLEN)=LVAR(3)
                WRITE(6,9111) XEM(I), YEM(I), STATNAMEM(I,SNAMLEN)
             ENDIF
          ELSE
             IF(NOUTM.NE.3)THEN
               READ(15,*) SLEM(I),SFEM(I)
             ELSE
               READ(15,'(A132)') STATLINE
               call parse(STATLINE, LVAR)
               SLEM(I)=a2f(LVAR(1))
               SFEM(I)=a2f(LVAR(2))
               STATNAMEM(I,SNAMLEN)=LVAR(3)
             ENDIF

CC    END OF ADDITION BY CF FOR NETCDF OPTION
               SLEM(I)=SLEM(I)*DEG2RAD
               SFEM(I)=SFEM(I)*DEG2RAD
               CALL CPP(XEM(I),YEM(I),SLEM(I),SFEM(I),SLAM0,SFEA0)
            ENDIF

          CALL CoordinateToElement(XEM(I), YEM(I),
     &         NNM(I), I, 'METEOROLOGICAL REC. STATION   ')

            IF(ICS.EQ.1) THEN
               WRITE(16,1880) I,NNM(I),XEM(I),YEM(I)
            ELSE
               WRITE(16,1883) I,NNM(I),SLEM(I)*RAD2DEG,
     &              SFEM(I)*RAD2DEG,XEM(I),YEM(I)
            ENDIF

C...  PRE-COMPUTE INFORMATION REQUIRED TO INTERPOLATE AT METEOROLOGICAL
C...  RECORDING STATIONS

            CALL ComputeInterpolatingFactors(XEM(I), YEM(I), NNM(I),
     &           STAIM1(I), STAIM2(I), STAIM3(I))

         END DO
      ENDIF
      IF (NSTAM.EQ.0) MNSTAM = 1

C...
C...  INPUT INFORMATION ABOUT GLOBAL ELEVATION DATA OUTPUT
C...

C...  READ IN NOUTGE,TOUTSGE,TOUTFGE,NSPOOLGE : IF NOUTGE<>0, GLOBAL ELEV.
C...  OUTPUT IS SPOOLED TO UNIT 63 EVERY NSPOOLGE TIME STEPS BETWEEN
C...  TIMES TOUTSGE AND TOUTFGE; IF ABS(NOUTGE)=2, OUTPUT WILL BE BINARY

      READ(15,*) NOUTGE,TOUTSGE,TOUTFGE,NSPOOLGE
      WRITE(16,3301) NOUTGE
 3301 FORMAT(////,1X,'GLOBAL NODAL ELEVATION INFORMATION OUTPUT: ',
     &     //,5X,'NOUTGE = ',I2)

C...  CHECK INPUT PARAMETER NOUTGE

      IF(ABS(NOUTGE).GT.3) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,3302)
         WRITE(16,3302)
 3302    FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',
     &        //,1X,'YOUR SELECTION OF THE UNIT 15 INPUT PARAMETER',
     &        ' NOUTGE',
     &        /,1X,'IS NOT AN ALLOWABLE VALUE.  CHECK YOUR INPUT!!')
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9973)
         WRITE(16,9973)
#ifdef CMPI
         call msg_fini()
#endif
         STOP
      ENDIF

C...  IF GLOBAL ELEVATION OUTPUT WILL NOT BE GENERATED

      IF(NOUTGE.EQ.0) THEN
         WRITE(16,3303)
 3303    FORMAT(///,5X,'NO GLOBAL ELEVATION OUTPUT WILL BE SPOOLED')
      ENDIF

C...  IF GLOBAL ELEVATION OUTPUT WILL BE GENERATED

      IF(NOUTGE.NE.0) THEN

C...  COMPUTE NTCYSGE, NTCYFGE, WHICH = TOUTSGE AND TOUTFGE IN TIMESTEPS
#ifdef IBM
         NTCYSGE=INT((TOUTSGE-STATIM)*(86400.D0/DTDP) + 0.5d0,
     &        KIND(0.0d0))
         NTCYFGE=INT((TOUTFGE-STATIM)*(86400.D0/DTDP) + 0.5d0, 
     &        KIND(0.0d0))
#else
         NTCYSGE=INT((TOUTSGE-STATIM)*(86400.D0/DTDP) + 0.5d0)
         NTCYFGE=INT((TOUTFGE-STATIM)*(86400.D0/DTDP) + 0.5d0)
#endif
         IF(NTCYFGE.GT.NT) NTCYFGE=NT

C...  CALCULATE NDSETSE = THE # OF DATA SETS TO BE SPOOLED TO UNIT 63

         IF(NSPOOLGE.EQ.0) NDSETSE=0
#ifdef IBM
         IF(NSPOOLGE.NE.0) NDSETSE=INT((NTCYFGE-NTCYSGE)/NSPOOLGE,
     &        KIND(0.0d0))
#else
         IF(NSPOOLGE.NE.0) NDSETSE=INT((NTCYFGE-NTCYSGE)/NSPOOLGE)
#endif

C...  WRITE NOUTGE,TOUTSGE,TOUTFGE,NTCYSGE,NTCYFGE,NSPOOLGE TO UNIT 16

         WRITE(16,3304) TOUTSGE,NTCYSGE,TOUTFGE,NTCYFGE,NSPOOLGE
 3304    FORMAT(/,5X,'DATA RECORDS WILL START AFTER TOUTSGE =',F8.3,
     &              ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &           I9,' TIME STEPS INTO THE SIMULATION',
     &        //,5X,'DATA RECORDS WILL STOP AFTER TOUTFGE =',F8.3,
     &              ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &           I9,' TIME STEPS INTO THE SIMULATION',
     &        //,5X,'INFORMATION WILL BE SPOOLED TO UNIT 63 EVERY ',
     &              'NSPOOLGE =',I8,' TIME STEPS')
         IF(ABS(NOUTGE).EQ.1) WRITE(16,3305)
 3305    FORMAT(/,5X,'UNIT 63 FORMAT WILL BE ASCII')
         IF(ABS(NOUTGE).EQ.2) WRITE(16,3306)
 3306    FORMAT(/,5X,'UNIT 63 FORMAT WILL BE BINARY')
         IF(ABS(NOUTGE).EQ.4) WRITE(16,3307)
 3307    FORMAT(/,5X,'UNIT 63 FORMAT WILL BE COMPACT ASCII')

      ENDIF

C...
C...  INPUT INFORMATION ABOUT GLOBAL VELOCITY DATA OUTPUT
C...

C...  READ IN NOUTGV,TOUTSGV,TOUTFGV,NSPOOLGV : IF NOUTGV<>0, GLOBAL VEL.
C...  OUTPUT IS SPOOLED TO UNIT 64 EVERY NSPOOLGV TIME STEPS BETWEEN
C...  TIMES TOUTSGV AND TOUTFGV; IF ABS(NOUTGV)=2, OUTPUT WILL BE BINARY

      READ(15,*) NOUTGV,TOUTSGV,TOUTFGV,NSPOOLGV
      WRITE(16,3351) NOUTGV
 3351 FORMAT(////,1X,'GLOBAL NODAL VELOCITY INFORMATION OUTPUT : ',
     &     //,5X,'NOUTGV = ',I2)

C...  CHECK INPUT PARAMETER NOUTGV

      IF(ABS(NOUTGV).GT.4) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,3352)
         WRITE(16,3352)
 3352    FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',
     &        //,1X,'YOUR SELECTION OF THE UNIT 15 INPUT PARAMETER',
     &        ' NOUTGV',
     &        /,1X,'IS NOT AN ALLOWABLE VALUE.  CHECK YOUR INPUT!!')
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9973)
         WRITE(16,9973)
#ifdef CMPI
         call msg_fini()
#endif
         STOP
      ENDIF

C...  IF GLOBAL VELOCITY OUTPUT WILL NOT BE GENERATED

      IF(NOUTGV.EQ.0) THEN
         WRITE(16,3353)
 3353    FORMAT(///,5X,'NO GLOBAL VELOCITY OUTPUT WILL BE SPOOLED')
      ENDIF

C...  IF GLOBAL VELOCITY OUTPUT WILL BE GENERATED

      IF(NOUTGV.NE.0) THEN

C...  COMPUTE NTCYSGV, NTCYFGV, WHICH = TOUTSGV AND TOUTFGV IN TIMESTEPS
#ifdef IBM
          NTCYSGV=INT((TOUTSGV-STATIM)*(86400.D0/DTDP) + 0.5d0, 
     &        KIND(0.0d0))
          NTCYFGV=INT((TOUTFGV-STATIM)*(86400.D0/DTDP) + 0.5d0, 
     &         KIND(0.0d0))
#else
          NTCYSGV=INT((TOUTSGV-STATIM)*(86400.D0/DTDP) + 0.5d0)
          NTCYFGV=INT((TOUTFGV-STATIM)*(86400.D0/DTDP) + 0.5d0)
#endif
          IF(NTCYFGV.GT.NT) NTCYFGV=NT

C...  CALCULATE NDSETSV = THE # OF DATA SETS TO BE SPOOLED TO UNIT 64

          IF(NSPOOLGV.EQ.0) NDSETSV=0
#ifdef IBM
          IF(NSPOOLGV.NE.0) NDSETSV=INT((NTCYFGV-NTCYSGV)/NSPOOLGV, 
     &         KIND(0.0d0))
#else
          IF(NSPOOLGV.NE.0) NDSETSV=INT((NTCYFGV-NTCYSGV)/NSPOOLGV)
#endif
C...  WRITE NOUTGV,TOUTSGV,TOUTFGV,NTCYSGV,NTCYFGV,NSPOOLGV TO UNIT 16

          WRITE(16,3354) TOUTSGV,NTCYSGV,TOUTFGV,NTCYFGV,NSPOOLGV
 3354     FORMAT(/,5X,'DATA RECORDS WILL START AFTER TOUTSGV =',F8.3,
     &         ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &         I9,' TIME STEPS INTO THE SIMULATION',
     &         //,5X,'DATA RECORDS WILL STOP AFTER TOUTFGV =',F8.3,
     &         ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &         I9,' TIME STEPS INTO THE SIMULATION',
     &         //,5X,'INFORMATION WILL BE SPOOLED TO UNIT 64 EVERY ',
     &         'NSPOOLGV =',I8,' TIME STEPS')
          IF(ABS(NOUTGV).EQ.1) WRITE(16,3355)
 3355     FORMAT(/,5X,'UNIT 64 FORMAT WILL BE ASCII')
          IF(ABS(NOUTGV).EQ.2) WRITE(16,3356)
 3356     FORMAT(/,5X,'UNIT 64 FORMAT WILL BE BINARY')
          IF(ABS(NOUTGV).EQ.4) WRITE(16,3357)
 3357     FORMAT(/,5X,'UNIT 64 FORMAT WILL BE COMPACT ASCII')

       ENDIF
       
C...
C...  IF TRANSPORT IS INCLUDED IN THE RUN, INPUT INFORMATION ABOUT GLOBAL
C...  CONCENTRATION DATA OUTPUT
C...
       NOUTGC=0
       IF(IM.EQ.10) THEN

C...  READ IN NOUTGC,TOUTSGC,TOUTFGC,NSPOOLGC : IF NOUTGC<>0, GLOBAL
C...  CONCENTRATION OUTPUT IS SPOOLED TO UNIT 73 EVERY NSPOOLGC TIME
C...  STEPS BETWEEN TIMES TOUTSGC AND TOUTFGC; IF ABS(NOUTGC)=2, OUTPUT
C...  WILL BE BINARY

          READ(15,*) NOUTGC,TOUTSGC,TOUTFGC,NSPOOLGC
          WRITE(16,3401) NOUTGC
 3401  FORMAT(////,1X,'GLOBAL NODAL CONCENTRATION INFORMATION OUTPUT:',
     &         //,5X,'NOUTGC = ',I2)

C...  CHECK INPUT PARAMETER NOUTGC

          IF(ABS(NOUTGC).GT.2) THEN
             IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,3402)
             WRITE(16,3402)
 3402     FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',
     &           //,1X,'YOUR SELECTION OF THE UNIT 15 INPUT PARAMETER',
     &                 ' NOUTGC',
     &            /,1X,'IS NOT AN ALLOWABLE VALUE.  CHECK YOUR INPUT!!')
             IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9973)
             WRITE(16,9973)
#ifdef CMPI
             call msg_fini()
#endif
             STOP
          ENDIF

C...  IF GLOBAL CONCENTRATION OUTPUT WILL NOT BE GENERATED

        IF(NOUTGC.EQ.0) THEN
           WRITE(16,3403)
 3403      FORMAT(///,5X,'NO GLOBAL CONCENTRATION OUTPUT WILL BE ',
     &          'SPOOLED')
        ENDIF

C...  IF GLOBAL CONCENTRATION OUTPUT WILL BE GENERATED

        IF(NOUTGC.NE.0) THEN

C...  COMPUTE NTCYSGC, NTCYFGC, WHICH = TOUTSGC AND TOUTFGC IN TIMESTEPS
#ifdef IBM
           NTCYSGC=INT((TOUTSGC-STATIM)*(86400.D0/DTDP) + 0.5d0, 
     &          KIND(0.0d0))
           NTCYFGC=INT((TOUTFGC-STATIM)*(86400.D0/DTDP) + 0.5d0,
     &          KIND(0.0d0))
#else
           NTCYSGC=INT((TOUTSGC-STATIM)*(86400.D0/DTDP) + 0.5d0)
           NTCYFGC=INT((TOUTFGC-STATIM)*(86400.D0/DTDP) + 0.5d0)
#endif
           IF(NTCYFGC.GT.NT) NTCYFGC=NT

C...  CALCULATE NDSETSC = THE # OF DATA SETS TO BE SPOOLED TO UNIT 73

           IF(NSPOOLGC.EQ.0) NDSETSC=0
#ifdef IBM
           IF(NSPOOLGC.NE.0) NDSETSC=INT((NTCYFGC-NTCYSGC)/NSPOOLGC, 
     &          KIND(0.0d0))
#else
           IF(NSPOOLGC.NE.0) NDSETSC=INT((NTCYFGC-NTCYSGC)/NSPOOLGC)
#endif

C...  WRITE NOUTGC,TOUTSGC,TOUTFGC,NTCYSGC,NTCYFGC,NSPOOLGC TO UNIT 16

           WRITE(16,3404) TOUTSGC,NTCYSGC,TOUTFGC,NTCYFGC,NSPOOLGC
 3404      FORMAT(/,5X,'DATA RECORDS WILL START AFTER TOUTSGC =',F8.3,
     &                ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &             I9,' TIME STEPS INTO THE SIMULATION',
     &          //,5X,'DATA RECORDS WILL STOP AFTER TOUTFGC =',F8.3,
     &                ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &             I9,' TIME STEPS INTO THE SIMULATION',
     &          //,5X,'INFORMATION WILL BE SPOOLED TO UNIT 73 EVERY ',
     &                'NSPOOLGC =',I8,' TIME STEPS')
           IF(ABS(NOUTGC).EQ.1) WRITE(16,3405)
 3405      FORMAT(/,5X,'UNIT 73 FORMAT WILL BE ASCII')
           IF(ABS(NOUTGC).EQ.2) WRITE(16,3406)
 3406      FORMAT(/,5X,'UNIT 73 FORMAT WILL BE BINARY')
        ENDIF

      ENDIF

C...
C...  IF NWS<>0   INPUT INFORMATION ABOUT GLOBAL WIND DATA OUTPUT
C...
      IF(NWS.NE.0) THEN

C...  READ IN NOUTGW,TOUTSGW,TOUTFGW,NSPOOLGW : IF NOUTGW<>0, GLOBAL
C...  WIND OUTPUT IS SPOOLED TO UNIT 74 EVERY NSPOOLGW TIME STEPS
C...  BETWEEN TIMES TOUTSGW AND TOUTFGW; IF ABS(NOUTGW)=2, OUTPUT WILL
C...  BE BINARY

         READ(15,*) NOUTGW,TOUTSGW,TOUTFGW,NSPOOLGW
         WRITE(16,3451) NOUTGW
 3451    FORMAT(////,1X,'GLOBAL WIND STRESS INFORMATION OUTPUT : ',
     &           //,5X,'NOUTGW = ',I2)

C...  CHECK INPUT PARAMETER NOUTGW

         IF(ABS(NOUTGW).GT.3) THEN
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,3452)
            WRITE(16,3452)
 3452       FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!',
     &           //,1X,'YOUR SELECTION OF THE UNIT 15 INPUT PARAMETER',
     &                 ' NOUTGW',
     &            /,1X,'IS NOT AN ALLOWABLE VALUE.  CHECK YOUR INPUT!!')
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,3453)
            WRITE(16,3453)
#ifdef CMPI
            call msg_fini()
#endif
            STOP
         ENDIF
          
C...  IF GLOBAL WIND STRESS OUTPUT WILL NOT BE GENERATED

         IF(NOUTGW.EQ.0) THEN
            WRITE(16,3453)
 3453    FORMAT(///,5X,'NO GLOBAL WIND STRESS OUTPUT WILL BE SPOOLED')
         ENDIF

C...  IF GLOBAL WIND STRESS OUTPUT WILL BE GENERATED

        IF(NOUTGW.NE.0) THEN

C...  COMPUTE NTCYSGW, NTCYFGW, WHICH = TOUTSGW AND TOUTFGW IN TIMESTEPS
#ifdef IBM
           NTCYSGW=INT((TOUTSGW-STATIM)*(86400.D0/DTDP) + 0.5d0, 
     &          KIND(0.0d0))
           NTCYFGW=INT((TOUTFGW-STATIM)*(86400.D0/DTDP) + 0.5d0, 
     &          KIND(0.0d0))
#else
           NTCYSGW=INT((TOUTSGW-STATIM)*(86400.D0/DTDP) + 0.5d0)
           NTCYFGW=INT((TOUTFGW-STATIM)*(86400.D0/DTDP) + 0.5d0)
#endif
           IF(NTCYFGW.GT.NT) NTCYFGW=NT

C...  CALCULATE NDSETSW = THE # OF DATA SETS TO BE SPOOLED TO UNIT 74
           
           IF(NSPOOLGW.EQ.0) NDSETSW=0
#ifdef IBM
           IF(NSPOOLGW.NE.0) NDSETSW=INT((NTCYFGW-NTCYSGW)/NSPOOLGW,
     &          KIND(0.0d0))
#else
           IF(NSPOOLGW.NE.0) NDSETSW=INT((NTCYFGW-NTCYSGW)/NSPOOLGW)
#endif

C...  WRITE NOUTGW,TOUTSGW,TOUTFGW,NTCYSGW,NTCYFGW,NSPOOLGW TO UNIT 16

           WRITE(16,3454) TOUTSGW,NTCYSGW,TOUTFGW,NTCYFGW,NSPOOLGW
 3454      FORMAT(/,5X,'DATA RECORDS WILL START AFTER TOUTSGW =',F8.3,
     &                ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &             I9,' TIME STEPS INTO THE SIMULATION',
     &          //,5X,'DATA RECORDS WILL STOP AFTER TOUTFGW =',F8.3,
     &                ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &             I9,' TIME STEPS INTO THE SIMULATION',
     &          //,5X,'INFORMATION WILL BE SPOOLED TO UNIT 74 EVERY ',
     &                'NSPOOLGW =',I8,' TIME STEPS')
           IF(ABS(NOUTGW).EQ.1) WRITE(16,3455)
 3455      FORMAT(/,5X,'UNIT 74 FORMAT WILL BE ASCII')
           IF(ABS(NOUTGW).EQ.2) WRITE(16,3456)
 3456      FORMAT(/,5X,'UNIT 74 FORMAT WILL BE BINARY')
           IF(ABS(NOUTGW).EQ.3) WRITE(16,3457)
 3457      FORMAT(/,5X,'UNIT 74 FORMAT WILL BE SPARSE ASCII')
        ENDIF

      ENDIF

C...
C...  READ AND CHECK INFORMATION ABOUT HARMONIC ANALYSIS OF MODEL RESULTS
C...  
      READ(15,*) NFREQ 
      WRITE(16,99392) NFREQ  
99392 FORMAT(////,1X,'HARMONIC ANALYSIS INFORMATION OUTPUT : ',
     &  //,5X,'HARMONIC ANALYSIS PERFORMED FOR ',I4,' CONSTITUENTS',/)
      MNHARF = NFREQ

      IF (NFREQ.EQ.0) MNHARF = 1

C     allocate harmonic analysis arrays

      IF (NFREQ.GT.0) THEN
         CALL ALLOC_HA()
         CALL ALLOC_MAIN14()
      ENDIF
      
      IF(NFREQ.LT.0) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,99391)
         WRITE(16,99391)
99391    FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',
     &        //,1X,'YOUR SELECTION OF NHARFR (A UNIT 15 '
     &        ,'INPUT PARAMETER) IS NOT AN ALLOWABLE VALUE',/,1X,
     &        'PLEASE CHECK YOUR INPUT',
     &        //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
#ifdef CMPI
         call msg_fini()
#endif
         STOP
      ENDIF
      IF(NFREQ.GT.0) WRITE(16,2330)
 2330 FORMAT(/,7X,'FREQUENCY',4X,'NODAL FACTOR',6X,'EQU.ARG(DEG)',
     +     1X,'CONSTITUENT',/)
      DO 1201 I=1,NFREQ  
         READ(15,'(A10)') NAMEFR(I)
         READ(15,*) HAFREQ(I),HAFF(I),HAFACE(I)
         WRITE(16,2331) HAFREQ(I),HAFF(I),HAFACE(I),NAMEFR(I)
 2331    FORMAT(4X,F15.12,2X,F10.7,5X,F10.3,7X,A10)
 1201 CONTINUE
      
C     read in interval information for harmonic analysis
C     compute thas and thaf in terms of the number of time steps

      READ(15,*) THAS,THAF,NHAINC,FMV
#ifdef IBM
      ITHAS=INT((THAS-STATIM)*(86400.D0/DTDP) + 0.5d0, KIND(0.0d0))
#else
      ITHAS=INT((THAS-STATIM)*(86400.D0/DTDP) + 0.5d0)
#endif
      THAS=ITHAS*DTDP/86400.D0 + STATIM
#ifdef IBM
      ITHAF=INT((THAF-STATIM)*(86400.D0/DTDP) + 0.5d0, KIND(0.0d0))
#else
      ITHAF=INT((THAF-STATIM)*(86400.D0/DTDP) + 0.5d0)
#endif
      THAF=ITHAF*DTDP/86400.D0 + STATIM
      ITMV = ITHAF - (ITHAF-ITHAS)*FMV
      IF(NFREQ.GT.0) THEN
         WRITE(16,34634) THAS,ITHAS,THAF,ITHAF,NHAINC
34634    FORMAT(/,5X,'HARMONIC ANALYSIS WILL START AFTER THAS =',F8.3,
     &        ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',I9,
     &        ' TIME STEPS INTO THE SIMULATION',
     &        //,5X,'HARMONIC ANALYSIS WILL STOP AFTER THAF =',F8.3,
     &        ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',I9,
     &        ' TIME STEPS INTO THE SIMULATION'
     &        ,//,5X,'INFORMATION WILL BE ANALYZED EVERY ',
     &        'NHAINC =',I8,' TIME STEPS.')
         WRITE(16,34639) FMV*100.,ITMV
34639    FORMAT(/,5X,'MEANS AND VARIANCES WILL BE COMPUTED FOR THE ',
     &        'FINAL ',F10.5,' %',/9X,'OF THE HARMONIC ANALYSIS ',
     &        'PERIOD OR AFTER ',I9,' TIME STEPS INTO THE ',
     &        'SIMULATION.',/9X,' RESULTS ARE WRITTEN TO UNIT 55.')
         
      ELSE
         WRITE(16,34645)
34645    FORMAT(///,5X,'NO HARMONIC ANALYSIS WILL BE DONE')
      ENDIF
      
      IF ((FMV.GT.0.).AND.(NFREQ.GT.0).AND.(C2DDI)) CHARMV = .TRUE.
      
C     read in and write out information on where harmonic analysis will
C     be done

      READ(15,*) NHASE,NHASV,NHAGE,NHAGV
      IF((NHASE.LT.0).OR.(NHASE.GT.1)) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,99661)
         WRITE(16,99661)
99661    FORMAT(////,1X,'!!!!!!!!!!  WARNING - NONFATAL ',
     &        'INPUT ERROR  !!!!!!!!!',//
     &        ,1X,'YOUR SELECTION OF NHASE (A UNIT 15 '
     &        ,'INPUT PARAMETER) IS NOT AN ALLOWABLE VALUE',/,1X,
     &        'PLEASE CHECK YOUR INPUT')
         IF(NFOVER.EQ.1) THEN
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,99671)
            WRITE(16,99671)
99671       FORMAT(/,1X,'PROGRAM WILL OVERRIDE SPECIFIED INPUT',
     &           ' AND SET NHASE EQUAL TO 0 ',
     &           //,1X,'!!!!!! EXECUTION WILL CONTINUE !!!!!!',//)
            NHASE=0
         ELSE
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9973)
            WRITE(16,9973)
#ifdef CMPI
            call msg_fini()
#endif
            STOP
         ENDIF
      ENDIF
      IF(NHASE.EQ.1) THEN
         WRITE(16,34641)
34641    FORMAT(///,5X,'STATION ELEVATION HARMONIC ANAL WILL BE ',
     &        'WRITTEN TO UNIT 51')
      ENDIF
      IF((NHASV.LT.0).OR.(NHASV.GT.1)) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,99662)
         WRITE(16,99662)
99662    FORMAT(////,1X,'!!!!!!!!!!  WARNING - NONFATAL ',
     &        'INPUT ERROR  !!!!!!!!!',//
     &        ,1X,'YOUR SELECTION OF NHASV (A UNIT 15 '
     &        ,'INPUT PARAMETER) IS NOT AN ALLOWABLE VALUE',/,1X,
     &        'PLEASE CHECK YOUR INPUT')
         IF(NFOVER.EQ.1) THEN
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,99672)
            WRITE(16,99672)
99672       FORMAT(/,1X,'PROGRAM WILL OVERRIDE SPECIFIED INPUT',
     &           ' AND SET NHASV EQUAL TO 0 ',
     &           //,1X,'!!!!!! EXECUTION WILL CONTINUE !!!!!!',//)
            NHASV=0
         ELSE
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9973)
            WRITE(16,9973)
#ifdef CMPI
            call msg_fini()
#endif
            STOP
         ENDIF
      ENDIF
      IF(NHASV.EQ.1) THEN
         WRITE(16,34642)
34642    FORMAT(///,5X,'STATION VELOCITY HARMONIC ANAL WILL BE ',
     &        'WRITTEN TO UNIT 52')
      ENDIF
      IF((NHAGE.LT.0).OR.(NHAGE.GT.1)) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,99663)
         WRITE(16,99663)
99663    FORMAT(////,1X,'!!!!!!!!!!  WARNING - NONFATAL ',
     &        'INPUT ERROR  !!!!!!!!!',//
     &        ,1X,'YOUR SELECTION OF NHAGE (A UNIT 15 '
     &        ,'INPUT PARAMETER) IS NOT AN ALLOWABLE VALUE',/,1X,
     &        'PLEASE CHECK YOUR INPUT')
         IF(NFOVER.EQ.1) THEN
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,99673)
            WRITE(16,99673)
99673       FORMAT(/,1X,'PROGRAM WILL OVERRIDE SPECIFIED INPUT',
     &           ' AND SET NHAGE EQUAL TO 0 ',
     &           //,1X,'!!!!!! EXECUTION WILL CONTINUE !!!!!!',//)
            NHAGE=0
         ELSE
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9973)
            WRITE(16,9973)
#ifdef CMPI
            call msg_fini()
#endif
            STOP
         ENDIF
      ENDIF
      IF(NHAGE.EQ.1) THEN
         WRITE(16,34643)
34643    FORMAT(///,5X,'GLOBAL ELEVATION HARMONIC ANAL WILL BE ',
     &        'WRITTEN TO UNIT 53')
      ENDIF
      IF((NHAGV.LT.0).OR.(NHAGV.GT.1)) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,99664)
         WRITE(16,99664)
99664    FORMAT(////,1X,'!!!!!!!!!!  WARNING - NONFATAL ',
     &        'INPUT ERROR  !!!!!!!!!',//
     &        ,1X,'YOUR SELECTION OF NHAGV (A UNIT 15 '
     &        ,'INPUT PARAMETER) IS NOT AN ALLOWABLE VALUE',/,1X,
     &        'PLEASE CHECK YOUR INPUT')
         IF(NFOVER.EQ.1) THEN
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,99674)
            WRITE(16,99674)
99674       FORMAT(/,1X,'PROGRAM WILL OVERRIDE SPECIFIED INPUT',
     &           ' AND SET NHAGV EQUAL TO 0 ',
     &           //,1X,'!!!!!! EXECUTION WILL CONTINUE !!!!!!',//)
            NHAGV=0
         ELSE
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9973)
            WRITE(16,9973)
#ifdef CMPI
            call msg_fini()
#endif
            STOP
         ENDIF
      ENDIF
      IF(NHAGV.EQ.1) THEN
         WRITE(16,34644)
34644    FORMAT(///,5X,'GLOBAL VELOCITY HARMONIC ANAL WILL BE ',
     &        'WRITTEN TO UNIT 54')
      ENDIF
        
C     compute flag telling whether any harmonic analysis will be done

      IHARIND=NFREQ*(NHASE+NHASV+NHAGE+NHAGV)
      IF(IHARIND.GT.0) IHARIND=1

C...
C...  Input information about hot start output
C...  
C     jgf45.07 added undocumented option to STOP after writing hot start file.
C     This option will be used in testing ADCIRC's hot start capabilities.
      READ(15,*) NHSTAR,NHSINC
      WRITE(16,99655)
99655 FORMAT(////,1X,'HOT START OUTPUT INFORMATION OUTPUT : ')
      IF((NHSTAR.NE.0).AND.(NHSTAR.NE.1).AND.
     &     (NHSTAR.NE.67).AND.(NHSTAR.NE.68) ) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,99665)
         WRITE(16,99665)
99665    FORMAT(////,1X,'!!!!!!!!!!  WARNING - NONFATAL ',
     &        'INPUT ERROR  !!!!!!!!!',//
     &        ,1X,'YOUR SELECTION OF NHSTAR (A UNIT 15 '
     &        ,'INPUT PARAMETER) IS NOT AN ALLOWABLE VALUE',/,1X,
     &        'PLEASE CHECK YOUR INPUT')
         IF(NFOVER.EQ.1) THEN
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,99675)
            WRITE(16,99675)
99675       FORMAT(/,1X,'PROGRAM WILL OVERRIDE SPECIFIED INPUT',
     &           ' AND SET NHSTAR EQUAL TO 0 ',
     &           //,1X,'!!!!!! EXECUTION WILL CONTINUE !!!!!!',//)
            NHSTAR=0
         ELSE
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9973)
            WRITE(16,9973)
#ifdef CMPI
            call msg_fini()
#endif
            STOP
         ENDIF
      ENDIF
      IF((NHSTAR.GT.0).AND.(NHSINC.EQ.0)) THEN
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,99666)
         WRITE(16,99666)
99666    FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',
     &        //,1X,'YOU HAVE SELECTED A HOT START INTERVAL OF ',
     &        'ZERO LENGTH.',
     &        /,1X,'NHSTAR = 1 AND NHSINC = 0 ARE INCOMPATIBLE ',
     &        'UNIT 15 INPUT PARAMETERS',
     &        /,1X,'PLEASE CHECK YOUR INPUT',
     &        //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
#ifdef CMPI
            call msg_fini()
#endif
         STOP
      ENDIF
      IF(NHSTAR.EQ.0) THEN
         WRITE(16,34646)
34646    FORMAT(///,5X,'NO HOT START OUTPUT WILL BE GENERATED')
      ENDIF
      IF(NHSTAR.GT.0) THEN
         WRITE(16,34636) NHSINC
34636    FORMAT(/,5X,'HOT START OUTPUT WILL BE WRITTEN TO UNIT',
     &        ' 67 OR 68 EVERY ',I5,' TIME STEPS')
      ENDIF
      IF(NHSTAR.EQ.67 .OR. NHSTAR.EQ.68 ) THEN
         WRITE(16,34626) NHSTAR
34626    FORMAT(/,5X,'ADCIRC will stop after writing to unit ',I2)
      ENDIF
      if (IHOT   <  0) NHSTAR = 0
      if (NHSINC <= 0) NHSINC = 1  ! rtm 46.xx NHSINC must have a
                                   ! reasonable value even when not
                                   ! generating hot start files.
C...
C...  Input information about GWCE solver
C...

C     read in and check matrix solver parameters        

      READ(15,*) ITITER,ISLDIA,CONVCR,ITMAX
        
      WRITE(16,99656)
99656 FORMAT(//,1X,'SOLVER INFORMATION OUTPUT : ')
      IF((ISLDIA.LT.0).OR.(ISLDIA.GT.5)) THEN
        IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9920)
        WRITE(16,9920)
 9920   FORMAT(////,1X,'!!!!!!!!!!  WARNING - NONFATAL INPUT ERROR ',
     &                 '!!!!!!!!!',
     &           //,1X,'ISLDIA (A UNIT 15 INPUT PARAMETER) MUST BE 0-5',
     &            /,1X,'PLEASE CHECK YOUR INPUT')
        IF(NFOVER.EQ.1) THEN
          IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9921)
          WRITE(16,9921)
 9921     FORMAT(/,1X,'PROGRAM WILL OVERRIDE SPECIFIED INPUT',
     &                ' AND SET ISLDIA EQUAL TO 0 ',
     &           //,1X,'!!!!!! EXECUTION WILL CONTINUE !!!!!!',//)
          ISLDIA=0
          ELSE
          IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,9973)
          WRITE(16,9973)
#ifdef CMPI
          call msg_fini()
#endif
          STOP
          ENDIF
        ENDIF

C       jgf47.12 Added section for reading netcdf data from fort.15. 
C
#ifdef NETCDF
! CF - Start addition for NETCDF I/O - 6/20/07

!    --------------------------------------------------------------------
!    IF output is in netCDF format read global attributes for netCDF file
!    --------------------------------------------------------------------

      IF( (noutge.EQ.3) .OR. (noutgv.EQ.3) .OR. (noutv.EQ.3) 
     &       .or. (noute.EQ.3) .or. (noutm.eq.3) )THEN

         READ(15,'(A80)') title
         WRITE(6,'(A80)') title
         READ(15,'(A80)') institution
         READ(15,'(A80)') source
         READ(15,'(A80)') history
         READ(15,'(A80)') references
         READ(15,'(A80)') comments
         READ(15,'(A80)') host
         READ(15,'(A80)') convention
         READ(15,'(A80)') contact
         READ(15,'(A80)') base_date
         WRITE(6,'(A80)') base_date
               WRITE(6,*)
     &                title, institution, source, history,
     &                references, comments, host, convention, contact,
     &                base_date

      ENDIF

! CF - Finish addition for NETCDF I/O
#endif

C     allocate arrays needed by GWCE matrix and iterative solver

      call alloc_main11()

C     initialize parameter arrays needed by iterative solver

      CALL DFAULT(IPARM,RPARM)
      IPARM(1)=ITMAX
      IPARM(2)=ISLDIA
      OPEN(33,FILE=TRIM(LOCALDIR)//'/'//'fort.33')
      IPARM(4)=33
      RPARM(1)=CONVCR
      NW = 4*NP + 4*ITMAX
C...
C...  Read input for 3D run
C...
      IF(C3D) THEN
        CALL READ_INPUT_3D(StaTim,NT)
c     ELSEIF(C3DDSS) THEN
c     CALL READ_INPUT_3DDSS(STATIM,NT)
      ENDIF

C...  INITIALIZE AVERAGING FOR INTERNAL BARRIER WATER LEVELS
C...  BARAVGWT=0.000 -> NO AVERAGING PERFORMED
cjjwm001 changed one line                
      BARAVGWT=0.000D0
      IBSTART=0
      DO I=1,NVEL
        RBARWL1AVG(I)=0.D0
        RBARWL2AVG(I)=0.D0
        END DO

C...  INITIALIZE NIBNODECODE(I)

      DO I=1,NP
         NIBNODECODE(I)=0
      END DO

C - allocate arrays dealing with wind forcing

      call alloc_main12()

      if (mnproc == 1) then
        NP_G    = NP               ! 
        NE_G    = NE               !
        NSTAE_G = NSTAE
        NSTAV_G = NSTAV
        NSTAM_G = NSTAM
        NSTAC_G = NSTAC
      endif

C     write table of ADCIRC parameter sizes
      WRITE(16,4010) MNPROC,MNWPROC,MNE,MNP,MNei,MNOPE,MNETA,MNBOU,MNVEL,
     & MNTIF,MNBFR,MNSTAE,MNSTAV,MNSTAC,MNSTAM,NWLAT,NWLON,MNHARF,MNFFR
      IF(NWS.EQ.0) WRITE(16,4011)
      IF(NWS.EQ.1) WRITE(16,4012)
      IF(ABS(NWS).EQ.2) WRITE(16,4013)
      IF(NWS.EQ.3) WRITE(16,4014)
      IF(ABS(NWS).EQ.4) WRITE(16,4015)
      IF(ABS(NWS).EQ.5) WRITE(16,4115)
      IF(ABS(NWS).EQ.7) WRITE(16,4013) !jgf46.00 Added NWS=7 (direct stress)
      IF(NWS.EQ.9) WRITE(16,4018)      !cf & cm added nws = 9 (jgf46.16 merged)
      IF(NWS.EQ.9) WRITE(16,4019)      !cf & cm added nws = 9 (jgf46.16 merged)
      IF(NWS.EQ.10) WRITE(16,4016)
      IF(NWS.EQ.11) WRITE(16,4017)
      IF(ABS(NWS).EQ.12) WRITE(16,4032) ! sb46.28sb01 Added NWS=12 (OWI format)
      IF((NFREQ.EQ.0).OR.(FMV.EQ.0.)) WRITE(16,4021)
      IF((NFREQ.GE.1).AND.(FMV.NE.0.)) WRITE(16,4022)
      IF(ILUMP.EQ.0) WRITE(16,4031)
      IF(ILUMP.EQ.1) WRITE(16,4032)
      IF(IM.EQ.0) WRITE(16,4101)
      IF(IM.EQ.10) WRITE(16,4109)
      IF(IM.EQ.1) WRITE(16,4102)
      IF(IM.EQ.2) WRITE(16,4103)
      WRITE(16,4105)
      WRITE(16,4108)
C
 4010 FORMAT(' *****************************************************',/,
     &       ' *   Based on information extracted from the ADCIRC  *',/,
     &       ' *   UNIT 14 and 15 (grid and horiz run info) files  *',/,
     &       ' *   the following paramter values will be set:      *',/,
     &       ' *                                                   *',/,
     &       ' *       MNPROC = ',I5,'     MWPROC = ',I5,'         *',/,
     &       ' *       MNE = ',I8,1X,'     MNP  = ',I8,1X,'        *',/,
     &       ' *       MNEI = ',I7,2X,'                            *',/,
     &       ' *       MNOPE = ',I6,3X,'   MNETA = ',I6,3X,'       *',/,
     &       ' *       MNBOU = ',I6,3X,'   MNVEL = ',I6,3X,'       *',/,
     &       ' *       MNTIF = ',I6,3X,'   MNBFR = ',I6,3X,'       *',/,
     &       ' *       MNSTAE = ',I5,4X,'  MNSTAV = ',I5,4X,'      *',/,
     &       ' *       MNSTAC = ',I5,4X,'  MNSTAM = ',I5,4X,'      *',/,
     &       ' *       MNWLAT = ',I5,4X,'  MNWLON = ',I5,4X,'      *',/,
     &       ' *       MNHARF = ',I5,4X,'  MNFFR = ',I6,3X,'       *',/,
     &       ' *                                                   *')
 4011 FORMAT(' *   Also, NO wind forcing will be used,             *')
 4012 FORMAT(' *   Also, Standard wind stress and pres will be used,*')
 4013 FORMAT(' *   Also, Semi-standard wind forcing will be used,  *')
 4014 FORMAT(' *   Also, Fleet numeric wind forcing will be used,  *')
 4015 FORMAT(' *   Also, PBL/JAG wind forcing will be used,        *')
 4115 FORMAT(' *   Also, Standard wind vel and pres will be used,  *')
 1236 FORMAT(' *   Also, surface stress forcing will be used,      *')
 4016 FORMAT(' *   Also, AVN wind & pressure forcing will be used, *')
 4017 FORMAT(' *   Also, ETA wind & pressure forcing will be used, *')
 4033 FORMAT(' *   Also, OWI format wind vel and pres will be used,*')
 4018 FORMAT(' *   Asymmetric hurricane wind and pressure forcing  *') 
 4019 FORMAT(' *              will be used,                        *') 
 4021 FORMAT(' *   means and variance calculation will NOT be made,*')
 4022 FORMAT(' *   means and variance calculation will be made,    *')
 4031 FORMAT(' *   the GWCE matrix will be left in consistent form *')
 4032 FORMAT(' *   the GWCE matrix will be LUMPED                  *')
 4101 FORMAT(' *   the model will be set up for a 2DDI run,        *')
 4109 FORMAT(' *   the model will be set up for a 2DDI run + transp*')
 4102 FORMAT(' *   the model will be set up for a 3D-VS run,       *')
 4103 FORMAT(' *   the model will be set up for a 3D-DSS run,      *')
 4105 FORMAT(' *   and an iterative solver will be used            *')
 4108 FORMAT(' *****************************************************',/)
C

C...
C...  CLOSE FILES FORT.14 and FORT.15                 
C...
      CLOSE(14)
      CLOSE(15)
C...  
      RETURN 
      END


c******************************************************************************
C   Subroutine to read in 3D portion of fort.15 file                          *
C                                                                             *
C  Note, initial conditions on density, temperature and/or salinity are read  *
C        in for a cold start in subroutine COLD_START_3D                      *
C                                                                             *
c******************************************************************************

      SUBROUTINE READ_INPUT_3D(StaTime,NT)
      USE SIZES
      USE GLOBAL_3DVS
#ifdef CMPI
      USE MESSENGER, ONLY : msg_fini
#endif
      IMPLICIT NONE !jgf45.09 added
C...  Declaration and definition of local variables used in this subroutine

      REAL(8)  :: StaTime                     !Model start time
      REAL(SZ) :: TO3DSDS, TO3DSDF            !3D Density station output starting & ending times
      REAL(SZ) :: TO3DSVS, TO3DSVF            !3D Velocity station output starting & ending times
      REAL(SZ) :: TO3DSTS, TO3DSTF            !3D Turbulence station output starting & ending times
      REAL(SZ) :: TO3DGDS, TO3DGDF            !3D Global density output starting & ending times
      REAL(SZ) :: TO3DGVS, TO3DGVF            !3D Global velocity output starting & ending times
      REAL(SZ) :: TO3DGTS, TO3DGTF            !3D Global turbulence output starting & ending times
      REAL(8)  :: X3DS,Y3DS,SF3DS,SL3DS       !3D station coordinates
      REAL(8)  :: X4,Y4
      REAL(8)  :: X1,X2,X3,Y1,Y2,Y3           !Used for computing interpolating factors for station output
      REAL(8)  :: A1,A2,A3,AA,AE,AEMin        !Used for computing interpolating factors for station output
      REAL(SZ) :: HH1                         !domain averaged depth used for some vertical FE grids
      REAL(SZ) :: Alp1,Alp2,Alp3              !time stepping coefficients

      INTEGER  :: NT                          !Total number of time steps in model run
      INTEGER  :: N1,N2,N3,KMin               !Used for computing interpolating factors for station output
      INTEGER  :: N,K,NN                      !loop counters
      INTEGER  :: NH                          !horizontal node loop counter
      INTEGER  :: IGC                         !vertical grid code

C     
C...  General invalid entry format statement used throughout
C     
 350  FORMAT(//,2X,'***** INVALID ENTRY IN THE 3D INPUT SECTION OF ',
     &     ' FILE (UNIT 15) ****',/,'****** RUN TERMINATED ******')

C...  
C...  BEGIN READING VERTICAL PARAMETER INFORMATION
C...  
      WRITE(16,300)
 300  FORMAT(//,1X,'3D SOLUTION INFORMATION',/)

C...  SPECIFY WHETHER A BAROTROPIC OR BAROCLINIC RUN

      READ(15,*) IDen
      WRITE(16,421) IDen
 421  FORMAT(/,5X,'IDen = ',I3)
      IF(ABS(IDen).GT.4) THEN
        IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
          WRITE(screenunit,350)
          WRITE(screenunit,423)
        ENDIF
        WRITE(16,350)
        WRITE(16,423)
 423    FORMAT(/,2X,'    IDen MUST = -4, -3, -2, -1, 0, 1, 2, 3, 4')
        STOP
      ENDIF

      IF(IDen.NE.0) CBaroclinic = .TRUE.
      IF(IDen.GT.0) C3D_BTrans  = .TRUE.

C...  READ IN THE TYPE OF BOTTOM BOUNDARY CONDITION AND THE SLIP COEFFICIENTS

      READ(15,*) ISlip,KP
      WRITE(16,355) ISlip,KP
 355  FORMAT(/,5X,'ISlip = ',I3,' KP = ',E10.5)
      IF((ISlip.LT.0).OR.(ISlip.GT.2)) THEN
        IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
          WRITE(ScreenUnit,350)
          WRITE(ScreenUnit,360)
        ENDIF
        WRITE(16,350)
        WRITE(16,360)
 360    FORMAT(/,2X,'    The Bottom Slip Code Must = 0,1,OR 2')
#ifdef CMPI
        call msg_fini()
#endif
        STOP
      ENDIF

C...  READ IN THE SURFACE AND BOTTOM ROUGHNESSES

      READ(15,*) Z0S, Z0B
      WRITE(16,380) Z0S,Z0B
 380  FORMAT(/,5X,'Z0S = ',E10.5,' Z0B = ',E10.5)

C...  READ IN THE TIME STEPPING COEFFICIENTS AND COMPUTE ASSOCIATED VARIABLES

      READ(15,*) Alp1,Alp2,Alp3
      WRITE(16,390) Alp1,Alp2,Alp3
 390  FORMAT(/,5X,'3D TIME STEPPING COEFFS Alp1 = ',E8.2,' Alp2 = ',
     &     E8.2,' Alp3 = ',E8.2)

      IDTAlp1 = I*DelT*Alp1
      IDT1MAlp1 = I*DelT*(1.-Alp1)
      DTAlp3 = DelT*Alp3
      DT1MAlp3 = DelT*(1-Alp3)
      DTAlp2 = DelT*Alp2
      DT1MAlp2 = DelT*(1.-Alp2)

C...  READ IN IGC & NFEN: F.E. GRID CODE & # NODES IN F.E. GRID

      READ(15,*) IGC,NFEN
      WRITE(16,400) IGC,NFEN
 400  FORMAT(/,5X,'Vertical grid code IGC = ',I3,
     &            '  Number of vertical nodes (NFEN) = ',I5)
      IF((IGC.LT.0).OR.(IGC.GT.6)) THEN
        IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
          WRITE(screenunit,350)
          WRITE(screenunit,408)
        ENDIF
        WRITE(16,350)
        WRITE(16,408)
 408    FORMAT(/,2X,'    IGC MUST BE 0, 1, 2, 3, 4, 5, 6')
        STOP
      ENDIF
      IF(NFEN.LT.0) THEN
        IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
          WRITE(screenunit,350)
          WRITE(screenunit,409)
        ENDIF
        WRITE(16,350)
        WRITE(16,409)
 409    FORMAT(/,2X,'    NFEN MUST BE > 0')
        STOP
      ENDIF

C...  SET MNFEN = NFEN

      MNFEN = NFEN

C...  ALLOCATE GENERAL 3D ARRAYS

      CALL ALLOC_3DVS()

C...  READ IN OR SET UP Vertical F.E. GRID

c     IGC = 0 - Read in grid from UNIT 15
     
      IF(IGC.EQ.0) then
        DO N=1,NFEN
          READ(15,*) Sigma(N)
        ENDDO
        IF(Sigma(1).NE.B) THEN
          IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
            WRITE(screenunit,350)
            WRITE(screenunit,1011)
            WRITE(screenunit,1012)
          ENDIF
          WRITE(16,350)
          WRITE(16,1011)
          WRITE(16,1012)
 1011     FORMAT(' Error reading in the vertical finite element grid')
 1012     FORMAT(' The first point in the finite element grid ',
     +           'must = b (-1) : run terminated'/)
          STOP
        ENDIF
        IF(Sigma(NFEN).NE.A) THEN
          IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
            WRITE(screenunit,350)
            WRITE(screenunit,1011)
            WRITE(screenunit,1013)
          ENDIF
          WRITE(16,350)
          WRITE(16,1011)
          WRITE(16,1013)
 1013     FORMAT(' The last point in the finite element grid ',
     +           'must = a (1) : run terminated'/)
          STOP
        ENDIF
      ENDIF

c     IGC <> 0 - Set up grid in subroutine FEGRIDS

      IF(IGC.NE.0) THEN
         HH1=0.d0
         DO NH=1,NP
            HH1=HH1+DP(NH)
         ENDDO 
         HH1=HH1/NP                        !domain averaged depth
         CALL FEGRIDS(IGC,HH1)
      ENDIF

c     write out the vertical grid in fort.16 file

      WRITE(16,1000)
 1000 FORMAT(//,5X,'Vertical Grid Information')
      WRITE(16,1001) 
 1001 FORMAT(/,5X,'V. Node #',5X,'V. Position',/)
      DO N = 1,NFEN
         WRITE(16,*) N,Sigma(N)
      ENDDO 


C...  SPECIFY TYPE OF EDDY VISCOSITY PROFILE

      READ(15,*) IEVC,EVMin,EVCon
      WRITE(16,410) IEVC,EVMin,EVCon
 410  FORMAT(/,5X,'IEVC = ',I3,2X,'EVMin = ',E15.8,2X,'EVCon = ',E15.8)
      IF((IEVC.NE.0 ).and.(IEVC.NE.1 ).and.
     &   (IEVC.NE.10).and.(IEVC.NE.11).and.
     &   (IEVC.NE.20).and.(IEVC.NE.21).and.
     &   (IEVC.NE.30).and.(IEVC.NE.31).and.(IEVC.NE.32).and.
     &   (IEVC.NE.33).and.
     &   (IEVC.NE.40).and.(IEVC.NE.41).and.(IEVC.NE.42).and.
     &   (IEVC.NE.43).and.
     &   (IEVC.NE.50)) THEN
        IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
          WRITE(screenunit,350)
          WRITE(screenunit,411)
        ENDIF
        WRITE(16,350)
        WRITE(16,411)
 411    FORMAT(/,2X,'    IEVC MUST BE 0,1,10,11,20,21,30,31,32,33,40,',
     &        '41,42,43,50')
        STOP
      ENDIF
      IF(IEVC.EQ.50) READ(15,*) Theta1,Theta2

C...  FOR IEVC=0, CONSTANT EDDY VISCOSITY, READ IN PROFILE

      IF(IEVC.EQ.0) THEN
        DO N=1,NFEN
          READ(15,*) EVTot(N)
        ENDDO
        WRITE(16,*) ' Vertical E.V. read in from UNIT 15'
      ENDIF

C
C...  READ IN 3D OUTPUT CONTROLS, COMPUTE NEEDED ANCILLARY PARAMETERS
C

C...  Format statements used for 3D Station output diagnostic information

 3108 FORMAT(/,7X,'STATION #   ELEMENT',9X,'X',13X,'Y',/)
 3109 FORMAT(/,5X,'STATION #   ELEMENT',3X,'LAMBDA(DEG)',
     &             4X,'FEA(DEG)',10X,'XCP',12X,'YCP',/)
 1880 FORMAT(8X,I6,5X,I9,2(2X,F14.2))
 1883 FORMAT(6X,I6,5X,I9,2(2X,F13.8),2X,2(1X,F13.2))
 9790 FORMAT(/,1X,'PROGRAM WILL ESTIMATE NEAREST ELEMENT',
     &       /,1X,'PROXIMITY INDEX FOR THIS STATION EQUALS ',E15.6,
     &      //,1X,'!!!!!! EXECUTION WILL CONTINUE !!!!!!',//)

C.... General variable definitions

C     NE       - total number of elements in grid
C     Areas(K) - 2x area of element K


C...3D Station Density Output (fort.41)
C
C     I3DSD    - flag for 3D station density output
C                If ABS(I3DSD)>0 densities are interpolated to stations and written out
C                If ABS(I3DSD)=1, output is in ascii. If ABS(I3DSD)=2, output is in binary.
C     TO3DSDS  - starting time in days for 3D station density output
C     TO3DSDF  - ending time in days for 3D station density output
C     NSpo3DSD - 3D station density output interval in timesteps
C     NTO3DSDS  - TO3DSDS converted to timesteps after StaTime
C     NTO3DSDF  - TO3DSDF converted to timesteps after StaTime
C     NDSet3DSD - total number of times that 3D density station data will be written
C     N3DSD    - counter to determine if it is time to write 3D density station data
C     NSta3DD  - number of 3D density stations 
C     MNSta3DD - larger of 1 or NSta3DD, used to dimension arrays
C     NE3DD(N) - element number containing 3D density station N
C     StaI3DD1(N),StaI3DD2(N),StaI3DD3(N) - interpolating factors used to compute output at 3D density station N
C     I3DSDRec - record pointer for ascii or binary output
C
      READ(15,*) I3DSD,TO3DSDS,TO3DSDF,NSpo3DSD

C.... CHECK INPUT PARAMETER I3DSD

      WRITE(16,501) I3DSD
 501  FORMAT(//,1X,'3D DENSITY STATION OUTPUT PARAMETER: I3DSD = ',I2)
      IF(ABS(I3DSD).GT.2) THEN
        IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
          WRITE(ScreenUnit,350)          
          WRITE(ScreenUnit,511)
        ENDIF
        WRITE(16,350)
        WRITE(16,511)
 511    FORMAT(/,1X,'YOUR SELECTION OF I3DSD IS NOT VALID. ',
     &              'CHECK YOUR INPUT!!')
#ifdef CMPI
        call msg_fini()
#endif
        STOP
      ENDIF
 
      IF((I3DSD.NE.0).AND.(NSpo3DSD.EQ.0)) THEN
        IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
          WRITE(ScreenUnit,350)
          WRITE(ScreenUnit,561) NSpo3DSD
        ENDIF
        WRITE(16,350)
        WRITE(16,561) NSpo3DSD
 561    FORMAT(2X,'    NSpo3DSD = ',I3,' IT MUST BE > 0')
#ifdef CMPI
        call msg_fini()
#endif
        STOP
      ENDIF

C.... IF NO STATION DENSITY OUTPUT WILL BE GENERATED

      IF(I3DSD.EQ.0) THEN
        WRITE(16,521)
 521    FORMAT(/,5X,'No 3D Station density will be output')

C.... IF STATION DENSITY OUTPUT WILL BE GENERATED

      ELSE

C....   COMPUTE NTO3DSDS, NTO3DSDF, WHICH = TO3DSDS AND TO3DSDF IN TIME STEPS
#ifdef IBM
        NTO3DSDS=INT((TO3DSDS-StaTime)*Day2Sec/DTDP+0.5d0,
     &        KIND(0.0d0))      !jgf45.11 was NINT
        NTO3DSDF=INT((TO3DSDF-StaTime)*Day2Sec/DTDP+0.5d0,
     &       KIND(0.0d0))       !jgf45.11 was NINT
#else
        NTO3DSDS=INT((TO3DSDS-StaTime)*Day2Sec/DTDP+0.5d0) !jgf45.11 was NINT
        NTO3DSDF=INT((TO3DSDF-StaTime)*Day2Sec/DTDP+0.5d0) !jgf45.11 was NINT
#endif
         IF(NTO3DSDS.LT.0) THEN
          IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
             WRITE(ScreenUnit,531) TO3DSDS
          ENDIF
          WRITE(16,531) TO3DSDS
 531      FORMAT(5X,'NOTE: TO3DSDS = ',E14.6,' < START TIME. ',
     &              'IT HAS BEEN RESET = START TIME.')
          TO3DSDS=StaTime
          NTO3DSDS=0
        ENDIF
        IF(NTO3DSDF.LT.NTO3DSDS) THEN
          IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
             WRITE(ScreenUnit,541) TO3DSDF
          ENDIF
          WRITE(16,541) TO3DSDF
 541      FORMAT(5X,'NOTE: TO3DSDF = ',E14.6,' < TO3DSDS. ',
     &              'IT HAS BEEN RESET = TO3DSDS.')
          TO3DSDF=TO3DSDS
          NTO3DSDF=NTO3DSDS
        ENDIF
        IF(NTO3DSDF.GT.NT) THEN
          IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
             WRITE(ScreenUnit,551) TO3DSDF
          ENDIF
          WRITE(16,551) TO3DSDF
 551      FORMAT(5X,'NOTE: TO3DSDF = ',E14.6,' > RNDAY. ',
     &              'IT HAS BEEN RESET = RNDAY.')
          NTO3DSDF=NT
        ENDIF

        NDSet3DSD = (NTO3DSDF-NTO3DSDS)/NSpo3DSD
        N3DSD=0
        I3DSDRec=0

C.....  WRITE Density station output information to UNIT 16

        WRITE(16,571) NSpo3DSD,TO3DSDS,(TO3DSDS+NSpo3DSD*DTDP*Sec2Day), 
     &                (NTO3DSDS+NSpo3DSD),TO3DSDF,NTO3DSDF,NSpo3DSD,
     &                NDSet3DSD
 571    FORMAT(/,5X,'3D Density station output will start',I9,
     &              ' timesteps after day TO3DSD = ',F9.3,
     &         /,9X,'This = ',F9.3,' days = ',I9,' timesteps relative',
     &              ' to the simulation Start Time.',
     &         /,9X,'Output will stop ',F9.3,' days = ',I9,' timesteps',
     &              ' relative to the simulation Start Time.',
     &         /,9X,'Data will be written to file fort.42 every ',I7,
     &              ' time steps at total of ',I7,' times.')
      ENDIF
 
C.... REGARDLESS OF WHETHER I3DSD=0, READ IN THE NUMBER OF 3D DENSITY
C.... RECORDING STATIONS

      READ(15,*) NSta3DD
      IF(I3DSD.NE.0) WRITE(16,581) NSta3DD
 581  FORMAT(/,5X,' Output will be written at ',I7,' stations.')

      IF (NSta3DD.EQ.0) THEN
        MNSta3DD=1
        ELSE
        MNSta3DD = NSta3DD
      ENDIF

C  Allocate arrays for station density output

      call alloc_3DSD()

C....Input the coordinates of the density output stations
C....and then compute the element # containing each station

      IF(NSta3DD.GT.0) THEN
        IF(ICS.EQ.1) WRITE(16,3108)
        IF(ICS.EQ.2) WRITE(16,3109)

        DO N=1,NSta3DD
          NE3DD(N)=0
          IF(ICS.EQ.1) THEN
            READ(15,*) X3DS,Y3DS
          ELSE
            READ(15,*) SL3DS,SF3DS
            SL3DS=SL3DS*Deg2Rad
            SF3DS=SF3DS*Deg2Rad
            CALL CPP(X3DS,Y3DS,SL3DS,SF3DS,SLAM0,SFEA0)
          ENDIF

          CALL CoordinateToElement(X3DS, Y3DS,
     &         NE3DD(N), N, '3D Density recording station  ')

C....     Write element and station

          IF(ICS.EQ.1) THEN
            WRITE(16,1880) N,NE3DD(N),X3DS,Y3DS
          ELSE
            WRITE(16,1883) N,NE3DD(N),SL3DS*Rad2Deg,
     &                     SF3DS*Rad2Deg,X3DS,Y3DS
          ENDIF

C....     COMPUTE INTERPOLATING FACTORS
          CALL ComputeInterpolatingFactors(X3DS, Y3DS, NE3DD(N),
     &         StaI3DD1(N), StaI3DD2(N), StaI3DD3(N)) 

        END DO

        WRITE(16,*) ' '
      ENDIF

C...3D Station Velocity Output (fort.42)
C
C     I3DSV    - flag for 3D station velocity output
C                If ABS(I3DSV)>0 velocities are interpolated to stations and written out
C                If ABS(I3DSV)=1, output is in ascii. If ABS(I3DSV)=2, output is in binary.
C     TO3DSVS  - starting time in days for 3D station velocity output
C     TO3DSVF  - ending time in days for 3D station velocity output
C     NSpo3DSV - 3D station velocity output interval in timesteps
C     NTO3DSVS  - TO3DSVS converted to timesteps after StaTime
C     NTO3DSVF  - TO3DSVF converted to timesteps after StaTime
C     NDSet3DSV - total number of times that 3D velocity station data will be written
C     N3DSV    - counter to determine if it is time to write 3D velocity station data
C     NSta3DV  - number of 3D velocity stations 
C     MNSta3DV - larger of 1 or NSta3DV, used to dimension arrays
C     NE3DV(N) - element number containing 3D velocity station N
C     StaI3DV1(N),StaI3DV2(N),StaI3DV3(N) - interpolating factors used to compute output at 3D velocity station N
C     I3DSVRec - record pointer for ascii or binary output
C
      READ(15,*) I3DSV,TO3DSVS,TO3DSVF,NSpo3DSV

C.... CHECK INPUT PARAMETER I3DSV

      WRITE(16,502) I3DSV
 502  FORMAT(//,1X,'3D VELOCITY STATION OUTPUT PARAMETER: I3DSV = ',I2)
      IF(ABS(I3DSV).GT.2) THEN
        IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
          WRITE(ScreenUnit,350)          
          WRITE(ScreenUnit,512)
        ENDIF
        WRITE(16,350)
        WRITE(16,512)
 512    FORMAT(/,1X,'YOUR SELECTION OF I3DSV IS NOT VALID. ',
     &              'CHECK YOUR INPUT!!')
#ifdef CMPI
        call msg_fini()
#endif
        STOP
      ENDIF
 
      IF((I3DSV.NE.0).AND.(NSpo3DSV.EQ.0)) THEN
        IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
          WRITE(ScreenUnit,350)
          WRITE(ScreenUnit,562) NSpo3DSV
        ENDIF
        WRITE(16,350)
        WRITE(16,562) NSpo3DSV
 562    FORMAT(2X,'    NSpo3DSV = ',I3,' IT MUST BE > 0')
#ifdef CMPI
        call msg_fini()
#endif
        STOP
      ENDIF

C.... IF NO STATION VELOCITY OUTPUT WILL BE GENERATED

      IF(I3DSV.EQ.0) THEN
        WRITE(16,522)
 522    FORMAT(/,5X,'No 3D Station velocity will be output')

C.... IF STATION VELOCITY OUTPUT WILL BE GENERATED

      ELSE

C....   COMPUTE NTO3DSVS, NTO3DSVF, WHICH = TO3DSVS AND TO3DSVF IN TIME STEPS
#ifdef IBM
        NTO3DSVS=INT((TO3DSVS-StaTime)*Day2Sec/DTDP+0.5d0,
     &        KIND(0.0d0))   !jgf45.11 was NINT
        NTO3DSVF=INT((TO3DSVF-StaTime)*Day2Sec/DTDP+0.5d0, 
     &       KIND(0.0d0))       !jgf45.11 was NINT
#else
        NTO3DSVS=INT((TO3DSVS-StaTime)*Day2Sec/DTDP+0.5d0) !jgf45.11 was NINT
        NTO3DSVF=INT((TO3DSVF-StaTime)*Day2Sec/DTDP+0.5d0) !jgf45.11 was NINT
#endif
        IF(NTO3DSVS.LT.0) THEN
          IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
             WRITE(ScreenUnit,532) TO3DSVS
          ENDIF
          WRITE(16,532) TO3DSVS
 532      FORMAT(5X,'NOTE: TO3DSVS = ',E14.6,' < START TIME. ',
     &              'IT HAS BEEN RESET = START TIME.')
          TO3DSVS=StaTime
          NTO3DSVS=0
        ENDIF
        IF(NTO3DSVF.LT.NTO3DSVS) THEN
          IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
             WRITE(ScreenUnit,542) TO3DSVF
          ENDIF
          WRITE(16,542) TO3DSVF
 542      FORMAT(5X,'NOTE: TO3DSVF = ',E14.6,' < TO3DSVS. ',
     &              'IT HAS BEEN RESET = TO3DSVS.')
          TO3DSVF=TO3DSVS
          NTO3DSVF=NTO3DSVS
        ENDIF
        IF(NTO3DSVF.GT.NT) THEN
          IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
             WRITE(ScreenUnit,552) TO3DSVF
          ENDIF
          WRITE(16,552) TO3DSVF
 552      FORMAT(5X,'NOTE: TO3DSVF = ',E14.6,' > RNDAY. ',
     &              'IT HAS BEEN RESET = RNDAY.')
          NTO3DSVF=NT
        ENDIF

        NDSet3DSV = (NTO3DSVF-NTO3DSVS)/NSpo3DSV
        N3DSV=0
        I3DSVRec=0

C.....  WRITE Velocity station output information to UNIT 16

        WRITE(16,572) NSpo3DSV,TO3DSVS,(TO3DSVS+NSpo3DSV*DTDP*Sec2Day), 
     &                (NTO3DSVS+NSpo3DSV),TO3DSVF,NTO3DSVF,NSpo3DSV,
     &                NDSet3DSV
 572    FORMAT(/,5X,'3D Velocity station output will start',I9,
     &              ' timesteps after day TO3DSV = ',F9.3,
     &         /,9X,'This = ',F9.3,' days = ',I9,' timesteps relative',
     &              ' to the simulation Start Time.',
     &         /,9X,'Output will stop ',F9.3,' days = ',I9,' timesteps',
     &              ' relative to the simulation Start Time.',
     &         /,9X,'Data will be written to file fort.42 every ',I7,
     &              ' time steps at total of ',I7,' times.')
      ENDIF
 
C.... REGARDLESS OF WHETHER I3DSV=0, READ IN THE NUMBER OF 3D VELOCITY
C.... RECORDING STATIONS

      READ(15,*) NSta3DV
      IF(I3DSV.NE.0) WRITE(16,582) NSta3DV
 582  FORMAT(/,5X,' Output will be written at ',I7,' stations.')

      IF (NSta3DV.EQ.0) THEN
        MNSta3DV=1
        ELSE
        MNSta3DV = NSta3DV
      ENDIF

C  Allocate arrays for station velocity output

      call alloc_3DSV()

C....Input the coordinates of the velocity output stations
C....and then compute the element # containing each station

      IF(NSta3DV.GT.0) THEN
        IF(ICS.EQ.1) WRITE(16,3108)
        IF(ICS.EQ.2) WRITE(16,3109)

        DO N=1,NSta3DV
          NE3DV(N)=0
          IF(ICS.EQ.1) THEN
            READ(15,*) X3DS,Y3DS
          ELSE
            READ(15,*) SL3DS,SF3DS
            SL3DS=SL3DS*Deg2Rad
            SF3DS=SF3DS*Deg2Rad
            CALL CPP(X3DS,Y3DS,SL3DS,SF3DS,SLAM0,SFEA0)
          ENDIF

          CALL CoordinateToElement(X3DS, Y3DS,
     &         NE3DV(N), N, '3D Velocity recording station ')

C....     Write element and station

          IF(ICS.EQ.1) THEN
            WRITE(16,1880) N,NE3DV(N),X3DS,Y3DS
          ELSE
            WRITE(16,1883) N,NE3DV(N),SL3DS*Rad2Deg,
     &                     SF3DS*Rad2Deg,X3DS,Y3DS
          ENDIF

C....     COMPUTE INTERPOLATING FACTORS
          CALL ComputeInterpolatingFactors(X3DS, Y3DS, NE3DV(N),
     &         StaI3DV1(N), StaI3DV2(N), StaI3DV3(N)) 

        END DO

        WRITE(16,*) ' '
      ENDIF

C...3D Station Turbulence Output  (fort.43)
C
C     I3DST    - flag for 3D station turbulence output
C                If ABS(I3DST)>0 turbulence is interpolated to stations and written out
C                If ABS(I3DST)=1, output is in ascii. If ABS(I3DST)=2, output is in binary.
C     TO3DSTS  - starting time in days for 3D station turbulence output
C     TO3DSTF  - ending time in days for 3D station turbulence output
C     NSpo3DST - 3D station turbulence output interval in timesteps
C     NTO3DSTS  - TO3DSTS converted to timesteps after StaTime
C     NTO3DSTF  - TO3DSTF converted to timesteps after StaTime
C     NDSet3DST - total number of times that 3D turbulence station data will be written
C     N3DST    - counter to determine if it is time to write 3D turbulence station data
C     NSta3DT  - number of 3D turbulence stations 
C     MNSta3DT - larger of 1 or NSta3DT, used to dimension arrays
C     NE3DT(N) - element number containing 3D turbulence station N
C     StaI3DT1(N),StaI3DT2(N),StaI3DT3(N) - interpolating factors used to compute output at 3D turbulence station N
C     I3DSTRec - record pointer for ascii or binary output
C
      READ(15,*) I3DST,TO3DSTS,TO3DSTF,NSpo3DST

C.... CHECK INPUT PARAMETER I3DST

      WRITE(16,503) I3DST
 503  FORMAT(//,1X,'3D Turbulence Station output parameter: I3DST = ',
     &                                                               I2)
      IF(ABS(I3DST).GT.2) THEN
        IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
          WRITE(ScreenUnit,350)          
          WRITE(ScreenUnit,513)
        ENDIF
        WRITE(16,350)
        WRITE(16,513)
 513    FORMAT(/,1X,'Your selection of I3DST is not valid. ',
     &              'CHECK YOUR INPUT!!')
#ifdef CMPI
        call msg_fini()
#endif
        STOP
      ENDIF

      IF((I3DST.NE.0).AND.(NSpo3DST.EQ.0)) THEN
        IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
          WRITE(ScreenUnit,350)
          WRITE(ScreenUnit,563) NSpo3DST
        ENDIF
        WRITE(16,350)
        WRITE(16,563) NSpo3DST
 563    FORMAT(2X,'    NSpo3DST = ',I3,' IT MUST BE > 0')
        STOP
      ENDIF

C.... IF NO 3D STATION TURBULENCE OUTPUT WILL BE GENERATED

      IF(I3DST.EQ.0) THEN
        WRITE(16,523)
 523    FORMAT(/,5X,'No 3D Station Turbulence will be output')

C.... IF 3D STATION TURBULENCE OUTPUT WILL BE GENERATED

      ELSE

C....   COMPUTE NTO3DSTS, NTO3DSTF, WHICH = TO3DSTS AND TO3DSTF IN TIME STEPS
#ifdef IBM
        NTO3DSTS=INT((TO3DSTS-StaTime)*Day2Sec/DTDP+0.5d0, 
     &        KIND(0.0d0))   !jgf45.11 was NINT
        NTO3DSTF=INT((TO3DSTF-StaTime)*Day2Sec/DTDP+0.5d0,
     &       KIND(0.0d0))       !jgf45.11 was NINT
#else
        NTO3DSTS=INT((TO3DSTS-StaTime)*Day2Sec/DTDP+0.5d0) !jgf45.11 was NINT
        NTO3DSTF=INT((TO3DSTF-StaTime)*Day2Sec/DTDP+0.5d0) !jgf45.11 was NINT
#endif
        IF(NTO3DSTS.LT.0) THEN
          IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
             WRITE(ScreenUnit,533) TO3DSTS
          ENDIF
          WRITE(16,533) TO3DSTS
 533      FORMAT(5X,'NOTE: TO3DSTS = ',E14.6,' < START TIME. ',
     &              'IT HAS BEEN RESET = START TIME.')
          TO3DSTS=StaTime
          NTO3DSTS=0
        ENDIF
        IF(NTO3DSTF.LT.NTO3DSTS) THEN
          IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
             WRITE(ScreenUnit,543) TO3DSTF
          ENDIF
          WRITE(16,543) TO3DSTF
 543      FORMAT(5X,'NOTE: TO3DSTF = ',E14.6,' < TO3DSTS. ',
     &              'IT HAS BEEN RESET = TO3DSTS.')
          TO3DSTF=TO3DSTS
          NTO3DSTF=NTO3DSTS
        ENDIF
        IF(NTO3DSTF.GT.NT) THEN
          IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
             WRITE(ScreenUnit,553) TO3DSTF
          ENDIF
          WRITE(16,553) TO3DSTF
 553      FORMAT(5X,'NOTE: TO3DSTF = ',E14.6,' > RNDay. ',
     &              'IT HAS BEEN RESET = RNDay.')
          NTO3DSTF=NT
        ENDIF

        NDSet3DST = (NTO3DSTF-NTO3DSTS)/NSpo3DST
        N3DST=0
        I3DSTRec=0

C.....  Write turbulence station output information to UNIT 16

        WRITE(16,573) NSpo3DST,TO3DSTS,(TO3DSTS+NSpo3DST*DTDP*Sec2Day), 
     &                (NTO3DSTS+NSpo3DST),TO3DSTF,NTO3DSTF,NSpo3DST,
     &                NDSet3DST
 573    FORMAT(/,5X,'3D Turbulence station output will start',I9,
     &              ' timesteps after day TO3DST = ',F9.3,
     &         /,9X,'This = ',F9.3,' days = ',I9,' timesteps relative',
     &              ' to the simulation Start Time.',
     &         /,9X,'Output will stop ',F9.3,' days = ',I9,' timesteps',
     &              ' relative to the simulation Start Time.',
     &         /,9X,'Data will be written to file fort.42 every ',I7,
     &              ' time steps at total of ',I7,' times.')
      ENDIF
 
C.... REGARDLESS OF WHETHER I3DST=0, READ IN THE NUMBER OF 3D Turbulence
C.... STATIONS

      READ(15,*) NSta3DT
      IF(I3DST.NE.0) WRITE(16,583) NSta3DT
 583  FORMAT(/,5X,'Output will be written at ',I7,' stations')

      IF (NSta3DT.EQ.0) THEN
         MNSta3DT=1
      ELSE
         MNSta3DT = NSta3DT
      ENDIF

C  Allocate arrays for station turbulence output

      call alloc_3DST()

C....Input the coordinates of the turbulence output stations
C....and then compute the element # containing each station

      IF(NSta3DT.GT.0) THEN
         IF(ICS.EQ.1) WRITE(16,3108)
         IF(ICS.EQ.2) WRITE(16,3109)

         DO N=1,NSta3DT
            NE3DT(N)=0
            IF(ICS.EQ.1) THEN
               READ(15,*) X3DS,Y3DS
            ELSE
               READ(15,*) SL3DS,SF3DS
               SL3DS=SL3DS*Deg2Rad
               SF3DS=SF3DS*Deg2Rad
               CALL CPP(X3DS,Y3DS,SL3DS,SF3DS,SLAM0,SFEA0)
            ENDIF

            CALL CoordinateToElement(X3DS, Y3DS,
     &           NE3DT(N), N, '3D Turbulence rec. station    ')
            
C.... Write element and station
            
            IF(ICS.EQ.1) THEN
               WRITE(16,1880) N,NE3DT(N),X3DS,Y3DS
            ELSE
               WRITE(16,1883) N,NE3DT(N),SL3DS*Rad2Deg,
     &              SF3DS*Rad2Deg,X3DS,Y3DS
            ENDIF
            
C.... COMPUTE INTERPOLATING FACTORS
            CALL ComputeInterpolatingFactors(X3DS, Y3DS, NE3DT(N),
     &           StaI3DT1(N), StaI3DT2(N), StaI3DT3(N)) 
            
         END DO
         
         WRITE(16,*) ' '
      ENDIF
      
C...3D Global Density Output (fort.44)
C
C     I3DGD    - flag for 3D global density output
C                If ABS(I3DGD)>0 global density is written out
C                If ABS(I3DGD)=1, output is in ascii. If ABS(I3DGD)=2, output is in binary.
C     TO3DGDS  - starting time in days for 3D global density output
C     TO3DGDF  - ending time in days for 3D global density output
C     NSpo3DGD - 3D global density output interval in timesteps
C     NTO3DGDS  - TO3DSDS converted to timesteps after StaTime
C     NTO3DGDF  - TO3DSDF converted to timesteps after StaTime
C     NDSet3DGD - total number of times that 3D global density data will be written
C     N3DGD    - counter to determine if it is time to write 3D global density data
C     I3DGDRec - record pointer for binary output
C
      READ(15,*) I3DGD,TO3DGDS,TO3DGDF,NSpo3DGD

C.... CHECK INPUT PARAMETER I3DGD

      WRITE(16,504) I3DGD
 504  FORMAT(//,1X,'3D DENSITY GLOBAL OUTPUT PARAMETER: I3DGD = ',I2)
      IF(ABS(I3DGD).GT.2) THEN
        IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
          WRITE(ScreenUnit,350)          
          WRITE(ScreenUnit,514)
        ENDIF
        WRITE(16,350)
        WRITE(16,514)
 514    FORMAT(/,1X,'YOUR SELECTION OF I3DGD IS NOT VALID. ',
     &              'CHECK YOUR INPUT!!')
        STOP
      ENDIF
 
      IF((I3DGD.NE.0).AND.(NSpo3DGD.EQ.0)) THEN
        IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
          WRITE(ScreenUnit,350)
          WRITE(ScreenUnit,564) NSpo3DGD
        ENDIF
        WRITE(16,350)
        WRITE(16,564) NSpo3DGD
 564    FORMAT(2X,'    NSpo3DGD = ',I3,' IT MUST BE > 0')
        STOP
      ENDIF

C.... IF NO GLOBAL DENSITY OUTPUT WILL BE GENERATED

      IF(I3DGD.EQ.0) THEN
        WRITE(16,524)
 524    FORMAT(/,5X,'No 3D global density will be output')

C.... IF GLOBAL DENSITY OUTPUT WILL BE GENERATED

      ELSE

C....   COMPUTE NTO3DGDS, NTO3DGDF, WHICH = TO3DGDS AND TO3DGDF IN TIME STEPS
#ifdef IBM
        NTO3DGDS=INT((TO3DGDS-StaTime)*Day2Sec/DTDP+0.5d0, 
     &        KIND(0.0d0))      !jgf45.11 was NINT
        NTO3DGDF=INT((TO3DGDF-StaTime)*Day2Sec/DTDP+0.5d0, 
     &       KIND(0.0d0))       !jgf45.11 was NINT
#else
        NTO3DGDS=INT((TO3DGDS-StaTime)*Day2Sec/DTDP+0.5d0) !jgf45.11 was NINT
        NTO3DGDF=INT((TO3DGDF-StaTime)*Day2Sec/DTDP+0.5d0) !jgf45.11 was NINT
#endif
        IF(NTO3DGDS.LT.0) THEN
          IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
             WRITE(ScreenUnit,534) TO3DGDS
          ENDIF
          WRITE(16,534) TO3DGDS
 534      FORMAT(5X,'NOTE: TO3DGDS = ',E14.6,' < START TIME. ',
     &              'IT HAS BEEN RESET = START TIME.')
          TO3DGDS=StaTime
          NTO3DGDS=0
        ENDIF
        IF(NTO3DGDF.LT.NTO3DGDS) THEN
          IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
             WRITE(ScreenUnit,544) TO3DGDF
          ENDIF
          WRITE(16,544) TO3DGDF
 544      FORMAT(5X,'NOTE: TO3DGDF = ',E14.6,' < TO3DGDS. ',
     &              'IT HAS BEEN RESET = TO3DGDS.')
          TO3DGDF=TO3DGDS
          NTO3DGDF=NTO3DGDS
        ENDIF
        IF(NTO3DGDF.GT.NT) THEN
          IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
             WRITE(ScreenUnit,554) TO3DGDF
          ENDIF
          WRITE(16,554) TO3DGDF
 554      FORMAT(5X,'NOTE: TO3DGDF = ',E14.6,' > RNDAY. ',
     &              'IT HAS BEEN RESET = RNDAY.')
          NTO3DGDF=NT
        ENDIF

        NDSet3DGD = (NTO3DGDF-NTO3DGDS)/NSpo3DGD
        N3DGD=0
        I3DGDRec=0

C.....  Write global 3D Density output information to UNIT 16

        WRITE(16,574) NSpo3DGD,TO3DGDS,(TO3DGDS+NSpo3DGD*DTDP*Sec2Day), 
     &                (NTO3DGDS+NSpo3DGD),TO3DGDF,NTO3DGDF,NSpo3DGD,
     &                NDSet3DGD
 574    FORMAT(/,5X,'3D Density global output will start',I9,
     &              ' timesteps after day TO3DGD = ',F9.3,
     &         /,9X,'This = ',F9.3,' days = ',I9,' timesteps relative',
     &              ' to the simulation Start Time.',
     &         /,9X,'Output will stop',F9.3,' days = ',I9,' timesteps',
     &              ' relative to the simulation Start Time.',
     &         /,9X,'Data will be written to file fort.42 every ',I7,
     &              ' time steps at total of ',I7,' times.')
 
        WRITE(16,*) ' '
      ENDIF

C...3D Global Velocity Output  (fort.45)
C
C     I3DGV    - flag for 3D global velocity output
C                If ABS(I3DGV)>0 global velocity is written out
C                If ABS(I3DGV)=1, output is in ascii. If ABS(I3DGV)=2, output is in binary.
C     TO3DGVS  - starting time in days for 3D global velocity output
C     TO3DGVF  - ending time in days for 3D global velocity output
C     NSpo3DGV - 3D global velocity output interval in timesteps
C     NTO3DGVS  - TO3DGVS converted to timesteps after StaTime
C     NTO3DGVF  - TO3DGVF converted to timesteps after StaTime
C     NDSet3DGV - total number of times that 3D velocity global data will be written
C     N3DGV    - counter to determine if it is time to write 3D velocity global data
C     I3DGVRec - record pointer for ascii or binary output
C
      READ(15,*) I3DGV,TO3DGVS,TO3DGVF,NSpo3DGV

C.... CHECK INPUT PARAMETER I3DGV

      WRITE(16,505) I3DGV
 505  FORMAT(//,1X,'3D Global Velocity Output Parameter: I3DGV = ',I2)
      IF(ABS(I3DGV).GT.2) THEN
        IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
          WRITE(ScreenUnit,350)          
          WRITE(ScreenUnit,515)
        ENDIF
        WRITE(16,350)
        WRITE(16,515)
 515    FORMAT(/,1X,'YOUR SELECTION OF I3DGV IS NOT VALID. ',
     &              'CHECK YOUR INPUT!!')
        STOP
      ENDIF
 
      IF((I3DGV.NE.0).AND.(NSpo3DGV.EQ.0)) THEN
        IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
          WRITE(ScreenUnit,350)
          WRITE(ScreenUnit,565) NSpo3DGV
        ENDIF
        WRITE(16,350)
        WRITE(16,565) NSpo3DGV
 565    FORMAT(2X,'    NSpo3DGV = ',I3,' IT MUST BE > 0')
        STOP
      ENDIF

C.... IF NO GLOBAL VELOCITY OUTPUT WILL BE GENERATED

      IF(I3DGV.EQ.0) THEN
        WRITE(16,525)
 525    FORMAT(/,5X,'No 3D Global velocity will be output')

C.... IF GLOBAL VELOCITY OUTPUT WILL BE GENERATED

      ELSE

C....   COMPUTE NTO3DGVS, NTO3DGVF, WHICH = TO3DGVS AND TO3DGVF IN TIME STEPS
#ifdef IBM
        NTO3DGVS=INT((TO3DGVS-StaTime)*Day2Sec/DTDP+0.5d0, 
     &        KIND(0.0d0))      !jgf45.11 was NINT
        NTO3DGVF=INT((TO3DGVF-StaTime)*Day2Sec/DTDP+0.5d0,
     &       KIND(0.0d0))       !jgf45.11 was NINT
#else
        NTO3DGVS=INT((TO3DGVS-StaTime)*Day2Sec/DTDP+0.5d0) !jgf45.11 was NINT
        NTO3DGVF=INT((TO3DGVF-StaTime)*Day2Sec/DTDP+0.5d0) !jgf45.11 was NINT
#endif

        IF(NTO3DGVS.LT.0) THEN
          IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
             WRITE(ScreenUnit,535) TO3DGVS
          ENDIF
          WRITE(16,535) TO3DGVS
 535      FORMAT(5X,'NOTE: TO3DGVS = ',E14.6,' < START TIME. ',
     &              'IT HAS BEEN RESET = START TIME.')
          TO3DGVS=StaTime
          NTO3DGVS=0
        ENDIF
        IF(NTO3DGVF.LT.NTO3DGVS) THEN
          IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
             WRITE(ScreenUnit,545) TO3DGVF
          ENDIF
          WRITE(16,545) TO3DGVF
 545      FORMAT(5X,'NOTE: TO3DGVF = ',E14.6,' < TO3DGVS. ',
     &              'IT HAS BEEN RESET = TO3DGVS.')
          TO3DGVF=TO3DGVS
          NTO3DGVF=NTO3DGVS
        ENDIF
        IF(NTO3DGVF.GT.NT) THEN
          IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
             WRITE(ScreenUnit,555) TO3DGVF
          ENDIF
          WRITE(16,555) TO3DGVF
 555      FORMAT(5X,'NOTE: TO3DGVF = ',E14.6,' > RNDAY. ',
     &              'IT HAS BEEN RESET = RNDAY.')
          NTO3DGVF=NT
        ENDIF

        NDSet3DGV = (NTO3DGVF-NTO3DGVS)/NSpo3DGV
        N3DGV=0
        I3DGVRec=0

C.....  Write global velocity output information to UNIT 16

        WRITE(16,575) NSpo3DGV,TO3DGVS,(TO3DGVS+NSpo3DGV*DTDP*Sec2Day), 
     &                (NTO3DGVS+NSpo3DGV),TO3DGVF,NTO3DGVF,NSpo3DGV,
     &                NDSet3DGV
 575    FORMAT(/,5X,'3D Global Velocity output will start',I9,
     &              ' timesteps after day TO3DGV = ',F9.3,
     &         /,9X,'This = ',F9.3,' days = ',I9,' timesteps relative',
     &              ' to the simulation Start Time.',
     &         /,9X,'Output will stop ',F9.3,' days = ',I9,' timesteps',
     &              ' relative to the simulation Start Time.',
     &         /,9X,'Data will be written to file fort.42 every ',I7,
     &              ' time steps at total of ',I7,' times.')

        WRITE(16,*) ' '
      ENDIF

C...3D Global Turbulence Output  (fort.46)
C
C     I3DGT    - flag for 3D global turbulence output
C                If ABS(I3DGT)>0 global turbulence is written out
C                If ABS(I3DGT)=1, output is in ascii. If ABS(I3DGT)=2, output is in binary.
C     TO3DGTS  - starting time in days for 3D global turbulence output
C     TO3DGTF  - ending time in days for 3D global turbulence output
C     NSpo3DGT - 3D global turbulence output interval in timesteps
C     NTO3DGTS  - TO3DGTS converted to timesteps after StaTime
C     NTO3DGTF  - TO3DGTF converted to timesteps after StaTime
C     NDSet3DGT - total number of times that 3D global turbulence data will be written
C     N3DGT    - counter to determine if it is time to write 3D global turbulence data
C     I3DGTRec - record pointer for ascii or binary output
C
      READ(15,*) I3DGT,TO3DGTS,TO3DGTF,NSpo3DGT

C.... CHECK INPUT PARAMETER I3DGT

      WRITE(16,506) I3DGT
 506  FORMAT(//,1X,'3D Global Turbulence output parameter: I3DGT = ',I2)
      IF(ABS(I3DGT).GT.2) THEN
        IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
          WRITE(ScreenUnit,350)          
          WRITE(ScreenUnit,516)
        ENDIF
        WRITE(16,350)
        WRITE(16,516)
 516    FORMAT(/,1X,'Your selection of I3DGT is not valid. ',
     &              'CHECK YOUR INPUT!!')
#ifdef CMPI
        call msg_fini()
#endif
        STOP
      ENDIF

      IF((I3DGT.NE.0).AND.(NSpo3DGT.EQ.0)) THEN
        IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
          WRITE(ScreenUnit,350)
          WRITE(ScreenUnit,566) NSpo3DGT
        ENDIF
        WRITE(16,350)
        WRITE(16,566) NSpo3DGT
 566    FORMAT(2X,'    NSpo3DGT = ',I3,' IT MUST BE > 0')
#ifdef CMPI
        call msg_fini()
#endif
        STOP
      ENDIF

C.... IF NO 3D GLOBAL TURBULENCE OUTPUT WILL BE GENERATED

      IF(I3DGT.EQ.0) THEN
        WRITE(16,526)
 526    FORMAT(/,5X,'No 3D Global Turbulence will be output')

C.... IF 3D GLOBAL TURBULENCE OUTPUT WILL BE GENERATED

      ELSE

C....   COMPUTE NTO3DGTS, NTO3DGTF, WHICH = TO3DGTS AND TO3DGTF IN TIME STEPS
#ifdef IBM
        NTO3DGTS=INT((TO3DGTS-StaTime)*Day2Sec/DTDP+0.5d0,
     &        KIND(0.0d0))      !jgf45.11 was NINT
        NTO3DGTF=INT((TO3DGTF-StaTime)*Day2Sec/DTDP+0.5d0, 
     &       KIND(0.0d0))       !jgf45.11 was NINT
#else
        NTO3DGTS=INT((TO3DGTS-StaTime)*Day2Sec/DTDP+0.5d0) !jgf45.11 was NINT
        NTO3DGTF=INT((TO3DGTF-StaTime)*Day2Sec/DTDP+0.5d0) !jgf45.11 was NINT
#endif

        IF(NTO3DGTS.LT.0) THEN
          IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
             WRITE(ScreenUnit,536) TO3DGTS
          ENDIF
          WRITE(16,536) TO3DGTS
 536      FORMAT(5X,'NOTE: TO3DGTS = ',E14.6,' < START TIME. ',
     &              'IT HAS BEEN RESET = START TIME.')
          TO3DGTS=StaTime
          NTO3DGTS=0
        ENDIF
        IF(NTO3DGTF.LT.NTO3DGTS) THEN
          IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
             WRITE(ScreenUnit,546) TO3DGTF
          ENDIF
          WRITE(16,546) TO3DGTF
 546      FORMAT(5X,'NOTE: TO3DGTF = ',E14.6,' < TO3DGTS. ',
     &              'IT HAS BEEN RESET = TO3DGTS.')
          TO3DGTF=TO3DGTS
          NTO3DGTF=NTO3DGTS
        ENDIF
        IF(NTO3DGTF.GT.NT) THEN
          IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
             WRITE(ScreenUnit,556) TO3DGTF
          ENDIF
          WRITE(16,556) TO3DGTF
 556      FORMAT(5X,'NOTE: TO3DGTF = ',E14.6,' > RNDay. ',
     &              'IT HAS BEEN RESET = RNDay.')
          NTO3DGTF=NT
        ENDIF

        NDSet3DGT = (NTO3DGTF-NTO3DGTS)/NSpo3DGT
        N3DGT=0
        I3DGTRec=0

C.....  Write global turbulence output information to UNIT 16

        WRITE(16,576) NSpo3DGT,TO3DGTS,(TO3DGTS+NSpo3DGT*DTDP*Sec2Day), 
     &                (NTO3DGTS+NSpo3DGT),TO3DGTF,NTO3DGTF,NSpo3DGT,
     &                NDSet3DGT
 576    FORMAT(/,5X,'3D Global Turbulence output will start',I9,
     &              ' timesteps after day TO3DGT = ',F9.3,
     &         /,9X,'This = ',F9.3,' days = ',I9,' timesteps relative',
     &              ' to the simulation Start Time.',
     &         /,9X,'Output will stop ',F9.3,' days = ',I9,' timesteps',
     &              ' relative to the simulation Start Time.',
     &         /,9X,'Data will be written to file fort.42 every ',I7,
     &              ' time steps at total of ',I7,' times.')

        WRITE(16,*) ' '
      ENDIF

!     Kendra45.12: Read in the new input values for the transport equation
C     jgf45.12: Made READs conditional on value of C3D_BTrans.

      if (C3D_Btrans) then
         
!     Kendra45.12: Must read in new values for lateral and vertical
!     diffusion
C...  READ IN NLSD, NLTD, NVTD & NVSD: Lateral and vertical diffusion
!     coefficients.
         READ(15,*) NLSD, NVSD
         WRITE(16,416) NLSD, NVSD
 416     FORMAT(/,5X,'Salinity Lateral Diffusion Coefficient = ',E8.2,
     &        'Salinity Vertical Diffusion Coefficient = ',E8.2)
         
         READ(15,*) NLTD, NVTD
         WRITE(16,417) NLTD, NVTD
 417     FORMAT(/,5X,'Temperature Lateral Diffusion Coefficient = ',
     &        E8.2,'Temperature Vertical Diffusion Coefficient = ',E8.2)
         
!     Kendra45.12: Read in the time stepping coefficient associated with the
!     transport equation terms.
         READ(15,*) ALP4
         WRITE(16,445) ALP4
 445     FORMAT(/,5X,'3D TIME STEPPING COEFFS ALP4 = ',E8.2)
         
         DTAlp4 = DelT*Alp4
         DT1MAlp4 = DelT*(1-Alp4)

!     Kendra45.12: Read in the temperature boundary condition file type
c     jgf45.12: Made READ conditional on dynamic temperature forcing.
         if ( IDEN .eq. 3 .or. IDEN .eq. 4 ) then
            READ(15,*) NTF
            WRITE(16,444) NTF
 444        FORMAT(/,5X,'Temperature flux conditions are ', I7)
         endif

      endif
      
      RETURN
C-----------------------------------------------------------------------
      END SUBROUTINE READ_INPUT_3D
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     S U B R O U T I N E   C O O R D I N A T E  T O  E L E M E N T 
C-----------------------------------------------------------------------
C
C     jgf45.12 Subroutine to take an X and Y cartesian coordinate and
C     find the corresponding element.
C
C-----------------------------------------------------------------------
      SUBROUTINE CoordinateToElement(InputXCoordinate, InputYCoordinate,
     &     OutputElement, StationNumber, Description)
      USE SIZES, ONLY : SZ, MyProc
      USE GLOBAL, ONLY: NE, NM, X, Y, Areas, NFOver, NScreen, ScreenUnit
#ifdef CMPI
      USE MESSENGER, ONLY : msg_fini
#endif
      IMPLICIT NONE

      REAL(8), intent(in) :: InputXCoordinate                  ! cartesian
      REAL(8), intent(in) :: InputYCoordinate                  ! cartesian
      INTEGER, intent(out) :: OutputElement                    
      INTEGER, intent(in) :: StationNumber                     ! for err. mesg.
      CHARACTER(len=30), intent(in) :: Description             ! for err. mesg.

      INTEGER Element         ! element loop counter
      INTEGER ClosestElement  ! element with closest match
      REAL(8) X1, X2, X3, X4, Y1, Y2, Y3, Y4       ! geometry
      REAL(8) A1, A2, A3, AE, AEMIN, AA            ! area
      LOGICAL ElementFound  ! .true. when a corresponding element is found
      REAL(8), PARAMETER :: Tolerance = 1.0E-5     ! area difference for match 

      ElementFound = .false.
      AEMIN=1.0E+25
      ClosestElement=0
      DO Element=1,NE
         X1=X(NM(Element,1))
         X2=X(NM(Element,2))
         X3=X(NM(Element,3))
         X4=InputXCoordinate
         Y1=Y(NM(Element,1))
         Y2=Y(NM(Element,2))
         Y3=Y(NM(Element,3))
         Y4=InputYCoordinate
         A1=(X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4)
         A2=(X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1)
         A3=(Y4-Y1)*(X2-X1)-(X4-X1)*(Y2-Y1)
         AA=ABS(A1)+ABS(A2)+ABS(A3)
         AE=ABS(AA-Areas(Element))/Areas(Element)
         IF (AE.LT.AEMIN) THEN
            AEMIN=AE
            ClosestElement=Element
         ENDIF
         IF (AE.LT.Tolerance) THEN
            ElementFound = .true.
            OutputElement=Element
         ENDIF
      ENDDO
      IF (.not. ElementFound ) THEN
         IF((NScreen.NE.0).AND.(MyProc.EQ.0)) THEN
            WRITE(ScreenUnit,593) StationNumber
         ENDIF
         WRITE(16,593) Description, StationNumber
         IF(NFOVER.EQ.1) THEN
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
               WRITE(ScreenUnit,9790) AEMIN
            ENDIF
            WRITE(16,9790) AEMIN
            OutputElement = ClosestElement
         ELSE
            IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
               WRITE(ScreenUnit,9791) AEMIN
            ENDIF
            WRITE(16,9791) AEMIN
#ifdef CMPI
            call msg_fini()
#endif
            STOP
         ENDIF
      ENDIF

 593  FORMAT(///,1X,'!!!!!!!!!!  WARNING - NONFATAL ',
     &     'INPUT ERROR  !!!!!!!!!',//
     &     ,1X,A30,1X,I6,' does ',
     &     'not lie within any element in the defined',
     &     /,1X,'computational domain.   PLEASE CHECK THE ',
     &     'INPUT COORDINATES FOR THIS STATION')
 9790 FORMAT(/,1X,'PROGRAM WILL ESTIMATE NEAREST ELEMENT',
     &     /,1X,'PROXIMITY INDEX FOR THIS STATION EQUALS ',E15.6,
     &     //,1X,'!!!!!! EXECUTION WILL CONTINUE !!!!!!',//)
 9791 FORMAT(/,1X,'PROGRAM WILL NOT CORRECT ERROR ',
     &     'SINCE NON-FATAL ERROR OVERIDE OPTION, NFOVER,',
     &     /,1X,'HAS BEEN SELECTED EQUAL TO 0',
     &     /,1X,'PROXIMITY INDEX FOR THIS STATION EQUALS ',E15.6,
     &     //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',
     &     //)

      RETURN
C-----------------------------------------------------------------------
      END SUBROUTINE CoordinateToElement
C-----------------------------------------------------------------------



C-----------------------------------------------------------------------
C     S U B R O U T I N E   
C     C O M P U T E  I N T E R P O L A T I N G  F A C T O R S 
C-----------------------------------------------------------------------
C
C     jgf45.12 Subroutine to pre-compute the interpolating factors for a
C     recording station.
C
C-----------------------------------------------------------------------
      SUBROUTINE ComputeInterpolatingFactors(InputXCoordinate,
     &     InputYCoordinate, InputElement, Factor1, Factor2, Factor3)

      USE GLOBAL, ONLY : NM, X, Y, Areas
      IMPLICIT NONE
      REAL(8), intent(in) :: InputXCoordinate                  ! cartesian
      REAL(8), intent(in) :: InputYCoordinate                  ! cartesian
      INTEGER, intent(in) :: InputElement                    
      REAL(8), intent(out):: Factor1, Factor2, Factor3         

      REAL(8) X1, X2, X3, X4, Y1, Y2, Y3, Y4                   ! geometry

      X1=X(NM(InputElement,1))
      X2=X(NM(InputElement,2))
      X3=X(NM(InputElement,3))
      X4=InputXCoordinate
      Y1=Y(NM(InputElement,1))
      Y2=Y(NM(InputElement,2))
      Y3=Y(NM(InputElement,3))
      Y4=InputYCoordinate

      Factor1=((X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4))/Areas(InputElement)
      Factor2=((X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1))/Areas(InputElement)
      Factor3=(-(X4-X1)*(Y2-Y1)+(Y4-Y1)*(X2-X1))/Areas(InputElement)
      
      RETURN
C-----------------------------------------------------------------------
      END SUBROUTINE ComputeInterpolatingFactors
C-----------------------------------------------------------------------

      !=================================================================
      !=================================================================
      !=================================================================
      !      =====                                            =====
      !      =====            SUBROUTINE PARSE                =====
      !      =====                                            =====
      !=================================================================
      !=================================================================
      !=================================================================

      !=================================================================
      ! This subroutines parse station names so they can passed on to
      ! the NetCDF I/O module of ADCIRC and be written out in the netCDF
      ! output files.
      !
      ! Revision history:
      !
      ! Date      Programmer                       Description of change
      ! ----      ----------                       ---------------------
      ! 12/16/07  Cristina Forbes, PSGS @ UNC-IMS  Wrote parser
      !                                            based on character
      !                                            manipulation subroutines
      !                                            developed at the
      !                                            University of Oslo, Norway
      !=================================================================
      subroutine parse (line, lvar)
      CHARACTER(132) LINE,L
      CHARACTER(50) LVAR(3)
      INTEGER LS, L1, LEN1

      LS=0
      L=TRIM(ADJUSTL(LINE))
      LEN1=LEN(L)

!     EXTRACT STATION DATA FROM LINE

      DO I=1,2
!        DO I=1,3
         L1=LEN(TRIM(ADJUSTL(L(1:INDEX(L," ")-1))))
         LVAR(I)=ADJUSTL(L(LS:L1))
         L=ADJUSTL(L(L1+1:LEN(L)))
      END DO
      LVAR(3)=L
      END SUBROUTINE parse


      FUNCTION front_trim(buf) RESULT(res)
      IMPLICIT NONE
      !// Arguments
      CHARACTER(LEN=*), INTENT(IN)   :: buf
      INTEGER                        :: res
      !// Local variables
      INTEGER                        :: i, lng
      lng = LEN_TRIM(buf)
      DO i = 1, lng
         IF (buf(i:i) .NE. ' ') THEN
            res = i
            RETURN
         END IF
      END DO
      END FUNCTION front_trim

 

      FUNCTION a2i(buf) RESULT(res)
      IMPLICIT NONE
      !// Arguments
      CHARACTER(LEN=*), INTENT(IN)   :: buf
      INTEGER                        :: res
      !// Local variables
      INTEGER                        :: i, foffs, lng
      INTEGER, EXTERNAL              :: front_trim
      LOGICAL                        :: neg
      neg = .FALSE.
      lng = LEN_TRIM(buf)
      foffs = front_trim(buf)
      res = 0
      DO i = foffs, lng
         IF(buf(i:i) .EQ. '-') THEN
            neg = .TRUE.
           CONTINUE
         END IF
         IF(buf(i:i) .GE. '0' .AND. buf(i:i) .LE. '9' .AND.
     &        buf(i:i) .NE. ' ') THEN
            res = res * 10
            res = res + (IACHAR(buf(i:i)) - 48)
         END IF
      END DO
      IF(neg) THEN
         res = res * (-1)
      END IF
      END FUNCTION a2i


      FUNCTION a2f(buf) RESULT(res)
      IMPLICIT NONE
      !// Arguments
      CHARACTER(LEN=*)               :: buf
      REAL                           :: res
      !// Local variables
      INTEGER                        :: i,j,k,itmp, foffs, lng
      INTEGER                        :: p, q
      REAL                           :: rtmp
      INTEGER, EXTERNAL              :: front_trim
      INTEGER, EXTERNAL              :: a2i
      LOGICAL                        :: neg
      LOGICAL                        :: exponential
      neg = .FALSE.
      exponential = .FALSE.
      lng = LEN_TRIM(buf)
      foffs = front_trim(buf)
      k = INDEX(buf,'.')
      p = INDEX(buf,'e')
      q = INDEX(buf,'E')
      IF(p /= -1 .AND. p /= 0) THEN
         exponential = .TRUE.
      END IF
      IF(q /= -1 .AND. q /= 0) THEN
         exponential = .TRUE.
         p = q
      END IF
      IF (k /= -1 .AND. k /= 0) THEN
         !// We have a floating point number
         itmp = 0
         !// Get the integer part of the number
         DO i = foffs, k - 1
            IF(buf(i:i) .EQ. '-') THEN
               neg = .TRUE.
               CONTINUE
            END IF
            IF(buf(i:i) .GE. '0' .AND. buf(i:i) .LE. '9' .AND.
     &           buf(i:i) .NE. ' ') THEN
               itmp = itmp * 10
               itmp = itmp + (IACHAR(buf(i:i)) - 48)
            END IF
         END DO
         res = DBLE(itmp)
         rtmp = 0.
         q = 0
         IF(.NOT. exponential) THEN
!// We do not have an exponential number
            DO i = LEN_TRIM(buf), k+1, -1
               itmp = (IACHAR(buf(i:i)) - 48)
               rtmp = rtmp + FLOAT(itmp)
               rtmp = rtmp / 10.
            END DO
         ELSE
!// We have an exponential number
            DO i = p-1, k+1, -1
               itmp = (IACHAR(buf(i:i)) - 48)
               rtmp = rtmp + FLOAT(itmp)
               rtmp = rtmp / 10.
            END DO
            q = a2i(buf(p+1:LEN_TRIM(buf)))
         END IF
         res = res + rtmp
         IF(exponential) THEN
            res = res * 10.**q
         END IF
         IF(neg) THEN
            res = res * (-1)
         END IF
      ELSE
         !// We have an integer
         res = REAL(a2i(buf))
      END IF
      END FUNCTION a2f


      FUNCTION a2d(buf) RESULT(res)
      IMPLICIT NONE
      !// Arguments
        CHARACTER(LEN=*)               :: buf
        DOUBLE PRECISION               :: res
!// Local variables
        INTEGER                        :: i,j,k,itmp, foffs, lng
        INTEGER                        :: p, q
        DOUBLE PRECISION               :: dtmp
        INTEGER, EXTERNAL              :: front_trim
        INTEGER, EXTERNAL              :: a2i
        LOGICAL                        :: neg
        LOGICAL                        :: exponential
        neg = .FALSE.
        exponential = .FALSE.
        lng = LEN_TRIM(buf)
        foffs = front_trim(buf)
        k = INDEX(buf,'.')
        p = INDEX(buf,'e')
        q = INDEX(buf,'E')
        IF(p /= -1 .AND. p /= 0) THEN
        exponential = .TRUE.
        END IF
        IF(q /= -1 .AND. q /= 0) THEN
        exponential = .TRUE.
        p = q
        END IF
        IF (k /= -1 .AND. k /= 0) THEN
!// We have a floating point number
        itmp = 0
!// Get the integer part of the number
        DO i = foffs, k - 1
        IF(buf(i:i) .EQ. '-') THEN
        neg = .TRUE.
        CONTINUE
        END IF
        IF(buf(i:i) .GE. '0' .AND. buf(i:i) .LE. '9' .AND.
     &  buf(i:i) .NE. ' ') THEN
        itmp = itmp * 10
        itmp = itmp + (IACHAR(buf(i:i)) - 48)
        END IF
        END DO
        res = DBLE(itmp)
        dtmp = 0.
        q = 0
        IF(.NOT. exponential) THEN
!// We do not have an exponential number
        DO i = LEN_TRIM(buf), k+1, -1
        itmp = (IACHAR(buf(i:i)) - 48)
        dtmp = dtmp + FLOAT(itmp)
        dtmp = dtmp / 10.
        END DO
        ELSE
!// We have an exponential number
        DO i = p-1, k+1, -1
        itmp = (IACHAR(buf(i:i)) - 48)
        dtmp = dtmp + FLOAT(itmp)
        dtmp = dtmp / 10.
        END DO
        q = a2i(buf(p+1:LEN_TRIM(buf)))
        END IF
        res = res + dtmp
        IF(exponential) THEN
        res = res * 10.**q
        END IF
        IF(neg) THEN
        res = res * (-1)
        END IF
        ELSE
!// We have an integer
        res = DBLE(a2i(buf))
        END IF
        END FUNCTION a2d



        SUBROUTINE i2a(int, res)

        IMPLICIT NONE
!// Arguments
        INTEGER, INTENT(IN)            :: int
        CHARACTER(LEN=*), INTENT(OUT)  :: res
!// Local variables
        INTEGER                        :: i,j,k
        CHARACTER(LEN=80)              :: sbuf
        sbuf = ' '
        res = ' '
        i = int
        k = 1
        DO
        j = MOD(i,10)
        sbuf(k:k) = ACHAR(j+48)
        i = i / 10
        IF (i <= 0) THEN
        EXIT
        END IF
        k = k + 1
        END DO
        k = LEN_TRIM(sbuf)
        IF (k>1) THEN
        j = k
        i = 1
        DO i = 1, k
         res(i:i) = sbuf(j:j)
         j = j - 1
        END DO
        ELSE
        res(1:1) = sbuf(1:1)
        END IF

        END SUBROUTINE i2a



