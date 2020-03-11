c***********************************************************************************
c  PROGRAM  P15
c
c  P15 will produce ADCIRC wind file for at most 100 track records which
c  are evenly spaced in time.  It will read your grid file an
c  automatically allocate the necessary amount of memory, so you do NOT
c  need to edit this file.  Also it will automatically read your track
c  file and make the necessary time interval adjustment.      
c  vjp 3/15/2006
c
c***********************************************************************************
c
c Modified to multiply wind speeds by 1.04 to convert 30-min winds to 10-min winds
c   -- Use this program with adcirc v46.xx or later
c   -- 1.04 was originally multiplied in interp.F of adcprep, but it will be 
c      removed in adcirc v46 series.  sb 8/15/2006
c
c***********************************************************************************
c
c  Program to read files with hurricane information in them from the website:  
c  http://wxp.eas.purdue.edu/hur_atlantic.  The program uses this data to prepare a 
c  wind file to be used with ADCIRC.  This program bypasses having to create a specific
c  input file to use with pblinput.f and then using pblinput and finally running a p 
c  (ex. p4.f) program to get the ADCIRC format.  Any time increment can be used and out-
c  put is every hour.
c                       
c  Modified from pblinput5.f and p4.f                       
c    p4.f      r.l. 7/11/96
c    p5_2.f    c.f. 3/30/98
c    p7.f      c.f. 10/29/98
c
c Modified from p7_04_29c.f by j.f. 12/19/01
c
c Modified from p08a by jjw         03/17/2003
c   Added wind reduction factor input file that asks for a grid file with 
c   wind reduction factors at each node.
c
c Modified from p08a by jjw         05/28/2003
c   Modified speed by 3x to allow 2 hourly input with faked times (still
c   6 hourly) so that will put out 20 min snaps (although PBL will
c   think of it as hourly)
c
c Modified further on 2003/05/28 so that better radius to max wind is computed
c   This has been done by adding Ed Thompson's subroutine JTrmax
c
c Eliminated any land reduction factor  - jjw 20030619
c modified from p08b_031703c_3xstormspeed_newrtmw.f
c
c p10 added hyperbolic tangent ramp
c p11 and p11 bug corrections on ramping function
c
c jgf47.01 Modified from v46.57 to change the name of MODULE SIZES
c to MODULE P15EXTENTS to avoid collision or confusion with the ADCIRC
c module named sizes.
C 
c***********************************************************************************

      MODULE P15EXTENTS
      use CkLunMod, only : CK_LUN
      INTEGER NBIG
!$omp threadprivate(NBIG)
      INTEGER NE, NBR
!$omp threadprivate(NE,NBR)
      END MODULE P15EXTENTS

      MODULE P1
      use CkLunMod, only : CK_LUN
      USE P15EXTENTS
      INTEGER ITOT
!$omp threadprivate(ITOT)
      REAL,ALLOCATABLE :: XXLAT(:),XXLON(:),XTIME(:),XXVEL(:),XXPRES(:)
!$omp threadprivate(XXLAT,XXLON,XTIME,XXVEL,XXPRES)
      END MODULE P1

      MODULE C0
      use CkLunMod, only : CK_LUN
      INTEGER  IMO,IDY,IYR,ISTART,NPTS
!$omp threadprivate(IMO,IDY,IYR,ISTART,NPTS)
      REAL IHR1,TIMEINC
!$omp threadprivate(IHR1,TIMEINC)
      END MODULE C0

      MODULE C1
      use CkLunMod, only : CK_LUN
      INTEGER  NSNAP,NZ,NM,IB
!$omp threadprivate(NSNAP,NZ,NM,IB)
      REAL AN1(600)
!$omp threadprivate(AN1)
      REAL DX,DT,F,SGW,UC,VC,UG,VG,CS,ST12
!$omp threadprivate(DX,DT,F,SGW,UC,VC,UG,VG,CS,ST12)
      END MODULE C1

      MODULE C2
      use CkLunMod, only : CK_LUN
      REAL,TARGET  :: AB(21,21,15),AC(21,7,21)                       
!$omp threadprivate(AB,AC)
      REAL,TARGET  ::  JA(15)
!$omp threadprivate(JA)
      REAL,POINTER ::  ITRACK    ! Map of JA
!$omp threadprivate(ITRACK)
      REAL,POINTER ::  EYELAT
!$omp threadprivate(EYELAT)
      REAL,POINTER ::  EYELONG 
!$omp threadprivate(EYELONG)
      REAL,POINTER ::  DIREC 
!$omp threadprivate(DIREC)
      REAL,POINTER ::  SPEED 
!$omp threadprivate(SPEED)
      REAL,POINTER ::  IQUAD 
!$omp threadprivate(IQUAD)
      REAL,POINTER ::  EYPRES 
!$omp threadprivate(EYPRES)
      REAL,POINTER,DIMENSION(:) ::  RADIUS, PFAR
!$omp threadprivate(RADIUS,PFAR)
      END MODULE C2

      MODULE C3
      use CkLunMod, only : CK_LUN
      REAL  U(21,21,5),V(21,21,5),UN(21,21,5),VN(21,21,5)
!$omp threadprivate(U,V,UN,VN)
      REAL  PX(21,21,5),PY(21,21,5),VTN(21,21,5),ANG(21,21,5)
!$omp threadprivate(PX,PY,VTN,ANG)
      REAL  LW(21,21,5)
!$omp threadprivate(LW)
      END MODULE C3

      MODULE C4
      use CkLunMod, only : CK_LUN
      REAL  CDR(100,2,2),UXV(100,2),TURN(100,2)    
!$omp threadprivate(CDR,UXV,TURN)
      END MODULE C4

      MODULE C5
      use CkLunMod, only : CK_LUN
      INTEGER K123,LS
!$omp threadprivate(K123,LS)
      REAL HH,HL,K2,K35
!$omp threadprivate(HH,HL,K2,K35)
      REAL FLAT,PTH,DTH(2),Z0LAND,VV(100),UX(3),UV(3)   
!$omp threadprivate(FLAT,PTH,DTH,Z0LAND,VV,UX,UV)
      REAL DUV(3),G,GA,DEN,VV2,Z0,ZLOG,AM,BM,CM,FF
!$omp threadprivate(DUV,G,GA,DEN,VV2,Z0,ZLOG,AM,BM,CM,FF)
      END MODULE C5

      MODULE C6
      use CkLunMod, only : CK_LUN
      USE P15EXTENTS
      REAL PFAR1
!$omp threadprivate(PFAR1)
      REAL RADIUS1(600)
!$omp threadprivate(RADIUS1)
      REAL,ALLOCATABLE :: ZLA(:),ZLO(:),ZPR(:),ZAN(:),ZSP(:)
!$omp threadprivate(ZLA,ZLO,ZPR,ZAN,ZSP)
      END MODULE C6

      MODULE C8
      use CkLunMod, only : CK_LUN
      USE P15EXTENTS
      INTEGER,ALLOCATABLE :: ILAD(:),ILAM(:),ILOD(:),ILOM(:)
!$omp threadprivate(ILAD,ILAM,ILOD,ILOM)
      INTEGER,ALLOCATABLE :: ISNAP1(:),ISNAP2(:),IROT(:)
!$omp threadprivate(ISNAP1,ISNAP2,IROT)
      REAL,ALLOCATABLE :: PCT1(:)
!$omp threadprivate(PCT1)
      END MODULE C8


      MODULE D13
      use CkLunMod, only : CK_LUN
      INTEGER LANSEA
!$omp threadprivate(LANSEA)
      REAL LA0,LO0,LA1,LO1
!$omp threadprivate(LA0,LO0,LA1,LO1)
      REAL ROT,DX,STHT,W1,TH1,D,AL
!$omp threadprivate(ROT,DX,STHT,W1,TH1,D,AL)
      REAL,TARGET :: UST   
!$omp threadprivate(UST)
      REAL XX(21,21,10)
!$omp threadprivate(XX)
      INTEGER NSNAP1(600),NSNAP2(600),IPI(600),NHT,INTVN,INTVI
!$omp threadprivate(NSNAP1,NSNAP2,IPI,NHT,INTVN,INTVI)
      REAL PCT(600)
!$omp threadprivate(PCT)
      END MODULE D13

      MODULE DJM
      use CkLunMod, only : CK_LUN
      INTEGER NTIMES
!$omp threadprivate(NTIMES)
      CHARACTER*4 SMF_NAM(100)
!$omp threadprivate(SMF_NAM)
      REAL SMF_UN(21,21,5,100), SMF_VN(21,21,5,100)
!$omp threadprivate(SMF_UN,SMF_VN)
      REAL ISMF_JA(15,100)
!$omp threadprivate(ISMF_JA)
      REAL SMF_DXX(100), SMF_SGW(100), SMF_AN1(100)
!$omp threadprivate(SMF_DXX,SMF_SGW,SMF_AN1)
      REAL SMF_ST1(100), SMF_HHH(100), SMF_GAR(100)
!$omp threadprivate(SMF_ST1,SMF_HHH,SMF_GAR)
      REAL SMF_PTH(100), SMF_K35(100), SMF_DTH1(100), SMF_DTH2(100)
!$omp threadprivate(SMF_PTH,SMF_K35,SMF_DTH1,SMF_DTH2)
      END MODULE DJM

      MODULE C57
      use CkLunMod, only : CK_LUN
      INTEGER K123
!$omp threadprivate(K123)
      REAL PTH,DTH,HH,ZCOEFF(3,5),LAKE,VV(100),UX(3),UV(3)     
!$omp threadprivate(PTH,DTH,HH,ZCOEFF,LAKE,VV,UX,UV)
      REAL DUV(3),K35,K2,G,GA,DEN,VV2,HL,Z0,ZLOG,AM,BM,CM
!$omp threadprivate(DUV,K35,K2,G,GA,DEN,VV2,HL,Z0,ZLOG,AM,BM,CM)
      REAL UXV(100,6),TURN(100,6),COST(100),SINT(100) 
!$omp threadprivate(UXV,TURN,COST,SINT)
      END MODULE C57

      MODULE P15GLOBAL
      use CkLunMod, only : CK_LUN
      USE P15EXTENTS
      ! WAS BLANK COMMON
      REAL U99,V99,PRESS
!$omp threadprivate(U99,V99,PRESS)
      REAL P0SNP(100),PFRSNP(100),RADSNP(100)
!$omp threadprivate(P0SNP,PFRSNP,RADSNP)
      INTEGER,ALLOCATABLE :: KSTA(:)
!$omp threadprivate(KSTA)
      REAL,ALLOCATABLE :: YLA(:),YLO(:),WINDRED(:,:)
!$omp threadprivate(YLA,YLO,WINDRED)
      END MODULE P15GLOBAL

      MODULE IHOUR
      use CkLunMod, only : CK_LUN
      USE P15EXTENTS
      INTEGER,ALLOCATABLE :: IHR(:)
!$omp threadprivate(IHR)
      END MODULE IHOUR


      PROGRAM P15        
      USE P15EXTENTS
      USE P1; USE C0; USE C1; USE C2; USE C3; USE C4; USE C5
      USE C6; USE C8; USE DJM; USE P15GLOBAL
      CHARACTER*30 TITLE,INPUTFILE,GRIDFILE
!$omp threadprivate(TITLE,INPUTFILE,GRIDFILE)
      CHARACTER*9  TSTAMP(200)
!$omp threadprivate(TSTAMP)
      REAL A(500,4)
!$omp threadprivate(A)
      INTEGER ISPD(500),IPRES(500)
!$omp threadprivate(ISPD,IPRES)
      REAL TIME(500),XLAT(500),XLON(500)
!$omp threadprivate(TIME,XLAT,XLON)
      REAL XPRES(500),XVEL(500)
!$omp threadprivate(XPRES,XVEL)
      REAL RZERO(600)
!$omp threadprivate(RZERO)
      REAL SPEEDUP
!$omp threadprivate(SPEEDUP)
      REAL,ALLOCATABLE :: ANGLE(:),SPED(:)
!$omp threadprivate(ANGLE,SPED)

      WRITE(*,*) 'Enter name of input file to use:  '
      read(*,'(A)') inputfile
      write(*,*) 'Enter name of grid file to use:  '
      read(*,'(A)') gridfile
      write(*,*) 'Enter the year of the storm (xxxx):  '
      read(*,*) iyr
      write(*,*) 'Enter the month of the storm (ex. 8 = August):  '
      read(*,*) imo
      write(*,*) 'Enter the day of the month storm starts:  '
      read(*,*) idy

      open(10+CK_LUN,file=gridfile,STATUS='unknown')
      READ(10+CK_LUN,'(A)') TITLE
      READ(10+CK_LUN,*)  NE, NBR
      NBIG = INT(1.25*REAL(NBR))

      CALL ALLOC1()
      ALLOCATE( ILAD(NBIG), ILAM(NBIG), ILOD(NBIG), ILOM(NBIG))
      ALLOCATE( ISNAP1(NBIG), ISNAP2(NBIG), IROT(NBIG))
      ALLOCATE( ANGLE(NBIG), SPED(NBIG) )

      IZERO=0
      IB=0
      HH=500
      NN=0

c   READ INFO FROM TRACK FILE AND GET TIMEINC FROM FIRST TWO RECORDS
      OPEN(11+CK_LUN,FILE=INPUTFILE,STATUS='UNKNOWN')
      READ(11+CK_LUN,*) 
      READ(11+CK_LUN,'(A)') title
      READ(11+CK_LUN,*)
      NPTS = 0
      DO I=1,100
         READ(11+CK_LUN,1000,END=1010) INUM,XLAT(I),XLON(I),TSTAMP(I),
     &                          ISPD(I),IPRES(I)
         print *, INUM,XLAT(I),XLON(I),TSTAMP(I),ISPD(I),IPRES(I)
         NPTS=NPTS+1
      ENDDO
1000  FORMAT(I3,2X,F5.2,2X,F6.2,1X,A9,2X,I3,2X,I4)
1010  continue

      READ(TSTAMP(1)(7:8),'(I2)') IHOUR1
      READ(TSTAMP(2)(7:8),'(I2)') IHOUR2

      IF (IHOUR2 > IHOUR1) THEN
        TIMEINC  = IHOUR2-IHOUR1
      ELSE
        TIMEINC  = 24+IHOUR2-IHOUR1
      ENDIF
      SPEEDUP = 6.0/TIMEINC

      print *, "1st hour = ", ihour1
      print *, "2nd hour = ", ihour2
      print *, "timeinc = ", timeinc
      print *, "speedup = ", speedup

      IHR1 = 0.    !  P15 Always thinks storm starts a hour = 0
      DO I=1,NPTS
         TIME(I)= IHR1+(I-1)*6.0     
         print *, i, time(i)
      ENDDO


      write(*,*) 'Start of the computations...  '
      DO N=1,NPTS
        XLON(N)=-1*XLON(N)
        IF(XLAT(N).EQ.0) GO TO 1025
        NN=NN+1
        XLAT(NN)=XLAT(N)
        XLON(NN)=XLON(N)
        XPRES(NN)=IPRES(N)
        XVEL(NN)=ISPD(N)
      end do

c...

 1025 CONTINUE
 
      IF(XPRES(1).LT.100.) XPRES(1)=1000.
      IF(XPRES(NPTS).LT.100.) XPRES(NPTS)=1000.
      IST=1

 1031 DO N=IST,1000
        IF(N.GE.NPTS) GO TO 1032
        IF(XPRES(N).EQ.0) GO TO 1027
        end do
 1027 DO NN=N+1,1000
        IF(XPRES(NN).GT.0) GO TO 1029
        end do
 1029 DP=XPRES(NN)-XPRES(N-1)
      NOS=NN-N+1
      DO NNN=1,NOS
        IST=N-1+NNN
        XMUL=REAL(NNN)/REAL(NOS)
        XPRES(IST)=XPRES(N-1)+DP*XMUL
        end do
      IST=IST+1
      GO TO 1031
 
 1032 CONTINUE
      
      IEND=1
      CALL CUBSPL(TIME,XLAT,NPTS,IEND,A)
      ITOT=(NPTS-1)*6 
      DO I=0,ITOT
        XX=FLOAT(I)
        CALL VALSPL(A,TIME,XX,NPTS,XXLAT(I+1))
        XTIME(I+1)=XX
        end do
      IEND=1
      CALL CUBSPL(TIME,XLON,NPTS,IEND,A)
      DO I=0,ITOT
        XX=FLOAT(I)
        CALL VALSPL(A,TIME,XX,NPTS,XXLON(I+1))
        end do
      IEND=1
      CALL CUBSPL(TIME,XPRES,NPTS,IEND,A)
      DO I=0,ITOT
        XX=FLOAT(I)
        CALL VALSPL(A,TIME,XX,NPTS,XXPRES(I+1))
        end do 
      IEND=1
      CALL CUBSPL(TIME,XVEL,NPTS,IEND,A)
      DO I=0,ITOT
        XX=FLOAT(I)
        CALL VALSPL(A,TIME,XX,NPTS,XXVEL(I+1))
        end do

      DELSNP = 6.0     
      NN=0
      HH=500.

c...

      IIYR=IYR-1900
      IF(IIYR.LT.0) IIYR=IIYR+100
      ITIME=IIYR*1000000+IMO*10000+IDY*100+IHR1
      NZ=NPTS
      PCTINC=1./DELSNP
      PCTF=0.

c      RADIUS1=10000.
      PI=3.141592653589796
      RAD=3950.
      dx=5.

      XXLAT(ITOT+2)=XXLAT(ITOT+1)+(XXLAT(ITOT+1)-XXLAT(ITOT))
      XXLON(ITOT+2)=XXLON(ITOT+1)+(XXLON(ITOT+1)-XXLON(ITOT))

      NUM = 0
      DO N=1,ITOT+1
        XLAT12=PI/180.*XXLAT(N)
        DX12=2.*PI*RAD*COS(XLAT12)/360.*(XXLON(N)-XXLON(N+1))    ! DX to East
        DY12=PI*RAD/180.*(XXLAT(N+1)-XXLAT(N))                   ! DY to North
        IF(DY12.EQ.0.) THEN
          ANGLE(n)=90.
          GO TO 1047
        ELSE
          ANGLE(n)=180./PI*(ATAN(DX12/DY12))
          IF(DX12.GT.0..AND.DY12.LT.0.) ANGLE(n)=180.+ANGLE(n)
          IF(DX12.LT.0..AND.DY12.LT.0.) ANGLE(n)=180.+ANGLE(n)
          IF(DX12.LT.0..AND.DY12.GT.0.) ANGLE(n)=360.+ANGLE(n)
        ENDIF

cjjw - modified speed by speedup=6/timeinc to allow timeinc hourly input with faked times (still 6 hours)
cjjw   so that will put out timeinc/6 min snaps (although PBL will think of it as hourly )

        XDUM=SPEEDUP*SQRT(DX12*DX12+DY12*DY12)
 1047   CONTINUE
        SPED(n)=XDUM*5280./6080.           !eye speed in nautical miles/hour
        ILAD(N)=XXLAT(N)
        ILAM(N)=(XXLAT(N)-ILAD(N))*60.+.5
        ILOD(N)=XXLON(N)
        ILOM(N)=(XXLON(N)-ILOD(N))*60.+.5
        PFAR1=1013.
        DP=PFAR1-XXPRES(N)
        IF(DP.LT.20.) DP=20.

cjjw modified this section as per Norm's recommendation on 7-10-2002
cjjw        VZERO=1.34*DP+38.
cjjw        RZERO(N)=-(XXVEL(N)*6080./5280.-VZERO)/.520
C ADDED BY ED THOMPSON
        xspd = xxvel(n)*6080./5280.     !convert max wind to mph
        call JTrmax(xspd,dp,rzero(n),SPED(n))   !added 11/15/99 eft
cjjw-end

        IF(RZERO(N).LT.10.) RZERO(N)=10.
        IF(RZERO(N).GT.50.) RZERO(N)=50.
        RZERO(N)=RZERO(N)*5280./6080.          !convert radius to maximum into nautical miles
        IQUAD=0
        IROT(N)=0
        XN=N
        PCT1(N)=(XN-1.)/DELSNP+.000001
        IPCT=PCT1(N)
        PCT1(N)=PCT1(N)-IPCT
        XSNAP=1.+(XN-1.)/DELSNP+.000001
        ISNAP1(N)=XSNAP
        ISNAP2(N)=ISNAP1(N)+1
        IF(PCT1(N).LT..001) ISNAP2(N)=0
        IF(ISNAP2(N).EQ.nz) THEN
          PCTF=PCTF+PCTINC
          PCT1(N)=PCTF
        ENDIF
        IHOUR=N-1
        IF(ISNAP2(N).EQ.0..OR.N.EQ.(ITOT+1)) THEN
          NUM=NUM+1
          ZLA(NUM)=XXLAT(N)
          ZLO(NUM)=XXLON(N)
          ZPR(NUM)=XXPRES(N)
          ZAN(NUM)=ANGLE(N)
          ZSP(NUM)=SPED(N)
          IHOUR=(N/timeinc)*timeinc
          SGW=6.                                   ! METERS PER SECOND
          AN1(NUM)=360.-ZAN(NUM)+90.               ! MEASURED COUNTERCLOCKWISE FROM EAST
          IF(AN1(NUM).GT.360.) AN1(NUM)=AN1(NUM)-360.
          RRR=RZERO(N)
          RADIUS1(NUM)=RZERO(N)
        ENDIF

      END DO

      ICNVRT=0
      NPRT=0
      ISTART=ITIME

      CALL ADCIRC()
      WRITE(*,*) '  ***  AFTER ADCIRC'

      CALL SNAP()
      WRITE(*,*) '  ***  AFTER SNAP  '

      CALL PCHIST()
      WRITE(*,*) '  ***  AFTER PCHIST  '

      STOP
      END 


      SUBROUTINE ALLOC1()
      USE P15EXTENTS
      USE P1; USE C0; USE C1; USE C2; USE C3; USE C4; USE C5
      USE C6; USE C8; USE D13; USE DJM; USE IHOUR; USE P15GLOBAL
      ALLOCATE( PCT1(NBIG) )
      ALLOCATE( YLA(NBIG), YLO(NBIG), KSTA(NBIG), WINDRED(NBIG,12))
      ITRACK  => JA(1)
      EYELAT  => JA(2)
      EYELONG => JA(3)
      DIREC   => JA(4)
      SPEED   => JA(5)
      IQUAD   => JA(6)
      EYPRES  => JA(7)
      RADIUS  => JA(8:11)
      PFAR    => JA(12:15)
      ALLOCATE( XXLAT(NBIG), XXLON(NBIG), XTIME(NBIG) )
      ALLOCATE( XXVEL(NBIG), XXPRES(NBIG))
      ALLOCATE( ZLA(NBIG),ZLO(NBIG),ZPR(NBIG),ZAN(NBIG),ZSP(NBIG))
      ALLOCATE( IHR(NBIG) )
      RETURN
      END SUBROUTINE ALLOC1


      SUBROUTINE CUBSPL(X,Y,N,IEND,A)
C THIS ROUTINE CALCULATES THE COEFFICIENTS TO CUBIC SPLINE POLYNOMIALS
C FOR TWO CORRESPONDING SETS OF VALUES X AND Y.  THE COEFFICIENTS ARE
C STORED IN THE MATRIX A. IEND REFERS TO THE END CONDITIONS OF THE
C SPLINE INTERPOLATION.
      INTEGER NR,N,IEND,NM1,NM2,J
      PARAMETER(NR=500)
      REAL X(NR),Y(NR),S(NR),A(NR,4)
      REAL DX1,DY1,DX2,DY2,DXN2,DXN1
      NM2=N-2
      NM1=N-1
      DX1=X(2)-X(1)
      DY1=(Y(2)-Y(1))/DX1*6.
C DETERMINE COEFFICIENTS IN MATRIX FOR SOLUTION OF THE SYSTEM
C OF EQUATIONS
      DO 10 I=1,NM2
      DX2=X(I+2)-X(I+1)
      DY2=(Y(I+2)-Y(I+1))/DX2*6.
      A(I,1)=DX1
      A(I,2)=2.*(DX1+DX2)
      A(I,3)=DX2
      A(I,4)=DY2-DY1
      DX1=DX2
      DY1=DY2
   10 CONTINUE
C ADJUST MATRIX ACCORDING TO END CONDITIONS
      IF(IEND.EQ.2) THEN
      A(1,2)=A(1,2)+X(2)-X(1)
      A(NM2,2)=A(NM2,2)+X(N)-X(NM1)
      ELSEIF(IEND.EQ.3) THEN
      DX1=X(2)-X(1)
      DX2=X(3)-X(2)
      A(1,2)=(DX1+DX2)*(DX1+2.*DX2)/DX2
      A(1,3)=(DX2*DX2-DX1*DX1)/DX2
      DXN2=X(NM1)-X(NM2)
      DXN1=X(N)-X(NM1)
      A(NM2,1)=(DXN2*DXN2-DXN1*DXN1)/DXN2
      A(NM2,2)=(DXN1+DXN2)*(DXN1+2.*DXN2)/DXN2
      ENDIF
C SOLVE SYSTEM OF EQUATIONS
      DO 110 I=2,NM2
      A(I,2)=A(I,2)-A(I,1)/A(I-1,2)*A(I-1,3)
      A(I,4)=A(I,4)-A(I,1)/A(I-1,2)*A(I-1,4)
  110 CONTINUE   
      A(NM2,4)=A(NM2,4)/A(NM2,2)
      DO 120 I=2,NM2
      J=NM1-I
      A(J,4)=(A(J,4)-A(J,3)*A(J+1,4))/A(J,2)
  120 CONTINUE   
      DO 130 I=1,NM2
      S(I+1)=A(I,4)
  130 CONTINUE   
C DETERMINE THE SECOND DERIVATIVE AT THE END POINTS
      IF(IEND.EQ.1) THEN
      S(1)=0
      S(N)=0
      ELSEIF(IEND.EQ.2) THEN
      S(1)=S(2)
      S(N)=S(N-1)
      ELSE
      S(1)=((DX1+DX2)*S(2)+DX1*S(3))/DX2
      S(N)=((DXN2+DXN1)*S(NM1)-DXN1*S(NM2))/DXN2
      ENDIF
C DETERMINE THE COEFFICIENTS FOR THE SPLINE POLYNOMIALS
      DO 11 I=1,N-1
      A(I,1)=(S(I+1)-S(I))/6./(X(I+1)-X(I))
      A(I,2)=S(I)/2
      A(I,3)=(Y(I+1)-Y(I))/(X(I+1)-X(I))-(X(I+1)-X(I))*(2*S(I)+
     1S(I+1))/6
      A(I,4)=Y(I)
   11 CONTINUE
      RETURN
      END



C THIS ROUTINE CALCULATES THE DEPTH AT A CERTAIN POINT (XX) ALONG
C THE PROFILE FROM THE CUBIC SPLINE POLYNOMIALS
      SUBROUTINE VALSPL(A,X,XX,NRV,HVAL)
      INTEGER I,NRV
      REAL A(500,4),X(500),XX,HVAL
      DO 10 I=1,NRV-1
      IF(X(I).LE.XX.AND.X(I+1).GE.XX) THEN
      HVAL=A(I,1)*(XX-X(I))**3+A(I,2)*(XX-X(I))**2+
     1A(I,3)*(XX-X(I))+A(I,4)
      ENDIF
   10 CONTINUE
      RETURN
      END


      SUBROUTINE ADCIRC ( )
      USE P15EXTENTS; USE P15GLOBAL
C THIS SUBROUTINE READS ALL NODES IN ADCIRC GRID
C      CHARACTER*20 TITLE
C--------------------------------------------------------------------------

      DO 10 N=1,NBR
      READ(10+CK_LUN,*) NODE,YLO(N),YLA(N),DEPP
      KSTA(N)=N
      YLO(N)=ABS(YLO(N))
   10 CONTINUE
      RETURN
      END


      SUBROUTINE SNAP
      USE C0; USE C1; USE C2; USE C3; USE C4; USE C5
      USE C6; USE DJM; USE IHOUR; USE P15GLOBAL
C     MAIN PROGRAM FOR SNAPSHOT WINDS ON NESTED GRID.                           
C     INPUT  IB IS SWITCH VARIABLE  IB = 0 TO SUPPRESS PRINTING                 
C                                          OF PRESSURE FIELD
C                                          AND INITIAL WIND FIELD               
C     NZ IS NUMBER OF WIND SNAPSHOTS TO PRODUCE.                         
C     FOR THE QUANTITIES EQUIVALENCED TO JA, SEE COMMENTS TO SUBROUTINE GRAD 
C     NM IS NUMBER OF TIMES TO CYCLE WIND COMPUTATION IN INNERMOST GRID NEST. 
C     SGW IS MAGNITUDE OF SURFACE GEOSTROPHIC WIND (MTRS/SEC)            
C     AN1 IS DIRECTION OF SGW, COUNTERCLOCKWISE FROM EAST (DEG)          
C     ST12 IS DIST (KM) FROM AXIS TO 1/2 MAGNITUDE OF SGW_) FROM A        
C     DX IS GRID DISTANCE OF INNERMOST NEST (KILOMETERS).                
C     FOR ALL ARRAYS IN C3 
C      1ST DIMENSION INCREASES FROM WEST TO EAST                          
C      2ND DIMENSION INCREASES FROM SOUTH TO NORTH                        
C      3RD DIMENSION(NEST NUMBER) INCREASES FROM                          
C      INNERMOST TO OUTERMOST                                   
C--------------------------------------------------------------------------
C*****   DIMENSIONED FOR 100 SNAPSHOTS
C------------------------------------------------------------------------------------
      PTH = 300.0          ! debug
      Z0LAND = 0.08        ! debug                                                  
      K35 = 0.35           ! debug                                                  
      G = 9.805            ! debug                                                 
      GARR = 0.0144        ! debug
      DTH(1) = 0.0         ! debug                                                          
      DTH(2) =-2.0         ! debug                                                          
C                                                                               
      K2=K35**2                                                 
      GA = GARR/G                                                      
C        
      DO 20 NSNAP= 1,NZ 
        write(*,*) ' ' 
        WRITE(*,*) '  TOP OF LOOP 20  ', NSNAP
        NM = 800                                                        
        ST12 = 0.                                                      
        ITRACK = 0                                 
        EYELAT=ZLA(NSNAP)
        EYLONG=ZLO(NSNAP)
        EYPRES=ZPR(NSNAP)
        DIREC=ZAN(NSNAP)
        SPEED=ZSP(NSNAP)
        RADIUS(1)=RADIUS1(NSNAP)
        PFAR(1)=PFAR1
        ihr(nsnap)=ihr1+((nsnap-1)*timeinc)

        write(*,*) DX,IQUAD,IHR(NSNAP)
        write(*,*) EYELAT,EYLONG,EYPRES
        write(*,*) DIREC,SPEED
        write(*,*) RADIUS
        write(*,*) PFAR
        write(*,*) SGW,AN1(NSNAP)
        write(*,*) NM,ST12,ITRACK
        write(*,*) DTH,HH,GARR,PTH,K35
        write(*,*) JA   

        IF(PFAR(1).LT.EYPRES) PFAR(1)=EYPRES
        P0SNP(NSNAP)=EYPRES
        PFRSNP(NSNAP)=PFAR(1)
        RADSNP(NSNAP)=RADIUS(1)
        DT = 10.0*DX   
        PHI = EYELAT/57.29578            !PHI IS LATITUDE IN RADIANS
        F=2.*7.29E-5*SIN(PHI)            !F IS CORIOLIS FORCE 
        FLAT=F      
        CALL CCROSS 
        write(*,*) ' '    
        write(*,*) '  AFTER CCROSS '  
        CALL BLOWUQ 
        write(*,*) ' '  
        write(*,*) '  AFTER BLOWUQ ' 

        DO K = 1,  5
        DO J = 1, 21
        DO I = 1, 21
           SMF_UN(I,J,K,NSNAP) = UN(I,J,K)
           SMF_VN(I,J,K,NSNAP) = VN(I,J,K)
        END DO
        END DO
        END DO
        DO K = 1, 15
           ISMF_JA(K,NSNAP) = JA(K)
        END DO
        SMF_DXX(NSNAP) = DX
        SMF_SGW(NSNAP) = SGW
        SMF_AN1(NSNAP) = AN1(NSNAP)
        SMF_ST1(NSNAP) = ST12
        SMF_HHH(NSNAP) = HH
        SMF_GAR(NSNAP) = GARR
        SMF_PTH(NSNAP) = PTH
        SMF_K35(NSNAP) = K35
        SMF_DTH1(NSNAP) = DTH(1)
        SMF_DTH2(NSNAP) = DTH(2)
        print *, "from snap:  ntimes = ", ntimes
        NTIMES=NZ
C------------------------------------------------------------------------------------
        write(*,*) ' '    
        write(*,*) '  BOTTOM OF LOOP 20  ', NSNAP
        write(*,*) DX,IQUAD,IHR(NSNAP)
        write(*,*) EYELAT,EYLONG,EYPRES
        write(*,*) DIREC,SPEED
        write(*,*) RADIUS
        write(*,*) PFAR
        write(*,*) SGW,AN1(NSNAP)
        write(*,*) NM,ST12,ITRACK
        write(*,*) DTH,HH,GARR,PTH,K35   
        write(*,*) JA   
 20   CONTINUE                                                          
      WRITE (6+CK_LUN,25)                                                     
   25 FORMAT('1 END OF SNAP/MAIN') 
      RETURN                                                            
      END 



      SUBROUTINE PCHIST()
      USE C0; USE C2; USE C8; USE DJM 
      USE C57; USE D13; USE P15GLOBAL
C-------------------------------------------------------------------------------
C     **** THIS IS SUBROUTINE HIST ****                                       
C    
C**   A.C.E. PROGRAM TO FIND WINDS AT STATIONS                                  
C                                                                               
C**   INPUT CARDS#                                                              
C       1. NAMELISTS NAME4,NAME5
C       2. HINDCAST LOCATIONS REQUESTED# LAT-DEG,LAT-MIN,LON-DEG,          
C          LON-MIN, TERRAIN CODE, STA HT IN METERS,                
C          STATION IDENTIFICATION NUMBER.                          
C          (5I4,F6.1,I4)                                           
C          WEST LONGITUDE IS POSITIVE.                             
C          1 CARD PER LOCATION. TERMINATED BY ?EOF.                
C       3. ONE CARD FOR EACH HOUR OF STORM.                                
C          LAT-DEG OF EYE,LAT-MIN OF EYE,LON-DEG OF                
C          EYE,LON-MIN OF EYE,SNAPSHOT 1 RECORD SEQUENCE           
C          NUMBER IN FILE,SNAP 2 REC SEQ NBR,                      
C          ROTATION ANGLE(DEG.CLOCKWISE                            
C          TO ROTATE WIND ON NESTED GRID),                         
C          INDICATOR - IF NONZERO WIND ON WAVE GRID IS PRIN        
C          (6I4,F8.4,2I4) .                                        
C          WEST LONGITUDE IS POSITIVE .  TERMINATED BY ?EOF        
C**   OUTPUT FILES#                                                             
C       1. FILE 20(LTROUT)# WIND DATA ON ICOSAHEDRAL GRID(R.A.F.)          
C                                                                               
C**   TEMPORARY FILES                                                           
C       1. FILE 10(LTEMP)  ALL SNAPSHOTS PLUS ALL INTERPOLATED WIND        
C                                                                               
C**   STANDARD UNITS# 5 IS CARD READER, 6(LP) IS PRINTER.                   
C                                                                               
C     TERRAIN CODE (LLAKE)  1-6,        LAKE  0-5.                              
C                                                                               
      PARAMETER (NHOURS=600) 
      PARAMETER (MAXI=11,MAXJ=11,IGRDHT=100) 
      CHARACTER*4 IZONE*3
      INTEGER LSTAB(MAXI,MAXJ)                                   
      INTEGER NAD(NHOURS), NAM(NHOURS), NOD(NHOURS), NOM(NHOURS) 
      INTEGER JSEQ(NHOURS), KDATE(NHOURS), KTIME(NHOURS)
      REAL DTHI(2),ZANG(MAXI,MAXJ)                             
      REAL ::  XY(21,21,10)
      REAL WIND(2,MAXI,MAXJ)                                       
      EQUIVALENCE (XY,WIND)                                    
C-------------------------------------------------------------------------------
      CON = 57.29578                                              
      LP  = 6                                                  
      KSNAP1 = 0    
      KSNAP2 = 0    
      LTEMP = 43                                                    
      OPEN(UNIT=LTEMP+CK_LUN, FILE='LTEMP.OUT',STATUS='UNKNOWN')
      OPEN(UNIT=69+CK_LUN, FILE='fort.22',STATUS='UNKNOWN')
C-----------------------------------------------------------------------------------
      IZONE='GMT'
      GRIDHT=FLOAT(IGRDHT)/10.                                          
      IMAX=MAXI                                                         
      JMAX=MAXJ
    4 DO 5 J=1,MAXJ                                                    
      DO 5 I=1,MAXI                                                     
      LSTAB(I,J)=0                                                      
    5 ZANG(I,J)=0.                                                      
      DTH=-2.                                                           
      LAKE=0                                                           
      G=9.806                                                           
      GARR=.035
      NTIM=(NTIMES-1)*6+1
C------------------------------------------------------------------------------------ 
      PTH= NHOURS
C------------------------------------------------------------------------------------
      K35=.35    
      DO J=1,5  
      DO I=1,3 
         ZCOEFF(I,J)=0. 
      END DO
      END DO
      REWIND LTEMP 
C------------------------------------------------------------------------------------
      N_DJM = 1
      DO K = 1,  5
      DO J = 1, 21
      DO I = 1, 21
        XX(I,J,K) = SMF_UN(I,J,K,N_DJM)
      END DO
      END DO
      END DO
      DO K = 6, 10
      DO J = 1, 21
      DO I = 1, 21
         L = K - 5
         XX(I,J,K) = SMF_VN(I,J,L,N_DJM)
      END DO
      END DO
      END DO
      DO K = 1, 15
         JA(K) = ISMF_JA(K,N_DJM)
      END DO
      DX   = SMF_DXX(N_DJM)
      SGW  = SMF_SGW(N_DJM)
      AN1  = SMF_AN1(N_DJM)
      ST12 = SMF_ST1(N_DJM)
      HH   = SMF_HHH(N_DJM)
      GARR = SMF_GAR(N_DJM)
      PTH  = SMF_PTH(N_DJM)
      K35  = SMF_K35(N_DJM)
      DTHI(1) = SMF_DTH1(N_DJM)
      DTHI(2) = SMF_DTH2(N_DJM)
C------------------------------------------------------------------------------------
 1233 FORMAT(A4)
      DTH=DTHI(2)
      K2=K35**2  
      GA=GARR/G 
      DT=10.0*DX
      CALL UXXV     
      CALL UPDOWN
C                                                                               
C     READ HOURLY INPUT CARDS .                                                 
      NHT=0    
   25 NHTP=NHT+1
      NAD(NHTP)=ILAD(NHTP)
      NAM(NHTP)=ILAM(NHTP)
      NOD(NHTP)=ILOD(NHTP)
      NOM(NHTP)=ILOM(NHTP)
      NSNAP1(NHTP)=ISNAP1(NHTP)
      NSNAP2(NHTP)=ISNAP2(NHTP)
      PCT(NHTP)=PCT1(NHTP)
      NHT=NHT+1
C------------------------------------------------------------------------------------
      IF (NHT.LT. NHOURS) GO TO 25
C------------------------------------------------------------------------------------
   30 CONTINUE 
      KY=ISTART/10**6+1900 
      KM=MOD((ISTART/10**4),100)
      KD=MOD((ISTART/100),100)
      KTIME(1)=MOD(ISTART,100)
      KJD=JULIAN(KM,KD,KY)
      KDATE(1)=ISTART/100                                           
      JSEQ(1)=1
      DO 35 J=2,NTIM 
      JSEQ(J)=JSEQ(J-1)
      IF(NSNAP1(J-1).NE.NSNAP1(J).OR.NSNAP2(J-1).NE.NSNAP2(J).OR.
     & PCT(J-1).NE.PCT(J))JSEQ(J)=JSEQ(J)+1
      KDATE(J)=KDATE(J-1)
C  *** REDUCED WINDS ARE COMPUTED AT 1 HR TIME STEPS ***                        
      KTIME(J)=KTIME(J-1)+1
      IF (KTIME(J).LT.24) GO TO 35
      KTIME(J)=0
      KJD=KJD+1
      CALL INVJD (KJD,KM,KD,KY)
      KDATE(J)=(KY-1900)*10**4+KM*100+KD
   35 CONTINUE 
      
      DO 80 KHR=1,NTIM      
C     
      IF (KSNAP1.EQ.NSNAP1(KHR)) GO TO 50 
      KSNAP1=NSNAP1(KHR)
                                                        
      DO K = 1,  5
      DO J = 1, 21
      DO I = 1, 21
        XX(I,J,K) = SMF_UN(I,J,K,KSNAP1)
      END DO
      END DO
      END DO

      DO K = 6, 10
      DO J = 1, 21
      DO I = 1, 21
        L = K - 5
        XX(I,J,K) = SMF_VN(I,J,L,KSNAP1)
      END DO
      END DO
      END DO
C-------------------------------------------------------------------------------------
   50 IF (NSNAP2(KHR).EQ.0) GO TO 60                                   
      IF (NSNAP2(KHR).EQ.KSNAP2) GO TO 60                          
      KSNAP2=NSNAP2(KHR)
      DO K = 1,  5
      DO J = 1, 21
      DO I = 1, 21
         XY(I,J,K) = SMF_UN(I,J,K,KSNAP2)
      END DO
      END DO
      END DO
      DO K = 6, 10
      DO J = 1, 21
      DO I = 1, 21
         L = K - 5
         XY(I,J,K) = SMF_VN(I,J,L,KSNAP2)
      END DO
      END DO
      END DO
C-------------------------------------------------------------------------------------
   60 CONTINUE
      IF (PCT(KHR).EQ.0.) GO TO 70                                     
      KSNAP1=0                                                         
      IF (PCT(KHR).LT.0..OR.PCT(KHR).GT.1.) STOP 444                   
      DO 65 K=1,10                                                      
      DO 65 J=1,21                                                      
      DO 65 I=1,21                                                      
      XX(I,J,K)=(1.-PCT(KHR))*XX(I,J,K)+PCT(KHR)*XY(I,J,K)              
   65 CONTINUE                                                          
   70 CONTINUE                                                         
C                                                                               
      IF (KHR.EQ.1) GO TO 75                                           
      IF (JSEQ(KHR).EQ.JSEQ(KHR-1)) GO TO 80
   75 WRITE(LTEMP+CK_LUN,*) XX
      WRITE(LTEMP+CK_LUN,*) JSEQ(KHR)
   80 CONTINUE
      IF (NBR.EQ.0) GO TO 106                                       
C                                                              
C     COMPUTE WIND FOR EACH ADCIRC GRID LOCATION AT THIS TIME STEP 
C                                                                               
      STHT=GRIDHT
      LANSEA=1
C--------------------------------------------------------------------------
      REWIND LTEMP 
C--------------------------------------------------------------------------
      DO 100 KHR = 1, NTIM
         ramp=tanh((khr-1)/18.d0)
      WRITE(69+CK_LUN,2222) KHR,NAD(KHR),NAM(KHR),NOD(KHR),NOM(KHR),
     1               PFRSNP(NSNAP1(KHR)),SGW
      WRITE(6+CK_LUN,2222) KHR,NAD(KHR),NAM(KHR),NOD(KHR),NOM(KHR),
     1              PFRSNP(NSNAP1(KHR)),SGW
 2222 FORMAT(' #',I10,' EYE AT (LAT,LON)',2I3,2X,2I3,
     $       ' PFAR =',F8.2,' SGW =',F8.2)
C-------------------------------------------------------------------------- 
      IF (KHR.EQ.1) GO TO 85                                           
      IF ( JSEQ(KHR) .EQ. JSEQ(KHR-1) ) GO TO 90
   85 READ (LTEMP+CK_LUN,*) XX
      READ(LTEMP+CK_LUN,*) KSEQ
 90   CONTINUE
C--------------------------------------------------------------------------
      iout=0
      DO 105 K = 1, NBR
      LA1=YLA(K)*3.14159/180.                                         
      LO1=YLO(K)*3.14159/180.
      XL=IABS(NAD(KHR))                                                 
      LA0=(XL+FLOAT(NAM(KHR))/60.0)/CON                                
      IF(NAD(KHR).LT.0)LA0=-LA0                                        
      XL=IABS(NOD(KHR))                                                
      LO0=(XL+FLOAT(NOM(KHR))/60.0)/CON  
      IF(NOD(KHR).LT.0)LO0=-LO0                                        
      ROT=IROT(KHR)
      CALL BREEZE
      IF( W1 .LT. 0.001 ) GO TO 105
      IXX=NSNAP1(KHR)
      IF(PCT(KHR).EQ.0.) THEN
            P0HR=P0SNP(IXX)
            PFHR=PFRSNP(IXX)
            RDHR=RADSNP(IXX)
         ELSE
            P0HR=(1.-PCT(KHR))*P0SNP(IXX)+PCT(KHR)*P0SNP(IXX+1)
            PFHR=(1.-PCT(KHR))*PFRSNP(IXX)+PCT(KHR)*PFRSNP(IXX+1)
            RDHR=(1.-PCT(KHR))*RADSNP(IXX)+PCT(KHR)*RADSNP(IXX+1)
      ENDIF
      IF(D.LE.0.1) D=0.1
      PRESS=P0HR+(PFHR-P0HR)*EXP(-RDHR*1.85/D)
      UST=UST*100.                                                      
C     UST IN CM/SEC
      MAD=YLA(K)
      MAM=(YLA(K)-MAD)*60.
      MMOD=YLO(K)
      MOM=(YLO(K)-MMOD)*60.
      TH1=TH1+180.
      U99=W1*SIN(3.14159/180.*TH1)
      V99=W1*COS(3.14159/180.*TH1)
      UST=UST*UST
      UTAU=UST*SIN(3.14159/180.*TH1)
      VTAU=UST*COS(3.14159/180.*TH1)
      IF(TH1.GT.360.) TH1=TH1-360.
      IF(AL.GT.360.) AL=AL-360.
C  CONVERT D IN KM TO NAUTICAL MILES
      D=D*0.5395
      iout=iout+1
c...convert marine wind speed from knots to m/s
      um=u99*1.04d0*0.5144d0
      vm=v99*1.04d0*0.5144d0
      uvm=(um**2.d0+vm**2.d0)**0.5d0
      fr=1.0000
Csb 08/15/2006
C-- Convert 30-min winds to 10-min winds
      U99 = U99 * 1.04d0
      V99 = V99 * 1.04d0
C--
      ramppress=1013.-ramp*(1013.-press)
      WRITE(69+CK_LUN,2223) KSTA(K), ramp*fr*U99, ramp*fr*V99, ramppress
 2223 FORMAT( I8, 5E13.5)
  105 CONTINUE
      if(iout.eq.0) WRITE(69+CK_LUN,2223) 1, 0., 0., 1013.
  100 CONTINUE                                             
  106 CONTINUE
  135 WRITE (LP+CK_LUN,205)                                                   
cvjp  CALL SYSTEM('rm LTEMP.OUT')
      STOP 999                   
  160 FORMAT (4X,I6,4X,I8,I3)                                       
  165 FORMAT (5I4,F6.1,1X,I3)                                          
  170 FORMAT(1X,2(I4,1X,I2),I3,F6.1,I3)                              
  175 FORMAT (//)                                                      
  180 FORMAT(6I4,F8.4,2I4)                                              
  185 FORMAT (1H0)                                                      
  190 FORMAT (1H1,///,T20,'STORM HISTORY   1ST HOUR IS  ',I8,1X,A3,//,  
     1  (1X,I4,6I4,F8.4,3I4,I8,I2))         
  195 FORMAT (1X,A4,I8,I2,'  W1=',F6.2,' TH1=',F5.1,'D =',F6.0,' AL=', 
     1  F5.1,' UST=',F6.2,'  STA',I4,'=',I3,1X,I2,I4,1X,I2,      
     2  '  EYE=',I3,1X,I2,I4,1X,I2,'  TERR=',I1,           
     3  ' HT=',F5.1,I4,I4)                                            
  200 FORMAT (/)                                                       
  205 FORMAT('0  END OF HIST/MAIN')                                     
      END                                                              



      SUBROUTINE ABCCC                           
      USE C57   
      Z0 = GA*UX(K123)**2/HH                                  
      ZLOG = ALOG(Z0)                                                
   53 HL = -HH*DEN/(UX(K123)**2*PTH*(ZLOG+CM))                         
      IF (HL+2.) 54,55,56                                           
   54 HLOG = ALOG(-HL)                                                
      AM = AMIN1(HLOG+1.5,-.875*ZLOG)
C                                                                               
C     CHECK FOR NEG. EXP. LIMIT ON HIS-6000                                     
C                                                                               
      IF (HL.LT.-400.) BM=0.                                           
      IF (HL.GE.-400.) BM=1.8*EXP(.2*HL)                               
      CM = AMIN1(HLOG+3.7,-.875*ZLOG)                               
      GO TO 60                                                         
   55 AM = 2.1931472                                                    
      BM = 1.2065761                                                   
      CM = 4.3931472                                                
      GO TO 60                                                        
   56 IF (HL-2.) 57,58,59                                              
   57 AM = 1.3865736-.4032868*HL                                     
      BM = 1.953288+.37335598*HL                                    
      CM = 2.5465736-HL*.9232868                                 
      GO TO 60                                                 
   58 AM = .58                                                     
      BM = 2.70                                                  
      CM = .7                                                       
      GO TO 60                                                    
   59 YLOG = ZLOG+4.7                                                  
      HL = .25*(YLOG+SQRT(YLOG**2+8.*HH*DEN/(UX(K123)**2*PTH)))  
      AM = AMAX1(2.5-.96*HL,-99.)                             
      BM = AMIN1(1.1+.8*HL,99.)                        
      CM = 4.7-2.*HL                                           
   60 UV(K123) = UX(K123)**2/K2*((ZLOG+AM)**2+BM**2)  
      DUV(K123) = UV(K123)-VV2         
      RETURN                                                   
      END      



      SUBROUTINE BREEZE                                                 
      USE C57; USE D13
      REAL LB,LL                                        
      REAL,TARGET ::  XY(2)
      REAL,POINTER :: U1,V1,UXX
      HAV(X) = (SIN(.5*X))**2                                          
      AHAV(X) = 2.*ASIN(SQRT(X))                             
      SQRU(X) = SQRT(ABS(X))     
C------------------------------------------------------------------------------
C     LA0,LO0 ARE LAT,LON OF EYE(PT0)            INPUT (RADIANS)                
C     LA1,LO1 ARE LAT,LON OF POINT AT WHICH WIND IS WANTED(PT1)-INPUT(RA        
C              LONGITUDE IS POSITIVE WEST.                                      
C     DX IS GRID SPACING OF INNERMOST NEST - INPUT(KILOMETERS)                  
C     W1 IS WIND SPEED AT PT1,20 METERS - OUTPUT(KNOTS)                         
C     D IS DISTANCE BETWEEN PT0 AND PT1  -  OUTPUT(KM)                          
C     ROT IS ANGLE FROM TRUE NORTH TO Y-AXIS OF NESTED RECTANGULAR              
C        GRID WIND FIELD (CLOCKWISE,DEGREES).                                   
C     TH1 IS DIRECTION TO WHICH WIND BLOWS, CLOCKWISE FROM SOUTH(DEG)           
C     AL IS BEARING OF POINT 1 FROM POINT 0, CLOCKWISE FROM NORTH(DEG)          
C     LANSEA IS TERRAIN CODE AT POINT 1                                         
C       1 IS OCEAN, 2 TO 6 ARE VARIOUS GROUNDS AND LAKES                        
C------------------------------------------------------------------------------
      U1 => XY(1)
      V1 => XY(2)
      UXX => UST
C------------------------------------------------------------------------------
    5 LL = LA0+LA1                                                    
      LB = LA0-LA1                                                     
      R = AHAV(HAV(LB)+COS(LA0)*COS(LA1)*HAV(LO0-LO1))
      IF(R.EQ.0.) GO TO 16                                    

C------------------------------------------------------------------------------
CCCC      IF (R.LT.1E-4) AB = ATAN2((LO0-LO1)*COS(.5+LL),-LB)   
      IF(R .LT. 1E-4) THEN
        DUMMY1 = (LO0-LO1) * COS(.5+LL)
        DUMMY2 = -LB
        IF( ABS( DUMMY2 ) .LE. 1.E-8 ) THEN
          AB = 1.570796                 ! PI / 2
        ELSE
          AB = ATAN2( DUMMY1, DUMMY2)
        END IF
      END IF
C------------------------------------------------------------------------------                    
      AB = AB*.5
C------------------------------------------------------------------------------                                                                
CCCC      IF(R.GE.1E-4) AB = ATAN2(SQRU(COS(.5*(LL+R))*SIN(.5*(R+LB))),             
CCCC     1  SQRU(COS(.5*(LL-R))*SIN(.5*(R-LB)))) 
      IF(R .GE. 1E-4) THEN
        DUMMY1 = SQRU( COS(.5*(LL+R)) * SIN(.5*(R+LB)) )
        DUMMY2 = SQRU( COS(.5*(LL-R)) * SIN(.5*(R-LB)) )
        IF( ABS( DUMMY2 ) .LE. 1.E-8 ) THEN
          AB = 1.570796                 ! PI / 2
        ELSE
          AB = ATAN2( DUMMY1, DUMMY2)
        END IF
      END IF                       
C------------------------------------------------------------------------------
      IF (LO1 .GT. LO0) AB = -AB                              
      AL = 2.*57.29578*AB                                          
      IF (AL .LT. 0.) AL = AL+360.                                 
      DE = (AL-ROT)/57.29578
      D = R*1.852*3437.7468
      AA = 10.0*DX 
      IZ = 1                                                          
      XY(1) = D*SIN(DE)                                              
      XY(2) = D*COS(DE)                                                
      AB = AMAX1(ABS(XY(1)),ABS(XY(2))) 
   10 IF (AB .LT. AA) GO TO 11                                 
      AA = AA+AA                                                  
      IZ = IZ+1
      GO TO 10  
   11 IF (IZ .LE. 5) GO TO 12
C------------------------------------------------------------------------------
      U1 = 0.   
      V1 = 0.   
      W1 = 0.                                              
  110 TH1 = 0.
C------------------------------------------------------------------------------
      UXX=0.
      RETURN
 12   CONTINUE
C-----------------------------------------------------------------------------
CCCC      WRITE(6,*) '    AA ', AA
C-----------------------------------------------------------------------------
      DO 13 IA = 1,2           
   13 XY(IA) = (AA+XY(IA))*10./AA
      I = XY(1)
      J = XY(2)                                                   
      F1 = XY(1)-I                                                 
      F2 = XY(2)-J  
      G11 = (1.-F1)*(1.-F2)                             
      G12 = (1.-F1)*F2                                            
      G21 = F1*(1.-F2)                                            
      G22 = F1*F2      
      DO 14 IA = 1,2                                                  
      IB = 5*IA+IZ-5                         
   14 XY(IA) = G11*XX(I+1,J+1,IB)+G12*XX(I+1,J+2,IB)+G21*   
     *XX(I+2,J+1,IB)+G22*XX(I+2,J+2,IB)           
   15 W1 = U1*U1+V1*V1 
      IF (W1.EQ.0.) GO TO 110                                     
      W1 = SQRT(W1)                                             
C     REDUCE W1 TO STATION HEIGHT                                               
 900  SPP = 1.25*W1                
      KPP = SPP             
      SPP = SPP-KPP                     
      IF (KPP .NE. 0) GO TO 910      
      UXX =W1*UXV(1,LANSEA)  
      TWIST = TURN(1,LANSEA)    
      GO TO 930                  
 910  IF (KPP .GE. 100) GO TO 920
      UXX = W1*((1.-SPP)*UXV(KPP,LANSEA)+SPP*UXV(KPP+1,LANSEA))
      TWIST = (1.-SPP)*TURN(KPP,LANSEA)+SPP*TURN(KPP+1,LANSEA)          
      GO TO 930                                          
 920  UXX = W1*UXV(100,LANSEA)                                 
      TWIST = TURN(100,LANSEA)                                     
  930 IF (LANSEA .EQ. 1) Z0 = GA*UST**2                              
      IF (LANSEA .NE. 1) Z0 = ZCOEFF(1,LANSEA-1)/UST+  
     $  ZCOEFF(2,LANSEA-1)*UST**2+ZCOEFF(3,LANSEA-1)
C-----------------------------------------------------------------------------
CCCC      WRITE(6,*) '    Z0 ', Z0
C-----------------------------------------------------------------------------
      Z20 = STHT/Z0                                                     
      W1 = (3600./1852.)*AMIN1                                          
     $  (UXX/K35*ALOG(Z20),W1)
C-----------------------------------------------------------------------------
      IF( ABS( V1 ) .LT. 1.E-8 ) THEN
         DUMMY1 = 1.570796
      ELSE
         DUMMY1 = ATAN2(-U1,-V1)
      END IF
      TH1 = (DUMMY1 + TWIST) * 57.29578 + ROT
C-----------------------------------------------------------------------------
      IF (TH1 .LT. 0.) TH1 = TH1+360.              
      RETURN                                                    
   16 U1=XX(11,11,1)                                         
      V1=XX(11,11,6)  
      D=0.                            
      AL=0.                          
      GO TO 15                 
      END                                                  



      SUBROUTINE INVJD (J,M,D,Y)
C     REVERSE OF FUNCTION 'JULIAN'.                                
C     COMPUTES INVERSE JULIAN DATE J FROM                                  
C     MONTH(M),DAY(D),AND YEAR(Y) INPUT.                                        
      INTEGER J,M,D,Y,TJ,TM,TD,TY,MTD   
      TJ=J-1721119                                                  
      TY=(4*TJ-1)/146097                                      
      TJ=4*TJ-1-146097*TY                                    
      TD=TJ/4                                                     
      MTD=4*TD+3                                                 
      TJ=MTD/1461                                               
      TD=MTD-1461*TJ                                             
      TD=(TD+4)/4                                                
      MTD=5*TD-3                                               
      TM=MTD/153                                                  
      TD=MTD-153*TM                                              
      D=(TD+5)/5                                                
      Y=100*TY+TJ                                                
      IF (TM.GE.10) GO TO 2                                      
    1 M=TM+3                                                      
      RETURN                                               
    2 M=TM-9                                                     
      Y=Y+1                                                           
      RETURN                                                          
      END 


      FUNCTION JULIAN(MO,DA,YR)
C     REVERSE OF SUBROUTINE 'INVJD'.                                            
C     COMPUTES JULIAN DATE.                                                     
      INTEGER D,Y,M,C,YA,MO,DA,YR      
      M=MO                  
      D=DA           
      Y=YR         
      IF (M.LE.2) GO TO 2     
    1 M=M-3          
      GO TO 3     
    2 M=M+9             
      Y=Y-1      
    3 C=Y/100                  
      YA = Y-100*C         
      JULIAN=(146097*C)/4+(1461*YA)/4+(153*M+2)/5+D+1721119     
      RETURN                                                         
      END 



      SUBROUTINE LPRINT(KHR,ISTART,IZONE,WIND,LSTAB) 
      PARAMETER (MAXI=11,MAXJ=11,NLN=2,NLT=2)
      COMMON /LGRID/ ILAT(MAXJ),ILONG(MAXI)                        
!$omp threadprivate(/LGRID/)
      DIMENSION WIND(2,MAXI,MAXJ),LSTAB(MAXI,MAXJ)        
      DIMENSION ILN(NLN),LLN(NLN),ILT(NLT),LLT(NLT),KODES(6)  
      DIMENSION LIST(24)                               
      DATA ILN/1,21/, LLN/20,36/                                  
      DATA ILT/26,12/, LLT/13,1/                                
C                                                                               
C     PAGE FORMAT IS 20 LON. BY 14 LAT.                                         
C                                                                               
      DATA KODES/3H   ,3H***,3H===,3H---,3H+++,3H$$$/       
      KODESP(N) = KODES(N)                          
      KT=(KHR-1)*3                                       
      DO 60 NPLAT=1,NLT                                        
      LLATSQ=LLT(NPLAT)                     
      ILATSQ=ILT(NPLAT)                                 
      DO 55 NPLONG=1,NLN                                      
      J1=ILN(NPLONG)                                          
      J2=LLN(NPLONG)                                            
      PRINT 20, KT,ISTART,IZONE,(ILONG(J),J=J1,J2)       
   20 FORMAT('1 STORM ',A6,'  GULF OF MEXICO WINDS AT HOUR ',I3, 
     1  8X,'1ST HOUR IS ',I10,1X,A3//5X,20I6)          
      DO 30 L=LLATSQ,ILATSQ                                    
      LAT=ILATSQ+LLATSQ-L                       
      KOUNT=0                                                 
      DO 22 I=J1,J2                                           
      KOUNT=KOUNT+1                                                 
      LIST(KOUNT)=KODESP(LSTAB(I,LAT))     
   22 CONTINUE                                          
      PRINT 23, ILAT(LAT),(WIND(1,J,LAT),J=J1,J2)             
      PRINT 24, ILAT(LAT),(WIND(2,J,LAT),J=J1,J2)               
      PRINT 25, (LIST(K),K=1,KOUNT)                               
   23 FORMAT('0',I4,20F6.1)                                    
   24 FORMAT(1X,I4,1X,20(1X,F4.0,1X))                            
   25 FORMAT(8X,20(A3,3X))                                        
   30 CONTINUE          
   55 CONTINUE  
   60 CONTINUE    
      RETURN    
      END 



      SUBROUTINE RDGRID(ZLA,ZLO,LSTAB)   
      PARAMETER (MAXI=11,MAXJ=11)                                
      COMMON /LGRID/ ILAT(MAXJ),ILONG(MAXI)                
!$omp threadprivate(/LGRID/)
      DIMENSION LSTAB(MAXI,MAXJ)                 
      DIMENSION ZLA(MAXI,MAXJ),ZLO(MAXI,MAXJ) 
C     
C     55 KM ATLANTIC HURRICANE GRID
C
      i=0
      do 10 ii=-14000,-15000,-100
        i=i+1
        j=0
        do 10 jj=1000,2000,100
        j=j+1
      ZLO(I,J)=ii*1.74532E-4
      ZLA(I,J)=jj*1.74532E-4
      LSTAB(I,J)=1
   10 CONTINUE
c       do 11 j=71,1,-1
c       read(30,300) (LSTAB(i,j),i=1,81)
c11     continue
 300    format(81I1)
      RETURN                                                          
      END                                                          



      SUBROUTINE UPDOWN 
      USE C57   
      REAL VW2(100),TOL(100),TARN(100)         
      IF (LAKE .EQ. 0) RETURN      
      BM = 1.95+K35         
      BM2 = BM**2             
      HLOG = 1.39-ALOG(HH)                                        
      DO 61 IA = 1,100     
        UM = VV(IA)*COST(IA)        
        VM = VV(IA)*SINT(IA)                                  
        VPLUS = -VV(IA)*UXV(IA,1)  
        VTOP = VM+VPLUS                
        VW2(IA) = UM**2+VTOP**2           
        TOL(IA) = 1E-4*VW2(IA)              
        TARN(IA) = UM*VPLUS/(UM**2+VM*VTOP)                       
   61 CONTINUE                                 
      DO 70 INCH = 1,LAKE                                  
        AZ = ZCOEFF(1,INCH)              
        BZ = ZCOEFF(2,INCH)        
        CZ = ZCOEFF(3,INCH)     
        UX(1) = 1.                         
        IF (AZ*BZ .NE. 0.) UX(1) = (.5*AZ/BZ)**(1./3.)       
        ZLOG = HLOG+ALOG(AZ/UX(1)+BZ*UX(1)**2+CZ)     
        DO 69 IA = 1,100                        
          UX(1) = K35*SQRT(VW2(IA)/(ZLOG**2+BM2))             
          ZLOG = HLOG+ALOG(AZ/UX(1)+BZ*UX(1)**2+CZ)   
          UV(1) = UX(1)**2/K2*(ZLOG**2+BM2)   
          DUV(1) = UV(1)-VW2(IA)         
          UX(2) = K35*SQRT(VW2(IA)/(ZLOG**2+BM2))        
          UV(2) = UX(2)**2/K2*(ZLOG**2+BM2)     
          DUV(2) = UV(2)-VW2(IA)          
          KSTOP = 0               
   62     IF (DUV(1).EQ.DUV(2)) GO TO 63                
          UX(3) = (UX(1)*DUV(2)-UX(2)*DUV(1))/(DUV(2)-DUV(1))      
          ZLOG = HLOG+ALOG(AZ/UX(3)+BZ*UX(3)**2+CZ)         
          UV(3) = UX(3)**2/K2*(ZLOG**2+BM2)         
          DUV(3) = UV(3)-VW2(IA)             
          IF (KSTOP .NE. 0) GO TO 63  
          IF (ABS(DUV(3)) .LT. TOL(IA)) KSTOP = 1 
          UX(1) = UX(2)     
          UX(2) = UX(3)                                   
          DUV(1) = DUV(2)      
          DUV(2) = DUV(3)        
          GO TO 62        
   63     UXV(IA,INCH+1) = K35*SQRT(VW2(IA)/(ZLOG**2+BM2))/VV(IA)
          TURN(IA,INCH+1) = (BM-ZLOG*TARN(IA))/(ZLOG+BM*TARN(IA))   
   69   CONTINUE              
   70 CONTINUE         
      RETURN    
      END 



      SUBROUTINE UXXV 
      USE C57   
      A125 = 1.25        
      DO 10 IA = 1,100 
        VV(IA) = FLOAT(IA)/A125 
   10 CONTINUE              
        DEN = G*K2*DTH     
        DO 40 IA = 1,100 
          IB = 101-IA    
          VV2 = VV(IB)**2      
          TOL = 1E-4*VV2       
          IF (IA .NE. 1) GO TO 20                                 
          CM = 2.55                       
          K123 = 1               
          UX(1) = 2.74          
          CALL ABCCC      
          K123 = 2             
          UX(2) = UX(1)*80./SQRT(UV(1))             
          GO TO 30                                                   
   20     UX(1) = UX(3)                                             
          UV(1) = UV(3)                                           
          DUV(1) = UV(1)-VV2                                       
          K123 = 2                         
          UX(2) = UX(1)*VV(IB)/VV(IB+1)  
   30     CALL ABCCC                         
          KSTOP = 0              
          K123 = 3                                      
   31     IF (DUV(1).EQ.DUV(2)) GO TO 32         
          UX(3) = AMAX1(.5*AMIN1(UX(1),UX(2)),         
     $      AMIN1(2.*AMAX1(UX(1),UX(2)),                  
     $      (UX(1)*DUV(2)-UX(2)*DUV(1))/(DUV(2)-DUV(1))))  
          CALL ABCCC                    
          IF (KSTOP .NE. 0) GO TO 32      
          IF (ABS(DUV(3)).LT. TOL) KSTOP = 1      
          UX(1) = UX(2)                 
          UX(2) = UX(3)        
          DUV(1) = DUV(2)         
          DUV(2) = DUV(3)                                  
          GO TO 31                     
   32     AB = SQRT((ZLOG+AM)**2+BM**2)    
          UXV(IB,1) = K35/AB             
          TURN(IB,1) = ATAN(BM/(ZLOG+AM))        
          COST(IB) = -(ZLOG+AM)/AB                
          SINT(IB) = BM/AB              
   40 CONTINUE                     
      RETURN             
      END                                                             



      SUBROUTINE AANGEL (I20)                           
      USE C3; USE C4; USE C5
C     CONVERTS U,V TO SPEED(KN), DIRECTION.                                     
C     COMPUTES SPEED AT 20 MTRS ALTITUDE IF I20 NE 0.                           
      DO 1000 NEST=1,5                                         
      DO 1000 I=1,21                                           
      DO 1000 J=1,21                                          
C COMPUTE VTN                                                                   
      VTN(I,J,NEST)=SQRT(UN(I,J,NEST)**2+VN(I,J,NEST)**2)             
      IF(VTN(I,J,NEST))  860,840,860                                    
 840  ANG(I,J,NEST)=0.0                                                 
      GO TO 1000                      
 860  AN=57.29578*ATAN2(VN(I,J,NEST),UN(I,J,NEST))                
      ANG(I,J,NEST)=AMOD(270.-AN, 360.)         
C     REDUCE SPEED TO 20 MTRS IF I20 NE 0                                       
C     REDUCE SPEED TO KNOTS IF I20 EQ 0                                         
      IF (I20 .NE. 0) GO TO 900                   
      VTN(I,J,NEST) = VTN(I,J,NEST)*(3600./1852.)             
      GO TO 1000                                 
 900  SPP = 1.25*VTN(I,J,NEST)                           
      KPP = SPP                              
      SPP = SPP-KPP                                
      LS = LW(I,J,NEST)        
      IF (KPP .NE. 0) GO TO 910                      
      UXX =VTN(I,J,NEST)*UXV(1,LS)         
      TWIST = TURN(1,LS)      
      GO TO 930                                           
 910  IF (KPP .GE. 100) GO TO 920            
      UXX = VTN(I,J,NEST)*((1.-SPP)*UXV(KPP,LS)+SPP*UXV(KPP+1,LS))
      TWIST = (1.-SPP)*TURN(KPP,LS)+SPP*TURN(KPP+1,LS)                
      GO TO 930                          
 920  UXX = VTN(I,J,NEST)*UXV(100,LS)                
      TWIST = TURN(100,LS)                                   
 930  Z20 = 10.0/Z0LAND                                   
      IF (LS .EQ. 2) Z20 = 10.0/(GA*UXX**2)                       
      VTN(I,J,NEST) = (3600./1852.)*AMIN1                           
     $  (UXX/K35*ALOG(Z20),VTN(I,J,NEST))                              
      ANG(I,J,NEST) = ANG(I,J,NEST)+57.29578*TWIST                     
      IF (ANG(I,J,NEST).LT. 0.) ANG(I,J,NEST) = ANG(I,J,NEST)+360.
 1000 CONTINUE                                     
      RETURN                                                       
      END                                                      



      SUBROUTINE ABCC 
      USE C4; USE C5
      IF (LS .EQ. 1) GO TO 53                                   
      Z0 = GA*UX(K123)**2/HH                                     
      ZLOG = ALOG(Z0)                                                
   53 HL = -HH*DEN/(UX(K123)**2*PTH*(ZLOG+CM))                        
      IF (HL+2.) 54,55,56                                            
   54 HLOG = ALOG(-1.*HL)                                                
      AM = AMIN1(HLOG+1.5,-.875*ZLOG)                                 
      BM = 1.8*EXP(.2*HL)                                               
      CM = AMIN1(HLOG+3.7,-.875*ZLOG)                     
      GO TO 60                           
   55 AM = 2.1931472                                        
      BM = 1.2065761               
      CM = 4.3931472                                              
      GO TO 60                                                   
   56 IF (HL-2.) 57,58,59                                   
   57 AM = 1.3865736-.4032868*HL                           
      BM = 1.953288+.37335598*HL                     
      CM = 2.5465736-HL*.9232868                                
      GO TO 60                                                    
   58 AM = .58                                                     
      BM = 2.70                                                       
      CM = .7                                                          
      GO TO 60                                                        
   59 YLOG = ZLOG+4.7                                             
      HL = .25*(YLOG+SQRT(YLOG**2+8.*HH*DEN/(UX(K123)**2*PTH)))  
      AM = AMAX1(2.5-.96*HL,-99.)                             
      BM = AMIN1(1.1+.8*HL,99.)                                        
      CM = 4.7-2.*HL                                                 
   60 UV(K123) = UX(K123)**2/K2*((ZLOG+AM)**2+BM**2)   
      DUV(K123) = UV(K123)-VV2              
      RETURN                
      END 



      SUBROUTINE BLOWUQ 
      USE C0; USE C1; USE C3; USE C4; USE C5
      INTEGER :: LEVEL(16) = (/1,5,1,2,1,3,1,2,1,4,1,2,1,3,1,2/)
      REAL CDH2(2)                                             
C     BLOWUP PRODUCES THE WIND FIELD FOR MOVING VORTEX IN THE PLANETARY         
C     BOUNDARY LAYER. ITS U(EASTWARD) AND V(NORTHWARD) COMPONENTS               
C     ARE IN ARRAYS UN,VN RESPECTIVELY UPON EXIT.                               
C   SGW IS MAGNITUDE OF SURFACE GEOSTROPHIC WIND                                
C   AN1 IS ANGLE BETWEEN SGW AND EAST, COUNTERCLOCKWISE FROM EAST               
C  ST12 IS DISTANCE IN KM FROM AXIS OF STORM TO HALF MAGNITUDE OF SGW           
C     (UG,VG) IS THE SFC GEOSTROPHIC WIND.                                      
C     CS IS SPEED OF STORM MOVEMENT.
      DEGG = AN1(NSNAP)/57.29578                                     
      UG=SGW*COS(DEGG)                                             
      VG=SGW*SIN(DEGG)       
      CALL SHORE                                                      
      CALL PXYM                                                     
      SPP = 1.25*SGW                                                   
      KPP = SPP                                                       
      SPP = SPP-KPP                                                 
      DO 1020 LS = 1,2                                       
        CDH2(LS) = (((1.-SPP)*UXV(KPP,LS)+SPP*UXV(KPP+1,LS))**2/HH)**2  
 1020 CONTINUE                                      
      CS=SQRT(UC**2+VC**2)                                     
      DO 1026 I = 1,21                                           
      DO 1026 J = 1,21                                          
      IF (I .EQ. 1) GO TO 1025                              
      IF (I .EQ. 21) GO TO 1025                                  
      IF (J .EQ. 1) GO TO 1025                                   
      IF (J .EQ. 21) GO TO 1025                               
      GO TO 1026                                                  
 1025 CONTINUE                                                    
      LS = LW(I,J,5)                                                   
      AA = CDH2(LS)*(PX(I,J,5)**2+PY(I,J,5)**2) 
      F2=F**2/2.0                                                   
      SK=AA/(F2+SQRT(F2**2+AA))                                        
      BK=SQRT(SK)                                                 
      V(I,J,5)=(F*PX(I,J,5)-BK*PY(I,J,5))/(F**2+SK)                  
      U(I,J,5)=-PY(I,J,5)/F-V(I,J,5)*BK/F             
      V(I,J,5)=V(I,J,5)-VC                                          
      U(I,J,5)=U(I,J,5)-UC 
                                                    
 1026 CONTINUE                                                
      DO 1028 I = 1,21                                          
      DO 1028 J = 1,21                                           
      UN(I,J,5)=U(I,J,5)                                            
      VN(I,J,5)=V(I,J,5)  
 1028 CONTINUE         

      DO 1030 NEST = 1,5                        
      DO 1030 I = 1,21                                    
      DO 1030 J = 1,21                                       
      PX(I,J,NEST)=PX(I,J,NEST)-F*VC                         
 1030 PY(I,J,NEST)=PY(I,J,NEST)+F*UC 
C                                                                               
      DO 1430 K=1,NM                                       
      KCK = MOD(K,16)                                            
      CALL COMQUT (LEVEL(KCK+1))                                  
 1430 CONTINUE                                                      
      CALL OUTFLO                                                     
C     SOLUTION OF WIND FIELD ON COORDINATE SYSTEM MOVING                        
C     WITH STORM IS NOW COMPLETE.                                               
C     COMPUTE SOLUTION WITH RESPECT TO SEA SURFACE.                             
      DO 1450 NEST=1,5                                           
      DO 1450 I=1,21                                                
      DO 1450 J=1,21                                                  
      UN(I,J,NEST)=UN(I,J,NEST)+UC                                   
      VN(I,J,NEST)=VN(I,J,NEST)+VC 
 1450 CONTINUE                                              
      RETURN                                                          
      END                                                            



      SUBROUTINE CCROSS                                               
      USE C4; USE C5
      FF = AMAX1(ABS(FLAT),1.832E-6)                             
      DO 10 IA = 1,100                                            
        VV(IA) = .8*FLOAT(IA)                                     
   10 CONTINUE                                                       
      DO 41 LS = 1,2                                                 
        DEN = G*K2*DTH(LS)                                              
        DO 40 IA = 1,100                                           
          IB = 101-IA                                                   
          VV2 = VV(IB)**2                                         
          IF (IA .NE. 1) GO TO 20                              
          CM = 2.55                                                
          K123 = 1                                             
          UX(1) = 2.74                                                
          IF (LS .EQ. 1) UX(1) = 3.38                                 
          Z0 = Z0LAND/HH                                            
          ZLOG = ALOG(Z0)                                           
          CALL ABCC                                                 
          K123 = 2                                                 
          UX(2) = UX(1)*80./SQRT(UV(1))                             
          GO TO 30                                                    
   20     UX(1) = UX(3)                                              
          UV(1) = UV(3)                                              
          DUV(1) = UV(1)-VV2                                       
          K123 = 2                                            
          UX(2) = UX(1)*VV(IB)/VV(IB+1)                             
   30     CALL ABCC                                                 
          KSTOP = 0                                               
          K123 = 3                                              
   31     IF (DUV(1).EQ.DUV(2)) GO TO 32                          
          UX(3) = AMAX1(.5*AMIN1(UX(1),UX(2)),                    
     $      AMIN1(2.*AMAX1(UX(1),UX(2)),                         
     $      (UX(1)*DUV(2)-UX(2)*DUV(1))/(DUV(2)-DUV(1))))         
          CALL ABCC                                              
          IF (KSTOP .NE. 0) GO TO 32                             
          IF (ABS(DUV(3)).LT. .0001*VV2) KSTOP = 1                   
          UX(1) = UX(2)                                   
          UX(2) = UX(3)                                           
          DUV(1) = DUV(2)                                         
          DUV(2) = DUV(3)                                         
          GO TO 31                                                 
   32     AB = SQRT((ZLOG+AM)**2+BM**2)                  
          UXV(IB,LS) = K35/AB                                     
          AB = UXV(IB,LS)**2/AB                                 
          CDR(IB,LS,1) = (ZLOG+AM)*AB                                   
          CDR(IB,LS,2) = BM*AB                                        
          TURN(IB,LS) = ATAN(BM/(ZLOG+AM))   
   40   CONTINUE                                                     
   41 CONTINUE                                                 
      RETURN                                                 
      END 



      SUBROUTINE COMQUT (LEVEL)
      USE C1; USE C3; USE C4; USE C5
      DIMENSION HKL(21,21),DRAG(2)
      DATA E1,E2/0.5,0.5/                                             
      DATA VONK/0.4/               
      SQQ(A,B) = SQRT(A*A+B*B)   ! vjp 3/10/06

      DO 800 NC = 1,LEVEL                                          
      NEST = LEVEL-NC+1                                          
      DO 720 I = 1,21                                               
      DO 720 J = 1,21                                               
      U(I,J,NEST) = UN(I,J,NEST)                                    
      V(I,J,NEST) = VN(I,J,NEST)                                    
  720 CONTINUE                                                        
      IF(NEST.NE.LEVEL) GO TO 721                                 
      IF(NEST.NE.5) CALL OUTBY2(NEST)                                
      GO TO 722                                                      
  721 CALL OUTBY1(NEST)                                             
  722 CONTINUE                                                     
      DXL = 2.0**(NEST-1)*DX*1000.0                                
      DXL2=DXL**2                                                 
      DTL=2.0**(NEST-1)*DT                                            
      FTL=F*DTL                                                    
      FTL1=FTL**2+1.0                                               
      FCL=2.0*VONK**2*(DXL/2.0)**2                                   
C  INNER BOUNDARY                                                               
      IF(NEST .EQ. 1) GO TO 730                           
      DO 724 J=1,10                                             
      U(7,J+6,NEST)=UN(3,2*J+1,NEST-1)                            
      V(7,J+6,NEST)=VN(3,2*J+1,NEST-1)                            
      V(15,J+6,NEST)=VN(19,2*J+1,NEST-1)                          
      U(15,J+6,NEST)=UN(19,2*J+1,NEST-1)                              
 724  CONTINUE                                                  
      DO 726 I=1,10                                                 
      U(I+6,7,NEST)=UN(2*I+1,3,NEST-1)                             
      V(I+6,7,NEST)=VN(2*I+1,3,NEST-1)                           
      U(I+6,15,NEST)=UN(2*I+1,19,NEST-1)                             
      V(I+6,15,NEST)=VN(2*I+1,19,NEST-1)                          
 726  CONTINUE                                                        
C  COMPUTATION OF INTERIOR POINT                                                
 730  DO 734 I=1,20                                            
      DO 734 J=1,20                                                  
      IF (NEST .EQ. 1) GO TO 733                                 
      IF(I .LE. 6 .OR. I .GE. 15) GO TO 733                       
      IF(J .LE. 6 .OR. J .GE. 15) GO TO 733                        
      HKL(I,J)=0.0                                                     
      GO TO 734                                                 
 733  D1=0.5*(U(I+1,J,NEST)-U(I,J,NEST)+U(I+1,J+1,NEST)-U(I,J+1,NEST)   
     1   -V(I,J+1,NEST)+V(I,J,NEST)-V(I+1,J+1,NEST)+V(I+1,J,NEST))/DXL
      D2=0.5*(V(I+1,J,NEST)-V(I,J,NEST)+V(I+1,J+1,NEST)-V(I,J+1,NEST) 
     1   +U(I,J+1,NEST)-U(I,J,NEST)+U(I+1,J+1,NEST)-U(I+1,J,NEST))/DXL  
      HKL(I,J) = FCL*SQQ(D1,D2)                                  
 734  CONTINUE                                                  
      DO 775 I=2,20                                                
      DO 775 J=2,20                                             
      IF( NEST .EQ. 1) GO TO 736                                     
      IF(I .LE. 6 .OR. I .GE. 16) GO TO 736                           
      IF(J .LE. 6 .OR. J .GE. 16) GO TO 736                        
      UN(I,J,NEST)=U(I,J,NEST)                                         
      VN(I,J,NEST)=V(I,J,NEST)       
      GO TO 775                                        
 736  U1=U(I,J,NEST)+UC                                           
      V1=V(I,J,NEST)+VC                                                
C DRAG(1) IS TANGENTIAL DRAG CORRECTION TERM                                    
C DRAG(2) IS NORMAL DRAG CORRECTION TERM                                        
      LS=LW(I,J,NEST)                        
      SPP=SQQ(U1,V1)                                           
      SPP1 = MAX(1.,MIN(99.99,1.25*SPP))    
      KPP=SPP1                   
      SPP1=SPP1-KPP       
      SPH = SPP/HH              
      DO 737 IA = 1,2    
        DRAG(IA) = SPH*(CDR(KPP,LS,IA)+                     
     $    SPP1*(CDR(KPP+1,LS,IA)-CDR(KPP,LS,IA)))         
  737 CONTINUE                                                    
      IF(NEST .NE. 5) GO TO 741                                   
      IF(I .NE. 2 .AND. I .NE. 20) GO TO 741                      
      IF(I-2) 738,738,740                                         
 738  IF(U(I,J,NEST)) 739,739,741          
 739  UKX=0.0                                                   
      VKX=0.0                                                    
      GO TO 742                                                  
 740  IF(U(I,J,NEST)) 741,739,739                            
 741  UKX=0.5*((HKL(I,J-1)+HKL(I,J))*(U(I+1,J,NEST)-U(I,J,NEST))     
     1   -(HKL(I-1,J-1)+HKL(I-1,J))*(U(I,J,NEST)-U(I-1,J,NEST)))/DXL2   
      VKX=0.5*((HKL(I,J-1)+HKL(I,J))*(V(I+1,J,NEST)-V(I,J,NEST))        
     1   -(HKL(I-1,J-1)+HKL(I-1,J))*(V(I,J,NEST)-V(I-1,J,NEST)))/DXL2  
      IF(NEST .NE. 5) GO TO 746                                       
 742  IF(J .NE. 2 .AND. J .NE. 20) GO TO 746                     
      IF(J-2) 743,743,745                                   
 743  IF(V(I,J,NEST)) 744,744,746                          
 744  UKY=0.0                                             
      VKY=0.0                                                  
      GO TO 747                                                
 745  IF(V(I,J,NEST)) 746,744,744                            
 746  UKY=0.5*((HKL(I-1,J)+HKL(I,J))*(U(I,J+1,NEST)-U(I,J,NEST))       
     1   -(HKL(I-1,J-1)+HKL(I,J-1)) *(U(I,J,NEST)-U(I,J-1,NEST)))/DXL2 
      VKY=0.5*((HKL(I-1,J)+HKL(I,J))*(V(I,J+1,NEST)-V(I,J,NEST))   
     1   -(HKL(I-1,J-1)+HKL(I,J-1)) *(V(I,J,NEST)-V(I,J-1,NEST)))/DXL2 
 747  IF(U(I,J,NEST)) 748,750,750                    
 748  WX=-U(I,J,NEST)*(U(I,J,NEST)-U(I+1,J,NEST))/DXL   
      VX=-U(I,J,NEST)*(V(I,J,NEST)-V(I+1,J,NEST))/DXL  
      GO TO 752                              
 750  WX= U(I,J,NEST)*(U(I,J,NEST)-U(I-1,J,NEST))/DXL 
      VX= U(I,J,NEST)*(V(I,J,NEST)-V(I-1,J,NEST))/DXL       
 752  IF(V(I,J,NEST)) 754,756,756                           
 754  UY=-V(I,J,NEST)*(U(I,J,NEST)-U(I,J+1,NEST))/DXL     
      VY=-V(I,J,NEST)*(V(I,J,NEST)-V(I,J+1,NEST))/DXL
      GO TO 758                                                  
 756  UY= V(I,J,NEST)*(U(I,J,NEST)-U(I,J-1,NEST))/DXL         
      VY= V(I,J,NEST)*(V(I,J,NEST)-V(I,J-1,NEST))/DXL                
 758  UR=0.5*(U(I,J,NEST)+V(I,J,NEST))             
      VR=0.5*(-U(I,J,NEST)+V(I,J,NEST))         
      IF(UR) 760,762,762                                     
 760  UX1=-UR*(U(I,J,NEST)-U(I+1,J+1,NEST))/DXL     
      VX1=-UR*(V(I,J,NEST)-V(I+1,J+1,NEST))/DXL      
      GO TO 764                                            
 762  UX1= UR*(U(I,J,NEST)-U(I-1,J-1,NEST))/DXL               
      VX1= UR*(V(I,J,NEST)-V(I-1,J-1,NEST))/DXL                    
 764  IF(VR) 766,768,768                    
 766  UY1=-VR*(U(I,J,NEST)-U(I-1,J+1,NEST))/DXL                 
      VY1=-VR*(V(I,J,NEST)-V(I-1,J+1,NEST))/DXL           
      GO TO 770                                         
 768  UY1= VR*(U(I,J,NEST)-U(I+1,J-1,NEST))/DXL   
      VY1= VR*(V(I,J,NEST)-V(I+1,J-1,NEST))/DXL           
 770  UXY=E1*(WX+UY)+E2*(UX1+UY1)                                
      VXY=E1*(VX+VY)+E2*(VX1+VY1)                              
      B1 = U(I,J,NEST)+DTL*                                     
     $  (UKX+UKY-PX(I,J,NEST)-UXY+U1*DRAG(1)+V1*DRAG(2))         
      B2 = V(I,J,NEST)+DTL*                                   
     $  (VKX+VKY-PY(I,J,NEST)-VXY+V1*DRAG(1)-U1*DRAG(2))           
      UN(I,J,NEST)=(B1+FTL*B2)/FTL1                             
      VN(I,J,NEST)=(B2-FTL*B1)/FTL1
 775  CONTINUE                                     
 800  CONTINUE                                                    
      RETURN                                                      
      END 



      SUBROUTINE GRAD                                           
      USE C1; USE C2
C     GRAD COMPUTES PRESSURE GRADIENT
      DIMENSION AD(21,4),AP(3,2),AT(2),BA(15),BC(4,2),BP(2,2)
      EQUIVALENCE (BC,BA(8))                                   
C                                                                               
C     IB IS SWITCH VARIABLE. IB = 0 TO SUPPRESS PRINTING.                       
C     ITRACK IS INDICATOR FOR CODING OF DIREC                                   
C     EYELAT IS NORTH LATITUDE OF EYE IN DEGREES                                
C              (SOUTH LATITUDE MUST HAVE MINUS SIGN)                            
C     EYLONG IS EAST LONGITUDE OF EYE IN DEGREES                                
C              (WEST LONGITUDE MUST HAVE MINUS SIGN)                            
C     DIREC IS DIRECTION OF TRACK OF HURRICANE, CLOCKWISE FROM NORTH            
C           ITRACK = 0, DIREC IN DEGREES                                        
C           ITRACK = 1, DIREC IN POINTS OF 11.25 DEG                            
C     SPEED IS FORWARD SPEED OF HURRICANE IN KNOTS                              
C     IQUAD IS INDICATOR FOR QUADRANTS OF PRESSURE FIELD                        
C           IQUAD = 0, CIRCULARLY SYMMETRIC PRESSURE FIELD                      
C           IQUAD = 1, FIRST QUADRANT IS RIGHT FRONT                            
C           IQUAD = 2, FIRST QUADRANT IS FORWARD                                
C           QUADRANTS FOLLOW CLOCKWISE FROM FIRST                               
C     EYPRES IS PRESSURE AT EYE IN MILLIBARS                                    
C     RADIUS(1,2,3,4) ARE R IN FOUR QUADRANTS IN UNITS OF 1.0 NM                
C                                                       (1852 METERS)           
C     PFAR(1,2,3,4) ARE FAR FIELD PRESSURE IN FOUR QUADRANTS,                   
C         IN MILLIBARS                                                          
C     IF IQUAD = 0, ENTER RADIUS(1) AND PFAR(1) ONLY                            
C     AB(I+11,J+11,K+1) IS DP/DX (EAST) IN MB/KM AT 5*I*2**K KM EAST OF         
C                                                 AND 5*J*2**K KM SOUTH         
C     AB(I+11,J+11,K+6) IS DP/DY (NORTH) IN MB/KM AT SAME POINT                 
C     AB(I+11,J+11,K+11) IS DP/DR (OUTWARD) IN MB/KM AT SAME POINT              
C     AC(I+11,1,J+11) IS COS OF ANGLE OF POINT I,J                              
C                                       (SAME FOR ALL GRID NESTS)               
C     AC(I+11,2,J+11) IS SIN OF ANGLE OF POINT I,J                              
C                                       (SAME FOR ALL GRID NESTS)               
C     AC(I+11,3...7,J+11) IS RADIUS OF POINT I,J IN NEST 1...5                  
C                                       RESPECTIVELY (METERS)                   
      DEG = .017453293   
      S45 = .70710678                                             
C     COMPUTE GRID ANGLES AND DISTANCES (FUNCTION OF DX ONLY)                   
      DO 10 IG = 1,21                                             
 10   AD(IG,1) = DX*FLOAT(IG-11)                                   
      DO 12 IC = 1,21                                                
      DO 12 ID = 1,21                                               
      AE = AD(IC,1)**2+AD(ID,1)**2                                 
      IF (AE.LE.0) GO TO 12                                         
      AC(IC,3,ID) = SQRT(AE)                                         
      AC(IC,1,ID) = AD(IC,1)/AC(IC,3,ID)                             
      AC(IC,2,ID) = -AD(ID,1)/AC(IC,3,ID)                             
      DO 11 IE = 4,7                                                   
 11   AC(IC,IE,ID) = AC(IC,IE-1,ID)+AC(IC,IE-1,ID)                      
 12   CONTINUE                                                     
C     PREVENT DIVISION BY ZERO LATER IN SUBROUTINE                              
      AC(11,3,11) = 0.0                                              
      BA(4) = DIREC*DEG                                                
      IF (JA(1).NE. 0) BA(4) = BA(4)*11.25                          
      BA(5) = SPEED*(1852./3600.)                                    
      IF (JA(6).EQ. 0) GO TO 280                                     
      DO 21 IC = 1,4                                                 
        BA(IC+7) = RADIUS(IC)*1.852                                
        BA(IC+11) = PFAR(IC)-EYPRES                                  
 21   CONTINUE                                                        
 24   DO 25 IC=1,2                                             
 25   AP(1,IC) = .25*(BC(1,IC)+BC(2,IC)+BC(3,IC)+BC(4,IC))              
      IF (JA(6).EQ. 2) GO TO 27                                      
      DO 26 IC = 1,2                                                 
      AP(2,IC) = .5*S45*(BC(1,IC)+BC(2,IC)-BC(3,IC)-BC(4,IC))       
 26   AP(3,IC) = .5*S45*(BC(1,IC)-BC(2,IC)-BC(3,IC)+BC(4,IC))      
      GO TO 29                                             
 27   DO 28 IC = 1,2                                               
      AP(2,IC) = .5*(BC(2,IC)-BC(4,IC))                     
      AP(3,IC) = .5*(BC(1,IC)-BC(3,IC))                           
 28   CONTINUE                                                     
      GO TO 29                                                     
 280  CONTINUE                                                    
      AP(1,1) = RADIUS(1)*1.852                                        
      AP(1,2) = PFAR(1)-EYPRES                                        
 29   CONTINUE                                                       
      AT(1) = COS(BA(4))                                          
      AT(2)= -SIN(BA(4))                                           
      CRR=AP(2,1)*AT(1)-AP(3,1)*AT(2)                              
      CRI=AP(2,1)*AT(2)+AP(3,1)*AT(1)                                 
      CDR=AP(2,2)*AT(1)-AP(3,2)*AT(2)                                  
      CDI=AP(2,2)*AT(2)+AP(3,2)*AT(1)                               
      AP(2,1)=CRR                                                   
      AP(3,1)=CRI                                                    
      AP(2,2)=CDR                                                    
      AP(3,2)=CDI                                                 
      UC=-AT(2)*BA(5)                                                
      VC = AT(1)*BA(5)                                               
 31   DO 40 IE = 1,5                                                 
      DO 38 ID = 1,21                                             
      DO 35 IC = 1,21                                              
      IF (JA(6).EQ. 0) GO TO 340                                       
      DO 32 IH = 1,2                              
      BP(1,IH) = AP(1,IH)+AP(2,IH)*AC(IC,1,ID)+AP(3,IH)*AC(IC,2,ID)   
   32 BP(2,IH) = -AP(2,IH)*AC(IC,2,ID)+AP(3,IH)*AC(IC,1,ID)            
      IF (AC(IC,3,ID).GT. 0) GO TO 33                        
      AD(IC,1) = 0.                                             
      AD(IC,2) = 0.                                        
      GO TO 34                                                   
   33 AE = EXP(-BP(1,1)/AC(IC,IE+2,ID))                              
      AF = AE*BP(1,2)                                                 
C     COMPUTE RADIAL PRESSURE GRADIENT                                          
      AD(IC,1) = AF*BP(1,1)/AC(IC,IE+2,ID)**2                     
      AD(IC,2) = (AE*BP(2,2)-AF/AC(IC,IE+2,ID)*BP(2,1))/AC(IC,IE+2,ID)  
C     TANGENTIAL PRESSURE GRADIENT                                              
   34 AD(IC,3) = AD(IC,1)*AC(IC,1,ID)-AD(IC,2)*AC(IC,2,ID)          
      AD(IC,4) = AD(IC,1)*AC(IC,2,ID)+AD(IC,2)*AC(IC,1,ID)        
      GO TO 341                                                        
  340 AD(IC,2) = 0.                                               
C     CIRCULARLY SYMMETRIC PRESSURE FIELD# ZERO TANGENTIAL GRADIENT             
      IF (AC(IC,3,ID) .GT. 0.0) GOTO 377                           
      AD(IC,1) = 0.0                                              
      AD(IC,3) = 0.0                                                   
      AD(IC,4) = 0.0                                                    
      GOTO 341                                                         
  377 AD(IC,1) = EXP(-AP(1,1)/AC(IC,IE+2,ID))*AP(1,2)*AP(1,1)/       
     &AC(IC,IE+2,ID)**2                                            
      AD(IC,3) = AD(IC,1)*AC(IC,1,ID)                             
      AD(IC,4) = AD(IC,1)*AC(IC,2,ID)                                
  341 CONTINUE                                                   
      AB(IC,ID,IE+10) = AD(IC,1)                                   
      AB(IC,ID,IE) = AD(IC,3)                                     
   35 AB(IC,ID,IE+5) = AD(IC,4)  
      IF (IB .EQ. 0) GO TO 38                                         
C     IF NO PROGRAM CHANGES, NORTH IS PRINTED AT TOP                            
C     OF PAGE, WEST AT LEFT, ETC.                                               
C     WRITE (LU,36) AD                                                          
   36 FORMAT (4(1X,1P21F6.2/))                                       
   38 CONTINUE                                                        
      IF (IB .EQ. 0) GO TO 40                                         
C     WRITE (LU,39)                                                             
   39 FORMAT (/////)                                                
   40 CONTINUE                                                       
      RETURN                                                       
      END



      SUBROUTINE OUTBY1(NEST)                                        
      USE C3
C-----------------------------------------------------------------------------
C  OUTER BOUNDARY ( NOT AT THE SAME TIME LEVEL)                                 
C-----------------------------------------------------------------------------
       UN(1,1,NEST)=0.5*(UN(6,6,NEST+1)+U(6,6,NEST+1))                 
       VN(1,1,NEST)=0.5*(VN(6,6,NEST+1)+V(6,6,NEST+1))                  
       DO 780 J=1,10                                                  
       UN(1,2*J+1,NEST)=0.5*(UN(6,J+6,NEST+1) +U(6,J+6,NEST+1))      
       VN(1,2*J+1,NEST)=0.5*(VN(6,J+6,NEST+1) +V(6,J+6,NEST+1))         
       UN(21,2*J+1,NEST)=0.5*(UN(16,J+6,NEST+1)+U(16,J+6,NEST+1))      
       VN(21,2*J+1,NEST)=0.5*(VN(16,J+6,NEST+1)+V(16,J+6,NEST+1))     
       UN(1 ,2*J,NEST)=0.250*(UN(6,J+5,NEST+1)+UN(6,J+6,NEST+1)         
     1               + U(6,J+5,NEST+1)+U(6,J+6,NEST+1))                 
       VN(1 ,2*J,NEST)=0.250*(VN(6,J+5,NEST+1)+VN(6,J+6,NEST+1)         
     1               + V(6,J+5,NEST+1)+V(6,J+6,NEST+1))                 
       UN(21,2*J,NEST)=0.250*(UN(16,J+5,NEST+1)+UN(16,J+6,NEST+1)       
     1               + U(16,J+5,NEST+1)+U(16,J+6,NEST+1))              
       VN(21,2*J,NEST)=0.250*(VN(16,J+5,NEST+1)+VN(16,J+6,NEST+1)      
     1               + V(16,J+5,NEST+1)+V(16,J+6,NEST+1))              
  780  CONTINUE                                                       
       DO 790 I=1,10                                                  
       UN(2*I+1,1 ,NEST)=0.5*(UN(I+6,6,NEST+1)+U(I+6,6,NEST+1))       
       VN(2*I+1,1 ,NEST)=0.5*(VN(I+6,6,NEST+1)+V(I+6,6,NEST+1))      
       VN(2*I+1,21,NEST)=0.5*(VN(I+6,16,NEST+1)+V(I+6,16,NEST+1))    
       UN(2*I+1,21,NEST)=0.5*(UN(I+6,16,NEST+1)+U(I+6,16,NEST+1))    
       UN(2*I,1 ,NEST)=0.250*(UN(I+5,6,NEST+1)+UN(I+6,6,NEST+1)    
     1               + U(I+5,6,NEST+1)+U(I+6,6,NEST+1))              
       VN(2*I,1 ,NEST)=0.250*(VN(I+5,6,NEST+1)+VN(I+6,6,NEST+1)        
     1               + V(I+5,6,NEST+1)+V(I+6,6,NEST+1))                
       UN(2*I,21,NEST)=0.250*(UN(I+5,16,NEST+1)+UN(I+6,16,NEST+1)     
     1               + U(I+5,16,NEST+1)+U(I+6,16,NEST+1))            
       VN(2*I,21,NEST)=0.250*(VN(I+5,16,NEST+1)+VN(I+6,16,NEST+1)     
     1               + V(I+5,16,NEST+1)+V(I+6,16,NEST+1))          
  790  CONTINUE                                                       
       RETURN                                                         
       END                                                            



      SUBROUTINE OUTBY2(NEST) 
      USE C3
C-----------------------------------------------------------------------------
C     OUTER BOUNDARY (AT SAME TIME LEVEL)                                       
C-----------------------------------------------------------------------------
       UN(1,1,NEST)=UN(6,6,NEST+1)                                 
       VN(1,1,NEST)=VN(6,6,NEST+1)                                   
       DO 820 J=1,10                                                
       UN(1 ,2*J+1,NEST)=UN(6,J+6,NEST+1)                          
       VN(1 ,2*J+1,NEST)=VN(6,J+6,NEST+1)                             
       VN(21,2*J+1,NEST)=VN(16,J+6,NEST+1)                           
       UN(21,2*J+1,NEST)=UN(16,J+6,NEST+1)                           
       UN(1 ,2*J,NEST)=0.5*(UN(6,J+5,NEST+1)+UN(6,J+6,NEST+1))       
       VN(1 ,2*J,NEST)=0.5*(VN(6,J+5,NEST+1)+VN(6,J+6,NEST+1))       
       VN(21,2*J,NEST)=0.5*(VN(16,J+5,NEST+1)+VN(16,J+6,NEST+1))       
       UN(21,2*J,NEST)=0.5*(UN(16,J+5,NEST+1)+UN(16,J+6,NEST+1))      
  820  CONTINUE                                                   
       DO 830 I=1,10                                                
       UN(2*I+1,1 ,NEST)=UN(I+6,6,NEST+1)                           
       VN(2*I+1,1 ,NEST)=VN(I+6,6,NEST+1)                           
       VN(2*I+1,21,NEST)=VN(I+6,16,NEST+1)                            
       UN(2*I+1,21,NEST)=UN(I+6,16,NEST+1)                            
       UN(2*I,1 ,NEST)=0.5*(UN(I+5,6,NEST+1)+UN(I+6,6,NEST+1))     
       VN(2*I,1 ,NEST)=0.5*(VN(I+5,6,NEST+1)+VN(I+6,6,NEST+1))      
       VN(2*I,21,NEST)=0.5*(VN(I+5,16,NEST+1)+VN(I+6,16,NEST+1))     
       UN(2*I,21,NEST)=0.5*(UN(I+5,16,NEST+1)+UN(I+6,16,NEST+1))     
  830  CONTINUE                                                      
       RETURN                                                         
       END                                                            



      SUBROUTINE OUTFLO                                           
      USE C3
      CO8 = .99026807                              
      SI8 = .1391731                              
      DO 10 NEST=1,5                                                   
      DO 10 IC=1,21                                                   
      DO 10 IR=1,21                                                     
        XX=UN(IR,IC,NEST)*CO8+VN(IR,IC,NEST)*SI8                       
        VN(IR,IC,NEST)=VN(IR,IC,NEST)*CO8-UN(IR,IC,NEST)*SI8          
        UN(IR,IC,NEST)=XX
   10 CONTINUE                                                          
      RETURN                                                           
      END                                                             



      SUBROUTINE PXYM                                              
      USE C1; USE C2; USE C3
      REAL,POINTER,DIMENSION(:,:,:) :: BA,BC
C-----------------------------------------------------------------------------
C     PXYM CALLS SUBROUTINE GRAD TO GET PRESSURE GRADIENT,                      
C     THEN REARRANGES THE PRESSURE GRADIENT, THEN PRODUCES                      
C     INITIAL GRADIENT WIND FIELD.   
C-----------------------------------------------------------------------------
      BA => AB
      BC => AC
C-----------------------------------------------------------------------------
      DO 50 NEST = 1,5                                                  
      GO TO (10,30,30,30,30),NEST                                       
   10 CALL GRAD                                                        
      AN2 = AN1(NSNAP)*3.14159265/180. 
      IF (ST12 .LE. 0) GO TO 30                                     
      DX2 = .5*DX/ST12                                               
      CO2 = COS(AN2)*DX2                                               
      SI2 = SIN(AN2)*DX2      
   30 DO 31 I = 1,21                                                    
      DO 31 J = 1,21                                                    
      MJ=22-J                                                        
      PX(I,J,NEST) = BA(I,MJ,NEST)*(1E-4/1.15E-3)                    
      PY(I,J,NEST) = BA(I,MJ,NEST+5)*(1E-4/1.15E-3) 
      AG = BA(I,MJ,NEST+10)*(1E-4/1.15E-3)                             
      AH=BC(I,NEST+2,MJ)*F*500.                                      
      AI=AG*BC(I,NEST+2,MJ)*1000.                                   
      IF(AI.NE.0.) AI = AI/(AH+SQRT(AH*AH+AI))                     
      U(I,J,NEST) = -AI*BC(I,2,MJ)                                 
      V(I,J,NEST) = AI*BC(I,1,MJ)                                   
      UN(I,J,NEST) = U(I,J,NEST)                                    
      VN(I,J,NEST) = V(I,J,NEST)       
   31 continue
C     SGW NON-ZERO FOR STEERING FLOW.                       
   32 IF (SGW.EQ.0.) GO TO 50                                      
      AG = F*UG                                                   
      AH = F*VG                                                    
      IF (ST12 .GT. 0) GO TO 34                                 
      DO 33 I = 1,21                                                
      DO 33 J = 1,21                                                
      PX(I,J,NEST) = PX(I,J,NEST)+AH                                    
      PY(I,J,NEST) = PY(I,J,NEST)-AG   
   33 continue
      GO TO 50                                          
   34 CO2 = CO2+CO2                                              
      SI2 = SI2+SI2                                               
      DO 36 I = 1,21                                            
      AI = (I-11)*SI2                                       
      DO 35 J = 1,21                                             
      AJ = (J-11)*CO2                                                
      FADE = -.69314718*(AJ-AI)**2                                 
      FADE = EXP(FADE)                                           
      IF (FADE .LE. 0.) GO TO 35                                     
      PX(I,J,NEST) = PX(I,J,NEST)+AH*FADE                            
      PY(I,J,NEST) = PY(I,J,NEST)-AG*FADE       
   35 CONTINUE                                         
   36 CONTINUE                                
   50 CONTINUE                                                   
   99 RETURN                          
      END 



      SUBROUTINE SHORE 
      USE C3
      DO 20 I=1,21                                             
      DO 20 J=1,21                                                 
      DO 20 K=1,5                                               
   20 LW(I,J,K)=2                                               
      RETURN                                                       
      END


      SUBROUTINE TVEL(VTN,ANG,IDENT,KSEQ,LV,I20)          
      USE C1
      CHARACTER*4 KREDUC, IDENT   ! debug
C     IF IN ARRAYS VTN,ANG 1ST DIMENSION INCREASES EASTWARD AND                 
C     2ND DIMENSION INCREASES NORTHWARD, SUBROUTINE TVEL PRINTS                 
C     WEST AT TOP OF PAGE,NORTH AT RIGHT OF PAGE, ETC.                          
C     NEST LV IS PRINTED AS INNER NEST, LV+1 AS OUTER NEST.
      DIMENSION VTN(21,21,5),ANG(21,21,5)
      DIMENSION KREDUC(2)                                         
      DATA KREDUC/'NOT ', '    '/ 
      DATA LU/6/                                         
      WRITE (LU+CK_LUN,100)                                               
  100 FORMAT(1H1)                                                  
      LVP=LV+1                                                       
      I20P=I20+1                                                    
      WRITE (LU+CK_LUN,200) IDENT,KSEQ,LVP,KREDUC(I20P)                     
  200 FORMAT(1H+,20X,A6,5X,A4,1X,I4,'   LEVEL',I3,5X,A4,'REDUCED')  
      DXI=2.0**(LV-1)*DS
      DXO=2.0*DXI
      WRITE(LU+CK_LUN,210) DXI,DXO
  210 FORMAT(21X,'( INNER-GRID STEP =',F6.1,' KM',6X,'OUTER-GRID STEP =
     1       ',F6.1,' KM',' )')
 1190 WRITE (LU+CK_LUN,1200) (J, J=1,12)                                      
 1200 FORMAT(1H0,11X,12(I2,8X))                                        
      DO 1210 I=1,5                                                     
      WRITE(LU+CK_LUN,1220) I,(VTN(I,J,LVP),J=1,12),(ANG(I,J,LVP),J=1,12)     
 1210 continue
 1220 FORMAT(1H0,6X,I2,12(1PF5.0,5X)/9X,12(0PF5.0,5X),/1H0,/)     
 1240 DO 1248 I=1,11                                              
      IL=I+5                                                      
      IO=2*I-1                                                          
      IE=2*I                                                           
      IF(I-11) 1242,1246,1246                                           
 1242 WRITE (LU+CK_LUN,1244) IL,(VTN(IL,J,LVP),J=1,5),(VTN(IO,N,LV),N=1,13)  
     1 ,(ANG(IL,J,LVP),J=1,5),(ANG(IO,N,LV),N=1,13)                     
     2 ,(VTN(IE,N,LV),N=1,13),(ANG(IE,N,LV),N=1,13)                  
 1244 FORMAT(1H0,6X,I2,5(1PF5.0,5X),13(1PF5.0)/9X,5(0PF5.0,5X),       
     1  13(0PF5.0)/1H0,58X,13(1PF5.0)/59X,13(0PF5.0))           
      GO TO 1248                                             
 1246 WRITE (LU+CK_LUN,1247) IL,(VTN(IL,J,LVP),J=1,5),(VTN(IO,N,LV),N=1,13)    
     1  ,(ANG(IL,J,LVP),J=1,5),(ANG(IO,N,LV),N=1,13)       
 1247 FORMAT(1H0,6X,I2,5(1PF5.0,5X),13(1PF5.0)/9X,5(0PF5.0,5X),   
     1  13(0PF5.0)/1H0,//)                            
 1248 CONTINUE                                          
      DO 1300 I=17,21                                          
      WRITE (LU+CK_LUN,1220) I,(VTN(I,J,LVP),J=1,12),(ANG(I,J,LVP),J=1,12)   
 1300 continue
      WRITE(LU+CK_LUN,100)                                                    
      WRITE (LU+CK_LUN,1340) (J,J=13,21)                                  
 1340 FORMAT(1H0,16X,9(I2,8X))                                        
      DO 1345 I=1,5                                                    
      WRITE (LU+CK_LUN,1350) (VTN(I,J,LVP),J=13,21),I,(ANG(I,J,LVP),J=13,21)  
 1345   continue
 1350 FORMAT(1H0,13X,9(1PF5.0,5X),I2/14X,9(0PF5.0,5X),/1H0,/)      
 1370 DO 1378 I=1,11                                              
      IL=I+5                                                      
      IO=2*I-1                                                      
      IE=2*I                                                          
      IF(I-11) 1372,1376,1376                                          
 1372 WRITE (LU+CK_LUN,1374) (VTN(IO,N,LV),N=14,21),(VTN(IL,J,LVP),J=17,21),IL 
     1  ,(ANG(IO,N,LV),N=14,21), (ANG(IL,J,LVP),J=17,21),               
     2  (VTN(IE,N,LV),N=14,21),(ANG(IE,N,LV),N=14,21)                  
 1374 FORMAT(1H0,8X,8(1PF5.0),5X,5(1PF5.0,5X),I2/9X,8(0PF5.0),5X,   
     1  5(0PF5.0,5X)/1H0,8X, 8(1PF5.0)/9X, 8(0PF5.0))                
      GO TO 1378                                                   
 1376 WRITE(LU+CK_LUN,1377) (VTN(IO,N,LV),N=14,21),(VTN(IL,J,LVP),J=17,21),IL, 
     1  (ANG(IO,N,LV),N=14,21), (ANG(IL,J,LVP),J=17,21)                
 1377 FORMAT(1H0,8X,8(1PF5.0),5X,5(1PF5.0,5X),I2/9X,8(0PF5.0),5X,     
     1  5(0PF5.0,5X)/1H0,//)                                       
 1378 CONTINUE                                                        
      DO 1400 I=17,21                                                  
      WRITE(LU+CK_LUN,1350)  (VTN(I,J,LVP),J=13,21),I,(ANG(I,J,LVP),J=13,21)   
 1400   continue
      RETURN                                                          
      END 
      
      

      SUBROUTINE BLEND
      PARAMETER (MAXI=11,MAXJ=11)
      DIMENSION WIND(2,MAXI,MAXJ),WIND7(2,MAXI,MAXJ),WIND8(2,MAXI,MAXJ)
C
C     BLEND STANDARD GRID WINDS WITH 25-KM GRID WINDS, AS NEEDED
C
C     FOR007 --> STANDARD GRID WINDS
C     FOR008 --> 25-KM GRID WINDS
C
      REWIND 7
      REWIND 8
      NW=0
      IG21=MAXI/21
      IR21=MOD(MAXI,21)
   10 READ(7+CK_LUN,'(I10)',END=140) NTOD
      READ(8+CK_LUN,'(I10)') NTOD
      NW=NW+1
      DO 60 NN=1,2
        IF (IG21.EQ.0) GOTO 40
        DO 30 IK=1,IG21
          IB=(IK-1)*21+1
          IE=IK*21
          DO 20 J=MAXJ,1,-1
            READ(7+CK_LUN,'(21F6.1)') (WIND7(NN,I,J),I=IB,IE)
            READ(8+CK_LUN,'(21F6.1)') (WIND8(NN,I,J),I=IB,IE)
   20     CONTINUE
   30   CONTINUE
   40   IF (IR21.EQ.0) GOTO 60
        IB=IG21*21+1
        IE=MAXI
        DO 50 J=MAXJ,1,-1
          READ(7+CK_LUN,'(21F6.1)') (WIND7(NN,I,J),I=IB,IE)
          READ(8+CK_LUN,'(21F6.1)') (WIND8(NN,I,J),I=IB,IE)
   50   CONTINUE
   60 CONTINUE
      NB=0
      DO 80 I=1,MAXI
        DO 70 J=1,MAXJ
          IF (WIND7(1,I,J).NE.0.) THEN
            WIND(1,I,J)=WIND7(1,I,J)
            WIND(2,I,J)=WIND7(2,I,J)
          ELSE
            NB=NB+1
            WIND(1,I,J)=WIND8(1,I,J)
            WIND(2,I,J)=WIND8(2,I,J)
          ENDIF
   70   CONTINUE
   80 CONTINUE
      WRITE(11+CK_LUN,'(I10)') NTOD
      DO 130 NN=1,2
        IF (IG21.EQ.0) GOTO 110
        DO 100 IK=1,IG21
          IB=(IK-1)*21+1
          IE=IK*21
          DO 90 J=MAXJ,1,-1
            WRITE(11+CK_LUN,'(21F6.1)') (WIND(NN,I,J),I=IB,IE)
   90     CONTINUE
  100   CONTINUE
  110   IF (IR21.EQ.0) GOTO 130
        IB=IG21*21+1
        IE=MAXI
        DO 120 J=MAXJ,1,-1
          WRITE(11+CK_LUN,'(21F6.1)') (WIND(NN,I,J),I=IB,IE)
  120   CONTINUE
  130 CONTINUE
      IF (NB.GT.0) THEN
        WRITE(6+CK_LUN,'(7H REC # ,I3,3X,I4,15H VALUES BLENDED)') NW,NB
      ENDIF
      GOTO 10
  140 STOP
      END


      subroutine JTrmax(ws,dp,rmax,spd)
      real ws,dp,rmax,ws30,slope
c*****************************************************************
cjjw - added to this version as per Norm's recommendation on 7-10-2002
c...  Computes radius to maximum wind as a function of pressure
c      difference & maximum sustained surface wind speed.
c      Calculation is based on relationships in Figure A1 of
c      Jelesnianski & Taylor (1973)

c        Complete reference:  Jelesnianski, C.P., and Taylor, A.D.,
c        "A Preliminary View of Storm Surges Before and After Storm
c        Modifications," NOAA Technical Memorandum ERL WMPO-3, 1973.
c*****************************************************************

c...  Compute max surface wind speed at rmax = 30 miles, in mph
c      (coefficients obtained by regression analysis of values
c      measured from Jelesnianski-Taylor (JT) curves)

      ws30 = 29.984 + 1.20035*dp - 0.002465642535*(dp**2)

c...  Compute slope of JT line of constant dp

      slope = 0.410604 + 0.00383248*dp

c...  Adjust best track max wind to be consistent with JT max wind
c       (developed for storm surge)
c       Subtract storm forward motion & convert from 1-min avg to 10-min avg.

      wsJT = (ws - spd*6080.0/5280.0)/1.05

c...  Compute rmax (in miles) for desired surface wind speed

      rmax = 30.0 - (wsJT - ws30)/slope

cvjp  write(20,100) dp,ws,wsJT,ws30,slope,rmax
  100 format(f5.0,3f6.1,f6.4,f6.1)

      return
      end
