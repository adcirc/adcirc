C
C  Version 1.2  vjp  12/07/99  All routines taken from ADCIRC
cjjw - fixes 053100
C

C***********************************************************************
C                                                                      *
C   THE FOLLOWING SUBROUTINES READ IN AND IN SOME CASES INTERPOLATE    *
C   ONTO THE ADCIRC GRID WIND AND PRESSURE FIELDS IN VARIOUS INPUT     *
C   FORMATS.                                                           *
C                                                                      *
C   ALL WIND SPEEDS ARE CONVERTED TO M/S AND ALL PRESSURES TO M OF H20 *
C   BEFORE THEY ARE RETURNED.                                          *
C                                                                      *
C***********************************************************************

C***********************************************************************
C                                                                      *
C   Convert time from year,month,day,hour,min,sec into seconds since   *
C   the beginning of the year.                                         *
C                                                                      *
C***********************************************************************

      SUBROUTINE TIMECONV(IYR,IMO,IDAY,IHR,IMIN,SEC,TIMESEC)
      USE PRESIZES, ONLY : SZ
      IMPLICIT NONE
      REAL(8) TIMESEC
C     jgf46.00 Explicitly declared the following variables
      INTEGER IDAY
      INTEGER IHR
      INTEGER IMIN
      REAL(SZ) SEC
      INTEGER IMO
      INTEGER ILEAP
      INTEGER IYR
C
      TIMESEC = (IDAY-1)*86400 + IHR*3600 + IMIN*60 + SEC
      IF(IMO.GE.2)  TIMESEC = TIMESEC + 31*86400
      ILEAP = (IYR/4)*4
      IF((ILEAP.EQ.IYR).AND.(IMO.GE.3)) TIMESEC = TIMESEC + 29*86400
      IF((ILEAP.NE.IYR).AND.(IMO.GE.3)) TIMESEC = TIMESEC + 28*86400
      IF(IMO.GE.4)  TIMESEC = TIMESEC + 31*86400
      IF(IMO.GE.5)  TIMESEC = TIMESEC + 30*86400
      IF(IMO.GE.6)  TIMESEC = TIMESEC + 31*86400
      IF(IMO.GE.7)  TIMESEC = TIMESEC + 30*86400
      IF(IMO.GE.8)  TIMESEC = TIMESEC + 31*86400
      IF(IMO.GE.9)  TIMESEC = TIMESEC + 31*86400
      IF(IMO.GE.10) TIMESEC = TIMESEC + 30*86400
      IF(IMO.GE.11) TIMESEC = TIMESEC + 31*86400
      IF(IMO.EQ.12) TIMESEC = TIMESEC + 30*86400
      IF(IMO.GT.12) THEN
        WRITE(6,*) 'FATAL ERROR IN SUBROUTINE TIMECONV - MONTH > 12 '
        WRITE(16,*) 'FATAL ERROR IN SUBROUTINE TIMECONV - MONTH > 12 '
        STOP
        ENDIF
      RETURN
      END

C***********************************************************************
C                                                                      *
C   READ IN AND INTERPOLATE ONTO THE ADCIRC GRID WIND FIELDS FROM U.S. *
C   NAVY FLEET NUMERIC WIND FILES.                                     *
C                                                                      *
C   NOTE: The ADCIRC grid information consists only of the Lon and Lat *
C   of the nodes.  THE LONS AND LATS MUST BE IN RADIANS!               *
C                                                                      *
C                                                                      *
C   MNWLAT = MAXIMUM NUMBER OF LATITUDES IN FLEET NUMERIC WIND FILE    *
C            SET = 1 IF FLEET NUMERIC WIND FILE NOT IN USE             *
C   MNWLON = MAXIMUM NUMBER OF LONGITUDES IN FLEET NUMERIC WIND FILE   *
C            SET = 1 IF FLEET NUMERIC WIND FILE NOT IN USE             *
C                                                                      *
C                           R.L. 4/17/96                               *
C***********************************************************************

      SUBROUTINE NWS3GET(X,Y,SLAM,SFEA,WVNX,WVNY,IWTIME,IWYR,WTIMED,NP,
     &    NWLON,NWLAT,WLATMAX,WLONMIN,WLATINC,WLONINC,ICS,ISTAT)
       USE PRESIZES
C
      IMPLICIT NONE
      INTEGER ISTAT
      REAL(8) WTIMED,RAD2DEG,PI,DEG2RAD
C      DIMENSION X(*),Y(*),SLAM(*),SFEA(*),WVNX(*),WVNY(*) !jgf46.00 comm. out
      DIMENSION WVNX(*),WVNY(*)
C RAWMET
      REAL(SZ),ALLOCATABLE :: WVXFN(:,:),WVYFN(:,:),PRN(:,:)
C     jgf46.00 Explicitly declared the following variables
      INTEGER IWTIME
      INTEGER IWYR
      INTEGER IWMO
      INTEGER IWDAY
      INTEGER IWHR
      INTEGER I
      INTEGER NWLAT
      INTEGER J
      INTEGER NWLON
      INTEGER NP
      INTEGER ICS
      INTEGER LATIND1
      INTEGER LATIND2
      INTEGER LONIND1
      INTEGER LONIND2
      REAL(SZ) WSPEED
      REAL(SZ) WDIR
      REAL(8) YCOOR
      REAL(8) XCOOR
      REAL(8), DIMENSION(NP) :: SFEA
      REAL(8), DIMENSION(NP) :: SLAM
      REAL(8), DIMENSION(NP) :: X
      REAL(8), DIMENSION(NP) :: Y
      REAL(8) WLATMAX
      REAL(8) WLATINC
      REAL(8) WLONMIN
      REAL(8) WLONINC
      REAL(8) WLONM
      REAL(8) WLATM
      REAL(8) XWRATIO
      REAL(8) YWRATIO
      REAL(8) WVNX
      REAL(8) WVNY
      INTRINSIC :: SIN, COS
C
      PI=3.141592653589793D0
      RAD2DEG = 180.D0/PI
      DEG2RAD = PI/180.D0
C
      ALLOCATE ( WVXFN(MNWLAT,MNWLON),WVYFN(MNWLAT,MNWLON),
     &           PRN(MNWLAT,MNWLON) )
C
      ISTAT = 0
      READ(22,*,END=9999) IWTIME
      IWYR = IWTIME/1000000
      IWMO = IWTIME/10000 - IWYR*100
      IWDAY = IWTIME/100 - IWYR*10000 - IWMO*100
      IWHR = IWTIME - IWYR*1000000 - IWMO*10000 - IWDAY*100
      CALL TIMECONV(IWYR,IWMO,IWDAY,IWHR,0,0.0d0,WTIMED)

      DO I=1,NWLAT
        READ(22,*) (WVXFN(I,J),J=1,NWLON)
      END DO
      DO I=1,NWLAT
        READ(22,*) (WVYFN(I,J),J=1,NWLON)
      END DO

      DO I=1,NWLAT                             !CONVERT TO X AND Y COMPONENTS
        DO J=1,NWLON
          WSPEED=WVXFN(I,J)
          WDIR=WVYFN(I,J)*DEG2RAD
          WVXFN(I,J)=-WSPEED*SIN(WDIR)
          WVYFN(I,J)=-WSPEED*COS(WDIR)
          END DO
        END DO

      DO I=1,NP                                !INTERPOLATE TO ADCIRC GRID
        IF(ICS.EQ.2) THEN
cvjp      YCOOR=SFEA(I)*RAD2DEG
cvjp      XCOOR=SLAM(I)*RAD2DEG
          YCOOR=SFEA(I)
          XCOOR=SLAM(I)
        ENDIF
        IF(ICS.EQ.1) THEN
          YCOOR=Y(I)
          XCOOR=X(I)
        ENDIF
        LATIND2=(WLATMAX-YCOOR)/WLATINC + 1
        IF(LATIND2.EQ.NWLAT) LATIND2=LATIND2-1
        LATIND1=LATIND2 + 1
        LONIND1=(XCOOR-WLONMIN)/WLONINC + 1
        IF(LONIND1.EQ.NWLON) LONIND1=LONIND1-1
        LONIND2=LONIND1+1
        WLONM = WLONMIN + (LONIND1-1)*WLONINC
        WLATM = WLATMAX - (LATIND1-1)*WLATINC
        XWRATIO=(XCOOR-WLONM)/WLONINC
        YWRATIO=(YCOOR-WLATM)/WLATINC
c       print *, LATIND1,LATIND2,LONIND1,LONIND2,XCOOR,YCOOR

        WVNX(I) = WVXFN(LATIND2,LONIND2)*XWRATIO*YWRATIO
     &           + WVXFN(LATIND2,LONIND1)*(1.-XWRATIO)*YWRATIO
     &           + WVXFN(LATIND1,LONIND2)*XWRATIO*(1.-YWRATIO)
     &           + WVXFN(LATIND1,LONIND1)*(1.-XWRATIO)*(1.-YWRATIO)
        WVNY(I) = WVYFN(LATIND2,LONIND2)*XWRATIO*YWRATIO
     &           + WVYFN(LATIND2,LONIND1)*(1.-XWRATIO)*YWRATIO
     &           + WVYFN(LATIND1,LONIND2)*XWRATIO*(1.-YWRATIO)
     &           + WVYFN(LATIND1,LONIND1)*(1.-XWRATIO)*(1.-YWRATIO)
c       print *, I,WVNX(I),WVNY(I)
      END DO
C
      GOTO 99
C
9999  CONTINUE
      ISTAT = 1
C
  99  RETURN
      END

C***********************************************************************
C                                                                      *
C   Read onto the ADCIRC grid wind fields from the PBL-JAG model       *
C                                                                      *
C   Output from this subroutine is U,V (M/S) and P (M H20) on the      *
C   ADCIRC grid.                                                       *
C                                                                      *
C   The background pressure is assumed to be 1013 Mbars                *
C                                                                      *
C                           R.L.11/06/96                               *
C***********************************************************************

cjjw - re-ordered the subroutine call
      SUBROUTINE NWS4GET(WVNX,WVNY,PRN,G,RHOWAT,NP,DONE)
      USE PRESIZES
      IMPLICIT NONE
      INTEGER NP, I, NHG
cjjw - fixed info here by adding rhowatg
      REAL(SZ) WVNX(*),WVNY(*),PRN(*),RHOWAT,RHOWATG
      REAL(8) G
      CHARACTER*80 PBLJAGF
      LOGICAL DONE

      DONE = .FALSE.
      RHOWATG=RHOWAT*G

cjjw - added following trace
cjjw      write(69,9091) g,rhowat,rhowatg,np
cjjw9091  format(/,' trace begin',50("+"),/,
cjjw     &  ' g = ' ,e30.20,/,
cjjw     &  ' rhowat = ' ,e30.20,/,
cjjw     &  ' rhowatg = ' ,e30.20,/,
cjjw     &  ' np = ' ,i30,/,
cjjw     &  ' trace end ',50("+"),/)
C
      DO I=1,NP
        WVNX(I)=0.d0
        WVNY(I)=0.d0
        PRN(I)=101300.d0/RHOWATG
      END DO
C
  170 READ(22,'(A80)',END=9999) PBLJAGF
      IF(PBLJAGF(2:2).EQ.'#') GOTO 170
  171 READ(PBLJAGF,'(I8,3E13.5)',END=9999) NHG,WVNX(NHG),WVNY(NHG),
     &                                     PRN(NHG)

C     jgf46.02 From now on, wind files must contain data which are
C     appropriate for the time span specified in the file: commented out
C     the following two lines.
C     WVNX(NHG)=WVNX(NHG)*1.04d0*0.5144d0  !CONVERT 30-MIN WINDS IN
C     WVNY(NHG)=WVNY(NHG)*1.04d0*0.5144d0  !KNOTS TO 10-MIN WIND IN M/S
C     jgf46.02 Added the following two lines.
      WVNX(NHG)=WVNX(NHG)*0.5144d0         !CONVERT KNOTS TO M/S
      WVNY(NHG)=WVNY(NHG)*0.5144d0
      PRN(NHG)=100.d0*PRN(NHG)/RHOWATG   !CONVERT MILLIBARS TO M OF WATER
      READ(22,'(A80)',END=9999) PBLJAGF
C
      IF(PBLJAGF(2:2).NE.'#') THEN
         GOTO 171
      ELSE
         GO TO 8888
      ENDIF
9999  DONE = .TRUE.
8888  RETURN
      END


C***********************************************************************
C                                                                      *
C   Read in and interpolate onto the ADCIRC grid wind fields from U.S. *
C   National Weather Service AVN model SFLUX meteorological files.     *
C                                                                      *
C   The input files are in binary and have been created by the GRIB    *
C   unpacking program unpkgrb1.f to extract only the U 10M, V 10M, and *
C   surface P fields.                                                  *
C                                                                      *
C   The SFLUX files utilize a global Gaussian Lon/Lat grid which is    *
C   constructed in these subroutines.                                  *
C                                                                      *
C   NOTE: The ADCIRC grid information consists only of the Lon and Lat *
C   of the nodes.  THE LONS AND LATS MUST BE IN RADIANS!               *
C                                                                      *
C   Output from this subroutine is U,V (M/S) and P (M H20) on the      *
C   ADCIRC grid.                                                       *
C                                                                      *
C   MNWLAT = LATB = 190    FOR GAUSSIAN GRID                           *
C   MNWLON = LONB = 384    FOR GAUSSIAN GRID                           *
C                                                                      *
C                           R.L. 4/18/96                               *
C***********************************************************************

      SUBROUTINE NWS10GET(NWSGGWI,FLON,FLAT,ULL,VLL,PLL,NP,RHOWAT,G,
     &                    LONB,LATB,GBL,FOUND)
      USE PRESIZES
C
      IMPLICIT NONE
      INTEGER,ALLOCATABLE :: N00(:),N10(:),N11(:),N01(:)
      REAL(8),ALLOCATABLE :: D00(:),D10(:),D11(:),D01(:) !jgf45.07 distances?
      REAL(SZ),ALLOCATABLE :: UG(:),VG(:),PG(:)
      REAL(SZ),ALLOCATABLE :: COLRAB(:),DUMMY(:),GCLAT(:),GCLON(:)
C      DIMENSION FLAT(*),FLON(*)      !jgf46.00 comm. out
C      DIMENSION ULL(*),VLL(*),PLL(*) !jgf46.00 comm. out
      INTEGER KGDS(200)
      CHARACTER*1 PDS(50),FNAME2(8)
      CHARACTER*8 FNAME1
C
Cvjp Return the Filename of the Global Wind Stress file
C
      CHARACTER*8 GBL
      EQUIVALENCE (FNAME1,FNAME2)
C
Cvjp  Modified Logic to return if file not found
C
      LOGICAL FOUND
C     jgf46.00 Explicitly declared the following variables:
      REAL(SZ), DIMENSION(*) :: PLL
      REAL(SZ), DIMENSION(*) :: ULL
      REAL(SZ), DIMENSION(*) :: VLL
      REAL(SZ) RHOWAT
      REAL(SZ) G
      INTEGER NWSGGWI
      INTEGER J
      INTEGER JJ
      INTEGER LATB
      INTEGER LONB
      INTEGER NP
      INTEGER IEXT
      INTEGER IDIG1
      INTEGER IDIG2
      INTEGER IDIG3
      INTEGER kerr
      INTEGER LENPDS
      INTEGER LENKGDS
      INTEGER NWORDS
      INTEGER N
      REAL(SZ) GDLON
      REAL(8), PARAMETER  ::  PI=3.141592653589793D0 !jgf46.00 was implicit??
      REAL(8), PARAMETER  ::  TWOPI=2.0D0*PI         !jgf46.00 was implicit??
      REAL(SZ), DIMENSION(*) :: FLAT
      REAL(SZ), DIMENSION(*) :: FLON
      REAL(SZ) P1
      REAL(SZ) P2
      REAL(SZ) P3
      REAL(SZ) P4
      REAL(SZ) U1
      REAL(SZ) U2
      REAL(SZ) U3
      REAL(SZ) U4
      REAL(SZ) V1
      REAL(SZ) V2
      REAL(SZ) V3
      REAL(SZ) V4
      REAL(SZ) RHOWATG
C
      RHOWATG=RHOWAT*G !jgf46.00 this implicitly declared and not initialized??
C
C INTERP10
      ALLOCATE ( N00(MNWP),N10(MNWP),N11(MNWP),N01(MNWP) )
      ALLOCATE ( D00(MNWP),D10(MNWP),D11(MNWP),D01(MNWP) )
C RAWMET
      ALLOCATE ( UG(MNWLAT*MNWLON),VG(MNWLAT*MNWLON),
     &  PG(MNWLAT*MNWLON) )
C
      ALLOCATE ( COLRAB(MNWLAT),DUMMY(MNWLAT),GCLAT(MNWLAT),
     &  GCLON(MNWLON) )

C
C...The first time the subroutine is called, setup the Gaussian grid and
C...determine the interpolating factors for the ADCIRC grid.
C
      IF (NWSGGWI.EQ.-1) THEN
        CALL GLATS(LATB/2,COLRAB,DUMMY,DUMMY,DUMMY)
        DO J=1,LATB/2
          GCLAT(J)=COLRAB(J)
          JJ=LATB-J+1
          GCLAT(JJ)=PI-COLRAB(J)
        ENDDO
        GDLON=TWOPI/LONB
        DO J=1,LONB
          GCLON(J)=GDLON*(J-1)
        END DO
        CALL G2RINI(GCLON,GCLAT,FLON,FLAT,N00,N10,N11,N01,D00,D10,D11,
     &              D01,NP,LONB,LATB)
        RETURN
       ENDIF

C...Figure out the GRIB data file name

      FNAME1='fort.   '
      IEXT=200 + NWSGGWI*6
      IDIG1=IEXT/100
      IDIG2=(IEXT-100*IDIG1)/10
      IDIG3=(IEXT-100*IDIG1-10*IDIG2)
      FNAME2(6)=CHAR(IDIG1+48)
      FNAME2(7)=CHAR(IDIG2+48)
      FNAME2(8)=CHAR(IDIG3+48)


C...Enter, locate and open the GRIB format data file

 1010 FORMAT(' File ',A8,' WAS NOT FOUND!  FATAL ERROR')
 1011 FORMAT(' File ',A8,' WAS FOUND!  Opening & Processing file')

      INQUIRE(FILE=FNAME1,EXIST=FOUND)
      IF(FOUND) GOTO 32
      GBL(1:8) = FNAME1
C
C--if file not found return logical FOUND to caller
C
      GBL(1:8) = FNAME1
      RETURN
C
  32  WRITE(*,1011) FNAME1
      GBL(1:8) = FNAME1
      OPEN(IEXT,FILE=FNAME1,status='old',access='sequential',
     &     form='unformatted',iostat=kerr)

C...Read the GRIB data file

      READ(IEXT,END=1100) LENPDS,LENKGDS,NWORDS
      IF(LENPDS.GT.0) READ(IEXT,END=1100) (pds(j),j=1,lenpds)
      IF(LENKGDS.GT.0) READ(IEXT,END=1100) (kgds(j),j=1,lenkgds)
      IF(NWORDS.GT.0) READ(IEXT,END=1100) (UG(J),J=1,NWORDS)

      READ(IEXT,END=1100) LENPDS,LENKGDS,NWORDS
      IF(LENPDS.GT.0) READ(IEXT,END=1100) (pds(j),j=1,lenpds)
      IF(LENKGDS.GT.0) READ(IEXT,END=1100) (kgds(j),j=1,lenkgds)
      IF(NWORDS.GT.0) READ(IEXT,END=1100) (VG(J),J=1,NWORDS)

      READ(IEXT,END=1100) LENPDS,LENKGDS,NWORDS
      IF(LENPDS.GT.0) READ(IEXT,END=1100) (pds(j),j=1,lenpds)
      IF(LENKGDS.GT.0) READ(IEXT,END=1100) (kgds(j),j=1,lenkgds)
      IF(NWORDS.GT.0) READ(IEXT,END=1100) (PG(J),J=1,NWORDS)

 1100 CLOSE(IEXT)


C.....Do from the Gaussian grid to the ADCIRC grid
C.....Convert pressure from N/M^2 to M of H20

      DO N=1,NP
        P1=PG(N00(N))
        P2=PG(N10(N))
        P3=PG(N11(N))
        P4=PG(N01(N))
        U1=UG(N00(N))
        U2=UG(N10(N))
        U3=UG(N11(N))
        U4=UG(N01(N))
        V1=VG(N00(N))
        V2=VG(N10(N))
        V3=VG(N11(N))
        V4=VG(N01(N))
        PLL(N)=P1*D00(N)+P2*D10(N)+P3*D11(N)+P4*D01(N)
        ULL(N)=U1*D00(N)+U2*D10(N)+U3*D11(N)+U4*D01(N)
        VLL(N)=V1*D00(N)+V2*D10(N)+V3*D11(N)+V4*D01(N)
        PLL(N)=PLL(N)/RHOWATG
        END DO

      RETURN
      end


C***********************************************************************
C  Subroutine to compute the latutudes in a Global Gaussian Lat/Lon    *
C  grid with T126 resolution (GRIB Grid type 126).                     *
C                                                                      *
C       modified from the original GLATS by R.L. 4/24/96               *
C***********************************************************************

      SUBROUTINE GLATS(LGGHAF,COLRAD,WGT,WGTCS,RCS2)
CCCC      %INCLUDE DBGLATS ;
CCCC  HALF PRECISION COLRAD,WGT,WGTCS,RCS2
C     DOUBLE PRECISION COLRAD,WGT,WGTCS,RCS2
      USE PRESIZES, ONLY : SZ
      IMPLICIT NONE

C      DIMENSION COLRAD(*),WGT(*),WGTCS(*),RCS2(*) !jgf46.00 comm. out
C     jgf46.00 Explicitly declare the following variables.
      REAL(SZ) EPS
      REAL(SZ) SI
      REAL(SZ) RL2
      REAL(SZ) SCALEVAL
      REAL(SZ) RAD
      REAL(SZ) PHI
      REAL(SZ) P1
      REAL(SZ) P2
      REAL(SZ) SN
      REAL(SZ) RC
      REAL(8) PI
      REAL(8) DRAD
      REAL(8) DRADZ
      REAL(8) X
      REAL(8) W
      REAL(SZ), DIMENSION(*) :: RCS2
      REAL(SZ), DIMENSION(*) :: COLRAD
      REAL(8), DIMENSION(*) :: WGT
      REAL(8), DIMENSION(*) :: WGTCS
      INTEGER K1
      INTEGER L2
      INTEGER LGGHAF
      INTEGER K
      INTEGER ITER
C
      EPS=1.E-6
C     EPS=1.E-12
C     PRINT 101
C101  FORMAT ('0 I   COLAT   COLRAD     WGT', 12X, 'WGTCS',
CCCC 1 10X, 'ITER  RES')
      SI = 1.0
      L2=2*LGGHAF
      RL2=L2
      SCALEVAL = 2.0/(RL2*RL2)
      K1=L2-1
      PI = ACOS(-1.)
      DRADZ = PI / 360.
      RAD = 0.0
      DO 1000 K=1,LGGHAF
      ITER=0
      DRAD=DRADZ
1     CALL POLY(L2,RAD,P2)
2     P1 =P2
      ITER=ITER+1
      RAD=RAD+DRAD
      CALL POLY(L2,RAD,P2)
      IF(SIGN(SI,P1).EQ.SIGN(SI,P2)) GO TO 2
      IF(DRAD.LT.EPS)GO TO 3
      RAD=RAD-DRAD
      DRAD = DRAD * 0.25
      GO TO 1
3     CONTINUE
      COLRAD(K)=RAD
      PHI = RAD * 180 / PI
      CALL POLY(K1,RAD,P1)
      X = COS(RAD)
      W = SCALEVAL * (1.0 - X*X)/ (P1*P1)
      WGT(K) = W
      SN = SIN(RAD)
      W=W/(SN*SN)
      WGTCS(K) = W
      RC=1./(SN*SN)
      RCS2(K) = RC
      CALL POLY(L2,RAD,P1)
C     PRINT 102,K,PHI,COLRAD(K),WGT(K),WGTCS(K),ITER,P1
C102  FORMAT(1H ,I2,2X,F6.2,2X,F10.7,2X,E13.7,2X,E13.7,2X,I4,2X,D13.7)
1000  CONTINUE
c     PRINT 100,LGGHAF
c100  FORMAT(1H ,'SHALOM FROM 0.0 GLATS FOR ',I3)
      RETURN
      END


C***********************************************************************
C  Subroutine used by GLATS.                                           *
C***********************************************************************

CFPP$ EXPAND(POLY)
      SUBROUTINE POLY(N,RAD,P)
CCCC      %INCLUDE DBPOLY ;
      IMPLICIT NONE
C     jgf46.00 Explicitly declared the following variables.
      REAL(8) X
      REAL(8) RAD
      REAL(8) Y1
      REAL(8) Y2
      REAL(8) Y3
      REAL(8) G
      REAL(8) P
      INTEGER I, N
      INTRINSIC :: COS, FLOAT
C
      X = COS(RAD)
      Y1 = 1.0
      Y2=X
      DO 1 I=2,N
      G=X*Y2
      Y3=G-Y1+G-(G-Y1)/FLOAT(I)
      Y1=Y2
      Y2=Y3
1     CONTINUE
      P=Y3
      RETURN
      END


C***********************************************************************
C  Subroutine to compute the factors to interpolate from a global      *
C  Gaussian Lat/Lon grid with T126 resolution (GRIB Grid type 126)     *
C  onto another grid.                                                  *
C                                                                      *
C  The new grid is a series of longitude and latitude points contained *
C  in the FLON and FLAT arrays with a total number of points NP        *
C                                                                      *
C       modified from the original G2RINI by R.L. 4/17/96              *
C***********************************************************************

      SUBROUTINE G2RINI(GCLON,GCLAT,FLON,FLAT,N00,N10,N11,N01,D00,D10,
     &                                          D11,D01,NP,LONB,LATB)
      USE PRESIZES, ONLY : SZ
      IMPLICIT NONE
C     jgf46.00 commented next three lines out
C      DIMENSION GCLAT(*),GCLON(*)
C      DIMENSION FLAT(*),FLON(*)
C      DIMENSION N00(*),N10(*),N11(*),N01(*),D00(*),D10(*),D11(*),D01(*)

      INTEGER, SAVE :: ICALL
      DATA ICALL/0/
C     jgf46.00 Explicitly declared the following variables.
      REAL(SZ), DIMENSION(*) :: GCLAT
      REAL(SZ), DIMENSION(*) :: GCLON
      REAL(SZ), DIMENSION(*) :: FLAT
      REAL(SZ), DIMENSION(*) :: FLON
      REAL(SZ), DIMENSION(*) :: D00
      REAL(SZ), DIMENSION(*) :: D10
      REAL(SZ), DIMENSION(*) :: D11
      REAL(SZ), DIMENSION(*) :: D01
      INTEGER, DIMENSION(*) :: N00
      INTEGER, DIMENSION(*) :: N10
      INTEGER, DIMENSION(*) :: N11
      INTEGER, DIMENSION(*) :: N01
      INTEGER NLAT
      INTEGER LATB
      INTEGER NLON
      INTEGER LONB
      INTEGER N
      INTEGER I
      INTEGER NP
      INTEGER LON
      INTEGER LONP1
      INTEGER LAT
      INTEGER LATP1
      REAL(SZ) DLAT
      REAL(SZ) DLON
      REAL(SZ) FLONWORK
      REAL(SZ) COLAT
      REAL(SZ) DDLAT
      REAL(SZ) XLAT
      REAL(SZ) DFLAT1
      REAL(SZ) DFLAT
      REAL(SZ) DDLON
      REAL(SZ) XLON
      REAL(SZ) DFLON
      REAL(SZ) DFLON1
      REAL(8) PI
      REAL(8) HFPI
      REAL(8) TWOPI
C
      IF( ICALL .EQ. 0 ) THEN
        ICALL = 1
c       PRINT 1234
c1234   FORMAT(' = IN ROUTINE G2RINI FOR HORIZONTAL INTERPOLATION = ')
        PI=ACOS(-1.)
        HFPI=PI/2.0
        TWOPI=PI*2.0

C...Compute estimated DLAT, true DLON for Gaussian grid

        NLAT=LATB
        NLON=LONB
        DLAT=PI/FLOAT(NLAT-1)
        DLON=TWOPI/FLOAT(NLON)
        N=0

C...Loop through all the nodes in the grid to be interpolated onto and
C.....compute the interpolating factors.

        DO I=1,NP

C.....Compute initial guess of which lon value FLON(I) is in the Gaussian file
C.......Check that this value is reasonable.

          FLONWORK=FLON(I)
          IF(FLONWORK.LT.0.) FLONWORK=FLONWORK+TWOPI
          LON=FLONWORK/DLON + 1
          LONP1=LON+1
          IF(LON.EQ.NLON) LONP1=1        !Circle condition
          IF((LON.LT.1).OR.(LON.GT.NLON)) THEN
            PRINT *,' ***** ERROR IN LON ****'
            PRINT *,' I ',I
            PRINT *,' LON ',LON
            PRINT *,' DLON ',DLON
            PRINT *,' FLON ',FLON(I)
            STOP
            endif

C.....Compute initial guess of which lat value FLAT(I) is in the Gaussian file
C.......Check that this value is reasonable.

          COLAT=HFPI-FLAT(I)
          LAT=COLAT/DLAT + 1
          IF(LAT.EQ.NLAT) LAT=LAT-1
          LATP1=LAT+1
          IF((LAT.LT.1).OR.(LAT.GT.NLAT)) THEN
            PRINT *,' ***** ERROR IN LAT ****'
            PRINT *,' I ',I
            PRINT *,' LAT ',LAT
            PRINT *,' DLAT ',DLAT
            PRINT *,' FLAT ',FLAT(I)
            STOP
          ENDIF

 5        CONTINUE
          IF((COLAT.GE.GCLAT(LAT)).AND.(COLAT.LE.GCLAT(LATP1))) GO TO 9
          IF(COLAT.LT.GCLAT(LAT)) THEN
            LATP1=LAT
            LAT=LAT-1
            IF(LAT.LE.0) THEN
              LAT=1
              LATP1=2
              GOTO 9
              ENDIF
            GOTO 5
            ENDIF
          IF(COLAT.GT.GCLAT(LATP1)) THEN
            LAT=LAT+1
            LATP1=LAT+1
            IF(LAT.GE.NLAT ) THEN
              LAT=NLAT-1
              LATP1=NLAT
              GOTO 9
              ENDIF
            GOTO 5
            ENDIF

 9        CONTINUE
          DDLAT=GCLAT(LATP1)-GCLAT(LAT)
          XLAT=GCLAT(LAT)
          DFLAT1=(COLAT-XLAT)/DDLAT
          IF(LAT.EQ.1) DFLAT1=MAX(0.0D0,DFLAT1)         !MODIFY THIS FOR POLAR POINTS
          IF(LATP1.EQ.NLAT) DFLAT1=MIN(1.0D0,DFLAT1)    !MODIFY THIS FOR POLAR POINTS
          DFLAT=1.-DFLAT1
          DDLON=DLON
          XLON=GCLON(LON)
          DFLON1=(FLONWORK-XLON)/DDLON
          DFLON=1.-DFLON1
          N=N+1
          D00(N)=DFLON*DFLAT
          D10(N)=DFLON1*DFLAT
          D11(N)=DFLON1*DFLAT1
          D01(N)=DFLON*DFLAT1
          N00(N)=LON+(LAT-1)*NLON
          N10(N)=LONP1+(LAT-1)*NLON
          N11(N)=LONP1+(LATP1-1)*NLON
          N01(N)=LON+(LATP1-1)*NLON

          END DO

c       WRITE(*,*) ' D00 TO D11 SHOULD BE ALL POSITIVE.'

        ELSE
c       WRITE(*,*) ' G2RINI ALREADY CALLED '
        ENDIF

      RETURN
      END

C***********************************************************************
C                                                                      *
C   Read in and interpolate onto the ADCIRC grid wind fields from U.S. *
C   National Weather Service ETA-29 model that have been stripped down *
C   and given to us by NOAA.                                           *
C                                                                      *
C   The input files are in binary and have been created by NOAA and    *
C   contain only the U 10M, V 10M, (M/S) and surface P fields (mbars). *
C                                                                      *
C   The ETA-29 model uses an E grid and therefore the U and V          *
C   components are not oriented along lines of constant latitute and   *
C   longitude. These must be converted to be useful in ADCIRC.         *
C                                                                      *
C   NOTE: The ADCIRC grid information consists only of the Lon and Lat *
C   of the nodes.  THE LONS AND LATS MUST BE IN RADIANS!               *
C                                                                      *
C   Output from this subroutine is U,V (M/S) and P (M H20) on the      *
C   ADCIRC grid.                                                       *
C                                                                      *
C   MNWLAT = LATB = 271    FOR ETA-29 GRID                             *
C   MNWLON = LONB = 181    FOR ETA-29 GRID                             *
C                                                                      *
C                           R.L. 1/11/97                               *
C***********************************************************************

      SUBROUTINE NWS11GET(NWSEGWI,IDSETFLG,FLON,FLAT,ULL,VLL,PLL,NP,
     &   NWLON,NWLAT,RHOWAT,G,GBL,FOUND)
      USE PRESIZES
C
      IMPLICIT NONE
      INTEGER NWLON,NWLAT,I,N
      INTEGER,ALLOCATABLE :: N1(:),N2(:),N3(:)
      REAL(SZ),ALLOCATABLE    :: D1(:),D2(:),D3(:),BETAU(:)
      REAL(SZ),ALLOCATABLE    :: UE(:),VE(:),PE(:)
C      DIMENSION FLAT(*),FLON(*)      !jgf46.00 commented out
C      DIMENSION ULL(*),VLL(*),PLL(*) !jgf46.00 commented out
      CHARACTER*1 FNAME2(8)
      CHARACTER*8 FNAME1,GBL
      EQUIVALENCE (FNAME1,FNAME2)
      LOGICAL FOUND
C     jgf46.00 Explicitly declared the following variables.
      REAL(SZ), DIMENSION(*) :: FLAT
      REAL(SZ), DIMENSION(*) :: FLON
      REAL(SZ), DIMENSION(*) :: ULL
      REAL(SZ), DIMENSION(*) :: VLL
      REAL(SZ), DIMENSION(*) :: PLL
      REAL(8) PI
      REAL(8) HFPI
      REAL(8) TWOPI
      REAL(8) DEG2RAD
      REAL(SZ) RHOWATG100
      REAL(SZ) RHOWAT
      REAL(SZ) G
      INTEGER NWSEGWI
      INTEGER IDSETFLG
      INTEGER NP
      REAL(SZ) FLONDEG
      REAL(SZ) FLATDEG
      INTEGER IEXT
      INTEGER IDIG1
      INTEGER IDIG2
      INTEGER IDIG3
      INTEGER kerr
      INTEGER IYEAR
      INTEGER IMONTH
      INTEGER IDAY
      INTEGER IHOUR
      REAL(SZ) P1
      REAL(SZ) P2
      REAL(SZ) P3
      REAL(SZ) U1
      REAL(SZ) U2
      REAL(SZ) U3
      REAL(SZ) DFLON1
      REAL(SZ) V1
      REAL(SZ) V2
      REAL(SZ) V3
      REAL(SZ) UE29
      REAL(SZ) VE29
      REAL(SZ) CBETAU
      REAL(SZ) SBETAU
C
C  INTERP11
      ALLOCATE ( N1(MNWP),N2(MNWP),N3(MNWP) )
      ALLOCATE ( D1(MNWP),D2(MNWP),D3(MNWP),BETAU(MNWP) )
C  RAWMET
      ALLOCATE ( UE(MNWLON*MNWLAT),VE(MNWLON*MNWLAT),PE(MNWLON*MNWLAT))

      PI=ACOS(-1.)
      TWOPI=PI*2.D0
      HFPI=PI/2.D0
      DEG2RAD=PI/180.D0
      RHOWATG100=RHOWAT*G*100.

C...The first time the subroutine is called, setup the interpolating factors
C...between the Eta-29 grid and the ADCIRC grid.

      IF((NWSEGWI.EQ.0).AND.(IDSETFLG.EQ.0)) THEN
        WRITE(*,*) '  '
        WRITE(*,*) 'Computing ETA29 met field interpolating factors'
C
        DO I=1,NP
          flondeg=flon(i)/deg2rad
          flatdeg=flat(i)/deg2rad
          CALL E29SEARCH(I,FLONDEG,FLATDEG,N1(I),N2(I),N3(I),
     &                                     D1(I),D2(I),D3(I),betau(i))
          END DO
        return
        ENDIF

C...Figure out the met data file name

      FNAME1='fort.   '
      IEXT=200 + NWSEGWI
      IDIG1=IEXT/100
      IDIG2=(IEXT-100*IDIG1)/10
      IDIG3=(IEXT-100*IDIG1-10*IDIG2)
      FNAME2(6)=CHAR(IDIG1+48)
      FNAME2(7)=CHAR(IDIG2+48)
      FNAME2(8)=CHAR(IDIG3+48)

C...If appropriate, enter, locate and open the met data file

 1010 FORMAT(' File ',A8,' WAS NOT FOUND!  FATAL ERROR')
 1011 FORMAT(' File ',A8,' WAS FOUND!  Opening & Processing file')

      INQUIRE(FILE=FNAME1,EXIST=FOUND)
      IF(FOUND) GOTO 32
C
C--if not found return logical FOUND to caller
C
      GBL(1:8) = FNAME1
      RETURN
c
  32  CONTINUE
      GBL(1:8) = FNAME1
      IF ((NWSEGWI.EQ.0).OR.(IDSETFLG.EQ.1)) THEN
        WRITE(*,1011) FNAME1
        OPEN(IEXT,FILE=FNAME1,status='old',access='sequential',
     &   form='unformatted',iostat=kerr)
      ENDIF

C...Read the met data file

      READ(IEXT,END=1100) IYEAR,IMONTH,IDAY,IHOUR
      READ(IEXT,END=1100) (UE(N),N=1,NWLON*NWLAT),
     &                    (VE(N),N=1,NWLON*NWLAT),
     &                    (PE(N),N=1,NWLON*NWLAT)

      IF(NWSEGWI.EQ.0) THEN  !If the first file, read until the end
        DO I=2,IDSETFLG
          READ(IEXT,END=1100) IYEAR,IMONTH,IDAY,IHOUR
          READ(IEXT,END=1100) (UE(N),N=1,NWLON*NWLAT),
     &                        (VE(N),N=1,NWLON*NWLAT),
     &                        (PE(N),N=1,NWLON*NWLAT)
          ENDDO
        ENDIF

 1100 IF(IDSETFLG.EQ.8) CLOSE(IEXT)

C.....Interpolate onto ADCIRC grid
C.....Convert velocity from the E grid reference to a lat/lon reference
C.....Convert pressure from millibars to N/M^2 to M of H20

      DO N=1,NP
        P1=PE(N1(N))
        P2=PE(N2(N))
        P3=PE(N3(N))
        U1=UE(N1(N))
        U2=UE(N2(N))
        U3=UE(N3(N))
        V1=VE(N1(N))
        V2=VE(N2(N))
        V3=VE(N3(N))
        UE29=U1*D1(N)+U2*D2(N)+U3*D3(N)
        VE29=V1*D1(N)+V2*D2(N)+V3*D3(N)
        CBETAU=COS(BETAU(N))
        SBETAU=SIN(BETAU(N))
        ULL(N)=UE29*CBETAU - VE29*SBETAU
        VLL(N)=UE29*SBETAU + VE29*CBETAU
        PLL(N)=P1*D1(N)+P2*D2(N)+P3*D3(N)
        PLL(N)=PLL(N)/RHOWATG100
        END DO

      RETURN
      end

C***********************************************************************
C  Subroutine to find where a given lon,lat falls in the Eta29 grid,   *
C     determine the interpolating factors to interpolate Eta29 fields  *
C     to that position, and finally to compute the angle to rotate the *
C     Eta29 velocity field to get to a lon, lat coordinated system.    *
C                                                                      *
C                    Written by R.L.       1/12/98                     *
C***********************************************************************

      subroutine e29search(node,FLON,FLAT,NN1,NN2,NN3,DD1,DD2,DD3,betau)
      use presizes, only : sz
      implicit none
      real(8) lamda0,phi0,rphi0,cphi0,sphi0,tphi0,dlamda,dphi,rdlamda,
     &       rdphi,rflat,tflat,sflat,cflat,a,rlamar,cphiicrlamda,phiarg,
     &       rphii,rlamda,ri1,ri2,rj,dgtora
      real(8) lamda,lamdaa,lamdab,lamdac,lamdad,lamdae,lamdag
      real(8) phi,phia,phib,phic,phid,phie,phig
C     jgf46.00 Explicitly declared the following variables.
      integer node, nn1, nn2, nn3
      real(sz) betau, dd1, dd2, dd3, flat, flon, ri
      real(8) X1, X2, X3, X4, Y1, Y2, Y3, Y4
      real(8) AEMIN, AREAS, A1, A2, A3, AA, AE
      integer icode, nwlon, nwlat, i, j
      integer jm2, im2, n, ia, ja, na, ib, jb, nb, ic, jc, nc
      integer ie, id, jd, nd, je, ne, ig, jg, ng, ifflag
      intrinsic :: mod

      icode=0

      nwlon=181
      nwlat=271
      dgtora=3.1415926d0/180.d0
      lamda0=-97.0d0
      phi0=41.0d0
      rphi0=dgtora*phi0
      cphi0=cos(rphi0)
      sphi0=sin(rphi0)
      tphi0=tan(rphi0)
      dlamda=7./36.d0
      dphi=5/27.d0
      rdlamda=dgtora*dlamda
      rdphi=dgtora*dphi

      rflat=flat*dgtora
        tflat=tan(rflat)
      sflat=sin(rflat)
      cflat=cos(rflat)

c     compute the position of the closest node in the E29 grid

        a=flon-lamda0
        rlamar=cos(a*dgtora)
        cphiicrlamda=(rlamar+tflat*tphi0)*cflat*cphi0
        phiarg=sflat
      rphii=asin((phiarg-sphi0*cphiicrlamda)/cphi0)
        rlamda=acos(cphiicrlamda/cos(rphii))
      if(flon.lt.lamda0) rlamda=-rlamda

        ri2=(rlamda/rdlamda+nwlon+1)/2.
      ri1=(rlamda/rdlamda+nwlon)/2.
      rj=rphii/rdphi+(nwlat+1)/2
      j=(rj+0.5)
      ri=ri1
      if(mod(j,2).eq.0) ri=ri2
      i=(ri+0.5)

c      write(*,*) "lamda, phi = ",flon,flat
c      write(*,*) "ri1, ri2, ri, rj = ",ri1,ri2,ri,rj
c      write(*,*) "i, j = ",i,j


        if((rj.lt.1).or.(rj.gt.nwlat)) then
c        write(333,*) 'ADCIRC grid node ',node,
c     &             ' falls outside of the ETA 29 grid'
        icode=1
        NN1=1
        NN2=1
        NN3=1
        DD1=0
        DD2=0
        DD3=0
        return
        endif

        if(mod(j,2).eq.0) then
        if((ri.lt.1).or.(ri.gt.(nwlon+0.5))) then
c          write(333,*) 'ADCIRC grid node ',node,
c     &                 ' falls outside of the ETA 29 grid'
          icode=1
          NN1=1
          NN2=1
          NN3=1
          DD1=0
          DD2=0
          DD3=0
          return
          endif
        endif

        if(mod(j,2).ne.0) then
        if((ri.lt.0.5).or.(ri.gt.nwlon)) then
c          write(333,*) 'ADCIRC grid node ',node,
c     &                 ' falls outside of the ETA 29 grid'
          icode=1
          NN1=1
          NN2=1
          NN3=1
          DD1=0
          DD2=0
          DD3=0
          return
          endif
        endif

c     compute the coordinates of the closest Eta29 grid node

      jm2=(nwlat+1)/2
      im2=nwlon*2
      call e29calc(i,j,lamda,phi,n)

c     compute the coordinates of neighbor node "a" (located SW of closest node)

      if((i.eq.1).and.(mod(j,2).eq.0)) then
        ia=i
        ja=j-2
        else
        ia=i
        if(mod(j,2).eq.0) ia=i-1
        ja=j-1
        endif
      if((ia.lt.1).or.(ja.lt.1)) then  !this neighbor lies outside of Eta29 grid
        na=0
        else
        call e29calc(ia,ja,lamdaa,phia,na)
        endif

c     compute the coordinates of neighbor node "b" (located W of closest node)

      ib=i-1
      jb=j
      if(ib.lt.1) then  !this neighbor lies outside of Eta29 grid
        nb=0
        else
        call e29calc(ib,jb,lamdab,phib,nb)
        endif

c     compute the coordinates of neighbor node "c" (located NW of closest node)

      if((i.eq.1).and.(mod(j,2).eq.0)) then
        ic=i
        jc=j+2
        else
        ic=ia
        jc=j+1
        endif
      if((ic.lt.1).or.(jc.gt.nwlat)) then  !this neighbor lies outside of Eta29 grid
        nc=0
        else
        call e29calc(ic,jc,lamdac,phic,nc)
        endif

c     compute the coordinates of neighbor node "d" (located NE of closest node)

      if((i.eq.181).and.(mod(j,2).ne.0)) then
        id=i
        jd=j+2
        else
        id=ic+1
        jd=j+1
        endif
      if((id.gt.nwlon).or.(jd.gt.nwlat)) then  !this neighbor lies outside of Eta29 grid
        nd=0
        else
        call e29calc(id,jd,lamdad,phid,nd)
        endif

c     compute the coordinates of neighbor node "e" (located E of closest node)

      ie=i+1
      je=j
      if(ie.gt.nwlon) then  !this neighbor lies outside of Eta29 grid
        ne=0
        else
        call e29calc(ie,je,lamdae,phie,ne)
        endif

c     compute the coordinates of neighbor node "g" (located SE of closest node)

      if((i.eq.181).and.(mod(j,2).ne.0)) then
        ig=i
        jg=j-2
        else
        ig=id
        jg=j-1
        endif
      if((ig.gt.nwlon).or.(jg.lt.1)) then  !this neighbor lies outside of Eta29 grid
        ng=0
        else
        call e29calc(ig,jg,lamdag,phig,ng)
        endif

c      write(*,*) 'closest E29 node i,j = ',n,i,j,lamda,phi
c      if(na.eq.0) write(*,*) 'point a falls outside of Eta29 grid'
c      if(na.ne.0) write(*,*) 'point a   = ',na,ia,ja,lamdaa,phia
c      if(nb.eq.0) write(*,*) 'point b falls outside of Eta29 grid'
c      if(nb.ne.0) write(*,*) 'point b   = ',nb,ib,jb,lamdab,phib
c      if(nc.eq.0) write(*,*) 'point c falls outside of Eta29 grid'
c      if(nc.ne.0) write(*,*) "point c   = ",nc,ic,jc,lamdac,phic
c      if(nd.eq.0) write(*,*) 'point d falls outside of Eta29 grid'
c      if(nd.ne.0) write(*,*) "point d   = ",nd,id,jd,lamdad,phid
c      if(ne.eq.0) write(*,*) 'point e falls outside of Eta29 grid'
c      if(ne.ne.0) write(*,*) "point e   = ",ne,ie,je,lamdae,phie
c      if(ng.eq.0) write(*,*) 'point g falls outside of Eta29 grid'
c      if(ng.ne.0) write(*,*) "point g   = ",ng,ig,jg,lamdag,phig

        NN1=1
      NN2=1
      NN3=1
      DD1=0
      DD2=0
      DD3=0
      X1=lamda
      X4=flon
      Y1=phi
      Y4=flat
      ifflag=0
      AEMIN=99999.

c     test if the point is in triangle ij - b - a

      if((na.ne.0).and.(nb.ne.0)) then
        X2=lamdab
        X3=lamdaa
        Y2=phib
        Y3=phia
        AREAS=ABS((X1-X3)*(Y2-Y3)+(X3-X2)*(Y1-Y3))
        A1=(X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4)
        A2=(X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1)
        A3=(Y4-Y1)*(X2-X1)-(X4-X1)*(Y2-Y1)
        AA=ABS(A1)+ABS(A2)+ABS(A3)
        AE=ABS(AA-AREAS)/AREAS
c        write(333,*) "AE = ",AE
        IF((AE.LT.1.0E-5).AND.(AE.LT.AEMIN)) THEN
          AEMIN=AE
          NN1=n
          NN2=nb
          NN3=na
          DD1=((X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4))/AREAS
          DD2=((X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1))/AREAS
          DD3=(-(X4-X1)*(Y2-Y1)+(Y4-Y1)*(X2-X1))/AREAS
          call betaucalc(i,j,DD1,ib,jb,DD2,ia,ja,DD3,betau)
          ifflag=ifflag+1
c          write(333,*) 'position found in triangle ij - b - a'
          ENDIF
        endif

c     if along the west boundary, test if the point is in triangle ij - c - a

      if((i.eq.1).and.(mod(j,2).ne.0)) then
        if((na.ne.0).and.(nc.ne.0)) then
          X2=lamdac
          X3=lamdaa
          Y2=phic
          Y3=phia
          AREAS=ABS((X1-X3)*(Y2-Y3)+(X3-X2)*(Y1-Y3))
          A1=(X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4)
          A2=(X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1)
          A3=(Y4-Y1)*(X2-X1)-(X4-X1)*(Y2-Y1)
          AA=ABS(A1)+ABS(A2)+ABS(A3)
          AE=ABS(AA-AREAS)/AREAS
c          write(333,*) "AE = ",AE
          IF((AE.LT.1.0E-5).AND.(AE.LT.AEMIN)) THEN
            NN1=n
            NN2=nc
            NN3=na
            DD1=((X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4))/AREAS
            DD2=((X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1))/AREAS
            DD3=(-(X4-X1)*(Y2-Y1)+(Y4-Y1)*(X2-X1))/AREAS
            call betaucalc(i,j,DD1,ic,jc,DD2,ia,ja,DD3,betau)
            ifflag=ifflag+1
c            write(333,*) 'position found in triangle ij - c - a'
            ENDIF
          endif
        endif

c     test if the point is in triangle ij - c - b

      if((nb.ne.0).and.(nc.ne.0)) then
        X2=lamdac
        X3=lamdab
        Y2=phic
        Y3=phib
        AREAS=ABS((X1-X3)*(Y2-Y3)+(X3-X2)*(Y1-Y3))
        A1=(X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4)
        A2=(X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1)
        A3=(Y4-Y1)*(X2-X1)-(X4-X1)*(Y2-Y1)
        AA=ABS(A1)+ABS(A2)+ABS(A3)
        AE=ABS(AA-AREAS)/AREAS
c        write(333,*) "AE = ",AE
        IF((AE.LT.1.0E-5).AND.(AE.LT.AEMIN)) THEN
          NN1=n
          NN2=nc
          NN3=nb
          DD1=((X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4))/AREAS
          DD2=((X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1))/AREAS
          DD3=(-(X4-X1)*(Y2-Y1)+(Y4-Y1)*(X2-X1))/AREAS
          call betaucalc(i,j,DD1,ic,jc,DD2,ib,jb,DD3,betau)
          ifflag=ifflag+1
c          write(333,*) 'position found in triangle ij - c - b'
          ENDIF
        endif

c     test if the point is in triangle ij - d - c

      if((nc.ne.0).and.(nd.ne.0)) then
        X2=lamdad
        X3=lamdac
        Y2=phid
        Y3=phic
        AREAS=ABS((X1-X3)*(Y2-Y3)+(X3-X2)*(Y1-Y3))
        A1=(X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4)
        A2=(X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1)
        A3=(Y4-Y1)*(X2-X1)-(X4-X1)*(Y2-Y1)
        AA=ABS(A1)+ABS(A2)+ABS(A3)
        AE=ABS(AA-AREAS)/AREAS
c        write(333,*) "AE = ",AE
        IF((AE.LT.1.0E-5).AND.(AE.LT.AEMIN)) THEN
          NN1=n
          NN2=nd
          NN3=nc
          DD1=((X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4))/AREAS
          DD2=((X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1))/AREAS
          DD3=(-(X4-X1)*(Y2-Y1)+(Y4-Y1)*(X2-X1))/AREAS
          call betaucalc(i,j,DD1,id,jd,DD2,ic,jc,DD3,betau)
          ifflag=ifflag+1
c          write(333,*) 'position found in triangle ij - d - c'
          ENDIF
        endif

c     if along the east boundary, test if the point is in triangle ij - g - d

      if((i.eq.181).and.(mod(j,2).eq.0)) then
        if((nd.ne.0).and.(ng.ne.0)) then
          X2=lamdag
          X3=lamdad
          Y2=phig
          Y3=phid
          AREAS=ABS((X1-X3)*(Y2-Y3)+(X3-X2)*(Y1-Y3))
          A1=(X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4)
          A2=(X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1)
          A3=(Y4-Y1)*(X2-X1)-(X4-X1)*(Y2-Y1)
          AA=ABS(A1)+ABS(A2)+ABS(A3)
          AE=ABS(AA-AREAS)/AREAS
c          write(333,*) "AE = ",AE
          IF((AE.LT.1.0E-5).AND.(AE.LT.AEMIN)) THEN
            NN1=n
            NN2=ng
            NN3=nd
            DD1=((X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4))/AREAS
            DD2=((X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1))/AREAS
            DD3=(-(X4-X1)*(Y2-Y1)+(Y4-Y1)*(X2-X1))/AREAS
            call betaucalc(i,j,DD1,ig,jg,DD2,id,jd,DD3,betau)
            ifflag=ifflag+1
c            write(333,*) 'position found in triangle ij - g - d'
            ENDIF
          endif
        endif

c     test if the point is in triangle ij - e - d

      if((nd.ne.0).and.(ne.ne.0)) then
        X2=lamdae
        X3=lamdad
        Y2=phie
        Y3=phid
        AREAS=ABS((X1-X3)*(Y2-Y3)+(X3-X2)*(Y1-Y3))
        A1=(X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4)
        A2=(X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1)
        A3=(Y4-Y1)*(X2-X1)-(X4-X1)*(Y2-Y1)
        AA=ABS(A1)+ABS(A2)+ABS(A3)
        AE=ABS(AA-AREAS)/AREAS
c        write(333,*) "AE = ",AE
        IF((AE.LT.1.0E-5).AND.(AE.LT.AEMIN)) THEN
          NN1=n
          NN2=ne
          NN3=nd
          DD1=((X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4))/AREAS
          DD2=((X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1))/AREAS
          DD3=(-(X4-X1)*(Y2-Y1)+(Y4-Y1)*(X2-X1))/AREAS
          call betaucalc(i,j,DD1,ie,je,DD2,id,jd,DD3,betau)
          ifflag=ifflag+1
c          write(333,*) 'position found in triangle ij - e - d'
          ENDIF
        endif

c     test if the point is in triangle ij - g - e

      if((ne.ne.0).and.(ng.ne.0)) then
        X2=lamdag
        X3=lamdae
        Y2=phig
        Y3=phie
        AREAS=ABS((X1-X3)*(Y2-Y3)+(X3-X2)*(Y1-Y3))
        A1=(X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4)
        A2=(X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1)
        A3=(Y4-Y1)*(X2-X1)-(X4-X1)*(Y2-Y1)
        AA=ABS(A1)+ABS(A2)+ABS(A3)
        AE=ABS(AA-AREAS)/AREAS
c        write(333,*) "AE = ",AE
        IF((AE.LT.1.0E-5).AND.(AE.LT.AEMIN)) THEN
          NN1=n
          NN2=ng
          NN3=ne
          DD1=((X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4))/AREAS
          DD2=((X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1))/AREAS
          DD3=(-(X4-X1)*(Y2-Y1)+(Y4-Y1)*(X2-X1))/AREAS
          call betaucalc(i,j,DD1,ig,jg,DD2,ie,je,DD3,betau)
          ifflag=ifflag+1
c          write(333,*) 'position found in triangle ij - g - e'
          ENDIF
        endif

c     test if the point is in triangle ij - a - g

      if((na.ne.0).and.(ng.ne.0)) then
        X2=lamdaa
        X3=lamdag
        Y2=phia
        Y3=phig
        AREAS=ABS((X1-X3)*(Y2-Y3)+(X3-X2)*(Y1-Y3))
        A1=(X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4)
        A2=(X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1)
        A3=(Y4-Y1)*(X2-X1)-(X4-X1)*(Y2-Y1)
        AA=ABS(A1)+ABS(A2)+ABS(A3)
        AE=ABS(AA-AREAS)/AREAS
c        write(333,*) "AE = ",AE
        IF((AE.LT.1.0E-5).AND.(AE.LT.AEMIN)) THEN
          NN1=n
          NN2=na
          NN3=ng
          DD1=((X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4))/AREAS
          DD2=((X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1))/AREAS
          DD3=(-(X4-X1)*(Y2-Y1)+(Y4-Y1)*(X2-X1))/AREAS
          call betaucalc(i,j,DD1,ia,ja,DD2,ig,jg,DD3,betau)
          ifflag=ifflag+1
c          write(333,*) 'position found in triangle ij - a - g'
          ENDIF
        endif

c      if(ifflag.eq.0) then
c        write(333,*) 'position not found'
c        write(*,*) 'position not found in subroutine E29SEARCH'
c        icode=3
c        else
c       write(*,*) 'i,j,NN1,NN2,NN3,DD1,DD2,DD3'
c        write(333,999) i,j,NN1,NN2,NN3,DD1,DD2,DD3,betau/dgtora
c 999    format(5I8,1x,3E13.6)
c        endif

      return
      end

C***********************************************************************
C  Subroutine to compute the longititude and latitude of a given i,j   *
C       position in the Eta29 grid.                                    *
C                                                                      *
C                    Written by R.L.       1/11/98                     *
C***********************************************************************

      subroutine e29calc(i,j,lamda,phi,n)
      use presizes, only : sz
      implicit none
      real(8) lamda0,phi0,rphi0,cphi0,sphi0,tphi0,dlamda,dphi,rdlamda,
     &       rdphi,a,rlamar,phiarg,rlamda,dgtora
      real(8) lamda,phi
C     jgf46.00 Explicitly declared the following variables
      integer nwlon, nwlat, jm2, im2, j, i1, i, i2, n, i1p1, i1m1, j1
      real(sz) phii, dlon, dlat, dlnt, arg, betau1
      intrinsic :: mod, sin, cos, sqrt, asin

      nwlon=181
      nwlat=271
      dgtora=3.1415926d0/180.d0
      lamda0=-97.0d0
      phi0=41.0d0
      rphi0=dgtora*phi0
      cphi0=cos(rphi0)
      sphi0=sin(rphi0)
      tphi0=tan(rphi0)
      dlamda=7.d0/36.d0
      dphi=5.d0/27.d0
      rdlamda=dgtora*dlamda
      rdphi=dgtora*dphi

      jm2=(nwlat+1)/2
      im2=nwlon*2

      phii=rdphi*float(j-jm2)
      i1=2*i-1
      i2=2*i
      if(mod(j,2).ne.0) then
        rlamda=rdlamda*float(i2-nwlon)
        else
        rlamda=rdlamda*float(i1-nwlon)
        endif
      phiarg= sin(phii)*cphi0+cos(phii)*sphi0*cos(rlamda)
      if(phiarg.gt.1.0) phiarg=1.0
      if(phiarg.lt.-1.0) phiarg=-1.0
      phi=asin(phiarg)
      rlamar= cos(phii)*cos(rlamda)/(cos(phi)*cphi0)-tan(phi)*tphi0
      if(rlamar.gt.1.0) rlamar=1.0
      if(rlamar.lt.-1.) rlamar=-1.
      a=acos(rlamar)/dgtora
      if(rlamda.le.0.) then
        lamda=lamda0-a
        else
        lamda=lamda0+a
        endif
      phi=phi/dgtora
      n=nwlon*(j-1)+i

      return
      end


C***********************************************************************
C  Subroutine to compute the conversion angle between the E29 velocity *
C       field and a lon,lat coordinate system.                         *
C                                                                      *
C                    Written by R.L.       1/12/98                     *
C***********************************************************************

      subroutine betaucalc(i1,j1,dd1,i2,j2,dd2,i3,j3,dd3,betau)
      use presizes, only : sz
      implicit none
      real(8) lamda,lamdap1,lamdam1,phi,phip1,phim1
C     jgf46.00 Explicitly declared the following variables
      integer i1, i1p1, i1m1, n, i2, i2p1, i2m1, j2, i3, i3p1
      integer i3m1, j3, j1
      real(sz) dgtora, dlon, dlat, dlnt, arg, betau1, betau2, betau3
      real(sz) betau, dd1, dd2, dd3
      intrinsic :: asin



      dgtora=3.1415926d0/180.d0

        if(i1.ne.181) then
        i1p1=i1+1
        else
        i1p1=i1
        endif
      if(i1.ne.1) then
        i1m1=i1-1
        else
        i1m1=i1
        endif
      call e29calc(i1,j1,lamda,phi,n)
      call e29calc(i1p1,j1,lamdap1,phip1,n)
      call e29calc(i1m1,j1,lamdam1,phim1,n)
      dlon=(lamdap1-lamdam1)*cos(phi*dgtora)
      dlat=phip1-phim1
      dlnt=sqrt(dlon*dlon+dlat*dlat)
      arg=dlat/dlnt
      if(arg.gt.1.) arg=1.
      if(arg.lt.-1.) arg=-1.
      betau1=asin(arg)

        if(i2.ne.181) then
        i2p1=i2+1
        else
        i2p1=i2
        endif
      if(i2.ne.1) then
        i2m1=i2-1
        else
        i2m1=i2
        endif
      call e29calc(i2,j2,lamda,phi,n)
      call e29calc(i2p1,j2,lamdap1,phip1,n)
      call e29calc(i2m1,j2,lamdam1,phim1,n)
      dlon=(lamdap1-lamdam1)*cos(phi*dgtora)
      dlat=phip1-phim1
      dlnt=sqrt(dlon*dlon+dlat*dlat)
      arg=dlat/dlnt
      if(arg.gt.1.) arg=1.
      if(arg.lt.-1.) arg=-1.
      betau2=asin(arg)

        if(i3.ne.181) then
        i3p1=i3+1
        else
        i3p1=i3
        endif
      if(i3.ne.1) then
        i3m1=i3-1
        else
        i3m1=i3
        endif
      call e29calc(i3,j3,lamda,phi,n)
      call e29calc(i3p1,j3,lamdap1,phip1,n)
      call e29calc(i3m1,j3,lamdam1,phim1,n)
      dlon=(lamdap1-lamdam1)*cos(phi*dgtora)
      dlat=phip1-phim1
      dlnt=sqrt(dlon*dlon+dlat*dlat)
      arg=dlat/dlnt
      if(arg.gt.1.) arg=1.
      if(arg.lt.-1.) arg=-1.
      betau3=asin(arg)

        betau=dd1*betau1+dd2*betau2+dd3*betau3

        return
      end

C***********************************************************************
C                                                                      *
C   End of subroutines to read wind and pressure fields.               *
C                                                                      *
C***********************************************************************
