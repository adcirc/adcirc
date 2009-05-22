C***********************************************************************
C PADCIRC VERSION 45.12 03/17/2006                                     *
C last changes in this file VERSION 45.05                              *
C                                                                      *
C***********************************************************************
C
       MODULE WIND
       USE SIZES, ONLY : SZ, LOCALDIR
       USE GLOBAL, ONLY : ScreenUnit
#ifdef CMPI
      use messenger, only : msg_fini
#endif

C
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
C
      REAL(8),PRIVATE,PARAMETER :: PI=3.141592653589793D0, 
     &                             TWOPI=PI*2.D0,
     &                             HFPI=PI/2.D0,
     &                             RAD2DEG = 180.D0/PI,
     &                             DEG2RAD = PI/180.D0
      INTEGER StormNumber !jgf46.28 For Holland Wind model wind multiplier
      REAL(8) WindRefTime !jgf46.29 seconds since beginning of year, this
                          !corresponds to time=0 of the simulation 
      REAL(8) BLAdj       !jgf46.32 boundary layer adjustment for Holland model

C------------------------end of data declarations-----------------------
       CONTAINS 
C-----------------------------------------------------------------------

C     ----------------------------------------------------------------
C      F U N C T I O N   W I N D   D R A G
C     ----------------------------------------------------------------
C
C     jgf46.00 Function to calculate wind drag coefficient from the
C     windspeed.
C       
C     ----------------------------------------------------------------
      REAL(SZ) FUNCTION WindDrag(WindSpeed, DragLimit, 
     &      DragLawString)
      REAL(SZ) WindSpeed         
      REAL(SZ) DragLimit
      CHARACTER(len=10) DragLawString
C
      SELECT CASE(TRIM(DragLawString))
C
      CASE("Garratt")
         WindDrag = 0.001d0*(0.75d0+0.067d0*WindSpeed)
C
      CASE DEFAULT
         WRITE(16,*) 'ERROR: Wind drag law not recognized:'
         WRITE(16,'(A10)') DragLawString
         WRITE(16,*) 'Execution will now be terminated.'
#ifdef CMPI
         call msg_fini()
#endif
         STOP
C
      END SELECT
      IF(WindDrag.gt.DragLimit) WindDrag=DragLimit
C
      RETURN
C     ----------------------------------------------------------------
      END FUNCTION WindDrag
C     ----------------------------------------------------------------


C***********************************************************************
C                                                                      *
C   Convert time from year,month,day,hour,min,sec into seconds since   *
C   the beginning of the year.                                         *
C                                                                      *
C***********************************************************************

      SUBROUTINE TIMECONV(IYR,IMO,IDAY,IHR,IMIN,SEC,TIMESEC, 
     &     MyProc, NScreen, ScreenUnit)
      IMPLICIT NONE
      INTEGER IYR,IMO,IDAY,IHR,IMIN,ILEAP
      INTEGER MyProc, NScreen, ScreenUnit
      REAL*8 TIMESEC,SEC
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
         IF (NScreen.ne.0.and.MyProc.eq.0) THEN
            WRITE(ScreenUnit,*) 
     &           'FATAL ERROR IN SUBROUTINE TIMECONV - MONTH > 12 '
         ENDIF
         WRITE(16,*) 'FATAL ERROR IN SUBROUTINE TIMECONV - MONTH > 12 '
#ifdef CMPI
         call msg_fini()
#endif
         STOP
      ENDIF
      RETURN
      END SUBROUTINE TIMECONV

C***********************************************************************
C                                                                      *
C   READ IN AND INTERPOLATE ONTO THE ADCIRC GRID WIND FIELDS FROM U.S. *
C   NAVY FLEET NUMERIC WIND FILES.                                     *
C                                                                      *
C   NOTE: The ADCIRC grid information consists only of the Lon and Lat *
C   of the nodes.  THE LONS AND LATS MUST BE IN RADIANS!               *
C                                                                      *
C                                                                      *
C   NWLAT = MAXIMUM NUMBER OF LATITUDES IN FLEET NUMERIC WIND FILE     *
C            SET = 1 IF FLEET NUMERIC WIND FILE NOT IN USE             *
C   NWLON = MAXIMUM NUMBER OF LONGITUDES IN FLEET NUMERIC WIND FILE    *
C            SET = 1 IF FLEET NUMERIC WIND FILE NOT IN USE             *
C                                                                      *
C                        R.L. 4/17/96                                  *
C                                                                      *
C   R.L. 4/2/01  changed MNWLAT,MNWLON in ALLOCATE statement to        *
C                NWLAT,NWLON                                           *
C***********************************************************************

      SUBROUTINE NWS3GET(X,Y,SLAM,SFEA,WVNX,WVNY,IWTIME,IWYR,WTIMED,NP,
     &     NWLON,NWLAT,WLATMAX,WLONMIN,WLATINC,WLONINC,ICS,
     &     NScreen,ScreenUnit)
      USE SIZES
      IMPLICIT NONE
      INTEGER, SAVE :: FIRSTCALL = 0
      INTEGER IWTIME,IWYR,IWMO,IWDAY,IWHR,NP,NWLON,NWLAT,ICS,I,J
      INTEGER NScreen,ScreenUnit
      REAL*8 WTIMED
      REAL*8 X(*),Y(*),SLAM(*),SFEA(*),XCOOR,YCOOR
      INTEGER  LATIND1,LATIND2,LONIND1,LONIND2
      REAL(SZ) WLATMAX,WLONMIN,WLATINC,WLONINC,WSPEED,WDIR
      REAL(SZ) WLATM,WLONM,XWRATIO,YWRATIO
      REAL(SZ),ALLOCATABLE,SAVE :: WVXFN(:,:),WVYFN(:,:),PRN(:,:)
      REAL(SZ) WVNX(*),WVNY(*)
C
      IF (FIRSTCALL.EQ.0) THEN
         FIRSTCALL = 1
         ALLOCATE ( WVXFN(NWLAT,NWLON),WVYFN(NWLAT,NWLON),
     &             PRN(NWLAT,NWLON) )
      ENDIF
C
      READ(22,*) IWTIME
      IWYR = IWTIME/1000000
      IWMO = IWTIME/10000 - IWYR*100
      IWDAY = IWTIME/100 - IWYR*10000 - IWMO*100
      IWHR = IWTIME - IWYR*1000000 - IWMO*10000 - IWDAY*100
      CALL TIMECONV(IWYR,IWMO,IWDAY,IWHR,0,0.0D0,WTIMED,
     &     MyProc,NScreen,ScreenUnit)
C
      DO I=1,NWLAT
         READ(22,*) (WVXFN(I,J),J=1,NWLON)
      END DO
      DO I=1,NWLAT
         READ(22,*) (WVYFN(I,J),J=1,NWLON)
      END DO
C
      DO I=1,NWLAT              !CONVERT TO X AND Y COMPONENTS
         DO J=1,NWLON
            WSPEED=WVXFN(I,J)
            WDIR=WVYFN(I,J)*DEG2RAD
            WVXFN(I,J)=-WSPEED*SIN(WDIR)
            WVYFN(I,J)=-WSPEED*COS(WDIR)
         END DO
      END DO
      
      DO I=1,NP                 !INTERPOLATE TO ADCIRC GRID
         IF(ICS.EQ.2) THEN
            YCOOR=SFEA(I)*RAD2DEG
            XCOOR=SLAM(I)*RAD2DEG
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
C     
         WVNX(I) = WVXFN(LATIND2,LONIND2)*XWRATIO*YWRATIO
     &        + WVXFN(LATIND2,LONIND1)*(1.d0-XWRATIO)*YWRATIO
     &        + WVXFN(LATIND1,LONIND2)*XWRATIO*(1.d0-YWRATIO)
     &        + WVXFN(LATIND1,LONIND1)*(1.d0-XWRATIO)*(1.d0-YWRATIO)
         WVNY(I) = WVYFN(LATIND2,LONIND2)*XWRATIO*YWRATIO
     &        + WVYFN(LATIND2,LONIND1)*(1.d0-XWRATIO)*YWRATIO
     &        + WVYFN(LATIND1,LONIND2)*XWRATIO*(1.d0-YWRATIO)
     &        + WVYFN(LATIND1,LONIND1)*(1.d0-XWRATIO)*(1.d0-YWRATIO)
      END DO
C     
      RETURN
      END SUBROUTINE NWS3GET



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
C   R.L.09/04/00 added RHOWAT0 to call                                 *
C   R.L. 4/2/01  changed MNP dimensions to *                           *
C   R.L. 3/15/03 accounted for PRN=0                                   *   !RAL0315+ OK
C***********************************************************************

      SUBROUTINE NWS4GET(WVNX,WVNY,PRN,NP,RHOWAT0,G)
      USE SIZES
      IMPLICIT NONE
      INTEGER   NP,I,NHG
      REAL(SZ)  WVNX(*),WVNY(*),PRN(*)
      REAL(SZ)  RHOWAT0,RHOWATG,G
      CHARACTER*80 PBLJAGF
C
      RHOWATG=RHOWAT0*G
      DO I=1,NP
        WVNX(I)=0.d0
        WVNY(I)=0.d0
        PRN(I)=101300.d0/RHOWATG
      END DO
 170  READ(22,'(A80)') PBLJAGF
      IF(PBLJAGF(2:2).EQ.'#') GOTO 170
 171  READ(PBLJAGF,'(I8,5E13.5)') NHG,WVNX(NHG),WVNY(NHG),
     &                            PRN(NHG)
C
C     jgf46.02 From now on, wind files must contain data that are
C     appropriate for the time increment listed within the
C     file. Therefore, the following two lines were commented out.
C     WVNX(NHG)=WVNX(NHG)*1.04d0*0.5144d0 !CONVERT 30-MIN WINDS IN
C     WVNY(NHG)=WVNY(NHG)*1.04d0*0.5144d0 !KNOTS TO 10-MIN WIND IN M/S
C     jgf46.02 Added the following two lines.
      WVNX(NHG)=WVNX(NHG)*0.5144d0 !CONVERT KNOTS TO  M/S
      WVNY(NHG)=WVNY(NHG)*0.5144d0 
      PRN(NHG)=100.d0*PRN(NHG)/RHOWATG !CONVERT MILLIBARS TO M OF WATER
      IF(PRN(NHG).EQ.0.) PRN(NHG)=101300.d0/RHOWATG                       !RAL0315+ OK
      READ(22,'(A80)') PBLJAGF
      IF(PBLJAGF(2:2).NE.'#') GOTO 171
      RETURN
      END SUBROUTINE


C***********************************************************************
C                                                                      *
C   Read in and interpolate onto the ADCIRC grid wind and pressure     *
C   fields from a meteorological file on a rectangular grid (either in *
C   Longitude, Latitude or Cartesian coordinates, consistent with the  *
C   ADCIRC grid coordinates).  If the ADCIRC grid is in Lon,Lat these  *
C   MUST BE IN RADIANS!                                                *
C                                                                      *
C   It is assumed that the meteorological grid is set up so that y     *
C   (e.g., latitude) varies from north (k=1) to south (k=NWLAT) and x  *
C   (e.g., longitude) varies from west (j=1) to east (j=NWLON).        *
C                                                                      *
C   The spatial extents of the meteorological grid must be consistent  *
C   with the ADCIRC model domain.  For example, if ADCIRC uses negative*
C   longitude values to indicate locations W of the Greenwich meridian,*
C   the meteorological file must be similarly organized.  Any grid that*
C   crosses the Greenwich Meridian should be organized so that the seam*
C   occurs @ 180 deg longitude. Therefore, the meteorological and      *
C   ADCIRC grids should use negative longitudes W of the Greenwich     *
C   Meridian and positive longitudes to the E.                         *
C                                                                      *
C                                                                      *
C   NOTE:  It is assumed that the met file data is oriented so that    *
C          the outer loop is on latitude and the inner loop is on      *
C          longitude.  For example:                                    *
C          line 1             lat 1,     lon 1                         *
C          line 2             lat 1,     lon 2                         *
C            .                                                         *
C          line nwlon         lat 1,     lon nwlon                     *
C          line nwlon+1       lat 2,     lon 1                         *
C          line nwlon+2       lat 2,     lon 2                         *
C            .                                                         *
C          line 2*nwlon       lat 2,     lon nwlon                     *
C          line 2*nwlon+1     lat 3,     lon 1                         *
C          line 2*nwlon+2     lat 3,     lon 2                         *
C            .                                                         *
C          line nwlon*nwlat   lat nwlat, lon nwlon                     *
C                                                                      *
C   NOTE:  It is assumed that he met file data is oriented so that     *
C          latitude varies from the northern most value (lat 1) to the *
C          southern most value (lat nwlat) and longitude varies in an  *
C          easterly direction (e.g. from 0 to 360 where positive       *
C          longitudes are angles measured easterly of the GM.          *
C                                                                      *
C   NOTE:  For the global AVN grid running from 0.5 - 359.5 deg lon    *
C          and 90 - -90 deg lat in 1 degree increments, NWLAT=181 and  *
C          NWLON=360 yielding a total number of entries in the file    *
C          of 65160.                                                   *    
C                                                                      *
C   NOTE:  It is assumed that wind velocity is in EAST,NORTH components*
C          in M/2 and pressure is in N/M^2                             *
C                                                                      *
C   NOTE:  WLATMAX,WLONMIN,WLATINC,WLONINC should be in deg.           *
C                                                                      *
C   NOTE:  This should wrap if XCOORD > WLONMIN+NWLON*WLONINC  or      *
C          XCOORD < WLONMIN                                            *
C                                                                      *
C                                                                      *
C   MNWLAT = MAXIMUM NUMBER OF LATITUDES IN WIND FILE                  *
C            SET = 1 IF FLEET NUMERIC WIND FILE NOT IN USE             *
C   MNWLON = MAXIMUM NUMBER OF LONGITUDES IN WIND FILE                 *
C            SET = 1 IF FLEET NUMERIC WIND FILE NOT IN USE             *
C                                                                      *
C                           R.L. 4/13/99                               *
C                           R.L.09/04/00 added RHOWAT0 to call         *
C   R.L.09/04/00 added RHOWAT0 to call                                 *
C   R.L. 4/2/01  changed MNWLAT,MNWLON in ALLOCATE statement to        *
C                NWLAT,NWLON                                           *
C   R.L. 8/10/05 eliminated adding 360 to negative longitudes to match *
C                AVN model grid setup.  User is now required to provide*
C                met and ADCIRC grid that correspond in space without  *
C                adjusting the longitude.  Also the input variable     *
C                order has been changed to U,V,P to be consistent with *
C                other NWS input formats.                              *
C***********************************************************************

      SUBROUTINE NWS6GET(X,Y,SLAM,SFEA,WVNX,WVNY,PRESS,NP,NWLON,NWLAT,
     &     WLATMAX,WLONMIN,WLATINC,WLONINC,ICS,RHOWAT0,G)
      USE SIZES

      IMPLICIT NONE
      INTEGER, SAVE :: FIRSTCALL = 0
      INTEGER NP,NWLON,NWLAT,I,J,ICS 
      REAL(SZ) RHOWAT0,RHOWATG,G
      INTEGER  LATIND1,LATIND2,LONIND1,LONIND2
      REAL(SZ) WLATMAX,WLONMIN,WLATINC,WLONINC,XWRATIO,YWRATIO
      REAL(SZ) WLATM,WLONM
      REAL*8 X(*),Y(*),SLAM(*),SFEA(*),XCOOR,YCOOR
      REAL(SZ) WVNX(*),WVNY(*),PRESS(*)
      REAL(SZ),SAVE,ALLOCATABLE :: WVXFN(:,:),WVYFN(:,:),PRN(:,:)
C     
      IF (FIRSTCALL.EQ.0) THEN
         FIRSTCALL = 1
         ALLOCATE ( WVXFN(NWLAT,NWLON),WVYFN(NWLAT,NWLON),
     &        PRN(NWLAT,NWLON) )
      ENDIF
C     
      RHOWATG=RHOWAT0*G
      DO I=1,NWLAT
         DO J=1,NWLON
            READ(22,*) WVXFN(I,J),WVYFN(I,J),PRN(I,J)
         END DO
      END DO

      DO I=1,NP                 !INTERPOLATE TO ADCIRC GRID
         IF(ICS.EQ.2) THEN
            YCOOR=SFEA(I)*RAD2DEG
            XCOOR=SLAM(I)*RAD2DEG
         ENDIF
         IF(ICS.EQ.1) THEN
            YCOOR=Y(I)
            XCOOR=X(I)
         ENDIF
         LATIND2=(WLATMAX-YCOOR)/WLATINC + 1
         IF(LATIND2.EQ.NWLAT) LATIND2=LATIND2-1
         LATIND1=LATIND2 + 1
         LONIND1=(XCOOR-WLONMIN)/WLONINC + 1
         LONIND2=LONIND1 + 1
C     
         WLONM = WLONMIN + (LONIND1-1)*WLONINC 
         WLATM = WLATMAX - (LATIND1-1)*WLATINC
         XWRATIO=(XCOOR-WLONM)/WLONINC
         YWRATIO=(YCOOR-WLATM)/WLATINC
C     
         IF(LONIND1.EQ.0) LONIND1=NWLON
         IF(LONIND1.EQ.NWLON) LONIND2=1
C     
         WVNX(I) = WVXFN(LATIND2,LONIND2)*XWRATIO*YWRATIO
     &        + WVXFN(LATIND2,LONIND1)*(1.d0-XWRATIO)*YWRATIO
     &        + WVXFN(LATIND1,LONIND2)*XWRATIO*(1.d0-YWRATIO)
     &        + WVXFN(LATIND1,LONIND1)*(1.d0-XWRATIO)*(1.d0-YWRATIO)
         WVNY(I) = WVYFN(LATIND2,LONIND2)*XWRATIO*YWRATIO
     &        + WVYFN(LATIND2,LONIND1)*(1.d0-XWRATIO)*YWRATIO
     &        + WVYFN(LATIND1,LONIND2)*XWRATIO*(1.d0-YWRATIO)
     &        + WVYFN(LATIND1,LONIND1)*(1.d0-XWRATIO)*(1.d0-YWRATIO)
         PRESS(I) = PRN(LATIND2,LONIND2)*XWRATIO*YWRATIO
     &        + PRN(LATIND2,LONIND1)*(1.d0-XWRATIO)*YWRATIO
     &        + PRN(LATIND1,LONIND2)*XWRATIO*(1.d0-YWRATIO)
     &        + PRN(LATIND1,LONIND1)*(1.d0-XWRATIO)*(1.d0-YWRATIO)
         PRESS(I) = PRESS(I)/RHOWATG

      END DO
C     
      RETURN
      END SUBROUTINE NWS6GET

C     ----------------------------------------------------------------
C      S U B R O U T I N E     N W S 7 G E T 
C     ----------------------------------------------------------------
C
C     jgf46.01 Subroutine to get a set of surface wind stresses and
C     barometric pressure on a rectangular grid (either in Longitude,
C     Latitude or Cartesian coordinates, consistent with the ADCIRC grid
C     coordinates) from unit 22 and interpolate them in space onto the
C     ADCIRC grid. If the ADCIRC grid is in Lon, Lat these must be in
C     radians.
C  
C     It is assumed that the meteorological grid is set up so that y
C     (e.g., latitude) varies from north (k=1) to south (k=NWLAT) and x
C     (e.g., longitude) varies from west (j=1) to east (j=NWLON).
C  
C     The spatial extents of the meteorological grid must be consistent
C     with the ADCIRC model domain.  For example, if ADCIRC uses
C     negative longitude values to indicate locations W of the Greenwich
C     meridian, the meteorological file must be similarly organized.
C     Any grid that crosses the Greenwich Meridian should be organized
C     so that the seam occurs @ 180 deg longitude. Therefore, the
C     meteorological and ADCIRC grids should use negative longitudes W
C     of the Greenwich Meridian and positive longitudes to the E.
C  
C   NOTE:  It is assumed that the met file data is oriented so that    
C          the outer loop is on latitude and the inner loop is on      
C          longitude.  For example:                                    
C          line 1             lat 1,     lon 1                         
C          line 2             lat 1,     lon 2                         
C            .                                                         
C          line nwlon         lat 1,     lon nwlon                     
C          line nwlon+1       lat 2,     lon 1                         
C          line nwlon+2       lat 2,     lon 2                         
C            .                                                         
C          line 2*nwlon       lat 2,     lon nwlon                     
C          line 2*nwlon+1     lat 3,     lon 1                         
C          line 2*nwlon+2     lat 3,     lon 2                         
C            .                                                         
C          line nwlon*nwlat   lat nwlat, lon nwlon                     
C                                                                      
C   NOTE:  It is assumed that he met file data is oriented so that     
C          latitude varies from the northern most value (lat 1) to the 
C          southern most value (lat nwlat) and longitude varies in an  
C          easterly direction (e.g. from 0 to 360 where positive       
C          longitudes are angles measured easterly of the GM.          
C
C   NOTE:  It is assumed that wind stress is in EAST, NORTH components
C          in m/s and pressure is in N/m^2 
C                                                                      
C   NOTE:  WLATMAX,WLONMIN,WLATINC,WLONINC should be in deg.           
C                                                                      
C   NOTE:  This should wrap if XCOORD > WLONMIN+NWLON*WLONINC  or      
C          XCOORD < WLONMIN                                            
C
C     ----------------------------------------------------------------
      SUBROUTINE NWS7GET(X,Y,SLAM,SFEA,WVNX,WVNY,PRESS,NP,NWLON,NWLAT,
     &     WLATMAX,WLONMIN,WLATINC,WLONINC,ICS,RHOWAT0,G)
      IMPLICIT NONE
      INTEGER, intent(in) :: NP     ! number of nodes
      REAL(8), intent(in), dimension(NP) :: X
      REAL(8), intent(in), dimension(NP) :: Y
      REAL(8), intent(in), dimension(NP) :: SLAM
      REAL(8), intent(in), dimension(NP) :: SFEA
      REAL(SZ), intent(out), dimension(NP) :: WVNX
      REAL(SZ), intent(out), dimension(NP) :: WVNY
      REAL(SZ), intent(out), dimension(NP) :: PRESS
      INTEGER, intent(in) :: NWLON  ! max. # of longitudes in stress file   
      INTEGER, intent(in) :: NWLAT  ! max. # of latitutes in stress file
      REAL(SZ), intent(in) :: WLATMAX
      REAL(SZ), intent(in) :: WLONMIN
      REAL(SZ), intent(in) :: WLATINC
      REAL(SZ), intent(in) :: WLONINC
      INTEGER, intent(in) :: ICS    ! coord.sys., 1=cartesian, 2=spherical
      REAL(SZ), intent(in):: RHOWAT0! ref. dens. of water
      REAL(SZ), intent(in):: G      ! gravitational constant

      INTEGER I                     ! node loop counter
      INTEGER K, J                  ! latitude, longitude loop counters
      REAL(SZ) RHOWATG             ! ref. dens. of water * grav. constant
      INTEGER  LATIND1,LATIND2,LONIND1,LONIND2
      REAL(SZ) WLATM,WLONM,XWRATIO,YWRATIO
      REAL(8) XCOOR,YCOOR
      REAL(SZ), SAVE, ALLOCATABLE :: WVXFN(:,:),WVYFN(:,:),PRN(:,:)
      LOGICAL, SAVE :: MemoryAllocated = .False.
C     
      IF (.not.MemoryAllocated) THEN
         ALLOCATE ( WVXFN(NWLAT,NWLON),WVYFN(NWLAT,NWLON),
     &        PRN(NWLAT,NWLON) )
         MemoryAllocated = .True.
      ENDIF
C     
      RHOWATG=RHOWAT0*G
      DO K=1,NWLAT
         DO J=1,NWLON
            READ(22,*) WVXFN(K,J),WVYFN(K,J),PRN(K,J)
         END DO
      END DO
C
      DO I=1,NP                 !INTERPOLATE TO ADCIRC GRID
         IF(ICS.EQ.2) THEN
            YCOOR=SFEA(I)*RAD2DEG
            XCOOR=SLAM(I)*RAD2DEG
         ENDIF
         IF(ICS.EQ.1) THEN
            YCOOR=Y(I)
            XCOOR=X(I)
         ENDIF
         LATIND2=(WLATMAX-YCOOR)/WLATINC + 1
         IF(LATIND2.EQ.NWLAT) LATIND2=LATIND2-1
         LATIND1=LATIND2 + 1
         LONIND1=(XCOOR-WLONMIN)/WLONINC + 1
         LONIND2=LONIND1 + 1
C     
         WLONM = WLONMIN + (LONIND1-1)*WLONINC 
         WLATM = WLATMAX - (LATIND1-1)*WLATINC
         XWRATIO=(XCOOR-WLONM)/WLONINC
         YWRATIO=(YCOOR-WLATM)/WLATINC
C     
         IF(LONIND1.EQ.0) LONIND1=NWLON
         IF(LONIND1.EQ.NWLON) LONIND2=1
C     
         WVNX(I) = WVXFN(LATIND2,LONIND2)*XWRATIO*YWRATIO
     &        + WVXFN(LATIND2,LONIND1)*(1.d0-XWRATIO)*YWRATIO
     &        + WVXFN(LATIND1,LONIND2)*XWRATIO*(1.d0-YWRATIO)
     &        + WVXFN(LATIND1,LONIND1)*(1.d0-XWRATIO)*(1.d0-YWRATIO)
         WVNY(I) = WVYFN(LATIND2,LONIND2)*XWRATIO*YWRATIO
     &        + WVYFN(LATIND2,LONIND1)*(1.d0-XWRATIO)*YWRATIO
     &        + WVYFN(LATIND1,LONIND2)*XWRATIO*(1.d0-YWRATIO)
     &        + WVYFN(LATIND1,LONIND1)*(1.d0-XWRATIO)*(1.d0-YWRATIO)
         PRESS(I) = PRN(LATIND2,LONIND2)*XWRATIO*YWRATIO
     &        + PRN(LATIND2,LONIND1)*(1.d0-XWRATIO)*YWRATIO
     &        + PRN(LATIND1,LONIND2)*XWRATIO*(1.d0-YWRATIO)
     &        + PRN(LATIND1,LONIND1)*(1.d0-XWRATIO)*(1.d0-YWRATIO)
         PRESS(I) = PRESS(I)/RHOWATG
      END DO
C
      RETURN
C     ----------------------------------------------------------------
      END SUBROUTINE NWS7GET
C     ----------------------------------------------------------------




C***********************************************************************
C                                                                      *
C   Read in and interpolate onto the ADCIRC grid wind fields from U.S. *
C   National Weather Service AVN model SFLUX meteorological files.     *
C                                                                      *
C   The input files are in binary and have been created by the GRIB    *
C   unpacking program unpkgrb1.f to extract only the U 10M, V 10M, and *
C   surface P fields.    THE BINARY INPUT HAS BEEN ELIMINATED!!!!      *
C   The input files are in ASCII and contain surface P, U 10M and V 10M*
C   fields.                                                            *
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
C                           R.L. 4/14/99                               *
C                           R.L.09/04/00 added RHOWAT0 to call         *
C   R.L. 4/2/01  changed MNWLAT,MNWLON in ALLOCATE statement to        *
C                LATB,LONB; elminiated MNWP as a dimension             *
C***********************************************************************

      SUBROUTINE NWS10GET(NWSGGWI,FLON,FLAT,ULL,VLL,PLL,NP,RHOWAT0,G,
     &     LONB,LATB,WTIMINC)
      USE SIZES
      IMPLICIT NONE
      INTEGER, SAVE :: FIRSTCALL = 0
      INTEGER N,NP,NWSGGWI,LONB,LATB,I,J,JJ,IEXT,IDIG1,IDIG2,
     &     IDIG3,KERR
      REAL*8 WTIMINC
      REAL*8 FLAT(*),FLON(*)
      REAL(SZ)  ULL(*),VLL(*),PLL(*)
      REAL(SZ) RHOWAT0,RHOWATG,G,GDLON,P1,P2,P3,P4,U1,U2,U3,U4,
     &     V1,V2,V3,V4
      INTEGER KGDS(200)
      INTEGER,SAVE,ALLOCATABLE ::  N00(:),N10(:),N11(:),N01(:)
      REAL(SZ),SAVE,ALLOCATABLE :: D00(:),D10(:),D11(:),D01(:)
      REAL(SZ),SAVE,ALLOCATABLE :: COLRAB(:),DUMMY(:),
     &     GCLAT(:),GCLON(:)
      REAL(SZ),SAVE,ALLOCATABLE ::  UG(:),VG(:),PG(:)
      CHARACTER*1 PDS(50),FNAME2(8)
      CHARACTER*8 FNAME1
      EQUIVALENCE (FNAME1,FNAME2)
      LOGICAL FOUND
C     
      IF (FIRSTCALL.EQ.0) THEN
         FIRSTCALL = 1
         ALLOCATE ( UG(LATB*LONB),VG(LATB*LONB),
     &        PG(LATB*LONB) )
         ALLOCATE ( N00(MNP),N10(MNP),N11(MNP),N01(MNP) )
         ALLOCATE ( D00(MNP),D10(MNP),D11(MNP),D01(MNP) )
         ALLOCATE ( COLRAB(LATB),DUMMY(LATB),GCLAT(LATB),
     &        GCLON(LONB) )
      ENDIF
C     
      RHOWATG=RHOWAT0*G
C     
C...  The first time the subroutine is called, setup the Gaussian grid and
C...  determine the interpolating factors for the ADCIRC grid.
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
     &        D01,NP,LONB,LATB)
         RETURN
      ENDIF

C...  Figure out the data file name

      FNAME1='fort.   '
      IEXT=200 + NWSGGWI*(WTIMINC/3600)
      IDIG1=IEXT/100
      IDIG2=(IEXT-100*IDIG1)/10
      IDIG3=(IEXT-100*IDIG1-10*IDIG2)
      FNAME2(6)=CHAR(IDIG1+48)
      FNAME2(7)=CHAR(IDIG2+48)
      FNAME2(8)=CHAR(IDIG3+48)


C...  Enter, locate and open the data file

 1010 FORMAT(' File ',A8,' WAS NOT FOUND!  FATAL ERROR',/)
 1011 FORMAT(' File ',A8,' WAS FOUND!  Opening & Processing file',/)

      if (myproc == 0) WRITE(screenunit,*) '  '
      INQUIRE(FILE=FNAME1,EXIST=FOUND)
      IF(FOUND) GOTO 32
      if (myproc == 0) WRITE(screenunit,1010) FNAME1
      WRITE(16,1010) FNAME1
#ifdef CMPI
      call msg_fini()
#endif
      STOP
 32   WRITE(screenunit,1011) FNAME1

C...Open and read the GRIB BINARY data file
c     OPEN(IEXT,FILE=FNAME1,status='old',access='sequential',
c    &     form='unformatted',iostat=kerr)
c     READ(IEXT,END=1100) LENPDS,LENKGDS,NWORDS
c     IF(LENPDS.GT.0) READ(IEXT,END=1100) (pds(j),j=1,lenpds)
c     IF(LENKGDS.GT.0) READ(IEXT,END=1100) (kgds(j),j=1,lenkgds)
c     IF(NWORDS.GT.0) READ(IEXT,END=1100) (UG(J),J=1,NWORDS)
c
c     READ(IEXT,END=1100) LENPDS,LENKGDS,NWORDS
c     IF(LENPDS.GT.0) READ(IEXT,END=1100) (pds(j),j=1,lenpds)
c     IF(LENKGDS.GT.0) READ(IEXT,END=1100) (kgds(j),j=1,lenkgds)
c     IF(NWORDS.GT.0) READ(IEXT,END=1100) (VG(J),J=1,NWORDS)
c
c     READ(IEXT,END=1100) LENPDS,LENKGDS,NWORDS
c     IF(LENPDS.GT.0) READ(IEXT,END=1100) (pds(j),j=1,lenpds)
c     IF(LENKGDS.GT.0) READ(IEXT,END=1100) (kgds(j),j=1,lenkgds)
c     IF(NWORDS.GT.0) READ(IEXT,END=1100) (PG(J),J=1,NWORDS)

C...Open and read the ASCII data file

      OPEN(IEXT,FILE=FNAME1,status='old',iostat=kerr)
      DO I=1,LONB*LATB
         READ(IEXT,*) PG(I),UG(I),VG(I)
      ENDDO

 1100 CLOSE(IEXT)


C.....Go from the Gaussian grid to the ADCIRC grid
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
C     
      RETURN
      END SUBROUTINE


C***********************************************************************
C  Subroutine to compute the latutudes in a Global Gaussian Lat/Lon    *
C  grid with T126 resolution (GRIB Grid type 126).                     *
C                                                                      *
C       modified from the original GLATS by R.L. 4/24/96               *
C***********************************************************************

      SUBROUTINE GLATS(LGGHAF,COLRAD,WGT,WGTCS,RCS2)
      USE SIZES
      IMPLICIT NONE
      REAL(SZ) COLRAD(*),WGT(*),WGTCS(*),RCS2(*)
      INTEGER LGGHAF,L2,K,K1,ITER
      REAL(SZ) SI,SCALE,RL2,DRAD,DRADZ,RAD,P1,P2,EPS,PHI,X,W,SN,RC
C     
      EPS=1.d-6
C     EPS=1.d-12
C     PRINT 101
C     101  FORMAT ('0 I   COLAT   COLRAD     WGT', 12X, 'WGTCS',
CCCC  1 10X, 'ITER  RES')
C     
      SI = 1.0d0
      L2=2*LGGHAF
      RL2=L2
      SCALE = 2.0d0/(RL2*RL2)
      K1=L2-1
      DRADZ = PI / 360.d0
      RAD = 0.0
      DO 1000 K=1,LGGHAF
         ITER=0
         DRAD=DRADZ
 1       CALL POLY(L2,RAD,P2)
 2       P1 =P2
         ITER=ITER+1
         RAD=RAD+DRAD
         CALL POLY(L2,RAD,P2)
         IF(SIGN(SI,P1).EQ.SIGN(SI,P2)) GO TO 2
         IF(DRAD.LT.EPS)GO TO 3
         RAD=RAD-DRAD
         DRAD = DRAD * 0.25d0
         GO TO 1
 3       CONTINUE
         COLRAD(K)=RAD
         PHI = RAD * 180.d0 / PI
         CALL POLY(K1,RAD,P1)
         X = COS(RAD)
         W = SCALE * (1.0d0 - X*X)/ (P1*P1)
         WGT(K) = W
         SN = SIN(RAD)
         W=W/(SN*SN)
         WGTCS(K) = W
         RC=1.d0/(SN*SN)
         RCS2(K) = RC
         CALL POLY(L2,RAD,P1)
C     PRINT 102,K,PHI,COLRAD(K),WGT(K),WGTCS(K),ITER,P1
C     102  FORMAT(1H ,I2,2X,F6.2,2X,F10.7,2X,E13.7,2X,E13.7,2X,I4,2X,D13.7)
 1000 CONTINUE
c     PRINT 100,LGGHAF
c     100  FORMAT(1H ,'SHALOM FROM 0.0 GLATS FOR ',I3)
      RETURN
      END SUBROUTINE


C***********************************************************************
C  Subroutine used by GLATS.                                           *
C***********************************************************************

      SUBROUTINE POLY(N,RAD,P)
      USE SIZES
      IMPLICIT NONE
      INTEGER N,I
      REAL(SZ) RAD,P,X,Y1,Y2,Y3,G
C     
      X = COS(RAD)
      Y1 = 1.0d0
      Y2=X
      DO 1 I=2,N
         G=X*Y2
         Y3=G-Y1+G-(G-Y1)/FLOAT(I)
         Y1=Y2
         Y2=Y3
 1    CONTINUE
      P=Y3
      RETURN
      END SUBROUTINE

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
     &     D11,D01,NP,LONB,LATB)
      USE SIZES
      IMPLICIT NONE
      INTEGER,SAVE :: ICALL = 0
      INTEGER NP,N,I,LONB,LATB,NLAT,NLON,LON,LONP1,LAT,LATP1
      REAL*8 DLAT,DLON,FLONWORK,COLAT,DDLAT,XLAT,DFLAT,DFLAT1,
     &     DDLON,XLON,DFLON,DFLON1
      REAL*8 FLAT(*),FLON(*)
      REAL(SZ) GCLAT(*),GCLON(*)
      INTEGER  N00(*),N10(*),N11(*),N01(*)
      REAL(SZ) D00(*),D10(*),D11(*),D01(*)
C     
      IF( ICALL .EQ. 0 ) THEN
         ICALL = 1
c       PRINT 1234
c1234   FORMAT(' = IN ROUTINE G2RINI FOR HORIZONTAL INTERPOLATION = ')

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
            IF(LON.EQ.NLON) LONP1=1 !Circle condition
            IF((LON.LT.1).OR.(LON.GT.NLON)) THEN
               PRINT *,' ***** ERROR IN LON ****'
               PRINT *,' I ',I
               PRINT *,' LON ',LON
               PRINT *,' DLON ',DLON
               PRINT *,' FLON ',FLON(I)
#ifdef CMPI
               call msg_fini()
#endif
               STOP
            ENDIF  
            
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
#ifdef CMPI
               call msg_fini()
#endif
               STOP
            ENDIF

 5          CONTINUE
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
            
 9          CONTINUE
            DDLAT=GCLAT(LATP1)-GCLAT(LAT)
            XLAT=GCLAT(LAT)
            DFLAT1=(COLAT-XLAT)/DDLAT
            IF(LAT.EQ.1) DFLAT1=MAX(0.d0,DFLAT1) !MODIFY THIS FOR POLAR POINTS
            IF(LATP1.EQ.NLAT) DFLAT1=MIN(1.d0,DFLAT1) !MODIFY THIS FOR POLAR POINTS
            DFLAT=1.d0-DFLAT1
            DDLON=DLON
            XLON=GCLON(LON)
            DFLON1=(FLONWORK-XLON)/DDLON
            DFLON=1.d0-DFLON1
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
c        if (myproc == 0) 
c    &     WRITE(screenunit,*) ' D00 TO D11 SHOULD BE ALL POSITIVE.'
         
      ELSE
c        if (myproc == 0) 
c          WRITE(screenunit,*) ' G2RINI ALREADY CALLED '
      ENDIF
      
      RETURN
      END SUBROUTINE


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
C   R.L.09/04/00 added RHOWAT0 to call                                 *
C   R.L. 4/02/01  elminiated MNWP as a dimension                       *
C   R.L. 9/14/01  changed MNWLAT,MNWLON in ALLOCATE statement to       *
C                271,181                                               *
C***********************************************************************

      SUBROUTINE NWS11GET(NWSEGWI,IDSETFLG,FLON,FLAT,ULL,VLL,PLL,NP,
     &     RHOWAT0,G)
      USE SIZES
      IMPLICIT NONE
      INTEGER,SAVE  ::  ICALL = 0
      INTEGER NWSEGWI,IDSETFLG,NP,I,IEXT,IDIG1,IDIG2,IDIG3,KERR,N
      INTEGER IYEAR,IMONTH,IDAY,IHOUR
      REAL*8 RHOWATG100,FLONDEG,FLATDEG
      REAL(SZ) P1,P2,P3,U1,U2,U3,V1,V2,V3,UE29,VE29,CBETAU,SBETAU,G
      REAL(SZ) RHOWAT0
      REAL(SZ) ULL(*),VLL(*),PLL(*)
      REAL*8 FLAT(*),FLON(*)
C     
      INTEGER,SAVE,ALLOCATABLE ::  N1(:),N2(:),N3(:)
      REAL(SZ),SAVE,ALLOCATABLE :: D1(:),D2(:),D3(:),BETAU(:)
      REAL(SZ),SAVE,ALLOCATABLE :: UE(:),VE(:),PE(:)
C     
      CHARACTER*1 FNAME2(8)
      CHARACTER*8 FNAME1
      EQUIVALENCE (FNAME1,FNAME2)
      LOGICAL FOUND
C     
      IF (ICALL.EQ.0) THEN
         ICALL = 1
         ALLOCATE ( N1(MNP),N2(MNP),N3(MNP) )
         ALLOCATE ( D1(MNP),D2(MNP),D3(MNP),BETAU(MNP) )
         ALLOCATE ( UE(181*271),VE(181*271),PE(181*271) )
      ENDIF
C     
      RHOWATG100=RHOWAT0*G*100.d0

C...  The first time the subroutine is called, setup the interpolating factors
C...  between the Eta-29 grid and the ADCIRC grid.

      IF((NWSEGWI.EQ.0).AND.(IDSETFLG.EQ.0)) THEN
         if (myproc == 0) then
          WRITE(screenunit,*) 'Computing ETA29 met field interp factors'
         endif
         DO I=1,NP
            flondeg=rad2deg*flon(i)
            flatdeg=rad2deg*flat(i)
            CALL E29SEARCH(I,FLONDEG,FLATDEG,N1(I),N2(I),N3(I),
     &           D1(I),D2(I),D3(I),betau(i))
         END DO
         RETURN
      ENDIF

C...  Figure out the met data file name

      FNAME1='fort.   '
      IEXT=200 + NWSEGWI
      IDIG1=IEXT/100
      IDIG2=(IEXT-100*IDIG1)/10
      IDIG3=(IEXT-100*IDIG1-10*IDIG2)
      FNAME2(6)=CHAR(IDIG1+48)
      FNAME2(7)=CHAR(IDIG2+48)
      FNAME2(8)=CHAR(IDIG3+48)

C...  If appropriate, enter, locate and open the met data file

 1010 FORMAT(' File ',A8,' WAS NOT FOUND!  FATAL ERROR',/)
 1011 FORMAT(' File ',A8,' WAS FOUND!  Opening & Processing file',/)

      if (myproc == 0) WRITE(screenunit,*) '  '
      INQUIRE(FILE=FNAME1,EXIST=FOUND)
      IF(FOUND) GOTO 32
      if (myproc == 0) WRITE(screenunit,1010) FNAME1
      WRITE(16,1010) FNAME1
#ifdef CMPI
      call msg_fini()
#endif
      STOP
 32   if (myproc == 0) WRITE(screenunit,1011) FNAME1
      IF((NWSEGWI.EQ.0).OR.(IDSETFLG.EQ.1)) OPEN(IEXT,FILE=FNAME1,
     &status='old',access='sequential',form='unformatted',iostat=kerr)

C...  Read the met data file

      READ(IEXT,END=1100) IYEAR,IMONTH,IDAY,IHOUR
      READ(IEXT,END=1100) UE,VE,PE

      IF(NWSEGWI.EQ.0) THEN     !If the first file, read until the end
         DO I=2,IDSETFLG
            READ(IEXT,END=1100) IYEAR,IMONTH,IDAY,IHOUR
            READ(IEXT,END=1100) UE,VE,PE
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
      END SUBROUTINE



C***********************************************************************
C  Subroutine to find where a given lon,lat falls in the Eta29 grid,   *
C     determine the interpolating factors to interpolate Eta29 fields  *
C     to that position, and finally to compute the angle to rotate the *
C     Eta29 velocity field to get to a lon, lat coordinated system.    *
C                                                                      *
C                    Written by R.L.       1/12/98                     *
C***********************************************************************

      subroutine e29search(node,FLON,FLAT,NN1,NN2,NN3,DD1,DD2,DD3,betau)
      implicit none
      integer nn1,nn2,nn3,node,icode,nwlon,nwlat,ifflag
      integer i,j,im2,jm2,n,ia,ja,na,ib,jb,nb,ic,jc,nc,id,jd,nd,
     &  ie,je,ne,ig,jg,ng,if
      real(sz) dd1,dd2,dd3,betau,ri,x1,x2,x3,x4,y1,y2,y3,y4
      real(sz) aemin,areas,a1,a2,a3,aa,ae,lambda
      real(8) lamda0,phi0,rphi0,cphi0,sphi0,tphi0,dlamda,dphi,rdlamda,
     &       rdphi,rflat,tflat,sflat,cflat,a,rlamar,cphiicrlamda,phiarg,
     &       rphii,rlamda,ri1,ri2,rj,dgtora,flon,flat
      real(sz) lamda,lamdaa,lamdab,lamdac,lamdad,lamdae,lamdag
      real(sz) phi,phia,phib,phic,phid,phie,phig
c
      icode=0
      nwlon=181
      nwlat=271
      dgtora=deg2rad
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
c
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
c
      ri2=(rlamda/rdlamda+nwlon+1)/2.
      ri1=(rlamda/rdlamda+nwlon)/2.
      rj=rphii/rdphi+(nwlat+1)/2
      j=(rj+0.5d0)
      ri=ri1
      if(mod(j,2).eq.0) ri=ri2
      i=(ri+0.5d0)

c     if (myproc == 0) then
c       write(screenunit,*) "lamda, phi = ",flon,flat
c       write(screenunit,*) "ri1, ri2, ri, rj = ",ri1,ri2,ri,rj
c       write(screenunit,*) "i, j = ",i,j
c     endif

      if ((rj.lt.1).or.(rj.gt.nwlat)) then
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

      if (mod(j,2).eq.0) then
         if ((ri.lt.1).or.(ri.gt.(nwlon+0.5d0))) then
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
      
      if (mod(j,2).ne.0) then
         if ((ri.lt.0.5).or.(ri.gt.nwlon)) then
c           write(333,*) 'ADCIRC grid node ',node,
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

      if ((i.eq.1).and.(mod(j,2).eq.0)) then
         ia=i
         ja=j-2
      else
         ia=i
         if(mod(j,2).eq.0) ia=i-1
         ja=j-1
      endif
c                                 this neighbor lies outside of Eta29 grid
      if ((ia.lt.1).or.(ja.lt.1)) then
         na=0
      else
         call e29calc(ia,ja,lamdaa,phia,na)
      endif

c     compute the coordinates of neighbor node "b" (located W of closest node)

      ib=i-1
      jb=j
      if (ib.lt.1) then         !this neighbor lies outside of Eta29 grid
         nb=0
      else
         call e29calc(ib,jb,lamdab,phib,nb)
      endif

c     compute the coordinates of neighbor node "c" (located NW of closest node)

      if ((i.eq.1).and.(mod(j,2).eq.0)) then
         ic=i
         jc=j+2
      else
         ic=ia
         jc=j+1
      endif
c                                    this neighbor lies outside of Eta29 grid
      if ((ic.lt.1).or.(jc.gt.nwlat)) then  
         nc=0
      else
         call e29calc(ic,jc,lamdac,phic,nc)
      endif

c     compute the coordinates of neighbor node "d" (located NE of closest node)

      if ((i.eq.181).and.(mod(j,2).ne.0)) then
         id=i
         jd=j+2
      else
         id=ic+1
         jd=j+1
      endif
c                                    this neighbor lies outside of Eta29 grid
      if ((id.gt.nwlon).or.(jd.gt.nwlat)) then  
         nd=0
      else
         call e29calc(id,jd,lamdad,phid,nd)
      endif

c     compute the coordinates of neighbor node "e" (located E of closest node)

      ie=i+1
      je=j
      if (ie.gt.nwlon) then     !this neighbor lies outside of Eta29 grid
         ne=0
      else
         call e29calc(ie,je,lamdae,phie,ne)
      endif
      
c     compute the coordinates of neighbor node "g" (located SE of closest node)

      if ((i.eq.181).and.(mod(j,2).ne.0)) then
         ig=i
         jg=j-2
      else
         ig=id
         jg=j-1
      endif
c                                    this neighbor lies outside of Eta29 grid
      if ((ig.gt.nwlon).or.(jg.lt.1)) then  
         ng=0
      else
         call e29calc(ig,jg,lamdag,phig,ng)
      endif

c     if (myproc == 0) then
c      write(screenunit,*) 'closest E29 node i,j = ',n,i,j,lamda,phi
c      if(na.eq.0) write(screenunit,*) 'point a falls outside of Eta29 grid'
c      if(na.ne.0) write(screenunit,*) 'point a   = ',na,ia,ja,lamdaa,phia
c      if(nb.eq.0) write(screenunit,*) 'point b falls outside of Eta29 grid'
c      if(nb.ne.0) write(screenunit,*) 'point b   = ',nb,ib,jb,lamdab,phib
c      if(nc.eq.0) write(screenunit,*) 'point c falls outside of Eta29 grid'
c      if(nc.ne.0) write(screenunit,*) "point c   = ",nc,ic,jc,lamdac,phic
c      if(nd.eq.0) write(screenunit,*) 'point d falls outside of Eta29 grid'
c      if(nd.ne.0) write(screenunit,*) "point d   = ",nd,id,jd,lamdad,phid
c      if(ne.eq.0) write(screenunit,*) 'point e falls outside of Eta29 grid'
c      if(ne.ne.0) write(screenunit,*) "point e   = ",ne,ie,je,lamdae,phie
c      if(ng.eq.0) write(screenunit,*) 'point g falls outside of Eta29 grid'
c      if(ng.ne.0) write(screenunit,*) "point g   = ",ng,ig,jg,lamdag,phig
c     endif

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
      AEMIN=99999.d0

c     test if the point is in triangle ij - b - a

      if ((na.ne.0).and.(nb.ne.0)) then
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
c     write(333,*) "AE = ",AE
         IF((AE.LT.1.0d-5).AND.(AE.LT.AEMIN)) THEN
            AEMIN=AE
            NN1=n
            NN2=nb
            NN3=na
            DD1=((X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4))/AREAS
            DD2=((X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1))/AREAS
            DD3=(-(X4-X1)*(Y2-Y1)+(Y4-Y1)*(X2-X1))/AREAS
            call betaucalc(i,j,DD1,ib,jb,DD2,ia,ja,DD3,betau)
            ifflag=ifflag+1
c     write(333,*) 'position found in triangle ij - b - a'
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
c     write(333,*) "AE = ",AE
            IF((AE.LT.1.0d-5).AND.(AE.LT.AEMIN)) THEN
               NN1=n
               NN2=nc
               NN3=na
               DD1=((X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4))/AREAS
               DD2=((X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1))/AREAS
               DD3=(-(X4-X1)*(Y2-Y1)+(Y4-Y1)*(X2-X1))/AREAS
               call betaucalc(i,j,DD1,ic,jc,DD2,ia,ja,DD3,betau)
               ifflag=ifflag+1
c     write(333,*) 'position found in triangle ij - c - a'
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
c     write(333,*) "AE = ",AE
         IF((AE.LT.1.0d-5).AND.(AE.LT.AEMIN)) THEN
            NN1=n
            NN2=nc
            NN3=nb
            DD1=((X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4))/AREAS
            DD2=((X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1))/AREAS
            DD3=(-(X4-X1)*(Y2-Y1)+(Y4-Y1)*(X2-X1))/AREAS
            call betaucalc(i,j,DD1,ic,jc,DD2,ib,jb,DD3,betau)
            ifflag=ifflag+1
c     write(333,*) 'position found in triangle ij - c - b'
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
c     write(333,*) "AE = ",AE
         IF((AE.LT.1.0d-5).AND.(AE.LT.AEMIN)) THEN
            NN1=n
            NN2=nd
            NN3=nc
            DD1=((X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4))/AREAS
            DD2=((X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1))/AREAS
            DD3=(-(X4-X1)*(Y2-Y1)+(Y4-Y1)*(X2-X1))/AREAS
            call betaucalc(i,j,DD1,id,jd,DD2,ic,jc,DD3,betau)
            ifflag=ifflag+1
c     write(333,*) 'position found in triangle ij - d - c'
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
c     write(333,*) "AE = ",AE
            IF((AE.LT.1.0d-5).AND.(AE.LT.AEMIN)) THEN
               NN1=n
               NN2=ng
               NN3=nd
               DD1=((X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4))/AREAS
               DD2=((X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1))/AREAS
               DD3=(-(X4-X1)*(Y2-Y1)+(Y4-Y1)*(X2-X1))/AREAS
               call betaucalc(i,j,DD1,ig,jg,DD2,id,jd,DD3,betau)
               ifflag=ifflag+1
c     write(333,*) 'position found in triangle ij - g - d'
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
c     write(333,*) "AE = ",AE
         IF((AE.LT.1.0d-5).AND.(AE.LT.AEMIN)) THEN
            NN1=n
            NN2=ne
            NN3=nd
            DD1=((X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4))/AREAS
            DD2=((X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1))/AREAS
            DD3=(-(X4-X1)*(Y2-Y1)+(Y4-Y1)*(X2-X1))/AREAS
            call betaucalc(i,j,DD1,ie,je,DD2,id,jd,DD3,betau)
            ifflag=ifflag+1
c     write(333,*) 'position found in triangle ij - e - d'
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
c     write(333,*) "AE = ",AE
         IF((AE.LT.1.0d-5).AND.(AE.LT.AEMIN)) THEN
            NN1=n
            NN2=ng
            NN3=ne
            DD1=((X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4))/AREAS
            DD2=((X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1))/AREAS
            DD3=(-(X4-X1)*(Y2-Y1)+(Y4-Y1)*(X2-X1))/AREAS
            call betaucalc(i,j,DD1,ig,jg,DD2,ie,je,DD3,betau)
            ifflag=ifflag+1
c     write(333,*) 'position found in triangle ij - g - e'
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
c     write(333,*) "AE = ",AE
         IF((AE.LT.1.0d-5).AND.(AE.LT.AEMIN)) THEN
            NN1=n
            NN2=na
            NN3=ng
            DD1=((X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4))/AREAS
            DD2=((X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1))/AREAS
            DD3=(-(X4-X1)*(Y2-Y1)+(Y4-Y1)*(X2-X1))/AREAS
            call betaucalc(i,j,DD1,ia,ja,DD2,ig,jg,DD3,betau)
            ifflag=ifflag+1
c     write(333,*) 'position found in triangle ij - a - g'
         ENDIF
      endif

c      if(ifflag.eq.0) then
c         write(333,*) 'position not found'
c        if (myproc == 0) then
c         write(screenunit,*) 'position not found in subroutine E29SEARCH'
c        endif
c        icode=3
c      else
c        if (myproc == 0) then
c         write(screenunit,*) 'i,j,NN1,NN2,NN3,DD1,DD2,DD3'
c        endif
c         write(333,999) i,j,NN1,NN2,NN3,DD1,DD2,DD3,betau/dgtora
c 999     format(5I8,1x,3E13.6)
c      endif

      return
      end subroutine



C***********************************************************************
C  Subroutine to compute the longititude and latitude of a given i,j   *
C       position in the Eta29 grid.                                    *
C                                                                      *
C                    Written by R.L.       1/11/98                     *
C***********************************************************************

      subroutine e29calc(i,j,lamda,phi,n)
      implicit none
      integer i,j,n,nwlon,nwlat,im2,jm2,i1,i2,i1p1,i1m1,i2p1,i2m1,
     &     i3p1,i3m1
      real(sz) lamda,phi,phii,dlon,dlat,dlnt,arg,betau1,betau2,betau3
      real(8) lamda0,phi0,rphi0,cphi0,sphi0,tphi0,dlamda,dphi,rdlamda,
     &     rdphi,a,rlamar,phiarg,rlamda,dgtora
c     
      nwlon=181
      nwlat=271
      dgtora=deg2rad
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
c     
      jm2=(nwlat+1)/2
      im2=nwlon*2
c     
      phii=rdphi*float(j-jm2)
      i1=2*i-1
      i2=2*i
      if(mod(j,2).ne.0) then
         rlamda=rdlamda*float(i2-nwlon)
      else
         rlamda=rdlamda*float(i1-nwlon)
      endif
      phiarg= sin(phii)*cphi0+cos(phii)*sphi0*cos(rlamda)
      if(phiarg.gt.1.0d0) phiarg=1.0d0
      if(phiarg.lt.-1.0d0) phiarg=-1.0d0
      phi=asin(phiarg)
      rlamar= cos(phii)*cos(rlamda)/(cos(phi)*cphi0)-tan(phi)*tphi0
      if(rlamar.gt.1.0d0) rlamar=1.0d0
      if(rlamar.lt.-1.d0) rlamar=-1.d0
      a=acos(rlamar)/dgtora
      if(rlamda.le.0.) then
         lamda=lamda0-a
      else
         lamda=lamda0+a
      endif
      phi=phi/dgtora
      n=nwlon*(j-1)+i
C     
      return
      end subroutine


C***********************************************************************
C  Subroutine to compute the conversion angle between the E29 velocity *
C       field and a lon,lat coordinate system.                         *
C                                                                      *
C                    Written by R.L.       1/12/98                     *
C***********************************************************************

      subroutine betaucalc(i1,j1,dd1,i2,j2,dd2,i3,j3,dd3,betau)
      implicit none
      integer i1,j1,i2,j2,i3,j3,n,i1p1,i1m1,i2p1,i2m1,i3p1,i3m1
      real(sz) dd1,dd2,dd3,betau
      real(sz) lamda,lamdap1,lamdam1,phi,phip1,phim1,dlon,dlat,
     &     dlnt,arg,betau1,betau2,betau3,dgtora
c     
      dgtora=deg2rad
c     
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
      if(arg.gt.1.d0) arg=1.d0
      if(arg.lt.-1.d0) arg=-1.d0
      betau1=asin(arg)
c     
      if(i2.ne.181) then
         i2p1=i2+1
      else
         i2p1=i2
      endif
c     
      if(i2.ne.1) then
         i2m1=i2-1
      else
         i2m1=i2
      endif
c     
      call e29calc(i2,j2,lamda,phi,n)
      call e29calc(i2p1,j2,lamdap1,phip1,n)
      call e29calc(i2m1,j2,lamdam1,phim1,n)
      dlon=(lamdap1-lamdam1)*cos(phi*dgtora)
      dlat=phip1-phim1
      dlnt=sqrt(dlon*dlon+dlat*dlat)
      arg=dlat/dlnt
      if(arg.gt.1.d0) arg=1.d0
      if(arg.lt.-1.d0) arg=-1.d0
      betau2=asin(arg)
c     
      if(i3.ne.181) then
         i3p1=i3+1
      else
         i3p1=i3
      endif
c     
      if(i3.ne.1) then
         i3m1=i3-1
      else
         i3m1=i3
      endif
c     
      call e29calc(i3,j3,lamda,phi,n)
      call e29calc(i3p1,j3,lamdap1,phip1,n)
      call e29calc(i3m1,j3,lamdam1,phim1,n)
      dlon=(lamdap1-lamdam1)*cos(phi*dgtora)
      dlat=phip1-phim1
      dlnt=sqrt(dlon*dlon+dlat*dlat)
      arg=dlat/dlnt
      if(arg.gt.1.d0) arg=1.d0
      if(arg.lt.-1.d0) arg=-1.d0
      betau3=asin(arg)
      betau=dd1*betau1+dd2*betau2+dd3*betau3
C     
      return
      end subroutine

C     ----------------------------------------------------------------
C      S U B R O U T I N E   H O L L A N D  G E T
C     ----------------------------------------------------------------
C
C     jgf46.05 Subroutine to calculate wind velocity at nodes from
C     the Holland Wind model.
C
C     The format statement takes into account whether the track data is
C     hindcast/nowcast (BEST) or forecast (OFCL).
C
C     The first line in the file MUST be a hindcast, since the central
C     pressure and the RMW are carried forward from hindcasts into
C     forecasts. So there needs to be at least one hindcast to carry the data
C     forward.
C
C     Assumes spherical coordinates (ICS=2 in fort.15 file).
C
C     Based on bob's NWS67GET (See below).
C
C     ----------------------------------------------------------------
      SUBROUTINE HollandGet(X,Y,SLAM,SFEA,WVNX,WVNY,PRESS,NP,
     &     ICS,RHOWAT0,G,TIME,NSCREEN,ScreenUnit)
      USE SIZES, ONLY : MyProc, LOCALDIR
      IMPLICIT NONE
      INTEGER, intent(in) :: NP,ICS
      REAL(SZ), intent(in) :: RHOWAT0,G,TIME
      REAL(8), intent(in), dimension(NP) :: X, Y, SLAM, SFEA
      REAL(SZ), intent(out), dimension(NP) :: WVNX,WVNY
      REAL(SZ), intent(out), dimension(NP) :: PRESS
      INTEGER, intent(in) :: NSCREEN
      INTEGER, intent(in) :: ScreenUnit
      INTEGER I, J
      REAL(SZ),SAVE,ALLOCATABLE :: RAD(:),DX(:),DY(:)
      REAL(SZ),SAVE,ALLOCATABLE :: XCOOR(:),YCOOR(:) !lat and lon of nodes
      REAL(SZ),SAVE,ALLOCATABLE :: V_r(:)            !wind sp. as fn. of dist
      REAL(SZ),SAVE,ALLOCATABLE :: THETA(:)
      REAL(SZ) :: TVX,TVY,RRP,RMW,A,B,WTRATIO
      REAL(SZ) :: TransSpdX, TransSpdY
      REAL(SZ) :: TM, cpress, lon, lat, spd
      REAL(SZ) :: ts        ! storm translation speed, m/s
      REAL(SZ) :: WindMultiplier         ! for storm 2 in LPFS ensemble
      REAL(SZ) :: centralPressureDeficit ! difference btw ambient and cpress
C
      REAL(8), PARAMETER :: Ambient_Pressure=101300.d0      ! Pascals!!
C      REAL(8), PARAMETER :: Coriolis=7.287588493541447D-05  ! 1/sec
      REAL(8), PARAMETER :: Rho_Air = 1.15d0                ! kg/m/m/m
      REAL(8), PARAMETER :: e=2.718281828459045d0
      REAL(8), PARAMETER :: SamplingTimeAdj=0.88d0 ! from 1 min to 10 min winds
      REAL(8), PARAMETER :: re=6378206.4d0
      REAL(8) :: pi,omega, coriolis
      REAL(8) :: mperdeg,DEG2RAD

      LOGICAL, SAVE :: FIRSTCALL = .True.
C
              pi = 4.0 * atan(1.0)
        mperdeg  = re * pi / 180.0
         DEG2RAD = pi / 180.0
           omega = 2.0d0*pi / 86164.2
C
      IF (FIRSTCALL) THEN
         FIRSTCALL = .False.
         ALLOCATE (RAD(NP),DX(NP),DY(NP),XCOOR(NP),YCOOR(NP))
         ALLOCATE (V_r(NP),THETA(NP))
C
C     The subroutine only works for ICS=2 (spherical coordinates)
         DO I=1,NP
            XCOOR(I)=SLAM(I)*RAD2DEG
            YCOOR(I)=SFEA(I)*RAD2DEG
         END DO
      ENDIF
C
C     Get data for this time step.
      CALL GetHollandStormData(lat,lon,cpress,spd,rrp,rmw,tvx,tvy,time,
     &     nscreen,screenunit)
C
C     Calculate and limit central pressure deficit; some track files
C     (e.g., Charley 2004) may have a central pressure greater than the
C     ambient pressure that this subroutine assumes
      centralPressureDeficit = Ambient_pressure - cpress
      IF ( centralPressureDeficit .lt. 100.d0 ) THEN
         centralPressureDeficit = 100.d0
      ENDIF
C
C     jgf46.29 Subtract the translational speed of the storm from the
C     observed max wind speed to avoid distortion in the Holland curve
C     fit. The translational speed will be added back later.
      ts=sqrt(tvx*tvx+tvy*tvy)
      spd=spd-ts
C
C     Convert wind speed from 10 meter altitude (which is what the
C     NHC forecast contains) to wind speed at the top of the atmospheric
C     boundary layer (which is what the Holland curve fit requires).
      spd=spd/BLAdj
C
C     Calculate Holland parameters and limit the result to its appropriate
C     range.
      B = Rho_Air*e*(spd**2.d0)/(centralPressureDeficit)
      IF (B.lt.1.0d0) B=1.0d0
      IF (B.gt.2.5d0) B=2.5d0
C
C     Calculate Holland A parameter. (jgf46.32jgf4 commented out)
C      A = (RMW*1000.d0)**B
C
#ifdef DEBUG_HOLLAND
      WRITE(16,4321) B
 4321 FORMAT(/,2x,'Holland B parameter is ',e16.8)
#endif
C     jgf46.28 If we are running storm 2 in the Lake Pontchartrain
C     Forecast System ensemble, the final wind speeds should be
C     multiplied by 1.2.
      IF (StormNumber.eq.2) THEN
         WindMultiplier=1.2d0
      ELSE
         WindMultiplier=1.0d0
      ENDIF
C
C     Calculate wind velocity and pressure at each node.
      DO I=1,NP

C         DX(I)=XCOOR(I)-lon
C         DY(I)=YCOOR(I)-lat

        DX(I)=(XCOOR(I)-lon)*DEG2RAD
        DY(I)=(YCOOR(I)-lat)*DEG2RAD

         THETA(I)=ATAN2(DY(I),DX(I))
C   RJW v48.45
C     compute the distances based on haversine formula for distance along a sphere
         Rad(i)=re*(2.0d0*ASIN(sqrt(sin(DY(I)/2.0d0)**2.0d0+
     & cos(lat*DEG2RAD)*cos(YCOOR(i)*DEG2RAD)*sin(DX(I)/2.0d0)**2.0d0)))
C calculate the coriolis at YCOOR
         coriolis = 2.0d0 * omega * sin(YCOOR(I)*DEG2RAD)

C        RAD(I)=SQRT(DX(I)*DX(I)+DY(I)*DY(I))*100.d0*1000.d0 ! into meters!!
C
C         PRESS(I)=(cpress+(Ambient_Pressure-cpress)*EXP(-A/RAD(I)**B))
C    &        / (RHOWAT0*G) ! jgf46.32jgf4 commented out
         PRESS(I)=(cpress+(centralPressureDeficit)*
     &        EXP(-(RMW*1000.d0/RAD(I))**B)) / (RHOWAT0*G)
C         V_r(I) = sqrt(
C     &        A*B*(Ambient_Pressure-cpress)*EXP(-A/RAD(I)**B)
C     &        / ( Rho_Air*RAD(I)**B )
C     &        + (RAD(I)**2.d0)*(CORIOLIS**2.d0)/4.d0
C     &        )
C     &        - RAD(I)*CORIOLIS/2.d0 ! jgf46.32jgf4 commented out
         V_r(I) = sqrt(
     &        (RMW*1000.d0/RAD(I))**B *
     &        EXP(1.d0-(RMW*1000.d0/RAD(I))**B)*spd**2.d0
     &        + (RAD(I)**2.d0)*(CORIOLIS**2.d0)/4.d0
     &        )
     &        - RAD(I)*CORIOLIS/2.d0
C
C     jgf46.31 Determine translation speed that should be added to final
C     storm wind speed. This is tapered to zero as the storm wind tapers
C     to zero toward the eye of the storm and at long distances from the
C     storm.
         TransSpdX = (abs(V_r(I))/spd)*TVX
         TransSpdY = (abs(V_r(I))/spd)*TVY
C     Apply mutliplier for Storm2 in LPFS ensemble.
         V_r(I) = V_r(I) * WindMultiplier
C
C     Find the velocity components.
         WVNX(I)=-V_r(I)*SIN(THETA(I))
         WVNY(I)= V_r(I)*COS(THETA(I))
C
C     jgf46.19 Convert wind velocity from top of atmospheric boundary
C     layer (which is what the Holland curve fit produces) to wind
C     velocity at 10m above the earth's surface (factor of
C     0.7).
         WVNX(I)=WVNX(I)*BLAdj
         WVNY(I)=WVNY(I)*BLAdj
C
C     jgf46.21 Also convert from 1 minute averaged winds to 10
C     minute averaged winds (0.88 factor).
         WVNX(I)=WVNX(I)*SamplingTimeAdj
         WVNY(I)=WVNY(I)*SamplingTimeAdj
C
C     jgf46.31 Add the storm translation speed.
         WVNX(I)=WVNX(I)+TransSpdX
         WVNY(I)=WVNY(I)+TransSpdY
C
C     jgf46.31 Set the wind velocities to zero outside the last closed
C     isobar.
C         IF (RAD(I).gt.rrp) THEN
C            WVNX(I)=0.0d0
C            WVNY(I)=0.0d0
C         ENDIF
C
      END DO
C
      RETURN
C     ----------------------------------------------------------------
      END SUBROUTINE HollandGet
C     ----------------------------------------------------------------

C     ----------------------------------------------------------------
C      S U B R O U T I N E   G E T  H O L L A N D  S T O R M  D A T A
C     ----------------------------------------------------------------
C
C     jgf46.08 Subroutine to support HollandGet. Gets the next line from
C     the file, skipping lines that are time repeats. Interpolates in
C     time if we are between wind data points. Does conversions to the
C     proper units. Uses old values of central pressure and RMW if the
C     line is a forecast, since forecasts do not have that data in them.
C     Assumes longitude is WEST longitude, latitude is NORTH latitude.
C
C     ----------------------------------------------------------------
C
      SUBROUTINE GetHollandStormData(LatOut,LonOut,CPressOut,SpdOut,
     &     RRPOut,RMWOut,TVXOut,TVYOut,Time,NScreen,ScreenUnit)
      USE SIZES, ONLY : MyProc
      IMPLICIT NONE
      REAL(SZ),intent(in) :: time
      INTEGER,intent(in) :: nscreen, screenunit
      REAL(SZ),intent(out) :: LatOut, LonOut, CPressOut
      REAL(SZ),intent(out) :: SpdOut, RRPOut, RMWOut, TVXOut, TVYOut
C
      INTEGER, ALLOCATABLE, SAVE :: iYear(:),iMth(:),iDay(:),iHr(:)
      INTEGER, ALLOCATABLE, SAVE :: iLat(:),iLon(:)
      INTEGER, ALLOCATABLE, SAVE :: iSpd(:),iCPress(:),iRRP(:),iRMW(:)
      REAL(SZ), ALLOCATABLE, SAVE :: Lat(:),Lon(:),Spd(:)
      REAL(SZ), ALLOCATABLE, SAVE :: CPress(:),RRP(:),RMW(:)
      REAL(SZ), ALLOCATABLE, SAVE :: TVX(:), TVY(:)      ! jgf46.32jgf9
      CHARACTER(len=4), ALLOCATABLE, SAVE :: CastType(:) !hindcast,forecast
      INTEGER, ALLOCATABLE, SAVE :: iFcstInc(:) ! hours
      REAL(SZ), ALLOCATABLE, SAVE :: FcstInc(:) ! seconds
      REAL(SZ),SAVE :: WTRATIO          !time ratio used for interpolation
      REAL(8), ALLOCATABLE, SAVE :: CastTime(:) ! seconds since start of year
      INTEGER,SAVE :: iNowcastCPress, iNowcastRRP, iNowcastRMW
      LOGICAL,SAVE :: FIRSTCALL = .True.
      INTEGER,SAVE :: i  ! Current array counter for fort.22 file
      INTEGER,SAVE :: nl ! Number of lines in the fort.22 file
      INTEGER,SAVE :: pl ! populated length of Holland Data array
      REAL(8), PARAMETER :: re=6378206.4d0
      REAL(8) :: pi
      REAL(8) :: mperdeg,DEG2RAD

        pi = 4.0 * atan(1.0)
        mperdeg  = re * pi / 180.0
        DEG2RAD = pi / 180.0

C     ------------------------------------------------------
C     BEGIN Code executed upon first call to this subroutine
C     ------------------------------------------------------
      IF (FIRSTCALL) THEN
C
C     Determine the number of lines in the file.
         nl=0
         OPEN(22,FILE=TRIM(LOCALDIR)//'/'//'fort.22')
         DO
            READ(UNIT=22,FMT='(A170)',END=8888)
            nl=nl+1
         ENDDO
 8888    CLOSE(22)
C
C     Dimension the arrays according to the number of lines in the file,
C     this will be greater than or equal to the size of the array we need
C     (probably greater because of the repeated lines that we throw away)
         ALLOCATE(iYear(nl),iMth(nl),iDay(nl),iHr(nl),iLat(nl),iLon(nl),
     &        iSpd(nl),iCpress(nl),iRRP(nl),iRMW(nl),iFcstInc(nl))
         ALLOCATE(Lat(nl),Lon(nl),Spd(nl),CPress(nl),RRP(nl),RMW(nl),
     &        FcstInc(nl),TVX(nl),TVY(nl))
         ALLOCATE(CastType(nl))
         ALLOCATE(CastTime(nl))
C
C     Now reopen the file and read the data into the arrays. The first
C     line must be a hindcast/nowcast.
         i=1
         OPEN(22,FILE=TRIM(LOCALDIR)//'/'//'fort.22')
C
         DO
C     Get another line of data from the file and check to see if the
C     line represents a new point in time, or is a repeated time
C     point. Repeated time points occur in hindcasts for the purpose of
C     describing winds in the quadrants of the storm. We don't use the
C     quadrant-by-quadrant wind data. Repeated time data occur in the
C     forecast because the time data is just the time that the forecast
C     was made. The important parameter in the forecast file is the
C     forecast increment.
            READ(UNIT=22,FMT=14,END=9999)
     &           iYear(i),iMth(i),iDay(i),iHr(i),
     &           CastType(i),iFcstInc(i),iLat(i),iLon(i),iSpd(i),
     &           iCPress(i),iRRP(i),iRMW(i)
C
C
            SELECT CASE(CastType(i))
C           ------------
            CASE("BEST")     ! nowcast/hindcast
C           ------------
C     Check to see if this is a repeated line. If so, go directly to the
C     next line without any processing.
               IF (i.gt.1.and.
     &              iYear(i).eq.iYear(i-1).and.iMth(i).eq.iMth(i-1).and.
     &              iDay(i).eq.iDay(i-1).and.iHr(i).eq.iHr(i-1)) THEN
                  CYCLE
               ENDIF
C
C     Save the central pressure, radius of last closed isobar, and
C     radius to max wind for use in forecasts
               iNowcastCPress=iCPress(i)
               iNowcastRMW=iRMW(i)
               iNowcastRRP=iRRP(i)
C
C     Determine the time of this hindcast in seconds since the beginning
C     of the year.
               CALL TimeConv(iYear(i),iMth(i),iDay(i),iHr(i),0,0.d0,
     &              CastTime(i),MyProc,NScreen,ScreenUnit)
C
C     Determine the CastTime in seconds since the beginning of the simulation.
               CastTime(i)=CastTime(i)-WindRefTime
               FcstInc(i)=iFcstInc(i)
C
C           ------------
            CASE("OFCL")        ! forecast
C           ------------
C     Check to see if this is a repeated line (i.e., a forecast that
C     coincides with the nowcast, or a repeated forecast). If so, go
C     directly to the next line without any processing.
               IF ( (iFcstInc(i).eq.0.and.
     &              (iYear(i).eq.iYear(i-1).and.iMth(i).eq.iMth(i-1)
     &              .and.iDay(i).eq.iDay(i-1).and.iHr(i).eq.iHr(i-1)))
     &              .or.
     &             (iFcstInc(i).ne.0.and.iFcstInc(i).eq.iFcstInc(i-1)))
     &              THEN
                  CYCLE
               ENDIF
               FcstInc(i) = iFcstInc(i)
C
C     Determine the time of this forecast in seconds since the beginning
C     of the year.
               IF ( iFcstInc(i).eq.0 ) THEN
                  CALL TimeConv(iYear(i),iMth(i),iDay(i),iHr(i),0,0.d0,
     &                 CastTime(i),MyProc,NScreen,ScreenUnit)
                  CastTime(i)=CastTime(i)-WindRefTime
               ELSE
                  FcstInc(i) = FcstInc(i) * 3600.d0 ! convert hours to seconds
                  CastTime(i) = CastTime(i-1) +
     &                 ( FcstInc(i) - FcstInc(i-1) )
               ENDIF
C
C     Set the central pressure and the radius to max wind to whatever
C     the nowcast values were.
               iCPress(i)=iNowcastCPress
               iRMW(i)=iNowcastRMW
               iRRP(i)=iNowcastRRP
C
            CASE DEFAULT        ! unrecognized
               WRITE(16,1000)   ! unit 22 Holland Storm File
               WRITE(16,1021) CastType(i),MyProc ! contains invalid name
               WRITE(16,1031)   ! describe valid input
               WRITE(16,1041)   ! tell which column failed
               IF (NScreen.ne.0.and.MyProc.eq.0) THEN
                  WRITE(ScreenUnit,1000) ! unit 22 Holland Storm File
                  WRITE(ScreenUnit,1025) CastType(i) ! contains invalid name
                  WRITE(ScreenUnit,1031) ! describe valid input
                  WRITE(ScreenUnit,1041) ! tell which column failed
               ENDIF
               STOP
            END SELECT
C
C     Convert integers to reals.
            Lat(i) = iLat(i)
            Lon(i) = iLon(i)
            Spd(i) = iSpd(i)
            CPress(i) = iCPress(i)
            RRP(i) = iRRP(i)
            RMW(i) = iRMW(i)
C
C     Convert units.
            Lat(i) = Lat(i) / 10.d0 ! convert 10ths of degs to degs
            Lon(i) = Lon(i) / 10.d0 ! convert 10ths of degs to degs
            Lon(i) = -1.d0 * Lon(i) ! negative b/c WEST longitude
            CPress(i) = CPress(i) * 100.d0 ! convert mbar to Pa
            RRP(i) = RRP(i) * 1.852000003180799d0 * 1000.0d0 ! convert nm to m
            RMW(i) = RMW(i) * 1.852000003180799d0 ! convert nm to km
            Spd(i) = Spd(i) * 0.51444444d0 ! convert kts to m/s
C
#ifdef DEBUG_HOLLAND
            WRITE(16,1244) CastTime(i),Lat(i),Lon(i),
     &           Spd(i),CPress(i),RRP(i),RMW(i),WindRefTime
            if ( i.gt.1 ) then
               write(16,2355) FcstInc(i-1)
            endif
            IF (NScreen.ne.0.and.MyProc.eq.0) THEN
               WRITE(ScreenUnit,1244) CastTime(i),Lat(i),Lon(i),
     &              Spd(i),CPress(i),RRP(i),RMW(i),WindRefTime
            ENDIF
 1244       FORMAT('CastTime ',e16.8,' Lat ',f6.2,' Lon ',f6.2,
     &           /,'Spd ',f8.2,' CPress ',f10.2,' RRP ',f12.2,
     &           ' RMW ',f8.2, ' WindRefTime ',e16.8)
 2355       format('FcstInc(i-1) ',e16.8)
#endif
C
C     Save the number of non-repeated lines from the fort.22 file, this
C     is the populated length of the array.
            pl=i
C
C     Increment array counter
            i=i+1
C
         ENDDO
 9999    CLOSE(22)
C
C     Calculate storm translation velocities based on change in position, then
C     convert degrees/time to m/s
         DO i=2, pl
C RJW  05.08.2009
C     use the formula haversine formula for distances along a sphere.
C         dist(i)=re*(2*ASIN(sqrt(sin(DY(I)/2)**2+
C     &  cos(lat*DEG2RAD)*cos(YCOOR(i)*DEG2RAD)*sin(DX(I)/2)**2) ) )

C     Calculate storm translation velocities based on change in position,

            TVX(i)=(re*(2.0d0*ASIN(sqrt(cos(LAT(I)*DEG2RAD)
     & * cos(LAT(I-1)*DEG2RAD)
     & * sin((LON(I)*DEG2RAD-LON(I-1)*DEG2RAD)/2.0d0)**2.0d0) ) ) )
     &              /(CastTime(i)-CastTime(i-1))

            TVY(i)=(re*(2.0d0*ASIN(sqrt(
     &  sin( (LAT(I)*DEG2RAD-LAT(I-1)*DEG2RAD)/2.0d0)**2.0d0))) )
     &              /(CastTime(i)-CastTime(i-1))


C            TVX(i)=(Lon(i)-Lon(i-1))*100.d0*1000.d0
C     &              /(CastTime(i)-CastTime(i-1))
C            TVY(i)=(Lat(i)-Lat(i-1))*100.d0*1000.d0
C     &              /(CastTime(i)-CastTime(i-1))
         ENDDO
         TVX(1)=TVX(2)
         TVY(1)=TVY(2)
C
C     Determine the correspondence between the current simulation time and
C     the fort.22 file.
         i=2
         DO
            IF (Time.ge.CastTime(i-1).and.Time.lt.CastTime(i)) THEN
               EXIT
            ELSE
               i=i+1
               IF (i.gt.pl) THEN
                  WRITE(16,1000) ! unit 22 Holland Storm File
                  WRITE(16,1051)
                  WRITE(16,1061)
 1051             FORMAT('Does not contain times/dates that correspond')
 1061             FORMAT('to the ADCIRC current model time.')
                  IF (NScreen.ne.0.and.MyProc.eq.0) THEN
                     WRITE(ScreenUnit,1000) ! unit 22 Holland Storm File
                     WRITE(ScreenUnit,1051)
                     WRITE(ScreenUnit,1061)
                  ENDIF
                  STOP
               ENDIF
            ENDIF
         ENDDO
C
         FIRSTCALL = .False.
      ENDIF
C     ----------------------------------------------------------
C     END Code executed only upon first call to this subroutine
C     ----------------------------------------------------------
C
C
C     ----------------------------------------------------------
C     BEGIN Code executed on every call to this subroutine
C     ----------------------------------------------------------
C
C     If time exceeds the next hindcast/nowcast/forecast time, increment the
C     array counter.
      IF (Time.ge.CastTime(i)) THEN
         i=i+1
      ENDIF
C
C     Interpolate w.r.t. time
      WTRATIO=(Time-CastTime(i-1))/(CastTime(i)-CastTime(i-1))
      CPressOut = CPress(i-1) + WTRATIO*(CPress(i)-CPress(i-1))
      LonOut = Lon(i-1) + WTRATIO*(Lon(i)-lon(i-1))
      LatOut = Lat(i-1) + WTRATIO*(Lat(i)-lat(i-1))
      SpdOut = Spd(i-1) + WTRATIO*(Spd(i)-Spd(i-1))
      RRPOut = RRP(i-1) + WTRATIO*(RRP(i)-RRP(i-1))
      RMWOut = RMW(i-1) + WTRATIO*(RMW(i)-RMW(i-1))
      TVXOut = TVX(i-1) + WTRATIO*(TVX(i)-TVX(i-1))
      TVYOut = TVY(i-1) + WTRATIO*(TVY(i)-TVY(i-1))
C
C     ----------------------------------------------------------
C     END Code executed on every call to this subroutine
C     ----------------------------------------------------------
 10   format(8x,i4,i2,i2,i2,6x,a4,7x,i3,4x,i3,3x,i3,2x,i4,52x,i3)
 12   format(8x,i4,i2,i2,i2,6x,a4,2x,i3,2x,i3,4x,i3,3x,i3,2x,i4,52x,i3)
 14   format(8x,i4,i2,i2,i2,6x,a4,2x,i3,1x,i4,3x,i4,3x,i3,2x,i4,47x,i3,
     &     2x,i3)
 1000 FORMAT('ERROR: The Storm Hindcast/Forecast Input File (unit 22)')
 1021 FORMAT('contains invalid TECH identifier: ',A4,' on PROC=',I4)
 1025 FORMAT('contains invalid TECH identifier: ',A4)
 1031 FORMAT('Valid TECH input is BEST (hindcast) or OFCL (forecast).')
 1041 FORMAT('Check the 5th column of the fort.22 file.')
C
      RETURN
C     ----------------------------------------------------------------
      END SUBROUTINE GetHollandStormData
C     ----------------------------------------------------------------


C     ----------------------------------------------------------------
C      S U B R O U T I N E   N W S  6 7  G E T
C     ----------------------------------------------------------------
C
C     jgf46.02 Subroutine written by Brian Blanton to calculate wind
C     velocity at nodes from the Holland Wind model. THIS SUBROUTINE IS
C     NOT CURRENTLY CALLED FROM ANYWHERE, AND IS INCLUDED HERE AS
C     REFERENCE MATERIAL. Code may be written in the future to call this
C     subroutine.
C
C     From Brian:
C
C     "This NWS67GET doesn't read a standard track file but rather it
C     reads a fort.22 file that contains along-track parameters.
C      
C     time      lon      lat      dp          du       dv     RMW    B
C                                [Pa]        [m/s]    [m/s]    km
C      0.000  -89.600   18.792 11500.000    0.000    8.200    56   1.40
C      6.000  -89.600   20.563 11500.000    0.000    8.200    56   1.40
C     12.000  -89.600   22.334 11500.000    0.000    8.200    56   1.40
C     18.000  -89.600   24.105 11500.000    0.000    8.200    56   1.40
C     24.000  -89.600   25.876 11500.000    0.000    8.200    56   1.40
C     30.000  -89.600   27.648 11500.000    0.000    8.200    56   1.40
C     36.000  -89.600   29.419 11500.000    0.000    8.200    56   1.40
C     42.000  -89.600   31.190 11500.000    0.000    8.200    56   1.40
C
C     time isn't actually used.
C     dp is press diff from ambient
C     du,dv is storm translation vel
C     not sure why RMW is in km, but it made sense at the time."
C
C     ----------------------------------------------------------------
C23456
      SUBROUTINE NWS67GET  (X,Y,SLAM,SFEA,WVNX,WVNY,PRESS,NP,ICS,
     &     RHOWAT0,G,TIME,WTIMINC,WTIME1,WTIME2)
      USE SIZES
      IMPLICIT NONE
      REAL(8), PARAMETER :: Ambient_Pressure=101300.d0      !   Pascals!!
      REAL(8), PARAMETER :: Coriolis=7.287588493541447D-05  ! 1/sec
      REAL(8), PARAMETER :: Rho_Air = 1.15d0                ! kg/ m/m/m
      INTEGER, SAVE :: FIRSTCALL = 0
      INTEGER NP,I,J,ICS
      REAL(SZ) RHOWAT0,G,WTIME1,WTIME2,TIME,WTIMINC
      REAL(SZ) X(*),Y(*),SLAM(*),SFEA(*)
      REAL(SZ) WVNX(*),WVNY(*),PRESS(*)
! LOCAL VARIABLES
      REAL(SZ),SAVE,ALLOCATABLE :: RAD(:),DX(:),DY(:),XCOOR (:),YCOOR(:)
      REAL(SZ),SAVE,ALLOCATABLE :: EARG(:)
      REAL(SZ),SAVE,ALLOCATABLE :: VG(:),THETA(:)
      REAL(SZ),SAVE :: XC1,YC1,DP1,TVX1,TVY1,RMW1,A1,B1
      REAL(SZ),SAVE :: XC2,YC2,DP2,TVX2,TVY2,RMW2,A2,B2
      REAL(SZ) :: XC,YC,DP,TVX,TVY,RMW,A,B,WTRATIO
      REAL(SZ) :: Central_Pressure,T1,T2,T3,D,N,TT,TM
      REAL(SZ),SAVE :: RHOWATG
      IF (FIRSTCALL.EQ.0) THEN
         FIRSTCALL = 1
         RHOWATG=RHOWAT0*G
         ALLOCATE (RAD(NP),DX(NP),DY(NP),XCOOR(NP),YCOOR(NP),EARG (NP))
         ALLOCATE (VG(NP),THETA(NP))
!         IF(ICS.EQ.2) THEN
            DO I=1,NP
            XCOOR(I)=SLAM(I)*RAD2DEG
            YCOOR(I)=SFEA(I)*RAD2DEG
            END DO
!         ENDIF
!         IF(ICS.EQ.1) THEN
!            DO I=1,NP
!               XCOOR(i)=X(i)
!                  YCOOR(i)=Y(i)
!            END DO
!         ENDIF
         READ(22,*) TT,XC1,YC1,DP1,TVX1,TVY1,RMW1,B1
         READ(22,*) TT,XC2,YC2,DP2,TVX2,TVY2,RMW2,B2
         OPEN(667,file='fort.22.adc')
      ENDIF
      IF(TIME.GT.WTIME2) THEN
         WTIME1=WTIME2
         WTIME2=WTIME2+WTIMINC
         XC1=XC2
         YC1=YC2
         DP1=DP2
         TVX1=TVX2
         TVY1=TVY2
         RMW1=RMW2
         B1=B2
         READ(22,*) TT,XC2,YC2,DP2,TVX2,TVY2,RMW2,B2
      END IF
      WTRATIO=(TIME-WTIME1)/WTIMINC
      DP  = DP1  + WTRATIO*(DP2-DP1)
      XC  = XC1  + WTRATIO*(XC2-XC1)
      YC  = YC1  + WTRATIO*(YC2-YC1)
      TVX = TVX1 + WTRATIO*(TVX2-TVX1)
      TVY = TVY1 + WTRATIO*(TVY2-TVY1)
      RMW = RMW1 + WTRATIO*(RMW2-RMW1)
      B   = B1   + WTRATIO*(B2-B1)

      Central_Pressure=Ambient_Pressure-DP
      A=(RMW*1000)**B
      DX=XCOOR-XC
      DY=YCOOR-YC
      THETA=ATAN2(DY,DX)
      RAD=SQRT(DX*DX+DY*DY)*100.d0*1000.d0  ! into meters!!
      EARG=EXP(-A/RAD**B)

      ! WVNX,WVNY,PRESS are returned!!
      DO I=1,NP
         PRESS(I)=(Central_Pressure + (DP)*EARG(I))/RHOWATG
         D=Rho_Air*RAD(I)**B
         N=A*B*DP*EARG(I)
         T1=N/D
         T2=((CORIOLIS**2.d0)/4.d0)*RAD(I)**2.d0
         T3=RAD(I)*CORIOLIS/2.d0
         VG(I)=sqrt(T1+T2)-T3
         TM=(RAD(I)/(RMW*1000))-5.
         IF (TM.LT.0.)THEN
            TM=0.
         ELSEIF (TM.GT.1.)THEN
            TM=1.
         END IF
         TM=-TM/2
         TM=COS(TM*PI)
         WVNX(I)=-VG(I)*SIN(THETA(I))+TVX*TM
         WVNY(I)= VG(I)*COS(THETA(I))+TVY*TM
      END DO

C      write(667,'(12(f12.4,1x))')TIME,XC,YC,TVX,TVY,DP,RMW,B
C     ----------------------------------------------------------------
      END SUBROUTINE NWS67GET
C     ----------------------------------------------------------------


      !=================================================================
      !=================================================================
      !=================================================================
      !      =====                                           =====
      !      =====            SUBROUTINE nws9get             =====
      !      =====                                           =====
      !=================================================================
      !=================================================================
      !=================================================================

      !=================================================================
      ! This subroutine takes a wind file in best track format and
      ! creates an asymmetric hurricane vortex, returning a wind field
      ! (wvnx,wvny) and atmospheric pressure-induced water surface
      ! elevation (press) on-the-fly for each model time step.
      !
      ! On input:
      !    slam     Longitude at model nodal points (radians)
      !    sfea     Latitude  at model nodal points (radians)
      !    np       Number of nodal points in model grid
      !    time     Model time (seconds)
      !    ics      Coordinate system selection parameter
      !                ics = 1: cartesian coordinates
      !                ics = 2: spherical coordinates
      !                Note: Subroutine valid only for ics = 2
      !
      ! On output:
      !    wvnx     x component of wind velocity at nodal points (m/s)
      !    wvny     y component of wind velocity at nodal points (m/s)
      !    press    Atmospheric pressure-induced water surface
      !             elevation at nodal points (m)
      !
      ! Revision history:
      !    Date         Programmer                 Description of change
      !    ----         ----------                 ---------------------
      !    06/19/06     Cristina Forbes, UNC-CEP   Wrote original code
      !    07/09/06     Craig  Mattocks, UNC-CEP   Code clean up
      !    07/20/06     Jason G. Fleming UNC-IMS   Add into v46.16
      !    07/20/06     Cristina Forbes  (wrote)   v46.16: comment out debug
      !                 Jason G. Fleming (merged)  msg, check for all 4 radii
      !    08/15/06     Cristina Forbes, UNC-CEP   if missing values or radii
      !                                            = 0 -> set winds = 0  and
      !                                            press = P0, comment out
      !                                            fort.300 output
      !    08/25/06     Cristina Forbes, UNC-CEP   Implemented hotstart
      !                                            capability
      !    08/25/06     Craig  Mattocks, UNC-CEP   Applied wind reduction:
      !                                            top of SFC layer --> SFC
      !    09/12/06     Craig  Mattocks, UNC-CEP   Subtracted translational
      !                                            wind speed from Vmax.
      !    07/05/07     Jason G. Fleming           Merged into v47.04, added
      !                                            ScreenUnit
      !    12/12/08     Cristina Forbes, UNC-IMS   Comment out debug statements
      !    04/21/08     Cristina Forbes, UNC-IMS   Added capability to compute
      !                                            winds in all hemispheres
      !    05/12/09     Robert Weaver, UNC-IMS     Modified damping of translational 
      !                                            velocity based on ratio of V/Vmax
      !    05/19/09     Cristina Forbes, UNC-IMS   Implemented old bug fix: units conversion 
      !    05/22/09     Cristina Forbes, UNC-IMS   Reverted damping of translational
      !                                            velocity back to the original formulation
      !                                            developed by Craig Mattocks
      !                                            
      !                                            
      ! Reference:
      !
      !     If you use this NWS=9 option capability please reference one of 
      !     the following publications:
      !
      !     Mattocks, C. and C. Forbes (2008): A real-time, event-triggered 
      !     storm surge forecasting system for the State of North Carolina, 
      !     Ocean Modelling, 25, pp. 95-119, 10.1016/j.ocemod.2008.06.008
      !
      !     Mattocks, C., C. Forbes and L. Ran (2006) Design and implementation 
      !     of a real-time storm surge and flood forecasting capability for the 
      !     State of North Carolina. UNC-CEP Technical Report, November 30, 2006, 103 pp.
      !
      !     http://www.unc.edu/~cmattock/ncfs/NCFSreport.pdf
      !=================================================================
      SUBROUTINE nws9get(slam,sfea, wvnx,wvny, press, np, time, ics)

         USE SIZES, ONLY : MyProc,mnproc
         !-------------------------------------------------------------
         ! Import custom precision types for cross-platform portability
         !-------------------------------------------------------------
         USE precision

         !------------------------
         ! Import global constants
         !------------------------
         USE constants

         !-----------------------------------------
         ! Import asymmetric hurricane vortex class
         !-----------------------------------------
         USE vortex

         !--------------------------------------------
         ! Force explicit declaration of all variables
         !--------------------------------------------
         IMPLICIT NONE

         INTEGER  :: np, ics
         REAL(sz) :: time, time1,time2, wtratio
         REAL( 8), SAVE, ALLOCATABLE :: xcoor(:),ycoor(:)
         REAL( 8) :: slam(np),sfea(np)
         REAL( 8) :: wvnx(np),wvny(np), press(np)

         INTEGER , SAVE :: iyear1, imth1, iday1, ihr1
         INTEGER , SAVE :: iyear2, imth2, iday2, ihr2
         INTEGER , SAVE :: ilat1,ilon1, ispd1, icpress1, irmw1
         INTEGER , SAVE :: ilat2,ilon2, ispd2, icpress2, irmw2
         REAL(sz), SAVE ::  lat1, lon1,  spd1,  cpress1,  rmw1
         REAL(sz), SAVE ::  lat2, lon2,  spd2,  cpress2,  rmw2
         REAL(sz), SAVE ::  dt, told, cLat1, cLon1
         CHARACTER(LEN=4) :: type1              ! hindcast/nowcast or forecast?
         CHARACTER(LEN=4) :: type2              ! hindcast/nowcast or forecast?
         INTEGER , SAVE :: iFcstInc1, iFcstInc2 ! hours between forecasts
         INTEGER , SAVE :: firstCall = 0
         LOGICAL , SAVE :: firstTime = .TRUE.

         INTEGER , SAVE :: adv
         REAL(sz), SAVE :: Vmax
         REAL(sz), SAVE :: Vr
         REAL(sz), SAVE, DIMENSION(4) :: r
         REAL(sz), SAVE :: Pn
         REAL(sz), SAVE :: Pc
         REAL(sz), SAVE :: cLat
         REAL(sz), SAVE :: cLon
         REAL(sz), SAVE :: offset

         INTEGER , SAVE :: ivr1, dir1, speed1
         INTEGER , SAVE :: ivr2, dir2, speed2
         INTEGER , DIMENSION(4), SAVE :: ir1,ir2
         INTEGER , SAVE :: ipn1,ipn2
         CHARACTER (LEN = 64), SAVE :: name
         REAL(sz), PARAMETER :: outFreq = thirtySixHundred
!        REAL(sz), PARAMETER :: outFreq = one
         INTEGER :: i, nwi

         REAL(sz),SAVE :: clatOld,clonOld, clatNew,clonNew,
     &                    timeOld,timeNew
         REAL(sz) :: uTrans = zero
         REAL(sz) :: vTrans = zero

         INTEGER , PARAMETER :: node2print = 6205
         CHARACTER*1 ew,E,W,ns,N,S
         E='E'
         W='W'
         N='N'
         S='S'

         timeNew = time

!         WRITE(*,*) 'SUBROUTINE nws9get: time', time

C     jgf46.16 commented out according to cf
!         IF (firstCall .EQ. 0 .OR. MOD(time,100.d0) .EQ. 0)
!     &      WRITE(*, '("SUBROUTINE nws9get: Time = ", F11.1)') time

         !----------------------------------------
         ! Read (lon,lat) at nodal points and
         ! transform them from radians --> degrees
         !----------------------------------------
         IF (firstTime) THEN
            firstTime = .FALSE.
            ALLOCATE (xcoor(np), ycoor(np))

            DO i=1,np
               xcoor(i) = slam(i) * rad2deg
               ycoor(i) = sfea(i) * rad2deg
            END DO
            offset=time
!         WRITE(*,*) 'SUBROUTINE nws9get: offset', offset
         END IF

         !---------------
         ! Initialization
         !---------------
         IF (firstCall .EQ. 0) THEN

            !----------------------------------------------
            ! Read parsed NHC advisory in best-track format
            !----------------------------------------------
            READ(22,22, END=999) adv, iyear1,imth1,iday1,ihr1, type1,
     &                           iFcstInc1, ilat1,ns,ilon1,ew, ispd1,
     &                           icpress1, ivr1,(ir1(i),i=1,4), ipn1, 
     &                           dir1,speed1, name

10          READ(22,22, END=999) adv, iyear2,imth2,iday2,ihr2, type2,
     &                           iFcstInc2, ilat2,ns,ilon2,ew, ispd2,
     &                           icpress2,ivr2,(ir2(i),i=1,4),ipn2, 
     &                           dir2, speed2, name

!
!         The logic for the code is set for W and N so if E or S
!         multiply lats and lons by -1
!
            IF(ew.EQ.E)THEN
               ilon1=(-1)*ilon1 
               ilon2=(-1)*ilon2
            ENDIF
            IF(ns.EQ.S)THEN
               ilat1=(-1)*ilat1
               ilat2=(-1)*ilat2
            ENDIF
            
            !---------------------------------
            ! Search for new NHC forecast time
            !---------------------------------
            IF (iFcstInc2 .EQ. iFcstInc1) GO TO 10

            !-----------------------------------------------
            ! Variables for calculating translation velocity
            ! at inital time = 0
            !-----------------------------------------------
!            timeOld =  iFcstInc1 * hour2sec
!            timeNew =  iFcstInc2 * hour2sec
            timeOld =  iFcstInc1 * hour2sec + offset
            timeNew =  iFcstInc2 * hour2sec + offset

!          WRITE(*,*) 'SUBROUTINE nws9get: firstCall=0, iFcstInc1,
!     &      iFcstInc2, timeOld, timeNew',
!     &      firstCall, iFcstInc1, iFcstInc2, timeOld, timeNew

            clatOld =  ilat1 * oneTenth
            clatNew =  ilat2 * oneTenth
            clonOld = -ilon1 * oneTenth
            clonNew = -ilon2 * oneTenth

         END IF   ! (firstCall .EQ. 0)

!         time1 = iFcstInc1 * hour2sec
!         time2 = iFcstInc2 * hour2sec

         time1 = iFcstInc1 * hour2sec + offset
         time2 = iFcstInc2 * hour2sec + offset
         dt=time-told

!          WRITE(*,*) 'SUBROUTINE nws9get:  time1,time2',
!     &      time1,time2

         !------------------------------------------
         ! If model reaches a new NHC forecast time,
         ! then save old values and advance forward.
         !------------------------------------------
         IF (time .EQ. time2) THEN

            iyear1    = iyear2
            imth1     = imth2
            iday1     = iday2
            ihr1      = ihr2
            type1     = type2
            iFcstInc1 = iFcstInc2
            ilat1     = ilat2
            ilon1     = ilon2
            ispd1     = ispd2
            icpress1  = icpress2
            ivr1      = ivr2
            ipn1      = ipn2
            dir1      = dir2
            speed1    = speed2
            lon1      = lon2
            lat1      = lat2
            spd1      = spd2
            cpress1   = cpress2

            DO i = 1,4
               ir1(i) = ir2(i)
            END DO

            !----------------------------------------------
            ! Read parsed NHC advisory in best-track format
            !----------------------------------------------
20          READ(22,22, END=999) adv, iyear2,imth2,iday2,ihr2, type2,
     &                           iFcstInc2,ilat2, ns,ilon2, ew,ispd2,
     &                           icpress2, ivr2,(ir2(i),i=1,4),ipn2,
     &                           dir2,speed2, name
!
!         The logic for the code is set for W and N so if E or S
!         multiply lats and lons by -1
!
            IF(ew.EQ.E)THEN
               ilon2=(-1)*ilon2
            ENDIF
            IF(ns.EQ.S)THEN
               ilat2=(-1)*ilat2
            ENDIF

            !---------------------------------
            ! Search for new NHC forecast time
            !---------------------------------
            IF (iFcstInc2 .EQ. iFcstInc1) GO TO 20

!            time1 = iFcstInc1 * hour2sec
!            time2 = iFcstInc2 * hour2sec

            time1 = iFcstInc1 * hour2sec + offset
            time2 = iFcstInc2 * hour2sec + offset

!          WRITE(*,*) 'SUBROUTINE nws9get:  iFcstInc2 .EQ. iFcstInc1
!     &      time1,time2',
!     &      time1,time2

         END IF   ! (time2 .EQ. time)

         !--------------------------------------------------------
         ! If radial velocity (and 4 wind radii) is missing in NHC
         ! advisory, then an asymmetric hurricane vortex cannot be
         ! constructed -- go to next time step.
         !--------------------------------------------------------
C     jgf46.16 added the .OR. condition according to cf
C     cf added if missing values set winds to zero and pressure to P0

         IF (ivr1 .EQ. missingInt .OR.
     &       (ir1(1)*ir1(2)*ir1(3)*ir1(4) .EQ. 0)) THEN

!            IF (firstCall .EQ. 0 .OR. MOD(time,outFreq) .EQ. 0)
!     &         WRITE(*, '("SUBROUTINE nws9get: Time = ", F8.1,
!     &                    " Missing radial velocity for asymmetric",
!     &                    " hurricane vortex calculation; continue",
!     &                    " w/next timestep ", I3)') time, firstCall

            DO i=1,np
               press(i) = Pzero
               wvnx(i) = zero
               wvny(i) = zero
            END DO

            firstCall = 2
            RETURN

         END IF   ! (ivr1 .EQ. missingInt)

         !-------------------------------------
         ! Interpolate NHC forecast interval in
         ! time to obtain values at model time.
         !-------------------------------------
         wtratio = (time-time1) / (time2-time1)

!         WRITE(*,*) 'SUBROUTINE nws9get:  wtratio
!     &      time1,time2, time',
!     &      time1,time2,time


!         IF (firstCall .EQ. 0 .OR. MOD(time,outFreq) .EQ. 0)
!     &      WRITE(ScreenUnit, 
!     &               '("SUBROUTINE nws9get: Time interpolation"   ,/,
!     &                 "   t1       t             t2"             ,/,
!     &                 "   +--------+-------------+"              ,/,
!     &                 "   time1 = ", F11.1                       ,/,
!     &                 "            time = ", F11.1               ,/,
!     &                 "                          time2 = ", F11.1,/,
!     &                 "   time fraction = ", F7.4)')
!     &           time1, time, time2, wtratio

         !--------------------------------------------------------------
         ! Perform time interpolation, transform variables from integers
         ! to real numbers for hurricane vortex calcualtions.
         !--------------------------------------------------------------
         Vmax =  one * ( ispd1    + wtratio * (ispd2    - ispd1   ) )
         Pn   =  one * ( ipn1     + wtratio * (ipn2     - ipn1    ) )
         Pc   =  one * ( icpress1 + wtratio * (icpress2 - icpress1) )
         cLat =  ilat1 * oneTenth + wtratio * (ilat2-ilat1) * oneTenth
         cLon = -ilon1 * oneTenth - wtratio * (ilon2-ilon1) * oneTenth

         Vr = one * ivr1
         IF (ivr1 .eq. ivr2) THEN
            !--------------------------------------------------------------
            ! Perform time interpolation of wind radii only when reported
            ! at the same wind intensity (34, 50, 64, 100 knots).
            !--------------------------------------------------------------
            DO i=1,4
               r(i) = one * ( ir1(i) + wtratio * (ir2(i) - ir1(i)) )
            END DO
         ELSE
            !--------------------------------------------------------------
            ! No time interpolation, just transform variables from integers
            ! to real numbers for hurricane vortex calcualtions.
            !--------------------------------------------------------------
            DO i=1,4
               r(i) = one * ir1(i)
            END DO
         END IF

!         IF (firstCall.EQ.0 .OR. MOD(time,outFreq) .EQ. 0) THEN
!            WRITE(ScreenUnit,*) 'SUBROUTINE nws9get: Time interpolation'
!            WRITE(ScreenUnit,*) 
!     &          '   ===== Vmax =====', ispd1   , Vmax, ispd2
!            WRITE(ScreenUnit,*) 
!     &          '   ===== Vr   =====', ivr1    , Vr  , ivr2
!            WRITE(ScreenUnit,*) 
!     &          '   ===== Pn   =====', ipn1    , Pn  , ipn2
!            WRITE(ScreenUnit,*) 
!     &          '   ===== Pc   =====', icpress1, Pc  , icpress2
!            WRITE(ScreenUnit,*) 
!     &          '   ===== cLat =====', ilat1   , cLat, ilat2
!            WRITE(ScreenUnit,*) 
!     &          '   ===== cLon =====', ilon1   , cLon, ilon2
!         END IF

!         IF (firstCall.EQ.0 .OR. MOD(time,outFreq) .EQ. 0)
!     &      WRITE(ScreenUnit,*) 
!     &                'SUBROUTINE nws9get: before = ', name, adv,
!     &                 iyear1,imth1,iday1,ihr1,iFcstInc1,
!     &                 Pn,Pc, cLat,cLon, Vmax, Vr,r

         !-------------------------------
         ! Calculate translation velocity
         !-------------------------------
         IF (firstCall .EQ. 1) THEN
            !-------------------------------------------------------
            ! If this is not the first time interval (not time = 0),
            ! uvTrans is calculated at smaller time intervals.
            !-------------------------------------------------------
            IF (time .NE. time1) THEN
                CALL uvtrans(clatOld, clonOld, cLat,cLon,
     &                       timeOld,timeNew, uTrans,vTrans)
                clatOld = cLat
                clonOld = cLon
                timeOld = time
            ELSE
               wtratio = (time+dt-time1) / (time2-time1)
               cLat1 =  ilat1 * oneTenth + wtratio * (ilat2-ilat1) * oneTenth
               cLon1 = -ilon1 * oneTenth - wtratio * (ilon2-ilon1) * oneTenth
               CALL uvTrans(clat, clon, cLat1,cLon1,
     &                      time1,time+dt, uTrans,vTrans)
               clatOld = cLat
               clonOld = cLon
               timeOld = time
            END IF

         ELSE IF (firstCall .EQ. 0) THEN
            !---------------------------------------------------
            ! If this is the first time interval (time = 0), use
            ! the entire time interval to calculate uvTrans.
            !---------------------------------------------------
            CALL uvtrans(clatOld,clonOld, clatNew,clonNew,
     &                   timeOld,timeNew, uTrans,vTrans)

         ELSE IF (firstCall .EQ. 2) THEN
            !---------------------------------------------------
            ! If this is the first time interval (time = 0), and
            ! data is missing at the start of the NHC advisory,
            ! use the entire time interval to calculate uvTrans.
            !---------------------------------------------------
!            timeOld =  iFcstInc1 * hour2sec
!            timeNew =  iFcstInc2 * hour2sec
            timeOld =  iFcstInc1 * hour2sec + offset
            timeNew =  iFcstInc2 * hour2sec + offset
            clatOld =  ilat1 * oneTenth
            clatNew =  ilat2 * oneTenth
            clonOld = -ilon1 * oneTenth
            clonNew = -ilon2 * oneTenth

!         WRITE(*,*) 'SUBROUTINE nws9get:  firstCall=2
!     &      timeOld,timeNew',
!     &      firstCall, timeOld,timeNew


            CALL uvtrans(clatOld,clonOld, clatNew,clonNew,
     &                   timeOld,timeNew, uTrans,vTrans)

         END IF   ! (firstCall)

         !-------------------------------------------------------
         ! Create a new asymmetric hurricane vortex.
         !
         ! Note: Subtract translational speed from Vmax, then
         ! scale (Vmax - Vt) and Vr up to the top of the surface,
         ! where the cylcostrophic wind balance is valid.
         !-------------------------------------------------------
! MCF implemented old bug fix: units conversion 5/19/09
         CALL newVortex(name, adv, iyear1,imth1,iday1,ihr1, iFcstInc1,
     &                  Pn,Pc, cLat,cLon,
     &                  ( (Vmax -  ms2kt*SQRT(uTrans*uTrans + vTrans*vTrans))
     &                    /windReduction), (Vr/windReduction), r)

         !---------------------------
         ! Print the hurricane vortex
         !---------------------------
!        IF (firstCall .EQ. 0 .OR. MOD(time,outFreq) .EQ. 0)
!     &     CALL printVortex()

         !-------------------------------------------------------------
         ! Calculate wind and pressure fields at model nodal points.
         !
         ! Note: the asymmetric vortex wind speed is reduced from the
         ! top of the surface layer to the surface, then converted from
         ! a 1-minute (max sustained) to a 10-minute average prior to
         ! adding the translational velocity in subroutine uvp.
         !-------------------------------------------------------------
         DO i=1,np
! RJW added to pass Vmax into uvp subroutine
! MCF implemented old bug fix: units conversion 5/19/09
!            CALL uvp(ycoor(i),xcoor(i), uTrans,vTrans,
!     &               wvnx(i),wvny(i), press(i),
!     &               ((Vmax -  ms2kt*SQRT(uTrans*uTrans + vTrans*vTrans))
!     &                 /windReduction ) )
            CALL uvp(ycoor(i),xcoor(i), uTrans,vTrans,
     &               wvnx(i),wvny(i), press(i))
!            WRITE (300,301) i, xcoor(i),ycoor(i),
!     &                      wvnx(i),wvny(i),press(i), uTrans,vTrans
301         FORMAT(i7, 7(1x,f12.4))

            !-------------------------------------------
            ! Convert atmospheric pressure (Pascals) to
            ! atmospheric pressure-induced water surface
            ! elevation (meters).
            !-------------------------------------------
            press(i) = press(i) / RhoWatG
         END DO

         !------------------------------------
         ! Write wind and pressure fields to a
         ! diagnostic output file (fort.300)
         !------------------------------------
!         IF (firstCall .EQ. 0 .OR. MOD(time,outFreq) .EQ. 0) THEN
!            nwi = time * sec2hour
!           nwi = time
!
!            DO i=1,np
!               WRITE (300,301) nwi,i, xcoor(i),ycoor(i),
!     &                         wvnx(i),wvny(i),press(i), uTrans,vTrans
!301            FORMAT(2(i7,1x), 7(1x,f9.4))
!            END DO
!         END IF

         !-----------------------------------
         ! NHC advisory best-track i/o format
         !-----------------------------------
22      FORMAT(3x, i3, 2x, i4, 3i2, 6x, a4, 2(2x,i3),a1, 3x, i3,a1, 2x, 
     &         i3,2x, i4, 6x, i3, 7x, 4(i4,2x), i4, 38x, 2(i3,2x), a10)

         firstCall = 1
         told=time

!        WRITE(*,*) 'SUBROUTINE nws9get: node #',node2print,' u,v,p = ',
!    &              wvnx(node2print),wvny(node2print), press(node2print)

999      RETURN

      END SUBROUTINE nws9get

C***********************************************************************
C                                                                      *
C   End of subroutines to read wind and pressure fields                * 
C                                                                      *
C***********************************************************************


C***********************************************************************
C                                                                      *
C   Read onto the ADCIRC grid radiation stress fields in the PBL-JAG   *
C   (hurricane) model format.                                          *
C                                                                      *
C                                                                      *
C                           R.L.05/12/99                               *
C***********************************************************************

      SUBROUTINE RSGET(RSNX,RSNY,NP)
      USE SIZES
      IMPLICIT NONE
      INTEGER NP,I,NHG
      REAL(SZ) RSNX(*),RSNY(*)
      CHARACTER*80 PBLJAGF
C     
      DO I=1,NP
         RSNX(I)=0.d0
         RSNY(I)=0.d0
      END DO
 170  READ(23,'(A80)') PBLJAGF
      IF(PBLJAGF(2:2).EQ.'#') GOTO 170
 171  READ(PBLJAGF,'(I8,5E13.5)') NHG,RSNX(NHG),RSNY(NHG)
      READ(23,'(A80)') PBLJAGF
      IF(PBLJAGF(2:2).NE.'#') GOTO 171
C     
      RETURN
      END SUBROUTINE


      END MODULE

