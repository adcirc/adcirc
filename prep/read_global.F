C----------------------------------------------------------------------------
C
C                           MODULE READ_GLOBAL
C
C----------------------------------------------------------------------------
C
C                  For use with ADCPREP Version 2.1 ( 03/17/2006 )
C
C                     current for ADCIRC v43.03   5/20/2003
C                 jgf Updated for ADCIRC v45.07  11/08/2005
C                 jgf Updated for ADCIRC v45.09  01/09/2006
C                 jgf Updated for ADCIRC v45.11  01/09/2006
C                 jgf Updated for ADCIRC v45.12  03/17/2006
C----------------------------------------------------------------------------

C
C  version 1.2 vjp 12/7/99
c  version 1.3 ral 10/10/01 revisions
C  version 1.6 meb 3/03 & ral 5/21/03
C  version 1.7 meb 4/04 added by jgf
C
      SUBROUTINE READ14()
      USE PRE_GLOBAL
      use memory_usage
      IMPLICIT NONE
C
C----------------------------------------------------------------------
C     Reads the global ADCIRC grid data file for data decomposition 
C     program ADCPP.  This version compatible with ADCIRC_v34.03
C     vjp  2/28/98
C---------------------------------------------------------------------
C
      integer :: nbytes = 0
      INTEGER I,J,JW,K,ITEMP,ITYPE
      INTEGER DISC,BBN,IBP,I1
      INTEGER MAXNEIGH
      INTEGER NUMLINES
c
c Dummy variables
c
      INTEGER DNBOU,DNVEL,DITYPE
      INTEGER,ALLOCATABLE :: DNVELL(:),DNBVV(:,:)
      INTEGER,ALLOCATABLE :: DIBTYPE(:)
      CHARACTER(80) DNBOUMSG,DNVELMSG
      CHARACTER(80),ALLOCATABLE ::  DNVELLMSG(:)
      CHARACTER(80) DUM

      ALLOCATE ( DNVELL(MNBOU),DNBVV(MNBOU,0:MNVEL),DIBTYPE(MNBOU) )
      ALLOCATE ( DNVELLMSG(MNBOU+1) )
      nbytes = 8*mnbou + 16*(mnbou+1) + 8*mnbou*(mnvel+1)
      call memory_alloc(nbytes)
C
C--The Grid file was opened in SIZEUP performed a rewind, and is ready here
C
C--Read Grid Title
C
      EXIST_FLUX = 0
      READ(14,80) AGRID
C
C--Read Total Number of Elements and Nodes
C
      READ(14,*) NELG,NNODG
c     READ(14,80) NSIZES
c     READ(NSIZES,*) NELG,NNODG
c     CALL GETMSG(NSIZES,SIZEMSG)
C
C--Read Nodal Coordinates and Bathymetry 
C  If ICS=2 Will Convert later in read15 
C
      NUMLINES=2
      DO I = 1,NNODG
         READ(14,*) J,X(J),Y(J),DP(J)
         IF (J.NE.I) THEN
            print *, I,J
            STOP 'Node Numbering not in Sequential Order'
         ENDIF
      ENDDO
C      
C--Read Element Connectivity Table
C
       DO I = 1,NELG
          READ(14,*) J,ITEMP,NNEG(1,J),NNEG(2,J),NNEG(3,J)
          IF (J.NE.I) THEN
            print *, I,J
            STOP 'Element Numbering not Sequential'
          ENDIF
       ENDDO
      NUMLINES=NUMLINES+NNODG+NELG
C...
C...Read Total Number of Open Boundary Segments
C...
      READ(14,80) NOPEMSG
      READ(NOPEMSG,*) NOPE
      NUMLINES=NUMLINES+1
C...
C...Read Total Number of Open Boundary Forcing Nodes
C...
      READ(14,80) NETAMSG
      READ(NETAMSG,*) NETA
      NUMLINES=NUMLINES+1
C...
C...Read Number of Nodes on Open Boundary Segment 
C...and Segment Nodes Numbers
C...
      J=0
      IBTYPEE(:) = 0 !jgf49.32 initialize all elements of the array to zero
      DO K=1,NOPE
         READ(14,80) NVDLLMSG(K) !jgf49.32: TODO: IBTYPEE should be read here
         NUMLINES=NUMLINES+1
         READ(NVDLLMSG(K),*) NVDLL(K)
         DO I=1,NVDLL(K)
            READ(14,*) NBDV(K,I)
            NUMLINES=NUMLINES+1
         ENDDO
         J=J+NVDLL(K)
      ENDDO
C     
      IF (NETA.NE.J) THEN
         print *, "Total Number of Boundary Nodes = ",J
         print *, "This exceeds NETA = ",NETA
         IF (NFOVER.EQ.1) THEN
            NETA = J
            print *, "ADCPP corrected this error"
         ELSE
            stop
         ENDIF
      ENDIF
C
C--Read Total Number of Land Boundary Segments
C     
      READ(14,80) NBOUMSG
      READ(NBOUMSG,*) NBOU
C
C--Read Total of Land Boundary Nodes
C
      READ(14,80) NVELMSG
      READ(NVELMSG,*) NVEL
C
C--Read Number of Nodes in the Land Boundary Segment and Boundary Type
C  and construct LBCODE array for read15 routine
C
      J=0
      NWEIR = 0
      DO K = 1,NBOU
C
         READ(14,80) NVELLMSG(K)
         READ(NVELLMSG(K),*) NVELL(K),IBTYPE(K)
         ITYPE = IBTYPE(K)
C
         DO I=1, NVELL(K)
            J = J+1
            LBCODE(J) = ITYPE
         ENDDO
C
C     jgf46.21 Added support for IBTYPE=52.
         IF   ((ITYPE.NE.0).AND.(ITYPE.NE.10).AND.(ITYPE.NE.20)
     &    .AND.(ITYPE.NE.1).AND.(ITYPE.NE.11).AND.(ITYPE.NE.21)      
     &    .AND.(ITYPE.NE.2).AND.(ITYPE.NE.12).AND.(ITYPE.NE.22)      
     &    .AND.(ITYPE.NE.52)
     &    .AND.(ITYPE.NE.3).AND.(ITYPE.NE.13).AND.(ITYPE.NE.23)      
     &    .AND.(ITYPE.NE.4).AND.                  (ITYPE.NE.24)
     &    .AND.(ITYPE.NE.5).AND.                  (ITYPE.NE.25)
     &    .AND.                 (ITYPE.NE.30).AND.(ITYPE.NE.32)
     &    .AND.                 (ITYPE.NE.40).AND.(ITYPE.NE.41)) THEN
            print *, "IBTYPE not set correctly for segment ",K
            stop
         ENDIF
C
         IF ((ITYPE.NE.3).AND.(ITYPE.NE.13).AND.(ITYPE.NE.23).AND.
     &     (ITYPE.NE.4).AND.(ITYPE.NE.24) .AND.
     &     (ITYPE.NE.52).AND.
     &     (ITYPE.NE.2).AND.(ITYPE.NE.12).AND.(ITYPE.NE.22)) THEN!jgf45.06
           DO I=1,NVELL(K)
              READ(14,*) NBVV(K,I)
              IBCONNR(K,I) = 0
           ENDDO
         ENDIF
C
         IF ((ITYPE.EQ.3).OR.(ITYPE.EQ.13).OR.(ITYPE.EQ.23)) THEN
           DO I=1,NVELL(K)
              READ(14,*) NBVV(K,I),BAR1(K,I),BAR2(K,I)
              IBCONNR(K,I) = 0
           ENDDO
         ENDIF
C
         IF ((ITYPE.EQ.4).OR.(ITYPE.EQ.24)) THEN
           DO I=1,NVELL(K)
              READ(14,*) NBVV(K,I),IBCONNR(K,I),
     &                   BAR1(K,I),BAR2(K,I),BAR3(K,I)
C
C--Construct List of WEIR nodes and their duals
C
              NWEIR = NWEIR+1
              WEIR(NWEIR) = NBVV(K,I)
              WEIRD(NWEIR) = IBCONNR(K,I)
           ENDDO
         ENDIF
         IF ((ITYPE.EQ.2).OR.(ITYPE.EQ.12).OR.(ITYPE.EQ.22).OR.
     &        (ITYPE.EQ.52)) THEN !jgf45.06
           DO I=1,NVELL(K)
              EXIST_FLUX = EXIST_FLUX + 1
              READ(14,*) NBVV(K,I)
              IBCONNR(K,I) = 0
           ENDDO
         ENDIF
C
      ENDDO

C     jgf45.06 added the following section
C MEB --------------------------------------------------
C MEB This entire section written 04/01/04 to accomodate
C MEB division of fort.20 files between PE directories.
C MEB This must be modified if fort.14 format changes.
C MEB --------------------------------------------------
      IF (EXIST_FLUX.GT.0) THEN 
         REWIND(14)
         ALLOCATE ( FLUX14_ARY(EXIST_FLUX) )
         EXIST_FLUX=0
c
C  SKIP OVER UNNEEDED LINES
c
         DO I=1,NUMLINES
            READ(14,*) DUM
         ENDDO
C
         READ(14,80) DNBOUMSG
         READ(NBOUMSG,*) DNBOU
C
C--Read Total of Land Boundary Nodes
C
         READ(14,80) DNVELMSG
         READ(NVELMSG,*) DNVEL
c
         DO K=1,DNBOU
            READ(14,80) DNVELLMSG(K)
            READ(DNVELLMSG(K),*) DNVELL(K), DIBTYPE(K)
            DITYPE=DIBTYPE(K)
C     jgf46.21 Added support for IBTYPE=52.
            IF ((DITYPE.EQ.2).OR.(DITYPE.EQ.12).OR.(DITYPE.EQ.22).OR.
     &           (DITYPE.EQ.52)) THEN
               DO I=1, DNVELL(K)
                  EXIST_FLUX = EXIST_FLUX + 1
                  READ(14,*) DNBVV(K,I) 
                  FLUX14_ARY(EXIST_FLUX) = DNBVV(K,I)
               ENDDO
            ELSE
               DO I=1, DNVELL(K)
                  READ(14,*) DNBVV(K,I)
               ENDDO
            ENDIF
         ENDDO
         DEALLOCATE ( DNVELL,DNBVV,DIBTYPE )
         DEALLOCATE ( DNVELLMSG )
         nbytes = 8*mnbou + 16*(mnbou+1) + 8*mnbou*(mnvel+1)
         call memory_dealloc(nbytes)
      ENDIF
c     jgf45.06 end of section written by MEB04/01/04
C
C--Close Global Grid file
C
      CLOSE(14)
      call memory_status()
C
 80   FORMAT(A80)
C
      RETURN
      END


      SUBROUTINE READ15() 
      USE PRE_GLOBAL 
      USE GLOBAL, ONLY : parse, a2f
      USE PRESIZES, ONLY : SZ
      use memory_usage
      IMPLICIT NONE

C    ===================================================================
C      MCF ADDITION FOR READING AND WRITING STATION INFO IN NETCDF
C
C      REAL, EXTERNAL              :: a2f
      CHARACTER(132) STATLINE
      CHARACTER(50) LVAR(3)
C
C    ===================================================================


C
C----------------------------------------------------------------------
C     Reads the global ADCIRC Input Data File for data decomposition 
C     program ADCPP.  This version compatible with ADCIRC_v34.03
C     vjp  3/28/98
C---------------------------------------------------------------------
C
      integer :: nbytes = 0
      REAL(SZ) RSTIMINC
      INTEGER N1, N2, N3, KMIN, JG, INDX
      INTEGER I,J,K,M
      INTEGER IG1,IG2,IG3,IL1,IL2,IL3
      REAL(8) X1, X2, X3, X4, Y1, Y2, Y3, Y4, A1, A2, A3
      REAL(8) SL1Sta,SF1Sta
      REAL(8) AE, AEMIN, AREASK, AA
      INTEGER NBV(MNVEL)
Casey 090825: Added dummy variables.
      INTEGER IDUM
      REAL(8) RDUM
C
C--The Run Info file was opened in SIZEUP, performed a rewind, and is ready here
C
C--Run Description and Run Identification
C

      READ(15,80) RUNDES
C
      READ(15,80) RUNID
C
      READ(15,80) OVERMSG
      READ(OVERMSG,*) NFOVER
      IF (NFOVER.EQ.1) THEN
c       print *, "Non-fatal errors will be corrected"
      ELSE
c      print *, "Non-fatal errors will stop execution"
      ENDIF
C
      READ(15,80) ABOUTMSG
      READ(ABOUTMSG,*) NABOUT
C
      READ(15,80) SCREENMSG
      READ(SCREENMSG,*) NSCREEN
C
      READ(15,80) HOTMSG
      READ(HOTMSG,*) IHOT
C
      READ(15,80) ICSMSG
      READ(ICSMSG,*) ICS
      IF ((ICS.NE.1).AND.(ICS.NE.2)) THEN
        print *, "ICS set incorrectly"
        STOP
        ENDIF
C
      READ(15,80) IMMSG
      READ(IMMSG,*) IM
      IF(IM.EQ.2) THEN
        PRINT *, "DSS Model type not presently supported"
        STOP
      ENDIF
C     jgf46.28 Read IDEN if necessary
      IF (CBaroclinic) READ(15,*) IDEN
C
      READ(15,80) IBFMSG
      READ(IBFMSG,*) NOLIBF
      IF((NOLIBF.LT.0).OR.(NOLIBF.GT.2)) THEN
         print *, "Value for NOLIBF not allowed"
         stop
      ENDIF
C
      READ(15,80) IFAMSG
      READ(IFAMSG,*) NOLIFA
      IF ((NOLIFA.LT.0).OR.(NOLIFA.GT.3)) THEN
         print *, "Value for NOLIFA not allowed"
         stop
      ENDIF
C
      READ(15,80) ICAMSG
      READ(ICAMSG,*) NOLICA
      IF ((NOLICA.LT.0).OR.(NOLICA.GT.1)) THEN
         print *, "Value for NOLICA not allowed"
         stop
      ENDIF
C
      READ(15,80) ICATMSG
      READ(ICATMSG,*) NOLICAT
      IF ((NOLICAT.LT.0).OR.(NOLICAT.GT.1)) THEN
         print *, "Value for NOLICAT not allowed"
         stop
      ENDIF
      IF ((NOLIFA.GE.1).AND.(NOLICAT.EQ.0)) THEN
         print *, "NOLIFA and NOLICAT are inconsistent"
         print *, "May lead to mass balance problems"
         IF(NFOVER.EQ.1) THEN
            print *, "Since NFOVER=1, Program will continue"
         ELSE
            stop
         ENDIF
      ENDIF
C
      READ(15,80) NWPMSG
      READ(NWPMSG,*) NWP
      IF (NWP.gt.0) THEN !jgf46.00 read nodal attribute labels
         ALLOCATE(NodalAttributes(NWP))
         nbytes = 4*nwp
         call memory_alloc(nbytes)
         DO I=1, NWP
            READ(15,80) NodalAttributes(I)
         ENDDO
      ENDIF
C
      READ(15,80) NCORMSG
      READ(NCORMSG,*) NCOR
      IF ((NCOR.NE.0).AND.(NCOR.NE.1)) THEN
         print *, "Value for NCOR not allowed"
         IF (NFOVER.EQ.1) THEN
           NCOR = 0
           print *, "NCOR has been reset to 0"
         ELSE
           stop
         ENDIF
      ENDIF
C
      IF ((ICS.EQ.1).AND.(NCOR.EQ.1)) THEN
        print *, "ICS=1 and NCOR=1 may lead to geometric distortions"
         IF(NFOVER.EQ.1) THEN
           print *, "Program will continue with these input values"
           print *, "for large domains it is recommended to use ICS=2"
         ELSE
           stop
         ENDIF
      ENDIF
C
      READ(15,80) NTIPMSG
      READ(NTIPMSG,*) NTIP
      IF ((NTIP.LT.0).OR.(NTIP.GT.2)) THEN
         print *, "Value for NTIP not allowed"
         IF(NFOVER.EQ.1) THEN
           NTIP = 0
           print *, "NTIP has been reset to 0"
         ELSE
           stop
         ENDIF
      ENDIF
C
      IF ((ICS.EQ.1).AND.(NTIP.GE.1)) THEN
         print *, "ICS=1 & NTIP >= 1 may lead to geometric distortions"
         print *, "for large domains it is recommended to use ICS=2"
         IF (NFOVER.EQ.1) THEN
            print *, "Program will continue with these input values"
         ELSE
            stop
         ENDIF
      ENDIF
C
      NRS=0
      READ(15,80) NWSMSG
      READ(NWSMSG,*) NWS
      IF(ABS(NWS/100).EQ.1) THEN ! sb46.28sb03
         NRS=1
         NWS=(ABS(NWS)-100)*(NWS/ABS(NWS))
      ENDIF
C     sb46.28sb03 Added NWS=2xx for STWAVE output direct read 09/xx/2006
      IF(ABS(NWS/100).EQ.2) THEN
         NRS=2
         NWS=(ABS(NWS)-200)*(NWS/ABS(NWS))
      ENDIF
Casey 090302: Added NWS=3xx for coupling to SWAN.
      IF(ABS(NWS/100).EQ.3) THEN
        NRS=3
        NWS=(ABS(NWS)-300)*(NWS/ABS(NWS))
      ENDIF
C     jgfdebug46.02 added NWS=45
C     jgf46.02 Added NWS=8.
C     jgf46.16 merged:
C     cf & cm added NWS=9: asymmetric hurricane wind model
C     sb46.28sb01 added NWS=12: OWI format
        IF((NWS.NE.0).AND.    (NWS.NE.1 ) .AND.(ABS(NWS).NE.2).AND.
     &       (NWS.NE.3).AND.(ABS(NWS).NE.4) .AND.(ABS(NWS).NE.5).AND.
     &       (ABS(NWS).NE.45).AND.(ABS(NWS).NE.6).AND.
     &       (ABS(NWS).NE.8).AND.(ABS(NWS).NE.9).AND.
     &       (ABS(NWS).NE.12) .AND. (ABS(NWS).NE.19))  THEN
           print *, "Value for NWS not supported by parallel code"
           stop
        ENDIF
C
C     jgf46.08 Modified to accomodate fine grained ramp functions.
      READ(15,80) RAMPMSG
      READ(RAMPMSG,*) NRAMP
      IF ((NRAMP.NE.0).AND.(NRAMP.GT.6)) THEN
         print *, "Value for NRAMP not allowed"
         IF (NFOVER.EQ.1) THEN
            print *, "Program will override and use NRAMP = 0"
            NRAMP = 0
         ELSE
            stop
         ENDIF
      ENDIF
C
      READ(15,80) GMSG
      READ(GMSG,*) G
      IF ((ICS.EQ.1).AND.(G.NE.9.81d0)) THEN
         IF ((NCOR.EQ.1).OR.(NTIP.EQ.1)) THEN
            print *, "G not consistent with ICS=1"
            print *, "in conjunction with NTIP=1 and/or NCOR=1"
            IF(NFOVER.EQ.1) THEN
               print *, "Program will override and set G=9.81"
               print *, "check to see that all input has SI units"
               G = 9.81d0
            ELSE
               stop
            ENDIF
         ENDIF
      ENDIF
C
      IF ((ICS.EQ.2).AND.(G.NE.9.81d0)) THEN
         print *, "G not consistent with ICS=2"
         IF(NFOVER.EQ.1) THEN
            print *, "Program will override and set G = 9.81 m/sec*sec"
            print *, "check to see that all input has SI units"
            print *, "execution will continue"
            G = 9.81d0
         ELSE
            stop
         ENDIF
      ENDIF
C     
C     jgf47.11 Allow the limits of time varying tau0 to be read in from
C     the fort.15 file
      READ(15,80) TAU0MSG
      READ(TAU0MSG,*) TAU0
      IF ( (TAU0.LE.-5.d0).AND.(TAU0.GT.-6.d0) ) THEN
         READ(15,80) TAU0LIMMSG
         READ(TAU0LIMMSG,*) Tau0FullDomainMin, Tau0FullDomainMax
      ENDIF
       
C     
      READ(15,80) DTMSG
      READ(DTMSG,*) DT
C
      READ(15,80) STATMSG
      READ(STATMSG,*) STATIM
C
      READ(15,80) REFTMSG
      READ(REFTMSG,*) REFTIM
C
C--If wind stress and surface pressures are applied process this.
C
      IF((NWS.EQ.0).AND.(NRS.GE.1)) READ(15,*) RSTIMMSG ! sb46.28sb03
      IF((NWS.EQ.1).AND.(NRS.GE.1)) READ(15,*) RSTIMMSG ! sb46.28sb03
C
      IF(NWS.EQ.3) THEN
         READ(15,80) WSMSG1
         READ(15,80) WSMSG2
         IF(NRS.EQ.0) THEN
            READ(WSMSG2,*) NWLAT,NWLON,WLATMAX,WLONMIN,WLATINC,
     &           WLONINC,WTIMINC
         ELSE
            READ(WSMSG2,*) NWLAT,NWLON,WLATMAX,WLONMIN,WLATINC,
     &           WLONINC,WTIMINC,RSTIMINC
         ENDIF
      ENDIF
C
C     jgfdebug46.02 Added NWS=45
C     jgf46.02 Added NWS=8
C     jgf46.16 Merged:
C     cf & cm Added NWS=9: asymmetric hurricane wind model
C     rjw Added NWS=19: asymmetric hurricane wind model v2.0
C     sb46.28sb01 added NWS=12: OWI format
      IF((ABS(NWS).EQ.2).OR.(ABS(NWS).EQ.4).OR.(ABS(NWS).EQ.5)
     &     .OR.(ABS(NWS).EQ.45).OR.(ABS(NWS).EQ.8).OR.
     &         (ABS(NWS).EQ.9).OR.(ABS(NWS).EQ.12)
     &     .OR.(ABS(NWS).EQ.19)) THEN
         READ(15,80) WSMSG1
         IF(NRS.EQ.0) READ(WSMSG1,*) WTIMINC
Casey 090825: Fix for NRS = 3.
C        IF(NRS.GE.1) READ(WSMSG1,*) WTIMINC,RSTIMINC ! sb46.28sb03
         IF((NRS.EQ.1).OR.(NRS.EQ.2)) READ(WSMSG1,*) WTIMINC,RSTIMINC ! sb46.28sb03
         IF(NRS.EQ.3)THEN
            IF(ABS(NWS).EQ.8)THEN
               READ(WSMSG1,*) IDUM,IDUM,IDUM,IDUM,IDUM,RDUM,RSTIMINC
            ELSE
               READ(WSMSG1,*) WTIMINC,RSTIMINC
            ENDIF
         ENDIF
      ENDIF
C
      IF(ABS(NWS).EQ.6) THEN
         READ(15,80) WSMSG1
         IF(NRS.EQ.0) THEN
            READ(WSMSG1,*) NWLAT,NWLON,WLATMAX,WLONMIN,WLATINC,
     &           WLONINC,WTIMINC
         ELSE
            READ(WSMSG1,*) NWLAT,NWLON,WLATMAX,WLONMIN,WLATINC,
     &           WLONINC,WTIMINC,RSTIMINC
         ENDIF
      ENDIF
C
      READ(15,80) RNDAYMSG
      READ(RNDAYMSG,*) RNDAY
C
      READ(15,80) DRAMPMSG
      READ(DRAMPMSG,*) DRAMP
C
      READ(15,80) COEFMSG
      READ(COEFMSG,*) A00,B00,C00
C
      READ(15,80) H0MSG
      IF (NOLIFA.NE.2) THEN
         READ(H0MSG,*) H0
      ELSE
         READ(H0MSG,*) H0,NODEDRYMIN,NODEWETRMP,VELMIN
      ENDIF
C
      READ(15,80) SLMSG
      READ(SLMSG,*) SLAM0,SFEA0
C
      SL0=DEG2RAD*SLAM0
      SF0=DEG2RAD*SFEA0

      DO I = 1,NNODG
         SLAM(I) = X(I)
         SFEA(I) = Y(I)
      ENDDO
      
!      WRITE(6,*)"SLAM0,SFEA0", SLAM0,SFEA0
C
C--If ICS = 2 then apply CPP projection
C
      IF (ICS.EQ.2) THEN
         DO I = 1,NNODG
            SL1(I) = DEG2RAD*SLAM(I)
            SF1(I) = DEG2RAD*SFEA(I)
            SLAM(I) = R*(SL1(I)-SL0)*COS(SF0)
            SFEA(I) = R*SF1(I)
         ENDDO
      ENDIF
      

C
      READ(15,80) TAUMSG
      IF (NOLIBF.EQ.0) THEN
         READ(TAUMSG,*) TAU
      ELSEIF (NOLIBF.EQ.1) THEN
         READ(TAUMSG,*) CF
      ELSEIF (NOLIBF.EQ.2) THEN
         READ(TAUMSG,*) CF,HBREAK,FTHETA,FGAMMA
      ENDIF
C
      READ(15,80) ESLMSG
      IF (IM.EQ.10) THEN
         READ(ESLMSG,*) ESLM,ESLC
      ELSE 
         READ(ESLMSG,*) ESLM
      ENDIF
C
      READ(15,80) CORIMSG
      READ(CORIMSG,*) CORI
C


      READ(15,80) NTIFMSG
      READ(NTIFMSG,*) NTIF
      IF (NTIF.GT.MNTIF) THEN
         print *, "NTIF = ",NTIF, " exceeds parameter MNTIF = ",MNTIF
         stop 
      ENDIF
C
      DO I=1,NTIF
        READ(15,80)  TIPOTAG(I)
        READ(15,80)  TPKMSG(I)
        READ(TPKMSG(I),*)  TPK(I),AMIGT(I),ETRF(I),FFT(I),FACET(I)
      ENDDO
C
      IF (((NTIP.EQ.0).AND.(NTIF.NE.0)).OR.((NTIP.NE.0).AND.
     &    (NTIF.EQ.0))) THEN
        print *, "NTIF and NTIP are not consistent"
        IF(NFOVER.EQ.1) THEN
          print *, "Program will reset NTIP = 0 and continue"
          NTIP = 0
        ELSE
          stop
        ENDIF
      ENDIF
C
      READ(15,80) NBFRMSG
      READ(NBFRMSG,*) NBFR
      IF (NBFR.GT.MNBFR) THEN
        print *, "NBFR = ",NBFR, " exceeds parameter MNBFR = ",MNBFR
        stop 
      ENDIF
C
      DO I=1,NBFR
        READ(15,80) BOUNTAG(I)
        READ(15,80) AMIGMSG(I)
        READ(AMIGMSG(I),*) AMIG(I),FF(I),FACE(I)
      ENDDO
C
      DO I=1,NBFR
        READ(15,80) ALPHA1(I)
        DO J=1,NETA
          READ(15,80) EMOMSG(I,J)
          READ(EMOMSG(I,J),*) EMO(I,J),EFA(I,J)
        ENDDO
      ENDDO
C
      READ(15,80) ANGMSG
      READ(ANGMSG,*) ANGINN
C
C--Determine whether there are any normal flow boundaries.
C
C     jgf46.21 Added support for IBTYPE=52.
      NFLUXF = 0
      DO K=1, NBOU
         IF ((IBTYPE(K).EQ.2).OR.(IBTYPE(K).EQ.12)
     &        .OR.(IBTYPE(K).EQ.22).OR.(IBTYPE(K).EQ.32)
     &        .OR.(IBTYPE(K).EQ.52)) THEN
            NFLUXF = 1
         ENDIF
      ENDDO
C     
C--If so, read the number of frequencies present and other info.
C  Use LBCODE array constructed in read14.
C
      IF (NFLUXF.EQ.1) THEN
         READ(15,80) NFFRMSG
         READ(NFFRMSG,*) NFFR
         IF (NFFR.GT.MNFFR) THEN
            print *, "NFFR = ",NFFR, " exceeds parameter MNFFR = ",MNFFR
            stop
         ENDIF         
C     
         IF (NFFR.NE.0) THEN
            DO I=1,NFFR
               READ(15,80) FBOUNTAG(I)
               READ(15,80) FREQMSG(I)
               READ(FREQMSG(I),*) FAMIG(I),FFF(I),FFACE(I)
            ENDDO
            DO I=1,NFFR
               READ(15,80) ALPHA2(I)
               DO J=1,NVEL
                  IF ((LBCODE(J).EQ.2).OR.(LBCODE(J).EQ.12)
     &                 .OR.(LBCODE(J).EQ.22).OR.(LBCODE(J).EQ.52)) THEN
                     READ(15,80) QNMSG(I,J)         
                     READ(QNMSG(I,J),*) QNAM(I,J),QNPH(I,J)
c     dbug             print *, "disc lbnode index = ",J
                  ENDIF
               ENDDO
            ENDDO        
            
Cvjp 5/1/99  added to help localize the flow boundary nodes
C  NBV    =  global node number of all boundary nodes
C  NFLBN  =  number of flow boundary nodes
C  FLBN   =  global node number of flow boundary nodes
C  FLBNX  =  index of flow boundary nodes as per NBV
C
            JG = 1
            DO K = 1,NBOU
               DO I=1, NVELL(K)
                  INDX = NBVV(K,I)
                  NBV(JG) = INDX
                  JG = JG + 1
               ENDDO
            ENDDO
C     
            NFLBN = 0
            DO J=1,NVEL
               IF ((LBCODE(J).EQ.2).OR.(LBCODE(J).EQ.12)
     &              .OR.(LBCODE(J).EQ.22).OR.(LBCODE(J).EQ.52)) THEN
                  NFLBN = NFLBN+1
                  FLBN(NFLBN) = NBV(J)
                  FLBNX(NFLBN) = J
               ENDIF
            ENDDO
         ELSE
            APERIODIC_FLOW_BC = .TRUE. !jgf45.09 need to break up a fort.20
         ENDIF
      ENDIF

C
C--Read Elevation Recording Stations
C...
      READ(15,80) STAEMSG
      READ(STAEMSG,*) NOUTE,TOUTSE,TOUTFE,NSPOOLE
      IF (ABS(NOUTE).GT.4) THEN
        print *, "The value of NOUTE is not allowed"
        stop
      ENDIF
C
      READ(15,80) NSTAEMSG
      READ(NSTAEMSG,*) NSTAE
      IF (NSTAE.GT.MNSTAE) THEN
         print *, "NSTAE = ",NSTAE, " exceeds parameter MNSTAE = ",
     &              MNSTAE
          stop   
      ENDIF
#ifdef NETCDF
C
C     ==================================================================
C     START ADDITION BY MCF FOR NETCDF OPTION
C
!      IF(.NOT. ALLOCATED(STATNUMB))ALLOCATE(STATNUMB(NSTAE,SNUMLEN))
      IF(.NOT. ALLOCATED(STATNAME))ALLOCATE(STATNAME(NSTAE))
C
C     END ADDITION BY MCF FOR NETCDF OPTION
C     ==================================================================
#endif
      DO I=1,NSTAE
         READ(15,'(A132)') STAELOC(I)
         IF(ICS.EQ.1) THEN

C     ==================================================================
C     START ADDITION BY MCF FOR NETCDF OPTION
C
              IF(NOUTE.NE.3)THEN
                 READ(STAELOC(I),*) XEL(I),YEL(I)
              ELSE
#ifdef NETCDF
!               READ(15,'(A132)') STATLINE
!               call parse(STATLINE, LVAR)
                call parse(STAELOC(I), LVAR)
                XEL(I)=a2f(LVAR(1))
                YEL(I)=a2f(LVAR(2))
                STATNAME(I)=LVAR(3)
                WRITE(6,9911) XEL(I), YEL(I), STATNAME(I)
#endif
              ENDIF
            
         ELSE
              IF(NOUTE.NE.3)THEN

                 READ(STAELOC(I),*) SLEL(I),SFEL(I)
              ELSE
#ifdef NETCDF
!                READ(15,'(A132)') STATLINE
!                call parse(STATLINE, LVAR)
                 call parse(STAELOC(I), LVAR)
                 SLEL(I)=a2f(LVAR(1))
                 SFEL(I)=a2f(LVAR(2))
                 STATNAME(I)=LVAR(3)
!                 WRITE(6,*)"parsed stations"
!                 WRITE(6,9911) SLEL(I),SFEL(I), STATNAME(I,SNAMLEN)
#endif
             ENDIF

            
            SL1Sta=DEG2RAD*SLEL(I)
            SF1Sta=DEG2RAD*SFEL(I)
c           CALL CPP(XEL(I),YEL(I),SL1,SF1,SL0,SF0)
            XEL(I) = R*(SL1Sta-SL0)*COS(SF0)
            YEL(I) = R*SF1Sta
         ENDIF
      ENDDO
 9911            FORMAT(F12.3,2X,F12.3, 6X, A50)

!      DO I=1,NSTAE
!      WRITE(6,*)" XEL ===", XEL(I),YEL(I)
!      END DO
!      DO I=1,NSTAE
!      WRITE(6,*)" SLEL(I),SFEL(I) ===",SLEL(I),SFEL(I)
!      END DO

C
C--For Each Elevation Station: 
C  Find the Global Index of the element it lies in.
C
      CALL CoordToEle(XEL,YEL,NNSEG,NSTAE,
     &     'Elevation recording station   ')
C
C--Read Velocity Recording Stations
C
      READ(15,80) STAVMSG
      READ(STAVMSG,*) NOUTV,TOUTSV,TOUTFV,NSPOOLV
      IF (ABS(NOUTV).GT.4) THEN
        print *, "Value for NOUTV is not allowable"
        stop
      ENDIF
C
      READ(15,80) NSTAVMSG
      READ(NSTAVMSG,*) NSTAV
      IF (NSTAV.GT.MNSTAV) THEN
        print *, "NSTAV = ",NSTAV, " exceeds parameter MNSTAV = ",
     &             MNSTAV
        stop 
      ENDIF
#ifdef NETCDF
C     ==================================================================
C     START ADDITION BY MCF FOR NETCDF OPTION
C
!      IF(.NOT. ALLOCATED(STATNUMBV))ALLOCATE(STATNUMBV(NSTAV,SNUMLEN))
      IF(.NOT. ALLOCATED(STATNAMEV))ALLOCATE(STATNAMEV(NSTAV))
C
C     END ADDITION BY MCF FOR NETCDF OPTION
C     ==================================================================
C
#endif
      IF (ICS.EQ.1) THEN
      

        DO I=1,NSTAV
           READ(15,'(A132)') STAVLOC(I)
C     ==================================================================
C     START ADDITION BY MCF FOR NETCDF OPTION
C
           IF(NOUTV.NE.3) THEN
              READ(STAVLOC(I),*) XEV(I),YEV(I)
           ELSE
#ifdef NETCDF
              call parse(STAVLOC(I), LVAR)
              XEV(I)=a2f(LVAR(1))
              YEV(I)=a2f(LVAR(2))
              STATNAMEV(I)=LVAR(3)
              WRITE(6,9911) XEV(I), YEV(I), STATNAMEV(I)
#endif
           ENDIF

        ENDDO
        
      ELSE
        DO I=1,NSTAV
           READ(15,'(A132)') STAVLOC(I)
           
           IF(NOUTV.NE.3) THEN
           
               READ(STAVLOC(I),*) SLEV(I),SFEV(I)

           ELSE
#ifdef NETCDF
              call parse(STAVLOC(I), LVAR)
              SLEV(I)=a2f(LVAR(1))
              SFEV(I)=a2f(LVAR(2))
              STATNAMEV(I)=LVAR(3)
!              WRITE(6,*)"parsed stations"
!              WRITE(6,9911) SLEV(I),SFEV(I), STATNAMEV(I,SNAMLEN)
!              WRITE(6,*)"STATION INFO - AFTER *************",
!     &           SLEV(I),SFEV(I), STATNAMEV(I,SNAMLEN)
#endif
         ENDIF
C
C     END ADDITION BY MCF FOR NETCDF OPTION
C     ==================================================================

           
           SL1Sta=DEG2RAD*SLEV(I)
           SF1Sta=DEG2RAD*SFEV(I)
C          CALL CPP(XEV(I),YEV(I),SL1,SF1,SL0,SF0)
           XEV(I) = R*(SL1Sta-SL0)*COS(SF0)
           YEV(I) = R*SF1Sta
        ENDDO
      ENDIF
C
C--For Each Velocity Station: 
C  Find the Global Index of the element it lies in.
C
      CALL CoordToEle(XEV,YEV,NNSVG,NSTAV,
     &     'Velocity recording station    ')
C
C     If Passive Transport is indicated, then read Concentration Station Info 
C
      IF (C2D_PTrans.or.C3D_PTrans) THEN
C
        READ(15,80) STACMSG
        READ(STACMSG,*) NOUTC,TOUTSC,TOUTFC,NSPOOLC
        IF (ABS(NOUTC).GT.2) THEN
          print *, "Value of NOUTC is not allowable"
          stop
        ENDIF
C
        READ(15,80) NSTACMSG
        READ(NSTACMSG,*) NSTAC
        IF (NSTAC.GT.MNSTAC) THEN
          print *, "NSTAC = ",NSTAC, " exceeds parameter MNSTAC = ",
     &              MNSTAC
          stop 
        ENDIF
C
        DO I=1,NSTAC
           IF (ICS.EQ.1) THEN
             READ(15,'(A132)') STACLOC(I)
             READ(STACLOC(I),*) XEC(I),YEC(I)
           ELSE
             READ(15,'(A132)') STACLOC(I)
             READ(STACLOC(I),*) SLEC(I),SFEC(I)
             SL1Sta=DEG2RAD*SLEC(I)
             SF1Sta=DEG2RAD*SFEC(I)
c            CALL CPP(XEC(I),YEC(I),SL1,SF1,SL0,SF0)
             XEC(I) = R*(SL1Sta-SL0)*COS(SF0)
             YEC(I) = R*SF1Sta
           ENDIF
        ENDDO
C
C--For Each Concentration Recording Station: 
C  Find the Global Index of the element it lies in.
C
         CALL CoordToEle(XEC,YEC,NNSCG,NSTAC,
     &        'Concentration station         ')
C
      ENDIF
C
C--If NWS <> 0 , then read Meteorlogical Station Info 
C
      NOUTM = 0
      IF (NWS.NE.0) THEN
C
         READ(15,80) STAMMSG
         READ(STAMMSG,*) NOUTM,TOUTSM,TOUTFM,NSPOOLM
         IF (ABS(NOUTM).GT.3) THEN
            print *, "Value of NOUTM is not allowable"
            stop
         ENDIF
C     
         READ(15,80) NSTAMMSG
         READ(NSTAMMSG,*) NSTAM
         IF (NSTAM.GT.MNSTAM) THEN
            print *, "NSTAM = ",NSTAM, " exceeds parameter MNSTAM = ",
     &           MNSTAM
            stop 
         ENDIF
#ifdef NETCDF
C     ==================================================================
C     START ADDITION BY MCF FOR NETCDF OPTION
C
!         IF(.NOT.ALLOCATED(STATNUMBM))ALLOCATE(STATNUMBM(NSTAM,SNUMLEN))
         IF(.NOT.ALLOCATED(STATNAMEM))ALLOCATE(STATNAMEM(NSTAM))
C
C     END ADDITION BY MCF FOR NETCDF OPTION
C     ==================================================================
C     
#endif
         DO I=1,NSTAM
            READ(15,'(A132)') STAMLOC(I)
!            READ(15,'(A132)') STAMLOC(I)
!            WRITE(6,*)"STAMLOC(I) ================",I,STAMLOC(I)
            IF (ICS.EQ.1) THEN
            
C     ==================================================================
C     START ADDITION BY MCF FOR NETCDF OPTION
C
             IF(NOUTM.NE.3)THEN
            
               READ(STAMLOC(I),*) XEM(I),YEM(I)
             ELSE
#ifdef NETCDF
                call parse(STAMLOC(I), LVAR)
                XEM(I)=a2f(LVAR(1))
                YEM(I)=a2f(LVAR(2))
                STATNAMEM(I)=LVAR(3)
                WRITE(6,9911) XEM(I), YEM(I), STATNAMEM(I)
#endif
             ENDIF
             
         ELSE

             IF(NOUTM.NE.3)THEN

!              READ(15,80) STAMLOC(I)
               READ(STAMLOC(I),*) SLEM(I),SFEM(I)
               
             ELSE
#ifdef NETCDF

               call parse(STAMLOC(I), LVAR)
               SLEM(I)=a2f(LVAR(1))
               SFEM(I)=a2f(LVAR(2))
               STATNAMEM(I)=LVAR(3)
!               WRITE(6,9911) SLEM(I), SFEM(I), STATNAMEM(I,SNAMLEN)
!               WRITE(6,*)"STATION INFO STAMLOC",
!     &             SLEM(I),SFEM(I),STATNAMEM(I,SNAMLEN)
!    &             SLEM(I),SFEM(I),STATNUMBM(I),STATNAMEM(I)
#endif
             ENDIF
C
C     END ADDITION BY MCF FOR NETCDF OPTION
C     ==================================================================

               SL1Sta=DEG2RAD*SLEM(I)
               SF1Sta=DEG2RAD*SFEM(I)
c     CALL CPP(XEM(I),YEM(I),SL1,SF1,SL0,SF0)
               XEM(I) = R*(SL1Sta-SL0)*COS(SF0)
               YEM(I) = R*SF1Sta
            ENDIF
         ENDDO
        
C
C--For Each Meterological Recording Station: 
C  Find the Global Index of the element it lies in.
C
         CALL CoordToEle(XEM,YEM,NNSMG,NSTAM,
     &        'Meteorological station        ')
C     
      ENDIF
C     
C
C--Read Global Elevation Data Output
C
      READ(15,80) OUTGEMSG
      READ(OUTGEMSG,*) NOUTGE,TOUTSGE,TOUTFGE,NSPOOLGE
      IF (ABS(NOUTGE).GT.4) THEN
         print *, "NOUTGE does not have an allowable value"
         stop
      ENDIF
C
C--Read Global Velocity Data Output
C
      READ(15,80) OUTGVMSG
      READ(OUTGVMSG,*) NOUTGV,TOUTSGV,TOUTFGV,NSPOOLGV
      IF (ABS(NOUTGV).GT.4) THEN
        print *, "NOUTGV does not have an allowable value"
        stop
      ENDIF
C
C     If Passive Transport is indicated, read Global Concentration Data Output
C
      IF (C2D_PTrans.or.C3D_PTrans) THEN
        READ(15,80) OUTGCMSG
        READ(OUTGCMSG,*) NOUTGC,TOUTSGC,TOUTFGC,NSPOOLGC
        IF (ABS(NOUTGC).GT.2) THEN
          print *, "NOUTGC does not have an allowable value"
          stop
        ENDIF
      ENDIF
C
      IF (NWS.NE.0) THEN
        READ(15,80) OUTGWMSG
        READ(OUTGWMSG,*) NOUTGW,TOUTSGW,TOUTFGW,NSPOOLGW
Casey 090302: Changed the following line for ABS(NOUTGW).EQ.4
        IF (ABS(NOUTGW).GT.3 .AND. ABS(NOUTGW).NE.4) THEN
          print *, "NOUTGW does not have an allowable value"
          stop
        ENDIF
      ENDIF
C
C--Read Harmonic Analysis Data
C
      READ(15,80) HARFRMSG
      READ(HARFRMSG,*) NHARFR
      IF (NHARFR.GT.MNHARF) THEN
        print *, "NHARFR exceeds parameter MNHARF"
        print *, "NHARFR = ",NHARFR, " exceeds parameter MNHARF = ",
     &            MNHARF
        stop 
      ENDIF
C
      DO I=1,NHARFR
          READ(15,80) HAFNAM(I)
          READ(15,80) HAFREMSG(I)
          READ(HAFREMSG(I),*) HAFREQ(I),HAFF(I),HAFACE(I)
      ENDDO
C
      READ(15,80) HARPARMSG
      READ(HARPARMSG,*) THAS,THAF,NHAINC,FMV
      READ(15,80) OUTHARMSG
      READ(OUTHARMSG,*) NHASE,NHASV,NHAGE,NHAGV
      IF ((NHASE.LT.0).OR.(NHASE.GT.1)) THEN
        print *, "NHASE does not have an allowable value"
        IF (NFOVER.EQ.1) THEN
          print *, "Program will override an reset NHASE=0 "
          NHASE = 0
        ELSE
          stop
        ENDIF
      ENDIF
C
      IF ((NHASV.LT.0).OR.(NHASV.GT.1)) THEN
         print *, "NHASV does not have an allowable value"
         IF (NFOVER.EQ.1) THEN
           print *, "Program will override an reset NHASV=0 "
           NHASV = 0
         ELSE
           stop
         ENDIF
      ENDIF
C
      IF ((NHAGE.LT.0).OR.(NHAGE.GT.1)) THEN
        print *, "NHAGE does not have an allowable value"
        IF (NFOVER.EQ.1) THEN
          print *, "Program will override an reset NHAGE=0 "
          NHAGE = 0
        ELSE
          stop
        ENDIF
      ENDIF
C
      IF ((NHAGV.LT.0).OR.(NHAGV.GT.1)) THEN
        print *, "NHAGV does not have an allowable value"
        IF (NFOVER.EQ.1) THEN
          print *, "Program will override an reset NHAGV=0 "
          NHAGV = 0
        ELSE
          stop
        ENDIF
      ENDIF
C
C--Read Hot Start Data
C
C     jgf45.07 added undocumented option to allow ADCIRC to stop after writing
C     hot start file. This is used to test hot starting capabilities.
      READ(15,80) HSTARMSG
      READ(HSTARMSG,*) NHSTAR,NHSINC
!      PRINT *, " NHSTAR,NHSINC ", NHSTAR,NHSINC
      IF ((NHSTAR.NE.0).AND.(NHSTAR.NE.1).AND.(NHSTAR.NE.67)
     &     .AND.(NHSTAR.NE.68).AND.(NHSTAR.NE.3)) THEN 
        print *, "NHSTAR does not have an allowable value"
        IF (NFOVER.EQ.1) THEN
          print *, "Program will override and reset NHSTAR=0 "
          NHSTAR = 0
        ELSE
          stop
        ENDIF
      ENDIF
!      PRINT *, " NHSTAR,NHSINC ", NHSTAR,NHSINC
C
C--Read Solver Data
C
      READ(15,80) SOLVMSG
      READ(SOLVMSG,*) ITITER,ISLDIA,CONVCR,ITMAX
      IF ((ITITER.NE.1).AND.(ITITER.NE.-1)) THEN
         print *, "WARNING: An invalid solver was selected, ITITER=",
     &            ITITER
         print *, "The JCG iterative solver (ITITER=1) and the explicit"
         print *, "direct solve (ITITER=-1) are the only valid choices."
         IF (NFOVER.EQ.1) THEN
           print *, "WARNING: ADCIRC will continue, using JCG."
           ITITER = 1
         ELSE
           print *, "ERROR: ITITER must be set to 1 or -1."
           stop
         ENDIF
      ENDIF
C
C
C--Read in 3D data
C
C      print *, " made it to the end of the 2D read of fort.15"
      IF(C3DVS) THEN
        CALL READ15_3DVS()
c     ELSEIF(C3DDSS) THEN
c       CALL READ15_3DDSS()
      ENDIF

#ifdef NETCDF
C     ==================================================================
C     START ADDITION BY MCF FOR NETCDF OPTION for NETCDF I/O - 6/20/07
C
!     --------------------------------------------------------------------
!     IF output is in netCDF format read global attributes for netCDF file
!     --------------------------------------------------------------------

      IF(useNetCDF.eqv..true.) THEN
#ifndef NETCDF
         print *, "NetCDF was selected but is not supported by your
     &            compilation flag.
     &
     &            You need to compile with:
     &
     &  gmake NETCDF=enable NETCDFHOME=the netcdf library location"
         stop
#endif
         READ(15,'(A80)') title
         READ(15,'(A80)') institution
         READ(15,'(A80)') source
         READ(15,'(A80)') history
         READ(15,'(A80)') references
         READ(15,'(A80)') comments
         READ(15,'(A80)') host
         READ(15,'(A80)') convention
         READ(15,'(A80)') contact
         READ(15,'(A80)') base_date

        ELSEIF(NOUTGC.eq.3.                                 ! fort.93
     &       or.NOUTC.eq.3.                                 ! fort.91
     &       or.NHASE.eq.3.                                 ! fort.51
     &       or.NHASV.eq.3.                                 ! fort.52
     &       or.NHAGE.eq.3.                                 ! fort.53
     &       or.NHAGV.eq.3.                                 ! fort.54
     &       or.I3DSD.eq.3.                                 ! fort.41
     &       or.I3DSV.eq.3.                                 ! fort.42
     &       or.I3DST.eq.3.                                 ! fort.43
     &       or.I3DGD.eq.3.                                 ! fort.44
     &       or.I3DGV.eq.3.                                 ! fort.45
     &       or.I3DGT.eq.3) THEN                            ! fort.46

        WRITE(6,*)
     &  "netCDF output not available yet for the options selected"

      ENDIF
#endif



C
C--Close Global Run Info file 
C
      CLOSE(15)
      call memory_status()
      RETURN
 80   FORMAT(A80)
      END


      SUBROUTINE READ15_3DVS() 
      USE PRE_GLOBAL 
      use memory_usage
C
C----------------------------------------------------------------------
C     Reads the 3DVS portion of the global ADCIRC Input Data File for
C     data decomposition program ADCPP.  This version compatible with
C     ADCIRC_v41.11a
C     tjc  6/24/02
C---------------------------------------------------------------------
C
      IMPLICIT NONE
      integer :: nbytes = 0
      REAL(8) SL1Sta,SF1Sta
      INTEGER I,J,K,M,N                                        ! loop counters
      
      REAL(8) :: STATIME
      REAL(8), PARAMETER :: day2sec=24.d0*3600.d0

      SL0=DEG2RAD*SLAM0 !jgf45.12
      SF0=DEG2RAD*SFEA0

C
 350  FORMAT(//,2X,'***** INVALID INPUT IN THE PRIMARY VERTICAL INPUT',
     &          ' FILE (UNIT 15) ****',/,'****** RUN TERMINATED ******')

C     jgf45.10 removed IDIAG
C... SPECIFY WHETHER A BAROTROPIC OR BAROCLINIC RUN

      READ(15,80) IDENMSG
      READ(IDENMSG,*) IDEN
      IF((IDEN.GT.4).OR.(IDEN.LT.-4)) THEN
         WRITE(*,*) "ERROR: IDEN=",IDEN
         WRITE(*,423) 
         WRITE(*,350)
 423     FORMAT(/,2X,' IDEN must be an integer between -4 and 4.')
         STOP
      ENDIF

C     jgf45.12 Set flag to read in unit 11 file (initial density) if
C     necessary.
      IF ( IDEN .ne. 0 ) THEN
         CBaroclinic = .true.
      ENDIF
C     jgf45.12 Set flag to read in additional parameters from unit 15
C     file that will be used when solving for transport of salinity,
C     temperature, etc in prognostic baroclinic simulation.
      IF ( IDEN .gt. 0 ) THEN
         C3D_BTrans = .true.
      ENDIF

C... READ IN THE TYPE OF BOTTOM BOUNDARY CONDITION AND THE SLIP COEFFICIENTS

      READ(15,80) SLIPMSG
      READ(SLIPMSG,*) ISLIP,KP
      IF((ISLIP.LT.0).OR.(ISLIP.GT.3)) THEN
         WRITE(6,350)
         WRITE(6,360)
 360     FORMAT(/,2X,'    THE SLIP CODE MUST = 0,1,2,OR 3.')
         STOP
      ENDIF

C... READ IN THE SURFACE AND BOTTOM ROUGHNESSES

      READ(15,80) Z0MSG
      READ(Z0MSG,*) Z0S, Z0B

C... READ IN THE TIME STEPPING COEFFICIENTS

      READ(15,80) ALPMSG
      READ(ALPMSG,*) ALP1,ALP2,ALP3

C... READ IN IGC & NFEN: F.E. GRID CODE & # NODES IN F.E. GRID

      READ(15,80) FEMSG
      READ(FEMSG,*) IGC,NFEN

C     jgf45.12 Add code to read in nondimensional thicknesses of
C     vertical layers.
      IF(IGC.EQ.0) THEN
         ALLOCATE ( Sigma(NFEN) )
         nbytes = 8*nfen
         call memory_alloc(nbytes)
         DO N=1,NFEN
            READ(15,*) Sigma(N)
         ENDDO
      ENDIF

C... SPECIFY TYPE OF EDDY VISCOSITY PROFILE

      READ(15,80) EVCMSG
      READ(EVCMSG,*) IEVC,EVMIN,EVCON
      IF((IEVC.NE.0 ).and.(IEVC.NE.1 ).and.
     &   (IEVC.NE.10).and.(IEVC.NE.11).and.
     &   (IEVC.NE.20).and.(IEVC.NE.21).and.
     &   (IEVC.NE.22).and.(IEVC.NE.23).and.
     &   (IEVC.NE.30).and.(IEVC.NE.31).and.(IEVC.NE.32).and.
     &   (IEVC.NE.33).and.
     &   (IEVC.NE.40).and.(IEVC.NE.41).and.(IEVC.NE.42).and.
     &   (IEVC.NE.43).and.
     &   (IEVC.NE.50)) THEN
         WRITE(*,350)
         WRITE(*,411)
 411     FORMAT(/,2X,'    IEVC MUST BE 0,1,10,11,20,21,22,23,',
     &              '30,31,32,33,40,41,42,43,50')
         STOP
      ENDIF
      IF(IEVC.EQ.50) THEN
         READ(15,80) THETAMSG
         READ(THETAMSG,*) THETA1,THETA2
      ENDIF
C     jgf45.12 Add code to read in vertical eddy viscosity profile.
      IF(IEVC.EQ.0) THEN
         ALLOCATE ( EVTot(NFEN) )
         nbytes = 8*nfen
         call memory_alloc(nbytes)
         DO N=1,NFEN
            READ(15,*) EVTot(N)
         ENDDO
      ENDIF
C     -----------------------------------------------------
C.... Station 3D Density, Temperature, Salinity output
C     jgf45.11 switched from node numbers to coordinates

      READ(15,80) DSDMSG
      READ(DSDMSG,*) I3DSD,TO3DSDS,TO3DSDF,NSPO3DSD
      READ(15,80) NSTA3DDMSG
      READ(NSTA3DDMSG,*) NSTA3DD
C  kmd48.33bc changed to -2 and I3DSD not equal to 0
      IF((I3DSD.LT.-2).OR.(I3DSD.GT.2)) THEN
         WRITE(*,350)
         WRITE(*,511)
 511     FORMAT(/,2X,'    I3DSD MUST BE -2,-1,0,1 OR 2')
         STOP
      ENDIF
      IF(I3DSD.NE.0) THEN
         ALLOCATE(NNSTA3DDP(MNPROC))
         nbytes = 4*mnproc
         call memory_alloc(nbytes)
         ALLOCATE(IMAP_STA3DD_LG(NSTA3DD,MNPROC))
         nbytes = 4*nsta3dd*mnproc
         call memory_alloc(nbytes)
         ALLOCATE(STA3DDLOC(NSTA3DD))
         nbytes = 8*nsta3dd
         call memory_alloc(nbytes)
         ALLOCATE(X3DD(NSTA3DD),Y3DD(NSTA3DD))
         nbytes = 16*nsta3dd
         call memory_alloc(nbytes)
         IF (ICS.EQ.1) THEN
            DO I=1,NSTA3DD
               READ(15,80) STA3DDLOC(I)
               READ(STA3DDLOC(I),*) X3DD(I),Y3DD(I)
            ENDDO
         ELSE
            ALLOCATE(SL3DD(NSTA3DD),SF3DD(NSTA3DD))
            nbytes = 16*nsta3dd
            call memory_alloc(nbytes)
            DO I=1,NSTA3DD
               READ(15,80) STA3DDLOC(I)
               READ(STA3DDLOC(I),*) SL3DD(I),SF3DD(I)
               SL1Sta=DEG2RAD*SL3DD(I)
               SF1Sta=DEG2RAD*SF3DD(I)
               X3DD(I) = R*(SL1Sta-SL0)*COS(SF0)
               Y3DD(I) = R*SF1Sta
            ENDDO
         ENDIF

C     jgf45.11 For each 3D density station, find the full-domain index of 
C     the corresponding element.

         ALLOCATE(NNS3DDG(NSTA3DD))
         nbytes = 4*nsta3dd
         call memory_alloc(nbytes)
         CALL CoordToEle(X3DD,Y3DD,NNS3DDG,NSTA3DD,
     &        '3D density recording station  ')

      ENDIF

C     -----------------------------------------------------
C.... Station 3D Velocity output
C     jgf45.11 switched from node numbers to coordinates

      READ(15,80) DSVMSG
      READ(DSVMSG,*) I3DSV,TO3DSVS,TO3DSVF,NSPO3DSV
      READ(15,80) NSTA3DVMSG
      READ(NSTA3DVMSG,*) NSTA3DV
C  kmd48.33 changed to -2
      IF((I3DSV.LT.-2).OR.(I3DSV.GT.2)) THEN
         WRITE(*,350)
         WRITE(*,512)
 512     FORMAT(/,2X,'    I3DSV MUST BE -2,-1,0,1 OR 2')
         STOP
      ENDIF
      IF(I3DSV.NE.0) THEN
         ALLOCATE(NNSTA3DVP(MNPROC))
         nbytes = 4*mnproc  
         call memory_alloc(nbytes)
         ALLOCATE(IMAP_STA3DV_LG(NSTA3DV,MNPROC))
         nbytes = 4*nsta3dv*mnproc  
         call memory_alloc(nbytes)
         ALLOCATE(STA3DVLOC(NSTA3DV))
         nbytes = 8*nsta3dv
         call memory_alloc(nbytes)
         ALLOCATE(X3DV(NSTA3DV),Y3DV(NSTA3DV))
         nbytes = 16*nsta3dv
         call memory_alloc(nbytes)
         IF (ICS.EQ.1) THEN
            DO I=1,NSTA3DV
               READ(15,80) STA3DVLOC(I)
               READ(STA3DVLOC(I),*) X3DV(I),Y3DV(I)
            ENDDO
         ELSE
            ALLOCATE(SL3DV(NSTA3DV),SF3DV(NSTA3DV))
            nbytes = 16*nsta3dv
            call memory_alloc(nbytes)
            DO I=1,NSTA3DV
               READ(15,80) STA3DVLOC(I)
               READ(STA3DVLOC(I),*) SL3DV(I),SF3DV(I)
               SL1Sta=DEG2RAD*SL3DV(I)
               SF1Sta=DEG2RAD*SF3DV(I)
               X3DV(I) = R*(SL1Sta-SL0)*COS(SF0)
               Y3DV(I) = R*SF1Sta
            ENDDO
         ENDIF

C     jgf45.11 For each 3D velocity station, find the full-domain index of 
C     the corresponding element.

         ALLOCATE(NNS3DVG(NSTA3DV))
         nbytes = 4*nsta3dv
         call memory_alloc(nbytes)
         CALL CoordToEle(X3DV,Y3DV,NNS3DVG,NSTA3DV,
     &        '3D velocity recording station ')

      ENDIF
C     -----------------------------------------------------
C.... Station 3D turbulence output
C     jgf45.11 switched from node numbers to coordinates

      READ(15,80) DSTMSG
      READ(DSTMSG,*) I3DST,TO3DSTS,TO3DSTF,NSPO3DST
      READ(15,80) NSTA3DTMSG
      READ(NSTA3DTMSG,*) NSTA3DT
C  kmd48.33 changed to -2
      IF((I3DST.LT.-2).OR.(I3DST.GT.2)) THEN
         WRITE(*,350)
         WRITE(*,513)
 513     FORMAT(/,2X,'    I3DST MUST BE -2,-1,0,1 OR 2')
         STOP
      ENDIF
      IF(I3DST.NE.0) THEN
         ALLOCATE(NNSTA3DTP(MNPROC))
         nbytes = 4*mnproc 
         call memory_alloc(nbytes)
         ALLOCATE(IMAP_STA3DT_LG(NSTA3DT,MNPROC))
         nbytes = 4*nsta3dt*mnproc 
         call memory_alloc(nbytes)
         ALLOCATE(STA3DTLOC(NSTA3DT))
         nbytes = 4*nsta3dt
         call memory_alloc(nbytes)
         ALLOCATE(X3DT(NSTA3DT),Y3DT(NSTA3DT))
         nbytes = 16*nsta3dt
         call memory_alloc(nbytes)
         IF (ICS.EQ.1) THEN
            DO I=1,NSTA3DT
               READ(15,80) STA3DTLOC(I)
               READ(STA3DTLOC(I),*) X3DT(I),Y3DT(I)
            ENDDO
         ELSE
            ALLOCATE(SL3DT(NSTA3DT),SF3DT(NSTA3DT))
            nbytes = 16*nsta3dt
            call memory_alloc(nbytes)
            DO I=1,NSTA3DT
               READ(15,80) STA3DTLOC(I)
               READ(STA3DTLOC(I),*) SL3DT(I),SF3DT(I)
               SL1Sta=DEG2RAD*SL3DT(I)
               SF1Sta=DEG2RAD*SF3DT(I)
               X3DT(I) = R*(SL1Sta-SL0)*COS(SF0)
               Y3DT(I) = R*SF1Sta
            ENDDO
         ENDIF

C     jgf45.11 For each 3D turbulence station, find the full-domain index of 
C     the corresponding element.

         ALLOCATE(NNS3DTG(NSTA3DT))
         nbytes = 4*nsta3dt
         call memory_alloc(nbytes)
         CALL CoordToEle(X3DT,Y3DT,NNS3DTG,NSTA3DT,
     &        '3D turbulence rec. station    ')

      ENDIF

C....  GLOBAL 3D DENSITY, TEMPERATURE, SALINITY OUTPUT

      READ(15,80) DGDMSG
      READ(DGDMSG,*) I3DGD,TO3DGDS,TO3DGDF,NSPO3DGD
C   kmd48.33bc changed to -2
      IF((I3DGD.LT.-2).OR.(I3DGD.GT.2)) THEN
         WRITE(*,350)
         WRITE(*,514)
 514     FORMAT(/,2X,'    I3DGD MUST BE -2,-1,0,1 OR 2')
         STOP
      ENDIF

C....  GLOBAL 3D VELOCITY OUTPUT

      READ(15,80) DGVMSG
      READ(DGVMSG,*) I3DGV,TO3DGVS,TO3DGVF,NSPO3DGV
C   kmd48.33bc changed to -2
      IF((I3DGV.LT.-2).OR.(I3DGV.GT.2)) THEN
         WRITE(*,350)
         WRITE(*,515)
 515     FORMAT(/,2X,'    I3DGV MUST BE -2,-1,0,1 OR 2')
         STOP
      ENDIF

C....  GLOBAL 3D TURBULENCE OUTPUT

      READ(15,80) DGTMSG
      READ(DGTMSG,*) I3DGT,TO3DGTS,TO3DGTF,NSPO3DGT
C   kmd48.33bc changed to -2
      IF((I3DGT.LT.-2).OR.(I3DGT.GT.2)) THEN
         WRITE(*,350)
         WRITE(16,516)
 516     FORMAT(/,2X,'    I3DGT MUST BE -2,-1,0,1 OR 2')
         STOP
      ENDIF

C   kmd48.33bc add in information for 3D boundary condition
      IF (CBAROCLINIC) THEN
        READ(15,80) RESBCFLAGMSG
        READ(RESBCFLAGMSG,*) RES_BC_FLAG, BCFLAG_LNM, BCFLAG_TEMP
        IF (RES_BC_FLAG.NE.0) THEN
           READ(15,80) BCTIMEMSG
           READ(15,80) BCSTATMSG
           IF (BCFLAG_TEMP.NE.0) THEN
              READ(15,80) TBCTIMEMSG
           END IF
        END IF
        IF ((RES_BC_FLAG.LT.0).OR.(RES_BC_FLAG.EQ.1))  THEN
           READ(BCTIMEMSG,*) RBCTIMEINC
           READ(BCSTATMSG,*) BCSTATIM
        ELSE IF (RES_BC_FLAG.EQ.2) THEN
           READ(BCTIMEMSG,*) RBCTIMEINC, SBCTIMEINC
           READ(BCSTATMSG,*) BCSTATIM, SBCSTATIM
        ELSE IF (RES_BC_FLAG.EQ.3) THEN
           READ(BCTIMEMSG,*) RBCTIMEINC, TBCTIMEINC
           READ(BCSTATMSG,*) BCSTATIM, TBCSTATIM
           IF (BCFLAG_TEMP.NE.0) THEN
              READ(TBCTIMEMSG,*) TTBCTIMEINC, TTBCSTATIM
           END IF
        ELSE IF (RES_BC_FLAG.EQ.4) THEN
           READ(BCTIMEMSG,*) RBCTIMEINC, SBCTIMEINC, TBCTIMEINC
           READ(BCSTATMSG,*) BCSTATIM, SBCSTATIM, TBCSTATIM
           IF (BCFLAG_TEMP.NE.0) THEN
              READ(TBCTIMEMSG,*) TTBCTIMEINC, TTBCSTATIM
           END IF
        END IF
      END IF

      IF (CBAROCLINIC) THEN
        READ(15,80) SPONGEDISTMSG
        READ(15,80) EqnstateMSG
      END IF

C
C     jgf45.12: Read in the parameters for the transport equation, if
C     necessary.
      IF (C3D_BTrans) THEN

C     Lateral and vertical diffusion coefficients.
         READ(15,*) NLSD, NVSD
         READ(15,*) NLTD, NVTD
C     Time stepping coefficient for the transport equation terms.
         READ(15,*) ALP4
C     Temperature boundary condition file type, if necessary
!         if ( IDEN .eq. 3 .or. IDEN .eq. 4 ) then
!            READ(15,*) NTF
!         endif
         endif

      call memory_status()
      RETURN
C
 80   FORMAT(A80)
      END SUBROUTINE



*******************************************************************************
C                                                                             *
C    Transform from lon,lat (lamda,phi) coordinates into CPP coordinates.     *
C    Lon,Lat must be in radians.                                              *
C                                                                             *
C******************************************************************************

      SUBROUTINE CPP(X,Y,RLAMBDA,PHI,RLAMBDA0,PHI0)
      IMPLICIT NONE
      REAL(8) X, Y, PHI, RLAMBDA
      REAL(8) RLAMBDA0, PHI0, R
C
      R=6378206.4
      X=R*(RLAMBDA-RLAMBDA0)*COS(PHI0)
      Y=PHI*R
      RETURN
      END


C******************************************************************************
C                                                                             *
C    Transform from CPP coordinates to lon,lat (lamda,phi) coordinates        *
C    Lon,Lat is in radians.                                                   *
C                                                                             *
C******************************************************************************

      SUBROUTINE INVCP(XXCP,YYCP,RLAMBDA,PHI,RLAMBDA0,PHI0)
      IMPLICIT NONE
      REAL(8) XXCP, YYCP, PHI, RLAMBDA
      REAL(8) RLAMBDA0, PHI0, R
C
      R=6378206.4
      RLAMBDA=RLAMBDA0+XXCP/(R*COS(PHI0))
      PHI=YYCP/R
      RETURN
      END


C-----------------------------------------------------------------------
C          S U B R O U T I N E   C O O R D  T O  E L E 
C-----------------------------------------------------------------------
C
C     jgf45.11 Subroutine to take a list of X and Y cartesian
C     coordinates and find the corresponding elements.
C
C-----------------------------------------------------------------------
      SUBROUTINE CoordToEle(XCoords,YCoords,FDEle,NumOfStations,Desc)
      USE PRE_GLOBAL
      IMPLICIT NONE
      INTEGER, intent(in) :: NumOfStations                     ! total
      REAL(8), intent(in), dimension(NumOfStations) :: XCoords ! cartesian
      REAL(8), intent(in), dimension(NumOfStations) :: YCoords ! cartesian
      INTEGER, intent(out), dimension(NumOfStations):: FDEle   ! FullDomain
      CHARACTER(len=30), intent(in) :: Desc                    ! description

      INTEGER I,J,K,M,N                                        ! loop counters
      INTEGER N1, N2, N3, KMIN
      REAL(8) X1, X2, X3, X4, Y1, Y2, Y3, Y4, A1, A2, A3       ! geometry
      REAL(8) AE, AEMIN, AREASK, AA                            ! area
      LOGICAL ElementFound  ! .true. when a corresponding element is found

      DO I = 1, NumOfStations
         ElementFound = .false.
         AEMIN=1.0E+25
         KMIN=0
         DO K=1,NELG
            N1=NNEG(1,K)
            N2=NNEG(2,K)
            N3=NNEG(3,K)
            X1=SLAM(N1)
            X2=SLAM(N2)
            X3=SLAM(N3)
            X4=XCoords(I)
            Y1=SFEA(N1)
            Y2=SFEA(N2)
            Y3=SFEA(N3)
            Y4=YCoords(I)
            A1=(X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4)
            A2=(X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1)
            A3=(Y4-Y1)*(X2-X1)-(X4-X1)*(Y2-Y1)
            AA=ABS(A1)+ABS(A2)+ABS(A3)
            AREASK=X2*Y3+X1*Y2+X3*Y1-Y1*X2-Y2*X3-Y3*X1
            AE=ABS(AA-AREASK)/AREASK
            IF (AE.LT.AEMIN) THEN
               AEMIN=AE
               KMIN=K
            ENDIF
            IF (AE.LT.1.0E-5) THEN
               ElementFound = .true.
               FDEle(I)=K
            ENDIF
         ENDDO
         IF ( ElementFound .eqv. .false.) THEN
            write(*,1234) Desc, I
            print *, "Please check the coordinates."
            IF (NFOVER.EQ.1) THEN
               print *, "The program will estimate nearest element."
               print *, "WARNING. Proximity index is ",AEMIN
               print *, " "
               FDEle(I) = KMIN
            ELSE
               print *, "ERROR. Proximity index is ",AEMIN
               stop
            ENDIF
         ENDIF
      ENDDO

 1234 format(A30,1x,I4,1x,'does not lie in the grid.')

      RETURN
C-----------------------------------------------------------------------
      END SUBROUTINE CoordToEle
C-----------------------------------------------------------------------



