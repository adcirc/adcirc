C
C  Version 1.1  vjp  5/04/99
C

      SUBROUTINE COMPARE61(DIR1,DIR2,TOL)
C
C---------------------------------------------------------------------------C
C                     (  Serial Version  3/28/98  )                         C 
C  Compare data at the elevation stations from two different directories.   C
C  This version is compatible with ADCIRC version 33.04                     C
C                                                                           C
C  Input Parameters:                                                        C
C                                                                           C
C    DIR1  =  pathname of first directory                                   C
C    DIR2  =  pathname of second directory                                  C
C    TOL   =  user supplied relative error tolerance                        C
C                                                                           C
C---------------------------------------------------------------------------C
C
      USE POST_GLOBAL
      INTEGER I,J,K,L,IPROC,IDUM
      INTEGER NTRSPE1,NTRSPE2
      INTEGER NSTEMP1,NSTEMP2
      INTEGER ITSE1,ITSE2
      INTEGER ITEMPE1,ITEMPE2
      INTEGER NSTAE1,NSTAE2 
      INTEGER NUMSTNS1, NUMSTNS2
      INTEGER NSPOOLE1,NSPOOLE2 
      INTEGER IREC1,IREC2
      INTEGER LEN1,LEN2
      REAL(SZ) DTE1,DTE2
      REAL(8) TIMEOUTSE1,TIMEOUTSE2
      REAL(SZ) TOL,ERR !jgf45.07 ERR is a fortran keyword...
      CHARACTER*80  DIR1,DIR2,FNAME1,FNAME2,FNAME3
      CHARACTER*85  INLINE
      CHARACTER*4 RDES4(8),RID4(6),AID4(6)
      CHARACTER*8 RDES8(8),RID8(6),AID8(6)
      LOGICAL FOUND1,FOUND2

      REAL(8),ALLOCATABLE ::  ETA1(:),ETA2(:)
      ALLOCATE ( ETA1(MNSTAE),ETA2(MNSTAE))
C
C--Construct File Names
C
      DO I=1, LEN(DIR1)
         IF (DIR1(I:I).EQ.' ') THEN
           LEN1 = I-1
           EXIT
         ENDIF
      ENDDO
      DO I=1, LEN(DIR2)
         IF (DIR2(I:I).EQ.' ') THEN
           LEN2 = I-1
           EXIT
         ENDIF
      ENDDO
C
      FNAME1 = DIR1(1:LEN1)//'/fort.61'
      FNAME2 = DIR2(1:LEN2)//'/fort.61'
      FNAME3 = 'diffs.61'
C
C--Determine whether Unit 61 is Sequential Formatted or Direct Access Binary
C
      IF (ABS(NOUTE).EQ.1) THEN
        GO TO 1000
      ELSE
        GO TO 2000
      ENDIF
C
1000  CONTINUE
C
C--Open Both Global Sequential Formatted fort.61 files
C
      INQUIRE(FILE=FNAME1,EXIST=FOUND1)
      INQUIRE(FILE=FNAME2,EXIST=FOUND2)
      IF (FOUND1.AND.FOUND2) THEN
        OPEN(UNIT=51,FILE=FNAME1)
        OPEN(UNIT=61,FILE=FNAME2)
        OPEN(UNIT=71,FILE=FNAME3)
      ELSE
        print *, "No fort.61 files found"
        RETURN
      ENDIF
C
      READ (51,'(A85)') INLINE
c     WRITE(71,'(A85)')  INLINE
      READ (61,'(A85)') INLINE
c     WRITE(71,'(A85)')  INLINE
C
      READ (51,3645) NTRSPE1,NSTAE1,DTE1,NSPOOLE1,ITEMPE1
      READ (61,3645) NTRSPE2,NSTAE2,DTE2,NSPOOLE2,ITEMPE2
C
      IF (NTRSPE1.NE.NTRSPE2) THEN
         WRITE(71,*) "NTRSPE1 not equal to NSTRSPE2"
         RETURN
      ENDIF
C
      IF (NSTAE1.NE.NSTAE2) THEN
        WRITE(71,*) "NSTAE1 not equal to NSTAE2"
        RETURN
      ENDIF
C
      IF (ABS(DTE1-DTE2).GT.1.0E-5) THEN
        WRITE(71,*) "DTE1 not equal to DTE2"
        RETURN
      ENDIF
C
      IF (NSPOOLE1.NE.NSPOOLE2) THEN
        WRITE(71,*) "NSPOOLE1 not equal to NSPOOLE2"
        RETURN
      ENDIF
      IF (ITEMPE1.NE.ITEMPE2) THEN
        WRITE(71,*) "ITEMPE1 not equal to ITEMPE2"
        RETURN
      ENDIF
C
c     WRITE(71,*) NTRSPE1,NSTAE1,DTE1,NSPOOLE1,ITEMPE1
C
      DO J=1,NTRSPE1
C
         READ(51,2120) TIMEOUTSE1,ITSE1
         READ(61,2120) TIMEOUTSE2,ITSE2
C
         IF (ABS(TIMEOUTSE1-TIMEOUTSE2).GT.1.0E-5) THEN
           WRITE(71,*) "TIMEOUTSE1 not equal to TIMEOUTSE2"
           RETURN
         ENDIF
C
        IF (ITSE1.NE.ITSE2) THEN
          WRITE(71,*) "ITSE1 not equal to ITSE2"
          RETURN
        ENDIF
C
C       WRITE(71,*) TIMEOUTSE1,ITSE1
C
         DO I=1, NSTAE1
            READ(51,*) IDUM,ETA1(I)
            READ(61,*) IDUM,ETA2(I)
            ERR = ABS(ETA1(I)-ETA2(I))
            IF (ERR.GT.TOL) WRITE(71,*) TIMEOUTSE1,I,ETA1(I),ETA2(I)
         ENDDO
C
      ENDDO                 
      GO TO 9999
C 
2000  CONTINUE
C
C--Open Both Global Direct Access Binary fort.61 files
C
      INQUIRE(FILE=FNAME1,EXIST=FOUND1)
      INQUIRE(FILE=FNAME2,EXIST=FOUND2)
      IF (FOUND1.AND.FOUND2) THEN
        OPEN(51,FILE=FNAME1,ACCESS='DIRECT',RECL=NBYTE)
        OPEN(61,FILE=FNAME2,ACCESS='DIRECT',RECL=NBYTE)
        OPEN(UNIT=71,FILE=FNAME3)
      ELSE
        print *, "No fort.61 files found"
        RETURN
      ENDIF
C
C--Read RUNDES RUNID and AGRID from 1st Global File
C
      IREC1 = 0
      IF (NBYTE.EQ.4) THEN
        DO I=1,8
           READ(51,REC=IREC1+I) RDES4(I)
        ENDDO
        IREC1=IREC1+8
        DO I=1,6
           READ(51,REC=IREC1+I) RID4(I)
        ENDDO
        IREC1=IREC1+6
        DO I=1,6
           READ(51,REC=IREC1+I) AID4(I)
        ENDDO
        IREC1=IREC1+6
      ENDIF
      IF (NBYTE.EQ.8) THEN
        DO I=1,4
           READ(51,REC=IREC1+I) RDES8(I)
        ENDDO
        IREC1=IREC1+4
        DO I=1,3
           READ(51,REC=IREC1+I) RID8(I)
        ENDDO
        IREC1=IREC1+3
        DO I=1,3
           READ(51,REC=IREC1+I) AID8(I)
        ENDDO
        IREC1=IREC1+3
      ENDIF
C
C--Read RUNDES RUNID and AGRID from 2nd Global File
C
      IREC2 = 0
      IF (NBYTE.EQ.4) THEN
        DO I=1,8
           READ(61,REC=IREC2+I) RDES4(I)
        ENDDO
        IREC2=IREC2+8
        DO I=1,6
           READ(61,REC=IREC2+I) RID4(I)
        ENDDO
        IREC2=IREC2+6
        DO I=1,6
           READ(61,REC=IREC2+I) AID4(I)
        ENDDO
        IREC2=IREC2+6
      ENDIF
      IF (NBYTE.EQ.8) THEN
        DO I=1,4
           READ(61,REC=IREC2+I) RDES8(I)
        ENDDO
        IREC2=IREC2+4
        DO I=1,3
           READ(61,REC=IREC2+I) RID8(I)
        ENDDO
        IREC2=IREC2+3
        DO I=1,3
           READ(61,REC=IREC2+I) AID8(I)
        ENDDO
        IREC2=IREC2+3
      ENDIF
C
c     IF (NBYTE.EQ.4) THEN
c       WRITE(71,*) RDES4
c       WRITE(71,*) RID4
c       WRITE(71,*) AID4
c     ELSE
c       WRITE(71,*) RDES8
c       WRITE(71,*) RID8
c       WRITE(71,*) AID8
c     ENDIF
C
C--Read NTRSPE, NSTAE, DT*NSPOOLE from both files
C  and then close both files to flush file buffers
C
      READ(51,REC=IREC1+1) NTRSPE1
      READ(51,REC=IREC1+2) NSTAE1
      READ(51,REC=IREC1+3) DTE1
      READ(51,REC=IREC1+4) NSPOOLE1
      READ(51,REC=IREC1+5) ITEMPE1
      IREC1 = IREC1+5
C
      READ(61,REC=IREC2+1) NTRSPE2
      READ(61,REC=IREC2+2) NSTAE2
      READ(61,REC=IREC2+3) DTE2
      READ(61,REC=IREC2+4) NSPOOLE2
      READ(61,REC=IREC2+5) ITEMPE2
      IREC2 = IREC2+5
C
      CLOSE(51)         ! Flush File Buffer for file 1
      CLOSE(61)         ! Flush File buffer for file 2
      OPEN(51,FILE=FNAME1,ACCESS='DIRECT',RECL=NBYTE)
      OPEN(61,FILE=FNAME2,ACCESS='DIRECT',RECL=NBYTE)
C
      DO J=1,NTRSPE1
C
         READ(51,REC=IREC1+1) TIMEOUTSE1
         READ(51,REC=IREC1+2) ITSE1
         IREC1 = IREC1+2
C
         READ(61,REC=IREC2+1) TIMEOUTSE2
         READ(61,REC=IREC2+2) ITSE2
         IREC2 = IREC2+2
C
         DO I=1, NSTAE1
            READ(51,REC=IREC1+I) ETA1(I)
            READ(61,REC=IREC2+I) ETA2(I)
            ERR = ABS(ETA1(I)-ETA2(I))
            IF (ERR.GT.TOL) WRITE(71,*) TIMEOUTSE1,I,ETA1(I),ETA2(I)
         ENDDO
         IREC1 = IREC1 + NSTAE1
         IREC2 = IREC2 + NSTAE2
C
      ENDDO                 
C
C--Close both Global fort.61 Files
C
9999  CONTINUE
      CLOSE(51)
      CLOSE(61)
      CLOSE(71)
C
   80 FORMAT(A40)
 2120 FORMAT(2X,E20.10,5X,I10)
 3645 FORMAT(1X,I10,1X,I10,1X,E15.7,1X,I5,1X,I5)
C
      RETURN
      END

      SUBROUTINE COMPARE62(DIR1,DIR2,TOL)
C
C---------------------------------------------------------------------------C
C                     (  Serial Version  3/28/98  )                         C 
C  Compare data at the velocity stations from two different directories.    C
C  This version is compatible with ADCIRC version 33.04                     C
C                                                                           C
C  Input Parameters:                                                        C
C                                                                           C
C    DIR1  =  pathname of first directory                                   C
C    DIR2  =  pathname of second directory                                  C
C    TOL   =  user supplied relative error tolerance                        C
C                                                                           C
C---------------------------------------------------------------------------C
C
      USE POST_GLOBAL
      INTEGER I,J,K,L,IPROC,IDUM
      INTEGER NTRSPV1,NTRSPV2
      INTEGER NSTEMP1,NSTEMP2
      INTEGER ITSV1,ITSV2
      INTEGER ITEMPV1,ITEMPV2
      INTEGER NSTAV1,NSTAV2 
      INTEGER NUMSTNS1, NUMSTNS2
      INTEGER NSPOOLV1,NSPOOLV2 
      INTEGER IREC1,IREC2
      INTEGER LEN1,LEN2
      REAL(SZ) DTV1,DTV2
      REAL(8) TIMEOUTSV1,TIMEOUTSV2
      REAL(SZ) TOL,ERR
      CHARACTER*80  DIR1,DIR2,FNAME1,FNAME2,FNAME3
      CHARACTER*85  INLINE
      CHARACTER*4 RDES4(8),RID4(6),AID4(6)
      CHARACTER*8 RDES8(8),RID8(6),AID8(6)
      LOGICAL FOUND1,FOUND2

      REAL(8),ALLOCATABLE  ::  UU1(:),VV1(:)
      REAL(8),ALLOCATABLE  ::  UU2(:),VV2(:)
      ALLOCATE ( UU1(MNSTAV),VV1(MNSTAV))
      ALLOCATE ( UU2(MNSTAV),VV2(MNSTAV))
C
C--Construct File Names
C
      DO I=1, LEN(DIR1)
         IF (DIR1(I:I).EQ.' ') THEN
           LEN1 = I-1
           EXIT
         ENDIF
      ENDDO
      DO I=1, LEN(DIR2)
         IF (DIR2(I:I).EQ.' ') THEN
           LEN2 = I-1
           EXIT
         ENDIF
      ENDDO
C
      FNAME1 = DIR1(1:LEN1)//'/fort.62'
      FNAME2 = DIR2(1:LEN2)//'/fort.62'
      FNAME3 = 'diffs.62'
C
C--Determine whether Unit 61 is Sequential Formatted or Direct Access Binary
C
      IF (ABS(NOUTV).EQ.1) THEN
        GO TO 1000
      ELSE
        GO TO 2000
      ENDIF
C
1000  CONTINUE
C
C--Open Both Global Sequential Formatted fort.62 files
C
      INQUIRE(FILE=FNAME1,EXIST=FOUND1)
      INQUIRE(FILE=FNAME2,EXIST=FOUND2)
      IF (FOUND1.AND.FOUND2) THEN
        OPEN(UNIT=52,FILE=FNAME1)
        OPEN(UNIT=62,FILE=FNAME2)
        OPEN(UNIT=72,FILE=FNAME3)
      ELSE
        print *, "No fort.62 files found"
        RETURN
      ENDIF
C
      READ (52,'(A85)') INLINE
c     WRITE(72,'(A85)')  INLINE
      READ (62,'(A85)') INLINE
c     WRITE(72,'(A85)')  INLINE
C
      READ (52,3645) NTRSPV1,NSTAV1,DTV1,NSPOOLV1,ITEMPV1
      READ (62,3645) NTRSPV2,NSTAV2,DTV2,NSPOOLV2,ITEMPV2
C
      IF (NTRSPV1.NE.NTRSPV2) THEN
         WRITE(72,*) "NTRSPV1 not equal to NSTRSPV2"
         RETURN
      ENDIF
C
      IF (NSTAV1.NE.NSTAV2) THEN
        WRITE(72,*) "NSTAV1 not equal to NSTAV2"
        RETURN
      ENDIF
C
      IF (ABS(DTV1-DTV2).GT.1.0E-5) THEN
        WRITE(72,*) "DTV1 not equal to DTV2"
        RETURN
      ENDIF
C
      IF (NSPOOLV1.NE.NSPOOLV2) THEN
        WRITE(72,*)  "NSPOOLV1 not equal to NSPOOLV2"
        RETURN
      ENDIF
C
      IF (ITEMPV1.NE.ITEMPV2) THEN
        WRITE(72,*) "ITEMPV1 not equal to ITEMPV2"
        RETURN
      ENDIF
C
c     WRITE(72,*) NTRSPV1,NSTAV1,DTV1,NSPOOLV1,ITEMPV1
C
      DO J=1,NTRSPV1
C
         READ(52,2120) TIMEOUTSV1,ITSV1
         READ(62,2120) TIMEOUTSV2,ITSV2
C
         IF (ABS(TIMEOUTSV1-TIMEOUTSV2).GT.1.0E-5) THEN
           WRITE(72,*) "TIMEOUTSV1 not equal to TIMEOUTSV2"
           RETURN
         ENDIF
C
        IF (ITSV1.NE.ITSV2) THEN
          WRITE(72,*) "ITSV1 not equal to ITSV2"
          RETURN
        ENDIF
C
C       WRITE(72,*) TIMEOUTSV1,ITSV1
C
         DO I=1, NSTAV1
            READ(52,*) IDUM,UU1(I),VV1(I)
            READ(62,*) IDUM,UU2(I),VV2(I)
            ERR = SQRT(((UU1(I)-UU2(I))**2+(VV1(I)-VV2(I))**2))
            IF (ERR.GT.TOL) WRITE(72,*) TIMEOUTSV1,
     &        I,UU1(I),UU2(I),VV1(I),VV2(I)
         ENDDO
C
      ENDDO                 
      GO TO 9999
C 
2000  CONTINUE
C
C--Open Both Global Direct Access Binary fort.62 files
C
      INQUIRE(FILE=FNAME1,EXIST=FOUND1)
      INQUIRE(FILE=FNAME2,EXIST=FOUND2)
      IF (FOUND1.AND.FOUND2) THEN
        OPEN(52,FILE=FNAME1,ACCESS='DIRECT',RECL=NBYTE)
        OPEN(62,FILE=FNAME2,ACCESS='DIRECT',RECL=NBYTE)
        OPEN(UNIT=72,FILE=FNAME3)
      ELSE
        print *, "No fort.61 files found"
        RETURN
      ENDIF
C
C--Read RUNDES RUNID and AGRID from 1st Global File
C
      IREC1 = 0
      IF (NBYTE.EQ.4) THEN
        DO I=1,8
           READ(52,REC=IREC1+I) RDES4(I)
        ENDDO
        IREC1=IREC1+8
        DO I=1,6
           READ(52,REC=IREC1+I) RID4(I)
        ENDDO
        IREC1=IREC1+6
        DO I=1,6
           READ(52,REC=IREC1+I) AID4(I)
        ENDDO
        IREC1=IREC1+6
      ENDIF
      IF (NBYTE.EQ.8) THEN
        DO I=1,4
           READ(52,REC=IREC1+I) RDES8(I)
        ENDDO
        IREC1=IREC1+4
        DO I=1,3
           READ(52,REC=IREC1+I) RID8(I)
        ENDDO
        IREC1=IREC1+3
        DO I=1,3
           READ(52,REC=IREC1+I) AID8(I)
        ENDDO
        IREC1=IREC1+3
      ENDIF
C
C--Read RUNDES RUNID and AGRID from 2nd Global File
C
      IREC2 = 0
      IF (NBYTE.EQ.4) THEN
        DO I=1,8
           READ(62,REC=IREC2+I) RDES4(I)
        ENDDO
        IREC2=IREC2+8
        DO I=1,6
           READ(62,REC=IREC2+I) RID4(I)
        ENDDO
        IREC2=IREC2+6
        DO I=1,6
           READ(62,REC=IREC2+I) AID4(I)
        ENDDO
        IREC2=IREC2+6
      ENDIF
      IF (NBYTE.EQ.8) THEN
        DO I=1,4
           READ(62,REC=IREC2+I) RDES8(I)
        ENDDO
        IREC2=IREC2+4
        DO I=1,3
           READ(62,REC=IREC2+I) RID8(I)
        ENDDO
        IREC2=IREC2+3
        DO I=1,3
           READ(62,REC=IREC2+I) AID8(I)
        ENDDO
        IREC2=IREC2+3
      ENDIF
C
c     IF (NBYTE.EQ.4) THEN
c       WRITE(72,*) RDES4
c       WRITE(72,*) RID4
c       WRITE(72,*) AID4
c     ELSE
c       WRITE(72,*) RDES8
c       WRITE(72,*) RID8
c       WRITE(72,*) AID8
c     ENDIF
C
C--Read NTRSPE, NSTAE, DT*NSPOOLE from both files
C  and then close both files to flush file buffers
C
      READ(52,REC=IREC1+1) NTRSPV1
      READ(52,REC=IREC1+2) NSTAV1
      READ(52,REC=IREC1+3) DTV1
      READ(52,REC=IREC1+4) NSPOOLV1
      READ(52,REC=IREC1+5) ITEMPV1
      IREC1 = IREC1+5
C
      READ(62,REC=IREC2+1) NTRSPV2
      READ(62,REC=IREC2+2) NSTAV2
      READ(62,REC=IREC2+3) DTV2
      READ(62,REC=IREC2+4) NSPOOLV2
      READ(62,REC=IREC2+5) ITEMPV2
      IREC2 = IREC2+5
C
      CLOSE(52)         ! Flush File Buffer for file 1
      CLOSE(62)         ! Flush File buffer for file 2
      OPEN(52,FILE=FNAME1,ACCESS='DIRECT',RECL=NBYTE)
      OPEN(62,FILE=FNAME2,ACCESS='DIRECT',RECL=NBYTE)
C
      DO J=1,NTRSPV1
C
         READ(52,REC=IREC1+1) TIMEOUTSV1
         READ(52,REC=IREC1+2) ITSV1
         IREC1 = IREC1+2
C
         READ(62,REC=IREC2+1) TIMEOUTSV2
         READ(62,REC=IREC2+2) ITSV2
         IREC2 = IREC2+2
C
         DO I=1, NSTAV1
            READ(52,REC=IREC1+2*I-1) UU1(I)
            READ(52,REC=IREC1+2*I)   VV1(I)
            READ(62,REC=IREC2+2*I-1) UU2(I)
            READ(62,REC=IREC2+2*I)   VV2(I)
            ERR = SQRT(((UU1(I)-UU2(I))**2+(VV1(I)-VV2(I))**2))
            IF (ERR.GT.TOL) WRITE(72,*) TIMEOUTSV1,
     &        I,UU1(I),UU2(I),VV1(I),VV2(I)
         ENDDO
         IREC1 = IREC1 + 2*NSTAV1
         IREC2 = IREC2 + 2*NSTAV2
C
      ENDDO                 
C
C--Close both Global fort.61 Files
C
9999  CONTINUE
      CLOSE(52)
      CLOSE(62)
      CLOSE(72)
C
   80 FORMAT(A40)
 2120 FORMAT(2X,E20.10,5X,I10)
 3645 FORMAT(1X,I10,1X,I10,1X,E15.7,1X,I5,1X,I5)
C
      RETURN
      END

      SUBROUTINE COMPARE63(DIR1,DIR2,TOL)
C
C---------------------------------------------------------------------------C
C                     (  Serial Version  3/28/98  )                         C 
C  Compare elevation data at all nodes from two different directories.      C
C  This version is compatible with ADCIRC version 33.04                     C
C                                                                           C
C  Input Parameters:                                                        C
C                                                                           C
C    DIR1  =  pathname of first directory                                   C
C    DIR2  =  pathname of second directory                                  C
C    TOL   =  user supplied relative error tolerance                        C
C                                                                           C
C---------------------------------------------------------------------------C
C
      USE POST_GLOBAL
      INTEGER I,J,K,L,IPROC,IDUM
      INTEGER NDSETGE1,NDSETGE2
      INTEGER NP1,NP2
      INTEGER NSTEMP1,NSTEMP2
      INTEGER ITEMPE1,ITEMPE2
      INTEGER ITE1,ITE2
      INTEGER IREC1,IREC2
      INTEGER LEN1,LEN2
      REAL(SZ) DTE1,DTE2
      REAL(SZ) TOL,ERR
      REAL(8)  TIMEOUTE1,TIMEOUTE2
      CHARACTER*80  DIR1,DIR2,FNAME1,FNAME2,FNAME3
      CHARACTER*85  INLINE
      CHARACTER*4 RDES4(8),RID4(6),AID4(6)
      CHARACTER*8 RDES8(8),RID8(6),AID8(6)
      LOGICAL FOUND1,FOUND2

      REAL(8),ALLOCATABLE  ::   ETA1(:),ETA2(:)
      ALLOCATE ( ETA1(MNP),ETA2(MNP))
C
C--Construct File Names
C
      DO I=1, LEN(DIR1)
         IF (DIR1(I:I).EQ.' ') THEN
           LEN1 = I-1
           EXIT
         ENDIF
      ENDDO
      DO I=1, LEN(DIR2)
         IF (DIR2(I:I).EQ.' ') THEN
           LEN2 = I-1
           EXIT
         ENDIF
      ENDDO
C
      FNAME1 = DIR1(1:LEN1)//'/fort.63'
      FNAME2 = DIR2(1:LEN2)//'/fort.63'
      FNAME3 = 'diffs.63'
C
C--Determine whether Unit 63 is Sequential Formatted or Direct Access Binary
C
      IF (ABS(NOUTGE).EQ.1) THEN
        GO TO 1000
      ELSE
        GO TO 2000
      ENDIF
C
1000  CONTINUE
C
C--Open Both Global Sequential Formatted fort.63 files
C
      INQUIRE(FILE=FNAME1,EXIST=FOUND1)
      INQUIRE(FILE=FNAME2,EXIST=FOUND2)
      IF (FOUND1.AND.FOUND2) THEN
        OPEN(UNIT=53,FILE=FNAME1)
        OPEN(UNIT=63,FILE=FNAME2)
        OPEN(UNIT=73,FILE=FNAME3)
      ELSE
        print *, "No fort.63 files found"
        RETURN
      ENDIF
C
      READ (53,'(A85)') INLINE
c     WRITE(73,'(A85)')  INLINE
      READ (63,'(A85)') INLINE
c     WRITE(73,'(A85)')  INLINE
C
      READ (53,3645) NDSETGE1,NP1,DTE1,NSTEMP1,ITEMPE1
      READ (63,3645) NDSETGE2,NP2,DTE2,NSTEMP2,ITEMPE2
C
      IF (NDSETGE1.NE.NDSETGE2) THEN
         WRITE(73,*) "NDSETGE1 not equal to NDSETGE2"
         RETURN
      ENDIF
C
      IF (NP1.NE.NP2) THEN
        WRITE(73,*) "NP1 not equal to NP2"
        RETURN
      ENDIF
C
      IF (ABS(DTE1-DTE2).GT.1.0E-5) THEN
        WRITE(73,*) "DTE1 not equal to DTE2"
        RETURN
      ENDIF
C
      IF (NSTEMP1.NE.NSTEMP2) THEN
        WRITE(73,*) "NSTEMP1 not equal to NSTEMP2"
        RETURN
      ENDIF
      IF (ITEMPE1.NE.ITEMPE2) THEN
        WRITE(73,*) "ITEMPE1 not equal to ITEMPE2"
        RETURN
      ENDIF
C
c     WRITE(73,*) NDSETGE1,NP1,DTE1,NSTEMP1,ITEMPE1
C
      DO J=1,NDSETGE1
C
         READ(53,2120) TIMEOUTE1,ITE1
         READ(63,2120) TIMEOUTE2,ITE2
C
         IF (ABS(TIMEOUTE1-TIMEOUTE2).GT.1.0E-5) THEN
           WRITE(73,*) "TIMEOUTE1 not equal to TIMEOUTE2"
           RETURN
         ENDIF
C
        IF (ITE1.NE.ITE2) THEN
          WRITE(73,*) "ITE1 not equal to ITE2"
          RETURN
        ENDIF
C
C       WRITE(73,*) TIMEOUTE1,ITE1
C
         DO I=1, NNODG
            READ(53,*) IDUM,ETA1(I)
            READ(63,*) IDUM,ETA2(I)
            ERR = ABS(ETA1(I)-ETA2(I))
            IF (ERR.GT.TOL) WRITE(73,*) TIMEOUTE1,I,ETA1(I),ETA2(I)
         ENDDO
C
      ENDDO                 
      GO TO 9999
C 
2000  CONTINUE
C
C--Open Both Global Direct Access Binary fort.63 files
C
      INQUIRE(FILE=FNAME1,EXIST=FOUND1)
      INQUIRE(FILE=FNAME2,EXIST=FOUND2)
      IF (FOUND1.AND.FOUND2) THEN
        OPEN(53,FILE=FNAME1,ACCESS='DIRECT',RECL=NBYTE)
        OPEN(63,FILE=FNAME2,ACCESS='DIRECT',RECL=NBYTE)
        OPEN(UNIT=73,FILE=FNAME3)
      ELSE
        print *, "No fort.63 files found"
        RETURN
      ENDIF
C
C--Read RUNDES RUNID and AGRID from 1st Global File
C
      IREC1 = 0
      IF (NBYTE.EQ.4) THEN
        DO I=1,8
           READ(53,REC=IREC1+I) RDES4(I)
        ENDDO
        IREC1=IREC1+8
        DO I=1,6
           READ(53,REC=IREC1+I) RID4(I)
        ENDDO
        IREC1=IREC1+6
        DO I=1,6
           READ(53,REC=IREC1+I) AID4(I)
        ENDDO
        IREC1=IREC1+6
      ENDIF
      IF (NBYTE.EQ.8) THEN
        DO I=1,4
           READ(53,REC=IREC1+I) RDES8(I)
        ENDDO
        IREC1=IREC1+4
        DO I=1,3
           READ(53,REC=IREC1+I) RID8(I)
        ENDDO
        IREC1=IREC1+3
        DO I=1,3
           READ(53,REC=IREC1+I) AID8(I)
        ENDDO
        IREC1=IREC1+3
      ENDIF
C
C--Read RUNDES RUNID and AGRID from 2nd Global File
C
      IREC2 = 0
      IF (NBYTE.EQ.4) THEN
        DO I=1,8
           READ(63,REC=IREC2+I) RDES4(I)
        ENDDO
        IREC2=IREC2+8
        DO I=1,6
           READ(63,REC=IREC2+I) RID4(I)
        ENDDO
        IREC2=IREC2+6
        DO I=1,6
           READ(63,REC=IREC2+I) AID4(I)
        ENDDO
        IREC2=IREC2+6
      ENDIF
      IF (NBYTE.EQ.8) THEN
        DO I=1,4
           READ(63,REC=IREC2+I) RDES8(I)
        ENDDO
        IREC2=IREC2+4
        DO I=1,3
           READ(63,REC=IREC2+I) RID8(I)
        ENDDO
        IREC2=IREC2+3
        DO I=1,3
           READ(63,REC=IREC2+I) AID8(I)
        ENDDO
        IREC2=IREC2+3
C
      ENDIF
C
c     IF (NBYTE.EQ.4) THEN
c       WRITE(73,*) RDES4
c       WRITE(73,*) RID4
c       WRITE(73,*) AID4
c     ELSE
c       WRITE(73,*) RDES8
c       WRITE(73,*) RID8
c       WRITE(73,*) AID8
c     ENDIF
C
C--Read NTRSPE, NSTAE, DT*NSPOOLE from both files
C  and then close both files to flush file buffers
C
      READ(53,REC=IREC1+1) NDSETGE1
      READ(53,REC=IREC1+2) NP1
      READ(53,REC=IREC1+3) DTE1
      READ(53,REC=IREC1+4) NSTEMP1
      READ(53,REC=IREC1+5) ITEMPE1
      IREC1 = IREC1+5
C
      READ(63,REC=IREC2+1) NDSETGE2
      READ(63,REC=IREC2+2) NP2
      READ(63,REC=IREC2+3) DTE2
      READ(63,REC=IREC2+4) NSTEMP2
      READ(63,REC=IREC2+5) ITEMPE2
      IREC2 = IREC2+5
C
      CLOSE(53)         ! Flush File Buffer for file 1
      CLOSE(63)         ! Flush File buffer for file 2
      OPEN(53,FILE=FNAME1,ACCESS='DIRECT',RECL=NBYTE)
      OPEN(63,FILE=FNAME2,ACCESS='DIRECT',RECL=NBYTE)
C
      DO J=1,NDSETGE1
C
         READ(53,REC=IREC1+1) TIMEOUTE1
         READ(53,REC=IREC1+2) ITE1
         IREC1 = IREC1+2
C
         READ(63,REC=IREC2+1) TIMEOUTE2
         READ(63,REC=IREC2+2) ITE2
         IREC2 = IREC2+2
C
         DO I=1, NNODG
            READ(53,REC=IREC1+I) ETA1(I)
            READ(63,REC=IREC2+I) ETA2(I)
            ERR = ABS(ETA1(I)-ETA2(I))
            IF (ERR.GT.TOL) WRITE(73,*) TIMEOUTE1,I,ETA1(I),ETA2(I)
         ENDDO
         IREC1 = IREC1 + NNODG
         IREC2 = IREC2 + NNODG
C
      ENDDO                 
C
C--Close both Global fort.63 Files
C
9999  CONTINUE
      CLOSE(53)
      CLOSE(63)
      CLOSE(73)
C
   80 FORMAT(A40)
 2120 FORMAT(2X,E20.10,5X,I10)
 3645 FORMAT(1X,I10,1X,I10,1X,E15.7,1X,I5,1X,I5)
C
      RETURN
      END

      SUBROUTINE COMPARE64(DIR1,DIR2,TOL)
C
C---------------------------------------------------------------------------C
C                     (  Serial Version  3/28/98  )                         C 
C  Compare velocity data at all nodes from two different directories.       C
C  This version is compatible with ADCIRC version 33.04                     C
C                                                                           C
C  Input Parameters:                                                        C
C                                                                           C
C    DIR1  =  pathname of first directory                                   C
C    DIR2  =  pathname of second directory                                  C
C    TOL   =  user supplied relative error tolerance                        C
C                                                                           C
C---------------------------------------------------------------------------C
C
      USE POST_GLOBAL
      INTEGER I,J,K,L,IPROC,IDUM
      INTEGER NDSETGV1,NDSETGV2
      INTEGER NSTEMP1,NSTEMP2
      INTEGER ITV1,ITV2
      INTEGER ITEMPV1,ITEMPV2
      INTEGER NP1,NP2 
      INTEGER IREC1,IREC2
      INTEGER LEN1,LEN2
      REAL(SZ) DTV1,DTV2
      REAL(8) TIMEOUTV1,TIMEOUTV2
      REAL(SZ) TOL,ERR !jgf45.07 ERR is a fortran keyword...
      CHARACTER*80  DIR1,DIR2,FNAME1,FNAME2,FNAME3
      CHARACTER*85  INLINE
      CHARACTER*4 RDES4(8),RID4(6),AID4(6)
      CHARACTER*8 RDES8(8),RID8(6),AID8(6)
      LOGICAL FOUND1,FOUND2

      REAL(8),ALLOCATABLE :: UU1(:),VV1(:)
      REAL(8),ALLOCATABLE :: UU2(:),VV2(:)

      ALLOCATE ( UU1(MNP),VV1(MNP))
      ALLOCATE ( UU2(MNP),VV2(MNP))
C
C--Construct File Names
C
      DO I=1, LEN(DIR1)
         IF (DIR1(I:I).EQ.' ') THEN
           LEN1 = I-1
           EXIT
         ENDIF
      ENDDO
      DO I=1, LEN(DIR2)
         IF (DIR2(I:I).EQ.' ') THEN
           LEN2 = I-1
           EXIT
         ENDIF
      ENDDO
C
      FNAME1 = DIR1(1:LEN1)//'/fort.64'
      FNAME2 = DIR2(1:LEN2)//'/fort.64'
      FNAME3 = 'diffs.64'
C
C--Determine whether Unit 64 is Sequential Formatted or Direct Access Binary
C
      IF (ABS(NOUTGV).EQ.1) THEN
        GO TO 1000
      ELSE
        GO TO 2000
      ENDIF
C
1000  CONTINUE
C
C--Open Both Global Sequential Formatted fort.64 files
C
      INQUIRE(FILE=FNAME1,EXIST=FOUND1)
      INQUIRE(FILE=FNAME2,EXIST=FOUND2)
      IF (FOUND1.AND.FOUND2) THEN
         OPEN(UNIT=54,FILE=FNAME1)
        OPEN(UNIT=64,FILE=FNAME2)
        OPEN(UNIT=74,FILE=FNAME3)
      ELSE
        print *, "No fort.64 files found"
        RETURN
      ENDIF
C
      READ (54,'(A85)') INLINE
c     WRITE(74,'(A85)')  INLINE
      READ (64,'(A85)') INLINE
c     WRITE(74,'(A85)')  INLINE
C
      READ (54,3645) NDSETGV1,NP1,DTV1,NSTEMP1,ITEMPV1
      READ (64,3645) NDSETGV2,NP2,DTV2,NSTEMP2,ITEMPV2
C
      IF (NDSETGV1.NE.NDSETGV2) THEN
         WRITE(74,*) "NDSETGV1 not equal to NDSETGV2"
         RETURN
      ENDIF
C
      IF (NP1.NE.NP2) THEN
        WRITE(74,*) "NP1 not equal to NP2"
        RETURN
      ENDIF
C
      IF (ABS(DTV1-DTV2).GT.1.0E-5) THEN
        WRITE(74,*)  "DTV1 not equal to DTV2"
        RETURN
      ENDIF
C
      IF (NSTEMP1.NE.NSTEMP2) THEN
        WRITE(74,*) "NSTEMP1 not equal to NSTEMP2"
        RETURN
      ENDIF
      IF (ITEMPV1.NE.ITEMPV2) THEN
        WRITE(74,*) "ITEMPV1 not equal to ITEMPV2"
        RETURN
      ENDIF
C
c     WRITE(74,*) NDSETGV1,NNODG,DTV1,NSTEMP1,ITEMPV1
C
      DO J=1,NDSETGV1
C
         READ(54,2120) TIMEOUTV1,ITV1
         READ(64,2120) TIMEOUTV2,ITV2
C
         IF (ABS(TIMEOUTV1-TIMEOUTV2).GT.1.0E-5) THEN
           WRITE(74,*) "TIMEOUTV1 not equal to TIMEOUTV2"
           RETURN
         ENDIF
C
        IF (ITV1.NE.ITV2) THEN
          WRITE(74,*) "ITV1 not equal to ITV2"
          RETURN
        ENDIF
C
C       WRITE(74,*) TIMEOUTV1,ITV1
C
         DO I=1, NNODG
            READ(54,*) IDUM,UU1(I),VV1(I)
            READ(64,*) IDUM,UU2(I),VV2(I)
            ERR = SQRT(((UU1(I)-UU2(I))**2+(VV1(I)-VV2(I))**2))
            IF (ERR.GT.TOL) WRITE(74,*) TIMEOUTV1,
     &        I,UU1(I),UU2(I),VV1(I),VV2(I)
         ENDDO
C
      ENDDO                 
      GO TO 9999
C 
2000  CONTINUE
C
C--Open Both Global Direct Access Binary fort.64 files
C
      INQUIRE(FILE=FNAME1,EXIST=FOUND1)
      INQUIRE(FILE=FNAME2,EXIST=FOUND2)
      IF (FOUND1.AND.FOUND2) THEN
        OPEN(54,FILE=FNAME1,ACCESS='DIRECT',RECL=NBYTE)
        OPEN(64,FILE=FNAME2,ACCESS='DIRECT',RECL=NBYTE)
        OPEN(UNIT=74,FILE=FNAME3)
      ELSE
        print *, "No fort.64 files found"
        RETURN
      ENDIF
C
C--Read RUNDES RUNID and AGRID from 1st Global File
C
      IREC1 = 0
      IF (NBYTE.EQ.4) THEN
        DO I=1,8
           READ(54,REC=IREC1+I) RDES4(I)
        ENDDO
        IREC1=IREC1+8
        DO I=1,6
           READ(54,REC=IREC1+I) RID4(I)
        ENDDO
        IREC1=IREC1+6
        DO I=1,6
           READ(54,REC=IREC1+I) AID4(I)
        ENDDO
        IREC1=IREC1+6
      ENDIF
      IF (NBYTE.EQ.8) THEN
        DO I=1,4
           READ(54,REC=IREC1+I) RDES8(I)
        ENDDO
        IREC1=IREC1+4
        DO I=1,3
           READ(54,REC=IREC1+I) RID8(I)
        ENDDO
        IREC1=IREC1+3
        DO I=1,3
           READ(54,REC=IREC1+I) AID8(I)
        ENDDO
        IREC1=IREC1+3
      ENDIF
C
C--Read RUNDES RUNID and AGRID from 2nd Global File
C
      IREC2 = 0
      IF (NBYTE.EQ.4) THEN
        DO I=1,8
           READ(64,REC=IREC2+I) RDES4(I)
        ENDDO
        IREC2=IREC2+8
        DO I=1,6
           READ(64,REC=IREC2+I) RID4(I)
        ENDDO
        IREC2=IREC2+6
        DO I=1,6
           READ(64,REC=IREC2+I) AID4(I)
        ENDDO
        IREC2=IREC2+6
      ENDIF
      IF (NBYTE.EQ.8) THEN
        DO I=1,4
           READ(64,REC=IREC2+I) RDES8(I)
        ENDDO
        IREC2=IREC2+4
        DO I=1,3
           READ(64,REC=IREC2+I) RID8(I)
        ENDDO
        IREC2=IREC2+3
        DO I=1,3
           READ(64,REC=IREC2+I) AID8(I)
        ENDDO
        IREC2=IREC2+3
      ENDIF
C
c     IF (NBYTE.EQ.4) THEN
c       WRITE(74,*) RDES4
c       WRITE(74,*) RID4
c       WRITE(74,*) AID4
c     ELSE
c       WRITE(74,*) RDES8
c       WRITE(74,*) RID8
c       WRITE(74,*) AID8
c     ENDIF
C
C--Read NTRSPE, NSTAE, DT*NSPOOLE from both files
C  and then close both files to flush file buffers
C
      READ(54,REC=IREC1+1) NDSETGV1
      READ(54,REC=IREC1+2) NP1
      READ(54,REC=IREC1+3) DTV1
      READ(54,REC=IREC1+4) NSTEMP1
      READ(54,REC=IREC1+5) ITEMPV1
      IREC1 = IREC1+5
C
      READ(64,REC=IREC2+1) NDSETGV2
      READ(64,REC=IREC2+2) NP2
      READ(64,REC=IREC2+3) DTV2
      READ(64,REC=IREC2+4) NSTEMP2
      READ(64,REC=IREC2+5) ITEMPV2
      IREC2 = IREC2+5
C
      CLOSE(54)         ! Flush File Buffer for file 1
      CLOSE(64)         ! Flush File buffer for file 2
      OPEN(54,FILE=FNAME1,ACCESS='DIRECT',RECL=NBYTE)
      OPEN(64,FILE=FNAME2,ACCESS='DIRECT',RECL=NBYTE)
C
      DO J=1,NDSETGV1
C
         READ(54,REC=IREC1+1) TIMEOUTV1
         READ(54,REC=IREC1+2) ITV1
         IREC1 = IREC1+2
C
         READ(64,REC=IREC2+1) TIMEOUTV2
         READ(64,REC=IREC2+2) ITV2
         IREC2 = IREC2+2
C
         DO I=1, NNODG
            READ(54,REC=IREC1+2*I-1) UU1(I)
            READ(54,REC=IREC1+2*I)   VV1(I)
            READ(64,REC=IREC2+2*I-1) UU2(I)
            READ(64,REC=IREC2+2*I)   VV2(I)
            ERR = SQRT(((UU1(I)-UU2(I))**2+(VV1(I)-VV2(I))**2))
            IF (ERR.GT.TOL) WRITE(74,*) TIMEOUTV1,
     &        I,UU1(I),UU2(I),VV1(I),VV2(I)
         ENDDO
         IREC1 = IREC1 + 2*NNODG
         IREC2 = IREC2 + 2*NNODG
C
      ENDDO                 
C
C--Close both Global fort.64 Files and Differences Log file
C
9999  CONTINUE
      CLOSE(54)
      CLOSE(64)
      CLOSE(74)
C
   80 FORMAT(A40)
 2120 FORMAT(2X,E20.10,5X,I10)
 3645 FORMAT(1X,I10,1X,I10,1X,E15.7,1X,I5,1X,I5)
C
      RETURN
      END

      SUBROUTINE COMPARE51(DIR1,DIR2,TOL)
C
C---------------------------------------------------------------------------C
C                     (  Serial Version  3/28/98  )                         C 
C  Compare the harmonic data at the elevation stations from two different   C
C  directories.                                                             C
C  This version is compatible with ADCIRC version 33.04                     C
C                                                                           C
C  Input Parameters:                                                        C
C                                                                           C
C    DIR1  =  pathname of first directory                                   C
C    DIR2  =  pathname of second directory                                  C
C    TOL   =  user supplied relative error tolerance                        C
C                                                                           C
C---------------------------------------------------------------------------C
C
      USE POST_GLOBAL
      INTEGER I,J,IPROC,IDUM
      INTEGER NFREQ1,NFREQ2,NUMSTNS
      INTEGER LEN1,LEN2
      REAL(SZ) TOL,ERR !jgf45.07 ERR is a fortran keyword...
      REAL(SZ) MAG1,MAG2,PHAS1,PHAS2
      CHARACTER*80  DIR1,DIR2,FNAME1,FNAME2,FNAME3
      CHARACTER*80  INLINE,HEADER1,HEADER2
      LOGICAL FOUND1,FOUND2

      CHARACTER*80,ALLOCATABLE  ::  HARDAT1(:),HARDAT2(:)
      ALLOCATE ( HARDAT1(MNSTAE),HARDAT2(MNSTAE))
C
C--Construct File Names
C
      DO I=1, LEN(DIR1)
         IF (DIR1(I:I).EQ.' ') THEN
           LEN1 = I-1
           EXIT
         ENDIF
      ENDDO
      DO I=1, LEN(DIR2)
         IF (DIR2(I:I).EQ.' ') THEN
           LEN2 = I-1
           EXIT
         ENDIF
      ENDDO
C
      FNAME1 = DIR1(1:LEN1)//'/fort.51'
      FNAME2 = DIR2(1:LEN2)//'/fort.51'
      FNAME3 = 'diffs.51'
C
C--Open Global Sequential Formatted fort.51 files and Differences Log file
C
      INQUIRE(FILE=FNAME1,EXIST=FOUND1)
      INQUIRE(FILE=FNAME2,EXIST=FOUND2)
      IF (FOUND1.AND.FOUND2) THEN
        OPEN(UNIT=51,FILE=FNAME1)
        OPEN(UNIT=61,FILE=FNAME2)
        OPEN(UNIT=71,FILE=FNAME3)
      ELSE
        print *, "No fort.51 files found"
        RETURN
      ENDIF
C
      READ(51,'(A80)')  INLINE
      READ(INLINE,*) NFREQ1
c     WRITE(71,'(A80)') INLINE
      READ(61,'(A80)')  INLINE
      READ(INLINE,*) NFREQ2
c     WRITE(71,'(A80)') INLINE
      IF (NFREQ1.NE.NFREQ2) THEN
        WRITE(71,*) "NFREQ different in the two fort.51 files"
        RETURN
      ELSE
c     WRITE(71,'(A80)') INLINE
      ENDIF
C
      DO I=1, NFREQ1
         READ(51,'(A80)') HEADER1
         READ(61,'(A80)') HEADER2
      ENDDO
      READ(51,*) NUMSTNS
      READ(61,*) NUMSTNS
C
      DO J=1,NSTAE
         READ(51,*) IDUM
         READ(61,*) IDUM
         DO I=1, NFREQ1
            READ(51,'(A80)') HARDAT1(J)
            READ(HARDAT1(J),*) MAG1,PHAS1
            READ(61,'(A80)') HARDAT2(J)
            READ(HARDAT2(J),*) MAG2,PHAS2
            IF ((ABS(MAG1-MAG2).GT.TOL)
     &        .OR.(ABS(PHAS1-PHAS2).GT.TOL)) THEN
              WRITE(71,*)  "difference at elevation station ",J
              WRITE(71,'(A80)') HARDAT1(J)
              WRITE(71,'(A80)') HARDAT2(J)
            ENDIF
         ENDDO
      ENDDO
C
C--Close the Global fort.51 Files and the Differences Log file
C
      CLOSE(51)
      CLOSE(61)
      CLOSE(71)
C
c6635 FORMAT(2X,E16.8,1X,F11.4)

C
      RETURN
      END

      SUBROUTINE COMPARE52(DIR1,DIR2,TOL)
C
C---------------------------------------------------------------------------C
C                     (  Serial Version  3/28/98  )                         C 
C  Compare the harmonic data at the velocity stations from two different    C
C  directories.                                                             C
C  This version is compatible with ADCIRC version 33.04                     C
C                                                                           C
C  Input Parameters:                                                        C
C                                                                           C
C    DIR1  =  pathname of first directory                                   C
C    DIR2  =  pathname of second directory                                  C
C    TOL   =  user supplied relative error tolerance                        C
C                                                                           C
C---------------------------------------------------------------------------C
C
      USE POST_GLOBAL
      INTEGER I,J,IPROC,IDUM
      INTEGER NFREQ1,NFREQ2,NUMSTNS
      INTEGER LEN1,LEN2
      REAL(SZ) TOL,ERR !jgf45.07 ERR is a fortran keyword...
      REAL(SZ) UMAG1,UMAG2,VMAG1,VMAG2
      REAL(SZ) UPH1,UPH2,VPH1,VPH2
      CHARACTER*80  DIR1,DIR2,FNAME1,FNAME2,FNAME3
      CHARACTER*80  INLINE,HEADER1,HEADER2
      LOGICAL FOUND1,FOUND2

      CHARACTER*80,ALLOCATABLE  ::  HARDAT1(:),HARDAT2(:)
      ALLOCATE ( HARDAT1(MNSTAE),HARDAT2(MNSTAE))
C
C--Construct File Names
C
      DO I=1, LEN(DIR1)
         IF (DIR1(I:I).EQ.' ') THEN
           LEN1 = I-1
           EXIT
         ENDIF
      ENDDO
      DO I=1, LEN(DIR2)
         IF (DIR2(I:I).EQ.' ') THEN
           LEN2 = I-1
           EXIT
         ENDIF
      ENDDO
C
      FNAME1 = DIR1(1:LEN1)//'/fort.52'
      FNAME2 = DIR2(1:LEN2)//'/fort.52'
      FNAME3 = 'diffs.52'
C
C--Open Global Sequential Formatted fort.52 files and Differences Log file
C
      INQUIRE(FILE=FNAME1,EXIST=FOUND1)
      INQUIRE(FILE=FNAME2,EXIST=FOUND2)
      IF (FOUND1.AND.FOUND2) THEN
        OPEN(UNIT=52,FILE=FNAME1)    
        OPEN(UNIT=62,FILE=FNAME2)    
        OPEN(UNIT=72,FILE=FNAME3)    
      ELSE
        print *, "No fort.52 files found"
        RETURN
      ENDIF
C
      READ(52,'(A80)')  INLINE
      READ(INLINE,*) NFREQ1
c     WRITE(72,'(A80)') INLINE
      READ(62,'(A80)')  INLINE
      READ(INLINE,*) NFREQ2
c     WRITE(72,'(A80)') INLINE
      IF (NFREQ1.NE.NFREQ2) THEN
        WRITE(72,*) "NFREQ different in the two fort.52 files"
        RETURN
      ELSE
c     WRITE(72,'(A80)') INLINE
      ENDIF
C
      DO I=1, NFREQ1
         READ(52,'(A80)') HEADER1
         READ(62,'(A80)') HEADER2
      ENDDO
C
      READ(52,*) NUMSTNS
      READ(62,*) NUMSTNS
C
      DO J=1,NSTAV
         READ(52,*) IDUM
         READ(62,*) IDUM
         DO I=1, NFREQ1
            READ(52,'(A80)') HARDAT1(J)
            READ(HARDAT1(J),*) UMAG1,UPH1,VMAG1,VPH1
            READ(62,'(A80)') HARDAT2(J)
            READ(HARDAT2(J),*) UMAG2,UPH2,VMAG2,VPH2
            IF ((ABS(UMAG1-UMAG2).GT.TOL)
     &       .OR.(ABS(UPH1-UPH2).GT.TOL)
     &       .OR.(ABS(VMAG1-VMAG2).GT.TOL)
     &       .OR.(ABS(VPH1-VPH2).GT.TOL)) THEN
             WRITE(72,*) "difference at velocity station ",J
             WRITE(72,'(A80)') HARDAT1(J)
             WRITE(72,'(A80)') HARDAT2(J)
            ENDIF
         ENDDO
      ENDDO
C
C--Close the Global fort.52 Files and the Differences Log file
C
      CLOSE(52)
      CLOSE(62)
      CLOSE(72)
C
c3679 FORMAT(1X,E20.10,1X,F10.7,1X,F12.8,1X,A10)
c6636 format(2x,e16.8,1x,f11.4,2x,e16.8,1x,f11.4)
C
      RETURN
      END

      SUBROUTINE COMPARE53(DIR1,DIR2,TOL)
C
C---------------------------------------------------------------------------C
C                     (  Serial Version  3/28/98  )                         C 
C  Compare the harmonic constituent elevations at all nodes from two        C
C  different directories.                                                   C
C  This version is compatible with ADCIRC version 33.04                     C
C                                                                           C
C  Input Parameters:                                                        C
C                                                                           C
C    DIR1  =  pathname of first directory                                   C
C    DIR2  =  pathname of second directory                                  C
C    TOL   =  user supplied relative error tolerance                        C
C                                                                           C
C---------------------------------------------------------------------------C
C
      USE POST_GLOBAL
      INTEGER I,J,IPROC,IDUM
      INTEGER NFREQ1,NFREQ2,NP1,NP2
      INTEGER LEN1,LEN2
      REAL(SZ) TOL,ERR !jgf45.07 ERR is a fortran keyword...
      REAL(SZ) MAG1,MAG2,PHAS1,PHAS2
      CHARACTER*80  DIR1,DIR2,FNAME1,FNAME2,FNAME3
      CHARACTER*80  INLINE,HEADER1,HEADER2
      LOGICAL FOUND1,FOUND2

      CHARACTER*80,ALLOCATABLE  ::   HARDAT1(:),HARDAT2(:)
      ALLOCATE (HARDAT1(MNP),HARDAT2(MNP))
C
C--Construct File Names
C
      DO I=1, LEN(DIR1)
         IF (DIR1(I:I).EQ.' ') THEN
           LEN1 = I-1
           EXIT
         ENDIF
      ENDDO
      DO I=1, LEN(DIR2)
         IF (DIR2(I:I).EQ.' ') THEN
           LEN2 = I-1
           EXIT
         ENDIF
      ENDDO
C
      FNAME1 = DIR1(1:LEN1)//'/fort.53'
      FNAME2 = DIR2(1:LEN2)//'/fort.53'
      FNAME3 = 'diffs.53'
C
C--Open Global Sequential Formatted fort.53 files and Differences Log file
C
      INQUIRE(FILE=FNAME1,EXIST=FOUND1)
      INQUIRE(FILE=FNAME2,EXIST=FOUND2)
      IF (FOUND1.AND.FOUND2) THEN
        OPEN(UNIT=53,FILE=FNAME1)    
        OPEN(UNIT=63,FILE=FNAME2)    
        OPEN(UNIT=73,FILE=FNAME3)    
      ELSE
        print *, "No fort.53 files found"
        RETURN
      ENDIF
C
      READ(53,'(A80)')  INLINE
      READ(INLINE,*) NFREQ1
c     WRITE(73,'(A80)') INLINE
      READ(63,'(A80)')  INLINE
      READ(INLINE,*) NFREQ2
c     WRITE(73,'(A80)') INLINE
      IF (NFREQ1.NE.NFREQ2) THEN
        WRITE(73,*) "NFREQ different in the two fort.53 files"
        RETURN
      ELSE
c     WRITE(73,'(A80)') INLINE
      ENDIF
C
      DO I=1, NFREQ1
         READ(53,'(A80)') HEADER1
         READ(63,'(A80)') HEADER2
      ENDDO
      READ(53,*) NP1     
      READ(63,*) NP2 
C
      DO J=1,NNODG
         READ(53,*) IDUM
         READ(63,*) IDUM
         DO I=1, NFREQ1
            READ(53,'(A80)') HARDAT1(J)
            READ(HARDAT1(J),*) MAG1,PHAS1
            READ(63,'(A80)') HARDAT2(J)
            READ(HARDAT2(J),*) MAG2,PHAS2
            IF ((ABS(MAG1-MAG2).GT.TOL)
     &        .OR.(ABS(PHAS1-PHAS2).GT.TOL)) THEN
              WRITE(73,*) "elevation difference at node ",J
              WRITE(73,'(A80)') HARDAT1(J)
              WRITE(73,'(A80)') HARDAT2(J)
            ENDIF
         ENDDO
      ENDDO
C
C--Close the Global fort.53 Files and the Differences Log file
C
      CLOSE(53)
      CLOSE(63)
      CLOSE(73)
C
c3679 FORMAT(1X,E20.10,1X,F10.7,1X,F12.8,1X,A10)
C
      RETURN
      END


      SUBROUTINE COMPARE54(DIR1,DIR2,TOL)
C
C---------------------------------------------------------------------------C
C                     (  Serial Version  3/28/98  )                         C 
C  Compare the harmonic constituent velocities at all nodes from two        C
C  different directories.                                                   C
C  This version is compatible with ADCIRC version 33.04                     C
C                                                                           C
C  Input Parameters:                                                        C
C                                                                           C
C    DIR1  =  pathname of first directory                                   C
C    DIR2  =  pathname of second directory                                  C
C    TOL   =  user supplied relative error tolerance                        C
C                                                                           C
C---------------------------------------------------------------------------C
C
      USE POST_GLOBAL
      INTEGER I,J,IPROC,IDUM
      INTEGER NFREQ1,NFREQ2,NP1,NP2
      INTEGER LEN1,LEN2
      REAL(SZ) TOL,ERR !jgf45.07 ERR is a fortran keyword...
      REAL(SZ) UMAG1,UMAG2,VMAG1,VMAG2
      REAL(SZ) UPH1,UPH2,VPH1,VPH2
      CHARACTER*80  DIR1,DIR2,FNAME1,FNAME2,FNAME3
      CHARACTER*80  INLINE,HEADER1,HEADER2
      LOGICAL FOUND1,FOUND2

      CHARACTER*80,ALLOCATABLE  ::  HARDAT1(:),HARDAT2(:)
      ALLOCATE ( HARDAT1(MNP),HARDAT2(MNP))
C
C--Construct File Names
C
      DO I=1, LEN(DIR1)
         IF (DIR1(I:I).EQ.' ') THEN
           LEN1 = I-1
           EXIT
         ENDIF
      ENDDO
      DO I=1, LEN(DIR2)
         IF (DIR2(I:I).EQ.' ') THEN
           LEN2 = I-1
           EXIT
         ENDIF
      ENDDO
C
      FNAME1 = DIR1(1:LEN1)//'/fort.54'
      FNAME2 = DIR2(1:LEN2)//'/fort.54'
      FNAME3 = 'diffs.54'
C
C--Open Global Sequential Formatted fort.54 files and Differences Log file
C
      INQUIRE(FILE=FNAME1,EXIST=FOUND1)
      INQUIRE(FILE=FNAME2,EXIST=FOUND2)
      IF (FOUND1.AND.FOUND2) THEN
        OPEN(UNIT=54,FILE=FNAME1)    
        OPEN(UNIT=64,FILE=FNAME2)    
        OPEN(UNIT=74,FILE=FNAME3)    
      ELSE
        print *, "No fort.54 files found"
        RETURN
      ENDIF
C
      READ(54,'(A80)')  INLINE
      READ(INLINE,*) NFREQ1
c     WRITE(74,'(A80)') INLINE
      READ(64,'(A80)')  INLINE
      READ(INLINE,*) NFREQ2
c     WRITE(74,'(A80)') INLINE
      IF (NFREQ1.NE.NFREQ2) THEN
        WRITE(74,*) "NFREQ different in the two fort.54 files"
        RETURN
      ELSE
c     WRITE(74,'(A80)') INLINE
      ENDIF
C
      DO I=1, NFREQ1
         READ(54,'(A80)') HEADER1
         READ(64,'(A80)') HEADER2
      ENDDO
C
      READ(54,*) NP1     
      READ(64,*) NP2 
C
      DO J=1,NNODG
         READ(54,*) IDUM
         READ(64,*) IDUM
         DO I=1, NFREQ1
            READ(54,'(A80)') HARDAT1(J)
            READ(HARDAT1(J),*) UMAG1,UPH1,VMAG1,VPH1
            READ(64,'(A80)') HARDAT2(J)
            READ(HARDAT2(J),*) UMAG2,UPH2,VMAG2,VPH2
            IF ((ABS(UMAG1-UMAG2).GT.TOL)
     &        .OR.(ABS(UPH1-UPH2).GT.TOL)
     &        .OR.(ABS(VMAG1-VMAG2).GT.TOL)
     &        .OR.(ABS(VPH1-VPH2).GT.TOL)) THEN
              WRITE(74,*) "velocity difference at node ",J
              WRITE(74,'(A80)') HARDAT1(J)
              WRITE(74,'(A80)') HARDAT2(J)
            ENDIF
         ENDDO
      ENDDO
C
C--Close the Global fort.54 Files and the Differences Log file
C
      CLOSE(54)
      CLOSE(64)
      CLOSE(74)
C
      RETURN
      END

      SUBROUTINE COMPARE55(DIR1,DIR2,TOL)
C
C---------------------------------------------------------------------------C
C                     (  Serial Version  4/13/98  )                         C 
C  Compare the harmonic constituent comparison files from two directories.  C
C  This version is compatible with ADCIRC version 33.04                     C
C                                                                           C
C  Input Parameters:                                                        C
C                                                                           C
C    DIR1  =  pathname of first directory                                   C
C    DIR2  =  pathname of second directory                                  C
C    TOL   =  user supplied relative error tolerance                        C
C                                                                           C
C---------------------------------------------------------------------------C
C
      USE POST_GLOBAL
      INTEGER I,J,IPROC,IDUM
      INTEGER LEN1,LEN2
      INTEGER NP1,NP2
      REAL(SZ) TOL,ERR !jgf45.07 ERR is a fortran keyword
      REAL(SZ) ELAV1,ELAV2
      REAL(SZ) EAV1,EAV2
      REAL(SZ) EAVDIF1,EAVDIF2
      REAL(SZ) ELVA1,ELVA2
      REAL(SZ) ESQ1,ESQ2
      REAL(SZ) EVADIF1,EVADIF2
      REAL(SZ) XVELAV1,XVELAV2
      REAL(SZ) UAV1,UAV2
      REAL(SZ) UAVDIF1,UAVDIF2
      REAL(SZ) XVELVA1,XVELVA2
      REAL(SZ) USQ1,USQ2
      REAL(SZ) UVADIF1,UVADIF2
      REAL(SZ) YVELAV1,YVELAV2
      REAL(SZ) VAV1,VAV2
      REAL(SZ) VAVDIF1,VAVDIF2
      REAL(SZ) YVELVA1,YVELVA2
      REAL(SZ) VSQ1,VSQ2
      REAL(SZ) VVADIF1,VVADIF2
      CHARACTER*80  DIR1,DIR2,FNAME1,FNAME2,FNAME3
      CHARACTER*114 ELEV1,ELEV2
      CHARACTER*114 UU1,VV1,UU2,VV2
      LOGICAL FOUND1,FOUND2
C
C--Construct File Names
C
      DO I=1, LEN(DIR1)
         IF (DIR1(I:I).EQ.' ') THEN
           LEN1 = I-1
           EXIT
         ENDIF
      ENDDO
      DO I=1, LEN(DIR2)
         IF (DIR2(I:I).EQ.' ') THEN
           LEN2 = I-1
           EXIT
         ENDIF
      ENDDO
C
      FNAME1 = DIR1(1:LEN1)//'/fort.55'
      FNAME2 = DIR2(1:LEN2)//'/fort.55'
      FNAME3 = 'diffs.55'
C
C--Open Global Sequential Formatted fort.55 files and Differences Log file
C
      INQUIRE(FILE=FNAME1,EXIST=FOUND1)
      INQUIRE(FILE=FNAME2,EXIST=FOUND2)
      IF (FOUND1.AND.FOUND2) THEN
        OPEN(UNIT=55,FILE=FNAME1)    
        OPEN(UNIT=65,FILE=FNAME2)    
        OPEN(UNIT=75,FILE=FNAME3)    
      ELSE
        print *, "No fort.55 files found"
        RETURN
      ENDIF
C
      READ(55,*) NP1     
      READ(65,*) NP2 
C
      DO J=1,NNODG
         READ(55,*) IDUM 
         READ(65,*) IDUM 
C
         READ(55,'(A114)') ELEV1
         READ(ELEV1,*) ELAV1,EAV1,EAVDIF1,ELVA1,ESQ1,EVADIF1
C
         READ(65,'(A114)') ELEV2
         READ(ELEV2,*) ELAV2,EAV2,EAVDIF2,ELVA2,ESQ2,EVADIF2
C
         IF ((ABS(ELAV1-ELAV2).GT.TOL)
     &   .OR.(ABS(EAV1-EAV2).GT.TOL)      
     &   .OR.(ABS(EAVDIF1-EAVDIF2).GT.TOL)
     &   .OR.(ABS(ELVA1-ELVA2).GT.TOL)
     &   .OR.(ABS(ESQ1-ESQ2).GT.TOL)
     &   .OR.(ABS(EVADIF1-EVADIF2).GT.TOL)) THEN
           WRITE(75,*) "elevation difference at node ",J
           WRITE(75,'(A114)') ELEV1       
           WRITE(75,'(A114)') ELEV2       
         ENDIF
      ENDDO
C
      DO J=1,NNODG
         READ(55,*) IDUM
         READ(55,'(A114)') UU1
         READ(55,'(A114)') VV1
         READ(UU1,*) XVELAV1,UAV1,UAVDIF1,XVELVA1,USQ1,UVADIF1
         READ(VV1,*) YVELAV1,VAV1,VAVDIF1,YVELVA1,VSQ1,VVADIF1
C
         READ(65,*) IDUM 
         READ(65,'(A114)') UU2
         READ(65,'(A114)') VV2
         READ(UU2,*) XVELAV2,UAV2,UAVDIF2,XVELVA2,USQ2,UVADIF2
         READ(VV2,*) YVELAV2,VAV2,VAVDIF2,YVELVA2,VSQ2,VVADIF2
C
         IF ((ABS(XVELAV1-XVELAV2).GT.TOL)
     &   .OR.(ABS(UAV1-UAV2).GT.TOL)      
     &   .OR.(ABS(UAVDIF1-UAVDIF2).GT.TOL)
     &   .OR.(ABS(XVELVA1-XVELVA2).GT.TOL)
     &   .OR.(ABS(USQ1-USQ2).GT.TOL)
     &   .OR.(ABS(UVADIF1-UVADIF2).GT.TOL)
C
     &   .OR.(ABS(YVELAV1-YVELAV2).GT.TOL)
     &   .OR.(ABS(VAV1-VAV2).GT.TOL)      
     &   .OR.(ABS(VAVDIF1-VAVDIF2).GT.TOL)
     &   .OR.(ABS(YVELVA1-YVELVA2).GT.TOL)
     &   .OR.(ABS(VSQ1-VSQ2).GT.TOL)
     &   .OR.(ABS(VVADIF1-VVADIF2).GT.TOL)) THEN
C
           WRITE(75,*) "velocity difference at node ",J
           WRITE(75,'(A114)') UU1         
           WRITE(75,'(A114)') VV1         
           WRITE(75,'(A114)') UU2         
           WRITE(75,'(A114)') VV2         
         ENDIF
      ENDDO
C
C
C--Close the Global fort.55 Files and the Differences Log file
C
      CLOSE(55)
      CLOSE(65)
      CLOSE(75)
C
      RETURN
      END

      SUBROUTINE COMPARE74(DIR1,DIR2,TOL)
C
C---------------------------------------------------------------------------C
C                     (  Serial Version  3/28/98  )                         C 
C  Compare Wind Stress data at all nodes from two different directories.    C
C  This version is compatible with ADCIRC version 33.04                     C
C                                                                           C
C  Input Parameters:                                                        C
C                                                                           C
C    DIR1  =  pathname of first directory                                   C
C    DIR2  =  pathname of second directory                                  C
C    TOL   =  user supplied relative error tolerance                        C
C                                                                           C
C---------------------------------------------------------------------------C
C
      USE POST_GLOBAL
      INTEGER I,J,K,L,IPROC,IDUM
      INTEGER NDSETGV1,NDSETGV2
      INTEGER NSTEMP1,NSTEMP2
      INTEGER ITV1,ITV2
      INTEGER ITEMPV1,ITEMPV2
      INTEGER NP1,NP2 
      INTEGER IREC1,IREC2
      INTEGER LEN1,LEN2
      REAL(SZ) DTV1,DTV2
      REAL(8) TIMEOUTV1,TIMEOUTV2
      REAL(SZ) TOL,ERR !jgf45.07 ERR is a fortran keyword...
      CHARACTER*80  DIR1,DIR2,FNAME1,FNAME2,FNAME3
      CHARACTER*85  INLINE
      CHARACTER*4 RDES4(8),RID4(6),AID4(6)
      CHARACTER*8 RDES8(8),RID8(6),AID8(6)
      LOGICAL FOUND1,FOUND2

      REAL(8),ALLOCATABLE :: UU1(:),VV1(:)
      REAL(8),ALLOCATABLE :: UU2(:),VV2(:)

      ALLOCATE ( UU1(MNP),VV1(MNP))
      ALLOCATE ( UU2(MNP),VV2(MNP))
C
C--Construct File Names
C
      DO I=1, LEN(DIR1)
         IF (DIR1(I:I).EQ.' ') THEN
           LEN1 = I-1
           EXIT
         ENDIF
      ENDDO
      DO I=1, LEN(DIR2)
         IF (DIR2(I:I).EQ.' ') THEN
           LEN2 = I-1
           EXIT
         ENDIF
      ENDDO
C
      FNAME1 = DIR1(1:LEN1)//'/fort.74'
      FNAME2 = DIR2(1:LEN2)//'/fort.74'
      FNAME3 = 'diffs.74'
C
C--Determine whether Unit 74 is Sequential Formatted or Direct Access Binary
C
      IF (ABS(NOUTGW).EQ.1) THEN
        GO TO 1000
      ELSE
        GO TO 2000
      ENDIF
C
1000  CONTINUE
C
C--Open Both Global Sequential Formatted fort.74 files
C
      INQUIRE(FILE=FNAME1,EXIST=FOUND1)
      INQUIRE(FILE=FNAME2,EXIST=FOUND2)
      IF (FOUND1.AND.FOUND2) THEN
         OPEN(UNIT=54,FILE=FNAME1)
        OPEN(UNIT=64,FILE=FNAME2)
        OPEN(UNIT=74,FILE=FNAME3)
      ELSE
        print *, "No fort.64 files found"
        RETURN
      ENDIF
C
      READ (54,'(A85)') INLINE
c     WRITE(74,'(A85)')  INLINE
      READ (64,'(A85)') INLINE
c     WRITE(74,'(A85)')  INLINE
C
      READ (54,3645) NDSETGV1,NP1,DTV1,NSTEMP1,ITEMPV1
      READ (64,3645) NDSETGV2,NP2,DTV2,NSTEMP2,ITEMPV2
C
      IF (NDSETGV1.NE.NDSETGV2) THEN
         WRITE(74,*) "NDSETGV1 not equal to NDSETGV2"
         RETURN
      ENDIF
C
      IF (NP1.NE.NP2) THEN
        WRITE(74,*) "NP1 not equal to NP2"
        RETURN
      ENDIF
C
      IF (ABS(DTV1-DTV2).GT.1.0E-5) THEN
        WRITE(74,*)  "DTV1 not equal to DTV2"
        RETURN
      ENDIF
C
      IF (NSTEMP1.NE.NSTEMP2) THEN
        WRITE(74,*) "NSTEMP1 not equal to NSTEMP2"
        RETURN
      ENDIF
      IF (ITEMPV1.NE.ITEMPV2) THEN
        WRITE(74,*) "ITEMPV1 not equal to ITEMPV2"
        RETURN
      ENDIF
C
c     WRITE(74,*) NDSETGV1,NNODG,DTV1,NSTEMP1,ITEMPV1
C
      DO J=1,NDSETGV1
C
         READ(54,2120) TIMEOUTV1,ITV1
         READ(64,2120) TIMEOUTV2,ITV2
C
         IF (ABS(TIMEOUTV1-TIMEOUTV2).GT.1.0E-5) THEN
           WRITE(74,*) "TIMEOUTV1 not equal to TIMEOUTV2"
           RETURN
         ENDIF
C
        IF (ITV1.NE.ITV2) THEN
          WRITE(74,*) "ITV1 not equal to ITV2"
          RETURN
        ENDIF
C
C       WRITE(74,*) TIMEOUTV1,ITV1
C
         DO I=1, NNODG
            READ(54,*) IDUM,UU1(I),VV1(I)
            READ(64,*) IDUM,UU2(I),VV2(I)
            ERR = SQRT(((UU1(I)-UU2(I))**2+(VV1(I)-VV2(I))**2))
            IF (ERR.GT.TOL) WRITE(74,*) I,UU1(I),UU2(I),VV1(I),VV2(I)
         ENDDO
C
      ENDDO                 
      GO TO 9999
C 
2000  CONTINUE
C
C--Open Both Global Direct Access Binary fort.64 files
C
      INQUIRE(FILE=FNAME1,EXIST=FOUND1)
      INQUIRE(FILE=FNAME2,EXIST=FOUND2)
      IF (FOUND1.AND.FOUND2) THEN
        OPEN(54,FILE=FNAME1,ACCESS='DIRECT',RECL=NBYTE)
        OPEN(64,FILE=FNAME2,ACCESS='DIRECT',RECL=NBYTE)
        OPEN(UNIT=74,FILE=FNAME3)
      ELSE
        print *, "No fort.64 files found"
        RETURN
      ENDIF
C
C--Read RUNDES RUNID and AGRID from 1st Global File
C
      IREC1 = 0
      IF (NBYTE.EQ.4) THEN
        DO I=1,8
           READ(54,REC=IREC1+I) RDES4(I)
        ENDDO
        IREC1=IREC1+8
        DO I=1,6
           READ(54,REC=IREC1+I) RID4(I)
        ENDDO
        IREC1=IREC1+6
        DO I=1,6
           READ(54,REC=IREC1+I) AID4(I)
        ENDDO
        IREC1=IREC1+6
      ENDIF
      IF (NBYTE.EQ.8) THEN
        DO I=1,4
           READ(54,REC=IREC1+I) RDES8(I)
        ENDDO
        IREC1=IREC1+4
        DO I=1,3
           READ(54,REC=IREC1+I) RID8(I)
        ENDDO
        IREC1=IREC1+3
        DO I=1,3
           READ(54,REC=IREC1+I) AID8(I)
        ENDDO
        IREC1=IREC1+3
      ENDIF
C
C--Read RUNDES RUNID and AGRID from 2nd Global File
C
      IREC2 = 0
      IF (NBYTE.EQ.4) THEN
        DO I=1,8
           READ(64,REC=IREC2+I) RDES4(I)
        ENDDO
        IREC2=IREC2+8
        DO I=1,6
           READ(64,REC=IREC2+I) RID4(I)
        ENDDO
        IREC2=IREC2+6
        DO I=1,6
           READ(64,REC=IREC2+I) AID4(I)
        ENDDO
        IREC2=IREC2+6
      ENDIF
      IF (NBYTE.EQ.8) THEN
        DO I=1,4
           READ(64,REC=IREC2+I) RDES8(I)
        ENDDO
        IREC2=IREC2+4
        DO I=1,3
           READ(64,REC=IREC2+I) RID8(I)
        ENDDO
        IREC2=IREC2+3
        DO I=1,3
           READ(64,REC=IREC2+I) AID8(I)
        ENDDO
        IREC2=IREC2+3
      ENDIF
C
c     IF (NBYTE.EQ.4) THEN
c       WRITE(74,*) RDES4
c       WRITE(74,*) RID4
c       WRITE(74,*) AID4
c     ELSE
c       WRITE(74,*) RDES8
c       WRITE(74,*) RID8
c       WRITE(74,*) AID8
c     ENDIF
C
C--Read NTRSPE, NSTAE, DT*NSPOOLE from both files
C  and then close both files to flush file buffers
C
      READ(54,REC=IREC1+1) NDSETGV1
      READ(54,REC=IREC1+2) NP1
      READ(54,REC=IREC1+3) DTV1
      READ(54,REC=IREC1+4) NSTEMP1
      READ(54,REC=IREC1+5) ITEMPV1
      IREC1 = IREC1+5
C
      READ(64,REC=IREC2+1) NDSETGV2
      READ(64,REC=IREC2+2) NP2
      READ(64,REC=IREC2+3) DTV2
      READ(64,REC=IREC2+4) NSTEMP2
      READ(64,REC=IREC2+5) ITEMPV2
      IREC2 = IREC2+5
C
      CLOSE(54)         ! Flush File Buffer for file 1
      CLOSE(64)         ! Flush File buffer for file 2
      OPEN(54,FILE=FNAME1,ACCESS='DIRECT',RECL=NBYTE)
      OPEN(64,FILE=FNAME2,ACCESS='DIRECT',RECL=NBYTE)
C
      DO J=1,NDSETGV1
C
         READ(54,REC=IREC1+1) TIMEOUTV1
         READ(54,REC=IREC1+2) ITV1
         IREC1 = IREC1+2
C
         READ(64,REC=IREC2+1) TIMEOUTV2
         READ(64,REC=IREC2+2) ITV2
         IREC2 = IREC2+2
C
         DO I=1, NNODG
            READ(54,REC=IREC1+2*I-1) UU1(I)
            READ(54,REC=IREC1+2*I)   VV1(I)
            READ(64,REC=IREC2+2*I-1) UU2(I)
            READ(64,REC=IREC2+2*I)   VV2(I)
            ERR = SQRT(((UU1(I)-UU2(I))**2+(VV1(I)-VV2(I))**2))
            IF (ERR.GT.TOL) WRITE(74,*) I,UU1(I),UU2(I),VV1(I),VV2(I)
         ENDDO
         IREC1 = IREC1 + 2*NNODG
         IREC2 = IREC2 + 2*NNODG
C
      ENDDO                 
C
C--Close both Global fort.64 Files and Differences Log file
C
9999  CONTINUE
      CLOSE(54)
      CLOSE(64)
      CLOSE(74)
C
   80 FORMAT(A40)
 2120 FORMAT(2X,E20.10,5X,I10)
 3645 FORMAT(1X,I10,1X,I10,1X,E15.7,1X,I5,1X,I5)
C
      RETURN
      END
