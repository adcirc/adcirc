
!   Version 1.1  vjp  5/04/99


    SUBROUTINE DIFFMERGE63(DIR1,DIR2,SCALEFAC)

!---------------------------------------------------------------------------C
!                     (  Serial Version  3/15/99  )                         C
!  Create a Difference File for elevation data at all nodes from two        C
!  different directories.                                                   C
!  This version is compatible with ADCIRC version 33.04                     C
!                                                                           C
!  Input Parameters:                                                        C
!                                                                           C
!    DIR1  =  pathname of first directory                                   C
!    DIR2  =  pathname of second directory                                  C
!    SCALEFAC =  user supplied scale factor for difference                  C
!                                                                           C
!---------------------------------------------------------------------------C

    USE POST_GLOBAL
    INTEGER :: I,J,K,L,IPROC,IDUM
    INTEGER :: NDSETGE1,NDSETGE2
    INTEGER :: NP1,NP2
    INTEGER :: NSTEMP1,NSTEMP2
    INTEGER :: ITEMPE1,ITEMPE2
    INTEGER :: ITE1,ITE2
    INTEGER :: IREC1,IREC2
    INTEGER :: LEN1,LEN2
    REAL(SZ) DTE1,DTE2,SCALEFAC
    REAL(8) TIMEOUTE1,TIMEOUTE2
    CHARACTER(80) ::  DIR1,DIR2,FNAME1,FNAME2,FNAME3
    CHARACTER(85) ::  INLINE
    CHARACTER(4) :: RDES4(8),RID4(6),AID4(6)
    CHARACTER(8) :: RDES8(8),RID8(6),AID8(6)
    LOGICAL :: FOUND1,FOUND2

    REAL(8),ALLOCATABLE  ::  ETA1(:),ETA2(:)
    ALLOCATE( ETA1(MNP),ETA2(MNP))

!--Construct File Names

    DO I=1, LEN(DIR1)
        IF (DIR1(I:I) == ' ') THEN
            LEN1 = I-1
            GO TO 100
        ENDIF
    ENDDO
    100 CONTINUE
    DO I=1, LEN(DIR2)
        IF (DIR2(I:I) == ' ') THEN
            LEN2 = I-1
            GO TO 200
        ENDIF
    ENDDO
    200 CONTINUE

    FNAME1 = DIR1(1:LEN1)//'/fort.63'
    FNAME2 = DIR2(1:LEN2)//'/fort.63'
    FNAME3 = 'diffmerge.63'

!--Determine whether Unit 63 is Sequential Formatted or Direct Access Binary

    IF (ABS(NOUTGE) == 1) THEN
        GO TO 1000
    ELSE
        GO TO 2000
    ENDIF

    1000 CONTINUE

!--Open Both Global Sequential Formatted fort.63 files

    INQUIRE(FILE=FNAME1,EXIST=FOUND1)
    INQUIRE(FILE=FNAME2,EXIST=FOUND2)
    IF (FOUND1 .AND. FOUND2) THEN
        OPEN(UNIT=53,FILE=FNAME1)
        OPEN(UNIT=63,FILE=FNAME2)
        OPEN(UNIT=73,FILE=FNAME3)
    ELSE
        print *, "No fort.63 files found"
        RETURN
    ENDIF

    READ (53,'(A85)') INLINE
    READ (63,'(A85)') INLINE
    WRITE(73,'(A85)')  INLINE

    READ (53,3645) NDSETGE1,NP1,DTE1,NSTEMP1,ITEMPE1
    READ (63,3645) NDSETGE2,NP2,DTE2,NSTEMP2,ITEMPE2

    IF (NDSETGE1 /= NDSETGE2) THEN
        print *, "NDSETGE1 not equal to NDSETGE2"
        RETURN
    ENDIF

    IF (NP1 /= NP2) THEN
        print *, "NP1 not equal to NP2"
        RETURN
    ENDIF

    IF (ABS(DTE1-DTE2) > 1.0E-5) THEN
        print *, "DTE1 not equal to DTE2"
        RETURN
    ENDIF

    IF (NSTEMP1 /= NSTEMP2) THEN
        print *, "NSTEMP1 not equal to NSTEMP2"
        RETURN
    ENDIF
    IF (ITEMPE1 /= ITEMPE2) THEN
        print *, "ITEMPE1 not equal to ITEMPE2"
        RETURN
    ENDIF

    WRITE(73,3645) NDSETGE1,NP1,DTE1,NSTEMP1,ITEMPE1

    DO J=1,NDSETGE1
    
        READ(53,2120) TIMEOUTE1,ITE1
        READ(63,2120) TIMEOUTE2,ITE2
    
        IF (ABS(TIMEOUTE1-TIMEOUTE2) > 1.0E-5) THEN
            print *, "TIMEOUTE1 not equal to TIMEOUTE2"
            RETURN
        ENDIF
    
        IF (ITE1 /= ITE2) THEN
            print *, "ITE1 not equal to ITE2"
            RETURN
        ENDIF
    
        WRITE(73,2120) TIMEOUTE1,ITE1
    
        DO I=1, NNODG
            READ(53,*) IDUM,ETA1(I)
            READ(63,*) IDUM,ETA2(I)
            WRITE(73,2453) I,SCALEFAC*(ETA1(I)-ETA2(I))
        ENDDO
    
    ENDDO
    GO TO 9999

    2000 CONTINUE

!--Open Both Global Direct Access Binary fort.63 files

    INQUIRE(FILE=FNAME1,EXIST=FOUND1)
    INQUIRE(FILE=FNAME2,EXIST=FOUND2)
    IF (FOUND1 .AND. FOUND2) THEN
        OPEN(53,FILE=FNAME1,ACCESS='DIRECT',RECL=NBYTE)
        OPEN(63,FILE=FNAME2,ACCESS='DIRECT',RECL=NBYTE)
        OPEN(73,FILE=FNAME3,ACCESS='DIRECT',RECL=NBYTE)
    ELSE
        print *, "No fort.63 files found"
        RETURN
    ENDIF

!--Read RUNDES RUNID and AGRID from 1st Global File

    IREC1 = 0
    IF (NBYTE == 4) THEN
        DO I=1,8
            READ(53,REC=IREC1+I) RDES4(I)
            WRITE(73,REC=IREC1+I) RDES4(I)
        ENDDO
        IREC1=IREC1+8
        DO I=1,6
            READ(53,REC=IREC1+I) RID4(I)
            WRITE(73,REC=IREC1+I) RID4(I)
        ENDDO
        IREC1=IREC1+6
        DO I=1,6
            READ(53,REC=IREC1+I) AID4(I)
            WRITE(73,REC=IREC1+I) AID4(I)
        ENDDO
        IREC1=IREC1+6
    ENDIF
    IF (NBYTE == 8) THEN
        DO I=1,4
            READ(53,REC=IREC1+I) RDES8(I)
            WRITE(73,REC=IREC1+I) RDES8(I)
        ENDDO
        IREC1=IREC1+4
        DO I=1,3
            READ(53,REC=IREC1+I) RID8(I)
            WRITE(73,REC=IREC1+I) RID8(I)
        ENDDO
        IREC1=IREC1+3
        DO I=1,3
            READ(53,REC=IREC1+I) AID8(I)
            WRITE(73,REC=IREC1+I) AID8(I)
        ENDDO
        IREC1=IREC1+3
    ENDIF

!--Read RUNDES RUNID and AGRID from 2nd Global File

    IREC2 = 0
    IF (NBYTE == 4) THEN
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
    IF (NBYTE == 8) THEN
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
    
    ENDIF


!--Read NTRSPE, NSTAE, DT*NSPOOLE from both files
!  and then close both files to flush file buffers

    READ(53,REC=IREC1+1) NDSETGE1
    WRITE(73,REC=IREC1+1) NDSETGE1
    READ(53,REC=IREC1+2) NP1
    WRITE(73,REC=IREC1+2) NP1
    READ(53,REC=IREC1+3) DTE1
    WRITE(73,REC=IREC1+3) DTE1
    READ(53,REC=IREC1+4) NSTEMP1
    WRITE(73,REC=IREC1+4) NSTEMP1
    READ(53,REC=IREC1+5) ITEMPE1
    WRITE(73,REC=IREC1+5) ITEMPE1
    IREC1 = IREC1+5

    READ(63,REC=IREC2+1) NDSETGE2
    READ(63,REC=IREC2+2) NP2
    READ(63,REC=IREC2+3) DTE2
    READ(63,REC=IREC2+4) NSTEMP2
    READ(63,REC=IREC2+5) ITEMPE2
    IREC2 = IREC2+5

    CLOSE(53)         ! Flush File Buffer for file 1
    CLOSE(63)         ! Flush File buffer for file 2
    CLOSE(73)         ! Flush File buffer for file 3
    OPEN(53,FILE=FNAME1,ACCESS='DIRECT',RECL=NBYTE)
    OPEN(63,FILE=FNAME2,ACCESS='DIRECT',RECL=NBYTE)
    OPEN(73,FILE=FNAME2,ACCESS='DIRECT',RECL=NBYTE)

    DO J=1,NDSETGE1
    
        READ(53,REC=IREC1+1) TIMEOUTE1
        WRITE(73,REC=IREC1+1) TIMEOUTE1
        READ(53,REC=IREC1+2) ITE1
        WRITE(73,REC=IREC1+2) ITE1
        IREC1 = IREC1+2
    
        READ(63,REC=IREC2+1) TIMEOUTE2
        READ(63,REC=IREC2+2) ITE2
        IREC2 = IREC2+2
    
        DO I=1, NNODG
            READ(53,REC=IREC1+I) ETA1(I)
            READ(63,REC=IREC2+I) ETA2(I)
            WRITE(73,REC=IREC1+I)  I,SCALEFAC*(ETA1(I)-ETA2(I))
        ENDDO
        IREC1 = IREC1 + NNODG
        IREC2 = IREC2 + NNODG
    
    ENDDO

!--Close both Global fort.63 Files

    9999 CONTINUE
    CLOSE(53)
    CLOSE(63)
    CLOSE(73)

    80 FORMAT(A40)
    2120 FORMAT(2X,E20.10,5X,I10)
    2453 FORMAT(2X,I8,2X,E15.8)
    3645 FORMAT(1X,I10,1X,I10,1X,E15.7,1X,I5,1X,I5)

    RETURN
    END SUBROUTINE DIFFMERGE63

    SUBROUTINE DIFFMERGE64(DIR1,DIR2,SCALEFAC)

!---------------------------------------------------------------------------C
!                     (  Serial Version  3/15/99  )                         C
!  Create a Difference file from velocity data at all nodes from two        C
!  different directories.                                                   C
!  This version is compatible with ADCIRC version 33.04                     C
!                                                                           C
!  Input Parameters:                                                        C
!                                                                           C
!    DIR1  =  pathname of first directory                                   C
!    DIR2  =  pathname of second directory                                  C
!    SCALEFAC =  user supplied scale factor for difference                  C
!                                                                           C
!---------------------------------------------------------------------------C

    USE POST_GLOBAL
    INTEGER :: I,J,K,L,IPROC,IDUM
    INTEGER :: NDSETGV1,NDSETGV2
    INTEGER :: NSTEMP1,NSTEMP2
    INTEGER :: ITV1,ITV2
    INTEGER :: ITEMPV1,ITEMPV2
    INTEGER :: NP1,NP2
    INTEGER :: IREC1,IREC2
    INTEGER :: LEN1,LEN2
    REAL(SZ) DTV1,DTV2,SCALEFAC
    REAL(8) TIMEOUTV1,TIMEOUTV2
    CHARACTER(80) ::  DIR1,DIR2,FNAME1,FNAME2,FNAME3
    CHARACTER(85) ::  INLINE
    CHARACTER(4) :: RDES4(8),RID4(6),AID4(6)
    CHARACTER(8) :: RDES8(8),RID8(6),AID8(6)
    LOGICAL :: FOUND1,FOUND2

    REAL(8),ALLOCATABLE  ::  UU1(:),VV1(:)
    REAL(8),ALLOCATABLE  ::  UU2(:),VV2(:)
    ALLOCATE ( UU1(MNP),VV1(MNP))
    ALLOCATE ( UU2(MNP),VV2(MNP))

!--Construct File Names

    DO I=1, LEN(DIR1)
        IF (DIR1(I:I) == ' ') THEN
            LEN1 = I-1
            GO TO 100
        ENDIF
    ENDDO
    100 CONTINUE
    DO I=1, LEN(DIR2)
        IF (DIR2(I:I) == ' ') THEN
            LEN2 = I-1
            GO TO 200
        ENDIF
    ENDDO
    200 CONTINUE

    FNAME1 = DIR1(1:LEN1)//'/fort.64'
    FNAME2 = DIR2(1:LEN2)//'/fort.64'
    FNAME3 = 'diffmerge.64'

!--Determine whether Unit 64 is Sequential Formatted or Direct Access Binary

    IF (ABS(NOUTGV) == 1) THEN
        GO TO 1000
    ELSE
        GO TO 2000
    ENDIF

    1000 CONTINUE

!--Open Both Global Sequential Formatted fort.64 files

    INQUIRE(FILE=FNAME1,EXIST=FOUND1)
    INQUIRE(FILE=FNAME2,EXIST=FOUND2)
    IF (FOUND1 .AND. FOUND2) THEN
        OPEN(UNIT=54,FILE=FNAME1)
        OPEN(UNIT=64,FILE=FNAME2)
        OPEN(UNIT=74,FILE=FNAME3)
    ELSE
        print *, "No fort.64 files found"
        RETURN
    ENDIF

    READ (54,'(A85)') INLINE
    READ (64,'(A85)') INLINE
    WRITE(74,'(A85)')  INLINE

    READ (54,3645) NDSETGV1,NP1,DTV1,NSTEMP1,ITEMPV1
    READ (64,3645) NDSETGV2,NP2,DTV2,NSTEMP2,ITEMPV2

    IF (NDSETGV1 /= NDSETGV2) THEN
        print *, "NDSETGV1 not equal to NDSETGV2"
        RETURN
    ENDIF

    IF (NP1 /= NP2) THEN
        print *, "NP1 not equal to NP2"
        RETURN
    ENDIF

    IF (ABS(DTV1-DTV2) > 1.0E-5) THEN
        print *,  "DTV1 not equal to DTV2"
        RETURN
    ENDIF

    IF (NSTEMP1 /= NSTEMP2) THEN
        print *, "NSTEMP1 not equal to NSTEMP2"
        RETURN
    ENDIF
    IF (ITEMPV1 /= ITEMPV2) THEN
        print *, "ITEMPV1 not equal to ITEMPV2"
        RETURN
    ENDIF

    WRITE(74,3645) NDSETGV1,NP1,DTV1,NSTEMP1,ITEMPV1

    DO J=1,NDSETGV1
    
        READ(54,2120) TIMEOUTV1,ITV1
        READ(64,2120) TIMEOUTV2,ITV2
    
        IF (ABS(TIMEOUTV1-TIMEOUTV2) > 1.0E-5) THEN
            print *, "TIMEOUTV1 not equal to TIMEOUTV2"
            RETURN
        ENDIF
    
        IF (ITV1 /= ITV2) THEN
            print *, "ITV1 not equal to ITV2"
            RETURN
        ENDIF
    
        WRITE(74,2120) TIMEOUTV1,ITV1
    
        DO I=1, NNODG
            READ(54,*) IDUM,UU1(I),VV1(I)
            READ(64,*) IDUM,UU2(I),VV2(I)
            WRITE(74,2454) I,SCALEFAC*(UU1(I)-UU2(I)), &
            SCALEFAC*(VV1(I)-VV2(I))
        ENDDO
    
    ENDDO
    GO TO 9999

    2000 CONTINUE

!--Open Both Global Direct Access Binary fort.64 files

    INQUIRE(FILE=FNAME1,EXIST=FOUND1)
    INQUIRE(FILE=FNAME2,EXIST=FOUND2)
    IF (FOUND1 .AND. FOUND2) THEN
        OPEN(54,FILE=FNAME1,ACCESS='DIRECT',RECL=NBYTE)
        OPEN(64,FILE=FNAME2,ACCESS='DIRECT',RECL=NBYTE)
        OPEN(74,FILE=FNAME3,ACCESS='DIRECT',RECL=NBYTE)
    ELSE
        print *, "No fort.64 files found"
        RETURN
    ENDIF

!--Read RUNDES RUNID and AGRID from 1st Global File

    IREC1 = 0
    IF (NBYTE == 4) THEN
        DO I=1,8
            READ(54,REC=IREC1+I) RDES4(I)
            WRITE(74,REC=IREC1+I) RDES4(I)
        ENDDO
        IREC1=IREC1+8
        DO I=1,6
            READ(54,REC=IREC1+I) RID4(I)
            WRITE(74,REC=IREC1+I) RID4(I)
        ENDDO
        IREC1=IREC1+6
        DO I=1,6
            READ(54,REC=IREC1+I) AID4(I)
            WRITE(74,REC=IREC1+I) AID4(I)
        ENDDO
        IREC1=IREC1+6
    ENDIF
    IF (NBYTE == 8) THEN
        DO I=1,4
            READ(54,REC=IREC1+I) RDES8(I)
            WRITE(74,REC=IREC1+I) RDES8(I)
        ENDDO
        IREC1=IREC1+4
        DO I=1,3
            READ(54,REC=IREC1+I) RID8(I)
            WRITE(74,REC=IREC1+I) RID8(I)
        ENDDO
        IREC1=IREC1+3
        DO I=1,3
            READ(54,REC=IREC1+I) AID8(I)
            WRITE(74,REC=IREC1+I) AID8(I)
        ENDDO
        IREC1=IREC1+3
    ENDIF

!--Read RUNDES RUNID and AGRID from 2nd Global File

    IREC2 = 0
    IF (NBYTE == 4) THEN
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
    IF (NBYTE == 8) THEN
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

!--Read NTRSPE, NSTAE, DT*NSPOOLE from both files
!  and then close both files to flush file buffers

    READ(54,REC=IREC1+1) NDSETGV1
    WRITE(74,REC=IREC1+1) NDSETGV1
    READ(54,REC=IREC1+2) NP1
    WRITE(74,REC=IREC1+2) NP1
    READ(54,REC=IREC1+3) DTV1
    WRITE(74,REC=IREC1+3) DTV1
    READ(54,REC=IREC1+4) NSTEMP1
    WRITE(74,REC=IREC1+4) NSTEMP1
    READ(54,REC=IREC1+5) ITEMPV1
    WRITE(74,REC=IREC1+5) ITEMPV1
    IREC1 = IREC1+5

    READ(64,REC=IREC2+1) NDSETGV2
    READ(64,REC=IREC2+2) NP2
    READ(64,REC=IREC2+3) DTV2
    READ(64,REC=IREC2+4) NSTEMP2
    READ(64,REC=IREC2+5) ITEMPV2
    IREC2 = IREC2+5

    CLOSE(54)         ! Flush File Buffer for file 1
    CLOSE(64)         ! Flush File buffer for file 2
    CLOSE(74)         ! Flush File buffer for file 3
    OPEN(54,FILE=FNAME1,ACCESS='DIRECT',RECL=NBYTE)
    OPEN(64,FILE=FNAME2,ACCESS='DIRECT',RECL=NBYTE)
    OPEN(74,FILE=FNAME3,ACCESS='DIRECT',RECL=NBYTE)

    DO J=1,NDSETGV1
    
        READ(54,REC=IREC1+1) TIMEOUTV1
        WRITE(74,REC=IREC1+1) TIMEOUTV1
        READ(54,REC=IREC1+2) ITV1
        WRITE(74,REC=IREC1+2) ITV1
        IREC1 = IREC1+2
    
        READ(64,REC=IREC2+1) TIMEOUTV2
        READ(64,REC=IREC2+2) ITV2
        IREC2 = IREC2+2
    
        DO I=1, NNODG
            READ(54,REC=IREC1+2*I-1) UU1(I)
            READ(54,REC=IREC1+2*I)   VV1(I)
            READ(64,REC=IREC2+2*I-1) UU2(I)
            READ(64,REC=IREC2+2*I)   VV2(I)
            WRITE(74,REC=IREC1+2*I-1) SCALEFAC*(UU1(I)-UU2(I))
            WRITE(74,REC=IREC1+2*I)   SCALEFAC*(VV1(I)-VV2(I))
        ENDDO
        IREC1 = IREC1 + 2*NNODG
        IREC2 = IREC2 + 2*NNODG
    
    ENDDO

!--Close both Global fort.64 Files and Differences Log file

    9999 CONTINUE
    CLOSE(54)
    CLOSE(64)
    CLOSE(74)

    80 FORMAT(A40)
    2120 FORMAT(2X,E20.10,5X,I10)
    2454 FORMAT(2X,I8,2(2X,E15.8))
    3645 FORMAT(1X,I10,1X,I10,1X,E15.7,1X,I5,1X,I5)

    RETURN
    END SUBROUTINE DIFFMERGE64
