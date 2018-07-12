C
C  version 1.1 vjp 5/04/99
C  updated by vjp 5/11/99 to handle incomplete files
C  version 1.2 vjp 12/12/99
C  updated by vjp  12/12/99 to allocate variables dynamically

C     jgf Updated for 45.06 09/07/2005 to reflect changes in hot start
C     file format and to add code for post processing fort.73 file.

C     jgf Updated for 45.07 11/07/2005 to handle fort.71 and fort.72;
C     made TIME double precision in order to fix a bug that was causing
C     TIME in merged output files to be different than subdomain
C     files. Added code to change the format of fort.80 file to handle
C     NOFF array in hotstart files.

C     jgf Updated for 45.11 Jan 2006 to handle new fort.41--46, new 3D
C     recording stations defined by coordinates rather than node
C     numbers, new formats of fort.41 and fort.44 based on value of
C     IDEN, and new fort.80 file format resulting from these changes.

C     jgf For use with ADCRIC v45.12 03/17/2006.

C     jgf Updated for 46.00 April 2006 to loop infinitely over data sets
C     in ascii files. This was done because hot started ascii files will
C     necessarily have the wrong number of data sets in the header.


C---------------------------------------------------------------------------
C                S U B R O U T I N E   P O S T _ I N I T 
C---------------------------------------------------------------------------
C                     (  Serial Version  3/28/98  )                          
C  This routine reads the domain decomposition information from file,       
C  "fort.80", which was written by the ADCIRC pre-processor ADCPREP.        
C                                                                           
C     jgf45.11 Updated to handle IDEN and 3D recording stations defined by
C     coordinates rather than node numbers.
C
C---------------------------------------------------------------------------
      SUBROUTINE POST_INIT()
      USE POST_GLOBAL
      IMPLICIT NONE
      INTEGER I,J,K
      CHARACTER DUMMY*80
C
C--Read the domain-decomposition information file created
C  by the Pre-processor
C
      OPEN(UNIT=80,FILE='fort.80')        
C
      READ(80,80) RUNDES
      READ(80,80) RUNID
      READ(80,80) AGRID
      READ(80,*)  NELG,NNODG
      READ(80,*)  NPROC
      READ(80,*)  MNPP
      READ(80,*)  MNEP !jgf45.07 max number of elements in any subdomain
      READ(80,*) IM    !jgf46.02 for post processing concentration sta. files
      READ(80,*) NWS   !jgf46.02 for post processing wind sta. files
      READ(80,*) NSTAE
      READ(80,*) NSTAV
      IF (IM.EQ.10) THEN
         READ(80,*) NSTAC  !jgf46.02 for post processing conc. sta. files
         MNSTAC = NSTAC
      ENDIF
      IF (NWS.NE.0) THEN
         READ(80,*) NSTAM  !jgf46.02 for post processing wind sta. files
         MNSTAM = NSTAM
      ENDIF
      READ(80,*) MNHARF   
      READ(80,*) MNWLAT,MNWLON
      MNSTAE = NSTAE
      MNSTAV = NSTAV
      MNP = NNODG
      MNE = NELG
      MNPROC = NPROC
C 
C  Allocate all global variables
      CALL ALLOC_MAIN1()
Casey 100301: Changed I8 to I12.
      DO I = 1,NPROC
         READ(80,*) J, NNODP(I), NOD_RES_TOT(I) 
         READ(80,'(9I12)') (IMAP_NOD_LG(K,I),K=1,NNODP(I))
      ENDDO
C
      READ(80,80) DUMMY
      DO I = 1,NNODG
         READ(80,1140) J, IMAP_NOD_GL(1,I), IMAP_NOD_GL(2,I)
         IMAP_NOD_GL(1,I) = IMAP_NOD_GL(1,I)+1
      ENDDO
C
C     jgf45.07 Add subdomain->fulldomain element mapping for NOFF processing
C     IMAP_EL_LG(I,PE) = Global Element Number of Local Element I on PE
Casey 100301: Changed I8 to I12.
      DO I = 1,NPROC
         READ(80,*) J, NELP(I)
         READ(80,'(9I12)') (IMAP_EL_LG(K,I),K=1,NELP(I))
      ENDDO

      READ(80,*) NOUTE,TOUTSE,TOUTFE,NSPOOLE 
C
C     jgf48.04 Added absolute value to station mappings.
      DO I = 1,NPROC
         READ(80,*) J,NSTAEP(I)
         DO K = 1,NSTAEP(I)
            READ(80,*) IMAP_STAE_LG(K,I)
            IMAP_STAE_LG(K,I) = ABS(IMAP_STAE_LG(K,I)) !jgf48.04
         ENDDO
      ENDDO
C
      READ(80,*) NOUTV,TOUTSV,TOUTFV,NSPOOLV
C
      DO I = 1,NPROC
         READ(80,*) J,NSTAVP(I)
         DO K = 1,NSTAVP(I)
            READ(80,*) IMAP_STAV_LG(K,I)
            IMAP_STAV_LG(K,I) = ABS(IMAP_STAV_LG(K,I)) !jgf48.04
         ENDDO
      ENDDO
C
      IF (IM.EQ.10) THEN        ! jgf46.02
         READ(80,*) NOUTC,TOUTSC,TOUTFC,NSPOOLC
         DO I = 1,NPROC
            READ(80,*) J,NSTACP(I)
            DO K = 1,NSTACP(I)
               READ(80,*) IMAP_STAC_LG(K,I)
               IMAP_STAC_LG(K,I) = ABS(IMAP_STAC_LG(K,I)) !jgf48.04
            ENDDO
         ENDDO
      ENDIF
C
      IF (NWS.NE.0) THEN  ! jgf46.02
         READ(80,*) NOUTM,TOUTSM,TOUTFM,NSPOOLM
         DO I = 1,NPROC
            READ(80,*) J,NSTAMP(I)
            DO K = 1,NSTAMP(I)
               READ(80,*) IMAP_STAM_LG(K,I)
               IMAP_STAM_LG(K,I) = ABS(IMAP_STAM_LG(K,I)) !jgf48.04
            ENDDO
         ENDDO
      ENDIF
C     
      READ(80,*) NOUTGE, TOUTSGE,TOUTFGE,NSPOOLGE
      READ(80,*) NOUTGV, TOUTSGV,TOUTFGV,NSPOOLGV
      READ(80,*) NOUTGC, TOUTSGC,TOUTFGC,NSPOOLGC
      READ(80,*) NOUTGW, TOUTSGW,TOUTFGW,NSPOOLGW
C
      READ(80,*) NHASE,NHASV,NHAGE,NHAGV
C     -------------------------------------------------------------
C     
C     S T A R T   3 D   D A T A 
C     
C     -------------------------------------------------------------

      READ(80,*) IDEN !jgf45.11 needed to post process the fort.44 file
C     -------------------------------------------------------------
C     jgf45.11 Read mappings for 3D density stations.
C     -------------------------------------------------------------
      READ(80,*) I3DSD, TO3DSDS, TO3DSDF, NSPO3DSD, NSTA3DD
      IF(NSTA3DD.GT.0)
     &  ALLOCATE ( NNSTA3DDP(MNPROC), IMAP_STA3DD_LG(NSTA3DD,MNPROC) )
      IF(I3DSD.NE.0) THEN
        DO I = 1, NPROC
           READ(80,*) J, NNSTA3DDP(I)
           DO K = 1, NNSTA3DDP(I)
              READ(80,*) IMAP_STA3DD_LG(K,I)
           ENDDO
        ENDDO
      ENDIF
C     -------------------------------------------------------------
C     jgf45.11 Read mappings for 3D velocity stations.
C     -------------------------------------------------------------
      READ(80,*) I3DSV, TO3DSVS, TO3DSVF, NSPO3DSV, NSTA3DV
      IF(NSTA3DV.GT.0)
     &  ALLOCATE ( NNSTA3DVP(MNPROC), IMAP_STA3DV_LG(NSTA3DV,MNPROC) )
      IF(I3DSV.NE.0) THEN
        DO I = 1, NPROC
           READ(80,*) J, NNSTA3DVP(I)
           DO K = 1, NNSTA3DVP(I)
              READ(80,*) IMAP_STA3DV_LG(K,I)
           ENDDO
        ENDDO
      ENDIF
C     -------------------------------------------------------------
C     jgf45.11 Read mappings for 3D turbulence stations.
C     -------------------------------------------------------------
      READ(80,*) I3DST, TO3DSTS, TO3DSTF, NSPO3DST, NSTA3DT
      IF(NSTA3DT.GT.0)
     &  ALLOCATE ( NNSTA3DTP(MNPROC), IMAP_STA3DT_LG(NSTA3DT,MNPROC) )
      IF(I3DST.NE.0) THEN
        DO I = 1, NPROC
           READ(80,*) J, NNSTA3DTP(I)
           write(*,*) j, i, nnsta3dtp(i)
           DO K = 1, NNSTA3DTP(I)
              READ(80,*) IMAP_STA3DT_LG(K,I)
           ENDDO
        ENDDO
      ENDIF
      READ(80,*) I3DGD,TO3DGDS,TO3DGDF,NSPO3DGD
      READ(80,*) I3DGV,TO3DGVS,TO3DGVF,NSPO3DGV
      READ(80,*) I3DGT,TO3DGTS,TO3DGTF,NSPO3DGT
C
C     End 3D data
C
      READ(80,*) NBYTE
C
      CLOSE(80)
C
  80  FORMAT(A80)
1130  FORMAT(8X,9I8)
1140  FORMAT(8X,3I8)
C
      RETURN
      END SUBROUTINE
C---------------------------------------------------------------------------
C     End of subroutine post_init
C---------------------------------------------------------------------------


      SUBROUTINE POST61()
C
C---------------------------------------------------------------------------C
C                     (  Serial Version  3/28/98  )                         C 
C  Globalize the elevation data at the elevation stations from the local    C
C  fort.61 files.                                                           C
C  This version is compatible with ADCIRC version 34.03                     C
C---------------------------------------------------------------------------C
C
      USE POST_GLOBAL
      IMPLICIT NONE
      INTEGER I,J,K,L,IPROC,IDUM
      INTEGER NTRSPE,NSTEMP,ITSE,ITEMPE,NUMSTNS
      INTEGER OREC
      REAL(SZ) DTE
      REAL(8) TIMEOUTSE
      CHARACTER*80  OUTMSG
      CHARACTER*85  INLINE
      CHARACTER*4 RDES4(8),RID4(6),AID4(6)
      CHARACTER*8 RDES8(8),RID8(6),AID8(6)
      LOGICAL FOUND

      INTEGER,ALLOCATABLE :: LOC2(:),IREC(:)
      REAL(SZ),ALLOCATABLE ::   ETBIN(:),ETBINP(:)
      CHARACTER*80,ALLOCATABLE ::  ETASC(:),ETASCP(:)
      CHARACTER*14,ALLOCATABLE ::  LOCNAME(:)

      ALLOCATE ( LOC2(MNPROC),IREC(MNPROC) )
      ALLOCATE ( ETBIN(MNSTAE),ETBINP(MNSTAE))
      ALLOCATE ( ETASC(MNSTAE),ETASCP(MNSTAE))
      ALLOCATE ( LOCNAME(MNPROC))

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
C--Open All Local Sequential Formatted fort.61 file
C
      DO IPROC = 1,NPROC
         LOCNAME(IPROC) = 'PE0000/fort.61'
         WRITE(LOCNAME(IPROC)(3:6),'(I4.4)') (IPROC-1)
         LOC2(IPROC) = 105 + (IPROC-1)
         INQUIRE(FILE=LOCNAME(IPROC),EXIST=FOUND)
         IF (FOUND) THEN
           OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC))
         ELSE
           print *, "No Local fort.61 files found"
           RETURN
         ENDIF
         READ(LOC2(IPROC),'(A85)')  INLINE
         READ(LOC2(IPROC),3645) NTRSPE,NUMSTNS,DTE,NSTEMP,ITEMPE
         IF (NUMSTNS.NE.NSTAEP(IPROC)) THEN
            STOP 'Inconsistency in Number of Elevation Stations'
         ENDIF
      ENDDO
C
C--Open Global Sequential Formatted fort.61 file
C
      OPEN(UNIT=61,FILE='fort.61')
C
      WRITE(61,'(A85)') INLINE
      WRITE(61,3645) NTRSPE,NSTAE,DTE,NSPOOLE,ITEMPE
C
C     jgf46.00 Made the loop over data sets an infinite loop because the
C     number of data sets at the top of the file will be wrong if the
C     run was hot started.
C
      DO ! loop infinitely over data sets
C
         DO IPROC=1,NPROC
            READ(LOC2(IPROC),2120,END=9999) TIMEOUTSE,ITSE
            IF (NSTAEP(IPROC).GT.0) THEN
              READ(LOC2(IPROC),80,END=9999) (ETASCP(K),K=1,NSTAEP(IPROC))
              DO K=1,NSTAEP(IPROC)
                 ETASC(IMAP_STAE_LG(K,IPROC)) = ETASCP(K)
              ENDDO
            ENDIF
         ENDDO            
C
         WRITE(61,2120) TIMEOUTSE,ITSE
         DO I=1, NSTAE
            CALL NEWINDEX(ETASC(I),OUTMSG,I)
            WRITE(61,*) TRIM(OUTMSG) !jgf46.00 TRIM to cut down file size 
         ENDDO
C
      ENDDO                 
      GO TO 9999
C 
2000  CONTINUE
C
C--Open All Local and the Global Direct Access Binary fort.61 file
C
      DO IPROC = 1,NPROC
         LOCNAME(IPROC) = 'PE0000/fort.61'
         WRITE(LOCNAME(IPROC)(3:6),'(I4.4)') (IPROC-1)
         LOC2(IPROC) = 105 + (IPROC-1)
         INQUIRE(FILE=LOCNAME(IPROC),EXIST=FOUND)
         IF (FOUND) THEN
           OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC),ACCESS='DIRECT',
     &       RECL=NBYTE)
         ELSE
           print *, "No Local fort.61 files found"
           RETURN
         ENDIF
      ENDDO
      OPEN(61,FILE='fort.61',ACCESS='DIRECT',RECL=NBYTE)
C
C
C--Read RUNDES RUNID and AGRID from each Local file
C
      DO IPROC = 1,NPROC
         IREC(IPROC) = 0
         IF (NBYTE.EQ.4) THEN
           DO I=1,8
              READ(LOC2(IPROC),REC=IREC(IPROC)+I) RDES4(I)
           ENDDO
           IREC(IPROC)=IREC(IPROC)+8
           DO I=1,6
              READ(LOC2(IPROC),REC=IREC(IPROC)+I) RID4(I)
           ENDDO
           IREC(IPROC)=IREC(IPROC)+6
           DO I=1,6
              READ(LOC2(IPROC),REC=IREC(IPROC)+I) AID4(I)
           ENDDO
           IREC(IPROC)=IREC(IPROC)+6
         ENDIF
         IF (NBYTE.EQ.8) THEN
           DO I=1,4
              READ(LOC2(IPROC),REC=IREC(IPROC)+I) RDES8(I)
           ENDDO
           IREC(IPROC)=IREC(IPROC)+4
           DO I=1,3
              READ(LOC2(IPROC),REC=IREC(IPROC)+I) RID8(I)
           ENDDO
           IREC(IPROC)=IREC(IPROC)+3
           DO I=1,3
              READ(LOC2(IPROC),REC=IREC(IPROC)+I) AID8(I)
           ENDDO
           IREC(IPROC)=IREC(IPROC)+3
         ENDIF
      ENDDO
C
C--Write RUNDES RUNID and AGRID to Global File
C
      OREC = 0
      IF (NBYTE.EQ.4) THEN
        DO I=1,8
           WRITE(61,REC=OREC+I) RDES4(I)
        ENDDO
        OREC=OREC+8
        DO I=1,6
           WRITE(61,REC=OREC+I) RID4(I)
        ENDDO
        OREC=OREC+6
        DO I=1,6
           WRITE(61,REC=OREC+I) AID4(I)
        ENDDO
        OREC=OREC+6
      ENDIF
      IF (NBYTE.EQ.8) THEN
        DO I=1,4
           WRITE(61,REC=OREC+I) RDES8(I)
        ENDDO
        OREC=OREC+4
        DO I=1,3
           WRITE(61,REC=OREC+I) RID8(I)
        ENDDO
        OREC=OREC+3
        DO I=1,3
           WRITE(61,REC=OREC+I) AID8(I)
        ENDDO
        OREC=OREC+3
      ENDIF
C
C--Read NTRSPE, NSTAE, DT*NSPOOLE from each Local file
C  and then close files to flush file buffers
C
      DO IPROC = 1,NPROC
         READ(LOC2(IPROC),REC=IREC(IPROC)+1) NTRSPE
         READ(LOC2(IPROC),REC=IREC(IPROC)+2) NUMSTNS        
         READ(LOC2(IPROC),REC=IREC(IPROC)+3) DTE
         READ(LOC2(IPROC),REC=IREC(IPROC)+4) NSTEMP 
         READ(LOC2(IPROC),REC=IREC(IPROC)+5) ITEMPE
         IREC(IPROC) = IREC(IPROC)+5
         CLOSE(LOC2(IPROC))    ! Flush the Write Buffer
         IF (NUMSTNS.NE.NSTAEP(IPROC)) THEN
            STOP 'Inconsistency in Number of Elevation Stations'
         ENDIF
         OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC),ACCESS='DIRECT',RECL=NBYTE)
      ENDDO
C
C--Write same info to Global file and close it to flush buffer also
C
      WRITE(61,REC=OREC+1) NTRSPE
      WRITE(61,REC=OREC+2) NSTAE
      WRITE(61,REC=OREC+3) DTE
      WRITE(61,REC=OREC+4) NSPOOLE
      WRITE(61,REC=OREC+5) ITEMPE
      OREC = OREC+5
      CLOSE(61)         ! Flush the Write Buffer
      OPEN(61,FILE='fort.61',ACCESS='DIRECT',RECL=NBYTE)
C
      DO J=1,NTRSPE
C
         DO IPROC=1,NPROC
            READ(LOC2(IPROC),REC=IREC(IPROC)+1) TIMEOUTSE
            READ(LOC2(IPROC),REC=IREC(IPROC)+2) ITSE
            IREC(IPROC) = IREC(IPROC) + 2
            IF (NSTAEP(IPROC).GT.0) THEN
              DO K=1, NSTAEP(IPROC)
                 READ(LOC2(IPROC),REC=IREC(IPROC)+K) ETBINP(K)
              ENDDO
              DO K=1,NSTAEP(IPROC)
                 ETBIN(IMAP_STAE_LG(K,IPROC)) = ETBINP(K)
              ENDDO
              IREC(IPROC) = IREC(IPROC) + NSTAEP(IPROC)
            ENDIF
         ENDDO            
C
         WRITE(61,REC=OREC+1) TIMEOUTSE
         WRITE(61,REC=OREC+2) ITSE
         OREC = OREC + 2
         DO K=1, NSTAE
            WRITE(61,REC=OREC+K) ETBIN(K)
         ENDDO
         OREC = OREC + NSTAE
C
      ENDDO                 
C
C--Close the Global and Local fort.61 Files
C
9999  CONTINUE
      CLOSE(61)
      DO IPROC = 1,NPROC
         CLOSE(LOC2(IPROC))
      ENDDO
C
   80 FORMAT(A40)
 2120 FORMAT(2X,E20.10,5X,I10)
 3645 FORMAT(1X,I10,1X,I10,1X,E15.7,1X,I5,1X,I5)
C
      RETURN
      END


      SUBROUTINE POST62()
C
C---------------------------------------------------------------------------C
C                     (  Serial Version  3/28/98  )                         C 
C  Globalize the elevation data at the velocity stations from the local     C
C  fort.62 files.                                                           C
C  This version is compatible with ADCIRC version 34.03                     C
C---------------------------------------------------------------------------C
C
      USE POST_GLOBAL
      IMPLICIT NONE
      INTEGER I,J,K,L,IPROC,IDUM
      INTEGER NTRSPV,NSTEMP,ITSV,ITEMPV,NUMSTNS
      INTEGER OREC
      REAL(SZ) DTV
      REAL(8) TIMEOUTSV 
      CHARACTER*80 OUTMSG
      CHARACTER*85  INLINE
      CHARACTER*4 RDES4(8),RID4(6),AID4(6)
      CHARACTER*8 RDES8(8),RID8(6),AID8(6)
      LOGICAL FOUND

      INTEGER,ALLOCATABLE :: LOC2(:),IREC(:)
      REAL(SZ),ALLOCATABLE ::  UUBIN(:),VVBIN(:)
      REAL(SZ),ALLOCATABLE ::  UUBINP(:),VVBINP(:)
      CHARACTER*80,ALLOCATABLE :: UUASC(:),UUASCP(:)
      CHARACTER*14,ALLOCATABLE  ::  LOCNAME(:)

      ALLOCATE (LOC2(MNPROC),IREC(MNPROC))
      ALLOCATE (UUBIN(MNSTAV),VVBIN(MNSTAV))
      ALLOCATE (UUBINP(MNSTAV),VVBINP(MNSTAV))
      ALLOCATE (UUASC(MNSTAV),UUASCP(MNSTAV))
      ALLOCATE (LOCNAME(MNPROC))
C
C--Determine whether Unit 62 is Sequential Formatted or Direct Access Binary
C
      IF (ABS(NOUTV).EQ.1) THEN
        GO TO 1000
      ELSE
        GO TO 2000
      ENDIF
C
1000  CONTINUE
C
      DO IPROC = 1,NPROC
         LOCNAME(IPROC) = 'PE0000/fort.62'
         WRITE(LOCNAME(IPROC)(3:6),'(I4.4)') (IPROC-1)
         LOC2(IPROC) = 105 + (IPROC-1)
         INQUIRE(FILE=LOCNAME(IPROC),EXIST=FOUND)
         IF (FOUND) THEN
           OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC))
         ELSE
           print *, "No Local fort.62 files found"
           RETURN
         ENDIF
         READ(LOC2(IPROC),'(A85)')  INLINE
         READ(LOC2(IPROC),3645) NTRSPV,NUMSTNS,DTV,NSTEMP,ITEMPV
         IF (NUMSTNS.NE.NSTAVP(IPROC)) THEN
            STOP 'Inconsistency in Number of Velocity Stations'
         ENDIF
      ENDDO
C
C--Open Global fort.62 file
C
      OPEN(UNIT=62,FILE='fort.62')
C
      WRITE(62,'(A85)') INLINE
      WRITE(62,3645) NTRSPV,NSTAV,DTV,NSPOOLV,ITEMPV
C
C     jgf46.00 Made the loop over data sets an infinite loop because the
C     number of data sets at the top of the file will be wrong if the
C     run was hot started.
C
      DO ! loop infinitely over data sets
C
         DO IPROC=1,NPROC
            READ(LOC2(IPROC),2120,END=9999) TIMEOUTSV,ITSV
            IF (NSTAVP(IPROC).GT.0) THEN
              READ(LOC2(IPROC),80,END=9999) (UUASCP(K),K=1,
     &                                     NSTAVP(IPROC))
              DO K=1,NSTAVP(IPROC)
                 UUASC(IMAP_STAV_LG(K,IPROC)) = UUASCP(K)
              ENDDO
            ENDIF
         ENDDO            
C
         WRITE(62,2120) TIMEOUTSV,ITSV
         DO I=1, NSTAV
            CALL NEWINDEX(UUASC(I),OUTMSG,I)
            WRITE(62,*) TRIM(OUTMSG) !jgf46.00 TRIM to cut down file size 
         ENDDO
C
      ENDDO                 
      GO TO 9999
C 
2000  CONTINUE
C
C--Open All Local and the Global Direct Access Binary fort.62 file
C
      DO IPROC = 1,NPROC
         LOCNAME(IPROC) = 'PE0000/fort.62'
         WRITE(LOCNAME(IPROC)(3:6),'(I4.4)') (IPROC-1)
         LOC2(IPROC) = 105 + (IPROC-1)
         INQUIRE(FILE=LOCNAME(IPROC),EXIST=FOUND)
         IF (FOUND) THEN
           OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC),ACCESS='DIRECT',
     &       RECL=NBYTE)
         ELSE
           print *, "No Local fort.61 files found"
           RETURN
         ENDIF
      ENDDO
      OPEN(62,FILE='fort.62',ACCESS='DIRECT',RECL=NBYTE)
C
C--Read RUNDES RUNID and AGRID from each Local file
C
      DO IPROC = 1,NPROC
         IREC(IPROC) = 0
         IF (NBYTE.EQ.4) THEN
           DO I=1,8
              READ(LOC2(IPROC),REC=IREC(IPROC)+I) RDES4(I)
           ENDDO
           IREC(IPROC)=IREC(IPROC)+8
           DO I=1,6
              READ(LOC2(IPROC),REC=IREC(IPROC)+I) RID4(I)
           ENDDO
           IREC(IPROC)=IREC(IPROC)+6
           DO I=1,6
              READ(LOC2(IPROC),REC=IREC(IPROC)+I) AID4(I)
           ENDDO
           IREC(IPROC)=IREC(IPROC)+6
         ENDIF
         IF (NBYTE.EQ.8) THEN
           DO I=1,4
              READ(LOC2(IPROC),REC=IREC(IPROC)+I) RDES8(I)
           ENDDO
           IREC(IPROC)=IREC(IPROC)+4
           DO I=1,3
              READ(LOC2(IPROC),REC=IREC(IPROC)+I) RID8(I)
           ENDDO
           IREC(IPROC)=IREC(IPROC)+3
           DO I=1,3
              READ(LOC2(IPROC),REC=IREC(IPROC)+I) AID8(I)
           ENDDO
           IREC(IPROC)=IREC(IPROC)+3
         ENDIF
      ENDDO
C
C--Write RUNDES RUNID and AGRID to Global File
C
      OREC = 0
      IF (NBYTE.EQ.4) THEN
        DO I=1,8
           WRITE(62,REC=OREC+I) RDES4(I)
        ENDDO
        OREC=OREC+8
        DO I=1,6
           WRITE(62,REC=OREC+I) RID4(I)
        ENDDO
        OREC=OREC+6
        DO I=1,6
           WRITE(62,REC=OREC+I) AID4(I)
        ENDDO
        OREC=OREC+6
      ENDIF
      IF (NBYTE.EQ.8) THEN
        DO I=1,4
           WRITE(62,REC=OREC+I) RDES8(I)
        ENDDO
        OREC=OREC+4
        DO I=1,3
           WRITE(62,REC=OREC+I) RID8(I)
        ENDDO
        OREC=OREC+3
        DO I=1,3
           WRITE(62,REC=OREC+I) AID8(I)
        ENDDO
        OREC=OREC+3
      ENDIF
C
C--Read NTRSPV, NSTAV, DT*NSPOOLV from each Local file
C  and then close files to flush file buffers
C
      DO IPROC = 1,NPROC
         READ(LOC2(IPROC),REC=IREC(IPROC)+1) NTRSPV
         READ(LOC2(IPROC),REC=IREC(IPROC)+2) NUMSTNS
         READ(LOC2(IPROC),REC=IREC(IPROC)+3) DTV
         READ(LOC2(IPROC),REC=IREC(IPROC)+4) NSTEMP
         READ(LOC2(IPROC),REC=IREC(IPROC)+5) ITEMPV
         IREC(IPROC) = IREC(IPROC)+5
         CLOSE(LOC2(IPROC))    ! Flush the Write Buffer
         IF (NUMSTNS.NE.NSTAVP(IPROC)) THEN
            STOP 'Inconsistency in Number of Velocity Stations'
         ENDIF
         OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC),ACCESS='DIRECT',RECL=NBYTE)
      ENDDO
C
C--Write same info to Global file and close it to flush buffer also
C
      WRITE(62,REC=OREC+1) NTRSPV
      WRITE(62,REC=OREC+2) NSTAV
      WRITE(62,REC=OREC+3) DTV
      WRITE(62,REC=OREC+4) NSPOOLV
      WRITE(62,REC=OREC+5) ITEMPV
      OREC = OREC+5
      CLOSE(62)         ! Flush the Write Buffer
      OPEN(62,FILE='fort.62',ACCESS='DIRECT',RECL=NBYTE)
C
      DO J=1,NTRSPV
C
         DO IPROC=1,NPROC
            READ(LOC2(IPROC),REC=IREC(IPROC)+1) TIMEOUTSV
            READ(LOC2(IPROC),REC=IREC(IPROC)+2) ITSV
            IREC(IPROC) = IREC(IPROC) + 2
            IF (NSTAVP(IPROC).GT.0) THEN
              DO K=1, NSTAVP(IPROC)
                 READ(LOC2(IPROC),REC=IREC(IPROC)+2*K-1) UUBINP(K)
                 READ(LOC2(IPROC),REC=IREC(IPROC)+2*K) VVBINP(K)
              ENDDO
              DO K=1,NSTAVP(IPROC)
                 UUBIN(IMAP_STAV_LG(K,IPROC)) = UUBINP(K)
                 VVBIN(IMAP_STAV_LG(K,IPROC)) = VVBINP(K)
              ENDDO
             IREC(IPROC) = IREC(IPROC) + 2*NSTAVP(IPROC)
            ENDIF
         ENDDO
C
         WRITE(62,REC=OREC+1) TIMEOUTSV
         WRITE(62,REC=OREC+2) ITSV
         OREC = OREC + 2
         DO K=1, NSTAV
            WRITE(62,REC=OREC+2*K-1) UUBIN(K)
            WRITE(62,REC=OREC+2*K) VVBIN(K)
         ENDDO
         OREC = OREC + 2*NSTAV
C
      ENDDO
C
 9999 CONTINUE
C
C--Close the Global and Local fort.62 Files
C
      CLOSE(62)
      DO IPROC = 1,NPROC
         CLOSE(LOC2(IPROC))
      ENDDO
C
 2120 FORMAT(2X,E20.10,5X,I10)
   80 FORMAT(A80)
c2454 FORMAT(2X,I8,2X,E15.8,2X,E15.8)
 3645 FORMAT(1X,I10,1X,I10,1X,E15.7,1X,I5,1X,I5)
C
      RETURN
      END


      SUBROUTINE POST63()
C
C---------------------------------------------------------------------------C
C                     (  Serial Version  3/28/98  )                         C 
C  Globalize the elevation data at all nodes from the local fort.63 files.  C
C  This version is compatible with ADCIRC version 34.03                     C
C---------------------------------------------------------------------------C
C
      USE POST_GLOBAL
      IMPLICIT NONE
      INTEGER I,J,K,L,IPROC,IDUM
      INTEGER NDSETGE,NSTEMP,NP,ITE,ITEMPE,INDX
      INTEGER OREC
      REAL(SZ) DTE
      REAL(SZ) ETA
      REAL(8) TIMEOUTE 
      CHARACTER*80 OUTMSG
      CHARACTER*85 INLINE
      CHARACTER*4 RDES4(8),RID4(6),AID4(6)
      CHARACTER*8 RDES8(8),RID8(6),AID8(6)
      LOGICAL FOUND

      INTEGER, ALLOCATABLE  ::  LOC2(:),IREC(:)
      REAL(SZ), ALLOCATABLE  ::  ETBIN(:),ETBINP(:)
      CHARACTER*80, ALLOCATABLE  ::  ETASC(:),ETASCP(:)
      CHARACTER*14, ALLOCATABLE  ::  LOCNAME(:)

      ALLOCATE (LOC2(MNPROC),IREC(MNPROC))
      ALLOCATE (ETBIN(MNP),ETBINP(MNPP))
      ALLOCATE (ETASC(MNP),ETASCP(MNPP))
      ALLOCATE (LOCNAME(MNPROC))

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
C--Open Global Sequential Formatted fort.63 file
C
      DO IPROC = 1,NPROC
         LOCNAME(IPROC) = 'PE0000/fort.63'
         WRITE(LOCNAME(IPROC)(3:6),'(I4.4)') (IPROC-1)
         LOC2(IPROC) = 105 + (IPROC-1)
         INQUIRE(FILE=LOCNAME(IPROC),EXIST=FOUND)
         IF (FOUND) THEN
           OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC))
         ELSE
           print *, "No Local fort.63 files found"
           RETURN
         ENDIF
         READ(LOC2(IPROC),'(A85)')  INLINE
         READ(LOC2(IPROC),3645)  NDSETGE,NP,DTE,NSTEMP,ITEMPE
         IF (NP.NE.NNODP(IPROC)) THEN
           print *, "NP = ",NP," NNODP = ",NNODP(IPROC)
           STOP 'Inconsistency in number of local nodes'
         ENDIF
      ENDDO
      OPEN(UNIT=63,FILE='fort.63')
C
      WRITE(63,'(A85)') INLINE
      WRITE(63,3645) NDSETGE,NNODG,DTE,NSPOOLGE,ITEMPE
C
C     jgf46.00 Made the loop over data sets an infinite loop because the
C     number of data sets at the top of the file will be wrong if the
C     run was hot started.
C
      DO ! loop infinitely over data sets
C
         DO IPROC=1,NPROC
            READ(LOC2(IPROC),FMT=2120,END=9999) TIMEOUTE,ITE
            IF (IPROC.EQ.1) WRITE(63,2120) TIMEOUTE,ITE
            DO K = 1,NNODP(IPROC)
               READ(LOC2(IPROC),80,END=9999) ETASCP(K)
               INDX = IMAP_NOD_LG(K,IPROC)
               IF (IMAP_NOD_GL(1,INDX).EQ.IPROC) THEN
                 ETASC(INDX) = ETASCP(K)
                ENDIF
            ENDDO  
         ENDDO 
c        DO I = 1,NNODG
c           READ(ETASC(I),*) IDUM, ETA
c           WRITE(63,2453) I, ETA
c        ENDDO
         DO I = 1,NNODG
            CALL NEWINDEX(ETASC(I),OUTMSG,I)
            WRITE(63,*) TRIM(OUTMSG) !jgf46.00 TRIM to cut down file size 
         ENDDO 
C              
      ENDDO
      GO TO 9999
C
2000  CONTINUE
C
C--Open All Local and the Global Direct Access Binary fort.63 file
C
      DO IPROC = 1,NPROC
         LOCNAME(IPROC) = 'PE0000/fort.63'
         WRITE(LOCNAME(IPROC)(3:6),'(I4.4)') (IPROC-1)
         LOC2(IPROC) = 105 + (IPROC-1)
         INQUIRE(FILE=LOCNAME(IPROC),EXIST=FOUND)
         IF (FOUND) THEN
           OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC),ACCESS='DIRECT',
     &       RECL=NBYTE)
         ELSE
           print *, "No Local fort.63 files found"
           RETURN
         ENDIF
      ENDDO
      OPEN(63,FILE='fort.63',ACCESS='DIRECT',RECL=NBYTE)
C
C--Read RUNDES RUNID and AGRID from each Local file
C
      DO IPROC = 1,NPROC
         IREC(IPROC) = 0
         IF (NBYTE.EQ.4) THEN
           DO I=1,8
              READ(LOC2(IPROC),REC=IREC(IPROC)+I) RDES4(I)
           ENDDO
           IREC(IPROC)=IREC(IPROC)+8
           DO I=1,6
              READ(LOC2(IPROC),REC=IREC(IPROC)+I) RID4(I)
           ENDDO
           IREC(IPROC)=IREC(IPROC)+6
           DO I=1,6
              READ(LOC2(IPROC),REC=IREC(IPROC)+I) AID4(I)
           ENDDO
           IREC(IPROC)=IREC(IPROC)+6
         ENDIF
         IF (NBYTE.EQ.8) THEN
           DO I=1,4
              READ(LOC2(IPROC),REC=IREC(IPROC)+I) RDES8(I)
           ENDDO
           IREC(IPROC)=IREC(IPROC)+4
           DO I=1,3
              READ(LOC2(IPROC),REC=IREC(IPROC)+I) RID8(I)
           ENDDO
           IREC(IPROC)=IREC(IPROC)+3
           DO I=1,3
              READ(LOC2(IPROC),REC=IREC(IPROC)+I) AID8(I)
           ENDDO
           IREC(IPROC)=IREC(IPROC)+3
         ENDIF
      ENDDO
C
C--Write RUNDES RUNID and AGRID to Global File
C
      OREC = 0
      IF (NBYTE.EQ.4) THEN
        DO I=1,8
           WRITE(63,REC=OREC+I) RDES4(I)
        ENDDO
        OREC=OREC+8
        DO I=1,6
           WRITE(63,REC=OREC+I) RID4(I)
        ENDDO
        OREC=OREC+6
        DO I=1,6
           WRITE(63,REC=OREC+I) AID4(I)
        ENDDO
        OREC=OREC+6
      ENDIF
      IF (NBYTE.EQ.8) THEN
        DO I=1,4
           WRITE(63,REC=OREC+I) RDES8(I)
        ENDDO
        OREC=OREC+4
        DO I=1,3
           WRITE(63,REC=OREC+I) RID8(I)
        ENDDO
        OREC=OREC+3
        DO I=1,3
           WRITE(63,REC=OREC+I) AID8(I)
        ENDDO
        OREC=OREC+3
      ENDIF
C
C--Read NDSETGE, NNODP, DTE, NSPOOLGE, ITEMPE from each Local file
C  and then close files to flush file buffers
C
      DO IPROC = 1,NPROC
         READ(LOC2(IPROC),REC=IREC(IPROC)+1) NDSETGE
         READ(LOC2(IPROC),REC=IREC(IPROC)+2) IDUM   
         READ(LOC2(IPROC),REC=IREC(IPROC)+3) DTE
         READ(LOC2(IPROC),REC=IREC(IPROC)+4) NSTEMP
         READ(LOC2(IPROC),REC=IREC(IPROC)+5) ITEMPE
         IREC(IPROC) = IREC(IPROC)+5
         CLOSE(LOC2(IPROC))    ! Flush the Write Buffer
         OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC),ACCESS='DIRECT',RECL=NBYTE)
      ENDDO
C
C--Write same info to Global file and close it to flush buffer also
C
      WRITE(63,REC=OREC+1) NDSETGE
      WRITE(63,REC=OREC+2) NNODG
      WRITE(63,REC=OREC+3) DTE
      WRITE(63,REC=OREC+4) NSPOOLGE
      WRITE(63,REC=OREC+5) ITEMPE
      OREC = OREC+5
      CLOSE(63)         ! Flush the Write Buffer
      OPEN(63,FILE='fort.63',ACCESS='DIRECT',RECL=NBYTE)
C
      DO J=1,NDSETGE
C
         DO IPROC=1,NPROC
            READ(LOC2(IPROC),REC=IREC(IPROC)+1) TIMEOUTE
            READ(LOC2(IPROC),REC=IREC(IPROC)+2) ITE
            IREC(IPROC) = IREC(IPROC) + 2
            DO K=1, NNODP(IPROC)
               READ(LOC2(IPROC),REC=IREC(IPROC)+K) ETBINP(K)
            ENDDO
            DO K=1,NNODP(IPROC)
               INDX = IMAP_NOD_LG(K,IPROC)
               IF (IMAP_NOD_GL(1,INDX).EQ.IPROC) THEN
                 ETBIN(INDX) = ETBINP(K)
               ENDIF
            ENDDO
            IREC(IPROC) = IREC(IPROC) + NNODP(IPROC)
         ENDDO
C
         WRITE(63,REC=OREC+1) TIMEOUTE
         WRITE(63,REC=OREC+2) ITE
         OREC = OREC + 2
         DO K=1, NNODG
            WRITE(63,REC=OREC+K) ETBIN(K)
         ENDDO
         OREC = OREC + NNODG
C
      ENDDO
C 
 9999 CONTINUE
C
C--Close the Global and Local fort.63 Files
C
      CLOSE(63)
      DO IPROC = 1,NPROC
         CLOSE(LOC2(IPROC))
      ENDDO
C
   80 FORMAT(A80)
 2120 FORMAT(2X,E20.10,5X,I10)
 2453 FORMAT(2X,I8,2X,E15.8)
 3645 FORMAT(1X,I10,1X,I10,1X,E15.7,1X,I5,1X,I5)
C
      RETURN
      END

      SUBROUTINE POST64()
C
C---------------------------------------------------------------------------C
C                     (  Serial Version  3/28/98  )                         C 
C  Globalize the Velocity Data at all nodes from the local fort.64 files.   C
C  This version is compatible with ADCIRC version 34.03                     C
C---------------------------------------------------------------------------C
C
      USE POST_GLOBAL
      IMPLICIT NONE
      INTEGER I,J,K,L,IPROC,IDUM
      INTEGER NDSETGV,NSTEMP,NP,ITV,ITEMPV,INDX
      INTEGER OREC
      REAL(SZ) DTV
      REAL(SZ) U1,V1
      REAL(8) TIMEOUTV
      CHARACTER*80  OUTMSG
      CHARACTER*85 INLINE
      CHARACTER*4 RDES4(8),RID4(6),AID4(6)
      CHARACTER*8 RDES8(8),RID8(6),AID8(6)
      LOGICAL FOUND

      INTEGER,ALLOCATABLE  :: LOC2(:),IREC(:)
      REAL(SZ),ALLOCATABLE  ::  UUBIN(:),VVBIN(:)
      REAL(SZ),ALLOCATABLE  ::  UUBINP(:),VVBINP(:)
      CHARACTER*80,ALLOCATABLE  ::  UU(:),UUP(:)
      CHARACTER*14,ALLOCATABLE  ::  LOCNAME(:)

      ALLOCATE (LOC2(MNPROC),IREC(MNPROC))
      ALLOCATE (UUBIN(MNP),VVBIN(MNP))
      ALLOCATE (UUBINP(MNPP),VVBINP(MNPP))
      ALLOCATE (UU(MNP),UUP(MNPP))
      ALLOCATE (LOCNAME(MNPROC))

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
C--Open Global fort.64 file
C
      DO IPROC = 1,NPROC
         LOCNAME(IPROC) = 'PE0000/fort.64'
         WRITE(LOCNAME(IPROC)(3:6),'(I4.4)') (IPROC-1)
         LOC2(IPROC) = 105 + (IPROC-1)
         INQUIRE(FILE=LOCNAME(IPROC),EXIST=FOUND)
         IF (FOUND) THEN
           OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC))
         ELSE
           print *, "No Local fort.64 files found"
           RETURN
         ENDIF
         READ(LOC2(IPROC),'(A85)')  INLINE
         READ(LOC2(IPROC),3645) NDSETGV,NP,DTV,NSTEMP,ITEMPV
         IF (NP.NE.NNODP(IPROC)) THEN
           STOP 'Inconsistency in number of local nodes'
         ENDIF
      ENDDO
      OPEN(UNIT=64,FILE='fort.64')
C
      WRITE(64,'(A85)') INLINE
      WRITE(64,3645) NDSETGV,NNODG,DTV,NSPOOLGV,ITEMPV
C
C     jgf46.00 Made the loop over data sets an infinite loop because the
C     number of data sets at the top of the file will be wrong if the
C     run was hot started.
C
      DO ! loop infinitely over data sets
C
         DO IPROC=1,NPROC
            READ(LOC2(IPROC),FMT=2120,END=9999) TIMEOUTV,ITV
            IF (IPROC.EQ.1) WRITE(64,2120) TIMEOUTV,ITV
            DO K = 1,NNODP(IPROC)
               READ(LOC2(IPROC),80,END=9999) UUP(K)
               INDX = IMAP_NOD_LG(K,IPROC)
               IF(IMAP_NOD_GL(1,INDX).EQ.IPROC) THEN
                UU(INDX) = UUP(K)
               ENDIF
            ENDDO
         ENDDO 
c        DO I = 1,NNODG
c           CALL NEWINDEX(UU(I),OUTMSG,I)
c           READ(UU(I),*)  IDUM, U1,V1
c           WRITE(64,2454) I, U1, V1
c        ENDDO
         DO I = 1,NNODG
            CALL NEWINDEX(UU(I),OUTMSG,I)
            WRITE(64,*) TRIM(OUTMSG)  !jgf46.00 TRIM to cut down file size 
         ENDDO
      ENDDO               
      GO TO 9999

C
2000  CONTINUE
C
C--Open All Local and the Global Direct Access Binary fort.64 file
C
      DO IPROC = 1,NPROC
         LOCNAME(IPROC) = 'PE0000/fort.64'
         WRITE(LOCNAME(IPROC)(3:6),'(I4.4)') (IPROC-1)
         LOC2(IPROC) = 105 + (IPROC-1)
         INQUIRE(FILE=LOCNAME(IPROC),EXIST=FOUND)
         IF (FOUND) THEN
           OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC),ACCESS='DIRECT',
     &       RECL=NBYTE)
         ELSE
           print *, "No Local fort.64 files found"
           RETURN
         ENDIF
      ENDDO
      OPEN(64,FILE='fort.64',ACCESS='DIRECT',RECL=NBYTE)
C
C--Read RUNDES RUNID and AGRID from each Local file
C
      DO IPROC = 1,NPROC
         IREC(IPROC) = 0
         IF (NBYTE.EQ.4) THEN
           DO I=1,8
              READ(LOC2(IPROC),REC=IREC(IPROC)+I) RDES4(I)
           ENDDO
           IREC(IPROC)=IREC(IPROC)+8
           DO I=1,6
              READ(LOC2(IPROC),REC=IREC(IPROC)+I) RID4(I)
           ENDDO
           IREC(IPROC)=IREC(IPROC)+6
           DO I=1,6
              READ(LOC2(IPROC),REC=IREC(IPROC)+I) AID4(I)
           ENDDO
           IREC(IPROC)=IREC(IPROC)+6
         ENDIF
         IF (NBYTE.EQ.8) THEN
           DO I=1,4
              READ(LOC2(IPROC),REC=IREC(IPROC)+I) RDES8(I)
           ENDDO
           IREC(IPROC)=IREC(IPROC)+4
           DO I=1,3
              READ(LOC2(IPROC),REC=IREC(IPROC)+I) RID8(I)
           ENDDO
           IREC(IPROC)=IREC(IPROC)+3
           DO I=1,3
              READ(LOC2(IPROC),REC=IREC(IPROC)+I) AID8(I)
           ENDDO
           IREC(IPROC)=IREC(IPROC)+3
         ENDIF
      ENDDO
C
C--Write RUNDES RUNID and AGRID to Global File
C
      OREC = 0
      IF (NBYTE.EQ.4) THEN
        DO I=1,8
           WRITE(64,REC=OREC+I) RDES4(I)
        ENDDO
        OREC=OREC+8
        DO I=1,6
           WRITE(64,REC=OREC+I) RID4(I)
        ENDDO
        OREC=OREC+6
        DO I=1,6
           WRITE(64,REC=OREC+I) AID4(I)
        ENDDO
        OREC=OREC+6
      ENDIF
      IF (NBYTE.EQ.8) THEN
        DO I=1,4
           WRITE(64,REC=OREC+I) RDES8(I)
        ENDDO
        OREC=OREC+4
        DO I=1,3
           WRITE(64,REC=OREC+I) RID8(I)
        ENDDO
        OREC=OREC+3
        DO I=1,3
           WRITE(64,REC=OREC+I) AID8(I)
        ENDDO
        OREC=OREC+3
      ENDIF
C
C--Read NTRSPV, NSTAV, DT*NSPOOLV from each Local file
C  and then close files to flush file buffers
C
      DO IPROC = 1,NPROC
         READ(LOC2(IPROC),REC=IREC(IPROC)+1) NDSETGV
         READ(LOC2(IPROC),REC=IREC(IPROC)+2) NP
         READ(LOC2(IPROC),REC=IREC(IPROC)+3) DTV
         READ(LOC2(IPROC),REC=IREC(IPROC)+4) NSTEMP
         READ(LOC2(IPROC),REC=IREC(IPROC)+5) ITEMPV
         IREC(IPROC) = IREC(IPROC)+5
         CLOSE(LOC2(IPROC))    ! Flush the Write Buffer
         IF (NP.NE.NNODP(IPROC)) THEN
           STOP 'Inconsistency in number of local nodes'
         ENDIF
         OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC),ACCESS='DIRECT',RECL=NBYTE)
      ENDDO
C
C--Write same info to Global file and close it to flush buffer also
C
      WRITE(64,REC=OREC+1) NDSETGV
      WRITE(64,REC=OREC+2) NNODG
      WRITE(64,REC=OREC+3) DTV
      WRITE(64,REC=OREC+4) NSPOOLGV
      WRITE(64,REC=OREC+5) ITEMPV
      OREC = OREC+5
      CLOSE(64)         ! Flush the Write Buffer
      OPEN(64,FILE='fort.64',ACCESS='DIRECT',RECL=NBYTE)
C
      DO J=1,NDSETGV
C
         DO IPROC=1,NPROC
            READ(LOC2(IPROC),REC=IREC(IPROC)+1) TIMEOUTV
            READ(LOC2(IPROC),REC=IREC(IPROC)+2) ITV
            IREC(IPROC) = IREC(IPROC) + 2
            DO K=1, NNODP(IPROC)
               READ(LOC2(IPROC),REC=IREC(IPROC)+2*K-1) UUBINP(K)
               READ(LOC2(IPROC),REC=IREC(IPROC)+2*K) VVBINP(K)
            ENDDO
            DO K=1,NNODP(IPROC)
               UUBIN(IMAP_NOD_LG(K,IPROC)) = UUBINP(K)
               VVBIN(IMAP_NOD_LG(K,IPROC)) = VVBINP(K)
            ENDDO
           IREC(IPROC) = IREC(IPROC) + 2*NNODP(IPROC)
         ENDDO
C
         WRITE(64,REC=OREC+1) TIMEOUTV
         WRITE(64,REC=OREC+2) ITV
         OREC = OREC + 2
         DO K=1, NNODG
            WRITE(64,REC=OREC+2*K-1) UUBIN(K)
            WRITE(64,REC=OREC+2*K)   VVBIN(K)
         ENDDO
         OREC = OREC + 2*NNODG
C
      ENDDO
C
9999  CONTINUE
C
C--Close the Global and Local fort.64 Files
C
      CLOSE(64)
      DO IPROC = 1,NPROC
         CLOSE(LOC2(IPROC))
      ENDDO
C
   80 FORMAT(A80)
 2120 FORMAT(2X,E20.10,5X,I10)
 2454 FORMAT(2X,I8,2(2X,E15.8))
 3645 FORMAT(1X,I10,1X,I10,1X,E15.7,1X,I5,1X,I5)
C
      RETURN
      END

      SUBROUTINE POST51()
C
C---------------------------------------------------------------------------C
C                     (  Serial Version  3/28/98  )                         C 
C  Globalize the harmonic data at the elevation stations from the local     C
C  fort.51 files.                                                           C
C  This version is compatible with ADCIRC version 34.03                     C
C  Fixed routine to handled multiple frequencies up to 32--9/18/98vjp       C
C---------------------------------------------------------------------------C
C
      USE POST_GLOBAL
      IMPLICIT NONE
      INTEGER I,J,IPROC,IDUM,IFREQ,INDX,IREC
      INTEGER NFREQ,NUMSTNS
      CHARACTER*80  NFREQMSG,INLINE
      LOGICAL FOUND


      INTEGER,ALLOCATABLE  :: LOC2(:)
      CHARACTER*14,ALLOCATABLE  ::  LOCNAME(:)
      CHARACTER*80,ALLOCATABLE  ::  HEADER(:),HARDAT(:)

      ALLOCATE ( LOC2(MNPROC))
      ALLOCATE ( LOCNAME(MNPROC))
      ALLOCATE ( HEADER(MNHARF),HARDAT(MNHARF))

C
C--If they exist open all Local & Global Sequential Formatted fort.51 file
C
      DO IPROC = 1,NPROC
         LOCNAME(IPROC) = 'PE0000/fort.51'
         WRITE(LOCNAME(IPROC)(3:6),'(I4.4)') (IPROC-1)
         LOC2(IPROC) = 105+ (IPROC-1)
         INQUIRE(FILE=LOCNAME(IPROC),EXIST=FOUND)
         IF (FOUND) THEN
           OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC))
         ELSE
           print *, "No Local unit 51 files found"
           RETURN
         ENDIF
      ENDDO
C
C   This is the global output file
C
      OPEN(UNIT=51,FILE='fort.51')
C
C   This is a DASD scratch file to save on memory.
C
      OPEN(UNIT=8,ACCESS='DIRECT',FORM='UNFORMATTED',RECL=80) 
C
      DO IPROC = 1,NPROC
         IF (IPROC.EQ.1) THEN
           READ(LOC2(IPROC),'(A80)')  NFREQMSG
           READ(NFREQMSG,*) NFREQ
           WRITE(51,'(A80)') NFREQMSG
           DO I=1, NFREQ  
             READ(LOC2(IPROC),'(A80)') HEADER(I)
             WRITE(51,*) trim(HEADER(I))
           ENDDO
           READ(LOC2(IPROC),*) NUMSTNS
           WRITE(51,*) NSTAE
         ELSE
           READ(LOC2(IPROC),'(A80)')  INLINE
           READ(INLINE,*) IDUM
           DO I=1, NFREQ  
             READ(LOC2(IPROC),'(A80)') INLINE
           ENDDO
           READ(LOC2(IPROC),*) IDUM
         ENDIF
      ENDDO
C
      DO IPROC = 1,NPROC
         DO J=1,NSTAEP(IPROC)
            READ(LOC2(IPROC),*) IDUM
            INDX = IMAP_STAE_LG(J,IPROC)
            IF (INDX.NE.0) THEN
              DO IFREQ=1, NFREQ
                 READ(LOC2(IPROC),'(A80)') HARDAT(IFREQ)
                 IREC = IFREQ + (INDX-1)*NFREQ
                 WRITE(8,REC=IREC) HARDAT(IFREQ)
              ENDDO
             ENDIF
          ENDDO
      ENDDO
C
C   Close the DASD file to flush the write buffer and re-open.
C
      CLOSE(8,STATUS='KEEP')
      OPEN(UNIT=8,ACCESS='DIRECT',FORM='UNFORMATTED',RECL=80) 
C
      DO J=1, NSTAE
          WRITE(51,*) J
          DO IFREQ=1, NFREQ
             IREC = IFREQ + (J-1)*NFREQ
             READ(8,REC=IREC) HARDAT(IFREQ)
             WRITE(51,*) trim(HARDAT(IFREQ))
          ENDDO
      ENDDO        
C
C--Close the Global and Local fort.51 Files
C
9999  CONTINUE
      CLOSE(51)
      CLOSE(8,STATUS='DELETE')
      DO IPROC = 1,NPROC
         CLOSE(LOC2(IPROC))
      ENDDO
C
c3679 FORMAT(1X,E20.10,1X,F10.7,1X,F12.8,1X,A10)
C
      RETURN
      END

      SUBROUTINE POST52()
C
C---------------------------------------------------------------------------C
C                     (  Serial Version  3/28/98  )                         C 
C  Globalize the harmonic data at the velocity stations from the local      C
C  fort.52 files.                                                           C
C  This version is compatible with ADCIRC version 34.03                     C
C  Fixed routine to handled multiple frequencies up to 32--9/18/98vjp       C
C---------------------------------------------------------------------------C
C
      USE POST_GLOBAL
      IMPLICIT NONE
      INTEGER I,J,IPROC,IDUM,IFREQ,INDX,IREC
      INTEGER NFREQ,NUMSTNS
      CHARACTER*80  NFREQMSG,INLINE
      LOGICAL FOUND


      INTEGER,ALLOCATABLE  :: LOC2(:)
      CHARACTER*14,ALLOCATABLE  ::  LOCNAME(:)
      CHARACTER*80,ALLOCATABLE  ::  HEADER(:),HARDAT(:)

      ALLOCATE ( LOC2(MNPROC))
      ALLOCATE ( LOCNAME(MNPROC))
      ALLOCATE ( HEADER(MNHARF),HARDAT(MNHARF))


C
C--If they exist open all Local & Global Sequential Formatted fort.52 file
C
      DO IPROC = 1,NPROC
         LOCNAME(IPROC) = 'PE0000/fort.52'
         WRITE(LOCNAME(IPROC)(3:6),'(I4.4)') (IPROC-1)
         LOC2(IPROC) = 105 + (IPROC-1)
         INQUIRE(FILE=LOCNAME(IPROC),EXIST=FOUND)
         IF (FOUND) THEN
           OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC))
         ELSE
           print *, "No Local unit 52 files found"
           RETURN
         ENDIF
      ENDDO
C
C   This is the global output file
C
      OPEN(UNIT=52,FILE='fort.52')
C
C   This is a DASD scratch file to save on memory.
C
      OPEN(UNIT=8,ACCESS='DIRECT',FORM='UNFORMATTED',RECL=80)
C
      DO IPROC = 1,NPROC
         IF (IPROC.EQ.1) THEN
           READ(LOC2(IPROC),'(A80)')  NFREQMSG
           READ(NFREQMSG,*) NFREQ
           WRITE(52,'(A80)') NFREQMSG
           DO I=1, NFREQ  
             READ(LOC2(IPROC),'(A80)') HEADER(I)
             WRITE(52,'(A80)') HEADER(I)
           ENDDO
           READ(LOC2(IPROC),*) NUMSTNS
           WRITE(52,*) NSTAV
         ELSE
           READ(LOC2(IPROC),'(A80)')  INLINE
           READ(INLINE,*) IDUM
           DO I=1, NFREQ  
             READ(LOC2(IPROC),'(A80)') INLINE
           ENDDO
           READ(LOC2(IPROC),*) IDUM
         ENDIF
      ENDDO
C
      DO IPROC = 1,NPROC
         DO J=1,NSTAVP(IPROC)
            READ(LOC2(IPROC),*) IDUM
            INDX = IMAP_STAV_LG(J,IPROC)
            IF (INDX.NE.0) THEN
              DO IFREQ=1, NFREQ
                 READ(LOC2(IPROC),'(A80)') HARDAT(IFREQ)
                 IREC = IFREQ + (INDX-1)*NFREQ
                 WRITE(8,REC=IREC) HARDAT(IFREQ)
              ENDDO
            ENDIF
           ENDDO
      ENDDO
C
C   Close the DASD file to flush the write buffer and re-open.
C
      CLOSE(8,STATUS='KEEP')
      OPEN(UNIT=8,ACCESS='DIRECT',FORM='UNFORMATTED',RECL=80) 
C
      DO J=1, NSTAV
          WRITE(52,*) J
          DO IFREQ=1, NFREQ
             IREC = IFREQ + (J-1)*NFREQ
             READ(8,REC=IREC) HARDAT(IFREQ)
             WRITE(52,'(A80)') HARDAT(IFREQ)
          ENDDO
      ENDDO              
C
C--Close the Global and Local fort.52 Files
C
9999  CONTINUE
      CLOSE(52)
      CLOSE(8,STATUS='DELETE')
      DO IPROC = 1,NPROC
         CLOSE(LOC2(IPROC))
      ENDDO
C
c3679 FORMAT(1X,E20.10,1X,F10.7,1X,F12.8,1X,A10)
C
      RETURN
      END

      SUBROUTINE POST53()
C
C---------------------------------------------------------------------------C
C                     (  Serial Version  3/28/98  )                         C 
C  Globalize the harmonic constituent elevations at all nodes from the      C
C  local fort.53 files.                                                     C
C  This version is compatible with ADCIRC version 34.03                     C
C  Fixed routine to handled arbitrary number of frequencies  1/30/99vjp     C
C---------------------------------------------------------------------------C
C
      USE POST_GLOBAL
      IMPLICIT NONE
      INTEGER I,J,IPROC,IDUM,INDX,NP1,IFREQ,IREC
      INTEGER NFREQ
      CHARACTER*80  NFREQMSG,INLINE
      LOGICAL FOUND


      INTEGER,ALLOCATABLE  :: LOC2(:)
      CHARACTER*14,ALLOCATABLE  ::  LOCNAME(:)
      CHARACTER*80,ALLOCATABLE  ::  HEADER(:),HARDAT(:)

      ALLOCATE ( LOC2(MNPROC))
      ALLOCATE ( LOCNAME(MNPROC))
      ALLOCATE (  HEADER(MNHARF),HARDAT(MNHARF))

C
C--If they exist open all Local & Global Sequential Formatted fort.53 file
C
      DO IPROC = 1,NPROC
         LOCNAME(IPROC) = 'PE0000/fort.53'
         WRITE(LOCNAME(IPROC)(3:6),'(I4.4)') (IPROC-1)
         LOC2(IPROC) = 105 + (IPROC-1)
         INQUIRE(FILE=LOCNAME(IPROC),EXIST=FOUND)
         IF (FOUND) THEN
           OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC))
         ELSE
           print *, "No Local unit 53 files found"
           RETURN
         ENDIF
      ENDDO
C
C   This is the global output file
C
      OPEN(UNIT=53,FILE='fort.53')
C
C   This is a DASD scratch file to reduce the memory requirement.
C
      OPEN(UNIT=8,ACCESS='DIRECT',FORM='UNFORMATTED',RECL=80)
C
      DO IPROC = 1,NPROC
         IF (IPROC.EQ.1) THEN
           READ(LOC2(IPROC),'(A80)')  NFREQMSG
           READ(NFREQMSG,*) NFREQ 
           WRITE(53,*) trim(NFREQMSG)
           DO I=1, NFREQ  
             READ(LOC2(IPROC),'(A80)') HEADER(I)
             WRITE(53,*) trim(HEADER(I))
           ENDDO
           READ(LOC2(IPROC),*) NP1
           WRITE(53,*) NNODG
         ELSE
           READ(LOC2(IPROC),'(A80)')  INLINE
           READ(INLINE,*) IDUM 
           DO I=1, NFREQ  
             READ(LOC2(IPROC),'(A80)') INLINE
           ENDDO
           READ(LOC2(IPROC),*) NP1
         ENDIF
      ENDDO
C
      DO IPROC = 1,NPROC
         DO J=1,NNODP(IPROC)
            READ(LOC2(IPROC),*) IDUM
            INDX = IMAP_NOD_LG(J,IPROC)
            IF (IMAP_NOD_GL(1,INDX).EQ.IPROC) THEN
              DO IFREQ=1, NFREQ 
                 READ(LOC2(IPROC),'(A80)') HARDAT(IFREQ)
                 IREC = IFREQ + (INDX-1)*NFREQ
                 WRITE(8,REC=IREC) HARDAT(IFREQ) 
              ENDDO
            ELSE
              DO I=1, NFREQ 
                 READ(LOC2(IPROC),'(A80)') INLINE
              ENDDO
            ENDIF
         ENDDO
      ENDDO
C
C   Close the DASD file to flush the write buffer and re-open.
C
      CLOSE(8,STATUS='KEEP')
      OPEN(UNIT=8,ACCESS='DIRECT',FORM='UNFORMATTED',RECL=80) 
C
      DO J=1, NNODG
         WRITE(53,*) J
         DO IFREQ=1, NFREQ 
            IREC = IFREQ + (J-1)*NFREQ
            READ(8,REC=IREC) HARDAT(IFREQ)    
            WRITE(53,*) trim(HARDAT(IFREQ))
         ENDDO
      ENDDO        
C
C--Close the Global and Local fort.53 Files
C
9999  CONTINUE
      CLOSE(53)
      CLOSE(8,STATUS='DELETE')
      DO IPROC = 1,NPROC
         CLOSE(LOC2(IPROC))
      ENDDO
C
c3679 FORMAT(1X,E20.10,1X,F10.7,1X,F12.8,1X,A10)
C
      RETURN
      END

      SUBROUTINE POST54()
C
C---------------------------------------------------------------------------C
C                     (  Serial Version  3/28/98  )                         C 
C  Globalize the harmonic constituent velocities at all nodes from the      C
C  local fort.54 files.                                                     C
C  This version is compatible with ADCIRC version 34.03                     C
C  Fixed routine to handled multiple frequencies up to 32--9/18/98vjp       C
C---------------------------------------------------------------------------C
C
      USE POST_GLOBAL
      IMPLICIT NONE
      INTEGER I,J,IPROC,IDUM,INDX,NP1,IFREQ,IREC
      INTEGER NFREQ
      CHARACTER*80  NFREQMSG,INLINE
      LOGICAL FOUND

      INTEGER,ALLOCATABLE  :: LOC2(:)
      CHARACTER*14,ALLOCATABLE  ::  LOCNAME(:)
      CHARACTER*80,ALLOCATABLE  ::  HEADER(:),HARDAT(:)

      ALLOCATE (LOC2(MNPROC))
      ALLOCATE (LOCNAME(MNPROC))
      ALLOCATE ( HEADER(MNHARF),HARDAT(MNHARF))

C
C--If they exist open all Local & Global Sequential Formatted fort.54 file
C
      DO IPROC = 1,NPROC
         LOCNAME(IPROC) = 'PE0000/fort.54'
         WRITE(LOCNAME(IPROC)(3:6),'(I4.4)') (IPROC-1)
         LOC2(IPROC) = 105 + (IPROC-1)
         INQUIRE(FILE=LOCNAME(IPROC),EXIST=FOUND)
         IF (FOUND) THEN
           OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC))
         ELSE
           print *, "No Local unit 54 files found"
           RETURN
         ENDIF
      ENDDO
C
C   This is the global output file
C
      OPEN(UNIT=54,FILE='fort.54')
C
C   This is a DASD scratch file to save on memory.
C
      OPEN(UNIT=8,ACCESS='DIRECT',FORM='UNFORMATTED',RECL=80)
C
      DO IPROC = 1,NPROC
         IF (IPROC.EQ.1) THEN
           READ(LOC2(IPROC),'(A80)')  NFREQMSG
           READ(NFREQMSG,*) NFREQ 
           WRITE(54,*) trim(NFREQMSG)
           DO I=1, NFREQ  
             READ(LOC2(IPROC),'(A80)') HEADER(I)
             WRITE(54,*) trim(HEADER(I))
           ENDDO
           READ(LOC2(IPROC),*) NP1
           WRITE(54,*) NNODG
         ELSE
           READ(LOC2(IPROC),'(A80)')  INLINE
           READ(INLINE,*) IDUM 
           DO I=1, NFREQ  
             READ(LOC2(IPROC),'(A80)') INLINE
           ENDDO
           READ(LOC2(IPROC),*) NP1
         ENDIF
      ENDDO
C
      DO IPROC = 1,NPROC
         DO J=1,NNODP(IPROC)
            READ(LOC2(IPROC),*) IDUM
            INDX = IMAP_NOD_LG(J,IPROC)
            IF (IMAP_NOD_GL(1,INDX).EQ.IPROC) THEN
              DO IFREQ=1, NFREQ 
                 READ(LOC2(IPROC),'(A80)') HARDAT(IFREQ)
                 IREC = IFREQ + (INDX-1)*NFREQ
                 WRITE(8,REC=IREC) HARDAT(IFREQ) 
              ENDDO
            ELSE
              DO I=1, NFREQ 
                 READ(LOC2(IPROC),'(A80)') INLINE
              ENDDO
            ENDIF
         ENDDO
      ENDDO
C
C   Close the DASD file to flush the write buffer and re-open.
C
      CLOSE(8,STATUS='KEEP')
      OPEN(UNIT=8,ACCESS='DIRECT',FORM='UNFORMATTED',RECL=80) 
C
      DO J=1, NNODG
         WRITE(54,*) J
         DO IFREQ=1, NFREQ 
            IREC = IFREQ + (J-1)*NFREQ
            READ(8,REC=IREC) HARDAT(IFREQ)    
            WRITE(54,*) trim(HARDAT(IFREQ))
         ENDDO
      ENDDO        
C
C--Close the Global and Local fort.54 Files
C
9999  CONTINUE
      CLOSE(54)
      CLOSE(8,STATUS='DELETE')
      DO IPROC = 1,NPROC
         CLOSE(LOC2(IPROC))
      ENDDO
C
c3679 FORMAT(1X,E20.10,1X,F10.7,1X,F12.8,1X,A10)
C
      RETURN
      END

      SUBROUTINE POST55()
C
C---------------------------------------------------------------------------C
C                     (  Serial Version  4/13/98  )                         C
C  Globalize the harmonic constituent comparison file from the local        C
C  fort.55 files.                                                           C
C  This version is compatible with ADCIRC version 34.03                     C
C---------------------------------------------------------------------------C
C
      USE POST_GLOBAL
      IMPLICIT NONE
      INTEGER I,J,IPROC,IDUM,INDX,NP1
      INTEGER NFREQ
      CHARACTER*114 INLINE1,INLINE2
      LOGICAL FOUND

      INTEGER,ALLOCATABLE  :: LOC2(:)
      CHARACTER*14,ALLOCATABLE  ::  LOCNAME(:)

      ALLOCATE ( LOC2(MNPROC))
      ALLOCATE (  LOCNAME(MNPROC))

C
C--If they exist open all Local & Global Sequential Formatted fort.55 file
C
      DO IPROC = 1,NPROC
         LOCNAME(IPROC) = 'PE0000/fort.55'
         WRITE(LOCNAME(IPROC)(3:6),'(I4.4)') (IPROC-1)
         LOC2(IPROC) = 105 + (IPROC-1)
         INQUIRE(FILE=LOCNAME(IPROC),EXIST=FOUND)
         IF (FOUND) THEN
           OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC))
         ELSE
           print *, "No Local unit 55 files found"
           RETURN
         ENDIF
      ENDDO
C
C   This is the global output file
C
      OPEN(UNIT=55,FILE='fort.55')
C
C   This is a DASD scratch file to save on memory.
C
      OPEN(UNIT=8,ACCESS='DIRECT',FORM='UNFORMATTED',RECL=114)       
C
C
      DO IPROC = 1,NPROC
         IF (IPROC.EQ.1) THEN
           READ(LOC2(IPROC),*) NP1
           WRITE(55,*) NNODG
         ELSE
           READ(LOC2(IPROC),*) NP1
         ENDIF
      ENDDO
C
      DO IPROC = 1,NPROC
         DO J=1,NNODP(IPROC)
            READ(LOC2(IPROC),*) IDUM
            READ(LOC2(IPROC),'(A114)') INLINE1
            INDX = IMAP_NOD_LG(J,IPROC)
            IF (IMAP_NOD_GL(1,INDX).EQ.IPROC) THEN
              WRITE(8,REC=INDX)  INLINE1
            ENDIF
         ENDDO
      ENDDO
C
C   Close the DASD scratch file to flush the write buffer and re-open.
C
      CLOSE(8,STATUS='KEEP')
      OPEN(UNIT=8,ACCESS='DIRECT',FORM='UNFORMATTED',RECL=114)       
C
      DO J=1, NNODG
          WRITE(55,*) J
          READ(8,REC=J) INLINE1
          WRITE(55,'(A114)') INLINE1
      ENDDO
      CLOSE(8,STATUS='DELETE')
C
      OPEN(UNIT=8,ACCESS='DIRECT',FORM='UNFORMATTED',RECL=114)    
      OPEN(UNIT=9,ACCESS='DIRECT',FORM='UNFORMATTED',RECL=114)    
C
      DO IPROC = 1,NPROC
         DO J=1,NNODP(IPROC)
            READ(LOC2(IPROC),*) IDUM
            READ(LOC2(IPROC),'(A114)') INLINE1
            READ(LOC2(IPROC),'(A114)') INLINE2
            INDX = IMAP_NOD_LG(J,IPROC)
            IF (IMAP_NOD_GL(1,INDX).EQ.IPROC) THEN
              WRITE(8,REC=INDX) INLINE1
              WRITE(9,REC=INDX) INLINE2
            ENDIF
         ENDDO
      ENDDO
C
C   Close both DASD scratch files to flush the write buffers and re-open.
C
      CLOSE(8,STATUS='KEEP')
      CLOSE(9,STATUS='KEEP')
      OPEN(UNIT=8,ACCESS='DIRECT',FORM='UNFORMATTED',RECL=114)    
      OPEN(UNIT=9,ACCESS='DIRECT',FORM='UNFORMATTED',RECL=114)    
C
      DO J=1, NNODG
          WRITE(55,*) J
          READ(8,REC=J) INLINE1
          READ(9,REC=J) INLINE2
          WRITE(55,'(A114)') INLINE1
          WRITE(55,'(A114)') INLINE2
      ENDDO
C
C--Close the Global and Local fort.55 Files
C
9999  CONTINUE
      CLOSE(55)
      CLOSE(8,STATUS='DELETE')
      CLOSE(9,STATUS='DELETE')
      DO IPROC = 1,NPROC
         CLOSE(LOC2(IPROC))
      ENDDO
C
      RETURN
      END

C---------------------------------------------------------------------------
C                  S U B R O U T I N E   P O S T 7 1 
C---------------------------------------------------------------------------
C     Merge the file containing atmospheric pressure data at the
C     atmospheric pressure stations (fort.71 file) from each subdomain
C     into a single, full domain fort.71 file. This version is
C     compatible with ADCIRC version 45.07. jgf 11/02/2005.
C
C     jgf46.02 Fixed fort.80 file to contain the station mapping, fixed
C     memory allocation bug in post_global.F (had MNSTAV instead of
C     MNSTAM).
C
C---------------------------------------------------------------------------
      SUBROUTINE POST71()
      USE POST_GLOBAL
      IMPLICIT NONE
c     counters
      INTEGER I,J,K,L ! loop counters 
      INTEGER IPROC   ! subdomain counter
      INTEGER OREC    ! record counter in binary output file
      INTEGER,ALLOCATABLE :: IREC(:) ! input record counters for binary files
c     containers
      INTEGER IDUM                         ! to read number of nodes in binary
      CHARACTER*80  OUTMSG                 ! text line of full domain data 
      CHARACTER*85  INLINE                 ! hold RUNDES, RUNID, AGRID in text
      CHARACTER*4 RDES4(8),RID4(6),AID4(6) ! hold RUNDES, RUNID, AGRID in bin.
      CHARACTER*8 RDES8(8),RID8(6),AID8(6) ! hold RUNDES, RUNID, AGRID in bin.
c     file stuff 
      LOGICAL FOUND                        ! TRUE if the file exists
      INTEGER,ALLOCATABLE::LOC2(:)          ! unit numbers of subdomains
      CHARACTER*14,ALLOCATABLE::LOCNAME(:) ! subdomain dir./file names
c     actual data 
      INTEGER NTRSPM  ! number of data sets to be written to fort.71
      INTEGER NUMSTNS ! number of subdomain meteorlogical recording stations
      INTEGER NSTEMP  ! time step interval at which fort.71 is written
      INTEGER ITEMPM  ! record type (1=elev.,2=vel.,3=3D vel.)
      REAL(SZ) DTM    ! total time: DT*NSPOOLM = step size * no. time steps
      INTEGER ITSM    ! model time step number
      REAL(8)TIMEOUTSM! model time
      REAL(SZ),ALLOCATABLE::AtmPresBIN(:)    ! full domain atm. pres. val.(bin)
      REAL(SZ),ALLOCATABLE::AtmPresBINP(:)   ! subdomain atm. pres. val. (bin)
      INTEGER,ALLOCATABLE::StaNum(:)         ! recording station number
      CHARACTER(len=80),ALLOCATABLE::AtmPresASC(:) ! atmospheric pressure data
C
      ALLOCATE ( LOC2(MNPROC),IREC(MNPROC) )
      ALLOCATE ( AtmPresBIN(MNSTAM),AtmPresBINP(MNSTAM))
      ALLOCATE ( AtmPresASC(MNSTAM))
      ALLOCATE ( LOCNAME(MNPROC))
C
C     Determine whether Unit 71 is Sequential Formatted or Direct Access Binary
C
      SELECT CASE(ABS(NOUTM))
C     -------      
      CASE(1)                   ! text (sequential formatted)
C     -------
C     Open All Local Sequential Formatted fort.71 files
         DO IPROC = 1,NPROC
            LOCNAME(IPROC) = 'PE0000/fort.71'
            WRITE(LOCNAME(IPROC)(3:6),'(I4.4)') (IPROC-1)
            LOC2(IPROC) = 105 + (IPROC-1)
            INQUIRE(FILE=LOCNAME(IPROC),EXIST=FOUND)
            IF (FOUND) THEN
               OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC))
            ELSE
               print *, "ERROR: No subdomain fort.71 files found."
               RETURN           ! *** EARLY RETURN ***
            ENDIF
            READ(LOC2(IPROC),'(A85)') INLINE
            READ(LOC2(IPROC),3645) NTRSPM,NUMSTNS,DTM,NSTEMP,ITEMPM
            IF (NUMSTNS.NE.NSTAMP(IPROC)) THEN
               STOP 'ERROR: Inconsistency in no. of atm. pres. sta.'
            ENDIF
         ENDDO
C     
C     Open Global Sequential Formatted fort.71 file
         OPEN(UNIT=71,FILE='fort.71')
C     
         WRITE(71,'(A85)') INLINE
         WRITE(71,3645) NTRSPM,NUMSTNS,DTM,NSPOOLM,ITEMPM
C     
         DO                     ! loop infinitely over data sets
C     
            DO IPROC=1,NPROC
               READ(LOC2(IPROC),2120,END=9999) TIMEOUTSM,ITSM
               IF (NSTAMP(IPROC).GT.0) THEN
                  READ(LOC2(IPROC),80,END=9999) 
     &                 (AtmPresASC(IMAP_STAM_LG(K,IPROC)),
     &                 K=1,NSTAMP(IPROC))
               ENDIF
            ENDDO            
C     
            WRITE(71,2120) TIMEOUTSM,ITSM
            DO I=1, NSTAM
               CALL NEWINDEX(AtmPresASC(I),OUTMSG,I)
               WRITE(71,*) TRIM(OUTMSG) !jgf46.00 TRIM to cut down file size 
            ENDDO
C     
         ENDDO                 
C     -------
      CASE(2)                   ! binary
C     -------
C     Open All Local and the Global Direct Access Binary fort.71 files
         DO IPROC = 1,NPROC
            LOCNAME(IPROC) = 'PE0000/fort.71'
            WRITE(LOCNAME(IPROC)(3:6),'(I4.4)') (IPROC-1)
            LOC2(IPROC) = 105 + (IPROC-1)
            INQUIRE(FILE=LOCNAME(IPROC),EXIST=FOUND)
            IF (FOUND) THEN
               OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC),ACCESS='DIRECT',
     &              RECL=NBYTE)
            ELSE
               print *, "ERROR: No subdomain fort.71 files found."
               RETURN
            ENDIF
         ENDDO
         OPEN(71,FILE='fort.71',ACCESS='DIRECT',RECL=NBYTE)
C     
C     
C--   Read RUNDES RUNID and AGRID from each subdomain file
C     
         DO IPROC = 1,NPROC
            IREC(IPROC) = 0
            IF (NBYTE.EQ.4) THEN
               DO I=1,8
                  READ(LOC2(IPROC),REC=IREC(IPROC)+I) RDES4(I)
               ENDDO
               IREC(IPROC)=IREC(IPROC)+8
               DO I=1,6
                  READ(LOC2(IPROC),REC=IREC(IPROC)+I) RID4(I)
               ENDDO
               IREC(IPROC)=IREC(IPROC)+6
               DO I=1,6
                  READ(LOC2(IPROC),REC=IREC(IPROC)+I) AID4(I)
               ENDDO
               IREC(IPROC)=IREC(IPROC)+6
            ENDIF
            IF (NBYTE.EQ.8) THEN
               DO I=1,4
                  READ(LOC2(IPROC),REC=IREC(IPROC)+I) RDES8(I)
               ENDDO
               IREC(IPROC)=IREC(IPROC)+4
               DO I=1,3
                  READ(LOC2(IPROC),REC=IREC(IPROC)+I) RID8(I)
               ENDDO
               IREC(IPROC)=IREC(IPROC)+3
               DO I=1,3
                  READ(LOC2(IPROC),REC=IREC(IPROC)+I) AID8(I)
               ENDDO
               IREC(IPROC)=IREC(IPROC)+3
            ENDIF
         ENDDO
C     
C--   Write RUNDES RUNID and AGRID to Global File
C     
         OREC = 0
         IF (NBYTE.EQ.4) THEN
            DO I=1,8
               WRITE(71,REC=OREC+I) RDES4(I)
            ENDDO
            OREC=OREC+8
            DO I=1,6
               WRITE(71,REC=OREC+I) RID4(I)
            ENDDO
            OREC=OREC+6
            DO I=1,6
               WRITE(71,REC=OREC+I) AID4(I)
            ENDDO
            OREC=OREC+6
         ENDIF
         IF (NBYTE.EQ.8) THEN
            DO I=1,4
               WRITE(71,REC=OREC+I) RDES8(I)
            ENDDO
            OREC=OREC+4
            DO I=1,3
               WRITE(71,REC=OREC+I) RID8(I)
            ENDDO
            OREC=OREC+3
            DO I=1,3
               WRITE(71,REC=OREC+I) AID8(I)
            ENDDO
            OREC=OREC+3
         ENDIF
C     
C--   Read NTRSPM, NSTAM, DT*NSPOOLM from each subdomain file
C     and then close files to flush file buffers
C     
         DO IPROC = 1,NPROC
            READ(LOC2(IPROC),REC=IREC(IPROC)+1) NTRSPM
            READ(LOC2(IPROC),REC=IREC(IPROC)+2) NUMSTNS        
            READ(LOC2(IPROC),REC=IREC(IPROC)+3) DTM
            READ(LOC2(IPROC),REC=IREC(IPROC)+4) NSPOOLM 
            READ(LOC2(IPROC),REC=IREC(IPROC)+5) ITEMPM
            IREC(IPROC) = IREC(IPROC)+5
            CLOSE(LOC2(IPROC))   ! Flush the Write Buffer
            IF (NSTAM.NE.NSTAMP(IPROC)) THEN
               STOP 'ERROR: Inconsistency in no. of atm. pres. sta.'
            ENDIF
            OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC),ACCESS='DIRECT',
     &           RECL=NBYTE)
         ENDDO
C     
C--   Write same info to full-domain file and close it to flush buffer also
C     
         WRITE(71,REC=OREC+1) NTRSPM
         WRITE(71,REC=OREC+2) NUMSTNS
         WRITE(71,REC=OREC+3) DTM
         WRITE(71,REC=OREC+4) NSPOOLM
         WRITE(71,REC=OREC+5) ITEMPM
         OREC = OREC+5
         CLOSE(71)              ! Flush the Write Buffer
         OPEN(71,FILE='fort.71',ACCESS='DIRECT',RECL=NBYTE)
C     
         DO J=1,NTRSPM
C     
            DO IPROC=1,NPROC
               READ(LOC2(IPROC),REC=IREC(IPROC)+1) TIMEOUTSM
               READ(LOC2(IPROC),REC=IREC(IPROC)+2) ITSM
               IREC(IPROC) = IREC(IPROC) + 2
               IF (NSTAMP(IPROC).GT.0) THEN
                  DO K=1, NSTAMP(IPROC)
                     READ(LOC2(IPROC),REC=IREC(IPROC)+K) AtmPresBINP(K)
                  ENDDO
                  DO K=1,NSTAMP(IPROC)
                     AtmPresBIN(IMAP_STAM_LG(K,IPROC)) = AtmPresBINP(K)
                  ENDDO
                  IREC(IPROC) = IREC(IPROC) + NSTAMP(IPROC)
               ENDIF
            ENDDO            
C     
            WRITE(71,REC=OREC+1) TIMEOUTSM
            WRITE(71,REC=OREC+2) ITSM
            OREC = OREC + 2
            DO K=1, NSTAM
               WRITE(71,REC=OREC+K) AtmPresBIN(K)
            ENDDO
            OREC = OREC + NSTAM
C     
         ENDDO                 
C
      CASE DEFAULT
C     
         WRITE(*,2125) NOUTM
C
      END SELECT
C     
C     Close the full domain and subdomain fort.71 Files
C
 9999 CONTINUE
      CLOSE(71)
      DO IPROC = 1,NPROC
         CLOSE(LOC2(IPROC))
      ENDDO
C     
 80   FORMAT(A40)
 2120 FORMAT(2X,E20.10,5X,I10)
 2125 FORMAT(2X,'ERROR: Incorrect value of NOUTM; NOUTM=',I3) 
 3645 FORMAT(1X,I10,1X,I10,1X,E15.7,1X,I5,1X,I5)
C     
      RETURN
C---------------------------------------------------------------------------
      END SUBROUTINE POST71
C---------------------------------------------------------------------------




C---------------------------------------------------------------------------
C                  S U B R O U T I N E   P O S T 7 2 
C---------------------------------------------------------------------------
C     Merge the file containing wind velocity data at the meteorological
C     recording stations (fort.72 file) from each subdomain into a
C     single, full domain fort.72 file. This version is compatible with
C     ADCIRC version 45.07. jgf 11/02/2005.
C
C     jgf46.02 Fixed memory allocation bug (had MNSTAV instead of MNSTAM).
C
C---------------------------------------------------------------------------
      SUBROUTINE POST72()
      USE POST_GLOBAL
      IMPLICIT NONE
c     counters
      INTEGER I,J,K,L             ! loop counters 
      INTEGER IPROC               ! subdomain counter
      INTEGER OREC                ! record counter in binary output file
      INTEGER,ALLOCATABLE::IREC(:)! input record counters for binary files
c     containers
      INTEGER IDUM                ! to read number of nodes in binary
      CHARACTER*80 OUTMSG         ! text line of full domain data 
      CHARACTER*85  INLINE        ! hold RUNDES, RUNID, AGRID in text
      CHARACTER*4 RDES4(8),RID4(6),AID4(6) ! hold RUNDES, RUNID, AGRID in bin.
      CHARACTER*8 RDES8(8),RID8(6),AID8(6) ! hold RUNDES, RUNID, AGRID in bin.
c     file stuff 
      LOGICAL FOUND               ! TRUE if the file exists
      INTEGER,ALLOCATABLE::LOC2(:) ! unit numbers of subdomains
      CHARACTER*14,ALLOCATABLE::LOCNAME(:) ! subdomain dir./file names
c     actual data 
      REAL(SZ) DTM    ! total time: DT*NSPOOLM = step size * no. time steps
      REAL(SZ)TIMEOUTSM! model time
      INTEGER NTRSPM  ! number of data sets to be written to fort.72
      INTEGER NSTEMP  ! time step interval at which fort.72 is written
      INTEGER ITSM    ! model time step number
      INTEGER ITEMPM  ! record type (1=elev.,2=vel.,3=3D vel.)
      INTEGER NUMSTNS ! number of subdomain meteorlogical recording stations
      REAL(SZ),ALLOCATABLE::UWindVelBIN(:)  ! fulldomain u wind velocity values
      REAL(SZ),ALLOCATABLE::VWindVelBIN(:)  ! fulldomain v wind velocity values
      REAL(SZ),ALLOCATABLE::UWindVelBINP(:) ! subdomain u wind velocity values
      REAL(SZ),ALLOCATABLE::VWindVelBINP(:) ! subdomain v wind velocity values
      CHARACTER*80,ALLOCATABLE::WindASC(:)  ! full domain wind vel.. val.
      CHARACTER*80,ALLOCATABLE::WindASCP(:) ! subdomain wind vel. val.
C
      ALLOCATE (LOC2(MNPROC),IREC(MNPROC))
      ALLOCATE (UWindVelBIN(MNSTAM),VWindVelBIN(MNSTAM))
      ALLOCATE (UWindVelBINP(MNSTAM),VWindVelBINP(MNSTAM))
      ALLOCATE (WindASC(MNSTAM))
      ALLOCATE (LOCNAME(MNPROC))
C
C--Determine whether Unit 72 is Sequential Formatted or Direct Access Binary
C
      SELECT CASE(ABS(NOUTM))

      CASE(1)
         DO IPROC = 1,NPROC
            LOCNAME(IPROC) = 'PE0000/fort.72'
            WRITE(LOCNAME(IPROC)(3:6),'(I4.4)') (IPROC-1)
            LOC2(IPROC) = 105 + (IPROC-1)
            INQUIRE(FILE=LOCNAME(IPROC),EXIST=FOUND)
            IF (FOUND) THEN
               OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC))
            ELSE
               print *, "ERROR: No subdomain fort.72 files found."
               RETURN
            ENDIF
            READ(LOC2(IPROC),'(A85)')  INLINE
            READ(LOC2(IPROC),3645) NTRSPM,NUMSTNS,DTM,NSTEMP,ITEMPM
            IF (NUMSTNS.NE.NSTAMP(IPROC)) THEN
               STOP 'ERROR: Inconsistency in number of vel. stations.'
            ENDIF
         ENDDO
C     
C--   Open full domain fort.72 file
C     
         OPEN(UNIT=72,FILE='fort.72')
C     
         WRITE(72,'(A85)') INLINE
         WRITE(72,3645) NTRSPM,MNSTAM,DTM,NSTEMP,ITEMPM
C     
C     jgf46.00 Made the loop over data sets an infinite loop because the
C     number of data sets at the top of the file will be wrong if the
C     run was hot started.
C     
         DO                     ! loop infinitely over data sets
C     
            DO IPROC=1,NPROC
               READ(LOC2(IPROC),2120,END=9999) TIMEOUTSM,ITSM
               IF (NSTAMP(IPROC).GT.0) THEN
                  READ(LOC2(IPROC),80,END=9999) 
     &                 (WindASC(IMAP_STAM_LG(K,IPROC)),
     &                 K=1,NSTAMP(IPROC))
               ENDIF
            ENDDO            
C     
            WRITE(72,2120) TIMEOUTSM,ITSM
            DO I=1, NSTAM
               CALL NEWINDEX(WindASC(I),OUTMSG,I)
               WRITE(72,*) TRIM(OUTMSG) !jgf46.00 TRIM to cut down file size 
            ENDDO
C     
         ENDDO                 

      CASE(2)
C     
C--   Open All subdomain and the full domain Direct Access Binary fort.72 file
C     
         DO IPROC = 1,NPROC
            LOCNAME(IPROC) = 'PE0000/fort.72'
            WRITE(LOCNAME(IPROC)(3:6),'(I4.4)') (IPROC-1)
            LOC2(IPROC) = 105 + (IPROC-1)
            INQUIRE(FILE=LOCNAME(IPROC),EXIST=FOUND)
            IF (FOUND) THEN
               OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC),ACCESS='DIRECT',
     &              RECL=NBYTE)
            ELSE
               print *, "ERROR: No subdomain fort.72 files found."
               RETURN
            ENDIF
         ENDDO
         OPEN(72,FILE='fort.72',ACCESS='DIRECT',RECL=NBYTE)
C     
C--   Read RUNDES RUNID and AGRID from each subdomain file
C     
         DO IPROC = 1,NPROC
            IREC(IPROC) = 0
            IF (NBYTE.EQ.4) THEN
               DO I=1,8
                  READ(LOC2(IPROC),REC=IREC(IPROC)+I) RDES4(I)
               ENDDO
               IREC(IPROC)=IREC(IPROC)+8
               DO I=1,6
                  READ(LOC2(IPROC),REC=IREC(IPROC)+I) RID4(I)
               ENDDO
               IREC(IPROC)=IREC(IPROC)+6
               DO I=1,6
                  READ(LOC2(IPROC),REC=IREC(IPROC)+I) AID4(I)
               ENDDO
               IREC(IPROC)=IREC(IPROC)+6
            ENDIF
            IF (NBYTE.EQ.8) THEN
               DO I=1,4
                  READ(LOC2(IPROC),REC=IREC(IPROC)+I) RDES8(I)
               ENDDO
               IREC(IPROC)=IREC(IPROC)+4
               DO I=1,3
                  READ(LOC2(IPROC),REC=IREC(IPROC)+I) RID8(I)
               ENDDO
               IREC(IPROC)=IREC(IPROC)+3
               DO I=1,3
                  READ(LOC2(IPROC),REC=IREC(IPROC)+I) AID8(I)
               ENDDO
               IREC(IPROC)=IREC(IPROC)+3
            ENDIF
         ENDDO
C     
C--   Write RUNDES RUNID and AGRID to full domain File
C     
         OREC = 0
         IF (NBYTE.EQ.4) THEN
            DO I=1,8
               WRITE(72,REC=OREC+I) RDES4(I)
            ENDDO
            OREC=OREC+8
            DO I=1,6
               WRITE(72,REC=OREC+I) RID4(I)
            ENDDO
            OREC=OREC+6
            DO I=1,6
               WRITE(72,REC=OREC+I) AID4(I)
            ENDDO
            OREC=OREC+6
         ENDIF
         IF (NBYTE.EQ.8) THEN
            DO I=1,4
               WRITE(72,REC=OREC+I) RDES8(I)
            ENDDO
            OREC=OREC+4
            DO I=1,3
               WRITE(72,REC=OREC+I) RID8(I)
            ENDDO
            OREC=OREC+3
            DO I=1,3
               WRITE(72,REC=OREC+I) AID8(I)
            ENDDO
            OREC=OREC+3
         ENDIF
C     
C--   Read NTRSPM, NSTAM, DT*NSPOOLM from each subdomain file
C     and then close files to flush file buffers
C     
         DO IPROC = 1,NPROC
            READ(LOC2(IPROC),REC=IREC(IPROC)+1) NTRSPM
            READ(LOC2(IPROC),REC=IREC(IPROC)+2) NUMSTNS
            READ(LOC2(IPROC),REC=IREC(IPROC)+3) DTM
            READ(LOC2(IPROC),REC=IREC(IPROC)+4) NSTEMP
            READ(LOC2(IPROC),REC=IREC(IPROC)+5) ITEMPM
            IREC(IPROC) = IREC(IPROC)+5
            CLOSE(LOC2(IPROC))   ! Flush the Write Buffer
            IF (NUMSTNS.NE.NSTAMP(IPROC)) THEN
               STOP 'ERROR: Inconsistency in number of vel. stations.'
            ENDIF
            OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC),
     &           ACCESS='DIRECT',RECL=NBYTE)
         ENDDO
C     
C--   Write same info to full domain file and close it to flush buffer also
C     
         WRITE(72,REC=OREC+1) NTRSPM
         WRITE(72,REC=OREC+2) NSTAM
         WRITE(72,REC=OREC+3) DTM
         WRITE(72,REC=OREC+4) NSPOOLM
         WRITE(72,REC=OREC+5) ITEMPM
         OREC = OREC+5
         CLOSE(72)              ! Flush the Write Buffer
         OPEN(72,FILE='fort.72',ACCESS='DIRECT',RECL=NBYTE)
C     
         DO J=1,NTRSPM
C     
            DO IPROC=1,NPROC
               READ(LOC2(IPROC),REC=IREC(IPROC)+1) TIMEOUTSM
               READ(LOC2(IPROC),REC=IREC(IPROC)+2) ITSM
               IREC(IPROC) = IREC(IPROC) + 2
               IF (NSTAMP(IPROC).GT.0) THEN
                  DO K=1, NSTAMP(IPROC)
                     READ(LOC2(IPROC),
     &                    REC=IREC(IPROC)+2*K-1) UWindVelBINP(K)
                     READ(LOC2(IPROC),
     &                    REC=IREC(IPROC)+2*K) VWindVelBINP(K)
                  ENDDO
                  DO K=1,NSTAMP(IPROC)
                     UWindVelBIN(IMAP_STAV_LG(K,IPROC)) 
     &                    = UWindVelBINP(K)
                     VWindVelBIN(IMAP_STAV_LG(K,IPROC)) 
     &                    = VWindVelBINP(K)
                  ENDDO
                  IREC(IPROC) = IREC(IPROC) + 2*NSTAMP(IPROC)
               ENDIF
            ENDDO
C     
            WRITE(72,REC=OREC+1) TIMEOUTSM
            WRITE(72,REC=OREC+2) ITSM
            OREC = OREC + 2
            DO K=1, NSTAM
               WRITE(72,REC=OREC+2*K-1) UWindVelBIN(K)
               WRITE(72,REC=OREC+2*K) VWindVelBIN(K)
            ENDDO
            OREC = OREC + 2*NSTAM
C     
         ENDDO
C     
         CASE DEFAULT
C
            WRITE(*,2125) NOUTM
C
      END SELECT
C
 9999 CONTINUE
C
C--Close the full domain and subdomain fort.72 Files
C
      CLOSE(72)
      DO IPROC = 1,NPROC
         CLOSE(LOC2(IPROC))
      ENDDO
C     
 2120 FORMAT(2X,E20.10,5X,I10)
 80   FORMAT(A80)
c     2454 FORMAT(2X,I8,2X,E15.8,2X,E15.8)
 2125 FORMAT(2X,'ERROR: Incorrect value of NOUTM; NOUTM=',I3) 
 3645 FORMAT(1X,I10,1X,I10,1X,E15.7,1X,I5,1X,I5)
C     
      RETURN
C---------------------------------------------------------------------------
      END SUBROUTINE POST72
C---------------------------------------------------------------------------





C---------------------------------------------------------------------------
C                  S U B R O U T I N E   P O S T 7 3 
C---------------------------------------------------------------------------
C     Merge the atmospheric pressure data (fort.73 files) from each
C     subdomain into a single, full domain fort.73 file. This version is
C     compatible with ADCIRC version 45.07. jgf 11/02/2005.
C---------------------------------------------------------------------------
      SUBROUTINE POST73()
      USE POST_GLOBAL
      IMPLICIT NONE
C     counters
      INTEGER I,J,K,L ! loop counters 
      INTEGER IPROC   ! subdomain counter
      INTEGER OREC    ! record counter in binary output file
      INTEGER, ALLOCATABLE::IREC(:) ! input record counters for binary files
      INTEGER INDX   ! full domain node number
C     containers
      CHARACTER*80 OUTMSG                  ! text line of full domain data 
      CHARACTER*85 INLINE                  ! hold RUNDES, RUNID, AGRID in text
      CHARACTER*4 RDES4(8),RID4(6),AID4(6) ! hold RUNDES, RUNID, AGRID in bin.
      CHARACTER*8 RDES8(8),RID8(6),AID8(6) ! hold RUNDES, RUNID, AGRID in bin.
      INTEGER IDUM                         ! to read number of nodes in binary
C     file stuff
      INTEGER, ALLOCATABLE::LOC2(:)         ! unit numbers of subdomains
      CHARACTER*14, ALLOCATABLE::LOCNAME(:)! subdomain dir./file names
      LOGICAL FOUND                        ! TRUE if the file exists
C     actual data
      INTEGER NDSETSGW! number of data sets to be written to fort.73
      INTEGER NP      ! number of nodes in the subdomain grid
      REAL(SZ) DTGW   ! total time: DT*NSPOOLGW = step size * no. time steps
      INTEGER NSTEMP  ! time step interval at which fort.73 is written
      INTEGER ITEMPGW ! record type (1=elev.,2=vel.,3=3D vel.)
C
      REAL(8)TIMEOUTGW! model time
      INTEGER ITGW    ! model time step
C
      CHARACTER*80,ALLOCATABLE::AtmPresSC(:) !full domain atm. pressure values
      CHARACTER*80,ALLOCATABLE::AtmPresSCP(:)!subdomain atm. pressure values
      REAL(SZ),ALLOCATABLE::AtmPresBIN(:)    ! full domain atm. pres. val.(bin)
      REAL(SZ),ALLOCATABLE::AtmPresBINP(:)   ! subdomain atm. pres. val. (bin)
C
      ALLOCATE (LOC2(MNPROC),IREC(MNPROC))
      ALLOCATE (AtmPresBIN(MNP),AtmPresBINP(MNPP))
      ALLOCATE (AtmPresSC(MNP),AtmPresSCP(MNPP))
      ALLOCATE ( LOCNAME(MNPROC))
C     
C--   Determine whether Unit 73 is Sequential Formatted or Direct Access Binary
C     
      IF (ABS(NOUTGW).EQ.1) THEN
         GO TO 1000 ! file is text
      ELSE
         GO TO 2000 ! file is binary(=2) or nothing (=0)
      ENDIF
C
C
 1000 CONTINUE !here begins the post-processing if the files are ASCII text
C
      DO IPROC = 1,NPROC
C     Open ASCII text fort.73 files in each subdirectory
         LOCNAME(IPROC) = 'PE0000/fort.73'
         WRITE(LOCNAME(IPROC)(3:6),'(I4.4)') (IPROC-1)
         LOC2(IPROC) = 105 + (IPROC-1)
         INQUIRE(FILE=LOCNAME(IPROC),EXIST=FOUND)
         IF (FOUND) THEN
            OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC))
         ELSE
            print *, "ERROR: No local fort.73 files found."
            RETURN !**** EARLY RETURN ****
         ENDIF
C     Read header information from each fort.73 file
         READ(LOC2(IPROC),'(A85)') INLINE ! read RUNDES, RUNID, AGRID
         READ(LOC2(IPROC),3645)  NDSETSGW,NP,DTGW,NSTEMP,ITEMPGW
         IF (NP.NE.NNODP(IPROC)) THEN
            print *, "NP = ",NP," NNODP(",IPROC,") = ",NNODP(IPROC)
            STOP 'ERROR: Inconsistency in expected subdomain nodes.'
         ENDIF
      ENDDO
C     Create fort.73 file for full domain
      OPEN(UNIT=73,FILE='fort.73')
      WRITE(73,'(A85)') INLINE ! write RUNDES, RUNID, AGRID
      WRITE(73,3645) NDSETSGW,NNODG,DTGW,NSPOOLGW,ITEMPGW 
C
C     jgf46.00 Made the loop over data sets an infinite loop because the
C     number of data sets at the top of the file will be wrong if the
C     run was hot started.
C
      DO ! loop infinitely over data sets
C
         DO IPROC=1,NPROC  ! loop over subdomains
            READ(LOC2(IPROC),FMT=2120,END=9999) TIMEOUTGW,ITGW
            IF (IPROC.EQ.1) WRITE(73,2120) TIMEOUTGW,ITGW    
            DO K = 1,NNODP(IPROC) ! loop over subdomain nodes
               READ(LOC2(IPROC),80,END=9999) AtmPresSCP(K)
               INDX = IMAP_NOD_LG(K,IPROC) ! sub --> full (lookup table)
               IF (IMAP_NOD_GL(1,INDX).EQ.IPROC) THEN
                  AtmPresSC(INDX) = AtmPresSCP(K)
               ENDIF
            ENDDO ! next subdomain node  
         ENDDO ! next subdomain
C
         DO I = 1,NNODG   ! loop over full domain nodes
            CALL NEWINDEX(AtmPresSC(I),OUTMSG,I)
            WRITE(73,*) TRIM(OUTMSG) !jgf46.00 TRIM to cut down file size 
         ENDDO 
      ENDDO ! next output time
      GO TO 9999 ! jump to end of subroutine where fort.73 file is closed
C
C     
 2000 CONTINUE !here begins the post-processing if files are binary
C--   Open all subdomain and the full domain Direct Access Binary fort.73 files
      DO IPROC = 1,NPROC
         LOCNAME(IPROC) = 'PE0000/fort.73'
         WRITE(LOCNAME(IPROC)(3:6),'(I4.4)') (IPROC-1)
         LOC2(IPROC) = 105 + (IPROC-1)
         INQUIRE(FILE=LOCNAME(IPROC),EXIST=FOUND)
         IF (FOUND) THEN
            OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC),ACCESS='DIRECT',
     &           RECL=NBYTE)
         ELSE
            print *, "ERROR: No local fort.73 files found."
            RETURN
         ENDIF
      ENDDO
      OPEN(73,FILE='fort.73',ACCESS='DIRECT',RECL=NBYTE)
C--   Read RUNDES RUNID and AGRID from each subdomain file
      DO IPROC = 1,NPROC
         IREC(IPROC) = 0
         IF (NBYTE.EQ.4) THEN
            DO I=1,8
               READ(LOC2(IPROC),REC=IREC(IPROC)+I) RDES4(I)
            ENDDO
            IREC(IPROC)=IREC(IPROC)+8
            DO I=1,6
               READ(LOC2(IPROC),REC=IREC(IPROC)+I) RID4(I)
            ENDDO
            IREC(IPROC)=IREC(IPROC)+6
            DO I=1,6
               READ(LOC2(IPROC),REC=IREC(IPROC)+I) AID4(I)
            ENDDO
            IREC(IPROC)=IREC(IPROC)+6
         ENDIF
         IF (NBYTE.EQ.8) THEN
            DO I=1,4
               READ(LOC2(IPROC),REC=IREC(IPROC)+I) RDES8(I)
            ENDDO
            IREC(IPROC)=IREC(IPROC)+4
            DO I=1,3
               READ(LOC2(IPROC),REC=IREC(IPROC)+I) RID8(I)
            ENDDO
            IREC(IPROC)=IREC(IPROC)+3
            DO I=1,3
               READ(LOC2(IPROC),REC=IREC(IPROC)+I) AID8(I)
            ENDDO
            IREC(IPROC)=IREC(IPROC)+3
         ENDIF
      ENDDO
C--   Write RUNDES RUNID and AGRID to full domain file
      OREC = 0
      IF (NBYTE.EQ.4) THEN
         DO I=1,8
            WRITE(73,REC=OREC+I) RDES4(I)
         ENDDO
         OREC=OREC+8
         DO I=1,6
            WRITE(73,REC=OREC+I) RID4(I)
         ENDDO
         OREC=OREC+6
         DO I=1,6
            WRITE(73,REC=OREC+I) AID4(I)
         ENDDO
         OREC=OREC+6
      ENDIF
      IF (NBYTE.EQ.8) THEN
         DO I=1,4
            WRITE(73,REC=OREC+I) RDES8(I)
         ENDDO
         OREC=OREC+4
         DO I=1,3
            WRITE(73,REC=OREC+I) RID8(I)
         ENDDO
         OREC=OREC+3
         DO I=1,3
            WRITE(73,REC=OREC+I) AID8(I)
         ENDDO
         OREC=OREC+3
      ENDIF
C     
C--   Read NDSETSGW, NNODP, DTGW, NSPOOLGW, ITEMPGW from each subdomain file
C     and then close files to flush file buffers
C     
      DO IPROC = 1,NPROC
         READ(LOC2(IPROC),REC=IREC(IPROC)+1) NDSETSGW
         READ(LOC2(IPROC),REC=IREC(IPROC)+2) IDUM   
         READ(LOC2(IPROC),REC=IREC(IPROC)+3) DTGW
         READ(LOC2(IPROC),REC=IREC(IPROC)+4) NSTEMP
         READ(LOC2(IPROC),REC=IREC(IPROC)+5) ITEMPGW
         IREC(IPROC) = IREC(IPROC)+5
         CLOSE(LOC2(IPROC))      ! flush the write buffer
         OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC),ACCESS='DIRECT',RECL=NBYTE)
      ENDDO
C     
C--   Write same info to full domain file and close it to flush buffer also
C     
      WRITE(73,REC=OREC+1) NDSETSGW
      WRITE(73,REC=OREC+2) NNODG
      WRITE(73,REC=OREC+3) DTGW
      WRITE(73,REC=OREC+4) NSPOOLGW
      WRITE(73,REC=OREC+5) ITEMPGW
      OREC = OREC+5
      CLOSE(73)                 ! flush the write buffer
      OPEN(73,FILE='fort.73',ACCESS='DIRECT',RECL=NBYTE)
C     
      DO J=1,NDSETSGW
C     
         DO IPROC=1,NPROC
            READ(LOC2(IPROC),REC=IREC(IPROC)+1) TIMEOUTGW
            READ(LOC2(IPROC),REC=IREC(IPROC)+2) ITGW
            IREC(IPROC) = IREC(IPROC) + 2
            DO K=1, NNODP(IPROC)
               READ(LOC2(IPROC),REC=IREC(IPROC)+K) AtmPresBINP(K)
            ENDDO
            DO K=1,NNODP(IPROC)
               INDX = IMAP_NOD_LG(K,IPROC)
               IF (IMAP_NOD_GL(1,INDX).EQ.IPROC) THEN
                  AtmPresBIN(INDX) = AtmPresBINP(K)
               ENDIF
            ENDDO
            IREC(IPROC) = IREC(IPROC) + NNODP(IPROC)
         ENDDO
C     
         WRITE(73,REC=OREC+1) TIMEOUTGW
         WRITE(73,REC=OREC+2) ITGW
         OREC = OREC + 2
         DO K=1, NNODG
            WRITE(73,REC=OREC+K) AtmPresBIN(K)
         ENDDO
         OREC = OREC + NNODG
C     
      ENDDO
C     
 9999 CONTINUE
C     
C--   Close the full domain and subdomain fort.73 Files
C     
      CLOSE(73)
      DO IPROC = 1,NPROC
         CLOSE(LOC2(IPROC))
      ENDDO
C     
 80   FORMAT(A80)
 2120 FORMAT(2X,E20.10,5X,I10)
 2453 FORMAT(2X,I8,2X,E15.8)
 3645 FORMAT(1X,I10,1X,I10,1X,E15.7,1X,I5,1X,I5)
C     
      RETURN
      END
C---------------------------------------------------------------------------
C     End of subroutine Post73 
C---------------------------------------------------------------------------



      SUBROUTINE POST74()
C
C---------------------------------------------------------------------------C
C                     (  Serial Version  3/28/98  )                         C 
C  Globalize the Wind Stress Data at all nodes from the local fort.74 files.C
C  This version is compatible with ADCIRC version 34.03                     C
C---------------------------------------------------------------------------C
C
      USE POST_GLOBAL
      IMPLICIT NONE
      INTEGER I,J,K,L,IPROC
      INTEGER NDSETGV,NSTEMP,NP,ITV,ITEMPV,INDX
      INTEGER OREC
      REAL(SZ) DTV
      REAL(8) TIMEOUTV
      CHARACTER*80  OUTMSG
      CHARACTER*85 INLINE
      CHARACTER*4 RDES4(8),RID4(6),AID4(6)
      CHARACTER*8 RDES8(8),RID8(6),AID8(6)
      LOGICAL FOUND

      INTEGER,ALLOCATABLE  :: LOC2(:),IREC(:)
      REAL(SZ),ALLOCATABLE  ::  UUBIN(:),VVBIN(:)
      REAL(SZ),ALLOCATABLE  ::  UUBINP(:),VVBINP(:)
      CHARACTER*80,ALLOCATABLE  ::  UU(:),UUP(:)
      CHARACTER*14,ALLOCATABLE  ::  LOCNAME(:)

      ALLOCATE ( LOC2(MNPROC),IREC(MNPROC))
      ALLOCATE (  UUBIN(MNP),VVBIN(MNP))
      ALLOCATE (  UUBINP(MNPP),VVBINP(MNPP))
      ALLOCATE ( UU(MNP),UUP(MNPP))
      ALLOCATE ( LOCNAME(MNPROC))

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
C--Open Global fort.74 file
C
      DO IPROC = 1,NPROC
         LOCNAME(IPROC) = 'PE0000/fort.74'
         WRITE(LOCNAME(IPROC)(3:6),'(I4.4)') (IPROC-1)
         LOC2(IPROC) = 105 + (IPROC-1)
         INQUIRE(FILE=LOCNAME(IPROC),EXIST=FOUND)
         IF (FOUND) THEN
           OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC))
         ELSE
           print *, "No Local fort.74 files found"
           RETURN
         ENDIF
         READ(LOC2(IPROC),'(A85)')  INLINE
         READ(LOC2(IPROC),3645) NDSETGV,NP,DTV,NSTEMP,ITEMPV
         IF (NP.NE.NNODP(IPROC)) THEN
           STOP 'Inconsistency in number of local nodes'
         ENDIF
      ENDDO
      OPEN(UNIT=74,FILE='fort.74')
C
      WRITE(74,'(A85)') INLINE
      WRITE(74,3645) NDSETGV,NNODG,DTV,NSPOOLGW,ITEMPV
C
C     jgf46.00 Made the loop over data sets an infinite loop because the
C     number of data sets at the top of the file will be wrong if the
C     run was hot started.
C
      DO ! loop infinitely over data sets
C
         DO IPROC=1,NPROC
            READ(LOC2(IPROC),FMT=2120,END=9999) TIMEOUTV,ITV
            IF (IPROC.EQ.1) WRITE(74,2120) TIMEOUTV,ITV
            DO K = 1,NNODP(IPROC)
               READ(LOC2(IPROC),80) UUP(K)
               INDX = IMAP_NOD_LG(K,IPROC)
               IF(IMAP_NOD_GL(1,INDX).EQ.IPROC) THEN
                UU(INDX) = UUP(K)
               ENDIF
            ENDDO
         ENDDO 
         DO I = 1,NNODG
            CALL NEWINDEX(UU(I),OUTMSG,I)
            WRITE(74,*) TRIM(OUTMSG)  !jgf46.00 TRIM to cut down file size
         ENDDO
      ENDDO               
      GO TO 9999
C
2000  CONTINUE
C
C--Open All Local and the Global Direct Access Binary fort.74 file
C
      DO IPROC = 1,NPROC
         LOCNAME(IPROC) = 'PE0000/fort.74'
         WRITE(LOCNAME(IPROC)(3:6),'(I4.4)') (IPROC-1)
         LOC2(IPROC) = 105 + (IPROC-1)
         INQUIRE(FILE=LOCNAME(IPROC),EXIST=FOUND)
         IF (FOUND) THEN
           OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC),ACCESS='DIRECT',
     &       RECL=NBYTE)
         ELSE
           print *, "No Local fort.74 files found"
           RETURN
         ENDIF
      ENDDO
      OPEN(74,FILE='fort.74',ACCESS='DIRECT',RECL=NBYTE)
C
C--Read RUNDES RUNID and AGRID from each Local file
C
      DO IPROC = 1,NPROC
         IREC(IPROC) = 0
         IF (NBYTE.EQ.4) THEN
           DO I=1,8
              READ(LOC2(IPROC),REC=IREC(IPROC)+I) RDES4(I)
           ENDDO
           IREC(IPROC)=IREC(IPROC)+8
           DO I=1,6
              READ(LOC2(IPROC),REC=IREC(IPROC)+I) RID4(I)
           ENDDO
           IREC(IPROC)=IREC(IPROC)+6
           DO I=1,6
              READ(LOC2(IPROC),REC=IREC(IPROC)+I) AID4(I)
           ENDDO
           IREC(IPROC)=IREC(IPROC)+6
         ENDIF
         IF (NBYTE.EQ.8) THEN
           DO I=1,4
              READ(LOC2(IPROC),REC=IREC(IPROC)+I) RDES8(I)
           ENDDO
           IREC(IPROC)=IREC(IPROC)+4
           DO I=1,3
              READ(LOC2(IPROC),REC=IREC(IPROC)+I) RID8(I)
           ENDDO
           IREC(IPROC)=IREC(IPROC)+3
           DO I=1,3
              READ(LOC2(IPROC),REC=IREC(IPROC)+I) AID8(I)
           ENDDO
           IREC(IPROC)=IREC(IPROC)+3
         ENDIF
      ENDDO
C
C--Write RUNDES RUNID and AGRID to Global File
C
      OREC = 0
      IF (NBYTE.EQ.4) THEN
        DO I=1,8
           WRITE(74,REC=OREC+I) RDES4(I)
        ENDDO
        OREC=OREC+8
        DO I=1,6
           WRITE(74,REC=OREC+I) RID4(I)
        ENDDO
        OREC=OREC+6
        DO I=1,6
           WRITE(74,REC=OREC+I) AID4(I)
        ENDDO
        OREC=OREC+6
      ENDIF
      IF (NBYTE.EQ.8) THEN
        DO I=1,4
           WRITE(74,REC=OREC+I) RDES8(I)
        ENDDO
        OREC=OREC+4
        DO I=1,3
           WRITE(74,REC=OREC+I) RID8(I)
        ENDDO
        OREC=OREC+3
        DO I=1,3
           WRITE(74,REC=OREC+I) AID8(I)
        ENDDO
        OREC=OREC+3
      ENDIF
C
C--Read NTRSPV, NSTAV, DT*NSPOOLV from each Local file
C  and then close files to flush file buffers
C
      DO IPROC = 1,NPROC
         READ(LOC2(IPROC),REC=IREC(IPROC)+1) NDSETGV
         READ(LOC2(IPROC),REC=IREC(IPROC)+2) NP
         READ(LOC2(IPROC),REC=IREC(IPROC)+3) DTV
         READ(LOC2(IPROC),REC=IREC(IPROC)+4) NSTEMP
         READ(LOC2(IPROC),REC=IREC(IPROC)+5) ITEMPV
         IREC(IPROC) = IREC(IPROC)+5
         CLOSE(LOC2(IPROC))    ! Flush the Write Buffer
         IF (NP.NE.NNODP(IPROC)) THEN
           STOP 'Inconsistency in number of local nodes'
         ENDIF
         OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC),ACCESS='DIRECT',RECL=NBYTE)
      ENDDO
C
C--Write same info to Global file and close it to flush buffer also
C
      WRITE(74,REC=OREC+1) NDSETGV
      WRITE(74,REC=OREC+2) NNODG
      WRITE(74,REC=OREC+3) DTV
      WRITE(74,REC=OREC+4) NSPOOLGW
      WRITE(74,REC=OREC+5) ITEMPV
      OREC = OREC+5
      CLOSE(74)         ! Flush the Write Buffer
      OPEN(74,FILE='fort.74',ACCESS='DIRECT',RECL=NBYTE)
C
      DO J=1,NDSETGV
C
         DO IPROC=1,NPROC
            READ(LOC2(IPROC),REC=IREC(IPROC)+1) TIMEOUTV
            READ(LOC2(IPROC),REC=IREC(IPROC)+2) ITV
            IREC(IPROC) = IREC(IPROC) + 2
            DO K=1, NNODP(IPROC)
               READ(LOC2(IPROC),REC=IREC(IPROC)+2*K-1) UUBINP(K)
               READ(LOC2(IPROC),REC=IREC(IPROC)+2*K) VVBINP(K)
            ENDDO
            DO K=1,NNODP(IPROC)
               UUBIN(IMAP_NOD_LG(K,IPROC)) = UUBINP(K)
               VVBIN(IMAP_NOD_LG(K,IPROC)) = VVBINP(K)
            ENDDO
           IREC(IPROC) = IREC(IPROC) + 2*NNODP(IPROC)
         ENDDO
C
         WRITE(74,REC=OREC+1) TIMEOUTV
         WRITE(74,REC=OREC+2) ITV
         OREC = OREC + 2
         DO K=1, NNODG
            WRITE(74,REC=OREC+2*K-1) UUBIN(K)
            WRITE(74,REC=OREC+2*K)   VVBIN(K)
         ENDDO
         OREC = OREC + 2*NNODG
C
      ENDDO
C
9999  CONTINUE
C
C--Close the Global and Local fort.74 Files
C
      CLOSE(74)
      DO IPROC = 1,NPROC
         CLOSE(LOC2(IPROC))
      ENDDO
C
   80 FORMAT(A80)
 2120 FORMAT(2X,E20.10,5X,I10)
c2454 FORMAT(2X,I8,2(2X,E15.8))
 3645 FORMAT(1X,I10,1X,I10,1X,E15.7,1X,I5,1X,I5)
C
      RETURN
      END


C---------------------------------------------------------------------------
C                  S U B R O U T I N E   P O S T 9 1 
C---------------------------------------------------------------------------
C     Merge the file containing ice concentration data at the ice 
C     concentration stations (fort.91 file) from each subdomain
C     into a single, full domain fort.91 file. This version is
C     compatible with ADCIRC version 45.07. jgf 11/02/2005.
C
C     jgf46.02 Fixed fort.80 file to contain the station mapping, fixed
C     memory allocation bug in post_global.F (had MNSTAV instead of
C     MNSTAM).
C
C     v49.64.01 tcm -- added this routine, modified from post71
C---------------------------------------------------------------------------
      SUBROUTINE POST91()
      USE POST_GLOBAL
      IMPLICIT NONE
c     counters
      INTEGER I,J,K,L ! loop counters 
      INTEGER IPROC   ! subdomain counter
      INTEGER OREC    ! record counter in binary output file
      INTEGER,ALLOCATABLE :: IREC(:) ! input record counters for binary files
c     containers
      INTEGER IDUM                         ! to read number of nodes in binary
      CHARACTER*80  OUTMSG                 ! text line of full domain data 
      CHARACTER*85  INLINE                 ! hold RUNDES, RUNID, AGRID in text
      CHARACTER*4 RDES4(8),RID4(6),AID4(6) ! hold RUNDES, RUNID, AGRID in bin.
      CHARACTER*8 RDES8(8),RID8(6),AID8(6) ! hold RUNDES, RUNID, AGRID in bin.
c     file stuff 
      LOGICAL FOUND                        ! TRUE if the file exists
      INTEGER,ALLOCATABLE::LOC2(:)          ! unit numbers of subdomains
      CHARACTER*14,ALLOCATABLE::LOCNAME(:) ! subdomain dir./file names
c     actual data 
      INTEGER NTRSPM  ! number of data sets to be written to fort.91
      INTEGER NUMSTNS ! number of subdomain meteorlogical recording stations
      INTEGER NSTEMP  ! time step interval at which fort.91 is written
      INTEGER ITEMPM  ! record type (1=elev.,2=vel.,3=3D vel.)
      REAL(SZ) DTM    ! total time: DT*NSPOOLM = step size * no. time steps
      INTEGER ITSM    ! model time step number
      REAL(8)TIMEOUTSM! model time
      REAL(SZ),ALLOCATABLE::CiceBIN(:)    ! full domain ice concent. val.(bin)
      REAL(SZ),ALLOCATABLE::CiceBINP(:)   ! subdomain ice concent. val. (bin)
      INTEGER,ALLOCATABLE::StaNum(:)         ! recording station number
      CHARACTER(len=80),ALLOCATABLE::CiceASC(:) ! ice concentration data
C
      ALLOCATE ( LOC2(MNPROC),IREC(MNPROC) )
      ALLOCATE ( CiceBIN(MNSTAM),CiceBINP(MNSTAM))
      ALLOCATE ( CiceASC(MNSTAM))
      ALLOCATE ( LOCNAME(MNPROC))
C
C     Determine whether Unit 91 is Sequential Formatted or Direct Access Binary
C
      SELECT CASE(ABS(NOUTM))
C     -------      
      CASE(1)                   ! text (sequential formatted)
C     -------
C     Open All Local Sequential Formatted fort.91 files
         DO IPROC = 1,NPROC
            LOCNAME(IPROC) = 'PE0000/fort.91'
            WRITE(LOCNAME(IPROC)(3:6),'(I4.4)') (IPROC-1)
            LOC2(IPROC) = 105 + (IPROC-1)
            INQUIRE(FILE=LOCNAME(IPROC),EXIST=FOUND)
            IF (FOUND) THEN
               OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC))
            ELSE
               print *, "ERROR: No subdomain fort.91 files found."
               RETURN           ! *** EARLY RETURN ***
            ENDIF
            READ(LOC2(IPROC),'(A85)') INLINE
            READ(LOC2(IPROC),3645) NTRSPM,NUMSTNS,DTM,NSTEMP,ITEMPM
            IF (NUMSTNS.NE.NSTAMP(IPROC)) THEN
               STOP 'ERROR: Inconsistency in no. of ice conc. sta.'
            ENDIF
         ENDDO
C     
C     Open Global Sequential Formatted fort.91 file
         OPEN(UNIT=91,FILE='fort.91')
C     
         WRITE(91,'(A85)') INLINE
         WRITE(91,3645) NTRSPM,NUMSTNS,DTM,NSPOOLM,ITEMPM
C     
         DO                     ! loop infinitely over data sets
C     
            DO IPROC=1,NPROC
               READ(LOC2(IPROC),2120,END=9999) TIMEOUTSM,ITSM
               IF (NSTAMP(IPROC).GT.0) THEN
                  READ(LOC2(IPROC),80,END=9999) 
     &                 (CiceASC(IMAP_STAM_LG(K,IPROC)),
     &                 K=1,NSTAMP(IPROC))
               ENDIF
            ENDDO
C     
            WRITE(91,2120) TIMEOUTSM,ITSM
            DO I=1, NSTAM
               CALL NEWINDEX(CiceASC(I),OUTMSG,I)
               WRITE(91,*) TRIM(OUTMSG) !jgf46.00 TRIM to cut down file size 
            ENDDO
C     
         ENDDO
C     -------
      CASE(2)                   ! binary
C     -------
C     Open All Local and the Global Direct Access Binary fort.91 files
         DO IPROC = 1,NPROC
            LOCNAME(IPROC) = 'PE0000/fort.91'
            WRITE(LOCNAME(IPROC)(3:6),'(I4.4)') (IPROC-1)
            LOC2(IPROC) = 105 + (IPROC-1)
            INQUIRE(FILE=LOCNAME(IPROC),EXIST=FOUND)
            IF (FOUND) THEN
               OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC),ACCESS='DIRECT',
     &              RECL=NBYTE)
            ELSE
               print *, "ERROR: No subdomain fort.91 files found."
               RETURN
            ENDIF
         ENDDO
         OPEN(91,FILE='fort.91',ACCESS='DIRECT',RECL=NBYTE)
C     
C     
C--   Read RUNDES RUNID and AGRID from each subdomain file
C     
         DO IPROC = 1,NPROC
            IREC(IPROC) = 0
            IF (NBYTE.EQ.4) THEN
               DO I=1,8
                  READ(LOC2(IPROC),REC=IREC(IPROC)+I) RDES4(I)
               ENDDO
               IREC(IPROC)=IREC(IPROC)+8
               DO I=1,6
                  READ(LOC2(IPROC),REC=IREC(IPROC)+I) RID4(I)
               ENDDO
               IREC(IPROC)=IREC(IPROC)+6
               DO I=1,6
                  READ(LOC2(IPROC),REC=IREC(IPROC)+I) AID4(I)
               ENDDO
               IREC(IPROC)=IREC(IPROC)+6
            ENDIF
            IF (NBYTE.EQ.8) THEN
               DO I=1,4
                  READ(LOC2(IPROC),REC=IREC(IPROC)+I) RDES8(I)
               ENDDO
               IREC(IPROC)=IREC(IPROC)+4
               DO I=1,3
                  READ(LOC2(IPROC),REC=IREC(IPROC)+I) RID8(I)
               ENDDO
               IREC(IPROC)=IREC(IPROC)+3
               DO I=1,3
                  READ(LOC2(IPROC),REC=IREC(IPROC)+I) AID8(I)
               ENDDO
               IREC(IPROC)=IREC(IPROC)+3
            ENDIF
         ENDDO
C     
C--   Write RUNDES RUNID and AGRID to Global File
C     
         OREC = 0
         IF (NBYTE.EQ.4) THEN
            DO I=1,8
               WRITE(91,REC=OREC+I) RDES4(I)
            ENDDO
            OREC=OREC+8
            DO I=1,6
               WRITE(91,REC=OREC+I) RID4(I)
            ENDDO
            OREC=OREC+6
            DO I=1,6
               WRITE(91,REC=OREC+I) AID4(I)
            ENDDO
            OREC=OREC+6
         ENDIF
         IF (NBYTE.EQ.8) THEN
            DO I=1,4
               WRITE(91,REC=OREC+I) RDES8(I)
            ENDDO
            OREC=OREC+4
            DO I=1,3
               WRITE(91,REC=OREC+I) RID8(I)
            ENDDO
            OREC=OREC+3
            DO I=1,3
               WRITE(91,REC=OREC+I) AID8(I)
            ENDDO
            OREC=OREC+3
         ENDIF
C     
C--   Read NTRSPM, NSTAM, DT*NSPOOLM from each subdomain file
C     and then close files to flush file buffers
C     
         DO IPROC = 1,NPROC
            READ(LOC2(IPROC),REC=IREC(IPROC)+1) NTRSPM
            READ(LOC2(IPROC),REC=IREC(IPROC)+2) NUMSTNS
            READ(LOC2(IPROC),REC=IREC(IPROC)+3) DTM
            READ(LOC2(IPROC),REC=IREC(IPROC)+4) NSPOOLM 
            READ(LOC2(IPROC),REC=IREC(IPROC)+5) ITEMPM
            IREC(IPROC) = IREC(IPROC)+5
            CLOSE(LOC2(IPROC))   ! Flush the Write Buffer
            IF (NSTAM.NE.NSTAMP(IPROC)) THEN
               STOP 'ERROR: Inconsistency in no. of ice conc. sta.'
            ENDIF
            OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC),ACCESS='DIRECT',
     &           RECL=NBYTE)
         ENDDO
C     
C--   Write same info to full-domain file and close it to flush buffer also
C     
         WRITE(91,REC=OREC+1) NTRSPM
         WRITE(91,REC=OREC+2) NUMSTNS
         WRITE(91,REC=OREC+3) DTM
         WRITE(91,REC=OREC+4) NSPOOLM
         WRITE(91,REC=OREC+5) ITEMPM
         OREC = OREC+5
         CLOSE(91)              ! Flush the Write Buffer
         OPEN(91,FILE='fort.91',ACCESS='DIRECT',RECL=NBYTE)
C     
         DO J=1,NTRSPM
C     
            DO IPROC=1,NPROC
               READ(LOC2(IPROC),REC=IREC(IPROC)+1) TIMEOUTSM
               READ(LOC2(IPROC),REC=IREC(IPROC)+2) ITSM
               IREC(IPROC) = IREC(IPROC) + 2
               IF (NSTAMP(IPROC).GT.0) THEN
                  DO K=1, NSTAMP(IPROC)
                     READ(LOC2(IPROC),REC=IREC(IPROC)+K) CiceBINP(K)
                  ENDDO
                  DO K=1,NSTAMP(IPROC)
                     CiceBIN(IMAP_STAM_LG(K,IPROC)) = CiceBINP(K)
                  ENDDO
                  IREC(IPROC) = IREC(IPROC) + NSTAMP(IPROC)
               ENDIF
            ENDDO
C     
            WRITE(91,REC=OREC+1) TIMEOUTSM
            WRITE(91,REC=OREC+2) ITSM
            OREC = OREC + 2
            DO K=1, NSTAM
               WRITE(91,REC=OREC+K) CiceBIN(K)
            ENDDO
            OREC = OREC + NSTAM
C     
         ENDDO
C
      CASE DEFAULT
C     
         WRITE(*,2125) NOUTM
C
      END SELECT
C     
C     Close the full domain and subdomain fort.91 Files
C
 9999 CONTINUE
      CLOSE(91)
      DO IPROC = 1,NPROC
         CLOSE(LOC2(IPROC))
      ENDDO
C     
 80   FORMAT(A40)
 2120 FORMAT(2X,E20.10,5X,I10)
 2125 FORMAT(2X,'ERROR: Incorrect value of NOUTM; NOUTM=',I3) 
 3645 FORMAT(1X,I10,1X,I10,1X,E15.7,1X,I5,1X,I5)
C     
      RETURN
C---------------------------------------------------------------------------
      END SUBROUTINE POST91
C---------------------------------------------------------------------------


C---------------------------------------------------------------------------
C                  S U B R O U T I N E   P O S T 9 3 
C---------------------------------------------------------------------------
C     Merge the ice concentration data (fort.93 files) from each
C     subdomain into a single, full domain fort.93 file. This version is
C     compatible with ADCIRC version 45.07. jgf 11/02/2005.
C
C     v49.64.01 tcm -- added
C---------------------------------------------------------------------------
      SUBROUTINE POST93()
      USE POST_GLOBAL
      IMPLICIT NONE
C     counters
      INTEGER I,J,K,L ! loop counters 
      INTEGER IPROC   ! subdomain counter
      INTEGER OREC    ! record counter in binary output file
      INTEGER, ALLOCATABLE::IREC(:) ! input record counters for binary files
      INTEGER INDX   ! full domain node number
C     containers
      CHARACTER*80 OUTMSG                  ! text line of full domain data 
      CHARACTER*85 INLINE                  ! hold RUNDES, RUNID, AGRID in text
      CHARACTER*4 RDES4(8),RID4(6),AID4(6) ! hold RUNDES, RUNID, AGRID in bin.
      CHARACTER*8 RDES8(8),RID8(6),AID8(6) ! hold RUNDES, RUNID, AGRID in bin.
      INTEGER IDUM                         ! to read number of nodes in binary
C     file stuff
      INTEGER, ALLOCATABLE::LOC2(:)         ! unit numbers of subdomains
      CHARACTER*14, ALLOCATABLE::LOCNAME(:)! subdomain dir./file names
      LOGICAL FOUND                        ! TRUE if the file exists
C     actual data
      INTEGER NDSETSGW! number of data sets to be written to fort.93
      INTEGER NP      ! number of nodes in the subdomain grid
      REAL(SZ) DTGW   ! total time: DT*NSPOOLGW = step size * no. time steps
      INTEGER NSTEMP  ! time step interval at which fort.93 is written
      INTEGER ITEMPGW ! record type (1=elev.,2=vel.,3=3D vel.)
C
      REAL(8)TIMEOUTGW! model time
      INTEGER ITGW    ! model time step
C
      CHARACTER*80,ALLOCATABLE::CiceSC(:) !full domain ice concentration values
      CHARACTER*80,ALLOCATABLE::CiceSCP(:)!subdomain ice concentration values
      REAL(SZ),ALLOCATABLE::CiceBIN(:)    ! full domain ice conc. val.(bin)
      REAL(SZ),ALLOCATABLE::CiceBINP(:)   ! subdomain ice conc. val. (bin)
C
      ALLOCATE (LOC2(MNPROC),IREC(MNPROC))
      ALLOCATE (CiceBIN(MNP),CiceBINP(MNPP))
      ALLOCATE (CiceSC(MNP),CiceSCP(MNPP))
      ALLOCATE ( LOCNAME(MNPROC))
C     
C--   Determine whether Unit 93 is Sequential Formatted or Direct Access Binary
C     
      IF (ABS(NOUTGW).EQ.1) THEN
         GO TO 1000 ! file is text
      ELSE
         GO TO 2000 ! file is binary(=2) or nothing (=0)
      ENDIF
C
C
 1000 CONTINUE !here begins the post-processing if the files are ASCII text
C
      DO IPROC = 1,NPROC
C     Open ASCII text fort.93 files in each subdirectory
         LOCNAME(IPROC) = 'PE0000/fort.93'
         WRITE(LOCNAME(IPROC)(3:6),'(I4.4)') (IPROC-1)
         LOC2(IPROC) = 105 + (IPROC-1)
         INQUIRE(FILE=LOCNAME(IPROC),EXIST=FOUND)
         IF (FOUND) THEN
            OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC))
         ELSE
            print *, "ERROR: No local fort.93 files found."
            RETURN !**** EARLY RETURN ****
         ENDIF
C     Read header information from each fort.93 file
         READ(LOC2(IPROC),'(A85)') INLINE ! read RUNDES, RUNID, AGRID
         READ(LOC2(IPROC),3645)  NDSETSGW,NP,DTGW,NSTEMP,ITEMPGW
         IF (NP.NE.NNODP(IPROC)) THEN
            print *, "NP = ",NP," NNODP(",IPROC,") = ",NNODP(IPROC)
            STOP 'ERROR: Inconsistency in expected subdomain nodes.'
         ENDIF
      ENDDO
C     Create fort.93 file for full domain
      OPEN(UNIT=93,FILE='fort.93')
      WRITE(93,'(A85)') INLINE ! write RUNDES, RUNID, AGRID
      WRITE(93,3645) NDSETSGW,NNODG,DTGW,NSPOOLGW,ITEMPGW 
C
C     jgf46.00 Made the loop over data sets an infinite loop because the
C     number of data sets at the top of the file will be wrong if the
C     run was hot started.
C
      DO ! loop infinitely over data sets
C
         DO IPROC=1,NPROC  ! loop over subdomains
            READ(LOC2(IPROC),FMT=2120,END=9999) TIMEOUTGW,ITGW
            IF (IPROC.EQ.1) WRITE(93,2120) TIMEOUTGW,ITGW    
            DO K = 1,NNODP(IPROC) ! loop over subdomain nodes
               READ(LOC2(IPROC),80,END=9999) CiceSCP(K)
               INDX = IMAP_NOD_LG(K,IPROC) ! sub --> full (lookup table)
               IF (IMAP_NOD_GL(1,INDX).EQ.IPROC) THEN
                  CiceSC(INDX) = CiceSCP(K)
               ENDIF
            ENDDO ! next subdomain node  
         ENDDO ! next subdomain
C
         DO I = 1,NNODG   ! loop over full domain nodes
            CALL NEWINDEX(CiceSC(I),OUTMSG,I)
            WRITE(93,*) TRIM(OUTMSG) !jgf46.00 TRIM to cut down file size 
         ENDDO 
      ENDDO ! next output time
      GO TO 9999 ! jump to end of subroutine where fort.93 file is closed
C
C     
 2000 CONTINUE !here begins the post-processing if files are binary
C--   Open all subdomain and the full domain Direct Access Binary fort.93 files
      DO IPROC = 1,NPROC
         LOCNAME(IPROC) = 'PE0000/fort.93'
         WRITE(LOCNAME(IPROC)(3:6),'(I4.4)') (IPROC-1)
         LOC2(IPROC) = 105 + (IPROC-1)
         INQUIRE(FILE=LOCNAME(IPROC),EXIST=FOUND)
         IF (FOUND) THEN
            OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC),ACCESS='DIRECT',
     &           RECL=NBYTE)
         ELSE
            print *, "ERROR: No local fort.93 files found."
            RETURN
         ENDIF
      ENDDO
      OPEN(93,FILE='fort.93',ACCESS='DIRECT',RECL=NBYTE)
C--   Read RUNDES RUNID and AGRID from each subdomain file
      DO IPROC = 1,NPROC
         IREC(IPROC) = 0
         IF (NBYTE.EQ.4) THEN
            DO I=1,8
               READ(LOC2(IPROC),REC=IREC(IPROC)+I) RDES4(I)
            ENDDO
            IREC(IPROC)=IREC(IPROC)+8
            DO I=1,6
               READ(LOC2(IPROC),REC=IREC(IPROC)+I) RID4(I)
            ENDDO
            IREC(IPROC)=IREC(IPROC)+6
            DO I=1,6
               READ(LOC2(IPROC),REC=IREC(IPROC)+I) AID4(I)
            ENDDO
            IREC(IPROC)=IREC(IPROC)+6
         ENDIF
         IF (NBYTE.EQ.8) THEN
            DO I=1,4
               READ(LOC2(IPROC),REC=IREC(IPROC)+I) RDES8(I)
            ENDDO
            IREC(IPROC)=IREC(IPROC)+4
            DO I=1,3
               READ(LOC2(IPROC),REC=IREC(IPROC)+I) RID8(I)
            ENDDO
            IREC(IPROC)=IREC(IPROC)+3
            DO I=1,3
               READ(LOC2(IPROC),REC=IREC(IPROC)+I) AID8(I)
            ENDDO
            IREC(IPROC)=IREC(IPROC)+3
         ENDIF
      ENDDO
C--   Write RUNDES RUNID and AGRID to full domain file
      OREC = 0
      IF (NBYTE.EQ.4) THEN
         DO I=1,8
            WRITE(93,REC=OREC+I) RDES4(I)
         ENDDO
         OREC=OREC+8
         DO I=1,6
            WRITE(93,REC=OREC+I) RID4(I)
         ENDDO
         OREC=OREC+6
         DO I=1,6
            WRITE(93,REC=OREC+I) AID4(I)
         ENDDO
         OREC=OREC+6
      ENDIF
      IF (NBYTE.EQ.8) THEN
         DO I=1,4
            WRITE(93,REC=OREC+I) RDES8(I)
         ENDDO
         OREC=OREC+4
         DO I=1,3
            WRITE(93,REC=OREC+I) RID8(I)
         ENDDO
         OREC=OREC+3
         DO I=1,3
            WRITE(93,REC=OREC+I) AID8(I)
         ENDDO
         OREC=OREC+3
      ENDIF
C     
C--   Read NDSETSGW, NNODP, DTGW, NSPOOLGW, ITEMPGW from each subdomain file
C     and then close files to flush file buffers
C     
      DO IPROC = 1,NPROC
         READ(LOC2(IPROC),REC=IREC(IPROC)+1) NDSETSGW
         READ(LOC2(IPROC),REC=IREC(IPROC)+2) IDUM   
         READ(LOC2(IPROC),REC=IREC(IPROC)+3) DTGW
         READ(LOC2(IPROC),REC=IREC(IPROC)+4) NSTEMP
         READ(LOC2(IPROC),REC=IREC(IPROC)+5) ITEMPGW
         IREC(IPROC) = IREC(IPROC)+5
         CLOSE(LOC2(IPROC))      ! flush the write buffer
         OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC),ACCESS='DIRECT',RECL=NBYTE)
      ENDDO
C     
C--   Write same info to full domain file and close it to flush buffer also
C     
      WRITE(93,REC=OREC+1) NDSETSGW
      WRITE(93,REC=OREC+2) NNODG
      WRITE(93,REC=OREC+3) DTGW
      WRITE(93,REC=OREC+4) NSPOOLGW
      WRITE(93,REC=OREC+5) ITEMPGW
      OREC = OREC+5
      CLOSE(93)                 ! flush the write buffer
      OPEN(93,FILE='fort.93',ACCESS='DIRECT',RECL=NBYTE)
C     
      DO J=1,NDSETSGW
C     
         DO IPROC=1,NPROC
            READ(LOC2(IPROC),REC=IREC(IPROC)+1) TIMEOUTGW
            READ(LOC2(IPROC),REC=IREC(IPROC)+2) ITGW
            IREC(IPROC) = IREC(IPROC) + 2
            DO K=1, NNODP(IPROC)
               READ(LOC2(IPROC),REC=IREC(IPROC)+K) CiceBINP(K)
            ENDDO
            DO K=1,NNODP(IPROC)
               INDX = IMAP_NOD_LG(K,IPROC)
               IF (IMAP_NOD_GL(1,INDX).EQ.IPROC) THEN
                  CiceBIN(INDX) = CiceBINP(K)
               ENDIF
            ENDDO
            IREC(IPROC) = IREC(IPROC) + NNODP(IPROC)
         ENDDO
C     
         WRITE(93,REC=OREC+1) TIMEOUTGW
         WRITE(93,REC=OREC+2) ITGW
         OREC = OREC + 2
         DO K=1, NNODG
            WRITE(93,REC=OREC+K) CiceBIN(K)
         ENDDO
         OREC = OREC + NNODG
C     
      ENDDO
C     
 9999 CONTINUE
C     
C--   Close the full domain and subdomain fort.93 Files
C     
      CLOSE(93)
      DO IPROC = 1,NPROC
         CLOSE(LOC2(IPROC))
      ENDDO
C     
 80   FORMAT(A80)
 2120 FORMAT(2X,E20.10,5X,I10)
 2453 FORMAT(2X,I8,2X,E15.8)
 3645 FORMAT(1X,I10,1X,I10,1X,E15.7,1X,I5,1X,I5)
C     
      RETURN
      END SUBROUTINE POST93
C---------------------------------------------------------------------------
C     End of subroutine Post93 
C---------------------------------------------------------------------------





      SUBROUTINE POST67_68()
      USE POST_GLOBAL 
C
C---------------------------------------------------------------------------C
C                     written 10/11/01 by RL                                C 
C                                                                           C  
C  This routine gathers the local hot start files (fort.67 and fort.68) and C
C  writes a global hotstart file                                            C
C  jgf Updated for v45.06 09/07/2005                                        C
C                                                                           C
C---------------------------------------------------------------------------C
C
      USE HARM, ONLY : FNAM8, HAFREQ, HAFF, HAFACE, HA, STAELV, STAULV,
     &                  STAVLV, GLOELV, GLOULV, GLOVLV, ELAV, ELVA, 
     &                  CHARMV, ICHA, IHARIND, ITHAS, MM, XELAV, YELAV,
     &                  XELVA, YELVA, ICALL, ITMV, ITUD, NF, NFREQ, NZ,
     &                  TIMEUD, NAMEFR, NTSTEPS
      IMPLICIT NONE
      LOGICAL FOUND
C     LOGICAL CHARMV

      INTEGER I,J,IPROC,LOCHSF,INDX,IHOTSTP
      INTEGER IMHSF,ITHSF,IHOT,NHSFILES
      INTEGER IESTP,NSCOUE,IVSTP,NSCOUV,ICSTP,NSCOUC,IPSTP,IWSTP,NSCOUM,
     &        IGEP,NSCOUGE,IGVP,NSCOUGV,IGCP,NSCOUGC,IGPP,IGWP,NSCOUGW
      CHARACTER FNAME*60,LOCFN*14
      INTEGER,ALLOCATABLE  :: LOC2(:),NODECODE(:),NOFF(:)
      REAL(SZ),ALLOCATABLE :: ETA1(:),ETA2(:),UU2(:),VV2(:),CH1(:)
C  kmd48.33bc updates for 3d output of hotstart information
      REAL(SZ),ALLOCATABLE :: EtaDisc(:)
      REAL(SZ),ALLOCATABLE :: DUU(:),DUV(:),DVV(:),UU(:),VV(:)
      REAL(SZ),ALLOCATABLE :: BSX(:),BSY(:)
      REAL(SZ),ALLOCATABLE :: REALQ(:,:),IMAGQ(:,:),WZ(:,:)
      REAL(SZ),ALLOCATABLE :: q20(:,:),l(:,:),SIGT(:,:)
      REAL(SZ),ALLOCATABLE :: SAL(:,:),TEMP(:,:)

      REAL(8) TIMEHSF
      INTEGER :: NSTAEX, NSTAVX

C  kmd48.33bc updates for 3d output of hotstart information
      INTEGER InputFileFmtVn,K
      INTEGER N3DSD,I3DSDRec,N3DSV,I3DSVRec,N3DST,I3DSTRec,N3DGD,
     &        I3DGDRec,N3DGV,I3DGVRec,N3DGT,I3DGTRec
      INTEGER NPX,NEX,NPG,NEG
      !
      INTEGER :: N
C
C Allocate local work arrays
C
      ALLOCATE ( LOC2(MNPROC) )
      ALLOCATE ( ETA1(MNP),ETA2(MNP),UU2(MNP),
     &     VV2(MNP),NODECODE(MNP),CH1(MNP) )
      ALLOCATE ( NOFF(MNE) )

C  kmd48.33bc updates for 3d output of hotstart information
      ALLOCATE ( DUU(MNP),DUV(MNP),DVV(MNP),UU(MNP),
     &           VV(MNP),BSX(MNP),BSY(MNP) )
      ALLOCATE ( EtaDisc(MNP) )


C--Determine which hotstart files to process

      WRITE(*,*) 'Globalize hotstart file 67,68 or both ( 67, 68, 135 )'
      READ(*,*) IHOT
      NHSFILES=1
      IF(IHOT.EQ.135) THEN
         NHSFILES=2
         IHOT=67
      ENDIF

C  kmd48.33bc add information for determining the number of layers
C             so the 3D portion can be evaluated
      WRITE(*,*) 'How many levels do we need to process?'
      READ(*,*) NFEN

C  kmd48.33bc allocate arrays that need NFEN

      ALLOCATE ( REALQ(MNP,NFEN),IMAGQ(MNP,NFEN),WZ(MNP,NFEN),
     &           q20(MNP,NFEN),l(MNP,NFEN) )
      ALLOCATE ( SIGT(MNP,NFEN),SAL(MNP,NFEN),TEMP(MNP,NFEN) )

      DO J=1,NHSFILES
         IHOT=IHOT+J-1  

C--Set up unit numbers for all Local Hot Start files 

         DO IPROC = 1,NPROC
            LOC2(IPROC) = 105 + (IPROC-1)
         ENDDO

C--Read common info from highest processor hotstart file

C  kmd48.33bc add information for hot start
         IF(IHOT.EQ.67) LOCFN(1:14) = 'PE0000/fort.67'
         IF(IHOT.EQ.68) LOCFN(1:14) = 'PE0000/fort.68' 
         LOCHSF=LOC2(NPROC)
         WRITE(LOCFN(3:6),'(I4.4)') (NPROC-1)
         OPEN (LOCHSF,FILE=LOCFN,ACCESS='DIRECT',RECL=8)
         IHOTSTP=1
         READ(LOCHSF,REC=IHOTSTP) InputFileFmtVn
         IHOTSTP=2
         READ(LOCHSF,REC=IHOTSTP) IMHSF
         IHOTSTP=3
         READ(LOCHSF,REC=IHOTSTP) TIMEHSF
         IHOTSTP=4
         READ(LOCHSF,REC=IHOTSTP) ITHSF
         IHOTSTP=5
         READ(LOCHSF,REC=IHOTSTP) NPG
         IHOTSTP=6
         READ(LOCHSF,REC=IHOTSTP) NEG
         IHOTSTP=7
         READ(LOCHSF,REC=IHOTSTP) NPX
         IHOTSTP=8
         READ(LOCHSF,REC=IHOTSTP) NEX
         CLOSE(LOCHSF)          !close to flush the buffer

C--Read stuff from each local file  
C  kmd48.33bc add in information for node and elements
         NPG=MNP
         NEG=MNE

C  kmd48.33bc update to utilize the global IO
         DO IPROC=1,NPROC
            IF(IHOT.EQ.67) LOCFN(1:14) = 'PE0000/fort.67'
            IF(IHOT.EQ.68) LOCFN(1:14) = 'PE0000/fort.68' 
            LOCHSF=LOC2(IPROC)
            WRITE(LOCFN(3:6),'(I4.4)') (IPROC-1)
            OPEN (LOCHSF,FILE=LOCFN,ACCESS='DIRECT',RECL=8)
            IHOTSTP=8
            DO I=1,NNODP(IPROC)
               INDX=IMAP_NOD_LG(I,IPROC)
               READ(LOCHSF,REC=IHOTSTP+1) ETA1(INDX)
               IHOTSTP=IHOTSTP+1
            END DO
            DO I=1,NNODP(IPROC)
               INDX=IMAP_NOD_LG(I,IPROC)
               READ(LOCHSF,REC=IHOTSTP+1) ETA2(INDX)
               IHOTSTP=IHOTSTP+1
            END DO
            DO I=1,NNODP(IPROC)
               INDX=IMAP_NOD_LG(I,IPROC)
               READ(LOCHSF,REC=IHOTSTP+1) EtaDisc(INDX)
               IHOTSTP=IHOTSTP+1
            END DO
            DO I=1,NNODP(IPROC)
               INDX=IMAP_NOD_LG(I,IPROC)
               READ(LOCHSF,REC=IHOTSTP+1) UU2(INDX)
               IHOTSTP=IHOTSTP+1
            END DO
            DO I=1,NNODP(IPROC)
               INDX=IMAP_NOD_LG(I,IPROC)
               READ(LOCHSF,REC=IHOTSTP+1) VV2(INDX)
               IHOTSTP=IHOTSTP+1
            END DO
            IF(IM.EQ.10) THEN
               DO I=1,NNODP(IPROC)
                 INDX=IMAP_NOD_LG(I,IPROC)
                 READ(LOCHSF,REC=IHOTSTP+1) CH1(INDX)
                 IHOTSTP=IHOTSTP+1
               END DO
            ENDIF
            DO I=1,NNODP(IPROC)
               INDX=IMAP_NOD_LG(I,IPROC)
               READ(LOCHSF,REC=IHOTSTP+1) NODECODE(INDX)
               IHOTSTP=IHOTSTP+1
            ENDDO
            DO I=1,NELP(IPROC)
               INDX=IMAP_EL_LG(I,IPROC)
               READ(LOCHSF,REC=IHOTSTP+1) NOFF(INDX)
               IHOTSTP=IHOTSTP+1
            END DO
           CLOSE(LOCHSF)
         END DO 

!            DO I=1,NNODP(IPROC)
!               INDX=IMAP_NOD_LG(I,IPROC)
!               READ(LOCHSF,REC=IHOTSTP+1) ETA1(INDX)
!               READ(LOCHSF,REC=IHOTSTP+2) ETA2(INDX)
!               READ(LOCHSF,REC=IHOTSTP+3) UU2(INDX)
!               READ(LOCHSF,REC=IHOTSTP+4) VV2(INDX)
!               IHOTSTP=IHOTSTP+4
!c              IF(IM.EQ.10) THEN
!c                 READ(LOCHSF,REC=IHOTSTP+1) CH1(INDX)
!c                 IHOTSTP=IHOTSTP+1
!c              ENDIF
!               READ(LOCHSF,REC=IHOTSTP+1) NODECODE(INDX)
!               IHOTSTP=IHOTSTP+1
!            ENDDO
!            DO I=1,NELP(IPROC)
!               INDX=IMAP_EL_LG(I,IPROC)
!               READ(LOCHSF,REC=IHOTSTP+1) NOFF(INDX)
!               IHOTSTP=IHOTSTP+1                                    
!            END DO

C--Read in more common info from higest processor hotstart file
         
         IF(IHOT.EQ.67) LOCFN(1:14) = 'PE0000/fort.67'
         IF(IHOT.EQ.68) LOCFN(1:14) = 'PE0000/fort.68' 
         LOCHSF=LOC2(NPROC)
         WRITE(LOCFN(3:6),'(I4.4)') (NPROC-1)
         OPEN (LOCHSF,FILE=LOCFN,ACCESS='DIRECT',RECL=8)
         READ(LOCHSF,REC=IHOTSTP+1 ) IESTP
         READ(LOCHSF,REC=IHOTSTP+2 ) NSCOUE
         READ(LOCHSF,REC=IHOTSTP+3 ) IVSTP
         READ(LOCHSF,REC=IHOTSTP+4 ) NSCOUV
         READ(LOCHSF,REC=IHOTSTP+5 ) ICSTP
         READ(LOCHSF,REC=IHOTSTP+6 ) NSCOUC
         READ(LOCHSF,REC=IHOTSTP+7 ) IPSTP
         READ(LOCHSF,REC=IHOTSTP+8 ) IWSTP
         READ(LOCHSF,REC=IHOTSTP+9 ) NSCOUM
         READ(LOCHSF,REC=IHOTSTP+10) IGEP
         READ(LOCHSF,REC=IHOTSTP+11) NSCOUGE
         READ(LOCHSF,REC=IHOTSTP+12) IGVP
         READ(LOCHSF,REC=IHOTSTP+13) NSCOUGV
         READ(LOCHSF,REC=IHOTSTP+14) IGCP
         READ(LOCHSF,REC=IHOTSTP+15) NSCOUGC
         READ(LOCHSF,REC=IHOTSTP+16) IGPP
         READ(LOCHSF,REC=IHOTSTP+17) IGWP
         READ(LOCHSF,REC=IHOTSTP+18) NSCOUGW

C  kmd48.33bc add in the 3D barotropic and baroclinic information
         IF ((IM.EQ.1).OR.(IM.EQ.2).OR.(IM.EQ.11)
     &          .OR.(IM.EQ.21).OR.(IM.EQ.31)) THEN
            READ(LOCHSF,REC=IHOTSTP+19) IDEN
            READ(LOCHSF,REC=IHOTSTP+20) N3DSD
            READ(LOCHSF,REC=IHOTSTP+21) I3DSDRec
            READ(LOCHSF,REC=IHOTSTP+22) N3DSV
            READ(LOCHSF,REC=IHOTSTP+23) I3DSVRec
            READ(LOCHSF,REC=IHOTSTP+24) N3DST
            READ(LOCHSF,REC=IHOTSTP+25) I3DSTRec
            READ(LOCHSF,REC=IHOTSTP+26) N3DGD
            READ(LOCHSF,REC=IHOTSTP+27) I3DGDRec
            READ(LOCHSF,REC=IHOTSTP+28) N3DGV
            READ(LOCHSF,REC=IHOTSTP+29) I3DGVRec
            READ(LOCHSF,REC=IHOTSTP+30) N3DGT
            READ(LOCHSF,REC=IHOTSTP+31) I3DGTRec
            IHOTSTP=IHOTSTP+31
         END IF

         CLOSE(LOCHSF)

C  kmd48.33bc add in the 3D barotropic and baroclinic information

         IF ((IM.EQ.1).OR.(IM.EQ.2).OR.(IM.EQ.11)
     &          .OR.(IM.EQ.21).OR.(IM.EQ.31)) THEN
            DO IPROC=1,NPROC
               IF(IHOT.EQ.67) LOCFN(1:14) = 'PE0000/fort.67'
               IF(IHOT.EQ.68) LOCFN(1:14) = 'PE0000/fort.68'
               LOCHSF=LOC2(IPROC)
               WRITE(LOCFN(3:6),'(I4.4)') (IPROC-1)
               OPEN (LOCHSF,FILE=LOCFN,ACCESS='DIRECT',RECL=8)
               IF (IM.EQ.10) THEN
                  IHOTSTP=8+(NNODP(IPROC)*7)+(NELP(IPROC)*1)+31
               ELSE
                  IHOTSTP=8+(NNODP(IPROC)*6)+(NELP(IPROC)*1)+31
               END IF
               DO I=1,NNODP(IPROC)
                 INDX=IMAP_NOD_LG(I,IPROC)
                 READ(LOCHSF,REC=IHOTSTP+1) DUU(INDX)
                 IHOTSTP=IHOTSTP+1
               END DO
               DO I=1,NNODP(IPROC)
                 INDX=IMAP_NOD_LG(I,IPROC)
                 READ(LOCHSF,REC=IHOTSTP+1) DUV(INDX)
                 IHOTSTP=IHOTSTP+1
               END DO
               DO I=1,NNODP(IPROC)
                 INDX=IMAP_NOD_LG(I,IPROC)
                 READ(LOCHSF,REC=IHOTSTP+1) DVV(INDX)
                 IHOTSTP=IHOTSTP+1
               END DO
               DO I=1,NNODP(IPROC)
                 INDX=IMAP_NOD_LG(I,IPROC)
                 READ(LOCHSF,REC=IHOTSTP+1) UU(INDX)
                 IHOTSTP=IHOTSTP+1
               END DO
               DO I=1,NNODP(IPROC)
                 INDX=IMAP_NOD_LG(I,IPROC)
                 READ(LOCHSF,REC=IHOTSTP+1) VV(INDX)
                 IHOTSTP=IHOTSTP+1
               END DO
               DO I=1,NNODP(IPROC)
                 INDX=IMAP_NOD_LG(I,IPROC)
                 READ(LOCHSF,REC=IHOTSTP+1) BSX(INDX)
                 IHOTSTP=IHOTSTP+1
               END DO
               DO I=1,NNODP(IPROC)
                 INDX=IMAP_NOD_LG(I,IPROC)
                 READ(LOCHSF,REC=IHOTSTP+1) BSY(INDX)
                 IHOTSTP=IHOTSTP+1
               END DO
               DO K=1,NFEN
                 DO I=1,NNODP(IPROC)
                   INDX=IMAP_NOD_LG(I,IPROC)
                   READ(LOCHSF,REC=IHOTSTP+1) REALQ(INDX,K)
                   IHOTSTP=IHOTSTP+1
                 END DO
               END DO
               DO K=1,NFEN
                 DO I=1,NNODP(IPROC)
                   INDX=IMAP_NOD_LG(I,IPROC)
                   READ(LOCHSF,REC=IHOTSTP+1) IMAGQ(INDX,K)
                   IHOTSTP=IHOTSTP+1
                 END DO
               END DO
               DO K=1,NFEN
                 DO I=1,NNODP(IPROC)
                   INDX=IMAP_NOD_LG(I,IPROC)
                   READ(LOCHSF,REC=IHOTSTP+1) WZ(INDX,K)
                   IHOTSTP=IHOTSTP+1
                 END DO
               END DO
               DO K=1,NFEN
                 DO I=1,NNODP(IPROC)
                   INDX=IMAP_NOD_LG(I,IPROC)
                   READ(LOCHSF,REC=IHOTSTP+1) q20(INDX,K)
                   IHOTSTP=IHOTSTP+1
                 END DO
               END DO
               DO K=1,NFEN
                 DO I=1,NNODP(IPROC)
                   INDX=IMAP_NOD_LG(I,IPROC)
                   READ(LOCHSF,REC=IHOTSTP+1) l(INDX,K)
                   IHOTSTP=IHOTSTP+1
                 END DO
               END DO
               IF (ABS(IDEN).EQ.1) THEN
                 DO K=1,NFEN
                   DO I=1,NNODP(IPROC)
                     INDX=IMAP_NOD_LG(I,IPROC)
                     READ(LOCHSF,REC=IHOTSTP+1) SIGT(INDX,K)
                     IHOTSTP=IHOTSTP+1
                   END DO
                 END DO
               ELSE IF (ABS(IDEN).EQ.2) THEN
                 DO K=1,NFEN
                   DO I=1,NNODP(IPROC)
                     INDX=IMAP_NOD_LG(I,IPROC)
                     READ(LOCHSF,REC=IHOTSTP+1) SAL(INDX,K)
                     IHOTSTP=IHOTSTP+1
                   END DO
                 END DO
               ELSE IF (ABS(IDEN).EQ.3) THEN
                 DO K=1,NFEN
                   DO I=1,NNODP(IPROC)
                     INDX=IMAP_NOD_LG(I,IPROC)
                     READ(LOCHSF,REC=IHOTSTP+1) TEMP(INDX,K)
                     IHOTSTP=IHOTSTP+1
                   END DO
                 END DO
               ELSE IF (ABS(IDEN).EQ.4) THEN
                 DO K=1,NFEN
                   DO I=1,NNODP(IPROC)
                     INDX=IMAP_NOD_LG(I,IPROC)
                     READ(LOCHSF,REC=IHOTSTP+1) SAL(INDX,K)
                     IHOTSTP=IHOTSTP+1
                   END DO
                 END DO
                 DO K=1,NFEN
                   DO I=1,NNODP(IPROC)
                     INDX=IMAP_NOD_LG(I,IPROC)
                     READ(LOCHSF,REC=IHOTSTP+1) TEMP(INDX,K)
                     IHOTSTP=IHOTSTP+1
                   END DO
                 END DO
               END IF
               CLOSE(LOCHSF)
            END DO
         END IF

c kmd 50 update to include harmonics in hotstart files

         IF ((IHARIND.EQ.1).AND.(ITHSF.GT.ITHAS)) THEN 
            DO IPROC=1,NPROC
               IF(IHOT.EQ.67) LOCFN(1:14) = 'PE0000/fort.67'
               IF(IHOT.EQ.68) LOCFN(1:14) = 'PE0000/fort.68'
               LOCHSF=LOC2(IPROC)
               WRITE(LOCFN(3:6),'(I4.4)') (IPROC-1)
               OPEN (LOCHSF,FILE=LOCFN,ACCESS='DIRECT',RECL=8)
               IF ((IM.EQ.10).AND.((IM.EQ.1).OR.(IM.EQ.2).OR. 
     &             (IM.EQ.11).OR.(IM.EQ.21).OR.(IM.EQ.31))) THEN
                  IHOTSTP=8+(NNODP(IPROC)*7)+(NELP(IPROC)*1)+31+
     &                    (NNODP(IPROC)*7)+(NNODP(IPROC)*NFEN*5)
               ELSE IF ((IM.EQ.10).AND.((IM.NE.1).OR.(IM.NE.2).OR.
     &                  (IM.NE.11).OR.(IM.NE.21).OR.(IM.NE.31))) THEN 
                  IHOTSTP=8+(NNODP(IPROC)*7)+(NELP(IPROC)*1)+31
               ELSE IF ((IM.NE.10).AND.((IM.EQ.1).OR.(IM.EQ.2).OR. 
     &             (IM.EQ.11).OR.(IM.EQ.21).OR.(IM.EQ.31))) THEN
                  IHOTSTP=8+(NNODP(IPROC)*6)+(NELP(IPROC)*1)+31+
     &                    (NNODP(IPROC)*7)+(NNODP(IPROC)*NFEN*5)
               ELSE IF ((IM.NE.10).AND.((IM.NE.1).OR.(IM.NE.2).OR.
     &                  (IM.NE.11).OR.(IM.NE.21).OR.(IM.NE.31))) THEN 
                  IHOTSTP=8+(NNODP(IPROC)*6)+(NELP(IPROC)*1)+31
               END IF
               IF ((IDEN.EQ.1).OR.(IDEN.EQ.2).OR.
     &             (IDEN.EQ.3)) THEN 
                  IHOTSTP=IHOTSTP+(NNODP(IPROC)*NFEN*1)
               ELSE IF (IDEN.EQ.4) THEN 
                  IHOTSTP=IHOTSTP+(NNODP(IPROC)*NFEN*2) 
               END IF 
               IHOTSTP=IHOTSTP+1
               READ(LOCHSF,REC=IHOTSTP) ICHA
               READ(LOCHSF,REC=IHOTSTP+1) NZ 
               READ(LOCHSF,REC=IHOTSTP+2) NF 
               READ(LOCHSF,REC=IHOTSTP+3) MM 
               READ(LOCHSF,REC=IHOTSTP+4) NPX 
               READ(LOCHSF,REC=IHOTSTP+5) NSTAEX 
               READ(LOCHSF,REC=IHOTSTP+6) NSTAVX
               READ(LOCHSF,REC=IHOTSTP+7) NHASE
               READ(LOCHSF,REC=IHOTSTP+8) NHASV 
               READ(LOCHSF,REC=IHOTSTP+9) NHAGE
               READ(LOCHSF,REC=IHOTSTP+10) NHAGV 
               READ(LOCHSF,REC=IHOTSTP+11) ICALL 
               READ(LOCHSF,REC=IHOTSTP+12) NFREQ 
               IHOTSTP = IHOTSTP + 12
            
               DO I=1,NFREQ+NF
                  FNAME=NAMEFR(I) 
                  READ(LOCHSF,REC=IHOTSTP+1) FNAM8(1)
                  READ(LOCHSF,REC=IHOTSTP+2) FNAM8(2)
                  IHOTSTP = IHOTSTP + 2 
                  READ(LOCHSF,REC=IHOTSTP+1) HAFREQ(I)
                  READ(LOCHSF,REC=IHOTSTP+2) HAFF(I)
                  READ(LOCHSF,REC=IHOTSTP+3) HAFACE(I)
                  IHOTSTP=IHOTSTP+3
               END DO 

C --       Write out time of most recent H.A. update

               READ(LOCHSF,REC=IHOTSTP+1) TIMEUD
               READ(LOCHSF,REC=IHOTSTP+2) ITUD
               IHOTSTP=IHOTSTP+2
 
C --       Write out LHS Matrix 

               IHOTSTP=IHOTSTP+1
               DO I=1,MM
                  DO K=1,MM
                     READ(LOCHSF,REC=IHOTSTP) HA(I,K)
                     IHOTSTP=IHOTSTP+1 
                  END DO 
               END DO 

C --      Read in load vectors to hotstart file.
               IF (NHASE.NE.0) THEN 
                 DO N=1,NSTAEX 
                    DO I=1,MM
                       READ(LOCHSF,REC=IHOTSTP) STAELV(I,N)
                       IHOTSTP=IHOTSTP+1
                    END DO 
                 END DO 
               END IF 

               IF (NHASV.NE.0) THEN 
                 DO N=1,NSTAVX 
                    DO I=1,MM
                       READ(LOCHSF,REC=IHOTSTP) STAULV(I,N)
                       IHOTSTP=IHOTSTP+1
                       READ(LOCHSF,REC=IHOTSTP) STAVLV(I,N)
                       IHOTSTP=IHOTSTP+1
                    END DO 
                 END DO 
               END IF 

               IF (NHAGE.NE.0) THEN 
                 DO N=1,NNODP(IPROC)
                    DO I=1,MM
                       INDX=IMAP_NOD_LG(I,IPROC) 
                       READ(LOCHSF,REC=IHOTSTP) GLOELV(I,INDX)
                       IHOTSTP=IHOTSTP+1
                    END DO 
                 END DO 
               END IF 

               IF (NHAGV.NE.0) THEN 
                 DO N=1,NNODP(IPROC) 
                    DO I=1,MM
                       INDX=IMAP_NOD_LG(I,IPROC)
                       READ(LOCHSF,REC=IHOTSTP) GLOULV(I,INDX)
                       IHOTSTP=IHOTSTP+1
                       READ(LOCHSF,REC=IHOTSTP) GLOVLV(I,INDX)
                       IHOTSTP=IHOTSTP+1
                    END DO 
                 END DO 
               END IF 

               IF (CHARMV) THEN 
                 IF ((IHARIND.EQ.1).AND.(ITHSF.GT.ITMV)) THEN 
                  READ(LOCHSF,REC=IHOTSTP) NTSTEPS
                  IHOTSTP=IHOTSTP+1
                  IF (NHAGE.NE.0) THEN 
                     DO I=1,NNODP(IPROC)
                        INDX=IMAP_NOD_LG(I,IPROC)
                        READ(LOCHSF,REC=IHOTSTP) ELAV(INDX)
                        IHOTSTP=IHOTSTP+1
                        READ(LOCHSF,REC=IHOTSTP) ELVA(INDX)
                        IHOTSTP=IHOTSTP+1
                     END DO 
                  END IF 
                  IF (NHAGV.NE.0) THEN 
                     DO I=1,NNODP(IPROC)
                        INDX=IMAP_NOD_LG(I,IPROC)
                        READ(LOCHSF,REC=IHOTSTP) XELAV(INDX)
                        IHOTSTP=IHOTSTP+1
                        READ(LOCHSF,REC=IHOTSTP) YELAV(INDX)
                        IHOTSTP=IHOTSTP+1
                        READ(LOCHSF,REC=IHOTSTP) XELVA(INDX)
                        IHOTSTP=IHOTSTP+1
                        READ(LOCHSF,REC=IHOTSTP) YELVA(INDX)
                        IHOTSTP=IHOTSTP+1
                     END DO 
                  END IF 
                 END IF 
               END IF 
               CLOSE(LOCHSF)
            END DO 
         END IF 


C--Open Appropriate Hot Start File based on the value of IHOT from the fort.15 file 

        IF(IHOT.EQ.67) FNAME='fort.67'
        IF(IHOT.EQ.68) FNAME='fort.68'
        OPEN(IHOT,FILE=FNAME,ACCESS='DIRECT',RECL=8)

C--Write out info to global file
C  kmd48.33bc update to new format

        IHOTSTP=1
        WRITE(IHOT,REC=IHOTSTP) InputFileFmtVn
        IHOTSTP=2
        WRITE(IHOT,REC=IHOTSTP) IMHSF
        IHOTSTP=3
        WRITE(IHOT,REC=IHOTSTP) TIMEHSF
        IHOTSTP=4
        WRITE(IHOT,REC=IHOTSTP) ITHSF
        IHOTSTP=5
        WRITE(IHOT,REC=IHOTSTP) NPG
        IHOTSTP=6
        WRITE(IHOT,REC=IHOTSTP) NEG
        IHOTSTP=7
        WRITE(IHOT,REC=IHOTSTP) NPX
        IHOTSTP=8
        WRITE(IHOT,REC=IHOTSTP) NEX
        DO I=1,MNP
          WRITE(IHOT,REC=IHOTSTP+1) ETA1(I)
          IHOTSTP=IHOTSTP+1
        END DO
        DO I=1,MNP
          WRITE(IHOT,REC=IHOTSTP+1) ETA2(I)
          IHOTSTP=IHOTSTP+1
        END DO
        DO I=1,MNP
          WRITE(IHOT,REC=IHOTSTP+1) EtaDisc(I)
          IHOTSTP=IHOTSTP+1
        END DO
        DO I=1,MNP
          WRITE(IHOT,REC=IHOTSTP+1) UU2(I)
          IHOTSTP=IHOTSTP+1
        END DO
        DO I=1,MNP
          WRITE(IHOT,REC=IHOTSTP+1) VV2(I)
          IHOTSTP=IHOTSTP+1
        END DO
        IF(IM.EQ.10) THEN
          DO I=1,MNP
            WRITE(IHOT,REC=IHOTSTP+1) CH1(I)
            IHOTSTP=IHOTSTP+1
          END DO
        END IF
        DO I=1,MNP
          WRITE(IHOT,REC=IHOTSTP+1) NODECODE(I)
          IHOTSTP=IHOTSTP+1
        END DO
        DO I=1,MNE
           WRITE(IHOT,REC=IHOTSTP+1) NOFF(I)
           IHOTSTP=IHOTSTP+1
        END DO
        WRITE(IHOT,REC=IHOTSTP+1 ) IESTP
        WRITE(IHOT,REC=IHOTSTP+2 ) NSCOUE
        WRITE(IHOT,REC=IHOTSTP+3 ) IVSTP
        WRITE(IHOT,REC=IHOTSTP+4 ) NSCOUV
        WRITE(IHOT,REC=IHOTSTP+5 ) ICSTP
        WRITE(IHOT,REC=IHOTSTP+6 ) NSCOUC
        WRITE(IHOT,REC=IHOTSTP+7 ) IPSTP
        WRITE(IHOT,REC=IHOTSTP+8 ) IWSTP
        WRITE(IHOT,REC=IHOTSTP+9 ) NSCOUM
        WRITE(IHOT,REC=IHOTSTP+10) IGEP
        WRITE(IHOT,REC=IHOTSTP+11) NSCOUGE
        WRITE(IHOT,REC=IHOTSTP+12) IGVP
        WRITE(IHOT,REC=IHOTSTP+13) NSCOUGV
        WRITE(IHOT,REC=IHOTSTP+14) IGCP
        WRITE(IHOT,REC=IHOTSTP+15) NSCOUGC
        WRITE(IHOT,REC=IHOTSTP+16) IGPP
        WRITE(IHOT,REC=IHOTSTP+17) IGWP
        WRITE(IHOT,REC=IHOTSTP+18) NSCOUGW
        IF ((IM.EQ.1).OR.(IM.EQ.2).OR.(IM.EQ.11)
     &        .OR.(IM.EQ.21).OR.(IM.EQ.31)) THEN
           WRITE(IHOT,REC=IHOTSTP+19) IDEN
           WRITE(IHOT,REC=IHOTSTP+20) N3DSD
           WRITE(IHOT,REC=IHOTSTP+21) I3DSDRec
           WRITE(IHOT,REC=IHOTSTP+22) N3DSV
           WRITE(IHOT,REC=IHOTSTP+23) I3DSVRec
           WRITE(IHOT,REC=IHOTSTP+24) N3DST
           WRITE(IHOT,REC=IHOTSTP+25) I3DSTRec
           WRITE(IHOT,REC=IHOTSTP+26) N3DGD
           WRITE(IHOT,REC=IHOTSTP+27) I3DGDRec
           WRITE(IHOT,REC=IHOTSTP+28) N3DGV
           WRITE(IHOT,REC=IHOTSTP+29) I3DGVRec
           WRITE(IHOT,REC=IHOTSTP+30) N3DGT
           WRITE(IHOT,REC=IHOTSTP+31) I3DGTRec
           IHOTSTP=IHOTSTP+31
           DO I=1,MNP
             WRITE(IHOT,REC=IHOTSTP+1) DUU(I)
             IHOTSTP=IHOTSTP+1
           END DO
           DO I=1,MNP
             WRITE(IHOT,REC=IHOTSTP+1) DUV(I)
             IHOTSTP=IHOTSTP+1
           END DO
           DO I=1,MNP
             WRITE(IHOT,REC=IHOTSTP+1) DVV(I)
             IHOTSTP=IHOTSTP+1
           END DO
           DO I=1,MNP
             WRITE(IHOT,REC=IHOTSTP+1) UU(I)
             IHOTSTP=IHOTSTP+1
           END DO
           DO I=1,MNP
             WRITE(IHOT,REC=IHOTSTP+1) VV(I)
             IHOTSTP=IHOTSTP+1
           END DO
           DO I=1,MNP
             WRITE(IHOT,REC=IHOTSTP+1) BSX(I)
             IHOTSTP=IHOTSTP+1
           END DO
           DO I=1,MNP
             WRITE(IHOT,REC=IHOTSTP+1) BSY(I)
             IHOTSTP=IHOTSTP+1
           END DO
           DO K=1,NFEN
             DO I=1,MNP
               WRITE(IHOT,REC=IHOTSTP+1) REALQ(I,K)
               IHOTSTP=IHOTSTP+1
             END DO
           END DO
           DO K=1,NFEN
             DO I=1,MNP
               WRITE(IHOT,REC=IHOTSTP+1) IMAGQ(I,K)
               IHOTSTP=IHOTSTP+1
             END DO
           END DO
           DO K=1,NFEN
             DO I=1,MNP
               WRITE(IHOT,REC=IHOTSTP+1) WZ(I,K)
               IHOTSTP=IHOTSTP+1
             END DO
           END DO
           DO K=1,NFEN
             DO I=1,MNP
               WRITE(IHOT,REC=IHOTSTP+1) q20(I,K)
               IHOTSTP=IHOTSTP+1
             END DO
           END DO
           DO K=1,NFEN
             DO I=1,MNP
               WRITE(IHOT,REC=IHOTSTP+1) l(I,K)
               IHOTSTP=IHOTSTP+1
             END DO
           END DO
           IF (ABS(IDEN).EQ.1) THEN
             DO K=1,NFEN
               DO I=1,MNP
                 WRITE(IHOT,REC=IHOTSTP+1) SIGT(I,K)
                 IHOTSTP=IHOTSTP+1
               END DO
             END DO
           ELSE IF (ABS(IDEN).EQ.2) THEN
             DO K=1,NFEN
               DO I=1,MNP
                  WRITE(IHOT,REC=IHOTSTP+1) SAL(I,K)
                  IHOTSTP=IHOTSTP+1
               END DO
             END DO
           ELSE IF (ABS(IDEN).EQ.3) THEN
             DO K=1,NFEN
               DO I=1,MNP
                  WRITE(IHOT,REC=IHOTSTP+1) TEMP(I,K)
                  IHOTSTP=IHOTSTP+1
               END DO
             END DO
           ELSE IF (ABS(IDEN).EQ.4) THEN
             DO K=1,NFEN
               DO I=1,MNP
                  WRITE(IHOT,REC=IHOTSTP+1) SAL(I,K)
                  IHOTSTP=IHOTSTP+1
               END DO
             END DO
             DO K=1,NFEN
               DO I=1,MNP
                  WRITE(IHOT,REC=IHOTSTP+1) TEMP(I,K)
                  IHOTSTP=IHOTSTP+1
               END DO
             END DO
           END IF
        END IF

        CLOSE(IHOT)

        ENDDO

!        IHOTSTP=1
!       WRITE(IHOT,REC=IHOTSTP) IMHSF
!       IHOTSTP=2
!       WRITE(IHOT,REC=IHOTSTP) TIMEHSF
!       IHOTSTP=3
!       WRITE(IHOT,REC=IHOTSTP) ITHSF
!       DO I=1,MNP
!         WRITE(IHOT,REC=IHOTSTP+1) ETA1(I)
!         WRITE(IHOT,REC=IHOTSTP+2) ETA2(I)
!         WRITE(IHOT,REC=IHOTSTP+3) UU2(I)
!         WRITE(IHOT,REC=IHOTSTP+4) VV2(I)
!         IHOTSTP = IHOTSTP + 4
!         IF(IM.EQ.10) THEN
!           WRITE(IHOT,REC=IHOTSTP+1) CH1(I)
!           IHOTSTP=IHOTSTP+1
!           ENDIF
!         WRITE(IHOT,REC=IHOTSTP+1) NODECODE(I)
!         IHOTSTP=IHOTSTP+1
!         END DO
!       DO I=1,MNE
!          WRITE(IHOT,REC=IHOTSTP+1) NOFF(I)                     
!          IHOTSTP=IHOTSTP+1                                    
!       END DO
!       WRITE(IHOT,REC=IHOTSTP+1 ) IESTP
!       WRITE(IHOT,REC=IHOTSTP+2 ) NSCOUE
!       WRITE(IHOT,REC=IHOTSTP+3 ) IVSTP
!       WRITE(IHOT,REC=IHOTSTP+4 ) NSCOUV
!       WRITE(IHOT,REC=IHOTSTP+5 ) ICSTP
!       WRITE(IHOT,REC=IHOTSTP+6 ) NSCOUC
!       WRITE(IHOT,REC=IHOTSTP+7 ) IPSTP
!       WRITE(IHOT,REC=IHOTSTP+8 ) IWSTP
!       WRITE(IHOT,REC=IHOTSTP+9 ) NSCOUM
!       WRITE(IHOT,REC=IHOTSTP+10) IGEP
!       WRITE(IHOT,REC=IHOTSTP+11) NSCOUGE
!       WRITE(IHOT,REC=IHOTSTP+12) IGVP
!       WRITE(IHOT,REC=IHOTSTP+13) NSCOUGV
!       WRITE(IHOT,REC=IHOTSTP+14) IGCP
!       WRITE(IHOT,REC=IHOTSTP+15) NSCOUGC
!       WRITE(IHOT,REC=IHOTSTP+16) IGPP
!       WRITE(IHOT,REC=IHOTSTP+17) IGWP
!       WRITE(IHOT,REC=IHOTSTP+18) NSCOUGW
!
!       CLOSE(IHOT)           

!       ENDDO
C
      RETURN
      END

C---------------------------------------------------------------------------
C 
C                  S U B R O U T I N E   P O S T 4 1 
C
C---------------------------------------------------------------------------
C                     (  Serial Version  3/28/98  )                          
C  Globalize the 3D density data at the density stations from the local     
C  fort.41 files.                                                           
C  This version is compatible with ADCIRC version 34.03                     
C
C     jgf45.11 Updated: added support for different formats based on
C     value of IDEN, added support for 3D density stations defined by 
C     coordinates rather than node numbers.
C     
C---------------------------------------------------------------------------
      SUBROUTINE POST41()
      USE POST_GLOBAL
      IMPLICIT NONE
      INTEGER IRTYPE ! record type (1=elevation,2=velocity,3=3D velocity)
      INTEGER N      ! loop counter
      INTEGER NSTA   ! number of stations reported in subdomain fort.41 file
      INTEGER I,J,K,L,IPROC,IDUM
      INTEGER NDSET,NSTEMP,NP,IT,ITEMP,IG
      INTEGER OREC
      REAL(SZ) DT
      REAL(8) TIMEOUT
      CHARACTER*80  OUTMSG
      CHARACTER*85 INLINE
      CHARACTER*4 RDES4(8),RID4(6),AID4(6)
      CHARACTER*8 RDES8(8),RID8(6),AID8(6)
      LOGICAL FOUND

      INTEGER,ALLOCATABLE  :: LOC2(:),IREC(:)
      REAL(SZ),ALLOCATABLE  ::  SIGMA(:,:)
      REAL(SZ),ALLOCATABLE  ::  SDVar(:,:,:)  ! subdomain station data
      REAL(SZ),ALLOCATABLE  ::  FDvar(:,:,:)  ! full domain station data
      INTEGER NVar  ! number of density-related vars in the run
      INTEGER M     ! counter for NVar
      CHARACTER*14,ALLOCATABLE  ::  LOCNAME(:)

      ALLOCATE (LOC2(MNPROC),IREC(MNPROC))
      ALLOCATE (LOCNAME(MNPROC))

C
C--Determine whether Unit 41 is Sequential Formatted or Direct Access Binary
C
      IF (ABS(I3DSD).EQ.1) THEN
        GO TO 1000
      ELSE
        GO TO 2000
      ENDIF
C
1000  CONTINUE
C
C     Open all local sequential formatted fort.41 files
C
      DO IPROC = 1,NPROC
         LOCNAME(IPROC) = 'PE0000/fort.41'
         WRITE(LOCNAME(IPROC)(3:6),'(I4.4)') (IPROC-1)
         LOC2(IPROC) = 105 + (IPROC-1)
         INQUIRE(FILE=LOCNAME(IPROC),EXIST=FOUND)
         IF (FOUND) THEN
           OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC))
         ELSE
           print *, "No Local fort.41 files found."
           RETURN
         ENDIF
         READ(LOC2(IPROC),'(A85)')  INLINE
         READ(LOC2(IPROC),3645) NDSET,NSTA,DT,NSPO3DSD,NFEN,IRTYPE
         IF (NSTA.NE.NNSTA3DDP(IPROC)) THEN
            STOP 'Inconsistency in number of 3D density stations.'
         ENDIF
      ENDDO
C
C     Open global sequential formatted fort.41 file
C
      OPEN(UNIT=41,FILE='fort.41')
C
      WRITE(41,'(A85)') INLINE
      WRITE(41,3645) NDSET,NSTA3DD,DT,NSPO3DSD,NFEN,IRTYPE
C
      SELECT CASE(ABS(IDEN)) ! (+) is prognostic, (-) is diagnostic

      CASE(1) ! baroclinic with density forcing
         NVar=1
      CASE(2,3) ! baroclinic with salinity (=2) or temperature (=3) forcing
         NVar=2
      CASE(4) ! baroclinic with salinity and temperature forcing
         NVar=3
      CASE DEFAULT  ! fall-through
         PRINT *, 'IDEN = ',IDEN
         STOP 'ERROR: This value of IDEN is not allowed.'

      END SELECT

      ALLOCATE(SIGMA(NVar,NFEN),SDVar(NSTA3DD,NVar,NFEN),
     &   FDVar(NSTA3DD,NVar,NFEN))
C
C     jgf46.00 Made the loop over data sets an infinite loop because the
C     number of data sets at the top of the file will be wrong if the
C     run was hot started.
      DO                     ! loop infinitely over data sets
         DO IPROC=1,NPROC
            READ(LOC2(IPROC),FMT=1100,END=9999) TIMEOUT,IT,
     &           ((SIGMA(M,N),M=1,NVar),N=1,NFEN-1),
     &           SIGMA(1,NFEN),SIGMA(2,NFEN)
            IF (NNSTA3DDP(IPROC).GT.0) THEN
               DO K = 1,NNSTA3DDP(IPROC)
                  READ(LOC2(IPROC),1104,END=9999) IDUM,
     &               ((SDVar(K,M,L),M=1,NVar),L=1,NFEN)
                  IG = IMAP_STA3DD_LG(K,IPROC) ! global station number
                  DO M=1,NVar
                     DO L=1,NFEN
                        FDVar(IG,M,L)=SDVar(K,M,L)
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
C        jgf46.32jgf11 The writing of the data to the full domain 
C        output file must be performed after all the stations have
C        been read and mapped for a particular data set. This is 
C        because stations that land in ghost elements will appear in
C        multiple subdomains, and we have to go through and map them 
C        all into the full domain, allowing duplicates to overwrite 
C        each other (with identical values) in the full domain array.
         WRITE(41,1100) TIMEOUT,IT,
     &        ((SIGMA(M,N),M=1,NVar),N=1,NFEN-1),
     &        SIGMA(1,NFEN),SIGMA(2,NFEN)
         DO I=1,NSTA3DD
             WRITE(41,1104) I,((FDVar(I,M,L),M=1,NVar),L=1,NFEN)
         ENDDO
      ENDDO
      DEALLOCATE(SIGMA,SDVar,FDVar)

      GO TO 9999
C 
2000  CONTINUE
C
C     Open all local and the global direct access binary fort.41 file
C
C !!! THIS NEEDS TO BE WRITTEN !!!
      PRINT *,'I3DSD = ',I3DSD
      PRINT *,'This value of I3DSD is not currently supported.'
C
C     Close the global and local fort.61 files
C
9999  CONTINUE
      CLOSE(41)
      DO IPROC = 1,NPROC
         CLOSE(LOC2(IPROC))
      ENDDO
C
   80 FORMAT(A40)
 1100 FORMAT(1X,E16.10,1X,I10,32000(2X,E12.6))
 1104 FORMAT(9X,I6,4X,32000(E12.6,2X))
 3645 FORMAT(1X,I10,1X,I10,1X,E15.7,1X,I10,1X,I10,1X,I10)
C
      RETURN
      END
C---------------------------------------------------------------------------
C     End of subroutine post41
C---------------------------------------------------------------------------



C---------------------------------------------------------------------------
C 
C                  S U B R O U T I N E   P O S T 4 2 
C
C---------------------------------------------------------------------------
C
C                     (  Serial Version  3/28/98  )                          
C  Globalize the 3D velocity data at the velocity stations from the local   
C  fort.42 files.                                                           
C  This version is compatible with ADCIRC version 34.03                     
C
C     jgf45.11 Updated to process 3D velocity recording stations defined
C     by coordinates rather than node numbers.
C
C---------------------------------------------------------------------------
      SUBROUTINE POST42()
      USE POST_GLOBAL
      IMPLICIT NONE
      INTEGER IRTYPE ! record type (1=elevation,2=velocity,3=3D velocity)
      INTEGER N      ! loop counter
      INTEGER NSTA   ! number of stations reported in subdomain fort.41 file
      INTEGER I,J,K,L,IPROC,IDUM
      INTEGER NDSET,NSTEMP,NP,IT,ITEMP,IG
      INTEGER OREC
      REAL(SZ) DT
      REAL(8) TIMEOUT
      CHARACTER*80  OUTMSG
      CHARACTER*85 INLINE
      CHARACTER*4 RDES4(8),RID4(6),AID4(6)
      CHARACTER*8 RDES8(8),RID8(6),AID8(6)
      LOGICAL FOUND

      INTEGER,ALLOCATABLE  :: LOC2(:),IREC(:)
      REAL(SZ),ALLOCATABLE  ::  SIGMA(:)
      REAL(SZ),ALLOCATABLE  ::  SDVar(:,:,:)  ! subdomain station data
      REAL(SZ),ALLOCATABLE  ::  FDvar(:,:,:)  ! full domain station data
      INTEGER NVar  ! number of density-related vars in the run
      INTEGER M     ! counter for NVar
      CHARACTER*14,ALLOCATABLE  ::  LOCNAME(:)

      ALLOCATE (LOC2(MNPROC),IREC(MNPROC))
      ALLOCATE (LOCNAME(MNPROC))
C
C     Determine whether unit 42 is sequential formatted or direct access binary
C
      IF (ABS(I3DSV).EQ.1) THEN
        GO TO 1000
      ELSE
        GO TO 2000
      ENDIF
C
1000  CONTINUE
C
C     Open all local sequential formatted fort.42 files
C
      DO IPROC = 1,NPROC
         LOCNAME(IPROC) = 'PE0000/fort.42'
         WRITE(LOCNAME(IPROC)(3:6),'(I4.4)') (IPROC-1)
         LOC2(IPROC) = 105 + (IPROC-1)
         INQUIRE(FILE=LOCNAME(IPROC),EXIST=FOUND)
         IF (FOUND) THEN
           OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC))
         ELSE
           print *, "No local fort.42 files found."
           RETURN
         ENDIF
         READ(LOC2(IPROC),'(A85)')  INLINE
         READ(LOC2(IPROC),3645) NDSET,NSTA,DT,NSPO3DSV,NFEN,IRTYPE
         IF (NSTA.NE.NNSTA3DVP(IPROC)) THEN
            STOP 'Inconsistency in number of 3D velocity stations.'
         ENDIF
      ENDDO
C
C     Open global sequential formatted fort.42 file
C
      OPEN(UNIT=42,FILE='fort.42')
C
      NVar=3 ! x, y, and z velocities
      ALLOCATE(SIGMA(NFEN),SDVar(NSTA3DV,NVar,NFEN),
     &   FDVar(NSTA3DV,NVar,NFEN))
C
      WRITE(42,'(A85)') INLINE
      WRITE(42,3645) NDSET,NSTA3DV,DT,NSPO3DSV,NFEN,IRTYPE
C
C     jgf46.00 Made the loop over data sets an infinite loop because the
C     number of data sets at the top of the file will be wrong if the
C     run was hot started.
C
      DO                        ! loop infinitely over data sets
         DO IPROC=1,NPROC
            READ(LOC2(IPROC),FMT=1100,END=9999) TIMEOUT,IT,
     &           (SIGMA(N),SIGMA(N),SIGMA(N),N=1,NFEN-1),
     &           SIGMA(NFEN),SIGMA(NFEN)
            IF (NNSTA3DVP(IPROC).GT.0) THEN
               DO K = 1,NNSTA3DVP(IPROC)
                  READ(LOC2(IPROC),1104,END=9999) IDUM,
     &               ((SDVar(K,M,L),M=1,NVar),L=1,NFEN)
                  IG = IMAP_STA3DV_LG(K,IPROC)
                  DO M=1,NVar
                     DO L=1,NFEN
                        FDVar(IG,M,L)=SDVar(K,M,L)
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
C        jgf46.32jgf11 The writing of the data to the full domain 
C        output file must be performed after all the stations have
C        been read and mapped for a particular data set. This is 
C        because stations that land in ghost elements will appear in
C        multiple subdomains, and we have to go through and map them 
C        all into the full domain, allowing duplicates to overwrite 
C        each other (with identical values) in the full domain array.
         WRITE(42,1100) TIMEOUT,IT,
     &      (SIGMA(N),SIGMA(N),SIGMA(N),N=1,NFEN-1),
     &      SIGMA(NFEN),SIGMA(NFEN)
         DO I=1,NSTA3DV
            WRITE(42,1104) I,((FDVar(I,M,L),M=1,NVar),L=1,NFEN)
         ENDDO
C     
      ENDDO

      GO TO 9999
C 
2000  CONTINUE
C
C     Open all local and the global direct access binary fort.42 file
C
C !!! THIS NEEDS TO BE WRITTEN !!!
      PRINT *,'I3DSV = ',I3DSV
      PRINT *,'This value of I3DSV is not currently supported.'
C
C     Close the global and local fort.42 files
C
9999  CONTINUE
      CLOSE(42)
      DO IPROC = 1,NPROC
         CLOSE(LOC2(IPROC))
      ENDDO
C
   80 FORMAT(A40)
 1100 FORMAT(1X,E16.10,1X,I10,32000(2X,E12.6))
 1104 FORMAT(9X,I6,4X,32000(E12.6,2X))
 3645 FORMAT(1X,I10,1X,I10,1X,E15.7,1X,I10,1X,I10,1X,I10)
C
      RETURN
      END
C---------------------------------------------------------------------------
C     End of subroutine post42
C---------------------------------------------------------------------------


C---------------------------------------------------------------------------
C 
C                  S U B R O U T I N E   P O S T 4 3 
C
C---------------------------------------------------------------------------
C                     (  Serial Version  3/28/98  )                          
C  Globalize the 3D turbulence data at the turbulence stations from the     
C  local fort.43 files.                                                     
C  This version is compatible with ADCIRC version 34.03                     
C
C     jgf45.11 Updated to process 3D turbulence recording stations defined
C     by coordinates rather than node numbers.
C
C---------------------------------------------------------------------------
      SUBROUTINE POST43()
      USE POST_GLOBAL
      IMPLICIT NONE
      INTEGER IRTYPE ! record type (1=elevation,2=velocity,3=3D velocity)
      INTEGER N      ! loop counter
      INTEGER NSTA   ! number of stations reported in subdomain fort.41 file
      INTEGER I,J,K,L,IPROC,IDUM
      INTEGER NDSET,NSTEMP,NP,IT,ITEMP,IG
      INTEGER OREC
      REAL(SZ) DT
      REAL(8) TIMEOUT
      CHARACTER*80  OUTMSG
      CHARACTER*85 INLINE
      CHARACTER*4 RDES4(8),RID4(6),AID4(6)
      CHARACTER*8 RDES8(8),RID8(6),AID8(6)
      LOGICAL FOUND
      INTEGER,ALLOCATABLE  :: LOC2(:),IREC(:)
      REAL(SZ),ALLOCATABLE  ::  SIGMA(:)
      REAL(SZ),ALLOCATABLE  ::  SDVar(:,:,:)  ! subdomain station data
      REAL(SZ),ALLOCATABLE  ::  FDvar(:,:,:)  ! full domain station data
      INTEGER NVar  ! number of density-related vars in the run
      INTEGER M     ! counter for NVar
      CHARACTER*14,ALLOCATABLE  ::  LOCNAME(:)

      ALLOCATE (LOC2(MNPROC),IREC(MNPROC))
      ALLOCATE (LOCNAME(MNPROC))
C
C     Determine whether unit 43 is sequential formatted or direct access binary
C
      IF (ABS(I3DST).EQ.1) THEN
        GO TO 1000
      ELSE
        GO TO 2000
      ENDIF
C
1000  CONTINUE
C
C     Open all local sequential formatted fort.43 files
C
      DO IPROC = 1,NPROC
         LOCNAME(IPROC) = 'PE0000/fort.43'
         WRITE(LOCNAME(IPROC)(3:6),'(I4.4)') (IPROC-1)
         LOC2(IPROC) = 105 + (IPROC-1)
         INQUIRE(FILE=LOCNAME(IPROC),EXIST=FOUND)
         IF (FOUND) THEN
           OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC))
         ELSE
           print *, "No local fort.43 files found."
           RETURN
         ENDIF
         READ(LOC2(IPROC),'(A85)')  INLINE
         READ(LOC2(IPROC),3645) NDSET,NSTA,DT,NSPO3DST,NFEN,IRTYPE
         IF (NSTA.NE.NNSTA3DTP(IPROC)) THEN
            STOP 'Inconsistency in number of 3D turbulence stations.'
         ENDIF
      ENDDO
C
C     Open global sequential formatted fort.43 file
C
      OPEN(UNIT=43,FILE='fort.43')
C
      NVar=3 ! turbulent kinetic energy, mixing length, horiz. eddy viscosity
      ALLOCATE(SIGMA(NFEN),SDVar(NSTA3DV,NVar,NFEN),
     &   FDVar(NSTA3DV,NVar,NFEN))
C
      WRITE(43,'(A85)') INLINE
      WRITE(43,3645) NDSET,NSTA,DT,NSPO3DST,NFEN,IRTYPE
C
C     jgf46.00 Made the loop over data sets an infinite loop because the
C     number of data sets at the top of the file will be wrong if the
C     run was hot started.
C
      DO                        ! loop infinitely over data sets
         DO IPROC=1,NPROC
            READ(LOC2(IPROC),FMT=1100,END=9999) TIMEOUT,IT,
     &           (SIGMA(N),SIGMA(N),SIGMA(N),N=1,NFEN-1),
     &           SIGMA(NFEN),SIGMA(NFEN)
            IF (NNSTA3DTP(IPROC).GT.0) THEN
               DO K = 1,NNSTA3DTP(IPROC)
                  READ(LOC2(IPROC),1104,END=9999) IDUM,
     &               ((SDVar(K,M,L),M=1,NVar),L=1,NFEN)
                  IG = IMAP_STA3DT_LG(K,IPROC)
                  DO M=1,NVar
                     DO L=1,NFEN
                        FDVar(IG,M,L)=SDVar(K,M,L)
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
C        jgf46.32jgf11 The writing of the data to the full domain 
C        output file must be performed after all the stations have
C        been read and mapped for a particular data set. This is 
C        because stations that land in ghost elements will appear in
C        multiple subdomains, and we have to go through and map them 
C        all into the full domain, allowing duplicates to overwrite 
C        each other (with identical values) in the full domain array.
         WRITE(43,1100) TIMEOUT,IT,
     &      (SIGMA(N),SIGMA(N),SIGMA(N),N=1,NFEN-1),
     &      SIGMA(NFEN),SIGMA(NFEN)
         DO I=1,NSTA3DT
            WRITE(43,1104) I,((FDVar(I,M,L),M=1,NVar),L=1,NFEN)
         ENDDO
C     
      ENDDO

      GO TO 9999
C 
2000  CONTINUE
C
C     Open all local and the global direct access binary fort.43 file
C
C !!! THIS NEEDS TO BE WRITTEN !!!
      PRINT *,'I3DST = ',I3DST
      PRINT *,'This value of I3DST is not currently supported.'
C
C     Close the global and local fort.43 files
C
9999  CONTINUE
      CLOSE(43)
      DO IPROC = 1,NPROC
         CLOSE(LOC2(IPROC))
      ENDDO
C
   80 FORMAT(A40)
 1100 FORMAT(1X,E16.10,1X,I10,32000(2X,E12.6))
 1104 FORMAT(9X,I6,4X,32000(E12.6,2X))
 3645 FORMAT(1X,I10,1X,I10,1X,E15.7,1X,I10,1X,I10,1X,I10)
C
      RETURN
      END SUBROUTINE
C---------------------------------------------------------------------------
C     End of subroutine post43
C---------------------------------------------------------------------------



C---------------------------------------------------------------------------
C 
C                  S U B R O U T I N E   P O S T 4 4 
C
C---------------------------------------------------------------------------
C                     (  Serial Version  3/28/98  )                          
C  Globalize the Density Data at all nodes from the local fort.44 files.    
C  This version is compatible with ADCIRC version 41.11a                    
C
C     jgf45.11 Updated to process 3D velocity recording stations defined
C     by coordinates rather than node numbers.
C
C---------------------------------------------------------------------------
!kmd48.33bc - added in the post processing of the fort.47 file for the top 
!       temperature boundary conditions. Since it writes out with the fort.44
!       file, I just decided not to make a special subroutine for it. We may 
!       want to make the file a global write with the changing of the 3D to 
!       write globally. 
      SUBROUTINE POST44()
      USE POST_GLOBAL
      IMPLICIT NONE
      INTEGER I,J,K,L,IPROC,IDUM
      INTEGER NDSET,NSTEMP,NP,IT,ITEMP,INDX
      INTEGER OREC
      INTEGER N      ! loop counter
      INTEGER IRTYPE ! record type, elevation=1, velocity=2, 3Dvelocity=3
      REAL(SZ) DT
      REAL(8) TIMEOUT
      CHARACTER*80  OUTMSG
      CHARACTER*85 INLINE
      CHARACTER*4 RDES4(8),RID4(6),AID4(6)
      CHARACTER*8 RDES8(8),RID8(6),AID8(6)
      LOGICAL FOUND

      INTEGER,ALLOCATABLE  :: LOC2(:),LOC1(:),IREC(:)
      REAL(SZ),ALLOCATABLE  ::  SIGMA(:),VAR1(:,:),VAR2(:,:),VAR3(:,:)
      REAL(SZ),ALLOCATABLE  ::  VAR1P(:,:),VAR2P(:,:),VAR3P(:,:)
      CHARACTER*14,ALLOCATABLE  ::  LOCNAME(:),LOCNAME1(:)

      ALLOCATE (LOC2(MNPROC),IREC(MNPROC),LOC1(MNPROC))
      ALLOCATE (LOCNAME(MNPROC),LOCNAME1(MNPROC))
C
C     jgf45.12 When post processing a diagnostic run, it does not make
C     sense to post process the fort.44 file, since the density field
C     does not change over time.
      IF (IDEN .LT. 0) THEN
         WRITE(*,*) "IDEN was negative; this was a diagnostic run."
         WRITE(*,*) "The fort.44 files will not be post processed."
         !RETURN !Note the early return
      ENDIF

C
C--Determine whether Unit 44 is Sequential Formatted or Direct Access Binary
C
      IF (ABS(I3DGD).EQ.1) THEN
         GO TO 1000
      ELSE
         GO TO 2000
      ENDIF
C
1000  CONTINUE
C
C--Open Global fort.44 file
C
      DO IPROC = 1,NPROC
         LOCNAME(IPROC) = 'PE0000/fort.44'
         WRITE(LOCNAME(IPROC)(3:6),'(I4.4)') (IPROC-1)
         LOC2(IPROC) = 105 + (IPROC-1)
         INQUIRE(FILE=LOCNAME(IPROC),EXIST=FOUND)
         IF (FOUND) THEN
            OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC))
         ELSE
            print *, "No Local fort.44 files found"
            RETURN
         ENDIF
         READ(LOC2(IPROC),'(A85)')  INLINE
         READ(LOC2(IPROC),3645) NDSET,NP,DT,NSPO3DGD,NFEN,IRTYPE
         IF (NP.NE.NNODP(IPROC)) THEN
            STOP 'Inconsistency in number of local nodes.'
         ENDIF
      ENDDO
      OPEN(UNIT=44,FILE='fort.44')
C     
      WRITE(44,'(A85)') INLINE
      WRITE(44,3645) NDSET,NNODG,DT,NSPO3DGD,NFEN,IRTYPE
C     
      SELECT CASE(ABS(IDEN)) ! (+) is prognostic, (-) is diagnostic     

      CASE(1) ! baroclinic model run, density forcing
         ALLOCATE(SIGMA(NFEN),VAR1(MNP,NFEN),VAR1P(MNPP,NFEN))
C
C     jgf46.00 Made the loop over data sets an infinite loop because the
C     number of data sets at the top of the file will be wrong if the
C     run was hot started.
C
         DO                     ! loop infinitely over data sets
            DO IPROC=1,NPROC
               READ(LOC2(IPROC),FMT=1100,END=9999) TIMEOUT,IT,
     &              (SIGMA(N),N=1,NFEN-1)

               IF (IPROC.EQ.1) WRITE(44,1100) TIMEOUT,IT,
     &              (SIGMA(N),N=1,NFEN-1)
               DO K = 1,NNODP(IPROC)
                  READ(LOC2(IPROC),1104,END=9999) IDUM,(VAR1P(K,L),
     &                 L=1,NFEN)
                  INDX = IMAP_NOD_LG(K,IPROC)
                  IF(IMAP_NOD_GL(1,INDX).EQ.IPROC) THEN
                     DO L=1,NFEN
                        VAR1(INDX,L) = VAR1P(K,L)
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO 
            DO I = 1,NNODG
               WRITE(44,1104) I,(VAR1(I,L),L=1,NFEN)
            ENDDO
         ENDDO 
         DEALLOCATE(SIGMA,VAR1,VAR1P)

      CASE(2,3) ! baroclinic with salinity(=2) or temperature(=3) forcing
         ALLOCATE(SIGMA(NFEN),VAR1(MNP,NFEN),VAR1P(MNPP,NFEN),
     &        VAR2(MNP,NFEN),VAR2P(MNPP,NFEN))
         DO                     ! loop infinitely over data sets
            DO IPROC=1,NPROC
               READ(LOC2(IPROC),FMT=1100,END=9999) TIMEOUT,IT,
     &              (SIGMA(N),SIGMA(N),N=1,NFEN-1),SIGMA(NFEN)
               IF (IPROC.EQ.1) WRITE(44,1100) TIMEOUT,IT,
     &              (SIGMA(N),SIGMA(N),N=1,NFEN-1),SIGMA(NFEN)
               DO K = 1,NNODP(IPROC)
                  READ(LOC2(IPROC),1104,END=9999) IDUM,(VAR1P(K,L),
     &                 VAR2P(K,L),L=1,NFEN)
                  INDX = IMAP_NOD_LG(K,IPROC)
                  IF(IMAP_NOD_GL(1,INDX).EQ.IPROC) THEN
                     DO L=1,NFEN
                        VAR1(INDX,L) = VAR1P(K,L)
                        VAR2(INDX,L) = VAR2P(K,L)
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO 
            DO I = 1,NNODG
               WRITE(44,1104) I,(VAR1(I,L),VAR2(I,L),L=1,NFEN)
            ENDDO
         ENDDO 
         DEALLOCATE(SIGMA,VAR1,VAR1P,VAR2,VAR2P)

      CASE(4) ! baroclinic with salinity and temperature forcing
         ALLOCATE(SIGMA(NFEN),VAR1(MNP,NFEN),VAR1P(MNPP,NFEN),
     &        VAR2(MNP,NFEN),VAR2P(MNPP,NFEN),
     &        VAR3(MNP,NFEN),VAR3P(MNPP,NFEN))
         DO                     ! loop infinitely over data sets
            DO IPROC=1,NPROC
               READ(LOC2(IPROC),FMT=1100,END=9999) TIMEOUT,IT,
     &              (SIGMA(N),SIGMA(N),SIGMA(N),N=1,NFEN-1),
     &              SIGMA(NFEN),SIGMA(NFEN)
               IF (IPROC.EQ.1) WRITE(44,1100) TIMEOUT,IT,
     &              (SIGMA(N),SIGMA(N),SIGMA(N),N=1,NFEN-1),
     &              SIGMA(NFEN),SIGMA(NFEN)
               DO K = 1,NNODP(IPROC)
                  READ(LOC2(IPROC),1104,END=9999) IDUM,(VAR1P(K,L),
     &             VAR2P(K,L),VAR3P(K,L),L=1,NFEN)
                  INDX = IMAP_NOD_LG(K,IPROC)
                  IF(IMAP_NOD_GL(1,INDX).EQ.IPROC) THEN
                     DO L=1,NFEN
                        VAR1(INDX,L) = VAR1P(K,L)
                        VAR2(INDX,L) = VAR2P(K,L)
                        VAR3(INDX,L) = VAR3P(K,L)
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO 
            DO I = 1,NNODG
               WRITE(44,1104) I,(VAR1(I,L),VAR2(I,L),VAR3(I,L),L=1,NFEN)
            ENDDO
         ENDDO 
         DEALLOCATE(SIGMA,VAR1,VAR1P,VAR2,VAR2P,VAR3,VAR3P)

      CASE DEFAULT  ! fall-through -> user can re-enter menu selection
         PRINT *, 'The value of IDEN is not supported by ADCPOST.'
         STOP

      END SELECT

      GO TO 9999
C
2000  CONTINUE
C
C--Open All Local and the Global Direct Access Binary fort.44 file
C
C !!!THIS NEEDS TO BE WRITTEN!!!
      PRINT *,'I3DGD = ',I3DGD
      PRINT *,'This value of I3DGD is not currently supported.'
C
9999  CONTINUE
C
C--Close the Global and Local fort.44 Files
C
      CLOSE(44)
      DO IPROC = 1,NPROC
         CLOSE(LOC2(IPROC))
      ENDDO
C
   80 FORMAT(A80)

 1100 FORMAT(1X,E16.10,1X,I10,32000(2X,E12.6))
 1104 FORMAT(9X,I6,4X,32000(E12.6,2X))
C 3645 FORMAT(1X,I10,1X,E15.7,I10,1X,I10,1X,I10//)
 3645 FORMAT(1X,I10,1X,I10,1X,E15.7,1X,I10,1X,I10,1X,I10)
C
      RETURN
      END SUBROUTINE
C---------------------------------------------------------------------------
C     End of subroutine post44
C---------------------------------------------------------------------------



      SUBROUTINE POST45()
C
C---------------------------------------------------------------------------C
C                     (  Serial Version  3/28/98  )                         C 
C  Globalize the Velocity Data at all nodes from the local fort.45 files.   C
C  This version is compatible with ADCIRC version 41.11a                    C
C
C     jgf45.11 Updated to reflect new 3D file format.
C
C---------------------------------------------------------------------------C
C
      USE POST_GLOBAL
      IMPLICIT NONE
      INTEGER I,J,K,L,IPROC,IDUM
      INTEGER NDSET,NSTEMP,NP,IT,ITEMP,INDX
      INTEGER OREC
      INTEGER N      ! loop counter
      INTEGER IRTYPE ! record type, elevation=1, velocity=2, 3Dvelocity=3
      REAL(SZ) DT
      REAL(8) TIMEOUT
      CHARACTER*80  OUTMSG
      CHARACTER*85 INLINE
      CHARACTER*1000 SUPERINLINE
      CHARACTER*4 RDES4(8),RID4(6),AID4(6)
      CHARACTER*8 RDES8(8),RID8(6),AID8(6)
      LOGICAL FOUND

      INTEGER,ALLOCATABLE  :: LOC2(:),IREC(:)
      REAL(SZ),ALLOCATABLE  ::  SIGMA(:),VAR1(:,:),VAR2(:,:),VAR3(:,:)
      REAL(SZ),ALLOCATABLE  ::  VAR1P(:,:),VAR2P(:,:),VAR3P(:,:)
      CHARACTER*14,ALLOCATABLE  ::  LOCNAME(:)

      ALLOCATE (LOC2(MNPROC),IREC(MNPROC))
      ALLOCATE (LOCNAME(MNPROC))
C
C--Determine whether Unit 45 is Sequential Formatted or Direct Access Binary
C
      IF (ABS(I3DGV).EQ.1) THEN
        GO TO 1000
      ELSE
        GO TO 2000
      ENDIF
C
1000  CONTINUE
C
C--Open Global fort.45 file
C
      DO IPROC = 1,NPROC
         LOCNAME(IPROC) = 'PE0000/fort.45'
         WRITE(LOCNAME(IPROC)(3:6),'(I4.4)') (IPROC-1)
         LOC2(IPROC) = 105 + (IPROC-1)
         INQUIRE(FILE=LOCNAME(IPROC),EXIST=FOUND)
         IF (FOUND) THEN
           OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC))
         ELSE
           print *, "No Local fort.45 files found"
           RETURN
         ENDIF
         READ(LOC2(IPROC),'(A85)')  INLINE
         READ(LOC2(IPROC),3645) NDSET,NP,DT,NSPO3DGV,NFEN,IRTYPE
         IF (NP.NE.NNODP(IPROC)) THEN
           STOP 'Inconsistency in number of local nodes'
         ENDIF
      ENDDO
      OPEN(UNIT=45,ACTION='WRITE',FILE='fort.45')
C
      ALLOCATE(SIGMA(NFEN))
      ALLOCATE(VAR1(MNP,NFEN),VAR2(MNP,NFEN),VAR3(MNP,NFEN))
      ALLOCATE(VAR1P(MNPP,NFEN),VAR2P(MNPP,NFEN),VAR3P(MNPP,NFEN))
C
      WRITE(45,'(A85)') INLINE
      WRITE(45,3645) NDSET,NNODG,DT,NSPO3DGV,NFEN,IRTYPE
C
 1100 FORMAT(1X,E16.10,1X,I10,32000(2X,E12.6))
C
C     jgf46.00 Made the loop over data sets an infinite loop because the
C     number of data sets at the top of the file will be wrong if the
C     run was hot started.
C
      DO                        ! loop infinitely over data sets
         DO IPROC=1,NPROC
            READ(LOC2(IPROC),FMT=1100,END=9999) TIMEOUT,IT,
     &                    (SIGMA(N),SIGMA(N),SIGMA(N),N=1,NFEN-1),
     &                     SIGMA(NFEN),SIGMA(NFEN)
            IF (IPROC.EQ.1) WRITE(45,1100) TIMEOUT,IT,
     &                    (SIGMA(N),SIGMA(N),SIGMA(N),N=1,NFEN-1),
     &                     SIGMA(NFEN),SIGMA(NFEN)
            DO K = 1,NNODP(IPROC)
               READ(LOC2(IPROC),1104,END=9999) IDUM,(VAR1P(K,L),
     &                                   VAR2P(K,L),VAR3P(K,L),L=1,NFEN)
               INDX = IMAP_NOD_LG(K,IPROC)
               IF(IMAP_NOD_GL(1,INDX).EQ.IPROC) THEN
                 DO L=1,NFEN
                   VAR1(INDX,L) = VAR1P(K,L)
                   VAR2(INDX,L) = VAR2P(K,L)
                   VAR3(INDX,L) = VAR3P(K,L)
                 ENDDO
               ENDIF
            ENDDO
         ENDDO 
         DO I = 1,NNODG
            WRITE(45,1104) I,(VAR1(I,L),VAR2(I,L),VAR3(I,L),L=1,NFEN)
         ENDDO
      ENDDO               
      GO TO 9999
C
2000  CONTINUE
C
C--Open All Local and the Global Direct Access Binary fort.45 file
C
C !!!THIS NEEDS TO BE WRITTEN!!!
      PRINT *,'I3DGV = ',I3DGV
      PRINT *,'This value of I3DGV is not currently supported.'
C
9999  CONTINUE
C
C--Close the Global and Local fort.45 Files
C
      CLOSE(45)
      DO IPROC = 1,NPROC
         CLOSE(LOC2(IPROC))
      ENDDO
C
      DEALLOCATE(SIGMA)
      DEALLOCATE(VAR1,VAR2,VAR3)
      DEALLOCATE(VAR1P,VAR2P,VAR3P)
C
   80 FORMAT(A80)
 1104 FORMAT(9X,I6,4X,32000(E12.6,2X))
 3645 FORMAT(1X,I10,1X,I10,1X,E15.7,1X,I10,1X,I10,1X,I10)
C
      RETURN
      END


      SUBROUTINE POST46()

C---------------------------------------------------------------------------C
C                     (  Serial Version  3/28/98  )                         C 
C  Globalize the Turbulence Data at all nodes from the local fort.46 files. C
C  This version is compatible with ADCIRC version 41.11a                    C
C
C     jgf45.11 Updated to reflect new 3D file format.
C
C---------------------------------------------------------------------------C

      USE POST_GLOBAL
      IMPLICIT NONE
      INTEGER I,J,K,L,IPROC,IDUM
      INTEGER NDSET,NSTEMP,NP,IT,ITEMP,INDX
      INTEGER OREC
      INTEGER N      ! loop counter
      INTEGER IRTYPE ! record type, elevation=1, velocity=2, 3Dvelocity=3
      REAL(SZ) DT
      REAL(8) TIMEOUT
      CHARACTER*80  OUTMSG
      CHARACTER*85 INLINE
      CHARACTER*4 RDES4(8),RID4(6),AID4(6)
      CHARACTER*8 RDES8(8),RID8(6),AID8(6)
      LOGICAL FOUND
      INTEGER,ALLOCATABLE  :: LOC2(:),IREC(:)
      REAL(SZ),ALLOCATABLE  ::  SIGMA(:),VAR1(:,:),VAR2(:,:),VAR3(:,:)
      REAL(SZ),ALLOCATABLE  ::  VAR1P(:,:),VAR2P(:,:),VAR3P(:,:)
      CHARACTER*14,ALLOCATABLE  ::  LOCNAME(:)

      ALLOCATE (LOC2(MNPROC),IREC(MNPROC))
      ALLOCATE (LOCNAME(MNPROC))

C
C--Determine whether Unit 46 is Sequential Formatted or Direct Access Binary
C
      IF (ABS(I3DGT).EQ.1) THEN
        GO TO 1000
      ELSE
        GO TO 2000
      ENDIF
C
1000  CONTINUE
C
C--Open Global fort.46 file
C
      DO IPROC = 1,NPROC
         LOCNAME(IPROC) = 'PE0000/fort.46'
         WRITE(LOCNAME(IPROC)(3:6),'(I4.4)') (IPROC-1)
         LOC2(IPROC) = 105 + (IPROC-1)
         INQUIRE(FILE=LOCNAME(IPROC),EXIST=FOUND)
         IF (FOUND) THEN
           OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC))
         ELSE
           print *, "No local fort.46 files found."
           RETURN
         ENDIF
         READ(LOC2(IPROC),'(A85)')  INLINE
         READ(LOC2(IPROC),3646) NDSET,NP,DT,NSPO3DGT,NFEN,IRTYPE
         IF (NP.NE.NNODP(IPROC)) THEN
           STOP 'Inconsistency in number of local nodes.'
         ENDIF
      ENDDO
      OPEN(UNIT=46,FILE='fort.46')
C
      ALLOCATE(SIGMA(NFEN))
      ALLOCATE(VAR1(MNP,NFEN),VAR2(MNP,NFEN),VAR3(MNP,NFEN))
      ALLOCATE(VAR1P(MNPP,NFEN),VAR2P(MNPP,NFEN),VAR3P(MNPP,NFEN))
C
      WRITE(46,'(A85)') INLINE
      WRITE(46,3646) NDSET,NNODG,DT,NSPO3DGT,NFEN,IRTYPE
C
C     jgf46.00 Made the loop over data sets an infinite loop because the
C     number of data sets at the top of the file will be wrong if the
C     run was hot started.
C
      DO                        ! loop infinitely over data sets
         DO IPROC=1,NPROC
            READ(LOC2(IPROC),FMT=1100,END=9999) TIMEOUT,IT,
     &                    (SIGMA(N),SIGMA(N),SIGMA(N),N=1,NFEN-1),
     &                     SIGMA(NFEN),SIGMA(NFEN)
            IF (IPROC.EQ.1) WRITE(46,1100) TIMEOUT,IT,
     &                    (SIGMA(N),SIGMA(N),SIGMA(N),N=1,NFEN-1),
     &                     SIGMA(NFEN),SIGMA(NFEN)
            DO K = 1,NNODP(IPROC)
               READ(LOC2(IPROC),1104,END=9999) IDUM,(VAR1P(K,L),
     &                                   VAR2P(K,L),VAR3P(K,L),L=1,NFEN)
               INDX = IMAP_NOD_LG(K,IPROC)
               IF(IMAP_NOD_GL(1,INDX).EQ.IPROC) THEN
                 DO L=1,NFEN
                   VAR1(INDX,L) = VAR1P(K,L)
                   VAR2(INDX,L) = VAR2P(K,L)
                   VAR3(INDX,L) = VAR3P(K,L)
                 ENDDO
               ENDIF
            ENDDO
         ENDDO 
         DO I = 1,NNODG
            WRITE(46,1104) I,(VAR1(I,L),VAR2(I,L),VAR3(I,L),L=1,NFEN)
         ENDDO
      ENDDO               
      GO TO 9999
C
2000  CONTINUE
C
C--Open All Local and the Global Direct Access Binary fort.46 file
C
C !!!THIS NEEDS TO BE WRITTEN!!!
      PRINT *,'I3DGT = ',I3DGT
      PRINT *,'This value of I3DGT is not currently supported.'
C
9999  CONTINUE
C
C--Close the Global and Local fort.46 Files
C
      CLOSE(46)
      DO IPROC = 1,NPROC
         CLOSE(LOC2(IPROC))
      ENDDO
C
      DEALLOCATE(SIGMA)
      DEALLOCATE(VAR1,VAR2,VAR3)
      DEALLOCATE(VAR1P,VAR2P,VAR3P)
C
   80 FORMAT(A80)
 1100 FORMAT(1X,E16.10,1X,I10,32000(2X,E12.6))
 1104 FORMAT(9X,I6,4X,32000(E12.6,2X))
C 3646 FORMAT(1X,I10,1X,E15.7,I10,1X,I10,1X,I10//)
 3646 FORMAT(1X,I10,1X,I10,1X,E15.7,1X,I10,1X,I10,1X,I10)
C
      RETURN
      END

C---------------------------------------------------------------------------
C 
C                  S U B R O U T I N E   P O S T 4 7 
C
C---------------------------------------------------------------------------
C                     (  Serial Version  3/28/98  )                          
C  Globalize the Density Data at all nodes from the local fort.47 files.    
C  This version is compatible with ADCIRC version 41.11a                    
C
C     jgf45.11 Updated to process 3D velocity recording stations defined
C     by coordinates rather than node numbers.
C
C---------------------------------------------------------------------------
!kmd48.33bc - added in the post processing of the fort.47 file for the top 
!       temperature boundary conditions. 
      SUBROUTINE POST47()
      USE POST_GLOBAL
      IMPLICIT NONE
      INTEGER I,J,K,L,IPROC,IDUM
      INTEGER NDSET,NSTEMP,NP,ITE,ITEMP,INDX
      INTEGER OREC
      INTEGER N      ! loop counter
      INTEGER IRTYPE ! record type, elevation=1, velocity=2, 3Dvelocity=3
      REAL(SZ) DT
      REAL(8) TIMEOUTE 
      CHARACTER*80  OUTMSG
      CHARACTER*85 INLINE
      CHARACTER*4 RDES4(8),RID4(6),AID4(6)
      CHARACTER*8 RDES8(8),RID8(6),AID8(6)
      LOGICAL FOUND

      INTEGER,ALLOCATABLE  :: LOC2(:),IREC(:)
      CHARACTER*80, ALLOCATABLE  ::  TOPTEMP(:),TOPTEMPP(:)
      CHARACTER*14,ALLOCATABLE  ::  LOCNAME(:)

      ALLOCATE (LOC2(MNPROC),IREC(MNPROC))
      ALLOCATE (TOPTEMP(MNP),TOPTEMPP(MNPP))
      ALLOCATE (LOCNAME(MNPROC))

C
C--Determine whether Unit 47 is Sequential Formatted or Direct Access Binary
C
      IF (ABS(I3DGD).EQ.1) THEN
         GO TO 1000
      ELSE
         GO TO 2000
      ENDIF
C
1000  CONTINUE
C
C--Open Global fort.47 file
C
      DO IPROC = 1,NPROC
         LOCNAME(IPROC) = 'PE0000/fort.47'
         WRITE(LOCNAME(IPROC)(3:6),'(I4.4)') (IPROC-1)
         LOC2(IPROC) = 105 + (IPROC-1)
         INQUIRE(FILE=LOCNAME(IPROC),EXIST=FOUND)
         IF (FOUND) THEN
            OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC))
         ELSE
            print *, "No Local fort.47 files found"
            RETURN
         ENDIF
         READ(LOC2(IPROC),'(A85)')  INLINE
         READ(LOC2(IPROC),3645) NDSET,NP,DT,NSPO3DGD,IRTYPE
         IF (NP.NE.NNODP(IPROC)) THEN
            STOP 'Inconsistency in number of local nodes.'
         ENDIF
      ENDDO
      OPEN(UNIT=47,FILE='fort.47')
C     
      WRITE(47,'(A85)') INLINE
      WRITE(47,3645) NDSET,NNODG,DT,NSPO3DGD,IRTYPE
C     
      DO ! loop infinitely over data sets
C
         DO IPROC=1,NPROC
            READ(LOC2(IPROC),FMT=2110,END=9999) TIMEOUTE,ITE
            IF (IPROC.EQ.1) WRITE(47,2110) TIMEOUTE,ITE
            DO K = 1,NNODP(IPROC)
               READ(LOC2(IPROC),80,END=9999) TOPTEMPP(K)
               INDX = IMAP_NOD_LG(K,IPROC)
               IF (IMAP_NOD_GL(1,INDX).EQ.IPROC) THEN
                 TOPTEMP(INDX) = TOPTEMPP(K)
                ENDIF
            ENDDO  
         ENDDO 
         DO I = 1,NNODG
            CALL NEWINDEX(TOPTEMP(I),OUTMSG,I)
            WRITE(47,*) TRIM(OUTMSG) !jgf46.00 TRIM to cut down file size 
            !WRITE(47,2413) I, TOPTEMP(I)
         ENDDO
C              
      ENDDO
      GO TO 9999
C
2000  CONTINUE
C
C--Open All Local and the Global Direct Access Binary fort.44 file
C
C !!!THIS NEEDS TO BE WRITTEN!!!
      PRINT *,'I3DGD = ',I3DGD
      PRINT *,'This value of I3DGD is not currently supported.'
C
9999  CONTINUE
C
C--Close the Global and Local fort.47 Files
C
      CLOSE(47)
      DO IPROC = 1,NPROC
         CLOSE(LOC2(IPROC))
      ENDDO
C
   80 FORMAT(A80)

 2110 FORMAT(2X,E20.10,5X,I10)
 2413 FORMAT(2X,I8,2X,E15.8)
C 3645 FORMAT(1X,I10,1X,E15.7,I10,1X,I10,1X,I10//)
C  kmd48.33bc Updated in code to get rid of return character
C           Updated all 3645 line to get rid of character
C           and also the 3646 line used in fort.46
 3645 FORMAT(1X,I10,1X,I10,1X,E15.7,1X,I10,1X,I10)
C
      RETURN
      END SUBROUTINE
C---------------------------------------------------------------------------
C     End of subroutine post47
C---------------------------------------------------------------------------

      SUBROUTINE POST81()
C
C----------------------------------------------------------------------------C
C                     (  Serial Version  3/28/98  )                          C 
C  Globalize the concentratiom data at the concentration stations from the   C
C  local fort.81 files.                                                      C
C  This version is compatible with ADCIRC version 34.03                      C
C  Added by mc67565 on June 19 2018, for definition of variables, similar to C
C  POST71                                                                    C
C----------------------------------------------------------------------------C
C
      USE POST_GLOBAL
      IMPLICIT NONE
      INTEGER I,J,K,L,IPROC,IDUM  
      INTEGER NTRSPC,NSTEMP,ITSC,ITEMPC,NUMSTNS
      INTEGER OREC
      REAL(SZ) DTC
      REAL(8) TIMEOUTSC
      CHARACTER*80  OUTMSG
      CHARACTER*85  INLINE
      CHARACTER*4 RDES4(8),RID4(6),AID4(6)
      CHARACTER*8 RDES8(8),RID8(6),AID8(6)
      LOGICAL FOUND

      INTEGER,ALLOCATABLE :: LOC2(:),IREC(:)
      REAL(SZ),ALLOCATABLE ::   ConcBIN(:),ConcBINP(:)
      CHARACTER*80,ALLOCATABLE ::  ConcASC(:),ConcASCP(:)
      CHARACTER*14,ALLOCATABLE ::  LOCNAME(:)

      ALLOCATE ( LOC2(MNPROC),IREC(MNPROC) )
      ALLOCATE ( ConcBIN(MNSTAC),ConcBINP(MNSTAC))
      ALLOCATE ( ConcASC(MNSTAC),ConcASCP(MNSTAC))
      ALLOCATE ( LOCNAME(MNPROC))

C
C--Determine whether Unit 81 is Sequential Formatted or Direct Access Binary
C
      IF (ABS(NOUTC).EQ.1) THEN
        GO TO 1000
      ELSE
        GO TO 2000
      ENDIF
C
1000  CONTINUE
C
C--Open All Local Sequential Formatted fort.81 file
C
      DO IPROC = 1,NPROC
         LOCNAME(IPROC) = 'PE0000/fort.81'
         WRITE(LOCNAME(IPROC)(3:6),'(I4.4)') (IPROC-1)
         LOC2(IPROC) = 105 + (IPROC-1)
         INQUIRE(FILE=LOCNAME(IPROC),EXIST=FOUND)
         IF (FOUND) THEN
           OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC))
         ELSE
           print *, "No Local fort.81 files found"
           RETURN
         ENDIF
         READ(LOC2(IPROC),'(A85)')  INLINE
         READ(LOC2(IPROC),3645) NTRSPC,NUMSTNS,DTC,NSTEMP,ITEMPC
         IF (NUMSTNS.NE.NSTACP(IPROC)) THEN
            STOP 'Inconsistency in Number of concentration Stations'
         ENDIF
      ENDDO
C
C--Open Global Sequential Formatted fort.81 file
C
      OPEN(UNIT=81,FILE='fort.81')
C
      WRITE(81,'(A85)') INLINE
      WRITE(81,3645) NTRSPC,NSTAC,DTC,NSPOOLC,ITEMPC
C
C     jgf46.00 Made the loop over data sets an infinite loop because the
C     number of data sets at the top of the file will be wrong if the
C     run was hot started.
C
      DO ! loop infinitely over data sets
C
         DO IPROC=1,NPROC
            READ(LOC2(IPROC),2120,END=9999) TIMEOUTSC,ITSC
            IF (NSTACP(IPROC).GT.0) THEN
              READ(LOC2(IPROC),80,END=9999) (ConcASCP(K),K=1,NSTACP(IPROC))
              DO K=1,NSTACP(IPROC)
                 ConcASC(IMAP_STAC_LG(K,IPROC)) = ConcASCP(K)
              ENDDO
            ENDIF
         ENDDO
C
         WRITE(81,2120) TIMEOUTSC,ITSC
         DO I=1, NSTAC
            CALL NEWINDEX(ConcASC(I),OUTMSG,I)
            WRITE(81,*) TRIM(OUTMSG) !jgf46.00 TRIM to cut down file size 
         ENDDO
C
      ENDDO
      GO TO 9999
C 
2000  CONTINUE
C
C--Open All Local and the Global Direct Access Binary fort.81 file
C
      DO IPROC = 1,NPROC
         LOCNAME(IPROC) = 'PE0000/fort.81'
         WRITE(LOCNAME(IPROC)(3:6),'(I4.4)') (IPROC-1)
         LOC2(IPROC) = 105 + (IPROC-1)
         INQUIRE(FILE=LOCNAME(IPROC),EXIST=FOUND)
         IF (FOUND) THEN
           OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC),ACCESS='DIRECT',
     &       RECL=NBYTE)
         ELSE
           print *, "No Local fort.81 files found"
          RETURN
         ENDIF
      ENDDO
      OPEN(81,FILE='fort.81',ACCESS='DIRECT',RECL=NBYTE)
C
C
C--Read RUNDES RUNID and AGRID from each Local file
C
      DO IPROC = 1,NPROC
         IREC(IPROC) = 0
         IF (NBYTE.EQ.4) THEN
           DO I=1,8
              READ(LOC2(IPROC),REC=IREC(IPROC)+I) RDES4(I)
           ENDDO
           IREC(IPROC)=IREC(IPROC)+8
           DO I=1,6
              READ(LOC2(IPROC),REC=IREC(IPROC)+I) RID4(I)
           ENDDO
           IREC(IPROC)=IREC(IPROC)+6
           DO I=1,6
              READ(LOC2(IPROC),REC=IREC(IPROC)+I) AID4(I)
           ENDDO
           IREC(IPROC)=IREC(IPROC)+6
         ENDIF
         IF (NBYTE.EQ.8) THEN
           DO I=1,4
              READ(LOC2(IPROC),REC=IREC(IPROC)+I) RDES8(I)
           ENDDO
           IREC(IPROC)=IREC(IPROC)+4
           DO I=1,3
              READ(LOC2(IPROC),REC=IREC(IPROC)+I) RID8(I)
           ENDDO
           IREC(IPROC)=IREC(IPROC)+3
           DO I=1,3
              READ(LOC2(IPROC),REC=IREC(IPROC)+I) AID8(I)
           ENDDO
           IREC(IPROC)=IREC(IPROC)+3
         ENDIF
      ENDDO
C
C--Write RUNDES RUNID and AGRID to Global File
C
      OREC = 0
      IF (NBYTE.EQ.4) THEN
        DO I=1,8
           WRITE(81,REC=OREC+I) RDES4(I)
        ENDDO
        OREC=OREC+8
        DO I=1,6
           WRITE(81,REC=OREC+I) RID4(I)
        ENDDO
        OREC=OREC+6
        DO I=1,6
           WRITE(81,REC=OREC+I) AID4(I)
        ENDDO
        OREC=OREC+6
      ENDIF
      IF (NBYTE.EQ.8) THEN
        DO I=1,4
           WRITE(81,REC=OREC+I) RDES8(I)
        ENDDO
        OREC=OREC+4
        DO I=1,3
           WRITE(81,REC=OREC+I) RID8(I)
        ENDDO
        OREC=OREC+3
        DO I=1,3
           WRITE(81,REC=OREC+I) AID8(I)
        ENDDO
        OREC=OREC+3
      ENDIF
C
C--Read NTRSPC, NSTAC, DT*NSPOOLC from each Local file
C  and then close files to flush file buffers
C
      DO IPROC = 1,NPROC
         READ(LOC2(IPROC),REC=IREC(IPROC)+1) NTRSPC
         READ(LOC2(IPROC),REC=IREC(IPROC)+2) NUMSTNS
         READ(LOC2(IPROC),REC=IREC(IPROC)+3) DTC
         READ(LOC2(IPROC),REC=IREC(IPROC)+4) NSTEMP
         READ(LOC2(IPROC),REC=IREC(IPROC)+5) ITEMPC
         IREC(IPROC) = IREC(IPROC)+5
         CLOSE(LOC2(IPROC))    ! Flush the Write Buffer
         IF (NUMSTNS.NE.NSTACP(IPROC)) THEN
            STOP 'Inconsistency in Number of concentration Stations'
         ENDIF
         OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC),ACCESS='DIRECT',RECL=NBYTE)
      ENDDO
C
C--Write same info to Global file and close it to flush buffer also
C
      WRITE(81,REC=OREC+1) NTRSPC
      WRITE(81,REC=OREC+2) NSTAC
      WRITE(81,REC=OREC+3) DTC
      WRITE(81,REC=OREC+4) NSPOOLC
      WRITE(81,REC=OREC+5) ITEMPC
      OREC = OREC+5
      CLOSE(81)         ! Flush the Write Buffer
      OPEN(81,FILE='fort.81',ACCESS='DIRECT',RECL=NBYTE)
C
      DO J=1,NTRSPC
C
         DO IPROC=1,NPROC
            READ(LOC2(IPROC),REC=IREC(IPROC)+1) TIMEOUTSC
            READ(LOC2(IPROC),REC=IREC(IPROC)+2) ITSC
            IREC(IPROC) = IREC(IPROC) + 2
            IF (NSTACP(IPROC).GT.0) THEN
              DO K=1, NSTACP(IPROC)
                 READ(LOC2(IPROC),REC=IREC(IPROC)+K) ConcBINP(K)
              ENDDO
              DO K=1,NSTACP(IPROC)
                 ConcBIN(IMAP_STAC_LG(K,IPROC)) = ConcBINP(K)
              ENDDO
              IREC(IPROC) = IREC(IPROC) + NSTACP(IPROC)
            ENDIF
         ENDDO
C
         WRITE(81,REC=OREC+1) TIMEOUTSC
         WRITE(81,REC=OREC+2) ITSC
         OREC = OREC + 2
         DO K=1, NSTAC
            WRITE(81,REC=OREC+K) ConcBIN(K)
         ENDDO
         OREC = OREC + NSTAC
C
      ENDDO
C
C--Close the Global and Local fort.81 Files
C
9999  CONTINUE
      CLOSE(81)
      DO IPROC = 1,NPROC
         CLOSE(LOC2(IPROC))
      ENDDO
C
   80 FORMAT(A40)
 2120 FORMAT(2X,E20.10,5X,I10)
 3645 FORMAT(1X,I10,1X,I10,1X,E15.7,1X,I5,1X,I5)
C
      RETURN
      END

C------------------------------------------------------------------------
C          END POST81
C------------------------------------------------------------------------

      SUBROUTINE POST83()
C
C---------------------------------------------------------------------------C
C                     (  Serial Version  3/28/98  )                         C 
C  Globalize the concentration data at all nodes from the local fort.83     C
C  files.                                                                   C
C  This version is compatible with ADCIRC version 34.03                     C
C  Added by mc67565 on June 19 2018, for definition of variables, similar   C
C  to POST73                                                                C
C---------------------------------------------------------------------------C
C
      USE POST_GLOBAL
      IMPLICIT NONE
      INTEGER I,J,K,L,IPROC,IDUM
      INTEGER NDSETGC,NSTEMP,NP,ITGC,ITEMPC,INDX
      INTEGER OREC
      REAL(SZ) DTGC
      REAL(SZ) ETA
      REAL(8) TIMEOUTSC
      CHARACTER*80 OUTMSG
      CHARACTER*85 INLINE
      CHARACTER*4 RDES4(8),RID4(6),AID4(6)
      CHARACTER*8 RDES8(8),RID8(6),AID8(6)
      LOGICAL FOUND

      INTEGER, ALLOCATABLE  ::  LOC2(:),IREC(:)
      REAL(SZ), ALLOCATABLE  ::  ConcBIN(:),ConcBINP(:)
      CHARACTER*80, ALLOCATABLE  ::  ConcSC(:),ConcSCP(:)
      CHARACTER*14, ALLOCATABLE  ::  LOCNAME(:)

      ALLOCATE (LOC2(MNPROC),IREC(MNPROC))
      ALLOCATE (ConcBIN(MNP),ConcBINP(MNPP))
      ALLOCATE (ConcSC(MNP),ConcSCP(MNPP))
      ALLOCATE (LOCNAME(MNPROC))

C
C--Determine whether Unit 83 is Sequential Formatted or Direct Access Binary
C
      IF (ABS(NOUTGC).EQ.1) THEN
        GO TO 1000
      ELSE
        GO TO 2000
      ENDIF
C
1000  CONTINUE
C
C--Open Global Sequential Formatted fort.83 file
C
      DO IPROC = 1,NPROC
         LOCNAME(IPROC) = 'PE0000/fort.83'
         WRITE(LOCNAME(IPROC)(3:6),'(I4.4)') (IPROC-1)
         LOC2(IPROC) = 105 + (IPROC-1)
         INQUIRE(FILE=LOCNAME(IPROC),EXIST=FOUND)
         IF (FOUND) THEN
           OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC))
         ELSE
           print *, "No Local fort.83 files found"
           RETURN
         ENDIF
         READ(LOC2(IPROC),'(A85)')  INLINE
         READ(LOC2(IPROC),3645)  NDSETGC,NP,DTGC,NSTEMP,ITEMPC
         IF (NP.NE.NNODP(IPROC)) THEN
           print *, "NP = ",NP," NNODP = ",NNODP(IPROC)
           STOP 'Inconsistency in number of local nodes'
         ENDIF
      ENDDO
      OPEN(UNIT=83,FILE='fort.83')
C
      WRITE(83,'(A85)') INLINE
      WRITE(83,3645) NDSETGC,NNODG,DTGC,NSPOOLGC,ITEMPC
C
C     jgf46.00 Made the loop over data sets an infinite loop because the
C     number of data sets at the top of the file will be wrong if the
C     run was hot started.
C
      DO ! loop infinitely over data sets
C
         DO IPROC=1,NPROC
            READ(LOC2(IPROC),FMT=2120,END=9999) TIMEOUTSC,ITGC
            IF (IPROC.EQ.1) WRITE(83,2120) TIMEOUTSC,ITGC
            DO K = 1,NNODP(IPROC)
               READ(LOC2(IPROC),80,END=9999) ConcSCP(K)
               INDX = IMAP_NOD_LG(K,IPROC)
               IF (IMAP_NOD_GL(1,INDX).EQ.IPROC) THEN
                 ConcSC(INDX) = ConcSCP(K)
                ENDIF
            ENDDO
         ENDDO
c        DO I = 1,NNODG
c           READ(ConcSC(I),*) IDUM, ETA
c           WRITE(83,2453) I, ETA
c        ENDDO
         DO I = 1,NNODG
            CALL NEWINDEX(ConcSC(I),OUTMSG,I)
            WRITE(83,*) TRIM(OUTMSG) !jgf46.00 TRIM to cut down file size 
         ENDDO
C              
      ENDDO
      GO TO 9999
C
2000  CONTINUE
C
C--Open All Local and the Global Direct Access Binary fort.83 file
C
      DO IPROC = 1,NPROC
         LOCNAME(IPROC) = 'PE0000/fort.83'
         WRITE(LOCNAME(IPROC)(3:6),'(I4.4)') (IPROC-1)
         LOC2(IPROC) = 105 + (IPROC-1)
         INQUIRE(FILE=LOCNAME(IPROC),EXIST=FOUND)
         IF (FOUND) THEN
           OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC),ACCESS='DIRECT',
     &       RECL=NBYTE)
         ELSE
           print *, "No Local fort.83 files found"
           RETURN
         ENDIF
      ENDDO
      OPEN(83,FILE='fort.83',ACCESS='DIRECT',RECL=NBYTE)
C
C--Read RUNDES RUNID and AGRID from each Local file
C
      DO IPROC = 1,NPROC
         IREC(IPROC) = 0
         IF (NBYTE.EQ.4) THEN
           DO I=1,8
              READ(LOC2(IPROC),REC=IREC(IPROC)+I) RDES4(I)
           ENDDO
           IREC(IPROC)=IREC(IPROC)+8
           DO I=1,6
              READ(LOC2(IPROC),REC=IREC(IPROC)+I) RID4(I)
           ENDDO
           IREC(IPROC)=IREC(IPROC)+6
           DO I=1,6
              READ(LOC2(IPROC),REC=IREC(IPROC)+I) AID4(I)
           ENDDO
           IREC(IPROC)=IREC(IPROC)+6
         ENDIF
         IF (NBYTE.EQ.8) THEN
           DO I=1,4
              READ(LOC2(IPROC),REC=IREC(IPROC)+I) RDES8(I)
           ENDDO
           IREC(IPROC)=IREC(IPROC)+4
           DO I=1,3
              READ(LOC2(IPROC),REC=IREC(IPROC)+I) RID8(I)
           ENDDO
           IREC(IPROC)=IREC(IPROC)+3
           DO I=1,3
              READ(LOC2(IPROC),REC=IREC(IPROC)+I) AID8(I)
           ENDDO
           IREC(IPROC)=IREC(IPROC)+3
         ENDIF
      ENDDO
C
C--Write RUNDES RUNID and AGRID to Global File
C
      OREC = 0
      IF (NBYTE.EQ.4) THEN
        DO I=1,8
           WRITE(83,REC=OREC+I) RDES4(I)
        ENDDO
        OREC=OREC+8
        DO I=1,6
           WRITE(83,REC=OREC+I) RID4(I)
        ENDDO
        OREC=OREC+6
        DO I=1,6
           WRITE(83,REC=OREC+I) AID4(I)
        ENDDO
        OREC=OREC+6
      ENDIF
      IF (NBYTE.EQ.8) THEN
        DO I=1,4
           WRITE(83,REC=OREC+I) RDES8(I)
        ENDDO
        OREC=OREC+4
        DO I=1,3
           WRITE(83,REC=OREC+I) RID8(I)
        ENDDO
        OREC=OREC+3
        DO I=1,3
           WRITE(83,REC=OREC+I) AID8(I)
        ENDDO
        OREC=OREC+3
      ENDIF
C
C--Read NDSETGC, NNODP, DTGC, NSPOOLGC, ITEMPC from each Local file
C  and then close files to flush file buffers
C
      DO IPROC = 1,NPROC
         READ(LOC2(IPROC),REC=IREC(IPROC)+1) NDSETGC
         READ(LOC2(IPROC),REC=IREC(IPROC)+2) IDUM
         READ(LOC2(IPROC),REC=IREC(IPROC)+3) DTGC
         READ(LOC2(IPROC),REC=IREC(IPROC)+4) NSTEMP
         READ(LOC2(IPROC),REC=IREC(IPROC)+5) ITEMPC
         IREC(IPROC) = IREC(IPROC)+5
         CLOSE(LOC2(IPROC))    ! Flush the Write Buffer
         OPEN(LOC2(IPROC),FILE=LOCNAME(IPROC),ACCESS='DIRECT',RECL=NBYTE)
      ENDDO
C
C--Write same info to Global file and close it to flush buffer also
C
      WRITE(83,REC=OREC+1) NDSETGC
      WRITE(83,REC=OREC+2) NNODG
      WRITE(83,REC=OREC+3) DTGC
      WRITE(83,REC=OREC+4) NSPOOLGC
      WRITE(83,REC=OREC+5) ITEMPC
      OREC = OREC+5
      CLOSE(83)         ! Flush the Write Buffer
      OPEN(83,FILE='fort.83',ACCESS='DIRECT',RECL=NBYTE)
C
      DO J=1,NDSETGC
C
         DO IPROC=1,NPROC
            READ(LOC2(IPROC),REC=IREC(IPROC)+1) TIMEOUTSC
            READ(LOC2(IPROC),REC=IREC(IPROC)+2) ITGC
            IREC(IPROC) = IREC(IPROC) + 2
            DO K=1, NNODP(IPROC)
               READ(LOC2(IPROC),REC=IREC(IPROC)+K) ConcBINP(K)
            ENDDO
            DO K=1,NNODP(IPROC)
               INDX = IMAP_NOD_LG(K,IPROC)
               IF (IMAP_NOD_GL(1,INDX).EQ.IPROC) THEN
                 ConcBIN(INDX) = ConcBINP(K)
               ENDIF
            ENDDO
            IREC(IPROC) = IREC(IPROC) + NNODP(IPROC)
         ENDDO
C
         WRITE(83,REC=OREC+1) TIMEOUTSC
         WRITE(83,REC=OREC+2) ITGC
         OREC = OREC + 2
         DO K=1, NNODG
            WRITE(83,REC=OREC+K) ConcBIN(K)
         ENDDO
         OREC = OREC + NNODG
C
      ENDDO
C 
 9999 CONTINUE
C
C--Close the Global and Local fort.83 Files
C
      CLOSE(83)
      DO IPROC = 1,NPROC
         CLOSE(LOC2(IPROC))
      ENDDO
C
   80 FORMAT(A80)
 2120 FORMAT(2X,E20.10,5X,I10)
 2453 FORMAT(2X,I8,2X,E15.8)
 3645 FORMAT(1X,I10,1X,I10,1X,E15.7,1X,I5,1X,I5)
C
      RETURN
      END
C------------------------------------------------------------------------
C       END POST83
C------------------------------------------------------------------------

      !------------------------------------------------------------------
      !  S U B R O U T I N E   M E R G E   S U B D O M A I N   F I L E S
      !------------------------------------------------------------------
      ! jgf52.04.03: Merge subdomain output files back into a fulldomain
      ! output file. Developed to handle nodecode and noff and to 
      ! analyze whether ghost and resident values really are the 
      ! same in all cases.
      !------------------------------------------------------------------
      subroutine mergeSubdomainFiles(filename)
      use post_global
      implicit none
      character(len=2048), intent(in) :: filename ! e.g. "fort.63" 
      integer, parameter :: NODE = 1         ! nodally based data
      integer, parameter :: ELEMENT = 2      ! elementally based data
      integer, parameter :: INTEGERTYPE = 1  ! for integer data
      integer, parameter :: REALTYPE = 2     ! for real(8) data
      integer :: dataType     ! set to INTEGERTYPE or REALTYPE
      integer :: dataLocation ! set to NODE or ELEMENT
  
      real(8) :: dataVal(3)  ! real values read from subdomain
      integer :: iRawValue  ! integer value read from subdomain
      integer :: ite        ! time step associated with dataset
      real(8) :: timeoute    ! time associated with dataset
      real(8) :: timeInc     ! time increment at which datasets were written 
      integer :: tsInc      ! time step increment at which datasets were written 
      integer :: numFullDomainValues
      integer :: numSubDomainValues
      integer :: numSets  ! number of time dependent datasets
      integer :: numVals  ! number of values in the subdomain file
      integer :: irtype   ! physical dimensionality of data (1, 2, or 3)
      real(8), allocatable :: fullDomainReal(:,:) ! dim (irtype,fdVals)
      integer, allocatable :: fullDomainInteger(:) ! dim (fdvals)
      integer :: fdLocation  ! fulldomain node or element number
      logical, allocatable :: setYet(:) ! .true. if the fulldomain value has been set
      logical, allocatable :: ambiguous(:) ! .true. if resident val != ghost val
      logical :: initialized ! .true. if data arrays have been allocated
      logical :: found       ! .true. if the subdomain file was found
      integer, allocatable  ::  loc2(:)
      character(len=2048), allocatable  ::  locname(:)
      character(len=2048) ::  outputFmt ! format string
      character(len=85) :: inline ! file metadata header
     
      integer :: i, j, k, n, iproc ! loop counters
      integer :: idum              ! toss the subdomain node number
      !
      ! allocate memory for arrays dimensioned by number of subdomains 
      allocate (loc2(mnproc))
      allocate (locname(mnproc))    
      ! set default datatype
      dataType = REALTYPE
      ! noff is the only data type that is elementally based
      if ((trim(filename).eq.'noff.100').or.
     &    (trim(filename).eq.'nodecode.63')) then
         dataType = INTEGERTYPE
      endif
      initialized = .false.
      ! 
      ! Loop over subdomains; find and open subdomain output files.
      do iproc = 1,nproc
         ! find and open subdomain output files.
         locname(iproc) = 'PE0000/'//trim(filename)
         write(locname(iproc)(3:6),'(i4.4)') (iproc-1)
         loc2(iproc) = 105 + (iproc-1)
         inquire(file=locname(iproc),exist=found)
         if (found) then
           open(loc2(iproc),file=locname(iproc))
         else
           write(*,'(a)') 'ERROR: Could not find "'//
     &        trim(locname(iproc))//
     &        '". Subdomain files of this type will not be merged.'
           return
         endif
      end do      
      ! 
      ! Loop over subdomains; read off header lines and initialize 
      ! arrays based on the characteristics of the file
      do iproc = 1,nproc
         ! read metadata header
         read(loc2(iproc),'(a85)')  inline
         ! read file characteristics header
         read(loc2(iproc),*)  numSets, numVals, timeInc, tsInc, irtype
         !
         ! Determine if this is a nodal or elemental quantity and 
         ! allocate fulldomain arrays to hold the data
         if ( (numVals.ne.nnodp(iproc)) .and. 
     &           (numVals.ne.nelp(iproc)) ) then
            write(*,'(a,i0,a,i0,a,i0,a)') 
     &      'ERROR: The number of values per data set in the file "'//
     &       trim(locname(iproc))//
     &      '" does not match the number of nodes or elements '//
     &      'in the subdomain; numVals=',numVals,
     &      ' number of nodes=',nnodp(iproc),
     &      ' number of elements=',nelp(iproc),'.'
            stop
         endif
         ! allocate data for fulldomain arrays
         if (initialized.eqv..false.) then
            if (numVals.eq.nnodp(iproc)) then
               dataLocation = NODE
               numFullDomainValues = nnodg
            endif
            if (numVals.eq.nelp(iproc)) then
               dataLocation = ELEMENT
               numFullDomainValues = nelg
            endif         
             if (dataType.eq.REALTYPE) then
               allocate(fulldomainReal(irtype,numFullDomainValues))
            endif
            if (dataType.eq.INTEGERTYPE) then
               allocate(fulldomainInteger(numFullDomainValues))
            endif
            allocate(setYet(numFullDomainValues))
            allocate(ambiguous(numFullDomainValues))
            setYet(:) = .false.
            ambiguous(:) = .false.
            initialized = .true.
         endif
      enddo
      !
      ! open fulldomain file and write header
      open(unit=100,file=trim(filename))
      write(100,'(a85)') inline
      write(100,3645) numSets,numFullDomainValues,timeInc,tsInc,irtype
3645  format(1x,i0,1x,i0,1x,e15.7,1x,i0,1x,i0)
C
C     jgf46.00 The loop over data sets is an infinite loop because the
C     number of data sets at the top of the file will be wrong if the
C     run was hot started.
C
      do ! loop infinitely over time dependent data sets
         !
         ! loop over subdomains
         do iproc=1,nproc
            ! read time and timestep
            read(loc2(iproc),fmt=2120,end=9999) timeoute,ite
            ! determine the number of lines in the dataset
            select case(dataLocation)          
            case(NODE) 
               numSubdomainValues = nnodp(iproc)
            case(ELEMENT)
               numSubdomainValues = nelp(iproc)
            case default
               write(*,'(a)') 'ERROR: Unrecognized data location.'
               stop
            end select
            !
            ! read all the datavalues in this subdomain, mapping them
            ! to their fulldomain positions 
            do k = 1,numSubdomainValues
               !
               ! Read the data according to its type
               if (dataType.eq.REALTYPE) then
                  read(loc2(iproc),fmt=*,end=9999) idum,
     &                 (dataVal(j),j=1,irtype)
               endif
               if (dataType.eq.INTEGERTYPE) then
                  read(loc2(iproc),fmt=*,end=9999) idum, iRawValue
               endif
               !
               ! Map the data according to its location.
               if (dataLocation.eq.NODE) then
                  ! look up fulldomain node number
                  fdLocation = imap_nod_lg(k,iproc)
               endif
               if (dataLocation.eq.ELEMENT) then
                  ! look up fulldomain element number, taking absolute
                  ! value so we can also process "ghost" elements
                  fdLocation = abs(imap_el_lg(k,iproc))
               endif
               !
               ! Check/set the data (resident and ghost values should match)
               if (setYet(fdLocation).eqv..true.) then
                  ! check to ensure that the value we've read from the 
                  ! subdomain for this fulldomain location matches the 
                  ! value we've already written to that fulldomain location
                  if (dataType.eq.REALTYPE) then
                     do j=1,irtype
                        ! this .ne. is subject to round off error
                        if (dataVal(j).ne.fullDomainReal(j,fdLocation)) then
                           write(*,*) 'ERROR: Differing values for '//
     &     'full domain location ',fdLocation,'. Subdomain ',iproc,
     &     ' has the value(s) ',(dataVal(n),n=1,irtype),' whereas '//
     &     'the value(s) ',(fullDomainReal(n,fdLocation),n=1,irtype),
     &     ' were already mapped.'
                           ambiguous(fdLocation) = .true.
                           exit ! break out of the irtype loop
                        endif 
                     end do
                  endif
                  if (dataType.eq.INTEGERTYPE) then
                     if ( iRawValue.ne.fullDomainInteger(fdLocation) ) then
                        outputFmt = '("ERROR: fd ",i0," Subdomain ",i0,"'//
     &                   ' sd value ",i0," vs  fd val ",i0).'
                        write(*,fmt=outputFmt) fdLocation, iproc, 
     &                   iRawValue, fullDomainInteger(fdLocation)
                        ambiguous(fdLocation) = .true.
                     endif
                  endif
               endif
               if (dataType.eq.REALTYPE) then
                  fullDomainReal(:,fdLocation) = dataVal(:)
               endif
               if (dataType.eq.INTEGERTYPE) then
                  fullDomainInteger(fdLocation) = iRawValue
               endif
               setYet(fdLocation) = .true.
            enddo ! next subdomain value

         enddo ! next subdomain
         ! all the subdomain values for this dataset have been mapped 
         ! 
         ! check to see if there were any fulldomain values that did not
         ! get set???
         do i=1, numFullDomainValues
            if (setYet(i).eqv..false.) then
               write(*,'(a,i0,a)') 'ERROR: Fulldomain location ',i,
     &            ' was never mapped.'
            endif
         end do
         !
         ! where the value at the resident location did not match the 
         ! value at the ghost location, set the ambiguous flag
         do i = 1,numFullDomainValues
            if (ambiguous(i).eqv..true.) then
               if (dataType.eq.REALTYPE) then
                  fullDomainReal(1:3,i) = -16384.d0
               endif
               if (dataType.eq.INTEGERTYPE) then
                  fullDomainInteger(i) = -1
               endif
            endif
         end do 
         !
         ! write this time dependent dataset header into the full domain file 
         write(100,2120) timeoute,ite
         ! write the mapped values for all subdomains into the fulldomain file
         ! for this dataset
         do i = 1,numFullDomainValues
            if (dataType.eq.REALTYPE) then
               write(100,'(2x,i8,3(2x,e15.8))') 
     &             i,(fullDomainReal(j,i),j=1,irtype)
            endif
            if (dataType.eq.INTEGERTYPE) then
               write(100,'(2x,i8,2x,i0)') i, fullDomainInteger(i)
            endif
         enddo 
         !              
         ! clear flags for the next dataset
         setYet(:) = .false. 
         ambiguous(:) = .false.
         fullDomainInteger(:) = -99  
      enddo ! next time dependent dataset

9999  close(100) ! jump to here when we run out of time dependent data
      do iproc = 1,nproc
         close(loc2(iproc))
      enddo
C
 2120 format(2x,e20.10,5x,i10)
 2453 format(2x,i8,3(2x,e15.8))
c
      return
      !------------------------------------------------------------------
      end subroutine mergeSubdomainFiles
      !------------------------------------------------------------------
