C----------------------------------------------------------------------------
C
C                           MODULE PREP
C
C----------------------------------------------------------------------------
C      current for ADCIRC v46.44   10/07/2006
C
C      Version 1.1 5/04/99 vjp
c      jjw fixes 053100
C      Revisions by rl 10/12/01, MEB 3/03, rl 3/03, rl 5/21/03, rl 5/18/04,
C                   vp 11/27/03 (by rl)
C      Revisions by MEB 4/04
C      jgf Created PREP11, rewrote PREP20 for ADCIRC v45.12 02/24/2006
C      jgf Created PREP13 for ADCRIC v46.00
C      vjp PREP13 updates for ADCIRC v46.44
C      vjp PREP67_68 localization updates for ADCIRC v46.44
C      vjp added Relocalize for relocalizing fort.13 and fort.15
C
C----------------------------------------------------------------------------


C---------------------------------------------------------------------------
C                S U B R O U T I N E   P R E P 1 0
C---------------------------------------------------------------------------
C
C     jgf46.28 from jgf45.16 This subroutine will break up the full
C     domain initial concentration file into subdomains. The fort.10
C     file may contain initial concentration for either 2D or 3D ADCIRC
C     runs.
C
C---------------------------------------------------------------------------
      SUBROUTINE PREP10()
C---------------------------------------------------------------------------
      USE PRE_GLOBAL
      use memory_usage
      IMPLICIT NONE
      integer :: nbytes = 0
      INTEGER I,J
      INTEGER iproc          ! subdomain index
      CHARACTER*80 header1, header2  ! header comments in unit 10 files
      CHARACTER*80 nvn_nvp           ! string representing nfen, np
      REAL(SZ) nvn           ! number of vertical nodes from unit 10 file
      REAL(SZ) nvp           ! number of horizontal nodes from unit 10 file
      INTEGER nhnn           ! horizontal nodes counter
      INTEGER nvnn           ! vertical nodes counter
      INTEGER sdu(nproc)     ! subdomain unit numbers for unit 10 files
      REAL(SZ), ALLOCATABLE :: fdData2D(:)   !(MNP)      full domain data
      REAL(SZ), ALLOCATABLE :: fdData3D(:,:) !(MNP,NFEN) full domain data
      REAL(SZ), ALLOCATABLE :: sdData2D(:)   !(MNP)      subdomain data
      REAL(SZ), ALLOCATABLE :: sdData3D(:,:) !(MNP,NFEN) subdomain data
      LOGICAL success        ! .true. if all files opened successfully

      CALL OpenPrepFiles(10, 'initial concentration         ',
     &     1, nproc, sdu, success)

      IF (.not.success) THEN
         WRITE(*,*) 'WARNING: Unit 10 files not preprocessed.'
         RETURN ! note early return
      ENDIF
C
C     Read header information from full domain unit 10 file
      READ(10,80) header1
      READ(10,80) header2
C
C     Transcribe header information.
      DO iproc = 1, nproc
         WRITE(sdu(iproc),80)  header1
         WRITE(sdu(iproc),80)  header2
      ENDDO
C
C     Check node number data for consistency (paranoia).
      READ(10,80) nvn_nvp
      READ(nvn_nvp,*) nvn, nvp
      IF ( nvn .ne. nfen .or. nvp .ne. nnodg ) then
         WRITE(*,*) 'ERROR: NVN or NVP not consistent with input data.'
         WRITE(*,*) 'NVN=',nvn,' although NFEN=',nfen
         WRITE(*,*) 'NVP=',nvp,' although NNODG=',nnodg
      ENDIF
C
C     Decompose concentration data
C
      IF (C2D_PTrans) THEN
c     read in the full domain data
         ALLOCATE ( fdData2D(MNP) )
         nbytes = 8*mnp
         call memory_alloc(nbytes)
         DO i=1, NNODG
            READ(10,*) nhnn, fdData2D(nhnn)
         ENDDO
c     write out subdomain data
         ALLOCATE ( sdData2D(MNP) )
         nbytes = 8*mnp
         call memory_alloc(nbytes)
         DO iproc = 1, nproc
            WRITE(sdu(iproc),*) nnodp(iproc)
            DO i=1, nnodp(iproc)
               sdData2D(i) = fdData2D(IMAP_NOD_LG(i,iproc))
               WRITE(sdu(iproc),*) i, sdData2D(i)
            ENDDO
         ENDDO
         DEALLOCATE ( fdData2D, sdData2D )
         nbytes = 16*mnp
         call memory_dealloc(nbytes)
      ENDIF
C
      IF (C3D_PTrans) THEN
c     read in the full domain data
         ALLOCATE ( fdData3D(MNP,NFEN) )
         nbytes = 8*mnp*nfen
         call memory_alloc(nbytes)
         DO i=1, NNODG
            DO j=1, nfen
               READ(10,*) nhnn, nvnn, fdData3D(nhnn,nvnn)
            ENDDO
         ENDDO
c     write out subdomain data
         ALLOCATE ( sdData3D(MNP,NFEN) )
         nbytes = 8*mnp*nfen
         call memory_alloc(nbytes)
         DO iproc = 1, nproc
            WRITE(sdu(iproc),*) nfen, nnodp(iproc)
            DO i=1, nnodp(iproc)
               DO j=1, nfen
                  sdData3D(i,j) = fdData3D(IMAP_NOD_LG(i,iproc),j)
                  WRITE(sdu(iproc),*) i, j, sdData3D(i,j)
               ENDDO
            ENDDO
         ENDDO
         DEALLOCATE ( fdData3D, sdData3D )
         nbytes = 16*mnp*nfen
         call memory_dealloc(nbytes)
      ENDIF
C
C     Close full domain file and all subdomain files
      CLOSE(10)
      DO iproc=1, nproc
         CLOSE(sdu(iproc))
      ENDDO
C
  80  FORMAT(A80)
C
      call memory_status()
      RETURN
C---------------------------------------------------------------------------
      END SUBROUTINE PREP10
C---------------------------------------------------------------------------


C---------------------------------------------------------------------------
C                S U B R O U T I N E   P R E P 1 1
C---------------------------------------------------------------------------
C
C     jgf45.12 This subroutine will break up the full domain initial
C     density forcing file into subdomains. The fort.11 file may
C     contain initial density, temperature, and/or salinity depending on
C     the value of IDEN in the fort.15 file.
C
C     jgf45.12 This subroutine is designed to work for baroclinic 3D
C     runs only, not 2D runs.
C
C---------------------------------------------------------------------------
      SUBROUTINE PREP11()
C---------------------------------------------------------------------------
      USE PRE_GLOBAL
      use memory_usage
      IMPLICIT NONE
      integer :: nbytes = 0
      INTEGER I,J
      INTEGER iproc          ! subdomain index
      CHARACTER*80 header1, header2  ! header comments in unit 11 files
      CHARACTER*80 nvn_nvp           ! string representing nfen, np
      REAL(SZ) nvn           ! number of vertical nodes from unit 11 file
      REAL(SZ) nvp           ! number of horizontal nodes from unit 11 file
      INTEGER nhnn           ! horizontal nodes counter
      INTEGER nvnn           ! vertical nodes counter
      INTEGER sdu(nproc)     ! subdomain unit numbers for unit 11 files
      REAL(SZ), ALLOCATABLE :: fdData1(:,:) !(MNP,NFEN) full domain data
      REAL(SZ), ALLOCATABLE :: fdData2(:,:) !(MNP,NFEN) full domain data
      REAL(SZ), ALLOCATABLE :: sdData1(:,:) !(MNP,NFEN) subdomain data
      REAL(SZ), ALLOCATABLE :: sdData2(:,:) !(MNP,NFEN) subdomain data
      LOGICAL success        ! .true. if all files opened successfully

      CALL OpenPrepFiles(11, 'initial density forcing       ',
     &     1, nproc, sdu, success)

      IF (.not.success) THEN
         WRITE(*,*) 'WARNING: Unit 11 files not preprocessed.'
         RETURN ! note early return
      ENDIF
C
C     Read header information from full domain unit 11 file
      READ(11,80) header1
      READ(11,80) header2
C
C     Transcribe header information.
      DO iproc = 1, nproc
         WRITE(sdu(iproc),80)  header1
         WRITE(sdu(iproc),80)  header2
      ENDDO
C
C     Check node number data for consistency (paranoia).
      READ(11,80) nvn_nvp
      READ(nvn_nvp,*) nvn, nvp
      IF ( nvn .ne. nfen .or. nvp .ne. nnodg ) then
         WRITE(*,*) 'ERROR: NVN or NVP not consistent with input data.'
         WRITE(*,*) 'NVN=',nvn,' although NFEN=',nfen
         WRITE(*,*) 'NVP=',nvp,' although NNODG=',nnodg
      ENDIF
C
C     Decompose density forcing data; format based on value of IDEN.
C     jgf45.12 This is designed to work for baroclinic 3D runs only, not
C     2D runs.

c     read in the full domain data
      SELECT CASE (ABS(IDEN))

      CASE(1,2,3)
         ALLOCATE ( fdData1(MNP,NFEN) )
         nbytes = 8*mnp*nfen
         call memory_alloc(nbytes)
         DO i=1, NNODG
            DO j=1, nfen
               READ(11,*) nhnn, nvnn, fdData1(nhnn,nvnn)
            ENDDO
         ENDDO

      CASE(4)
         ALLOCATE ( fdData1(MNP,NFEN) )
         ALLOCATE ( fdData2(MNP,NFEN) )
         nbytes = 16*mnp*nfen
         call memory_alloc(nbytes)
         DO i=1, NNODG
            DO j=1, nfen
               READ(11,*) nhnn, nvnn,
     &              fdData1(nhnn,nvnn),fdData2(nhnn,nvnn)
            ENDDO
         ENDDO

      END SELECT

c     write out subdomain data
      SELECT CASE (ABS(IDEN))

      CASE(1,2,3)
         ALLOCATE ( sdData1(MNP,NFEN) )
         nbytes = 8*mnp*nfen
         call memory_alloc(nbytes)
         DO iproc = 1, nproc
            WRITE(sdu(iproc),*) nfen, nnodp(iproc)
            DO i=1, nnodp(iproc)
               DO j=1, nfen
                  sdData1(i,j) = fdData1(IMAP_NOD_LG(i,iproc),j)
                  WRITE(sdu(iproc),*) i, j, sdData1(i,j)
               ENDDO
            ENDDO
         ENDDO
         DEALLOCATE ( fdData1, sdData1 )
         nbytes = 16*mnp*nfen
         call memory_dealloc(nbytes)

      CASE(4)
         ALLOCATE ( sdData1(MNP,NFEN) )
         ALLOCATE ( sdData2(MNP,NFEN) )
         nbytes = 16*mnp*nfen
         call memory_alloc(nbytes)
          DO iproc = 1, nproc
            WRITE(sdu(iproc),*) nfen, nnodp(iproc)
            DO i=1, nnodp(iproc)
               DO j=1, nfen
                  sdData1(i,j) = fdData1(IMAP_NOD_LG(i,iproc),j)
                  sdData2(i,j) = fdData2(IMAP_NOD_LG(i,iproc),j)
                  WRITE(sdu(iproc),*) i, j, sdData1(i,j), sdData2(i,j)
               ENDDO
            ENDDO
         ENDDO
         DEALLOCATE ( fdData1, fdData2 )
         DEALLOCATE ( sdData1, sdData2 )
         nbytes = 32*mnp*nfen
         call memory_dealloc(nbytes)

      END SELECT
C
C     Close full domain file and all subdomain files
      CLOSE(11)
      DO iproc=1, nproc
         CLOSE(sdu(iproc))
      ENDDO
C
  80  FORMAT(A80)
C
      call memory_status()
      RETURN
C---------------------------------------------------------------------------
      END SUBROUTINE PREP11
C---------------------------------------------------------------------------




C---------------------------------------------------------------------------
C                S U B R O U T I N E   R E L O C A L I Z E
C---------------------------------------------------------------------------
C
C     This routine allows re-localizing selected files after a prepall
C     operation.  vjp 10/2006
C     tcm v50.21 20110610 -- changed I8 to I12 format specifications
C
C---------------------------------------------------------------------------
      SUBROUTINE RELOCALIZE()
C---------------------------------------------------------------------------
      USE PRE_GLOBAL
      use memory_usage
      IMPLICIT NONE
      integer :: nbytes = 0
      INTEGER I, J, IPROC, INDX, ITEMP, idumy
      CHARACTER(14) LOCFN
      CHARACTER(80) skipped
      INTEGER,ALLOCATABLE :: LUNP(:)
C
      print *, "entering relocalize"
!      print *, "nproc = ", nproc

      allocate( lunp(nproc) )  ! logical unit number for each subdomain
      do iproc=1, nproc
         lunp(iproc) = 105 + (iproc-1)
      enddo

      if (.not.allocated(nnodp)) then
        ALLOCATE(NNODP(NPROC))
        nbytes = 4*nproc
        call memory_alloc(nbytes)
      endif

      if (.not.allocated(nelp)) then
        ALLOCATE(NELP(NPROC))
        nbytes = 4*nproc
        call memory_alloc(nbytes)
      endif

      if (.not.allocated(netap)) then
        ALLOCATE(NETAP(NPROC))
        nbytes = 4*nproc
        call memory_alloc(nbytes)
      endif

      if (.not.allocated(nfluxfp)) then
        ALLOCATE(NFLUXFP(NPROC))
        nbytes = 4*nproc
        call memory_alloc(nbytes)
      endif

      DO IPROC = 1,NPROC
         LOCFN(1:14) = 'PE0000/fort.18'
         CALL IWRITE(LOCFN,3,6,IPROC-1)
         OPEN (LUNP(IPROC),FILE=LOCFN)
      ENDDO


      print *, "from relocalize: reading local-to-global element maps"
      DO IPROC = 1,NPROC
         READ(LUNP(IPROC),'(A)') skipped    !  read past fileFmt header
Casey 100208: Changed I8 to I12.
         READ(LUNP(IPROC),'(8X,3I12)') NELG, MNEP, NELP(IPROC)
      ENDDO
!      print *, "nelg = ", nelg

      if (.not.allocated(imap_el_lg)) then
        ALLOCATE ( IMAP_EL_LG(MNEP, NPROC) )
        nbytes = 4*nproc*mnep
        call memory_alloc(nbytes)
      endif

      DO IPROC = 1,NPROC
      DO I=1, NELP(IPROC)
Casey 100208: Changed I8 to I12.
         READ(LUNP(IPROC),'(I12)') idumy
         IMAP_EL_LG(I,IPROC) = abs(idumy)
      ENDDO
      ENDDO

      print *, "from relocalize: reading local-to-global node maps"
      DO IPROC = 1,NPROC
         READ(LUNP(IPROC),'(8X,3I12)') NNODG, MNPP, NNODP(IPROC)   !tcm v50.21
      ENDDO
!      print *, "nnodg = ", nnodg

      if (.not.allocated(imap_nod_lg)) then
         ALLOCATE ( IMAP_NOD_LG(MNPP, NPROC) )
         nbytes = 4*nproc*mnpp
         call memory_alloc(nbytes)
      endif

      DO IPROC = 1,NPROC
      DO I=1, NNODP(IPROC)
         READ(LUNP(IPROC),'(I12)') idumy        !tcm v50.21
         IMAP_NOD_LG(I,IPROC) = abs(idumy)
      ENDDO
      ENDDO

C  This section for prep15
      IF ((PREP_15.eqv..true.).or.(PREP_20.eqv..true.)) THEN
        print *, "from relocalize: reading nfluxf for each subdomain"
        DO IPROC = 1,NPROC
           READ(LUNP(IPROC),'(8X,I12)') NFLUXFP(IPROC)        !tcm v50.21
        ENDDO

        print *, "from relocalize: reading neta for each subdomain"
        DO IPROC = 1,NPROC
           READ(LUNP(IPROC),'(8X,3I12)') idumy, NETA_MAX, NETAP(IPROC)   !tcm v50.21
        ENDDO

        if (.not.allocated(obnode_lg)) then
           ALLOCATE ( OBNODE_LG(NETA_MAX, NPROC) )
           nbytes = 4*nproc*neta_max
           call memory_alloc(nbytes)
        endif

        print *, "from relocalize: reading open boundary table"
        DO IPROC = 1,NPROC
           DO I=1, NETAP(IPROC)
              READ(LUNP(IPROC),'(I12)') OBNODE_LG(I,IPROC)    !tcm v50.21
           ENDDO
        ENDDO
      ENDIF

! Build Global-to-Local Node Map
      if (.not.allocated(itotproc)) then
        ALLOCATE ( ITOTPROC(NNODG) )
        nbytes = 4*nnodg
        call memory_alloc(nbytes)
      endif
      DO I = 1,NNODG
         ITOTPROC(I) = 0
      ENDDO
      DO IPROC = 1,NPROC
         DO I = 1,NNODP(IPROC)
            INDX = IMAP_NOD_LG(I,IPROC)
            ITOTPROC(INDX) = ITOTPROC(INDX) + 1
         ENDDO
      ENDDO
      MNEI = 0
      DO I = 1,NNODG
         IF (ITOTPROC(I) .gt. MNEI) MNEI = ITOTPROC(I)
      ENDDO
      print *, "MNEI = ", MNEI

      if (.not.allocated(imap_nod_gl2)) then
        ALLOCATE( IMAP_NOD_GL2( 2*MNEI, NNODG) )
        nbytes = 8*mnei*nnodg
        call memory_alloc(nbytes)
      endif
      print *, "allocated imap_nod_GL2"

      DO I = 1,NNODG
         ITOTPROC(I) = 0
      ENDDO
      DO IPROC = 1,NPROC
         DO J = 1,NNODP(IPROC)
            INDX = IMAP_NOD_LG(J,IPROC)
            ITOTPROC(INDX) = ITOTPROC(INDX) + 1
            ITEMP = (ITOTPROC(INDX)-1)*2 + 1
            IMAP_NOD_GL2(ITEMP,INDX) = IPROC
            IMAP_NOD_GL2(ITEMP+1,INDX) = J
         ENDDO
      ENDDO

      !jgf50.35: Need this for --prep13 option.
      if (.not.allocated(imap_nod_gl)) then
         ALLOCATE (IMAP_NOD_GL(2,NNODG))
         nbytes = nbytes + 8*mnp
         call memory_alloc(nbytes)
      endif
      print *, "allocated imap_nod_GL"
      ! jgf50.35: Formulate the global-to-local
      ! mapping for resident nodes
      DO IPROC=1,NPROC
         DO J=1,NNODP(IPROC)
            INDX = IMAP_NOD_LG(J,IPROC)
            IF (ITOTPROC(INDX).eq.1) THEN
               IMAP_NOD_GL(1,INDX) = IPROC
               IMAP_NOD_GL(2,INDX) = J
            ENDIF
         ENDDO
      ENDDO
      call memory_status()
      print *, "leaving relocalize"
      RETURN
C---------------------------------------------------------------------------
      END SUBROUTINE RELOCALIZE
C---------------------------------------------------------------------------


C---------------------------------------------------------------------------
C                S U B R O U T I N E   P R E P 1 3
C---------------------------------------------------------------------------
C
C     jgf46.00 This subroutine will break up the full domain nodal
C     attributes file into subdomains.
C
C     jgf48.47 Rewritten to only open a smaller number of subdomains
C     at a time (256 by default)
C
C     jgf48.50 Fixed bug (remove ALLOCATABLE from vars that don't need
C     it). Also adding following documentation:
C
C     BY DEFAULT, ONLY 256 SUBDOMAINS WILL BE PREPPED AT A TIME TO AVOID
C     OPENING TOO MANY FILES ON CERTAIN PLATFORMS. THIS NUMBER CAN BE
C     CONTROLLED BY THE PARAMETER MAXOPENFILES.
C
C---------------------------------------------------------------------------
      SUBROUTINE PREP13()
C---------------------------------------------------------------------------
      USE PRE_GLOBAL
      use memory_usage
      IMPLICIT NONE
C
      integer :: nbytes = 0
      INTEGER i              ! node loop counter
      INTEGER k              ! attribute loop counter
      INTEGER ll             ! line loop counter
      INTEGER m              ! attribute default value counter
      INTEGER iproc          ! subdomain loop counter
      INTEGER sdu(nproc)     ! subdomain unit number for unit 13 files
      INTEGER NumOfNodes     ! must match the no. of nodes in grid (unit 14)
      INTEGER NAttr          ! number of nodal attributes in the file
      INTEGER NumNotDefault  ! number of nodes specified in the file
      CHARACTER(len=80) header   ! header comments in unit 13 files
      CHARACTER(len=80) AttrName ! label for attribute
      CHARACTER(len=80) Units    ! label for physical units
      CHARACTER(len=80) Skipped  ! data we want to skip over
      REAL(SZ) DefaultVal(12) ! default value of attribute
      INTEGER NoOfVals       ! at each node for an attribute
      INTEGER Mode           !=0 to count, =1 to write
      LOGICAL success        ! .true. if all files opened successfully
      INTEGER, ALLOCATABLE :: SDNumND(:,:)  ! subdomain # of nodes not default
C     jgf48.47 Do the decomposition for a max of 256 subdomains at a
C     time ... some platforms/compilers limit the number of files that
C     can be open at any one time.
      INTEGER, PARAMETER :: maxOpenFiles = 256
      INTEGER startProc
      INTEGER endProc
      INTEGER deltaProc
C
C     Perform decomposition over range of subdomains.
      startProc = 1
      DO WHILE ( startProc .lt. nproc )
         deltaProc = nproc - startProc
         IF ( deltaProc .gt. maxOpenFiles ) deltaProc = maxOpenFiles
         endProc = startProc + deltaProc
C
C        Open full domain and subdomain fort.13 files.
         CALL OpenPrepFiles(13, 'nodal attributes              ',
     &       startProc, endProc, sdu, success)
C
         IF (.not.success) THEN
            WRITE(*,*) 'WARNING: Unit 13 files not preprocessed.'
            RETURN ! note early return
         ENDIF
C
C        Read header information from full domain unit 13 file
         READ(13,'(A80)') header
         READ(13,*) NumOfNodes     ! number of nodes according to unit 13
         READ(13,*) NAttr          ! number of attributes in the unit 13 file
C
C        Check to make sure that the number of nodes in the nodal
C        attributes file is the same as in the grid file (unit 14).
         IF (NumOfNodes.NE.NNODG) THEN
            WRITE(6,9900)
 9900       FORMAT(////,1X,'!!!!!!!!!!  FATAL ERROR  !!!!!!!!!',
     &        //,1X,'The number of nodes in the grid file (unit 14) and'
     &        /,1X,'the nodal attributes file (unit 13) must match.',
     &        //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
            STOP                   ! We're toast.
         ENDIF
C
C        Transcribe header information into subdomain unit 13 files
         DO iproc = startProc, endProc
            WRITE(sdu(iproc),'(A80)') header
            WRITE(sdu(iproc),*) NNODP(iproc)
            WRITE(sdu(iproc),*) NAttr
         ENDDO
C
C        Transcribe attribute names from full domain file to subdomains.
         DO k=1, NAttr
            READ(13,'(A80)') AttrName
            READ(13,'(A80)') Units
            READ(13,*) NoOfVals
            READ(13,*) (DefaultVal(m),m=1,NoOfVals)
            DO iproc=startProc, endProc
               WRITE(sdu(IPROC),'(A80)') AttrName
               WRITE(sdu(IPROC),'(A80)') Units
               WRITE(sdu(IPROC),*) NoOfVals
               WRITE(sdu(IPROC),'(12(1x,e16.9))')
     &          (DefaultVal(m),m=1,NoOfVals)
            END DO
         END DO
C
C        Allocate and initialize the matrix for the number of Non Default
C        nodes in each SubDomain for each nodal attribute
         ALLOCATE(SDNumND(nproc,NAttr))
         nbytes = 8*nproc*nattr
         call memory_alloc(nbytes)
         DO iproc=startProc,endProc
            DO k=1, NAttr
               SDNumND(iproc,k)=0
            END DO
         END DO
C        We need to figure out how many nodes go into each subdomain
C        for each attribute.
         CALL processNodalAttr(NAttr, 0, sdu, SDNumND,
     &      startProc, endProc)
C
C     Now rewind and advance to the beginning of the data again
         REWIND(13)
         DO ll=1, 3
            READ(13,*) skipped    ! skip header, NumOfNodes, NAttr
         END DO
         DO k=1, NAttr
            DO ll=1, 4
               READ(13,*) skipped ! skip AttrName,Units,NoOfVals,Default
            END DO
         END DO
C
C        Now read each of the nodal attributes and transcribe them to the
C        appropriate subdomain.
         CALL processNodalAttr(NAttr, 1, sdu, SDNumND,
     &      startProc, endProc)
         DEALLOCATE(SDNumND)
         nbytes = 8*nproc*nattr
         call memory_dealloc(nbytes)
C
C        Close full domain and subdomain files
         CLOSE (13)
         DO iproc=startProc, endProc
            CLOSE(sdu(iproc))
         ENDDO
         startProc = endProc + 1
      END DO
C
      if (allocated(NodalAttributes)) then
         DEALLOCATE(NodalAttributes)
         nbytes = 4*nwp
         call memory_dealloc(nbytes)
      endif
C
      call memory_status()
      RETURN
C---------------------------------------------------------------------------
      END SUBROUTINE PREP13
C---------------------------------------------------------------------------

C     ----------------------------------------------------------------
C        S U B R O U T I N E   P R O C E S S   N O D A L   A T T R
C     ----------------------------------------------------------------
C
C     jgf46.00 Subroutine to support PREP13. This subroutine is called
C     twice; once to determine the number of nodes with non-default
C     values going into each subdomain, then a second time to actually
C     place the data in the subdomain files.
C
C     This is necessary because the attributes in the subdomain files
C     must each have the number of non-default values at the top, and
C     this information cannot be known until we have processed the
C     entire fulldomain file.
C
C     ----------------------------------------------------------------
      SUBROUTINE ProcessNodalAttr(NAttr, mode, sdu, SDNumND,
     &   startProc, endProc)
      USE PRE_GLOBAL
      IMPLICIT NONE
C
      INTEGER,intent(in) :: NAttr  ! number of attributes in the file
      INTEGER,intent(in) :: Mode   !=0 to count and return, =1 to write
      INTEGER,intent(in),dimension(nproc) :: sdu !i/o unit number array
      INTEGER,intent(inout),dimension(nproc,NAttr) :: SDNumND
      INTEGER,intent(in) :: startProc
      INTEGER,intent(in) :: endProc
      INTEGER NumNotDefault      ! number of nodes specified in the file
      INTEGER NumCol             ! number of values per node for an attr
      INTEGER NodeNum            ! full domain node number
      INTEGER SDNode             ! subdomain node number
      INTEGER i                  ! node loop counter
      INTEGER j                  ! column loop counter
      INTEGER k                  ! attribute loop counter
      INTEGER m                  ! mapping loop counter
      INTEGER iproc              ! subdomain loop counter
      INTEGER iproc2             ! mapped subdomain
      CHARACTER(len=80) AttrName ! label for attribute
      REAL(SZ), ALLOCATABLE :: AttrData(:)  ! attribute data
      CHARACTER(len=80) Skipped  ! data we want to skip over
C
      DO k=1, NAttr
         READ(13,*) AttrName
         IF (Mode.eq.1) THEN
            DO iproc=startProc,endProc
               WRITE(sdu(iproc),'(A80)') AttrName
            END DO
         ENDIF
         SELECT CASE (AttrName)
         CASE("primitive_weighting_in_continuity_equation")
            NumCol=1
         CASE("surface_submergence_state")
            NumCol=1
         CASE("quadratic_friction_coefficient_at_sea_floor")
            NumCol=1
         CASE("mannings_n_at_sea_floor")
            NumCol=1
         CASE("bottom_roughness_length")
            NumCol=1
         CASE("chezy_friction_coefficient_at_sea_floor")
            NumCol=1
         CASE("sea_surface_height_above_geoid")
            NumCol=1
         CASE("surface_directional_effective_roughness_length")
            NumCol=12
         CASE("surface_canopy_coefficient")
            NumCol=1
         CASE("bridge_pilings_friction_parameters")
            NumCol=4
         CASE
     &   ("average_horizontal_eddy_viscosity_in_sea_water_wrt_depth")
            NumCol=1
         CASE
     &   ("average_horizontal_eddy_diffusivity_in_sea_water_wrt_depth")
            NumCol=1
         CASE
     &   ("min_and_max_primitive_weighting_in_continuity_equation")
            NumCol=2
         CASE
     &   ("subtract_from_depth_to_equal_mean_sea_level")
            NumCol=1
Casey 100210: Allow SWAN to handle wave refraction as a nodal attribute.
         CASE("wave_refraction_in_swan")
            NumCol=1
         CASE("elemental_slope_limiter")
            NumCol=1
Corbitt 120328: Allow local advection as a nodal attribute.
         CASE("advection_state")
            NumCol=1
         CASE DEFAULT
            NumCol=1
            IF (Mode.eq.0) THEN
               WRITE(6,1001)    ! Nodal Attributes file
               WRITE(6,1021) AttrName ! contains invalid name
            ELSE
               WRITE(6,1031) AttrName ! Process 1st column only
            ENDIF
         END SELECT
         ALLOCATE(AttrData(NumCol))
         READ(13,*) NumNotDefault
         IF (Mode.eq.1) THEN
            DO iproc=startProc,endProc
               WRITE(sdu(iproc),*) SDNumND(iproc,k)
            END DO
         ENDIF
         DO i=1, NumNotDefault
            READ(13,*) NodeNum, (AttrData(j),j=1,NumCol)
            IF (ITOTPROC(NodeNum).eq.1) THEN
               iproc = IMAP_NOD_GL(1,NodeNum)
               IF ( (iproc.lt.startProc) .or. (iproc.gt.endProc) ) THEN
                  CYCLE ! skip it if it does not map to our range of procs
               ENDIF
               IF (Mode.eq.0) SDNumND(iproc,k) = SDNumND(iproc,k)+1
               IF (Mode.eq.1) THEN
                  SDNode = IMAP_NOD_GL(2,NodeNum)
                  WRITE(sdu(iproc),1100) SDNode,(AttrData(j),j=1,NumCol)
               ENDIF
            ELSE
               DO m=1, ITOTPROC(NodeNum)
                  iproc2 = IMAP_NOD_GL2(2*(m-1)+1,NodeNum)
                  DO iproc=startProc, endProc
                     IF (iproc.EQ.iproc2) THEN !f.d. node maps to this s.d.
                        IF (Mode.eq.0) THEN
                           SDNumND(iproc,k)=SDNumND(iproc,k)+1
                        ENDIF
                        IF (Mode.eq.1) THEN
                           SDNode = IMAP_NOD_GL2(2*(m-1)+2,NodeNum)
                           WRITE(sdu(iproc),1100) SDNode,
     &                          (AttrData(j),j=1,NumCol)
                        ENDIF
                     ENDIF
                  END DO
               END DO
            END IF
         END DO
         DEALLOCATE(AttrData)
         IF (Mode.eq.1) THEN
            WRITE(6,'(A25,A80)') '     Finished processing ', AttrName
            WRITE(6,*) 'for processor range ',startProc,' to ',endProc
         ENDIF
      END DO
C
C
 1001 FORMAT('ERROR: The Nodal Attributes File (unit 13)')
 1021 FORMAT('contains invalid name: ',A80)
 1031 FORMAT('WARNING: Processed only one column of unrecognized ',A80)
 1100 FORMAT(I10,32000(2X,E16.8))
C
      RETURN
C     ----------------------------------------------------------------
      END SUBROUTINE ProcessNodalAttr
C     ----------------------------------------------------------------


      SUBROUTINE PREP14()
      USE PRE_GLOBAL
C
C---------------------------------------------------------------------------C
C                     (  Serial Version  2/28/98  )                         C
C  This routine writes a Local Grid file "fort.14" file for each subdomain  C
C  using the domain decomposition of the ADCIRC grid created by the routine C
C  DECOMP.                                                                  C
C                                                                           C
C  The Decomposition Variables are defined in the include file adcprep.inc  C
C  This version is compatible with ADCIRC version 34.03                     C
C                                                                           C
C---------------------------------------------------------------------------C
C
      IMPLICIT NONE
      INTEGER I,I1,J,K,M,ETYPE,ITEMP,ITEMP2,ILNODE,ILNODE2,ILNODE3
      INTEGER JD,JG,JP,IPROC,IPROC2,IPROC3,DISC,BBN,IBP
      INTEGER INDX,INDEX2,ITOT,ITYPE,NUMS(10)
      CHARACTER LOCFN*14,PE*6
      CHARACTER*80 OUTMSG
C
      ETYPE = 3   ! The only Element-Type supported by ADCIRC is 3.
C
C--------------------------------------------------------------------------
C--MAIN LOOP:   Write a Local Grid File ( fort.14 ) for each PE
C--------------------------------------------------------------------------
C
      NETA_MAX = 0   ! max number of open boundary nodes on any subdomain

      DO 1000 IPROC = 1,NPROC
C
         LOCFN(1:14) = 'PE0000/fort.14'
         CALL IWRITE(LOCFN,3,6,IPROC-1)
         OPEN (14,FILE=LOCFN)
C
C--------------------------------------------------------------------------
C--OPEN BOUNDARY NODES PROCESSING BEGINS HERE
C--------------------------------------------------------------------------
C
C--Partition the open boundary nodes between various processors
C
         NETAP(IPROC) = 0
         DO K=1, NOPE
            NVDLLP(K) = 0
            DO J=1, NETA
               OBNODE_LG(J,IPROC) = 0
               NBDVP(K,J) = 0
            ENDDO
         ENDDO
C
         ITOT = 0
         DO K = 1,NOPE
            DO I = 1,NVDLL(K)
               ITOT = ITOT + 1
               INDX = NBDV(K,I)
               DO J = 1,ITOTPROC(INDX)
                  ITEMP = (J-1)*2+1
                  IPROC2 = IMAP_NOD_GL2(ITEMP,INDX)
                  ILNODE = IMAP_NOD_GL2(ITEMP+1,INDX)
                  IF (IPROC.EQ.IPROC2) THEN
                    NETAP(IPROC) = NETAP(IPROC)+1
                    NVDLLP(K) = NVDLLP(K) + 1
                    NBDVP(K,NVDLLP(K)) = ILNODE
                    OBNODE_LG(NETAP(IPROC),IPROC)=ITOT
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
         IF (NETAP(IPROC) > NETA_MAX) NETA_MAX = NETAP(IPROC)
C
         NOPEP(IPROC) = 0
         DO K = 1,NOPE
            IF (NVDLLP(K).NE.0) THEN
              NOPEP(IPROC) = NOPEP(IPROC) + 1
            ENDIF
         ENDDO
C
C
C--------------------------------------------------------------------------
C--LAND BOUNDARY NODES PROCESSING BEGINS HERE
C--------------------------------------------------------------------------
C
C--Partition Land Boundary Segments between PEs
C
         NVELP(IPROC) = 0
         DO K = 1,NBOU
            NVELLP(K) = 0
            IBTYPEP(K,IPROC) = IBTYPE(K)
         ENDDO
C
         DO K = 1,NBOU
C
C--Weir Land Boundary Node-Pair Case
Cmod vjp 3/8/99
C  mod to allow that each of Weir-node pair might be ghosts nodes
C
            IF ((IBTYPE(K).EQ.4).OR.(IBTYPE(K).EQ.24)) THEN
               DO I = 1,NVELL(K)
                  INDX = NBVV(K,I)
                  INDEX2 = IBCONNR(K,I)
                  DO J = 1,ITOTPROC(INDX)
                     ITEMP = (J-1)*2 + 1
                     IPROC2  =  IMAP_NOD_GL2(ITEMP,INDX)
                     ILNODE2 =  IMAP_NOD_GL2(ITEMP+1,INDX)
                     IF (IPROC.EQ.IPROC2) THEN
                        DO JD = 1, ITOTPROC(INDEX2)
                           ITEMP2 = (JD-1)*2 + 1
                           IPROC3  = IMAP_NOD_GL2(ITEMP2,INDEX2)
                           ILNODE3 = IMAP_NOD_GL2(ITEMP2+1,INDEX2)
                           IF (IPROC.EQ.IPROC3) THEN
                              NVELP(IPROC) = NVELP(IPROC) + 1
                              NVELLP(K) = NVELLP(K) + 1
                              LBINDEX_LG(K,NVELLP(K)) = I
                              NBVVP(K,NVELLP(K))   = ILNODE2
                              IBCONNRP(K,NVELLP(K)) = ILNODE3
                           ENDIF
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
C
C--All Other Land Boundary Node types
C
            ELSE
C
               DO I = 1,NVELL(K)
                  INDX = NBVV(K,I)
                  DO J = 1,ITOTPROC(INDX)
                     ITEMP = (J-1)*2 + 1
                     IPROC2 =  IMAP_NOD_GL2(ITEMP,INDX)
                     ILNODE =  IMAP_NOD_GL2(ITEMP+1,INDX)
                     IF (IPROC.EQ.IPROC2) THEN
                        NVELP(IPROC) = NVELP(IPROC) + 1
                        NVELLP(K) = NVELLP(K) + 1
                        LBINDEX_LG(K,NVELLP(K)) = I
                        NBVVP(K,NVELLP(K)) = ILNODE
                     ENDIF
                  ENDDO
               ENDDO
C
            ENDIF
C
         ENDDO
C
Cmod 05/18/2004 rl -- I don't think this next part is the correct
c  way to handle islands.  Rather, if an island is split by a domain, it
c  should remain an island.  This will ensure that the boundary is
c  closed.  The only error would occur in ghost node space, which is
c  not a problem since the answers are not used there anyway.

Cmod 12/18/98 vjp --this section re-written
C--If a PE has only part of a closed internal land boundary
C  modify its local IBTYPE to be an external land boundary segment
C  of the same type by decrementing its IBTYPE.
C  and remove a closing loop node if present

C
c        DO K=1, NBOU
c          IF (NVELLP(K).LT.NVELL(K)) THEN
c            IF (  (IBTYPEP(K,IPROC).EQ.1)
c    &         .OR.(IBTYPEP(K,IPROC).EQ.11)
c    &         .OR.(IBTYPEP(K,IPROC).EQ.21)) THEN
c decrement ibtype
c              IBTYPEP(K,IPROC) = IBTYPEP(K,IPROC)-1
C remove loop closing node
c              IF (NVELLP(K).GT.1.AND.
c    &           NBVVP(K,NVELLP(K)).EQ.NBVVP(K,1)) THEN
c                NVELLP(K) = NVELLP(K)-1
c              ENDIF
c            ENDIF
c          ENDIF
c        ENDDO

C If a segment contains only one node, remove the segment from the list
C (NOTE: rl 5/18/04 I don't see how this could possibly happen, including
C  ghost nodes)

         DO K=1, NBOU
           IF (NVELLP(K).EQ.1) NVELLP(K) = 0
         ENDDO

C
C--Count the number of land boundary segments on PE IPROC.
C
         NBOUP(IPROC) = 0
         DO K = 1,NBOU
            IF (NVELLP(K).NE.0) THEN
              NBOUP(IPROC) = NBOUP(IPROC) + 1
            ENDIF
         ENDDO
C
C--Count to check correctness of NVELP
C
         DISC=0  ! LB Nodes with non-zero normal discharge
         BBN=0   ! Mainland Barrier Boundary Nodes
         IBP=0   ! Internal Barrier Boundary Pairs
         ITEMP = 0
C
C     jgf46.21 Added support for IBTYPE=52.
         DO 400 K=1,NBOU
            IF (NVELLP(K).EQ.0) GOTO 400
            ITYPE = IBTYPEP(K,IPROC)
! kmd - added for rivers in baroclinic simulation
            IF (ABS(ITYPE/100).EQ.1) THEN
               ITYPE = (ABS(ITYPE)-100)*(ITYPE/ABS(ITYPE))
            END IF
            ! jgf50.21: Added support for IBTYPE=32 and replaced
            ! if/then statements with a select statement.
            select case(ITYPE)
            case(2,12,22,32,52)
              DISC = DISC + NVELLP(K)
            case(3,13,23)
              BBN = BBN + NVELLP(K)
            case(4,24)
              IBP = IBP + NVELLP(K)
            case default
               ITEMP = ITEMP + NVELLP(K)
            end select
            I1 = 0
            DO I=1,NVELLP(K)
               IF ((ITYPE.EQ.1).OR.(ITYPE.EQ.11).OR.
     &             (ITYPE.EQ.21)) THEN
                 IF ((I.EQ.NVELLP(K)).AND.(NBVVP(K,I).NE.I1)) THEN
                   ITEMP = ITEMP + 1
                 ENDIF
               ENDIF
               IF (I.EQ.1) I1 = NBVVP(K,I)
            ENDDO
 400     CONTINUE
C
c        print *, IPROC-1,ITEMP,DISC,BBN,2*IBP
         ITEMP  = ITEMP + DISC + BBN + 2*IBP
         IF (ITEMP.NE.NVELP(IPROC)) THEN
c          print *, "changed value from ",NVELP(IPROC)," to ",ITEMP
           NVELP(IPROC) = ITEMP
         ENDIF
         IF (NVELP(IPROC)+1.GT.MNVEL) THEN
           print *, "NVEL exceeds parameter value MNVEL on PE",IPROC
           print *, "local NVEL value = ",ITEMP
           stop
         ENDIF
C
C--Construct a LBCODE for each Land Boundary Node of this PE
C
         JP=0
         DO K = 1,NBOU
            DO I=1, NVELLP(K)
               JP = JP+1
               LBCODEP(JP,IPROC) = IBTYPEP(K,IPROC)
            ENDDO
         ENDDO
C
C--Determine whether there are any normal flow boundaries local to PE
C
! kmd - changed for rivers in baroclinic simulations
         NFLUXFP(IPROC) = 0
         DO K=1, NBOU
            IF (NVELLP(K).GT.0) THEN
              ITYPE=IBTYPE(K)
              IF (ABS(ITYPE/100).EQ.1) THEN
                 ITYPE = (ABS(ITYPE)-100)*(ITYPE/ABS(ITYPE))
              END IF
              IF ((ITYPE.EQ.2).OR.(ITYPE.EQ.12)
     &            .OR.(ITYPE.EQ.32)
     &            .OR.(ITYPE.EQ.22).OR.(ITYPE.EQ.52)) THEN
                NFLUXFP(IPROC) = 1
              ENDIF
            ENDIF
         ENDDO
C
C--------------------------------------------------------------------------
C--BEGIN WRITING LOCAL GRID ( fort.14 ) FILE HERE
C--------------------------------------------------------------------------
C
C--Write Mesh Data
C
          WRITE(14,80) AGRID
C
          NUMS(1) = NELP(IPROC)
          NUMS(2) = NNODP(IPROC)
c
cjgf45.06    CALL INSERT(SIZEMSG,OUTMSG,NUMS,2)
cjgf45.06    WRITE(14,80) OUTMSG
          WRITE(14,43) NELP(IPROC),NNODP(IPROC) !jgf45.06
C
          DO J = 1,NNODP(IPROC)
             INDX = IMAP_NOD_LG(J,IPROC)
             WRITE(14,44) J,X(INDX),Y(INDX),DP(INDX)
          ENDDO
c
          DO J = 1,NELP(IPROC)
             WRITE(14,45) J,ETYPE,NNEP(1,J,IPROC),NNEP(2,J,IPROC),
     &                           NNEP(3,J,IPROC)
          ENDDO
 43       FORMAT(2I8)
 44       FORMAT(I8,3(E24.12))
 45       FORMAT(5I8)
C
C--Write Open Boundary Data
C
          CALL NEWINDEX(NOPEMSG,OUTMSG,NOPEP(IPROC))
          WRITE(14,80) OUTMSG
C
          CALL NEWINDEX(NETAMSG,OUTMSG,NETAP(IPROC))
          WRITE(14,80) OUTMSG
C
          ITOT = 0
          DO K = 1,NOPE
             IF (NVDLLP(K).GT.0)THEN
                ITOT = ITOT + 1
Casey 090304: Added the following section.  If we are coupling to SWAN,
C             then we also want to give the global number of each
C             boundary segment.
#ifndef ADCSWAN
                CALL NEWINDEX(NVDLLMSG(K),OUTMSG,NVDLLP(K))
#else
                NUMS(1) = NVDLLP(K)
                NUMS(2) = K
                CALL INSERT(NVDLLMSG(K),OUTMSG,NUMS,2)
#endif
                WRITE(14,80) OUTMSG
                DO I = 1,NVDLLP(K)
                   WRITE(14,*) NBDVP(K,I)
                ENDDO
             ENDIF
          ENDDO
C
C--Write Land Boundary Data
C
          CALL NEWINDEX(NBOUMSG,OUTMSG,NBOUP(IPROC))
          WRITE(14,80) OUTMSG
C
          CALL NEWINDEX(NVELMSG,OUTMSG,NVELP(IPROC))
          WRITE(14,80) OUTMSG
C
          DO K = 1,NBOU
             IF(NVELLP(K).GT.0)THEN
                ITYPE = IBTYPEP(K,IPROC)
                NUMS(1) = NVELLP(K)
                NUMS(2) = ITYPE
Casey 090304: Added the following section.  If we are coupling to SWAN,
C             then we also want to give the global number of each
C             boundary segment.
#ifndef ADCSWAN
                CALL INSERT(NVELLMSG(K),OUTMSG,NUMS,2)
#else
                NUMS(3) = NOPE + K
                CALL INSERT(NVELLMSG(K),OUTMSG,NUMS,3)
#endif
                WRITE(14,80) OUTMSG
C
                IF ((ITYPE.NE.3).AND.(ITYPE.NE.13).AND.
     &             (ITYPE.NE.23).AND.(ITYPE.NE.4).AND.
     &             (ITYPE.NE.24)) THEN
                  DO I = 1,NVELLP(K)
                     WRITE(14,'(I8)') NBVVP(K,I)
                  ENDDO
                ELSEIF ((ITYPE.EQ.3).OR.(ITYPE.EQ.13).OR.
     &                  (ITYPE.EQ.23)) THEN
                  DO I = 1,NVELLP(K)
                     INDX = LBINDEX_LG(K,I)
                     WRITE(14,81) NBVVP(K,I),BAR1(K,INDX),BAR2(K,INDX)
                  ENDDO
C
                ELSEIF ((ITYPE.EQ.4).OR.(ITYPE.EQ.24)) THEN
                  DO I = 1,NVELLP(K)
                     INDX = LBINDEX_LG(K,I)
                     WRITE(14,82) NBVVP(K,I),IBCONNRP(K,I),
     &                        BAR1(K,INDX),BAR2(K,INDX),BAR3(K,INDX)
                  ENDDO
                ENDIF
             ENDIF
          ENDDO
C
          CLOSE(14)
C
1000  CONTINUE
C
C--Print Summary of Boundary Node Decomposition
C
      print *, " "
      print *, "Boundary Node Decomposition Data"
      print *, "DOMAIN      NOPE    NETA    NBOU  NVEL    NWEIR"
      WRITE(*,90)  "GLOBAL",NOPE, NETA, NBOU, NVEL, NWEIR
      DO IPROC=1, NPROC
         PE(1:6) = 'PE0000'
         CALL IWRITE(PE,3,6,IPROC-1)
         WRITE(*,90)  PE,NOPEP(IPROC),NETAP(IPROC),
     &                NBOUP(IPROC),NVELP(IPROC),NWEIRP(IPROC)
      ENDDO
C
  80  FORMAT(A80)
  81  FORMAT(I8,2X,E13.6,2X,E13.6)
  82  FORMAT(I8,2X,I8,2X,E13.6,2X,E13.6,2X,E13.6)
  90  FORMAT(1X,A6,5I8)
C
      RETURN
      END

C---------------------------------------------------------------------------C
C                     (  Serial Version  2/28/98  )                         C
C  This routine writes a Local Input file "fort.15" file for each subdomain C
C  using the domain decomposition of the ADCIRC grid created by the routine C
C  DECOMP.                                                                  C
C                                                                           C
C  The Decomposition Variables are defined in the include file adcprep.inc  C
C  This version is compatible with ADCIRC version 34.03                     C
C                                                                           C
C           Modifications by RL on 10/9/01 to accomodate NWS = -2           C
C---------------------------------------------------------------------------C
      SUBROUTINE PREP15()
      USE PRE_GLOBAL
      use memory_usage
      USE HARM, ONLY : NAMEFR
      IMPLICIT NONE
      integer :: nbytes = 0
      INTEGER I,J,K,M,JG,JP,KK
      INTEGER INDX,ITOT,ILNODE,IPROC,IPROC2,ITYPE,NUMS(10)
      CHARACTER LOCFN*14,PE*6
      CHARACTER*80 OUTMSG
C
C--Write a Local Input file ( fort.15 ) for each PE
C
      ! max number of stations in any subdomain
      NSTAE_MAX = 0; NSTAV_MAX = 0; NSTAM_MAX = 0; NSTAC_MAX = 0
C
      DO 1000 IPROC = 1,NPROC
C
         LOCFN(1:14) = 'PE0000/fort.15'
         CALL IWRITE(LOCFN,3,6,IPROC-1)
         OPEN (15,FILE=LOCFN)
C
         WRITE(15,80) RUNDES
         WRITE(15,80) RUNID
         WRITE(15,80) OVERMSG
         WRITE(15,80) ABOUTMSG
         WRITE(15,80) SCREENMSG
         WRITE(15,80) HOTMSG
         WRITE(15,80) ICSMSG
         WRITE(15,80) IMMSG
         IF (CBaroclinic) THEN  !jgf46.28
            WRITE(15,80) IDENMSG
         ENDIF
         WRITE(15,80) IBFMSG
         WRITE(15,80) IFAMSG
         WRITE(15,80) ICAMSG
         WRITE(15,80) ICATMSG
         WRITE(15,80) NWPMSG
         IF (NWP.gt.0) THEN !jgf46.00 write nodal attributes
            DO I=1, NWP
               WRITE(15,80) NodalAttributes(I)
            ENDDO
         ENDIF
         WRITE(15,80) NCORMSG
         WRITE(15,80) NTIPMSG
C     jgfdebug46.02 Added check for NWS=45 to write NWS=5
         IF (NWS.EQ.45) THEN
            WRITE(15,'(A1)') "5"
         ELSE
            WRITE(15,80) NWSMSG
         ENDIF
         WRITE(15,80) RAMPMSG
         WRITE(15,80) GMSG
         WRITE(15,80) TAU0MSG
C        jgf47.11 Added writing of limits for time varying tau0
         IF ( (TAU0.le.-5.d0).AND.(TAU0.gt.-6.d0) ) THEN
             WRITE(15,80) TAU0LIMMSG
         ENDIF
         WRITE(15,80) DTMSG
         WRITE(15,80) STATMSG
         WRITE(15,80) REFTMSG

!       tcm v49.64.01 No changes needed here for the use of ICE
         IF((NWS.EQ.0).AND.(NRS.GE.1)) WRITE(15,80) RSTIMMSG  ! sb46.28sb03
         IF((NWS.EQ.1).AND.(NRS.GE.1)) WRITE(15,80) RSTIMMSG  ! sb46.28sb03
C     jgfdebug46.02 Added check for NWS=45.
C     jgf46.02 Added NWS=8.
C     jgf46.16 Merged:
C     rjw added NWS=19: asymmetric hurricane wind model v2.0
C     sb46.28sb01 added NWS=12: OWI format
C     jgf50.38.05: added NWS=15: HWind format
         IF ((ABS(NWS).EQ.2).OR.(ABS(NWS).EQ.4).OR.(ABS(NWS).EQ.45).OR.
     &        (ABS(NWS).EQ.5).OR.(ABS(NWS).EQ.6).OR.(ABS(NWS).EQ.8)
     &        .OR.(ABS(NWS).EQ.12).OR.(ABS(NWS).EQ.15)
     &         .OR.(ABS(NWS).EQ.19).OR.(ABS(NWS).EQ.29))THEN
            WRITE(15,80) WSMSG1
         ENDIF
         IF (NWS.EQ.3) THEN
            WRITE(15,80) WSMSG1
            WRITE(15,80) WSMSG2
         ENDIF

         WRITE(15,80) RNDAYMSG
         WRITE(15,80) DRAMPMSG
         WRITE(15,80) COEFMSG
         WRITE(15,80) H0MSG
         WRITE(15,80) SLMSG
         WRITE(15,80) TAUMSG
         WRITE(15,80) ESLMSG
         WRITE(15,80) CORIMSG
         WRITE(15,80) NTIFMSG
         DO I=1,NTIF
            WRITE(15,80)  TIPOTAG(I)
            WRITE(15,80)  TPKMSG(I)
         ENDDO

         WRITE(15,80) NBFRMSG
         DO I=1,NBFR
            WRITE(15,80) BOUNTAG(I)
            WRITE(15,80) AMIGMSG(I)
         ENDDO
         DO I=1,NBFR
            WRITE(15,80) ALPHA1(I)
            DO J=1,NETAP(IPROC)
               WRITE(15,80) EMOMSG(I,OBNODE_LG(J,IPROC))
            ENDDO
         ENDDO

         WRITE(15,80) ANGMSG
C
C--If there were any normal flow boundaries local to PE, process them
C
C         PRINT *, NFFRMSG
C         PRINT *, "NFLUXFP(",IPROC,") = ", NFLUXFP(IPROC)

         IF (NFLUXFP(IPROC).EQ.1) THEN
C
            NFLBNP = 0
            DO I=1, NFLBN
               INDX = FLBN(I)
               DO J=1, ITOTPROC(INDX)
                  IPROC2 = IMAP_NOD_GL2(2*(J-1)+1,INDX)
                  IF (IPROC.EQ.IPROC2) THEN
                     NFLBNP = NFLBNP + 1
                     FLBNXP(NFLBNP) = FLBNX(I)
                  ENDIF
               ENDDO
            ENDDO
C
            WRITE(15,80) NFFRMSG
            IF ((NFFR.NE.0).AND.(NFFR.NE.-1)) THEN
               DO I=1,NFFR
                  WRITE(15,80) FBOUNTAG(I)
                  WRITE(15,80) FREQMSG(I)
               ENDDO
               DO I=1,NFFR
                  WRITE(15,80) ALPHA2(I)
                  DO J=1,NFLBNP
                     WRITE(15,80) QNMSG(I,FLBNXP(J))
cdebug               print *, "PE=",IPROC," FLUXNODE=",FLBNXP(J)
                  ENDDO
               ENDDO
            ENDIF
C
         ENDIF
C
C--Write Local Elevation Station Info:
C--Create Local-to-Global element "ownership" of an elevation station
C
C     WRITE(15,80) STAEMSG !jgf45.07 we may have changed NOUTE in adcprep
         WRITE(15,*) NOUTE,TOUTSE,TOUTFE,NSPOOLE
C
         NSTAEP(IPROC) = 0
         DO K = 1,NSTAE
            DO J=1,NELP(IPROC)
               INDX = abs(IMAP_EL_LG(J,IPROC))
               IF (INDX.EQ.NNSEG(K)) THEN
                  NSTAEP(IPROC) = NSTAEP(IPROC) + 1
                  KK = K
                  if (STAE_SHARE(K) > -1) KK = -K
                  IMAP_STAE_LG(NSTAEP(IPROC),IPROC) = KK
                  STAE_SHARE(K) = IPROC
               ENDIF
            ENDDO
         ENDDO
         NSTAE_MAX = MAX(NSTAEP(IPROC),NSTAE_MAX)
C
         CALL INSERT(NSTAEMSG,OUTMSG,NSTAEP(IPROC),1)
         WRITE(15,80) OUTMSG
C
         DO K=1,NSTAEP(IPROC)
            INDX = abs(IMAP_STAE_LG(K,IPROC))
            WRITE(15,80) STAELOC(INDX)
         ENDDO
C
C--Write Local Velocity Station Info:
C--Create Local-to-Global element "ownership" of an velocity station
C
         WRITE(15,*) NOUTV,TOUTSV,TOUTFV,NSPOOLV
C
         NSTAVP(IPROC) = 0
         DO K = 1,NSTAV
            DO J=1,NELP(IPROC)
               INDX = abs(IMAP_EL_LG(J,IPROC))
               IF (INDX.EQ.NNSVG(K)) THEN
                  NSTAVP(IPROC) = NSTAVP(IPROC) + 1
                  KK = K
                  if (STAV_SHARE(K) > -1) KK = -K
                  IMAP_STAV_LG(NSTAVP(IPROC),IPROC) = KK
                  STAV_SHARE(K) = IPROC
               ENDIF
            ENDDO
         ENDDO
         NSTAV_MAX = MAX(NSTAVP(IPROC),NSTAV_MAX)
C
         CALL INSERT(NSTAVMSG,OUTMSG,NSTAVP(IPROC),1)
         WRITE(15,80) OUTMSG
C
         DO K=1,NSTAVP(IPROC)
            INDX = abs(IMAP_STAV_LG(K,IPROC))
            WRITE(15,80) STAVLOC(INDX)
         ENDDO
C
C--If IM=10 Write Concentration Station Info:
C--Create Local-to-Global element "ownership" of an concentration station
C
         NSTACP(IPROC) = 0
         IF (C2D_PTrans.or.C3D_PTrans) THEN !jgf46.28
C
C     WRITE(15,80) STACMSG   !jgf45.07 we may have changed NOUTC in adcprep
            WRITE(15,*) NOUTC,TOUTSC,TOUTFC,NSPOOLC
C
            DO K = 1,NSTAC
               DO J=1,NELP(IPROC)
                  INDX = abs(IMAP_EL_LG(J,IPROC))
                  IF (INDX.EQ.NNSCG(K)) THEN
                     NSTACP(IPROC) = NSTACP(IPROC) + 1
                     KK = K
                     if (STAC_SHARE(K) > -1) KK = -K
                     IMAP_STAC_LG(NSTACP(IPROC),IPROC) = KK
                     STAC_SHARE(K) = IPROC
                  ENDIF
               ENDDO
            ENDDO
            NSTAC_MAX = MAX(NSTACP(IPROC),NSTAC_MAX)
C
            CALL INSERT(NSTACMSG,OUTMSG,NSTACP(IPROC),1)
            WRITE(15,80) OUTMSG
C
            DO K=1,NSTACP(IPROC)
               INDX = abs(IMAP_STAC_LG(K,IPROC))
               WRITE(15,80) STACLOC(INDX)
            ENDDO
C
         ENDIF
C
C--Write Local Meterological Station Info:
C--Create Local-to-Global element "ownership" of an elevation station
C
         NSTAMP(IPROC) = 0
         IF (NWS.NE.0) THEN
            WRITE(15,*) NOUTM,TOUTSM,TOUTFM,NSPOOLM
            DO K = 1,NSTAM
               DO J=1,NELP(IPROC)
                  INDX = abs(IMAP_EL_LG(J,IPROC))
                  IF (INDX.EQ.NNSMG(K)) THEN
                     NSTAMP(IPROC) = NSTAMP(IPROC) + 1
                     KK = K
                     if (STAM_SHARE(K) > -1) KK = -K
                     IMAP_STAM_LG(NSTAMP(IPROC),IPROC) = KK
                     STAM_SHARE(K) = IPROC
                  ENDIF
               ENDDO
            ENDDO
            NSTAM_MAX = MAX(NSTAMP(IPROC),NSTAM_MAX)
C
            CALL INSERT(NSTAMMSG,OUTMSG,NSTAMP(IPROC),1)
            WRITE(15,80) OUTMSG
C
            DO K=1,NSTAMP(IPROC)
               INDX = abs(IMAP_STAM_LG(K,IPROC))
               WRITE(15,80) STAMLOC(INDX)
            ENDDO
         ENDIF

C
C--Write Local Elevation Data Output Info
C
C      WRITE(15,80) OUTGEMSG !jgf45.07 we may have changed NOUTGE in adcprep
         WRITE(15,*) NOUTGE,TOUTSGE,TOUTFGE,NSPOOLGE
C
C--Write Local Velocity Data Output Info
C
C     WRITE(15,80) OUTGVMSG !jgf45.07 we may have changed NOUTGV in adcprep
         WRITE(15,*) NOUTGV,TOUTSGV,TOUTFGV,NSPOOLGV
C
C     jgf45.07 write subdomain concentration data output info if necessary
C
         IF (IM.EQ.10) WRITE(15,*) NOUTGC,TOUTSGC,TOUTFGC,NSPOOLGC
C
C--Write Local Wind Velocity Data Output Info ( added 4/16/98 vjp )
C
C     jgf45.07 we may have changed NOUTGW in adcprep
C     IF (NWS.NE.0) WRITE(15,80) OUTGWMSG
         IF (NWS.NE.0) WRITE(15,*) NOUTGW,TOUTSGW,TOUTFGW,NSPOOLGW
C
C--Write Harmonic Analysis Data
C
         WRITE(15,80) HARFRMSG
         DO I=1,NHARFR
            WRITE(15,'(A10)') NAMEFR(I)
            WRITE(15,80) HAFREMSG(I)
c           WRITE(15,*) HAFREQ(I),HAFF(I),HAFACE(I)
         ENDDO
C
         WRITE(15,80) HARPARMSG
         WRITE(15,80) OUTHARMSG
C
C--Write Hot Start Info
C
         WRITE(15,80) HSTARMSG
C
C--Write Solver Info
C
         WRITE(15,80) SOLVMSG
C
C--Write 3DVS Info
C
         IF(C3DVS) THEN
            CALL PREP15_3DVS(IPROC)
c        ELSEIF(C3DDSS) THEN
c           CALL PREP15_3DDSS(IPROC)
         ENDIF
C
C     jgf48.03 Write netCDF metadata, if necessary
         IF (useNetCDF.eqv..true.) THEN
            !write(*,*) 'writing netcdf metadata to fort.15' ! jgfdebug
            WRITE(15,*) trim(adjustl(title))
            WRITE(15,*) trim(adjustl(institution))
            WRITE(15,*) trim(adjustl(source))
            WRITE(15,*) trim(adjustl(history))
            WRITE(15,*) trim(adjustl(references))
            WRITE(15,*) trim(adjustl(comments))
            WRITE(15,*) trim(adjustl(host))
            WRITE(15,*) trim(adjustl(convention))
            WRITE(15,*) trim(adjustl(contact))
            WRITE(15,*) trim(adjustl(base_date))
         ELSE
            !write(*,*) 'not writing netcdf metadata' ! jgfdebug
         ENDIF

c...     tcm v50.66.02 additions for time varying bathymetry
         IF (FOUND_TBC_NML) then  !If there was a namelist in the original fort.15 put it in the decomp 15's
!         IF (NDDT.NE.0) THEN
            write(15,*) '! -- Begin Time Varying Bathymetry Inputs --'
            write(15,TimeBathyControl)
            write(15,*) '! -- End Time Varying Bathymetry Inputs --'
         ENDIF
C
#if defined CSWAN || defined ADCSWAN
         write(15,*) '! -- Begin SWAN Output Control Namelist --'
         write(15,SWANOutputControl)
         write(15,*) '! -- End SWAN Output Control Namelist --'
#endif

! tcm v50.79 Changed so that metControl namelist is only written if it was found in the
! original fort.15 file.  Also changed the single line write, which is missing some commas
! to a multiple line write.  The single line write was causing problems on
! some compilers because the character DragLawString could end up being written
! on multiple lines and this caused issues.  This section
! should only be written if there was a namelist in the original fort.15.
!
         if (found_metCon_nml) then  !metControl namelist was found so write it in the parsed files
            write(15,*) '! -- Begin Met Control Namelist --'
!         write(15,*) "&metControl WindDragLimit=",WindDragLimit,
!     &         " DragLawString='",DragLawString,"' rhoAir=",rhoAir," /"

            write(15,*) "&metControl "
            write(15,*) "    WindDragLimit=",Winddraglimit,","
            write(15,*) "    DragLawString='",Draglawstring,"',"
            write(15,*) "    rhoAir=",rhoAir
            write(15,*) "/"
            write(15,*) '! -- End Met Control Namelist --'
         endif

         CLOSE(15)
C
 1000 CONTINUE

      IF(C3DVS.and.(IEVC.EQ.0)) THEN
         DEALLOCATE ( EVTot )
         nbytes = 8*nfen
         call memory_dealloc(nbytes)
      ENDIF

C
C--Print Summary of Stations
C
      print *, " "
      print *, "Station Data"
      print *, "DOMAIN      NSTAE   NSTAV    NSTAC    NSTAM"
      WRITE(*,92)  "GLOBAL",NSTAE,NSTAV,NSTAC,NSTAM
      DO IPROC=1, NPROC
         PE(1:6) = 'PE0000'
         CALL IWRITE(PE,3,6,IPROC-1)
         WRITE(*,92)  PE,NSTAEP(IPROC),NSTAVP(IPROC),
     .        NSTACP(IPROC),NSTAMP(IPROC)
      ENDDO
C
      RETURN
 80   FORMAT(A80)
 92   FORMAT(1X,A6,4I8)
      END SUBROUTINE PREP15

C
C---------------------------------------------------------------------------C
C                     (  Serial Version  6/24/02  )                         C
C  This routine writes the 3DVS info in the Local Input file "fort.15" file C
C  for each subdomain using the domain decomposition of the ADCIRC grid     C
C  created by the routine DECOMP.                                           C
C                                                                           C
C  The Decomposition Variables are defined in the include file adcprep.inc  C
C  This version is compatible with ADCIRC version 41.11a                    C
C
C     jgf45.11 Updated to handle new format of 3D input files with stations
C     defined by coordinates rather than node numbers.
C---------------------------------------------------------------------------C
      SUBROUTINE PREP15_3DVS(IPROC)
      USE PRE_GLOBAL
      USE GLOBAL_3DVS, ONLY : SIGMA
      IMPLICIT NONE
      INTEGER N          ! vertical grid layer counter
      INTEGER IPROC      ! subdomain counter
      INTEGER SDStation  ! subdomain station
      INTEGER FDStation  ! full domain station
      INTEGER SDEle      ! subdomain station element
      INTEGER FDEle      ! full domain element

c     jgf45.10 removed IDIAG
      WRITE(15,80) IDENMSG
      WRITE(15,80) SLIPMSG
      WRITE(15,80) Z0MSG
      WRITE(15,80) ALPMSG
      WRITE(15,80) FEMSG
C     jgf45.12 Added code to record thicknesses of vertical grid layers,
C     if necessary.
      IF(IGC.EQ.0) THEN
         DO N=1,NFEN
            WRITE(15,*) Sigma(N)
         ENDDO
      ENDIF
      WRITE(15,80) EVCMSG
C     jgf45.12 Add code to record vertical eddy viscosity profile.
      IF(IEVC.EQ.0) THEN
         DO N=1,NFEN
            WRITE(15,*) EVTot(N)
         ENDDO
      ENDIF
      IF((IEVC.EQ.50).OR.(IEVC.EQ.51)) WRITE(15,80) THETAMSG
C     -------------------------------------------------------------
C     jgf45.11 Create mapping from full domain 3D density station
C     elements to corresponding elements in subdomains. Write out
C     subdomain station locations to fort.15 file.
C     -------------------------------------------------------------
      WRITE(15,*) I3DSD,TO3DSDS,TO3DSDF,NSPO3DSD
C   kmd48.33bc changed
      IF(NSTA3DD.NE.0) THEN
         NNSTA3DDP(IPROC) = 0
         DO FDStation = 1, NSTA3DD
            DO SDEle = 1, NELP(IPROC)
               FDEle = abs(IMAP_EL_LG(SDEle,IPROC))
               IF ( FDEle .eq. NNS3DDG(FDStation) ) THEN
                  NNSTA3DDP(IPROC) = NNSTA3DDP(IPROC) + 1
                  IMAP_STA3DD_LG(NNSTA3DDP(IPROC),IPROC) = FDStation
               ENDIF
            END DO
         END DO
         WRITE(15,*) NNSTA3DDP(IPROC)
         DO SDStation = 1, NNSTA3DDP(IPROC)
            FDStation = IMAP_STA3DD_LG(SDStation,IPROC)
            WRITE(15,80) STA3DDLOC(FDStation)
         ENDDO
      ELSE
         WRITE(15,80) NSTA3DDMSG
      ENDIF
C     -------------------------------------------------------------
C     jgf45.11 Create mapping from full domain 3D velocity station
C     elements to corresponding elements in subdomains. Write out
C     velocity subdomain station locations to fort.15 file.
C     -------------------------------------------------------------
      WRITE(15,*) I3DSV,TO3DSVS,TO3DSVF,NSPO3DSV
C   kmd48.33bc changed
      IF(NSTA3DV.NE.0) THEN
         NNSTA3DVP(IPROC) = 0
         DO FDStation = 1, NSTA3DV
            DO SDEle = 1, NELP(IPROC)
               FDEle = abs(IMAP_EL_LG(SDEle,IPROC))
               IF ( FDEle .eq. NNS3DVG(FDStation) ) THEN
                  NNSTA3DVP(IPROC) = NNSTA3DVP(IPROC) + 1
                  IMAP_STA3DV_LG(NNSTA3DVP(IPROC),IPROC) = FDStation
               ENDIF
            END DO
         END DO
         WRITE(15,*) NNSTA3DVP(IPROC)
         DO SDStation = 1, NNSTA3DVP(IPROC)
            FDStation = IMAP_STA3DV_LG(SDStation,IPROC)
            WRITE(15,80) STA3DVLOC(FDStation)
         ENDDO
      ELSE
         WRITE(15,80) NSTA3DVMSG
      ENDIF
C     -------------------------------------------------------------
C     jgf45.11 Create mapping from full domain 3D turbulence station
C     elements to corresponding elements in subdomains. Write out
C     turbulence subdomain station locations to fort.15 file.
C     -------------------------------------------------------------
      WRITE(15,*) I3DST,TO3DSTS,TO3DSTF,NSPO3DST
C   kmd48.33bc changed
      IF(NSTA3DT.NE.0) THEN
         NNSTA3DTP(IPROC) = 0
         DO FDStation = 1, NSTA3DT
            DO SDEle = 1, NELP(IPROC)
               FDEle = abs(IMAP_EL_LG(SDEle,IPROC))
               IF ( FDEle .eq. NNS3DTG(FDStation) ) THEN
                  NNSTA3DTP(IPROC) = NNSTA3DTP(IPROC) + 1
                  IMAP_STA3DT_LG(NNSTA3DTP(IPROC),IPROC) = FDStation
               ENDIF
            END DO
         END DO
         WRITE(15,*) NNSTA3DTP(IPROC)
         DO SDStation = 1, NNSTA3DTP(IPROC)
            FDStation = IMAP_STA3DT_LG(SDStation,IPROC)
            WRITE(15,80) STA3DTLOC(FDStation)
         ENDDO
      ELSE
         WRITE(15,80) NSTA3DTMSG
      ENDIF

      WRITE(15,80) DGDMSG
      WRITE(15,80) DGVMSG
      WRITE(15,80) DGTMSG

C    kmd48.33bc add 3D boundary condition information
      IF (CBAROCLINIC) THEN
        WRITE(15,80) RESBCFLAGMSG
        IF (RES_BC_FLAG.NE.0) THEN
          IF (NOPEP(IPROC).GT.0) THEN
             WRITE(15,80) BCTIMEMSG
             WRITE(15,80) BCSTATMSG
          END IF
          IF (BCFLAG_TEMP.NE.0) THEN
             WRITE(15,80) TBCTIMEMSG
          END IF
        END IF
      END IF

      IF (CBAROCLINIC) THEN
        WRITE(15,80) SPONGEDISTMSG
        WRITE(15,80) EqnstateMSG
      END IF

C
C     jgf45.12: Write out the parameters for the transport equation, if
C     necessary.
      IF (C3D_BTrans) THEN
C     Lateral and vertical diffusion coefficients.
         WRITE(15,*) NLSD, NVSD
         WRITE(15,*) NLTD, NVTD
C     Time stepping coefficient for the transport equation terms.
         WRITE(15,*) ALP4
C   kmd48.33 took out as it is no longer needed with new heat flux boundary conditions
C     Temperature boundary condition file type, if necessary
!         IF ( IDEN .eq. 3 .or. IDEN .eq. 4 ) THEN
!            WRITE(15,*) NTF
!         ENDIF
      ENDIF

      RETURN
  80  FORMAT(A80)
  81  FORMAT(I8,2E15.8,2I8,A32)
  82  FORMAT(500I8)
      END
C-----------------------------------------------------------------------
C     End of subroutine PREP15_3DVS
C-----------------------------------------------------------------------


      SUBROUTINE PREP18()
      USE PRE_GLOBAL
      use memory_usage
C
C---------------------------------------------------------------------------C
C                     (  Serial Version  6/10/2011  )                       C
C  This Routine writes a message-passing file "fort.18" for each subdomain  C
C  of the domain decomposition created by DECOMP.                           C
C                                                                           C
C  The Decomposition Variables are defined in the include file adcprep.inc  C
C  This version is compatible with ADCIRC version 50.21                     C
C                                                                           C
C  tcm V50.21 -- Changed all I8 formats to I12                              C
C---------------------------------------------------------------------------C
C
      IMPLICIT NONE
      integer :: nbytes = 0
      INTEGER N1, N2, N3, KMIN
      INTEGER I,J,K,M,ITEMP,IPR,IPR1
      INTEGER INDX,ITOT,IEL,IELG,ILNODE,IPROC,ITYPE
      INTEGER,ALLOCATABLE :: RES_NODE(:)
      CHARACTER LOCFN*14,PE*6
C
C  Allocate local arrays
C
      ALLOCATE ( RES_NODE(MNPP) )
      nbytes = 4*mnpp
      call memory_alloc(nbytes)
C
C--Write Message-Passing File for each PE
C
      DO 1000 I = 1,NPROC
C
         LOCFN(1:14) = 'PE0000/fort.18'
         CALL IWRITE(LOCFN,3,6,I-1)
         OPEN (18,FILE=LOCFN)

         write(18,3050) FileFmtVersion, 0, 0

!vjp 9/17/06
C--Write the Global indexes of all local elements in local element order

Casey 100209: Changed I8 to I12.
         WRITE(18,3000) NELG, MNEP, NELP(I)  ! number of Global elements
         DO J = 1,NELP(I)
            INDX = IMAP_EL_LG(J,I)
            WRITE(18,'(I12)') INDX        ! Global index of local element
         ENDDO

C--Write the Global indexes of all local nodes in local node order
C  write global index as positive if a resident node and negative
C  if a ghost node

         WRITE(18,3001) NNODG, MNPP, NNODP(I)   ! number of Global nodes
         ITOT = 0
         DO J = 1,NNODP(I)
            INDX = IMAP_NOD_LG(J,I)
            IPR = IMAP_NOD_GL(1,INDX)
            IF (IPR.EQ.I)THEN
               ITOT = ITOT + 1
               RES_NODE(ITOT) = J
               WRITE(18,'(I12)') INDX        ! Global index of resident node
            ELSE
               WRITE(18,'(I12)') -1*INDX     ! Global index of ghost node
            ENDIF
         ENDDO
         IF (ITOT.NE.NOD_RES_TOT(I)) STOP 'ERROR IN # OF RES. NODES'

C--Write local normal flow boundary flag
C--vjp This info is used only for relocalizing fort.15
         WRITE(18,3002) NFLUXFP(I) ! normal flow b.c. flag for subdomain

C--Write global and local total number of elevation boundary nodes
C--vjp This info is used only for relocalizing fort.15
         WRITE(18,3003) NETA, NETA_MAX, NETAP(I) ! number of global elevation b.c. nodes
         DO J = 1,NETAP(I)
            INDX = OBNODE_LG(J,I)
            WRITE(18,'(I12)') INDX           ! Global open boundary node index
         ENDDO

C--Write the Global indexes of all Elevation Stations in local node order
C  write global index as positive if a resident node and negative
C  if a ghost node

         WRITE(18,3004) NSTAE, NSTAE_MAX, NSTAEP(I) ! number of Global Elevation Stations
         DO J = 1,NSTAEP(I)
            INDX = IMAP_STAE_LG(J,I)
            WRITE(18,'(I12)') INDX           ! Global station number
         ENDDO

C--Write the Global indexes of all Velocity Stations in local node order
C  write global index as positive if a resident node and negative
C  if a ghost node

         WRITE(18,3005) NSTAV, NSTAV_MAX, NSTAVP(I) ! number of Global Velocity Stations
         DO J = 1,NSTAVP(I)
            INDX = IMAP_STAV_LG(J,I)
            WRITE(18,'(I12)') INDX           ! Global station number
         ENDDO

C--Write the Global indexes of all Elevation Stations in local node order
C  write global index as positive if a resident node and negative
C  if a ghost node

         WRITE(18,3006) NSTAM, NSTAM_MAX, NSTAMP(I) ! number of Global Meteorlogical Stations
         DO J = 1,NSTAMP(I)
            INDX = IMAP_STAM_LG(J,I)
            WRITE(18,'(I12)') INDX           ! Global station number
         ENDDO

C--Write the Global indexes of all Concentration Stations in local node order
C  write global index as positive if a resident node and negative
C  if a ghost node

         WRITE(18,3007) NSTAC, NSTAC_MAX, NSTACP(I) ! number of Global Concentration Stations
         DO J = 1,NSTACP(I)
            INDX = IMAP_STAC_LG(J,I)
            WRITE(18,'(I12)') INDX           ! Global station number
         ENDDO

C---------------------------------------------------------------------------------
C---------------------------------------------------------------------------------
C---------------------------------------------------------------------------------
C
C--Write the Resident Node List
C
         WRITE(18,3010) (I-1),NOD_RES_TOT(I)
         WRITE(18,1130) (RES_NODE(J),J=1,ITOT)
C
C--Write the Number of Communicating PEs
C
         WRITE(18,3020) NUM_COMM_PE(I)
C
C--Write the Receive List
C
         DO J = 1,NUM_COMM_PE(I)
            IPR = COMM_PE_NUM(J,I)
            IRECV_TOT(J,I) = 0
            DO K = 1,NNODP(I)
               INDX = IMAP_NOD_LG(K,I)
               IF (IMAP_NOD_GL(1,INDX).EQ.IPR) THEN
                  IRECV_TOT(J,I) = IRECV_TOT(J,I) + 1
                  IRECV(IRECV_TOT(J,I)) = K
c uncomment next line and comment preceding line for debugging
c                 IRECV(IRECV_TOT(J,I)) = INDX
               ENDIF
            ENDDO
            WRITE(18,3030) (IPR-1), IRECV_TOT(J,I)
            WRITE(18,1130) (IRECV(K),K=1,IRECV_TOT(J,I))
         ENDDO
C
C--write the send list
C
         DO J = 1,NUM_COMM_PE(I)
            IPR = COMM_PE_NUM(J,I)
            ISEND_TOT(J,I) = 0
            DO K = 1,NNODP(IPR)
               INDX = IMAP_NOD_LG(K,IPR)
               IF (IMAP_NOD_GL(1,INDX).EQ.I) THEN
                  ISEND_TOT(J,I) = ISEND_TOT(J,I) + 1
                  ISEND(ISEND_TOT(J,I)) = IMAP_NOD_GL(2,INDX)
c uncomment next line and comment preceding line for debugging
c                 ISEND(ISEND_TOT(J,I)) = INDX
               ENDIF
            ENDDO
            WRITE(18,3040)  IPR-1, ISEND_TOT(J,I)
            WRITE(18,1130) (ISEND(K),K=1,ISEND_TOT(J,I))
         ENDDO
C
         IF (C3D.eqv..true.) THEN
C           jgf49.43.18: Add 3D station mappings from subdomain to fulldomain
C           to accomodate globalio.
C
C           Write the fulldomain station numbers of all 3D density stations
C           in local node order; write the fulldomain station number as positive
C           for resident stations and negative for ghost stations.
            WRITE(18,3060) NSTA3DD, MAXVAL(NNSTA3DDP), NNSTA3DDP(I)
            DO J=1,NNSTA3DDP(I)
               WRITE(18,1131) IMAP_STA3DD_LG(J,I)
            END DO
C           3D velocity stations
            WRITE(18,3061) NSTA3DV, MAXVAL(NNSTA3DVP), NNSTA3DVP(I)
            DO J=1,NNSTA3DVP(I)
               WRITE(18,'(I12)') IMAP_STA3DV_LG(J,I)
            END DO
C           3D turbulence stations
            WRITE(18,3062) NSTA3DT, MAXVAL(NNSTA3DTP), NNSTA3DTP(I)
            DO J=1,NNSTA3DTP(I)
               WRITE(18,'(I12)') IMAP_STA3DT_LG(J,I)
            ENDDO
         ENDIF
C
         CLOSE(18)
C
1000  CONTINUE
C
C--Compute the surface to volume ratio (in %)
C
      DO I = 1,NPROC
         ITOT = 0
         DO J = 1,NUM_COMM_PE(I)
            ITOT = ITOT + IRECV_TOT(J,I)
         ENDDO
         PROC_SV(I) = (ITOT/REAL(NOD_RES_TOT(I)))*100.0
c        WRITE(6,*) I-1,PROC_SV(I)
      ENDDO
C
      print *, " "
      print *, "Communication Data"
      print *, "DOMAIN  COMM_PE  %(SURF/VOL)"
      print *, "------  -------  -----------"
      DO I=1, NPROC
         PE(1:6) = 'PE0000'
         CALL IWRITE(PE,3,6,I-1)
         WRITE(6,92) PE, NUM_COMM_PE(I),PROC_SV(I)
      ENDDO

      deallocate( res_node )
      nbytes = 4*mnpp
      call memory_dealloc(nbytes)
      call memory_status()
      RETURN
C
  92  FORMAT(1X,A6,2X,I7,2X,F8.2)
1130  FORMAT(8X,6I12) !(8X,9I8)
1131  FORMAT(:,I12)
Casey 100209: Changed I8 to I12 through this section.
3000  FORMAT('NELG    ',3I12)
3001  FORMAT('NNODG   ',3I12)
3002  FORMAT('NFLUXF  ',I12)
3003  FORMAT('NETA    ',3I12)
3004  FORMAT('NSTAE   ',3I12)
3005  FORMAT('NSTAV   ',3I12)
3006  FORMAT('NSTAM   ',3I12)
3007  FORMAT('NSTAC   ',3I12)
3010  FORMAT('RES NODE',2I12)
3020  FORMAT('COMM PE ',2I12)
3030  FORMAT('RECV PE ',2I12)
3040  FORMAT('SEND PE ',2I12)
3050  FORMAT('FileFmt ',3I12)
3060  FORMAT('NSTA3DD ',3I12)
3061  FORMAT('NSTA3DV ',3I12)
3062  FORMAT('NSTA3DT ',3I12)
      END SUBROUTINE PREP18


      SUBROUTINE PREP19()
      USE PRE_GLOBAL
      use memory_usage
C
C---------------------------------------------------------------------------C
C                     (  Serial Version  2/28/98  )                         C
C  This routine writes a Local "Aperiodic Elevation Boundary Condtions"     C
C  (fort.19) file for each subdomain using the domain decomposition of      C
C  the ADCIRC grid created by the routine DECOMP.                           C
C                                                                           C
C  The Decomposition Variables are defined in the include file adcprep.inc  C
C  This version is compatible with ADCIRC version 34.03                     C
C
C     jgf45.12 Added subroutine call to open files.
C                                                                           C
C---------------------------------------------------------------------------C
C
      IMPLICIT NONE
      integer :: nbytes = 0
      INTEGER I,J,IPROC
      INTEGER SDU(NPROC) ! subdomain unit numbers
      LOGICAL Success    ! .true. if files opened without errors
      CHARACTER*40  ETIMINC,ESBINP
      CHARACTER*40,ALLOCATABLE :: ESBIN(:)

C
C--Enter, Locate, Open, and Read the ADCIRC UNIT 19
C  Global Aperiodic Elevation Boundary Conditions file
C
C     Open full domain and subdomain fort.19 files
      CALL OpenPrepFiles(19, 'aperiodic elevation boundary  ',
     &     1, nproc, SDU, Success)
      IF (.not.Success) THEN
         WRITE(*,*) 'WARNING: Unit 19 files not preprocessed.'
         RETURN ! note early return
      ENDIF
C
C--Allocate local arrays
C
      ALLOCATE ( ESBIN(MNETA) )
      nbytes = 8*mneta
      call memory_alloc(nbytes)
C
      READ(19,40) ETIMINC
      DO IPROC = 1,NPROC
         WRITE(SDU(IPROC),40)  ETIMINC
      ENDDO
C
C--While ( NOT EOF ) Read NETA BCs from Global File
C
1000  CONTINUE
      DO I=1, NETA
         READ(19,40,END=9999)  ESBIN(I)
      ENDDO
C
      DO IPROC= 1,NPROC
         DO I=1, NETAP(IPROC)
            ESBINP = ESBIN(OBNODE_LG(I,IPROC))
            WRITE(SDU(IPROC),40) ESBINP
         ENDDO
      ENDDO
C
      GO TO 1000
C
C--Close Global file and all the Local Files
C
 9999 CLOSE (19)
      DO IPROC=1, NPROC
         CLOSE (SDU(IPROC))
      ENDDO

      deallocate( esbin )
      nbytes = 8*mneta
      call memory_dealloc(nbytes)
      call memory_status()
      RETURN
  40  FORMAT(A40)
      END


C---------------------------------------------------------------------------
C                S U B R O U T I N E   P R E P 2 0
C---------------------------------------------------------------------------
C
C     jgf45.12 This subroutine will break up the full domain aperiodic
C     flux boundaries into subdomains using the domain decomposition of
C     the ADCIRC grid created by the routine DECOMP.
C
C     -Written by MEB 04/01/04
C     -Added by jgf to 45.06 10/07/2005
C     -jgf45.12 Rewritten to correct bugs in subdomain fort.20
C     formatting as well as the erroneous use of the GL mapping instead
C     of GL2. Also added subroutine call to open files.
C
C---------------------------------------------------------------------------
      SUBROUTINE PREP20()
C---------------------------------------------------------------------------
      USE PRE_GLOBAL
      use memory_usage
      IMPLICIT NONE
      integer :: nbytes = 0
      INTEGER IPROC
      INTEGER INDEX14, I
      REAL(SZ)  FLUX_INC, FLUX_VAL
      INTEGER SDU(NPROC)  ! subdomain unit numbers
      LOGICAL Success     ! .true. if all files open without error
      INTEGER INDX ! full domain node number for a flow boundary node
      INTEGER J     ! counter for subdomains that corrsp. to a single f.d. node
      INTEGER IPROC2! PE of a subdomain that matches a single full domain node

C     Open full domain and subdomain fort.20 files
      CALL OpenPrepFiles(20, 'aperiodic flux boundary       ',
     &     1, nproc, SDU, Success)
      IF (.not.Success) THEN
         WRITE(*,*) 'WARNING: Unit 20 files not preprocessed.'
         RETURN ! note early return
      ENDIF
c
c     Write Increment into all flux files
c
      READ(20,*) FLUX_INC
      DO IPROC=1,NPROC
         WRITE(SDU(IPROC),*) FLUX_INC
      ENDDO
C
C     jgf45.12 Write each full domain nodal flux value into each of the
C     subdomains that that full domain node maps to. The full domain
C     node may map to more than one subdomain node if it falls on a
C     boundary between subdomains (ghost nodes).
C
 33   DO I=1, EXIST_FLUX            ! loop through full domain flow nodes
         READ(20,*,END=40) FLUX_VAL ! get a flo val for this f.d. flow node
         INDX = FLUX14_ARY(I)      ! get full domain flow node number
         DO J=1, ITOTPROC(INDX)    ! loop over subdomains for 1 f.d. node
            IPROC2 = IMAP_NOD_GL2(2*(J-1)+1,INDX) ! find next subdomain
            DO IPROC=1, NPROC
               IF (IPROC.EQ.IPROC2) THEN ! full domain node maps to this s.d.
                  WRITE(SDU(IPROC),50) FLUX_VAL
               ENDIF
            END DO
         END DO
      END DO
      GOTO 33
 40   CLOSE (20)
      DO IPROC=1, NPROC
         CLOSE (SDU(IPROC))
      ENDDO

      IF (allocated(FLUX14_ARY)) then
        DEALLOCATE (FLUX14_ARY)
        nbytes = 4*exist_flux
        call memory_dealloc(nbytes)
      ENDIF
      call memory_status()
      return
 50   FORMAT (F16.8,1x,I6,1x,I6,1x,I6)
c----------------------------------------------------------------------------
      END SUBROUTINE PREP20
c----------------------------------------------------------------------------

C---------------------------------------------------------------------------
C                S U B R O U T I N E   P R E P 8 8
C---------------------------------------------------------------------------
C
C     kmd49 This subroutine will break up the full domain elevation
C     changes due to the river boundary information being above mean
C     sea level. It writes the fort.88 file into each subdomain using
C     the ADCIRC grid created by the routine DECOMP.
C
C     - added as part of Evan's changes for rivers above MSL.
C
C     TCM v50.99.14 -- Added the decomposition for a max of 256 subdomains
C     at a time ... some platforms/compilers limit the number of files that
C     can be open at any one time.
C
C---------------------------------------------------------------------------

      SUBROUTINE PREP88()
      USE PRE_GLOBAL
      use memory_usage
C
      IMPLICIT NONE
      integer :: nbytes = 0
      INTEGER I,J,IPROC
      INTEGER SDU(NPROC) ! subdomain unit numbers
      LOGICAL Success    ! .true. if files opened without errors
      INTEGER :: NODP
      CHARACTER*80 :: et_tempsWSE
      CHARACTER*80,ALLOCATABLE :: et_SWSE(:)
      INTEGER, PARAMETER :: maxOpenFiles = 256
      INTEGER startProc
      INTEGER endProc
      INTEGER deltaProc
C
C     Perform decomposition over range of subdomains.
      startProc = 1
      DO WHILE ( startProc .lt. nproc )
         deltaProc = nproc - startProc
         IF ( deltaProc .gt. maxOpenFiles ) deltaProc = maxOpenFiles
         endProc = startProc + deltaProc

C        Open full domain and subdomain fort.88 files.      
         CALL OpenPrepFiles(88, '  river elevation data  ',
     &     startProc, endProc, SDU, Success)
     
         IF (.not.Success) THEN
            WRITE(*,*) 'WARNING: Unit 88 files not preprocessed.'
            RETURN ! note early return
         ENDIF

         ALLOCATE(et_SWSE(NNODG))
         DO I=1, NNODG
            READ(88,80,END=9999) et_SWSE(I)
         END DO

         DO IPROC=startProc,endProc
            DO I=1, NNODP(IPROC)
               NODP=IMAP_NOD_LG(I,IPROC)
               et_tempsWSE=et_SWSE(NODP)
               WRITE(SDU(IPROC),80) et_tempsWSE
            END DO
         END DO

C        Close full domain and subdomain files
         CLOSE (88)
         DO iproc=startProc, endProc
            CLOSE(sdu(iproc))
         ENDDO
         startProc = endProc + 1
         
         DEALLOCATE(et_SWSE)         
         
         WRITE(6,'(A25,A80)') '     Finished processing ',
     &    'river elevation data'
	 WRITE(6,*) 'for processor range ',startProc,' to ',endProc

      END DO  !Loop over Procs
 80   FORMAT(A80)
9999  CLOSE(88)
      END SUBROUTINE PREP88
C  End SUBROUTINE PREP88


C---------------------------------------------------------------------------
C                S U B R O U T I N E   P R E P 2 2
C---------------------------------------------------------------------------
C
C                     (  Serial Version  2/28/98  )                         C
C  This routine reads a global external meteorology file when NWS=1,+-2,    C
C  +-4,+-5.  In each case it wites a local meteorology file of the same     C
C  format for each subdomain using the domain decomposition of the ADCIRC   C
C  grid created by the routine DECOMP.                                      C
C                                                                           C
C  The Decomposition Variables are defined in the include file adcprep.inc  C
C  This version is compatible with ADCIRC version 43.03                     C
C
C     jgf46.02 Added subroutine call to open prep files; this provides
C     the user with the ability to skip the prepping of wind data files.
C
C     jgfdebug46.02 Added NWS=45 to imitate the behavior of the v42 (IPET)
C     code.
C
C     jgf46.02 Added NWS=8 to copy the wind files for the Holland model
C     into the subdomains.
C
C     tcm_v49.04 Removed NWS=3 and NWS=6 to correspond with the use of a
C      global file rather than local.
C
C---------------------------------------------------------------------------
      SUBROUTINE PREP22()
C---------------------------------------------------------------------------
      USE PRE_GLOBAL
      use memory_usage
      IMPLICIT NONE
      integer :: nbytes = 0
      LOGICAL FOUND,DONE
      INTEGER I,J,IPROC,IPROC2,ILNODE,INDX,NHG,LINDEX
      CHARACTER*80 PBLJAGF
C      CHARACTER FNAME*60,LOCFN*14,CMD1*63,CMD2*7,CMD*70,INLINE*80
      CHARACTER FNAME*60,CMD1*63,CMD2*7,CMD*70
      CHARACTER*170 Line ! line of data from NWS=8 (Holland) file
      CHARACTER*270 Line19 ! line of data from NWS=19 (AsymmHollandv2.0) file
      INTEGER SDU(NPROC)  ! subdomain unit numbers
      LOGICAL Success     ! .true. if all files open without error
      INTEGER,ALLOCATABLE  :: NG(:)
      REAL(SZ),ALLOCATABLE :: WVNXG(:),WVNYG(:),PRG(:)
      REAL(SZ),ALLOCATABLE :: WVNXL(:),WVNYL(:),PRL(:)
      REAL(SZ) U,V,PR
      REAL(SZ) RHOWAT     !jgfdebug46.02
C     jgf48.47 Do the decomposition for a max of 256 subdomains at a
C     time ... some platforms/compilers limit the number of files that
C     can be open at any one time.
      INTEGER, PARAMETER :: maxOpenFiles = 256
      INTEGER startProc
      INTEGER endProc
      INTEGER deltaProc
C
C     Allocate local work arrays
C
      ALLOCATE ( NG(MNWP) )
      nbytes = 4*mnwp
      call memory_alloc(nbytes)
      ALLOCATE ( WVNXG(MNWP),WVNYG(MNWP),PRG(MNWP) )
      nbytes = 24*mnwp
      call memory_alloc(nbytes)
      ALLOCATE ( WVNXL(MNWP),WVNYL(MNWP),PRL(MNWP) )
      nbytes = 24*mnwp
      call memory_alloc(nbytes)
C
C     Perform decomposition over a range of subdomains.
      startProc = 1
      DO WHILE ( startProc .lt. nproc )
         deltaProc = nproc - startProc
         IF ( deltaProc .gt. maxOpenFiles ) deltaProc = maxOpenFiles
         endProc = startProc + deltaProc

C        Open full domain and all subdomain fort.22 files
         CALL OpenPrepFiles(22, 'wind information              ',
     &      startProc, endProc, sdu, success)

         IF (.not.success) THEN
           WRITE(*,*) 'WARNING: Unit 22 files not preprocessed.'
           RETURN ! note early return
         ENDIF
C
C--Branch to Appropriate Code
C
         SELECT CASE(ABS(NWS))
C        -------------
         CASE(1,2,5,7)
C        -------------
C
C     MAIN LOOP FOR NWS = 1, +-2,+-5,+-7
C     (1)  Read a record from Global Wind Stress File
C     (2)  Use Decomp arrarys to Localize record to a subdomain
C     (3)  Write Local Wind Stress record in same format
           DO                     ! loop forever (or until file ends)
              READ(22,*,END=9999)
     &           (NG(I),WVNXG(I),WVNYG(I),PRG(I),I=1,NNODG)
              DO IPROC = STARTPROC, ENDPROC
                 DO I=1, NNODP(IPROC)
                    INDX = IMAP_NOD_LG(I,IPROC)
                    WVNXL(I) = WVNXG(INDX)
                    WVNYL(I) = WVNYG(INDX)
                    PRL(I) = PRG(INDX)
                 ENDDO
                 DO I=1, NNODP(IPROC)
                    WRITE(SDU(IPROC),1100)  I,WVNXL(I),WVNYL(I),PRL(I)
                 ENDDO
              ENDDO
           ENDDO
C
C        -------
         CASE(4)
C        -------
C        MAIN LOOP FOR NWS = +- 4  ( PBL Format )
C        (1)  Read a record from Global Wind Stress File
C        (2)  Use Decomp arrarys to Localize record to a subdomain
C        (3)  Write out in PBL Format on subdomain
C
C--Read a wind field record from the global input file
C
         DO
            READ(22,'(A80)',END=9999) PBLJAGF
            IF(PBLJAGF(2:2).EQ.'#') THEN
               DO IPROC =  STARTPROC,ENDPROC
                  WRITE(SDU(IPROC),1101)
                  WRITE(SDU(IPROC),1100) 1,0.0,0.0,0.0 !victor didn't like this line 27/11/03
               ENDDO
            ELSE
!     vjp 27/11/03
!     rewrote this section to handle ghost-nodes
C              READ(PBLJAGF,'(I8,3E13.5)',END=9999) NHG,U,V,PR
               READ(PBLJAGF,*,END=9999) NHG,U,V,PR
               DO J=1, ITOTPROC(NHG)
                  IPROC  = IMAP_NOD_GL2(2*(J-1)+1,NHG)
                  LINDEX = IMAP_NOD_GL2(2*(J-1)+2,NHG)
                  WRITE(SDU(IPROC),1100) LINDEX,U,V,PR
               ENDDO
            ENDIF
         END DO
C
C        --------
         CASE(45)
C        --------
C        jgf46.02 Convert NWS=4 winds to NWS=5 winds to imitate the Katrina
C        (IPET) version of the code.
C
C-- Read a wind field record from the global input file
C
         DO
            RHOWAT=1000.0d0
            CALL NWS4GET(WVNXG,WVNYG,PRG,G,RHOWAT,NNODG,DONE)
C
            DO IPROC = STARTPROC,ENDPROC
               DO I=1, NNODP(IPROC)
                  INDX = IMAP_NOD_LG(I,IPROC)
                  WVNXL(I) = WVNXG(INDX)
                  WVNYL(I) = WVNYG(INDX)
                  PRL(I) = PRG(INDX)
               ENDDO
               DO I=1, NNODP(IPROC)
                  WRITE(SDU(IPROC),1100)  I,WVNXL(I),WVNYL(I),PRL(I)
               ENDDO
            ENDDO
C--   If reached EOF in NWS4GET last time go close files and return
C
            IF (DONE) GOTO 9999
         ENDDO

C        ------------
         CASE DEFAULT
C        ------------
         print *, "NWS=",NWS," has incorrect value in PREP22"
         RETURN

         END SELECT

C
C--Close Global file and all the Local Files
C
 9999    CLOSE (22)
         DO IPROC=STARTPROC, ENDPROC
            CLOSE (SDU(IPROC))
         ENDDO
         startProc=endProc+1
      ENDDO


      DEALLOCATE ( NG,  WVNXG, WVNYG, PRG )
      DEALLOCATE ( WVNXL, WVNYL, PRL )
      nbytes = 52*mnwp
      call memory_dealloc(nbytes)
      call memory_status()
      RETURN
  60  FORMAT(A60)
 170  FORMAT(A170)
 270  FORMAT(A270)
 1010 FORMAT(' File ',A60,/,' WAS NOT FOUND!  Try again',/)
 1011 FORMAT(' File ',A60,/,' WAS FOUND!  Opening & Processing file',/)
 1100 FORMAT(I8,3E13.5)
 1101 FORMAT(' #')
c----------------------------------------------------------------------------
      END SUBROUTINE PREP22
c----------------------------------------------------------------------------



      SUBROUTINE PREP23()
      USE PRE_GLOBAL
C
C---------------------------------------------------------------------------C
C                           (  add MEB 03/04/03  )                          C
C  This routine writes a Local Input file "fort.23" file for each subdomain C
C  using the domain decomposition of the ADCIRC grid created by the routine C
C  DECOMP.                                                                  C
C                                                                           C
C  The Decomposition Variables are defined in the include file adcprep.inc  C
C  This version is compatible with ADCIRC version 34.03                     C
C                                                                           C
C---------------------------------------------------------------------------C
C
      IMPLICIT NONE
      INTEGER IPROC, NHG, J, LINDEX
      CHARACTER*80 PBLJAGF
      INTEGER SDU(NPROC)  ! subdomain unit numbers
      LOGICAL Success     ! .true. if all files open without error
      REAL(SZ)                U,V
C
C--Open Global Wave Stress File ( UNIT 23 )
C
C     Open full domain and subdomain fort.23 files
      CALL OpenPrepFiles(23, 'wave stress                   ',
     &     1, nproc, SDU, Success)
      IF (.not.Success) THEN
         WRITE(*,*) 'WARNING: Unit 23 files not preprocessed.'
         RETURN ! note early return
      ENDIF
C--------------------------------------------------------------------------
C--MAIN LOOP
C   (1)  Read a record from Global Wave Stress File
C   (2)  Use Decomp arrays to Localize record to a subdomain
C   (3)  Write Local Wave Stress record in standard PBL format
C--------------------------------------------------------------------------
C
C--Read a wave field record from the global input file
C--and write out to respective local fort.23 file.
C
 170  READ(23,'(A80)',END=9999) PBLJAGF
      IF(PBLJAGF(2:2).EQ.'#') THEN
         DO IPROC = 1,NPROC
            WRITE(SDU(IPROC),1101)
            WRITE(SDU(IPROC),1100) 1,0.0,0.0 !victor didn't like this line 27/11/03
         ENDDO
      ELSE
! vjp 27/11/03
! rewrote this section to handle ghost-nodes
! and changed if test from "and" to "or"
         READ(PBLJAGF,'(I8,2E13.5)',END=9999) NHG,U,V
         IF ((U.NE.0.).OR.(V.NE.0.)) THEN
            DO J=1, ITOTPROC(NHG)
               IPROC  = IMAP_NOD_GL2(2*(J-1)+1,NHG)
               LINDEX = IMAP_NOD_GL2(2*(J-1)+2,NHG)
               WRITE(SDU(IPROC),1100) LINDEX,U,V
            ENDDO
         ENDIF
      ENDIF

      GOTO 170

 9999 CLOSE(23)
      DO IPROC=1,NPROC
         CLOSE(SDU(IPROC))
      ENDDO

 1100 FORMAT(I8,2E13.5)
 1101 FORMAT (' #')

 99   RETURN
      END


C---------------------------------------------------------------------------
C                S U B R O U T I N E   P R E P 141
C---------------------------------------------------------------------------
C
C                     (  Serial Version  4/13/12  )                         C
C  This routine reads a global external bathymetry file when NDDT=+-1,+-2.  C
C  In each case it wites a local bathymetry file of the same                C
C  format for each subdomain using the domain decomposition of the ADCIRC   C
C  grid created by the routine DECOMP.                                      C
C                                                                           C
C  The Decomposition Variables are defined in the include file adcprep.inc  C
C  This version is compatible with ADCIRC version 50.66                     C
C                                                                           C
C  TCM -v 50.66.03 Addition for time varying Bathymetry                     C
C          This routine adopted/modified from the prep22 subroutine.        C
C                                                                           C
C---------------------------------------------------------------------------
      SUBROUTINE PREP141()
C---------------------------------------------------------------------------
      USE PRE_GLOBAL
      use memory_usage
      IMPLICIT NONE
      integer :: nbytes = 0
      LOGICAL FOUND,DONE
      INTEGER I,J,IPROC,IPROC2,ILNODE,INDX,NHG,LINDEX
      CHARACTER*80 PBLJAGF
      INTEGER SDU(NPROC)  ! subdomain unit numbers
      LOGICAL Success     ! .true. if all files open without error
      INTEGER,ALLOCATABLE  :: NG(:)
      REAL(SZ),ALLOCATABLE :: DPG(:)  !global array
      REAL(SZ),ALLOCATABLE :: DPL(:)  !local array
      REAL(SZ) DPTMP
C     jgf48.47 Do the decomposition for a max of 256 subdomains at a
C     time ... some platforms/compilers limit the number of files that
C     can be open at any one time.
      INTEGER, PARAMETER :: maxOpenFiles = 256
      INTEGER startProc
      INTEGER endProc
      INTEGER deltaProc
C
C     Allocate local work arrays
C
      ALLOCATE ( NG(MNP) )
      nbytes = 4*mnp
      call memory_alloc(nbytes)
      ALLOCATE ( DPG(MNP) )   !global
      nbytes = 8*mnp
      call memory_alloc(nbytes)
      ALLOCATE ( DPL(MNP) )   !local
      nbytes = 8*mnp
      call memory_alloc(nbytes)
C
C     Perform decomposition over a range of subdomains.
      startProc = 1
      DO WHILE ( startProc .lt. nproc )
         deltaProc = nproc - startProc
         IF ( deltaProc .gt. maxOpenFiles ) deltaProc = maxOpenFiles
         endProc = startProc + deltaProc

C        Open full domain and all subdomain fort.141 files
         CALL OpenPrepFiles(141, 'bathymetry information        ',
     &      startProc, endProc, sdu, success)

         IF (.not.success) THEN
           WRITE(*,*) 'WARNING: Unit 141 files not preprocessed.'
           RETURN ! note early return
         ENDIF
C
C--Branch to Appropriate Code
C
         SELECT CASE(ABS(NDDT))
C        -------------
         CASE(1)
C        -------------
C
C     MAIN LOOP FOR NWS = +-1
C     (1)  Read a record from Global Bathymetry File
C     (2)  Use Decomp arrarys to Localize record to a subdomain
C     (3)  Write Local Bathymetry record in same format
           DO                     ! loop forever (or until file ends)
              READ(141,*,END=9999)
     &           (NG(I),DPG(I),I=1,NNODG)
              DO IPROC = STARTPROC, ENDPROC
                 DO I=1, NNODP(IPROC)
                    INDX = IMAP_NOD_LG(I,IPROC)
                    DPL(I) = DPG(INDX)
                 ENDDO
                 DO I=1, NNODP(IPROC)
                    WRITE(SDU(IPROC),*)  I,DPL(I)
                 ENDDO
              ENDDO
           ENDDO
C
C        -------
         CASE(2)
C        -------
C        MAIN LOOP FOR NWS = +- 2  ( PBL Format )
C        (1)  Read a record from Global Bathymetry File
C        (2)  Use Decomp arrarys to Localize record to a subdomain
C        (3)  Write out in PBL Format on subdomain
C
C--Read a bathymetry field record from the global input file
C--- during the decomp phase, after each time record indicator is written (#)
C--- we write a single entry (1,-99999.d0) to ensure that there will be no
C--- empty records.  When this file is read by ADCIRC using nddt2get, the
C--- extra entry (1,-99999.d0) will be ignored, and if node 1 actually is
C--- changed then it will be read regardless if it appears twice.
C
         DO
            READ(141,'(A80)',END=9999) PBLJAGF
            IF(PBLJAGF(2:2).EQ.'#') THEN
               DO IPROC =  STARTPROC,ENDPROC
                  WRITE(SDU(IPROC),1101)
!  write a default value to ensure that no empty records
!  are produced during the decomp phase (default values will be ignored by ADCIRC)
                  WRITE(SDU(IPROC),1100) 1,-99999.d0
               ENDDO
            ELSE

               READ(PBLJAGF,*,END=9999) NHG,DPTMP
               DO J=1, ITOTPROC(NHG)
                  IPROC  = IMAP_NOD_GL2(2*(J-1)+1,NHG)
                  LINDEX = IMAP_NOD_GL2(2*(J-1)+2,NHG)
                  WRITE(SDU(IPROC),1100) LINDEX,DPTMP
               ENDDO
            ENDIF
         END DO
C

C        ------------
         CASE DEFAULT
C        ------------
         WRITE(*,*) "NDDT = ",NDDT," has incorrect value in PREP141"
         RETURN

         END SELECT

C
C--Close Global file and all the Local Files
C
 9999    CLOSE (141)
         DO IPROC=STARTPROC, ENDPROC
            CLOSE (SDU(IPROC))
         ENDDO
         startProc=endProc+1
      ENDDO


      DEALLOCATE ( NG,  DPG )
      DEALLOCATE ( DPL )
      nbytes = 20*mnp
      call memory_dealloc(nbytes)
      call memory_status()
      RETURN
!  60  FORMAT(A60)
! 170  FORMAT(A170)
 !270  FORMAT(A270)
 !1010 FORMAT(' File ',A60,/,' WAS NOT FOUND!  Try again',/)
 !1011 FORMAT(' File ',A60,/,' WAS FOUND!  Opening & Processing file',/)
 1100 FORMAT(I8,E13.5)
 1101 FORMAT(' #')
c----------------------------------------------------------------------------
      END SUBROUTINE PREP141
c----------------------------------------------------------------------------








C   kmd48.33bc add in prep subroutines for 3D boundary condition files
      SUBROUTINE PREP35()
      USE PRE_GLOBAL
      use memory_usage
C
C---------------------------------------------------------------------------C
C                                                                           C
C  This routine writes a Local "Residual Boundary Condtions Baroclinic"     C
C  (fort.35) file for each subdomain using the domain decomposition of     C
C  the ADCIRC grid created by the routine DECOMP.                           C
C                                                                           C
C                   Added by Kendra Dresback (Aug. 18, 2007)                C
C---------------------------------------------------------------------------C
C
      IMPLICIT NONE
      INTEGER :: nbytes = 0
      INTEGER I,J,IPROC
      INTEGER SDU(NPROC) ! subdomain unit numbers
      LOGICAL Success    ! .true. if files opened without errors
      CHARACTER*40  ETIMINC,RESBCBINP,GRIDINC
      CHARACTER*40,ALLOCATABLE :: RESBCBIN(:)

C
C--Enter, Locate, Open, and Read the ADCIRC UNIT 35
C  Global Level of No Motion Boundary Conditions file for baroclinic
C
C     Open full domain and subdomain fort.35 files
C      Print *, "Made it to prepping the files"
      CALL OpenPrepFiles(35, 'level of no motion boundary  ',
     &     1, nproc, SDU, Success)
C      Print *, "Made it out of prepping the files"
      IF (.not.Success) THEN
         WRITE(*,*) 'WARNING: Unit 35 files not preprocessed.'
         RETURN ! note early return
      ENDIF
C
C--Allocate local arrays
C
      ALLOCATE ( RESBCBIN(MNETA) )
      nbytes = 8*mneta
      call memory_alloc(nbytes)
C
C--While ( NOT EOF ) Read NETA BCs from Global File
C
C      PRINT *, "Made it to the reading in of the 35 file"
      DO  ! loop until end of file
      READ(35,40,END=9999) ETIMINC
      DO IPROC = 1,NPROC
         WRITE(SDU(IPROC),40)  ETIMINC
      ENDDO
      DO I=1, NETA
         READ(35,40,END=9999)  RESBCBIN(I)
      ENDDO
C
      DO IPROC= 1,NPROC
         DO I=1, NETAP(IPROC)
            RESBCBINP = RESBCBIN(OBNODE_LG(I,IPROC))
            WRITE(SDU(IPROC),40) RESBCBINP
         ENDDO
      ENDDO
      END DO
C
C
C--Close Global file and all the Local Files
C
 9999 CLOSE (35)
      DO IPROC=1, NPROC
         CLOSE (SDU(IPROC))
      ENDDO

      deallocate(resbcbin)
      nbytes = 8*mneta
      call memory_dealloc(nbytes)
      call memory_status()
C
  40  FORMAT(A40)
C
      RETURN
      END

      SUBROUTINE PREP36()
      USE PRE_GLOBAL
      use memory_usage
C
C---------------------------------------------------------------------------C
C                                                                           C
C  This routine writes a Local "Salinity Boundary Conditions Values"        C
C  (fort.36) file for each subdomain using the domain decomposition of      C
C  the ADCIRC grid created by the routine DECOMP.                           C
C                                                                           C
C                Added by Kendra Dresback (January 15, 2008)                C
C---------------------------------------------------------------------------C
C
      IMPLICIT NONE
      integer :: nbytes = 0
      INTEGER I,J,IPROC
      INTEGER SDU(NPROC) ! subdomain unit numbers
      LOGICAL Success    ! .true. if files opened without errors
      CHARACTER*40  ETIMINC,GRIDINC
      INTEGER :: NODP, M
      INTEGER,ALLOCATABLE :: NOD(:)
      REAL(SZ),ALLOCATABLE :: SalBC(:,:)
      REAL(SZ),ALLOCATABLE :: RESBCBINP(:)


C
C--Enter, Locate, Open, and Read the ADCIRC UNIT 36
C  Global Salinity Boundary Conditions file for baroclinic
C
C     Open full domain and subdomain fort.36 files
      CALL OpenPrepFiles(36, 'salinity boundary             ',
     &     1, nproc, SDU, Success)
      IF (.not.Success) THEN
         WRITE(*,*) 'WARNING: Unit 36 files not preprocessed.'
         RETURN ! note early return
      ENDIF
C
C--Allocate local arrays
C
      ALLOCATE ( NOD(MNETA) )
      ALLOCATE ( RESBCBINP(NFEN) )
      ALLOCATE ( SalBC(MNETA,NFEN) )
      nbytes = 8*mneta
      call memory_alloc(nbytes)
C
C--While ( NOT EOF ) Read NETA BCs from Global File
C
      DO  ! loop until end of file
      READ(36,40,END=9999) ETIMINC
      DO IPROC = 1,NPROC
         WRITE(SDU(IPROC),40)  ETIMINC
      ENDDO

      DO I=1, NETA
         READ(36,*,END=9999)  NOD(I), (SalBC(I,M),M=1,NFEN)
      ENDDO
C
      DO IPROC= 1,NPROC
         DO I=1, NETAP(IPROC)
            NODP = NOD(OBNODE_LG(I,IPROC))
            DO M=1,NFEN
              RESBCBINP(M) = SalBC(OBNODE_LG(I,IPROC),M)
            END DO
              WRITE(SDU(IPROC),80) NODP, (RESBCBINP(M),M=1,NFEN)
         ENDDO
      ENDDO
      END DO
C
C
C--Close Global file and all the Local Files
C
 9999 CLOSE (36)
      DO IPROC=1, NPROC
         CLOSE (SDU(IPROC))
      ENDDO
C
      deallocate(salbc)
      nbytes = 8*mneta
      call memory_dealloc(nbytes)
      call memory_status()
  40  FORMAT(A40)
  80  FORMAT(1X,I6,1X,32000(F11.7,2X))
C
      RETURN
      END

      SUBROUTINE PREP37()
      USE PRE_GLOBAL
      use memory_usage
C
C---------------------------------------------------------------------------C
C                                                                           C
C  This routine writes a Local "Temperature Boundary Conditions Values"     C
C  (fort.37) file for each subdomain using the domain decomposition of      C
C  the ADCIRC grid created by the routine DECOMP.                           C
C                                                                           C
C                Added by Kendra Dresback (January 15, 2008)                C
C---------------------------------------------------------------------------C
C
      IMPLICIT NONE
      integer :: nbytes = 0
      INTEGER I,J,IPROC
      INTEGER SDU(NPROC) ! subdomain unit numbers
      LOGICAL Success    ! .true. if files opened without errors
      CHARACTER*40  ETIMINC
      INTEGER :: NODP, M
      INTEGER,ALLOCATABLE :: NOD(:)
      REAL(SZ),ALLOCATABLE :: TempBC(:,:)
      REAL(SZ),ALLOCATABLE :: RESBCBINP(:)

C
C--Enter, Locate, Open, and Read the ADCIRC UNIT 37
C  Global Temperature Boundary Conditions file for baroclinic
C
C     Open full domain and subdomain fort.37 files
      CALL OpenPrepFiles(37, 'temperature boundary          ',
     &     1, nproc, SDU, Success)
      IF (.not.Success) THEN
         WRITE(*,*) 'WARNING: Unit 37 files not preprocessed.'
         RETURN ! note early return
      ENDIF
C
C--Allocate local arrays
C
      ALLOCATE ( NOD(MNETA) )
      ALLOCATE ( RESBCBINP(NFEN) )
      ALLOCATE ( TempBC(MNETA,NFEN) )
      nbytes = 8*mneta
      call memory_alloc(nbytes)
C
C
C--While ( NOT EOF ) Read NETA BCs from Global File
C
      DO  ! loop around until the end of the file
      READ(37,40,END=9999) ETIMINC
      DO IPROC = 1,NPROC
         WRITE(SDU(IPROC),40)  ETIMINC
      ENDDO

      DO I=1, NETA
         READ(37,*,END=9999)  NOD(I), (TempBC(I,M),M=1,NFEN)
      ENDDO
C
      DO IPROC= 1,NPROC
         DO I=1, NETAP(IPROC)
            NODP = NOD(OBNODE_LG(I,IPROC))
            DO M=1,NFEN
              RESBCBINP(M) = TempBC(OBNODE_LG(I,IPROC),M)
            END DO
              WRITE(SDU(IPROC),80) NODP, (RESBCBINP(M),M=1,NFEN)
         ENDDO
      ENDDO
      END DO
C
C
C--Close Global file and all the Local Files
C
 9999 CLOSE (37)
      DO IPROC=1, NPROC
         CLOSE (SDU(IPROC))
      ENDDO
C
      deallocate(TempBC)
      nbytes = 8*mneta
      call memory_dealloc(nbytes)
      call memory_status()
  40  FORMAT(A40)
  80  FORMAT(1X,I6,1X,32000(F11.7,2X))
C
      RETURN
      END

      SUBROUTINE PREP38()
      USE PRE_GLOBAL
      use memory_usage
C
C---------------------------------------------------------------------------C
C                                                                           C
C  This routine writes a Local "Temperature Boundary Conditions Values      C
C  for the surface" (fort.38) file for each subdomain using the domain     C
C  decomposition of the ADCIRC grid created by the routine DECOMP.          C
C                                                                           C
C                Added by Kendra Dresback (October 15, 2008)                C
C---------------------------------------------------------------------------C
C
      IMPLICIT NONE
      integer :: nbytes = 0
      INTEGER I,J,IPROC
      INTEGER SDU(NPROC) ! subdomain unit numbers
      LOGICAL Success    ! .true. if files opened without errors
      CHARACTER*40  ETIMINC,GRIDINC
      INTEGER :: NODP, M, NFLUX
      INTEGER,ALLOCATABLE :: NOD(:)
      REAL(SZ),ALLOCATABLE :: TopTempBC(:,:)
      REAL(SZ),ALLOCATABLE :: RESBCBINP(:,:)


C
C--Enter, Locate, Open, and Read the ADCIRC UNIT 38
C  Global Salinity Boundary Conditions file for baroclinic
C
C     Open full domain and subdomain fort.38 files
      CALL OpenPrepFiles(38, 'top temperature boundary      ',
     &     1, nproc, SDU, Success)
      IF (.not.Success) THEN
         WRITE(*,*) 'WARNING: Unit 38 files not preprocessed.'
         RETURN ! note early return
      ENDIF
C
C  Determine how many values are in the top temperature boundary
C  condition

      IF (BCFLAG_TEMP.EQ.1) THEN
         NFLUX = 1
      ELSE IF (BCFLAG_TEMP.EQ.2) THEN
         NFLUX = 6
      ELSE IF (BCFLAG_TEMP.EQ.3) THEN
         NFLUX = 4
      END IF

      MNP=nnodg
C
C--Allocate local arrays
C
      ALLOCATE ( NOD(MNP) )
      ALLOCATE ( RESBCBINP(MNP,NFLUX) )
      ALLOCATE ( TopTempBC(MNP,NFLUX) )
      nbytes = 24*mnp
      call memory_alloc(nbytes)
C
C--While ( NOT EOF ) Read NETA BCs from Global File
C
      DO  ! loop until end of file

       READ(38,*,END=9999)  (NOD(I),(TopTempBC(I,M),M=1,NFLUX),I=1,NNODG)
C
      DO IPROC= 1,NPROC
         DO I=1, NNODP(IPROC)
            NODP = IMAP_NOD_LG(I,IPROC)
            DO M=1,NFLUX
              RESBCBINP(I,M) = TopTempBC(NODP,M)
            END DO
              WRITE(SDU(IPROC),80) I, (RESBCBINP(I,M),M=1,NFLUX)
         ENDDO
      ENDDO
      END DO
C
C
C--Close Global file and all the Local Files
C
 9999 CLOSE (38)
      DO IPROC=1, NPROC
         CLOSE (SDU(IPROC))
      ENDDO
C
      deallocate(toptempbc)
      nbytes = 24*mnp
      call memory_dealloc(nbytes)
      call memory_status()
  40  FORMAT(A40)
  80  FORMAT(1X,I8,1X,32(F12.6,2X))
C
      RETURN
      END

      SUBROUTINE PREP39()
C---------------------------------------------------------------------------
C                                                                           C
C  This routine writes a Local river boundary file for the baroclnic        C
C  simulation (fort.39) for each subdomain using the domain                 C
C  decomposition of the ADCIRC grid created by the routine DECOMP.          C
C                                                                           C
C                Added by Kendra Dresback (January 14, 2010)                C
C---------------------------------------------------------------------------C
C
      USE PRE_GLOBAL
      use memory_usage
      IMPLICIT NONE
      integer :: nbytes = 0
      INTEGER IPROC
      INTEGER INDEX14, I
      REAL(SZ) :: FLUX_INC
      REAL(SZ),ALLOCATABLE ::  FLUX_VAL(:,:)
      INTEGER SDU(NPROC)  ! subdomain unit numbers
      LOGICAL Success     ! .true. if all files open without error
      INTEGER INDX ! full domain node number for a flow boundary node
      INTEGER J,M     ! counter for subdomains that corrsp. to a single f.d. node
      INTEGER IPROC2! PE of a subdomain that matches a single full domain node
      INTEGER, ALLOCATABLE :: NOD(:)
      REAL(SZ),ALLOCATABLE ::  RESBCBINP(:)

C     Open full domain and subdomain fort.20 files
      CALL OpenPrepFiles(39, 'aperiodic river temp and salinity   ',
     &     1, nproc, SDU, Success)
      IF (.not.Success) THEN
         WRITE(*,*) 'WARNING: Unit 39 files not preprocessed.'
         RETURN ! note early return
      ENDIF

      ALLOCATE ( NOD(MNVEL) )
      ALLOCATE ( RESBCBINP(NFEN) )
      ALLOCATE ( FLUX_VAL(MNVEL,NFEN) )
c
c     Write Increment into all flux files
c
      READ(39,*) FLUX_INC
      DO IPROC=1,NPROC
         WRITE(SDU(IPROC),*) FLUX_INC
      ENDDO
C
C     jgf45.12 Write each full domain nodal flux value into each of the
C     subdomains that that full domain node maps to. The full domain
C     node may map to more than one subdomain node if it falls on a
C     boundary between subdomains (ghost nodes).
C

      DO  ! continue to loop over file until you reach the end of the file

      DO I=1, EXIST_BC_TS      ! loop through full domain flow nodes
         INDX=BCTS14_ARY(I)
         READ(39,*,END=40)  (FLUX_VAL(INDX,M),M=1,NFEN)
      END DO

      DO I=1, EXIST_BC_TS
         INDX = BCTS14_ARY(I)      ! get full domain flow node number
         DO J=1, ITOTPROC(INDX)    ! loop over subdomains for 1 f.d. node
            IPROC2 = IMAP_NOD_GL2(2*(J-1)+1,INDX) ! find next subdomain
            DO IPROC=1, NPROC
               IF (IPROC.EQ.IPROC2) THEN ! full domain node maps to this s.d.
                 DO M=1,NFEN
                    RESBCBINP(M) = FLUX_VAL(INDX,M)
                 END DO
                 WRITE(SDU(IPROC),80) (RESBCBINP(M),M=1,NFEN)
               ENDIF
            END DO
         END DO
      END DO
      END DO

 40   CLOSE (39)
      DO IPROC=1, NPROC
         CLOSE (SDU(IPROC))
      ENDDO

      IF (allocated(BCTS14_ARY)) then
        DEALLOCATE (BCTS14_ARY)
        nbytes = 4*exist_bc_ts
        call memory_dealloc(nbytes)
      ENDIF
      call memory_status()
      return
 80   FORMAT(1X,32000(F11.7,2X))
c----------------------------------------------------------------------------
      END SUBROUTINE PREP39
c----------------------------------------------------------------------------


C   kmd48.33bc add information for initial condition file
      SUBROUTINE HOTINITCOND()
      USE PRE_GLOBAL
      use presizes; use memory_usage
C
C---------------------------------------------------------------------------C
C                     written 10/11/01 by RL                                C
C             started mods for harmonic analysis and 3D RL 5/22/03          C
C         jgf Updated for v45.06 09/07/2005 not incl. harmonic or 3D        C
C         kmd Updated for v48.33 07/07/2008 to bring in initial conditions  C
C                                                                           C
C  This routine reads the global initial condition file (fort.17)           C
C  and writes local hot start files of the same format.                     C
C                                                                           C
C---------------------------------------------------------------------------C
C
      IMPLICIT NONE
      integer :: nbytes = 0
      LOGICAL FOUND
      INTEGER I,J,IPROC,IINDX,IHOTSTP, not_active
      INTEGER IMHSF,ITHSF
      CHARACTER FNAME*60,LOCFN*14
      CHARACTER*16 FNAME1
      CHARACTER*8 FNAM8(2)
      EQUIVALENCE (FNAM8(1),FNAME1)

      INTEGER,ALLOCATABLE  :: LOC2(:),NOFF(:), domA(:)
      REAL(SZ),ALLOCATABLE :: ETA1(:),ETA2(:),EtaDisc(:),
     &    UU2(:),VV2(:),CH1(:)
      REAL(8) TIMEHSF
      integer :: InputFileFmtVn, NP_G_IN, NE_G_IN, NP_A_IN, NE_A_IN
      CHARACTER*60 FileFmtVn
      INTEGER SDU(NPROC)  ! subdomain unit numbers
      LOGICAL Success     ! .true. if all files open without error

#if 0
      ! vjp 2006/9/30 not supporting harmonic analysis or C3D yet
      INTEGER INZ,INF,IMM,INP,INSTAE,INSTAV,IISTAE,IISTAV,IIGLOE,IIGLOV,
     &                                       IICALL,INFREQ,ITUD,NTSTEPS
      INTEGER ITHAS,ITHAF,ITMV,IHABEG,ICHA
      CHARACTER*10,ALLOCATABLE     ::  INAMEFR(:)
      REAL(8)  TIMEUD
      REAL(SZ),ALLOCATABLE ::  HA(:,:)
      REAL(SZ),ALLOCATABLE ::  ELAV(:),ELVA(:),XVELAV(:),XVELVA(:),
     &                                         YVELAV(:),YVELVA(:)
      REAL(SZ),ALLOCATABLE ::  IFREQ(:),IFF(:),IFACE(:)
      REAL(SZ),ALLOCATABLE ::  GLOELV(:,:)
      REAL(SZ),ALLOCATABLE ::  GLOULV(:,:),GLOVLV(:,:)
      REAL(SZ),ALLOCATABLE ::  STAELV(:,:)
      REAL(SZ),ALLOCATABLE ::  STAULV(:,:),STAVLV(:,:)
#endif
C
C--   Open the Initial Condition Start File based on the value of IHOT from
C--   the fort.15 file
C
C     Open full domain and subdomain fort.17 files
      Print *, "Made it to prepping the files"
      CALL OpenPrepFiles(17, 'initial condition file  ',
     &     1, nproc, SDU, Success)
      Print *, "Made it out of prepping the files"
      IF (.not.Success) THEN
         WRITE(*,*) 'WARNING: Unit 17 files not preprocessed.'
         RETURN ! note early return
      ENDIF
      IHOT=17

C--   Read in info from global initial condition file

      READ(IHOT,*) FileFmtVn

      READ(IHOT,*) IMHSF
      READ(IHOT,*) TIMEHSF
      READ(IHOT,*) ITHSF
      READ(IHOT,*) NP_G_IN
      READ(IHOT,*) NE_G_IN
      READ(IHOT,*) NP_A_IN
      READ(IHOT,*) NE_A_IN
      if (nnodg == np_g_in) then
        MNP = nnodg
      else
        print *, "number global nodes does not match hotstart file"
        write(*,'(A,I8)') "expected value   = ", nnodg
        write(*,'(A,I8)') "hotstart value = ", np_g_in
        stop
      endif
      if (nelg ==  ne_g_in) then
        MNE = nelg
      else
        print *, "number global elements does not match hotstart file"
        write(*,'(A,I8)') "expected value   = ", nelg
        write(*,'(A,I8)') "hotstart value = ", ne_g_in
        stop
      endif
C
C Allocate local work arrays
C
      MNP = nnodg
      MNE = nelg
      nbytes = 4*nproc
      call memory_alloc(nbytes)
      ALLOCATE ( ETA1(MNP),ETA2(MNP),EtaDisc(MNP),UU2(MNP),
     &           VV2(MNP),NODECODE(MNP),CH1(MNP) )
      nbytes = 7*mnp
      call memory_alloc(nbytes)
      ALLOCATE ( NOFF(MNE) )
      nbytes = 4*mne
      call memory_alloc(nbytes)

#if 0
      ! vjp 2006/9/30 not supporting harmonic analysis or C3D yet
      ALLOCATE ( HA(2*MNHARF,2*MNHARF) )
      nbytes = 32*mnharf
      call memory_alloc(nbytes)
      ALLOCATE ( GLOELV(2*MNHARF,MNP) )
      nbytes = 16*mnharf*mnp
      call memory_alloc(nbytes)
      ALLOCATE ( GLOULV(2*MNHARF,MNP),GLOVLV(2*MNHARF,MNP) )
      nbytes = 32*mnharf*mnp
      call memory_alloc(nbytes)
      ALLOCATE ( STAELV(2*MNHARF,MNSTAE) )
      nbytes = 16*mnharf*mnstae
      call memory_alloc(nbytes)
      ALLOCATE ( STAULV(2*MNHARF,MNSTAV),STAVLV(2*MNHARF,MNSTAV) )
      nbytes = 16*mnharf*mnstav
      call memory_alloc(nbytes)
      ALLOCATE ( ELAV(MNP),ELVA(MNP) )
      nbytes = 16*mnp
      call memory_alloc(nbytes)
      ALLOCATE ( XVELAV(MNP),XVELVA(MNP),YVELAV(MNP),YVELVA(MNP) )
      nbytes = 32*mnp
      call memory_alloc(nbytes)
      ALLOCATE ( IFREQ(MNHARF),IFF(MNHARF),IFACE(MNHARF) )
      nbytes = 12*mnharf
      call memory_alloc(nbytes)
      ALLOCATE ( INAMEFR(MNHARF) )
      nbytes = 4*mnharf
      call memory_alloc(nbytes)
#endif
!  Continue reading global initial condition file
      print *, "continuing to read global initial condition file"
      write(*,*) "enter number of layers: "
      read(*,*) NFEN

      DO I=1,MNP
         READ(IHOT,*) ETA1(I)
      END DO
      DO I=1,MNP
         READ(IHOT,*) ETA2(I)
      END DO
      DO I=1,MNP
         READ(IHOT,*) UU2(I)
      END DO
      DO I=1,MNP
         READ(IHOT,*) VV2(I)
      END DO
      IF(IM.EQ.10) THEN
        DO I=1,MNP
           IHOTSTP=IHOTSTP+1
           READ(IHOT,REC=IHOTSTP) CH1(I)
        END DO
      ENDIF
      DO I=1,MNP
         READ(IHOT,*) NODECODE(I)
      END DO
      DO I=1,MNE
         READ(IHOT,*) NOFF(I)
      END DO

      PRINT *, "Made it through the 2D values"

       ! vjp 2006/9/30 not supporting harmonic analysis or C3D yet
C     jgf46.02 Read in 3D hotstart data if appropriate
      IF (IMHSF.GT.10) THEN
         PRINT *, "set to go into 3D read"
         PRINT *, "NFEN = ", NFEN
         CALL ReadInitCond3D(IHOT)
      ENDIF
#if 0
C
C.....DETERMINE HARMONIC ANALYSIS PARAMETERS

      IHARIND=NHARFR*(NHASE+NHASV+NHAGE+NHAGV)
      IF(IHARIND.GT.0) IHARIND=1

C.....IF HARMONIC ANALYSIS IS INCLUDED IN THE RUN, PROCESS HOT START INFORMATION FOR
C.....IN PROGRESS HARMONIC ANALYSIS

      IF(IHARIND.EQ.1) THEN
         ITHAS=INT((THAS-STATIM)*(86400.D0/DT) + 0.5d0)
         ITHAF=INT((THAF-STATIM)*(86400.D0/DT) + 0.5d0)
         ITMV = ITHAF - (ITHAF-ITHAS)*FMV
         IHABEG=ITHAS+NHAINC

C.......IF HARMONIC ANALYSIS HAS ALREADY BEGUN, READ IN HOT START
C........HARMONIC ANALYSIS, MEAN AND SQUARE INFO

         IF(ITHSF.GT.ITHAS) THEN
            IHOTSTP=IHOTSTP+1
            READ(IHOT,REC=IHOTSTP) ICHA
         ENDIF

         IF(ITHSF.GE.IHABEG) THEN
            READ(IHOT,REC=IHOTSTP+1) INZ
            READ(IHOT,REC=IHOTSTP+2) INF
            READ(IHOT,REC=IHOTSTP+3) IMM
            READ(IHOT,REC=IHOTSTP+4) INP
            READ(IHOT,REC=IHOTSTP+5) INSTAE
            READ(IHOT,REC=IHOTSTP+6) INSTAV
            READ(IHOT,REC=IHOTSTP+7) IISTAE
            READ(IHOT,REC=IHOTSTP+8) IISTAV
            READ(IHOT,REC=IHOTSTP+9) IIGLOE
            READ(IHOT,REC=IHOTSTP+10) IIGLOV
            READ(IHOT,REC=IHOTSTP+11) IICALL
            READ(IHOT,REC=IHOTSTP+12) INFREQ
            IHOTSTP = IHOTSTP+12

            DO I=1,INFREQ+INF
               READ(IHOT,REC=IHOTSTP+1) FNAM8(1)
               READ(IHOT,REC=IHOTSTP+2) FNAM8(2)
               IHOTSTP = IHOTSTP + 2
               INAMEFR(I) = FNAME1
               READ(IHOT,REC=IHOTSTP+1) IFREQ(I)
               READ(IHOT,REC=IHOTSTP+2) IFF(I)
               READ(IHOT,REC=IHOTSTP+3) IFACE(I)
               IHOTSTP = IHOTSTP + 3
            ENDDO

            READ(IHOT,REC=IHOTSTP+1) TIMEUD
            READ(IHOT,REC=IHOTSTP+2) ITUD
            IHOTSTP = IHOTSTP + 2

            DO I=1,IMM
               DO J=1,IMM
                  IHOTSTP = IHOTSTP + 1
                  READ(IHOT,REC=IHOTSTP) HA(I,J)
               ENDDO
            ENDDO

            IF(NHASE.EQ.1) THEN
               DO J=1,INSTAE
                  DO I=1,IMM
                     IHOTSTP=IHOTSTP+1
                     READ(IHOT,REC=IHOTSTP) STAELV(I,J)
                  ENDDO
               ENDDO
            ENDIF

            IF(NHASV.EQ.1) THEN
               DO J=1,INSTAV
                  DO I=1,IMM
                     READ(IHOT,REC=IHOTSTP+1) STAULV(I,J)
                     READ(IHOT,REC=IHOTSTP+2) STAVLV(I,J)
                     IHOTSTP = IHOTSTP + 2
                  ENDDO
               ENDDO
            ENDIF

            IF(NHAGE.EQ.1) THEN
               DO J=1,INP
                  DO I=1,IMM
                     IHOTSTP=IHOTSTP+1
                     READ(IHOT,REC=IHOTSTP) GLOELV(I,J)
                  ENDDO
               ENDDO
            ENDIF

            IF(NHAGV.EQ.1) THEN
               DO J=1,INP
                  DO I=1,IMM
                     READ(IHOT,REC=IHOTSTP+1) GLOULV(I,J)
                     READ(IHOT,REC=IHOTSTP+2) GLOVLV(I,J)
                     IHOTSTP = IHOTSTP + 2
                  ENDDO
               ENDDO
            ENDIF

         ENDIF

         IF((FMV.GT.0.).AND.(INFREQ.GT.0).AND.(IM.EQ.0)) THEN !include means and variances
            IF(ITHSF.GT.ITMV) THEN
               IHOTSTP=IHOTSTP+1
               READ(IHOT,REC=IHOTSTP) NTSTEPS
               IF(NHAGE.EQ.1) THEN
                  DO I=1,INP
                     READ(IHOT,REC=IHOTSTP+1) ELAV(I)
                     READ(IHOT,REC=IHOTSTP+2) ELVA(I)
                     IHOTSTP=IHOTSTP+2
                  ENDDO
               ENDIF
               IF(NHAGV.EQ.1) THEN
                  DO I=1,INP
                     READ(IHOT,REC=IHOTSTP+1) XVELAV(I)
                     READ(IHOT,REC=IHOTSTP+2) YVELAV(I)
                     READ(IHOT,REC=IHOTSTP+3) XVELVA(I)
                     READ(IHOT,REC=IHOTSTP+4) YVELVA(I)
                     IHOTSTP=IHOTSTP+4
                  ENDDO
               ENDIF
            ENDIF
         ENDIF    ! charmv
      ENDIF     ! HARIND
#endif

C
C--Open All Local Hot Start files
C
      ALLOCATE ( LOC2(NPROC) )
      DO IPROC = 1,NPROC
         LOC2(IPROC) = 105 + (IPROC-1)
         LOCFN(1:14) = 'PE0000/'//FNAME(1:7)
         CALL IWRITE(LOCFN,3,6,IPROC-1)
         OPEN (LOC2(IPROC),FILE=LOCFN)
      ENDDO
C
C--Write out info to local hot start files
C
      DO IPROC = 1,NPROC
         WRITE(LOC2(IPROC),*) FileFmtVn
         WRITE(LOC2(IPROC),*) IMHSF
         WRITE(LOC2(IPROC),*) TIMEHSF
         WRITE(LOC2(IPROC),*) ITHSF
         WRITE(LOC2(IPROC),*) NNODP(IPROC)
         WRITE(LOC2(IPROC),*) NELP(IPROC)
         WRITE(LOC2(IPROC),*) NNODP(IPROC)
         WRITE(LOC2(IPROC),*) NELP(IPROC)

         DO I=1, NNODP(IPROC)
            IINDX = ABS(IMAP_NOD_LG(I,IPROC))
            WRITE(LOC2(IPROC),*) ETA1(IINDX)
         END DO

         DO I=1, NNODP(IPROC)
            IINDX = ABS(IMAP_NOD_LG(I,IPROC))
            WRITE(LOC2(IPROC),*) ETA2(IINDX)
         END DO

         DO I=1, NNODP(IPROC)
            IINDX = ABS(IMAP_NOD_LG(I,IPROC))
            WRITE(LOC2(IPROC),*) UU2(IINDX)
         END DO

         DO I=1, NNODP(IPROC)
            IINDX = ABS(IMAP_NOD_LG(I,IPROC))
            WRITE(LOC2(IPROC),*) VV2(IINDX)
         END DO

         IF(IM.EQ.10) THEN
            DO I=1, NNODP(IPROC)
               IINDX = ABS(IMAP_NOD_LG(I,IPROC))
               WRITE(LOC2(IPROC),*) CH1(IINDX)
            END DO
         ENDIF

         DO I=1, NNODP(IPROC)
            IINDX = ABS(IMAP_NOD_LG(I,IPROC))
            WRITE(LOC2(IPROC),*) NODECODE(IINDX)
         END DO

         DO I=1,NELP(IPROC)
            IINDX=ABS(IMAP_EL_LG(I,IPROC))
            WRITE(LOC2(IPROC),*) NOFF(IINDX)
         END DO

C
C     jgf46.02 Write out 3D hotstart data if appropriate
         IF (IMHSF.GT.10) THEN
            CALL WriteInitCond3D(LOC2(IPROC),IPROC)
         ENDIF
#if 0
C
C....IF APPROPRIATE, WRITE OUT HOT START INFORMATION FOR IN PROGRESS HARMONIC ANALYSIS

c       IF((IHARIND.EQ.1).AND.(ITHSF.GT.ITHAS)) THEN
c         WRITE(LOC2(IPROC),REC=IHOTSTP+1) ICHA
c         IHOTSTP = IHOTSTP + 1
c         CALL HAHOUT(NP,NSTAE,NSTAV,NHASE,NHASV,NHAGE,NHAGV,
c    &                LOC2(IPROC),IHOTSTP)
c
c         IF(NHASE.EQ.1) CALL HAHOUTES(NSTAE,LOC2(IPROC),IHOTSTP)
c         IF(NHASV.EQ.1) CALL HAHOUTVS(NSTAV,LOC2(IPROC),IHOTSTP)
c         IF(NHAGE.EQ.1) CALL HAHOUTEG(MNP,LOC2(IPROC),IHOTSTP)
c         IF(NHAGV.EQ.1) CALL HAHOUTVG(MNP,LOC2(IPROC),IHOTSTP)
c         ENDIF
c
c       if(CHARMV) then
c         IF((IHARIND.EQ.1).AND.(ITHSF.GT.ITMV)) THEN
c           IHOTSTP=IHOTSTP+1
c           WRITE(LOC2(IPROC),REC=IHOTSTP) NTSTEPS
c           IF(NHAGE.EQ.1) THEN
c             DO I=1, NNODP(IPROC)
c               IINDX = IMAP_NOD_LG(I,IPROC)
c               DO I=1,MNP
c                 WRITE(LOC2(IPROC),REC=IHOTSTP+1) ELAV(IINDX)
c                 WRITE(LOC2(IPROC),REC=IHOTSTP+2) ELVA(IINDX)
c                 IHOTSTP=IHOTSTP+2
c                 END DO
c             ENDIF
c           IF(NHAGV.EQ.1) THEN
c             DO I=1,NNODP(IPROC)
c               WRITE(LOC2(IPROC),REC=IHOTSTP+1) XVELAV(IINDX)
c               WRITE(LOC2(IPROC),REC=IHOTSTP+2) YVELAV(IINDX)
c               WRITE(LOC2(IPROC),REC=IHOTSTP+3) XVELVA(IINDX)
c               WRITE(LOC2(IPROC),REC=IHOTSTP+4) YVELVA(IINDX)
c               IHOTSTP=IHOTSTP+4
c               END DO
c             ENDIF
c           ENDIF
c         ENDIF
#endif

      ENDDO
C
C--Close Global file and all the Local Files
C
      CLOSE (IHOT)
      DO IPROC=1, NPROC
         CLOSE (LOC2(IPROC))
      ENDDO
C
      DEALLOCATE ( LOC2 )
      nbytes = 4*nproc
      call memory_dealloc(nbytes)
      DEALLOCATE ( ETA1, ETA2, EtaDisc, UU2, VV2, NODECODE, CH1 )
      nbytes = 7*mnp
      call memory_dealloc(nbytes)
      DEALLOCATE ( NOFF )
      nbytes = 6*mne
      call memory_dealloc(nbytes)
      call memory_status()
C
      RETURN
 1001 FORMAT('ERROR: The hot start file')
 1010 FORMAT(' File ',A60,/,' WAS NOT FOUND!  ADCPrep Terminated!!!',/)
 1011 FORMAT(' File ',A60,/,' WAS FOUND!  Opening & Processing file',/)
 1012 FORMAT('was a nonmatching version')
 1005 FORMAT('exists but cannot be opened.')
 9973 FORMAT(/,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
      END SUBROUTINE HOTINITCOND


      SUBROUTINE HOTLOCALIZE()
      USE VERSION
      USE PRE_GLOBAL
      use presizes; use memory_usage
C
C---------------------------------------------------------------------------C
C                     written 10/11/01 by RL                                C
C             started mods for harmonic analysis and 3D RL 5/22/03          C
C         jgf Updated for v45.06 09/07/2005 not incl. harmonic or 3D        C
C         kmd48.33bc updated with 3D information                            C
C                                                                           C
C  This routine reads the global hot start file (either fort.67 or fort.68) C
C  and writes local hot start files of the same format.                     C
C                                                                           C
C---------------------------------------------------------------------------C
C
      IMPLICIT NONE
      integer :: nbytes = 0
      LOGICAL FOUND
      INTEGER I,J,IPROC,INDX,IHOTSTP, not_active
      INTEGER IMHSF,ITHSF, NH, N
      CHARACTER FNAME*60,LOCFN*14
      CHARACTER*16 FNAME1
      CHARACTER*8 FNAM8(2)
      EQUIVALENCE (FNAM8(1),FNAME1)

      INTEGER,ALLOCATABLE  :: LOC2(:),NOFF(:), domA(:)
      REAL(SZ),ALLOCATABLE :: ETA1(:),ETA2(:),EtaDisc(:),
     &    UU2(:),VV2(:),CH1(:)
      REAL(8) TIMEHSF
      integer :: InputFileFmtVn, NP_G_IN, NE_G_IN, NP_A_IN, NE_A_IN

      INTEGER INZ,INF,IMM,INP,INSTAE,INSTAV,IISTAE,IISTAV,IIGLOE,IIGLOV,
     &                                       IICALL,INFREQ,ITUD,NTSTEPS
      INTEGER IHARIND,ITHAS,ITHAF,ITMV,IHABEG,ICHA
      CHARACTER*10,ALLOCATABLE     ::  INAMEFR(:)
      REAL(8)  TIMEUD
      REAL(SZ),ALLOCATABLE ::  HA(:,:)
      REAL(SZ),ALLOCATABLE ::  ELAV(:),ELVA(:),XVELAV(:),XVELVA(:),
     &                                         YVELAV(:),YVELVA(:)
      REAL(SZ),ALLOCATABLE ::  IFREQ(:),IFF(:),IFACE(:)
      REAL(SZ),ALLOCATABLE ::  GLOELV(:,:)
      REAL(SZ),ALLOCATABLE ::  GLOULV(:,:),GLOVLV(:,:)
      REAL(SZ),ALLOCATABLE ::  STAELV(:,:)
      REAL(SZ),ALLOCATABLE ::  STAULV(:,:),STAVLV(:,:)
      REAL(SZ) TIME

      REAL(SZ) DUMMY
      INTEGER IDUMMY
      INTEGER LUN
      INTEGER NHS
C
C--   Open Appropriate Hot Start File based on the value of IHOT from
C--   the fort.15 file
C
      write(*,*) "enter IHOT: "
      read(*,*) IHOT
      SELECT CASE (IHOT)
      CASE(67)
         FNAME='fort.67'
      CASE(68)
         FNAME='fort.68'
      CASE(367,368)
         write(*,*) "INFO: IHOT=",IHOT,
     &      " means parallel ADCIRC should read a NetCDF hotstart file."
         write(*,*)
     &      "INFO: NetCDF hotstart files do not require decomposition."
         RETURN
      CASE DEFAULT
         write(*,*) "ERROR: The IHOT value ",IHOT,
     &      " is not a valid option."
         write(*,*) "INFO: 67 and 68 are the only valid options."
         RETURN
      END SELECT
C
      INQUIRE(FILE=FNAME,EXIST=FOUND)
      IF (FOUND) THEN
         WRITE(*,1011) FNAME
         IF(IHOT.EQ.67.OR.IHOT.EQ.68)
     &      OPEN(IHOT,FILE=FNAME,ACCESS='DIRECT',RECL=8)
      ELSE
         WRITE(*,1010) FNAME
         STOP
      ENDIF

C--   Read in info from global hot start files

      IHOTSTP=1
      READ(IHOT,REC=IHOTSTP) InputFileFmtVn ; IHOTSTP = IHOTSTP + 1

      if (.not. CMP_VERSION_NUMBERS(InputFileFmtVn, FileFmtVersion)) then
        write(*, 1001)
        write(*, 1012)
        write(*, 9973)
        !stop
      endif

      READ(IHOT,REC=IHOTSTP) IMHSF        ; IHOTSTP = IHOTSTP + 1
      READ(IHOT,REC=IHOTSTP) TIMEHSF      ; IHOTSTP = IHOTSTP + 1
      READ(IHOT,REC=IHOTSTP) ITHSF        ; IHOTSTP = IHOTSTP + 1
      READ(IHOT,REC=IHOTSTP) NP_G_IN      ; IHOTSTP = IHOTSTP + 1
      READ(IHOT,REC=IHOTSTP) NE_G_IN      ; IHOTSTP = IHOTSTP + 1
      READ(IHOT,REC=IHOTSTP) NP_A_IN      ; IHOTSTP = IHOTSTP + 1
      READ(IHOT,REC=IHOTSTP) NE_A_IN      ; IHOTSTP = IHOTSTP + 1
      if (nnodg == np_g_in) then
        MNP = nnodg
        print *, "MNP = ", MNP
      else
        print *, "number global nodes does not match hotstart file"
        write(*,'(A,I8)') "expected value   = ", nnodg
        write(*,'(A,I8)') "hotstart value = ", np_g_in
        stop
      endif
      if (nelg ==  ne_g_in) then
        MNE = nelg
      else
        print *, "number global elements does not match hotstart file"
        write(*,'(A,I8)') "expected value   = ", nelg
        write(*,'(A,I8)') "hotstart value = ", ne_g_in
        stop
      endif
      PRINT *, "IMHSF ", IMHSF
C
C Allocate local work arrays
C
      nbytes = 4*nproc
      call memory_alloc(nbytes)
      ALLOCATE ( ETA1(MNP),ETA2(MNP),EtaDisc(MNP),UU2(MNP),
     &           VV2(MNP),NODECODE(MNP),CH1(MNP) )
      nbytes = 7*mnp*8
      call memory_alloc(nbytes)
      ALLOCATE ( NOFF(MNE) )
      nbytes = 4*mne
      call memory_alloc(nbytes)

!  Continue reading global hot start file
      print *, "continuing to read global hotstart file"

      DO I=1,MNP
         READ(IHOT,REC=IHOTSTP) ETA1(I) ; IHOTSTP = IHOTSTP + 1
      END DO
      DO I=1,MNP
         READ(IHOT,REC=IHOTSTP) ETA2(I) ; IHOTSTP = IHOTSTP + 1
      END DO
      DO I=1,MNP
         READ(IHOT,REC=IHOTSTP) EtaDisc(I) ; IHOTSTP = IHOTSTP + 1
      END DO
      DO I=1,MNP
         READ(IHOT,REC=IHOTSTP) UU2(I) ; IHOTSTP = IHOTSTP + 1
      END DO
      DO I=1,MNP
         READ(IHOT,REC=IHOTSTP) VV2(I) ; IHOTSTP = IHOTSTP + 1
      END DO
      IF(IMHSF.EQ.10) THEN
        DO I=1,MNP
           READ(IHOT,REC=IHOTSTP) CH1(I) ; IHOTSTP = IHOTSTP + 1
        END DO
      ENDIF
      DO I=1,MNP
         READ(IHOT,REC=IHOTSTP) NODECODE(I) ; IHOTSTP = IHOTSTP + 1
      END DO

      DO I=1,MNE
         READ(IHOT,REC=IHOTSTP) NOFF(I)  ; IHOTSTP = IHOTSTP + 1
      END DO

      READ(IHOT,REC=IHOTSTP) IESTP ; IHOTSTP = IHOTSTP + 1
      READ(IHOT,REC=IHOTSTP) NSCOUE ; IHOTSTP = IHOTSTP + 1

      READ(IHOT,REC=IHOTSTP) IVSTP ; IHOTSTP = IHOTSTP + 1
      READ(IHOT,REC=IHOTSTP) NSCOUV ; IHOTSTP = IHOTSTP + 1

      READ(IHOT,REC=IHOTSTP) ICSTP ; IHOTSTP = IHOTSTP + 1
      READ(IHOT,REC=IHOTSTP) NSCOUC ; IHOTSTP = IHOTSTP + 1

      READ(IHOT,REC=IHOTSTP) IPSTP ; IHOTSTP = IHOTSTP + 1
      READ(IHOT,REC=IHOTSTP) IWSTP ; IHOTSTP = IHOTSTP + 1
      READ(IHOT,REC=IHOTSTP) NSCOUM ; IHOTSTP = IHOTSTP + 1

      READ(IHOT,REC=IHOTSTP) IGEP ; IHOTSTP = IHOTSTP + 1
      READ(IHOT,REC=IHOTSTP) NSCOUGE ; IHOTSTP = IHOTSTP + 1

      READ(IHOT,REC=IHOTSTP) IGVP ; IHOTSTP = IHOTSTP + 1
      READ(IHOT,REC=IHOTSTP) NSCOUGV ; IHOTSTP = IHOTSTP + 1

      READ(IHOT,REC=IHOTSTP) IGCP ; IHOTSTP = IHOTSTP + 1
      READ(IHOT,REC=IHOTSTP) NSCOUGC ; IHOTSTP = IHOTSTP + 1

      READ(IHOT,REC=IHOTSTP) IGPP ; IHOTSTP = IHOTSTP + 1
      READ(IHOT,REC=IHOTSTP) IGWP ; IHOTSTP = IHOTSTP + 1
      READ(IHOT,REC=IHOTSTP) NSCOUGW ; IHOTSTP = IHOTSTP + 1

C kmd48.33 moved 3D hot start information to subroutine
C          and took out other lines
C jgf49.17 refined check of IMHSF so that it picks up only
C     IM values that indicate 3D (and so we can use six integer IM values).
      IF ((IMHSF.EQ.1).OR.(IMHSF.EQ.11).OR.
     &    (IMHSF.EQ.21).OR.(IMHSF.EQ.31)) THEN
         CALL ReadHotStart3D(IHOT,IHOTSTP)
      ENDIF

C     jgf48.03 harmonic analysis not supported yet
#if 0
C
C....DETERMINE HARMONIC ANALYSIS PARAMETERS

      IHARIND=NHARFR*(NHASE+NHASV+NHAGE+NHAGV)
      IF(IHARIND.GT.0) IHARIND=1

C.....IF HARMONIC ANALYSIS IS INCLUDED IN THE RUN, PROCESS HOT START
C     INFORMATION FOR IN PROGRESS HARMONIC ANALYSIS

      IF(IHARIND.EQ.1) THEN
         ITHAS=INT((THAS-STATIM)*(86400.D0/DT) + 0.5d0)
         ITHAF=INT((THAF-STATIM)*(86400.D0/DT) + 0.5d0)
         ITMV = ITHAF - (ITHAF-ITHAS)*FMV
         IHABEG=ITHAS+NHAINC

C.......IF HARMONIC ANALYSIS HAS ALREADY BEGUN, READ IN HOT START
C........HARMONIC ANALYSIS, MEAN AND SQUARE INFO

         IF(ITHSF.GT.ITHAS) THEN
            READ(IHOT,REC=IHOTSTP) ICHA
            IHOTSTP=IHOTSTP+1
         ENDIF

         IF(ITHSF.GE.IHABEG) THEN
            READ(IHOT,REC=IHOTSTP) INZ ; IHOTSTP = IHOTSTP + 1
            READ(IHOT,REC=IHOTSTP) INF ; IHOTSTP = IHOTSTP + 1
            READ(IHOT,REC=IHOTSTP) IMM ; IHOTSTP = IHOTSTP + 1
            READ(IHOT,REC=IHOTSTP) INP ; IHOTSTP = IHOTSTP + 1
            READ(IHOT,REC=IHOTSTP) INSTAE ; IHOTSTP = IHOTSTP + 1
            READ(IHOT,REC=IHOTSTP) INSTAV ; IHOTSTP = IHOTSTP + 1
            READ(IHOT,REC=IHOTSTP) IISTAE ; IHOTSTP = IHOTSTP + 1
            READ(IHOT,REC=IHOTSTP) IISTAV ; IHOTSTP = IHOTSTP + 1
            READ(IHOT,REC=IHOTSTP) IIGLOE ; IHOTSTP = IHOTSTP + 1
            READ(IHOT,REC=IHOTSTP) IIGLOV ; IHOTSTP = IHOTSTP + 1
            READ(IHOT,REC=IHOTSTP) IICALL ; IHOTSTP = IHOTSTP + 1
            READ(IHOT,REC=IHOTSTP) INFREQ ; IHOTSTP = IHOTSTP + 1

            DO I=1,INFREQ+INF
               READ(IHOT,REC=IHOTSTP) FNAM8(1) ; IHOTSTP = IHOTSTP + 1
               READ(IHOT,REC=IHOTSTP) FNAM8(2) ; IHOTSTP = IHOTSTP + 1

               INAMEFR(I) = FNAME1
               READ(IHOT,REC=IHOTSTP) IFREQ(I) ; IHOTSTP = IHOTSTP + 1
               READ(IHOT,REC=IHOTSTP) IFF(I) ; IHOTSTP = IHOTSTP + 1
               READ(IHOT,REC=IHOTSTP) IFACE(I) ; IHOTSTP = IHOTSTP + 1
            ENDDO

            READ(IHOT,REC=IHOTSTP) TIMEUD ; IHOTSTP = IHOTSTP + 1
            READ(IHOT,REC=IHOTSTP) ITUD ; IHOTSTP = IHOTSTP + 1

            DO I=1,IMM
               DO J=1,IMM
                  READ(IHOT,REC=IHOTSTP) HA(I,J) ; IHOTSTP = IHOTSTP + 1
               ENDDO
            ENDDO

            IF(NHASE.EQ.1) THEN
               DO J=1,INSTAE
                  DO I=1,IMM
                     READ(IHOT,REC=IHOTSTP) STAELV(I,J)
                     IHOTSTP=IHOTSTP+1
                  ENDDO
               ENDDO
            ENDIF

            IF(NHASV.EQ.1) THEN
               DO J=1,INSTAV
                  DO I=1,IMM
                     READ(IHOT,REC=IHOTSTP) STAULV(I,J)
                     IHOTSTP = IHOTSTP + 1
                     READ(IHOT,REC=IHOTSTP) STAVLV(I,J)
                     IHOTSTP = IHOTSTP + 1
                  ENDDO
               ENDDO
            ENDIF

            IF(NHAGE.EQ.1) THEN
               DO J=1,INP
                  DO I=1,IMM
                     READ(IHOT,REC=IHOTSTP) GLOELV(I,J)
                     IHOTSTP=IHOTSTP+1
                  ENDDO
               ENDDO
            ENDIF

            IF(NHAGV.EQ.1) THEN
               DO J=1,INP
                  DO I=1,IMM
                     READ(IHOT,REC=IHOTSTP) GLOULV(I,J)
                     IHOTSTP = IHOTSTP + 1
                     READ(IHOT,REC=IHOTSTP) GLOVLV(I,J)
                     IHOTSTP = IHOTSTP + 1
                  ENDDO
               ENDDO
            ENDIF

         ENDIF

         IF((FMV.GT.0.).AND.(INFREQ.GT.0).AND.(IM.EQ.0)) THEN !include means and variances
            IF(ITHSF.GT.ITMV) THEN
               READ(IHOT,REC=IHOTSTP) NTSTEPS
               IHOTSTP=IHOTSTP+1
               IF(NHAGE.EQ.1) THEN
                  DO I=1,INP
                     READ(IHOT,REC=IHOTSTP) ELAV(I)
                     IHOTSTP=IHOTSTP+1
                     READ(IHOT,REC=IHOTSTP) ELVA(I)
                     IHOTSTP=IHOTSTP+1
                  ENDDO
               ENDIF
               IF(NHAGV.EQ.1) THEN
                  DO I=1,INP
                     READ(IHOT,REC=IHOTSTP) XVELAV(I)
                     IHOTSTP=IHOTSTP+1
                     READ(IHOT,REC=IHOTSTP) YVELAV(I)
                     IHOTSTP=IHOTSTP+1
                     READ(IHOT,REC=IHOTSTP) XVELVA(I)
                     IHOTSTP=IHOTSTP+1
                     READ(IHOT,REC=IHOTSTP) YVELVA(I)
                     IHOTSTP=IHOTSTP+1
                  ENDDO
               ENDIF
            ENDIF
         ENDIF    ! charmv
      ENDIF     ! HARIND
#endif
C
C--Open All Local Hot Start files
C

      ALLOCATE ( LOC2(NPROC) )
      DO IPROC = 1,NPROC
         LOC2(IPROC) = 105 + (IPROC-1)
         LOCFN(1:14) = 'PE0000/'//FNAME(1:7)
         CALL IWRITE(LOCFN,3,6,IPROC-1)
         OPEN (LOC2(IPROC),FILE=LOCFN,ACCESS='DIRECT',RECL=8)
      ENDDO
C
C--Write out info to local hot start files
C
      DO IPROC = 1,NPROC
         IHOTSTP=1
         WRITE(LOC2(IPROC),REC=IHOTSTP) InputFileFmtVn ; IHOTSTP = IHOTSTP + 1
         WRITE(LOC2(IPROC),REC=IHOTSTP) IMHSF          ; IHOTSTP = IHOTSTP + 1
         WRITE(LOC2(IPROC),REC=IHOTSTP) TIMEHSF        ; IHOTSTP = IHOTSTP + 1
         WRITE(LOC2(IPROC),REC=IHOTSTP) ITHSF          ; IHOTSTP = IHOTSTP + 1
         WRITE(LOC2(IPROC),REC=IHOTSTP) NNODP(IPROC)   ; IHOTSTP = IHOTSTP + 1
         WRITE(LOC2(IPROC),REC=IHOTSTP) NELP(IPROC)    ; IHOTSTP = IHOTSTP + 1
         WRITE(LOC2(IPROC),REC=IHOTSTP) NNODP(IPROC)   ; IHOTSTP = IHOTSTP + 1
         WRITE(LOC2(IPROC),REC=IHOTSTP) NELP(IPROC)    ; IHOTSTP = IHOTSTP + 1

         DO I=1, NNODP(IPROC)
            INDX = ABS(IMAP_NOD_LG(I,IPROC))
            WRITE(LOC2(IPROC),REC=IHOTSTP) ETA1(INDX)
            IHOTSTP=IHOTSTP+1
         END DO

         DO I=1, NNODP(IPROC)
            INDX = ABS(IMAP_NOD_LG(I,IPROC))
            WRITE(LOC2(IPROC),REC=IHOTSTP) ETA2(INDX)
            IHOTSTP=IHOTSTP+1
         END DO

         DO I=1, NNODP(IPROC)
            INDX = ABS(IMAP_NOD_LG(I,IPROC))
            WRITE(LOC2(IPROC),REC=IHOTSTP) EtaDisc(INDX)
            IHOTSTP=IHOTSTP+1
         END DO

         DO I=1, NNODP(IPROC)
            INDX = ABS(IMAP_NOD_LG(I,IPROC))
            WRITE(LOC2(IPROC),REC=IHOTSTP) UU2(INDX)
            IHOTSTP=IHOTSTP+1
         END DO

         DO I=1, NNODP(IPROC)
            INDX = ABS(IMAP_NOD_LG(I,IPROC))
            WRITE(LOC2(IPROC),REC=IHOTSTP) VV2(INDX)
            IHOTSTP=IHOTSTP+1
         END DO

         IF(IM.EQ.10) THEN
            DO I=1, NNODP(IPROC)
               INDX = ABS(IMAP_NOD_LG(I,IPROC))
               WRITE(LOC2(IPROC),REC=IHOTSTP) CH1(INDX)
               IHOTSTP=IHOTSTP+1
            END DO
         ENDIF

         DO I=1, NNODP(IPROC)
            INDX = ABS(IMAP_NOD_LG(I,IPROC))
            WRITE(LOC2(IPROC),REC=IHOTSTP) NODECODE(INDX)
            IHOTSTP=IHOTSTP+1
         END DO

         DO I=1,NELP(IPROC)
            INDX=ABS(IMAP_EL_LG(I,IPROC))
            WRITE(LOC2(IPROC),REC=IHOTSTP) NOFF(INDX)
            IHOTSTP=IHOTSTP+1
         END DO

         WRITE(LOC2(IPROC),REC=IHOTSTP) IESTP  ; IHOTSTP = IHOTSTP + 1
         WRITE(LOC2(IPROC),REC=IHOTSTP) NSCOUE  ; IHOTSTP = IHOTSTP + 1

         WRITE(LOC2(IPROC),REC=IHOTSTP) IVSTP  ; IHOTSTP = IHOTSTP + 1
         WRITE(LOC2(IPROC),REC=IHOTSTP) NSCOUV  ; IHOTSTP = IHOTSTP + 1

         WRITE(LOC2(IPROC),REC=IHOTSTP) ICSTP  ; IHOTSTP = IHOTSTP + 1
         WRITE(LOC2(IPROC),REC=IHOTSTP) NSCOUC  ; IHOTSTP = IHOTSTP + 1

         WRITE(LOC2(IPROC),REC=IHOTSTP) IPSTP  ; IHOTSTP = IHOTSTP + 1
         WRITE(LOC2(IPROC),REC=IHOTSTP) IWSTP  ; IHOTSTP = IHOTSTP + 1
         WRITE(LOC2(IPROC),REC=IHOTSTP) NSCOUM  ; IHOTSTP = IHOTSTP + 1

         WRITE(LOC2(IPROC),REC=IHOTSTP) IGEP  ; IHOTSTP = IHOTSTP + 1
         WRITE(LOC2(IPROC),REC=IHOTSTP) NSCOUGE  ; IHOTSTP = IHOTSTP + 1

         WRITE(LOC2(IPROC),REC=IHOTSTP) IGVP  ; IHOTSTP = IHOTSTP + 1
         WRITE(LOC2(IPROC),REC=IHOTSTP) NSCOUGV  ; IHOTSTP = IHOTSTP + 1

         WRITE(LOC2(IPROC),REC=IHOTSTP) IGCP  ; IHOTSTP = IHOTSTP + 1
         WRITE(LOC2(IPROC),REC=IHOTSTP) NSCOUGC  ; IHOTSTP = IHOTSTP + 1

         WRITE(LOC2(IPROC),REC=IHOTSTP) IGPP  ; IHOTSTP = IHOTSTP + 1
         WRITE(LOC2(IPROC),REC=IHOTSTP) IGWP  ; IHOTSTP = IHOTSTP + 1
         WRITE(LOC2(IPROC),REC=IHOTSTP) NSCOUGW  ; IHOTSTP = IHOTSTP + 1

C  kmd48.33bc moved 3D hot start information to subroutine
C jgf49.43 refined check of IMHSF so that it picks up only
C     IM values that indicate 3D (and so we can use six integer IM values).
      IF ((IMHSF.EQ.1).OR.(IMHSF.EQ.11).OR.
     &    (IMHSF.EQ.21).OR.(IMHSF.EQ.31)) THEN
            CALL WriteHotStart3D(LOC2(IPROC),IHOTSTP,IPROC)
      ENDIF

#if 0
C
C....IF APPROPRIATE, WRITE OUT HOT START INFORMATION FOR IN PROGRESS HARMONIC ANALYSIS

c       IF((IHARIND.EQ.1).AND.(ITHSF.GT.ITHAS)) THEN
c         WRITE(LOC2(IPROC),REC=IHOTSTP+1) ICHA
c         IHOTSTP = IHOTSTP + 1
c         CALL HAHOUT(NP,NSTAE,NSTAV,NHASE,NHASV,NHAGE,NHAGV,
c    &                LOC2(IPROC),IHOTSTP)
c
c         IF(NHASE.EQ.1) CALL HAHOUTES(NSTAE,LOC2(IPROC),IHOTSTP)
c         IF(NHASV.EQ.1) CALL HAHOUTVS(NSTAV,LOC2(IPROC),IHOTSTP)
c         IF(NHAGE.EQ.1) CALL HAHOUTEG(MNP,LOC2(IPROC),IHOTSTP)
c         IF(NHAGV.EQ.1) CALL HAHOUTVG(MNP,LOC2(IPROC),IHOTSTP)
c         ENDIF
c
c       if(CHARMV) then
c         IF((IHARIND.EQ.1).AND.(ITHSF.GT.ITMV)) THEN
c           IHOTSTP=IHOTSTP+1
c           WRITE(LOC2(IPROC),REC=IHOTSTP) NTSTEPS
c           IF(NHAGE.EQ.1) THEN
c             DO I=1, NNODP(IPROC)
c               INDX = IMAP_NOD_LG(I,IPROC)
c               DO I=1,MNP
c                 WRITE(LOC2(IPROC),REC=IHOTSTP+1) ELAV(INDX)
c                 WRITE(LOC2(IPROC),REC=IHOTSTP+2) ELVA(INDX)
c                 IHOTSTP=IHOTSTP+2
c                 END DO
c             ENDIF
c           IF(NHAGV.EQ.1) THEN
c             DO I=1,NNODP(IPROC)
c               WRITE(LOC2(IPROC),REC=IHOTSTP+1) XVELAV(INDX)
c               WRITE(LOC2(IPROC),REC=IHOTSTP+2) YVELAV(INDX)
c               WRITE(LOC2(IPROC),REC=IHOTSTP+3) XVELVA(INDX)
c               WRITE(LOC2(IPROC),REC=IHOTSTP+4) YVELVA(INDX)
c               IHOTSTP=IHOTSTP+4
c               END DO
c             ENDIF
c           ENDIF
c         ENDIF
#endif

      ENDDO
C
C--Close Global file and all the Local Files
C
      CLOSE (IHOT)
      DO IPROC=1, NPROC
         CLOSE (LOC2(IPROC))
      ENDDO
C
      IF(ALLOCATED(LOC2)) DEALLOCATE ( LOC2 )
      nbytes = 4*nproc
      call memory_dealloc(nbytes)
      IF(ALLOCATED( ETA1 ))DEALLOCATE ( ETA1  )
      IF(ALLOCATED( ETA2 ))DEALLOCATE ( ETA2  )
      IF(ALLOCATED( EtaDisc ))DEALLOCATE ( EtaDisc )
      IF(ALLOCATED( UU2 ))DEALLOCATE ( UU2 )
      IF(ALLOCATED( VV2 ))DEALLOCATE ( VV2 )
      IF(ALLOCATED( NODECODE ))DEALLOCATE ( NODECODE )
      IF(ALLOCATED( CH1 ))DEALLOCATE ( CH1 )
      nbytes = 7*mnp*8
      call memory_dealloc(nbytes)
      IF(ALLOCATED(NOFF))DEALLOCATE ( NOFF )
      nbytes = 6*mne
      call memory_dealloc(nbytes)
      IF(ALLOCATED( DUU ))DEALLOCATE ( DUU  )
      IF(ALLOCATED( DUV ))DEALLOCATE ( DUV  )
      IF(ALLOCATED( DVV ))DEALLOCATE ( DVV )
      nbytes = 3*mnp*8
      call memory_dealloc(nbytes)
      IF(ALLOCATED( UU )) DEALLOCATE ( UU )
      IF(ALLOCATED( VV )) DEALLOCATE ( VV )
      nbytes = 2*mnp*8
      call memory_dealloc(nbytes)
      IF(ALLOCATED  ( BSX )) DEALLOCATE ( BSX )
      IF(ALLOCATED  ( BSY )) DEALLOCATE ( BSY )
      nbytes = 2*mnp*8
      call memory_dealloc(nbytes)
      IF(ALLOCATED  ( WZ )) DEALLOCATE ( WZ )
      IF(ALLOCATED  ( q20 )) DEALLOCATE (q20 )
      nbytes = (mnp*nfen*8) + (mnp*nfen*8)
      call memory_dealloc(nbytes)
      IF(ALLOCATED  ( RealQ ))  DEALLOCATE ( RealQ)
      IF(ALLOCATED  ( ImagQ ))  DEALLOCATE ( ImagQ)
      nbytes = (mnp*nfen*8) + (mnp*nfen*8)
      call memory_dealloc(nbytes)
      IF(ALLOCATED  ( l )) DEALLOCATE ( l )
      IF(ALLOCATED  ( SigT )) DEALLOCATE ( SigT )
      nbytes = (mnp*nfen*8) + (mnp*nfen*8)
      call memory_dealloc(nbytes)
      IF(ALLOCATED ( Sal)) DEALLOCATE ( Sal  )
      IF(ALLOCATED ( Temp )) DEALLOCATE (  Temp )
      nbytes = (mnp*nfen*8) + (mnp*nfen*8)
      call memory_dealloc(nbytes)
      call memory_status()
C
      RETURN
 1001 FORMAT('ERROR: The hot start file')
 1010 FORMAT(' File ',A60,/,' WAS NOT FOUND!  ADCPrep Terminated!!!',/)
 1011 FORMAT(' File ',A60,/,' WAS FOUND!  Opening & Processing file',/)
 1012 FORMAT('was a nonmatching version')
 1005 FORMAT('exists but cannot be opened.')
 9973 FORMAT(/,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
      END SUBROUTINE HOTLOCALIZE


      SUBROUTINE HOTGLOBALIZE()
      USE PRE_GLOBAL
      use presizes; use memory_usage
C
C---------------------------------------------------------------------------C
C                     written 10/11/01 by RL                                C
C             started mods for harmonic analysis and 3D RL 5/22/03          C
C         jgf Updated for v45.06 09/07/2005 not incl. harmonic or 3D        C
C         kmd48.33bc updated with 3D hot start                              C
C                                                                           C
C  This routine reads the global hot start file (either fort.67 or fort.68) C
C  and writes local hot start files of the same format.                     C
C                                                                           C
C---------------------------------------------------------------------------C
C
      IMPLICIT NONE
      integer :: nbytes = 0
      LOGICAL FOUND
      INTEGER I,J,IPROC,INDX,IHOTSTP, not_active
      INTEGER IMHSF,ITHSF,IVALUE,IDUMY, NH, N
      CHARACTER FNAME*60,LOCFN*14
      CHARACTER*16 FNAME1
      CHARACTER*8 FNAM8(2)
      EQUIVALENCE (FNAM8(1),FNAME1)

      INTEGER,ALLOCATABLE  :: LOC2(:),NOFF(:), domA(:)
      REAL(SZ),ALLOCATABLE :: ETA1(:),ETA2(:),EtaDisc(:),
     &    UU2(:),VV2(:),CH1(:)
      REAL(8) TIMEHSF, RVALUE
      integer :: InputFileFmtVn, NP_G_IN, NE_G_IN, NP_A_IN, NE_A_IN

#if 0
      ! vjp 2006/9/30 not supporting harmonic analysis or C3D yet
      INTEGER INZ,INF,IMM,INP,INSTAE,INSTAV,IISTAE,IISTAV,IIGLOE,IIGLOV,
     &                                       IICALL,INFREQ,ITUD,NTSTEPS
      INTEGER IHARIND,ITHAS,ITHAF,ITMV,IHABEG,ICHA
      CHARACTER*10,ALLOCATABLE     ::  INAMEFR(:)
      REAL(8)  TIMEUD
      REAL(SZ),ALLOCATABLE ::  HA(:,:)
      REAL(SZ),ALLOCATABLE ::  ELAV(:),ELVA(:),XVELAV(:),XVELVA(:),
     &                                         YVELAV(:),YVELVA(:)
      REAL(SZ),ALLOCATABLE ::  IFREQ(:),IFF(:),IFACE(:)
      REAL(SZ),ALLOCATABLE ::  GLOELV(:,:)
      REAL(SZ),ALLOCATABLE ::  GLOULV(:,:),GLOVLV(:,:)
      REAL(SZ),ALLOCATABLE ::  STAELV(:,:)
      REAL(SZ),ALLOCATABLE ::  STAULV(:,:),STAVLV(:,:)
#endif
C
C--   Open Appropriate Hot Start File based on the value of IHOT from
C--   the fort.15 file
C
      write(*,*) "enter IHOT: "
      read(*,*) IHOT
      IF(IHOT.EQ.67) FNAME='fort.67'
      IF(IHOT.EQ.68) FNAME='fort.68'

C
C--Open All Local Hot Start files
C
      ALLOCATE ( LOC2(NPROC) )
      DO IPROC = 1,NPROC
         LOC2(IPROC) = 105 + (IPROC-1)
         LOCFN(1:14) = 'PE0000/'//FNAME(1:7)
         CALL IWRITE(LOCFN,3,6,IPROC-1)
         INQUIRE(FILE=LOCFN,EXIST=FOUND)
         IF (FOUND) THEN
           WRITE(*,1011) LOCFN
           OPEN (LOC2(IPROC),FILE=LOCFN,ACCESS='DIRECT',RECL=8)
         ELSE
          WRITE(*,1010) FNAME
          STOP
         ENDIF
      ENDDO
C
C Allocate local work arrays
C
      MNP  =  nnodg    !  global number of nodes    ( read from fort.18 )
      print *, "MNP =", MNP

      nbytes = 4*nproc
      call memory_alloc(nbytes)
      ALLOCATE ( ETA1(MNP),ETA2(MNP),EtaDisc(MNP),UU2(MNP),
     &           VV2(MNP),NODECODE(MNP),CH1(MNP) )
      nbytes = 7*mnp

      MNE  =  nelg     !  global number of elements ( read from fort.18 )
      print *, "MNE =", MNE

      call memory_alloc(nbytes)
      ALLOCATE ( NOFF(MNE) )
      nbytes = 4*mne
      call memory_alloc(nbytes)

#if HA
      ! vjp 2006/9/30 not supporting harmonic analysis or C3D yet
      ALLOCATE ( HA(2*MNHARF,2*MNHARF) )
      nbytes = 32*mnharf
      call memory_alloc(nbytes)
      ALLOCATE ( GLOELV(2*MNHARF,MNP) )
      nbytes = 16*mnharf*mnp
      call memory_alloc(nbytes)
      ALLOCATE ( GLOULV(2*MNHARF,MNP),GLOVLV(2*MNHARF,MNP) )
      nbytes = 32*mnharf*mnp
      call memory_alloc(nbytes)
      ALLOCATE ( STAELV(2*MNHARF,MNSTAE) )
      nbytes = 16*mnharf*mnstae
      call memory_alloc(nbytes)
      ALLOCATE ( STAULV(2*MNHARF,MNSTAV),STAVLV(2*MNHARF,MNSTAV) )
      nbytes = 16*mnharf*mnstav
      call memory_alloc(nbytes)
      ALLOCATE ( ELAV(MNP),ELVA(MNP) )
      nbytes = 16*mnp
      call memory_alloc(nbytes)
      ALLOCATE ( XVELAV(MNP),XVELVA(MNP),YVELAV(MNP),YVELVA(MNP) )
      nbytes = 32*mnp
      call memory_alloc(nbytes)
      ALLOCATE ( IFREQ(MNHARF),IFF(MNHARF),IFACE(MNHARF) )
      nbytes = 12*mnharf
      call memory_alloc(nbytes)
      ALLOCATE ( INAMEFR(MNHARF) )
      nbytes = 4*mnharf
      call memory_alloc(nbytes)
#endif

C
C--Read info from local hot start files
C
      DO IPROC = 1,NPROC
         IHOTSTP=1
         READ(LOC2(IPROC),REC=IHOTSTP) InputFileFmtVn ; IHOTSTP = IHOTSTP + 1
         READ(LOC2(IPROC),REC=IHOTSTP) IMHSF          ; IHOTSTP = IHOTSTP + 1
         READ(LOC2(IPROC),REC=IHOTSTP) TIMEHSF        ; IHOTSTP = IHOTSTP + 1
         READ(LOC2(IPROC),REC=IHOTSTP) ITHSF          ; IHOTSTP = IHOTSTP + 1
         READ(LOC2(IPROC),REC=IHOTSTP) IDUMY          ; IHOTSTP = IHOTSTP + 1
         READ(LOC2(IPROC),REC=IHOTSTP) IDUMY          ; IHOTSTP = IHOTSTP + 1
         READ(LOC2(IPROC),REC=IHOTSTP) IDUMY          ; IHOTSTP = IHOTSTP + 1
         READ(LOC2(IPROC),REC=IHOTSTP) IDUMY          ; IHOTSTP = IHOTSTP + 1

         DO I=1, NNODP(IPROC)
            INDX = IMAP_NOD_LG(I,IPROC)
            READ(LOC2(IPROC),REC=IHOTSTP) RVALUE
            IHOTSTP=IHOTSTP+1
            IF (INDX > 0) ETA1(INDX) = RVALUE
         END DO

         DO I=1, NNODP(IPROC)
            INDX = IMAP_NOD_LG(I,IPROC)
            READ(LOC2(IPROC),REC=IHOTSTP) RVALUE
            IHOTSTP=IHOTSTP+1
            IF (INDX > 0) ETA2(INDX) = RVALUE
         END DO

         DO I=1, NNODP(IPROC)
            INDX = IMAP_NOD_LG(I,IPROC)
            READ(LOC2(IPROC),REC=IHOTSTP) RVALUE
            IHOTSTP=IHOTSTP+1
            IF (INDX > 0) EtaDisc(INDX) = RVALUE
         END DO

         DO I=1, NNODP(IPROC)
            INDX = IMAP_NOD_LG(I,IPROC)
            READ(LOC2(IPROC),REC=IHOTSTP) RVALUE
            IHOTSTP=IHOTSTP+1
            IF (INDX > 0) UU2(INDX) = RVALUE
         END DO

         DO I=1, NNODP(IPROC)
            INDX = IMAP_NOD_LG(I,IPROC)
            READ(LOC2(IPROC),REC=IHOTSTP) RVALUE
            IHOTSTP=IHOTSTP+1
            IF (INDX > 0) VV2(INDX) = RVALUE
         END DO

         IF(IM.EQ.10) THEN
            DO I=1, NNODP(IPROC)
               INDX = IMAP_NOD_LG(I,IPROC)
               READ(LOC2(IPROC),REC=IHOTSTP) RVALUE
               IHOTSTP=IHOTSTP+1
               IF (INDX > 0) CH1(INDX) = RVALUE
            END DO
         ENDIF

         DO I=1, NNODP(IPROC)
            INDX = IMAP_NOD_LG(I,IPROC)
            READ(LOC2(IPROC),REC=IHOTSTP) IVALUE
            IHOTSTP=IHOTSTP+1
            IF (INDX > 0) NODECODE(INDX) = IVALUE
         END DO

         DO I=1,NELP(IPROC)
            INDX = IMAP_EL_LG(I,IPROC)
            READ(LOC2(IPROC),REC=IHOTSTP) IVALUE
            IHOTSTP=IHOTSTP+1
            IF (INDX > 0) NOFF(INDX) = IVALUE
         END DO

         READ(LOC2(IPROC),REC=IHOTSTP) IESTP ; IHOTSTP = IHOTSTP + 1
         READ(LOC2(IPROC),REC=IHOTSTP) NSCOUE ; IHOTSTP = IHOTSTP + 1

         READ(LOC2(IPROC),REC=IHOTSTP) IVSTP ; IHOTSTP = IHOTSTP + 1
         READ(LOC2(IPROC),REC=IHOTSTP) NSCOUV ; IHOTSTP = IHOTSTP + 1

         READ(LOC2(IPROC),REC=IHOTSTP) ICSTP ; IHOTSTP = IHOTSTP + 1
         READ(LOC2(IPROC),REC=IHOTSTP) NSCOUC ; IHOTSTP = IHOTSTP + 1

         READ(LOC2(IPROC),REC=IHOTSTP) IPSTP ; IHOTSTP = IHOTSTP + 1
         READ(LOC2(IPROC),REC=IHOTSTP) IWSTP ; IHOTSTP = IHOTSTP + 1
         READ(LOC2(IPROC),REC=IHOTSTP) NSCOUM ; IHOTSTP = IHOTSTP + 1

         READ(LOC2(IPROC),REC=IHOTSTP) IGEP ; IHOTSTP = IHOTSTP + 1
         READ(LOC2(IPROC),REC=IHOTSTP) NSCOUGE ; IHOTSTP = IHOTSTP + 1

         READ(LOC2(IPROC),REC=IHOTSTP) IGVP ; IHOTSTP = IHOTSTP + 1
         READ(LOC2(IPROC),REC=IHOTSTP) NSCOUGV ; IHOTSTP = IHOTSTP + 1

         READ(LOC2(IPROC),REC=IHOTSTP) IGCP ; IHOTSTP = IHOTSTP + 1
         READ(LOC2(IPROC),REC=IHOTSTP) NSCOUGC ; IHOTSTP = IHOTSTP + 1

         READ(LOC2(IPROC),REC=IHOTSTP) IGPP ; IHOTSTP = IHOTSTP + 1
         READ(LOC2(IPROC),REC=IHOTSTP) IGWP ; IHOTSTP = IHOTSTP + 1
         READ(LOC2(IPROC),REC=IHOTSTP) NSCOUGW ; IHOTSTP = IHOTSTP + 1
C
C   kmd48.33bc add information for 3D hot start
C     jgf46.02 read in 3D hotstart data if appropriate
         IF (C3D) THEN
            CALL ReadHotStart3DGlobal(LOC2(IPROC),IHOTSTP,IPROC)
         ENDIF
#if HA
C
C....IF APPROPRIATE, WRITE OUT HOT START INFORMATION FOR IN PROGRESS HARMONIC ANALYSIS

c       IF((IHARIND.EQ.1).AND.(ITHSF.GT.ITHAS)) THEN
c         READ(LOC2(IPROC),REC=IHOTSTP+1) ICHA
c         IHOTSTP = IHOTSTP + 1
c         CALL HAHOUT(NP,NSTAE,NSTAV,NHASE,NHASV,NHAGE,NHAGV,
c    &                LOC2(IPROC),IHOTSTP)
c
c         IF(NHASE.EQ.1) CALL HAHOUTES(NSTAE,LOC2(IPROC),IHOTSTP)
c         IF(NHASV.EQ.1) CALL HAHOUTVS(NSTAV,LOC2(IPROC),IHOTSTP)
c         IF(NHAGE.EQ.1) CALL HAHOUTEG(MNP,LOC2(IPROC),IHOTSTP)
c         IF(NHAGV.EQ.1) CALL HAHOUTVG(MNP,LOC2(IPROC),IHOTSTP)
c         ENDIF
c
c       if(CHARMV) then
c         IF((IHARIND.EQ.1).AND.(ITHSF.GT.ITMV)) THEN
c           IHOTSTP=IHOTSTP+1
c           READ(LOC2(IPROC),REC=IHOTSTP) NTSTEPS
c           IF(NHAGE.EQ.1) THEN
c             DO I=1, NNODP(IPROC)
c               INDX = IMAP_NOD_LG(I,IPROC)
c               DO I=1,MNP
c                 READ(LOC2(IPROC),REC=IHOTSTP+1) ELAV(INDX)
c                 READ(LOC2(IPROC),REC=IHOTSTP+2) ELVA(INDX)
c                 IHOTSTP=IHOTSTP+2
c                 END DO
c             ENDIF
c           IF(NHAGV.EQ.1) THEN
c             DO I=1,NNODP(IPROC)
c               READ(LOC2(IPROC),REC=IHOTSTP+1) XVELAV(INDX)
c               READ(LOC2(IPROC),REC=IHOTSTP+2) YVELAV(INDX)
c               READ(LOC2(IPROC),REC=IHOTSTP+3) XVELVA(INDX)
c               READ(LOC2(IPROC),REC=IHOTSTP+4) YVELVA(INDX)
c               IHOTSTP=IHOTSTP+4
c               END DO
c             ENDIF
c           ENDIF
c         ENDIF
#endif
          CLOSE (LOC2(IPROC))

      ENDDO

C-----------------------------------------------------------------------
C--   Write info to global hot start files
C-----------------------------------------------------------------------

      OPEN(IHOT,FILE=trim(FNAME),ACCESS='DIRECT',RECL=8)
      print *, "opening global hotstart file"

      IHOTSTP=1
      WRITE(IHOT,REC=IHOTSTP) InputFileFmtVn ; IHOTSTP = IHOTSTP + 1
      WRITE(IHOT,REC=IHOTSTP) IMHSF        ; IHOTSTP = IHOTSTP + 1
      WRITE(IHOT,REC=IHOTSTP) TIMEHSF      ; IHOTSTP = IHOTSTP + 1
      WRITE(IHOT,REC=IHOTSTP) ITHSF        ; IHOTSTP = IHOTSTP + 1
      WRITE(IHOT,REC=IHOTSTP) MNP          ; IHOTSTP = IHOTSTP + 1
      WRITE(IHOT,REC=IHOTSTP) MNE          ; IHOTSTP = IHOTSTP + 1
      WRITE(IHOT,REC=IHOTSTP) MNP          ; IHOTSTP = IHOTSTP + 1
      WRITE(IHOT,REC=IHOTSTP) MNE          ; IHOTSTP = IHOTSTP + 1

      DO I=1,MNP
         WRITE(IHOT,REC=IHOTSTP) ETA1(I)
         IHOTSTP = IHOTSTP + 1
      END DO
      DO I=1,MNP
         WRITE(IHOT,REC=IHOTSTP) ETA2(I)
         IHOTSTP = IHOTSTP + 1
      END DO
      DO I=1,MNP
         WRITE(IHOT,REC=IHOTSTP) EtaDisc(I)
         IHOTSTP = IHOTSTP + 1
      END DO
      DO I=1,MNP
         WRITE(IHOT,REC=IHOTSTP) UU2(I)
         IHOTSTP = IHOTSTP + 1
      END DO
      DO I=1,MNP
         WRITE(IHOT,REC=IHOTSTP) VV2(I)
         IHOTSTP = IHOTSTP + 1
      END DO
      IF(IM.EQ.10) THEN
        DO I=1,MNP
           WRITE(IHOT,REC=IHOTSTP) CH1(I)
           IHOTSTP=IHOTSTP+1
        END DO
      ENDIF
      DO I=1,MNP
         WRITE(IHOT,REC=IHOTSTP) NODECODE(I)
         IHOTSTP=IHOTSTP+1
      END DO
      DO I=1,MNE
         WRITE(IHOT,REC=IHOTSTP) NOFF(I)
         IHOTSTP=IHOTSTP+1
      END DO

      WRITE(IHOT,REC=IHOTSTP) IESTP ; IHOTSTP = IHOTSTP + 1
      WRITE(IHOT,REC=IHOTSTP) NSCOUE ; IHOTSTP = IHOTSTP + 1

      WRITE(IHOT,REC=IHOTSTP) IVSTP ; IHOTSTP = IHOTSTP + 1
      WRITE(IHOT,REC=IHOTSTP) NSCOUV ; IHOTSTP = IHOTSTP + 1

      WRITE(IHOT,REC=IHOTSTP) ICSTP ; IHOTSTP = IHOTSTP + 1
      WRITE(IHOT,REC=IHOTSTP) NSCOUC ; IHOTSTP = IHOTSTP + 1

      WRITE(IHOT,REC=IHOTSTP) IPSTP ; IHOTSTP = IHOTSTP + 1
      WRITE(IHOT,REC=IHOTSTP) IWSTP ; IHOTSTP = IHOTSTP + 1
      WRITE(IHOT,REC=IHOTSTP) NSCOUM ; IHOTSTP = IHOTSTP + 1

      WRITE(IHOT,REC=IHOTSTP) IGEP ; IHOTSTP = IHOTSTP + 1
      WRITE(IHOT,REC=IHOTSTP) NSCOUGE ; IHOTSTP = IHOTSTP + 1

      WRITE(IHOT,REC=IHOTSTP) IGVP ; IHOTSTP = IHOTSTP + 1
      WRITE(IHOT,REC=IHOTSTP) NSCOUGV ; IHOTSTP = IHOTSTP + 1

      WRITE(IHOT,REC=IHOTSTP) IGCP ; IHOTSTP = IHOTSTP + 1
      WRITE(IHOT,REC=IHOTSTP) NSCOUGC ; IHOTSTP = IHOTSTP + 1

      WRITE(IHOT,REC=IHOTSTP) IGPP ; IHOTSTP = IHOTSTP + 1
      WRITE(IHOT,REC=IHOTSTP) IGWP ; IHOTSTP = IHOTSTP + 1
      WRITE(IHOT,REC=IHOTSTP) NSCOUGW ; IHOTSTP = IHOTSTP + 1
C

       ! vjp 2006/9/30 not supporting harmonic analysis or C3D yet
C   kmd48.33bc add information for 3D hot start
C     jgf46.02 Write in 3D hotstart data if appropriate
      IF (C3D) THEN
         CALL WriteHotStart3DGlobal(IHOT,IHOTSTP,IPROC)
      ENDIF
C
#if 0
C.....DETERMINE HARMONIC ANALYSIS PARAMETERS

      IHARIND=NHARFR*(NHASE+NHASV+NHAGE+NHAGV)
      IF(IHARIND.GT.0) IHARIND=1

C.....IF HARMONIC ANALYSIS IS INCLUDED IN THE RUN, PROCESS HOT START INFORMATION FOR
C.....IN PROGRESS HARMONIC ANALYSIS

      IF(IHARIND.EQ.1) THEN
         ITHAS=INT((THAS-STATIM)*(86400.D0/DT) + 0.5d0)
         ITHAF=INT((THAF-STATIM)*(86400.D0/DT) + 0.5d0)
         ITMV = ITHAF - (ITHAF-ITHAS)*FMV
         IHABEG=ITHAS+NHAINC

C.......IF HARMONIC ANALYSIS HAS ALREADY BEGUN, READ IN HOT START
C........HARMONIC ANALYSIS, MEAN AND SQUARE INFO

         IF(ITHSF.GT.ITHAS) THEN
            WRITE(IHOT,REC=IHOTSTP) ICHA
            IHOTSTP=IHOTSTP+1
         ENDIF

         IF(ITHSF.GE.IHABEG) THEN
            WRITE(IHOT,REC=IHOTSTP) INZ ; IHOTSTP = IHOTSTP + 1
            WRITE(IHOT,REC=IHOTSTP) INF ;  IHOTSTP = IHOTSTP + 1
            WRITE(IHOT,REC=IHOTSTP) IMM ;  IHOTSTP = IHOTSTP + 1
            WRITE(IHOT,REC=IHOTSTP) INP ;  IHOTSTP = IHOTSTP + 1
            WRITE(IHOT,REC=IHOTSTP) INSTAE ;  IHOTSTP = IHOTSTP + 1
            WRITE(IHOT,REC=IHOTSTP) INSTAV ;  IHOTSTP = IHOTSTP + 1
            WRITE(IHOT,REC=IHOTSTP) IISTAE ;  IHOTSTP = IHOTSTP + 1
            WRITE(IHOT,REC=IHOTSTP) IISTAV ;  IHOTSTP = IHOTSTP + 1
            WRITE(IHOT,REC=IHOTSTP) IIGLOE ;  IHOTSTP = IHOTSTP + 1
            WRITE(IHOT,REC=IHOTSTP) IIGLOV ;  IHOTSTP = IHOTSTP + 1
            WRITE(IHOT,REC=IHOTSTP) IICALL ;  IHOTSTP = IHOTSTP + 1
            WRITE(IHOT,REC=IHOTSTP) INFREQ ;  IHOTSTP = IHOTSTP + 1

            DO I=1,INFREQ+INF
               WRITE(IHOT,REC=IHOTSTP) FNAM8(1) ;  IHOTSTP = IHOTSTP + 1
               WRITE(IHOT,REC=IHOTSTP) FNAM8(2) ;  IHOTSTP = IHOTSTP + 1
               INAMEFR(I) = FNAME1
               WRITE(IHOT,REC=IHOTSTP) IFREQ(I) ;  IHOTSTP = IHOTSTP + 1
               WRITE(IHOT,REC=IHOTSTP) IFF(I) ;  IHOTSTP = IHOTSTP + 1
               WRITE(IHOT,REC=IHOTSTP) IFACE(I) ;  IHOTSTP = IHOTSTP + 1
            ENDDO

            WRITE(IHOT,REC=IHOTSTP) TIMEUD ;  IHOTSTP = IHOTSTP + 1
            WRITE(IHOT,REC=IHOTSTP) ITUD ;  IHOTSTP = IHOTSTP + 1

            DO I=1,IMM
               DO J=1,IMM
                  WRITE(IHOT,REC=IHOTSTP) HA(I,J)
                  IHOTSTP = IHOTSTP + 1
               ENDDO
            ENDDO

            IF(NHASE.EQ.1) THEN
               DO J=1,INSTAE
                  DO I=1,IMM
                     WRITE(IHOT,REC=IHOTSTP) STAELV(I,J)
                     IHOTSTP=IHOTSTP+1
                  ENDDO
               ENDDO
            ENDIF

            IF(NHASV.EQ.1) THEN
               DO J=1,INSTAV
                  DO I=1,IMM
                     WRITE(IHOT,REC=IHOTSTP) STAULV(I,J)
                     IHOTSTP=IHOTSTP+1
                     WRITE(IHOT,REC=IHOTSTP) STAVLV(I,J)
                     IHOTSTP=IHOTSTP+1
                  ENDDO
               ENDDO
            ENDIF

            IF(NHAGE.EQ.1) THEN
               DO J=1,INP
                  DO I=1,IMM
                     WRITE(IHOT,REC=IHOTSTP) GLOELV(I,J)
                     IHOTSTP=IHOTSTP+1
                  ENDDO
               ENDDO
            ENDIF

            IF(NHAGV.EQ.1) THEN
               DO J=1,INP
                  DO I=1,IMM
                     WRITE(IHOT,REC=IHOTSTP) GLOULV(I,J)
                     IHOTSTP = IHOTSTP + 1
                     WRITE(IHOT,REC=IHOTSTP) GLOVLV(I,J)
                     IHOTSTP = IHOTSTP + 1
                  ENDDO
               ENDDO
            ENDIF

         ENDIF

         IF((FMV.GT.0.).AND.(INFREQ.GT.0).AND.(IM.EQ.0)) THEN !include means and variances
            IF(ITHSF.GT.ITMV) THEN
               WRITE(IHOT,REC=IHOTSTP) NTSTEPS
               IHOTSTP=IHOTSTP+1
               IF(NHAGE.EQ.1) THEN
                  DO I=1,INP
                     WRITE(IHOT,REC=IHOTSTP) ELAV(I)
                     IHOTSTP=IHOTSTP+1
                     WRITE(IHOT,REC=IHOTSTP) ELVA(I)
                     IHOTSTP=IHOTSTP+1
                  ENDDO
               ENDIF
               IF(NHAGV.EQ.1) THEN
                  DO I=1,INP
                     WRITE(IHOT,REC=IHOTSTP) XVELAV(I)
                     IHOTSTP=IHOTSTP+1
                     WRITE(IHOT,REC=IHOTSTP) YVELAV(I)
                     IHOTSTP=IHOTSTP+1
                     WRITE(IHOT,REC=IHOTSTP) XVELVA(I)
                     IHOTSTP=IHOTSTP+1
                     WRITE(IHOT,REC=IHOTSTP) YVELVA(I)
                     IHOTSTP=IHOTSTP+1
                  ENDDO
               ENDIF
            ENDIF
         ENDIF    ! charmv
      ENDIF     ! HARIND
#endif

C
C--Close Global file and all the Local Files
C
      CLOSE (IHOT)
C
      DEALLOCATE ( LOC2 )
      nbytes = 4*nproc
      call memory_dealloc(nbytes)
      DEALLOCATE ( ETA1, ETA2, EtaDisc, UU2, VV2, NODECODE, CH1 )
      nbytes = 7*mnp
      call memory_dealloc(nbytes)
      DEALLOCATE ( NOFF )
      nbytes = 6*mne
      call memory_dealloc(nbytes)
      call memory_status()
C
      RETURN
 1001 FORMAT('ERROR: The hot start file')
 1010 FORMAT(' File ',A60,/,' WAS NOT FOUND!  ADCPrep Terminated!!!',/)
 1011 FORMAT(' File ',A60,/,' WAS FOUND!  Opening & Processing file',/)
 1012 FORMAT('was a nonmatching version')
 1005 FORMAT('exists but cannot be opened.')
 9973 FORMAT(/,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
      END SUBROUTINE HOTGLOBALIZE






C     ----------------------------------------------------------------------
C          S U B R O U T I N E     R E A D  H O T  S T A R T  3 D
C     ----------------------------------------------------------------------
C
C     jgf46.02 This subroutine supports PREP67_68. It reads in the 3D
C     section of the full domain hot start file.
C
C     ----------------------------------------------------------------------
      SUBROUTINE ReadHotStart3D(UnitNumber,FilePosition)
C     ----------------------------------------------------------------------
      USE PRE_GLOBAL
      IMPLICIT NONE
      INTEGER, intent(in) :: UnitNumber ! i/o unit of full domain file
      INTEGER, intent(inout) :: FilePosition ! position in binary file
      INTEGER IHOTSTP,NH,K
C
C     Start reading in the data

!kmd Added information for 3D hotstart
      ALLOCATE ( DUU(MNP),DUV(MNP),DVV(MNP))
      ALLOCATE ( UU(MNP),VV(MNP))
      ALLOCATE ( BSX(MNP),BSY(MNP))
!kmd end of additions

        PRINT *, "NFEN = ", NFEN
        IHOT=UnitNumber
        IHOTSTP=FilePosition

        PRINT *, "How many layers need to be evaluated:"
        READ *, NFEN

        PRINT *, "made it to 3D portion of code"
        PRINT *, "MNP = ", MNP
        PRINT *, "NFEN = ", NFEN
        PRINT *, "IHOTSTP = ", IHOTSTP

!kmd Added information for 3D hotstart
      ALLOCATE ( WZ(MNP,NFEN), q20(MNP,NFEN))
      ALLOCATE ( RealQ(MNP,NFEN), ImagQ(MNP,NFEN))
      ALLOCATE ( l(MNP,NFEN), SigT(MNP,NFEN))
      ALLOCATE ( Sal(MNP,NFEN), Temp(MNP,NFEN))
!kmd end of additions

        READ(IHOT,REC=IHOTSTP) IDEN  ;  IHOTSTP = IHOTSTP + 1
        READ(IHOT,REC=IHOTSTP) N3DSD ;  IHOTSTP = IHOTSTP + 1
        READ(IHOT,REC=IHOTSTP) I3DSDRec ;  IHOTSTP = IHOTSTP + 1
        READ(IHOT,REC=IHOTSTP) N3DSV ;  IHOTSTP = IHOTSTP + 1
        READ(IHOT,REC=IHOTSTP) I3DSVRec ;  IHOTSTP = IHOTSTP + 1
        READ(IHOT,REC=IHOTSTP) N3DST ;  IHOTSTP = IHOTSTP + 1
        READ(IHOT,REC=IHOTSTP) I3DSTRec ;  IHOTSTP = IHOTSTP + 1
        READ(IHOT,REC=IHOTSTP) N3DGD ;  IHOTSTP = IHOTSTP + 1
        READ(IHOT,REC=IHOTSTP) I3DGDRec ;  IHOTSTP = IHOTSTP + 1
        READ(IHOT,REC=IHOTSTP) N3DGV ;  IHOTSTP = IHOTSTP + 1
        READ(IHOT,REC=IHOTSTP) I3DGVRec ;  IHOTSTP = IHOTSTP + 1
        READ(IHOT,REC=IHOTSTP) N3DGT ;  IHOTSTP = IHOTSTP + 1
        READ(IHOT,REC=IHOTSTP) I3DGTRec ;  IHOTSTP = IHOTSTP + 1

        DO NH=1,MNP
          READ(IHOT,REC=IHOTSTP) DUU(NH)
          IHOTSTP=IHOTSTP+1
        END DO
        DO NH=1,MNP
          READ(IHOT,REC=IHOTSTP) DUV(NH)
          IHOTSTP=IHOTSTP+1
        END DO
        DO NH=1,MNP
          READ(IHOT,REC=IHOTSTP) DVV(NH)
          IHOTSTP=IHOTSTP+1
        END DO
        DO NH=1,MNP
          READ(IHOT,REC=IHOTSTP) UU(NH)
          IHOTSTP=IHOTSTP+1
        END DO
        DO NH=1,MNP
          READ(IHOT,REC=IHOTSTP) VV(NH)
          IHOTSTP=IHOTSTP+1
        END DO
        DO NH=1,MNP
          READ(IHOT,REC=IHOTSTP) BSX(NH)
          IHOTSTP=IHOTSTP+1
        END DO
        DO NH=1,MNP
          READ(IHOT,REC=IHOTSTP) BSY(NH)
          IHOTSTP=IHOTSTP+1
        ENDDO

        DO K=1,NFEN
          DO NH=1,MNP
            READ(IHOT,REC=IHOTSTP) RealQ(NH,K)
            IHOTSTP=IHOTSTP+1
          END DO
        END DO
        DO K=1,NFEN
          DO NH=1,MNP
            READ(IHOT,REC=IHOTSTP) ImagQ(NH,K)
            IHOTSTP=IHOTSTP+1
          END DO
        END DO
        DO K=1,NFEN
          DO NH=1,MNP
            READ(IHOT,REC=IHOTSTP) WZ(NH,K)
            IHOTSTP=IHOTSTP+1
          END DO
        END DO
        DO K=1,NFEN
          DO NH=1,MNP
            READ(IHOT,REC=IHOTSTP) q20(NH,K)
            IHOTSTP=IHOTSTP+1
          END DO
        END DO
        DO K=1,NFEN
          DO NH=1,MNP
            READ(IHOT,REC=IHOTSTP) l(NH,K)
            IHOTSTP=IHOTSTP+1
          END DO
        END DO
        IF (ABS(IDEN).EQ.1) THEN
          DO K=1,NFEN
            DO NH=1,MNP
              READ(IHOT,REC=IHOTSTP) SigT(NH,K)
              IHOTSTP=IHOTSTP+1
            END DO
          END DO
        END IF
        IF(ABS(IDen).EQ.2) THEN
          DO K=1,NFEN
            DO NH=1,MNP
              READ(IHOT,REC=IHOTSTP) Sal(NH,K)
              IHOTSTP=IHOTSTP+1
            END DO
          END DO
        ENDIF
        IF(ABS(IDen).EQ.3) THEN
          DO K=1,NFEN
            DO NH=1,MNP
              READ(IHOT,REC=IHOTSTP) Temp(NH,K)
              IHOTSTP=IHOTSTP+1
            END DO
          END DO
        ENDIF
        IF(ABS(IDen).EQ.4) THEN
          DO K=1,NFEN
            DO NH=1,MNP
              READ(IHOT,REC=IHOTSTP) Sal(NH,K)
              IHOTSTP=IHOTSTP+1
            END DO
          END DO
          DO K=1,NFEN
            DO NH=1,MNP
              READ(IHOT,REC=IHOTSTP) Temp(NH,K)
              IHOTSTP=IHOTSTP+1
            ENDDO
          ENDDO
        END IF

      RETURN
C     ----------------------------------------------------------------------
      END SUBROUTINE ReadHotStart3D
C     ----------------------------------------------------------------------


C     ----------------------------------------------------------------------
C          S U B R O U T I N E     W R I T E   H O T  S T A R T  3 D
C     ----------------------------------------------------------------------
C
C     jgf46.02 This subroutine supports PREP67_68. It writes out the 3D
C     section of the full domain hot start file.
C
C     ----------------------------------------------------------------------
      SUBROUTINE WriteHotStart3D(UnitNumber,FilePosition,IPROC)
C     ----------------------------------------------------------------------
      USE PRE_GLOBAL
      IMPLICIT NONE
      INTEGER, intent(in) :: UnitNumber ! i/o unit of subdomain file
      INTEGER, intent(inout) :: FilePosition ! position in binary file
      INTEGER IHOTSTP, LOCHSF, I, N, IINDX
      INTEGER, intent(in) :: IPROC

C     Start writing out the 3D hotstart information

            LOCHSF=UnitNumber
            IHOTSTP=FilePosition

            WRITE(LOCHSF,REC=IHOTSTP) IDEN ; IHOTSTP = IHOTSTP + 1
            WRITE(LOCHSF,REC=IHOTSTP) N3DSD ; IHOTSTP = IHOTSTP + 1
            WRITE(LOCHSF,REC=IHOTSTP) I3DSDRec ; IHOTSTP = IHOTSTP + 1
            WRITE(LOCHSF,REC=IHOTSTP) N3DSV ; IHOTSTP = IHOTSTP + 1
            WRITE(LOCHSF,REC=IHOTSTP) I3DSVRec ; IHOTSTP = IHOTSTP + 1
            WRITE(LOCHSF,REC=IHOTSTP) N3DST ; IHOTSTP = IHOTSTP + 1
            WRITE(LOCHSF,REC=IHOTSTP) I3DSTRec ; IHOTSTP = IHOTSTP + 1
            WRITE(LOCHSF,REC=IHOTSTP) N3DGD ; IHOTSTP = IHOTSTP + 1
            WRITE(LOCHSF,REC=IHOTSTP) I3DGDRec ; IHOTSTP = IHOTSTP + 1
            WRITE(LOCHSF,REC=IHOTSTP) N3DGV ; IHOTSTP = IHOTSTP + 1
            WRITE(LOCHSF,REC=IHOTSTP) I3DGVRec ; IHOTSTP = IHOTSTP + 1
            WRITE(LOCHSF,REC=IHOTSTP) N3DGT ; IHOTSTP = IHOTSTP + 1
            WRITE(LOCHSF,REC=IHOTSTP) I3DGTRec ; IHOTSTP = IHOTSTP + 1

            DO I=1,NNODP(IPROC)
               IINDX=ABS(IMAP_NOD_LG(I,IPROC))
               WRITE(LOCHSF,REC=IHOTSTP) DUU(IINDX)
               IHOTSTP=IHOTSTP+1
            END DO
            DO I=1,NNODP(IPROC)
               IINDX=ABS(IMAP_NOD_LG(I,IPROC))
               WRITE(LOCHSF,REC=IHOTSTP) DUV(IINDX)
               IHOTSTP=IHOTSTP+1
            END DO
            DO I=1,NNODP(IPROC)
               IINDX=ABS(IMAP_NOD_LG(I,IPROC))
               WRITE(LOCHSF,REC=IHOTSTP) DVV(IINDX)
               IHOTSTP=IHOTSTP+1
            END DO
            DO I=1,NNODP(IPROC)
               IINDX=ABS(IMAP_NOD_LG(I,IPROC))
               WRITE(LOCHSF,REC=IHOTSTP) UU(IINDX)
               IHOTSTP=IHOTSTP+1
            END DO
            DO I=1,NNODP(IPROC)
               IINDX=ABS(IMAP_NOD_LG(I,IPROC))
               WRITE(LOCHSF,REC=IHOTSTP) VV(IINDX)
               IHOTSTP=IHOTSTP+1
            END DO
            DO I=1,NNODP(IPROC)
               IINDX=ABS(IMAP_NOD_LG(I,IPROC))
               WRITE(LOCHSF,REC=IHOTSTP) BSX(IINDX)
               IHOTSTP=IHOTSTP+1
            END DO
            DO I=1,NNODP(IPROC)
               IINDX=ABS(IMAP_NOD_LG(I,IPROC))
               WRITE(LOCHSF,REC=IHOTSTP) BSY(IINDX)
               IHOTSTP=IHOTSTP+1
             ENDDO

            DO N=1,NFEN
              DO I=1,NNODP(IPROC)
                IINDX=ABS(IMAP_NOD_LG(I,IPROC))
                WRITE(LOCHSF,REC=IHOTSTP) RealQ(IINDX,N)
                IHOTSTP=IHOTSTP+1
              END DO
            END DO

            DO N=1,NFEN
              DO I=1,NNODP(IPROC)
                IINDX=ABS(IMAP_NOD_LG(I,IPROC))
                WRITE(LOCHSF,REC=IHOTSTP) ImagQ(IINDX,N)
                IHOTSTP=IHOTSTP+1
              END DO
            END DO

            DO N=1,NFEN
              DO I=1,NNODP(IPROC)
                IINDX=ABS(IMAP_NOD_LG(I,IPROC))
                WRITE(LOCHSF,REC=IHOTSTP) WZ(IINDX,N)
                IHOTSTP=IHOTSTP+1
              END DO
            END DO

            DO N=1,NFEN
              DO I=1,NNODP(IPROC)
                IINDX=ABS(IMAP_NOD_LG(I,IPROC))
                WRITE(LOCHSF,REC=IHOTSTP) q20(IINDX,N)
                IHOTSTP=IHOTSTP+1
              END DO
            END DO

            DO N=1,NFEN
              DO I=1,NNODP(IPROC)
                IINDX=ABS(IMAP_NOD_LG(I,IPROC))
                WRITE(LOCHSF,REC=IHOTSTP) l(IINDX,N)
                IHOTSTP=IHOTSTP+1
              END DO
            END DO

            IF (ABS(IDEN).EQ.1) THEN
              DO N=1,NFEN
                DO I=1,NNODP(IPROC)
                  IINDX=ABS(IMAP_NOD_LG(I,IPROC))
                  WRITE(LOCHSF,REC=IHOTSTP) SigT(IINDX,N)
                  IHOTSTP=IHOTSTP+1
                END DO
              END DO
            ENDIF
            IF(ABS(IDen).EQ.2) THEN
              DO N=1,NFEN
                DO I=1,NNODP(IPROC)
                  IINDX=ABS(IMAP_NOD_LG(I,IPROC))
                  WRITE(LOCHSF,REC=IHOTSTP) Sal(IINDX,N)
                  IHOTSTP=IHOTSTP+1
                END DO
              END DO
            ENDIF
            IF(ABS(IDen).EQ.3) THEN
              DO N=1,NFEN
                DO I=1,NNODP(IPROC)
                  IINDX=ABS(IMAP_NOD_LG(I,IPROC))
                  WRITE(LOCHSF,REC=IHOTSTP) Temp(IINDX,N)
                  IHOTSTP=IHOTSTP+1
                END DO
              END DO
            ENDIF
            IF(ABS(IDen).EQ.4) THEN
              DO N=1,NFEN
                DO I=1,NNODP(IPROC)
                  IINDX=ABS(IMAP_NOD_LG(I,IPROC))
                  WRITE(LOCHSF,REC=IHOTSTP) Sal(IINDX,N)
                  IHOTSTP=IHOTSTP+1
                END DO
              END DO
              DO N=1,NFEN
                DO I=1,NNODP(IPROC)
                  IINDX=ABS(IMAP_NOD_LG(I,IPROC))
                  WRITE(LOCHSF,REC=IHOTSTP) Temp(IINDX,N)
                  IHOTSTP=IHOTSTP+1
                END DO
              END DO
            END IF

      RETURN
C     ----------------------------------------------------------------------
      END SUBROUTINE WriteHotStart3D
C     ----------------------------------------------------------------------

C   kmd48.33bc add in 3D global hot start files
C     ----------------------------------------------------------------------
C          S U B R O U T I N E   R E A D  H O T  S T A R T  3 D  G L O B A L
C     ----------------------------------------------------------------------
C
C     This subroutine supports PREP67_68. It reads in the 3D
C     section of the full domain hot start file.
C
C     ----------------------------------------------------------------------
      SUBROUTINE ReadHotStart3DGlobal(UnitNumber,FilePosition,ProcessNo)
C     ----------------------------------------------------------------------
      USE PRE_GLOBAL
      IMPLICIT NONE
      INTEGER, intent(in) :: UnitNumber ! i/o unit of full domain file
      INTEGER, intent(inout) :: FilePosition ! position in binary file
      INTEGER, intent(in) :: ProcessNo ! i/o unit of full domain file
      INTEGER IHOTSTP,NH,K,I,IPROC
      INTEGER IINDX
      REAL(8) RVALUE
C
C     Start reading in the data

!kmd Added information for 3D hotstart
      ALLOCATE ( DUU(MNP),DUV(MNP),DVV(MNP))
      ALLOCATE ( UU(MNP),VV(MNP))
      ALLOCATE ( BSX(MNP),BSY(MNP))
!kmd end of additions

        IHOT=UnitNumber
        IHOTSTP=FilePosition
        IPROC=ProcessNo

        PRINT*, "How many layers need to be evaluated:"
        READ *, NFEN

C        PRINT *, "Made it to the 3d portion of this"
        PRINT *, "NFEN = ", NFEN

!kmd Added information for 3D hotstart
      ALLOCATE ( WZ(MNP,NFEN), q20(MNP,NFEN))
      ALLOCATE ( RealQ(MNP,NFEN), ImagQ(MNP,NFEN))
      ALLOCATE ( l(MNP,NFEN), SigT(MNP,NFEN))
      ALLOCATE ( Sal(MNP,NFEN), Temp(MNP,NFEN))
!kmd end of additions

        READ(IHOT,REC=IHOTSTP) IDEN ; IHOTSTP = IHOTSTP + 1
        READ(IHOT,REC=IHOTSTP) N3DSD ; IHOTSTP = IHOTSTP + 1
        READ(IHOT,REC=IHOTSTP) I3DSDRec ; IHOTSTP = IHOTSTP + 1
        READ(IHOT,REC=IHOTSTP) N3DSV ; IHOTSTP = IHOTSTP + 1
        READ(IHOT,REC=IHOTSTP) I3DSVRec ; IHOTSTP = IHOTSTP + 1
        READ(IHOT,REC=IHOTSTP) N3DST ; IHOTSTP = IHOTSTP + 1
        READ(IHOT,REC=IHOTSTP) I3DSTRec ; IHOTSTP = IHOTSTP + 1
        READ(IHOT,REC=IHOTSTP) N3DGD ; IHOTSTP = IHOTSTP + 1
        READ(IHOT,REC=IHOTSTP) I3DGDRec ; IHOTSTP = IHOTSTP + 1
        READ(IHOT,REC=IHOTSTP) N3DGV ; IHOTSTP = IHOTSTP + 1
        READ(IHOT,REC=IHOTSTP) I3DGVRec ; IHOTSTP = IHOTSTP + 1
        READ(IHOT,REC=IHOTSTP) N3DGT ; IHOTSTP = IHOTSTP + 1
        READ(IHOT,REC=IHOTSTP) I3DGTRec ; IHOTSTP = IHOTSTP + 1

        DO I=1,NNODP(IPROC)
          IINDX=IMAP_NOD_LG(I,IPROC)
          READ(IHOT,REC=IHOTSTP) RVALUE
          IF (IINDX.GT.0) DUU(IINDX) = RVALUE
          IHOTSTP=IHOTSTP+1
        END DO
        DO I=1,NNODP(IPROC)
          IINDX=IMAP_NOD_LG(I,IPROC)
          READ(IHOT,REC=IHOTSTP) RVALUE
          IF (IINDX.GT.0) DUV(IINDX) = RVALUE
          IHOTSTP=IHOTSTP+1
        END DO
        DO I=1,NNODP(IPROC)
          IINDX=IMAP_NOD_LG(I,IPROC)
          READ(IHOT,REC=IHOTSTP) RVALUE
          IF (IINDX.GT.0) DVV(IINDX) = RVALUE
          IHOTSTP=IHOTSTP+1
        END DO
        DO I=1,NNODP(IPROC)
          IINDX=IMAP_NOD_LG(I,IPROC)
          READ(IHOT,REC=IHOTSTP) RVALUE
          IF (IINDX.GT.0) UU(IINDX) = RVALUE
          IHOTSTP=IHOTSTP+1
        END DO
        DO I=1,NNODP(IPROC)
          IINDX=IMAP_NOD_LG(I,IPROC)
          READ(IHOT,REC=IHOTSTP) RVALUE
          IF (IINDX.GT.0) VV(IINDX) = RVALUE
          IHOTSTP=IHOTSTP+1
        END DO
        DO I=1,NNODP(IPROC)
          IINDX=IMAP_NOD_LG(I,IPROC)
          READ(IHOT,REC=IHOTSTP) RVALUE
          IF (IINDX.GT.0) BSX(IINDX) = RVALUE
          IHOTSTP=IHOTSTP+1
        END DO
        DO I=1,NNODP(IPROC)
          IINDX=IMAP_NOD_LG(I,IPROC)
          READ(IHOT,REC=IHOTSTP) RVALUE
          IF (IINDX.GT.0) BSY(IINDX) = RVALUE
          IHOTSTP=IHOTSTP+1
        END DO

        DO K=1,NFEN
          DO I=1,NNODP(IPROC)
            IINDX=IMAP_NOD_LG(I,IPROC)
            READ(IHOT,REC=IHOTSTP) RVALUE
            IF (IINDX.GT.0) REALQ(IINDX,K) = RVALUE
            IHOTSTP=IHOTSTP+1
          END DO
        END DO
        DO K=1,NFEN
          DO I=1,NNODP(IPROC)
            IINDX=IMAP_NOD_LG(I,IPROC)
            READ(IHOT,REC=IHOTSTP) RVALUE
            IF (IINDX.GT.0) ImagQ(IINDX,K) = RVALUE
            IHOTSTP=IHOTSTP+1
          END DO
        END DO
        DO K=1,NFEN
          DO I=1,NNODP(IPROC)
            IINDX=IMAP_NOD_LG(I,IPROC)
            READ(IHOT,REC=IHOTSTP) RVALUE
            IF (IINDX.GT.0) WZ(IINDX,K) = RVALUE
            IHOTSTP=IHOTSTP+1
          END DO
        END DO
        DO K=1,NFEN
          DO I=1,NNODP(IPROC)
            IINDX=IMAP_NOD_LG(I,IPROC)
            READ(IHOT,REC=IHOTSTP) RVALUE
            IF (IINDX.GT.0) q20(IINDX,K) = RVALUE
            IHOTSTP=IHOTSTP+1
          END DO
        END DO
        DO K=1,NFEN
          DO I=1,NNODP(IPROC)
            IINDX=IMAP_NOD_LG(I,IPROC)
            READ(IHOT,REC=IHOTSTP) RVALUE
            IF (IINDX.GT.0) l(IINDX,K) = RVALUE
            IHOTSTP=IHOTSTP+1
          END DO
        END DO
        IF (ABS(IDEN).EQ.1) THEN
          DO K=1,NFEN
            DO I=1,NNODP(IPROC)
              IINDX=IMAP_NOD_LG(I,IPROC)
              READ(IHOT,REC=IHOTSTP) RVALUE
              IF (IINDX.GT.0) SigT(IINDX,K) = RVALUE
              IHOTSTP=IHOTSTP+1
            END DO
          END DO
        END IF
        IF(ABS(IDen).EQ.2) THEN
          DO K=1,NFEN
            DO I=1,NNODP(IPROC)
              IINDX=IMAP_NOD_LG(I,IPROC)
              READ(IHOT,REC=IHOTSTP) RVALUE
              IF (IINDX.GT.0) Sal(IINDX,K) = RVALUE
              IHOTSTP=IHOTSTP+1
            END DO
          END DO
        ENDIF
        IF(ABS(IDen).EQ.3) THEN
          DO K=1,NFEN
            DO I=1,NNODP(IPROC)
              IINDX=IMAP_NOD_LG(I,IPROC)
              READ(IHOT,REC=IHOTSTP) RVALUE
              IF (IINDX.GT.0) Temp(IINDX,K) = RVALUE
              IHOTSTP=IHOTSTP+1
            END DO
          END DO
        ENDIF
        IF(ABS(IDen).EQ.4) THEN
          DO K=1,NFEN
            DO I=1,NNODP(IPROC)
              IINDX=IMAP_NOD_LG(I,IPROC)
              READ(IHOT,REC=IHOTSTP) RVALUE
              IF (IINDX.GT.0) Sal(IINDX,K) = RVALUE
              IHOTSTP=IHOTSTP+1
            END DO
          END DO
          DO K=1,NFEN
            DO I=1,NNODP(IPROC)
              IINDX=IMAP_NOD_LG(I,IPROC)
              READ(IHOT,REC=IHOTSTP) RVALUE
              IF (IINDX.GT.0) TEMP(IINDX,K) = RVALUE
              IHOTSTP=IHOTSTP+1
            END DO
          ENDDO
        END IF

      RETURN
C     ----------------------------------------------------------------------
      END SUBROUTINE ReadHotStart3DGlobal
C     ----------------------------------------------------------------------

C   kmd48.33bc added for global 3D hot start
C     ----------------------------------------------------------------------
C       S U B R O U T I N E   W R I T E   H O T  S T A R T  3 D  G L O B A L
C     ----------------------------------------------------------------------
C
C     This subroutine supports PREP67_68. It writes out the 3D
C     section of the full domain hot start file.
C
C     ----------------------------------------------------------------------
      SUBROUTINE WriteHotStart3DGlobal(UnitNumber,FilePosition,IPROC)
C     ----------------------------------------------------------------------
      USE PRE_GLOBAL
      IMPLICIT NONE
      INTEGER, intent(in) :: UnitNumber ! i/o unit of subdomain file
      INTEGER, intent(inout) :: FilePosition ! position in binary file
      INTEGER, intent(in) :: IPROC
      INTEGER IHOTSTP, LOCHSF, I, N, IINDX

C     Start writing out the 3D hotstart information

            LOCHSF=UnitNumber
            IHOTSTP=FilePosition

            WRITE(LOCHSF,REC=IHOTSTP) IDEN ; IHOTSTP = IHOTSTP + 1
            WRITE(LOCHSF,REC=IHOTSTP) N3DSD ; IHOTSTP = IHOTSTP + 1
            WRITE(LOCHSF,REC=IHOTSTP) I3DSDRec ; IHOTSTP = IHOTSTP + 1
            WRITE(LOCHSF,REC=IHOTSTP) N3DSV ; IHOTSTP = IHOTSTP + 1
            WRITE(LOCHSF,REC=IHOTSTP) I3DSVRec ; IHOTSTP = IHOTSTP + 1
            WRITE(LOCHSF,REC=IHOTSTP) N3DST ; IHOTSTP = IHOTSTP + 1
            WRITE(LOCHSF,REC=IHOTSTP) I3DSTRec ; IHOTSTP = IHOTSTP + 1
            WRITE(LOCHSF,REC=IHOTSTP) N3DGD ; IHOTSTP = IHOTSTP + 1
            WRITE(LOCHSF,REC=IHOTSTP) I3DGDRec ; IHOTSTP = IHOTSTP + 1
            WRITE(LOCHSF,REC=IHOTSTP) N3DGV ; IHOTSTP = IHOTSTP + 1
            WRITE(LOCHSF,REC=IHOTSTP) I3DGVRec ; IHOTSTP = IHOTSTP + 1
            WRITE(LOCHSF,REC=IHOTSTP) N3DGT ; IHOTSTP = IHOTSTP + 1
            WRITE(LOCHSF,REC=IHOTSTP) I3DGTRec ; IHOTSTP = IHOTSTP + 1

            DO I=1,MNP
               IINDX=ABS(IMAP_NOD_LG(I,IPROC))
               WRITE(LOCHSF,REC=IHOTSTP) DUU(IINDX)
               IHOTSTP=IHOTSTP+1
            END DO
            DO I=1,NNODP(IPROC)
               IINDX=ABS(IMAP_NOD_LG(I,IPROC))
               WRITE(LOCHSF,REC=IHOTSTP) DUV(IINDX)
               IHOTSTP=IHOTSTP+1
            END DO
            DO I=1,NNODP(IPROC)
               IINDX=ABS(IMAP_NOD_LG(I,IPROC))
               WRITE(LOCHSF,REC=IHOTSTP) DVV(IINDX)
               IHOTSTP=IHOTSTP+1
            END DO
            DO I=1,NNODP(IPROC)
               IINDX=ABS(IMAP_NOD_LG(I,IPROC))
               WRITE(LOCHSF,REC=IHOTSTP) UU(IINDX)
               IHOTSTP=IHOTSTP+1
            END DO
            DO I=1,NNODP(IPROC)
               IINDX=ABS(IMAP_NOD_LG(I,IPROC))
               WRITE(LOCHSF,REC=IHOTSTP) VV(IINDX)
               IHOTSTP=IHOTSTP+1
            END DO
            DO I=1,NNODP(IPROC)
               IINDX=ABS(IMAP_NOD_LG(I,IPROC))
               WRITE(LOCHSF,REC=IHOTSTP) BSX(IINDX)
               IHOTSTP=IHOTSTP+1
            END DO
            DO I=1,NNODP(IPROC)
               IINDX=ABS(IMAP_NOD_LG(I,IPROC))
               WRITE(LOCHSF,REC=IHOTSTP) BSY(IINDX)
               IHOTSTP=IHOTSTP+1
             ENDDO

            DO N=1,NFEN
              DO I=1,NNODP(IPROC)
                IINDX=ABS(IMAP_NOD_LG(I,IPROC))
                WRITE(LOCHSF,REC=IHOTSTP) RealQ(IINDX,N)
                IHOTSTP=IHOTSTP+1
              END DO
            END DO
            DO N=1,NFEN
              DO I=1,NNODP(IPROC)
                IINDX=ABS(IMAP_NOD_LG(I,IPROC))
                WRITE(LOCHSF,REC=IHOTSTP) ImagQ(IINDX,N)
                IHOTSTP=IHOTSTP+1
              END DO
            END DO
            DO N=1,NFEN
              DO I=1,NNODP(IPROC)
                IINDX=ABS(IMAP_NOD_LG(I,IPROC))
                WRITE(LOCHSF,REC=IHOTSTP) WZ(IINDX,N)
                IHOTSTP=IHOTSTP+1
              END DO
            END DO
            DO N=1,NFEN
              DO I=1,NNODP(IPROC)
                IINDX=ABS(IMAP_NOD_LG(I,IPROC))
                WRITE(LOCHSF,REC=IHOTSTP) q20(IINDX,N)
                IHOTSTP=IHOTSTP+1
              END DO
            END DO
            DO N=1,NFEN
              DO I=1,NNODP(IPROC)
                IINDX=ABS(IMAP_NOD_LG(I,IPROC))
                WRITE(LOCHSF,REC=IHOTSTP) l(IINDX,N)
                IHOTSTP=IHOTSTP+1
              END DO
            END DO
            IF (ABS(IDEN).EQ.1) THEN
              DO N=1,NFEN
                DO I=1,NNODP(IPROC)
                  IINDX=ABS(IMAP_NOD_LG(I,IPROC))
                  WRITE(LOCHSF,REC=IHOTSTP) SigT(IINDX,N)
                  IHOTSTP=IHOTSTP+1
                END DO
              END DO
            ENDIF
            IF(ABS(IDen).EQ.2) THEN
              DO N=1,NFEN
                DO I=1,NNODP(IPROC)
                  IINDX=ABS(IMAP_NOD_LG(I,IPROC))
                  WRITE(LOCHSF,REC=IHOTSTP) Sal(IINDX,N)
                  IHOTSTP=IHOTSTP+1
                END DO
              END DO
            ENDIF
            IF(ABS(IDen).EQ.3) THEN
              DO N=1,NFEN
                DO I=1,NNODP(IPROC)
                  IINDX=ABS(IMAP_NOD_LG(I,IPROC))
                  WRITE(LOCHSF,REC=IHOTSTP) Temp(IINDX,N)
                  IHOTSTP=IHOTSTP+1
                END DO
              END DO
            ENDIF
            IF(ABS(IDen).EQ.4) THEN
              DO N=1,NFEN
                DO I=1,NNODP(IPROC)
                  IINDX=ABS(IMAP_NOD_LG(I,IPROC))
                  WRITE(LOCHSF,REC=IHOTSTP) Sal(IINDX,N)
                  IHOTSTP=IHOTSTP+1
                END DO
              END DO
              DO N=1,NFEN
                DO I=1,NNODP(IPROC)
                  IINDX=ABS(IMAP_NOD_LG(I,IPROC))
                  WRITE(LOCHSF,REC=IHOTSTP) Temp(IINDX,N)
                  IHOTSTP=IHOTSTP+1
                END DO
              END DO
            END IF

C
      RETURN
C     ----------------------------------------------------------------------
      END SUBROUTINE WriteHotStart3DGLOBAL
C     ----------------------------------------------------------------------

C   kmd48.33bc add read information for the initial condition file
C     ----------------------------------------------------------------------
C          S U B R O U T I N E     R E A D  I N I T C O N D  3 D
C     ----------------------------------------------------------------------
C
C     kmd47.22 reads in the 3D information from an initial condition
C     file
C
C     ----------------------------------------------------------------------
      SUBROUTINE ReadInitCond3D(UnitNumber)
C     ----------------------------------------------------------------------
      USE PRE_GLOBAL
      IMPLICIT NONE
      INTEGER, intent(in) :: UnitNumber ! i/o unit of full domain file
      INTEGER NH,K
C
C     Start reading in the data

!kmd Added information for 3D hotstart
      ALLOCATE ( DUU(MNP),DUV(MNP),DVV(MNP))
      ALLOCATE ( UU(MNP),VV(MNP))
      ALLOCATE ( BSX(MNP),BSY(MNP))
      ALLOCATE ( WZ(MNP,NFEN), q20(MNP,NFEN))
      ALLOCATE ( RealQ(MNP,NFEN), ImagQ(MNP,NFEN))
      ALLOCATE ( l(MNP,NFEN), SigT(MNP,NFEN))
      ALLOCATE ( Sal(MNP,NFEN), Temp(MNP,NFEN))
!kmd end of additions

        IHOT=UnitNumber

        READ(IHOT,*) IDEN

        DO NH=1,MNP
          READ(IHOT,*) BSX(NH)
        END DO

        DO NH=1,MNP
          READ(IHOT,*) BSY(NH)
        END DO

        DO K=1,NFEN
          DO NH=1,MNP
            READ(IHOT,*) RealQ(NH,K)
          END DO
        END DO
        DO K=1,NFEN
          DO NH=1,MNP
            READ(IHOT,*) ImagQ(NH,K)
          END DO
        END DO
        DO K=1,NFEN
          DO NH=1,MNP
            READ(IHOT,*) WZ(NH,K)
          END DO
        END DO
        DO K=1,NFEN
          DO NH=1,MNP
            READ(IHOT,*) q20(NH,K)
          END DO
        END DO
        DO K=1,NFEN
          DO NH=1,MNP
            READ(IHOT,*) l(NH,K)
          END DO
        END DO
        IF (ABS(IDEN).EQ.1) THEN
          DO K=1,NFEN
            DO NH=1,MNP
              READ(IHOT,*) SigT(NH,K)
            END DO
          END DO
        END IF
        IF(ABS(IDen).EQ.2) THEN
          DO K=1,NFEN
            DO NH=1,MNP
              READ(IHOT,*) Sal(NH,K)
            END DO
          END DO
        ENDIF
        IF(ABS(IDen).EQ.3) THEN
          DO K=1,NFEN
            DO NH=1,MNP
              READ(IHOT,*) Temp(NH,K)
            END DO
          END DO
        ENDIF
        IF(ABS(IDen).EQ.4) THEN
          DO K=1,NFEN
            DO NH=1,MNP
              READ(IHOT,*) Sal(NH,K)
            END DO
          END DO
          DO K=1,NFEN
            DO NH=1,MNP
              READ(IHOT,*) Temp(NH,K)
            ENDDO
          ENDDO
        END IF
      RETURN

C     ----------------------------------------------------------------------
      END SUBROUTINE ReadInitCond3D
C     ----------------------------------------------------------------------


C     ----------------------------------------------------------------------
C          S U B R O U T I N E     W R I T E   I N I T C O N D  3 D
C     ----------------------------------------------------------------------
C
C     kmd47.22 This subroutine writes out the 3D section of the
C     full domain initial condtion file.
C
C     ----------------------------------------------------------------------
      SUBROUTINE WriteInitCond3D(UnitNumber,IPROC)
C     ----------------------------------------------------------------------
      USE PRE_GLOBAL
      IMPLICIT NONE
      INTEGER, intent(in) :: UnitNumber ! i/o unit of subdomain file
      INTEGER LOCHSF, I, N, IINDX
      INTEGER, intent(in) :: IPROC

C     Start writing out the 3D initial condition information

            LOCHSF=UnitNumber

            WRITE(LOCHSF,*) IDEN

            DO I=1,NNODP(IPROC)
              IINDX=ABS(IMAP_NOD_LG(I,IPROC))
              WRITE(LOCHSF,*) BSX(IINDX)
            END DO

            DO I=1,NNODP(IPROC)
              IINDX=ABS(IMAP_NOD_LG(I,IPROC))
              WRITE(LOCHSF,*) BSY(IINDX)
            END DO

            DO N=1,NFEN
              DO I=1,NNODP(IPROC)
                IINDX=ABS(IMAP_NOD_LG(I,IPROC))
                WRITE(LOCHSF,*) RealQ(IINDX,N)
              END DO
            END DO
            DO N=1,NFEN
              DO I=1,NNODP(IPROC)
                IINDX=ABS(IMAP_NOD_LG(I,IPROC))
                WRITE(LOCHSF,*) ImagQ(IINDX,N)
              END DO
            END DO
            DO N=1,NFEN
              DO I=1,NNODP(IPROC)
                IINDX=ABS(IMAP_NOD_LG(I,IPROC))
                WRITE(LOCHSF,*) WZ(IINDX,N)
              END DO
            END DO
            DO N=1,NFEN
              DO I=1,NNODP(IPROC)
                IINDX=ABS(IMAP_NOD_LG(I,IPROC))
                WRITE(LOCHSF,*) q20(IINDX,N)
              END DO
            END DO
            DO N=1,NFEN
              DO I=1,NNODP(IPROC)
                IINDX=ABS(IMAP_NOD_LG(I,IPROC))
                WRITE(LOCHSF,*) l(IINDX,N)
              END DO
            END DO
            IF (ABS(IDEN).EQ.1) THEN
              DO N=1,NFEN
                DO I=1,NNODP(IPROC)
                  IINDX=ABS(IMAP_NOD_LG(I,IPROC))
                  WRITE(LOCHSF,*) SigT(IINDX,N)
                END DO
              END DO
            ENDIF
            IF(ABS(IDen).EQ.2) THEN
              DO N=1,NFEN
                DO I=1,NNODP(IPROC)
                  IINDX=ABS(IMAP_NOD_LG(I,IPROC))
                  WRITE(LOCHSF,*) Sal(IINDX,N)
                END DO
              END DO
            ENDIF
            IF(ABS(IDen).EQ.3) THEN
              DO N=1,NFEN
                DO I=1,NNODP(IPROC)
                  IINDX=ABS(IMAP_NOD_LG(I,IPROC))
                  WRITE(LOCHSF,*) Temp(IINDX,N)
                END DO
              END DO
            ENDIF
            IF(ABS(IDen).EQ.4) THEN
              DO N=1,NFEN
                DO I=1,NNODP(IPROC)
                  IINDX=ABS(IMAP_NOD_LG(I,IPROC))
                  WRITE(LOCHSF,*) Sal(IINDX,N)
                END DO
              END DO
              DO N=1,NFEN
                DO I=1,NNODP(IPROC)
                  IINDX=ABS(IMAP_NOD_LG(I,IPROC))
                  WRITE(LOCHSF,*) Temp(IINDX,N)
                END DO
              END DO
            END IF
C
      RETURN
C     ----------------------------------------------------------------------
      END SUBROUTINE WriteInitCond3D
C     ----------------------------------------------------------------------


C---------------------------------------------------------------------------
C     S U B R O U T I N E    O P E N  F U L L  D O M A I N  F I L E
C---------------------------------------------------------------------------
C
C     jgf47.02 This subroutine will open the full domain file
C
C---------------------------------------------------------------------------
      SUBROUTINE OpenFullDomainFile(UnitNumber, Description, Success)
C---------------------------------------------------------------------------
      USE PRE_GLOBAL
      IMPLICIT NONE
      INTEGER, intent(in) :: UnitNumber     ! i/o unit number to open
      CHARACTER(len=30), intent(in) :: Description ! description of file
      LOGICAL, intent(out):: Success     ! .true. if file opened w/o errors
      LOGICAL Found               !.true. if the full domain file exists
      CHARACTER(len=80) FileName   ! name of full domain file
      CHARACTER(len=7) DefaultName! default name of full domain file
      INTEGER ErrorIO             ! zero if file opened successfully
      CHARACTER(len=4) skipstring ! indicates user wants to skip this file

      Found = .false.
      Success = .false.
      ErrorIO = 1
      skipstring = 'skip'

      DefaultName(1:5) = 'fort.'
      WRITE(DefaultName(6:7),2) UnitNumber
C
C     Determine the name of the file; if found, open it
 31   IF (USE_DEFAULT) THEN
         FileName = DefaultName
      ELSE
         WRITE(*,850) ! type skip to bypass
         WRITE(*,900) Description
         WRITE(*,910) UnitNumber
         READ(*,'(A)') FileName
         FileName = trim(FILENAME)
      ENDIF
C
C     Determine if full domain file exists
      INQUIRE(FILE=FileName,EXIST=FOUND)
C
C     If it does exist, open it
      IF ( FOUND ) THEN
         WRITE(*,1011) FileName !found
         OPEN(UNIT=UnitNumber, FILE=FileName, IOSTAT=ErrorIO)
         Success = .true.
         IF ( ErrorIO .GT. 0 ) THEN
            WRITE(*,*) "ERROR: Full domain file exists but"
            WRITE(*,*) "cannot be opened."
            Success = .false.
         ELSE
            WRITE(*,*) "Successfully opened full domain file."
         ENDIF
      ELSE
C     Give the user a chance to opt out of prepping this file.
         IF (FileName .eq. skipstring) RETURN ! note the early RETURN
         WRITE(*,1010) FileName !not found
         GOTO 31
      ENDIF

 2    FORMAT(I2)
 30   FORMAT(A30)
 850  FORMAT(/,'Type ''skip'' to bypass preprocessing or')
 900  FORMAT('Enter the name of the ',A30)
 910  FORMAT('file (unit ',I3,'): ')
 1010 FORMAT('File ',A7,/,' WAS NOT FOUND! Try again or type "skip"',/)
 1011 FORMAT('File ',A7,/,' WAS FOUND!  Opening & Processing file.',/)
      RETURN
C---------------------------------------------------------------------------
      END SUBROUTINE OpenFullDomainFile
C---------------------------------------------------------------------------



C---------------------------------------------------------------------------
C     S U B R O U T I N E   O P E N  S U B D O M A I N  F I L E
C---------------------------------------------------------------------------
C
C     jgf47.02 This subroutine will open a single subdomain file
C
C---------------------------------------------------------------------------
      SUBROUTINE OpenSubDomainFile(UnitNumber, IProc, sdu, Success)
C---------------------------------------------------------------------------
      USE PRE_GLOBAL
      IMPLICIT NONE
      INTEGER, intent(in) :: UnitNumber     ! i/o unit number of full dom file
      INTEGER, intent(in) :: iproc          ! subdomain number
      INTEGER, intent(out) :: sdu  ! i/o unit nunber that was opened
      LOGICAL, intent(out):: Success     ! .true. if files opened w/o errors
      LOGICAL Found               !.true. if the full domain file exists
      CHARACTER(len=80) FileName   ! name of full domain file
      CHARACTER(len=7) DefaultName! default name of full domain file
      INTEGER ErrorIO             ! zero if file opened successfully
      CHARACTER*14 sdFileName     ! subdomain file name
      CHARACTER(len=4) skipstring ! indicates user wants to skip this file

      Found = .false.
      Success = .false.
      ErrorIO = 1
      skipstring = 'skip'

      DefaultName(1:5) = 'fort.'
      WRITE(DefaultName(6:7),2) UnitNumber

C     Open subdomain file
      sdu = 105 + (iproc-1)
      sdFileName(1:7) = 'PE0000/'
      sdFileName(8:14) = DefaultName
      CALL IWRITE(sdFileName, 3, 6, iproc-1)
      OPEN (UNIT=sdu, FILE=sdFileName, IOSTAT=ErrorIO)
      Success = .true.
      IF ( ErrorIO .GT. 0 ) THEN
         WRITE(*,*) "ERROR: Subdomain file cannot be opened."
         print *,sdu
         print *,sdFileName
         Success = .false.
      ENDIF
 2    FORMAT(I2)
      RETURN

C---------------------------------------------------------------------------
      END SUBROUTINE OpenSubDomainFile
C---------------------------------------------------------------------------



C---------------------------------------------------------------------------
C           S U B R O U T I N E   O P E N  P R E P  F I L E S
C---------------------------------------------------------------------------
C
C     jgf45.12 This subroutine will open the full domain file and each
C     of the subdomain files in the range of subdomains provided in the
C     arguments. It assumes that all unit numbers are between 10 and 99.
C     
C     tcm v50.66.03 -- incrased unit numbers to include 10 - 999.
C---------------------------------------------------------------------------
      SUBROUTINE OpenPrepFiles(UnitNumber, Description,
     &     startProc, endProc, SDU, Success)
C---------------------------------------------------------------------------
      USE PRE_GLOBAL
      IMPLICIT NONE
      INTEGER, intent(in) :: UnitNumber     ! i/o unit number to open
      CHARACTER(*), intent(in) :: Description ! description of file
      INTEGER, intent(in) :: startProc        ! subdomains to start with
      INTEGER, intent(in) :: endProc          ! subdomain to end on
      INTEGER, intent(out), dimension(nproc) :: SDU ! Subdomain unit numbers
      LOGICAL, intent(out):: Success     ! .true. if files opened w/o errors
      LOGICAL Found               !.true. if the full domain file exists
      CHARACTER(len=80) FileName   ! name of full domain file
      CHARACTER(len=8) DefaultName! default name of full domain file  !increased from 7 to 8 tcm v50.66.03
      INTEGER :: dnlen  !length of defaultname (7 or 8)
      INTEGER ErrorIO             ! zero if file opened successfully
      INTEGER iproc               ! subdomain index
      CHARACTER(len=15) sdFileName     ! subdomain file name   !increased from 14 to 15 tcm v50.66.03
      CHARACTER(len=4) skipstring ! indicates user wants to skip this file

      Found = .false.
      Success = .false.
      ErrorIO = 1
      skipstring = 'skip'
      DefaultName(:) = ' '  !initialize to all blanks !tcm v50.77
  
      DefaultName(1:5) = 'fort.'

c.....tcm v50.66.03 increased unit number to 100's places
      if (UnitNumber < 100) then
         WRITE(DefaultName(6:7),2) UnitNumber
         dnlen = 7
      else
         WRITE(DefaultName(6:8),3) UnitNumber
         dnlen = 8
      endif

C
C     Determine the name of the file; if found, open it
 31   IF (USE_DEFAULT) THEN
         FileName = trim(DefaultName(1:dnlen))   !tcm v50.66.03 added trim !tcm v50.77 added 1:dnlen
Casey 120402: Avoid an endless loop.  If the default file does not exist,
C             then give the user a chance to specify the file name or skip.a
C     ELSE
         GOTO 33
      ENDIF
 32   CONTINUE
         WRITE(*,850) ! type skip to bypass
         WRITE(*,900) Description
         WRITE(*,910) UnitNumber
         READ(*,'(A)') FileName
         FileName = trim(FILENAME)
C     ENDIF
 33   CONTINUE

C
C     Determine if full domain file exists
      INQUIRE(FILE=FileName,EXIST=FOUND)
C
C     If it does exist, open it
      IF ( FOUND ) THEN
         WRITE(*,1011) FileName !found
         OPEN(UNIT=UnitNumber, FILE=FileName, IOSTAT=ErrorIO)
         Success = .true.
         IF ( ErrorIO .GT. 0 ) THEN
            WRITE(*,*) "ERROR: Full domain file exists but"
            WRITE(*,*) "cannot be opened."
            Success = .false.
         ENDIF
      ELSE
C     Give the user a chance to opt out of prepping this file.
         IF (FileName .eq. skipstring) RETURN ! note the early RETURN
         WRITE(*,1010) FileName !not found
Casey 120402: Avoid an endless loop.  If the default file does not exist,
C             then give the user a chance to specify the file name or skip.
C        GOTO 31
         GOTO 32
      ENDIF
C
      If (.not.Success) RETURN ! failed to open full domain file
C
C     Open each of the subdomain files
      DO iproc = startProc, endProc
         sdu(iproc) = 105 + (iproc-1)
         sdFileName(1:7) = 'PE0000/'
c........tcm v50.66.03 increased unit number to 100's places
         if (UnitNumber < 100 ) then
            sdFileName(8:14) = DefaultName(1:7)
         else
            sdFileName(8:15) = DefaultName
         endif
#ifdef ADCSWAN
         sdFileName = 'PE0000/'//FileName
#endif
         CALL IWRITE(sdFileName, 3, 6, iproc-1)
         OPEN (UNIT=SDU(iproc), FILE=sdFileName, IOSTAT=ErrorIO)
         Success = .true.
         IF ( ErrorIO .GT. 0 ) THEN
            WRITE(*,*) "ERROR: Subdomain file cannot be opened."
            Success = .false.
            RETURN ! failed to open at least one subdomain file
         ENDIF
      ENDDO

 2    FORMAT(I2)
 3    FORMAT(I3)
 30   FORMAT(A30)
 850  FORMAT(/,'Type ''skip'' to bypass preprocessing or')
 900  FORMAT('Enter the name of the ',A30)
 910  FORMAT('file (unit ',I3,'): ')
 1010 FORMAT('File ',A8,/,' WAS NOT FOUND! Try again or type "skip"',/)  !increased A7 to A8  tcm v50.66.03
 1011 FORMAT('File ',A8,/,' WAS FOUND!  Opening & Processing file.',/)   !increased A7 to A8
      RETURN
C---------------------------------------------------------------------------
      END SUBROUTINE OpenPrepFiles
C---------------------------------------------------------------------------



c***********************************************************************
c   Subroutine to write out to the hotstart file (UNITS 67 and 68)     *
c   header information and the LHS matrix for the harmonic analysis    *
c                                                                      *
c                        R.L.  11/8/95                                 *
c***********************************************************************
c
c     SUBROUTINE HAHOUT(NP,NSTAE,NSTAV,ISTAE,ISTAV,IGLOE,IGLOV,
c    &  IOUNIT,IHOTSTP)
c     implicit none
c     INTEGER NP,NSTAE,NSTAV,ISTAE,AE,ISTAV
c     INTEGER IGLOE,IGLOV,IOUNIT,IHOTSTP,I,J
c     CHARACTER*16 FNAME
c     CHARACTER*8 FNAM8(2)
c     EQUIVALENCE (FNAM8(1),FNAME)c
c
c
c***** Write Out various parameter values
c
c     WRITE(IOUNIT,REC=IHOTSTP+1) NZ
c     WRITE(IOUNIT,REC=IHOTSTP+2) NF
c     WRITE(IOUNIT,REC=IHOTSTP+3) MM
c     WRITE(IOUNIT,REC=IHOTSTP+4) NP
c     WRITE(IOUNIT,REC=IHOTSTP+5) NSTAE
c     WRITE(IOUNIT,REC=IHOTSTP+6) NSTAV
c     WRITE(IOUNIT,REC=IHOTSTP+7) ISTAE
c     WRITE(IOUNIT,REC=IHOTSTP+8) ISTAV
c     WRITE(IOUNIT,REC=IHOTSTP+9) IGLOE
c     WRITE(IOUNIT,REC=IHOTSTP+10) IGLOV
c     WRITE(IOUNIT,REC=IHOTSTP+11) ICALL
c     WRITE(IOUNIT,REC=IHOTSTP+12) NFREQ
c     IHOTSTP = IHOTSTP+12
c
c     do i=1,nfreq+nf
c        FNAME=NAMEFR(I)
c        WRITE(IOUNIT,REC=IHOTSTP+1) FNAM8(1)
c        WRITE(IOUNIT,REC=IHOTSTP+2) FNAM8(2)
c        IHOTSTP=IHOTSTP+2
c        WRITE(IOUNIT,REC=IHOTSTP+1) hafreq(i)
c        WRITE(IOUNIT,REC=IHOTSTP+2) HAFF(i)
c        WRITE(IOUNIT,REC=IHOTSTP+3) HAFACE(i)
c        IHOTSTP=IHOTSTP+3
c     end do
c
c
c***** Write Out time of most recent H.A. update
c
c     WRITE(IOUNIT,REC=IHOTSTP+1) TIMEUD
c     WRITE(IOUNIT,REC=IHOTSTP+2) ITUD
c     IHOTSTP=IHOTSTP+2
c
c***** Write Out LHS Matrix
c
c     do i=1,mm
c        do j=1,mm
c           IHOTSTP = IHOTSTP + 1
c           WRITE(IOUNIT,REC=IHOTSTP) HA(I,J)
c        END DO
c     END DO
c
c     return
c     end subroutine
c
c***********************************************************************
c   Subroutine to write global elevation harmonic analysis RHS load    *
c   vector to a hot start file (UNITS 67 and 68)                       *
c                                                                      *
c                        R.L.  11/8/95                                 *
c***********************************************************************
c
c     SUBROUTINE HAHOUTEG(NP,IOUNIT,IHOTSTP)
c     implicit none
c     INTEGER IOUNIT
c     INTEGER NP,IHOTSTP,N,I
c
c***** Write Out Global Elevation RHS load vector
c
c     do n=1,np
c        do i=1,mm
c           IHOTSTP=IHOTSTP+1
c           WRITE(IOUNIT,REC=IHOTSTP) GLOELV(I,N)
c        end do
c     end do
c
c     return
c     end subroutine

c***********************************************************************
c   Subroutine to write elevation station harmonic analysis RHS load   *
c   vector to a hot start file (UNITS 67 and 68)                       *
c                                                                      *
c                        R.L.  11/8/95                                 *
c***********************************************************************
c
c     SUBROUTINE HAHOUTES(NSTAE,IOUNIT,IHOTSTP)
c     implicit none
c     INTEGER NSTAE,IOUNIT,IHOTSTP,N,I
c
c***** Write Out Station Elevation RHS load vector
c
c     do n=1,NSTAE
c        do i=1,mm
c           IHOTSTP=IHOTSTP+1
c           WRITE(IOUNIT,REC=IHOTSTP) STAELV(I,N)
c        end do
c     end do
c
c     return
c     end subroutine
c
c***********************************************************************
c   Subroutine to write global velocity harmonic analysis RHS load     *
c   vector to a hot start file (UNITS 67 and 68)                       *
c                                                                      *
c                        R.L.  11/8/95                                 *
c***********************************************************************
c
c     SUBROUTINE HAHOUTVG(NP,IOUNIT,IHOTSTP)
c     implicit none
c     INTEGER NP,IOUNIT,IHOTSTP,N,I
c
c***** Write Out Global Velocity RHS load vector
c
c     do n=1,np
c        do i=1,mm
c           IHOTSTP=IHOTSTP+1
c           WRITE(IOUNIT,REC=IHOTSTP) GLOULV(I,N)
c           IHOTSTP=IHOTSTP+1
c           WRITE(IOUNIT,REC=IHOTSTP) GLOVLV(I,N)
c        end do
c     end do
c
c     return
c     end subroutine
c
c***********************************************************************
c   Subroutine to write velocity station harmonic analysis RHS load    *
c   vector to a hot start file (UNITS 67 and 68)                       *
c                                                                      *
c                        R.L.  11/8/95                                 *
c***********************************************************************
c
c     SUBROUTINE HAHOUTVS(NSTAV,IOUNIT,IHOTSTP)
c     implicit none
c     INTEGER NSTAV,IOUNIT,IHOTSTP,N,I
c
c***** Write Out Station Velocity LHS load vector
c
c     do N=1,NSTAV
c        do i=1,mm
c           IHOTSTP=IHOTSTP+1
c           WRITE(IOUNIT,REC=IHOTSTP) STAULV(I,N)
c           IHOTSTP=IHOTSTP+1
c           WRITE(IOUNIT,REC=IHOTSTP) STAVLV(I,N)
c        end do
c     end do
c
c     return
c     end subroutine

C
C---------------------------------------------------------------------------C
C                     (  Serial Version  2/28/98  )                         C
C  This routine writes the domain decomposition information into a file,    C
C  "fort.80".  This file is used by the ADCIRC post-processor ADCPOST.      C
C  This version is compatible with ADCIRC version 34.03                     C
C                                                                           C
C     jgf45.07 Added subdomain->fulldomain element mapping to handle the    C
C     processing of the NOFF array.
C     jgf45.11 Added IDEN for processing fort.44 file, added additional
C     arrays to handle processing of 3D recording stations.
C---------------------------------------------------------------------------C
CC    Addition by CF   8/2007
      SUBROUTINE PREP80()
      USE PRE_GLOBAL
      USE SIZES, ONLY : MNHARF
      USE HARM, ONLY : NHASE, NHASV, NHAGE, NHAGV
      IMPLICIT NONE
      INTEGER I,K
C
      OPEN(UNIT=80,FILE='fort.80')              ! output for ADCPOST
C
C--Write out the domain decomposition information into a file
C  which will later be used in post-processing the results
C
      WRITE(80,80) RUNDES
      WRITE(80,80) RUNID
      WRITE(80,80) AGRID
      WRITE(80,'(2I8,16x,A)') NELG,NNODG,'! Total # elements & nodes'
      WRITE(80,'(I8,24x,A)') NPROC,' ! Number of processors'
      WRITE(80,'(I8,24x,A)') MNPP,'! Max nodes on any processor'
      WRITE(80,'(I8,24x,A)') MNEP,'! Max elements on any processor'!jgf45.07
      WRITE(80,'(I8,24x,A)') IM,'! IM, run type'         !jgf46.02
      WRITE(80,'(I8,24x,A)') NWS,'! NWS, wind data type' !jgf46.02
      WRITE(80,'(I8,24x,A)') NSTAE,'! NSTAE'
      WRITE(80,'(I8,24x,A)') NSTAV,'! NSTAV'
      IF (IM.EQ.10) THEN
         WRITE(80,'(I8,24x,A)') NSTAC,' ! NSTAC' !jgf46.02
      ENDIF
      IF (NWS.NE.0) THEN
         WRITE(80,'(I8,24x,A)') NSTAM,'! NSTAM' !jgf46.02
      ENDIF
      WRITE(80,'(I8,24x,A)') MNHARF,'! MNHARF'
      WRITE(80,'(2I8,16x,A)') MNWLAT,MNWLON,'! NWLON, NWLAT'
C
Casey 100301: Changed I8 to I12.
      DO I = 1,NPROC
         WRITE(80,'(3I8,A33)') I-1, NNODP(I), NOD_RES_TOT(I),
     &        '  ! PE, NNODP(PE), NOD_RES_TOT(PE)'
         WRITE(80,'(9I12)') (IMAP_NOD_LG(K,I),K=1,NNODP(I))
      ENDDO
C
      WRITE(80,*) "GLOBAL   PE     LOCAL   ( Global-to-Local Nodes )"
      DO I = 1,NNODG
         WRITE(80,1140) I, IMAP_NOD_GL(1,I)-1, IMAP_NOD_GL(2,I)
      ENDDO
C
C     jgf45.07 Add subdomain->fulldomain mapping to handle NOFF processing
C     IMAP_EL_LG(I,PE) = Global Element Number of Local Element I on PE
Casey 100301: Changed I8 to I12.
      DO I = 1,NPROC
         WRITE(80,'(2I8,A33)') I-1, NELP(I), '  ! PE, NELP(PE)'
         WRITE(80,'(9I12)') (IMAP_EL_LG(K,I),K=1,NELP(I))
      ENDDO
C
      WRITE(80,'(I8,2E15.8,I8,A32)') NOUTE,TOUTSE,TOUTFE,NSPOOLE,
     &    '   ! NOUTE,TOUTSE,TOUTFE,NSPOOLE'
C
      DO I = 1,NPROC
         WRITE(80,*) I,NSTAEP(I)
         DO K = 1,NSTAEP(I)
            WRITE(80,*) IMAP_STAE_LG(K,I)
         ENDDO
      ENDDO
C
      WRITE(80,'(I8,2E15.8,I8,A32)') NOUTV,TOUTSV,TOUTFV,NSPOOLV,
     &    '   ! NOUTV,TOUTSV,TOUTFV,NSPOOLV'
C
      DO I = 1,NPROC
         WRITE(80,*) I,NSTAVP(I)
         DO K = 1,NSTAVP(I)
            WRITE(80,*) IMAP_STAV_LG(K,I)
         ENDDO
      ENDDO
C
      IF (IM.EQ.10) THEN ! jgf46.02
         WRITE(80,'(I8,2E15.8,I8,A32)') NOUTC,TOUTSC,TOUTFC,NSPOOLC,
     &        '   ! NOUTC,TOUTSC,TOUTFC,NSPOOLC'
         DO I = 1,NPROC
            WRITE(80,*) I,NSTACP(I)
            DO K = 1,NSTACP(I)
               WRITE(80,*) IMAP_STAC_LG(K,I)
            ENDDO
         ENDDO
      ENDIF
C
      IF (NWS.NE.0) THEN ! jgf46.02
         WRITE(80,'(I8,2E15.8,I8,A32)') NOUTM,TOUTSM,TOUTFM,NSPOOLM,
     &        '   ! NOUTM,TOUTSM,TOUTFM,NSPOOLM'
         DO I = 1,NPROC
            WRITE(80,*) I,NSTAMP(I)
            DO K = 1,NSTAMP(I)
               WRITE(80,*) IMAP_STAM_LG(K,I)
            ENDDO
         ENDDO
      ENDIF
C
      WRITE(80,'(I8,2E15.8,I8,A32)') NOUTGE, TOUTSGE,TOUTFGE,NSPOOLGE,
     &    '   ! NOUTGE, TOUTSGE, TOUTFGE, NSPOOLGE'
C
      WRITE(80,'(I8,2E15.8,I8,A32)') NOUTGV, TOUTSGV,TOUTFGV,NSPOOLGV,
     &    '   ! NOUTGV, TOUTSGV, TOUTFGV, NSPOOLGV'
C
      WRITE(80,'(I8,2E15.8,I8,A32)') NOUTGC, TOUTSGC,TOUTFGC,NSPOOLGC,
     &    '   ! NOUTGC, TOUTSGC, TOUTFGC, NSPOOLGC'
C
      WRITE(80,'(I8,2E15.8,I8,A32)') NOUTGW, TOUTSGW,TOUTFGW,NSPOOLGW,
     &    '   ! NOUTGW, TOUTSGW, TOUTFGW, NSPOOLGW'
C
      WRITE(80,'(4I4,A32)') NHASE,NHASV,NHAGE,NHAGV,
     &    '   ! NHASE, NHASV, NHAGE, NHAGV'

C     -------------------------------------------------------------
C
C     S T A R T   3 D   D A T A
C
C     -------------------------------------------------------------

      WRITE(80,*) IDEN !jgf45.11 needed to post process the fort.44 file

C     -------------------------------------------------------------
C     jgf45.11 Write mappings for 3D density stations.
C     -------------------------------------------------------------
      WRITE(80,81) I3DSD, TO3DSDS, TO3DSDF, NSPO3DSD, NSTA3DD,
     &       '   ! I3DSD, TO3DSDS, TO3DSDF, NSPO3DSD, NSTA3DD'
      IF(I3DSD.NE.0) THEN
         DO I = 1, NPROC
            WRITE(80,*) I, NNSTA3DDP(I)
            DO K = 1, NNSTA3DDP(I)
               WRITE(80,*) IMAP_STA3DD_LG(K,I)
            ENDDO
         ENDDO
      ENDIF
C     -------------------------------------------------------------
C     jgf45.11 Write mappings for 3D velocity stations.
C     -------------------------------------------------------------
      WRITE(80,81) I3DSV, TO3DSVS, TO3DSVF, NSPO3DSV, NSTA3DV,
     &       '   ! I3DSV, TO3DSVS, TO3DSVF, NSPO3DSV, NSTA3DV'
      IF(I3DSV.NE.0) THEN
         DO I = 1, NPROC
            WRITE(80,*) I, NNSTA3DVP(I)
            DO K = 1, NNSTA3DVP(I)
               WRITE(80,*) IMAP_STA3DV_LG(K,I)
            ENDDO
         ENDDO
      ENDIF
C     -------------------------------------------------------------
C     jgf45.11 Write mappings for 3D turbulence stations.
C     -------------------------------------------------------------
      WRITE(80,81) I3DST, TO3DSTS, TO3DSTF, NSPO3DST, NSTA3DT,
     &       '   ! I3DST, TO3DSST, TO3DFST, NSPO3DST, NSTA3DT'
      IF(I3DST.NE.0) THEN
         DO I = 1, NPROC
            WRITE(80,*) I, NNSTA3DTP(I)
            DO K = 1, NNSTA3DTP(I)
               WRITE(80,*) IMAP_STA3DT_LG(K,I)
            ENDDO
         ENDDO
      ENDIF
      WRITE(80,82) I3DGD, TO3DGDS, TO3DGDF, NSPO3DGD,
     &       '   ! I3DGD, TO3DGDS, TO3DGDF, NSPO3DGD'
      WRITE(80,82) I3DGV, TO3DGVS, TO3DGVF, NSPO3DGV,
     &       '   ! I3DGV, TO3DGVS, TO3DGVF, NSPO3DGV'
      WRITE(80,82) I3DGT, TO3DGTS, TO3DGTF, NSPO3DGT,
     &       '   ! I3DGT, TO3DGTS, TO3DGTF, NSPO3DGT'
C
C     End 3D data
C
      WRITE(80,*) NBYTE
C
      CLOSE(80)
C
  80  FORMAT(A80)
  81  FORMAT(I8,2E15.8,2I8,A32)
  82  FORMAT(I8,2E15.8,I8,A32)
1130  FORMAT(8X,9I8)
1140  FORMAT(8X,3I8)
C
      RETURN
      END


C---------------------------------------------------------------------------
C                S U B R O U T I N E   P R E P   U N S W A N
C---------------------------------------------------------------------------
C
C     jgf48.17 This subroutine just copies the UnSWAN input file to each
C     subdomain. Currently, it is only called if the compiler directive
C     ADCSWAN is defined.
C
C     jgf50.15: After discussion with Casey, it was decided that changes
C     should be made such that adcprep looks for the swaninit file and
C     parses it to find the name of the swan control file (conventionally
C     called fort.26 in the adcirc community). It should then attempt to
C     copy that file to the subdomains.
C---------------------------------------------------------------------------
      SUBROUTINE  prepUnSWAN
C---------------------------------------------------------------------------
      USE PRE_GLOBAL
      use memory_usage
      IMPLICIT NONE
      INTEGER I,J,IPROC,IPROC2,ILNODE,INDX,NHG,LINDEX
      CHARACTER*170 Line ! line of data from UnSWAN file
      INTEGER SDU(NPROC)  ! subdomain unit numbers
      LOGICAL Success     ! .true. if all files open without error
      LOGICAL swanInitFound !.true. if the swaninit file was found
      LOGICAL origUseDefault ! there is no default UnSWAN input file name
      INTEGER swiLUN ! logical unit number of swaninit file
      INTEGER ErrorIO ! zero if the file was opened successfully
      LOGICAL readError ! .true. if there was an error reading a file
      CHARACTER(36) swanComFile ! name of swan command file
      LOGICAL swanComFound ! .true. if swan command file was found
      CHARACTER(len=43) sdFileName     ! subdomain file name
C
C
Casey 110627: Set the unit number.
      swiLUN = 26

Casey 120402: Changes through this section.  The procedure here is:
C             1. Check for the swaninit file.  If it exists, then read the name
C                of the SWAN control file from the fourth line of swaninit.
C             2. If the swaninit file does not exist, then the user may be using
C                the default control file of INPUT.  Check for the INPUT file.
C                If it exists, then use it as the SWAN control file.
C             3. Otherwise, notify the user of the error.

      swanInitFound = .false.
      swanComFound = .false.
      readError = .true.
      INQUIRE(FILE='swaninit',EXIST=swanInitFound)
      IF (swanInitFound.eqv..true.) THEN
         WRITE(*,*) "INFO: The swaninit file was found."
         OPEN(swiLUN,FILE='swaninit',ACTION='READ',
     &       ACCESS='SEQUENTIAL',IOSTAT=ErrorIO,STATUS='OLD')
         IF (ErrorIO.ne.0) THEN
            WRITE(*,*) "ERROR: The swaninit file could not be opened."
         ELSE
            WRITE(*,*)
     &      "INFO: Parsing swaninit file for swan command file name."
            ! skip down to the 4th line
            DO I=1,4
               READ(swiLUN,*,ERR=4321,END=4321,IOSTAT=ErrorIO) Line
            ENDDO
            READ(Line,*,ERR=4321,END=4321,IOSTAT=ErrorIO) swanComFile
Casey 110627: Close the swaninit file.
            CLOSE(UNIT=swiLUN,STATUS='KEEP')
            WRITE(*,*) "INFO: The swan command file is '",
     &         trim(swanComFile),"'."
            readError = .false.
Casey 120402: Changes to handle the default INPUT control file.
            ! If swan was started without a swaninit file, it will create
            ! its own swaninit, and it will use the string INPUT to
            ! represent the name of the swan command file.
C           IF (TRIM(swanComFile).eq.'INPUT') THEN
C              WRITE(*,*) "ERROR: 'INPUT' is an invalid name."
C           ELSE
               ! check to see if the swan command file is present
               INQUIRE(FILE=trim(swanComFile),EXIST=swanComFound)
               IF (swanComFound.eqv..true.) THEN
                  WRITE(*,*) "INFO: The swan command file '",
     &               trim(swanComFile),"' was found."
                  OPEN(26,FILE=trim(swanComFile),ACTION='READ',
     &            ACCESS='SEQUENTIAL',IOSTAT=ErrorIO,STATUS='OLD')
                  IF (ErrorIO.ne.0) THEN
                     WRITE(*,*) "ERROR: The swan command file '",
     &                  trim(swanComFile),"' could not be opened."
                  ENDIF
               ELSE
                   WRITE(*,*) "ERROR: The swan command file '",
     &               trim(swanComFile),"' was not found."
               ENDIF
C           ENDIF
         ENDIF
Casey 120402: Changes to handle the default INPUT control file.
      ELSE
         INQUIRE(FILE='INPUT',EXIST=swanComFound)
         IF(swanComFound.eqv..true.)THEN
            WRITE(swanComFile,'(A)') "INPUT"
            readError = .false.
            ErrorIO = 0
            swanInitFound = .true.
            swanComFound = .true.
            WRITE(*,*) "INFO: The swan command file is '",
     &         trim(swanComFile),"'."
            ! Casey 120402: Changes to handle the default INPUT control file.
            OPEN(UNIT=26,FILE=TRIM(swanComFile),ACTION='READ')
         ENDIF
      ENDIF
4321  IF ((readError.eqv..true.).or.(ErrorIO.ne.0).or.
     &     (swanInitFound.eqv..false.).or.
     &     (swanComFound.eqv..false.)) THEN
         WRITE(*,*) 'ERROR: There was an error reading swan files.'
         WRITE(*,*) 'WARNING: swan files not preprocessed.'
         RETURN ! note early return
      ENDIF
C
C     Open each of the subdomain files
      DO iproc = 1, nproc
         sdu(iproc) = 105 + (iproc-1)
         sdFileName = 'PE0000/'//swanComFile
         CALL IWRITE(sdFileName, 3, 6, iproc-1)
         OPEN (UNIT=SDU(iproc), FILE=TRIM(sdFileName), IOSTAT=ErrorIO)
         IF ( ErrorIO .GT. 0 ) THEN
            WRITE(*,*) "ERROR: Subdomain file '",TRIM(sdFileName),
     &         " cannot be opened."
            RETURN ! failed to open at least one subdomain file
         ENDIF
      ENDDO
C
Casey 090304: Changed the formatting through the next section.

      DO
         READ(26,'(A)',END=9999) Line
         DO IPROC = 1,NPROC
            WRITE(SDU(IPROC),'(A)') trim(Line)
         ENDDO
      ENDDO
C
C--Close fulldomain file and all the subdomain files
C
 9999 CLOSE (26)
      DO IPROC=1, NPROC
         CLOSE (SDU(IPROC))
      ENDDO

      RETURN
  60  FORMAT(A60)
 170  FORMAT(A170)
 1100 FORMAT(I8,3E13.5)
 1101 FORMAT(' #')
c----------------------------------------------------------------------------
      END SUBROUTINE prepUnSWAN
c----------------------------------------------------------------------------

C---------------------------------------------------------------------------
C                S U B R O U T I N E   P R E P   N E T C D F
C---------------------------------------------------------------------------
C
C     jgf49.31 This subroutine will initialize the fulldomain netcdf
C     output files for the parallel run. These files are initialized
C     in the adcprep phase because they contain fulldomain data which
C     are not available to the subdomains during a parallel run.
C
C---------------------------------------------------------------------------
      SUBROUTINE prepNetCDF()
#ifdef ADCNETCDF
      USE PRESIZES, ONLY: NM, MNTIF, MNE, MNP
      USE HARM, ONLY : IHARIND, CHARMV
      USE VERSION
      USE PRE_GLOBAL, ONLY: NSTAE, NOUTE, RUNDES, RUNID, AGRID,
     &    DTDP, TAU0, STATIM, REFTIM, RNDAY, DRAMP, A00, B00, C00,
     &    H0, SL0, SF0, CF, ESLM, CORI, DP, NBDV, NBVV,
     &    SLEL, SFEL, ICS, NT, NTRSPE, IHOT, NOLIBF, NOLIFA, IM,
     &    XEV,YEV,XEM,YEM,XEL,YEL,X,Y,NELG,DT, NTIF, NODECODE,
     &    NOLICA, NOLICAT, NWP, NCOR, NTIP, NWS, NRS, NRAMP, NBFR,
     &    NHY, NOPE, NETA, NBOU, NVEL, RAD2DEG, G, SL1, SF1,
     &    NVDLL, NVELL, IBTYPE, IBTYPEE,
     &    NSTAV, NOUTV, SLEV, SFEV, NTRSPV,
     &    NSTAM, NOUTM, SLEM, SFEM, NTRSPM,
     &    NOUTGE, SLAM, SFEA, NDSETSE,
     &    NOUTGV, NDSETSV,
     &    NOUTGW, NDSETSW,
     &    IESTP,NSCOUE,IVSTP,NSCOUV,ICSTP,NSCOUC,IPSTP,IWSTP,NSCOUM,
     &    IGEP,NSCOUGE,IGVP,NSCOUGV,IGCP,NSCOUGC,IGPP,IGWP,NSCOUGW,
     &    nhstar,nhsinc,
     &    title, institution,
     &    source, history, references, comments, host, convention,
     &    contact, base_date, NABOUT, NSCREEN,
     &    I3DSD, NFEN, NSTA3DD, I3DSV, NSTA3DV, I3DST, NSTA3DT,
     &    I3DGD, I3DGV, I3DGT, IDEN, islip, kp, z0s, z0b, theta1,
     &    theta2, ievc, evmin, evcon, alp1, alp2, alp3, igc, nlsd,
     &    nvsd, nltd, nvtd, alp4, C3D, X3DD, Y3DD, SL3DD, SF3DD,
     &    X3DV, Y3DV, SL3DV, SF3DV, X3DT, Y3DT, SL3DT, SF3DT, C3D
      USE NETCDFIO, ONLY : setADCIRCParameters, initNetCDFOutputFile,
     &    initNetCDFHotstart3D, initNetCDFHotstart,
     &    initNetCDFHotstartHarmonic,
     &    initNetCDFHotstartHarmonicMeansVariances,
     &    freeNetCDFCoord
      USE GLOBAL, ONLY : OutputDataDescript_t, screenUnit, myProc
#ifdef ADCSWAN
      USE GLOBAL, ONLY : OutputDataDescript_t, screenUnit, myProc,
     &                   SWAN_OutputHS,SWAN_OutputDIR,SWAN_OutputTM01,
     &                   SWAN_OutputTPS,SWAN_OutputWIND,SWAN_OutputTM02,
     &                   SWAN_OutputTMM10
#else
      USE GLOBAL, ONLY : OutputDataDescript_t, screenUnit, myProc
#endif
C
      IMPLICIT NONE
      LOGICAL reterr
C
C     The initialization of the output data descriptors for each ADCIRC
C     output file had to be cut and pasted from write_output.F. At some
C     point in the future, adcprep will be part of adcirc, making this
C     unfortunate cut-and-paste duplication unnecessary.
      type(OutputDataDescript_t), SAVE :: ElevStaDescript    ! fort.61
      type(OutputDataDescript_t), SAVE :: VelStaDescript     ! fort.62
      type(OutputDataDescript_t), SAVE :: ElevDescript       ! fort.63
      type(OutputDataDescript_t), SAVE :: VelDescript        ! fort.64
      type(OutputDataDescript_t), SAVE :: PrStaDescript      ! fort.71
      type(OutputDataDescript_t), SAVE :: WindVelStaDescript ! fort.72
      type(OutputDataDescript_t), SAVE :: PrDescript         ! fort.73
      type(OutputDataDescript_t), SAVE :: WindVelDescript    ! fort.74
      type(OutputDataDescript_t), SAVE :: EtaMaxDescript     ! maxele.63
      type(OutputDataDescript_t), SAVE :: UMaxDescript       ! maxvel.63
      type(OutputDataDescript_t), SAVE :: PrMinDescript      ! minpr.63
      type(OutputDataDescript_t), SAVE :: WVMaxDescript      ! maxwvel.63
      type(OutputDataDescript_t), SAVE :: RSMaxDescript    

!      tcm v50.75 moved RSDescript outside of the ifdef adcswan for use with 
!      nrs = 3 or nrs= 4
      type(OutputDataDescript_t), SAVE :: RSDescript

#ifdef ADCSWAN
Cobell 20120510: SWAN Output Data
      type(OutputDataDescript_t), SAVE :: SwanHSDescript
      type(OutputDataDescript_t), SAVE :: SwanDIRDescript
      type(OutputDataDescript_t), SAVE :: SwanTM01Descript
      type(OutputDataDescript_t), SAVE :: SwanTPSDescript
      type(OutputDataDescript_t), SAVE :: SwanWindDescript
      type(OutputDataDescript_t), SAVE :: SwanTM02Descript
      type(OutputDataDescript_t), SAVE :: SwanTMM10Descript
!      type(OutputDataDescript_t), SAVE :: RSDescript   ! tcm v50.75 moved
      type(OutputDataDescript_t), SAVE :: SwanHSMaxDescript
      type(OutputDataDescript_t), SAVE :: SwanDIRMaxDescript
      type(OutputDataDescript_t), SAVE :: SwanTM01MaxDescript
      type(OutputDataDescript_t), SAVE :: SwanTPSMaxDescript
      type(OutputDataDescript_t), SAVE :: SwanWindMaxDescript
      type(OutputDataDescript_t), SAVE :: SwanTM02MaxDescript
      type(OutputDataDescript_t), SAVE :: SwanTMM10MaxDescript
#endif
C     3D output data
      type(OutputDataDescript_t), save :: SigTStaDescript    ! fort.41
      type(OutputDataDescript_t), save :: SalStaDescript
      type(OutputDataDescript_t), save :: TempStaDescript
      type(OutputDataDescript_t), save :: QSurfKp1Descript
      type(OutputDataDescript_t), save :: RealQStaDescript   ! fort.42
      type(OutputDataDescript_t), save :: ImaginaryQStaDescript
      type(OutputDataDescript_t), save :: WZStaDescript
      type(OutputDataDescript_t), save :: Q20StaDescript     ! fort.43
      type(OutputDataDescript_t), save :: LStaDescript
      type(OutputDataDescript_t), save :: EVStaDescript
      type(OutputDataDescript_t), save :: SigTDescript       ! fort.44
      type(OutputDataDescript_t), save :: SalDescript
      type(OutputDataDescript_t), save :: TempDescript
      type(OutputDataDescript_t), save :: RealQDescript      ! fort.45
      type(OutputDataDescript_t), save :: ImaginaryQDescript
      type(OutputDataDescript_t), save :: WZDescript
      type(OutputDataDescript_t), save :: Q20Descript        ! fort.46
      type(OutputDataDescript_t), save :: LDescript
      type(OutputDataDescript_t), save :: EVDescript
C     For hotstart files:
      type(OutputDataDescript_t), SAVE :: Elev1Descript
      type(OutputDataDescript_t), SAVE :: Elev2Descript
      type(OutputDataDescript_t), SAVE :: CH1Descript
      type(OutputDataDescript_t), SAVE :: EtaDiscDescript
      type(OutputDataDescript_t), SAVE :: NodeCodeDescript
      type(OutputDataDescript_t), SAVE :: NOFFDescript
C     for hotstart 3D data
      type(OutputDataDescript_t),SAVE :: Duudescript
      type(OutputDataDescript_t),SAVE :: Duvdescript
      type(OutputDataDescript_t),SAVE :: Dvvdescript
      type(OutputDataDescript_t),SAVE :: Uudescript
      type(OutputDataDescript_t),SAVE :: Vvdescript
      type(OutputDataDescript_t),SAVE :: Bsxdescript
      type(OutputDataDescript_t),SAVE :: Bsydescript
C     for hotstart harmonic analysis
      type(OutputDataDescript_t), SAVE :: HarmElevFDLVDescript
      type(OutputDataDescript_t), SAVE :: HarmElevSLVDescript
      type(OutputDataDescript_t), SAVE :: HarmUVelFDLVDescript
      type(OutputDataDescript_t), SAVE :: HarmVVelFDLVDescript
      type(OutputDataDescript_t), SAVE :: HarmUVelSLVDescript
      type(OutputDataDescript_t), SAVE :: HarmVVelSLVDescript
C     for hotstart harmoinc analysis means and variance calculations
      type(OutputDataDescript_t), SAVE :: ELAVDescript
      type(OutputDataDescript_t), SAVE :: ELVADescript
      type(OutputDataDescript_t), SAVE :: XVELAVDescript
      type(OutputDataDescript_t), SAVE :: YVELAVDescript
      type(OutputDataDescript_t), SAVE :: XVELVADescript
      type(OutputDataDescript_t), SAVE :: YVELVADescript
C
      INTEGER numHotstartWrites ! number writes to hot start files
      INTEGER nextLun           ! next LUN to write to, after initial write
C
C     fort.61
      ElevStaDescript % specifier            = NOUTE
      ElevStaDescript % initial_value        = 0.0
      ElevStaDescript % num_items_per_record = 1
      ElevStaDescript % num_fd_records       = NSTAE
      ElevStaDescript % num_records_this     = NSTAE
      ElevStaDescript % ConsiderWetDry       = .TRUE.
      ElevStaDescript % alternate_value      = -99999.0
      ElevStaDescript % field_name           = 'ElevSta'
      IF (ICS.eq.2) THEN
         ElevStaDescript % x_coord           => SLEL
         ElevStaDescript % y_coord           => SFEL
      ELSE
         ElevStaDescript % x_coord           => XEL
         ElevStaDescript % y_coord           => YEL
      ENDIF
C
C     fort.62
      VelStaDescript % specifier            = NOUTV
      VelStaDescript % initial_value        = 0.0
      VelStaDescript % num_items_per_record = 2
      VelStaDescript % num_fd_records       = NSTAV
      VelStaDescript % num_records_this     = NSTAV
      VelStaDescript % ConsiderWetDry       = .FALSE.
      VelStaDescript % alternate_value      = 0.0
      VelStaDescript % field_name           = 'VelSta'
      IF (ICS.eq.2) THEN
         VelStaDescript % x_coord           => SLEV
         VelStaDescript % y_coord           => SFEV
      ELSE
         VelStaDescript % x_coord           => XEV
         VelStaDescript % y_coord           => YEV
      ENDIF
C
C     fort.63
      ElevDescript % specifier            = NOUTGE
      ElevDescript % initial_value        = 0.0
      ElevDescript % num_items_per_record = 1
      ElevDescript % num_fd_records       = MNP
      ElevDescript % num_records_this     = MNP
      ElevDescript % ConsiderWetDry       = .TRUE.
      ElevDescript % alternate_value      = -99999.0
      ElevDescript % field_name           = 'Elev'
C
C     fort.64
      VelDescript % specifier            = NOUTGV
      VelDescript % initial_value        = 0.0
      VelDescript % num_items_per_record = 2
      VelDescript % num_fd_records       = MNP
      VelDescript % num_records_this     = MNP
      VelDescript % ConsiderWetDry       = .FALSE.
      VelDescript % alternate_value      = 0.0
      VelDescript % field_name           = 'Vel'
C
C     maxele.63
      EtaMaxDescript % specifier            = NOUTGE
      EtaMaxDescript % initial_value        = 0.0
      EtaMaxDescript % num_items_per_record = 1
      EtaMaxDescript % num_fd_records       = MNP
      EtaMaxDescript % num_records_this     = MNP
      EtaMaxDescript % ConsiderWetDry       = .FALSE.
      EtaMaxDescript % alternate_value      = 0.0
      EtaMaxDescript % field_name           = 'EtaMax'
C
C     maxvel.63
      UMaxDescript % specifier            = NOUTGV
      UMaxDescript % initial_value        = 0.0
      UMaxDescript % num_items_per_record = 2
      UMaxDescript % num_fd_records       = MNP
      UMaxDescript % num_records_this     = MNP
      UMaxDescript % ConsiderWetDry       = .FALSE.
      UMaxDescript % alternate_value      = 0.0
      UMaxDescript % field_name           = 'UMax'      
C
C     fort.71
      PrStaDescript % specifier            = NOUTM
      PrStaDescript % initial_value        = 0.0
      PrStaDescript % num_items_per_record = 1
      PrStaDescript % num_fd_records       = NSTAM
      PrStaDescript % num_records_this     = NSTAM
      PrStaDescript % ConsiderWetDry       = .FALSE.
      PrStaDescript % alternate_value      = 0.0
      PrStaDescript % field_name           = 'PrSta'
      IF (ICS.eq.2) THEN
         PrStaDescript % x_coord           => SLEM
         PrStaDescript % y_coord           => SFEM
      ELSE
         PrStaDescript % x_coord           => XEM
         PrStaDescript % y_coord           => YEM
      ENDIF
C
C     fort.72
      WindVelStaDescript % specifier            = NOUTM
      WindVelStaDescript % initial_value        = 0.0
      WindVelStaDescript % num_items_per_record = 2
      WindVelStaDescript % num_fd_records       = NSTAM
      WindVelStaDescript % num_records_this     = NSTAM
      WindVelStaDescript % ConsiderWetDry       = .FALSE.
      WindVelStaDescript % alternate_value      = 0.0
      WindVelStaDescript % field_name           = 'WindVelSta'
      IF (ICS.eq.2) THEN
         WindVelStaDescript % x_coord           => SLEM
         WindVelStaDescript % y_coord           => SFEM
      ELSE
         WindVelStaDescript % x_coord           => XEM
         WindVelStaDescript % y_coord           => YEM
      ENDIF
C
C     fort.73
      PrDescript % specifier            = NOUTGW
      PrDescript % initial_value        = 0.0
      PrDescript % num_items_per_record = 1
      PrDescript % num_fd_records       = MNP
      PrDescript % num_records_this     = MNP
      PrDescript % ConsiderWetDry       = .FALSE.
      PrDescript % alternate_value      = 0.0
      PrDescript % field_name           = 'Pr'
C
C     fort.74
      WindVelDescript % specifier            = NOUTGW
      WindVelDescript % initial_value        = 0.0
      WindVelDescript % num_items_per_record = 2
      WindVelDescript % num_fd_records       = MNP
      WindVelDescript % num_records_this     = MNP
      WindVelDescript % ConsiderWetDry       = .FALSE.
      WindVelDescript % alternate_value      = 0.0
      WindVelDescript % field_name           = 'WindVel'
C
C     minpr.63
      PrMinDescript % specifier            = NOUTGW
      PrMinDescript % initial_value        = 0.0
      PrMinDescript % num_items_per_record = 1
      PrMinDescript % num_fd_records       = MNP
      PrMinDescript % num_records_this     = MNP
      PrMinDescript % ConsiderWetDry       = .FALSE.
      PrMinDescript % alternate_value      = 0.0
      PrMinDescript % field_name           = 'PrMin'
C
C     maxwvel.63
      WVMaxDescript % specifier            = NOUTGW
      WVMaxDescript % initial_value        = 0.0
      WVMaxDescript % num_items_per_record = 2
      WVMaxDescript % num_fd_records       = MNP
      WVMaxDescript % num_records_this     = MNP
      WVMaxDescript % ConsiderWetDry       = .FALSE.
      WVMaxDescript % alternate_value      = 0.0
      WVMaxDescript % field_name           = 'WVMax'

      RSMaxDescript % specifier            = NOUTGW
      RSMaxDescript % initial_value        = 0.0
      RSMaxDescript % num_items_per_record = 2
      RSMaxDescript % num_fd_records       = MNP
      RSMaxDescript % num_records_this     = MNP
      RSMaxDescript % ConsiderWetDry       = .FALSE.
      RSMaxDescript % alternate_value      = 0.0
      RSMaxDescript % field_name           = "RSMax"
      RSMaxDescript % file_name            = "maxrs.63"         

!  tcm v50.75 removed ifdef adcswan to allow for use whenever nrs=3 or nrs=4
!#ifdef ADCSWAN
Cobell 20120510: SWAN Output Data
C........Radiation Stress
         RSDescript % specifier            = NOUTGW
         RSDescript % initial_value        = 0.0
         RSDescript % num_items_per_record = 2
         RSDescript % num_fd_records       = MNP
         RSDescript % num_records_this     = MNP
         RSDescript % ConsiderWetDry       = .FALSE.
         RSDescript % alternate_value      = -99999.0
         RSDescript % field_name           = "rads"
         RSDescript % file_name            = "rads.64"

! tcm v50.75 moved the ifdef adcswan down past the RSDescript
#ifdef ADCSWAN
C........Significant Wave Height (HS)
         SwanHSDescript % specifier            = NOUTGW
         SwanHSDescript % initial_value        = 0.0
         SwanHSDescript % num_items_per_record = 1
         SwanHSDescript % num_fd_records       = MNP
         SwanHSDescript % num_records_this     = MNP
         SwanHSDescript % ConsiderWetDry       = .FALSE.
         SwanHSDescript % alternate_value      = -99999.0
         SwanHSDescript % field_name           = "swan_HS"
         SwanHSDescript % file_name            = "swan_HS.63"

         SwanHSMaxDescript % specifier            = NOUTGW
         SwanHSMaxDescript % initial_value        = 0.0
         SwanHSMaxDescript % num_items_per_record = 1
         SwanHSMaxDescript % num_fd_records       = MNP
         SwanHSMaxDescript % num_records_this     = MNP
         SwanHSMaxDescript % ConsiderWetDry       = .FALSE.
         SwanHSMaxDescript % alternate_value      = -99999.0
         SwanHSMaxDescript % field_name           = "swan_HS_max"
         SwanHSMaxDescript % file_name            = "swan_HS_max.63"

C........Mean Wave Direction (DIR)
         SwanDIRDescript % specifier            = NOUTGW
         SwanDIRDescript % initial_value        = 0.0
         SwanDIRDescript % num_items_per_record = 1
         SwanDIRDescript % num_fd_records       = MNP
         SwanDIRDescript % num_records_this     = MNP
         SwanDIRDescript % ConsiderWetDry       = .FALSE.
         SwanDIRDescript % alternate_value      = -99999.0
         SwanDIRDescript % field_name           = "swan_DIR"
         SwanDIRDescript % file_name            = "swan_DIR.63"

         SwanDIRMaxDescript % specifier            = NOUTGW
         SwanDIRMaxDescript % initial_value        = 0.0
         SwanDIRMaxDescript % num_items_per_record = 1
         SwanDIRMaxDescript % num_fd_records       = MNP
         SwanDIRMaxDescript % num_records_this     = MNP
         SwanDIRMaxDescript % ConsiderWetDry       = .FALSE.
         SwanDIRMaxDescript % alternate_value      = -99999.0
         SwanDIRMaxDescript % field_name           = "swan_DIR_max"
         SwanDIRMaxDescript % file_name            = "swan_DIR_max.63"

C........Mean Wave Period (TM01)
         SwanTM01Descript % specifier            = NOUTGW
         SwanTM01Descript % initial_value        = 0.0
         SwanTM01Descript % num_items_per_record = 1
         SwanTM01Descript % num_fd_records       = MNP
         SwanTM01Descript % num_records_this     = MNP
         SwanTM01Descript % ConsiderWetDry       = .FALSE.
         SwanTM01Descript % alternate_value      = -99999.0
         SwanTM01Descript % field_name           = "swan_TM01"
         SwanTM01Descript % file_name            = "swan_TM01.63"

         SwanTM01MaxDescript % specifier            = NOUTGW
         SwanTM01MaxDescript % initial_value        = 0.0
         SwanTM01MaxDescript % num_items_per_record = 1
         SwanTM01MaxDescript % num_fd_records       = MNP
         SwanTM01MaxDescript % num_records_this     = MNP
         SwanTM01MaxDescript % ConsiderWetDry       = .FALSE.
         SwanTM01MaxDescript % alternate_value      = -99999.0
         SwanTM01MaxDescript % field_name           = "swan_TM01_max"
         SwanTM01MaxDescript % file_name            = "swan_TM01_max.63"

C........Peak Wave Period (TPS)
         SwanTPSDescript % specifier            = NOUTGW
         SwanTPSDescript % initial_value        = 0.0
         SwanTPSDescript % num_items_per_record = 1
         SwanTPSDescript % num_fd_records       = MNP
         SwanTPSDescript % num_records_this     = MNP
         SwanTPSDescript % ConsiderWetDry       = .FALSE.
         SwanTPSDescript % alternate_value      = -99999.0
         SwanTPSDescript % field_name           = "Swan_TPS"
         SwanTPSDescript % file_name            = "swan_TPS.63"

         SwanTPSMaxDescript % specifier            = NOUTGW
         SwanTPSMaxDescript % initial_value        = 0.0
         SwanTPSMaxDescript % num_items_per_record = 1
         SwanTPSMaxDescript % num_fd_records       = MNP
         SwanTPSMaxDescript % num_records_this     = MNP
         SwanTPSMaxDescript % ConsiderWetDry       = .FALSE.
         SwanTPSMaxDescript % alternate_value      = -99999.0
         SwanTPSMaxDescript % field_name           = "Swan_TPS_max"
         SwanTPSMaxDescript % file_name            = "swan_TPS_max.63"

C........SWAN Wind Values (WINDX,WINDY)
         SwanWindDescript % specifier            = NOUTGW
         SwanWindDescript % initial_value        = 0.0
         SwanWindDescript % num_items_per_record = 2
         SwanWindDescript % num_fd_records       = MNP
         SwanWindDescript % num_records_this     = MNP
         SwanWindDescript % ConsiderWetDry       = .FALSE.
         SwanWindDescript % alternate_value      = -99999.0
         SwanWindDescript % field_name           = "swan_WIND"
         SwanWindDescript % file_name            = "swan_WIND.64"

         SwanWindMaxDescript % specifier            = NOUTGW
         SwanWindMaxDescript % initial_value        = 0.0
         SwanWindMaxDescript % num_items_per_record = 1
         SwanWindMaxDescript % num_fd_records       = MNP
         SwanWindMaxDescript % num_records_this     = MNP
         SwanWindMaxDescript % ConsiderWetDry       = .FALSE.
         SwanWindMaxDescript % alternate_value      = -99999.0
         SwanWindMaxDescript % field_name           = "swan_WIND_max"
         SwanWindMaxDescript % file_name            = "swan_WIND_max.63"

C........Mean Wave Period (TM02)
         SwanTM02Descript % specifier            = NOUTGW
         SwanTM02Descript % initial_value        = 0.0
         SwanTM02Descript % num_items_per_record = 1
         SwanTM02Descript % num_fd_records       = MNP
         SwanTM02Descript % num_records_this     = MNP
         SwanTM02Descript % ConsiderWetDry       = .FALSE.
         SwanTM02Descript % alternate_value      = -99999.0
         SwanTM02Descript % field_name           = "swan_TM02"
         SwanTM02Descript % file_name            = "swan_TM02.63"

         SwanTM02MaxDescript % specifier            = NOUTGW
         SwanTM02MaxDescript % initial_value        = 0.0
         SwanTM02MaxDescript % num_items_per_record = 1
         SwanTM02MaxDescript % num_fd_records       = MNP
         SwanTM02MaxDescript % num_records_this     = MNP
         SwanTM02MaxDescript % ConsiderWetDry       = .FALSE.
         SwanTM02MaxDescript % alternate_value      = -99999.0
         SwanTM02MaxDescript % field_name           = "swan_TM02_max"
         SwanTM02MaxDescript % file_name            = "swan_TM02_max.63"

C........Mean Wave Period (TMM10)
         SwanTMM10Descript % specifier            = NOUTGW
         SwanTMM10Descript % initial_value        = 0.0
         SwanTMM10Descript % num_items_per_record = 1
         SwanTMM10Descript % num_fd_records       = MNP
         SwanTMM10Descript % num_records_this     = MNP
         SwanTMM10Descript % ConsiderWetDry       = .FALSE.
         SwanTMM10Descript % alternate_value      = -99999.0
         SwanTMM10Descript % field_name           = "swan_TMM10"
         SwanTMM10Descript % file_name            = "swan_TMM10.63"

         SwanTMM10MaxDescript % specifier            = NOUTGW
         SwanTMM10MaxDescript % initial_value        = 0.0
         SwanTMM10MaxDescript % num_items_per_record = 1
         SwanTMM10MaxDescript % num_fd_records       = MNP
         SwanTMM10MaxDescript % num_records_this     = MNP
         SwanTMM10MaxDescript % ConsiderWetDry       = .FALSE.
         SwanTMM10MaxDescript % alternate_value      = -99999.0
         SwanTMM10MaxDescript % field_name           = "swan_TMM10_max"
         SwanTMM10MaxDescript % file_name          = "swan_TMM10_max.63"

#endif
C     3D data
C     fort.41
      SigTStaDescript % specifier            = I3DSD
      SigTStaDescript % initial_value        = 0.0
      SigTStaDescript % num_items_per_record = NFEN
      SigTStaDescript % num_fd_records       = NSTA3DD
      SigTStaDescript % num_records_this     = NSTA3DD
      SigTStaDescript % field_name           = 'SigmaTStations'
      IF (ICS.eq.2) THEN
         SigTStaDescript % x_coord           => SL3DD
         SigTStaDescript % y_coord           => SF3DD
      ELSE
         SigTStaDescript % x_coord           => X3DD
         SigTStaDescript % y_coord           => Y3DD
      ENDIF

      SalStaDescript % specifier            = I3DSD
      SalStaDescript % initial_value        = 0.0
      SalStaDescript % num_items_per_record = NFEN
      SalStaDescript % num_fd_records       = NSTA3DD
      SalStaDescript % num_records_this     = NSTA3DD
      SalStaDescript % field_name           = 'SalinityStations'
      IF (ICS.eq.2) THEN
         SalStaDescript % x_coord           => SL3DD
         SalStaDescript % y_coord           => SF3DD
      ELSE
         SalStaDescript % x_coord           => X3DD
         SalStaDescript % y_coord           => Y3DD
      ENDIF

      TempStaDescript % specifier            = I3DSD
      TempStaDescript % initial_value        = 0.0
      TempStaDescript % num_items_per_record = NFEN
      TempStaDescript % num_fd_records       = NSTA3DD
      TempStaDescript % num_records_this     = NSTA3DD
      TempStaDescript % field_name           = 'TemperatureStations'
      IF (ICS.eq.2) THEN
         TempStaDescript % x_coord           => SL3DD
         TempStaDescript % y_coord           => SF3DD
      ELSE
         TempStaDescript % x_coord           => X3DD
         TempStaDescript % y_coord           => Y3DD
      ENDIF
C
C     fort.42
      RealQStaDescript % specifier            =  I3DSV
      RealQStaDescript % initial_value        =  0.0
      RealQStaDescript % num_items_per_record =  NFEN
      RealQStaDescript % num_fd_records       =  NSta3DV
      RealQStaDescript % num_records_this     =  NSta3DV
      RealQStaDescript % field_name           = 'RealQStations'
      IF (ICS.eq.2) THEN
         RealQStaDescript % x_coord           => SL3DV
         RealQStaDescript % y_coord           => SF3DV
      ELSE
         RealQStaDescript % x_coord           => X3DV
         RealQStaDescript % y_coord           => Y3DV
      ENDIF
      ImaginaryQStaDescript % specifier            =  I3DSV
      ImaginaryQStaDescript % initial_value        =  0.0
      ImaginaryQStaDescript % num_items_per_record =  NFEN
      ImaginaryQStaDescript % num_fd_records       =  NSTA3DV
      ImaginaryQStaDescript % num_records_this     =  NSTA3DV
      ImaginaryQStaDescript % field_name      ='ImaginaryQStations'
      IF (ICS.eq.2) THEN
         ImaginaryQStaDescript % x_coord           => SL3DV
         ImaginaryQStaDescript % y_coord           => SF3DV
      ELSE
         ImaginaryQStaDescript % x_coord           => X3DV
         ImaginaryQStaDescript % y_coord           => Y3DV
      ENDIF
      WZStaDescript % specifier            =  I3DSV
      WZStaDescript % initial_value        =  0.0
      WZStaDescript % num_items_per_record =  NFEN
      WZStaDescript % num_fd_records       =  NSTA3DV
      WZStaDescript % num_records_this     =  NSTA3DV
      WZStaDescript % field_name           = 'WZStations'
      IF (ICS.eq.2) THEN
         WZStaDescript % x_coord           => SL3DV
         WZStaDescript % y_coord           => SF3DV
      ELSE
         WZStaDescript % x_coord           => X3DV
         WZStaDescript % y_coord           => Y3DV
      ENDIF
C
C     fort.43
      Q20StaDescript % specifier            =  I3DST
      Q20StaDescript % initial_value        =  0.0
      Q20StaDescript % num_items_per_record =  NFEN
      Q20StaDescript % num_fd_records       =  NSta3DT
      Q20StaDescript % num_records_this     =  NSta3DT
      Q20StaDescript % field_name           = 'q20Stations'
      IF (ICS.eq.2) THEN
         Q20StaDescript % x_coord           => SL3DT
         Q20StaDescript % y_coord           => SF3DT
      ELSE
         Q20StaDescript % x_coord           => X3DT
         Q20StaDescript % y_coord           => Y3DT
      ENDIF

      LStaDescript % specifier        =  I3DST
      LStaDescript % initial_value        =  0.0
      LStaDescript % num_items_per_record =  NFEN
      LStaDescript % num_fd_records       =  NSTA3DT
      LStaDescript % num_records_this     =  NSTA3DT
      LStaDescript % field_name           = 'LStations'
      IF (ICS.eq.2) THEN
         LStaDescript % x_coord           => SL3DT
         LStaDescript % y_coord           => SF3DT
      ELSE
         LStaDescript % x_coord           => X3DT
         LStaDescript % y_coord           => Y3DT
      ENDIF

      EVStaDescript % specifier            =  I3DST
      EVStaDescript % initial_value        =  0.0
      EVStaDescript % num_items_per_record =  NFEN
      EVStaDescript % num_fd_records       =  NSTA3DT
      EVStaDescript % num_records_this     =  NSTA3DT
      EVStaDescript % field_name           = 'EVStations'
      IF (ICS.eq.2) THEN
         EVStaDescript % x_coord           => SL3DT
         EVStaDescript % y_coord           => SF3DT
      ELSE
         EVStaDescript % x_coord           => X3DT
         EVStaDescript % y_coord           => Y3DT
      ENDIF
C
C     fort.44
      SigTDescript % specifier            =  I3DGD
      SigTDescript % initial_value        =  0.0
      SigTDescript % num_items_per_record =  NFEN
      SigTDescript % num_fd_records        =  MNP
      SigTDescript % num_records_this     =  MNP
      SigTDescript % field_name           = 'SigmaT'
      SalDescript % specifier            =  I3DGD
      SalDescript % initial_value        =  0.0
      SalDescript % num_items_per_record =  NFEN
      SalDescript % num_fd_records       =  MNP
      SalDescript % num_records_this     =  MNP
      SalDescript % field_name           = 'Salinity'
      TempDescript % specifier            =  I3DGD
      TempDescript % initial_value        =  0.0
      TempDescript % num_items_per_record =  NFEN
      TempDescript % num_fd_records       =  MNP
      TempDescript % num_records_this     =  MNP
      TempDescript % field_name           = 'Temperature'
C
C     fort.45
      RealQDescript % specifier            =  I3DGV
      RealQdescript % initial_value        =  0.0
      RealQDescript % num_items_per_record =  NFEN
      RealQDescript % num_fd_records       =  MNP
      RealQDescript % num_records_this     =  MNP
      RealQDescript % field_name           = 'RealQ'
      ImaginaryQDescript % specifier            =  I3DGV
      ImaginaryQDescript % initial_value        =  0.0
      ImaginaryQDescript % num_items_per_record =  NFEN
      ImaginaryQDescript % num_fd_records       =  MNP
      ImaginaryQDescript % num_records_this     =  MNP
      ImaginaryQDescript % field_name           = 'ImaginaryQ'
      WZDescript % specifier            =  I3DGV
      WZDescript % initial_value        =  0.0
      WZDescript % num_items_per_record =  NFEN
      WZDescript % num_fd_records       =  MNP
      WZDescript % num_records_this     =  MNP
      WZDescript % field_name           = 'WZ'
C
C     fort.46
      Q20Descript % specifier            =  I3DGT
      Q20Descript % initial_value        =  0.0
      Q20Descript % num_items_per_record =  NFEN
      Q20Descript % num_fd_records       =  MNP
      Q20Descript % num_records_this     =  MNP
      Q20Descript % field_name           = 'q20'
      LDescript % specifier            =  I3DGT
      LDescript % initial_value        =  0.0
      LDescript % num_items_per_record =  NFEN
      LDescript % num_fd_records       =  MNP
      LDescript % num_records_this     =  MNP
      LDescript % field_name           = 'L'
      EVDescript % specifier            =  I3DGT
      EVDescript % initial_value        =  0.0
      EVDescript % num_items_per_record =  NFEN
      EVDescript % num_fd_records       =  MNP
      EVDescript % num_records_this     =  MNP
      EVDescript % field_name           = 'EV'
C
C     fort.47
      QSurfKp1Descript % specifier            =  I3DGD
      QSurfKp1Descript % initial_value        =  0.0
      QSurfKp1Descript % num_items_per_record =  1
      QSurfKp1Descript % num_fd_records     =  MNP
      QSurfKp1Descript % num_records_this   =  MNP
      QSurfKp1Descript % field_name           = 'qsurfkp1'
C
C     fort.67 and fort.68
      Elev1Descript % specifier            = NHSTAR
      Elev1Descript % initial_value        = 0.0
      Elev2Descript % specifier            = NHSTAR
      Elev2Descript % initial_value        = 0.0
      CH1Descript % specifier            = NHSTAR
      CH1Descript % initial_value        = 0.0
      EtaDiscDescript % specifier            = NHSTAR
      EtaDiscDescript % initial_value        = 0.0
      NodeCodeDescript % specifier            = NHSTAR
      NodeCodeDescript % int_initial_value        = 0
      NOFFDescript % specifier            = NHSTAR
      NOFFDescript % int_initial_value        = 0
C     hotstart 3D
      Duudescript % specifier            =  NHSTAR
      Duudescript % initial_value        =  0.0
      Duudescript % num_items_per_record =  1
      Duudescript % num_fd_records       =  MNP
      Duudescript % num_records_this     =  MNP
      Duvdescript % specifier            =  NHSTAR
      Duvdescript % initial_value        =  0.0
      Duvdescript % num_items_per_record =  1
      Duvdescript % num_fd_records       =  MNP
      Duvdescript % num_records_this     =  MNP
      Dvvdescript % specifier            =  NHSTAR
      Dvvdescript % initial_value        =  0.0
      Dvvdescript % num_items_per_record =  1
      Dvvdescript % num_fd_records       =  MNP
      Dvvdescript % num_records_this     =  MNP
      Uudescript % specifier            =  NHSTAR
      Uudescript % initial_value        =  0.0
      Uudescript % num_items_per_record =  1
      Uudescript % num_fd_records       =  MNP
      Uudescript % num_records_this     =  MNP
      Vvdescript % specifier            =  NHSTAR
      Vvdescript % initial_value        =  0.0
      Vvdescript % num_items_per_record =  1
      Vvdescript % num_fd_records       =  MNP
      Vvdescript % num_records_this     =  MNP
      Bsxdescript % specifier            =  NHSTAR
      Bsxdescript % initial_value        =  0.0
      Bsxdescript % num_items_per_record =  1
      Bsxdescript % num_fd_records       =  MNP
      Bsxdescript % num_records_this     =  MNP
      Bsydescript % specifier            =  NHSTAR
      Bsydescript % initial_value        =  0.0
      Bsydescript % num_items_per_record =  1
      Bsydescript % num_fd_records       =  MNP
      Bsydescript % num_records_this     =  MNP
C     hotstart harmonic analysis
      HarmElevFDLVDescript % specifier            = NHSTAR
      HarmElevFDLVDescript % initial_value        = 0.0
      HarmElevFDLVDescript % num_fd_records       = MNP
      HarmElevSLVDescript % specifier            = NHSTAR
      HarmElevSLVDescript % initial_value        = 0.0
      HarmElevSLVDescript % num_fd_records       = NSTAE
      HarmUVelFDLVDescript % specifier            = NHSTAR
      HarmUVelFDLVDescript % initial_value        = 0.0
      HarmUVelFDLVDescript % num_fd_records       = MNP
      HarmVVelFDLVDescript % specifier            = NHSTAR
      HarmVVelFDLVDescript % initial_value        = 0.0
      HarmVVelFDLVDescript % num_fd_records       = MNP
      HarmUvelSLVDescript % specifier            = NHSTAR
      HarmUVelSLVDescript % initial_value        = 0.0
      HarmUVelSLVDescript % num_fd_records       = NSTAV
      HarmVVelSLVDescript % specifier            = NHSTAR
      HarmVVelSLVDescript % initial_value        = 0.0
      HarmVVelSLVDescript % num_fd_records       = NSTAV
C     hotstart means and variance calculations
      ELAVDescript % specifier            = NHSTAR
      ELAVDescript % initial_value        = 0.0
      ELAVDescript % num_fd_records       = MNP
      ELVADescript % specifier            = NHSTAR
      ELVADescript % initial_value        = 0.0
      ELVADescript % num_fd_records       = MNP
      XVELAVDescript % specifier            = NHSTAR
      XVELAVDescript % initial_value        = 0.0
      XVELAVDescript % num_fd_records       = MNP
      YVELAVDescript % specifier            = NHSTAR
      YVELAVDescript % initial_value        = 0.0
      YVELAVDescript % num_fd_records       = MNP
      XVELVADescript % specifier            = NHSTAR
      XVELVADescript % initial_value        = 0.0
      XVELVADescript % num_fd_records       = MNP
      YVELVADescript % specifier            = NHSTAR
      YVELVADescript % initial_value        = 0.0
      YVELVADescript % num_fd_records       = MNP
C
C     Need to populate the global and nodal attributes modules with
C     these parameters, since the netcdfio module relies on those
C     modules, rather than the pre_global module. Some day, adcprep
C     will be integrated with ADCIRC and this subroutine call will
C     not be needed.
      IF (.not.ALLOCATED(NODECODE)) THEN
         ALLOCATE(NODECODE(MNP))
      ENDIF
C
C     jgf49.44: Set parameters in global module based on the data
C     we collected in read_global.F and stored in pre_global.F.
      CALL setADCIRCParameters(
     &   base_date, MNE, NBOU,
     &   NVEL, NOPE, MNP, SL0, SF0, NBVV, NVDLL, NBDV, NVELL, X, Y,
     &   IBTYPE, IBTYPEE, SL1, SF1, NODECODE, G, FileFmtRev,
     &   FileFmtMinor, FileFmtMajor, im, iestp, nscoue, ivstp, nscouv,
     &   icstp, nscouc, ipstp, iwstp, nscoum, igep, nscouge, igvp,
     &   nscougv, igcp, nscougc, igpp, igwp, nscougw, NM,
     &   DP, RUNDES, AGRID, title, institution, source, history,
     &   references, comments, host, convention, contact, DT, ihot,
     &   ics, nolifa, nolica, nolicat, ncor, ntip, nws, nramp, statim,
     &   reftim, rnday, dramp, a00, b00, c00, h0, cori, ntif, nbfr,
     &   myProc, screenUnit, nolibf, nwp, tau0, cf, eslm,
     &   nstae, nstav, nstam, neta, nabout, nscreen,
     &   nfen, iden, islip, kp, z0s, z0b, theta1, theta2,
     &   ievc, evmin, evcon, alp1, alp2, alp3, igc, nlsd, nvsd, nltd,
     &   nvtd, alp4, C3D)
C
C     Create NetCDF output files for those output files where NetCDF
C     was specified.
      reterr = .false.
      CALL initNetCDFOutputFile(61, ElevStaDescript, reterr)
      CALL initNetCDFOutputFile(62, VelStaDescript, reterr)
      CALL initNetCDFOutputFile(63, ElevDescript, reterr)
      CALL initNetCDFOutputFile(64, VelDescript, reterr)
      CALL initNetCDFOutputFile(71, PrStaDescript, reterr)
      CALL initNetCDFOutputFile(72, WindVelStaDescript, reterr)
      CALL initNetCDFOutputFile(73, PrDescript, reterr)
      CALL initNetCDFOutputFile(74, WindVelDescript, reterr)
      CALL initNetCDFOutputFile(311, EtaMaxDescript, reterr)
      CALL initNetCDFOutputFile(312, UMaxDescript, reterr)
      CALL initNetCDFOutputFile(313, PrMinDescript, reterr)
      CALL initNetCDFOutputFile(314, WVMaxDescript, reterr)      
      CALL initNetCDFOutputFile(315,RSMaxDescript,reterr)
 
! tcm v50.75 moved ifdef adcswan below RSDescript only to allow
!  for use whenever nrs=3 or nrs=4     
!#ifdef ADCSWAN
Cobell 20120510: Added for SWAN NetCDF
      IF ((NRS.EQ.3).OR.(NRS.EQ.4)) THEN
        CALL initNetCDFOutputFile(164,RSDescript,reterr)
      ENDIF
! tcm v50.75 moved ifdef adcswan to here
#ifdef ADCSWAN
Cobell 20120510: Added for SWAN NetCDF
      IF(NRS.EQ.3)THEN
        IF(SWAN_OutputHS)THEN
            CALL initNetCDFOutputFile(301,SwanHSDescript,reterr)
            CALL initNetCDFOutputFile(316,SwanHSMaxDescript,reterr)
        ENDIF
        IF(SWAN_OutputDIR)THEN
            CALL initNetCDFOutputFile(302,SwanDIRDescript,reterr)
            CALL initNetCDFOutputFile(317,SwanDIRMaxDescript,reterr)
        ENDIF
        IF(SWAN_OutputTM01)THEN
            CALL initNetCDFOutputFile(303,SwanTM01Descript,reterr)
            CALL initNetCDFOutputFile(318,SwanTM01MaxDescript,reterr)
        ENDIF
        IF(SWAN_OutputTPS)THEN
            CALL initNetCDFOutputFile(304,SwanTPSDescript,reterr)
            CALL initNetCDFOutputFile(319,SwanTPSMaxDescript,reterr)
        ENDIF
        IF(SWAN_OutputWIND)THEN
            CALL initNetCDFOutputFile(305,SwanWINDDescript,reterr)
            CALL initNetCDFOutputFile(320,SwanWINDMaxDescript,reterr)
        ENDIF
        IF(SWAN_OutputTM02)THEN
            CALL initNetCDFOutputFile(306,SwanTM02Descript,reterr)
            CALL initNetCDFOutputFile(321,SwanTM02MaxDescript,reterr)
        ENDIF
        IF(SWAN_OutputTMM10)THEN
            CALL initNetCDFOutputFile(307,SwanTMM10Descript,reterr)
            CALL initNetCDFOutputFile(322,SwanTMM10MaxDescript,reterr)
        ENDIF
      ENDIF
#endif

C
      IF (C3D.eqv..true.) THEN
         CALL initNetCDFOutputFile(41, SigTStaDescript, reterr,
     &       SalStaDescript, TempStaDescript)
         CALL initNetCDFOutputFile(42, RealQStaDescript, reterr,
     &      ImaginaryQStaDescript, WZStaDescript)
         CALL initNetCDFOutputFile(43, Q20StaDescript, reterr,
     &      LStaDescript, EVStaDescript)
         CALL initNetCDFOutputFile(44, SigTDescript, reterr,
     &      SalDescript, TempDescript)
         CALL initNetCDFOutputFile(45, RealQDescript, reterr,
     &      ImaginaryQDescript, WZDescript)
         CALL initNetCDFOutputFile(46, Q20Descript, reterr,
     &      LDescript, EVDescript)
         CALL initNetCDFOutputFile(47, QSurfKp1Descript, reterr)
      ENDIF
C
C     Create NetCDF hotstart files if NetCDF was specified.
      IF ((NHSTAR.EQ.3).OR.(NHSTAR.EQ.367).OR.(NHSTAR.EQ.368).OR.
     & (NHSTAR.EQ.5).OR.(NHSTAR.EQ.567).OR.(NHSTAR.EQ.568)) THEN
         ! must init both hotstart files ... there is not enough information
         ! available to determine if they will both be needed
         IF ((IHOT.ne.367).and.(IHOT.ne.567)) THEN
            CALL initNetCDFHotstart(67, Elev1Descript,
     &         Elev2Descript, VelDescript, CH1Descript, EtaDiscDescript,
     &         NodeCodeDescript, NOFFDescript, reterr)
            IF (C3D.eqv..true.) THEN
               CALL initNetCDFHotstart3D(67,NHSTAR)
            ENDIF
         ENDIF
         IF ((IHOT.ne.368).and.(IHOT.ne.568)) THEN
            CALL initNetCDFHotstart(68, Elev1Descript,
     &         Elev2Descript, VelDescript, CH1Descript, EtaDiscDescript,
     &         NodeCodeDescript, NOFFDescript, reterr)
            IF (C3D.eqv..true.) THEN
               CALL initNetCDFHotstart3D(68,NHSTAR)
            ENDIF
         ENDIF
         IF (IHARIND.eq.1) THEN
            IF ((IHOT.ne.367).and.(IHOT.ne.567)) THEN
               CALL initNetCDFHotstartHarmonic(67,
     &            HarmElevFDLVDescript, HarmElevSLVDescript,
     &            HarmUVelFDLVDescript, HarmVVelFDLVDescript,
     &            HarmUVelSLVDescript, HarmVVelSLVDescript, reterr)
            ENDIF
            IF (IHOT.ne.368) THEN
               CALL initNetCDFHotstartHarmonic(68,
     &            HarmElevFDLVDescript, HarmElevSLVDescript,
     &            HarmUVelFDLVDescript, HarmVVelFDLVDescript,
     &            HarmUVelSLVDescript, HarmVVelSLVDescript, reterr)
            ENDIF
            IF (CHARMV.eqv..true.) THEN
               IF ((IHOT.ne.367).and.(IHOT.ne.567)) THEN
                  CALL initNetCDFHotstartHarmonicMeansVariances(
     &               67, ELAVDescript, ELVADescript,
     &               XVELAVDescript, YVELAVDescript, XVELVADescript,
     &               YVELVADescript, reterr)
               ENDIF
               IF ((IHOT.ne.368).and.(IHOT.ne.568)) THEN
                  CALL initNetCDFHotstartHarmonicMeansVariances(
     &               68, ELAVDescript, ELVADescript,
     &               XVELAVDescript, YVELAVDescript, XVELVADescript,
     &               YVELVADescript, reterr)
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      ! free up memory allocated for mesh and boundaries
      CALL freeNetCDFCoord()
#endif
c----------------------------------------------------------------------------
      END SUBROUTINE prepNetCDF
c----------------------------------------------------------------------------
