!----------------------------------------------------------------------------

!                           MODULE PREP

!----------------------------------------------------------------------------
!      current for ADCIRC v46.44   10/07/2006

!      Version 1.1 5/04/99 vjp
!      jjw fixes 053100
!      Revisions by rl 10/12/01, MEB 3/03, rl 3/03, rl 5/21/03, rl 5/18/04,
!                   vp 11/27/03 (by rl)
!      Revisions by MEB 4/04
!      jgf Created PREP11, rewrote PREP20 for ADCIRC v45.12 02/24/2006
!      jgf Created PREP13 for ADCRIC v46.00
!      vjp PREP13 updates for ADCIRC v46.44
!      vjp PREP67_68 localization updates for ADCIRC v46.44
!      vjp added Relocalize for relocalizing fort.13 and fort.15

!----------------------------------------------------------------------------


!---------------------------------------------------------------------------
!                S U B R O U T I N E   P R E P 1 0
!---------------------------------------------------------------------------

!     jgf46.28 from jgf45.16 This subroutine will break up the full
!     domain initial concentration file into subdomains. The fort.10
!     file may contain initial concentration for either 2D or 3D ADCIRC
!     runs.

!---------------------------------------------------------------------------
    SUBROUTINE PREP10()
!---------------------------------------------------------------------------
    USE PRE_GLOBAL
    use memory_usage
    IMPLICIT NONE
    integer :: nbytes = 0
    INTEGER :: I,J
    INTEGER :: iproc          ! subdomain index
    CHARACTER(80) :: header1, header2  ! header comments in unit 10 files
    CHARACTER(80) :: nvn_nvp           ! string representing nfen, np
    REAL(SZ) nvn           ! number of vertical nodes from unit 10 file
    REAL(SZ) nvp           ! number of horizontal nodes from unit 10 file
    INTEGER :: nhnn           ! horizontal nodes counter
    INTEGER :: nvnn           ! vertical nodes counter
    INTEGER :: sdu(nproc)     ! subdomain unit numbers for unit 10 files
    REAL(SZ), ALLOCATABLE :: fdData2D(:)   !(MNP)      full domain data
    REAL(SZ), ALLOCATABLE :: fdData3D(:,:) !(MNP,NFEN) full domain data
    REAL(SZ), ALLOCATABLE :: sdData2D(:)   !(MNP)      subdomain data
    REAL(SZ), ALLOCATABLE :: sdData3D(:,:) !(MNP,NFEN) subdomain data
    LOGICAL :: success        ! .TRUE. if all files opened successfully

    CALL OpenPrepFiles(10, 'initial concentration         ', &
    &      1, nproc, sdu, success)

    IF ( .NOT. success) THEN
        WRITE(*,*) 'WARNING: Unit 10 files not preprocessed.'
        RETURN ! note early return
    ENDIF

!     Read header information from full domain unit 10 file
    READ(10,80) header1
    READ(10,80) header2

!     Transcribe header information.
    DO iproc = 1, nproc
        WRITE(sdu(iproc),80)  header1
        WRITE(sdu(iproc),80)  header2
    ENDDO

!     Check node number data for consistency (paranoia).
    READ(10,80) nvn_nvp
    READ(nvn_nvp,*) nvn, nvp
    IF ( nvn /= nfen .OR. nvp /= nnodg ) then
        WRITE(*,*) 'ERROR: NVN or NVP not consistent with input data.'
        WRITE(*,*) 'NVN=',nvn,' although NFEN=',nfen
        WRITE(*,*) 'NVP=',nvp,' although NNODG=',nnodg
    ENDIF

!     Decompose concentration data

    IF (C2D_PTrans) THEN
    !     read in the full domain data
        ALLOCATE ( fdData2D(MNP) )
        nbytes = 8*mnp
        call memory_alloc(nbytes)
        DO i=1, NNODG
            READ(10,*) nhnn, fdData2D(nhnn)
        ENDDO
    !     write out subdomain data
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

    IF (C3D_PTrans) THEN
    !     read in the full domain data
        ALLOCATE ( fdData3D(MNP,NFEN) )
        nbytes = 8*mnp*nfen
        call memory_alloc(nbytes)
        DO i=1, NNODG
            DO j=1, nfen
                READ(10,*) nhnn, nvnn, fdData3D(nhnn,nvnn)
            ENDDO
        ENDDO
    !     write out subdomain data
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

!     Close full domain file and all subdomain files
    CLOSE(10)
    DO iproc=1, nproc
        CLOSE(sdu(iproc))
    ENDDO

    80 FORMAT(A80)

    call memory_status()
    RETURN
!---------------------------------------------------------------------------
    END SUBROUTINE PREP10
!---------------------------------------------------------------------------


!---------------------------------------------------------------------------
!                S U B R O U T I N E   P R E P 1 1
!---------------------------------------------------------------------------

!     jgf45.12 This subroutine will break up the full domain initial
!     density forcing file into subdomains. The fort.11 file may
!     contain initial density, temperature, and/or salinity depending on
!     the value of IDEN in the fort.15 file.

!     jgf45.12 This subroutine is designed to work for baroclinic 3D
!     runs only, not 2D runs.

!---------------------------------------------------------------------------
    SUBROUTINE PREP11()
!---------------------------------------------------------------------------
    USE PRE_GLOBAL
    use memory_usage
    IMPLICIT NONE
    integer :: nbytes = 0
    INTEGER :: I,J
    INTEGER :: iproc          ! subdomain index
    CHARACTER(80) :: header1, header2  ! header comments in unit 11 files
    CHARACTER(80) :: nvn_nvp           ! string representing nfen, np
    REAL(SZ) nvn           ! number of vertical nodes from unit 11 file
    REAL(SZ) nvp           ! number of horizontal nodes from unit 11 file
    INTEGER :: nhnn           ! horizontal nodes counter
    INTEGER :: nvnn           ! vertical nodes counter
    INTEGER :: sdu(nproc)     ! subdomain unit numbers for unit 11 files
    REAL(SZ), ALLOCATABLE :: fdData1(:,:) !(MNP,NFEN) full domain data
    REAL(SZ), ALLOCATABLE :: fdData2(:,:) !(MNP,NFEN) full domain data
    REAL(SZ), ALLOCATABLE :: sdData1(:,:) !(MNP,NFEN) subdomain data
    REAL(SZ), ALLOCATABLE :: sdData2(:,:) !(MNP,NFEN) subdomain data
    LOGICAL :: success        ! .TRUE. if all files opened successfully

    CALL OpenPrepFiles(11, 'initial density forcing       ', &
    &      1, nproc, sdu, success)

    IF ( .NOT. success) THEN
        WRITE(*,*) 'WARNING: Unit 11 files not preprocessed.'
        RETURN ! note early return
    ENDIF

!     Read header information from full domain unit 11 file
    READ(11,80) header1
    READ(11,80) header2

!     Transcribe header information.
    DO iproc = 1, nproc
        WRITE(sdu(iproc),80)  header1
        WRITE(sdu(iproc),80)  header2
    ENDDO

!     Check node number data for consistency (paranoia).
    READ(11,80) nvn_nvp
    READ(nvn_nvp,*) nvn, nvp
    IF ( nvn /= nfen .OR. nvp /= nnodg ) then
        WRITE(*,*) 'ERROR: NVN or NVP not consistent with input data.'
        WRITE(*,*) 'NVN=',nvn,' although NFEN=',nfen
        WRITE(*,*) 'NVP=',nvp,' although NNODG=',nnodg
    ENDIF

!     Decompose density forcing data; format based on value of IDEN.
!     jgf45.12 This is designed to work for baroclinic 3D runs only, not
!     2D runs.

!     read in the full domain data
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
            READ(11,*) nhnn, nvnn, &
            fdData1(nhnn,nvnn),fdData2(nhnn,nvnn)
        ENDDO
    ENDDO

    END SELECT

!     write out subdomain data
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

!     Close full domain file and all subdomain files
    CLOSE(11)
    DO iproc=1, nproc
        CLOSE(sdu(iproc))
    ENDDO

    80 FORMAT(A80)

    call memory_status()
    RETURN
!---------------------------------------------------------------------------
    END SUBROUTINE PREP11
!---------------------------------------------------------------------------




!---------------------------------------------------------------------------
!                S U B R O U T I N E   R E L O C A L I Z E
!---------------------------------------------------------------------------

!     This routine allows re-localizing selected files after a prepall
!     operation.  vjp 10/2006
!     tcm v50.21 20110610 -- changed I8 to I12 format specifications

!---------------------------------------------------------------------------
    SUBROUTINE RELOCALIZE()
!---------------------------------------------------------------------------
    USE PRE_GLOBAL
    USE PREP_WEIR,ONLY:NWEIRBNDRY
    use memory_usage
    IMPLICIT NONE
    integer :: nbytes = 0
    INTEGER :: I, J, K, IPROC, INDX, ITEMP, idumy
    CHARACTER(14) LOCFN
    CHARACTER(80) skipped
    INTEGER,ALLOCATABLE :: LUNP(:)
    LOGICAL,ALLOCATABLE :: ISWEIR(:)

    print *, "entering relocalize"
!      print *, "nproc = ", nproc

    allocate( lunp(nproc) )  ! logical unit number for each subdomain
    do iproc=1, nproc
        lunp(iproc) = 105 + (iproc-1)
    enddo

    if ( .NOT. allocated(nnodp)) then
        ALLOCATE(NNODP(NPROC))
        nbytes = 4*nproc
        call memory_alloc(nbytes)
    endif

    if ( .NOT. allocated(nelp)) then
        ALLOCATE(NELP(NPROC))
        nbytes = 4*nproc
        call memory_alloc(nbytes)
    endif

    if ( .NOT. allocated(netap)) then
        ALLOCATE(NETAP(NPROC))
        nbytes = 4*nproc
        call memory_alloc(nbytes)
    endif

    if ( .NOT. allocated(nfluxfp)) then
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
    ! Casey 100208: Changed I8 to I12.
        READ(LUNP(IPROC),'(8X,3I12)') NELG, MNEP, NELP(IPROC)
    ENDDO
!      print *, "nelg = ", nelg

    if ( .NOT. allocated(imap_el_lg)) then
        ALLOCATE ( IMAP_EL_LG(MNEP, NPROC) )
        nbytes = 4*nproc*mnep
        call memory_alloc(nbytes)
    endif

    DO IPROC = 1,NPROC
        DO I=1, NELP(IPROC)
        ! Casey 100208: Changed I8 to I12.
            READ(LUNP(IPROC),'(I12)') idumy
            IMAP_EL_LG(I,IPROC) = abs(idumy)
        ENDDO
    ENDDO

    print *, "from relocalize: reading local-to-global node maps"
    DO IPROC = 1,NPROC
        READ(LUNP(IPROC),'(8X,3I12)') NNODG, MNPP, NNODP(IPROC)   !tcm v50.21
    ENDDO
!      print *, "nnodg = ", nnodg

    if ( .NOT. allocated(imap_nod_lg)) then
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

!  This section for prep15
    IF ((PREP_15.eqv. .TRUE. ) .OR. (PREP_20.eqv. .TRUE. )) THEN
        print *, "from relocalize: reading nfluxf for each subdomain"
        DO IPROC = 1,NPROC
            READ(LUNP(IPROC),'(8X,I12)') NFLUXFP(IPROC)        !tcm v50.21
        ENDDO

        print *, "from relocalize: reading neta for each subdomain"
        DO IPROC = 1,NPROC
            READ(LUNP(IPROC),'(8X,3I12)') idumy, NETA_MAX, NETAP(IPROC)   !tcm v50.21
        ENDDO

        if ( .NOT. allocated(obnode_lg)) then
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
    if ( .NOT. allocated(itotproc)) then
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
        IF (ITOTPROC(I) > MNEI) MNEI = ITOTPROC(I)
    ENDDO
    print *, "MNEI = ", MNEI

    if ( .NOT. allocated(imap_nod_gl2)) then
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

! jgf50.35: Need this for --prep13 option.
    if ( .NOT. allocated(imap_nod_gl)) then
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
            IF (ITOTPROC(INDX) == 1) THEN
                IMAP_NOD_GL(1,INDX) = IPROC
                IMAP_NOD_GL(2,INDX) = J
            ENDIF
        ENDDO
    ENDDO
! c...If we're using time varying weirs, we need to
!     find out which processors have weir (4,24) boundaries
!     in the case that we're not passing through PREP14
    IF(USE_TVW)THEN
        ALLOCATE(ISWEIR(1:NNODG))
        ISWEIR(:) = .FALSE. 
        DO K=1,NBOU
            SELECT CASE(IBTYPE(K))
            CASE(4,24,5,25)
            DO J=1,NVELL(K)
                ISWEIR(NBVV(K,J))= .TRUE. 
                ISWEIR(IBCONNR(K,J))= .TRUE. 
            ENDDO
            CASE(3,13,23)
            DO J=1,NVELL(K)
                ISWEIR(NBVV(K,J))= .TRUE. 
            ENDDO
            CASE DEFAULT
            I=I+NVELL(K)
            END SELECT
        ENDDO
        DO IPROC=1,NPROC
            DO J=1,NNODP(IPROC)
                INDX=IMAP_NOD_LG(J,IPROC)
                IF(ISWEIR(INDX))NWEIRBNDRY(IPROC)=1
            ENDDO
        ENDDO
    ENDIF
    call memory_status()
    print *, "leaving relocalize"
    RETURN
!---------------------------------------------------------------------------
    END SUBROUTINE RELOCALIZE
!---------------------------------------------------------------------------


!---------------------------------------------------------------------------
!                S U B R O U T I N E   P R E P 1 3
!---------------------------------------------------------------------------

!     jgf46.00 This subroutine will break up the full domain nodal
!     attributes file into subdomains.

!     jgf48.47 Rewritten to only open a smaller number of subdomains
!     at a time (256 by default)

!     jgf48.50 Fixed bug (remove ALLOCATABLE from vars that don't need
!     it). Also adding following documentation:

!     BY DEFAULT, ONLY 256 SUBDOMAINS WILL BE PREPPED AT A TIME TO AVOID
!     OPENING TOO MANY FILES ON CERTAIN PLATFORMS. THIS NUMBER CAN BE
!     CONTROLLED BY THE PARAMETER MAXOPENFILES.

!---------------------------------------------------------------------------
    subroutine prep13()
!---------------------------------------------------------------------------
    use pre_global, only : useNodalAttrNames, nnodg, nnodp, nproc
    use memory_usage
    use nodalattributes
    IMPLICIT NONE

    integer :: nbytes = 0
    INTEGER :: ll             ! line loop counter
    INTEGER :: m              ! attribute default value counter
    INTEGER :: iproc          ! subdomain loop counter
    INTEGER :: sdu(nproc)     ! subdomain unit number for unit 13 files
    INTEGER :: NumNotDefault  ! number of nodes specified in the file
    CHARACTER(len=80) header   ! header comments in unit 13 files
    CHARACTER(len=80) AttrName ! label for attribute
    CHARACTER(len=80) Units    ! label for physical units
    CHARACTER(len=80) Skipped  ! data we want to skip over
    REAL(SZ) DefaultVal(12) ! default value of attribute
    INTEGER :: NoOfVals       ! at each node for an attribute
    INTEGER :: Mode           !=0 to count, =1 to write
    LOGICAL :: success        ! .TRUE. if all files opened successfully
    INTEGER, ALLOCATABLE :: SDNumND(:,:)  ! subdomain# of nodes not default
!     jgf48.47 Do the decomposition for a max of 256 subdomains at a
!     time ... some platforms/compilers limit the number of files that
!     can be open at any one time.
    INTEGER, PARAMETER :: maxOpenFiles = 256
    INTEGER :: startProc
    INTEGER :: endProc
    INTEGER :: deltaProc

!     Perform decomposition over range of subdomains.
    startProc = 1
    DO WHILE ( startProc < nproc )
        deltaProc = nproc - startProc
        IF ( deltaProc > maxOpenFiles ) deltaProc = maxOpenFiles
        endProc = startProc + deltaProc
    
    !        Open full domain and subdomain fort.13 files.
        CALL OpenPrepFiles(13, 'nodal attributes              ', &
        startProc, endProc, sdu, success)
    
        IF ( .NOT. success) THEN
            WRITE(*,*) 'WARNING: Unit 13 files not preprocessed.'
            RETURN ! note early return
        ENDIF
    
    !        Read header information from full domain unit 13 file
        READ(13,'(A80)') header
        READ(13,*) NumOfNodes     ! number of nodes according to unit 13
        READ(13,*) NAttr          ! number of attributes in the unit 13 file
    
    !        Check to make sure that the number of nodes in the nodal
    !        attributes file is the same as in the grid file (unit 14).
        IF (NumOfNodes /= NNODG) THEN
            WRITE(6,9900)
            9900 FORMAT(////,1X,'!!!!!!!!!!  FATAL ERROR  !!!!!!!!!', &
            //,1X,'The number of nodes in the grid file (unit 14) and' &
            /,1X,'the nodal attributes file (unit 13) must match.', &
            //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
            STOP                   ! We're toast.
        ENDIF
    
    !        Transcribe header information into subdomain unit 13 files
        DO iproc = startProc, endProc
            WRITE(sdu(iproc),'(A80)') header
            WRITE(sdu(iproc),*) NNODP(iproc)
            WRITE(sdu(iproc),*) NAttr
        ENDDO
    
    !        Transcribe attribute names from full domain file to subdomains.
        DO k=1, NAttr
            READ(13,'(A80)') AttrName
            READ(13,'(A80)') Units
            READ(13,*) NoOfVals
            READ(13,*) (DefaultVal(m),m=1,NoOfVals)
            DO iproc=startProc, endProc
                WRITE(sdu(IPROC),'(A80)') AttrName
                WRITE(sdu(IPROC),'(A80)') Units
                WRITE(sdu(IPROC),*) NoOfVals
                WRITE(sdu(IPROC),'(12(1x,e16.9))') &
                (DefaultVal(m),m=1,NoOfVals)
            END DO
        END DO
    
    !        Allocate and initialize the matrix for the number of Non Default
    !        nodes in each SubDomain for each nodal attribute
        ALLOCATE(SDNumND(nproc,NAttr))
        nbytes = 8*nproc*nattr
        call memory_alloc(nbytes)
        DO iproc=startProc,endProc
            DO k=1, NAttr
                SDNumND(iproc,k)=0
            END DO
        END DO
    !        We need to figure out how many nodes go into each subdomain
    !        for each attribute.
        CALL processNodalAttr(NAttr, 0, sdu, SDNumND, &
        startProc, endProc, naType)
    
    !     Now rewind and advance to the beginning of the data again
        REWIND(13)
        DO ll=1, 3
            READ(13,*) skipped    ! skip header, NumOfNodes, NAttr
        END DO
        DO k=1, NAttr
            DO ll=1, 4
                READ(13,*) skipped ! skip AttrName,Units,NoOfVals,Default
            END DO
        END DO
    
    !        Now read each of the nodal attributes and transcribe them to the
    !        appropriate subdomain.
        CALL processNodalAttr(NAttr, 1, sdu, SDNumND, &
        startProc, endProc, naType)
        DEALLOCATE(SDNumND)
        nbytes = 8*nproc*nattr
        call memory_dealloc(nbytes)
    
    !        Close full domain and subdomain files
        CLOSE (13)
        DO iproc=startProc, endProc
            CLOSE(sdu(iproc))
        ENDDO
        startProc = endProc + 1
    END DO

    if (allocated(useNodalAttrNames)) then
        DEALLOCATE(useNodalAttrNames)
        nbytes = 4*nwp
        call memory_dealloc(nbytes)
    endif

    call memory_status()
    RETURN
!---------------------------------------------------------------------------
    END SUBROUTINE PREP13
!---------------------------------------------------------------------------

!     ----------------------------------------------------------------
!        S U B R O U T I N E   P R O C E S S   N O D A L   A T T R
!     ----------------------------------------------------------------

!     jgf46.00 Subroutine to support PREP13. This subroutine is called
!     twice; once to determine the number of nodes with non-default
!     values going into each subdomain, then a second time to actually
!     place the data in the subdomain files.

!     This is necessary because the attributes in the subdomain files
!     must each have the number of non-default values at the top, and
!     this information cannot be known until we have processed the
!     entire fulldomain file.

!     ----------------------------------------------------------------
    subroutine processNodalAttr(NAttr, mode, sdu, SDNumND, &
    startProc, endProc, naType)
    use pre_global
    use sizes, only : ASCII, XDMF
    use nodalattributes, only : na
    implicit none

    INTEGER,intent(in) :: NAttr  ! number of attributes in the file
    INTEGER,intent(in) :: Mode   !=0 to count and return, =1 to write
    INTEGER,intent(in),dimension(nproc) :: sdu !i/o unit number array
    INTEGER,intent(inout),dimension(nproc,NAttr) :: SDNumND
    INTEGER,intent(in) :: startProc
    INTEGER,intent(in) :: endProc
    INTEGER,intent(in) :: naType ! ascii and xdmf are supported
    INTEGER :: NumNotDefault      ! number of nodes specified in the file
    INTEGER :: NumCol             ! number of values per node for an attr
    INTEGER :: NodeNum            ! full domain node number
    INTEGER :: SDNode             ! subdomain node number
    INTEGER :: i                  ! node loop counter
    INTEGER :: j                  ! column loop counter
    INTEGER :: k                  ! attribute loop counter
    INTEGER :: m                  ! mapping loop counter
    INTEGER :: iproc              ! subdomain loop counter
    INTEGER :: iproc2             ! mapped subdomain
    CHARACTER(len=80) AttrName ! label for attribute
    REAL(SZ), ALLOCATABLE :: AttrData(:)  ! attribute data

    CHARACTER(len=80) Skipped  ! data we want to skip over

    DO k=1, NAttr
        select case(naType)
        case(ASCII)
        read(13,*) AttrName
        case(XDMF)
        attrName = trim(adjustl(na(k)%attrName))
        case default
        write(6,'(a,i0,a)') 'ERROR: Nodal attribute file format ',naType, &
        ' is not supported by adcprep.'
        stop
        end select
        IF (Mode == 1) THEN
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
        CASE("initial_river_elevation")
        NumCol=1
        CASE &
        ("average_horizontal_eddy_viscosity_in_sea_water_wrt_depth")
        NumCol=1
        CASE &
        ("average_horizontal_eddy_diffusivity_in_sea_water_wrt_depth")
        NumCol=1
        CASE &
        ("min_and_max_primitive_weighting_in_continuity_equation")
        NumCol=2
        CASE &
        ("subtract_from_depth_to_equal_mean_sea_level")
        NumCol=1
    ! Casey 100210: Allow SWAN to handle wave refraction as a nodal attribute.
        CASE("wave_refraction_in_swan")
        NumCol=1
        CASE("elemental_slope_limiter")
        NumCol=1
    ! Corbitt 120328: Allow local advection as a nodal attribute.
        CASE("advection_state")
        NumCol=1
        CASE DEFAULT
        NumCol=1
        IF (Mode == 0) THEN
            WRITE(6,1001)    ! Nodal Attributes file
            WRITE(6,1021) AttrName ! contains invalid name
        ELSE
            WRITE(6,1031) AttrName ! Process 1st column only
        ENDIF
        END SELECT
        ALLOCATE(attrData(numCol))
        select case(naType)
        case(ASCII)
        read(13,*) numNotDefault
        case(XDMF)
        numNotDefault = na(k) % numNodesNotDefault
        case default
        write(6,'(a,i0,a)') 'ERROR: Nodal attribute file format ',naType, &
        ' is not supported by adcprep.'
        stop
        end select
        IF (Mode == 1) THEN
            DO iproc=startProc,endProc
                WRITE(sdu(iproc),*) SDNumND(iproc,k)
            END DO
        ENDIF
        DO i=1, NumNotDefault
            select case(naType)
            case(ASCII)
            READ(13,*) nodeNum, (AttrData(j),j=1,NumCol)
            case(XDMF)
            nodeNum = na(k) % nonDefaultNodes(i)
            attrData(:) = na(k) % nonDefaultVals(:,nodeNum)
            case default
            write(6,'(a,i0,a)') 'ERROR: Nodal attribute file format ',naType, &
            ' is not supported by adcprep.'
            stop
            end select
            IF (ITOTPROC(NodeNum) == 1) THEN
                iproc = IMAP_NOD_GL(1,NodeNum)
                IF ( (iproc < startProc) .OR. (iproc > endProc) ) THEN
                    CYCLE ! skip it if it does not map to our range of procs
                ENDIF
                IF (Mode == 0) SDNumND(iproc,k) = SDNumND(iproc,k)+1
                IF (Mode == 1) THEN
                    SDNode = IMAP_NOD_GL(2,NodeNum)
                    WRITE(sdu(iproc),1100) SDNode,(AttrData(j),j=1,NumCol)
                ENDIF
            ELSE
                DO m=1, ITOTPROC(NodeNum)
                    iproc2 = IMAP_NOD_GL2(2*(m-1)+1,NodeNum)
                    DO iproc=startProc, endProc
                        IF (iproc == iproc2) THEN !f.d. node maps to this s.d.
                            IF (Mode == 0) THEN
                                SDNumND(iproc,k)=SDNumND(iproc,k)+1
                            ENDIF
                            IF (Mode == 1) THEN
                                SDNode = IMAP_NOD_GL2(2*(m-1)+2,NodeNum)
                                WRITE(sdu(iproc),1100) SDNode, &
                                (AttrData(j),j=1,NumCol)
                            ENDIF
                        ENDIF
                    END DO
                END DO
            END IF
        END DO
        DEALLOCATE(AttrData)
        IF (Mode == 1) THEN
            WRITE(6,'(A25,A80)') '     Finished processing ', AttrName
            WRITE(6,*) 'for processor range ',startProc,' to ',endProc
        ENDIF
    END DO


    1001 FORMAT('ERROR: The Nodal Attributes File (unit 13)')
    1021 FORMAT('contains invalid name: ',A80)
    1031 FORMAT('WARNING: Processed only one column of unrecognized ',A80)
    1100 FORMAT(I10,32000(2X,E16.8))

    RETURN
!     ----------------------------------------------------------------
    END SUBROUTINE ProcessNodalAttr
!     ----------------------------------------------------------------


!---------------------------------------------------------------------------
    subroutine prep13XDMF()
!---------------------------------------------------------------------------
    use pre_global, only : useNodalAttrNames, nnodg, nnodp, nproc
    use memory_usage
    use nodalattributes
    use sizes, only : naFileName
    implicit none

    integer :: nbytes = 0
    integer :: iproc          ! subdomain loop counter
    integer :: sdu(nproc)     ! subdomain unit number for unit 13 files
    integer :: Mode           !=0 to count, =1 to write
    logical :: success        ! .TRUE. if all files opened successfully
    character(len=15) sdFileName     ! subdomain file name   !increased from 14 to 15 tcm v50.66.03
    integer, allocatable :: SDNumND(:,:)  ! subdomain# of nodes not default
!     jgf48.47 Do the decomposition for a max of 256 subdomains at a
!     time ... some platforms/compilers limit the number of files that
!     can be open at any one time.
    integer, parameter :: maxOpenFiles = 256
    integer :: startProc
    integer :: endProc
    integer :: deltaProc
    integer :: errorIO
    real(sz), allocatable :: diff(:) ! difference between nodal value and default value
    logical, allocatable :: areDefaultValues(:)
    integer :: nonDefaultCount

    call readNodalAttrXDMF()

! we need to compute the number of nondefault values and
! create a list of nodes with nondefault values
    do i=1,nAttr
        if (na(i)%numVals == 1) then
        ! machine precision prevents us from simply checking whether the
        ! value .ne. the default value
            diff = abs(na(i)%xdmfArray - na(i)%defaultVals(1))
            na(i)%numNodesNotDefault = count(diff > 1.e-6)
        ! now allocate space for the non default values and populate them
            allocate(na(i)%nonDefaultVals(1,na(i)%numNodesNotDefault))
            allocate(na(i)%nonDefaultNodes(na(i)%numNodesNotDefault))
        ! now record the node number and value where the values are not
        ! the default
            nonDefaultCount = 1
            do j=1,nnodg
                if (diff(j) > 1.e-6) then
                    na(i)%nonDefaultNodes(nonDefaultCount) = j
                    na(i)%nonDefaultVals(1,nonDefaultCount) = na(i)%xdmfArray(j)
                    nonDefaultCount = nonDefaultCount + 1
                endif
            end do
        else
        ! determine the number of nondefault values
            areDefaultValues = .TRUE. 
            do j=1,nnodg
                do k=1,na(i)%numVals
                    if (abs(na(i)%xdmfMatrix(k,j)-na(i)%defaultVals(k)) > 1.e-6) then
                        areDefaultValues(j) = .FALSE. 
                    endif
                enddo
            enddo
        ! now allocate space for the non default values and populate them
            na(i)%numNodesNotDefault = count(areDefaultValues.eqv. .FALSE. )
            allocate(na(i)%nonDefaultVals(na(i)%numVals,na(i)%numNodesNotDefault))
            allocate(na(i)%nonDefaultNodes(na(i)%numNodesNotDefault))
            nonDefaultCount = 1
            do j=1,nnodg
                if (areDefaultValues(j).eqv. .FALSE. ) then
                    na(i)%nonDefaultNodes(nonDefaultCount) = j
                    do k=1,na(i)%numVals
                        na(i)%nonDefaultVals(k,nonDefaultCount) = &
                        na(i)%xdmfMatrix(k,j)
                    end do
                    nonDefaultCount = nonDefaultCount + 1
                endif
            end do
        endif
    end do

!     Perform decomposition over range of subdomains.
    startProc = 1
    do while ( startProc < nproc )
        deltaProc = nproc - startProc
        if ( deltaProc > maxOpenFiles ) deltaProc = maxOpenFiles
        endProc = startProc + deltaProc
    ! Open each of the subdomain files
        do iproc = startProc, endProc
            sdu(iproc) = 105 + (iproc-1)
            sdFileName = 'PE0000/fort.13'
            call iwrite(sdFileName, 3, 6, iproc-1)
            open(unit=sdu(iproc), file=sdFileName, iostat=ErrorIO)
            Success = .TRUE. 
            IF ( ErrorIO > 0 ) THEN
                write(6,'(a,a,a)') "ERROR: Subdomain file " &
                // trim(sdFileName) // " cannot be opened."
                Success = .FALSE. 
                stop
            endif
        enddo
    
    !        Transcribe header information into subdomain unit 13 files
        do iproc = startProc, endProc
            write(sdu(iproc),'(a)') trim(adjustl(nodalAttributesComment))
            write(sdu(iproc),*) NNODP(iproc)
            write(sdu(iproc),*) NAttr
            do i=1,nAttr
                write(sdu(iproc),'(a)') trim(adjustl(na(i)%attrName))
                write(sdu(iproc),'(a)') trim(adjustl(na(i)%units))
                write(sdu(iproc),'(99(i0))') na(i)%numVals
                write(sdu(iproc),'(99(F15.7))') (na(i)%defaultVals(j), j=1,na(i)%numVals)
            end do
        end do
    
    !        Allocate and initialize the matrix for the number of Non Default
    !        nodes in each SubDomain for each nodal attribute
        ALLOCATE(SDNumND(nproc,nAttr))
        nbytes = 8*nproc*nattr
        call memory_alloc(nbytes)
        SDNumND(:,:)=0
    !        We need to figure out how many nodes go into each subdomain
    !        for each attribute.
        CALL processNodalAttr(NAttr, 0, sdu, SDNumND, &
        startProc, endProc, naType)
    
    !        Now read each of the nodal attributes and transcribe them to the
    !        appropriate subdomain.
        CALL processNodalAttr(NAttr, 1, sdu, SDNumND, &
        startProc, endProc, naType)
        DEALLOCATE(SDNumND)
        nbytes = 8*nproc*nattr
        call memory_dealloc(nbytes)
    
    !        Close subdomain files
        DO iproc=startProc, endProc
            CLOSE(sdu(iproc))
        ENDDO
        startProc = endProc + 1
    END DO

    if (allocated(useNodalAttrNames)) then
        DEALLOCATE(useNodalAttrNames)
        nbytes = 4*nwp
        call memory_dealloc(nbytes)
    endif

    call memory_status()
!---------------------------------------------------------------------------
    end subroutine prep13XDMF
!---------------------------------------------------------------------------



    SUBROUTINE PREP14()
    USE PRE_GLOBAL
    USE PREP_WEIR,ONLY:NWEIRBNDRY

!---------------------------------------------------------------------------C
!                     (  Serial Version  2/28/98  )                         C
!  This routine writes a Local Grid file "fort.14" file for each subdomain  C
!  using the domain decomposition of the ADCIRC grid created by the routine C
!  DECOMP.                                                                  C
!                                                                           C
!  The Decomposition Variables are defined in the include file adcprep.inc  C
!  This version is compatible with ADCIRC version 34.03                     C
!                                                                           C
!---------------------------------------------------------------------------C

    IMPLICIT NONE
    INTEGER :: I,I1,J,K,M,ETYPE,ITEMP,ITEMP2,ILNODE,ILNODE2,ILNODE3
    INTEGER :: JD,JG,JP,IPROC,IPROC2,IPROC3,DISC,BBN,IBP
    INTEGER :: INDX,INDEX2,ITOT,ITYPE,NUMS(10)
    CHARACTER LOCFN*14,PE*6
    CHARACTER(80) :: OUTMSG

    ETYPE = 3   ! The only Element-Type supported by ADCIRC is 3.

!--------------------------------------------------------------------------
!--MAIN LOOP:   Write a Local Grid File ( fort.14 ) for each PE
!--------------------------------------------------------------------------

    NETA_MAX = 0   ! max number of open boundary nodes on any subdomain

    DO 1000 IPROC = 1,NPROC
    
        LOCFN(1:14) = 'PE0000/fort.14'
        CALL IWRITE(LOCFN,3,6,IPROC-1)
        OPEN (14,FILE=LOCFN)
    
    !--------------------------------------------------------------------------
    !--OPEN BOUNDARY NODES PROCESSING BEGINS HERE
    !--------------------------------------------------------------------------
    
    !--Partition the open boundary nodes between various processors
    
        NETAP(IPROC) = 0
        DO K=1, NOPE
            NVDLLP(K) = 0
            DO J=1, NETA
                OBNODE_LG(J,IPROC) = 0
                NBDVP(K,J) = 0
            ENDDO
        ENDDO
    
        ITOT = 0
        DO K = 1,NOPE
            DO I = 1,NVDLL(K)
                ITOT = ITOT + 1
                INDX = NBDV(K,I)
                DO J = 1,ITOTPROC(INDX)
                    ITEMP = (J-1)*2+1
                    IPROC2 = IMAP_NOD_GL2(ITEMP,INDX)
                    ILNODE = IMAP_NOD_GL2(ITEMP+1,INDX)
                    IF (IPROC == IPROC2) THEN
                        NETAP(IPROC) = NETAP(IPROC)+1
                        NVDLLP(K) = NVDLLP(K) + 1
                        NBDVP(K,NVDLLP(K)) = ILNODE
                        OBNODE_LG(NETAP(IPROC),IPROC)=ITOT
                    ENDIF
                ENDDO
            ENDDO
        ENDDO
        IF (NETAP(IPROC) > NETA_MAX) NETA_MAX = NETAP(IPROC)
    
        NOPEP(IPROC) = 0
        DO K = 1,NOPE
            IF (NVDLLP(K) /= 0) THEN
                NOPEP(IPROC) = NOPEP(IPROC) + 1
            ENDIF
        ENDDO
    
    
    !--------------------------------------------------------------------------
    !--LAND BOUNDARY NODES PROCESSING BEGINS HERE
    !--------------------------------------------------------------------------
    
    !--Partition Land Boundary Segments between PEs
    
        NVELP(IPROC) = 0
        DO K = 1,NBOU
            NVELLP(K) = 0
            IBTYPEP(K,IPROC) = IBTYPE(K)
        ENDDO
    
        DO K = 1,NBOU
        
        !--Weir Land Boundary Node-Pair Case
        ! od vjp 3/8/99
        !  mod to allow that each of Weir-node pair might be ghosts nodes
        
            IF ((IBTYPE(K) == 4) .OR. (IBTYPE(K) == 24) .OR. &
            (IBTYPE(K) == 5) .OR. (IBTYPE(K) == 25)) THEN
                DO I = 1,NVELL(K)
                    INDX = NBVV(K,I)
                    INDEX2 = IBCONNR(K,I)
                    DO J = 1,ITOTPROC(INDX)
                        ITEMP = (J-1)*2 + 1
                        IPROC2  =  IMAP_NOD_GL2(ITEMP,INDX)
                        ILNODE2 =  IMAP_NOD_GL2(ITEMP+1,INDX)
                        IF (IPROC == IPROC2) THEN
                            DO JD = 1, ITOTPROC(INDEX2)
                                ITEMP2 = (JD-1)*2 + 1
                                IPROC3  = IMAP_NOD_GL2(ITEMP2,INDEX2)
                                ILNODE3 = IMAP_NOD_GL2(ITEMP2+1,INDEX2)
                                IF (IPROC == IPROC3) THEN
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
            
            !--All Other Land Boundary Node types
            
            ELSE
            
                DO I = 1,NVELL(K)
                    INDX = NBVV(K,I)
                    DO J = 1,ITOTPROC(INDX)
                        ITEMP = (J-1)*2 + 1
                        IPROC2 =  IMAP_NOD_GL2(ITEMP,INDX)
                        ILNODE =  IMAP_NOD_GL2(ITEMP+1,INDX)
                        IF (IPROC == IPROC2) THEN
                            NVELP(IPROC) = NVELP(IPROC) + 1
                            NVELLP(K) = NVELLP(K) + 1
                            LBINDEX_LG(K,NVELLP(K)) = I
                            NBVVP(K,NVELLP(K)) = ILNODE
                        ENDIF
                    ENDDO
                ENDDO
            
            ENDIF
        
        ENDDO
    
    ! od 05/18/2004 rl -- I don't think this next part is the correct
    !  way to handle islands.  Rather, if an island is split by a domain, it
    !  should remain an island.  This will ensure that the boundary is
    !  closed.  The only error would occur in ghost node space, which is
    !  not a problem since the answers are not used there anyway.

    ! od 12/18/98 vjp --this section re-written
    !--If a PE has only part of a closed internal land boundary
    !  modify its local IBTYPE to be an external land boundary segment
    !  of the same type by decrementing its IBTYPE.
    !  and remove a closing loop node if present

    
    !        DO K=1, NBOU
    !          IF (NVELLP(K).LT.NVELL(K)) THEN
    !            IF (  (IBTYPEP(K,IPROC).EQ.1)
    !    &         .OR.(IBTYPEP(K,IPROC).EQ.11)
    !    &         .OR.(IBTYPEP(K,IPROC).EQ.21)) THEN
    ! decrement ibtype
    !              IBTYPEP(K,IPROC) = IBTYPEP(K,IPROC)-1
    ! remove loop closing node
    !              IF (NVELLP(K).GT.1.AND.
    !    &           NBVVP(K,NVELLP(K)).EQ.NBVVP(K,1)) THEN
    !                NVELLP(K) = NVELLP(K)-1
    !              ENDIF
    !            ENDIF
    !          ENDIF
    !        ENDDO

    ! If a segment contains only one node, remove the segment from the list
    ! (NOTE: rl 5/18/04 I don't see how this could possibly happen, including
    !  ghost nodes)

        DO K=1, NBOU
            IF (NVELLP(K) == 1) NVELLP(K) = 0
        ENDDO

    
    !--Count the number of land boundary segments on PE IPROC.
    
        NBOUP(IPROC) = 0
        DO K = 1,NBOU
            IF (NVELLP(K) /= 0) THEN
                NBOUP(IPROC) = NBOUP(IPROC) + 1
            ENDIF
        ENDDO
    
    !--Count to check correctness of NVELP
    
        DISC=0  ! LB Nodes with non-zero normal discharge
        BBN=0   ! Mainland Barrier Boundary Nodes
        IBP=0   ! Internal Barrier Boundary Pairs
        ITEMP = 0
    
    !     jgf46.21 Added support for IBTYPE=52.
        DO 400 K=1,NBOU
            IF (NVELLP(K) == 0) GOTO 400
            ITYPE = IBTYPEP(K,IPROC)
        ! kmd - added for rivers in baroclinic simulation
            IF (ABS(ITYPE/100) == 1) THEN
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
                IF ((ITYPE == 1) .OR. (ITYPE == 11) .OR. &
                (ITYPE == 21)) THEN
                    IF ((I == NVELLP(K)) .AND. (NBVVP(K,I) /= I1)) THEN
                        ITEMP = ITEMP + 1
                    ENDIF
                ENDIF
                IF (I == 1) I1 = NBVVP(K,I)
            ENDDO
        400 END DO
    
    !        print *, IPROC-1,ITEMP,DISC,BBN,2*IBP
        ITEMP  = ITEMP + DISC + BBN + 2*IBP
        IF (ITEMP /= NVELP(IPROC)) THEN
        !          print *, "changed value from ",NVELP(IPROC)," to ",ITEMP
            NVELP(IPROC) = ITEMP
        ENDIF
        IF (NVELP(IPROC)+1 > MNVEL) THEN
            print *, "NVEL exceeds parameter value MNVEL on PE",IPROC
            print *, "local NVEL value = ",ITEMP
            stop
        ENDIF
    
    !--Construct a LBCODE for each Land Boundary Node of this PE
    
        JP=0
        DO K = 1,NBOU
            DO I=1, NVELLP(K)
                JP = JP+1
                LBCODEP(JP,IPROC) = IBTYPEP(K,IPROC)
            ENDDO
        ENDDO
    
    !--Determine whether there are any normal flow boundaries local to PE
    
    ! kmd - changed for rivers in baroclinic simulations
        NFLUXFP(IPROC) = 0
        DO K=1, NBOU
            IF (NVELLP(K) > 0) THEN
                ITYPE=IBTYPE(K)
                IF (ABS(ITYPE/100) == 1) THEN
                    ITYPE = (ABS(ITYPE)-100)*(ITYPE/ABS(ITYPE))
                END IF
                IF ((ITYPE == 2) .OR. (ITYPE == 12) &
                 .OR. (ITYPE == 32) &
                 .OR. (ITYPE == 22) .OR. (ITYPE == 52)) THEN
                    NFLUXFP(IPROC) = 1
                ENDIF
            ENDIF
        ENDDO
    
    !--------------------------------------------------------------------------
    !--BEGIN WRITING LOCAL GRID ( fort.14 ) FILE HERE
    !--------------------------------------------------------------------------
    
    !--Write Mesh Data
    
        WRITE(14,80) AGRID
    
        NUMS(1) = NELP(IPROC)
        NUMS(2) = NNODP(IPROC)
    
    ! jgf45.06    CALL INSERT(SIZEMSG,OUTMSG,NUMS,2)
    ! jgf45.06    WRITE(14,80) OUTMSG
        WRITE(14,43) NELP(IPROC),NNODP(IPROC) !jgf45.06
    
        DO J = 1,NNODP(IPROC)
            INDX = IMAP_NOD_LG(J,IPROC)
            WRITE(14,44) J,X(INDX),Y(INDX),DP(INDX)
        ENDDO
    
        DO J = 1,NELP(IPROC)
            WRITE(14,45) J,ETYPE,NNEP(1,J,IPROC),NNEP(2,J,IPROC), &
            NNEP(3,J,IPROC)
        ENDDO
        43 FORMAT(2I8)
        44 FORMAT(I8,3(E24.12))
        45 FORMAT(5I8)
    
    !--Write Open Boundary Data
    
        CALL NEWINDEX(NOPEMSG,OUTMSG,NOPEP(IPROC))
        WRITE(14,80) OUTMSG
    
        CALL NEWINDEX(NETAMSG,OUTMSG,NETAP(IPROC))
        WRITE(14,80) OUTMSG
    
        ITOT = 0
        DO K = 1,NOPE
            IF (NVDLLP(K) > 0)THEN
                ITOT = ITOT + 1
            ! Casey 090304: Added the following section.  If we are coupling to SWAN,
            !             then we also want to give the global number of each
            !             boundary segment.
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
    
    !--Write Land Boundary Data
    
        CALL NEWINDEX(NBOUMSG,OUTMSG,NBOUP(IPROC))
        WRITE(14,80) OUTMSG
    
        CALL NEWINDEX(NVELMSG,OUTMSG,NVELP(IPROC))
        WRITE(14,80) OUTMSG
    
        DO K = 1,NBOU
            IF(NVELLP(K) > 0)THEN
                ITYPE = IBTYPEP(K,IPROC)
                NUMS(1) = NVELLP(K)
                NUMS(2) = ITYPE
            ! Casey 090304: Added the following section.  If we are coupling to SWAN,
            !             then we also want to give the global number of each
            !             boundary segment.
#ifndef ADCSWAN
                CALL INSERT(NVELLMSG(K),OUTMSG,NUMS,2)
#else
                NUMS(3) = NOPE + K
                CALL INSERT(NVELLMSG(K),OUTMSG,NUMS,3)
#endif
                WRITE(14,80) OUTMSG
            
                IF ((ITYPE /= 3) .AND. (ITYPE /= 13) .AND. &
                (ITYPE /= 23) .AND. (ITYPE /= 4) .AND. &
                (ITYPE /= 24) .AND. (ITYPE /= 5) .AND. &
                (ITYPE /= 25)) THEN
                    DO I = 1,NVELLP(K)
                        WRITE(14,'(I8)') NBVVP(K,I)
                    ENDDO
                ELSEIF ((ITYPE == 3) .OR. (ITYPE == 13) .OR. &
                    (ITYPE == 23)) THEN
                    IF(USE_TVW)NWEIRBNDRY(IPROC) = 1
                    DO I = 1,NVELLP(K)
                        INDX = LBINDEX_LG(K,I)
                        WRITE(14,81) NBVVP(K,I),BAR1(K,INDX),BAR2(K,INDX)
                    ENDDO
                
                ELSEIF ((ITYPE == 4) .OR. (ITYPE == 24)) THEN
                    IF(USE_TVW)NWEIRBNDRY(IPROC) = 1
                    DO I = 1,NVELLP(K)
                        INDX = LBINDEX_LG(K,I)
                        WRITE(14,82) NBVVP(K,I),IBCONNRP(K,I), &
                        BAR1(K,INDX),BAR2(K,INDX),BAR3(K,INDX)
                    ENDDO
                ELSEIF ((ITYPE == 5) .OR. (ITYPE == 25)) THEN
                    DO I = 1,NVELLP(K)
                        INDX = LBINDEX_LG(K,I)
                        WRITE(14,83) NBVVP(K,I),IBCONNRP(K,I), &
                        BAR1(K,INDX),BAR2(K,INDX),BAR3(K,INDX), &
                        BAR4(K,INDX),BAR5(K,INDX),BAR6(K,INDX)
                    ENDDO
                ENDIF
            ENDIF
        ENDDO
    
        CLOSE(14)
    
    1000 END DO

!--Print Summary of Boundary Node Decomposition

    print *, " "
    print *, "Boundary Node Decomposition Data"
    print *, "DOMAIN      NOPE    NETA    NBOU  NVEL    NWEIR"
    WRITE(*,90)  "GLOBAL",NOPE, NETA, NBOU, NVEL, NWEIR
    DO IPROC=1, NPROC
        PE(1:6) = 'PE0000'
        CALL IWRITE(PE,3,6,IPROC-1)
        WRITE(*,90)  PE,NOPEP(IPROC),NETAP(IPROC), &
        NBOUP(IPROC),NVELP(IPROC),NWEIRP(IPROC)
    ENDDO

    80 FORMAT(A80)
    81 FORMAT(I8,2X,E13.6,2X,E13.6)
    82 FORMAT(I8,2X,I8,2X,E13.6,2X,E13.6,2X,E13.6)
    83 FORMAT(I8,2X,I8,6(2X,E13.6))
    90 FORMAT(1X,A6,5I8)

    RETURN
    END SUBROUTINE PREP14

!---------------------------------------------------------------------------C
!                     (  Serial Version  2/28/98  )                         C
!  This routine writes a Local Input file "fort.15" file for each subdomain C
!  using the domain decomposition of the ADCIRC grid created by the routine C
!  DECOMP.                                                                  C
!                                                                           C
!  The Decomposition Variables are defined in the include file adcprep.inc  C
!  This version is compatible with ADCIRC version 34.03                     C
!                                                                           C
!           Modifications by RL on 10/9/01 to accomodate NWS = -2           C
!---------------------------------------------------------------------------C
    SUBROUTINE PREP15()
    USE PRE_GLOBAL
    use memory_usage
    USE PREP_WEIR
    USE HARM, ONLY : NAMEFR
    use subprep, only : subdomainOn, found_sm_nml ! NCSU Subdomain Modeling
    use nodalattributes, only : outputTau0
    IMPLICIT NONE
    integer :: nbytes = 0
    INTEGER :: I,J,K,M,JG,JP,KK, ios_stations
    INTEGER :: INDX,ITOT,ILNODE,IPROC,IPROC2,ITYPE,NUMS(10)
    CHARACTER LOCFN*80,PE*6,LOCSTATFN*20
    INTEGER,ALLOCATABLE :: NTIMEVARYINGWEIRP(:)
    CHARACTER(80) :: OUTMSG

!--Write a Local Input file ( fort.15 ) for each PE

! max number of stations in any subdomain
    NSTAE_MAX = 0; NSTAV_MAX = 0; NSTAM_MAX = 0; NSTAC_MAX = 0


    ALLOCATE(NTIMEVARYINGWEIRP(NPROC))
    NTIMEVARYINGWEIRP(:)=0

    DO 1000 IPROC = 1,NPROC
    

        LOCFN = 'PE0000/fort.15'
        CALL IWRITE(LOCFN,3,6,IPROC-1)
        OPEN (15,FILE=TRIM(LOCFN))
    
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
        IF (NWP > 0) THEN !jgf46.00 write nodal attributes
            DO I=1, NWP
                WRITE(15,80) useNodalAttrNames(I)
            ENDDO
        ENDIF
        WRITE(15,80) NCORMSG
        WRITE(15,80) NTIPMSG
    !     jgfdebug46.02 Added check for NWS=45 to write NWS=5
        IF (NWS == 45) THEN
            WRITE(15,'(A1)') "5"
        ELSE
            WRITE(15,80) NWSMSG
        ENDIF
        WRITE(15,80) RAMPMSG
        WRITE(15,80) GMSG
        WRITE(15,80) TAU0MSG

    !        jgf47.11 Added writing of limits for time varying tau0
        IF ( (TAU0 <= -5.d0) .AND. (TAU0 > -6.d0) ) THEN
            WRITE(15,80) TAU0LIMMSG
        ENDIF

        WRITE(15,80) DTMSG
        WRITE(15,80) STATMSG
        WRITE(15,80) REFTMSG

    !       tcm v49.64.01 No changes needed here for the use of ICE
        IF((NWS == 0) .AND. (NRS >= 1)) WRITE(15,80) RSTIMMSG  ! sb46.28sb03
        IF((NWS == 1) .AND. (NRS >= 1)) WRITE(15,80) RSTIMMSG  ! sb46.28sb03
    !     jgfdebug46.02 Added check for NWS=45.
    !     jgf46.02 Added NWS=8.
    !     jgf46.16 Merged:
    !     rjw added NWS=19: asymmetric hurricane wind model v2.0
    !     jie added NWS=20: generalized asymmetric vortex model
    !     sb46.28sb01 added NWS=12: OWI format
    !     jgf50.38.05: added NWS=15: HWind format
    !     tcm v51.06.02 added NWS=16: GFDL Met Data
        IF ((ABS(NWS) == 2) .OR. (ABS(NWS) == 4) .OR. (ABS(NWS) == 45) .OR. &
        (ABS(NWS) == 5) .OR. (ABS(NWS) == 6) .OR. (ABS(NWS) == 8) &
         .OR. (ABS(NWS) == 12) .OR. (ABS(NWS) == 15) &
         .OR. (ABS(NWS) == 16) &
         .OR. (ABS(NWS) == 19) .OR. (ABS(NWS) == 29) .OR. &
        (ABS(NWS) == 20))THEN
            WRITE(15,80) WSMSG1
        ENDIF
        IF (NWS == 3) THEN
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
    
    !--If there were any normal flow boundaries local to PE, process them
    
    !         PRINT *, NFFRMSG
    !         PRINT *, "NFLUXFP(",IPROC,") = ", NFLUXFP(IPROC)

        IF (NFLUXFP(IPROC) == 1) THEN
        
            NFLBNP = 0
            DO I=1, NFLBN
                INDX = FLBN(I)
                DO J=1, ITOTPROC(INDX)
                    IPROC2 = IMAP_NOD_GL2(2*(J-1)+1,INDX)
                    IF (IPROC == IPROC2) THEN
                        NFLBNP = NFLBNP + 1
                        FLBNXP(NFLBNP) = FLBNX(I)
                    ENDIF
                ENDDO
            ENDDO
        
            WRITE(15,80) NFFRMSG
            IF ((NFFR /= 0) .AND. (NFFR /= -1)) THEN
                DO I=1,NFFR
                    WRITE(15,80) FBOUNTAG(I)
                    WRITE(15,80) FREQMSG(I)
                ENDDO
                DO I=1,NFFR
                    WRITE(15,80) ALPHA2(I)
                    DO J=1,NFLBNP
                        WRITE(15,80) QNMSG(I,FLBNXP(J))
                    ! ebug               print *, "PE=",IPROC," FLUXNODE=",FLBNXP(J)
                    ENDDO
                ENDDO
            ENDIF
        
        ENDIF

    ! bell
    !--IF THERE ARE INTERNAL/EXTERNAL TIME VARYING FLUX BOUNDARIES, PARSE TO PEs
    
        IF(USE_TVW)THEN
            IF(ALLOCATED(NWEIRBNDRY) .AND. NTIMEVARYINGWEIR > 0)THEN
                IF(NWEIRBNDRY(IPROC) == 1)THEN
                !...If this processor has Time varying weirs, count them
                    NTIMEVARYINGWEIRP(IPROC) = 0
                    DO I=1,NNODP(IPROC)
                        INDX = IMAP_NOD_LG(I,IPROC)
                        DO J = 1,NTIMEVARYINGWEIR
                            IF(INDX == NODES_TVW(J))THEN
                                NTIMEVARYINGWEIRP(IPROC) = &
                                NTIMEVARYINGWEIRP(IPROC) + 1
                            ENDIF
                        ENDDO
                    ENDDO
                ENDIF
            !...Open the local PE fort.tvw
                LOCFN = 'PE0000/'//TRIM(tvw_file)
                CALL IWRITE(LOCFN,3,6,IPROC-1)
                OPEN(FILE=TRIM(LOCFN),UNIT=98,ACTION="WRITE")
                IF(NTIMEVARYINGWEIRP(IPROC) > 0)THEN
                !...Write fort.tvw namelist file for the
                !   local PE
                    WRITE(98,*) NTIMEVARYINGWEIRP(IPROC)
                    DO I=1,NNODP(IPROC)
                        INDX = IMAP_NOD_LG(I,IPROC)
                        DO J = 1,NTIMEVARYINGWEIR
                            IF(INDX == NODES_TVW(J))THEN
                                WRITE(98,'(A)') &
                                TRIM(ADJUSTL(TIMEVARYINGWEIRMSSG(J)))
                            ENDIF
                        ENDDO
                    ENDDO
                ELSE
                !...Write an empty fort.tvw to avoid old files from
                !   previous decompositions
                    WRITE(98,*) 0
                ENDIF
                CLOSE(98)
            ENDIF
        ENDIF
    
    !--Write Local Elevation Station Info:
    !--Create Local-to-Global element "ownership" of an elevation station
    
    !     WRITE(15,80) STAEMSG !jgf45.07 we may have changed NOUTE in adcprep
        WRITE(15,*) NOUTE,TOUTSE,TOUTFE,NSPOOLE
    
        NSTAEP(IPROC) = 0
        DO K = 1,abs(NSTAE)                          !tcm -- added the comments below
            DO J=1,NELP(IPROC)                        !nelp(iproc) lists the number of elements from processor iproc
                INDX = abs(IMAP_EL_LG(J,IPROC))        ! global element number
                IF (INDX == NNSEG(K)) THEN             !nnseg(k) contains the element number station k resides in
                    NSTAEP(IPROC) = NSTAEP(IPROC) + 1
                    KK = K
                    if (STAE_SHARE(K) > -1) KK = -K
                    IMAP_STAE_LG(NSTAEP(IPROC),IPROC) = KK
                    STAE_SHARE(K) = IPROC
                ! tcm v51.20.03 once found exit the element loop
                    exit
                ENDIF
            ENDDO
        ENDDO
        NSTAE_MAX = MAX(NSTAEP(IPROC),NSTAE_MAX)
    
    !...     update the number of stations for this proc's domain

        if (use_elev_stat_file ) then  !tcm v51.20.03
            CALL INSERT(NSTAEMSG,OUTMSG,-NSTAEP(IPROC),1)  !keep the negative sign for fort.15
            write(15,80) OUTMSG
            CALL INSERT(NSTAEMSG,OUTMSG,NSTAEP(IPROC),1)
            LOCSTATFN(1:20) = 'PE0000/elev_stat.151'
            CALL IWRITE(LOCstatFN,3,6,IPROC-1)
            ios_stations = 0
            open(unit=151,file=locstatfn, &
            status='unknown',iostat=ios_stations)
            write(151,80) OUTMSG
        else
            CALL INSERT(NSTAEMSG,OUTMSG,NSTAEP(IPROC),1)
            write(15,80) OUTMSG
        endif
    
    !...     write the stations located in this proc's domain
        DO K=1,NSTAEP(IPROC)
            INDX = abs(IMAP_STAE_LG(K,IPROC))
            if (use_elev_stat_file) then
                write(151,80) STAELOC(INDX)
            else
                WRITE(15,80) STAELOC(INDX)
            endif
        ENDDO
        if (use_elev_stat_file ) close(151)
    
    !--Write Local Velocity Station Info:
    !--Create Local-to-Global element "ownership" of an velocity station
    
        WRITE(15,*) NOUTV,TOUTSV,TOUTFV,NSPOOLV
    
        NSTAVP(IPROC) = 0
        DO K = 1,abs(NSTAV)
            DO J=1,NELP(IPROC)
                INDX = abs(IMAP_EL_LG(J,IPROC))
                IF (INDX == NNSVG(K)) THEN
                    NSTAVP(IPROC) = NSTAVP(IPROC) + 1
                    KK = K
                    if (STAV_SHARE(K) > -1) KK = -K
                    IMAP_STAV_LG(NSTAVP(IPROC),IPROC) = KK
                    STAV_SHARE(K) = IPROC
                ! tcm v51.20.03 once found exit the element loop
                    exit
                ENDIF
            ENDDO
        ENDDO
        NSTAV_MAX = MAX(NSTAVP(IPROC),NSTAV_MAX)
    
        if (use_vel_stat_file ) then  !tcm v51.20.03
            CALL INSERT(NSTAVMSG,OUTMSG,-NSTAVP(IPROC),1)  !keep the negative sign for fort.15
            write(15,80) OUTMSG
            CALL INSERT(NSTAVMSG,OUTMSG,NSTAVP(IPROC),1)
            LOCSTATFN(1:19) = 'PE0000/vel_stat.151'
            CALL IWRITE(LOCstatFN,3,6,IPROC-1)
            ios_stations = 0
            open(unit=151,file=locstatfn(1:19), &
            status='unknown',iostat=ios_stations)
            write(151,80) OUTMSG
        else
            CALL INSERT(NSTAVMSG,OUTMSG,NSTAVP(IPROC),1)
            write(15,80) OUTMSG
        endif

    
        DO K=1,NSTAVP(IPROC)
            INDX = abs(IMAP_STAV_LG(K,IPROC))
            if (use_vel_stat_file ) then
                WRITE(151,80) STAVLOC(INDX)
            else
                WRITE(15,80) STAVLOC(INDX)
            endif
        ENDDO
        if (use_vel_stat_file  ) close(151)
    
    !--If IM=10 Write Concentration Station Info:
    !--Create Local-to-Global element "ownership" of an concentration station
    
        NSTACP(IPROC) = 0
        IF (C2D_PTrans .OR. C3D_PTrans) THEN !jgf46.28
        
        !     WRITE(15,80) STACMSG   !jgf45.07 we may have changed NOUTC in adcprep
            WRITE(15,*) NOUTC,TOUTSC,TOUTFC,NSPOOLC
        
            DO K = 1,abs(NSTAC)
                DO J=1,NELP(IPROC)
                    INDX = abs(IMAP_EL_LG(J,IPROC))
                    IF (INDX == NNSCG(K)) THEN
                        NSTACP(IPROC) = NSTACP(IPROC) + 1
                        KK = K
                        if (STAC_SHARE(K) > -1) KK = -K
                        IMAP_STAC_LG(NSTACP(IPROC),IPROC) = KK
                        STAC_SHARE(K) = IPROC
                    ! tcm v51.20.03 once found exit the element loop
                        exit
                    ENDIF
                ENDDO
            ENDDO
            NSTAC_MAX = MAX(NSTACP(IPROC),NSTAC_MAX)
        
        !...     update the number of stations for this proc's domain
            if (use_conc_stat_file ) then  !tcm v51.20.03
                CALL INSERT(NSTACMSG,OUTMSG,-NSTACP(IPROC),1)  !keep the negative sign for fort.15
                write(15,80) OUTMSG
                CALL INSERT(NSTACMSG,OUTMSG,NSTACP(IPROC),1)
                LOCSTATFN(1:20) = 'PE0000/conc_stat.151'
                CALL IWRITE(LOCstatFN,3,6,IPROC-1)
                ios_stations = 0
                open(unit=151,file=locstatfn, &
                status='unknown',iostat=ios_stations)
                write(151,80) OUTMSG
            else
                CALL INSERT(NSTACMSG,OUTMSG,NSTACP(IPROC),1)
                write(15,80) OUTMSG
            endif
        ! ... write the stations located in this proc's domain
            DO K=1,NSTACP(IPROC)
                INDX = abs(IMAP_STAC_LG(K,IPROC))
                IF (use_conc_stat_file) then
                    write(151,80) STACLOC(INDX)
                ELSE
                    WRITE(15,80) STACLOC(INDX)
                ENDIF
            ENDDO
            IF (use_conc_stat_file) close(151)
        
        ENDIF
    
    !--Write Local Meterological Station Info:
    !--Create Local-to-Global element "ownership" of an elevation station
    
        NSTAMP(IPROC) = 0
        IF (NWS /= 0) THEN
            WRITE(15,*) NOUTM,TOUTSM,TOUTFM,NSPOOLM
            DO K = 1,abs(NSTAM)
                DO J=1,NELP(IPROC)
                    INDX = abs(IMAP_EL_LG(J,IPROC))
                    IF (INDX == NNSMG(K)) THEN
                        NSTAMP(IPROC) = NSTAMP(IPROC) + 1
                        KK = K
                        if (STAM_SHARE(K) > -1) KK = -K
                        IMAP_STAM_LG(NSTAMP(IPROC),IPROC) = KK
                        STAM_SHARE(K) = IPROC
                    ! tcm v51.20.03 once found exit the element loop
                        exit
                    ENDIF
                ENDDO
            ENDDO
            NSTAM_MAX = MAX(NSTAMP(IPROC),NSTAM_MAX)
        
        !...        update the number of stations for this proc's domain
            if (use_met_stat_file ) then  !tcm v51.20.03
                CALL INSERT(NSTAMMSG,OUTMSG,-NSTAMP(IPROC),1)  !keep the negative sign for fort.15
                write(15,80) OUTMSG
                CALL INSERT(NSTAMMSG,OUTMSG,NSTAMP(IPROC),1)
                LOCSTATFN(1:19) = 'PE0000/met_stat.151'
                CALL IWRITE(LOCstatFN,3,6,IPROC-1)
                ios_stations = 0
                open(unit=151,file=locstatfn(1:19), &
                status='unknown',iostat=ios_stations)
                write(151,80) OUTMSG
            else
                CALL INSERT(NSTAMMSG,OUTMSG,NSTAMP(IPROC),1)
                write(15,80) OUTMSG
            endif
        
        !...        write the stations located in this proc's domain
            DO K=1,NSTAMP(IPROC)
                INDX = abs(IMAP_STAM_LG(K,IPROC))
                IF (use_met_stat_file) then
                    WRITE(151,80) STAMLOC(INDX)
                ELSE
                    WRITE(15,80) STAMLOC(INDX)
                ENDIF
            ENDDO
            IF (use_met_stat_file) CLOSE(151)
        ENDIF

    
    !--Write Local Elevation Data Output Info
    
    !      WRITE(15,80) OUTGEMSG !jgf45.07 we may have changed NOUTGE in adcprep
        WRITE(15,*) NOUTGE,TOUTSGE,TOUTFGE,NSPOOLGE
    
    !--Write Local Velocity Data Output Info
    
    !     WRITE(15,80) OUTGVMSG !jgf45.07 we may have changed NOUTGV in adcprep
        WRITE(15,*) NOUTGV,TOUTSGV,TOUTFGV,NSPOOLGV
    
    !     jgf45.07 write subdomain concentration data output info if necessary
    
        IF (IM == 10) WRITE(15,*) NOUTGC,TOUTSGC,TOUTFGC,NSPOOLGC
    
    !--Write Local Wind Velocity Data Output Info ( added 4/16/98 vjp )
    
    !     jgf45.07 we may have changed NOUTGW in adcprep
    !     IF (NWS.NE.0) WRITE(15,80) OUTGWMSG
        IF (NWS /= 0) WRITE(15,*) NOUTGW,TOUTSGW,TOUTFGW,NSPOOLGW
    
    !--Write Harmonic Analysis Data
    
        WRITE(15,80) HARFRMSG
        DO I=1,NHARFR
            WRITE(15,'(A10)') NAMEFR(I)
            WRITE(15,80) HAFREMSG(I)
        !           WRITE(15,*) HAFREQ(I),HAFF(I),HAFACE(I)
        ENDDO
    
        WRITE(15,80) HARPARMSG
        WRITE(15,80) OUTHARMSG
    
    !--Write Hot Start Info
    
        WRITE(15,80) HSTARMSG
    
    !--Write Solver Info
    
        WRITE(15,80) SOLVMSG
    
    !--Write 3DVS Info
    
        IF(C3DVS) THEN
            CALL PREP15_3DVS(IPROC)
        !        ELSEIF(C3DDSS) THEN
        !           CALL PREP15_3DDSS(IPROC)
        ENDIF
    
    !     jgf48.03 Write netCDF metadata, if necessary
        IF (useNetCDF.eqv. .TRUE. ) THEN
        ! rite(*,*) 'writing netcdf metadata to fort.15' ! jgfdebug
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
        ! rite(*,*) 'not writing netcdf metadata' ! jgfdebug
        ENDIF

    !...     tcm v50.66.02 additions for time varying bathymetry
        IF (FOUND_TBC_NML) then  !If there was a namelist in the original fort.15 put it in the decomp 15's
        !         IF (NDDT.NE.0) THEN
            write(15,*) '! -- Begin Time Varying Bathymetry Inputs --'
            write(15,TimeBathyControl)
            write(15,*) '! -- End Time Varying Bathymetry Inputs --'
        ENDIF
    
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
    
        if (found_metCon_nml) then  !metControl namelist was found so write it in the parsed files
            write(15,*) '! -- Begin Met Control Namelist --'
        !         write(15,*) "&metControl WindDragLimit=",WindDragLimit,
        !     &         " DragLawString='",DragLawString,"' rhoAir=",rhoAir," /"

            write(15,*) "&metControl "
            write(15,*) "    WindDragLimit=",Winddraglimit,","
            write(15,*) "    DragLawString='",Draglawstring,"',"
            write(15,*) "    rhoAir=",rhoAir,","
            write(15,*) "    invertedBarometerOnElevationBoundary=", &
            invertedBarometerOnElevationBoundary
            write(15,*) "/"
            write(15,*) '! -- End Met Control Namelist --'
        endif

        if (found_tvw_nml) then
            write(15,*) '! -- Begin TVW Control Namelist --'
            write(15,*) "&TVWControl "
            write(15,*) "    USE_TVW=",USE_TVW,","
            write(15,*) "    TVW_FILE='",TRIM(ADJUSTL(TVW_FILE)),"',"
            write(15,*) "    NOUT_TVW=",NOUT_TVW,","
            write(15,*) "    TOUTS_TVW=",TOUTS_TVW,","
            write(15,*) "    TOUTF_TVW=",TOUTF_TVW,","
            write(15,*) "    NSPOOL_TVW=",NSPOOL_TVW
            write(15,*) "/"
            write(15,*) '! -- End TVW Control Namelist --'
        endif

    ! Casey 121019: Added multiplication factor to be used before sending winds to coupled wave models.
        IF(FOUND_WC_NML)THEN
            WRITE(15,*) '! -- Begin Wave Coupling Namelist --'
            WRITE(15,waveCoupling)
            WRITE(15,*) '! -- End Wave Coupling Namelist --'
        ENDIF

    ! NCSU Subdomain Modeling
        if (FOUND_SM_NML) then
            WRITE(15,*) "&subdomainModeling subdomainOn=",subdomainOn," /"
        endif

    ! wet dry control
        if (foundWetDryControlNameList) then
            WRITE(15,'(a,l,a,l,a,l,a,l)') &
            "&wetDryControl outputNodeCode=",outputNodeCode, &
            " outputNOFF=",outputNOFF," noffActive=",noffActive," /"
        endif

    ! jgf52.08.02: inundation output control
        if (foundInundationOutputControlNamelist) then
            write(15,'(a,l,a,f15.8,a)') &
            '&inundationOutputControl inundationOutput=', &
            inundationOutput,' inunThresh=',inunThresh," /"
        endif

        CLOSE(15)
    
    1000 END DO

    IF(C3DVS .AND. (IEVC == 0)) THEN
        DEALLOCATE ( EVTot )
        nbytes = 8*nfen
        call memory_dealloc(nbytes)
    ENDIF


!--Print Summary of Stations

    print *, " "
    print *, "Station Data"
    print *, "DOMAIN      NSTAE   NSTAV    NSTAC    NSTAM"
    WRITE(*,92)  "GLOBAL",abs(NSTAE),abs(NSTAV), &
    abs(NSTAC),abs(NSTAM)
    DO IPROC=1, NPROC
        PE(1:6) = 'PE0000'
        CALL IWRITE(PE,3,6,IPROC-1)
        WRITE(*,92)  PE,NSTAEP(IPROC),NSTAVP(IPROC), &
        NSTACP(IPROC),NSTAMP(IPROC)
    ENDDO

    RETURN
    80 FORMAT(A80)
    92 FORMAT(1X,A6,4I8)
    END SUBROUTINE PREP15


!---------------------------------------------------------------------------C
!                     (  Serial Version  6/24/02  )                         C
!  This routine writes the 3DVS info in the Local Input file "fort.15" file C
!  for each subdomain using the domain decomposition of the ADCIRC grid     C
!  created by the routine DECOMP.                                           C
!                                                                           C
!  The Decomposition Variables are defined in the include file adcprep.inc  C
!  This version is compatible with ADCIRC version 41.11a                    C

!     jgf45.11 Updated to handle new format of 3D input files with stations
!     defined by coordinates rather than node numbers.
!---------------------------------------------------------------------------C
    SUBROUTINE PREP15_3DVS(IPROC)
    USE PRE_GLOBAL
    USE GLOBAL_3DVS, ONLY : SIGMA
    IMPLICIT NONE
    INTEGER :: N          ! vertical grid layer counter
    INTEGER :: IPROC      ! subdomain counter
    INTEGER :: SDStation  ! subdomain station
    INTEGER :: FDStation  ! full domain station
    INTEGER :: SDEle      ! subdomain station element
    INTEGER :: FDEle      ! full domain element

!     jgf45.10 removed IDIAG
    WRITE(15,80) IDENMSG
    WRITE(15,80) SLIPMSG
    WRITE(15,80) Z0MSG
    WRITE(15,80) ALPMSG
    WRITE(15,80) FEMSG
!     jgf45.12 Added code to record thicknesses of vertical grid layers,
!     if necessary.
    IF(IGC == 0) THEN
        DO N=1,NFEN
            WRITE(15,*) Sigma(N)
        ENDDO
    ENDIF
    WRITE(15,80) EVCMSG
!     jgf45.12 Add code to record vertical eddy viscosity profile.
    IF(IEVC == 0) THEN
        DO N=1,NFEN
            WRITE(15,*) EVTot(N)
        ENDDO
    ENDIF
    IF((IEVC == 50) .OR. (IEVC == 51)) WRITE(15,80) THETAMSG
!     -------------------------------------------------------------
!     jgf45.11 Create mapping from full domain 3D density station
!     elements to corresponding elements in subdomains. Write out
!     subdomain station locations to fort.15 file.
!     -------------------------------------------------------------
    WRITE(15,*) I3DSD,TO3DSDS,TO3DSDF,NSPO3DSD
!   kmd48.33bc changed
    IF(NSTA3DD /= 0) THEN
        NNSTA3DDP(IPROC) = 0
        DO FDStation = 1, NSTA3DD
            DO SDEle = 1, NELP(IPROC)
                FDEle = abs(IMAP_EL_LG(SDEle,IPROC))
                IF ( FDEle == NNS3DDG(FDStation) ) THEN
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
!     -------------------------------------------------------------
!     jgf45.11 Create mapping from full domain 3D velocity station
!     elements to corresponding elements in subdomains. Write out
!     velocity subdomain station locations to fort.15 file.
!     -------------------------------------------------------------
    WRITE(15,*) I3DSV,TO3DSVS,TO3DSVF,NSPO3DSV
!   kmd48.33bc changed
    IF(NSTA3DV /= 0) THEN
        NNSTA3DVP(IPROC) = 0
        DO FDStation = 1, NSTA3DV
            DO SDEle = 1, NELP(IPROC)
                FDEle = abs(IMAP_EL_LG(SDEle,IPROC))
                IF ( FDEle == NNS3DVG(FDStation) ) THEN
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
!     -------------------------------------------------------------
!     jgf45.11 Create mapping from full domain 3D turbulence station
!     elements to corresponding elements in subdomains. Write out
!     turbulence subdomain station locations to fort.15 file.
!     -------------------------------------------------------------
    WRITE(15,*) I3DST,TO3DSTS,TO3DSTF,NSPO3DST
!   kmd48.33bc changed
    IF(NSTA3DT /= 0) THEN
        NNSTA3DTP(IPROC) = 0
        DO FDStation = 1, NSTA3DT
            DO SDEle = 1, NELP(IPROC)
                FDEle = abs(IMAP_EL_LG(SDEle,IPROC))
                IF ( FDEle == NNS3DTG(FDStation) ) THEN
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

!    kmd48.33bc add 3D boundary condition information
    IF (CBAROCLINIC) THEN
        WRITE(15,80) RESBCFLAGMSG
        IF (RES_BC_FLAG /= 0) THEN
            IF (NOPEP(IPROC) > 0) THEN
                WRITE(15,80) BCTIMEMSG
                WRITE(15,80) BCSTATMSG
            END IF
            IF (BCFLAG_TEMP /= 0) THEN
                WRITE(15,80) TBCTIMEMSG
            END IF
        END IF
    END IF

    IF (CBAROCLINIC) THEN
        WRITE(15,80) SPONGEDISTMSG
        WRITE(15,80) EqnstateMSG
    END IF


!     jgf45.12: Write out the parameters for the transport equation, if
!     necessary.
    IF (C3D_BTrans) THEN
    !     Lateral and vertical diffusion coefficients.
        WRITE(15,*) NLSD, NVSD
        WRITE(15,*) NLTD, NVTD
    !     Time stepping coefficient for the transport equation terms.
        WRITE(15,*) ALP4
    !   kmd48.33 took out as it is no longer needed with new heat flux boundary conditions
    !     Temperature boundary condition file type, if necessary
    !         IF ( IDEN .eq. 3 .or. IDEN .eq. 4 ) THEN
    !            WRITE(15,*) NTF
    !         ENDIF
    ENDIF

    RETURN
    80 FORMAT(A80)
    81 FORMAT(I8,2E15.8,2I8,A32)
    82 FORMAT(500I8)
    END SUBROUTINE PREP15_3DVS
!-----------------------------------------------------------------------
!     End of subroutine PREP15_3DVS
!-----------------------------------------------------------------------


    SUBROUTINE PREP18()
    USE PRE_GLOBAL
    use memory_usage

!---------------------------------------------------------------------------C
!                     (  Serial Version  6/10/2011  )                       C
!  This Routine writes a message-passing file "fort.18" for each subdomain  C
!  of the domain decomposition created by DECOMP.                           C
!                                                                           C
!  The Decomposition Variables are defined in the include file adcprep.inc  C
!  This version is compatible with ADCIRC version 50.21                     C
!                                                                           C
!  tcm V50.21 -- Changed all I8 formats to I12                              C
!---------------------------------------------------------------------------C

    IMPLICIT NONE
    integer :: nbytes = 0
    INTEGER :: N1, N2, N3, KMIN
    INTEGER :: I,J,K,M,ITEMP,IPR,IPR1
    INTEGER :: INDX,ITOT,IEL,IELG,ILNODE,IPROC,ITYPE
    INTEGER,ALLOCATABLE :: RES_NODE(:)
    CHARACTER LOCFN*14,PE*6

!  Allocate local arrays

    ALLOCATE ( RES_NODE(MNPP) )
    nbytes = 4*mnpp
    call memory_alloc(nbytes)

!--Write Message-Passing File for each PE

    DO 1000 I = 1,NPROC
    
        LOCFN(1:14) = 'PE0000/fort.18'
        CALL IWRITE(LOCFN,3,6,I-1)
        OPEN (18,FILE=LOCFN)

        write(18,3050) FileFmtVersion, 0, 0

    ! jp 9/17/06
    !--Write the Global indexes of all local elements in local element order

    ! Casey 100209: Changed I8 to I12.
        WRITE(18,3000) NELG, MNEP, NELP(I)  ! number of Global elements
        DO J = 1,NELP(I)
            INDX = IMAP_EL_LG(J,I)
            WRITE(18,'(I12)') INDX        ! Global index of local element
        ENDDO

    !--Write the Global indexes of all local nodes in local node order
    !  write global index as positive if a resident node and negative
    !  if a ghost node

        WRITE(18,3001) NNODG, MNPP, NNODP(I)   ! number of Global nodes
        ITOT = 0
        DO J = 1,NNODP(I)
            INDX = IMAP_NOD_LG(J,I)
            IPR = IMAP_NOD_GL(1,INDX)
            IF (IPR == I)THEN
                ITOT = ITOT + 1
                RES_NODE(ITOT) = J
                WRITE(18,'(I12)') INDX        ! Global index of resident node
            ELSE
                WRITE(18,'(I12)') -1*INDX     ! Global index of ghost node
            ENDIF
        ENDDO
        IF (ITOT /= NOD_RES_TOT(I)) STOP 'ERROR IN# OF RES. NODES'

    !--Write local normal flow boundary flag
    !--vjp This info is used only for relocalizing fort.15
        WRITE(18,3002) NFLUXFP(I) ! normal flow b.c. flag for subdomain

    !--Write global and local total number of elevation boundary nodes
    !--vjp This info is used only for relocalizing fort.15
        WRITE(18,3003) NETA, NETA_MAX, NETAP(I) ! number of global elevation b.c. nodes
        DO J = 1,NETAP(I)
            INDX = OBNODE_LG(J,I)
            WRITE(18,'(I12)') INDX           ! Global open boundary node index
        ENDDO

    !--Write the Global indexes of all Elevation Stations in local node order
    !  write global index as positive if a resident node and negative
    !  if a ghost node

        WRITE(18,3004) abs(NSTAE), NSTAE_MAX, NSTAEP(I) ! number of Global Elevation Stations
        DO J = 1,NSTAEP(I)
            INDX = IMAP_STAE_LG(J,I)
            WRITE(18,'(I12)') INDX           ! Global station number
        ENDDO

    !--Write the Global indexes of all Velocity Stations in local node order
    !  write global index as positive if a resident node and negative
    !  if a ghost node

        WRITE(18,3005) abs(NSTAV), NSTAV_MAX, NSTAVP(I) ! number of Global Velocity Stations
        DO J = 1,NSTAVP(I)
            INDX = IMAP_STAV_LG(J,I)
            WRITE(18,'(I12)') INDX           ! Global station number
        ENDDO

    !--Write the Global indexes of all Elevation Stations in local node order
    !  write global index as positive if a resident node and negative
    !  if a ghost node

        WRITE(18,3006) abs(NSTAM), NSTAM_MAX, NSTAMP(I) ! number of Global Meteorlogical Stations
        DO J = 1,NSTAMP(I)
            INDX = IMAP_STAM_LG(J,I)
            WRITE(18,'(I12)') INDX           ! Global station number
        ENDDO

    !--Write the Global indexes of all Concentration Stations in local node order
    !  write global index as positive if a resident node and negative
    !  if a ghost node

        WRITE(18,3007) abs(NSTAC), NSTAC_MAX, NSTACP(I) ! number of Global Concentration Stations
        DO J = 1,NSTACP(I)
            INDX = IMAP_STAC_LG(J,I)
            WRITE(18,'(I12)') INDX           ! Global station number
        ENDDO

    !---------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------
    
    !--Write the Resident Node List
    
        WRITE(18,3010) (I-1),NOD_RES_TOT(I)
        WRITE(18,1130) (RES_NODE(J),J=1,ITOT)
    
    !--Write the Number of Communicating PEs
    
        WRITE(18,3020) NUM_COMM_PE(I)
    
    !--Write the Receive List
    
    ! loop over this subdomain's neighbors
        DO J = 1,NUM_COMM_PE(I)
        ! get the subdomain number of the Jth neighbor
        ! of this subdomain
            IPR = COMM_PE_NUM(J,I)
        ! zero out the total number of ghost nodes on this
        ! subdomain that are residents on the Jth neighbor
        ! of this subdomain
            IRECV_TOT(J,I) = 0
        ! loop over nodes in this subdomain
            DO K = 1,NNODP(I)
            ! get corresponding fulldomain node number
                INDX = IMAP_NOD_LG(K,I)
            ! if the fulldomain node is a ghost and is a resident
            ! of the Jth neighbor of this subdomain
                IF (IMAP_NOD_GL(1,INDX) == IPR) THEN
                ! increment the total number of ghost nodes that
                ! are residents of the Jth neighbor subdomain
                    IRECV_TOT(J,I) = IRECV_TOT(J,I) + 1
                ! record the local node number of this subdomain
                ! that will receive data
                    IRECV(IRECV_TOT(J,I)) = K
                ! uncomment next line and comment preceding line for debugging
                !                 IRECV(IRECV_TOT(J,I)) = INDX
                ENDIF
            ENDDO
            WRITE(18,3030) (IPR-1), IRECV_TOT(J,I)
            WRITE(18,1130) (IRECV(K),K=1,IRECV_TOT(J,I))
        ENDDO
    
    !--write the send list
    
    ! loop over this subdomain's neighbors
        DO J = 1,NUM_COMM_PE(I)
        ! get the subdomain number of the Jth neighbor
        ! of this subdomain
            IPR = COMM_PE_NUM(J,I)
            ISEND_TOT(J,I) = 0
        ! loop over nodes in the Jth neighbor subdomain
            DO K = 1,NNODP(IPR)
            ! get the fulldomain node number in the Jth neighbor subdomain
                INDX = IMAP_NOD_LG(K,IPR)
            ! if the node is a resident of this subdomain
                IF (IMAP_NOD_GL(1,INDX) == I) THEN
                ! increment the total number of nodes on this subdomain
                ! that are ghosts on the Jth neighbor
                    ISEND_TOT(J,I) = ISEND_TOT(J,I) + 1
                ! record the local node number on this subdomain
                ! that will send data
                    ISEND(ISEND_TOT(J,I)) = IMAP_NOD_GL(2,INDX)
                ! uncomment next line and comment preceding line for debugging
                !                 ISEND(ISEND_TOT(J,I)) = INDX
                ENDIF
            ENDDO
            WRITE(18,3040)  IPR-1, ISEND_TOT(J,I)
            WRITE(18,1130) (ISEND(K),K=1,ISEND_TOT(J,I))
        ENDDO
    
        IF (C3D.eqv. .TRUE. ) THEN
        !           jgf49.43.18: Add 3D station mappings from subdomain to fulldomain
        !           to accomodate globalio.
        
        !           Write the fulldomain station numbers of all 3D density stations
        !           in local node order; write the fulldomain station number as positive
        !           for resident stations and negative for ghost stations.
            WRITE(18,3060) NSTA3DD, MAXVAL(NNSTA3DDP), NNSTA3DDP(I)
            DO J=1,NNSTA3DDP(I)
                WRITE(18,1131) IMAP_STA3DD_LG(J,I)
            END DO
        !           3D velocity stations
            WRITE(18,3061) NSTA3DV, MAXVAL(NNSTA3DVP), NNSTA3DVP(I)
            DO J=1,NNSTA3DVP(I)
                WRITE(18,'(I12)') IMAP_STA3DV_LG(J,I)
            END DO
        !           3D turbulence stations
            WRITE(18,3062) NSTA3DT, MAXVAL(NNSTA3DTP), NNSTA3DTP(I)
            DO J=1,NNSTA3DTP(I)
                WRITE(18,'(I12)') IMAP_STA3DT_LG(J,I)
            ENDDO
        ENDIF
    
        CLOSE(18)
    
    1000 END DO

!--Compute the surface to volume ratio (in %)

    DO I = 1,NPROC
        ITOT = 0
        DO J = 1,NUM_COMM_PE(I)
            ITOT = ITOT + IRECV_TOT(J,I)
        ENDDO
        PROC_SV(I) = (ITOT/REAL(NOD_RES_TOT(I)))*100.0
    !        WRITE(6,*) I-1,PROC_SV(I)
    ENDDO

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

    92 FORMAT(1X,A6,2X,I7,2X,F8.2)
    1130 FORMAT(8X,6I12) !(8X,9I8)
    1131 FORMAT(:,I12)
! Casey 100209: Changed I8 to I12 through this section.
    3000 FORMAT('NELG    ',3I12)
    3001 FORMAT('NNODG   ',3I12)
    3002 FORMAT('NFLUXF  ',I12)
    3003 FORMAT('NETA    ',3I12)
    3004 FORMAT('NSTAE   ',3I12)
    3005 FORMAT('NSTAV   ',3I12)
    3006 FORMAT('NSTAM   ',3I12)
    3007 FORMAT('NSTAC   ',3I12)
    3010 FORMAT('RES NODE',2I12)
    3020 FORMAT('COMM PE ',2I12)
    3030 FORMAT('RECV PE ',2I12)
    3040 FORMAT('SEND PE ',2I12)
    3050 FORMAT('FileFmt ',3I12)
    3060 FORMAT('NSTA3DD ',3I12)
    3061 FORMAT('NSTA3DV ',3I12)
    3062 FORMAT('NSTA3DT ',3I12)
    END SUBROUTINE PREP18


    SUBROUTINE PREP19()
    USE PRE_GLOBAL
    use memory_usage

!---------------------------------------------------------------------------C
!                     (  Serial Version  2/28/98  )                         C
!  This routine writes a Local "Aperiodic Elevation Boundary Condtions"     C
!  (fort.19) file for each subdomain using the domain decomposition of      C
!  the ADCIRC grid created by the routine DECOMP.                           C
!                                                                           C
!  The Decomposition Variables are defined in the include file adcprep.inc  C
!  This version is compatible with ADCIRC version 34.03                     C

!     jgf45.12 Added subroutine call to open files.
!                                                                           C
!---------------------------------------------------------------------------C

    IMPLICIT NONE
    integer :: nbytes = 0
    INTEGER :: I,J,IPROC
    INTEGER :: SDU(NPROC) ! subdomain unit numbers
    LOGICAL :: Success    ! .TRUE. if files opened without errors
    CHARACTER(40) ::  ETIMINC,ESBINP
    CHARACTER*40,ALLOCATABLE :: ESBIN(:)


!--Enter, Locate, Open, and Read the ADCIRC UNIT 19
!  Global Aperiodic Elevation Boundary Conditions file

!     Open full domain and subdomain fort.19 files
    CALL OpenPrepFiles(19, 'aperiodic elevation boundary  ', &
    &      1, nproc, SDU, Success)
    IF ( .NOT. Success) THEN
        WRITE(*,*) 'WARNING: Unit 19 files not preprocessed.'
        RETURN ! note early return
    ENDIF

!--Allocate local arrays

    ALLOCATE ( ESBIN(MNETA) )
    nbytes = 8*mneta
    call memory_alloc(nbytes)

    READ(19,40) ETIMINC
    DO IPROC = 1,NPROC
        WRITE(SDU(IPROC),40)  ETIMINC
    ENDDO

!--While ( NOT EOF ) Read NETA BCs from Global File

    1000 CONTINUE
    DO I=1, NETA
        READ(19,40,END=9999)  ESBIN(I)
    ENDDO

    DO IPROC= 1,NPROC
        DO I=1, NETAP(IPROC)
            ESBINP = ESBIN(OBNODE_LG(I,IPROC))
            WRITE(SDU(IPROC),40) ESBINP
        ENDDO
    ENDDO

    GO TO 1000

!--Close Global file and all the Local Files

    9999 CLOSE (19)
    DO IPROC=1, NPROC
        CLOSE (SDU(IPROC))
    ENDDO

    deallocate( esbin )
    nbytes = 8*mneta
    call memory_dealloc(nbytes)
    call memory_status()
    RETURN
    40 FORMAT(A40)
    END SUBROUTINE PREP19


!---------------------------------------------------------------------------
!                S U B R O U T I N E   P R E P 2 0
!---------------------------------------------------------------------------

!     jgf45.12 This subroutine will break up the full domain aperiodic
!     flux boundaries into subdomains using the domain decomposition of
!     the ADCIRC grid created by the routine DECOMP.

!     -Written by MEB 04/01/04
!     -Added by jgf to 45.06 10/07/2005
!     -jgf45.12 Rewritten to correct bugs in subdomain fort.20
!     formatting as well as the erroneous use of the GL mapping instead
!     of GL2. Also added subroutine call to open files.

!---------------------------------------------------------------------------
    SUBROUTINE PREP20()
!---------------------------------------------------------------------------
    USE PRE_GLOBAL
    use memory_usage
    IMPLICIT NONE
    integer :: nbytes = 0
    INTEGER :: IPROC
    INTEGER :: INDEX14, I
    REAL(SZ)  FLUX_INC, FLUX_VAL
    INTEGER :: SDU(NPROC)  ! subdomain unit numbers
    LOGICAL :: Success     ! .TRUE. if all files open without error
    INTEGER :: INDX ! full domain node number for a flow boundary node
    INTEGER :: J     ! counter for subdomains that corrsp. to a single f.d. node
    INTEGER :: IPROC2! PE of a subdomain that matches a single full domain node

!     Open full domain and subdomain fort.20 files
    CALL OpenPrepFiles(20, 'aperiodic flux boundary       ', &
    &      1, nproc, SDU, Success)
    IF ( .NOT. Success) THEN
        WRITE(*,*) 'WARNING: Unit 20 files not preprocessed.'
        RETURN ! note early return
    ENDIF

!     Write Increment into all flux files

    READ(20,*) FLUX_INC
    DO IPROC=1,NPROC
        WRITE(SDU(IPROC),*) FLUX_INC
    ENDDO

!     jgf45.12 Write each full domain nodal flux value into each of the
!     subdomains that that full domain node maps to. The full domain
!     node may map to more than one subdomain node if it falls on a
!     boundary between subdomains (ghost nodes).

    33 DO I=1, EXIST_FLUX            ! loop through full domain flow nodes
        READ(20,*,END=40) FLUX_VAL ! get a flo val for this f.d. flow node
        INDX = FLUX14_ARY(I)      ! get full domain flow node number
        DO J=1, ITOTPROC(INDX)    ! loop over subdomains for 1 f.d. node
            IPROC2 = IMAP_NOD_GL2(2*(J-1)+1,INDX) ! find next subdomain
            DO IPROC=1, NPROC
                IF (IPROC == IPROC2) THEN ! full domain node maps to this s.d.
                    WRITE(SDU(IPROC),50) FLUX_VAL
                ENDIF
            END DO
        END DO
    END DO
    GOTO 33
    40 CLOSE (20)
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
    50 FORMAT (F16.8,1x,I6,1x,I6,1x,I6)
!----------------------------------------------------------------------------
    END SUBROUTINE PREP20
!----------------------------------------------------------------------------

!---------------------------------------------------------------------------
!                S U B R O U T I N E   P R E P 8 8
!---------------------------------------------------------------------------

!     kmd49 This subroutine will break up the full domain elevation
!     changes due to the river boundary information being above mean
!     sea level. It writes the fort.88 file into each subdomain using
!     the ADCIRC grid created by the routine DECOMP.

!     - added as part of Evan's changes for rivers above MSL.

!     TCM v51.24 -- Added the decomposition for a max of 256 subdomains
!     at a time ... some platforms/compilers limit the number of files that
!     can be open at any one time.
!     TCM v51.27 -- Commented out, as the fort.88 river elevation has
!         now been made a nodal attritube
!---------------------------------------------------------------------------

!      SUBROUTINE PREP88()
!      USE PRE_GLOBAL
!      use memory_usage
!C
!      IMPLICIT NONE
!      integer :: nbytes = 0
!      INTEGER I,J,IPROC
!      INTEGER SDU(NPROC) ! subdomain unit numbers
!      LOGICAL Success    ! .true. if files opened without errors
!      INTEGER :: NODP
!      CHARACTER*80 :: et_tempsWSE
!      CHARACTER*80,ALLOCATABLE :: et_SWSE(:)
!      INTEGER, PARAMETER :: maxOpenFiles = 256
!      INTEGER startProc
!      INTEGER endProc
!      INTEGER deltaProc
!C
!C     Perform decomposition over range of subdomains.
!      startProc = 1
!      DO WHILE ( startProc .lt. nproc )
!         deltaProc = nproc - startProc
!         IF ( deltaProc .gt. maxOpenFiles ) deltaProc = maxOpenFiles
!         endProc = startProc + deltaProc

!C        Open full domain and subdomain fort.88 files.
!         CALL OpenPrepFiles(88, '  river elevation data  ',
!     &     startProc, endProc, SDU, Success)

!         IF (.not.Success) THEN
!            WRITE(*,*) 'WARNING: Unit 88 files not preprocessed.'
!         RETURN ! note early return
!      ENDIF

!      ALLOCATE(et_SWSE(NNODG))
!      DO I=1, NNODG
!            READ(88,80,END=9999) et_SWSE(I)
!      END DO

!         DO IPROC=startProc,endProc
!            DO I=1, NNODP(IPROC)
!               NODP=IMAP_NOD_LG(I,IPROC)
!               et_tempsWSE=et_SWSE(NODP)
!               WRITE(SDU(IPROC),80) et_tempsWSE
!            END DO
!         END DO

!C        Close full domain and subdomain files
!         CLOSE (88)
!         DO iproc=startProc, endProc
!            CLOSE(sdu(iproc))
!         ENDDO
!         startProc = endProc + 1

!         DEALLOCATE(et_SWSE)

!         WRITE(6,'(A25,A80)') '     Finished processing ',
!     &                        'river elevation data'
!         WRITE(6,*) 'for processor range ',startProc,' to ',endProc

!      END DO  !Loop over Procs
! 80   FORMAT(A80)
! 999  CLOSE(88)
!      END SUBROUTINE PREP88
!C  End SUBROUTINE PREP88


!---------------------------------------------------------------------------
!                S U B R O U T I N E   P R E P 2 2
!---------------------------------------------------------------------------

!                     (  Serial Version  2/28/98  )                         C
!  This routine reads a global external meteorology file when NWS=1,+-2,    C
!  +-4,+-5.  In each case it wites a local meteorology file of the same     C
!  format for each subdomain using the domain decomposition of the ADCIRC   C
!  grid created by the routine DECOMP.                                      C
!                                                                           C
!  The Decomposition Variables are defined in the include file adcprep.inc  C
!  This version is compatible with ADCIRC version 43.03                     C

!     jgf46.02 Added subroutine call to open prep files; this provides
!     the user with the ability to skip the prepping of wind data files.

!     jgfdebug46.02 Added NWS=45 to imitate the behavior of the v42 (IPET)
!     code.

!     jgf46.02 Added NWS=8 to copy the wind files for the Holland model
!     into the subdomains.

!     tcm_v49.04 Removed NWS=3 and NWS=6 to correspond with the use of a
!      global file rather than local.

!---------------------------------------------------------------------------
    SUBROUTINE PREP22()
!---------------------------------------------------------------------------
    USE PRE_GLOBAL
    use memory_usage
    IMPLICIT NONE
    integer :: nbytes = 0
    LOGICAL :: FOUND,DONE
    INTEGER :: I,J,IPROC,IPROC2,ILNODE,INDX,NHG,LINDEX
    CHARACTER(80) :: PBLJAGF
!      CHARACTER FNAME*60,LOCFN*14,CMD1*63,CMD2*7,CMD*70,INLINE*80
    CHARACTER FNAME*60,CMD1*63,CMD2*7,CMD*70
    CHARACTER(170) :: Line ! line of data from NWS=8 (Holland) file
    CHARACTER(270) :: Line19 ! line of data from NWS=19 (AsymmHollandv2.0) file
    INTEGER :: SDU(NPROC)  ! subdomain unit numbers
    LOGICAL :: Success     ! .TRUE. if all files open without error
    INTEGER,ALLOCATABLE  :: NG(:)
    REAL(SZ),ALLOCATABLE :: WVNXG(:),WVNYG(:),PRG(:)
    REAL(SZ),ALLOCATABLE :: WVNXL(:),WVNYL(:),PRL(:)
    REAL(SZ) U,V,PR
    REAL(SZ) RHOWAT     !jgfdebug46.02
!     jgf48.47 Do the decomposition for a max of 256 subdomains at a
!     time ... some platforms/compilers limit the number of files that
!     can be open at any one time.
    INTEGER, PARAMETER :: maxOpenFiles = 256
    INTEGER :: startProc
    INTEGER :: endProc
    INTEGER :: deltaProc

!     Allocate local work arrays

    ALLOCATE ( NG(MNWP) )
    nbytes = 4*mnwp
    call memory_alloc(nbytes)
    ALLOCATE ( WVNXG(MNWP),WVNYG(MNWP),PRG(MNWP) )
    nbytes = 24*mnwp
    call memory_alloc(nbytes)
    ALLOCATE ( WVNXL(MNWP),WVNYL(MNWP),PRL(MNWP) )
    nbytes = 24*mnwp
    call memory_alloc(nbytes)

!     Perform decomposition over a range of subdomains.
    startProc = 1
    DO WHILE ( startProc < nproc )
        deltaProc = nproc - startProc
        IF ( deltaProc > maxOpenFiles ) deltaProc = maxOpenFiles
        endProc = startProc + deltaProc

    !        Open full domain and all subdomain fort.22 files
        CALL OpenPrepFiles(22, 'wind information              ', &
        startProc, endProc, sdu, success)

        IF ( .NOT. success) THEN
            WRITE(*,*) 'WARNING: Unit 22 files not preprocessed.'
            RETURN ! note early return
        ENDIF
    
    !--Branch to Appropriate Code
    
        SELECT CASE(ABS(NWS))
    !        -------------
        CASE(1,2,5,7)
    !        -------------
    
    !     MAIN LOOP FOR NWS = 1, +-2,+-5,+-7
    !     (1)  Read a record from Global Wind Stress File
    !     (2)  Use Decomp arrarys to Localize record to a subdomain
    !     (3)  Write Local Wind Stress record in same format
        DO                     ! loop forever (or until file ends)
            READ(22,*,END=9999) &
            (NG(I),WVNXG(I),WVNYG(I),PRG(I),I=1,NNODG)
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
    
    !        -------
        CASE(4)
    !        -------
    !        MAIN LOOP FOR NWS = +- 4  ( PBL Format )
    !        (1)  Read a record from Global Wind Stress File
    !        (2)  Use Decomp arrarys to Localize record to a subdomain
    !        (3)  Write out in PBL Format on subdomain
    
    !--Read a wind field record from the global input file
    
        DO
        READ(22,'(A80)',END=9999) PBLJAGF
        IF(PBLJAGF(2:2) == '#') THEN
            DO IPROC =  STARTPROC,ENDPROC
                WRITE(SDU(IPROC),1101)
                WRITE(SDU(IPROC),1100) 1,0.0,0.0,0.0 !victor didn't like this line 27/11/03
            ENDDO
        ELSE
        !     vjp 27/11/03
        !     rewrote this section to handle ghost-nodes
        !              READ(PBLJAGF,'(I8,3E13.5)',END=9999) NHG,U,V,PR
            READ(PBLJAGF,*,END=9999) NHG,U,V,PR
            DO J=1, ITOTPROC(NHG)
                IPROC  = IMAP_NOD_GL2(2*(J-1)+1,NHG)
                LINDEX = IMAP_NOD_GL2(2*(J-1)+2,NHG)
                WRITE(SDU(IPROC),1100) LINDEX,U,V,PR
            ENDDO
        ENDIF
    END DO

!        --------
    CASE(45)
!        --------
!        jgf46.02 Convert NWS=4 winds to NWS=5 winds to imitate the Katrina
!        (IPET) version of the code.

!-- Read a wind field record from the global input file

    DO
    RHOWAT=1000.0d0
    CALL NWS4GET(WVNXG,WVNYG,PRG,G,RHOWAT,NNODG,DONE)

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
!--   If reached EOF in NWS4GET last time go close files and return

    IF (DONE) GOTO 9999
    ENDDO

!        ------------
    CASE DEFAULT
!        ------------
    print *, "NWS=",NWS," has incorrect value in PREP22"
    RETURN

    END SELECT


!--Close Global file and all the Local Files

    9999 CLOSE (22)
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
    60 FORMAT(A60)
    170 FORMAT(A170)
    270 FORMAT(A270)
    1010 FORMAT(' File ',A60,/,' WAS NOT FOUND!  Try again',/)
    1011 FORMAT(' File ',A60,/,' WAS FOUND!  Opening & Processing file',/)
    1100 FORMAT(I8,3E13.5)
    1101 FORMAT('#')
!----------------------------------------------------------------------------
    END SUBROUTINE PREP22
!----------------------------------------------------------------------------



    SUBROUTINE PREP23()
    USE PRE_GLOBAL

!---------------------------------------------------------------------------C
!                           (  add MEB 03/04/03  )                          C
!  This routine writes a Local Input file "fort.23" file for each subdomain C
!  using the domain decomposition of the ADCIRC grid created by the routine C
!  DECOMP.                                                                  C
!                                                                           C
!  The Decomposition Variables are defined in the include file adcprep.inc  C
!  This version is compatible with ADCIRC version 34.03                     C
!                                                                           C
!---------------------------------------------------------------------------C

    IMPLICIT NONE
    INTEGER :: IPROC, NHG, J, LINDEX
    CHARACTER(80) :: PBLJAGF
    INTEGER :: SDU(NPROC)  ! subdomain unit numbers
    LOGICAL :: Success     ! .TRUE. if all files open without error
    REAL(SZ)                U,V

!--Open Global Wave Stress File ( UNIT 23 )

!     Open full domain and subdomain fort.23 files
    CALL OpenPrepFiles(23, 'wave stress                   ', &
    &      1, nproc, SDU, Success)
    IF ( .NOT. Success) THEN
        WRITE(*,*) 'WARNING: Unit 23 files not preprocessed.'
        RETURN ! note early return
    ENDIF
!--------------------------------------------------------------------------
!--MAIN LOOP
!   (1)  Read a record from Global Wave Stress File
!   (2)  Use Decomp arrays to Localize record to a subdomain
!   (3)  Write Local Wave Stress record in standard PBL format
!--------------------------------------------------------------------------

!--Read a wave field record from the global input file
!--and write out to respective local fort.23 file.

    170 READ(23,'(A80)',END=9999) PBLJAGF
    IF(PBLJAGF(2:2) == '#') THEN
        DO IPROC = 1,NPROC
            WRITE(SDU(IPROC),1101)
            WRITE(SDU(IPROC),1100) 1,0.0,0.0 !victor didn't like this line 27/11/03
        ENDDO
    ELSE
    ! vjp 27/11/03
    ! rewrote this section to handle ghost-nodes
    ! and changed if test from "and" to "or"
        READ(PBLJAGF,'(I8,2E13.5)',END=9999) NHG,U,V
        IF ((U /= 0.) .OR. (V /= 0.)) THEN
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
    1101 FORMAT ('#')

    99 RETURN
    END SUBROUTINE PREP23


!---------------------------------------------------------------------------
!                S U B R O U T I N E   P R E P 141
!---------------------------------------------------------------------------

!                     (  Serial Version  4/13/12  )                         C
!  This routine reads a global external bathymetry file when NDDT=+-1,+-2.  C
!  In each case it wites a local bathymetry file of the same                C
!  format for each subdomain using the domain decomposition of the ADCIRC   C
!  grid created by the routine DECOMP.                                      C
!                                                                           C
!  The Decomposition Variables are defined in the include file adcprep.inc  C
!  This version is compatible with ADCIRC version 50.66                     C
!                                                                           C
!  TCM -v 50.66.03 Addition for time varying Bathymetry                     C
!          This routine adopted/modified from the prep22 subroutine.        C
!                                                                           C
!---------------------------------------------------------------------------
    SUBROUTINE PREP141()
!---------------------------------------------------------------------------
    USE PRE_GLOBAL
    use memory_usage
    IMPLICIT NONE
    integer :: nbytes = 0
    LOGICAL :: FOUND,DONE
    INTEGER :: I,J,IPROC,IPROC2,ILNODE,INDX,NHG,LINDEX
    CHARACTER(80) :: PBLJAGF
    INTEGER :: SDU(NPROC)  ! subdomain unit numbers
    LOGICAL :: Success     ! .TRUE. if all files open without error
    INTEGER,ALLOCATABLE  :: NG(:)
    REAL(SZ),ALLOCATABLE :: DPG(:)  !global array
    REAL(SZ),ALLOCATABLE :: DPL(:)  !local array
    REAL(SZ) DPTMP
!     jgf48.47 Do the decomposition for a max of 256 subdomains at a
!     time ... some platforms/compilers limit the number of files that
!     can be open at any one time.
    INTEGER, PARAMETER :: maxOpenFiles = 256
    INTEGER :: startProc
    INTEGER :: endProc
    INTEGER :: deltaProc

!     Allocate local work arrays

    ALLOCATE ( NG(MNP) )
    nbytes = 4*mnp
    call memory_alloc(nbytes)
    ALLOCATE ( DPG(MNP) )   !global
    nbytes = 8*mnp
    call memory_alloc(nbytes)
    ALLOCATE ( DPL(MNP) )   !local
    nbytes = 8*mnp
    call memory_alloc(nbytes)

!     Perform decomposition over a range of subdomains.
    startProc = 1
    DO WHILE ( startProc < nproc )
        deltaProc = nproc - startProc
        IF ( deltaProc > maxOpenFiles ) deltaProc = maxOpenFiles
        endProc = startProc + deltaProc

    !        Open full domain and all subdomain fort.141 files
        CALL OpenPrepFiles(141, 'bathymetry information        ', &
        startProc, endProc, sdu, success)

        IF ( .NOT. success) THEN
            WRITE(*,*) 'WARNING: Unit 141 files not preprocessed.'
            RETURN ! note early return
        ENDIF
    
    !--Branch to Appropriate Code
    
        SELECT CASE(ABS(NDDT))
    !        -------------
        CASE(1)
    !        -------------
    
    !     MAIN LOOP FOR NWS = +-1
    !     (1)  Read a record from Global Bathymetry File
    !     (2)  Use Decomp arrarys to Localize record to a subdomain
    !     (3)  Write Local Bathymetry record in same format
        DO                     ! loop forever (or until file ends)
            READ(141,*,END=9999) &
            (NG(I),DPG(I),I=1,NNODG)
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
    
    !        -------
        CASE(2)
    !        -------
    !        MAIN LOOP FOR NWS = +- 2  ( PBL Format )
    !        (1)  Read a record from Global Bathymetry File
    !        (2)  Use Decomp arrarys to Localize record to a subdomain
    !        (3)  Write out in PBL Format on subdomain
    
    !--Read a bathymetry field record from the global input file
    !--- during the decomp phase, after each time record indicator is written (#)
    !--- we write a single entry (1,-99999.d0) to ensure that there will be no
    !--- empty records.  When this file is read by ADCIRC using nddt2get, the
    !--- extra entry (1,-99999.d0) will be ignored, and if node 1 actually is
    !--- changed then it will be read regardless if it appears twice.
    
        DO
            PBLJAGF(:) = ' '
            READ(141,'(A80)',END=9999) PBLJAGF
            IF(PBLJAGF(2:2) == '#') THEN
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
                    IF ( (IPROC >= STARTPROC) .AND. &
                         (IPROC <= ENDPROC) ) THEN
                        WRITE(SDU(IPROC),1100) LINDEX,DPTMP
                    ENDIF
                ENDDO
            ENDIF
        ENDDO


!        ------------
    CASE DEFAULT
!        ------------
        WRITE(*,*) "NDDT = ",NDDT," has incorrect value in PREP141"
        RETURN

    END SELECT


!--Close Global file and all the Local Files
9999    CLOSE (141)
        DO IPROC=STARTPROC, ENDPROC
            CLOSE (SDU(IPROC))
        ENDDO
        write(*,*) "     Finished processing fort.141 file"
        write(*,*) "for processor range ",startproc," to ",endproc
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
! 70  FORMAT(A270)
! 010 FORMAT(' File ',A60,/,' WAS NOT FOUND!  Try again',/)
! 011 FORMAT(' File ',A60,/,' WAS FOUND!  Opening & Processing file',/)
    1100 FORMAT(I8,E13.5)
    1101 FORMAT('#')
!----------------------------------------------------------------------------
    END SUBROUTINE PREP141
!----------------------------------------------------------------------------








!   kmd48.33bc add in prep subroutines for 3D boundary condition files
    SUBROUTINE PREP35()
    USE PRE_GLOBAL
    use memory_usage

!---------------------------------------------------------------------------C
!                                                                           C
!  This routine writes a Local "Residual Boundary Condtions Baroclinic"     C
!  (fort.35) file for each subdomain using the domain decomposition of     C
!  the ADCIRC grid created by the routine DECOMP.                           C
!                                                                           C
!                   Added by Kendra Dresback (Aug. 18, 2007)                C
!---------------------------------------------------------------------------C

    IMPLICIT NONE
    INTEGER :: nbytes = 0
    INTEGER :: I,J,IPROC
    INTEGER :: SDU(NPROC) ! subdomain unit numbers
    LOGICAL :: Success    ! .TRUE. if files opened without errors
    CHARACTER(40) ::  ETIMINC,RESBCBINP,GRIDINC
    CHARACTER*40,ALLOCATABLE :: RESBCBIN(:)


!--Enter, Locate, Open, and Read the ADCIRC UNIT 35
!  Global Level of No Motion Boundary Conditions file for baroclinic

!     Open full domain and subdomain fort.35 files
!      Print *, "Made it to prepping the files"
    CALL OpenPrepFiles(35, 'level of no motion boundary  ', &
    &      1, nproc, SDU, Success)
!      Print *, "Made it out of prepping the files"
    IF ( .NOT. Success) THEN
        WRITE(*,*) 'WARNING: Unit 35 files not preprocessed.'
        RETURN ! note early return
    ENDIF

!--Allocate local arrays

    ALLOCATE ( RESBCBIN(MNETA) )
    nbytes = 8*mneta
    call memory_alloc(nbytes)

!--While ( NOT EOF ) Read NETA BCs from Global File

!      PRINT *, "Made it to the reading in of the 35 file"
    DO  ! loop until end of file
        READ(35,40,END=9999) ETIMINC
        DO IPROC = 1,NPROC
            WRITE(SDU(IPROC),40)  ETIMINC
        ENDDO
        DO I=1, NETA
            READ(35,40,END=9999)  RESBCBIN(I)
        ENDDO
    
        DO IPROC= 1,NPROC
            DO I=1, NETAP(IPROC)
                RESBCBINP = RESBCBIN(OBNODE_LG(I,IPROC))
                WRITE(SDU(IPROC),40) RESBCBINP
            ENDDO
        ENDDO
    END DO


!--Close Global file and all the Local Files

    9999 CLOSE (35)
    DO IPROC=1, NPROC
        CLOSE (SDU(IPROC))
    ENDDO

    deallocate(resbcbin)
    nbytes = 8*mneta
    call memory_dealloc(nbytes)
    call memory_status()

    40 FORMAT(A40)

    RETURN
    END SUBROUTINE PREP35

    SUBROUTINE PREP36()
    USE PRE_GLOBAL
    use memory_usage

!---------------------------------------------------------------------------C
!                                                                           C
!  This routine writes a Local "Salinity Boundary Conditions Values"        C
!  (fort.36) file for each subdomain using the domain decomposition of      C
!  the ADCIRC grid created by the routine DECOMP.                           C
!                                                                           C
!                Added by Kendra Dresback (January 15, 2008)                C
!---------------------------------------------------------------------------C

    IMPLICIT NONE
    integer :: nbytes = 0
    INTEGER :: I,J,IPROC
    INTEGER :: SDU(NPROC) ! subdomain unit numbers
    LOGICAL :: Success    ! .TRUE. if files opened without errors
    CHARACTER(40) ::  ETIMINC,GRIDINC
    INTEGER :: NODP, M
    INTEGER,ALLOCATABLE :: NOD(:)
    REAL(SZ),ALLOCATABLE :: SalBC(:,:)
    REAL(SZ),ALLOCATABLE :: RESBCBINP(:)



!--Enter, Locate, Open, and Read the ADCIRC UNIT 36
!  Global Salinity Boundary Conditions file for baroclinic

!     Open full domain and subdomain fort.36 files
    CALL OpenPrepFiles(36, 'salinity boundary             ', &
    &      1, nproc, SDU, Success)
    IF ( .NOT. Success) THEN
        WRITE(*,*) 'WARNING: Unit 36 files not preprocessed.'
        RETURN ! note early return
    ENDIF

!--Allocate local arrays

    ALLOCATE ( NOD(MNETA) )
    ALLOCATE ( RESBCBINP(NFEN) )
    ALLOCATE ( SalBC(MNETA,NFEN) )
    nbytes = 8*mneta
    call memory_alloc(nbytes)

!--While ( NOT EOF ) Read NETA BCs from Global File

    DO  ! loop until end of file
        READ(36,40,END=9999) ETIMINC
        DO IPROC = 1,NPROC
            WRITE(SDU(IPROC),40)  ETIMINC
        ENDDO

        DO I=1, NETA
            READ(36,*,END=9999)  NOD(I), (SalBC(I,M),M=1,NFEN)
        ENDDO
    
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


!--Close Global file and all the Local Files

    9999 CLOSE (36)
    DO IPROC=1, NPROC
        CLOSE (SDU(IPROC))
    ENDDO

    deallocate(salbc)
    nbytes = 8*mneta
    call memory_dealloc(nbytes)
    call memory_status()
    40 FORMAT(A40)
    80 FORMAT(1X,I6,1X,32000(F11.7,2X))

    RETURN
    END SUBROUTINE PREP36

    SUBROUTINE PREP37()
    USE PRE_GLOBAL
    use memory_usage

!---------------------------------------------------------------------------C
!                                                                           C
!  This routine writes a Local "Temperature Boundary Conditions Values"     C
!  (fort.37) file for each subdomain using the domain decomposition of      C
!  the ADCIRC grid created by the routine DECOMP.                           C
!                                                                           C
!                Added by Kendra Dresback (January 15, 2008)                C
!---------------------------------------------------------------------------C

    IMPLICIT NONE
    integer :: nbytes = 0
    INTEGER :: I,J,IPROC
    INTEGER :: SDU(NPROC) ! subdomain unit numbers
    LOGICAL :: Success    ! .TRUE. if files opened without errors
    CHARACTER(40) ::  ETIMINC
    INTEGER :: NODP, M
    INTEGER,ALLOCATABLE :: NOD(:)
    REAL(SZ),ALLOCATABLE :: TempBC(:,:)
    REAL(SZ),ALLOCATABLE :: RESBCBINP(:)


!--Enter, Locate, Open, and Read the ADCIRC UNIT 37
!  Global Temperature Boundary Conditions file for baroclinic

!     Open full domain and subdomain fort.37 files
    CALL OpenPrepFiles(37, 'temperature boundary          ', &
    &      1, nproc, SDU, Success)
    IF ( .NOT. Success) THEN
        WRITE(*,*) 'WARNING: Unit 37 files not preprocessed.'
        RETURN ! note early return
    ENDIF

!--Allocate local arrays

    ALLOCATE ( NOD(MNETA) )
    ALLOCATE ( RESBCBINP(NFEN) )
    ALLOCATE ( TempBC(MNETA,NFEN) )
    nbytes = 8*mneta
    call memory_alloc(nbytes)


!--While ( NOT EOF ) Read NETA BCs from Global File

    DO  ! loop around until the end of the file
        READ(37,40,END=9999) ETIMINC
        DO IPROC = 1,NPROC
            WRITE(SDU(IPROC),40)  ETIMINC
        ENDDO

        DO I=1, NETA
            READ(37,*,END=9999)  NOD(I), (TempBC(I,M),M=1,NFEN)
        ENDDO
    
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


!--Close Global file and all the Local Files

    9999 CLOSE (37)
    DO IPROC=1, NPROC
        CLOSE (SDU(IPROC))
    ENDDO

    deallocate(TempBC)
    nbytes = 8*mneta
    call memory_dealloc(nbytes)
    call memory_status()
    40 FORMAT(A40)
    80 FORMAT(1X,I6,1X,32000(F11.7,2X))

    RETURN
    END SUBROUTINE PREP37

    SUBROUTINE PREP38()
    USE PRE_GLOBAL
    use memory_usage

!---------------------------------------------------------------------------C
!                                                                           C
!  This routine writes a Local "Temperature Boundary Conditions Values      C
!  for the surface" (fort.38) file for each subdomain using the domain     C
!  decomposition of the ADCIRC grid created by the routine DECOMP.          C
!                                                                           C
!                Added by Kendra Dresback (October 15, 2008)                C
!---------------------------------------------------------------------------C

    IMPLICIT NONE
    integer :: nbytes = 0
    INTEGER :: I,J,IPROC
    INTEGER :: SDU(NPROC) ! subdomain unit numbers
    LOGICAL :: Success    ! .TRUE. if files opened without errors
    CHARACTER(40) ::  ETIMINC,GRIDINC
    INTEGER :: NODP, M, NFLUX
    INTEGER,ALLOCATABLE :: NOD(:)
    REAL(SZ),ALLOCATABLE :: TopTempBC(:,:)
    REAL(SZ),ALLOCATABLE :: RESBCBINP(:,:)



!--Enter, Locate, Open, and Read the ADCIRC UNIT 38
!  Global Salinity Boundary Conditions file for baroclinic

!     Open full domain and subdomain fort.38 files
    CALL OpenPrepFiles(38, 'top temperature boundary      ', &
    &      1, nproc, SDU, Success)
    IF ( .NOT. Success) THEN
        WRITE(*,*) 'WARNING: Unit 38 files not preprocessed.'
        RETURN ! note early return
    ENDIF

!  Determine how many values are in the top temperature boundary
!  condition

    IF (BCFLAG_TEMP == 1) THEN
        NFLUX = 1
    ELSE IF (BCFLAG_TEMP == 2) THEN
        NFLUX = 6
    ELSE IF (BCFLAG_TEMP == 3) THEN
        NFLUX = 4
    END IF

    MNP=nnodg

!--Allocate local arrays

    ALLOCATE ( NOD(MNP) )
    ALLOCATE ( RESBCBINP(MNP,NFLUX) )
    ALLOCATE ( TopTempBC(MNP,NFLUX) )
    nbytes = 24*mnp
    call memory_alloc(nbytes)

!--While ( NOT EOF ) Read NETA BCs from Global File

    DO  ! loop until end of file

        READ(38,*,END=9999)  (NOD(I),(TopTempBC(I,M),M=1,NFLUX),I=1,NNODG)
    
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


!--Close Global file and all the Local Files

    9999 CLOSE (38)
    DO IPROC=1, NPROC
        CLOSE (SDU(IPROC))
    ENDDO

    deallocate(toptempbc)
    nbytes = 24*mnp
    call memory_dealloc(nbytes)
    call memory_status()
    40 FORMAT(A40)
    80 FORMAT(1X,I8,1X,32(F12.6,2X))

    RETURN
    END SUBROUTINE PREP38

    SUBROUTINE PREP39()
!---------------------------------------------------------------------------
!                                                                           C
!  This routine writes a Local river boundary file for the baroclnic        C
!  simulation (fort.39) for each subdomain using the domain                 C
!  decomposition of the ADCIRC grid created by the routine DECOMP.          C
!                                                                           C
!                Added by Kendra Dresback (January 14, 2010)                C
!---------------------------------------------------------------------------C

    USE PRE_GLOBAL
    use memory_usage
    IMPLICIT NONE
    integer :: nbytes = 0
    INTEGER :: IPROC
    INTEGER :: INDEX14, I
    REAL(SZ) :: FLUX_INC
    REAL(SZ),ALLOCATABLE ::  FLUX_VAL(:,:)
    INTEGER :: SDU(NPROC)  ! subdomain unit numbers
    LOGICAL :: Success     ! .TRUE. if all files open without error
    INTEGER :: INDX ! full domain node number for a flow boundary node
    INTEGER :: J,M     ! counter for subdomains that corrsp. to a single f.d. node
    INTEGER :: IPROC2! PE of a subdomain that matches a single full domain node
    INTEGER, ALLOCATABLE :: NOD(:)
    REAL(SZ),ALLOCATABLE ::  RESBCBINP(:)

!     Open full domain and subdomain fort.20 files
    CALL OpenPrepFiles(39, 'aperiodic river temp and salinity   ', &
    &      1, nproc, SDU, Success)
    IF ( .NOT. Success) THEN
        WRITE(*,*) 'WARNING: Unit 39 files not preprocessed.'
        RETURN ! note early return
    ENDIF

    ALLOCATE ( NOD(MNVEL) )
    ALLOCATE ( RESBCBINP(NFEN) )
    ALLOCATE ( FLUX_VAL(MNVEL,NFEN) )

!     Write Increment into all flux files

    READ(39,*) FLUX_INC
    DO IPROC=1,NPROC
        WRITE(SDU(IPROC),*) FLUX_INC
    ENDDO

!     jgf45.12 Write each full domain nodal flux value into each of the
!     subdomains that that full domain node maps to. The full domain
!     node may map to more than one subdomain node if it falls on a
!     boundary between subdomains (ghost nodes).


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
                    IF (IPROC == IPROC2) THEN ! full domain node maps to this s.d.
                        DO M=1,NFEN
                            RESBCBINP(M) = FLUX_VAL(INDX,M)
                        END DO
                        WRITE(SDU(IPROC),80) (RESBCBINP(M),M=1,NFEN)
                    ENDIF
                END DO
            END DO
        END DO
    END DO

    40 CLOSE (39)
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
    80 FORMAT(1X,32000(F11.7,2X))
!----------------------------------------------------------------------------
    END SUBROUTINE PREP39
!----------------------------------------------------------------------------


!   kmd48.33bc add information for initial condition file
    SUBROUTINE HOTINITCOND()
    USE PRE_GLOBAL
    use presizes; use memory_usage

!---------------------------------------------------------------------------C
!                     written 10/11/01 by RL                                C
!             started mods for harmonic analysis and 3D RL 5/22/03          C
!         jgf Updated for v45.06 09/07/2005 not incl. harmonic or 3D        C
!         kmd Updated for v48.33 07/07/2008 to bring in initial conditions  C
!                                                                           C
!  This routine reads the global initial condition file (fort.17)           C
!  and writes local hot start files of the same format.                     C
!                                                                           C
!---------------------------------------------------------------------------C

    IMPLICIT NONE
    integer :: nbytes = 0
    LOGICAL :: FOUND
    INTEGER :: I,J,IPROC,IINDX,IHOTSTP, not_active
    INTEGER :: IMHSF,ITHSF
    CHARACTER FNAME*60,LOCFN*14
    CHARACTER(16) :: FNAME1
    CHARACTER(8) :: FNAM8(2)
    EQUIVALENCE (FNAM8(1),FNAME1)

    INTEGER,ALLOCATABLE  :: LOC2(:),NOFF(:), domA(:)
    REAL(SZ),ALLOCATABLE :: ETA1(:),ETA2(:),EtaDisc(:), &
    UU2(:),VV2(:),CH1(:)
    REAL(8) TIMEHSF
    integer :: InputFileFmtVn, NP_G_IN, NE_G_IN, NP_A_IN, NE_A_IN
    CHARACTER(60) :: FileFmtVn
    INTEGER :: SDU(NPROC)  ! subdomain unit numbers
    LOGICAL :: Success     ! .TRUE. if all files open without error

#if 0
! vjp 2006/9/30 not supporting harmonic analysis or C3D yet
    INTEGER :: INZ,INF,IMM,INP,INSTAE,INSTAV,IISTAE,IISTAV,IIGLOE,IIGLOV, &
    IICALL,INFREQ,ITUD,NTSTEPS
    INTEGER :: ITHAS,ITHAF,ITMV,IHABEG,ICHA
    CHARACTER*10,ALLOCATABLE     ::  INAMEFR(:)
    REAL(8)  TIMEUD
    REAL(SZ),ALLOCATABLE ::  HA(:,:)
    REAL(SZ),ALLOCATABLE ::  ELAV(:),ELVA(:),XVELAV(:),XVELVA(:), &
    YVELAV(:),YVELVA(:)
    REAL(SZ),ALLOCATABLE ::  IFREQ(:),IFF(:),IFACE(:)
    REAL(SZ),ALLOCATABLE ::  GLOELV(:,:)
    REAL(SZ),ALLOCATABLE ::  GLOULV(:,:),GLOVLV(:,:)
    REAL(SZ),ALLOCATABLE ::  STAELV(:,:)
    REAL(SZ),ALLOCATABLE ::  STAULV(:,:),STAVLV(:,:)
#endif

!--   Open the Initial Condition Start File based on the value of IHOT from
!--   the fort.15 file

!     Open full domain and subdomain fort.17 files
    Print *, "Made it to prepping the files"
    CALL OpenPrepFiles(17, 'initial condition file  ', &
    &      1, nproc, SDU, Success)
    Print *, "Made it out of prepping the files"
    IF ( .NOT. Success) THEN
        WRITE(*,*) 'WARNING: Unit 17 files not preprocessed.'
        RETURN ! note early return
    ENDIF
    IHOT=17

!--   Read in info from global initial condition file

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

! Allocate local work arrays

    MNP = nnodg
    MNE = nelg
    nbytes = 4*nproc
    call memory_alloc(nbytes)
    ALLOCATE ( ETA1(MNP),ETA2(MNP),EtaDisc(MNP),UU2(MNP), &
    VV2(MNP),NODECODE(MNP),CH1(MNP) )
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
    IF(IM == 10) THEN
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
!     jgf46.02 Read in 3D hotstart data if appropriate
    IF (IMHSF > 10) THEN
        PRINT *, "set to go into 3D read"
        PRINT *, "NFEN = ", NFEN
        CALL ReadInitCond3D(IHOT)
    ENDIF
#if 0

!.....DETERMINE HARMONIC ANALYSIS PARAMETERS

    IHARIND=NHARFR*(NHASE+NHASV+NHAGE+NHAGV)
    IF(IHARIND > 0) IHARIND=1

!.....IF HARMONIC ANALYSIS IS INCLUDED IN THE RUN, PROCESS HOT START INFORMATION FOR
!.....IN PROGRESS HARMONIC ANALYSIS

    IF(IHARIND == 1) THEN
        ITHAS=INT((THAS-STATIM)*(86400.D0/DT) + 0.5d0)
        ITHAF=INT((THAF-STATIM)*(86400.D0/DT) + 0.5d0)
        ITMV = ITHAF - (ITHAF-ITHAS)*FMV
        IHABEG=ITHAS+NHAINC

    !.......IF HARMONIC ANALYSIS HAS ALREADY BEGUN, READ IN HOT START
    !........HARMONIC ANALYSIS, MEAN AND SQUARE INFO

        IF(ITHSF > ITHAS) THEN
            IHOTSTP=IHOTSTP+1
            READ(IHOT,REC=IHOTSTP) ICHA
        ENDIF

        IF(ITHSF >= IHABEG) THEN
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

            IF(NHASE == 1) THEN
                DO J=1,INSTAE
                    DO I=1,IMM
                        IHOTSTP=IHOTSTP+1
                        READ(IHOT,REC=IHOTSTP) STAELV(I,J)
                    ENDDO
                ENDDO
            ENDIF

            IF(NHASV == 1) THEN
                DO J=1,INSTAV
                    DO I=1,IMM
                        READ(IHOT,REC=IHOTSTP+1) STAULV(I,J)
                        READ(IHOT,REC=IHOTSTP+2) STAVLV(I,J)
                        IHOTSTP = IHOTSTP + 2
                    ENDDO
                ENDDO
            ENDIF

            IF(NHAGE == 1) THEN
                DO J=1,INP
                    DO I=1,IMM
                        IHOTSTP=IHOTSTP+1
                        READ(IHOT,REC=IHOTSTP) GLOELV(I,J)
                    ENDDO
                ENDDO
            ENDIF

            IF(NHAGV == 1) THEN
                DO J=1,INP
                    DO I=1,IMM
                        READ(IHOT,REC=IHOTSTP+1) GLOULV(I,J)
                        READ(IHOT,REC=IHOTSTP+2) GLOVLV(I,J)
                        IHOTSTP = IHOTSTP + 2
                    ENDDO
                ENDDO
            ENDIF

        ENDIF

        IF((FMV > 0.) .AND. (INFREQ > 0) .AND. (IM == 0)) THEN !include means and variances
            IF(ITHSF > ITMV) THEN
                IHOTSTP=IHOTSTP+1
                READ(IHOT,REC=IHOTSTP) NTSTEPS
                IF(NHAGE == 1) THEN
                    DO I=1,INP
                        READ(IHOT,REC=IHOTSTP+1) ELAV(I)
                        READ(IHOT,REC=IHOTSTP+2) ELVA(I)
                        IHOTSTP=IHOTSTP+2
                    ENDDO
                ENDIF
                IF(NHAGV == 1) THEN
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


!--Open All Local Hot Start files

    ALLOCATE ( LOC2(NPROC) )
    DO IPROC = 1,NPROC
        LOC2(IPROC) = 105 + (IPROC-1)
        LOCFN(1:14) = 'PE0000/'//FNAME(1:7)
        CALL IWRITE(LOCFN,3,6,IPROC-1)
        OPEN (LOC2(IPROC),FILE=LOCFN)
    ENDDO

!--Write out info to local hot start files

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

        IF(IM == 10) THEN
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

    
    !     jgf46.02 Write out 3D hotstart data if appropriate
        IF (IMHSF > 10) THEN
            CALL WriteInitCond3D(LOC2(IPROC),IPROC)
        ENDIF
#if 0
    
    !....IF APPROPRIATE, WRITE OUT HOT START INFORMATION FOR IN PROGRESS HARMONIC ANALYSIS

    !       IF((IHARIND.EQ.1).AND.(ITHSF.GT.ITHAS)) THEN
    !         WRITE(LOC2(IPROC),REC=IHOTSTP+1) ICHA
    !         IHOTSTP = IHOTSTP + 1
    !         CALL HAHOUT(NP,NSTAE,NSTAV,NHASE,NHASV,NHAGE,NHAGV,
    !    &                LOC2(IPROC),IHOTSTP)
    
    !         IF(NHASE.EQ.1) CALL HAHOUTES(NSTAE,LOC2(IPROC),IHOTSTP)
    !         IF(NHASV.EQ.1) CALL HAHOUTVS(NSTAV,LOC2(IPROC),IHOTSTP)
    !         IF(NHAGE.EQ.1) CALL HAHOUTEG(MNP,LOC2(IPROC),IHOTSTP)
    !         IF(NHAGV.EQ.1) CALL HAHOUTVG(MNP,LOC2(IPROC),IHOTSTP)
    !         ENDIF
    
    !       if(CHARMV) then
    !         IF((IHARIND.EQ.1).AND.(ITHSF.GT.ITMV)) THEN
    !           IHOTSTP=IHOTSTP+1
    !           WRITE(LOC2(IPROC),REC=IHOTSTP) NTSTEPS
    !           IF(NHAGE.EQ.1) THEN
    !             DO I=1, NNODP(IPROC)
    !               IINDX = IMAP_NOD_LG(I,IPROC)
    !               DO I=1,MNP
    !                 WRITE(LOC2(IPROC),REC=IHOTSTP+1) ELAV(IINDX)
    !                 WRITE(LOC2(IPROC),REC=IHOTSTP+2) ELVA(IINDX)
    !                 IHOTSTP=IHOTSTP+2
    !                 END DO
    !             ENDIF
    !           IF(NHAGV.EQ.1) THEN
    !             DO I=1,NNODP(IPROC)
    !               WRITE(LOC2(IPROC),REC=IHOTSTP+1) XVELAV(IINDX)
    !               WRITE(LOC2(IPROC),REC=IHOTSTP+2) YVELAV(IINDX)
    !               WRITE(LOC2(IPROC),REC=IHOTSTP+3) XVELVA(IINDX)
    !               WRITE(LOC2(IPROC),REC=IHOTSTP+4) YVELVA(IINDX)
    !               IHOTSTP=IHOTSTP+4
    !               END DO
    !             ENDIF
    !           ENDIF
    !         ENDIF
#endif

    ENDDO

!--Close Global file and all the Local Files

    CLOSE (IHOT)
    DO IPROC=1, NPROC
        CLOSE (LOC2(IPROC))
    ENDDO

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

!---------------------------------------------------------------------------C
!                     written 10/11/01 by RL                                C
!             started mods for harmonic analysis and 3D RL 5/22/03          C
!         jgf Updated for v45.06 09/07/2005 not incl. harmonic or 3D        C
!         kmd48.33bc updated with 3D information                            C
!                                                                           C
!  This routine reads the global hot start file (either fort.67 or fort.68) C
!  and writes local hot start files of the same format.                     C
!                                                                           C
!---------------------------------------------------------------------------C

    IMPLICIT NONE
    integer :: nbytes = 0
    LOGICAL :: FOUND
    INTEGER :: I,J,IPROC,INDX,IHOTSTP, not_active
    INTEGER :: IMHSF,ITHSF, NH, N
    CHARACTER FNAME*60,LOCFN*14
    CHARACTER(16) :: FNAME1
    CHARACTER(8) :: FNAM8(2)
    EQUIVALENCE (FNAM8(1),FNAME1)

    INTEGER,ALLOCATABLE  :: LOC2(:),NOFF(:), domA(:)
    REAL(SZ),ALLOCATABLE :: ETA1(:),ETA2(:),EtaDisc(:), &
    UU2(:),VV2(:),CH1(:)
    REAL(8) TIMEHSF
    integer :: InputFileFmtVn, NP_G_IN, NE_G_IN, NP_A_IN, NE_A_IN

    INTEGER :: INZ,INF,IMM,INP,INSTAE,INSTAV,IISTAE,IISTAV,IIGLOE,IIGLOV, &
    IICALL,INFREQ,ITUD,NTSTEPS
    INTEGER :: IHARIND,ITHAS,ITHAF,ITMV,IHABEG,ICHA
    CHARACTER*10,ALLOCATABLE     ::  INAMEFR(:)
    REAL(8)  TIMEUD
    REAL(SZ),ALLOCATABLE ::  HA(:,:)
    REAL(SZ),ALLOCATABLE ::  ELAV(:),ELVA(:),XVELAV(:),XVELVA(:), &
    YVELAV(:),YVELVA(:)
    REAL(SZ),ALLOCATABLE ::  IFREQ(:),IFF(:),IFACE(:)
    REAL(SZ),ALLOCATABLE ::  GLOELV(:,:)
    REAL(SZ),ALLOCATABLE ::  GLOULV(:,:),GLOVLV(:,:)
    REAL(SZ),ALLOCATABLE ::  STAELV(:,:)
    REAL(SZ),ALLOCATABLE ::  STAULV(:,:),STAVLV(:,:)
    REAL(SZ) TIME

    REAL(SZ) DUMMY
    INTEGER :: IDUMMY
    INTEGER :: LUN
    INTEGER :: NHS

!--   Open Appropriate Hot Start File based on the value of IHOT from
!--   the fort.15 file

    write(*,*) "enter IHOT: "
    read(*,*) IHOT
    SELECT CASE (IHOT)
    CASE(67)
    FNAME='fort.67'
    CASE(68)
    FNAME='fort.68'
    CASE(367,368)
    write(*,*) "INFO: IHOT=",IHOT, &
    " means parallel ADCIRC should read a NetCDF hotstart file."
    write(*,*) &
    "INFO: NetCDF hotstart files do not require decomposition."
    RETURN
    CASE DEFAULT
    write(*,*) "ERROR: The IHOT value ",IHOT, &
    " is not a valid option."
    write(*,*) "INFO: 67 and 68 are the only valid options."
    RETURN
    END SELECT

    INQUIRE(FILE=FNAME,EXIST=FOUND)
    IF (FOUND) THEN
        WRITE(*,1011) FNAME
        IF(IHOT == 67 .OR. IHOT == 68) &
        OPEN(IHOT,FILE=FNAME,ACCESS='DIRECT',RECL=8)
    ELSE
        WRITE(*,1010) FNAME
        STOP
    ENDIF

!--   Read in info from global hot start files

    IHOTSTP=1
    READ(IHOT,REC=IHOTSTP) InputFileFmtVn ; IHOTSTP = IHOTSTP + 1

    if ( .NOT. CMP_VERSION_NUMBERS(InputFileFmtVn, FileFmtVersion)) then
        write(*, 1001)
        write(*, 1012)
        write(*, 9973)
    ! top
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

! Allocate local work arrays

    nbytes = 4*nproc
    call memory_alloc(nbytes)
    ALLOCATE ( ETA1(MNP),ETA2(MNP),EtaDisc(MNP),UU2(MNP), &
    VV2(MNP),NODECODE(MNP),CH1(MNP) )
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
    IF(IMHSF == 10) THEN
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

! kmd48.33 moved 3D hot start information to subroutine
!          and took out other lines
! jgf49.17 refined check of IMHSF so that it picks up only
!     IM values that indicate 3D (and so we can use six integer IM values).
    IF ((IMHSF == 1) .OR. (IMHSF == 11) .OR. &
    (IMHSF == 21) .OR. (IMHSF == 31)) THEN
        CALL ReadHotStart3D(IHOT,IHOTSTP)
    ENDIF

!     jgf48.03 harmonic analysis not supported yet
#if 0

!....DETERMINE HARMONIC ANALYSIS PARAMETERS

    IHARIND=NHARFR*(NHASE+NHASV+NHAGE+NHAGV)
    IF(IHARIND > 0) IHARIND=1

!.....IF HARMONIC ANALYSIS IS INCLUDED IN THE RUN, PROCESS HOT START
!     INFORMATION FOR IN PROGRESS HARMONIC ANALYSIS

    IF(IHARIND == 1) THEN
        ITHAS=INT((THAS-STATIM)*(86400.D0/DT) + 0.5d0)
        ITHAF=INT((THAF-STATIM)*(86400.D0/DT) + 0.5d0)
        ITMV = ITHAF - (ITHAF-ITHAS)*FMV
        IHABEG=ITHAS+NHAINC

    !.......IF HARMONIC ANALYSIS HAS ALREADY BEGUN, READ IN HOT START
    !........HARMONIC ANALYSIS, MEAN AND SQUARE INFO

        IF(ITHSF > ITHAS) THEN
            READ(IHOT,REC=IHOTSTP) ICHA
            IHOTSTP=IHOTSTP+1
        ENDIF

        IF(ITHSF >= IHABEG) THEN
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

            IF(NHASE == 1) THEN
                DO J=1,INSTAE
                    DO I=1,IMM
                        READ(IHOT,REC=IHOTSTP) STAELV(I,J)
                        IHOTSTP=IHOTSTP+1
                    ENDDO
                ENDDO
            ENDIF

            IF(NHASV == 1) THEN
                DO J=1,INSTAV
                    DO I=1,IMM
                        READ(IHOT,REC=IHOTSTP) STAULV(I,J)
                        IHOTSTP = IHOTSTP + 1
                        READ(IHOT,REC=IHOTSTP) STAVLV(I,J)
                        IHOTSTP = IHOTSTP + 1
                    ENDDO
                ENDDO
            ENDIF

            IF(NHAGE == 1) THEN
                DO J=1,INP
                    DO I=1,IMM
                        READ(IHOT,REC=IHOTSTP) GLOELV(I,J)
                        IHOTSTP=IHOTSTP+1
                    ENDDO
                ENDDO
            ENDIF

            IF(NHAGV == 1) THEN
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

        IF((FMV > 0.) .AND. (INFREQ > 0) .AND. (IM == 0)) THEN !include means and variances
            IF(ITHSF > ITMV) THEN
                READ(IHOT,REC=IHOTSTP) NTSTEPS
                IHOTSTP=IHOTSTP+1
                IF(NHAGE == 1) THEN
                    DO I=1,INP
                        READ(IHOT,REC=IHOTSTP) ELAV(I)
                        IHOTSTP=IHOTSTP+1
                        READ(IHOT,REC=IHOTSTP) ELVA(I)
                        IHOTSTP=IHOTSTP+1
                    ENDDO
                ENDIF
                IF(NHAGV == 1) THEN
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

!--Open All Local Hot Start files


    ALLOCATE ( LOC2(NPROC) )
    DO IPROC = 1,NPROC
        LOC2(IPROC) = 105 + (IPROC-1)
        LOCFN(1:14) = 'PE0000/'//FNAME(1:7)
        CALL IWRITE(LOCFN,3,6,IPROC-1)
        OPEN (LOC2(IPROC),FILE=LOCFN,ACCESS='DIRECT',RECL=8)
    ENDDO

!--Write out info to local hot start files

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

        IF(IM == 10) THEN
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

    !  kmd48.33bc moved 3D hot start information to subroutine
    ! jgf49.43 refined check of IMHSF so that it picks up only
    !     IM values that indicate 3D (and so we can use six integer IM values).
        IF ((IMHSF == 1) .OR. (IMHSF == 11) .OR. &
        (IMHSF == 21) .OR. (IMHSF == 31)) THEN
            CALL WriteHotStart3D(LOC2(IPROC),IHOTSTP,IPROC)
        ENDIF

#if 0
    
    !....IF APPROPRIATE, WRITE OUT HOT START INFORMATION FOR IN PROGRESS HARMONIC ANALYSIS

    !       IF((IHARIND.EQ.1).AND.(ITHSF.GT.ITHAS)) THEN
    !         WRITE(LOC2(IPROC),REC=IHOTSTP+1) ICHA
    !         IHOTSTP = IHOTSTP + 1
    !         CALL HAHOUT(NP,NSTAE,NSTAV,NHASE,NHASV,NHAGE,NHAGV,
    !    &                LOC2(IPROC),IHOTSTP)
    
    !         IF(NHASE.EQ.1) CALL HAHOUTES(NSTAE,LOC2(IPROC),IHOTSTP)
    !         IF(NHASV.EQ.1) CALL HAHOUTVS(NSTAV,LOC2(IPROC),IHOTSTP)
    !         IF(NHAGE.EQ.1) CALL HAHOUTEG(MNP,LOC2(IPROC),IHOTSTP)
    !         IF(NHAGV.EQ.1) CALL HAHOUTVG(MNP,LOC2(IPROC),IHOTSTP)
    !         ENDIF
    
    !       if(CHARMV) then
    !         IF((IHARIND.EQ.1).AND.(ITHSF.GT.ITMV)) THEN
    !           IHOTSTP=IHOTSTP+1
    !           WRITE(LOC2(IPROC),REC=IHOTSTP) NTSTEPS
    !           IF(NHAGE.EQ.1) THEN
    !             DO I=1, NNODP(IPROC)
    !               INDX = IMAP_NOD_LG(I,IPROC)
    !               DO I=1,MNP
    !                 WRITE(LOC2(IPROC),REC=IHOTSTP+1) ELAV(INDX)
    !                 WRITE(LOC2(IPROC),REC=IHOTSTP+2) ELVA(INDX)
    !                 IHOTSTP=IHOTSTP+2
    !                 END DO
    !             ENDIF
    !           IF(NHAGV.EQ.1) THEN
    !             DO I=1,NNODP(IPROC)
    !               WRITE(LOC2(IPROC),REC=IHOTSTP+1) XVELAV(INDX)
    !               WRITE(LOC2(IPROC),REC=IHOTSTP+2) YVELAV(INDX)
    !               WRITE(LOC2(IPROC),REC=IHOTSTP+3) XVELVA(INDX)
    !               WRITE(LOC2(IPROC),REC=IHOTSTP+4) YVELVA(INDX)
    !               IHOTSTP=IHOTSTP+4
    !               END DO
    !             ENDIF
    !           ENDIF
    !         ENDIF
#endif

    ENDDO

!--Close Global file and all the Local Files

    CLOSE (IHOT)
    DO IPROC=1, NPROC
        CLOSE (LOC2(IPROC))
    ENDDO

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

!---------------------------------------------------------------------------C
!                     written 10/11/01 by RL                                C
!             started mods for harmonic analysis and 3D RL 5/22/03          C
!         jgf Updated for v45.06 09/07/2005 not incl. harmonic or 3D        C
!         kmd48.33bc updated with 3D hot start                              C
!                                                                           C
!  This routine reads the global hot start file (either fort.67 or fort.68) C
!  and writes local hot start files of the same format.                     C
!                                                                           C
!---------------------------------------------------------------------------C

    IMPLICIT NONE
    integer :: nbytes = 0
    LOGICAL :: FOUND
    INTEGER :: I,J,IPROC,INDX,IHOTSTP, not_active
    INTEGER :: IMHSF,ITHSF,IVALUE,IDUMY, NH, N
    CHARACTER FNAME*60,LOCFN*14
    CHARACTER(16) :: FNAME1
    CHARACTER(8) :: FNAM8(2)
    EQUIVALENCE (FNAM8(1),FNAME1)

    INTEGER,ALLOCATABLE  :: LOC2(:),NOFF(:), domA(:)
    REAL(SZ),ALLOCATABLE :: ETA1(:),ETA2(:),EtaDisc(:), &
    UU2(:),VV2(:),CH1(:)
    REAL(8) TIMEHSF, RVALUE
    integer :: InputFileFmtVn, NP_G_IN, NE_G_IN, NP_A_IN, NE_A_IN

#if 0
! vjp 2006/9/30 not supporting harmonic analysis or C3D yet
    INTEGER :: INZ,INF,IMM,INP,INSTAE,INSTAV,IISTAE,IISTAV,IIGLOE,IIGLOV, &
    IICALL,INFREQ,ITUD,NTSTEPS
    INTEGER :: IHARIND,ITHAS,ITHAF,ITMV,IHABEG,ICHA
    CHARACTER*10,ALLOCATABLE     ::  INAMEFR(:)
    REAL(8)  TIMEUD
    REAL(SZ),ALLOCATABLE ::  HA(:,:)
    REAL(SZ),ALLOCATABLE ::  ELAV(:),ELVA(:),XVELAV(:),XVELVA(:), &
    YVELAV(:),YVELVA(:)
    REAL(SZ),ALLOCATABLE ::  IFREQ(:),IFF(:),IFACE(:)
    REAL(SZ),ALLOCATABLE ::  GLOELV(:,:)
    REAL(SZ),ALLOCATABLE ::  GLOULV(:,:),GLOVLV(:,:)
    REAL(SZ),ALLOCATABLE ::  STAELV(:,:)
    REAL(SZ),ALLOCATABLE ::  STAULV(:,:),STAVLV(:,:)
#endif

!--   Open Appropriate Hot Start File based on the value of IHOT from
!--   the fort.15 file

    write(*,*) "enter IHOT: "
    read(*,*) IHOT
    IF(IHOT == 67) FNAME='fort.67'
    IF(IHOT == 68) FNAME='fort.68'


!--Open All Local Hot Start files

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

! Allocate local work arrays

    MNP  =  nnodg    !  global number of nodes    ( read from fort.18 )
    print *, "MNP =", MNP

    nbytes = 4*nproc
    call memory_alloc(nbytes)
    ALLOCATE ( ETA1(MNP),ETA2(MNP),EtaDisc(MNP),UU2(MNP), &
    VV2(MNP),NODECODE(MNP),CH1(MNP) )
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


!--Read info from local hot start files

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

        IF(IM == 10) THEN
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
    
    !   kmd48.33bc add information for 3D hot start
    !     jgf46.02 read in 3D hotstart data if appropriate
        IF (C3D) THEN
            CALL ReadHotStart3DGlobal(LOC2(IPROC),IHOTSTP,IPROC)
        ENDIF
#if HA
    
    !....IF APPROPRIATE, WRITE OUT HOT START INFORMATION FOR IN PROGRESS HARMONIC ANALYSIS

    !       IF((IHARIND.EQ.1).AND.(ITHSF.GT.ITHAS)) THEN
    !         READ(LOC2(IPROC),REC=IHOTSTP+1) ICHA
    !         IHOTSTP = IHOTSTP + 1
    !         CALL HAHOUT(NP,NSTAE,NSTAV,NHASE,NHASV,NHAGE,NHAGV,
    !    &                LOC2(IPROC),IHOTSTP)
    
    !         IF(NHASE.EQ.1) CALL HAHOUTES(NSTAE,LOC2(IPROC),IHOTSTP)
    !         IF(NHASV.EQ.1) CALL HAHOUTVS(NSTAV,LOC2(IPROC),IHOTSTP)
    !         IF(NHAGE.EQ.1) CALL HAHOUTEG(MNP,LOC2(IPROC),IHOTSTP)
    !         IF(NHAGV.EQ.1) CALL HAHOUTVG(MNP,LOC2(IPROC),IHOTSTP)
    !         ENDIF
    
    !       if(CHARMV) then
    !         IF((IHARIND.EQ.1).AND.(ITHSF.GT.ITMV)) THEN
    !           IHOTSTP=IHOTSTP+1
    !           READ(LOC2(IPROC),REC=IHOTSTP) NTSTEPS
    !           IF(NHAGE.EQ.1) THEN
    !             DO I=1, NNODP(IPROC)
    !               INDX = IMAP_NOD_LG(I,IPROC)
    !               DO I=1,MNP
    !                 READ(LOC2(IPROC),REC=IHOTSTP+1) ELAV(INDX)
    !                 READ(LOC2(IPROC),REC=IHOTSTP+2) ELVA(INDX)
    !                 IHOTSTP=IHOTSTP+2
    !                 END DO
    !             ENDIF
    !           IF(NHAGV.EQ.1) THEN
    !             DO I=1,NNODP(IPROC)
    !               READ(LOC2(IPROC),REC=IHOTSTP+1) XVELAV(INDX)
    !               READ(LOC2(IPROC),REC=IHOTSTP+2) YVELAV(INDX)
    !               READ(LOC2(IPROC),REC=IHOTSTP+3) XVELVA(INDX)
    !               READ(LOC2(IPROC),REC=IHOTSTP+4) YVELVA(INDX)
    !               IHOTSTP=IHOTSTP+4
    !               END DO
    !             ENDIF
    !           ENDIF
    !         ENDIF
#endif
        CLOSE (LOC2(IPROC))

    ENDDO

!-----------------------------------------------------------------------
!--   Write info to global hot start files
!-----------------------------------------------------------------------

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
    IF(IM == 10) THEN
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


! vjp 2006/9/30 not supporting harmonic analysis or C3D yet
!   kmd48.33bc add information for 3D hot start
!     jgf46.02 Write in 3D hotstart data if appropriate
    IF (C3D) THEN
        CALL WriteHotStart3DGlobal(IHOT,IHOTSTP,IPROC)
    ENDIF

#if 0
!.....DETERMINE HARMONIC ANALYSIS PARAMETERS

    IHARIND=NHARFR*(NHASE+NHASV+NHAGE+NHAGV)
    IF(IHARIND > 0) IHARIND=1

!.....IF HARMONIC ANALYSIS IS INCLUDED IN THE RUN, PROCESS HOT START INFORMATION FOR
!.....IN PROGRESS HARMONIC ANALYSIS

    IF(IHARIND == 1) THEN
        ITHAS=INT((THAS-STATIM)*(86400.D0/DT) + 0.5d0)
        ITHAF=INT((THAF-STATIM)*(86400.D0/DT) + 0.5d0)
        ITMV = ITHAF - (ITHAF-ITHAS)*FMV
        IHABEG=ITHAS+NHAINC

    !.......IF HARMONIC ANALYSIS HAS ALREADY BEGUN, READ IN HOT START
    !........HARMONIC ANALYSIS, MEAN AND SQUARE INFO

        IF(ITHSF > ITHAS) THEN
            WRITE(IHOT,REC=IHOTSTP) ICHA
            IHOTSTP=IHOTSTP+1
        ENDIF

        IF(ITHSF >= IHABEG) THEN
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

            IF(NHASE == 1) THEN
                DO J=1,INSTAE
                    DO I=1,IMM
                        WRITE(IHOT,REC=IHOTSTP) STAELV(I,J)
                        IHOTSTP=IHOTSTP+1
                    ENDDO
                ENDDO
            ENDIF

            IF(NHASV == 1) THEN
                DO J=1,INSTAV
                    DO I=1,IMM
                        WRITE(IHOT,REC=IHOTSTP) STAULV(I,J)
                        IHOTSTP=IHOTSTP+1
                        WRITE(IHOT,REC=IHOTSTP) STAVLV(I,J)
                        IHOTSTP=IHOTSTP+1
                    ENDDO
                ENDDO
            ENDIF

            IF(NHAGE == 1) THEN
                DO J=1,INP
                    DO I=1,IMM
                        WRITE(IHOT,REC=IHOTSTP) GLOELV(I,J)
                        IHOTSTP=IHOTSTP+1
                    ENDDO
                ENDDO
            ENDIF

            IF(NHAGV == 1) THEN
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

        IF((FMV > 0.) .AND. (INFREQ > 0) .AND. (IM == 0)) THEN !include means and variances
            IF(ITHSF > ITMV) THEN
                WRITE(IHOT,REC=IHOTSTP) NTSTEPS
                IHOTSTP=IHOTSTP+1
                IF(NHAGE == 1) THEN
                    DO I=1,INP
                        WRITE(IHOT,REC=IHOTSTP) ELAV(I)
                        IHOTSTP=IHOTSTP+1
                        WRITE(IHOT,REC=IHOTSTP) ELVA(I)
                        IHOTSTP=IHOTSTP+1
                    ENDDO
                ENDIF
                IF(NHAGV == 1) THEN
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


!--Close Global file and all the Local Files

    CLOSE (IHOT)

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

    RETURN
    1001 FORMAT('ERROR: The hot start file')
    1010 FORMAT(' File ',A60,/,' WAS NOT FOUND!  ADCPrep Terminated!!!',/)
    1011 FORMAT(' File ',A60,/,' WAS FOUND!  Opening & Processing file',/)
    1012 FORMAT('was a nonmatching version')
    1005 FORMAT('exists but cannot be opened.')
    9973 FORMAT(/,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
    END SUBROUTINE HOTGLOBALIZE






!     ----------------------------------------------------------------------
!          S U B R O U T I N E     R E A D  H O T  S T A R T  3 D
!     ----------------------------------------------------------------------

!     jgf46.02 This subroutine supports PREP67_68. It reads in the 3D
!     section of the full domain hot start file.

!     ----------------------------------------------------------------------
    SUBROUTINE ReadHotStart3D(UnitNumber,FilePosition)
!     ----------------------------------------------------------------------
    USE PRE_GLOBAL
    IMPLICIT NONE
    INTEGER, intent(in) :: UnitNumber ! i/o unit of full domain file
    INTEGER, intent(inout) :: FilePosition ! position in binary file
    INTEGER :: IHOTSTP,NH,K

!     Start reading in the data

! md Added information for 3D hotstart
    ALLOCATE ( DUU(MNP),DUV(MNP),DVV(MNP))
    ALLOCATE ( UU(MNP),VV(MNP))
    ALLOCATE ( BSX(MNP),BSY(MNP))
! md end of additions

    PRINT *, "NFEN = ", NFEN
    IHOT=UnitNumber
    IHOTSTP=FilePosition

    PRINT *, "How many layers need to be evaluated:"
    READ *, NFEN

    PRINT *, "made it to 3D portion of code"
    PRINT *, "MNP = ", MNP
    PRINT *, "NFEN = ", NFEN
    PRINT *, "IHOTSTP = ", IHOTSTP

! md Added information for 3D hotstart
    ALLOCATE ( WZ(MNP,NFEN), q20(MNP,NFEN))
    ALLOCATE ( RealQ(MNP,NFEN), ImagQ(MNP,NFEN))
    ALLOCATE ( l(MNP,NFEN), SigT(MNP,NFEN))
    ALLOCATE ( Sal(MNP,NFEN), Temp(MNP,NFEN))
! md end of additions

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
    IF (ABS(IDEN) == 1) THEN
        DO K=1,NFEN
            DO NH=1,MNP
                READ(IHOT,REC=IHOTSTP) SigT(NH,K)
                IHOTSTP=IHOTSTP+1
            END DO
        END DO
    END IF
    IF(ABS(IDen) == 2) THEN
        DO K=1,NFEN
            DO NH=1,MNP
                READ(IHOT,REC=IHOTSTP) Sal(NH,K)
                IHOTSTP=IHOTSTP+1
            END DO
        END DO
    ENDIF
    IF(ABS(IDen) == 3) THEN
        DO K=1,NFEN
            DO NH=1,MNP
                READ(IHOT,REC=IHOTSTP) Temp(NH,K)
                IHOTSTP=IHOTSTP+1
            END DO
        END DO
    ENDIF
    IF(ABS(IDen) == 4) THEN
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
!     ----------------------------------------------------------------------
    END SUBROUTINE ReadHotStart3D
!     ----------------------------------------------------------------------


!     ----------------------------------------------------------------------
!          S U B R O U T I N E     W R I T E   H O T  S T A R T  3 D
!     ----------------------------------------------------------------------

!     jgf46.02 This subroutine supports PREP67_68. It writes out the 3D
!     section of the full domain hot start file.

!     ----------------------------------------------------------------------
    SUBROUTINE WriteHotStart3D(UnitNumber,FilePosition,IPROC)
!     ----------------------------------------------------------------------
    USE PRE_GLOBAL
    IMPLICIT NONE
    INTEGER, intent(in) :: UnitNumber ! i/o unit of subdomain file
    INTEGER, intent(inout) :: FilePosition ! position in binary file
    INTEGER :: IHOTSTP, LOCHSF, I, N, IINDX
    INTEGER, intent(in) :: IPROC

!     Start writing out the 3D hotstart information

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

    IF (ABS(IDEN) == 1) THEN
        DO N=1,NFEN
            DO I=1,NNODP(IPROC)
                IINDX=ABS(IMAP_NOD_LG(I,IPROC))
                WRITE(LOCHSF,REC=IHOTSTP) SigT(IINDX,N)
                IHOTSTP=IHOTSTP+1
            END DO
        END DO
    ENDIF
    IF(ABS(IDen) == 2) THEN
        DO N=1,NFEN
            DO I=1,NNODP(IPROC)
                IINDX=ABS(IMAP_NOD_LG(I,IPROC))
                WRITE(LOCHSF,REC=IHOTSTP) Sal(IINDX,N)
                IHOTSTP=IHOTSTP+1
            END DO
        END DO
    ENDIF
    IF(ABS(IDen) == 3) THEN
        DO N=1,NFEN
            DO I=1,NNODP(IPROC)
                IINDX=ABS(IMAP_NOD_LG(I,IPROC))
                WRITE(LOCHSF,REC=IHOTSTP) Temp(IINDX,N)
                IHOTSTP=IHOTSTP+1
            END DO
        END DO
    ENDIF
    IF(ABS(IDen) == 4) THEN
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
!     ----------------------------------------------------------------------
    END SUBROUTINE WriteHotStart3D
!     ----------------------------------------------------------------------

!   kmd48.33bc add in 3D global hot start files
!     ----------------------------------------------------------------------
!          S U B R O U T I N E   R E A D  H O T  S T A R T  3 D  G L O B A L
!     ----------------------------------------------------------------------

!     This subroutine supports PREP67_68. It reads in the 3D
!     section of the full domain hot start file.

!     ----------------------------------------------------------------------
    SUBROUTINE ReadHotStart3DGlobal(UnitNumber,FilePosition,ProcessNo)
!     ----------------------------------------------------------------------
    USE PRE_GLOBAL
    IMPLICIT NONE
    INTEGER, intent(in) :: UnitNumber ! i/o unit of full domain file
    INTEGER, intent(inout) :: FilePosition ! position in binary file
    INTEGER, intent(in) :: ProcessNo ! i/o unit of full domain file
    INTEGER :: IHOTSTP,NH,K,I,IPROC
    INTEGER :: IINDX
    REAL(8) RVALUE

!     Start reading in the data

! md Added information for 3D hotstart
    ALLOCATE ( DUU(MNP),DUV(MNP),DVV(MNP))
    ALLOCATE ( UU(MNP),VV(MNP))
    ALLOCATE ( BSX(MNP),BSY(MNP))
! md end of additions

    IHOT=UnitNumber
    IHOTSTP=FilePosition
    IPROC=ProcessNo

    PRINT*, "How many layers need to be evaluated:"
    READ *, NFEN

!        PRINT *, "Made it to the 3d portion of this"
    PRINT *, "NFEN = ", NFEN

! md Added information for 3D hotstart
    ALLOCATE ( WZ(MNP,NFEN), q20(MNP,NFEN))
    ALLOCATE ( RealQ(MNP,NFEN), ImagQ(MNP,NFEN))
    ALLOCATE ( l(MNP,NFEN), SigT(MNP,NFEN))
    ALLOCATE ( Sal(MNP,NFEN), Temp(MNP,NFEN))
! md end of additions

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
        IF (IINDX > 0) DUU(IINDX) = RVALUE
        IHOTSTP=IHOTSTP+1
    END DO
    DO I=1,NNODP(IPROC)
        IINDX=IMAP_NOD_LG(I,IPROC)
        READ(IHOT,REC=IHOTSTP) RVALUE
        IF (IINDX > 0) DUV(IINDX) = RVALUE
        IHOTSTP=IHOTSTP+1
    END DO
    DO I=1,NNODP(IPROC)
        IINDX=IMAP_NOD_LG(I,IPROC)
        READ(IHOT,REC=IHOTSTP) RVALUE
        IF (IINDX > 0) DVV(IINDX) = RVALUE
        IHOTSTP=IHOTSTP+1
    END DO
    DO I=1,NNODP(IPROC)
        IINDX=IMAP_NOD_LG(I,IPROC)
        READ(IHOT,REC=IHOTSTP) RVALUE
        IF (IINDX > 0) UU(IINDX) = RVALUE
        IHOTSTP=IHOTSTP+1
    END DO
    DO I=1,NNODP(IPROC)
        IINDX=IMAP_NOD_LG(I,IPROC)
        READ(IHOT,REC=IHOTSTP) RVALUE
        IF (IINDX > 0) VV(IINDX) = RVALUE
        IHOTSTP=IHOTSTP+1
    END DO
    DO I=1,NNODP(IPROC)
        IINDX=IMAP_NOD_LG(I,IPROC)
        READ(IHOT,REC=IHOTSTP) RVALUE
        IF (IINDX > 0) BSX(IINDX) = RVALUE
        IHOTSTP=IHOTSTP+1
    END DO
    DO I=1,NNODP(IPROC)
        IINDX=IMAP_NOD_LG(I,IPROC)
        READ(IHOT,REC=IHOTSTP) RVALUE
        IF (IINDX > 0) BSY(IINDX) = RVALUE
        IHOTSTP=IHOTSTP+1
    END DO

    DO K=1,NFEN
        DO I=1,NNODP(IPROC)
            IINDX=IMAP_NOD_LG(I,IPROC)
            READ(IHOT,REC=IHOTSTP) RVALUE
            IF (IINDX > 0) REALQ(IINDX,K) = RVALUE
            IHOTSTP=IHOTSTP+1
        END DO
    END DO
    DO K=1,NFEN
        DO I=1,NNODP(IPROC)
            IINDX=IMAP_NOD_LG(I,IPROC)
            READ(IHOT,REC=IHOTSTP) RVALUE
            IF (IINDX > 0) ImagQ(IINDX,K) = RVALUE
            IHOTSTP=IHOTSTP+1
        END DO
    END DO
    DO K=1,NFEN
        DO I=1,NNODP(IPROC)
            IINDX=IMAP_NOD_LG(I,IPROC)
            READ(IHOT,REC=IHOTSTP) RVALUE
            IF (IINDX > 0) WZ(IINDX,K) = RVALUE
            IHOTSTP=IHOTSTP+1
        END DO
    END DO
    DO K=1,NFEN
        DO I=1,NNODP(IPROC)
            IINDX=IMAP_NOD_LG(I,IPROC)
            READ(IHOT,REC=IHOTSTP) RVALUE
            IF (IINDX > 0) q20(IINDX,K) = RVALUE
            IHOTSTP=IHOTSTP+1
        END DO
    END DO
    DO K=1,NFEN
        DO I=1,NNODP(IPROC)
            IINDX=IMAP_NOD_LG(I,IPROC)
            READ(IHOT,REC=IHOTSTP) RVALUE
            IF (IINDX > 0) l(IINDX,K) = RVALUE
            IHOTSTP=IHOTSTP+1
        END DO
    END DO
    IF (ABS(IDEN) == 1) THEN
        DO K=1,NFEN
            DO I=1,NNODP(IPROC)
                IINDX=IMAP_NOD_LG(I,IPROC)
                READ(IHOT,REC=IHOTSTP) RVALUE
                IF (IINDX > 0) SigT(IINDX,K) = RVALUE
                IHOTSTP=IHOTSTP+1
            END DO
        END DO
    END IF
    IF(ABS(IDen) == 2) THEN
        DO K=1,NFEN
            DO I=1,NNODP(IPROC)
                IINDX=IMAP_NOD_LG(I,IPROC)
                READ(IHOT,REC=IHOTSTP) RVALUE
                IF (IINDX > 0) Sal(IINDX,K) = RVALUE
                IHOTSTP=IHOTSTP+1
            END DO
        END DO
    ENDIF
    IF(ABS(IDen) == 3) THEN
        DO K=1,NFEN
            DO I=1,NNODP(IPROC)
                IINDX=IMAP_NOD_LG(I,IPROC)
                READ(IHOT,REC=IHOTSTP) RVALUE
                IF (IINDX > 0) Temp(IINDX,K) = RVALUE
                IHOTSTP=IHOTSTP+1
            END DO
        END DO
    ENDIF
    IF(ABS(IDen) == 4) THEN
        DO K=1,NFEN
            DO I=1,NNODP(IPROC)
                IINDX=IMAP_NOD_LG(I,IPROC)
                READ(IHOT,REC=IHOTSTP) RVALUE
                IF (IINDX > 0) Sal(IINDX,K) = RVALUE
                IHOTSTP=IHOTSTP+1
            END DO
        END DO
        DO K=1,NFEN
            DO I=1,NNODP(IPROC)
                IINDX=IMAP_NOD_LG(I,IPROC)
                READ(IHOT,REC=IHOTSTP) RVALUE
                IF (IINDX > 0) TEMP(IINDX,K) = RVALUE
                IHOTSTP=IHOTSTP+1
            END DO
        ENDDO
    END IF

    RETURN
!     ----------------------------------------------------------------------
    END SUBROUTINE ReadHotStart3DGlobal
!     ----------------------------------------------------------------------

!   kmd48.33bc added for global 3D hot start
!     ----------------------------------------------------------------------
!       S U B R O U T I N E   W R I T E   H O T  S T A R T  3 D  G L O B A L
!     ----------------------------------------------------------------------

!     This subroutine supports PREP67_68. It writes out the 3D
!     section of the full domain hot start file.

!     ----------------------------------------------------------------------
    SUBROUTINE WriteHotStart3DGlobal(UnitNumber,FilePosition,IPROC)
!     ----------------------------------------------------------------------
    USE PRE_GLOBAL
    IMPLICIT NONE
    INTEGER, intent(in) :: UnitNumber ! i/o unit of subdomain file
    INTEGER, intent(inout) :: FilePosition ! position in binary file
    INTEGER, intent(in) :: IPROC
    INTEGER :: IHOTSTP, LOCHSF, I, N, IINDX

!     Start writing out the 3D hotstart information

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
    IF (ABS(IDEN) == 1) THEN
        DO N=1,NFEN
            DO I=1,NNODP(IPROC)
                IINDX=ABS(IMAP_NOD_LG(I,IPROC))
                WRITE(LOCHSF,REC=IHOTSTP) SigT(IINDX,N)
                IHOTSTP=IHOTSTP+1
            END DO
        END DO
    ENDIF
    IF(ABS(IDen) == 2) THEN
        DO N=1,NFEN
            DO I=1,NNODP(IPROC)
                IINDX=ABS(IMAP_NOD_LG(I,IPROC))
                WRITE(LOCHSF,REC=IHOTSTP) Sal(IINDX,N)
                IHOTSTP=IHOTSTP+1
            END DO
        END DO
    ENDIF
    IF(ABS(IDen) == 3) THEN
        DO N=1,NFEN
            DO I=1,NNODP(IPROC)
                IINDX=ABS(IMAP_NOD_LG(I,IPROC))
                WRITE(LOCHSF,REC=IHOTSTP) Temp(IINDX,N)
                IHOTSTP=IHOTSTP+1
            END DO
        END DO
    ENDIF
    IF(ABS(IDen) == 4) THEN
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
!     ----------------------------------------------------------------------
    END SUBROUTINE WriteHotStart3DGLOBAL
!     ----------------------------------------------------------------------

!   kmd48.33bc add read information for the initial condition file
!     ----------------------------------------------------------------------
!          S U B R O U T I N E     R E A D  I N I T C O N D  3 D
!     ----------------------------------------------------------------------

!     kmd47.22 reads in the 3D information from an initial condition
!     file

!     ----------------------------------------------------------------------
    SUBROUTINE ReadInitCond3D(UnitNumber)
!     ----------------------------------------------------------------------
    USE PRE_GLOBAL
    IMPLICIT NONE
    INTEGER, intent(in) :: UnitNumber ! i/o unit of full domain file
    INTEGER :: NH,K

!     Start reading in the data

! md Added information for 3D hotstart
    ALLOCATE ( DUU(MNP),DUV(MNP),DVV(MNP))
    ALLOCATE ( UU(MNP),VV(MNP))
    ALLOCATE ( BSX(MNP),BSY(MNP))
    ALLOCATE ( WZ(MNP,NFEN), q20(MNP,NFEN))
    ALLOCATE ( RealQ(MNP,NFEN), ImagQ(MNP,NFEN))
    ALLOCATE ( l(MNP,NFEN), SigT(MNP,NFEN))
    ALLOCATE ( Sal(MNP,NFEN), Temp(MNP,NFEN))
! md end of additions

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
    IF (ABS(IDEN) == 1) THEN
        DO K=1,NFEN
            DO NH=1,MNP
                READ(IHOT,*) SigT(NH,K)
            END DO
        END DO
    END IF
    IF(ABS(IDen) == 2) THEN
        DO K=1,NFEN
            DO NH=1,MNP
                READ(IHOT,*) Sal(NH,K)
            END DO
        END DO
    ENDIF
    IF(ABS(IDen) == 3) THEN
        DO K=1,NFEN
            DO NH=1,MNP
                READ(IHOT,*) Temp(NH,K)
            END DO
        END DO
    ENDIF
    IF(ABS(IDen) == 4) THEN
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

!     ----------------------------------------------------------------------
    END SUBROUTINE ReadInitCond3D
!     ----------------------------------------------------------------------


!     ----------------------------------------------------------------------
!          S U B R O U T I N E     W R I T E   I N I T C O N D  3 D
!     ----------------------------------------------------------------------

!     kmd47.22 This subroutine writes out the 3D section of the
!     full domain initial condtion file.

!     ----------------------------------------------------------------------
    SUBROUTINE WriteInitCond3D(UnitNumber,IPROC)
!     ----------------------------------------------------------------------
    USE PRE_GLOBAL
    IMPLICIT NONE
    INTEGER, intent(in) :: UnitNumber ! i/o unit of subdomain file
    INTEGER :: LOCHSF, I, N, IINDX
    INTEGER, intent(in) :: IPROC

!     Start writing out the 3D initial condition information

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
    IF (ABS(IDEN) == 1) THEN
        DO N=1,NFEN
            DO I=1,NNODP(IPROC)
                IINDX=ABS(IMAP_NOD_LG(I,IPROC))
                WRITE(LOCHSF,*) SigT(IINDX,N)
            END DO
        END DO
    ENDIF
    IF(ABS(IDen) == 2) THEN
        DO N=1,NFEN
            DO I=1,NNODP(IPROC)
                IINDX=ABS(IMAP_NOD_LG(I,IPROC))
                WRITE(LOCHSF,*) Sal(IINDX,N)
            END DO
        END DO
    ENDIF
    IF(ABS(IDen) == 3) THEN
        DO N=1,NFEN
            DO I=1,NNODP(IPROC)
                IINDX=ABS(IMAP_NOD_LG(I,IPROC))
                WRITE(LOCHSF,*) Temp(IINDX,N)
            END DO
        END DO
    ENDIF
    IF(ABS(IDen) == 4) THEN
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

    RETURN
!     ----------------------------------------------------------------------
    END SUBROUTINE WriteInitCond3D
!     ----------------------------------------------------------------------


!---------------------------------------------------------------------------
!     S U B R O U T I N E    O P E N  F U L L  D O M A I N  F I L E
!---------------------------------------------------------------------------

!     jgf47.02 This subroutine will open the full domain file

!---------------------------------------------------------------------------
    SUBROUTINE OpenFullDomainFile(UnitNumber, Description, Success)
!---------------------------------------------------------------------------
    USE PRE_GLOBAL
    IMPLICIT NONE
    INTEGER, intent(in) :: UnitNumber     ! i/o unit number to open
    CHARACTER(len=30), intent(in) :: Description ! description of file
    LOGICAL, intent(out):: Success     ! .TRUE. if file opened w/o errors
    LOGICAL :: Found               ! .TRUE. if the full domain file exists
    CHARACTER(len=80) FileName   ! name of full domain file
    CHARACTER(len=7) DefaultName! default name of full domain file
    INTEGER :: ErrorIO             ! zero if file opened successfully
    CHARACTER(len=4) skipstring ! indicates user wants to skip this file

    Found = .FALSE. 
    Success = .FALSE. 
    ErrorIO = 1
    skipstring = 'skip'

    DefaultName(1:5) = 'fort.'
    WRITE(DefaultName(6:7),2) UnitNumber

!     Determine the name of the file; if found, open it
    31 IF (USE_DEFAULT) THEN
        FileName = DefaultName
    ELSE
        WRITE(*,850) ! type skip to bypass
        WRITE(*,900) Description
        WRITE(*,910) UnitNumber
        READ(*,'(A)') FileName
        FileName = trim(FILENAME)
    ENDIF

!     Determine if full domain file exists
    INQUIRE(FILE=FileName,EXIST=FOUND)

!     If it does exist, open it
    IF ( FOUND ) THEN
        WRITE(*,1011) FileName !found
        OPEN(UNIT=UnitNumber, FILE=FileName, IOSTAT=ErrorIO)
        Success = .TRUE. 
        IF ( ErrorIO > 0 ) THEN
            WRITE(*,*) "ERROR: Full domain file exists but"
            WRITE(*,*) "cannot be opened."
            Success = .FALSE. 
        ELSE
            WRITE(*,*) "Successfully opened full domain file."
        ENDIF
    ELSE
    !     Give the user a chance to opt out of prepping this file.
        IF (FileName == skipstring) RETURN ! note the early RETURN
        WRITE(*,1010) FileName !not found
        GOTO 31
    ENDIF

    2 FORMAT(I2)
    30 FORMAT(A30)
    850 FORMAT(/,'Type ''skip'' to bypass preprocessing or')
    900 FORMAT('Enter the name of the ',A30)
    910 FORMAT('file (unit ',I3,'): ')
    1010 FORMAT('File ',A7,/,' WAS NOT FOUND! Try again or type "skip"',/)
    1011 FORMAT('File ',A7,/,' WAS FOUND!  Opening & Processing file.',/)
    RETURN
!---------------------------------------------------------------------------
    END SUBROUTINE OpenFullDomainFile
!---------------------------------------------------------------------------



!---------------------------------------------------------------------------
!     S U B R O U T I N E   O P E N  S U B D O M A I N  F I L E
!---------------------------------------------------------------------------

!     jgf47.02 This subroutine will open a single subdomain file

!---------------------------------------------------------------------------
    SUBROUTINE OpenSubDomainFile(UnitNumber, IProc, sdu, Success)
!---------------------------------------------------------------------------
    USE PRE_GLOBAL
    IMPLICIT NONE
    INTEGER, intent(in) :: UnitNumber     ! i/o unit number of full dom file
    INTEGER, intent(in) :: iproc          ! subdomain number
    INTEGER, intent(out) :: sdu  ! i/o unit nunber that was opened
    LOGICAL, intent(out):: Success     ! .TRUE. if files opened w/o errors
    LOGICAL :: Found               ! .TRUE. if the full domain file exists
    CHARACTER(len=80) FileName   ! name of full domain file
    CHARACTER(len=7) DefaultName! default name of full domain file
    INTEGER :: ErrorIO             ! zero if file opened successfully
    CHARACTER(14) :: sdFileName     ! subdomain file name
    CHARACTER(len=4) skipstring ! indicates user wants to skip this file

    Found = .FALSE. 
    Success = .FALSE. 
    ErrorIO = 1
    skipstring = 'skip'

    DefaultName(1:5) = 'fort.'
    WRITE(DefaultName(6:7),2) UnitNumber

!     Open subdomain file
    sdu = 105 + (iproc-1)
    sdFileName(1:7) = 'PE0000/'
    sdFileName(8:14) = DefaultName
    CALL IWRITE(sdFileName, 3, 6, iproc-1)
    OPEN (UNIT=sdu, FILE=sdFileName, IOSTAT=ErrorIO)
    Success = .TRUE. 
    IF ( ErrorIO > 0 ) THEN
        WRITE(*,*) "ERROR: Subdomain file cannot be opened."
        print *,sdu
        print *,sdFileName
        Success = .FALSE. 
    ENDIF
    2 FORMAT(I2)
    RETURN

!---------------------------------------------------------------------------
    END SUBROUTINE OpenSubDomainFile
!---------------------------------------------------------------------------



!---------------------------------------------------------------------------
!           S U B R O U T I N E   O P E N  P R E P  F I L E S
!---------------------------------------------------------------------------

!     jgf45.12 This subroutine will open the full domain file and each
!     of the subdomain files in the range of subdomains provided in the
!     arguments. It assumes that all unit numbers are between 10 and 99.

!     tcm v50.66.03 -- incrased unit numbers to include 10 - 999.
!---------------------------------------------------------------------------
    SUBROUTINE OpenPrepFiles(UnitNumber, Description, &
    startProc, endProc, SDU, Success)
!---------------------------------------------------------------------------
    USE PRE_GLOBAL
    IMPLICIT NONE
    INTEGER, intent(in) :: UnitNumber     ! i/o unit number to open
    CHARACTER(*), intent(in) :: Description ! description of file
    INTEGER, intent(in) :: startProc        ! subdomains to start with
    INTEGER, intent(in) :: endProc          ! subdomain to end on
    INTEGER, intent(out), dimension(nproc) :: SDU ! Subdomain unit numbers
    LOGICAL, intent(out):: Success     ! .TRUE. if files opened w/o errors
    LOGICAL :: Found               ! .TRUE. if the full domain file exists
    CHARACTER(len=80) FileName   ! name of full domain file
    CHARACTER(len=8) DefaultName! default name of full domain file  !increased from 7 to 8 tcm v50.66.03
    INTEGER :: dnlen  !length of defaultname (7 or 8)
    INTEGER :: ErrorIO             ! zero if file opened successfully
    INTEGER :: iproc               ! subdomain index
    CHARACTER(len=15) sdFileName     ! subdomain file name   !increased from 14 to 15 tcm v50.66.03
    CHARACTER(len=4) skipstring ! indicates user wants to skip this file

    Found = .FALSE. 
    Success = .FALSE. 
    ErrorIO = 1
    skipstring = 'skip'
    DefaultName(:) = ' '  !initialize to all blanks !tcm v50.77
      
    DefaultName(1:5) = 'fort.'

!.....tcm v50.66.03 increased unit number to 100's places
    if (UnitNumber < 100) then
        WRITE(DefaultName(6:7),2) UnitNumber
        dnlen = 7
    else
        WRITE(DefaultName(6:8),3) UnitNumber
        dnlen = 8
    endif


!     Determine the name of the file; if found, open it
    31 IF (USE_DEFAULT) THEN
        FileName = trim(DefaultName(1:dnlen))   !tcm v50.66.03 added trim !tcm v50.77 added 1:dnlen
    ! Casey 120402: Avoid an endless loop.  If the default file does not exist,
    !             then give the user a chance to specify the file name or skip.a
    !     ELSE
        GOTO 33
    ENDIF
    32 CONTINUE
    WRITE(*,850) ! type skip to bypass
    WRITE(*,900) Description
    WRITE(*,910) UnitNumber
    READ(*,'(A)') FileName
    FileName = trim(FILENAME)
!     ENDIF
    33 CONTINUE


!     Determine if full domain file exists
    INQUIRE(FILE=FileName,EXIST=FOUND)

!     If it does exist, open it
    IF ( FOUND ) THEN
        WRITE(*,1011) FileName !found
        OPEN(UNIT=UnitNumber, FILE=FileName, IOSTAT=ErrorIO)
        Success = .TRUE. 
        IF ( ErrorIO > 0 ) THEN
            WRITE(*,*) "ERROR: Full domain file exists but"
            WRITE(*,*) "cannot be opened."
            Success = .FALSE. 
        ENDIF
    ELSE
    !     Give the user a chance to opt out of prepping this file.
        IF (FileName == skipstring) RETURN ! note the early RETURN
        WRITE(*,1010) FileName !not found
    ! Casey 120402: Avoid an endless loop.  If the default file does not exist,
    !             then give the user a chance to specify the file name or skip.
    !        GOTO 31
        GOTO 32
    ENDIF

    If ( .NOT. Success) RETURN ! failed to open full domain file

!     Open each of the subdomain files
    DO iproc = startProc, endProc
        sdFileName(:) = ' '
        sdu(iproc) = 505 + (iproc-1)  !tcm v51.31 changed 105 to 505 to avoid conflicts with fort.141 file
        sdFileName(1:7) = 'PE0000/'
    !........tcm v50.66.03 increased unit number to 100's places
        if (UnitNumber < 100 ) then
            sdFileName(8:14) = DefaultName(1:7)
        else
            sdFileName(8:15) = DefaultName
        endif
#ifdef ADCSWAN
        sdFileName = 'PE0000/'//FileName
#endif
        CALL IWRITE(sdFileName, 3, 6, iproc-1)
        OPEN (UNIT=SDU(iproc), FILE=trim(sdFileName), IOSTAT=ErrorIO)
        Success = .TRUE. 
        IF ( ErrorIO > 0 ) THEN
            WRITE(*,*) "ERROR: Subdomain file cannot be opened."
            Success = .FALSE. 
            RETURN ! failed to open at least one subdomain file
        ENDIF
    ENDDO

    2 FORMAT(I2)
    3 FORMAT(I3)
    30 FORMAT(A30)
    850 FORMAT(/,'Type ''skip'' to bypass preprocessing or')
    900 FORMAT('Enter the name of the ',A30)
    910 FORMAT('file (unit ',I3,'): ')
    1010 FORMAT('File ',A8,/,' WAS NOT FOUND! Try again or type "skip"',/)  !increased A7 to A8  tcm v50.66.03
    1011 FORMAT('File ',A8,/,' WAS FOUND!  Opening & Processing file.',/)   !increased A7 to A8
    RETURN
!---------------------------------------------------------------------------
    END SUBROUTINE OpenPrepFiles
!---------------------------------------------------------------------------



!***********************************************************************
!   Subroutine to write out to the hotstart file (UNITS 67 and 68)     *
!   header information and the LHS matrix for the harmonic analysis    *
!                                                                      *
!                        R.L.  11/8/95                                 *
!***********************************************************************

!     SUBROUTINE HAHOUT(NP,NSTAE,NSTAV,ISTAE,ISTAV,IGLOE,IGLOV,
!    &  IOUNIT,IHOTSTP)
!     implicit none
!     INTEGER NP,NSTAE,NSTAV,ISTAE,AE,ISTAV
!     INTEGER IGLOE,IGLOV,IOUNIT,IHOTSTP,I,J
!     CHARACTER*16 FNAME
!     CHARACTER*8 FNAM8(2)
!     EQUIVALENCE (FNAM8(1),FNAME)c


!***** Write Out various parameter values

!     WRITE(IOUNIT,REC=IHOTSTP+1) NZ
!     WRITE(IOUNIT,REC=IHOTSTP+2) NF
!     WRITE(IOUNIT,REC=IHOTSTP+3) MM
!     WRITE(IOUNIT,REC=IHOTSTP+4) NP
!     WRITE(IOUNIT,REC=IHOTSTP+5) NSTAE
!     WRITE(IOUNIT,REC=IHOTSTP+6) NSTAV
!     WRITE(IOUNIT,REC=IHOTSTP+7) ISTAE
!     WRITE(IOUNIT,REC=IHOTSTP+8) ISTAV
!     WRITE(IOUNIT,REC=IHOTSTP+9) IGLOE
!     WRITE(IOUNIT,REC=IHOTSTP+10) IGLOV
!     WRITE(IOUNIT,REC=IHOTSTP+11) ICALL
!     WRITE(IOUNIT,REC=IHOTSTP+12) NFREQ
!     IHOTSTP = IHOTSTP+12

!     do i=1,nfreq+nf
!        FNAME=NAMEFR(I)
!        WRITE(IOUNIT,REC=IHOTSTP+1) FNAM8(1)
!        WRITE(IOUNIT,REC=IHOTSTP+2) FNAM8(2)
!        IHOTSTP=IHOTSTP+2
!        WRITE(IOUNIT,REC=IHOTSTP+1) hafreq(i)
!        WRITE(IOUNIT,REC=IHOTSTP+2) HAFF(i)
!        WRITE(IOUNIT,REC=IHOTSTP+3) HAFACE(i)
!        IHOTSTP=IHOTSTP+3
!     end do


!***** Write Out time of most recent H.A. update

!     WRITE(IOUNIT,REC=IHOTSTP+1) TIMEUD
!     WRITE(IOUNIT,REC=IHOTSTP+2) ITUD
!     IHOTSTP=IHOTSTP+2

!***** Write Out LHS Matrix

!     do i=1,mm
!        do j=1,mm
!           IHOTSTP = IHOTSTP + 1
!           WRITE(IOUNIT,REC=IHOTSTP) HA(I,J)
!        END DO
!     END DO

!     return
!     end subroutine

!***********************************************************************
!   Subroutine to write global elevation harmonic analysis RHS load    *
!   vector to a hot start file (UNITS 67 and 68)                       *
!                                                                      *
!                        R.L.  11/8/95                                 *
!***********************************************************************

!     SUBROUTINE HAHOUTEG(NP,IOUNIT,IHOTSTP)
!     implicit none
!     INTEGER IOUNIT
!     INTEGER NP,IHOTSTP,N,I

!***** Write Out Global Elevation RHS load vector

!     do n=1,np
!        do i=1,mm
!           IHOTSTP=IHOTSTP+1
!           WRITE(IOUNIT,REC=IHOTSTP) GLOELV(I,N)
!        end do
!     end do

!     return
!     end subroutine

!***********************************************************************
!   Subroutine to write elevation station harmonic analysis RHS load   *
!   vector to a hot start file (UNITS 67 and 68)                       *
!                                                                      *
!                        R.L.  11/8/95                                 *
!***********************************************************************

!     SUBROUTINE HAHOUTES(NSTAE,IOUNIT,IHOTSTP)
!     implicit none
!     INTEGER NSTAE,IOUNIT,IHOTSTP,N,I

!***** Write Out Station Elevation RHS load vector

!     do n=1,NSTAE
!        do i=1,mm
!           IHOTSTP=IHOTSTP+1
!           WRITE(IOUNIT,REC=IHOTSTP) STAELV(I,N)
!        end do
!     end do

!     return
!     end subroutine

!***********************************************************************
!   Subroutine to write global velocity harmonic analysis RHS load     *
!   vector to a hot start file (UNITS 67 and 68)                       *
!                                                                      *
!                        R.L.  11/8/95                                 *
!***********************************************************************

!     SUBROUTINE HAHOUTVG(NP,IOUNIT,IHOTSTP)
!     implicit none
!     INTEGER NP,IOUNIT,IHOTSTP,N,I

!***** Write Out Global Velocity RHS load vector

!     do n=1,np
!        do i=1,mm
!           IHOTSTP=IHOTSTP+1
!           WRITE(IOUNIT,REC=IHOTSTP) GLOULV(I,N)
!           IHOTSTP=IHOTSTP+1
!           WRITE(IOUNIT,REC=IHOTSTP) GLOVLV(I,N)
!        end do
!     end do

!     return
!     end subroutine

!***********************************************************************
!   Subroutine to write velocity station harmonic analysis RHS load    *
!   vector to a hot start file (UNITS 67 and 68)                       *
!                                                                      *
!                        R.L.  11/8/95                                 *
!***********************************************************************

!     SUBROUTINE HAHOUTVS(NSTAV,IOUNIT,IHOTSTP)
!     implicit none
!     INTEGER NSTAV,IOUNIT,IHOTSTP,N,I

!***** Write Out Station Velocity LHS load vector

!     do N=1,NSTAV
!        do i=1,mm
!           IHOTSTP=IHOTSTP+1
!           WRITE(IOUNIT,REC=IHOTSTP) STAULV(I,N)
!           IHOTSTP=IHOTSTP+1
!           WRITE(IOUNIT,REC=IHOTSTP) STAVLV(I,N)
!        end do
!     end do

!     return
!     end subroutine


!---------------------------------------------------------------------------C
!                     (  Serial Version  2/28/98  )                         C
!  This routine writes the domain decomposition information into a file,    C
!  "fort.80".  This file is used by the ADCIRC post-processor ADCPOST.      C
!  This version is compatible with ADCIRC version 34.03                     C
!                                                                           C
!     jgf45.07 Added subdomain->fulldomain element mapping to handle the    C
!     processing of the NOFF array.
!     jgf45.11 Added IDEN for processing fort.44 file, added additional
!     arrays to handle processing of 3D recording stations.
!---------------------------------------------------------------------------C
!C    Addition by CF   8/2007
    SUBROUTINE PREP80()
    USE PRE_GLOBAL
    USE SIZES, ONLY : MNHARF
    USE HARM, ONLY : NHASE, NHASV, NHAGE, NHAGV
    IMPLICIT NONE
    INTEGER :: I,K

    OPEN(UNIT=80,FILE='fort.80')              ! output for ADCPOST

!--Write out the domain decomposition information into a file
!  which will later be used in post-processing the results

    WRITE(80,80) RUNDES
    WRITE(80,80) RUNID
    WRITE(80,80) AGRID
    WRITE(80,'(2I8,16x,A)') NELG,NNODG,'! Total# elements & nodes'
    WRITE(80,'(I8,24x,A)') NPROC,' ! Number of processors'
    WRITE(80,'(I8,24x,A)') MNPP,'! Max nodes on any processor'
    WRITE(80,'(I8,24x,A)') MNEP,'! Max elements on any processor'!jgf45.07
    WRITE(80,'(I8,24x,A)') IM,'! IM, run type'         !jgf46.02
    WRITE(80,'(I8,24x,A)') NWS,'! NWS, wind data type' !jgf46.02
    WRITE(80,'(I8,24x,A)') abs(NSTAE),'! NSTAE'
    WRITE(80,'(I8,24x,A)') abs(NSTAV),'! NSTAV'
    IF (IM == 10) THEN
        WRITE(80,'(I8,24x,A)') abs(NSTAC),' ! NSTAC' !jgf46.02
    ENDIF
    IF (NWS /= 0) THEN
        WRITE(80,'(I8,24x,A)') abs(NSTAM),'! NSTAM' !jgf46.02
    ENDIF
    WRITE(80,'(I8,24x,A)') MNHARF,'! MNHARF'
    WRITE(80,'(2I8,16x,A)') MNWLAT,MNWLON,'! NWLON, NWLAT'

! Casey 100301: Changed I8 to I12.
    DO I = 1,NPROC
        WRITE(80,'(3I8,A33)') I-1, NNODP(I), NOD_RES_TOT(I), &
        '  ! PE, NNODP(PE), NOD_RES_TOT(PE)'
        WRITE(80,'(9I12)') (IMAP_NOD_LG(K,I),K=1,NNODP(I))
    ENDDO

    WRITE(80,*) "GLOBAL   PE     LOCAL   ( Global-to-Local Nodes )"
    DO I = 1,NNODG
        WRITE(80,1140) I, IMAP_NOD_GL(1,I)-1, IMAP_NOD_GL(2,I)
    ENDDO

!     jgf45.07 Add subdomain->fulldomain mapping to handle NOFF processing
!     IMAP_EL_LG(I,PE) = Global Element Number of Local Element I on PE
! Casey 100301: Changed I8 to I12.
    DO I = 1,NPROC
        WRITE(80,'(2I8,A33)') I-1, NELP(I), '  ! PE, NELP(PE)'
        WRITE(80,'(9I12)') (IMAP_EL_LG(K,I),K=1,NELP(I))
    ENDDO

    WRITE(80,'(I8,2E15.8,I8,A32)') NOUTE,TOUTSE,TOUTFE,NSPOOLE, &
    '   ! NOUTE,TOUTSE,TOUTFE,NSPOOLE'

    DO I = 1,NPROC
        WRITE(80,*) I,NSTAEP(I)
        DO K = 1,NSTAEP(I)
            WRITE(80,*) IMAP_STAE_LG(K,I)
        ENDDO
    ENDDO

    WRITE(80,'(I8,2E15.8,I8,A32)') NOUTV,TOUTSV,TOUTFV,NSPOOLV, &
    '   ! NOUTV,TOUTSV,TOUTFV,NSPOOLV'

    DO I = 1,NPROC
        WRITE(80,*) I,NSTAVP(I)
        DO K = 1,NSTAVP(I)
            WRITE(80,*) IMAP_STAV_LG(K,I)
        ENDDO
    ENDDO

    IF (IM == 10) THEN ! jgf46.02
        WRITE(80,'(I8,2E15.8,I8,A32)') NOUTC,TOUTSC,TOUTFC,NSPOOLC, &
        '   ! NOUTC,TOUTSC,TOUTFC,NSPOOLC'
        DO I = 1,NPROC
            WRITE(80,*) I,NSTACP(I)
            DO K = 1,NSTACP(I)
                WRITE(80,*) IMAP_STAC_LG(K,I)
            ENDDO
        ENDDO
    ENDIF

    IF (NWS /= 0) THEN ! jgf46.02
        WRITE(80,'(I8,2E15.8,I8,A32)') NOUTM,TOUTSM,TOUTFM,NSPOOLM, &
        '   ! NOUTM,TOUTSM,TOUTFM,NSPOOLM'
        DO I = 1,NPROC
            WRITE(80,*) I,NSTAMP(I)
            DO K = 1,NSTAMP(I)
                WRITE(80,*) IMAP_STAM_LG(K,I)
            ENDDO
        ENDDO
    ENDIF

    WRITE(80,'(I8,2E15.8,I8,A32)') NOUTGE, TOUTSGE,TOUTFGE,NSPOOLGE, &
    '   ! NOUTGE, TOUTSGE, TOUTFGE, NSPOOLGE'

    WRITE(80,'(I8,2E15.8,I8,A32)') NOUTGV, TOUTSGV,TOUTFGV,NSPOOLGV, &
    '   ! NOUTGV, TOUTSGV, TOUTFGV, NSPOOLGV'

    WRITE(80,'(I8,2E15.8,I8,A32)') NOUTGC, TOUTSGC,TOUTFGC,NSPOOLGC, &
    '   ! NOUTGC, TOUTSGC, TOUTFGC, NSPOOLGC'

    WRITE(80,'(I8,2E15.8,I8,A32)') NOUTGW, TOUTSGW,TOUTFGW,NSPOOLGW, &
    '   ! NOUTGW, TOUTSGW, TOUTFGW, NSPOOLGW'

    WRITE(80,'(4I4,A32)') NHASE,NHASV,NHAGE,NHAGV, &
    '   ! NHASE, NHASV, NHAGE, NHAGV'

!     -------------------------------------------------------------

!     S T A R T   3 D   D A T A

!     -------------------------------------------------------------

    WRITE(80,*) IDEN !jgf45.11 needed to post process the fort.44 file

!     -------------------------------------------------------------
!     jgf45.11 Write mappings for 3D density stations.
!     -------------------------------------------------------------
    WRITE(80,81) I3DSD, TO3DSDS, TO3DSDF, NSPO3DSD, NSTA3DD, &
    '   ! I3DSD, TO3DSDS, TO3DSDF, NSPO3DSD, NSTA3DD'
    IF(I3DSD /= 0) THEN
        DO I = 1, NPROC
            WRITE(80,*) I, NNSTA3DDP(I)
            DO K = 1, NNSTA3DDP(I)
                WRITE(80,*) IMAP_STA3DD_LG(K,I)
            ENDDO
        ENDDO
    ENDIF
!     -------------------------------------------------------------
!     jgf45.11 Write mappings for 3D velocity stations.
!     -------------------------------------------------------------
    WRITE(80,81) I3DSV, TO3DSVS, TO3DSVF, NSPO3DSV, NSTA3DV, &
    '   ! I3DSV, TO3DSVS, TO3DSVF, NSPO3DSV, NSTA3DV'
    IF(I3DSV /= 0) THEN
        DO I = 1, NPROC
            WRITE(80,*) I, NNSTA3DVP(I)
            DO K = 1, NNSTA3DVP(I)
                WRITE(80,*) IMAP_STA3DV_LG(K,I)
            ENDDO
        ENDDO
    ENDIF
!     -------------------------------------------------------------
!     jgf45.11 Write mappings for 3D turbulence stations.
!     -------------------------------------------------------------
    WRITE(80,81) I3DST, TO3DSTS, TO3DSTF, NSPO3DST, NSTA3DT, &
    '   ! I3DST, TO3DSST, TO3DFST, NSPO3DST, NSTA3DT'
    IF(I3DST /= 0) THEN
        DO I = 1, NPROC
            WRITE(80,*) I, NNSTA3DTP(I)
            DO K = 1, NNSTA3DTP(I)
                WRITE(80,*) IMAP_STA3DT_LG(K,I)
            ENDDO
        ENDDO
    ENDIF
    WRITE(80,82) I3DGD, TO3DGDS, TO3DGDF, NSPO3DGD, &
    '   ! I3DGD, TO3DGDS, TO3DGDF, NSPO3DGD'
    WRITE(80,82) I3DGV, TO3DGVS, TO3DGVF, NSPO3DGV, &
    '   ! I3DGV, TO3DGVS, TO3DGVF, NSPO3DGV'
    WRITE(80,82) I3DGT, TO3DGTS, TO3DGTF, NSPO3DGT, &
    '   ! I3DGT, TO3DGTS, TO3DGTF, NSPO3DGT'

!     End 3D data

    WRITE(80,*) NBYTE

    CLOSE(80)

    80 FORMAT(A80)
    81 FORMAT(I8,2E15.8,2I8,A32)
    82 FORMAT(I8,2E15.8,I8,A32)
    1130 FORMAT(8X,9I8)
    1140 FORMAT(8X,3I8)

    RETURN
    END SUBROUTINE PREP80


!---------------------------------------------------------------------------
!                S U B R O U T I N E   P R E P   U N S W A N
!---------------------------------------------------------------------------

!     jgf48.17 This subroutine just copies the UnSWAN input file to each
!     subdomain. Currently, it is only called if the compiler directive
!     ADCSWAN is defined.

!     jgf50.15: After discussion with Casey, it was decided that changes
!     should be made such that adcprep looks for the swaninit file and
!     parses it to find the name of the swan control file (conventionally
!     called fort.26 in the adcirc community). It should then attempt to
!     copy that file to the subdomains.
!---------------------------------------------------------------------------
    SUBROUTINE  prepUnSWAN
!---------------------------------------------------------------------------
    USE PRE_GLOBAL
    use memory_usage
    IMPLICIT NONE
    INTEGER :: I,J,IPROC,IPROC2,ILNODE,INDX,NHG,LINDEX
    CHARACTER(170) :: Line ! line of data from UnSWAN file
    INTEGER :: SDU(NPROC)  ! subdomain unit numbers
    LOGICAL :: Success     ! .TRUE. if all files open without error
    LOGICAL :: swanInitFound ! .TRUE. if the swaninit file was found
    LOGICAL :: origUseDefault ! there is no default UnSWAN input file name
    INTEGER :: swiLUN ! logical unit number of swaninit file
    INTEGER :: ErrorIO ! zero if the file was opened successfully
    LOGICAL :: readError ! .TRUE. if there was an error reading a file
    CHARACTER(36) swanComFile ! name of swan command file
    LOGICAL :: swanComFound ! .TRUE. if swan command file was found
    CHARACTER(len=43) sdFileName     ! subdomain file name


! Casey 110627: Set the unit number.
    swiLUN = 26

! Casey 120402: Changes through this section.  The procedure here is:
!             1. Check for the swaninit file.  If it exists, then read the name
!                of the SWAN control file from the fourth line of swaninit.
!             2. If the swaninit file does not exist, then the user may be using
!                the default control file of INPUT.  Check for the INPUT file.
!                If it exists, then use it as the SWAN control file.
!             3. Otherwise, notify the user of the error.

    swanInitFound = .FALSE. 
    swanComFound = .FALSE. 
    readError = .TRUE. 
    INQUIRE(FILE='swaninit',EXIST=swanInitFound)
    IF (swanInitFound.eqv. .TRUE. ) THEN
        WRITE(*,*) "INFO: The swaninit file was found."
        OPEN(swiLUN,FILE='swaninit',ACTION='READ', &
        ACCESS='SEQUENTIAL',IOSTAT=ErrorIO,STATUS='OLD')
        IF (ErrorIO /= 0) THEN
            WRITE(*,*) "ERROR: The swaninit file could not be opened."
        ELSE
            WRITE(*,*) &
            "INFO: Parsing swaninit file for swan command file name."
        ! skip down to the 4th line
            DO I=1,4
                READ(swiLUN,*,ERR=4321,END=4321,IOSTAT=ErrorIO) Line
            ENDDO
            READ(Line,*,ERR=4321,END=4321,IOSTAT=ErrorIO) swanComFile
        ! Casey 110627: Close the swaninit file.
            CLOSE(UNIT=swiLUN,STATUS='KEEP')
            WRITE(*,*) "INFO: The swan command file is '", &
            trim(swanComFile),"'."
            readError = .FALSE. 
        ! Casey 120402: Changes to handle the default INPUT control file.
        ! If swan was started without a swaninit file, it will create
        ! its own swaninit, and it will use the string INPUT to
        ! represent the name of the swan command file.
        !           IF (TRIM(swanComFile).eq.'INPUT') THEN
        !              WRITE(*,*) "ERROR: 'INPUT' is an invalid name."
        !           ELSE
        ! check to see if the swan command file is present
            INQUIRE(FILE=trim(swanComFile),EXIST=swanComFound)
            IF (swanComFound.eqv. .TRUE. ) THEN
                WRITE(*,*) "INFO: The swan command file '", &
                trim(swanComFile),"' was found."
                OPEN(26,FILE=trim(swanComFile),ACTION='READ', &
                ACCESS='SEQUENTIAL',IOSTAT=ErrorIO,STATUS='OLD')
                IF (ErrorIO /= 0) THEN
                    WRITE(*,*) "ERROR: The swan command file '", &
                    trim(swanComFile),"' could not be opened."
                ENDIF
            ELSE
                WRITE(*,*) "ERROR: The swan command file '", &
                trim(swanComFile),"' was not found."
            ENDIF
        !           ENDIF
        ENDIF
    ! Casey 120402: Changes to handle the default INPUT control file.
    ELSE
        INQUIRE(FILE='INPUT',EXIST=swanComFound)
        IF(swanComFound.eqv. .TRUE. )THEN
            WRITE(swanComFile,'(A)') "INPUT"
            readError = .FALSE. 
            ErrorIO = 0
            swanInitFound = .TRUE. 
            swanComFound = .TRUE. 
            WRITE(*,*) "INFO: The swan command file is '", &
            trim(swanComFile),"'."
        ! Casey 120402: Changes to handle the default INPUT control file.
            OPEN(UNIT=26,FILE=TRIM(swanComFile),ACTION='READ')
        ENDIF
    ENDIF
    4321 IF ((readError.eqv. .TRUE. ) .OR. (ErrorIO /= 0) .OR. &
    (swanInitFound.eqv. .FALSE. ) .OR. &
    (swanComFound.eqv. .FALSE. )) THEN
        WRITE(*,*) 'ERROR: There was an error reading swan files.'
        WRITE(*,*) 'WARNING: swan files not preprocessed.'
        RETURN ! note early return
    ENDIF

!     Open each of the subdomain files
    DO iproc = 1, nproc
        sdu(iproc) = 105 + (iproc-1)
        sdFileName = 'PE0000/'//swanComFile
        CALL IWRITE(sdFileName, 3, 6, iproc-1)
        OPEN (UNIT=SDU(iproc), FILE=TRIM(sdFileName), IOSTAT=ErrorIO)
        IF ( ErrorIO > 0 ) THEN
            WRITE(*,*) "ERROR: Subdomain file '",TRIM(sdFileName), &
            " cannot be opened."
            RETURN ! failed to open at least one subdomain file
        ENDIF
    ENDDO

! Casey 090304: Changed the formatting through the next section.

    DO
        READ(26,'(A)',END=9999) Line
        DO IPROC = 1,NPROC
            WRITE(SDU(IPROC),'(A)') trim(Line)
        ENDDO
    ENDDO

!--Close fulldomain file and all the subdomain files

    9999 CLOSE (26)
    DO IPROC=1, NPROC
        CLOSE (SDU(IPROC))
    ENDDO

    RETURN
    60 FORMAT(A60)
    170 FORMAT(A170)
    1100 FORMAT(I8,3E13.5)
    1101 FORMAT('#')
!----------------------------------------------------------------------------
    END SUBROUTINE prepUnSWAN
!----------------------------------------------------------------------------

!---------------------------------------------------------------------------
!                S U B R O U T I N E   P R E P   N E T C D F
!---------------------------------------------------------------------------

!     jgf49.31 This subroutine will initialize the fulldomain netcdf
!     output files for the parallel run. These files are initialized
!     in the adcprep phase because they contain fulldomain data which
!     are not available to the subdomains during a parallel run.

!---------------------------------------------------------------------------
    SUBROUTINE prepNetCDF()
#ifdef ADCNETCDF
    use sizes, only : ASCII, XDMF
    USE PRESIZES, ONLY: NM, MNTIF, MNE, MNP
    USE HARM, ONLY : IHARIND, CHARMV
    USE VERSION
    USE PRE_GLOBAL, ONLY: NSTAE, NOUTE, RUNDES, RUNID, AGRID, &
    DTDP, TAU0, STATIM, REFTIM, RNDAY, DRAMP, A00, B00, C00, &
    H0, SL0, SF0, CF, ESLM, CORI, DP, NBDV, NBVV, &
    SLEL, SFEL, ICS, NT, NTRSPE, IHOT, NOLIBF, NOLIFA, IM, &
    XEV,YEV,XEM,YEM,XEL,YEL,X,Y,NELG,DT, NTIF, NODECODE, &
    NOLICA, NOLICAT, NWP, NCOR, NTIP, NWS, NRS, NRAMP, NBFR, &
    NHY, NOPE, NETA, NBOU, NVEL, G, SL1, SF1, &
    NVDLL, NVELL, IBTYPE, IBTYPEE, &
    NSTAV, NOUTV, SLEV, SFEV, NTRSPV, &
    NSTAM, NOUTM, SLEM, SFEM, NTRSPM, &
    NOUTGE, SLAM, SFEA, NDSETSE, &
    NOUTGV, NDSETSV, &
    NOUTGW, NDSETSW, &
    IESTP,NSCOUE,IVSTP,NSCOUV,ICSTP,NSCOUC,IPSTP,IWSTP,NSCOUM, &
    IGEP,NSCOUGE,IGVP,NSCOUGV,IGCP,NSCOUGC,IGPP,IGWP,NSCOUGW, &
    nhstar,nhsinc, &
    title, institution, &
    source, history, references, comments, host, convention, &
    contact, base_date, NABOUT, NSCREEN, inundationOutput, &
    I3DSD, NFEN, NSTA3DD, I3DSV, NSTA3DV, I3DST, NSTA3DT, &
    I3DGD, I3DGV, I3DGT, IDEN, islip, kp, z0s, z0b, theta1, &
    theta2, ievc, evmin, evcon, alp1, alp2, alp3, igc, nlsd, &
    nvsd, nltd, nvtd, alp4, C3D, X3DD, Y3DD, SL3DD, SF3DD, &
    X3DV, Y3DV, SL3DV, SF3DV, X3DT, Y3DT, SL3DT, SF3DT, C3D
    USE NODALATTRIBUTES, ONLY : outputTau0
    USE NETCDFIO, ONLY : setADCIRCParameters, initNetCDFOutputFile, &
    initNetCDFHotstart3D, initNetCDFHotstart, &
    initNetCDFHotstartHarmonic, &
    initNetCDFHotstartHarmonicMeansVariances, &
    freeNetCDFCoord
#ifdef ADCSWAN
    USE GLOBAL, ONLY : OutputDataDescript_t, screenUnit, myProc, &
    SWAN_OutputHS,SWAN_OutputDIR,SWAN_OutputTM01, &
    SWAN_OutputTPS,SWAN_OutputWIND,SWAN_OutputTM02, &
    SWAN_OutputTMM10,NOUT_TVW
#else
    USE GLOBAL, ONLY : OutputDataDescript_t, screenUnit, myProc, &
    NOUT_TVW
#endif

    IMPLICIT NONE
    LOGICAL :: reterr

!     The initialization of the output data descriptors for each ADCIRC
!     output file had to be cut and pasted from write_output.F. At some
!     point in the future, adcprep will be part of adcirc, making this
!     unfortunate cut-and-paste duplication unnecessary.
    type(OutputDataDescript_t), SAVE :: ElevStaDescript    ! fort.61
    type(OutputDataDescript_t), SAVE :: VelStaDescript     ! fort.62
    type(OutputDataDescript_t), SAVE :: ElevDescript       ! fort.63
    type(OutputDataDescript_t), SAVE :: Tau0Descript       ! fort.90
    type(OutputDataDescript_t), SAVE :: VelDescript        ! fort.64
    type(OutputDataDescript_t), SAVE :: PrStaDescript      ! fort.71
    type(OutputDataDescript_t), SAVE :: WindVelStaDescript ! fort.72
    type(OutputDataDescript_t), SAVE :: PrDescript         ! fort.73
    type(OutputDataDescript_t), SAVE :: WindVelDescript    ! fort.74
    type(OutputDataDescript_t), SAVE :: weirElevDescript   ! fort.77
    type(OutputDataDescript_t), SAVE :: EtaMaxDescript     ! maxele.63
    type(OutputDataDescript_t), SAVE :: UMaxDescript       ! maxvel.63
    type(OutputDataDescript_t), SAVE :: PrMinDescript      ! minpr.63
    type(OutputDataDescript_t), SAVE :: WVMaxDescript      ! maxwvel.63
    type(OutputDataDescript_t), SAVE :: RSMaxDescript
    type(OutputDataDescript_t), SAVE :: InundationTimeDescript ! inundationtime.63
    type(OutputDataDescript_t), SAVE :: MaxInunDepthDescript   ! maxinundepth.63
    type(OutputDataDescript_t), SAVE :: InitiallyDryDescript   ! initiallydry.63
    type(OutputDataDescript_t), SAVE :: EndRisingInunDescript  ! endrisinginun.63
    type(OutputDataDescript_t), SAVE :: EverDriedDescript      ! everdried.63

!      tcm v50.75 moved RSDescript outside of the ifdef adcswan for use with
!      nrs = 3 or nrs= 4
    type(OutputDataDescript_t), SAVE :: RSDescript

#ifdef ADCSWAN
! bell 20120510: SWAN Output Data
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
!     3D output data
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
!     For hotstart files:
    type(OutputDataDescript_t), SAVE :: Elev1Descript
    type(OutputDataDescript_t), SAVE :: Elev2Descript
    type(OutputDataDescript_t), SAVE :: CH1Descript
    type(OutputDataDescript_t), SAVE :: EtaDiscDescript
    type(OutputDataDescript_t), SAVE :: NodeCodeDescript
    type(OutputDataDescript_t), SAVE :: NOFFDescript
!     for hotstart 3D data
    type(OutputDataDescript_t),SAVE :: Duudescript
    type(OutputDataDescript_t),SAVE :: Duvdescript
    type(OutputDataDescript_t),SAVE :: Dvvdescript
    type(OutputDataDescript_t),SAVE :: Uudescript
    type(OutputDataDescript_t),SAVE :: Vvdescript
    type(OutputDataDescript_t),SAVE :: Bsxdescript
    type(OutputDataDescript_t),SAVE :: Bsydescript
!     for hotstart harmonic analysis
    type(OutputDataDescript_t), SAVE :: HarmElevFDLVDescript
    type(OutputDataDescript_t), SAVE :: HarmElevSLVDescript
    type(OutputDataDescript_t), SAVE :: HarmUVelFDLVDescript
    type(OutputDataDescript_t), SAVE :: HarmVVelFDLVDescript
    type(OutputDataDescript_t), SAVE :: HarmUVelSLVDescript
    type(OutputDataDescript_t), SAVE :: HarmVVelSLVDescript
!     for hotstart harmoinc analysis means and variance calculations
    type(OutputDataDescript_t), SAVE :: ELAVDescript
    type(OutputDataDescript_t), SAVE :: ELVADescript
    type(OutputDataDescript_t), SAVE :: XVELAVDescript
    type(OutputDataDescript_t), SAVE :: YVELAVDescript
    type(OutputDataDescript_t), SAVE :: XVELVADescript
    type(OutputDataDescript_t), SAVE :: YVELVADescript

    INTEGER :: numHotstartWrites ! number writes to hot start files
    INTEGER :: nextLun           ! next LUN to write to, after initial write

!     fort.61
    ElevStaDescript % specifier            = NOUTE
    ElevStaDescript % lun                  = 61
    ElevStaDescript % initial_value        = 0.0
    ElevStaDescript % num_items_per_record = 1
    ElevStaDescript % num_fd_records       = abs(NSTAE)
    ElevStaDescript % num_records_this     = abs(NSTAE)
    ElevStaDescript % ConsiderWetDry       = .TRUE. 
    ElevStaDescript % alternate_value      = -99999.0
    ElevStaDescript % field_name           = 'ElevSta'
    IF (ICS == 2) THEN
        ElevStaDescript % x_coord           => SLEL
        ElevStaDescript % y_coord           => SFEL
    ELSE
        ElevStaDescript % x_coord           => XEL
        ElevStaDescript % y_coord           => YEL
    ENDIF
    ElevStaDescript % file_extension       = 61
    ElevStaDescript % file_basename        = 'fort'
    ElevStaDescript % readMaxMin           = .FALSE. 
    call makeFileName(ElevStaDescript)

!     fort.62
    VelStaDescript % specifier            = NOUTV
    VelStaDescript % lun                  = 62
    VelStaDescript % initial_value        = 0.0
    VelStaDescript % num_items_per_record = 2
    VelStaDescript % num_fd_records       = abs(NSTAV)
    VelStaDescript % num_records_this     = abs(NSTAV)
    VelStaDescript % ConsiderWetDry       = .FALSE. 
    VelStaDescript % alternate_value      = 0.0
    VelStaDescript % field_name           = 'VelSta'
    IF (ICS == 2) THEN
        VelStaDescript % x_coord           => SLEV
        VelStaDescript % y_coord           => SFEV
    ELSE
        VelStaDescript % x_coord           => XEV
        VelStaDescript % y_coord           => YEV
    ENDIF
    VelStaDescript % file_extension       = 62
    VelStaDescript % file_basename        = 'fort'
    VelStaDescript % readMaxMin           = .FALSE. 
    call makeFileName(VelStaDescript)

!     fort.63
    ElevDescript % specifier            = NOUTGE
    ElevDescript % lun                  = 63
    ElevDescript % initial_value        = 0.0
    ElevDescript % num_items_per_record = 1
    ElevDescript % num_fd_records       = MNP
    ElevDescript % num_records_this     = MNP
    ElevDescript % ConsiderWetDry       = .TRUE. 
    ElevDescript % alternate_value      = -99999.0
    ElevDescript % field_name           = 'Elev'
    ElevDescript % file_extension       = 63
    ElevDescript % file_basename        = 'fort'
    ElevDescript % readMaxMin           = .FALSE. 
    call makeFileName(ElevDescript)

! fort.90 (tau0)
    Tau0Descript % lun                  = 90
    Tau0Descript % specifier            = NOUTGE
    Tau0Descript % initial_value        = 0.d0
    Tau0Descript % num_fd_records       = MNP
    Tau0Descript % num_records_this     = MNP
    Tau0Descript % ConsiderWetDry       = .FALSE. 
    Tau0Descript % alternate_value      = -99999.0
    Tau0Descript % field_name           = 'Tau0'
    Tau0Descript % file_extension       = 90
    Tau0Descript % file_basename        = 'fort'
    Tau0Descript % readMaxMin            = .FALSE. 
    call makeFileName(Tau0Descript)

!     fort.64
    VelDescript % specifier            = NOUTGV
    VelDescript % lun                  = 64
    VelDescript % initial_value        = 0.0
    VelDescript % num_items_per_record = 2
    VelDescript % num_fd_records       = MNP
    VelDescript % num_records_this     = MNP
    VelDescript % ConsiderWetDry       = .FALSE. 
    VelDescript % alternate_value      = 0.0
    VelDescript % field_name           = 'Vel'
    VelDescript % file_extension       = 64
    VelDescript % file_basename        = 'fort'
    VelDescript % readMaxMin           = .FALSE. 
    call makeFileName(VelDescript)

!     maxele.63
!     jgf52.08.11: Take absolute value of specifier to avoid overwriting
!     (reinitializing) existing min/max file, if any.
    EtaMaxDescript % specifier            = abs(NOUTGE)
    EtaMaxDescript % lun                  = 311
    EtaMaxDescript % initial_value        = -99999.d0
    EtaMaxDescript % num_items_per_record = 1
    EtaMaxDescript % num_fd_records       = MNP
    EtaMaxDescript % num_records_this     = MNP
    EtaMaxDescript % ConsiderWetDry       = .FALSE. 
    EtaMaxDescript % alternate_value      = -99999.d0
    EtaMaxDescript % field_name           = 'EtaMax'
    EtaMaxDescript % file_extension       = 63
    EtaMaxDescript % file_basename        = 'maxele'
    EtaMaxDescript % readMaxMin           = .TRUE. 
    call makeFileName(EtaMaxDescript)

!     maxvel.63
    UMaxDescript % specifier            = abs(NOUTGV)
    UMaxDescript % lun                  = 312
    UMaxDescript % initial_value        = 0.0
    UMaxDescript % num_items_per_record = 2
    UMaxDescript % num_fd_records       = MNP
    UMaxDescript % num_records_this     = MNP
    UMaxDescript % ConsiderWetDry       = .FALSE. 
    UMaxDescript % alternate_value      = 0.0
    UMaxDescript % field_name           = 'UMax'
    UMaxDescript % file_extension       = 63
    UMaxDescript % file_basename        = 'maxvel'
    UMaxDescript % readMaxMin           = .TRUE. 
    call makeFileName(UMaxDescript)

!     fort.71
    PrStaDescript % specifier            = NOUTM
    PrStaDescript % lun                  = 71
    PrStaDescript % initial_value        = 0.0
    PrStaDescript % num_items_per_record = 1
    PrStaDescript % num_fd_records       = abs(NSTAM)
    PrStaDescript % num_records_this     = abs(NSTAM)
    PrStaDescript % ConsiderWetDry       = .FALSE. 
    PrStaDescript % alternate_value      = 0.0
    PrStaDescript % field_name           = 'PrSta'
    IF (ICS == 2) THEN
        PrStaDescript % x_coord           => SLEM
        PrStaDescript % y_coord           => SFEM
    ELSE
        PrStaDescript % x_coord           => XEM
        PrStaDescript % y_coord           => YEM
    ENDIF
    PrStaDescript % file_extension       = 71
    PrStaDescript % file_basename        = 'fort'
    PrStaDescript % readMaxMin           = .FALSE. 
    call makeFileName(PrStaDescript)

!     fort.72
    WindVelStaDescript % specifier            = NOUTM
    WindVelStaDescript % lun                  = 72
    WindVelStaDescript % initial_value        = 0.0
    WindVelStaDescript % num_items_per_record = 2
    WindVelStaDescript % num_fd_records       = abs(NSTAM)
    WindVelStaDescript % num_records_this     = abs(NSTAM)
    WindVelStaDescript % ConsiderWetDry       = .FALSE. 
    WindVelStaDescript % alternate_value      = 0.0
    WindVelStaDescript % field_name           = 'WindVelSta'
    IF (ICS == 2) THEN
        WindVelStaDescript % x_coord           => SLEM
        WindVelStaDescript % y_coord           => SFEM
    ELSE
        WindVelStaDescript % x_coord           => XEM
        WindVelStaDescript % y_coord           => YEM
    ENDIF
    WindVelStaDescript % file_extension       = 72
    WindVelStaDescript % file_basename        = 'fort'
    WindVelStaDescript % readMaxMin           = .FALSE. 
    call makeFileName(WindVelStaDescript)

!     fort.73
    PrDescript % specifier            = NOUTGW
    PrDescript % lun                  = 73
    PrDescript % initial_value        = 0.0
    PrDescript % num_items_per_record = 1
    PrDescript % num_fd_records       = MNP
    PrDescript % num_records_this     = MNP
    PrDescript % ConsiderWetDry       = .FALSE. 
    PrDescript % alternate_value      = 0.0
    PrDescript % field_name           = 'Pr'
    PrDescript % file_extension       = 73
    PrDescript % file_basename        = 'fort'
    PrDescript % readMaxMin           = .FALSE. 
    call makeFileName(PrDescript)

!     fort.74
    WindVelDescript % specifier            = NOUTGW
    WindVelDescript % lun                  = 74
    WindVelDescript % initial_value        = 0.0
    WindVelDescript % num_items_per_record = 2
    WindVelDescript % num_fd_records       = MNP
    WindVelDescript % num_records_this     = MNP
    WindVelDescript % ConsiderWetDry       = .FALSE. 
    WindVelDescript % alternate_value      = 0.0
    WindVelDescript % field_name           = 'WindVel'
    WindVelDescript % file_extension       = 74
    WindVelDescript % file_basename        = 'fort'
    WindVelDescript % readMaxMin           = .FALSE. 
    call makeFileName(WindVelDescript)

!     fort.77
    WeirElevDescript % specifier            = NOUT_TVW
    WeirElevDescript % lun                  = 77
    WeirElevDescript % initial_value        = 0.0
    WeirElevDescript % num_items_per_record = 1
    WeirElevDescript % num_fd_records       = MNP
    WeirElevDescript % num_records_this     = MNP
    WeirElevDescript % ConsiderWetDry       = .FALSE. 
    WeirElevDescript % alternate_value      = 0.0
    WeirElevDescript % field_name           = 'weir_dz'
    WeirElevDescript % file_extension       = 77
    WeirElevDescript % file_basename        = 'fort'
    call makeFileName(weirElevDescript)


!     minpr.63
    PrMinDescript % specifier            = abs(NOUTGW)
    PrMinDescript % lun                  = 313
    PrMinDescript % initial_value        = 99999.d0
    PrMinDescript % num_items_per_record = 1
    PrMinDescript % num_fd_records       = MNP
    PrMinDescript % num_records_this     = MNP
    PrMinDescript % ConsiderWetDry       = .FALSE. 
    PrMinDescript % alternate_value      = 0.0
    PrMinDescript % field_name           = 'PrMin'
    PrMinDescript % file_extension       = 63
    PrMinDescript % file_basename        = 'minpr'
    PrMinDescript % readMaxMin           = .TRUE. 
    call makeFileName(PrMinDescript)

!     maxwvel.63
    WVMaxDescript % specifier            = abs(NOUTGW)
    WVMaxDescript % lun                  = 314
    WVMaxDescript % initial_value        = 0.0
    WVMaxDescript % num_items_per_record = 2
    WVMaxDescript % num_fd_records       = MNP
    WVMaxDescript % num_records_this     = MNP
    WVMaxDescript % ConsiderWetDry       = .FALSE. 
    WVMaxDescript % alternate_value      = 0.0
    WVMaxDescript % field_name           = 'WVMax'
    WVMaxDescript % file_extension       = 63
    WVMaxDescript % file_basename        = 'maxwvel'
    WVMaxDescript % readMaxMin           = .TRUE. 
    call makeFileName(WVMaxDescript)

    RSMaxDescript % specifier            = abs(NOUTGW)
    RSMaxDescript % lun                  = 315
    RSMaxDescript % initial_value        = 0.0
    RSMaxDescript % num_items_per_record = 2
    RSMaxDescript % num_fd_records       = MNP
    RSMaxDescript % num_records_this     = MNP
    RSMaxDescript % ConsiderWetDry       = .FALSE. 
    RSMaxDescript % alternate_value      = 0.0
    RSMaxDescript % field_name           = "RSMax"
    RSMaxDescript % file_name            = "maxrs.63"
    RSMaxDescript % file_extension       = 63
    RSMaxDescript % file_basename        = 'maxrs'
    RSMaxDescript % readMaxMin           = .TRUE. 
    call makeFileName(RSMaxDescript)

!  tcm v50.75 removed ifdef adcswan to allow for use whenever nrs=3 or nrs=4
!#ifdef ADCSWAN
! Cobell 20120510: SWAN Output Data
!........Radiation Stress
      RSDescript % specifier            = NOUTGW
      RSDescript % lun                  = 164
      RSDescript % initial_value        = 0.0
      RSDescript % num_items_per_record = 2
      RSDescript % num_fd_records       = MNP
      RSDescript % num_records_this     = MNP
      RSDescript % ConsiderWetDry       = .FALSE.
      RSDescript % alternate_value      = -99999.0
      RSDescript % field_name           = "rads"
      RSDescript % file_name            = "rads.64"
      RSDescript % file_extension       = 64
      RSDescript % file_basename        = 'rads'
      RSDescript % readMaxMin           = .false.
      call makeFileName(RSDescript)
      !
      !  D E T A I L E D    I N U N D A T I O N    O U T P U T
      !
      ! inundationtime.63 (works like a min/max file)
      InundationTimeDescript % lun                  = 400
      InundationTimeDescript % specifier            = abs(NOUTGE)
      InundationTimeDescript % num_items_per_record = 1
      InundationTimeDescript % initial_value        = 0.0      
      InundationTimeDescript % num_fd_records       = MNP
      InundationTimeDescript % num_records_this     = MNP  
      InundationTimeDescript % ConsiderWetDry       = .false.
      InundationTimeDescript % alternate_value      = 0.0
      InundationTimeDescript % field_name           = 'inundationTime'
      InundationTimeDescript % file_basename        = 'inundationtime'
      InundationTimeDescript % file_extension       = 63
      InundationTimeDescript % readMaxMin           = .true.
      if ( InundationTimeDescript % specifier .eq. XDMF ) then
         InundationTimeDescript % specifier = ASCII
      endif
      call makeFileName(InundationTimeDescript)
      !
      ! maxinundepth.63 (works like a min/max file)
      MaxInunDepthDescript % lun                  = 401
      MaxInunDepthDescript % specifier            = abs(NOUTGE)
      MaxInunDepthDescript % initial_value        = 0.0
      MaxInunDepthDescript % num_fd_records       = MNP
      MaxInunDepthDescript % num_records_this     = MNP 
      MaxInunDepthDescript % ConsiderWetDry       = .false.
      MaxInunDepthDescript % alternate_value      = 0.0 
      MaxInunDepthDescript % num_items_per_record = 1   
      MaxInunDepthDescript % field_name           = 'maxInunDepth'
      MaxInunDepthDescript % file_basename        = 'maxinundepth'
      MaxInunDepthDescript % file_extension       = 63
      MaxInunDepthDescript % readMaxMin           = .true.      
      if ( MaxInunDepthDescript % specifier .eq. XDMF ) then
         MaxInunDepthDescript % specifier = ASCII
      endif
      call makeFileName(MaxInunDepthDescript)      
      !
      ! initiallydry.63 (nodes that are deemed dry upon cold start)
      InitiallyDryDescript % lun                  = 402
      InitiallyDryDescript % specifier            = NOUTGE
      InitiallyDryDescript % initial_value        = 0.0
      InitiallyDryDescript % num_fd_records       = MNP
      InitiallyDryDescript % num_records_this     = MNP 
      InitiallyDryDescript % ConsiderWetDry       = .false.
      InitiallyDryDescript % alternate_value      = 0.0 
      InitiallyDryDescript % num_items_per_record = 1   
      InitiallyDryDescript % field_name           = 'initiallyDry'
      InitiallyDryDescript % file_basename        = 'initiallydry'
      InitiallyDryDescript % file_extension       = 63
      InitiallyDryDescript % isInteger            = .true.
      if ( InitiallyDryDescript % specifier .eq. XDMF ) then
         InitiallyDryDescript % specifier = ASCII
      endif
      call makeFileName(InitiallyDryDescript)      
      !
      ! endrisinginun.63 (1 if water is rising at end of simulation)
      EndRisingInunDescript % lun                  = 403
      EndRisingInunDescript % specifier            = NOUTGE
      EndRisingInunDescript % initial_value        = 0.0
      EndRisingInunDescript % num_fd_records       = MNP
      EndRisingInunDescript % num_records_this     = MNP 
      EndRisingInunDescript % ConsiderWetDry       = .false.
      EndRisingInunDescript % alternate_value      = 0.0 
      EndRisingInunDescript % num_items_per_record = 1        
      EndRisingInunDescript % field_name           = 'endRisingInun'
      EndRisingInunDescript % file_basename        = 'endrisinginun'
      EndRisingInunDescript % file_extension       = 63
      EndRisingInunDescript % isInteger            = .true.
      if ( EndRisingInunDescript % specifier .eq. XDMF ) then
         EndRisingInunDescript % specifier = ASCII
      endif
      call makeFileName(EndRisingInunDescript)
      ! 
      ! everdried.63 (works like a min/max file)
      EverDriedDescript % lun                  = 404
      EverDriedDescript % specifier            = abs(NOUTGE)
      EverDriedDescript % num_items_per_record = 1
      EverDriedDescript % initial_value        = 0.0      
      EverDriedDescript % num_fd_records       = MNP
      EverDriedDescript % num_records_this     = MNP  
      EverDriedDescript % ConsiderWetDry       = .false.
      EverDriedDescript % alternate_value      = -99999.d0
      EverDriedDescript % field_name           = 'EverDried'
      EverDriedDescript % file_basename        = 'everdried'
      EverDriedDescript % file_extension       = 63
      EverDriedDescript % readMaxMin           = .true.
      if ( EverDriedDescript % specifier .eq. XDMF ) then
         EverDriedDescript % specifier = ASCII
      endif
      call makeFileName(EverDriedDescript)

         weirElevDescript % specifier            = NOUT_TVW
         weirElevDescript % initial_value        = 0D0
         weirElevDescript % num_items_per_record = 1
         weirElevDescript % num_fd_records       = MNP
         weirElevDescript % num_records_this     = MNP
         weirElevDescript % ConsiderWetDry       = .FALSE.
         weirElevDescript % alternate_value      = 0D0
         weirElevDescript % field_name           = "weir_dz"
         weirElevDescript % file_name            = "fort.77"

! tcm v50.75 moved the ifdef adcswan down past the RSDescript
#ifdef ADCSWAN
!........Significant Wave Height (HS)
    SwanHSDescript % specifier            = NOUTGW
    SwanHSDescript % lun                  = 301
    SwanHSDescript % initial_value        = 0.0
    SwanHSDescript % num_items_per_record = 1
    SwanHSDescript % num_fd_records       = MNP
    SwanHSDescript % num_records_this     = MNP
    SwanHSDescript % ConsiderWetDry       = .FALSE. 
    SwanHSDescript % alternate_value      = -99999.0
    SwanHSDescript % field_name           = "swan_HS"
    SwanHSDescript % file_name            = "swan_HS.63"
    SwanHSDescript % file_extension       = 63
    SwanHSDescript % file_basename        = 'swan_HS'
    call makeFileName(SwanHSDescript)

    SwanHSMaxDescript % specifier            = NOUTGW
    SwanHSMaxDescript % lun                  = 316
    SwanHSMaxDescript % initial_value        = 0.0
    SwanHSMaxDescript % num_items_per_record = 1
    SwanHSMaxDescript % num_fd_records       = MNP
    SwanHSMaxDescript % num_records_this     = MNP
    SwanHSMaxDescript % ConsiderWetDry       = .FALSE. 
    SwanHSMaxDescript % alternate_value      = -99999.0
    SwanHSMaxDescript % field_name           = "swan_HS_max"
    SwanHSMaxDescript % file_name            = "swan_HS_max.63"
    SwanHSMaxDescript % file_extension       = 63
    SwanHSMaxDescript % file_basename        = 'swan_HS_max'
    call makeFileName(SwanHSMaxDescript)

!........Mean Wave Direction (DIR)
    SwanDIRDescript % specifier            = NOUTGW
    SwanDIRDescript % lun                  = 302
    SwanDIRDescript % initial_value        = 0.0
    SwanDIRDescript % num_items_per_record = 1
    SwanDIRDescript % num_fd_records       = MNP
    SwanDIRDescript % num_records_this     = MNP
    SwanDIRDescript % ConsiderWetDry       = .FALSE. 
    SwanDIRDescript % alternate_value      = -99999.0
    SwanDIRDescript % field_name           = "swan_DIR"
    SwanDIRDescript % file_name            = "swan_DIR.63"
    SwanDIRDescript % file_extension       = 63
    SwanDIRDescript % file_basename        = 'swan_DIR'
    call makeFileName(SwanDIRDescript)

    SwanDIRMaxDescript % specifier            = NOUTGW
    SwanDIRMaxDescript % lun                  = 317
    SwanDIRMaxDescript % initial_value        = 0.0
    SwanDIRMaxDescript % num_items_per_record = 1
    SwanDIRMaxDescript % num_fd_records       = MNP
    SwanDIRMaxDescript % num_records_this     = MNP
    SwanDIRMaxDescript % ConsiderWetDry       = .FALSE. 
    SwanDIRMaxDescript % alternate_value      = -99999.0
    SwanDIRMaxDescript % field_name           = "swan_DIR_max"
    SwanDIRMaxDescript % file_name            = "swan_DIR_max.63"
    SwanDIRMaxDescript % file_extension       = 63
    SwanDIRMaxDescript % file_basename        = 'swan_DIR_max'
    call makeFileName(SwanDIRMaxDescript)

!........Mean Wave Period (TM01)
    SwanTM01Descript % specifier            = NOUTGW
    SwanTM01Descript % lun                  = 303
    SwanTM01Descript % initial_value        = 0.0
    SwanTM01Descript % num_items_per_record = 1
    SwanTM01Descript % num_fd_records       = MNP
    SwanTM01Descript % num_records_this     = MNP
    SwanTM01Descript % ConsiderWetDry       = .FALSE. 
    SwanTM01Descript % alternate_value      = -99999.0
    SwanTM01Descript % field_name           = "swan_TM01"
    SwanTM01Descript % file_name            = "swan_TM01.63"
    SwanTM01Descript % file_extension       = 63
    SwanTM01Descript % file_basename        = 'swan_TM01'
    call makeFileName(SwanTM01Descript)

    SwanTM01MaxDescript % specifier            = NOUTGW
    SwanTM01MaxDescript % lun                  = 318
    SwanTM01MaxDescript % initial_value        = 0.0
    SwanTM01MaxDescript % num_items_per_record = 1
    SwanTM01MaxDescript % num_fd_records       = MNP
    SwanTM01MaxDescript % num_records_this     = MNP
    SwanTM01MaxDescript % ConsiderWetDry       = .FALSE. 
    SwanTM01MaxDescript % alternate_value      = -99999.0
    SwanTM01MaxDescript % field_name           = "swan_TM01_max"
    SwanTM01MaxDescript % file_name            = "swan_TM01_max.63"
    SwanTM01MaxDescript % file_extension       = 63
    SwanTM01MaxDescript % file_basename        = 'swan_TM01_max'
    call makeFileName(SwanTM01MaxDescript)

!........Peak Wave Period (TPS)
    SwanTPSDescript % specifier            = NOUTGW
    SwanTPSDescript % lun                  = 304
    SwanTPSDescript % initial_value        = 0.0
    SwanTPSDescript % num_items_per_record = 1
    SwanTPSDescript % num_fd_records       = MNP
    SwanTPSDescript % num_records_this     = MNP
    SwanTPSDescript % ConsiderWetDry       = .FALSE. 
    SwanTPSDescript % alternate_value      = -99999.0
    SwanTPSDescript % field_name           = "Swan_TPS"
    SwanTPSDescript % file_name            = "swan_TPS.63"
    SwanTPSDescript % file_extension       = 63
    SwanTPSDescript % file_basename        = 'swan_TPS'
    call makeFileName(SwanTPSDescript)

    SwanTPSMaxDescript % specifier            = NOUTGW
    SwanTPSMaxDescript % lun                  = 319
    SwanTPSMaxDescript % initial_value        = 0.0
    SwanTPSMaxDescript % num_items_per_record = 1
    SwanTPSMaxDescript % num_fd_records       = MNP
    SwanTPSMaxDescript % num_records_this     = MNP
    SwanTPSMaxDescript % ConsiderWetDry       = .FALSE. 
    SwanTPSMaxDescript % alternate_value      = -99999.0
    SwanTPSMaxDescript % field_name           = "Swan_TPS_max"
    SwanTPSMaxDescript % file_name            = "swan_TPS_max.63"
    SwanTPSMaxDescript % file_extension       = 63
    SwanTPSMaxDescript % file_basename        = 'swan_TPS_max'
    call makeFileName(SwanTPSMaxDescript)

!........SWAN Wind Values (WINDX,WINDY)
    SwanWindDescript % specifier            = NOUTGW
    SwanWindDescript % lun                  = 305
    SwanWindDescript % initial_value        = 0.0
    SwanWindDescript % num_items_per_record = 2
    SwanWindDescript % num_fd_records       = MNP
    SwanWindDescript % num_records_this     = MNP
    SwanWindDescript % ConsiderWetDry       = .FALSE. 
    SwanWindDescript % alternate_value      = -99999.0
    SwanWindDescript % field_name           = "swan_WIND"
    SwanWindDescript % file_name            = "swan_WIND.64"
    SwanWindDescript % file_extension       = 64
    SwanWindDescript % file_basename        = 'swan_WIND'
    call makeFileName(SwanWindDescript)

    SwanWindMaxDescript % specifier            = NOUTGW
    SwanWindMaxDescript % lun                  = 320
    SwanWindMaxDescript % initial_value        = 0.0
    SwanWindMaxDescript % num_items_per_record = 1
    SwanWindMaxDescript % num_fd_records       = MNP
    SwanWindMaxDescript % num_records_this     = MNP
    SwanWindMaxDescript % ConsiderWetDry       = .FALSE. 
    SwanWindMaxDescript % alternate_value      = -99999.0
    SwanWindMaxDescript % field_name           = "swan_WIND_max"
    SwanWindMaxDescript % file_name            = "swan_WIND_max.63"
    SwanWindMaxDescript % file_extension       = 63
    SwanWindMaxDescript % file_basename        = 'swan_WIND_max'
    call makeFileName(SwanWindMaxDescript)

!........Mean Wave Period (TM02)
    SwanTM02Descript % specifier            = NOUTGW
    SwanTM02Descript % lun                  = 306
    SwanTM02Descript % initial_value        = 0.0
    SwanTM02Descript % num_items_per_record = 1
    SwanTM02Descript % num_fd_records       = MNP
    SwanTM02Descript % num_records_this     = MNP
    SwanTM02Descript % ConsiderWetDry       = .FALSE. 
    SwanTM02Descript % alternate_value      = -99999.0
    SwanTM02Descript % field_name           = "swan_TM02"
    SwanTM02Descript % file_name            = "swan_TM02.63"
    SwanTM02Descript % file_extension       = 63
    SwanTM02Descript % file_basename        = 'swan_TM02'
    call makeFileName(SwanTM02Descript)

    SwanTM02MaxDescript % specifier            = NOUTGW
    SwanTM02MaxDescript % lun                  = 321
    SwanTM02MaxDescript % initial_value        = 0.0
    SwanTM02MaxDescript % num_items_per_record = 1
    SwanTM02MaxDescript % num_fd_records       = MNP
    SwanTM02MaxDescript % num_records_this     = MNP
    SwanTM02MaxDescript % ConsiderWetDry       = .FALSE. 
    SwanTM02MaxDescript % alternate_value      = -99999.0
    SwanTM02MaxDescript % field_name           = "swan_TM02_max"
    SwanTM02MaxDescript % file_name            = "swan_TM02_max.63"
    SwanTM02MaxDescript % file_extension       = 63
    SwanTM02MaxDescript % file_basename        = 'swan_TM02_max'
    call makeFileName(SwanTM02MaxDescript)

!........Mean Wave Period (TMM10)
    SwanTMM10Descript % specifier            = NOUTGW
    SwanTMM10Descript % lun                  = 307
    SwanTMM10Descript % initial_value        = 0.0
    SwanTMM10Descript % num_items_per_record = 1
    SwanTMM10Descript % num_fd_records       = MNP
    SwanTMM10Descript % num_records_this     = MNP
    SwanTMM10Descript % ConsiderWetDry       = .FALSE. 
    SwanTMM10Descript % alternate_value      = -99999.0
    SwanTMM10Descript % field_name           = "swan_TMM10"
    SwanTMM10Descript % file_name            = "swan_TMM10.63"
    SwanTMM10Descript % file_extension       = 63
    SwanTMM10Descript % file_basename        = 'swan_TMM10'
    call makeFileName(SwanTMM10Descript)

    SwanTMM10MaxDescript % specifier            = NOUTGW
    SwanTMM10MaxDescript % lun                  = 322
    SwanTMM10MaxDescript % initial_value        = 0.0
    SwanTMM10MaxDescript % num_items_per_record = 1
    SwanTMM10MaxDescript % num_fd_records       = MNP
    SwanTMM10MaxDescript % num_records_this     = MNP
    SwanTMM10MaxDescript % ConsiderWetDry       = .FALSE. 
    SwanTMM10MaxDescript % alternate_value      = -99999.0
    SwanTMM10MaxDescript % field_name           = "swan_TMM10_max"
    SwanTMM10MaxDescript % file_name          = "swan_TMM10_max.63"
    SwanTMM10MaxDescript % file_extension       = 63
    SwanTMM10MaxDescript % file_basename        = 'swan_TMM10_max'
    call makeFileName(SwanTMM10MaxDescript)

#endif
!     3D data
!     fort.41
    SigTStaDescript % specifier            = I3DSD
    SigTStaDescript % lun                  = 41
    SigTStaDescript % initial_value        = 0.0
    SigTStaDescript % num_items_per_record = NFEN
    SigTStaDescript % num_fd_records       = NSTA3DD
    SigTStaDescript % num_records_this     = NSTA3DD
    SigTStaDescript % field_name           = 'SigmaTStations'
    IF (ICS == 2) THEN
        SigTStaDescript % x_coord           => SL3DD
        SigTStaDescript % y_coord           => SF3DD
    ELSE
        SigTStaDescript % x_coord           => X3DD
        SigTStaDescript % y_coord           => Y3DD
    ENDIF
    SigTStaDescript % file_extension       = 41
    SigTStaDescript % file_basename        = 'fort'
    call makeFileName(SigTStaDescript)

    SalStaDescript % specifier            = I3DSD
    SalStaDescript % lun                  = 41
    SalStaDescript % initial_value        = 0.0
    SalStaDescript % num_items_per_record = NFEN
    SalStaDescript % num_fd_records       = NSTA3DD
    SalStaDescript % num_records_this     = NSTA3DD
    SalStaDescript % field_name           = 'SalinityStations'
    IF (ICS == 2) THEN
        SalStaDescript % x_coord           => SL3DD
        SalStaDescript % y_coord           => SF3DD
    ELSE
        SalStaDescript % x_coord           => X3DD
        SalStaDescript % y_coord           => Y3DD
    ENDIF

    TempStaDescript % specifier            = I3DSD
    TempStaDescript % lun                  = 41
    TempStaDescript % initial_value        = 0.0
    TempStaDescript % num_items_per_record = NFEN
    TempStaDescript % num_fd_records       = NSTA3DD
    TempStaDescript % num_records_this     = NSTA3DD
    TempStaDescript % field_name           = 'TemperatureStations'
    IF (ICS == 2) THEN
        TempStaDescript % x_coord           => SL3DD
        TempStaDescript % y_coord           => SF3DD
    ELSE
        TempStaDescript % x_coord           => X3DD
        TempStaDescript % y_coord           => Y3DD
    ENDIF

!     fort.42
    RealQStaDescript % specifier            =  I3DSV
    RealQStaDescript % lun                  =  42
    RealQStaDescript % initial_value        =  0.0
    RealQStaDescript % num_items_per_record =  NFEN
    RealQStaDescript % num_fd_records       =  NSta3DV
    RealQStaDescript % num_records_this     =  NSta3DV
    RealQStaDescript % field_name           = 'RealQStations'
    IF (ICS == 2) THEN
        RealQStaDescript % x_coord           => SL3DV
        RealQStaDescript % y_coord           => SF3DV
    ELSE
        RealQStaDescript % x_coord           => X3DV
        RealQStaDescript % y_coord           => Y3DV
    ENDIF
    RealQStaDescript % file_extension       = 42
    RealQStaDescript % file_basename        = 'fort'
    call makeFileName(RealQStaDescript)

    ImaginaryQStaDescript % specifier            =  I3DSV
    ImaginaryQStaDescript % lun                  =  42
    ImaginaryQStaDescript % initial_value        =  0.0
    ImaginaryQStaDescript % num_items_per_record =  NFEN
    ImaginaryQStaDescript % num_fd_records       =  NSTA3DV
    ImaginaryQStaDescript % num_records_this     =  NSTA3DV
    ImaginaryQStaDescript % field_name      ='ImaginaryQStations'
    IF (ICS == 2) THEN
        ImaginaryQStaDescript % x_coord           => SL3DV
        ImaginaryQStaDescript % y_coord           => SF3DV
    ELSE
        ImaginaryQStaDescript % x_coord           => X3DV
        ImaginaryQStaDescript % y_coord           => Y3DV
    ENDIF

    WZStaDescript % specifier            =  I3DSV
    WZStaDescript % lun                  =  42
    WZStaDescript % initial_value        =  0.0
    WZStaDescript % num_items_per_record =  NFEN
    WZStaDescript % num_fd_records       =  NSTA3DV
    WZStaDescript % num_records_this     =  NSTA3DV
    WZStaDescript % field_name           = 'WZStations'
    IF (ICS == 2) THEN
        WZStaDescript % x_coord           => SL3DV
        WZStaDescript % y_coord           => SF3DV
    ELSE
        WZStaDescript % x_coord           => X3DV
        WZStaDescript % y_coord           => Y3DV
    ENDIF

!     fort.43
    Q20StaDescript % specifier            =  I3DST
    Q20StaDescript % lun                  =  43
    Q20StaDescript % initial_value        =  0.0
    Q20StaDescript % num_items_per_record =  NFEN
    Q20StaDescript % num_fd_records       =  NSta3DT
    Q20StaDescript % num_records_this     =  NSta3DT
    Q20StaDescript % field_name           = 'q20Stations'
    IF (ICS == 2) THEN
        Q20StaDescript % x_coord           => SL3DT
        Q20StaDescript % y_coord           => SF3DT
    ELSE
        Q20StaDescript % x_coord           => X3DT
        Q20StaDescript % y_coord           => Y3DT
    ENDIF
    Q20StaDescript % file_extension       = 43
    Q20StaDescript % file_basename        = 'fort'
    call makeFileName(Q20StaDescript)

    LStaDescript % specifier            =  I3DST
    LStaDescript % lun                  =  43
    LStaDescript % initial_value        =  0.0
    LStaDescript % num_items_per_record =  NFEN
    LStaDescript % num_fd_records       =  NSTA3DT
    LStaDescript % num_records_this     =  NSTA3DT
    LStaDescript % field_name           = 'LStations'
    IF (ICS == 2) THEN
        LStaDescript % x_coord           => SL3DT
        LStaDescript % y_coord           => SF3DT
    ELSE
        LStaDescript % x_coord           => X3DT
        LStaDescript % y_coord           => Y3DT
    ENDIF

    EVStaDescript % specifier            =  I3DST
    EVStaDescript % lun                  =  43
    EVStaDescript % initial_value        =  0.0
    EVStaDescript % num_items_per_record =  NFEN
    EVStaDescript % num_fd_records       =  NSTA3DT
    EVStaDescript % num_records_this     =  NSTA3DT
    EVStaDescript % field_name           = 'EVStations'
    IF (ICS == 2) THEN
        EVStaDescript % x_coord           => SL3DT
        EVStaDescript % y_coord           => SF3DT
    ELSE
        EVStaDescript % x_coord           => X3DT
        EVStaDescript % y_coord           => Y3DT
    ENDIF

!     fort.44
    SigTDescript % specifier            =  I3DGD
    SigTDescript % lun                  =  44
    SigTDescript % initial_value        =  0.0
    SigTDescript % num_items_per_record =  NFEN
    SigTDescript % num_fd_records        =  MNP
    SigTDescript % num_records_this     =  MNP
    SigTDescript % field_name           = 'SigmaT'
    SigTDescript % file_extension       = 44
    SigTDescript % file_basename        = 'fort'
    call makeFileName(SigTDescript)

    SalDescript % specifier            =  I3DGD
    SalDescript % lun                  =  44
    SalDescript % initial_value        =  0.0
    SalDescript % num_items_per_record =  NFEN
    SalDescript % num_fd_records       =  MNP
    SalDescript % num_records_this     =  MNP
    SalDescript % field_name           = 'Salinity'
    TempDescript % specifier            =  I3DGD
    TempDescript % lun                  =  44
    TempDescript % initial_value        =  0.0
    TempDescript % num_items_per_record =  NFEN
    TempDescript % num_fd_records       =  MNP
    TempDescript % num_records_this     =  MNP
    TempDescript % field_name           = 'Temperature'

!     fort.45
    RealQDescript % specifier            =  I3DGV
    RealQDescript % lun                  =  45
    RealQdescript % initial_value        =  0.0
    RealQDescript % num_items_per_record =  NFEN
    RealQDescript % num_fd_records       =  MNP
    RealQDescript % num_records_this     =  MNP
    RealQDescript % field_name           = 'RealQ'
    RealQDescript % file_extension       = 45
    RealQDescript % file_basename        = 'fort'
    call makeFileName(RealQDescript)

    ImaginaryQDescript % specifier            =  I3DGV
    ImaginaryQDescript % lun                  =  45
    ImaginaryQDescript % initial_value        =  0.0
    ImaginaryQDescript % num_items_per_record =  NFEN
    ImaginaryQDescript % num_fd_records       =  MNP
    ImaginaryQDescript % num_records_this     =  MNP
    ImaginaryQDescript % field_name           = 'ImaginaryQ'
    WZDescript % specifier            =  I3DGV
    WZDescript % lun                  =  45
    WZDescript % initial_value        =  0.0
    WZDescript % num_items_per_record =  NFEN
    WZDescript % num_fd_records       =  MNP
    WZDescript % num_records_this     =  MNP
    WZDescript % field_name           = 'WZ'

!     fort.46
    Q20Descript % specifier            =  I3DGT
    Q20Descript % lun                  =  46
    Q20Descript % initial_value        =  0.0
    Q20Descript % num_items_per_record =  NFEN
    Q20Descript % num_fd_records       =  MNP
    Q20Descript % num_records_this     =  MNP
    Q20Descript % field_name           = 'q20'
    Q20Descript % file_extension       = 46
    Q20Descript % file_basename        = 'fort'
    call makeFileName(Q20Descript)

    LDescript % specifier            =  I3DGT
    LDescript % lun                  =  46
    LDescript % initial_value        =  0.0
    LDescript % num_items_per_record =  NFEN
    LDescript % num_fd_records       =  MNP
    LDescript % num_records_this     =  MNP
    LDescript % field_name           = 'L'
    EVDescript % specifier            =  I3DGT
    EVDescript % lun                  =  46
    EVDescript % initial_value        =  0.0
    EVDescript % num_items_per_record =  NFEN
    EVDescript % num_fd_records       =  MNP
    EVDescript % num_records_this     =  MNP
    EVDescript % field_name           = 'EV'

!     fort.47
    QSurfKp1Descript % specifier            =  I3DGD
    QSurfKp1Descript % lun                  =  47
    QSurfKp1Descript % initial_value        =  0.0
    QSurfKp1Descript % num_items_per_record =  1
    QSurfKp1Descript % num_fd_records     =  MNP
    QSurfKp1Descript % num_records_this   =  MNP
    QSurfKp1Descript % field_name           = 'qsurfkp1'
    call makeFileName(QSurfKp1Descript)

!     fort.67 and fort.68
    Elev1Descript % specifier            = NHSTAR
    Elev1Descript % initial_value        = 0.0
    Elev1Descript % file_basename        = 'fort'
    Elev1Descript % file_extension       = 67
    call makeFileName(Elev1Descript)
    Elev2Descript % specifier            = NHSTAR
    Elev2Descript % initial_value        = 0.0
    CH1Descript % specifier            = NHSTAR
    CH1Descript % initial_value        = 0.0
    EtaDiscDescript % specifier            = NHSTAR
    EtaDiscDescript % initial_value        = 0.0
    NodeCodeDescript % specifier            = NHSTAR
    NodeCodeDescript % initial_value        = 0.d0
    NOFFDescript % specifier            = NHSTAR
    NOFFDescript % initial_value        = 0.d0
!     hotstart 3D
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
!     hotstart harmonic analysis
    HarmElevFDLVDescript % specifier            = NHSTAR
    HarmElevFDLVDescript % initial_value        = 0.0
    HarmElevFDLVDescript % num_fd_records       = MNP
    HarmElevSLVDescript % specifier            = NHSTAR
    HarmElevSLVDescript % initial_value        = 0.0
    HarmElevSLVDescript % num_fd_records       = abs(NSTAE)
    HarmUVelFDLVDescript % specifier            = NHSTAR
    HarmUVelFDLVDescript % initial_value        = 0.0
    HarmUVelFDLVDescript % num_fd_records       = MNP
    HarmVVelFDLVDescript % specifier            = NHSTAR
    HarmVVelFDLVDescript % initial_value        = 0.0
    HarmVVelFDLVDescript % num_fd_records       = MNP
    HarmUvelSLVDescript % specifier            = NHSTAR
    HarmUVelSLVDescript % initial_value        = 0.0
    HarmUVelSLVDescript % num_fd_records       = abs(NSTAV)
    HarmVVelSLVDescript % specifier            = NHSTAR
    HarmVVelSLVDescript % initial_value        = 0.0
    HarmVVelSLVDescript % num_fd_records       = abs(NSTAV)
!     hotstart means and variance calculations
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

!     Need to populate the global and nodal attributes modules with
!     these parameters, since the netcdfio module relies on those
!     modules, rather than the pre_global module. Some day, adcprep
!     will be integrated with ADCIRC and this subroutine call will
!     not be needed.
    IF ( .NOT. ALLOCATED(NODECODE)) THEN
        ALLOCATE(NODECODE(MNP))
    ENDIF

!     jgf49.44: Set parameters in global module based on the data
!     we collected in read_global.F and stored in pre_global.F.
    CALL setADCIRCParameters( &
    base_date, MNE, NBOU, &
    NVEL, NOPE, MNP, SL0, SF0, NBVV, NVDLL, NBDV, NVELL, X, Y, &
    IBTYPE, IBTYPEE, SL1, SF1, NODECODE, G, FileFmtRev, &
    FileFmtMinor, FileFmtMajor, im, iestp, nscoue, ivstp, nscouv, &
    icstp, nscouc, ipstp, iwstp, nscoum, igep, nscouge, igvp, &
    nscougv, igcp, nscougc, igpp, igwp, nscougw, NM, &
    DP, RUNDES, AGRID, title, institution, source, history, &
    references, comments, host, convention, contact, DT, ihot, &
    ics, nolifa, nolica, nolicat, ncor, ntip, nws, nramp, statim, &
    reftim, rnday, dramp, a00, b00, c00, h0, cori, ntif, nbfr, &
    myProc, screenUnit, nolibf, nwp, tau0, cf, eslm, &
    abs(nstae), abs(nstav), abs(nstam), neta, nabout, nscreen, &
    nfen, iden, islip, kp, z0s, z0b, theta1, theta2, &
    ievc, evmin, evcon, alp1, alp2, alp3, igc, nlsd, nvsd, nltd, &
    nvtd, alp4, C3D, runid)

!     Create NetCDF output files for those output files where NetCDF
!     was specified.
    reterr = .FALSE. 
    CALL initNetCDFOutputFile(ElevStaDescript, reterr)
    CALL initNetCDFOutputFile(VelStaDescript, reterr)
    CALL initNetCDFOutputFile(ElevDescript, reterr)
    if (outputTau0.eqv. .TRUE. ) then
        CALL initNetCDFOutputFile(Tau0Descript, reterr)
    endif
    CALL initNetCDFOutputFile(VelDescript, reterr)
    CALL initNetCDFOutputFile(PrStaDescript, reterr)
    CALL initNetCDFOutputFile(WindVelStaDescript, reterr)
    CALL initNetCDFOutputFile(PrDescript, reterr)
    CALL initNetCDFOutputFile(WindVelDescript, reterr)
    CALL initNetCDFOutputFile(WeirElevDescript, reterr)
    CALL initNetCDFOutputFile(EtaMaxDescript, reterr)
    CALL initNetCDFOutputFile(UMaxDescript, reterr)
    CALL initNetCDFOutputFile(PrMinDescript, reterr)
    CALL initNetCDFOutputFile(WVMaxDescript, reterr)
    CALL initNetCDFOutputFile(RSMaxDescript,reterr)
    if (inundationOutput.eqv. .TRUE. ) then
        CALL initNetCDFOutputFile(InundationTimeDescript,reterr)
        CALL initNetCDFOutputFile(MaxInunDepthDescript,reterr)
        CALL initNetCDFOutputFile(InitiallyDryDescript,reterr)
        CALL initNetCDFOutputFile(EndRisingInunDescript,reterr)
        CALL initNetCDFOutputFile(EverDriedDescript,reterr)
    endif

! tcm v50.75 moved ifdef adcswan below RSDescript only to allow
! for use whenever nrs=3 or nrs=4
! Cobell 20120510: Added for SWAN NetCDF
    IF ((NRS == 3) .OR. (NRS == 4)) THEN
        CALL initNetCDFOutputFile(RSDescript,reterr)
    ENDIF
! tcm v50.75 moved ifdef adcswan to here
#ifdef ADCSWAN
! Cobell 20120510: Added for SWAN NetCDF
    IF(NRS == 3)THEN
        IF(SWAN_OutputHS)THEN
            CALL initNetCDFOutputFile(SwanHSDescript,reterr)
            CALL initNetCDFOutputFile(SwanHSMaxDescript,reterr)
        ENDIF
        IF(SWAN_OutputDIR)THEN
            CALL initNetCDFOutputFile(SwanDIRDescript,reterr)
            CALL initNetCDFOutputFile(SwanDIRMaxDescript,reterr)
        ENDIF
        IF(SWAN_OutputTM01)THEN
            CALL initNetCDFOutputFile(SwanTM01Descript,reterr)
            CALL initNetCDFOutputFile(SwanTM01MaxDescript,reterr)
        ENDIF
        IF(SWAN_OutputTPS)THEN
            CALL initNetCDFOutputFile(SwanTPSDescript,reterr)
            CALL initNetCDFOutputFile(SwanTPSMaxDescript,reterr)
        ENDIF
        IF(SWAN_OutputWIND)THEN
            CALL initNetCDFOutputFile(SwanWINDDescript,reterr)
            CALL initNetCDFOutputFile(SwanWINDMaxDescript,reterr)
        ENDIF
        IF(SWAN_OutputTM02)THEN
            CALL initNetCDFOutputFile(SwanTM02Descript,reterr)
            CALL initNetCDFOutputFile(SwanTM02MaxDescript,reterr)
        ENDIF
        IF(SWAN_OutputTMM10)THEN
            CALL initNetCDFOutputFile(SwanTMM10Descript,reterr)
            CALL initNetCDFOutputFile(SwanTMM10MaxDescript,reterr)
        ENDIF
    ENDIF
#endif

    IF (C3D.eqv. .TRUE. ) THEN
        CALL initNetCDFOutputFile(SigTStaDescript, reterr, &
        SalStaDescript, TempStaDescript)
        CALL initNetCDFOutputFile(RealQStaDescript, reterr, &
        ImaginaryQStaDescript, WZStaDescript)
        CALL initNetCDFOutputFile(Q20StaDescript, reterr, &
        LStaDescript, EVStaDescript)
        CALL initNetCDFOutputFile(SigTDescript, reterr, &
        SalDescript, TempDescript)
        CALL initNetCDFOutputFile(RealQDescript, reterr, &
        ImaginaryQDescript, WZDescript)
        CALL initNetCDFOutputFile(Q20Descript, reterr, &
        LDescript, EVDescript)
        CALL initNetCDFOutputFile(QSurfKp1Descript, reterr)
    ENDIF

!     Create NetCDF hotstart files if NetCDF was specified.
    IF ((NHSTAR == 3) .OR. (NHSTAR == 367) .OR. (NHSTAR == 368) .OR. &
    (NHSTAR == 5) .OR. (NHSTAR == 567) .OR. (NHSTAR == 568)) THEN
    ! must init both hotstart files ... there is not enough information
    ! available to determine if they will both be needed
        IF ((IHOT /= 367) .AND. (IHOT /= 567)) THEN
            Elev1Descript % file_name = 'fort.67'
            CALL initNetCDFHotstart(67, Elev1Descript, &
            Elev2Descript, VelDescript, CH1Descript, EtaDiscDescript, &
            NodeCodeDescript, NOFFDescript, reterr)
            IF (C3D.eqv. .TRUE. ) THEN
                CALL initNetCDFHotstart3D(67,NHSTAR)
            ENDIF
        ENDIF
        IF ((IHOT /= 368) .AND. (IHOT /= 568)) THEN
            Elev1Descript % file_name = 'fort.68'
            CALL initNetCDFHotstart(68, Elev1Descript, &
            Elev2Descript, VelDescript, CH1Descript, EtaDiscDescript, &
            NodeCodeDescript, NOFFDescript, reterr)
            IF (C3D.eqv. .TRUE. ) THEN
                CALL initNetCDFHotstart3D(68,NHSTAR)
            ENDIF
        ENDIF
        IF (IHARIND == 1) THEN
            IF ((IHOT /= 367) .AND. (IHOT /= 567)) THEN
                CALL initNetCDFHotstartHarmonic(67, &
                HarmElevFDLVDescript, HarmElevSLVDescript, &
                HarmUVelFDLVDescript, HarmVVelFDLVDescript, &
                HarmUVelSLVDescript, HarmVVelSLVDescript, reterr)
            ENDIF
            IF ((IHOT /= 368) .AND. (IHOT /= 568)) THEN
                CALL initNetCDFHotstartHarmonic(68, &
                HarmElevFDLVDescript, HarmElevSLVDescript, &
                HarmUVelFDLVDescript, HarmVVelFDLVDescript, &
                HarmUVelSLVDescript, HarmVVelSLVDescript, reterr)
            ENDIF
            IF (CHARMV.eqv. .TRUE. ) THEN
                IF ((IHOT /= 367) .AND. (IHOT /= 567)) THEN
                    CALL initNetCDFHotstartHarmonicMeansVariances( &
                    &                67, ELAVDescript, ELVADescript, &
                    XVELAVDescript, YVELAVDescript, XVELVADescript, &
                    YVELVADescript, reterr)
                ENDIF
                IF ((IHOT /= 368) .AND. (IHOT /= 568)) THEN
                    CALL initNetCDFHotstartHarmonicMeansVariances( &
                    &                68, ELAVDescript, ELVADescript, &
                    XVELAVDescript, YVELAVDescript, XVELVADescript, &
                    YVELVADescript, reterr)
                ENDIF
            ENDIF
        ENDIF
    ENDIF
! free up memory allocated for mesh and boundaries
    CALL freeNetCDFCoord()
#endif
!----------------------------------------------------------------------------
    END SUBROUTINE prepNetCDF
!----------------------------------------------------------------------------

!---------------------------------------------------------------------------
!                S U B R O U T I N E   M A K E   F I L E   N A M E
!---------------------------------------------------------------------------
!     jgf51.21.41 A little subroutine to make the file name from the
!     base name and the file extension. When the write_output module
!     is integrated into adcprep, this subroutine will be redundant.
!---------------------------------------------------------------------------
    subroutine makeFileName(descript)
    use global, only : OutputDataDescript_t
    implicit none
    type(OutputDataDescript_t), intent(inout) :: descript
    character(len=10) :: extString
          
    write(extString,'(i0)') descript % file_extension
    descript % file_name = trim(descript % file_basename) // &
    '.' // trim(extString)
!----------------------------------------------------------------------------
    end subroutine makeFileName
!----------------------------------------------------------------------------
         
          
