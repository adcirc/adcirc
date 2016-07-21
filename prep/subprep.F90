
!                ADCIRC - SUBDOMAIN MODELING PREPROCESSING MODULE

!    ========================================================================
!    |                                                                      |
!    |   This file contains the subroutines required by Subdomain Modeling, |
!    |   an approach to reduce the total runtime of a series of hurricane   |
!    |   storm surge simulations on a smaller grid, a subdomain grid,       |
!    |   extracted from the original grid.                                  |
!    |                                                                      |
!    |   Written by Alper Altuntas, aaltunt@ncsu.edu                        |
!    |   North Carolina State University,                                   |
!    |   2013                                                               |
!    |                                                                      |
!    ========================================================================

    module subprep
    use sizes, only : sz
    use global, only : DEBUG, ECHO, INFO, WARNING, ERROR, &
    openFileForRead, scratchMessage, logMessage, allMessage, &
    setMessageSource, unsetMessageSource

!   NCSU Subdomain Modeling variables:
    logical :: subdomainOn
    integer :: enforceBN
    integer :: psbtiminc, pncbn, pnobn, pnibn
    real(sz) :: pebn,pubn,pvbn
    integer,allocatable :: pcbn(:),pncbnp(:)
    integer,allocatable :: pobn(:),pnobnp(:)
    integer,allocatable :: pibn(:),pnibnp(:)
    integer :: pwdbn
    logical :: found_sm_nml

    contains


    SUBROUTINE readFort015prep()

    implicit none
    integer :: dummy
    integer :: ioerror  ! zero if the fort.15 file was opened
    namelist /subdomainModeling/ subdomainOn

    FOUND_SM_NML = .FALSE. 

    call setMessageSource("readFort015prep")
#if defined(SUBPREP_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Enter.")
#endif

    ioerror = 0

! Look for the subdomainModeling namelist in the fort.15 file
! and if it is found, read the value of subdomainOn.
    call openFileForRead(15,'fort.15',ioerror)
    if (ioerror /= 0) then
        call allMessage(ERROR,'Failed to open fort.15 file.')
        stop
    endif

! jgf51.42: Add a namelist for the user to control subdomain
! modeling.
    READ(UNIT=15,NML=subdomainModeling,IOSTAT=ioerror)
    IF (ioerror > 0) THEN
        call logMessage(INFO, &
        'The subdomainModeling namelist was not found.')
    else if(ioerror == 0) then
        FOUND_SM_NML = .TRUE. 
    endif
    if ((ioerror > 0) .OR. (subdomainOn.eqv. .FALSE. )) then
        call logMessage(INFO, &
        'The subdomainModeling capability is not active.')
    endif
    write(scratchMessage,*) "subdomainOn=",subdomainOn
    call logMessage(ECHO,trim(scratchMessage))
    close(15)
                    
    if (subdomainOn) then
        call openFileForRead(1015,'fort.015',ioerror)
        if (ioerror == 0) then
            print *, "Subdomain Active"
            open(1015, file='fort.015')
            read(1015,*) dummy
            read(1015,*) dummy
            read(1015,*) enforceBN
            close(1015)
        else
            call allMessage(ERROR,'The fort.015 file was not found.')
            stop
        endif
    endif

#if defined(SUBPREP_TRACE) || defined(ALL_TRACE)
    call allMessage(DEBUG,"Return.")
#endif
    call unsetMessageSource()

    END SUBROUTINE readFort015prep




    SUBROUTINE openFort019prep()

    Use PRE_GLOBAL, only : nproc, itotproc, imap_nod_gl2

    implicit none
    integer :: i,j,n,node,iproc
    integer :: gn,p,totalProc,itemp,procTemp
     
    open(1019, file='fort.019')
    read(1019,*)
    read(1019,*) psbtiminc, pncbn
    allocate( pcbn(pncbn) )
    allocate( pncbnp(nproc) )
    do i=1,pncbn
        read(1019,*) pcbn(i)
    enddo

    do i=1,nproc
        pncbnp(i) = 0
    enddo

    do i=1,pncbn
        gn = pcbn(i)
        totalProc = itotproc(gn)
        do p=1,totalProc
            itemp = (p-1)*2+1
            procTemp = IMAP_NOD_GL2(itemp,gn)
            pncbnp(procTemp) = pncbnp(procTemp) + 1
        enddo
    enddo

    END SUBROUTINE openFort019prep
         


     

    SUBROUTINE subPrep019()

    Use PRE_GLOBAL, only : nproc, itotproc, imap_nod_gl2

    implicit none
    integer :: iproc,i,n,tstep,gnode,lnode
    integer :: totalProc,p,itemp,proctemp,gn,lntemp
    integer :: sdu(nproc)
    logical :: success

    call openFort019prep()
    call subOpenPrep(19,'subdomain boundary conditions ', &
    &         1,nproc, SDU, Success)


    do iproc=1,nproc
        write(sdu(iproc),*) 'subdomain boundary conditions '
        write(sdu(iproc),*) psbtiminc,pncbnp(iproc)
    enddo

    do i=1,pncbn
        gnode = pcbn(i)
        totalProc = itotproc(gnode)
        do p=1,totalProc
            itemp = (p-1)*2+1
            procTemp = IMAP_NOD_GL2(itemp,gnode)
            lnTemp = IMAP_NOD_GL2(itemp+1,gnode)
            write(sdu(procTemp),*) lnTemp
        enddo
    enddo

    1000 CONTINUE

    read(1019,*,end=1905) tstep
    do iproc=1,nproc
        write(sdu(iproc),*) tstep
    enddo
       

    do i=1,pncbn
        read(1019,*) gnode,pebn,pubn
        read(1019,*) pvbn,pwdbn
        totalProc = itotproc(gnode)
        do p=1,totalProc
            itemp = (p-1)*2+1
            procTemp = IMAP_NOD_GL2(itemp,gnode)
            lnTemp = IMAP_NOD_GL2(itemp+1,gnode)
            write(sdu(procTemp),*) lnTemp,pebn,pubn
            write(sdu(procTemp),*) pvbn,pwdbn
        enddo
    enddo
                
    go to 1000

    1905 close(1019)
    do iproc=1,nproc
        close(sdu(iproc))
    enddo

    END SUBROUTINE subPrep019





    SUBROUTINE openFort020prep()

    Use PRE_GLOBAL, only : nproc, itotproc, imap_nod_gl2

    implicit none
    integer :: i,j,n,node,iproc
    integer :: gn,p, totalProc, itemp,procTemp

    open(1020, file='fort.020')
    read(1020,*)
    read(1020,*) psbtiminc, pnobn
    allocate( pobn(pnobn) )
    allocate( pnobnp(nproc) )
    do i=1,pnobn
        read(1020,*) pobn(i)
    enddo


    do i=1,nproc
        pnobnp(i) = 0
    enddo


    do i=1,pnobn
        gn = pobn(i)
        totalProc = itotproc(gn)
        do p=1,totalProc
            itemp = (p-1)*2+1
            procTemp = IMAP_NOD_GL2(itemp,gn)
        ! nTemp = IMAP_NOD_GL2(itemp+1,gn)
            pnobnp(procTemp) = pnobnp(procTemp) + 1
        enddo
    enddo

                 
    do i=1,nproc
        print *, i, pnobnp(i)
    enddo
         
    END SUBROUTINE openFort020prep




    SUBROUTINE subPrep020()

    use pre_global, only : nproc, itotproc, imap_nod_gl2
    implicit none
    integer :: iproc,i,n,tstep,gnode,lnode
    integer :: totalProc,p,itemp,proctemp,gn,lntemp
    integer :: sdu(nproc)
    logical :: success

    call openFort020prep()
    call subOpenPrep(20,'subdomain boundary conditions ', &
    &         1,nproc, SDU, Success)


    do iproc=1,nproc
        write(sdu(iproc),*) 'subdomain boundary conditions '
        write(sdu(iproc),*) psbtiminc,pnobnp(iproc)
    enddo
     
    do i=1,pnobn
        gnode = pobn(i)
        totalProc = itotproc(gnode)
        do p=1,totalProc
            itemp = (p-1)*2+1
            procTemp = IMAP_NOD_GL2(itemp,gnode)
            lnTemp = IMAP_NOD_GL2(itemp+1,gnode)
            write(sdu(procTemp),*) lnTemp
        enddo
    enddo

    1000 CONTINUE

    read(1020,*,end=1905) tstep
    do iproc=1,nproc
        write(sdu(iproc),*) tstep
    enddo

!        do i=1,pnobn
!            read(1020,*) gnode,pebn,pubn
!            read(1020,*) pvbn,pwdbn
!            iproc = IMAP_NOD_GL(1,gnode)
!            lnode = IMAP_NOD_GL(2,gnode)
!            write(sdu(iproc),*) lnode,pebn,pubn
!            write(sdu(iproc),*) pvbn,pwdbn
!        enddo


    do i=1,pnobn
        read(1020,*) gnode,pebn,pubn
        read(1020,*) pvbn,pwdbn
        totalProc = itotproc(gnode)
        do p=1,totalProc
            itemp = (p-1)*2+1
            procTemp = IMAP_NOD_GL2(itemp,gnode)
            lnTemp = IMAP_NOD_GL2(itemp+1,gnode)
            write(sdu(procTemp),*) lnTemp,pebn,pubn
            write(sdu(procTemp),*) pvbn,pwdbn
        enddo
    enddo


    go to 1000

    1905 close(1020)
    do iproc=1,nproc
        close(sdu(iproc))
    enddo

    END SUBROUTINE subPrep020




    SUBROUTINE openFort021prep()

    Use PRE_GLOBAL, only : nproc, itotproc, imap_nod_gl2

    implicit none
    integer :: i,j,n,node,iproc
    integer :: gn,p, totalProc, itemp,procTemp


    open(1021, file='fort.021')
    read(1021,*)
    read(1021,*) psbtiminc, pnibn
    allocate( pibn(pnibn) )
    allocate( pnibnp(nproc) )
    do i=1,pnibn
        read(1021,*) pibn(i)
    enddo




    do i=1,nproc
        pnibnp(i) = 0
    enddo


    do i=1,pnibn
        gn = pibn(i)
        totalProc = itotproc(gn)
        do p=1,totalProc
            itemp = (p-1)*2+1
            procTemp = IMAP_NOD_GL2(itemp,gn)
        ! nTemp = IMAP_NOD_GL2(itemp+1,gn)
            pnibnp(procTemp) = pnibnp(procTemp) + 1
        enddo
    enddo

    END SUBROUTINE openFort021prep




    SUBROUTINE subPrep021()

    use PRE_GLOBAL, only : nproc, itotproc, imap_nod_gl2

    implicit none
    integer :: iproc,i,n,tstep,gnode,lnode
    integer :: totalProc,p,itemp,proctemp,gn,lntemp
    integer :: sdu(nproc)
    logical :: success

    call openFort021prep()
    call subOpenPrep(21,'subdomain boundary conditions ', &
    &         1,nproc, SDU, Success)


    do iproc=1,nproc
        write(sdu(iproc),*) 'subdomain boundary conditions '
        write(sdu(iproc),*) psbtiminc,pnibnp(iproc)
    enddo

    do i=1,pnibn
        gnode = pibn(i)
        totalProc = itotproc(gnode)
        do p=1,totalProc
            itemp = (p-1)*2+1
            procTemp = IMAP_NOD_GL2(itemp,gnode)
            lnTemp = IMAP_NOD_GL2(itemp+1,gnode)
            write(sdu(procTemp),*) lnTemp
        enddo
    enddo


    1000 CONTINUE

    read(1021,*,end=1905) tstep
    do iproc=1,nproc
        write(sdu(iproc),*) tstep
    enddo

!        do i=1,pnibn
!            read(1021,*) gnode,pebn
!            iproc = IMAP_NOD_GL(1,gnode)
!            lnode = IMAP_NOD_GL(2,gnode)
!            write(sdu(iproc),*) lnode,pebn
!        enddo


    do i=1,pnibn
        read(1021,*) gnode,pebn
        totalProc = itotproc(gnode)
        do p=1,totalProc
            itemp = (p-1)*2+1
            procTemp = IMAP_NOD_GL2(itemp,gnode)
            lnTemp = IMAP_NOD_GL2(itemp+1,gnode)
            write(sdu(procTemp),*) lnTemp,pebn
        enddo
    enddo

    go to 1000

    1905 close(1021)
    do iproc=1,nproc
        close(sdu(iproc))
    enddo

    END SUBROUTINE subPrep021








!     This subroutine is copied from prep.F (OpenPrepFiles) and is modified to
!     open subdomain modeling b.c. files

    SUBROUTINE subOpenPrep(UnitNumber, Description, &
    startProc, endProc, SDU, Success)
!     NCSU SUBDOMAIN: modified to handle "fort.019"
!---------------------------------------------------------------------------
    USE PRE_GLOBAL, only : nproc
    IMPLICIT NONE
    INTEGER :: UnitNumber     ! i/o unit number to open
    CHARACTER(len=30), intent(in) :: Description ! description of file
    INTEGER, intent(in) :: startProc        ! subdomains to start with
    INTEGER, intent(in) :: endProc          ! subdomain to end on
    INTEGER, intent(out), dimension(nproc) :: SDU ! Subdomain unit numbers
    LOGICAL, intent(out):: Success     ! .TRUE. if files opened w/o errors
    LOGICAL :: Found               ! .TRUE. if the full domain file exists
    CHARACTER(len=80) FileName   ! name of full domain file
    CHARACTER(len=8) DefaultName! default name of full domain file  !NCSU SUBDOMAIN
    INTEGER :: ErrorIO             ! zero if file opened successfully
    INTEGER :: iproc               ! subdomain index
    CHARACTER(len=15) sdFileName     ! subdomain file name  !NCSU SUBDOMAIN

    Found = .FALSE. 
    Success = .FALSE. 
    ErrorIO = 1

! NCSU Subdomain:
    if (UnitNumber == 19) then
        DefaultName= 'fort.019'
        FileName = 'fort.019'
        UnitNumber = 1019
    else if (UnitNumber == 20) then
        DefaultName= 'fort.020'
        FileName = 'fort.020'
        UnitNumber = 1020
    else if (UnitNumber == 21) then
        DefaultName= 'fort.021'
        FileName = 'fort.021'
        UnitNumber = 1021
    endif


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
    ENDIF

    If ( .NOT. Success) RETURN ! failed to open full domain file

    DO iproc = startProc, endProc
        sdu(iproc) = 105 + (iproc-1)
        sdFileName(1:7) = 'PE0000/'
        sdFileName(8:15) = DefaultName   !NCSU SUBDOMAIN
#ifdef ADCSWAN
        sdFileName = 'PE0000/'//FileName
#endif
        CALL IWRITE(sdFileName, 3, 6, iproc-1)
        OPEN (UNIT=SDU(iproc), FILE=sdFileName, IOSTAT=ErrorIO)
        Success = .TRUE. 
        IF ( ErrorIO > 0 ) THEN
            WRITE(*,*) "ERROR: Subdomain file cannot be opened."
            Success = .FALSE. 
            RETURN ! failed to open at least one subdomain file
        ENDIF
    ENDDO

    2 FORMAT(I2)
    30 FORMAT(A30)
    1010 FORMAT('File ',A8,/,' WAS NOT FOUND! Try again or type "skip"',/)
    1011 FORMAT('File ',A8,/,' WAS FOUND!  Opening & Processing file.',/)
    RETURN
!---------------------------------------------------------------------------
    END SUBROUTINE subOpenPrep
!---------------------------------------------------------------------------

    end module subprep








































