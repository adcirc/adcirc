!******************************************************************************
! PADCIRC VERSION 45.12 03/17/2006                                            *
!  last changes in this file VERSION 45.56                                    *
!                                                                             *
!******************************************************************************

    MODULE SIZES
    IMPLICIT NONE

!...SET NUMBER OF BYTES "SZ" IN REAL(SZ) DECLARATIONS
!...SET "NBYTE" FOR PROCESSING INPUT DATA RECORD LENGTH

#ifdef REAL4
    INTEGER, PARAMETER :: SZ = 4
    INTEGER, PARAMETER :: NBYTE=4
#endif
#ifdef REAL8
    INTEGER, PARAMETER :: SZ = 8
    INTEGER, PARAMETER :: NBYTE=8
#endif

!...SET MAX OF DIGITS OF PRECISION "NPREC" THE GRID CAN BE EXPECTED TO HAVE
!...NOTE: IF THE GRID WAS BUILT ON A 32 BIT COMPUTER, IT SHOULD BE
!   ACCURATE TO ABOUT 7 DIGITS.  THUS, IF THE NODAL SPACING REQUIRES MORE
!   THAN 5 DIGITS OF PRECISION, THE MODEL RESULTS MAY NOT BE TRUSTWORTHY.

    INTEGER, PARAMETER ::  NPREC=7

    INTEGER ::  MNE,MNP,MNEI,MNOPE,MNETA,MNBOU,MNVEL, &
    MNTIF,MNBFR,MNFFR,MNSTAE,MNSTAV,MNSTAC,MNSTAM,MNHARF
    INTEGER :: MNPROC         ! Number of compute processors
    INTEGER :: MNWPROC        ! Number of writer processors
    INTEGER :: MNWPROH        ! Number of Hwriter processors !st3 100711
    LOGICAL :: ISPLIT         ! .TRUE. if splitting output for writer !st3 100708:
    INTEGER :: MNALLPROC      ! Number of all processors (= MNPROC + MNWPROC)

!     Dimension of vertical FE mesh (To interface 2D & 3D)
!     INTEGER :: MNODES                                        !changed to MNFEN 08/16/2005
    INTEGER :: MNFEN

!     LOGICAL C2DDI,C3D,C3DDSS,C3DVS,CLUMP,CTIP,CHARMV         !moved to GLOBAL.F  04/09/2004

! For Definition of Working Directory

    INTEGER :: MYPROC
    LOGICAL      :: WRITE_LOCAL_FILES
    LOGICAL      :: WRITE_LOCAL_HOT_START_FILES
! jgf49.44: Make it an option to read local binary hotstart files
    LOGICAL      :: READ_LOCAL_HOT_START_FILES
! kmd51: Trying to make harmonic files write to local directory
    LOGICAL      :: WRITE_LOCAL_HARM_FILES

    CHARACTER(256) :: ROOTDIR
    CHARACTER(2048), TARGET :: INPUTDIR  ! directory containing input files
    CHARACTER(2048), TARGET :: GLOBALDIR ! directory for fulldomain files
    CHARACTER(2048), TARGET :: LOCALDIR  ! directory for subdomain files
    CHARACTER(2048), TARGET :: HARMDIR   ! directory where harmonic analysis file are written
    CHARACTER(2048) :: GBLINPUTDIR

! file formats
    integer, parameter :: numFormats = 7
    integer, parameter :: OFF = 0
    integer, parameter :: ASCII = 1
    integer, parameter :: BINARY = 2
    integer, parameter :: NETCDF3 = 3
    integer, parameter :: SPARSE_ASCII = 4
    integer, parameter :: NETCDF4 = 5
    integer, parameter :: XDMF = 7
    character(len=20) :: formatString(0:7) = (/ 'OFF         ', &
    'ASCII       ', 'BINARY      ', 'NETCDF3     ', 'SPARSE_ASCII', &
    'NETCDF4     ', 'UNDEF       ', 'XDMF        ' /)
               

! parameters related to input files

! mesh (fort.14) file
    INTEGER :: meshType
    LOGICAL :: meshFileNameSpecified ! .TRUE. if the user explicitly gave the fn
    CHARACTER(2048) :: meshFileName ! if different from the default fort.14
! name of full domain mesh file in parallel (only used by XDMF in parallel)
    CHARACTER(2048) :: meshFileName_g

! control (fort.15) file
    INTEGER :: controlType
    LOGICAL :: controlFileNameSpecified ! .TRUE. if the user explicitly gave the fn
    CHARACTER(2048) :: controlFileName ! if different from the default (fort.15)
! name of full domain mesh file in parallel (only used by XDMF in parallel)
    CHARACTER(2048) :: controlFileName_g

! nodal attributes (fort.13) file
    INTEGER :: naType
    LOGICAL :: naFileNameSpecified ! .TRUE. if the user explicitly gave the fn
    CHARACTER(2048) :: naFileName ! if different from the default fort.14
      

!---------------------end of data declarations--------------------------------C


    CONTAINS

!-----------------------------------------------------------------------
!     S U B R O U T I N E   M A K E _ D I R N A M E
!-----------------------------------------------------------------------
!     jgf: Determine which directories (paths) should be used
!     to locate input files and write output files; process command
!     line options related to file i/o.

!     jgf51.21.11: Added capability to allow user to specify name and
!     file format of the mesh file.
!-----------------------------------------------------------------------
    SUBROUTINE MAKE_DIRNAME()
! cm v50.85 added for windows builds
#ifdef CMPI
#ifdef WINDOWS
    USE IFPORT
#endif
#endif
    IMPLICIT NONE
    INTEGER :: LNAME, IARGC, ARGCOUNT, I, iprefix, res
    CHARACTER(2048) :: CMDLINEARG
    CHARACTER(2000) :: wholeCommandLine ! used to launch this executable, for log message
    CHARACTER(8)    :: PREFIX(2) = (/ '/PE0000 ', '/DOM0000' /)
    logical         :: fileFound
    logical(4) :: dir_result

    INPUTDIR  = ""
    GLOBALDIR = ""
    LOCALDIR  = ""
    ARGCOUNT  = IARGC()
    WRITE_LOCAL_HOT_START_FILES = .FALSE. 
    READ_LOCAL_HOT_START_FILES = .FALSE. 
    WRITE_LOCAL_FILES = .FALSE. 
    WRITE_LOCAL_HARM_FILES = .FALSE. 

! defaults
    meshType = ASCII
    meshFileNameSpecified = .FALSE. 
    meshFileName = 'fort.14'
    meshFileName_g = 'fort.14' ! used by XDMF writer in parallel to read mesh
          
    controlType = ASCII
    controlFileNameSpecified = .FALSE. 
    controlFileName = 'fort.15'
    controlFileName_g = 'fort.15' ! used by XDMF writer in parallel to read metadata
          
    naType = ASCII
    naFileNameSpecified = .FALSE. 
    naFileName = 'fort.13'

! Casey 090527: Debug.
    ROOTDIR = "."

! bell 110518: Add compiler flag for local hot start.
#ifdef LOCALHOT
    WRITE_LOCAL_HOT_START_FILES = .TRUE. 
#else
    WRITE_LOCAL_HOT_START_FILES = .FALSE. 
#endif

    if (ARGCOUNT > 0) then
    
    ! jgf51.21.41: Echo the entire command line that was used to
    ! launch this executable.
        i=0
        wholeCommandLine(:) = ' '
        do while (i < argcount)
            call getarg(i, cmdlinearg)
            wholeCommandLine = ' ' // trim(adjustl(wholeCommandLine)) &
            // ' ' // trim(adjustl(cmdlinearg))
            i=i+1
        end do
        if (myproc == 0) then
            write(*,'(a,a)') 'ECHO: Application launched with: ' &
            // trim(adjustl(wholeCommandLine))
        endif
    
    ! Now process the command line arguments.
        i = 0
        do while (i < ARGCOUNT)
            i = i + 1
            call getarg(i, CMDLINEARG)
            select case(cmdlinearg(1:2))
            case("-I","-i")
            i = i + 1
            call getarg(i,INPUTDIR)
            if (myProc == 0) then
                write(*,*) "INFO: Processing '",cmdlinearg(1:2)," ", &
                trim(INPUTDIR),"'."
            endif
            case("-O","-o")
            i = i + 1
            call getarg(i,GLOBALDIR)
            if (myProc == 0) then
                write(*,*) "INFO: Processing '",cmdlinearg(1:2)," ", &
                trim(GLOBALDIR),"'."
            endif
            case("-L","-l")
            WRITE_LOCAL_FILES = .TRUE. 
            WRITE_LOCAL_HOT_START_FILES = .TRUE. 
            if (myProc == 0) then
                write(*,*) "INFO: Processing '",cmdlinearg(1:2),"."
            endif
            case("-U","-u")
            WRITE_LOCAL_FILES = .TRUE. 
            if (myProc == 0) then
                write(*,*) "INFO: Processing '",cmdlinearg(1:2),"."
            endif
            case("-S","-s")
            WRITE_LOCAL_HOT_START_FILES = .TRUE. 
            if (myProc == 0) then
                write(*,*) "INFO: Processing '",cmdlinearg(1:2),"."
            endif
            case("-R","-r")
            READ_LOCAL_HOT_START_FILES = .TRUE. 
            if (myProc == 0) then
                write(*,*) "INFO: Processing '",cmdlinearg(1:2),"."
            endif
            case("-W","-H") !tcm v51.09 added test for -H (Hot start writers)
            i = i + 1
        ! this is the number of writer processors, and will
        ! be parsed in a later call to GET_NUMWRITERS
        !  ... could also be -Ws for round-robin writers
            case("-M","-m")
        ! jgf51.21.11: Added options related to the mesh file
            select case(cmdlinearg(1:4))
            case("-MFT","-mft") ! mesh file type
            i = i + 1
            call getarg(i, CMDLINEARG)
            select case(trim(cmdlinearg))
            case("ASCII","Ascii","ascii")
            meshType = ASCII ! this is the default anyway
            case("XDMF","Xdmf","xdmf")
            meshType = XDMF
            case default
            if (myProc == 0) then
                write(*,*) "WARNING: The mesh file type '", &
                trim(cmdlinearg), &
                "' is not valid. Valid types are ascii and xdmf."
            endif
            end select
            case("-MFN","-mfn") ! mesh file name
            meshFileNameSpecified = .TRUE. 
            i = i + 1
            call getarg(i,meshFileName)
            case default
        ! jgf51.44: -m by itself triggers subdomain harmonic
        ! analysis output files.
            if ( len(trim(cmdlinearg)) == 2 ) then
                WRITE_LOCAL_HARM_FILES = .TRUE. 
                if (myProc == 0) then
                    write(*,'(a,a,a)') "INFO: Processing " &
                    // cmdlinearg(1:2) // "."
                endif
            endif
            if (myProc == 0) then
                write(*,*) "WARNING: The command line option '", &
                cmdlinearg(1:4),"' is not valid and will be ignored."
            endif
            end select
            case("-N","-n")
        ! jgf51.21.11: Added options related to the nodal attributes file
            select case(cmdlinearg(1:4))
            case("-NFT","-nft") ! nodal attributes file type
            i = i + 1
            call getarg(i, CMDLINEARG)
            select case(trim(cmdlinearg))
            case("ASCII","Ascii","ascii")
            naType = ASCII ! this is the default anyway
            case("XDMF","Xdmf","xdmf")
            naType = XDMF
            case default
            if (myProc == 0) then
                write(*,*) &
                "WARNING: The nodal attributes file type '", &
                trim(cmdlinearg), &
                "' is not valid. Valid types are ascii and xdmf."
            endif
            end select
            case("-NFN","-nfn") ! nodal attributes file name
            naFileNameSpecified = .TRUE. 
            i = i + 1
            call getarg(i,naFileName)
            case default
            if (myProc == 0) then
                write(*,*) "WARNING: The command line option '", &
                cmdlinearg(1:4),"' is not valid and will be ignored."
            endif
            end select
            case default
            if (myProc == 0) then
                write(*,*) "WARNING: The command line option '", &
                cmdlinearg(1:2),"' is not valid and will be ignored."
            endif
            end select
        end do
    end if

!.....Default root working directory

    if (len_trim(INPUTDIR) /= 0) then
        ROOTDIR = INPUTDIR
    endif

    GBLINPUTDIR = ROOTDIR
#ifdef CMPI
    iprefix = 0
    if (myProc == 0) then
        write(*,*) "INFO: Searching for ADCIRC subdomain directories:"
    endif
    do i = 1, 2
        if (myProc == 0) then
            write(*,*) "INFO: Looking for '",trim(ROOTDIR), &
            trim(PREFIX(i)),"/fort.14' ..."
        endif
        INQUIRE(file=TRIM(ROOTDIR)//TRIM(PREFIX(i))//'/'//'fort.14', &
        exist=fileFound)
        if (fileFound.eqv. .TRUE. ) then
            if (myProc == 0) then
                write(*,*) "INFO: File '",trim(ROOTDIR), &
                trim(PREFIX(i)),"/fort.14' was found!"
                write(*,*) "INFO: The search for the subdomain ", &
                "directory was completed successfully."
            endif
            iprefix = i
            exit
        else
            write(*,*) "ERROR: Processor ",myProc,": File '", &
            trim(ROOTDIR),trim(PREFIX(i)), &
            "/fort.14' was not found."

            print *, "ERROR: ADCIRC stopping."
            call msg_abort()
        end if
    end do

    WRITE(INPUTDIR,'(2A)') TRIM(ROOTDIR),PREFIX(iprefix)
    LNAME = LEN_TRIM(INPUTDIR)
    WRITE(INPUTDIR(LNAME-3:LNAME),'(I4.4)') MYPROC
#else
    WRITE(INPUTDIR,'(A)') TRIM(ROOTDIR)
#endif

    if (len_trim(GLOBALDIR) /= 0) then
        ROOTDIR = GLOBALDIR
    endif

    WRITE(GLOBALDIR,'(A)') TRIM(ROOTDIR)


#ifdef CMPI
    WRITE(LOCALDIR,'(2A)') TRIM(ROOTDIR),TRIM(PREFIX(iprefix))
    LNAME = LEN_TRIM(LOCALDIR)
    WRITE(LOCALDIR(LNAME-3:LNAME),'(I4.4)') MYPROC

!... tcm v49.67 -- additions for Windows based PC compilations
#ifdef WINDOWS
! all system('mkdir '//trim(LOCALDIR))  !tcm v50.85 this no longer works for INTEL MPI
    dir_result = makedirqq(trim(localdir)) !tcm v50.85 added
#else
    call MAKEDIR(trim(LOCALDIR))
#endif

#else
    WRITE(LOCALDIR,'(A)') TRIM(ROOTDIR)
#endif

! jgf51.21.27: Form mesh file name and nodal attributes file name.
    if (meshFileNameSpecified.eqv. .FALSE. ) then
        meshFileName = trim(LOCALDIR)//'/fort.14'
    endif
    if (naFileNameSpecified.eqv. .FALSE. ) then
        naFileName = trim(LOCALDIR)//'/fort.13'
    endif

!     jgf49.17.01 Summarize and log results.
    if (myProc == 0) then
        write(*,*) "INFO: The ROOTDIR is  '",trim(ROOTDIR),"'."
        write(*,*) "INFO: The INPUTDIR is '",trim(INPUTDIR),"'."
        write(*,*) "INFO: The GBLINPUTDIR is '",trim(GBLINPUTDIR),"'."
        write(*,*) "INFO: The GLOBALDIR is '",trim(GLOBALDIR),"'."
        write(*,*) "INFO: The LOCALDIR is '",trim(LOCALDIR),"'."
        if ((WRITE_LOCAL_FILES.eqv. .TRUE. ) .AND. (MNPROC > 1)) then
            write(*,*) "INFO: ADCIRC will write subdomain output files."
        endif
        if ((WRITE_LOCAL_HOT_START_FILES.eqv. .TRUE. ) &
         .AND. (MNPROC > 1)) then
            write(*,*) &
            "INFO: ADCIRC will write subdomain hotstart files."
        endif
        if ((WRITE_LOCAL_HARM_FILES.eqv. .TRUE. ) .AND. (MNPROC > 1)) then
            write(*,*) 'INFO: ADCIRC will write ' &
            //  'subdomain harmonic analysis files.'
        endif
    endif
    RETURN
!-----------------------------------------------------------------------
    END SUBROUTINE MAKE_DIRNAME
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
!     S U B R O U T I N E   G E T _ N U M  W R I T E R S
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    SUBROUTINE GET_NUMWRITERS
    IMPLICIT NONE
    INTEGER :: IARGC, ARGCOUNT, I
    CHARACTER(2048) :: CMDLINEARG
    ARGCOUNT  = IARGC()
    MNWPROC = 0
    isplit = .FALSE. !st3 100708:split file

    if (ARGCOUNT > 0) then
        i = 0
        do while (i < ARGCOUNT)
            i = i + 1
            call getarg(i, CMDLINEARG)
            if (cmdlinearg(1:2) == "-W") then
                if( cmdlinearg(1:3) == "-Ws" ) then  !st3 100708: split file
                    isplit = .TRUE.                   !st3 100708: split file
                endif                                !st3 100708: split file
                i = i + 1
                call getarg(i,cmdlinearg)
                read(cmdlinearg,*) MNWPROC
            endif
        end do
    end if

    MNWPROH = 0  !st3 100711 for hsfile
    if (ARGCOUNT > 0) then
        i = 0
        do while (i < ARGCOUNT)
            i = i + 1
            call getarg(i, CMDLINEARG)
            if (cmdlinearg(1:2) == "-H") then
                i = i + 1
                call getarg(i,cmdlinearg)
                read(cmdlinearg,*) MNWPROH
            endif
        end do
    end if

!.....Default root working directory

    RETURN
!-----------------------------------------------------------------------
    END SUBROUTINE GET_NUMWRITERS
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    END MODULE SIZES
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
