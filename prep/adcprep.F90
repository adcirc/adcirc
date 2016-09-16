!----------------------------------------------------------------------------

!                        ADCPREP Version 2.2 (  09/20/2006 )

!----------------------------------------------------------------------------

!  Program Development History
!  ---------------------------
!   Written for ADCIRC_v24.05    ( S. Chippada 1996 )
!   Updated for ADCIRC_v24.05    ( M. Martinez 1997 )
!   Hilbert Space Filling Curve  ( C. Edwards 1997 )
!   Updated for ADCIRC_v34.04    ( V. Parr 1998 )
!   METIS 4.0 Interface          ( V. Parr 1999 )
!   Added UNIT 12 for polders    ( V. Parr 1999 )
!   hotstarts for Windstress     ( V. Parr 1999 )
!   Updated for ADCIRC_v34.08    ( V. Parr 1999 )
!   Modified for PADCIRC90       ( V. Parr & J. Romo 1999 )
!   Parallel version for 2 PEs   ( V. Parr & J. Romo 2000 )
!   Add hotprep1 command         ( V. Parr & J. Romo 2001 )
!   Add hotprep2 command         ( R. Luettich, 10/2001 )
!   Expanded to include 3D       ( T. Campbell 2002)
!   Bug fix                      ( J. Fleming/R. Luettich 2/2002)
!   Added wave stress capability ( M. Brown 3/2003)
!   NWS=3,6 fixes                ( M. Brown 3/2003)
!   NWS=+-4 revised to not convert to NWS = 5  ( R. Luettich 3/2003)
!   Wave stress routines revise  ( R. Luettich 3/2003)
!   Time-varying Flux files      ( M. Brown 4/2004 )
!   Added default filename capability  ( M. Brown 4/2004 )
!   Islands kept as islands      ( R. Luettich 5/2004)
!   Vic Parr's bug fix for NWS=+-4, Wave stress  ( R. Luettich 5/2004)
!   Removed NWS +/- conversion   ( R. Luettich/jgf 11/2005 )
!   Add PREP14_15 option         ( J. Fleming 11/2005 )
!   3D Update: Recording Stations( J. Fleming 02/2006 )
!   Wrote PREP11,rewrote PREP20  ( J. Fleming 03/2006 )

!   Add partmesh command         ( V. Parr 09/2006 )
!   Add routine relocalize       ( V. Parr 09/2006 )
!   Add prep13 ONLY command      ( V. Parr 09/2006 )
!   Add prep15 ONLY command      ( V. Parr 10/2006 )
!   Add prephot ONLY command     ( V. Parr 10/2006 )
!   Added NWS=19 asymmetric holland v2.0 ( R. Weaver 07/2009)
!   Removed NWS=3 and NWS=6      (T.C. Massey 11/2009)
!   Add NWS=20 generalized asymmetric hurricane vortex (Jie 5/2013)

!----------------------------------------------------------------------------

!  ADCPREP performs 7 operations:
!    (1) partmesh  =  partition the mesh with metis
!    (2) prepall   =  Localize grid, input, and hot start files using default
!                     file names
!    (3) prepspec  =  same as prepall but allows custom file names and
!                     skipping over files.
!    (4) prep15    =  Localizes ONLY a RunInfo file ( fort.15)
!                     after a prepall has been run
!    (5) prep13    =  Localizes ONLY a nodal attributes file ( fort.13)
!                     after a prepall has been run
!    (6) hotlocal  =  Localize  ONLY a hotstart file  (fort.67 or fort.68)
!    (7) hotglobal =  Globalize ONLY a hotstart files (fort.67 or fort.68)

!  Pre-processor
!  -------------
!  The Pre-processor "prep" performs data decomposition of the ADCIRC grid
!  and localizes global input files to subdomains of the decomposition for
!  parallel runs of ADCIRC.

!  Subdirectories "PExxxx" of the working directory ( where "xxxx" is a
!  4-digit integer representing the MPI_rank of a Processing Element )
!  are created and "localized" input files are placed in the appropriate
!  subdirectory with the same filename as the global file in the working
!  directory.  A unit 14 and 15 file, called respectively, setup.14 and
!  setup.15  are written to the working directory to allow setup to
!  prepare the source code for a parallel run.

!  The global input files currently localized by ADCPREP are:

!   fort.13    Global Nodal Attributes File
!   fort.14    Global Mesh File
!   fort.15    Global Input File
!   fort.19    Global Time-Variant Elevation Boundary Conditions
!   fort.20    Global Time-Variant Normal Flow Boundary Conditions
!   fort.22    Global Wind Stress and Atmospheric Pressure
!   fort.23    Global Wave Stress Forcing
!   fort.67    Global Hot Start file
!   fort.68    Global Hot Start file
!   fort.150   METIS input information
!   fort.151   METIS output information
!   fort.141   Global Time-Variant Bathymetry

!  When a parallel ADCIRC job is run, each PExxxx takes its input from its
!  subdirectory PExxxx and writes local output files to the same directory.

!----------------------------------------------------------------------------

    PROGRAM ADCPREP
    USE GLOBAL, ONLY : screenUnit, myProc, initLogging, Rearth, deg2rad, &
    rad2deg
    USE SIZES, ONLY : meshFileName, meshFileNameSpecified, meshType, &
    ASCII, XDMF, naFileNameSpecified, naFileName, naType, controlFileName, &
    controlType, controlFileNameSpecified, globaldir
    USE SUBPREP, ONLY : readFort015prep
    USE PRE_GLOBAL
    USE VERSION
    use memory_usage
    IMPLICIT NONE

    CHARACTER JOB*80 ! user input for menu selection
    INTEGER :: ARGCOUNT ! number of command line arguments
    INTEGER :: IARGC    ! function to return command line arguments
    INTEGER :: I        ! loop counter for command line arguments
    CHARACTER(2048) :: CMDLINEARG ! content of cmd line arg

! set the names of the log levels for use in log messages
    call initLogging()
! set the present working directory as the directory containing
! the fulldomain input files as the
    globaldir = '.'

! jgf51.42: NCSU Subdomain Modeling
    call readFort015prep()

!  Define File Format Version

    FileFmtVersion = VERSION_NUMBER(FileFmtMajor, &
    FileFmtMinor, &
    FileFmtRev)


!     Initialize flags.

    USE_DEFAULT = .FALSE. 
    APERIODIC_FLOW_BC = .FALSE. 
    PARTITION  = .FALSE. 
    PREP_ALL   = .FALSE. 
    PREP_13    = .FALSE. 
    PREP_15    = .FALSE. 
    HOT_LOCAL   = .FALSE. 
    HOT_GLOBAL  = .FALSE. 
    PREP_20 = .FALSE. 
! REP_88 = .FALSE.
    MNPROC = 0

! jgf51.21.27: Add defaults for file types now that we support the
! loading of different file types.
    naType = ASCII
    naFileNameSpecified = .FALSE. 
    naFileName = 'fort.13'
    meshType = ASCII
    meshFileNameSpecified = .FALSE. 
    meshFileName = 'fort.14'
    controlType = ASCII
    controlFileNameSpecified = .FALSE. 
    controlFileName = 'fort.15'

!--Init Memory Usage

    call memory_init()

!     jgf50.11: Parse command line arguments. This is needed initially
!     so that fort.20 and fort.88 can be prepped individually (without
!     having to do a prep all and without having to add menu items).
    ARGCOUNT = IARGC() ! count up command line options
    IF (ARGCOUNT > 0) THEN
        I=0
        DO WHILE (I < ARGCOUNT)
            I = I + 1
            CALL GETARG(I, CMDLINEARG)
            WRITE(*,*) "INFO: Processing ",TRIM(CMDLINEARG)
            SELECT CASE(TRIM(CMDLINEARG))
            CASE("--partmesh")
            PARTITION = .TRUE. 
            CASE("--prepall")
            PREP_ALL = .TRUE. 
            CASE("--prep13")
            PREP_13 = .TRUE. 
            CASE("--prep15")
            PREP_15 = .TRUE. 
            CASE("--prep20")
            PREP_20 = .TRUE. 
        !            CASE("--prep88")
        !               PREP_88 = .TRUE.
            CASE("--np")
            I = I + 1
            CALL GETARG(I,CMDLINEARG)
            READ(CMDLINEARG,*) MNPROC
        ! jgf51.21.27: Added options related to the mesh file
            case("-MFT","-mft","--MFT","--mft") ! mesh file type
            i = i + 1
            call getarg(i, CMDLINEARG)
            select case(trim(cmdlinearg))
            case("ASCII","Ascii","ascii")
            meshType = ASCII ! this is the default anyway
            case("XDMF","Xdmf","xdmf")
            meshType = XDMF
            case default
            write(*,*) "WARNING: The mesh file type '", &
            trim(cmdlinearg), &
            "' is not valid. Valid types are ascii and xdmf."
            end select
            case("-MFN","-mfn","--MFN","--mfn") ! mesh file name
            meshFileNameSpecified = .TRUE. 
            i = i + 1
            call getarg(i,meshFileName)
        ! jgf51.21.27: Added options related to the nodal attributes file
            case("-NFT","-nft","--NFT","--nft") ! nodal attributes file type
            i = i + 1
            call getarg(i, CMDLINEARG)
            select case(trim(cmdlinearg))
            case("ASCII","Ascii","ascii")
            naType = ASCII ! this is the default anyway
            case("XDMF","Xdmf","xdmf")
            naType = XDMF
            case default
            write(*,*) &
            "WARNING: The nodal attributes file type '", &
            trim(cmdlinearg), &
            "' is not valid. Valid types are ascii and xdmf."
            end select
            case("--NFN","--nfn") ! nodal attributes file name
            naFileNameSpecified = .TRUE. 
            i = i + 1
            call getarg(i,naFileName)
            case default
            write(*,*) "WARNING: The command line option '", &
            trim(cmdlinearg),"' is not valid and will be ignored."
            end select
        END DO
    ! Conk out if MNPROC was not specified
        IF (MNPROC == 0) THEN
            WRITE(*,*) "ERROR: MNPROC was not specified."
            STOP
        ELSE
            NPROC=MNPROC
        ENDIF
    ! The command-line and menu-driven options are mutually
    ! exclusive.
        USE_DEFAULT = .TRUE. 
        CALL PREPINPUT()
        call memory_status()
        STOP
    ENDIF


!-- Say Hello Gracie

    print *," *****************************************"
    print *," ADCPREP Fortran90 Version 2.3  10/18/2006"
    print *," Serial version of ADCIRC Pre-processor   "
    print *," *****************************************"
    print *, " "

!-- Prompt for user input

    print *, "Input number of processors for parallel ADCIRC run:"
    READ(*,*) MNPROC

!-- Copy MNPROC to NPROC since they are to be the same

    NPROC=MNPROC

    print *, '-------------------------------------------------------'
    print *, 'Preparing input files for subdomains.'
    print *, 'Select number or action:'
    print *, ' 1. partmesh'
    print *, ' - partition mesh using metis ( perform this first)'
    print *, ''
    print *, ' 2. prepall'
    print *, ' - Full pre-process using default names (i.e., fort.14)'
    print *, ''
    print *, ' 3. prepspec'
    print *, ' - Full pre-process except user may specify the names'
    print *, '   of input files. This option also allows the user'
    print *, '   to skip the preprocessing of certain files.'
    print *, ''
    print *, ' 4. prep15'
    print *, ' - Localizes RunInfo (fort.15) file ONLY'
    print *, '   Assumes a prepall has been run previously'
    print *, ''
    print *, ' 5. prep13'
    print *, ' - Localizes NodalAttributes (fort.13) file ONLY'
    print *, '   Assumes a prepall has been run previously'
    print *, ''
    print *, ' 6. hotLocalize'
    print *, ' - Localizes global hotstart file ONLY'
    print *, ''
    print *, ' 7. hotGlobalize'
    print *, ' - Globalizes local hotstart files ONLY'
    print *, '-------------------------------------------------------'

    9999 READ(*,*)  JOB

    SELECT CASE(JOB)

    CASE("1","partmesh","PARTMESH")
    PARTITION   = .TRUE. 

    CASE("2","prepall","PREPALL")
    USE_DEFAULT = .TRUE. 
    PREP_ALL = .TRUE. 
    PRINT*,''
    PRINT*,'Using default filenames.'
    PRINT*,''

    CASE("3","prepspec","PREPSPEC")
    PREP_ALL = .TRUE. 

    CASE("4","prep15","PREP15")
    PREP_15  = .TRUE. 

    CASE("5","prep13","PREP13")
    PREP_13  = .TRUE. 

    CASE("6","hotlocal","HOTLOCAL")
    HOT_LOCAL  = .TRUE. 

    CASE("7","hotglobal","HOTGLOBAL")
    HOT_GLOBAL  = .TRUE. 

    CASE DEFAULT  ! fall-through -> user can re-enter menu selection
    PRINT *, 'Input was misunderstood.'
    PRINT *, 'Please select number or action.'
    GO TO 9999

    END SELECT

    print *, "calling: prepinput"
    print *, "use_default = ", use_default
    print *, "partition = ", partition
    print *, "prep_all  = ", prep_all
    print *, "prep_15   = ", prep_15
    print *, "prep_13   = ", prep_13
    print *, "hot_local  = ", hot_local
    print *, "hot_global  = ", hot_global
    CALL PREPINPUT()
    print *, " "
    call memory_status()

    STOP
!---------------------------------------------------------------------------
    END PROGRAM ADCPREP
!---------------------------------------------------------------------------




!---------------------------------------------------------------------------
!               S U B R O U T I N E   P R E P I N P U T
!---------------------------------------------------------------------------
!     Read in the full domain input files, perform domain decomposition,
!     and write new subdomain input files.
!     This version is compatible with ADCIRC version 46.45. vjp 10/8/2006.
!---------------------------------------------------------------------------

    SUBROUTINE PREPINPUT()
    USE PRE_GLOBAL
    use sizes, only : naType, meshType, controlType, ASCII, XDMF, &
    naFileName, meshFileName, controlFileName, formatString
    use mesh, only : readMeshXDMF
    use subprep, only : subdomainOn, enforceBN, subPrep020, subPrep021, &
    subPrep019
! cm v50.85
#ifdef WINDOWS
    USE IFPORT
#endif
    IMPLICIT NONE

    INTEGER :: I           ! loop counter for nodes, elements
    INTEGER :: PE          ! loop counter for MPI processors
    CHARACTER CMD*6     ! string to hold shell command
    CHARACTER PENUM*6   ! string to hold directory name
    CHARACTER DIRCMD*72 ! string to hold complete shell command line
#ifdef WINDOWS
    logical(4) :: dir_result  !tcm v50.85
    INTEGER :: dir_error   !tcm v50.85
#endif

!-- If ONLY localizing Nodal-Attributes File ( fort.13 )
    if (PREP_13) then
        write(*,'(/)')
        write(*,*) 'Re-Writing subdomain Nodal Attributes (unit 13)'
        write(*,*) 'file for each PE.'
        CALL RELOCALIZE()
        CALL PREP13()
        GO TO 9999
    endif

!-- If ONLY processing hotstart File(s) ( fort.67-68 )
    if (HOT_LOCAL .OR. HOT_GLOBAL) then
        write(*,'(/)')
        write(*,*) 'Writing subdomain hotstart file for each subdomain'
        CALL RELOCALIZE()
        if (HOT_LOCAL)  CALL HOTLOCALIZE()
        if (HOT_GLOBAL) CALL HOTGLOBALIZE()
        GO TO 9999
    endif

!--   Read the Global Grid File ( Unit 14 )

    select case(meshType)
    case(ASCII)
    call sizeup14()
    if ( (partition.eqv. .FALSE. ) .AND. &
    (hot_local.eqv. .FALSE. ) .AND. &
    (hot_global.eqv. .FALSE. ) ) then
        call sizeup15()
    endif
    call alloc_main1()
    call prepRead14()  ! use the existing read14 from prep/read_global.f
    case(XDMF)
    call readMeshXDMF() ! from src/mesh.F
! copy mesh parameters from mesh module into pre_global module
    call getMeshParametersXDMF() ! defined in prep/read_global.F
    call sizeup15()
    case default
    write(6,'(a,i0,a)') 'ERROR: The mesh type ',meshType, &
    ' is invalid'
    stop
    end select
    print *, 'Global Grid file read successfully.'

! if XDMF output format was indicated for any file, then write
! the names of the fulldomain input files for later use by
! XDMF output routines
    if (useXDMF.eqv. .TRUE. ) then
        open(unit=30,file='fulldomainInputFiles', &
        status='replace',action='write')
        write(30,'(a)') 'naType'
        write(30,'(a)') trim(formatString(naType))
        write(30,'(a)') 'naFileName'
        write(30,'(a)') trim(naFileName)
        write(30,'(a)') 'meshType'
        write(30,'(a)') trim(formatString(meshType))
        write(30,'(a)') 'meshFileName'
        write(30,'(a)') trim(meshFileName)
        write(30,'(a)') 'controlType'
        write(30,'(a)') trim(formatString(controlType))
        write(30,'(a)') 'controlFileName'
        write(30,'(a)') trim(controlFileName)
        close(30)
    endif

!--   If ONLY localizing Run Info File ( fort.15 )
    if (PREP_15) then
        CALL READ15()
        print *, 'Global Run Info file read successfully.'
        CALL RELOCALIZE()
        CALL PREP15()
    ! jgf48.07 Initialize netCDF output files
        IF (useNetCDF) THEN
            write(*,*) "INFO: Initializing the netCDF output files."
        ! call netCDF module
            CALL PREPNETCDF()
        ENDIF
        GO TO 9999
    endif

!   If ONLY decomposing a fort.20 (variable flux) file.
!   -- Write a Local Flux File ( fort.20 ) for each PE if needed
!   -- EXIST_FLUX and other variables are set when the fort.14 is
!   -- initially read in within the READ14() subroutine.
    IF (PREP_20.eqv. .TRUE. ) THEN
        print *, 'INFO: Start reading in fort.15.'
        CALL READ15() ! need to determine value of APERIODIC_FLOW_BC
        print *, 'INFO: Fulldomain fort.15 file read successfully.'
        if (EXIST_FLUX /= 0 .AND. APERIODIC_FLOW_BC) then
            WRITE(*,*) &
            "INFO: Decomposing fort.20 (variable flux) file."
            CALL RELOCALIZE()
            CALL PREP20()
            GOTO 9999
        endif
    ENDIF

!-- Partition Nodes with METIS 4.0 Graph Partition Library

    IF (PARTITION) THEN
        CALL METIS()
        print *, " "
        print *, 'INFO: METIS has partitioned nodes'
        GO TO 9999
    ELSE
        print *, 'INFO: Opening file partmesh.txt'
        OPEN(990,FILE='partmesh.txt')
        DO I=1, MNP
            READ(990,*) PROC(I)
        !           print *," PROC for ",I," = ", PROC(I)
        ENDDO
        CLOSE(990)
        print *, 'INFO: Closed partmesh.txt file.'
    ENDIF

!-- Read the Global Input File ( Unit 15 )

    print *, 'INFO: Start reading in fort.15.'
    CALL READ15()
    print *, 'INFO: Fulldomain fort.15 file read successfully.'

!     jgf45.07 In the case of starting a parallel job from a full domain
!     hot start file, we must reverse the signs of NOUTM, NOUTV, etc.
!     This is because ADCIRC wouldn't be able to find preexisting
!     output files in the subdomains to append to.

    IF (PREP_ALL .AND. (IHOT == 67 .OR. IHOT == 68)) THEN
        print *, 'INFORMATION: Subdomain output files will be started'
        print *, '             anew rather than appended.'
        IF (NOUTE > 0 .AND. ABS(NOUTE) /= 3) NOUTE =-NOUTE
        IF (NOUTV > 0 .AND. ABS(NOUTV) /= 3) NOUTV =-NOUTV
        IF (NOUTC > 0.) NOUTC=-NOUTC
        IF (NOUTM > 0 .AND. ABS(NOUTM) /= 3) NOUTM =-NOUTM
        IF (NOUTGE > 0 .AND. ABS(NOUTGE) /= 3) NOUTGE=-NOUTGE
        IF (NOUTGV > 0 .AND. ABS(NOUTGV) /= 3) NOUTGV=-NOUTGV
        IF (NOUTGC > 0.                    ) NOUTGC=-NOUTGC
        IF (NOUTGW > 0 .AND. ABS(NOUTGW) /= 3) NOUTGW=-NOUTGW
    !         WRITE(6,*)'==NOUT==',NOUTE,NOUTV,NOUTM,NOUTGE,NOUTGV
    !   kmd48.33bc add information for 3D
        IF (C3DVS) THEN
            IF (I3DSD > 0) I3DSD=-I3DSD
            IF (I3DSV > 0) I3DSV=-I3DSV
            IF (I3DST > 0) I3DST=-I3DST
            IF (I3DGD > 0) I3DGD=-I3DGD
            IF (I3DGV > 0) I3DGV=-I3DGV
            IF (I3DGT > 0) I3DGT=-I3DGT
        ENDIF

    ENDIF

! jgf48.07 Initialize netCDF output files
    IF ((useNetCDF.eqv. .TRUE. ) .AND. (PREP_ALL.eqv. .TRUE. )) THEN
        write(*,*) "INFO: Initializing the netCDF output files."
    ! call netCDF module
        CALL PREPNETCDF()
    ENDIF

!-- Decompose the ADCIRC grid into MNPROC subdomains

    print *, " "
    print *, "Determine the parameters MNPP and MNEP"
    CALL DOMSIZE()

    print *, "Allocate arrays dimensioned by MNPP and MNEP"
    CALL ALLOC_MAIN2()

    print *, " "
    print *, "Decomposition of grid begins"
    CALL DECOMP()
    print *, "Decomposition successful"

!-- Create MNPROC sub-directories of the working directory

    DO PE=0, MNPROC-1
        PENUM  = 'PE0000'
        CALL IWRITE(PENUM,3,6,PE)
    ! cm 20110510 v50.05 adding support for windows
#ifdef WINDOWS
        dir_result = makedirqq(PENUM) !tcm v50.85 added
        if (dir_result.eqv. .FALSE. ) then
            dir_error = getlasterrorqq()
            SELECT CASE (dir_error)
            CASE(2)
            print *, "Error Creating ",PENUM, &
            ". The file or path specified was not found."
            stop
            CASE(13)
            print *, "Error Creating ",PENUM,".  Access Denied."
            stop
            CASE(17)
            print *, PENUM," already exists."
            CASE DEFAULT
            END SELECT
        ENDIF
#elif PC_DIG_FORT
        call  system('mkdir '//PENUM)  !tcm v50.85 this no longer works for INTEL MPI
#else
        CALL MAKEDIR(PENUM)
#endif
    ENDDO

!-- Write a Local Grid File ( fort.14 ) for each PE

    print *, "Writing Local UNIT 14 (Grid) File for each PE"
    CALL PREP14()

!-- Write a Local Input file ( fort.15 ) for each PE

    print *, "Writing Local UNIT 15 (Run Info) File for each PE"
    CALL PREP15()

!     jgf48.17 If compiling for ADCSWAN, need to copy UnSWAN file
! Casey 090304: Moved this call from earlier in this file.
#ifdef ADCSWAN
    CALL PREPUNSWAN()
#endif

!-- Write Message-Passing File for each PE

    print *, "Writing Message-Passing Info Files for each PE"
    IF (PREP_ALL) CALL PREP18()

!     Write a subdomain initial concentration file ( fort.10 ) for each PE

    IF (PREP_ALL) THEN
        IF (C2D_PTrans .OR. C3D_PTrans) THEN
            PRINT *, "Writing subdomain UNIT 10 file for each PE."
            CALL PREP10()
        ENDIF
    ENDIF

!     Write a subdomain initial density file ( fort.11 ) for each PE

!   kmd48.33bc changed the if statement to include IHOT parameter
    IF (PREP_ALL) THEN
        IF ((CBaroclinic) .AND. (IHOT == 0)) THEN
            PRINT *, "Writing subdomain UNIT 11 file for each PE."
            CALL PREP11()
        ENDIF
    ENDIF

!   kmd48.33bc added writing of subdomain initial condition file

!     jgf50.05.01: Commented out b/c we now (as of v49release) map
!     this file to subdomains on the fly during padcirc execution.
!     Write a subdomain initial condtion file (fort.17) for each PE

!      IF (PREP_ALL) THEN
!         IF ((CBaroclinic).AND.(IHOT.EQ.17)) THEN
!            PRINT *, "Writing subdomain UNIT 17 file for each PE."
!            CALL HOTINITCOND()
!         ENDIF
!      ENDIF


!     Write a subdomain nodal attributes file (fort.13) for each PE

    if (NWP > 0) then
        write(*,'(/)')
        write(*,*) 'Writing subdomain Nodal Attributes (unit 13)'
        write(*,*) 'file for each PE.'
        CALL PREP13()
    endif

!-- If required write a Local fort.19 file for each PE


    if (subdomainOn) then                     ! NCSU Subdomain
        if (enforceBN == 1) call subprep019() ! NCSU Subdomain
        if (enforceBN == 2) call subprep020() ! NCSU Subdomain
        if (enforceBN == 2) call subprep021() ! NCSU Subdomain
    else                                       ! NCSU Subdomain
        IF (PREP_ALL) THEN
            IF ((NOPE > 0) .AND. (NBFR == 0)) THEN
                print *, "Ready to Write Local UNIT 19 File for each PE"
                CALL PREP19()
            ENDIF
        ENDIF
    endif ! NCSU Subdomain

!-- Write a Local Flux File ( fort.20 ) for each PE if needed
!   -- EXIST_FLUX and other variables are set when the fort.14 is
!   -- initially read in within the READ14() subroutine.

    IF (PREP_ALL) THEN
        if (EXIST_FLUX /= 0 .AND. APERIODIC_FLOW_BC) then
            CALL PREP20()
        endif
    ENDIF

!-- If required write a Local Wind Stress file for each PE

    IF (PREP_ALL) THEN
    !     sb46.28sb01 NWS=12 doesn't need fort.22 decomposition
    !     tcm_v49.04 NWS=3 and NWS=6 no longer need fort.22 decomposition
    !     jgf49.0804 NWS=29 does not need fort.22 decomposition.
    !     tcm v51.06.02 added NWS=16 GFDL Met Data
    !         IF (NWS.NE.0.AND.ABS(NWS).NE.12)  THEN
        IF ((NWS /= 0) .AND. (ABS(NWS) /= 3) .AND. &
        (ABS(NWS) /= 6) .AND. (ABS(NWS) /= 12) &
         .AND. (ABS(NWS) /= 16) &
         .AND. (ABS(NWS) /= 29))  THEN
            print *, "Ready to Write Local UNIT 22 File for each PE"
        !     jgf46.00 Added NWS=7
        !     jgfdebug46.02 Added NWS=45
        !     jgf46.02 Added NWS=8
        !     jgf46.16 Merged:
        !     cf & cm added NWS=9: asymmetric hurricane wind model
        !     rjw added NWS=19: asymmetric hurricane wind model v2.0
        !     jie added NWS=20: generalized asymmetric vortex model
        !     tcm_v49.04 NWS=3 and NWS=6 no longer need fort.22 decomposition
            IF((NWS == 1) .OR. (ABS(NWS) == 2) .OR. (ABS(NWS) == 4) .OR. &
            (ABS(NWS) == 5) .OR. (NWS == 7) .OR. &
            (ABS(NWS) == 45)) THEN
                CALL PREP22()
            !           ELSEIF (NWS.EQ.10.OR.NWS.EQ.11) THEN
            !              CALL PREP200()
            ENDIF
        ENDIF
    ENDIF

!-- If required write a Local Wave Stress file for each PE

    IF (PREP_ALL) THEN
        IF(NRS == 1) CALL PREP23()
    ENDIF

!... tcm v50.66.03 Added for time varying bathymetry (fort.141)

!-- If required write a Local Time Bathymetry file for each PE

    IF (PREP_ALL) THEN
        IF ((ABS(NDDT) == 1) .OR. (ABS(NDDT) == 2)) THEN
            CALL PREP141()
        ENDIF
    ENDIF


!  kmd48.33bc add information for 3D baroclinic simulations
    IF (PREP_ALL) THEN
        IF (CBAROCLINIC) THEN
            IF ((RES_BC_FLAG >= 1) .OR. (RES_BC_FLAG <= -1)) THEN
                IF (NOPE > 0) THEN
                    IF (BCFLAG_LNM == 1) THEN
                        print *, "Ready to Write Local UNIT 35 File for each PE"
                        CALL PREP35()
                    END IF
                END IF
            END IF
            IF (RES_BC_FLAG == 2) THEN
                IF (NOPE > 0) THEN
                    print *, "Ready to Write Local UNIT 36 File for each PE"
                    CALL PREP36()
                END IF
            ELSE IF (RES_BC_FLAG == 3) THEN
                IF (NOPE > 0) THEN
                    print *, "Ready to Write Local UNIT 37 File for each PE"
                    CALL PREP37()
                END IF
                IF (BCFLAG_TEMP /= 0) THEN
                    print *, "Ready to Write Local UNIT 38 File for each PE"
                    CALL PREP38()
                END IF
            ELSE IF (RES_BC_FLAG == 4) THEN
                IF (NOPE > 0) THEN
                    print *, "Ready to Write Local UNIT 36 File for each PE"
                    CALL PREP36()
                    print *, "Ready to Write Local UNIT 37 File for each PE"
                    CALL PREP37()
                END IF
                IF (BCFLAG_TEMP /= 0) THEN
                    print *, "Ready to Write Local UNIT 38 File for each PE"
                    CALL PREP38()
                END IF
            END IF
        END IF
    END IF

! kmd - added for rivers in a baroclinic simulation
    IF (PREP_ALL) THEN
        IF ((EXIST_BC_TS /= 0) .AND. (APERIODIC_BC_TS)) THEN
            print *, "Ready to Write Local UNIT 39 File for each PE"
            CALL PREP39()
        END IF
    END IF


!-- Save domain-decomposition information for post-processor

!C -- Start Addition by CF  8/2007
    IF (PREP_ALL) THEN
        print *, "Writing domain-decomposition file for post-processor"
        CALL PREP80()
    ENDIF
!C -- Finish Addition by CF

    9999 CONTINUE
    print *, ""
    print *, "INFO: Finished pre-processing input files."
    RETURN
    END SUBROUTINE
!---------------------------------------------------------------------------
!     End of subroutine prepinput
!---------------------------------------------------------------------------




    SUBROUTINE GETMSG( STRING, MSG )
    INTEGER :: I, I1
    CHARACTER  STRING*(*),MSG*(*), TARGET

    I1 = 0
    TARGET = "!"

!-- Find beginning of message

    DO I=1, 80
        IF (STRING(I:I) == TARGET) THEN
            I1 = I
            GOTO 100
        ENDIF
    ENDDO

    100 CONTINUE

!--Copy message to ouput string

    DO I=1, I1-1
        MSG(I:I) = " "
    ENDDO
    MSG(I1:80)  = STRING(I1:80)

    RETURN
    END SUBROUTINE GETMSG



    SUBROUTINE NEWINDEX ( ISTRING, OSTRING, INDX )
    INTEGER :: I,I1,I2,I3,I4,INDX
    CHARACTER  ISTRING*(*),OSTRING*(*),TARGET
    CHARACTER TEMP1*80, TEMP2*100

    I1 = 0
    I2 = 0
    I3 = 0
    I4 = 0
    TARGET = " "

!-- Find first non-blank character of String

    DO I=1, 80
        IF (ISTRING(I:I) /= TARGET) THEN
            I1 = I
            GOTO 100
        ENDIF
    ENDDO

!-- Find next blank character of String

    100 CONTINUE
    DO I=I1+1,80
        IF (ISTRING(I:I) == TARGET) THEN
            I2 = I
            GOTO 200
        ENDIF
    ENDDO

    200 CONTINUE

!-- Create a temporary string containing new index

    WRITE(TEMP1(1:80),'(I8)') INDX

!-- Find first non-blank character of String

    DO I=1, 80
        IF (TEMP1(I:I) /= TARGET) THEN
            I3 = I
            GOTO 300
        ENDIF
    ENDDO

!-- Find next blank character of String

    300 CONTINUE
    DO I=I3+1,80
        IF (TEMP1(I:I) == TARGET) THEN
            I4 = I
            GOTO 400
        ENDIF
    ENDDO

    400 CONTINUE

! ebug print *, "i1 i2 i3 i4 ",I1, I2, I3 , I4
    TEMP2(1:100) = TEMP1(I3-1:I4-1)//ISTRING(I2:80)

!-- Write out first 80 characters of concatenated strings

    OSTRING(1:80) = TEMP2(1:80)

    RETURN
    END SUBROUTINE NEWINDEX


    SUBROUTINE INSERT( ISTRING, OSTRING, NUMS, N )
    INTEGER :: I,J,I1,N,NUMS(N)
    CHARACTER  ISTRING*80,OSTRING*80,BLANK
    CHARACTER  TEMP1*80, TEMP2*160

    I1 = 0
    BLANK = " "

!-- Create Tempoarary String TEMP1 containing NUMS

    IF (N == 1) THEN
        WRITE(TEMP1(1:80),'(I8)') NUMS(1)
    ! Casey 090304: Changed this section to allow N = 3.
    ELSEIF (N == 2) THEN
        WRITE(TEMP1(1:80),'(2I8)') NUMS(1),NUMS(2)
    ELSE
        WRITE(TEMP1(1:80),'(3I8)') NUMS(1),NUMS(2),NUMS(3)
    ENDIF

!-- Find length of TEMP1 string

    DO I=80,1,-1
        IF (TEMP1(I:I) /= BLANK) THEN
            LEN1 = I
            GOTO 10
        ENDIF
    ENDDO
    10 CONTINUE

!-- Scan input string for character after old number list

    I = 1
    DO NUM=1, N
        DO J=I,80
            IF (ISTRING(J:J) /= BLANK) THEN
                I = J
                GOTO 100
            ENDIF
        ENDDO
        100 CONTINUE
        DO J=I,80
            IF (ISTRING(J:J) == BLANK) THEN
                I = J
                GOTO 200
            ENDIF
        ENDDO
        200 CONTINUE
    ENDDO
    I1 = MAX(0,I)

!-- Insert Integer List into Message

    IF (I1 /= 0) THEN
    !-- if there is a message
        TEMP2(1:160) = TEMP1(1:LEN1+1)//ISTRING(I1:80)
    ELSE
        TEMP2(1:160) = TEMP1(1:LEN1+1)
    ENDIF

!-- Write out first 80 characters of concatenated string

    OSTRING(1:80) = TEMP2(1:80)

    RETURN
    END SUBROUTINE INSERT
