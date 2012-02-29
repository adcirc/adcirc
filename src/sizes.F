C******************************************************************************
C PADCIRC VERSION 45.12 03/17/2006                                            *
C  last changes in this file VERSION 45.56                                    *
C                                                                             *
C******************************************************************************
C
      MODULE SIZES
      IMPLICIT NONE
C
C...SET NUMBER OF BYTES "SZ" IN REAL(SZ) DECLARATIONS
C...SET "NBYTE" FOR PROCESSING INPUT DATA RECORD LENGTH

#ifdef REAL4
      INTEGER, PARAMETER :: SZ = 4
      INTEGER, PARAMETER :: NBYTE=4
#endif
#ifdef REAL8
      INTEGER, PARAMETER :: SZ = 8
      INTEGER, PARAMETER :: NBYTE=8
#endif

C...SET MAX OF DIGITS OF PRECISION "NPREC" THE GRID CAN BE EXPECTED TO HAVE
C...NOTE: IF THE GRID WAS BUILT ON A 32 BIT COMPUTER, IT SHOULD BE
C   ACCURATE TO ABOUT 7 DIGITS.  THUS, IF THE NODAL SPACING REQUIRES MORE
C   THAN 5 DIGITS OF PRECISION, THE MODEL RESULTS MAY NOT BE TRUSTWORTHY.

      INTEGER, PARAMETER ::  NPREC=7
C
      INTEGER ::  MNE,MNP,MNEI,MNOPE,MNETA,MNBOU,MNVEL,
     &  MNTIF,MNBFR,MNFFR,MNSTAE,MNSTAV,MNSTAC,MNSTAM,MNHARF
      INTEGER :: MNPROC         ! Number of compute processors
      INTEGER :: MNWPROC        ! Number of writer processors
      INTEGER :: MNWPROH        ! Number of Hwriter processors !st3 100711
      INTEGER :: ISPLIT         ! ON(1)/OFF(0) parameter of splitting output for writer !st3 100708:
      INTEGER :: MNALLPROC      ! Number of all processors (= MNPROC + MNWPROC)
C
C     Dimension of vertical FE mesh (To interface 2D & 3D)
c     INTEGER :: MNODES                                        !changed to MNFEN 08/16/2005
      INTEGER :: MNFEN

c     LOGICAL C2DDI,C3D,C3DDSS,C3DVS,CLUMP,CTIP,CHARMV         !moved to GLOBAL.F  04/09/2004
C
C For Definition of Working Directory
C
      INTEGER,SAVE :: MYPROC
      LOGICAL      :: WRITE_LOCAL_FILES
      LOGICAL      :: WRITE_LOCAL_HOT_START_FILES
      ! jgf49.44: Make it an option to read local binary hotstart files
      LOGICAL      :: READ_LOCAL_HOT_START_FILES

      CHARACTER(256),  SAVE :: ROOTDIR
      CHARACTER(2048), TARGET, SAVE :: INPUTDIR, GLOBALDIR, LOCALDIR
      CHARACTER(2048), SAVE :: GBLINPUTDIR, HOTSTARTDIR

C---------------------end of data declarations--------------------------------C


      CONTAINS


      SUBROUTINE MAKE_DIRNAME()
      IMPLICIT NONE
      INTEGER :: LNAME, IARGC, ARGCOUNT, I, iprefix
      CHARACTER(2048) :: CMDLINEARG
      CHARACTER(8)    :: PREFIX(2) = (/ '/PE0000 ', '/DOM0000' /)
      logical         :: fileFound

      INPUTDIR  = ""
      GLOBALDIR = ""
      LOCALDIR  = ""
      ARGCOUNT  = IARGC()
      WRITE_LOCAL_HOT_START_FILES = .FALSE.
      READ_LOCAL_HOT_START_FILES = .FALSE.
      WRITE_LOCAL_FILES = .FALSE.

Casey 090527: Debug.
      ROOTDIR = "."

Cobell 110518: Add compiler flag for local hot start.
#ifdef LOCALHOT
      WRITE_LOCAL_HOT_START_FILES = .TRUE.
#else
      WRITE_LOCAL_HOT_START_FILES = .FALSE.
#endif

      IF ( MNPROC.EQ.1 ) THEN
         WRITE_LOCAL_FILES = .TRUE.
         WRITE_LOCAL_HOT_START_FILES = .TRUE.
         READ_LOCAL_HOT_START_FILES = .TRUE.
      ENDIF

      if (ARGCOUNT > 0) then
C         write(*,*) "DEBUG: ",myProc,
C     &       ": Processing command line options."
         i = 0
         do while (i < ARGCOUNT)
            i = i + 1
            call getarg(i, CMDLINEARG)
            select case(cmdlinearg(1:2))
            case("-I")
               i = i + 1
               call getarg(i,INPUTDIR)
               if (myProc.eq.0) then
                  write(*,*) "INFO: Processing '-I ",trim(INPUTDIR),"'."
               endif
            case("-O")
               i = i + 1
               call getarg(i,GLOBALDIR)
               if (myProc.eq.0) then
                 write(*,*) "INFO: Processing '-O ",trim(GLOBALDIR),"'."
               endif
            case("-L")
               WRITE_LOCAL_FILES = .TRUE.
               WRITE_LOCAL_HOT_START_FILES = .TRUE.
               if (myProc.eq.0) then
                  write(*,*) "INFO: Processing '-L'."
               endif
            case("-S")
               WRITE_LOCAL_HOT_START_FILES = .TRUE.
               if (myProc.eq.0) then
                  write(*,*) "INFO: Processing '-S'."
               endif
            case("-R")
               READ_LOCAL_HOT_START_FILES = .TRUE.
               if (myProc.eq.0) then
                  write(*,*) "INFO: Processing '-R'."
               endif
            case("-W")
               i = i + 1
               ! this is the number of writer processors, and will
               ! be parsed in a later call to GET_NUMWRITERS
               !  ... could also be -Ws for round-robin writers
            case default
               if (myProc.eq.0) then
                  write(*,*) "WARNING: The command line option '",
     &             cmdlinearg(1:2),"' is not valid and will be ignored."
               endif
            end select
         end do
      end if

C.....Default root working directory

      if (len_trim(INPUTDIR) /= 0) then
        ROOTDIR = INPUTDIR
      endif

      GBLINPUTDIR = ROOTDIR
#ifdef CMPI
      iprefix = 0
      if (myProc.eq.0) then
         write(*,*) "INFO: Searching for ADCIRC subdomain directories:"
      endif
      do i = 1, 2
         if (myProc.eq.0) then
            write(*,*) "INFO: Looking for '",trim(ROOTDIR),
     &         trim(PREFIX(i)),"/fort.14' ..."
         endif
         INQUIRE(file=TRIM(ROOTDIR)//TRIM(PREFIX(i))//'/'//'fort.14',
     &      exist=fileFound)
         if (fileFound.eqv..TRUE.) then
            if (myProc.eq.0) then
               write(*,*) "INFO: File '",trim(ROOTDIR),
     &            trim(PREFIX(i)),"/fort.14' was found!"
               write(*,*) "INFO: The search for the subdomain ",
     &            "directory was completed successfully."
            endif
            iprefix = i
            exit
         else
            write(*,*) "ERROR: Processor ",myProc,": File '",
     &         trim(ROOTDIR),trim(PREFIX(i)),
     &         "/fort.14' was not found."

            print *, "ERROR: ADCIRC stopping."
            call msg_abort()
         end if
      end do
C
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

C... tcm v49.67 -- additions for Windows based PC compilations
#ifdef WINDOWS
      call system('mkdir '//trim(LOCALDIR))
#else
      call MAKEDIR(trim(LOCALDIR))
#endif

#else
      WRITE(LOCALDIR,'(A)') TRIM(ROOTDIR)
#endif
      if (WRITE_LOCAL_FILES) GLOBALDIR = LOCALDIR
      HOTSTARTDIR = LOCALDIR
C
C     jgf49.17.01 Summarize and log results.
      if (myProc.eq.0) then
         write(*,*) "INFO: The ROOTDIR is  '",trim(ROOTDIR),"'."
         write(*,*) "INFO: The INPUTDIR is '",trim(INPUTDIR),"'."
         write(*,*) "INFO: The GBLINPUTDIR is '",trim(GBLINPUTDIR),"'."
         write(*,*) "INFO: The GLOBALDIR is '",trim(GLOBALDIR),"'."
         write(*,*) "INFO: The LOCALDIR is '",trim(LOCALDIR),"'."
         write(*,*) "INFO: The HOTSTARTDIR is '",trim(HOTSTARTDIR),"'."
         if ((WRITE_LOCAL_FILES.eqv..TRUE.).AND.(MNPROC.GT.1)) then
            write(*,*) "INFO: ADCIRC will write subdomain output files."
         endif
         if ((WRITE_LOCAL_HOT_START_FILES.eqv..TRUE.)
     &                            .AND.(MNPROC.GT.1)) then
            write(*,*)
     &         "INFO: ADCIRC will write subdomain hotstart files."
         endif
      endif
      RETURN
      END SUBROUTINE MAKE_DIRNAME

      SUBROUTINE GET_NUMWRITERS
      IMPLICIT NONE
      INTEGER :: IARGC, ARGCOUNT, I
      CHARACTER(2048) :: CMDLINEARG
      ARGCOUNT  = IARGC()
      MNWPROC = 0
      isplit = 0 !st3 100708:split file

      if (ARGCOUNT > 0) then
        i = 0
        do while (i < ARGCOUNT)
          i = i + 1
          call getarg(i, CMDLINEARG)
          if (cmdlinearg(1:2) == "-W") then
            if( cmdlinearg(1:3) == "-Ws" ) then  !st3 100708: split file
               isplit = 1                        !st3 100708: split file
            endif                                !st3 100708: split file
            i = i + 1
            call getarg(i,cmdlinearg)
            read(cmdlinearg,*) MNWPROC
          endif
        end do
      end if
!
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

C.....Default root working directory

      RETURN
      END SUBROUTINE GET_NUMWRITERS

      END MODULE SIZES
