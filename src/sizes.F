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

      CHARACTER(256),  SAVE :: ROOTDIR
      CHARACTER(2048), SAVE :: INPUTDIR, GLOBALDIR, LOCALDIR
      CHARACTER(2048), SAVE :: GBLINPUTDIR, HOTSTARTDIR
      
C---------------------end of data declarations--------------------------------C


      CONTAINS


      SUBROUTINE MAKE_DIRNAME
      INTEGER :: LNAME, IARGC, ARGCOUNT, I, iprefix
      CHARACTER(2048) :: CMDLINEARG 
      CHARACTER(8)    :: PREFIX(2) = (/ '/PE0000 ', '/DOM0000' /)
      logical         :: fileFound
      
      INPUTDIR  = ""
      GLOBALDIR = ""
      LOCALDIR  = ""
      ARGCOUNT  = IARGC()

      IF ( MNPROC.EQ.1 ) THEN
         WRITE_LOCAL_FILES = .TRUE.
         WRITE_LOCAL_HOT_START_FILES = .TRUE.
      ENDIF

      if (ARGCOUNT > 0) then
        i = 0
        do while (i < ARGCOUNT)
          i = i + 1
          call getarg(i, CMDLINEARG)
          if (cmdlinearg(1:2) == "-I") then
            i = i + 1
            call getarg(i,INPUTDIR)
          else if (cmdlinearg(1:2) == "-O") then
            i = i + 1
            call getarg(i,GLOBALDIR)
          else if (cmdlinearg(1:2) == "-L") then
            WRITE_LOCAL_FILES = .TRUE.
            WRITE_LOCAL_HOT_START_FILES = .TRUE.
          endif
        end do
      end if

C.....Default root working directory

      if (INPUTDIR == "") then
        ROOTDIR = '.'
      else
        ROOTDIR = INPUTDIR
      endif

      GBLINPUTDIR = ROOTDIR
#ifdef CMPI
      iprefix = 0
      do i = 1, 2
        INQUIRE(file=TRIM(ROOTDIR)//TRIM(PREFIX(i))//'/'//'fort.14',
     &    exist=fileFound)
        if (fileFound) then
          iprefix = i
          exit
        end if
      end do
      if (.not. fileFound) then
        print *, "Failed to find prefix directory"
        call msg_abort()
      end if

      WRITE(INPUTDIR,'(2A)') TRIM(ROOTDIR),PREFIX(iprefix)
      LNAME = LEN_TRIM(INPUTDIR)
      WRITE(INPUTDIR(LNAME-3:LNAME),'(I4.4)') MYPROC
#else
      WRITE(INPUTDIR,'(A)') TRIM(ROOTDIR)
#endif

      if (GLOBALDIR == "") then
        ROOTDIR = '.'
      else
        ROOTDIR = GLOBALDIR
      endif

      WRITE(GLOBALDIR,'(A)') TRIM(ROOTDIR)


#ifdef CMPI
      WRITE(LOCALDIR,'(2A)') TRIM(ROOTDIR),TRIM(PREFIX(iprefix))
      LNAME = LEN_TRIM(LOCALDIR)
      WRITE(LOCALDIR(LNAME-3:LNAME),'(I4.4)') MYPROC
      call MAKEDIR(trim(LOCALDIR))
#else
      WRITE(LOCALDIR,'(A)') TRIM(ROOTDIR)
#endif
      if (WRITE_LOCAL_FILES) GLOBALDIR = LOCALDIR
      HOTSTARTDIR = LOCALDIR

      
      RETURN
      END SUBROUTINE MAKE_DIRNAME

      SUBROUTINE GET_NUMWRITERS
      INTEGER :: IARGC, ARGCOUNT, I
      CHARACTER(2048) :: CMDLINEARG 
      
      ARGCOUNT  = IARGC()

      MNWPROC = 0

      if (ARGCOUNT > 0) then
        i = 0
        do while (i < ARGCOUNT)
          i = i + 1
          call getarg(i, CMDLINEARG)
          if (cmdlinearg(1:2) == "-W") then
            i = i + 1
            call getarg(i,cmdlinearg)
            read(cmdlinearg,*) MNWPROC
          endif
        end do
      end if

C.....Default root working directory

      RETURN
      END SUBROUTINE GET_NUMWRITERS

      END MODULE SIZES
