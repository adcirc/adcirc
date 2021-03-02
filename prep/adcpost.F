C----------------------------------------------------------------------------
C
C                 ADCPOST:  Parallel ADCIRC Post-Processor
C
C
C                   current to ADCIRC v45.12   03/17/2006
C----------------------------------------------------------------------------
C
C
C                      Serial Version 1.1 ( vjp 5/04/99 )
C                      Serial Version 1.3 ( meb 4/01/04 ) 
C                      Serial Version 1.4 ( jgf 11/08/2005 ) 
C                      Serial Version 1.5 ( jgf 02/02/2006 ) 
C
C  Program Development History
C  ---------------------------
C   Written for ADCIRC_v24.05   ( S. Chippada 1996 )
C   Updated for ADCIRC_v24.05   ( M. Martinez 1997 )
C   Hilbert Space Filling Curve ( C. Edwards 1997 ) 
C   Updated for ADCIRC_v33.04   ( V. Parr 1998 ) 
C   Updated for ADCIRC_v34.08   ( V. Parr 1999 ) 
C   Updated to gather hotstart files    ( R. Luettich 10/01)
C   Added processing of ind. files   ( M. Brown    04/04 )
C   Updated for ADCIRC v45.07   ( J. Fleming 11/05 )
C   3D Updated for ADCIRC v45.11   ( J. Fleming 02/2006 )
C
C----------------------------------------------------------------------------
C
C  ADCPOST performs 3 operations:
C
C    post      =   globalize output files
C    compare   =   compare global output files in two directories
C    diffmerge =   merge difference of outputs in two directories
C
C  Post-processor
C  --------------
C  The Post-processor "post" reads the local output files and globalizes them,
C  placing them in the working directory producing the same  global output
C  files as a normal non-parallel run.
C
C  Compare
C  -------
C  Compares global output files from two different directories using
C  a user supplied tolerance. A log of the differences for Unit xx is
C  written to the current directory in file "diffs.xx".   This utility
C  is useful for validation of PADCIRC on a new computer platform, 
C  or comparing the results from a serial and a parallel run.
C
C  Diffmerge
C  ---------
C  Computes the scaled difference of fort.63 and fort.64 from two different 
C  directories and writes the difference in the same format as those
C  files and calls them diffmerge.63 and diffmerge.64, respectively.
C  The scale factor is supplied by the user.
C
C----------------------------------------------------------------------------
C
      PROGRAM ADCPOST
      USE POST_GLOBAL
C
      INTEGER UNIT   ! user input for file to post-process
      CHARACTER DUM  ! variable used in waiting for menu choice  
      INTEGER I,PE,PRECIS,ISP,NWTYPE

      REAL(SZ) TOL,SCALEVAL
      LOGICAL FOUND
      CHARACTER DIRCMD*72,CMD*6,PENUM*6,PHASE*9,ISET*6,DIR1*80,DIR2*80,
     &     CONVDIR*6,HOTANS*1

      integer :: argcount ! number of command line arguments
      integer :: iargc    ! function to return command line arguments
      character(len=2048) :: cmdlinearg ! content of cmd line arg
      character(len=2048) :: filename   ! e.g. "fort.63" 
C
C--Initialize constants for case ICS = 2 inputs
C
      R  =  6378206.4D0
      DEG2RAD = 3.141592653589793D0/180.0D0
      RAD2DEG = 180.0D0/3.141592653589793D0
C
      print *, "initializing post-processor"
      CALL POST_INIT()
      !
      ! parse command line option
      argcount = iargc() ! count up command line options
      if (argcount.gt.0) then
         i=0
         do while (i.lt.argcount)
            i = i + 1
            call getarg(i, cmdlinearg)
            select case(trim(cmdlinearg))
            case("--filename")
               i = i + 1
               call getarg(i,filename)
               write(6,'(a)') 'INFO: Processing "'//trim(cmdlinearg)//
     &           ' '//trim(filename)//'".'
               call mergeSubdomainFiles(filename) 
            case default
               write(*,*) 'WARNING: The command line option "',
     &         trim(cmdlinearg),'" is not valid and will be ignored.'
            end select
         end do
         CALL EXIT(1)
      endif
C
C--Say Hello Gracie
C
      print *," ***********************************"
      print *," ADCPOST version 1.5 (03/17/2006)"
      print *," Parallel ADCIRC Post-processor"
      print *," ***********************************"
99    print *, " "
      print *," Select operation "
      print *, " "
      print *," post      = Globalize output files"
      print *," compare   = Compare outputs in two directories"
      print *," diffmerge = Merge difference in two directories"
      print *," quit      = Exit program"
C
      read(*,'(A7)') PHASE
C
C--set NPROC = MNPROC
C
      NPROC = MNPROC
C
C---------------------------------------------------------------------------
C  Parallel ADCIRC Post-Processor Starts here
C---------------------------------------------------------------------------
C
      UNIT = 0
 100  IF (PHASE(1:4).EQ.'post') THEN
C     jgf45.06 BEGIN block of code written by MEB04/01/04 - - - - - - - - - 
C         call ISHELL("clear")
         print*, ""
         print*, "-----------------------"
         print*, "Post Processing choices"
         print*, "-----------------------"
         print*, "Enter Unit number to Post Process"
         if (UNIT .NE. 100) then
            print*, "Enter 100 to process all applicable unit numbers"
         endif
         print*, "Enter 555 to see list of Unit Numbers"
         print*, "Enter 999 to End Processing"
         read*, UNIT
C     
         if (UNIT == 555) then
            print*,'-------------------------------------------------'
            print*,'                 2D output files'
            print*,''
            print*,'51 Harm. Station Elevation  61 Station Elevation'
            print*,'52 Harm. Station Velocity   62 Station Velocity'
            print*,'53 Harm. Global  Elevation  63 Global  Elevation'
            print*,'54 Harm. Global  Velocity   64 Global  Velocity'
            print*,''
            print*,'71 Station Atm. Pres.       73 Global Atm. Pres.'
            print*,'72 Station Wind Vel.        74 Global Wind Vel.'
            print*,''
            print*,'91 Ice Conc. Station      93 Global Ice Conc.'
            print*,''
            print*,'                 3D output files'
            print*,''
            print*,'41 Station Density          44 Global Density'
            print*,'42 Station 3d Velocity      45 Global 3d Velocity'
            print*,'43 Station Turbulence       46 Global Turbulence'
            print*,'-------------------------------------------------'
            print*,'Press return to continue'
            read(*,'(A)') DUM
            goto 100
         endif

         if ((UNIT == 100) .or. (UNIT == 51)) then
            IF (NHASE.NE.0) THEN
               print *, "Post-Processing Unit 51"
               CALL POST51()
            ENDIF
         endif
C     
         if ((UNIT == 100) .or. (UNIT == 52)) then
            IF (NHASV.NE.0) THEN
               print *, "Post-Processing Unit 52"
               CALL POST52()
            ENDIF
         endif
C     
         if ((UNIT == 100) .or. (UNIT == 53)) then
            IF (NHAGE.NE.0) THEN
               print *, "Post-Processing Unit 53"
               CALL POST53()
            ENDIF
         endif
C     
         if ((UNIT == 100) .or. (UNIT == 54)) then
            IF (NHAGV.NE.0) THEN
               print *, "Post-Processing Unit 54"
               CALL POST54()
            ENDIF
         endif
C     
         if ((UNIT == 100) .or. (UNIT == 55)) then
            print *, "Post-Processing Unit 55"
            CALL POST55()
         endif
C     
         if ((UNIT == 100) .or. (UNIT == 61)) then
            IF (NOUTE.NE.0) THEN
               print *, "Post-Processing Unit 61"
               CALL POST61()
            ENDIF
         endif
C     
         if ((UNIT == 100) .or. (UNIT == 62)) then
            IF (NOUTV.NE.0) THEN
               print *, "Post-Processing Unit 62"
               CALL POST62()
            ENDIF
         endif
C     
         if ((UNIT == 100) .or. (UNIT == 63)) then
            IF (NOUTGE.NE.0) THEN
               print *, "Post-Processing Unit 63"
               CALL POST63()
            ENDIF
         endif
C     
         if ((UNIT == 100) .or. (UNIT == 64)) then
            IF (NOUTGV.NE.0) THEN
               print *, "Post-Processing Unit 64"
               CALL POST64()
            ENDIF
         endif
C     jgf45.07 wrote code to handle fort.71, fort.72, fort.73
         if ((UNIT == 100) .or. (UNIT == 71)) then
            IF (NOUTM.NE.0) THEN
               print *, "Post-Processing Unit 71"
               CALL POST71()
            ENDIF
         endif
C
         if ((UNIT == 100) .or. (UNIT == 72)) then
            IF (NOUTM.NE.0) THEN
               print *, "Post-Processing Unit 72"
               CALL POST72()
            ENDIF
         endif
C
         if ((UNIT == 100) .or. (UNIT == 73)) then
            IF (NOUTGW.NE.0) THEN
               print *, "Post-Processing Unit 73"
               CALL POST73()
            ENDIF
         endif
C     
         if ((UNIT == 100) .or. (UNIT == 74)) then
            IF (NOUTGW.NE.0) THEN
               print *, "Post-Processing Unit 74"
               CALL POST74()
            ENDIF
         endif
C     
C
!         if ((UNIT == 100) .or. (UNIT == 91)) then
         IF (UNIT == 91) THEN
            IF (NOUTM.NE.0) THEN
               print *, "Post-Processing Unit 91"
               CALL POST91()
            ENDIF
         endif
C
!         if ((UNIT == 100) .or. (UNIT == 93)) then
         IF (UNIT == 93) THEN
            IF (NOUTGW.NE.0) THEN
               print *, "Post-Processing Unit 93"
               CALL POST93()
            ENDIF
         endif

C ***POST PROCESSING FOR 3D FILES***
C
         if ((UNIT == 100) .or. (UNIT == 41)) then
            IF (I3DSD.NE.0) THEN
               print *, "Post-Processing Unit 41"
               CALL POST41()
            ENDIF
         endif
C     
         if ((UNIT == 100) .or. (UNIT == 42)) then
            IF (I3DSV.NE.0) THEN
               print *, "Post-Processing Unit 42"
               CALL POST42()
            ENDIF
         endif
C     
         if ((UNIT == 100) .or. (UNIT == 43)) then
            IF (I3DST.NE.0) THEN
               print *, "Post-Processing Unit 43"
               CALL POST43()
            ENDIF
         endif
C     
         if ((UNIT == 100) .or. (UNIT == 44)) then
            IF (I3DGD.NE.0) THEN
               print *, "Post-Processing Unit 44"
               CALL POST44()
               IF ((IDEN.EQ.3).OR.(IDEN.EQ.4)) THEN 
                  print *, "Post-Processing Unit 47"
                  CALL POST47()
               END IF 
            ENDIF
         endif
C     
         if ((UNIT == 100) .or. (UNIT == 45)) then
            IF (I3DGV.NE.0) THEN
               print *, "Post-Processing Unit 45"
               CALL POST45()
            ENDIF
         endif
C     
         if ((UNIT == 100) .or. (UNIT == 46)) then
            IF (I3DGT.NE.0) THEN
               print *, "Post-Processing Unit 46"
               CALL POST46()
            ENDIF
         endif
         
         if ((UNIT .NE. 100) .and. (UNIT .NE. 999)) goto 100
C     jgf45.06 END block of code written by MEB04/01/04 - - - - - - - - - 
         PRINT *,"Do you want to Post-Process hotstart files? ( Y/N )"
         READ(*,'(A1)') HOTANS
         IF((HOTANS.EQ.'Y').OR.(HOTANS.EQ.'y')) THEN
            print *, "Post-Processing hot start files"
            CALL POST67_68()
         ENDIF

C
         GO TO 99
C
C---------------------------------------------------------------------------
C  Parallel ADCIRC COMPARE Starts here
C---------------------------------------------------------------------------
C
      ELSEIF (PHASE(1:7).EQ.'compare') THEN
C
         print *, "enter pathame of first directory"
         READ(*,'(A80)') DIR1
         print *, "enter pathame of second directory"
         READ(*,'(A80)') DIR2
         print *, "enter error tolerance"
         READ(*,*) TOL
C     
        IF (NHASE.NE.0) THEN
           print *, "Comparing Unit 51 files"
           CALL COMPARE51(DIR1,DIR2,TOL)
        ENDIF
C     
        IF (NHASV.NE.0) THEN
           print *, "Comparing Unit 52 files"
           CALL COMPARE52(DIR1,DIR2,TOL)
        ENDIF
C     
        IF (NHAGE.NE.0) THEN
           print *, "Comparing Unit 53 files"
           CALL COMPARE53(DIR1,DIR2,TOL)
        ENDIF
C     
        IF (NHAGV.NE.0) THEN
           print *, "Comparing Unit 54 files"
           CALL COMPARE54(DIR1,DIR2,TOL)
        ENDIF
C     
        print *, "Comparing Unit 55 files"
        CALL COMPARE55(DIR1,DIR2,TOL)
C     
        IF (NOUTE.NE.0) THEN
           print *, "Comparing Unit 61 files"
           CALL COMPARE61(DIR1,DIR2,TOL)
        ENDIF
C     
        IF (NOUTV.NE.0) THEN
           print *, "Comparing Unit 62 files"
           CALL COMPARE62(DIR1,DIR2,TOL)
        ENDIF
C     
        IF (NOUTGE.NE.0) THEN
           print *, "Comparing Unit 63 files"
           CALL COMPARE63(DIR1,DIR2,TOL)
        ENDIF
C     
        IF (NOUTGV.NE.0) THEN
           print *, "Comparing Unit 64 files"
           CALL COMPARE64(DIR1,DIR2,TOL)
        ENDIF
C     
        IF (NOUTGW.NE.0) THEN
           print *, "Comparing Unit 74 files"
           CALL COMPARE74(DIR1,DIR2,TOL)
        ENDIF
C     
        GO TO 99
C
C
C---------------------------------------------------------------------------
C  Parallel ADCIRC DIFFMERGE Starts here
C---------------------------------------------------------------------------
C
      ELSEIF (PHASE(1:4).EQ.'diff') THEN
C     
         print *, "enter pathame of first directory"
         READ(*,'(A80)') DIR1
         print *, "enter pathame of second directory"
         READ(*,'(A80)') DIR2
         print *, "enter scale factor"    
         READ(*,*) SCALEVAL
C     
         IF (NOUTGE.NE.0) THEN
            print *, "DiffMerge of Unit 63 files"
            CALL DIFFMERGE63(DIR1,DIR2,SCALEVAL)
         ENDIF
C     
         IF (NOUTGV.NE.0) THEN
            print *, "DiffMerge of Unit 64 files"
            CALL DIFFMERGE64(DIR1,DIR2,SCALEVAL)
         ENDIF
C     
         GO TO 99
C     
      ELSEIF (PHASE(1:4).EQ.'quit') THEN
C     
         CALL EXIT(1)
C     
      ELSE
C     
         print *, "Not a Valid Operation ( Try Again )"
         GO TO 99
C     
      ENDIF
C     
      END
      
      
      SUBROUTINE GETMSG( STRING, MSG )
      INTEGER I, I1
      CHARACTER  STRING*(*),MSG*(*), TARGET
c
      I1 = 0
      TARGET = "!"
c
c--Find beginning of message
c
      DO I=1, 80
         IF (STRING(I:I).EQ.TARGET) THEN
            I1 = I
            GOTO 100
         ENDIF
      ENDDO
c     
 100  CONTINUE
c     
c--   Copy message to ouput string
c     
      DO I=1, I1-1
         MSG(I:I) = " "
      ENDDO
      MSG(I1:80)  = STRING(I1:80)
c     
      RETURN
      END


      SUBROUTINE NEWINDEX ( ISTRING, OSTRING, INDX )
      INTEGER I,I1,I2,I3,I4,INDX
      CHARACTER  ISTRING*(*),OSTRING*(*),TARGET
      CHARACTER TEMP1*80, TEMP2*100
c
      I1 = 0
      I2 = 0
      I3 = 0
      I4 = 0
      TARGET = " "
c
c--Find first non-blank character of String
c
      DO I=1, 80
         IF (ISTRING(I:I).NE.TARGET) THEN
            I1 = I
            GOTO 100
         ENDIF
      ENDDO
c     
c--   Find next blank character of String
c     
 100  CONTINUE
      DO I=I1+1,80
         IF (ISTRING(I:I).EQ.TARGET) THEN
            I2 = I
            GOTO 200
         ENDIF
      ENDDO
c     
 200  CONTINUE
c     
c--   Create a temporary string containing new index
c     
      WRITE(TEMP1(1:80),'(I8)') INDX
c     
c--   Find first non-blank character of String
c     
      DO I=1, 80
         IF (TEMP1(I:I).NE.TARGET) THEN
            I3 = I
            GOTO 300
         ENDIF
      ENDDO
c     
c--   Find next blank character of String
c     
 300  CONTINUE
      DO I=I3+1,80
         IF (TEMP1(I:I).EQ.TARGET) THEN
            I4 = I
            GOTO 400
         ENDIF
      ENDDO
c     
 400  CONTINUE
c     
Cdebug print *, "i1 i2 i3 i4 ",I1, I2, I3 , I4
      TEMP2(1:100) = TEMP1(I3-1:I4-1)//ISTRING(I2:80)
c
c--Write out first 80 characters of concatenated strings
c
      OSTRING(1:80) = TEMP2(1:80)
c
      RETURN
      END

      SUBROUTINE INSERT( ISTRING, OSTRING, NUMS, N )
      INTEGER I,J,I1,N,NUMS(N)
      CHARACTER  ISTRING*80,OSTRING*80,BLANK
      CHARACTER  TEMP1*80, TEMP2*160
c
      I1 = 0
      BLANK = " "
c
c--Create Tempoarary String TEMP1 containing NUMS
c
      WRITE(TEMP1(1:80),*) (NUMS(I),I=1,N)
c
c--Find length of TEMP1 string
c
      DO I=80,1,-1
         IF (TEMP1(I:I).NE.BLANK) THEN
           LEN1 = I
           GOTO 10
         ENDIF
      ENDDO
 10   CONTINUE  
c
c--Scan input string for character after old number list
c
      I = 1
      DO NUM=1, N
         DO J=I,80
            IF (ISTRING(J:J).NE.BLANK) THEN
              I = J
              GOTO 100
            ENDIF
         ENDDO
 100     CONTINUE
         DO J=I,80
           IF (ISTRING(J:J).EQ.BLANK) THEN
             I = J
             GOTO 200
           ENDIF
         ENDDO
 200     CONTINUE
      ENDDO
      I1 = MAX(0,I)
c
c--Insert Integer List into Message
c
      IF (I1.NE.0) THEN
C--if there is a message
        TEMP2(1:160) = TEMP1(1:LEN1+1)//ISTRING(I1:80)
      ELSE
        TEMP2(1:160) = TEMP1(1:LEN1+1)
      ENDIF
c
c--Write out first 80 characters of concatenated string
c
      OSTRING(1:80) = TEMP2(1:80)
c
      RETURN
      END
