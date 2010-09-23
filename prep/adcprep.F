C----------------------------------------------------------------------------
C
C                        ADCPREP Version 2.2 (  09/20/2006 )
C
C----------------------------------------------------------------------------
C
C  Program Development History
C  ---------------------------
C   Written for ADCIRC_v24.05    ( S. Chippada 1996 )
C   Updated for ADCIRC_v24.05    ( M. Martinez 1997 )
C   Hilbert Space Filling Curve  ( C. Edwards 1997 )
C   Updated for ADCIRC_v34.04    ( V. Parr 1998 )
C   METIS 4.0 Interface          ( V. Parr 1999 )
C   Added UNIT 12 for polders    ( V. Parr 1999 )
C   hotstarts for Windstress     ( V. Parr 1999 )
C   Updated for ADCIRC_v34.08    ( V. Parr 1999 )
C   Modified for PADCIRC90       ( V. Parr & J. Romo 1999 )
C   Parallel version for 2 PEs   ( V. Parr & J. Romo 2000 )
C   Add hotprep1 command         ( V. Parr & J. Romo 2001 )
C   Add hotprep2 command         ( R. Luettich, 10/2001 )
C   Expanded to include 3D       ( T. Campbell 2002)
C   Bug fix                      ( J. Fleming/R. Luettich 2/2002)
C   Added wave stress capability ( M. Brown 3/2003)
C   NWS=3,6 fixes                ( M. Brown 3/2003)
C   NWS=+-4 revised to not convert to NWS = 5  ( R. Luettich 3/2003)
C   Wave stress routines revise  ( R. Luettich 3/2003)
C   Time-varying Flux files      ( M. Brown 4/2004 )
C   Added default filename capability  ( M. Brown 4/2004 )
C   Islands kept as islands      ( R. Luettich 5/2004)
C   Vic Parr's bug fix for NWS=+-4, Wave stress  ( R. Luettich 5/2004)
C   Removed NWS +/- conversion   ( R. Luettich/jgf 11/2005 )
C   Add PREP14_15 option         ( J. Fleming 11/2005 )
C   3D Update: Recording Stations( J. Fleming 02/2006 )
C   Wrote PREP11,rewrote PREP20  ( J. Fleming 03/2006 )
C
C   Add partmesh command         ( V. Parr 09/2006 )
C   Add routine relocalize       ( V. Parr 09/2006 )
C   Add prep13 ONLY command      ( V. Parr 09/2006 )
C   Add prep15 ONLY command      ( V. Parr 10/2006 )
C   Add prephot ONLY command     ( V. Parr 10/2006 )
C   Added NWS=19 asymmetric holland v2.0 ( R. Weaver 07/2009)
c   Removed NWS=3 and NWS=6      (T.C. Massey 11/2009) 
C
C----------------------------------------------------------------------------
C
C  ADCPREP performs 7 operations:
C    (1) partmesh  =  partition the mesh with metis                        
C    (2) prepall   =  Localize grid, input, and hot start files using default
C                     file names
C    (3) prepspec  =  same as prepall but allows custom file names and 
C                     skipping over files.
C    (4) prep15    =  Localizes ONLY a RunInfo file ( fort.15)         
C                     after a prepall has been run
C    (5) prep13    =  Localizes ONLY a nodal attributes file ( fort.13)         
C                     after a prepall has been run
C    (6) hotlocal  =  Localize  ONLY a hotstart file  (fort.67 or fort.68)       
C    (7) hotglobal =  Globalize ONLY a hotstart files (fort.67 or fort.68)       
C
C  Pre-processor
C  -------------
C  The Pre-processor "prep" performs data decomposition of the ADCIRC grid 
C  and localizes global input files to subdomains of the decomposition for 
C  parallel runs of ADCIRC.
C   
C  Subdirectories "PExxxx" of the working directory ( where "xxxx" is a 
C  4-digit integer representing the MPI_rank of a Processing Element ) 
C  are created and "localized" input files are placed in the appropriate 
C  subdirectory with the same filename as the global file in the working 
C  directory.  A unit 14 and 15 file, called respectively, setup.14 and
C  setup.15  are written to the working directory to allow setup to
C  prepare the source code for a parallel run.
C
C  The global input files currently localized by ADCPREP are:
C
C   fort.13    Global Nodal Attributes File
C   fort.14    Global Mesh File
C   fort.15    Global Input File
C   fort.19    Global Time-Variant Elevation Boundary Conditions
C   fort.20    Global Time-Variant Normal Flow Boundary Conditions
C   fort.22    Global Wind Stress and Atmospheric Pressure
C   fort.23    Global Wave Stress Forcing
C   fort.67    Global Hot Start file
C   fort.68    Global Hot Start file
C   fort.150   METIS input information
C   fort.151   METIS output information
C 
C  When a parallel ADCIRC job is run, each PExxxx takes its input from its 
C  subdirectory PExxxx and writes local output files to the same directory.
C
C----------------------------------------------------------------------------
C
      PROGRAM ADCPREP
      USE PRE_GLOBAL 
      USE VERSION
      use memory_usage
C
      CHARACTER JOB*80 ! user input for menu selection

C  Define File Format Version

      FileFmtVersion = VERSION_NUMBER(FileFmtMajor, 
     &                                FileFmtMinor,
     &                                FileFmtRev)

C
C     Initialize constants for case ICS = 2 inputs.
C
      R  =  6378206.4D0
      DEG2RAD = 3.141592653589793D0/180.0D0
      RAD2DEG = 180.0D0/3.141592653589793D0
C
C     Initialize flags.
C
      USE_DEFAULT = .FALSE.      
      APERIODIC_FLOW_BC = .FALSE. 
      PARTITION  = .FALSE.
      PREP_ALL   = .FALSE.
      PREP_13    = .FALSE.
      PREP_15    = .FALSE.
      HOT_LOCAL   = .FALSE.
      HOT_GLOBAL  = .FALSE.
C
C-- Say Hello Gracie
C
      print *," *****************************************"
      print *," ADCPREP Fortran90 Version 2.3  10/18/2006"
      print *," Serial version of ADCIRC Pre-processor   "
      print *," *****************************************"
      print *, " "
C
C--Init Memory Usage
C
      call memory_init()
C
C-- Prompt for user input
C
      print *, "Input number of processors for parallel ADCIRC run:"
      READ(*,*) MNPROC
C
C-- Copy MNPROC to NPROC since they are to be the same
C
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
C
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
C---------------------------------------------------------------------------
      END PROGRAM ADCPREP
C---------------------------------------------------------------------------




C---------------------------------------------------------------------------
C               S U B R O U T I N E   P R E P I N P U T 
C---------------------------------------------------------------------------
C     Read in the full domain input files, perform domain decomposition,
C     and write new subdomain input files.
C     This version is compatible with ADCIRC version 46.45. vjp 10/8/2006.
C---------------------------------------------------------------------------
C
      SUBROUTINE PREPINPUT()
      USE PRE_GLOBAL
      IMPLICIT NONE

      INTEGER I           ! loop counter for nodes, elements
      INTEGER PE          ! loop counter for MPI processors
      CHARACTER CMD*6     ! string to hold shell command
      CHARACTER PENUM*6   ! string to hold directory name
      CHARACTER DIRCMD*72 ! string to hold complete shell command line

C     jgf48.17 If compiling for ADCSWAN, need to copy UnSWAN file
Casey 090304: Moved this call to later, after the PE directories have been made.
#ifdef ADCSWAN
C     CALL PREPUNSWAN()
#endif

 
C-- If ONLY localizing Nodal-Attributes File ( fort.13 )
      if (PREP_13) then
        write(*,'(/)')
        write(*,*) 'Re-Writing subdomain Nodal Attributes (unit 13)'
        write(*,*) 'file for each PE.'
        CALL RELOCALIZE()
        CALL PREP13()
        GO TO 9999
      endif

C-- If ONLY processing hotstart File(s) ( fort.67-68 )
      if (HOT_LOCAL .or. HOT_GLOBAL) then
        write(*,'(/)')
        write(*,*) 'Writing subdomain hotstart file for each subdomain'
        CALL RELOCALIZE()
        if (HOT_LOCAL)  CALL HOTLOCALIZE()
        if (HOT_GLOBAL) CALL HOTGLOBALIZE()
        GO TO 9999
      endif

      CALL MAKE_DIRNAME()
      CALL SIZEUP()
      CALL ALLOC_MAIN1()
C
C-- Read the Global Grid File ( Unit 14 )
C
      CALL READ14() 
      print *, 'Global Grid file read successfully.'

C-- If ONLY localizing Run Info File ( fort.15 )
      if (PREP_15) then
        CALL READ15()
        print *, 'Global Run Info file read successfully.'
        CALL RELOCALIZE()
        CALL PREP15()
        GO TO 9999
      endif

C
C-- Partition Nodes with METIS 4.0 Graph Partition Library
C
      IF (PARTITION) THEN
        CALL METIS()
        print *, " "
        print *, 'METIS has partitioned nodes'
        GO TO 9999
      ELSE
        print *, 'Opening file partmesh.txt'
        OPEN(990,FILE='partmesh.txt')
C        print *, 'made it this part and MNP = ', MNP
        DO I=1, MNP
           READ(990,*) PROC(I)
C           print *," PROC for ",I," = ", PROC(I)
        ENDDO
        CLOSE(990)           
        print *, 'closed file'
      ENDIF
C
C-- Read the Global Input File ( Unit 15 )
C
      print *, 'start reading in fort.15'
      CALL READ15()
      print *, 'Global Input file read successfully.'
C
C     jgf45.07 In the case of starting a parallel job from a full domain
C     hot start file, we must reverse the signs of NOUTM, NOUTV, etc.
C     This is because ADCIRC wouldn't be able to find preexisting
C     output files in the subdomains to append to.
C
      IF (PREP_ALL.AND.(IHOT.EQ.67.OR.IHOT.EQ.68)) THEN
         print *, 'INFORMATION: Subdomain output files will be started'
         print *, '             anew rather than appended.'
         IF (NOUTE .GT.0.AND.ABS(NOUTE) .NE.3) NOUTE =-NOUTE
         IF (NOUTV .GT.0.AND.ABS(NOUTV) .NE.3) NOUTV =-NOUTV
         IF (NOUTC .GT.0.) NOUTC=-NOUTC
         IF (NOUTM .GT.0.AND.ABS(NOUTM) .NE.3) NOUTM =-NOUTM
         IF (NOUTGE.GT.0.AND.ABS(NOUTGE).NE.3) NOUTGE=-NOUTGE
         IF (NOUTGV.GT.0.AND.ABS(NOUTGV).NE.3) NOUTGV=-NOUTGV
         IF (NOUTGC.GT.0.                    ) NOUTGC=-NOUTGC
         IF (NOUTGW.GT.0.AND.ABS(NOUTGW).NE.3) NOUTGW=-NOUTGW
!         WRITE(6,*)'==NOUT==',NOUTE,NOUTV,NOUTM,NOUTGE,NOUTGV
C   kmd48.33bc add information for 3D
         IF (C3DVS) THEN
            IF (I3DSD.GT.0) I3DSD=-I3DSD
            IF (I3DSV.GT.0) I3DSV=-I3DSV
            IF (I3DST.GT.0) I3DST=-I3DST
            IF (I3DGD.GT.0) I3DGD=-I3DGD
            IF (I3DGV.GT.0) I3DGV=-I3DGV
            IF (I3DGT.GT.0) I3DGT=-I3DGT
      ENDIF

      ENDIF
C
C-- Decompose the ADCIRC grid into MNPROC subdomains
C
      print *, " "
      print *, "Determine the parameters MNPP and MNEP"
      CALL DOMSIZE()
C
      print *, "Allocate arrays dimensioned by MNPP and MNEP"
      CALL ALLOC_MAIN2()

      print *, " "
      print *, "Decomposition of grid begins"
      CALL DECOMP()
      print *, "Decomposition successful"
C
C-- Create MNPROC sub-directories of the working directory
C
      DO PE=0, MNPROC-1
         PENUM  = 'PE0000'
         CALL IWRITE(PENUM,3,6,PE)
         CALL MAKEDIR(PENUM)
      ENDDO
C
C-- Write a Local Grid File ( fort.14 ) for each PE
C
      print *, "Writing Local UNIT 14 (Grid) File for each PE"
      CALL PREP14()
C
C-- Write a Local Input file ( fort.15 ) for each PE
C
      print *, "Writing Local UNIT 15 (Run Info) File for each PE"
      CALL PREP15()
C
C     jgf48.17 If compiling for ADCSWAN, need to copy UnSWAN file
Casey 090304: Moved this call from earlier in this file.
#ifdef ADCSWAN
      CALL PREPUNSWAN()
#endif
C
C-- Write Message-Passing File for each PE
C
      print *, "Writing Message-Passing Info Files for each PE"
      IF (PREP_ALL) CALL PREP18()
C
C     Write a subdomain initial concentration file ( fort.10 ) for each PE
C
      IF (PREP_ALL) THEN
         IF (C2D_PTrans.or.C3D_PTrans) THEN
            PRINT *, "Writing subdomain UNIT 10 file for each PE."
            CALL PREP10()
         ENDIF
      ENDIF
C
C     Write a subdomain initial density file ( fort.11 ) for each PE
C
C   kmd48.33bc changed the if statement to include IHOT parameter
      IF (PREP_ALL) THEN
         IF ((CBaroclinic).AND.(IHOT.EQ.0)) THEN
            PRINT *, "Writing subdomain UNIT 11 file for each PE."
            CALL PREP11()
         ENDIF
      ENDIF

C   kmd48.33bc added writing of subdomain initial condition file
C
C     Write a subdomain initial condtion file (fort.17) for each PE
C
      IF (PREP_ALL) THEN
         IF ((CBaroclinic).AND.(IHOT.EQ.17)) THEN
            PRINT *, "Writing subdomain UNIT 17 file for each PE."
            CALL HOTINITCOND()
         ENDIF
      ENDIF

C
C     Write a subdomain nodal attributes file (fort.13) for each PE
C
      if (NWP.gt.0) then
        write(*,'(/)')
        write(*,*) 'Writing subdomain Nodal Attributes (unit 13)'
        write(*,*) 'file for each PE.'
        CALL PREP13()
      endif

C
C-- If required write a Local fort.19 file for each PE
C
      IF (PREP_ALL) THEN
         IF ((NOPE.GT.0).AND.(NBFR.EQ.0)) THEN
            print *, "Ready to Write Local UNIT 19 File for each PE"
            CALL PREP19()
         ENDIF
      ENDIF
C
C-- Write a Local Flux File ( fort.20 ) for each PE if needed
C   -- EXIST_FLUX and other variables are set when the fort.14 is 
C   -- initially read in within the READ14() subroutine.
C
      IF (PREP_ALL) THEN
         if (EXIST_FLUX .ne. 0 .and. APERIODIC_FLOW_BC) then
            CALL PREP20()
         endif
      ENDIF
C
C-- If required write a Local Wind Stress file for each PE
C
      IF (PREP_ALL) THEN
C     sb46.28sb01 NWS=12 doesn't need fort.22 decomposition
C     tcm_v49.04 NWS=3 and NWS=6 no longer need fort.22 decomposition
!         IF (NWS.NE.0.AND.ABS(NWS).NE.12)  THEN
         IF ((NWS.NE.0).AND.(ABS(NWS).NE.3).AND.
     &      (ABS(NWS).NE.6).AND.(ABS(NWS).NE.12))  THEN
            print *, "Ready to Write Local UNIT 22 File for each PE"
C     jgf46.00 Added NWS=7
C     jgfdebug46.02 Added NWS=45
C     jgf46.02 Added NWS=8
C     jgf46.16 Merged: 
C     cf & cm added NWS=9: asymmetric hurricane wind model
C     rjw added NWS=9: asymmetric hurricane wind model v2.0
c     tcm_v49.04 NWS=3 and NWS=6 no longer need fort.22 decomposition
            IF((NWS.EQ.1).OR.(ABS(NWS).EQ.2).OR.(ABS(NWS).EQ.4).OR.    
     &           (ABS(NWS).EQ.5).OR.(NWS.EQ.7).OR.
     &           (ABS(NWS).EQ.45).OR.(ABS(NWS).EQ.8)
     &           .OR.(ABS(NWS).EQ.9) .OR.(ABS(NWS).EQ.19)) THEN
               CALL PREP22()
c           ELSEIF (NWS.EQ.10.OR.NWS.EQ.11) THEN
c              CALL PREP200()
            ENDIF
         ENDIF
      ENDIF
C
C-- If required write a Local Wave Stress file for each PE
C
      IF (PREP_ALL) THEN
         IF(NRS.EQ.1) CALL PREP23()
      ENDIF

C  kmd48.33bc add information for 3D baroclinic simulations
      IF (PREP_ALL) THEN
        IF (CBAROCLINIC) THEN
           IF ((RES_BC_FLAG.GE.1).OR.(RES_BC_FLAG.LE.-1)) THEN
             IF (BCFLAG_LNM.EQ.1) THEN
               print *, "Ready to Write Local UNIT 35 File for each PE"
               CALL PREP35()
             END IF
           END IF
           IF (RES_BC_FLAG.EQ.2) THEN
             print *, "Ready to Write Local UNIT 36 File for each PE"
             CALL PREP36()
           ELSE IF (RES_BC_FLAG.EQ.3) THEN
             print *, "Ready to Write Local UNIT 37 File for each PE"
             CALL PREP37()
             IF (BCFLAG_TEMP.NE.0) THEN
               print *, "Ready to Write Local UNIT 38 File for each PE"
               CALL PREP38()
             END IF
           ELSE IF (RES_BC_FLAG.EQ.4) THEN
             print *, "Ready to Write Local UNIT 36 File for each PE"
             CALL PREP36()
             print *, "Ready to Write Local UNIT 37 File for each PE"
             CALL PREP37()
             IF (BCFLAG_TEMP.NE.0) THEN
               print *, "Ready to Write Local UNIT 38 File for each PE"
               CALL PREP38()
             END IF
           END IF
        END IF
      END IF

C
C-- Save domain-decomposition information for post-processor
C     
CC -- Start Addition by CF  8/2007    
      IF (PREP_ALL) THEN
         print *, "Writing domain-decomposition file for post-processor"
         CALL PREP80()
      ENDIF
CC -- Finish Addition by CF      

      !jgf48.07 Initialize netCDF output files
      IF (useNetCDF) THEN
         write(*,*) "INFO: Initializing the netCDF output file(s)."
         ! call netCDF module
         CALL PREPNETCDF()
      ENDIF

9999  CONTINUE
      print *, ""
      print *, "Finished pre-processing input files"
      RETURN
      END SUBROUTINE
C---------------------------------------------------------------------------
C     End of subroutine prepinput
C---------------------------------------------------------------------------




      SUBROUTINE GETMSG( STRING, MSG )
      INTEGER I, I1
      CHARACTER  STRING*(*),MSG*(*), TARGET
c
      I1 = 0
      TARGET = "!"
c
c-- Find beginning of message
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
c--Copy message to ouput string
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
c-- Find first non-blank character of String
c
      DO I=1, 80
         IF (ISTRING(I:I).NE.TARGET) THEN
           I1 = I
           GOTO 100
         ENDIF
      ENDDO
c
c-- Find next blank character of String
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
c-- Create a temporary string containing new index
c
      WRITE(TEMP1(1:80),'(I8)') INDX
c
c-- Find first non-blank character of String
c
      DO I=1, 80
         IF (TEMP1(I:I).NE.TARGET) THEN
           I3 = I
           GOTO 300
         ENDIF
      ENDDO
c
c-- Find next blank character of String
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
c-- Write out first 80 characters of concatenated strings
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
c-- Create Tempoarary String TEMP1 containing NUMS
c
      IF (N.EQ.1) THEN
        WRITE(TEMP1(1:80),'(I8)') NUMS(1)
Casey 090304: Changed this section to allow N = 3.
      ELSEIF (N.EQ.2) THEN
        WRITE(TEMP1(1:80),'(2I8)') NUMS(1),NUMS(2)
      ELSE
        WRITE(TEMP1(1:80),'(3I8)') NUMS(1),NUMS(2),NUMS(3)
      ENDIF
c
c-- Find length of TEMP1 string
c
      DO I=80,1,-1
         IF (TEMP1(I:I).NE.BLANK) THEN
           LEN1 = I
           GOTO 10
         ENDIF
      ENDDO
 10   CONTINUE  
c
c-- Scan input string for character after old number list
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
c-- Insert Integer List into Message
c
      IF (I1.NE.0) THEN
C-- if there is a message
        TEMP2(1:160) = TEMP1(1:LEN1+1)//ISTRING(I1:80)
      ELSE
        TEMP2(1:160) = TEMP1(1:LEN1+1)
      ENDIF
c
c-- Write out first 80 characters of concatenated string
c
      OSTRING(1:80) = TEMP2(1:80)
c
      RETURN
      END
