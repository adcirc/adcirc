C***************************************************************************
C  PROGRAM TO READ A BINARY HOTSTART FILE AND PRINT TIME TO STDOUT         *
C                                                                          *
C     Written by Jason G. Fleming 12 Oct 2006                                 *
C     Updated for greater flexibility and support of netcdf hotstart
C        files jgf20110513     
C***************************************************************************

      PROGRAM hstime
      IMPLICIT NONE
#ifdef ADCNETCDF
      include 'netcdf.inc'
      INTEGER :: iret     ! return status from netcdf calls
!$omp threadprivate(iret)
      INTEGER :: ncid     ! id of netcdf file
!$omp threadprivate(ncid)
      INTEGER :: dimid    ! id of time dimension
!$omp threadprivate(dimid)
      INTEGER :: timesLen ! size of times array
!$omp threadprivate(timesLen)
      REAL(8), ALLOCATABLE :: times(:) ! array of times in hs file
!$omp threadprivate(times)
      INTEGER :: varid    ! id of times variable
!$omp threadprivate(varid)
#endif
      INTEGER :: i
!$omp threadprivate(i)
      INTEGER :: ihotstp  ! record counter for binary hotstart file
!$omp threadprivate(ihotstp)
      INTEGER :: iArgC    ! function returns the number of command line options 
!$omp threadprivate(iArgC)
      INTEGER :: argCount ! number of command line options
!$omp threadprivate(argCount)
      CHARACTER(2048) :: cmdLineArg  ! the command line arguments
!$omp threadprivate(cmdLineArg)
      LOGICAL :: fileFound           ! .true. if the file was found
!$omp threadprivate(fileFound)
      INTEGER :: ErrorIO             ! nonzero if there was an error
!$omp threadprivate(ErrorIO)
      CHARACTER(2048) :: fileName    ! full path name of hotstart file
!$omp threadprivate(fileName)
      LOGICAL :: netCDFFormat        ! .true. if hs file is in netcdf format
!$omp threadprivate(netCDFFormat)
      REAL(8) timeHSF
!$omp threadprivate(timeHSF)
C
      netCDFFormat = .false.
      fileFound = .false.
C
      i=0
      argCount = iArgC()
      DO WHILE (i < argCount) 
         i = i + 1
         CALL getArg(i,cmdLineArg)
         SELECT CASE(cmdLineArg(1:2))
         CASE("-n") 
            netCDFFormat = .true.
#ifndef ADCNETCDF
            WRITE(0+CK_LUN,*) "ERROR: hstime: NetCDF hotstart format was " //
     &     " specified but hstime was not compiled with NetCDF support." //
     &     " You must recompile hstime with NetCDF support." 
            write(6+CK_LUN,*) 'null'
            STOP
#endif
         CASE("-f")
            i = i + 1
            CALL getArg(i,fileName)
         CASE DEFAULT
            write(0+CK_LUN,*) "ERROR: hstime: command line option '",
     &         cmdLineArg(1:2),"' was not recognized."
            write(6+CK_LUN,*) 'null'
            stop
         END SELECT
      END DO
C
      INQUIRE(FILE=TRIM(fileName),EXIST=fileFound)
      IF (fileFound.eqv..false.) THEN
         WRITE(0+CK_LUN,*) "ERROR: hstime: the hotstart file '",
     &      TRIM(fileName),"' was not found."
         write(6+CK_LUN,*) 'null'
         STOP
      ENDIF         
C
      IF (netCDFFormat.eqv..true.) THEN
#ifdef ADCNETCDF
         iret = nf_open(TRIM(fileName),NF_NOWRITE, ncid)
         call check_err(iret)
         iret = nf_inq_unlimdim(ncid,dimid)
         call check_err(iret)
         iret = nf_inq_dimlen(ncid,dimid,timesLen)
         call check_err(iret)        
         !jgf51.41: Check to make sure the file actually has data.
         if (timesLen.eq.0) then
            write(6+CK_LUN,*) 'null' ! return nothing to stdout 
            close(99+CK_LUN)
            ! write an error message to stderr for analyst
            write(0+CK_LUN,*) 'ERROR: hstime: Hotstart file contains no data.'
            stop
         endif
         allocate(times(timesLen))
         iret = nf_inq_varid(ncid,"time",varid)
         call check_err(iret)    
         iret = nf_get_vara_double(ncid,varid,1,timesLen,times)        
         call check_err(iret)    
         write(*,*) times(timesLen) ! write the last if theres more than 1
         iret = nf_close(ncid)
         call check_err(iret)        
#endif
      ELSE
         OPEN(99+CK_LUN,FILE=TRIM(fileName),ACCESS='DIRECT',RECL=8,
     &      ACTION='READ',IOSTAT=errorIO,STATUS='OLD')
         IF (errorIO.gt.0) THEN
            WRITE(0+CK_LUN,*) "ERROR: hstime: I/O error when opening existing
     &         hotstart file '",TRIM(fileName),"'."
            write(6+CK_LUN,*) 'null'
            STOP
         ENDIF
         IHOTSTP=3
         READ(99+CK_LUN,REC=IHOTSTP) TIMEHSF
         WRITE(*,*) TIMEHSF
         CLOSE(99+CK_LUN)           
      ENDIF
      END PROGRAM hstime
C
#ifdef ADCNETCDF
      subroutine check_err(iret)
      IMPLICIT NONE
      include 'netcdf.inc'
      INTEGER, intent(in) :: iret
      if (iret .ne. NF_NOERR) then
         WRITE(0+CK_LUN,*) "ERROR: hstime: ",nf_strerror(iret)
         write(6+CK_LUN,*) 'null'
         stop
      endif
      end subroutine check_err
#endif
