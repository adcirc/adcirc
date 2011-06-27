C***************************************************************************
C  PROGRAM TO READ A BINARY HOTSTART FILE AND PRINT TIME TO STDOUT         *
C                                                                          *
C     Written by Jason G. Fleming 12 Oct 2006                                 *
C     Updated for greater flexibility and support of netcdf hotstart
C        files jgf20110513     
C***************************************************************************

      PROGRAM hstime
      IMPLICIT NONE
#ifdef NETCDF
      include 'netcdf.inc'
      INTEGER :: iret     ! return status from netcdf calls
      INTEGER :: ncid     ! id of netcdf file
      INTEGER :: dimid    ! id of time dimension
      INTEGER :: timesLen ! size of times array
      REAL(8), ALLOCATABLE :: times(:) ! array of times in hs file
      INTEGER :: varid    ! id of times variable
#endif
      INTEGER :: i
      INTEGER :: ihotstp  ! record counter for binary hotstart file
      INTEGER :: iArgC    ! function returns the number of command line options 
      INTEGER: : argCount ! number of command line options
      CHARACTER(2048) :: cmdLineArg  ! the command line arguments
      LOGICAL :: fileFound           ! .true. if the file was found
      INTEGER :: ErrorIO             ! nonzero if there was an error
      CHARACTER(2048) :: fileName    ! full path name of hotstart file
      LOGICAL :: netCDFFormat        ! .true. if hs file is in netcdf format
      REAL(8) timeHSF
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
#ifndef NETCDF
            WRITE(*,*) "ERROR: hstime: NetCDF hotstart format was ",
     &     " specified but hstime was not compiled with NetCDF support."
     &     " You must recompile hstime with NetCDF support." 
            STOP
#endif
         CASE("-f")
            i = i + 1
            CALL getArg(i,fileName)
         CASE DEFAULT
            write (*,*) "ERROR: hstime: command line option '",
     &         cmdLineArg(1:2),"' was not recognized."
         END SELECT
      END DO
C
      INQUIRE(FILE=TRIM(fileName),EXIST=fileFound)
      IF (fileFound.eqv..false.) THEN
         WRITE(*,*) "ERROR: hstime: the hotstart file '",
     &      TRIM(fileName),"' was not found."
         STOP
      ENDIF         
C
      IF (netCDFFormat.eqv..true.) THEN
#ifdef NETCDF
         iret = nf_open(TRIM(fileName),NF_NOWRITE, ncid)
         call check_err(iret)
         iret = nf_inq_unlimdim(ncid,dimid)
         call check_err(iret)
         iret = nf_inq_dimlen(ncid,dimid,timesLen)
         call check_err(iret)        
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
         OPEN(99,FILE=TRIM(fileName),ACCESS='DIRECT',RECL=8,
     &      ACTION='READ',IOSTAT=errorIO,STATUS='OLD')
         IF (errorIO.gt.0) THEN
            WRITE(*,*) "ERROR: hstime: I/O error when opening existing
     &         hotstart file '",TRIM(fileName),"'."
            STOP
         ENDIF
         IHOTSTP=3
         READ(99,REC=IHOTSTP) TIMEHSF
         WRITE(*,*) TIMEHSF
         CLOSE(99)           
      ENDIF
      END PROGRAM hstime
C
#ifdef NETCDF
      subroutine check_err(iret)
      IMPLICIT NONE
      include 'netcdf.inc'
      INTEGER, intent(in) :: iret
      if (iret .ne. NF_NOERR) then
         WRITE(*,*) "ERROR: hstime: ",nf_strerror(iret)
         stop
      endif
      end subroutine check_err
#endif
