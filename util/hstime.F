!-------------------------------------------------------------------------------!
!
! ADCIRC - The ADvanced CIRCulation model
! Copyright (C) 1994-2023 R.A. Luettich, Jr., J.J. Westerink
! 
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU Lesser General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!-------------------------------------------------------------------------------!
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
      INTEGER :: ncid     ! id of netcdf file
      INTEGER :: dimid    ! id of time dimension
      INTEGER :: timesLen ! size of times array
      REAL(8), ALLOCATABLE :: times(:) ! array of times in hs file
      INTEGER :: varid    ! id of times variable
#endif
      INTEGER :: i
      INTEGER :: ihotstp  ! record counter for binary hotstart file
      INTEGER :: iArgC    ! function returns the number of command line options 
      INTEGER :: argCount ! number of command line options
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
#ifndef ADCNETCDF
            WRITE(0,*) "ERROR: hstime: NetCDF hotstart format was " //
     &     " specified but hstime was not compiled with NetCDF support." //
     &     " You must recompile hstime with NetCDF support." 
            write(6,*) 'null'
            STOP
#endif
         CASE("-f")
            i = i + 1
            CALL getArg(i,fileName)
         CASE DEFAULT
            write(0,*) "ERROR: hstime: command line option '",
     &         cmdLineArg(1:2),"' was not recognized."
            write(6,*) 'null'
            stop
         END SELECT
      END DO
C
      INQUIRE(FILE=TRIM(fileName),EXIST=fileFound)
      IF (fileFound.eqv..false.) THEN
         WRITE(0,*) "ERROR: hstime: the hotstart file '",
     &      TRIM(fileName),"' was not found."
         write(6,*) 'null'
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
            write(6,*) 'null' ! return nothing to stdout 
            close(99)
            ! write an error message to stderr for analyst
            write(0,*) 'ERROR: hstime: Hotstart file contains no data.'
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
         OPEN(99,FILE=TRIM(fileName),ACCESS='DIRECT',RECL=8,
     &      ACTION='READ',IOSTAT=errorIO,STATUS='OLD')
         IF (errorIO.gt.0) THEN
            WRITE(0,*) "ERROR: hstime: I/O error when opening existing
     &         hotstart file '",TRIM(fileName),"'."
            write(6,*) 'null'
            STOP
         ENDIF
         IHOTSTP=3
         READ(99,REC=IHOTSTP) TIMEHSF
         WRITE(*,*) TIMEHSF
         CLOSE(99)           
      ENDIF
      END PROGRAM hstime
C
#ifdef ADCNETCDF
      subroutine check_err(iret)
      IMPLICIT NONE
      include 'netcdf.inc'
      INTEGER, intent(in) :: iret
      if (iret .ne. NF_NOERR) then
         WRITE(0,*) "ERROR: hstime: ",nf_strerror(iret)
         write(6,*) 'null'
         stop
      endif
      end subroutine check_err
#endif
