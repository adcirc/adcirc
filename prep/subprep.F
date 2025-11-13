!-------------------------------------------------------------------------------!
!
! ADCIRC - The ADvanced CIRCulation model
! Copyright (C) 1994-2025 R.A. Luettich, Jr., J.J. Westerink
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
C
C                ADCIRC - SUBDOMAIN MODELING PREPROCESSING MODULE
C       
C    ========================================================================
C    |                                                                      |
C    |   This file contains the subroutines required by Subdomain Modeling, |
C    |   an approach to reduce the total runtime of a series of hurricane   |
C    |   storm surge simulations on a smaller grid, a subdomain grid,       |
C    |   extracted from the original grid.                                  |
C    |                                                                      |
C    |   Written by Alper Altuntas, aaltunt@ncsu.edu                        |
C    |   North Carolina State University,                                   |
C    |   2013                                                               |
C    |                                                                      |
C    ========================================================================

      module subprep
      use mod_io, only: openFileForRead
      use mod_logging, only : logMessage, allMessage,
     &  setMessageSource, unsetMessageSource, DEBUG, ECHO, INFO, WARNING, ERROR

      implicit none

C   NCSU Subdomain Modeling variables:
       logical :: subdomainOn
       integer :: enforceBN
       integer :: psbtiminc, pncbn, pnobn, pnibn
       real(8) :: pebn,pubn,pvbn
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
         character(1024) :: scratchMessage
         namelist /subdomainModeling/ subdomainOn

         FOUND_SM_NML = .false.

      call setMessageSource("readFort015prep")
#if defined(SUBPREP_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

         ioerror = 0
         !
         ! Look for the subdomainModeling namelist in the fort.15 file
         ! and if it is found, read the value of subdomainOn. 
         call openFileForRead(15,'fort.15',ioerror)
         if (ioerror.ne.0) then
            call allMessage(ERROR,'Failed to open fort.15 file.')
            CALL EXIT(1)
         endif
         !
         ! jgf51.42: Add a namelist for the user to control subdomain
         ! modeling. 
         READ(UNIT=15,NML=subdomainModeling,IOSTAT=ioerror)
         IF (ioerror.gt.0) THEN
            call logMessage(INFO,
     &         'The subdomainModeling namelist was not found.')
         else if(ioerror.eq.0) then
               FOUND_SM_NML = .true.
         endif
         if ((ioerror.gt.0).or.(subdomainOn.eqv..false.)) then
            call logMessage(INFO,
     &         'The subdomainModeling capability is not active.')     
         endif
         write(scratchMessage,*) "subdomainOn=",subdomainOn
         call logMessage(ECHO,trim(scratchMessage))
         close(15)
                
         if (subdomainOn) then
            call openFileForRead(1015,'fort.015',ioerror)
            if (ioerror.eq.0) then
               print *, "Subdomain Active"
               open(1015, file='fort.015')
               read(1015,*) dummy
               read(1015,*) dummy
               read(1015,*) enforceBN
               close(1015)
            else
               call allMessage(ERROR,'The fort.015 file was not found.')
               CALL EXIT(1)
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
         integer sdu(nproc)
         logical success

         call openFort019prep()
         call subOpenPrep(1019,'subdomain boundary conditions ',
     &        1,nproc, SDU, Success)


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

        do while(.true.)

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
            
        enddo

1905   close(1019)
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
                !lnTemp = IMAP_NOD_GL2(itemp+1,gn)
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
         integer sdu(nproc)
         logical success

         call openFort020prep()
         call subOpenPrep(1020,'subdomain boundary conditions ',
     &        1,nproc, SDU, Success)


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

        do while(.true.)

            read(1020,*,end=1905) tstep
            do iproc=1,nproc
                write(sdu(iproc),*) tstep
            enddo

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

        enddo

1905   close(1020)
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
                !lnTemp = IMAP_NOD_GL2(itemp+1,gn)
                 pnibnp(procTemp) = pnibnp(procTemp) + 1
             enddo
         enddo

      END SUBROUTINE openFort021prep




      SUBROUTINE subPrep021()

         use PRE_GLOBAL, only : nproc, itotproc, imap_nod_gl2

         implicit none
         integer :: iproc,i,n,tstep,gnode,lnode
         integer :: totalProc,p,itemp,proctemp,gn,lntemp
         integer sdu(nproc)
         logical success

         call openFort021prep()
         call subOpenPrep(1021,'subdomain boundary conditions ',
     &        1,nproc, SDU, Success)


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


        do while(.true.)

            read(1021,*,end=1905) tstep
            do iproc=1,nproc
                write(sdu(iproc),*) tstep
            enddo

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

        enddo

1905   close(1021)
       do iproc=1,nproc
           close(sdu(iproc))
       enddo

      END SUBROUTINE subPrep021








C     This subroutine is copied from prep.F (OpenPrepFiles) and is modified to
C     open subdomain modeling b.c. files

      SUBROUTINE subOpenPrep(UnitNumber, Description,
     &     startProc, endProc, SDU, Success)
C     NCSU SUBDOMAIN: modified to handle "fort.019"
C---------------------------------------------------------------------------
      USE PRE_GLOBAL, only : nproc
      IMPLICIT NONE
      INTEGER, intent(in) :: UnitNumber     ! i/o unit number to open
      CHARACTER(len=30), intent(in) :: Description ! description of file
      INTEGER, intent(in) :: startProc        ! subdomains to start with
      INTEGER, intent(in) :: endProc          ! subdomain to end on
      INTEGER, intent(out), dimension(nproc) :: SDU ! Subdomain unit numbers
      LOGICAL, intent(out):: Success     ! .true. if files opened w/o errors

      LOGICAL Found               !.true. if the full domain file exists
      CHARACTER(len=80) FileName  ! name of full domain file
      CHARACTER(len=8) DefaultName! default name of full domain file  !NCSU SUBDOMAIN
      INTEGER ErrorIO             ! zero if file opened successfully
      INTEGER iproc               ! subdomain index
      INTEGER UnitNumberIn
      CHARACTER(len=15) sdFileName     ! subdomain file name  !NCSU SUBDOMAIN

      Found = .false.
      Success = .false.
      ErrorIO = 1
      UnitNumberIn = UnitNumber

      ! NCSU Subdomain: 
      if (UnitNumberIn.eq.1019) then
          DefaultName= 'fort.019'    
          FileName = 'fort.019'      
      else if (UnitNumberIn.eq.1020) then
          DefaultName= 'fort.020'    
          FileName = 'fort.020' 
      else if (UnitNumberIn.eq.1021) then
          DefaultName= 'fort.021'
          FileName = 'fort.021'
      endif

C
C     Determine if full domain file exists
      INQUIRE(FILE=FileName,EXIST=FOUND)
C
C     If it does exist, open it
      IF ( FOUND ) THEN
         WRITE(*,1011) FileName !found
         OPEN(UNIT=UnitNumberIn, FILE=FileName, IOSTAT=ErrorIO)
         Success = .true.
         IF ( ErrorIO .GT. 0 ) THEN
            WRITE(*,*) "ERROR: Full domain file exists but"
            WRITE(*,*) "cannot be opened."
            Success = .false.
         ENDIF
      ENDIF
C
      If (.not.Success) RETURN ! failed to open full domain file
C
      DO iproc = startProc, endProc
         sdu(iproc) = 105 + (iproc-1)
         sdFileName = 'PE0000/'//FileName
         WRITE(sdFileName(3:6), '(I4.4)') iproc-1
         OPEN (UNIT=SDU(iproc), FILE=sdFileName, IOSTAT=ErrorIO)
         Success = .true.
         IF ( ErrorIO .GT. 0 ) THEN
            WRITE(*,*) "ERROR: Subdomain file cannot be opened."
            Success = .false.
            RETURN ! failed to open at least one subdomain file
         ENDIF
      ENDDO

 2    FORMAT(I2)
 30   FORMAT(A30)
 1010 FORMAT('File ',A8,/,' WAS NOT FOUND! Try again or type "skip"',/)
 1011 FORMAT('File ',A8,/,' WAS FOUND!  Opening & Processing file.',/)
      RETURN
C---------------------------------------------------------------------------
      END SUBROUTINE subOpenPrep
C---------------------------------------------------------------------------

      end module subprep








































