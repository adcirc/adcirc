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
      module soln
         implicit none
         integer, parameter :: PATH = 2048, VS=30
         type Soln_t
            integer          :: maxsize
            integer          :: nsize
            integer          :: rectype
            integer          :: it
            integer          :: itype     ! is 1 when ROW, 2 for natural order
            real(8)          :: time
            real(8), pointer :: a(:), b(:)
            integer, pointer :: nc(:)
         end type Soln_t



         type SolnFile_t
            character(PATH) :: filename
            character(VS)   :: varname(2)
            integer         :: rectype
            integer         :: lun
            logical         :: etaType
            real(8)         :: rel_error(2)
            real(8)         :: abs_error(2)
            real(8)         :: max_error(2)
         end type SolnFile_t

         integer, pointer   :: ncList(:)
         integer,save       :: icount = 0
         integer,save       :: isz    = 0


         contains

         subroutine initialize_soln(x)
         type (Soln_t) :: x 
         x % maxsize = 0

         end subroutine 

         subroutine build_solnfile_entry(sf, filename,varname,
     &     rectype, etaType, lun)
         implicit none

         character(*)     :: filename
         character(VS)    :: varname(2)
         integer          :: rectype, lun
         logical          :: etaType
         type(SolnFile_t) :: sf

         sf % filename     = filename
         sf % varname(:)   = varname(:)
         sf % rectype      = rectype
         sf % etaType      = etaType
         sf % lun          = lun
         sf % rel_error(:) = -1.0
         sf % abs_error(:) = -1.0
         sf % max_error(:) = -1.0

         end subroutine

         subroutine fill_error_values(nfiles, fileA, varname,
     &                                str_rel_err, str_abs_err)

         implicit none

         integer          :: nfiles
         type(SolnFile_t) :: fileA(nfiles)
         character(VS)    :: varname, str_rel_err, str_abs_err

         integer          :: i,j
         real(8)          :: rel_err, abs_err

         read(str_rel_err,'(E20.10)') rel_err
         read(str_abs_err,'(E20.10)') abs_err

         do i = 1, nfiles
           do j = 1, fileA(i) % rectype
             if (varname == fileA(i) % varname(j)) then
               fileA(i) % rel_error(j) = rel_err
               fileA(i) % abs_error(j) = abs_err
             endif
           enddo
         enddo
         end subroutine

      end module soln

!-------------------------------------------------------------------------------
!  Main program
!-------------------------------------------------------------------------------

      program adccmp

      use soln
      implicit none

      type(SolnFile_t), allocatable :: fileA(:)
      type(Soln_t) :: trial, gold
      integer      :: count, iargc, nvar, ivar, lunGold, lunTrial
      integer      :: iarg , nstepsGold, nstepsTrial, istep
      integer      :: nfiles, i, k
      logical      :: foundGold, foundTrial, compare, results, passed
      logical      :: eof_gold, eof_trial, readarray1, readarray2
      logical      :: etaType
      character(PATH) :: gold_dir, trial_dir, fnameGold, fnameTrial
      character(VS)   :: varname, str_rel_err, str_abs_err

      character(80) :: usage, flag

      common /report_cmn/ reportFile
      logical reportFile


      usage = "Usage: adccmp gold_dir trial_dir varname1 "//
     &    "rel_err abs_err ..."

      count = COMMAND_ARGUMENT_COUNT()
      if (count < 5) then
        print *, usage
        stop
      endif

      print *,"ADCCMP Version 1.1"

      iarg = 1

      reportFile = .FALSE.
      call GET_COMMAND_ARGUMENT(iarg, flag)
      if (flag(1:2) == "-r") then
        reportFile = .TRUE.
        iarg = iarg + 1
      endif

      ! save away list of nodes where the nodecodes differ.
      isz = 10
      allocate(ncList(isz))

      nfiles = 9
      allocate(fileA(nfiles))

      call build_fileA(nfiles, fileA)

      call GET_COMMAND_ARGUMENT(iarg,gold_dir);  iarg = iarg + 1
      call GET_COMMAND_ARGUMENT(iarg,trial_dir); iarg = iarg + 1

      if (mod((count - iarg + 1),3) /= 0) then
        print *, usage
        stop
      endif

      nvar = (count - iarg + 1)/3
      do ivar = 1,nvar
        call GET_COMMAND_ARGUMENT(iarg, varname);     iarg = iarg + 1
        call GET_COMMAND_ARGUMENT(iarg, str_rel_err); iarg = iarg + 1
        call GET_COMMAND_ARGUMENT(iarg, str_abs_err); iarg = iarg + 1

        call fill_error_values(nfiles, fileA, varname,
     &    str_rel_err, str_abs_err)
      end do

      passed = .TRUE.

      print *, " VarName   Relative_Error    Absolute_Error "

      do i = 1, nfiles

        do ivar = 1, 2
          if (fileA(i) % rel_error(ivar) < 0) cycle
          print 2000, trim(fileA(i) % varname(ivar)),
     &      fileA(i) % rel_error(ivar),
     &      fileA(i) % abs_error(ivar)
        enddo
      enddo

      call initialize_soln(gold)
      call initialize_soln(trial)

      do i = 1, nfiles

        etaType = fileA(i) % etaType
        if (fileA(i) % rel_error(1) < 0) cycle

        fnameGold  = trim(gold_dir) //'/'// fileA(i) % filename
        fnameTrial = trim(trial_dir) //'/'// fileA(i) % filename

        inquire (file=fnameGold,  exist = foundGold)
        inquire (file=fnameTrial, exist = foundTrial)

        if (.not. foundGold) cycle

        if (foundGold .neqv. foundTrial) then
          call report("failed", "File '"//trim(fnameTrial)
     &      //"' does not exist", .TRUE.)
        endif

        lunGold  = fileA(i) % lun + 200
        lunTrial = fileA(i) % lun
        call readheader(lunGold,  fnameGold,  gold)
        call readheader(lunTrial, fnameTrial, trial)

        print 1000, trim(fileA(i) % filename)

        do   ! loop over time steps
          if (fileA(i) % rectype == 1) then
            eof_gold  = readarray1(lunGold,  etaType, gold)
            eof_trial = readarray1(lunTrial, etaType, trial)
          else
            eof_gold  = readarray2(lunGold,  gold)
            eof_trial = readarray2(lunTrial, trial)
          endif

          if (eof_gold .and. eof_trial) exit

          if (eof_gold .or. eof_trial) then
            call report("failed", "Number of time steps do not match",
     $      .TRUE.)
          endif

          results =  compare(fileA(i), gold, trial)

          if (.not. results) then
            passed = .FALSE.
            if (icount > 0) then
              print *,""
              print *,"Number of nodes where nodecode differ: ",icount
              print 1010, (nclist(k),k= 1,icount)
            endif
          endif
          


        enddo  ! timesteps
        close(lunGold)
        close(lunTrial)

      end do   ! files

      if (passed) then
        call report("passed","",.FALSE.)
      else
        call report("diff","", .FALSE.)
      endif


 1000 format(/,1x,'!-----------------------------------------',
     &       /,1x,'#  ',a,
     &       /,1x,'!-----------------------------------------')
 1010 format(1x, 10i6)
 2000 format(1x,a5,5x, 1pg12.5, 6x,1pg12.5)

      end program

!-------------------------------------------------------------------------------
!  Read Header
!-------------------------------------------------------------------------------

      subroutine readheader(lun, fname, x)
      use soln
      implicit none

      type(Soln_t)  :: x
      character(*)  :: fname
      character(80) :: title
      integer       :: nsteps, nsize, rectype
      real (8)      :: dt
      integer       :: lun, idummy

      call check_format(lun, fname, x % itype)

      open (lun, file=fname)


      read(lun, '(a)') title
      read(lun, *) nsteps, nsize, dt, idummy, rectype

      if (nsize < 1) then
        call report("failed", "Array size is zero in: "//fname,
     $      .TRUE.)
      endif
        
      call resize(nsize, rectype, x)

      end subroutine

!-------------------------------------------------------------------------------
!  Read array 1
!-------------------------------------------------------------------------------

      function readarray1(lun, etaType, x) result(eof)
      use soln
      implicit none

      integer :: lun, ndomains, idom, kdom, nitems, i
      integer :: iglobal
      logical :: eof, etaType
      real(8) :: x1, x2
      type (Soln_t) :: x

      eof = .FALSE.

      if (x % itype == 1) then
        read(lun, 1100, end = 200) x % time, x % it, nitems
        x % a(1: x % nsize) = 0.0
        if (etatype) x % nc(1:x % nsize) = 0       ! vjp
        do i = 1, nitems
          read(lun, *, err = 100, end = 150) iglobal, x1
          x % a(iglobal)  = x1
          if (etatype) x % nc(iglobal) = 1         ! vjp
        end do
        return

      else if ( x % itype == 2) then
        read(lun, 1100, end = 200) x % time, x % it

        x % nc(1:x % nsize) = 1
        do i = 1, x % nsize
          read(lun, *, end = 200) iglobal, x1
          if (etaType .and. x1 == -99999.0) then
             x1 = 0.0
             x % nc(iglobal) = 0
          endif
          x % a(iglobal) = x1
        end do

      else
        eof = .TRUE.
      endif

      return
  100 continue
      backspace lun
  150 continue 
      return
  200 continue
      eof = .TRUE.

 1100 FORMAT(2X,E20.10,5X,I10, 5x, I10)
 1200 format(9x,i8,12x, i8)
 1300 format(2x, i8, 2x, E20.10)
      end function
      
!-------------------------------------------------------------------------------
!  Read array 2
!-------------------------------------------------------------------------------

      function readarray2(lun, x) result(eof)
      use soln
      implicit none

      integer :: lun, ndomains, idom, kdom, nitems, i
      integer :: iglobal
      logical :: eof
      real(8) :: x1, x2
      type (Soln_t) :: x

      eof = .FALSE.

      if (x % itype == 1) then
        read(lun, 1100, end = 200) x % time, x % it, nitems
        x % a(1: x % nsize) = 0.0
        x % b(1: x % nsize) = 0.0
        do i = 1, nitems    
          read(lun, *, end = 150, err = 100) iglobal, x1, x2
          x % a(iglobal)  = x1
          x % b(iglobal)  = x2
        end do

      else if ( x % itype == 2) then
        read(lun, 1100, end = 200) x % time, x % it
        do i = 1, x % nsize
          read(lun, *, end = 200) iglobal, x1, x2
          x % a(iglobal) = x1
          x % b(iglobal) = x2
        end do

      else
        eof = .TRUE.
      endif

      return
  100 continue
      backspace lun
  150 continue
      return 
  200 continue
      eof = .TRUE.

 1100 FORMAT(2X,E20.10,5X,I10, 5x, I10)
 1200 format(9x,i8,12x, i8)
 1300 format(2x, i8, 2x, E20.10, E20.10)
      end function

!-----------------------------------------------------------------------
! Resize solution array to be big enough
!-----------------------------------------------------------------------
      subroutine resize(nsize, rectype, x)
      use soln
      implicit none
      type(Soln_t) :: x
      integer :: nsize, rectype

      if (x % maxsize == 0) then
        x % maxsize = max(1,nsize)
        allocate( x % a(nsize))
        allocate( x % b(nsize))
        allocate( x % nc(nsize))   
        print *, 'allocating nsize:', nsize
      else if (x % maxsize < nsize) then
        deallocate (x % a)
        deallocate (x % b)
        deallocate (x % nc) 
        allocate( x % a(nsize))
        allocate( x % b(nsize))
        allocate( x % nc(nsize))   
        x % maxsize = max(1,nsize)
        print *, 'reallocating nsize:', nsize
      end if

      x % nsize = nsize
      x % rectype = rectype

      end subroutine

      !-----------------------------------------------------------------------
      ! List of files to test
      !-----------------------------------------------------------------------

      subroutine build_fileA(nfiles, fileA)

      use soln
      implicit none

      integer          :: nfiles
      type(SolnFile_t) :: fileA(nfiles)
      character(VS)    :: varnames(2)
      logical          :: t,f
      t = .TRUE.
      f = .FALSE.

      varnames(1) = 'ET00'; varnames(2) = ""
      call build_solnfile_entry(fileA(1),'fort.61', varnames, 1, t, 61)

      varnames(1) = 'UU00'; varnames(2) = "VV00"
      call build_solnfile_entry(fileA(2),'fort.62', varnames, 2, f, 62)

      varnames(1) = 'ETA2'; varnames(2) = ""
      call build_solnfile_entry(fileA(3),'fort.63', varnames, 1, t, 63)

      varnames(1) = 'UU2'; varnames(2) = "VV2"
      call build_solnfile_entry(fileA(4),'fort.64', varnames, 2, f, 64)

      varnames(1) = 'RMP00'; varnames(2) = ""
      call build_solnfile_entry(fileA(5),'fort.71', varnames, 1, f, 71)

      varnames(1) = 'RMU00'; varnames(2) = "RMV00"
      call build_solnfile_entry(fileA(6),'fort.72', varnames, 2, f, 72)

      varnames(1) = 'PR2'; varnames(2) = ""
      call build_solnfile_entry(fileA(7),'fort.73', varnames, 1, f, 73)

      varnames(1) = 'WVNXOut'; varnames(2) = "WVNYOut"
      call build_solnfile_entry(fileA(8),'fort.74', varnames, 2, f, 74)

      varnames(1) = 'CC00'; varnames(2) = ""
      call build_solnfile_entry(fileA(9),'fort.81', varnames, 1, f, 81)
      end subroutine

!-------------------------------------------------------------------------------
!     Compare
!-------------------------------------------------------------------------------
      function compare(file, gold, trial) result(r)

      use soln
      implicit none

      type (Soln_t)     :: gold, trial
      type (SolnFile_t) :: file
      logical           :: r
      integer           :: i, ivar, ierror = 0, iszOld, nodDiff(2)
      logical           :: relFlag, absFlag, ncFlag

      real(8)           :: sumDiff(2), diff(2), varMax(2)
      real(8)           :: ref(2), norm(2), diffSln(2), relDiff
      character(30)     :: myResults(2)
      integer,pointer   :: itmp(:)

      r = .TRUE.
      if (gold % nsize /= trial % nsize) then
        call report("failed", " Array sizes do not match", .TRUE.)
      end if

      if (gold % rectype /= trial % rectype) then
        call report("failed", "incompable rectype", .TRUE.)
      end if

      print *, "time = ", gold % time, trial % time
      if (abs (gold % time - trial % time) > 1.e-6) then
        call report("failed", "time differs", .TRUE.)
      end if


      sumDiff(:) = 0.0
      diffSln(:) = 0.0
      varMax(:)  = 0.0
      norm(:)    = 0.0
      nodDiff(:) = 0
      icount     = 0   ! zero the list of nodecode problem children.


      ncFlag     = .FALSE.
      if (gold % rectype == 1) then
        do i = 1, gold % nsize
          ref(1)     = gold % a(i)
          varMax(1)  = max(varMax(1), abs(gold % a(i)))
          diff(1)    = abs(gold % a(i) - trial % a(i))

          norm(1)    = norm(1) + ref(1) * ref(1)
          sumDiff(1) = sumDiff(1) + diff(1)* diff(1)
          if (diff(1) > diffSln(1) ) then
            diffSln(1) = diff(1)
            nodDiff(1) = i
          endif

          if (file%etatype .and. gold % nc(i) /= trial % nc(i)) then   ! vjp
            ncFlag = .TRUE.
            icount = icount + 1
            if (icount > isz) then
              itmp => ncList
              iszOld = isz
              isz    = 2*isz
              allocate(ncList(isz))
              ncList(1:iszOld) = itmp(:)
              deallocate(itmp)
            end if
            ncList(icount) = i
          endif
        end do

      else
        do i = 1, gold % nsize
          ref(1)     = gold % a(i)
          ref(2)     = gold % b(i)
          varMax(1)  = max(varMax(1), abs(gold % a(i)))
          varMax(2)  = max(varMax(2), abs(gold % b(i)))
          diff(1)    = abs(gold % a(i) - trial % a(i))
          diff(2)    = abs(gold % b(i) - trial % b(i))

          norm(1)    = norm(1) + ref(1) * ref(1)
          norm(2)    = norm(2) + ref(2) * ref(2)
          sumDiff(1) = sumDiff(1) + diff(1)* diff(1)
          sumDiff(2) = sumDiff(2) + diff(2)* diff(2)
          if (diff(1) > diffSln(1) ) then
            diffSln(1) = diff(1)
            nodDiff(1) = i
          endif
          if (diff(2) > diffSln(2) ) then
            diffSln(2) = diff(2)
            nodDiff(2) = i
          endif
! vjp
!         if (gold % nc(i) /= trial % nc(i)) then
!           ncFlag = .TRUE.
!           icount = icount + 1
!           if (icount > isz) then
!             itmp => ncList
!             iszOld = isz
!             isz    = 2*isz
!             allocate(ncList(isz))
!             ncList(1:iszOld) = itmp(:)
!             deallocate(itmp)
!           end if
!           ncList(icount) = i
!         endif
! vjp
        end do

      endif

      print 1020,  gold % time, gold % it

      do ivar = 1, gold % rectype
        relFlag = .FALSE.
        absFlag = .FALSE.

        sumDiff(ivar) = sqrt(sumDiff(ivar))
        norm(ivar)    = sqrt(norm(ivar))

        relDiff = 0.0
        if (norm(ivar) > 1.0e-20) then
          relDiff = sumDiff(ivar)/norm(ivar)
          if (relDiff > file % rel_error(ivar)) then
            relFlag = .TRUE.
          endif
        endif

        if (sumDiff(ivar) > file % abs_error(ivar)) then
          absFlag = .TRUE.
        endif

        myResults(ivar) = "passed"
        if (relFlag .and. absFlag) then
          myResults(ivar) = "failed"
          ierror = ierror + 1
        endif

        if (ncFlag) then
          myResults(ivar) = "failed : nodecode"
          ierror = ierror + 1
        endif

        if (abs(varMax(ivar)) < 1.0e-20) varMax(ivar) = 1.0

        write(*,1000) trim(file % varname(ivar)), norm(ivar),
     &    sumDiff(ivar), relDiff, diffSln(ivar),
     &    diffSln(ivar)/varMax(ivar), nodDiff(ivar), varMax(ivar),
     &    trim(myResults(ivar))
      enddo


      r = ierror == 0
 1000 format(1x,a6,2x,5(1pg12.5,2x),i10,4x,1pg12.5,2x,a)
 1020 format(/,1x,"Time: ", 1pg12.5, 2x, "iteration: ", i10,
     &  /,1x,"Name      Norm         L2_diff      Rel_L2_diff",
     & "     L-inf       Rel_L-inf      MaxDiffNode     MaxVal  ",
     &  "   Result")
      end function

!-------------------------------------------------------------------------------
!     Report
!-------------------------------------------------------------------------------
      subroutine report(flag, msg, abort)
      implicit none

      logical       :: abort
      character(*)  :: flag, msg
      character(8)  :: day
      character(10) :: tm
      common /report_cmn/ reportFile
      logical reportFile

      if (.not. reportFile) then
        print *, trim(msg)
      else

        call date_and_time(day,tm)

        open(27,file="results.lua")
        write(27,1000) trim(day), trim(tm), trim(flag), trim(msg)

      endif
      if (abort) stop

 1000 format("-- -*- lua -*-",/,
     &    "-- ",a,1x,a,/,
     &    'myTbl.result = "',a,'"'/,
     &    'myTbl.reason = "',a,'"')
      end subroutine

!-------------------------------------------------------------------------------
!     Check format
!-------------------------------------------------------------------------------
      subroutine check_format(lun, fname, itype)
      implicit none

      character(*)  :: fname
      character(80) :: line
      integer       :: lun, itype, i(2), isum, icount
      real(8)       :: time

      open (lun, file=fname)

      itype = 0
      read(lun,'(A)',end=100) line
      read(lun,'(A)',end=100) line
      read(lun,'(A)',end=100) line

      i(:) = 0
      read(line,*,end=10) time,i
!
10    isum = icount(i)

      if (isum == 2) then
         itype = 1
      else
         itype = 2
      endif

      close(lun)
      print *, "file format = ", itype
      return
  100 continue
      close(lun)
      end subroutine check_format

      function icount(i) result (isum)
      integer i(2), isum, j
      isum = 0
      do j=1, 2
         if (i(j) .ne.0) isum = isum + 1
      enddo
      return
      end function icount


