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
!******************************************************************************
! PADCIRC VERSION 45.12 03/17/2006                                            *
!  last changes in this file VERSION 45.56                                    *
!                                                                             *
!******************************************************************************
!
module MOD_SIZES
   implicit none

   public

   !...SET MAX OF DIGITS OF PRECISION "NPREC" THE GRID CAN BE EXPECTED TO HAVE
   !...NOTE: IF THE GRID WAS BUILT ON A 32 BIT COMPUTER, IT SHOULD BE
   !   ACCURATE TO ABOUT 7 DIGITS.  THUS, IF THE NODAL SPACING REQUIRES MORE
   !   THAN 5 DIGITS OF PRECISION, THE MODEL RESULTS MAY NOT BE TRUSTWORTHY.

   integer, parameter :: NBYTE = 8
   integer, parameter :: NPREC = 7
   !
   integer :: MNE, MNP, MNEI, MNOPE, MNETA, MNBOU, MNVEL, &
              MNTIF, MNBFR, MNFFR, MNSTAE, MNSTAV, MNSTAC, MNSTAM, MNHARF
   integer :: MNPROC ! Number of compute processors
   integer :: MNWPROC ! Number of writer processors
   logical :: hotstartWritersRequested = .false. ! -H (no longer supported) was given
   logical :: ISPLIT ! .true. if splitting output for writer !st3 100708:
   integer :: MNALLPROC ! Number of all processors (= MNPROC + MNWPROC)
   !
   !     Dimension of vertical FE mesh (To interface 2D & 3D)
   !     INTEGER :: MNODES                                        !changed to MNFEN 08/16/2005
   integer :: MNFEN

   !     LOGICAL C2DDI,C3D,C3DDSS,C3DVS,CLUMP,CTIP,CHARMV         !moved to GLOBAL.F  04/09/2004
   !
   ! For Definition of Working Directory
   !
   integer :: MYPROC
   logical :: WRITE_LOCAL_FILES
   logical :: WRITE_LOCAL_HOT_START_FILES
   ! jgf49.44: Make it an option to read local binary hotstart files
   logical :: READ_LOCAL_HOT_START_FILES
   ! kmd51: Trying to make harmonic files write to local directory
   logical :: WRITE_LOCAL_HARM_FILES

   character(2048) :: ROOTDIR
   character(2048), target :: INPUTDIR ! directory containing input files
   character(2048), target :: GLOBALDIR ! directory for fulldomain files
   character(2048), target :: LOCALDIR ! directory for subdomain files
   character(2048), target :: HARMDIR ! directory where harmonic analysis file are written
   character(2048) :: GBLINPUTDIR
   !
   ! file formats
   integer, parameter :: numFormats = 7
   integer, parameter :: OFF = 0
   integer, parameter :: ASCII = 1
   integer, parameter :: NETCDF3 = 3
   integer, parameter :: SPARSE_ASCII = 4
   integer, parameter :: NETCDF4 = 5
   integer, parameter :: XDMF = 7
   character(len=20) :: formatString(0:7) = ['OFF         ', &
                                             'ASCII       ', 'UNDEF       ', 'NETCDF3     ', 'SPARSE_ASCII', &
                                             'NETCDF4     ', 'UNDEF       ', 'XDMF        ']

   !
   ! parameters related to input files
   !
   ! mesh (fort.14) file
   integer :: meshType
   logical :: meshFileNameSpecified ! .true. if the user explicitly gave the fn
   character(2048) :: meshFileName ! if different from the default fort.14
   ! name of full domain mesh file in parallel (only used by XDMF in parallel)
   character(2048) :: meshFileName_g
   !
   ! control (fort.15) file
   integer :: controlType
   logical :: controlFileNameSpecified ! .true. if the user explicitly gave the fn
   character(2048) :: controlFileName ! if different from the default (fort.15)
   ! name of full domain mesh file in parallel (only used by XDMF in parallel)
   character(2048) :: controlFileName_g
   !
   ! nodal attributes (fort.13) file
   integer :: naType
   logical :: naFileNameSpecified ! .true. if the user explicitly gave the fn
   character(2048) :: naFileName ! if different from the default fort.14

   !---------------------end of data declarations--------------------------------C

contains

   !-----------------------------------------------------------------------
   !     S U B R O U T I N E   M A K E _ D I R N A M E
   !-----------------------------------------------------------------------
   !     jgf: Determine which directories (paths) should be used
   !     to locate input files and write output files; process command
   !     line options related to file i/o.
   !
   !     jgf51.21.11: Added capability to allow user to specify name and
   !     file format of the mesh file.
   !-----------------------------------------------------------------------
   subroutine MAKE_DIRNAME()
      !tcm v50.85 added for windows builds
#ifdef CMPI
      use mpi_f08, only: MPI_Abort, MPI_COMM_WORLD
#ifdef WINDOWS
      use IFPORT, only: makedirqq
#endif
#endif
      implicit none
#if defined(CMPI) && !defined(WINDOWS)
      interface
         ! C routine in prep/mkdir.c; takes a Fortran string with hidden length
         subroutine makedir(path)
            implicit none
            character(len=*), intent(in) :: path
         end subroutine makedir
      end interface
#endif
      integer :: ARGCOUNT, I
      character(2048) :: CMDLINEARG
      character(2000) :: wholeCommandLine ! used to launch this executable, for log message
#ifdef CMPI
      integer :: LNAME, iprefix, ierr
      character(8), parameter :: PREFIX(2) = ['/PE0000 ', '/DOM0000']
      logical :: fileFound
#ifdef WINDOWS
      logical(4) :: dir_result
#endif
#endif

      INPUTDIR = ""
      GLOBALDIR = ""
      LOCALDIR = ""
      ARGCOUNT = command_argument_count()
      WRITE_LOCAL_HOT_START_FILES = .false.
      READ_LOCAL_HOT_START_FILES = .false.
      WRITE_LOCAL_FILES = .false.
      WRITE_LOCAL_HARM_FILES = .false.
      !
      ! defaults
      meshType = ASCII
      meshFileNameSpecified = .false.
      meshFileName = 'fort.14'
      meshFileName_g = 'fort.14' ! used by XDMF writer in parallel to read mesh

      controlType = ASCII
      controlFileNameSpecified = .false.
      controlFileName = 'fort.15'
      controlFileName_g = 'fort.15' ! used by XDMF writer in parallel to read metadata

      naType = ASCII
      naFileNameSpecified = .false.
      naFileName = 'fort.13'

      ! Casey 090527: Debug.
      ROOTDIR = "."

      ! Cobell 110518: Add compiler flag for local hot start.
#ifdef LOCALHOT
      WRITE_LOCAL_HOT_START_FILES = .true.
#else
      WRITE_LOCAL_HOT_START_FILES = .false.
#endif

      if (ARGCOUNT > 0) then
         !
         ! jgf51.21.41: Echo the entire command line that was used to
         ! launch this executable.
         i = 0
         wholeCommandLine(:) = ' '
         do while (i < argcount)
            call get_command_argument(i, cmdlinearg)
            wholeCommandLine = ' '//trim(adjustl(wholeCommandLine)) &
                               //' '//trim(adjustl(cmdlinearg))
            i = i + 1
         end do
         if (myproc == 0) then
            write (*, '(a,a)') 'ECHO: Application launched with: ' &
               //trim(adjustl(wholeCommandLine))
         end if
         !
         ! Now process the command line arguments.
         i = 0
         do while (i < ARGCOUNT)
            i = i + 1
            call get_command_argument(i, CMDLINEARG)
            select case (cmdlinearg(1:2))
            case ("-I", "-i")
               i = i + 1
               call get_command_argument(i, INPUTDIR)
               if (myProc == 0) then
                  write (*, *) "INFO: Processing '", cmdlinearg(1:2), " ", &
                     trim(INPUTDIR), "'."
               end if
            case ("-O", "-o")
               i = i + 1
               call get_command_argument(i, GLOBALDIR)
               if (myProc == 0) then
                  write (*, *) "INFO: Processing '", cmdlinearg(1:2), " ", &
                     trim(GLOBALDIR), "'."
               end if
            case ("-L", "-l")
               WRITE_LOCAL_FILES = .true.
               WRITE_LOCAL_HOT_START_FILES = .true.
               if (myProc == 0) then
                  write (*, *) "INFO: Processing '", cmdlinearg(1:2), "."
               end if
            case ("-U", "-u")
               WRITE_LOCAL_FILES = .true.
               if (myProc == 0) then
                  write (*, *) "INFO: Processing '", cmdlinearg(1:2), "."
               end if
            case ("-S", "-s")
               WRITE_LOCAL_HOT_START_FILES = .true.
               if (myProc == 0) then
                  write (*, *) "INFO: Processing '", cmdlinearg(1:2), "."
               end if
            case ("-R", "-r")
               READ_LOCAL_HOT_START_FILES = .true.
               if (myProc == 0) then
                  write (*, *) "INFO: Processing '", cmdlinearg(1:2), "."
               end if
            case ("-W", "-H")
               i = i + 1
               ! this is the number of writer processors, and will
               ! be parsed in a later call to GET_NUMWRITERS
               !  ... could also be -Ws for round-robin writers.
               ! -H (hotstart writers) is no longer supported;
               ! ADCIRC_Init stops with an error if it is given
            case ("-M", "-m")
               ! jgf51.21.11: Added options related to the mesh file
               select case (cmdlinearg(1:4))
               case ("-MFT", "-mft") ! mesh file type
                  i = i + 1
                  call get_command_argument(i, CMDLINEARG)
                  select case (trim(cmdlinearg))
                  case ("ASCII", "Ascii", "ascii")
                     meshType = ASCII ! this is the default anyway
                  case ("XDMF", "Xdmf", "xdmf")
                     meshType = XDMF
                  case default
                     if (myProc == 0) then
                        write (*, *) "WARNING: The mesh file type '", &
                           trim(cmdlinearg), &
                           "' is not valid. Valid types are ascii and xdmf."
                     end if
                  end select
               case ("-MFN", "-mfn") ! mesh file name
                  meshFileNameSpecified = .true.
                  i = i + 1
                  call get_command_argument(i, meshFileName)
               case default
                  ! jgf51.44: -m by itself triggers subdomain harmonic
                  ! analysis output files.
                  if (len(trim(cmdlinearg)) == 2) then
                     WRITE_LOCAL_HARM_FILES = .true.
                     if (myProc == 0) then
                        write (*, '(a,a,a)') "INFO: Processing " &
                           //cmdlinearg(1:2)//"."
                     end if
                  end if
                  if (myProc == 0) then
                     write (*, *) "WARNING: The command line option '", &
                        cmdlinearg(1:4), "' is not valid and will be ignored."
                  end if
               end select
            case ("-N", "-n")
               ! jgf51.21.11: Added options related to the nodal attributes file
               select case (cmdlinearg(1:4))
               case ("-NFT", "-nft") ! nodal attributes file type
                  i = i + 1
                  call get_command_argument(i, CMDLINEARG)
                  select case (trim(cmdlinearg))
                  case ("ASCII", "Ascii", "ascii")
                     naType = ASCII ! this is the default anyway
                  case ("XDMF", "Xdmf", "xdmf")
                     naType = XDMF
                  case default
                     if (myProc == 0) then
                        write (*, *) &
                           "WARNING: The nodal attributes file type '", &
                           trim(cmdlinearg), &
                           "' is not valid. Valid types are ascii and xdmf."
                     end if
                  end select
               case ("-NFN", "-nfn") ! nodal attributes file name
                  naFileNameSpecified = .true.
                  i = i + 1
                  call get_command_argument(i, naFileName)
               case default
                  if (myProc == 0) then
                     write (*, *) "WARNING: The command line option '", &
                        cmdlinearg(1:4), "' is not valid and will be ignored."
                  end if
               end select
            case default
               if (myProc == 0) then
                  write (*, *) "WARNING: The command line option '", &
                     cmdlinearg(1:2), "' is not valid and will be ignored."
               end if
            end select
         end do
      end if

      !.....Default root working directory

      if (len_trim(INPUTDIR) /= 0) then
         ROOTDIR = INPUTDIR
      end if

      GBLINPUTDIR = ROOTDIR
#ifdef CMPI
      iprefix = 0
      if (myProc == 0) then
         write (*, *) "INFO: Searching for ADCIRC subdomain directories:"
      end if
      do i = 1, 2
         if (myProc == 0) then
            write (*, *) "INFO: Looking for '", trim(ROOTDIR), &
               trim(PREFIX(i)), "/fort.14' ..."
         end if
         inquire (file=trim(ROOTDIR)//trim(PREFIX(i))//'/'//'fort.14', &
                  exist=fileFound)
         if (fileFound .eqv. .true.) then
            if (myProc == 0) then
               write (*, *) "INFO: File '", trim(ROOTDIR), &
                  trim(PREFIX(i)), "/fort.14' was found!"
               write (*, *) "INFO: The search for the subdomain ", &
                  "directory was completed successfully."
            end if
            iprefix = i
            exit
         else
            write (*, *) "ERROR: Processor ", myProc, ": File '", &
               trim(ROOTDIR), trim(PREFIX(i)), &
               "/fort.14' was not found."

            print *, "ERROR: ADCIRC stopping."
            call mpi_abort(mpi_comm_world, 1, ierr)
         end if
      end do
!
      write (INPUTDIR, '(2A)') trim(ROOTDIR), PREFIX(iprefix)
      LNAME = len_trim(INPUTDIR)
      write (INPUTDIR(LNAME - 3:LNAME), '(I4.4)') MYPROC
#else
      write (INPUTDIR, '(A)') trim(ROOTDIR)
#endif

      if (len_trim(GLOBALDIR) /= 0) then
         ROOTDIR = GLOBALDIR
      end if

      write (GLOBALDIR, '(A)') trim(ROOTDIR)

#ifdef CMPI
      write (LOCALDIR, '(2A)') trim(ROOTDIR), trim(PREFIX(iprefix))
      LNAME = len_trim(LOCALDIR)
      write (LOCALDIR(LNAME - 3:LNAME), '(I4.4)') MYPROC

!... tcm v49.67 -- additions for Windows based PC compilations
#ifdef WINDOWS
!call system('mkdir '//trim(LOCALDIR))  !tcm v50.85 this no longer works for INTEL MPI
      dir_result = makedirqq(trim(localdir)) !tcm v50.85 added
#else
      call MAKEDIR(trim(LOCALDIR))
#endif

#else
      write (LOCALDIR, '(A)') trim(ROOTDIR)
#endif
      !
      ! jgf51.21.27: Form mesh file name and nodal attributes file name.
      if (meshFileNameSpecified .eqv. .false.) then
         meshFileName = trim(LOCALDIR)//'/fort.14'
      end if
      if (naFileNameSpecified .eqv. .false.) then
         naFileName = trim(LOCALDIR)//'/fort.13'
      end if
      !
      !     jgf49.17.01 Summarize and log results.
      if (myProc == 0) then
         write (*, *) "INFO: The ROOTDIR is  '", trim(ROOTDIR), "'."
         write (*, *) "INFO: The INPUTDIR is '", trim(INPUTDIR), "'."
         write (*, *) "INFO: The GBLINPUTDIR is '", trim(GBLINPUTDIR), "'."
         write (*, *) "INFO: The GLOBALDIR is '", trim(GLOBALDIR), "'."
         write (*, *) "INFO: The LOCALDIR is '", trim(LOCALDIR), "'."
         if ((WRITE_LOCAL_FILES .eqv. .true.) .and. (MNPROC > 1)) then
            write (*, *) "INFO: ADCIRC will write subdomain output files."
         end if
         if ((WRITE_LOCAL_HOT_START_FILES .eqv. .true.) &
             .and. (MNPROC > 1)) then
            write (*, *) &
               "INFO: ADCIRC will write subdomain hotstart files."
         end if
         if ((WRITE_LOCAL_HARM_FILES .eqv. .true.) .and. (MNPROC > 1)) then
            write (*, *) 'INFO: ADCIRC will write ' &
               //'subdomain harmonic analysis files.'
         end if
      end if
      return
      !-----------------------------------------------------------------------
   end subroutine MAKE_DIRNAME
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !     S U B R O U T I N E   G E T _ N U M  W R I T E R S
   !-----------------------------------------------------------------------
   !-----------------------------------------------------------------------
   subroutine GET_NUMWRITERS
      implicit none
      integer :: ARGCOUNT, I
      character(2048) :: CMDLINEARG
      ARGCOUNT = command_argument_count()
      MNWPROC = 0
      isplit = .false. !st3 100708:split file

      if (ARGCOUNT > 0) then
         i = 0
         do while (i < ARGCOUNT)
            i = i + 1
            call get_command_argument(i, CMDLINEARG)
            if (cmdlinearg(1:2) == "-W") then
               if (cmdlinearg(1:3) == "-Ws") then !st3 100708: split file
                  isplit = .true. !st3 100708: split file
               end if !st3 100708: split file
               i = i + 1
               call get_command_argument(i, cmdlinearg)
               read (cmdlinearg, *) MNWPROC
            end if
         end do
      end if
      !
      ! Dedicated hotstart writer processors (-H) are no longer
      ! supported; note that the option was given so that ADCIRC_Init
      ! can stop with an error once MPI is running
      hotstartWritersRequested = .false.
      if (ARGCOUNT > 0) then
         do i = 1, ARGCOUNT
            call get_command_argument(i, CMDLINEARG)
            if (cmdlinearg(1:2) == "-H") then
               hotstartWritersRequested = .true.
            end if
         end do
      end if

      !.....Default root working directory

      return
      !-----------------------------------------------------------------------
   end subroutine GET_NUMWRITERS
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !-----------------------------------------------------------------------
end module MOD_SIZES
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
