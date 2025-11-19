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
! PADCIRC VERSION 47.xx 10/13/2006                                             *
!******************************************************************************
!
module MESSENGER

   use mpi_f08, only: MPI_Comm, MPI_Group, MPI_Request, MPI_Status, MPI_Datatype

   implicit none

   private

   !--------------------------------------------------------------------------
   !  This module supplies the MPI Message-Passing Interface for PADCIRC.
   !  Uses asynchronous communication with buffer packing as performance
   !  enhancement for "cluster" architectures.
   !--------------------------------------------------------------------------

   !  Message-Passing Array space
   type(MPI_Comm) :: MPI_COMM_ADCIRC ! Local communicator
   integer :: NEIGHPROC ! number of subdomains neighboring this one
   integer :: RDIM ! total send/recieve requests when passing info
   ! between subdomain neighbors, =2*neighproc
   integer :: IERR ! error status of mpi subroutine call
   integer, parameter :: TAG = 100
   type(MPI_Comm) :: COMM_COMP ! COMMUNICATOR FOR COMPUTATION
   type(MPI_Group) :: GROUP_WORLD, GROUP_COMP !, GROUP_WRITER_ONLY
   type(MPI_Group), allocatable :: GROUP_WRITER(:) ! GROUPS FOR GLOBAL FILE WRITING
   type(MPI_Group), allocatable :: GROUP_WRITEH(:) ! GROUPS FOR HOTSTART FILE WRITING !st3 100711
   type(MPI_Group), allocatable :: GROUP_HSLEEP(:) ! GROUPS FOR HOTSTART FILE WRITING !st3 100711
   logical, allocatable :: RESNODE(:)

   integer, allocatable :: IPROC(:), NNODELOC(:), NNODSEND(:), &
                           NNODRECV(:), IBELONGTO(:), ISENDLOC(:, :), IRECVLOC(:, :), &
                           ISENDBUF(:, :), IRECVBUF(:, :)

   type(MPI_Request), allocatable :: REQ_I1(:), REQ_I2(:)
   type(MPI_Status), allocatable :: STAT_I1(:), STAT_I2(:)
   type(MPI_Request), allocatable :: REQ_R1(:), REQ_R2(:), REQ_R3(:)
   type(MPI_Status), allocatable :: STAT_R1(:), STAT_R2(:), STAT_R3(:)
   type(MPI_Request), allocatable :: REQ_M4R(:)
   type(MPI_Status), allocatable :: STAT_M4R(:) !sb 10/13/2022
   type(MPI_Request), allocatable :: REQ_R3D(:)
   type(MPI_Status), allocatable :: STAT_R3D(:)
   type(MPI_Request), allocatable :: REQ_C3D(:)
   type(MPI_Status), allocatable :: STAT_C3D(:)
   integer, allocatable :: INDX(:)
   real(8), allocatable :: SENDBUF(:, :), RECVBUF(:, :)
   complex(8), allocatable :: SENDBUF_COMPLEX(:, :), RECVBUF_COMPLEX(:, :)
   !jgf50.82: Create a flag for unrecoverable issue on a subdomain
   logical :: subdomainFatalError ! true if mpi_abort should be called
   logical :: writers_active = .false.
   logical :: hs_writers_active = .false.

   public :: MSG_INIT, MSG_FINI, MSG_TABLE, MSG_START, UPDATEI, UPDATER, &
             MSG_RScalar_Reduce, MSG_BARRIER, MSG_RBCASTD, MSG_RBCAST, &
             MSG_CBCAST, MSG_LBCAST, MSG_IBCAST, EarlyTermSum, WarnElevSum, &
             WetDrySum, MapToSubdomainIntMpi, AllNodes, ps3dots, ps2dots, &
             msg_imax, updatec3d, updater_w_perbc, mapToSubdomainrealmpi, &
             updater3d, updatem4r, psdot, subdomainFatalError, tag, &
             mpi_comm_adcirc, writers_active, hs_writers_active, resnode, &
             msg_abort, neighproc, nnodrecv, irecvloc

contains

   !---------------------end of data declarations--------------------------------C

   !----------------------------------------------------------------------
   !   S U B R O U T I N E    M S G _ I N I T
   !----------------------------------------------------------------------
   !  Routine performs following steps:
   !   (1)  define mpi data types to be used
   !   (2)  initialize MPI,
   !   (3)  get number of processors,
   !   (4)  get MPI rank of processor
   !   (5)  initialize adcirc COMM communicator MPI_COMM_WORLD
   !  vjp  10/3/2006
   !----------------------------------------------------------------------
   subroutine MSG_INIT(MPI_COMM_IN)
      use mpi_f08, only: MPI_COMM_WORLD, MPI_Comm_dup, &
                         MPI_Comm_size, MPI_Comm_rank, MPI_Group_incl, MPI_Comm_create, &
                         MPI_Comm_group, MPI_Init, MPI_Comm
      use SIZES, only: MNALLPROC, MNPROC, MYPROC, MNWPROH, MNWPROC
      use GLOBAL, only: COMM_WRITER, COMM_WRITEH, COMM_HSLEEP, setMessageSource, &
                        screenMessage, allMessage, unsetMessageSource, WRITER_ID, COMM
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      use GLOBAL, only: DEBUG
#endif
      implicit none
      type(MPI_Comm), optional, intent(in) :: MPI_COMM_IN
      integer :: I
      integer, allocatable :: RANKS(:) ! array of mpi ranks for compute processors
      integer :: IRANK_SLEEP(2) !st3 100711  for hsfile
      call setMessageSource("msg_init")
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call screenMessage(DEBUG, "Enter.") ! log to screen; don't have log dirname
#endif
      subdomainFatalError = .false.
      if (present(MPI_COMM_IN)) then
         !.......Duplicate communicator passed from outside
         call MPI_Comm_dup(MPI_COMM_IN, MPI_COMM_ADCIRC, IERR)
      else
         !.......Initialize MPI
         call MPI_Init(IERR)
         !.......Duplicate communicator
         call MPI_Comm_dup(MPI_COMM_WORLD, MPI_COMM_ADCIRC, IERR)
      end if

      call MPI_COMM_SIZE(MPI_COMM_ADCIRC, MNALLPROC, IERR) ! Get number of procs
      call MPI_COMM_RANK(MPI_COMM_ADCIRC, MYPROC, IERR) ! Get MPI rank

      MNPROC = MNALLPROC - MNWPROC & ! MNALLPROC = MNPROC + MNWPROC
               - MNWPROH !             + MNWPROH  !st3 100711 for hsfile

      allocate (RANKS(MNPROC + 1))

      !...  Create a communicator for computation
      do I = 1, MNPROC
         RANKS(I) = I - 1
      end do
      ! jgf51.21.27: Create group and communicator consisting
      ! only of compute processors
      call MPI_COMM_GROUP(MPI_COMM_ADCIRC, GROUP_WORLD, IERR)
      call MPI_GROUP_INCL(GROUP_WORLD, MNPROC, RANKS, GROUP_COMP, IERR)
      call MPI_COMM_CREATE(MPI_COMM_ADCIRC, GROUP_COMP, COMM_COMP, IERR)
      WRITER_ID = 0

      ! if we have dedicated writer processors, then for each one
      ! create a group and communicator that includes
      ! that one dedicated writer processor and all the compute
      ! processors ... so if we have 4 dedicated writer processors
      ! there will be 4 writer communicators and 4 writer groups
      ! ... in each one of these the writer is the highest rank processor
      ! in the group (the writer will be of rank mnproc)
      if (MNWPROC > 0) then
         !...     Allocate memory for groups and communicators
         allocate (GROUP_WRITER(MNWPROC), COMM_WRITER(MNWPROC))
         !...     Create communicators for global file writings
         do I = 1, MNWPROC
            RANKS(MNPROC + 1) = MNPROC - 1 + I
            call MPI_GROUP_INCL(GROUP_WORLD, MNPROC + 1, RANKS, &
                                GROUP_WRITER(I), IERR)
            call MPI_COMM_CREATE(MPI_COMM_ADCIRC, GROUP_WRITER(I), &
                                 COMM_WRITER(I), IERR)

            if (MYPROC == MNPROC - 1 + I) then
               WRITER_ID = I
            end if
         end do
      end if
      !st3 for hsfile 05.14.2010
#ifdef ADCNETCDF
      MNWPROH = 0 !jgfdebug20120215: why is this here?
#endif
#ifndef CMPI
      MNWPROH = 0
#endif
      if (MNWPROH > 0) then
         IRANK_SLEEP(1) = 0
         allocate (GROUP_WRITEH(MNWPROH), COMM_WRITEH(MNWPROH))
         allocate (GROUP_HSLEEP(MNWPROH), COMM_HSLEEP(MNWPROH))
         do I = 1, MNWPROH
            RANKS(MNPROC + 1) = MNPROC + MNWPROC - 1 + I

            call MPI_GROUP_INCL(GROUP_WORLD, MNPROC + 1, RANKS, &
                                GROUP_WRITEH(I), IERR)
            call MPI_COMM_CREATE(MPI_COMM_ADCIRC, GROUP_WRITEH(I), &
                                 COMM_WRITEH(I), IERR)

            IRANK_SLEEP(2) = MNPROC + MNWPROC - 1 + I
            call MPI_GROUP_INCL(GROUP_WORLD, 2, IRANK_SLEEP, &
                                GROUP_HSLEEP(I), IERR)
            call MPI_COMM_CREATE(MPI_COMM_ADCIRC, GROUP_HSLEEP(I), &
                                 COMM_HSLEEP(I), IERR)

            if (MYPROC == MNPROC + MNWPROC - 1 + I) then
               WRITER_ID = I
            end if
         end do
      end if
      deallocate (RANKS)
      COMM = COMM_COMP

#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
      !---------------------------------------------------------------------
   end subroutine MSG_INIT
   !---------------------------------------------------------------------

   !---------------------------------------------------------------------
   !      S U B R O U T I N E   M S G _ F I N I
   !---------------------------------------------------------------------
   !  Delete MPI resources and Shutdown MPI library.
   !  vjp  8/29/1999
   !---------------------------------------------------------------------
   subroutine MSG_FINI(NO_MPI_FINALIZE)
      use mpi_f08, only: MPI_Send, MPI_INTEGER, MPI_Barrier, MPI_Abort, MPI_Finalize, &
                         MPI_COMM_WORLD

      use SIZES, only: MNWPROC, MYPROC, MNPROC, MNWPROH
      use GLOBAL, only: SIG_TERM, COMM_WRITER, COMM_WRITEH, COMM_HSLEEP, &
                        CPL2STWAVE, Flag_ElevError, Flag_VelError, setMessageSource, &
                        allMessage, unsetMessageSource
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      use GLOBAL, only: DEBUG
#endif
      implicit none
      logical, optional, intent(in) :: NO_MPI_FINALIZE
      integer :: I

      call setMessageSource("msg_fini")
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      if (MNWPROC > 0 .and. writers_active) then
         if (MYPROC == 0) then
            do I = 1, MNWPROC
               write (16, *) 'PROC ', MYPROC, ' IS SENDING SIG_TERM TO WRITER ', &
                  I
               call MPI_SEND(SIG_TERM, 1, MPI_INTEGER, MNPROC, &
                             TAG, COMM_WRITER(I), IERR)
            end do
         end if
      end if

      if (MNWPROH > 0 .and. hs_writers_active) then !st3 100711 for hsfile
         if (MYPROC == 0) then
            do I = 1, MNWPROH
               write (16, *) 'PROC ', MYPROC, ' IS SENDING SIG_TERM TO HSWRITER', &
                  I
               call MPI_BARRIER(COMM_HSLEEP(I), IERR)
               call MPI_SEND(SIG_TERM, 1, MPI_INTEGER, MNPROC, &
                             TAG, COMM_WRITEH(I), IERR)
            end do
         end if
      end if

      if (subdomainFatalError .eqv. .true.) then
         ! jgf50.82: Return the rank of the offending processor
         ! as the error code
         call MPI_ABORT(MPI_COMM_ADCIRC, MYPROC, IERR)
      end if

      ! tcm v51.32  added a "go nuclear" option for killing
      ! all MPI processes when Elevation Greater than Error
      ! elevation is exceeded, not just the mpi processes owned
      ! by ADCIRC when coupled with STWAVE via CSTORM-MS
      ! Note: mpi_comm_world is defined in mpich.f
      if (CPL2STWAVE .eqv. .true.) then
         if (Flag_ElevError .eqv. .true.) then
            call mpi_abort(mpi_comm_world, myproc, ierr)
         end if
         ! DMW 202401 Check for velocity exceeding error level
         if (Flag_VelError .eqv. .true.) then
            call mpi_abort(mpi_comm_world, myproc, ierr)
         end if
      end if

      if (present(NO_MPI_FINALIZE)) then
         if (.not. NO_MPI_FINALIZE) then
            call MPI_FINALIZE(IERR)
            if (MYPROC == 0) &
               print *, "MPI terminated with Status = ", IERR
         end if
      else
         call MPI_FINALIZE(IERR)
         if (MYPROC == 0) &
            print *, "MPI terminated with Status = ", IERR
      end if

#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
      return
      !---------------------------------------------------------------------
   end subroutine MSG_FINI
   !---------------------------------------------------------------------

   !---------------------------------------------------------------------
   !     S U B R O U T I N E   M S G _ T A B L E
   !---------------------------------------------------------------------
   !  Routine preforms following steps:
   !
   !   (1) Read Message-Passing Information from file "fort.18"
   !   (2) Determine resident nodes: RESNODE(I) is true  if I is resident node
   !   (3) Determine ghost nodes:    RESNODE(I) is false if I is ghost node
   !   (4) Determine number of neighbor subdomains
   !   (5) MPI rank of each neighbor and number of ghosts nodes to receive
   !   (6) Read Message-Passing Receive List
   !   (7) MPI rank of each neighbor and number of ghosts nodes to send
   !   (8) Read Message-Passing Send List
   !  vjp  10/13/2006
   !
   !  tcm v50.21 20110610 -- Changed I8 to I12 formats
   !---------------------------------------------------------------------
   subroutine MSG_TABLE()
      use SIZES, only: MNE, MNP, INPUTDIR, MYPROC
      use GLOBAL, only: IMAP_EL_LG, NODES_LG, FileFmtVersion, NP_G, NE_G, &
                        NSTAE, NSTAV, NSTAM, NSTAC, NSTAE_G, NSTAV_G, NSTAM_G, NSTAC_G, &
                        C3D, CMP_VERSION_NUMBERS, IMAP_STAE_LG, IMAP_STAV_LG, &
                        IMAP_STAM_LG, IMAP_STAC_LG, ERROR, setMessageSource, &
                        allMessage, unsetMessageSource
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      use GLOBAL, only: DEBUG
#endif
      use GLOBAL_3DVS, only: NSTA3DD, NSTA3DV, NSTA3DT, &
                             NSTA3DD_G, NSTA3DV_G, NSTA3DT_G, &
                             IMAP_STA3DD_LG, IMAP_STA3DV_LG, IMAP_STA3DT_LG
      implicit none
      integer :: IDPROC, NLOCAL, I, J, jdumy_loc
      integer :: jdumy, jdumy_G, jdumy_max, inputFileFmtVn
      character(10) :: BlkName
      logical :: FileFound

      call setMessageSource("msg_table")
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      FileFound = .false.
      inquire (FILE=trim(INPUTDIR)//'/'//'fort.18', EXIST=FileFound)
      if (FileFound .eqv. .false.) then
         write (*, *) "ERROR: ", myProc, ": The file ", &
            trim(INPUTDIR)//'/'//'fort.18 was not found.'
         write (*, 9973) ! execution terminated
         call MSG_FINI()
         call exit(1)
      end if

      open (18, FILE=trim(INPUTDIR)//'/'//'fort.18', STATUS='OLD')

      read (18, 3020) BlkName, inputFileFmtVn
      if (trim(BlkName) /= 'FileFmt' .and. &
          .not. CMP_VERSION_NUMBERS(FileFmtVersion, inputFileFmtVn)) then
         write (16, *) 'File Format of Fort.18 does not match aborting'
         call exit(1)
      end if

      ! Read Global number of elements and Local-to_Global element map ( used by module global_io )
      ! Casey 100209: Changed from FMT=3015.
      read (18, '(8X,3I12)') NE_G
      allocate (IMAP_EL_LG(MNE))
      do I = 1, MNE
         read (18, *) IMAP_EL_LG(I)
      end do

      ! Read Global number of nodes and Local-to_Global node map ( used by module global_io )
      read (18, 3015) NP_G
      allocate (NODES_LG(MNP))
      do I = 1, MNP
         read (18, *) NODES_LG(I)
      end do

      !  This information is provided for relocalizing fort.15
      !  Just read past it.
      read (18, '(8X,I12)') jdumy ! nfluxf for subdomain
      read (18, '(8X,3I12)') jdumy_g, jdumy_max, jdumy_loc ! neta for subdomain
      do I = 1, jdumy_loc
         read (18, '(I12)') jdumy ! obnode_lg table
      end do

      ! Read Global indexes of Elevation Station nodes ( used by module global_io )
      read (18, 3015) NSTAE_G, jdumy_max, jdumy_loc !tcm v51.20.05
      if (NSTAE /= jdumy_loc) then
         call allMessage(ERROR, "Elevation Station Dimensioning Error in fort.18")
         call MSG_FINI()
         call exit(1)
      end if
      if (NSTAE > 0) then
         allocate (IMAP_STAE_LG(NSTAE))
         do I = 1, NSTAE
            read (18, '(I12)') IMAP_STAE_LG(I)
         end do
      end if

      ! Read Global indexes of Velocity  Station nodes ( used by module global_io )
      read (18, 3015) NSTAV_G, jdumy_max, jdumy_loc !tcm v51.20.05
      if (NSTAV /= jdumy_loc) then
         call allMessage(ERROR, "Velocity Station Dimensioning Error in fort.18")
         call MSG_FINI()
         call exit(1)
      end if
      if (NSTAV > 0) then
         allocate (IMAP_STAV_LG(NSTAV))
         do I = 1, NSTAV
            read (18, '(I12)') IMAP_STAV_LG(I)
         end do
      end if

      ! Read Global indexes of Meteorlogical Station nodes ( used by module global_io )
      read (18, 3015) NSTAM_G, jdumy_max, jdumy_loc !tcm v51.20.05
      if (NSTAM /= jdumy_loc) then
         call allMessage(ERROR, "Met. Station Dimensioning Error in fort.18")
         call MSG_FINI()
         call exit(1)
      end if
      if (NSTAM > 0) then
         allocate (IMAP_STAM_LG(NSTAM))
         do I = 1, NSTAM
            read (18, '(I12)') IMAP_STAM_LG(I)
         end do
      end if

      ! Read Global indexes of Concentration Station nodes ( used by module global_io )
      read (18, 3015) NSTAC_G, jdumy_max, jdumy_loc !tcm v51.20.05
      if (NSTAC /= jdumy_loc) then
         call allMessage(ERROR, "Conc. Station Dimensioning Error in fort.18")
         call MSG_FINI()
         call exit(1)
      end if
      if (NSTAC > 0) then
         allocate (IMAP_STAC_LG(NSTAC))
         do I = 1, NSTAC
            read (18, '(I12)') IMAP_STAC_LG(I)
         end do
      end if

      !---------------------------------------------------------------------------------
      !--Message-Passing tables start here
      !---------------------------------------------------------------------------------

      read (18, 3010) IDPROC, NLOCAL
      allocate (NNODELOC(NLOCAL))
      read (18, 1130) (NNODELOC(I), I=1, NLOCAL)
      allocate (IBELONGTO(MNP), RESNODE(MNP))

      do I = 1, MNP
         IBELONGTO(I) = 0
      end do
      do I = 1, NLOCAL
         IBELONGTO(NNODELOC(I)) = IDPROC + 1
      end do
      do I = 1, MNP
         if (IBELONGTO(I) - 1 == MYPROC) then
            RESNODE(I) = .true.
         else
            RESNODE(I) = .false.
         end if
      end do

      read (18, 3015) NEIGHPROC

      RDIM = 2*NEIGHPROC
      allocate (INDX(RDIM))
      allocate (IPROC(NEIGHPROC), NNODRECV(NEIGHPROC))
      allocate (IRECVLOC(MNP, NEIGHPROC))

      do J = 1, NEIGHPROC
         read (18, 3010) IPROC(J), NNODRECV(J)
         read (18, 1130) (IRECVLOC(I, J), I=1, NNODRECV(J))
      end do

      allocate (NNODSEND(NEIGHPROC))
      allocate (ISENDLOC(MNP, NEIGHPROC))

      do J = 1, NEIGHPROC
         read (18, 3010) IPROC(J), NNODSEND(J)
         read (18, 1130) (ISENDLOC(I, J), I=1, NNODSEND(J))
      end do

      !     jgf49.43.18: Add 3D station mappings if appropriate. Used by globalio.
      if (C3D) then
         read (18, 3015) NSTA3DD_G
         if (NSTA3DD > 0) then
            allocate (IMAP_STA3DD_LG(NSTA3DD))
            do I = 1, NSTA3DD
               read (18, '(I12)') IMAP_STA3DD_LG(I)
            end do
         end if

         read (18, 3015) NSTA3DV_G
         if (NSTA3DV > 0) then
            allocate (IMAP_STA3DV_LG(NSTA3DV))
            do I = 1, NSTA3DV
               read (18, '(I12)') IMAP_STA3DV_LG(I)
            end do
         end if

         read (18, 3015) NSTA3DT_G
         if (NSTA3DT > 0) then
            allocate (IMAP_STA3DT_LG(NSTA3DT))
            do I = 1, NSTA3DT
               read (18, '(I12)') IMAP_STA3DT_LG(I)
            end do
         end if
      end if

      close (18)
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
      return
1130  format(8x, 6i12)
3010  format(8x, 2i12)
3015  format(8x, 3i12)
3020  format(a8, I12)
9973  format(/, 1x, '!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!', //)
      !---------------------------------------------------------------------
   end subroutine MSG_TABLE
   !---------------------------------------------------------------------

   !---------------------------------------------------------------------
   !     S U B R O U T I N E   M S G _ S T A R T
   !---------------------------------------------------------------------
   !  Routine allocates message-passing space
   !  vjp  10/01/2006
   !---------------------------------------------------------------------
   subroutine MSG_START()
      use SIZES, only: MNP, MNFEN
      use GLOBAL, only: C3D, setMessageSource, allMessage, unsetMessageSource
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      use GLOBAL, only: DEBUG
#endif
      implicit none
      call setMessageSource("msg_start")
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      allocate (ISENDBUF(MNP, NEIGHPROC), IRECVBUF(MNP, NEIGHPROC))
      allocate (REQ_I1(RDIM), REQ_I2(RDIM))
      allocate (REQ_R1(RDIM), REQ_R2(RDIM), REQ_R3(RDIM))
      allocate (STAT_I1(RDIM), STAT_I2(RDIM))
      allocate (STAT_R1(RDIM), STAT_R2(RDIM), STAT_R3(RDIM))
      allocate (REQ_M4R(RDIM)) !sb 10/13/2022
      allocate (STAT_M4R(RDIM)) !sb 10/13/2022

      if (C3D) then
         allocate (SENDBUF(2*MNP*MNFEN, NEIGHPROC))
         allocate (RECVBUF(2*MNP*MNFEN, NEIGHPROC))
         allocate (SENDBUF_COMPLEX(MNP*MNFEN, NEIGHPROC))
         allocate (RECVBUF_COMPLEX(MNP*MNFEN, NEIGHPROC))
         allocate (REQ_R3D(RDIM))
         allocate (STAT_R3D(RDIM))
         allocate (REQ_C3D(RDIM))
         allocate (STAT_C3D(RDIM))
      else
         allocate (SENDBUF(MNP, NEIGHPROC))
         allocate (RECVBUF(MNP, NEIGHPROC))
      end if
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
      return
      !---------------------------------------------------------------------
   end subroutine MSG_START
   !---------------------------------------------------------------------

   !---------------------------------------------------------------------
   !     S U B R O U T I N E   U P D A T E  I
   !---------------------------------------------------------------------
   !  Update 1 or 2 Integer Arrays's Ghost Cells using asynchronous
   !  and persistent message-passing.
   !
   !  vjp  8/06/1999
   !---------------------------------------------------------------------
   subroutine UPDATEI(IVEC1, IVEC2, NMSG)
      use mpi_f08, only: MPI_IRECV, MPI_ISEND, MPI_INTEGER, MPI_WAITSOME
      use GLOBAL, only: setMessageSource, allMessage, unsetMessageSource, COMM
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      use GLOBAL, only: DEBUG
#endif
      implicit none
      integer, intent(IN) :: NMSG
      integer, intent(INOUT) :: IVEC1(NMSG), IVEC2(NMSG)
      integer :: N, I, J, NCOUNT, NFINI, TOT

      call setMessageSource("updatei")
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      !..Pack 1 or 2 Messages
      do J = 1, NEIGHPROC
         NCOUNT = 0
         do I = 1, NNODSEND(J)
            NCOUNT = NCOUNT + 1
            ISENDBUF(NCOUNT, J) = IVEC1(ISENDLOC(I, J))
         end do
         if (NMSG > 1) then
            do I = 1, NNODSEND(J)
               NCOUNT = NCOUNT + 1
               ISENDBUF(NCOUNT, J) = IVEC2(ISENDLOC(I, J))
            end do
         end if
      end do
      ! Send/receive messages to/from all neighbors
      if (NMSG == 1) then
         do J = 1, NEIGHPROC
            call MPI_IRECV(IRECVBUF(1, J), NNODRECV(J), &
                           MPI_INTEGER, IPROC(J), TAG, COMM, REQ_I1(J), IERR)
            call MPI_ISEND(ISENDBUF(1, J), NNODSEND(J), &
                           MPI_INTEGER, IPROC(J), TAG, COMM, REQ_I1(J + NEIGHPROC), IERR)
         end do
      else
         do J = 1, NEIGHPROC
            call MPI_IRECV(IRECVBUF(1, J), 2*NNODRECV(J), &
                           MPI_INTEGER, IPROC(J), TAG, COMM, REQ_I2(J), IERR)
            call MPI_ISEND(ISENDBUF(1, J), 2*NNODSEND(J), &
                           MPI_INTEGER, IPROC(J), TAG, COMM, REQ_I2(J + NEIGHPROC), IERR)
         end do
      end if
      !..Unpack Received messages as they arrive
      if (NMSG == 1) then
         TOT = 0
         do while (TOT < RDIM)
            do N = 1, RDIM
               INDX(N) = 0
            end do
            call MPI_WAITSOME(RDIM, REQ_I1, NFINI, INDX, STAT_I1, IERR)
            TOT = TOT + NFINI
            do N = 1, NFINI
               if (INDX(N) > 0 .and. INDX(N) <= RDIM) then
                  if (INDX(N) <= NEIGHPROC) then
                     J = INDX(N)
                     NCOUNT = 0
                     do I = 1, NNODRECV(J)
                        NCOUNT = NCOUNT + 1
                        IVEC1(IRECVLOC(I, J)) = IRECVBUF(NCOUNT, J)
                     end do
                  end if
               end if
            end do
         end do
      else
         TOT = 0
         do while (TOT < RDIM)
            do N = 1, RDIM
               INDX(N) = 0
            end do
            call MPI_WAITSOME(RDIM, REQ_I2, NFINI, INDX, STAT_I2, IERR)
            TOT = TOT + NFINI
            do N = 1, NFINI
               if (INDX(N) > 0 .and. INDX(N) <= RDIM) then
                  if (INDX(N) <= NEIGHPROC) then
                     J = INDX(N)
                     NCOUNT = 0
                     do I = 1, NNODRECV(J)
                        NCOUNT = NCOUNT + 1
                        IVEC1(IRECVLOC(I, J)) = IRECVBUF(NCOUNT, J)
                     end do
                     do I = 1, NNODRECV(J)
                        NCOUNT = NCOUNT + 1
                        IVEC2(IRECVLOC(I, J)) = IRECVBUF(NCOUNT, J)
                     end do
                  end if
               end if
            end do
         end do
      end if

#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
      return
      !---------------------------------------------------------------------
   end subroutine UPDATEI
   !---------------------------------------------------------------------

   !---------------------------------------------------------------------
   !     S U B R O U T I N E   U P D A T E R
   !---------------------------------------------------------------------
   !  Update 1, 2, or 3 Real Arrays's Ghost Cells using asynchronous
   !  and persistent message-passing.
   !
   !  vjp  8/06/1999
   !---------------------------------------------------------------------
   subroutine UPDATER(VEC1, VEC2, VEC3, NMSG)
      use mpi_f08, only: MPI_IRECV, MPI_ISEND, MPI_WAITSOME, MPI_DOUBLE_PRECISION
      use GLOBAL, only: setMessageSource, allMessage, unsetMessageSource, COMM
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      use GLOBAL, only: DEBUG
#endif
      implicit none
      integer, intent(IN) :: NMSG ! number of arrays to pass
      real(8), intent(INOUT) :: VEC1(NMSG), VEC2(NMSG), VEC3(NMSG)
      integer :: j ! loop counter for neighboring subdomains
      integer :: i ! loop counter for nodes shared with a neighboring subdomain
      integer :: nfini ! number of just-completed requests
      integer :: n ! loop counter for just-completed requests
      integer :: tot ! total number of completed requests
      integer :: ncount ! num values passed to a particular subdomain neighbor

      call setMessageSource("updater")
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      !
      ! loop over neighboring subdomains
      do J = 1, NEIGHPROC
         NCOUNT = 0
         ! loop over nodes shared as ghost nodes with this neighboring subdomain
         do I = 1, NNODSEND(J)
            NCOUNT = NCOUNT + 1
            ! store the array values in the send buffer
            SENDBUF(NCOUNT, J) = VEC1(ISENDLOC(I, J))
         end do
         if (NMSG > 1) then
            do I = 1, NNODSEND(J)
               NCOUNT = NCOUNT + 1
               SENDBUF(NCOUNT, J) = VEC2(ISENDLOC(I, J))
            end do
         end if
         if (NMSG > 2) then
            do I = 1, NNODSEND(J)
               NCOUNT = NCOUNT + 1
               SENDBUF(NCOUNT, J) = VEC3(ISENDLOC(I, J))
            end do
         end if
      end do

      ! receive and send data to each neighboring subdomain, populating
      ! the request handler array
      if (NMSG == 1) then
         do J = 1, NEIGHPROC
            call MPI_IRECV(RECVBUF(1, J), NNODRECV(J), &
                           MPI_DOUBLE_PRECISION, IPROC(J), TAG, COMM, REQ_R1(J), IERR)
            call MPI_ISEND(SENDBUF(1, J), NNODSEND(J), &
                           MPI_DOUBLE_PRECISION, IPROC(J), TAG, COMM, REQ_R1(J + NEIGHPROC), IERR)
         end do
      elseif (NMSG == 2) then
         do J = 1, NEIGHPROC
            call MPI_IRECV(RECVBUF(1, J), 2*NNODRECV(J), &
                           MPI_DOUBLE_PRECISION, IPROC(J), TAG, COMM, REQ_R2(J), IERR)
            call MPI_ISEND(SENDBUF(1, J), 2*NNODSEND(J), &
                           MPI_DOUBLE_PRECISION, IPROC(J), TAG, COMM, REQ_R2(J + NEIGHPROC), IERR)
         end do
      else
         do J = 1, NEIGHPROC
            call MPI_IRECV(RECVBUF(1, J), 3*NNODRECV(J), &
                           MPI_DOUBLE_PRECISION, IPROC(J), TAG, COMM, REQ_R3(J), IERR)
            call MPI_ISEND(SENDBUF(1, J), 3*NNODSEND(J), &
                           MPI_DOUBLE_PRECISION, IPROC(J), TAG, COMM, REQ_R3(J + NEIGHPROC), IERR)
         end do
      end if

      select case (NMSG)
      case (1)
         TOT = 0
         ! keep looping until the total number of completed communications
         ! equals the total number that were requested
         do while (TOT < RDIM)
            do N = 1, RDIM
               INDX(N) = 0 ! zero out the array of just-completed requests
            end do
            call MPI_WAITSOME(RDIM, REQ_R1, NFINI, INDX, STAT_R1, IERR)
            ! add the number of just-completed requests to the total
            TOT = TOT + NFINI
            ! loop over the just-completed requests
            do N = 1, NFINI
               ! if the request has just completed (TODO: I don't see how
               ! indx(n) can be greater than rdim
               if (INDX(N) > 0 .and. INDX(N) <= RDIM) then
                  ! if the request that just completed is a receive
                  if (INDX(N) <= NEIGHPROC) then
                     ! j = subdomain neighbor number that sent this data
                     J = INDX(N)
                     NCOUNT = 0
                     do I = 1, NNODRECV(J)
                        NCOUNT = NCOUNT + 1
                        ! update this subdomain's ghost value with the real
                        ! value from the other subdomain
                        VEC1(IRECVLOC(I, J)) = RECVBUF(NCOUNT, J)
                     end do
                  end if
               end if
            end do
         end do
      case (2)
         TOT = 0
         do while (TOT < RDIM)
            do N = 1, RDIM
               INDX(N) = 0
            end do
            call MPI_WAITSOME(RDIM, REQ_R2, NFINI, INDX, STAT_R2, IERR)
            TOT = TOT + NFINI
            do N = 1, NFINI
               if (INDX(N) > 0 .and. INDX(N) <= RDIM) then
                  if (INDX(N) <= NEIGHPROC) then
                     J = INDX(N)
                     NCOUNT = 0
                     do I = 1, NNODRECV(J)
                        NCOUNT = NCOUNT + 1
                        VEC1(IRECVLOC(I, J)) = RECVBUF(NCOUNT, J)
                     end do
                     do I = 1, NNODRECV(J)
                        NCOUNT = NCOUNT + 1
                        VEC2(IRECVLOC(I, J)) = RECVBUF(NCOUNT, J)
                     end do
                  end if
               end if
            end do
         end do
      case DEFAULT
         TOT = 0
         do while (TOT < RDIM)
            do N = 1, RDIM
               INDX(N) = 0
            end do
            call MPI_WAITSOME(RDIM, REQ_R3, NFINI, INDX, STAT_R3, IERR)
            TOT = TOT + NFINI
            !debug     print *, myproc, tot,nfini,INDX(1),INDX(2)
            do N = 1, NFINI
               if (INDX(N) > 0 .and. INDX(N) <= RDIM) then
                  if (INDX(N) <= NEIGHPROC) then
                     J = INDX(N)
                     NCOUNT = 0
                     do I = 1, NNODRECV(J)
                        NCOUNT = NCOUNT + 1
                        VEC1(IRECVLOC(I, J)) = RECVBUF(NCOUNT, J)
                     end do
                     do I = 1, NNODRECV(J)
                        NCOUNT = NCOUNT + 1
                        VEC2(IRECVLOC(I, J)) = RECVBUF(NCOUNT, J)
                     end do
                     do I = 1, NNODRECV(J)
                        NCOUNT = NCOUNT + 1
                        VEC3(IRECVLOC(I, J)) = RECVBUF(NCOUNT, J)
                     end do
                  end if
               end if
            end do
         end do
      end select

#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
      return
      !---------------------------------------------------------------------
   end subroutine UPDATER
   !---------------------------------------------------------------------

   !---------------------------------------------------------------------
   !     S U B R O U T I N E   U P D A T E M A T 4 R
   !---------------------------------------------------------------------
   !  Update a 4-row real matrix's ghost cells using asynchronous
   !  and persistent message-passing.
   !
   !  sb  10/13/2022
   !---------------------------------------------------------------------
   subroutine UPDATEM4R(M4R)
      use mpi_f08, only: MPI_Irecv, MPI_Isend, MPI_Waitsome, MPI_DOUBLE_PRECISION
      use SIZES, only: MNP
      use GLOBAL, only: setMessageSource, allMessage, unsetMessageSource, COMM
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      use GLOBAL, only: DEBUG
#endif
      implicit none
      real(8), intent(INOUT) :: M4R(:, :)
      integer :: j ! loop counter for neighboring subdomains
      integer :: i ! loop counter for nodes shared with a neighboring subdomain
      integer :: nfini ! number of just-completed requests
      integer :: n ! loop counter for just-completed requests
      integer :: tot ! total number of completed requests
      integer :: ncount ! num values passed to a particular subdomain neighbor
      integer :: ir ! loop counter for rows
      !
      call setMessageSource("updater")
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      NCOUNT = 0

      ! loop over neighboring subdomains
      do J = 1, NEIGHPROC
         NCOUNT = 0
         do IR = 1, 4
            ! loop over nodes shared as ghost nodes with this neighboring subdomain
            do I = 1, NNODSEND(J)
               NCOUNT = NCOUNT + 1
               ! store the matrix values in the send buffer
               SENDBUF(NCOUNT, J) = M4R(IR, ISENDLOC(I, J))
            end do
         end do
      end do
      if (NCOUNT > MNP) then
         write (*, '(A)') &
            'INSUFFICIENT MEMORY ALLOCATION FOR SENDBUF. TERMINATING.'
         call MSG_FINI()
         stop
      end if

      ! receive and send data to each neighboring subdomain, populating
      ! the request handler array
      do J = 1, NEIGHPROC
         call MPI_IRECV(RECVBUF(1, J), 4*NNODRECV(J), &
                        MPI_DOUBLE_PRECISION, IPROC(J), TAG, COMM, REQ_M4R(J), IERR)
         call MPI_ISEND(SENDBUF(1, J), 4*NNODSEND(J), &
                        MPI_DOUBLE_PRECISION, IPROC(J), TAG, COMM, REQ_M4R(J + NEIGHPROC), IERR)
      end do

      TOT = 0
      ! keep looping until the total number of completed communications
      ! equals the total number that were requested
      do while (TOT < RDIM)
         do N = 1, RDIM
            INDX(N) = 0 ! zero out the array of just-completed requests
         end do
         call MPI_WAITSOME(RDIM, REQ_M4R, NFINI, INDX, STAT_M4R, IERR)
         ! add the number of just-completed requests to the total
         TOT = TOT + NFINI
         ! loop over the just-completed requests
         do N = 1, NFINI
            ! if the request has just completed (TODO: I don't see how
            ! indx(n) can be greater than rdim
            if (INDX(N) > 0 .and. INDX(N) <= RDIM) then
               ! if the request that just completed is a receive
               if (INDX(N) <= NEIGHPROC) then
                  ! j = subdomain neighbor number that sent this data
                  J = INDX(N)
                  NCOUNT = 0
                  do IR = 1, 4
                     do I = 1, NNODRECV(J)
                        NCOUNT = NCOUNT + 1
                        ! update this subdomain's ghost value with the real
                        ! value from the other subdomain
                        M4R(IR, IRECVLOC(I, J)) = RECVBUF(NCOUNT, J)
                     end do
                  end do
               end if
            end if
         end do
      end do

#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
      return
      !---------------------------------------------------------------------
   end subroutine UPDATEM4R
   !---------------------------------------------------------------------

   !..... DW:, for periodic bcs (this subroutine looks a bit ugly and probably does belongs here)
   !        - experiment periodic boundary conditions
   !        - use in updating the pseudo residual of slave nodes when they are ghost nodes
   !          0. Assign  the residual of slave nodes to a temporay vector
   !          1. Owner PE assigns the residual of the periodic-primary nodes to theirs slave nodes
   !          2. Update ghost nodes through CALL UPDATER()
   !          3. Restore the residual of slave nodes
   subroutine UPDATER_W_PERBC(VEC1, VEC2, VEC3, NMSG)
      use global, only: NPERSEG, NNPERBC, IPERCONN

      implicit none

      integer, intent(IN) :: NMSG
      real(8), intent(INOUT) :: VEC1(NMSG), VEC2(NMSG), VEC3(NMSG)

      real(8) :: VECTMP(NNPERBC, NMSG); 
      if (NPERSEG > 0) then
         VECTMP(:, 1) = VEC1(IPERCONN(1:NNPERBC, 2)); 
         VEC1(IPERCONN(1:NNPERBC, 2)) = VEC1(IPERCONN(1:NNPERBC, 1)); 
         if (NMSG == 2) then
            VECTMP(:, 2) = VEC2(IPERCONN(1:NNPERBC, 2)); 
            VEC2(IPERCONN(1:NNPERBC, 2)) = VEC2(IPERCONN(1:NNPERBC, 1)); 
         elseif (NMSG == 3) then
            VECTMP(:, 2) = VEC2(IPERCONN(1:NNPERBC, 2)); 
            VEC2(IPERCONN(1:NNPERBC, 2)) = VEC2(IPERCONN(1:NNPERBC, 1)); 
            VECTMP(:, 3) = VEC3(IPERCONN(1:NNPERBC, 2)); 
            VEC3(IPERCONN(1:NNPERBC, 2)) = VEC3(IPERCONN(1:NNPERBC, 1)); 
         end if
      end if

      call UPDATER(VEC1, VEC2, VEC3, NMSG); 
      if (NPERSEG > 0) then
         VEC1(IPERCONN(1:NNPERBC, 2)) = VECTMP(:, 1)
         if (NMSG == 2) then
            VEC2(IPERCONN(1:NNPERBC, 2)) = VECTMP(:, 2); 
         elseif (NMSG == 3) then
            VEC2(IPERCONN(1:NNPERBC, 2)) = VECTMP(:, 2); 
            VEC3(IPERCONN(1:NNPERBC, 2)) = VECTMP(:, 3); 
         end if
      end if

      return; 
   end subroutine UPDATER_W_PERBC
   !..... DW
   !
   !---------------------------------------------------------------------
   !     S U B R O U T I N E   U P D A T E R 3 D
   !---------------------------------------------------------------------
   !  Update 1 Three-dimensional Real Arrays's Ghost Cells using asynchronous
   !  and persistent message-passing.
   !
   !  tjc  6/24/2002
   !---------------------------------------------------------------------
   subroutine UPDATER3D(VEC)
      use mpi_f08, only: MPI_Irecv, MPI_Isend, MPI_Waitsome, MPI_DOUBLE_PRECISION
      use SIZES, only: MNP, MNFEN
      use GLOBAL, only: setMessageSource, allMessage, unsetMessageSource, COMM
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      use GLOBAL, only: DEBUG
#endif
      implicit none
      real(8), intent(INOUT) :: VEC(MNP, MNFEN)
      integer :: N, I, J, K, NCOUNT, NFINI, TOT

      call setMessageSource("updater3d")
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      !..Pack Messages
      do J = 1, NEIGHPROC
         NCOUNT = 0
         do I = 1, NNODSEND(J)
            do K = 1, MNFEN
               NCOUNT = NCOUNT + 1
               SENDBUF(NCOUNT, J) = VEC(ISENDLOC(I, J), K)
            end do
         end do
      end do

      ! Send/receive messages to/from all neighbors
      do J = 1, NEIGHPROC
         call MPI_IRECV(RECVBUF(1, J), MNFEN*NNODRECV(J), &
                        MPI_DOUBLE_PRECISION, IPROC(J), TAG, COMM, REQ_R3D(J), IERR)
         call MPI_ISEND(SENDBUF(1, J), MNFEN*NNODSEND(J), &
                        MPI_DOUBLE_PRECISION, IPROC(J), TAG, COMM, REQ_R3D(J + NEIGHPROC), IERR)
      end do

      !..Unpack Received messages as they arrive
      TOT = 0
      do while (TOT < RDIM)
         do N = 1, RDIM
            INDX(N) = 0
         end do
         call MPI_WAITSOME(RDIM, REQ_R3D, NFINI, INDX, STAT_R3D, IERR)
         TOT = TOT + NFINI
         do N = 1, NFINI
            if (INDX(N) > 0 .and. INDX(N) <= RDIM) then
               if (INDX(N) <= NEIGHPROC) then
                  J = INDX(N)
                  NCOUNT = 0
                  do I = 1, NNODRECV(J)
                     do K = 1, MNFEN
                        NCOUNT = NCOUNT + 1
                        VEC(IRECVLOC(I, J), K) = RECVBUF(NCOUNT, J)
                     end do
                  end do
               end if
            end if
         end do
      end do
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
      return
      !---------------------------------------------------------------------
   end subroutine UPDATER3D
   !---------------------------------------------------------------------

   !---------------------------------------------------------------------
   !     S U B R O U T I N E   U P D A T E C  3 D
   !---------------------------------------------------------------------
   !  Update 1 Three-dimensional Complex Arrays' Ghost Cells using asynchronous
   !  and persistent message-passing.
   !  tjc  6/24/2002
   !---------------------------------------------------------------------
   subroutine UPDATEC3D(VEC)
      use mpi_f08, only: MPI_Irecv, MPI_Isend, MPI_Waitsome, MPI_DOUBLE_COMPLEX
      use SIZES, only: MNP, MNFEN
      use GLOBAL, only: setMessageSource, allMessage, unsetMessageSource, COMM
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      use GLOBAL, only: DEBUG
#endif
      implicit none
      complex(8), intent(INOUT) :: VEC(MNP, MNFEN)
      integer :: N, I, J, K, NCOUNT, NFINI, TOT

      call setMessageSource("updatec3d")
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      !..Pack Messages
      do J = 1, NEIGHPROC
         NCOUNT = 0
         do I = 1, NNODSEND(J)
            do K = 1, MNFEN
               NCOUNT = NCOUNT + 1
               SENDBUF_COMPLEX(NCOUNT, J) = VEC(ISENDLOC(I, J), K)
            end do
         end do
      end do

      ! Send/receive messages to/from all neighbors
      do J = 1, NEIGHPROC
         call MPI_IRECV(RECVBUF_COMPLEX(1, J), MNFEN*NNODRECV(J), &
                        MPI_DOUBLE_COMPLEX, IPROC(J), TAG, COMM, REQ_C3D(J), IERR)
         call MPI_ISEND(SENDBUF_COMPLEX(1, J), MNFEN*NNODSEND(J), &
                        MPI_DOUBLE_COMPLEX, IPROC(J), TAG, COMM, REQ_C3D(J + NEIGHPROC), IERR)
      end do

      !..Unpack Received messages as they arrive
      TOT = 0
      do while (TOT < RDIM)
         do N = 1, RDIM
            INDX(N) = 0
         end do
         call MPI_WAITSOME(RDIM, REQ_C3D, NFINI, INDX, STAT_C3D, IERR)
         TOT = TOT + NFINI
         do N = 1, NFINI
            if (INDX(N) > 0 .and. INDX(N) <= RDIM) then
               if (INDX(N) <= NEIGHPROC) then
                  J = INDX(N)
                  NCOUNT = 0
                  do I = 1, NNODRECV(J)
                     do K = 1, MNFEN
                        VEC(IRECVLOC(I, J), K) = RECVBUF_COMPLEX(NCOUNT + 1, J)
                        NCOUNT = NCOUNT + 1
                     end do
                  end do
               end if
            end if
         end do
      end do

#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
   end subroutine UPDATEC3D
   !---------------------------------------------------------------------

   !---------------------------------------------------------------------
   !     F U N C T I O N   M S G _ M A X
   !---------------------------------------------------------------------
   !  Find a maximum integer value
   !---------------------------------------------------------------------
   function msg_imax(v) result(vmax)
      use mpi_f08, only: MPI_ALLREDUCE, MPI_INTEGER, MPI_MAX
      use GLOBAL, only: setMessageSource, allMessage, unsetMessageSource, COMM
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      use GLOBAL, only: DEBUG
#endif
      implicit none
      integer, intent(in) :: v
      integer :: kount
      integer :: vmax

      call setMessageSource("psdot")
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      kount = 1
      call MPI_ALLREDUCE(v, vmax, kount, MPI_INTEGER, &
                         MPI_MAX, COMM, ierr)

#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
      return
      !---------------------------------------------------------------------
   end function msg_imax
   !---------------------------------------------------------------------

   !---------------------------------------------------------------------
   !     F U N C T I O N   P S D O T
   !---------------------------------------------------------------------
   !  Parallel version of SDOT for ITPACKV module
   !---------------------------------------------------------------------
   real(8) function psdot(n, sx, sy) result(gsum)
      use mpi_f08, only: MPI_ALLREDUCE, MPI_SUM, MPI_DOUBLE_PRECISION
      use GLOBAL, only: setMessageSource, allMessage, unsetMessageSource, COMM
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      use GLOBAL, only: DEBUG
#endif
      implicit none
      integer, intent(in) :: n
      real(8), intent(in) :: sx(n), sy(n)
      integer :: i
      real(8) :: lsum
      integer :: kount !jgf46.00 added

      call setMessageSource("psdot")
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      gsum = 0.0d0
      if (n <= 0) then
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG, "Return.")
#endif
         call unsetMessageSource()
         return
      end if

      lsum = 0.0d0
      do i = 1, n
         if (.not. RESNODE(i)) cycle
         lsum = lsum + sx(i)*sy(i)
      end do

      kount = 1
      call MPI_ALLREDUCE(lsum, gsum, kount, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, COMM, ierr)

#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
      return
      !---------------------------------------------------------------------
   end function psdot
   !---------------------------------------------------------------------

   !---------------------------------------------------------------------
   !     S U B R O U T I N E   P S 2 D O T S
   !---------------------------------------------------------------------
   !  Parallel version of 2-SDOTs for ITPACKV module
   !  jbr  6/17/00
   !  zc   5/01/12 performance enhancements
   !---------------------------------------------------------------------
   subroutine ps2dots(n, sd, sdt, dot3rray)
      use mpi_f08, only: MPI_ALLREDUCE, MPI_SUM, MPI_DOUBLE_PRECISION
      use GLOBAL, only: setMessageSource, allMessage, unsetMessageSource, COMM
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      use GLOBAL, only: DEBUG
#endif
      implicit none
      integer, intent(in) :: n
      real(8), intent(in) :: sd(n), sdt(n)
      real(8), intent(out) :: dot3rray(3)
      integer :: i
      real(8) :: lsum(2), gsum(2)
      integer :: kount !jgf46.00 added

      call setMessageSource("ps2dots")
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      dot3rray(1:3) = 0.0d0
      if (n <= 0) then
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG, "Return.")
#endif
         call unsetMessageSource()
         return
      end if

      lsum(1:2) = 0.0d0

      do i = 1, n
         if (.not. RESNODE(I)) cycle
         lsum(1) = lsum(1) + sd(i)*sd(i)
         lsum(2) = lsum(2) + sd(i)*sdt(i)
      end do

      kount = 2
      call MPI_ALLREDUCE(lsum, gsum, kount, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, COMM, ierr)

      dot3rray(1:2) = gsum(1:2)
      dot3rray(3) = 1.0d0
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
      return
      !---------------------------------------------------------------------
   end subroutine ps2dots
   !---------------------------------------------------------------------

   !---------------------------------------------------------------------
   !     S U B R O U T I N E   P S 3 D O T S
   !---------------------------------------------------------------------
   !     Parallel version of 3-SDOTs for ITPACKV module
   !     jbr  6/17/00
   !     zc   5/01/12 performance enhancements
   !---------------------------------------------------------------------
   subroutine ps3dots(n, sd, sdt, su, dot3rray)
      use mpi_f08, only: MPI_ALLREDUCE, MPI_SUM, MPI_DOUBLE_PRECISION
      use GLOBAL, only: setMessageSource, allMessage, unsetMessageSource, COMM
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      use GLOBAL, only: DEBUG
#endif
      implicit none
      integer, intent(in) :: n
      real(8), intent(in) :: sd(n), sdt(n), su(n)
      real(8), intent(out) :: dot3rray(3)
      real(8) :: gsum(3)
      integer :: i
      real(8) :: lsum(3)
      integer :: kount ! jgf46.00 added

      call setMessageSource("ps3dots")
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      dot3rray(1:3) = 0d0
      if (n <= 0) then
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG, "Return.")
#endif
         call unsetMessageSource()
         return
      end if

      lsum(1:3) = 0d0

      do i = 1, n
         if (.not. RESNODE(i)) cycle
         lsum(1) = lsum(1) + sd(i)*sd(i)
         lsum(2) = lsum(2) + sd(i)*sdt(i)
         lsum(3) = lsum(3) + su(i)*su(i)
      end do

      kount = 3
      call MPI_ALLREDUCE(lsum, gsum, kount, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, COMM, ierr)
      dot3rray(1:3) = gsum(1:3)

#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
      return
      !---------------------------------------------------------------------
   end subroutine ps3dots
   !---------------------------------------------------------------------

   !---------------------------------------------------------------------
   !     S U B R O U T I N E   A L L N O D E S
   !---------------------------------------------------------------------
   !  Compute Number of nodes in entire domain.
   !---------------------------------------------------------------------
   subroutine ALLNODES(TOTNODES)
      use mpi_f08, only: MPI_Allreduce, MPI_INTEGER, MPI_SUM
      use SIZES, only: MNP
      use GLOBAL, only: setMessageSource, allMessage, unsetMessageSource, COMM
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      use GLOBAL, only: DEBUG
#endif
      implicit none
      integer, intent(out) :: TOTNODES
      integer :: I
      integer :: LNODES
      integer :: kount !jgf46.00 added

      call setMessageSource("allnodes")
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      LNODES = 0
      do I = 1, MNP
         if (RESNODE(I)) LNODES = LNODES + 1
      end do

      kount = 1
      call MPI_ALLREDUCE(LNODES, TOTNODES, kount, MPI_INTEGER, MPI_SUM, &
                         COMM, IERR)

#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
      !---------------------------------------------------------------------
   end subroutine ALLNODES
   !---------------------------------------------------------------------

   !---------------------------------------------------------------------
   subroutine mapToSubdomainRealMPI(np_global, np_local, local_array, &
                                    mapping, global_array)
      !---------------------------------------------------------------------
      ! Mapping to a subdomain for parallel code to avoid all processors
      ! reading the same file at the same time
      !---------------------------------------------------------------------
      use mpi_f08, only: MPI_Bcast, MPI_DOUBLE_PRECISION
      use SIZES, only: MYPROC
      use GLOBAL, only: setMessageSource, allMessage, unsetMessageSource, COMM
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      use GLOBAL, only: DEBUG
#endif
      implicit none
      real(8), intent(in) :: global_array(:)
      integer, intent(in) :: mapping(:)
      integer, intent(in) :: np_global
      integer, intent(in) :: np_local
      real(8) :: global_array_local(np_global)
      real(8), intent(out) :: local_array(:)
      integer :: sd_node_number
      integer :: ierr
      call setMessageSource("mapToSubdomainRealMPI")
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif
      if (myproc == 0) then
         global_array_local(1:np_global) = global_array(1:np_global)
      end if
      call mpi_bcast(global_array_local, np_global, MPI_DOUBLE_PRECISION, 0, &
                     COMM, ierr)
      do sd_node_number = 1, np_local
         local_array(sd_node_number) = &
            global_array_local(abs(mapping(sd_node_number)))
      end do
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
      !---------------------------------------------------------------------
   end subroutine mapToSubdomainRealMPI
   !---------------------------------------------------------------------

   !---------------------------------------------------------------------
   subroutine mapToSubdomainIntMPI(np_global, np_local, local_array, &
                                   global_array, mapping)
      !---------------------------------------------------------------------
      ! Mapping to a subdomain for parallel code to avoid all processors
      ! reading the same file at the same time
      !---------------------------------------------------------------------
      use mpi_f08, only: MPI_Bcast, MPI_INTEGER
      use SIZES, only: MYPROC
      use GLOBAL, only: setMessageSource, allMessage, unsetMessageSource, COMM
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      use GLOBAL, only: DEBUG
#endif
      implicit none
      integer, intent(in) :: global_array(:)
      integer, intent(in) :: mapping(:)
      integer, intent(in) :: np_global
      integer, intent(in) :: np_local
      integer :: global_array_local(np_global)
      integer, intent(out) :: local_array(:)
      integer :: sd_node_number
      integer :: ierr
      call setMessageSource("mapToSubdomainIntMPI")
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif
      if (myproc == 0) then
         global_array_local(1:np_global) = global_array(1:np_global)
      end if
      call mpi_bcast(global_array_local, np_global, mpi_integer, 0, &
                     COMM, ierr)
      do sd_node_number = 1, np_local
         local_array(sd_node_number) = &
            global_array_local(abs(mapping(sd_node_number)))
      end do
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
      !---------------------------------------------------------------------
   end subroutine mapToSubdomainIntMPI
   !---------------------------------------------------------------------

   !------------------------------------------------------------------------------
   !               S U B R O U T I N E   W E T D R Y S U M
   !------------------------------------------------------------------------------
   !     jgf45.06 The GWCE left hand side must be reset whenever wetting or
   !     drying has caused the grid to change. In parallel execution,
   !     wetting or drying may occur in one subdomain but not another
   !     during any particular time step. In this case, the NCChange flag
   !     will be 1 in one subdomain but 0 in another, causing the
   !     subdomains to get out of sync with each other's MPI calls. PADCIRC
   !     will necessarily hang under these circumstances. This subroutine
   !     sums the NCChange flags from all subdomains. If the sum comes back
   !     greater than 0, they must all recompute the GWCE lhs, thus
   !     preventing desynchronization of MPI communications.
   !------------------------------------------------------------------------------
   !
   subroutine WetDrySum(NCCHANGE)
      use mpi_f08, only: MPI_ALLREDUCE, MPI_INTEGER, MPI_SUM
      use GLOBAL, only: setMessageSource, allMessage, unsetMessageSource, COMM
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      use GLOBAL, only: DEBUG
#endif
      implicit none
      integer, intent(inout) :: NCCHANGE !input flag,=1 if this subdomain has wetted or dried
      integer :: SumNCChange !sum total of all flags from all subdomains
      integer :: kount !jgf46.00 to avoid compiler bug on certain platforms

      call setMessageSource("wetdrysum")
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      SumNCChange = 0
      kount = 1
      call MPI_ALLREDUCE(NCCHANGE, SumNCChange, kount, MPI_INTEGER, &
                         MPI_SUM, COMM, ierr)
      NCCHANGE = SumNCChange !resets GWCE for all subdomains if any s.d. resets

#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
      return
      !---------------------------------------------------------------------
   end subroutine WetDrySum
   !---------------------------------------------------------------------

   !------------------------------------------------------------------------------
   !               S U B R O U T I N E   W A R N  E L E V  S U M
   !------------------------------------------------------------------------------
   !     jgf46.11 If the warning elevation was exceeded in one subdomain,
   !     that information is propagated to the other subdomains so that the
   !     velocities can be dumped to a file for debugging.
   !------------------------------------------------------------------------------
   subroutine WarnElevSum(WarnElevExceeded)
      use mpi_f08, only: MPI_ALLREDUCE, MPI_INTEGER, MPI_SUM
      use GLOBAL, only: setMessageSource, allMessage, unsetMessageSource, COMM
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      use GLOBAL, only: DEBUG
#endif
      implicit none
      integer, intent(inout) :: WarnElevExceeded !=1 if this subdomain has exceeded warning elev
      integer :: SumWarnElevExceeded !sum total of all flags from all subdomains
      integer :: kount ! to avoid compiler bug on certain platforms

      call setMessageSource("warnelevsum")
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      SumWarnElevExceeded = 0
      kount = 1
      call MPI_ALLREDUCE(WarnElevExceeded, SumWarnElevExceeded, kount, &
                         MPI_INTEGER, MPI_SUM, COMM, ierr)
      WarnElevExceeded = SumWarnElevExceeded

#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
      !------------------------------------------------------------------------------
   end subroutine WarnElevSum
   !------------------------------------------------------------------------------

   !------------------------------------------------------------------------------
   !               S U B R O U T I N E   E A R L Y T E R M S U M
   !------------------------------------------------------------------------------
   !     DMW 202401 If a NaN was detected in the velocity soln in one subdomain,
   !     that information is propagated to the other subdomains so that the
   !     solns can be output for debugging before exiting the model.
   !------------------------------------------------------------------------------
   subroutine EarlyTermSum(earlyterminate)
      use mpi_f08, only: MPI_ALLREDUCE, MPI_INTEGER, MPI_SUM
      use GLOBAL, only: setMessageSource, allMessage, unsetMessageSource, COMM
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      use GLOBAL, only: DEBUG
#endif
      implicit none
      integer, intent(inout) :: earlyterminate !=1 if this subdomain has detected a NaN in the velocity soln and 20 timesteps have passed since then
      integer :: SumEarlyterminate !sum total of all flags from all subdomains
      integer :: kount ! to avoid compiler bug on certain platforms

      call setMessageSource("earlytermsum")
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      SumEarlyterminate = 0
      kount = 1
      call MPI_ALLREDUCE(earlyterminate, SumEarlyterminate, kount, &
                         MPI_INTEGER, MPI_SUM, COMM, ierr)
      earlyterminate = SumEarlyterminate

#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
      !------------------------------------------------------------------------------
   end subroutine EarlyTermSum
   !------------------------------------------------------------------------------

   !---------------------------------------------------------------------
   !     S U B R O U T I N E   M S G _ I B C A S T
   !---------------------------------------------------------------------
   !  Broadcast integer array from processor 0.
   !  vjp  9/26/2006
   !---------------------------------------------------------------------
   subroutine MSG_IBCAST(array, n)
      use mpi_f08, only: MPI_BCAST, MPI_INTEGER
      use GLOBAL, only: setMessageSource, allMessage, unsetMessageSource, COMM
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      use GLOBAL, only: DEBUG
#endif
      implicit none
      integer, intent(in) :: n
      integer, intent(inout) :: array(n)

      call setMessageSource("msg_ibcast")
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      call MPI_BCAST(array, n, mpi_integer, 0, comm, ierr)

#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
      return
      !---------------------------------------------------------------------
   end subroutine MSG_IBCAST
   !---------------------------------------------------------------------

   !---------------------------------------------------------------------
   !     S U B R O U T I N E   M S G _ L B C A S T
   !---------------------------------------------------------------------
   !---------------------------------------------------------------------
   subroutine MSG_LBCAST(array, n)
      use mpi_f08, only: MPI_BCAST, MPI_LOGICAL
      use GLOBAL, only: setMessageSource, allMessage, unsetMessageSource, COMM
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      use GLOBAL, only: DEBUG
#endif
      implicit none
      integer, intent(in) :: n
      logical, intent(inout) :: array(n)

      call setMessageSource("msg_lbcast")
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      call MPI_BCAST(array, n, mpi_logical, 0, comm, ierr)

#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
      return
      !---------------------------------------------------------------------
   end subroutine MSG_LBCAST
   !---------------------------------------------------------------------

   !---------------------------------------------------------------------
   !     S U B R O U T I N E   M S G _ C B C A S T
   !---------------------------------------------------------------------
   !  Broadcast integer array from processor 0.
   !  vjp  9/26/2006
   !---------------------------------------------------------------------
   subroutine MSG_CBCAST(msg, n)
      use mpi_f08, only: MPI_BCAST, MPI_CHARACTER
      use GLOBAL, only: setMessageSource, allMessage, unsetMessageSource, COMM
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      use GLOBAL, only: DEBUG
#endif
      implicit none
      integer, intent(in) :: n
      character(n), intent(inout) :: msg

      call setMessageSource("msg_cbcast")
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      call MPI_BCAST(msg, n, mpi_character, 0, comm, ierr)

#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
      return
      !---------------------------------------------------------------------
   end subroutine MSG_CBCAST
   !---------------------------------------------------------------------

   !---------------------------------------------------------------------
   !     S U B R O U T I N E   M S G _ R B C A S T
   !---------------------------------------------------------------------
   !  Broadcast integer array from processor 0.
   !  vjp  9/26/2006
   !---------------------------------------------------------------------
   subroutine MSG_RBCAST(array, n)
      use mpi_f08, only: MPI_BCAST, MPI_DOUBLE_PRECISION
      use GLOBAL, only: setMessageSource, allMessage, unsetMessageSource, COMM
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      use GLOBAL, only: DEBUG
#endif
      implicit none
      integer, intent(in) :: n
      real(8), intent(inout) :: array(:, :, :)

      call setMessageSource("msg_rbcast")
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      call MPI_BCAST(array, n, MPI_DOUBLE_PRECISION, 0, comm, ierr)

#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
      !---------------------------------------------------------------------
   end subroutine MSG_RBCAST
   !---------------------------------------------------------------------

   !---------------------------------------------------------------------
   !     S U B R O U T I N E   M S G _ R B C A S T D
   !---------------------------------------------------------------------
   !asey 090327: Implement Seizo's buffering of radiation stress broadcasts.
   !=====Seizo Dividing Messenger
   !---------------------------------------------------------------------
   subroutine MSG_RBCASTD(array, in, jn, kn)
      use mpi_f08, only: MPI_BCAST, MPI_DOUBLE_PRECISION
      use GLOBAL, only: setMessageSource, allMessage, unsetMessageSource, COMM
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      use GLOBAL, only: DEBUG
#endif
      implicit none
      integer, intent(IN) :: in, jn, kn
      real(8), intent(INOUT) :: array(in, jn, kn)
      integer :: i, j, k, icount, ntimes, nlast
      real(8), allocatable :: buffer(:)
      integer, allocatable :: nbox(:)
      integer, parameter :: limit_buff = 1000000

      call setMessageSource("msg_rbcastd")
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      allocate (buffer(in*jn*kn))

      ntimes = in*jn*kn/limit_buff
      nlast = in*jn*kn - limit_buff*ntimes
      allocate (nbox(ntimes + 1))
      do i = 1, ntimes
         nbox(i) = limit_buff
      end do
      if (nlast /= 0) then
         ntimes = ntimes + 1
         nbox(ntimes) = nlast
      end if

      icount = 0
      do i = 1, in
         do j = 1, jn
            do k = 1, kn
               icount = icount + 1
               buffer(icount) = array(i, j, k)
            end do
         end do
      end do

      do i = 1, ntimes
         icount = limit_buff*(i - 1) + 1
         call MPI_BCAST(buffer(icount), nbox(i), MPI_DOUBLE_PRECISION, 0, comm, ierr)
      end do

      icount = 0
      do i = 1, in
         do j = 1, jn
            do k = 1, kn
               icount = icount + 1
               array(i, j, k) = buffer(icount)
            end do
         end do
      end do
      deallocate (buffer, nbox)

#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
      return
      !---------------------------------------------------------------------
   end subroutine MSG_RBCASTD
   !---------------------------------------------------------------------

   !---------------------------------------------------------------------
   !---------------------------------------------------------------------
   !---------------------------------------------------------------------
   subroutine MSG_BARRIER()
      use mpi_f08, only: MPI_Barrier
      use GLOBAL, only: COMM
      implicit none
      call MPI_Barrier(comm)
   end subroutine MSG_BARRIER
   !---------------------------------------------------------------------
   !=============================================================================

   !=============================================================================
   !---------------------------------------------------------------------
   ! DW, 2024
   subroutine MSG_RScalar_Reduce(OP, vali, valo, loco)
      use mpi_f08, only: MPI_REDUCE, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, &
                         MPI_MAX, MPI_MINLOC, MPI_MIN, MPI_SUM, MPI_DOUBLE_PRECISION
      use GLOBAL, only: setMessageSource, allMessage, unsetMessageSource, COMM
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      use GLOBAL, only: DEBUG
#endif
      implicit none

      character(LEN=*), intent(in) :: OP
      real(8), intent(in) :: vali
      real(8), intent(out) :: valo
      real(8), intent(inout), optional :: loco !/ in - rank, out - max/min rank

      real(8) :: tempi(2), tempo(2)

      call setMessageSource("msg_rvecreduce")
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      select case (OP)
      case ('max', 'MAX')
         if (present(loco)) then
            tempi = [vali, loco]
            call MPI_Reduce(tempi, tempo, 1, MPI_2DOUBLE_PRECISION, &
                            MPI_MAXLOC, 0, comm, ierr)
            valo = tempo(1)
            loco = tempo(2)
         else
            call MPI_Reduce(vali, valo, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, &
                            comm, ierr)
         end if
      case ('min', 'MIN')
         if (present(loco)) then
            tempi = [vali, loco]
            call MPI_Reduce(tempi, tempo, 1, MPI_2DOUBLE_PRECISION, &
                            MPI_MINLOC, 0, comm, ierr)
            valo = tempo(1)
            loco = tempo(2)
         else
            call MPI_Reduce(vali, valo, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, &
                            comm, ierr)
         end if
      case ('sum', 'SUM')
         call MPI_Reduce(vali, valo, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
                         comm, ierr)
      end select

#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
   end subroutine MSG_RScalar_Reduce
   !---------------------------------------------------------------------
   !=============================================================================

   subroutine MSG_ABORT()
      use mpi_f08, only: MPI_Abort
      use global, only: comm
      implicit none
      integer :: myerrcode
      call MPI_Abort(comm, myerrcode)
   end subroutine MSG_ABORT

end module MESSENGER
