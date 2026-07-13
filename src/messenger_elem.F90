module MESSENGER_ELEM

   use SIZES, only: myproc, mne, inputdir, mnproc
   use GLOBAL, only: C3D
   use mpi_f08, only: MPI_COMM_SIZE, MPI_Comm, MPI_Group, MPI_Request, MPI_Status, &
                      MPI_Comm_World, MPI_Integer, MPI_Sum, Mpi_Double_Precision, MPI_Waitsome, &
                      MPI_Send_Init, MPI_Recv_init, MPI_startall, MPI_Abort, MPI_Finalize, &
                      MPI_Comm_Create, MPI_Group_Incl, MPI_Comm_Group, MPI_Comm_rank, &
                      MPI_Allreduce, MPI_start
   implicit none

   private

!
!--------------------------------------------------------------------------
!  This module supplies the MPI Message-Passing Interface
!
!  Uses asynchronous and persistent communication with buffer packing
!  as performance enhancements for "cluster" architectures.
!
!  vjp  8/29/1999
!--------------------------------------------------------------------------
!
!

!
!  Message-Passing Array space
!
!sb-PDG1

   integer, parameter, private :: sz = 8, dofh = 3
   integer, private ::  NEIGHPROC_R, NEIGHPROC_S, RDIM, IERR
   integer, private ::  TAG = 101
   logical, allocatable :: RESELEM(:)
!
   integer, private, allocatable ::IPROC_R(:), IPROC_S(:), NELEMLOC(:), &
                                    NELEMSEND(:), NELEMRECV(:), ISENDLOC(:, :), IBELONGTO(:), &
                                    IRECVLOC(:, :), ISENDBUF(:, :), IRECVBUF(:, :)
!
   type(MPI_Request), private, allocatable :: REQ_I1(:), REQ_I2(:)
   type(MPI_Status), private, allocatable :: STAT_I1(:), STAT_I2(:)
   type(MPI_Request), private, allocatable :: REQ_R1(:), REQ_R2(:), REQ_R3(:), &
                                              REQ_LZ(:)
   type(MPI_Status), private, allocatable :: STAT_R1(:), STAT_R2(:), &
                                             STAT_R3(:), STAT_LZ(:)
   integer, private, allocatable :: index(:)
   real(8), private, allocatable :: SENDBUF(:, :), RECVBUF(:, :)
!--
!

   public :: updater_elem_mod, message_start_elem, msg_table_elem
!---------------------end of data declarations--------------------------------C

contains

   subroutine MSG_TABLE_ELEM()
!
!--------------------------------------------------------------------------
!  Routine preforms following steps:
!
!   (1) Read Message-Passing Information from file "DG.18"
!   (2) Determine resident nodes: RESNODE(I) is true  if I is resident node
!   (3) Determine ghost nodes:    RESNODE(I) is false if I is ghost node
!   (4) Determine number of neighbor subdomains
!   (5) MPI rank of each neighbor and number of ghosts nodes to receive
!   (6) Read Message-Passing Receive List
!   (7) MPI rank of each neighbor and number of ghosts nodes to send
!   (8) Read Message-Passing Send List
!
!  vjp  8/06/1999
!--------------------------------------------------------------------------
!
      implicit none
      integer :: IDPROC, NLOCAL, I, J, JJ
      open (18, FILE=trim(inputdir)//'/'//'DG.18')
      read (18, 3010) IDPROC, NLOCAL
      allocate (NELEMLOC(NLOCAL))
      read (18, 1130) (NELEMLOC(I), I=1, NLOCAL)
      allocate (IBELONGTO(MNE), RESELEM(MNE))
      do I = 1, MNE
         IBELONGTO(I) = 0
      end do
      do I = 1, NLOCAL
         IBELONGTO(NELEMLOC(I)) = IDPROC + 1
      end do
      do I = 1, MNE
         if (IBELONGTO(I) - 1 == MYPROC) then
            RESELEM(I) = .true.
         else
            RESELEM(I) = .false.
         end if
      end do
      read (18, 3010) NEIGHPROC_R, NEIGHPROC_S
      RDIM = NEIGHPROC_R + NEIGHPROC_S
      allocate (index(RDIM))
      allocate (IPROC_R(NEIGHPROC_R), NELEMRECV(NEIGHPROC_R))
      allocate (IRECVLOC(MNE, NEIGHPROC_R))
      do JJ = 1, NEIGHPROC_R
         J = mod(JJ - 1 + MYPROC, NEIGHPROC_R) + 1
         read (18, 3010) IPROC_R(J), NELEMRECV(J)
         read (18, 1130) (IRECVLOC(I, J), I=1, NELEMRECV(J))
      end do
!
      allocate (IPROC_S(NEIGHPROC_S), NELEMSEND(NEIGHPROC_S))
      allocate (ISENDLOC(MNE, NEIGHPROC_S))
!
      do JJ = 1, NEIGHPROC_S
         J = mod(JJ - 1 + MYPROC, NEIGHPROC_S) + 1
         read (18, 3010) IPROC_S(J), NELEMSEND(J)
         read (18, 1130) (ISENDLOC(I, J), I=1, NELEMSEND(J))
      end do
!
      close (18)
      return
!
1130  format(8x, 9i8)
3010  format(8x, 2i8)
   end subroutine MSG_TABLE_ELEM

   subroutine MESSAGE_START_ELEM()
!
!--------------------------------------------------------------------------
!  Routine preforms following steps:
!   (1)  allocate message-passing space
!   (2)  setup MPI data structures for "persistent" message-passing.
!
!  vjp  8/06/1999
!--------------------------------------------------------------------------
!
      implicit none
      integer :: J, NCOMMELEM_S, NCOMMELEM_R
!
      NCOMMELEM_R = 0
      NCOMMELEM_S = 0
      do J = 1, NEIGHPROC_R
         NCOMMELEM_R = NCOMMELEM_R + NELEMRECV(J)
      end do
      do J = 1, NEIGHPROC_S
         NCOMMELEM_S = NCOMMELEM_S + NELEMSEND(J)
      end do
!
      allocate (ISENDBUF(NCOMMELEM_S*2, NEIGHPROC_S))
      allocate (IRECVBUF(NCOMMELEM_R*2, NEIGHPROC_R))
!
      if (C3D) then
!         ALLOCATE ( SENDBUF(2*MNP*MNODES,NEIGHPROC) )
!         ALLOCATE ( RECVBUF(2*MNP*MNODES,NEIGHPROC) )
         stop 'NOT UPDATED'
      else
         allocate (SENDBUF(NCOMMELEM_S*DOFH*4, NEIGHPROC_S))
         allocate (RECVBUF(NCOMMELEM_R*DOFH*4, NEIGHPROC_R))
      end if
!
      allocate (REQ_I1(RDIM), REQ_I2(RDIM))
      allocate (REQ_R1(RDIM), REQ_R2(RDIM), REQ_R3(RDIM), REQ_LZ(RDIM))
!
      allocate (STAT_I1(RDIM), &
                STAT_I2(RDIM))

      allocate (STAT_R1(RDIM), &
                STAT_R2(RDIM), &
                STAT_R3(RDIM), &
                STAT_LZ(RDIM))
!
      if (C3D) then
!         ALLOCATE ( REQ_R3D(RDIM) )
!         ALLOCATE ( STAT_R3D(RDIM) )
!         ALLOCATE ( REQ_C3D(RDIM) )
!         ALLOCATE ( STAT_C3D(RDIM) )
         stop 'NOT UPDATED'
      end if
!
      !  Setup persistent structures for integer arrays
!
      do J = 1, NEIGHPROC_R
         call MPI_RECV_INIT(IRECVBUF(1, J), NELEMRECV(J), &
                            MPI_INTEGER, IPROC_R(J), TAG, MPI_COMM_WORLD, &
                            REQ_I1(J), IERR)
      end do
      do J = 1, NEIGHPROC_S
         call MPI_SEND_INIT(ISENDBUF(1, J), NELEMSEND(J), &
                            MPI_INTEGER, IPROC_S(J), TAG, MPI_COMM_WORLD, &
                            REQ_I1(J + NEIGHPROC_R), IERR)
      end do
!
!
      do J = 1, NEIGHPROC_R
         call MPI_RECV_INIT(IRECVBUF(1, J), 2*NELEMRECV(J), &
                            MPI_INTEGER, IPROC_R(J), TAG, MPI_COMM_WORLD, &
                            REQ_I2(J), IERR)
      end do
      do J = 1, NEIGHPROC_S
         call MPI_SEND_INIT(ISENDBUF(1, J), 2*NELEMSEND(J), &
                            MPI_INTEGER, IPROC_S(J), TAG, MPI_COMM_WORLD, &
                            REQ_I2(J + NEIGHPROC_R), IERR)
      end do
!
      !  Setup persistent structures for real arrays
!
      do J = 1, NEIGHPROC_R
         call MPI_RECV_INIT(RECVBUF(1, J), DOFH*NELEMRECV(J), &
                            MPI_DOUBLE_PRECISION, IPROC_R(J), TAG, MPI_COMM_WORLD, &
                            REQ_R1(J), IERR)
      end do
      do J = 1, NEIGHPROC_S
         call MPI_SEND_INIT(SENDBUF(1, J), DOFH*NELEMSEND(J), &
                            MPI_DOUBLE_PRECISION, IPROC_S(J), TAG, MPI_COMM_WORLD, &
                            REQ_R1(J + NEIGHPROC_R), IERR)
      end do
!
      do J = 1, NEIGHPROC_R
         call MPI_RECV_INIT(RECVBUF(1, J), 2*DOFH*NELEMRECV(J), &
                            MPI_DOUBLE_PRECISION, IPROC_R(J), TAG, MPI_COMM_WORLD, &
                            REQ_R2(J), IERR)
      end do
      do J = 1, NEIGHPROC_S
         call MPI_SEND_INIT(SENDBUF(1, J), 2*DOFH*NELEMSEND(J), &
                            MPI_DOUBLE_PRECISION, IPROC_S(J), TAG, MPI_COMM_WORLD, &
                            REQ_R2(J + NEIGHPROC_R), IERR)
      end do
!
      do J = 1, NEIGHPROC_R
         call MPI_RECV_INIT(RECVBUF(1, J), 3*DOFH*NELEMRECV(J), &
                            MPI_DOUBLE_PRECISION, IPROC_R(J), TAG, MPI_COMM_WORLD, &
                            REQ_R3(J), IERR)
      end do
      do J = 1, NEIGHPROC_S
         call MPI_SEND_INIT(SENDBUF(1, J), 3*DOFH*NELEMSEND(J), &
                            MPI_DOUBLE_PRECISION, IPROC_S(J), TAG, MPI_COMM_WORLD, &
                            REQ_R3(J + NEIGHPROC_R), IERR)
      end do
!
      do J = 1, NEIGHPROC_R
         call MPI_RECV_INIT(RECVBUF(1, J), 4*DOFH*NELEMRECV(J), &
                            MPI_DOUBLE_PRECISION, IPROC_R(J), TAG, MPI_COMM_WORLD, &
                            REQ_LZ(J), IERR)
      end do
      do J = 1, NEIGHPROC_S
         call MPI_SEND_INIT(SENDBUF(1, J), 4*DOFH*NELEMSEND(J), &
                            MPI_DOUBLE_PRECISION, IPROC_S(J), TAG, MPI_COMM_WORLD, &
                            REQ_LZ(J + NEIGHPROC_R), IERR)
      end do
      if (C3D) then
         stop 'NOT UPDATED'
      end if

   end subroutine MESSAGE_START_ELEM

   subroutine UPDATER_ELEM_MOD(VEC1, VEC2, VEC3, IRK, NMSG)
!
!--------------------------------------------------------------------------
!  Update 1, 2, or 3 Real Arrays's Ghost Cells using asynchronous
!  and persistent message-passing.
!
!  vjp  8/06/1999
!--------------------------------------------------------------------------
!
      implicit none
!
      integer, intent(IN) ::  IRK, NMSG
      real(SZ), intent(INOUT) ::  VEC1(:, :, :), VEC2(:, :, :), VEC3(:, :, :)
      integer :: N, I, J, K, NCOUNT, NFINI, TOT
!
      !..Pack 1, 2, or 3 Messages
!      DO JJ=1,NEIGHPROC_S
!         J = MOD(JJ-1+MYPROC,NEIGHPROC_S)+1
      do J = 1, NEIGHPROC_S
         NCOUNT = 0
         do I = 1, NELEMSEND(J)
            do K = 1, DOFH
               NCOUNT = NCOUNT + 1
               SENDBUF(NCOUNT, J) = VEC1(K, ISENDLOC(I, J), IRK)
            end do
         end do
         if (NMSG > 1) then
            do I = 1, NELEMSEND(J)
               do K = 1, DOFH
                  NCOUNT = NCOUNT + 1
                  SENDBUF(NCOUNT, J) = VEC2(K, ISENDLOC(I, J), IRK)
               end do
            end do
         end if
         if (NMSG > 2) then
            do I = 1, NELEMSEND(J)
               do K = 1, DOFH
                  NCOUNT = NCOUNT + 1
                  SENDBUF(NCOUNT, J) = VEC3(K, ISENDLOC(I, J), IRK)
               end do
            end do
         end if

         ! Start sending a message
         if (NMSG == 1) then
            call MPI_START(REQ_R1(J + NEIGHPROC_R), IERR)
         elseif (NMSG == 2) then
            call MPI_START(REQ_R2(J + NEIGHPROC_R), IERR)
         else
            call MPI_START(REQ_R3(J + NEIGHPROC_R), IERR)
         end if
      end do
!
      ! Send/receive messages to/from all neighbors
!
      if (NMSG == 1) then
         call MPI_STARTALL(RDIM - NEIGHPROC_S, REQ_R1, IERR)
      elseif (NMSG == 2) then
         call MPI_STARTALL(RDIM - NEIGHPROC_S, REQ_R2, IERR)
      else
         call MPI_STARTALL(RDIM - NEIGHPROC_S, REQ_R3, IERR)
      end if
!
      !..Unpack Received messages as they arrive
!
      if (NMSG == 1) then
         TOT = 0
         do while (TOT < RDIM)
            do N = 1, RDIM
               index(N) = 0
            end do
            call MPI_WAITSOME(RDIM, REQ_R1, NFINI, INDEX, STAT_R1, IERR)
            TOT = TOT + NFINI
            do N = 1, NFINI
               if (index(N) > 0 .and. index(N) <= RDIM) then
                  if (index(N) <= NEIGHPROC_R) then
                     J = index(N)
                     NCOUNT = 0
                     do I = 1, NELEMRECV(J)
                        do K = 1, DOFH
                           NCOUNT = NCOUNT + 1
                           VEC1(K, IRECVLOC(I, J), IRK) = RECVBUF(NCOUNT, J)
                        end do
                     end do
                  end if
               end if
            end do
         end do
         goto 999
      elseif (NMSG == 2) then
         TOT = 0
         do while (TOT < RDIM)
            do N = 1, RDIM
               index(N) = 0
            end do
            call MPI_WAITSOME(RDIM, REQ_R2, NFINI, INDEX, STAT_R2, IERR)
            TOT = TOT + NFINI
            do N = 1, NFINI
               if (index(N) > 0 .and. index(N) <= RDIM) then
                  if (index(N) <= NEIGHPROC_R) then
                     J = index(N)
                     NCOUNT = 0
                     do I = 1, NELEMRECV(J)
                        do K = 1, DOFH
                           NCOUNT = NCOUNT + 1
                           VEC1(K, IRECVLOC(I, J), IRK) = RECVBUF(NCOUNT, J)
                        end do
                     end do
                     do I = 1, NELEMRECV(J)
                        do K = 1, DOFH
                           NCOUNT = NCOUNT + 1
                           VEC2(K, IRECVLOC(I, J), IRK) = RECVBUF(NCOUNT, J)
                        end do
                     end do
                  end if
               end if
            end do
         end do
         goto 999
      else
         TOT = 0
         do while (TOT < RDIM)
            do N = 1, RDIM
               index(N) = 0
            end do
            call MPI_WAITSOME(RDIM, REQ_R3, NFINI, INDEX, STAT_R3, IERR)
            TOT = TOT + NFINI
!debug     print *, myproc, tot,nfini,index(1),index(2)
            do N = 1, NFINI
               if (index(N) > 0 .and. index(N) <= RDIM) then
                  if (index(N) <= NEIGHPROC_R) then
                     J = index(N)
                     NCOUNT = 0
                     do I = 1, NELEMRECV(J)
                        do K = 1, DOFH
                           NCOUNT = NCOUNT + 1
                           VEC1(K, IRECVLOC(I, J), IRK) = RECVBUF(NCOUNT, J)
                        end do
                     end do
                     do I = 1, NELEMRECV(J)
                        do K = 1, DOFH
                           NCOUNT = NCOUNT + 1
                           VEC2(K, IRECVLOC(I, J), IRK) = RECVBUF(NCOUNT, J)
                        end do
                     end do
                     do I = 1, NELEMRECV(J)
                        do K = 1, DOFH
                           NCOUNT = NCOUNT + 1
                           VEC3(K, IRECVLOC(I, J), IRK) = RECVBUF(NCOUNT, J)
                        end do
                     end do
                  end if
               end if
            end do
         end do
         goto 999
      end if
!
999   continue
      return
   end subroutine UPDATER_ELEM_MOD

end module MESSENGER_ELEM
