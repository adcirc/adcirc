      SUBROUTINE METIS()
      USE PRE_GLOBAL
      USE METIS_INTERFACE, ONLY: METIS_SETDEFAULTOPTIONS,
     &   METIS_PARTGRAPHKWAY, METIS_NOPTIONS,
     &   METIS_OPTION_NUMBERING, idx_t, METIS_OK
      IMPLICIT NONE
C----------------------------------------------------------------------
C  INTERFACE ROUTINE FOR PADCIRC TO USE THE METIS LIBRARY
C  A GRAPH PARTITION TOOL, FOR DOMAIN DECOMPOSITION
C  Version 2.1  vjp  6/7/2006
C  Redimensioned CO_NODES and IDUALS and added bounds test
C  Added check that adjacency matrix is symmetric
C  KJR, August 25, 2021 - support for metis 5.1.0
C----------------------------------------------------------------------
C
        LOGICAL FOUND, SYMMETRIC
        INTEGER MNED, MNEDLOC, IDUMP, IERR
        INTEGER I, J, K, IEL, INODE, JNODE, NCOUNT, ITOT, NEDGETOT
        INTEGER NCON
        INTEGER MAXDUALS
        INTEGER  WEIGHTFLAG, NUMFLAG, NPARTS
        INTEGER  EDGECUT, OPTYPE, NBYTES
C
        INTEGER,ALLOCATABLE :: IDUALS(:,:),ITVECT(:),ITVECT2(:)
        INTEGER,ALLOCATABLE :: XADJ(:), ADJNCY(:)
        INTEGER,ALLOCATABLE :: VWGTS(:), EWGTS(:)
        INTEGER,ALLOCATABLE :: CO_NODES(:,:),NEDGES(:), NEDLOC(:)
        INTEGER,ALLOCATABLE :: NUMDUALS(:)
        INTEGER(idx_t)      :: OPTIONS(0:METIS_NOPTIONS-1)
C...... DW
        INTEGER:: ID1, ID2
        INTEGER, ALLOCATABLE:: NNEG_TMP(:,:), PERBC_IDN_MAP(:)
C...... DW

        integer :: contiguous ! number of nodes on a boundary within a particular subdomain
        integer :: subdomainNumber ! subdomain that a boundary node falls within

        ALLOCATE ( ITVECT(MNP),ITVECT2(MNP) )
        ALLOCATE ( NUMDUALS(MNP) )
        ALLOCATE ( XADJ(MNP+1), VWGTS(MNP), NEDGES(MNP) )
        ALLOCATE ( NEDLOC(MNP) )
C
C--COMPUTE INDEX OF WEIR DUALS WHICH IS ZERO IF NOT A WEIR NODE
C
        DO INODE=1, MNP
           NUMDUALS(INODE) = 0
        ENDDO
        DO J=1, NWEIR
           NUMDUALS(WEIR(J)) = NUMDUALS(WEIR(J))+1
           NUMDUALS(WEIRD(J)) = NUMDUALS(WEIRD(J))+1
        ENDDO

        MAXDUALS = 0
        DO J=1, MNP
           IF (NUMDUALS(J) .ge. MAXDUALS) MAXDUALS = NUMDUALS(J)
        ENDDO

        write(*,'(a,i0,a)')
     &   'INFO: This mesh has ',NWEIR,' weir node pairs.'
        write(*,'(a,i0,a)')
     &  'INFO: Maximum number of duals for any weir node is ',maxduals,
     &  '.'

        ALLOCATE ( IDUALS(MAXDUALS,MNP) )


        DO INODE=1, MNP
        DO K=1, MAXDUALS
           IDUALS(K,INODE) = 0
        ENDDO
        ENDDO

        DO J=1, NWEIR
        DO K=1, MAXDUALS
           IF (IDUALS(K,WEIR(J)) == 0) THEN
             IDUALS(K,WEIR(J)) = WEIRD(J)
             EXIT
           ENDIF
        ENDDO
        ENDDO

        DO J=1, NWEIR
        DO K=1, MAXDUALS
           IF (IDUALS(K,WEIRD(J)) == 0) THEN
             IDUALS(K,WEIRD(J)) = WEIR(J)
             EXIT
           ENDIF
        ENDDO
        ENDDO
C......
C...... DW, periodic boundary conditions
C...... - a workaround approach based on a temporary connectivity table
C......   that uses the periodic node indices
        IF ( NPERBC > 0 ) THEN
           ALLOCATE( NNEG_TMP(3,MNE) ) ;
           NNEG_TMP = NNEG ; ! dummy arr holding the original element table

           ALLOCATE( PERBC_IDN_MAP(MNP) ) ; ! create mapping between the original node indice and perbc node indices
           DO J = 1,MNP
              PERBC_IDN_MAP(J) = J;
           ENDDO
           PERBC_IDN_MAP(IPERCONN(1:NPERBC,2)) = IPERCONN(1:NPERBC,1) ;

           DO I = 1, MNE
               NNEG(:,I) = PERBC_IDN_MAP( NNEG(:,I) ) ;  ! swicth to the element with the node-per bc indices !
           END DO
        END IF
C...... DW
C......

C-------------------------------------------------------------
C  COMPUTES THE TOTAL NUMBER OF EDGES        -->    MNED
C  AND THE MAX NUMBER OF EDGES FOR ANY NODE  -->    MNEDLOC
C  BOTH COUNTS INCLUDE WEIR-NODE PAIRS
C-------------------------------------------------------------

        MNED = 0
        DO INODE = 1,MNP
           NEDLOC(INODE) = 0
        ENDDO

        DO J=1, 3
           DO IEL=1, MNE
              INODE = NNEG(J,IEL)
              NCOUNT = NEDLOC(INODE) + 2
              MNED = MNED + 2
              DO K=1, MAXDUALS
                 IF (IDUALS(K,INODE).NE.0) THEN
                   NCOUNT = NCOUNT + 1
                   MNED = MNED + 1
                 ENDIF
              ENDDO
              NEDLOC(INODE) = NCOUNT
           ENDDO
        ENDDO

        MNEDLOC = 0
        DO INODE=1, MNP
           IF (NEDLOC(INODE) .ge. MNEDLOC) MNEDLOC = NEDLOC(INODE)
        ENDDO

c       print *, "total number of edges = ", MNED
        print *, "maximum co-nodes for any node = ", MNEDLOC
        print *, MNED, MNEDLOC, MNP

        ALLOCATE ( ADJNCY(MNED), EWGTS(MNED) )
        ALLOCATE ( CO_NODES(MNEDLOC,MNP) )

C
C--COMPUTE CO_NODES LISTS AND NUMBER OF EDGES CONTAINING A NODE
C
        DO INODE = 1,MNP
           NEDGES(INODE) = 0
        ENDDO
C
        DO IEL=1, MNE
           INODE = NNEG(1,IEL)
           CO_NODES(NEDGES(INODE)+1,INODE) = NNEG(2,IEL)
           CO_NODES(NEDGES(INODE)+2,INODE) = NNEG(3,IEL)
           NCOUNT = NEDGES(INODE) + 2
           DO K=1, MAXDUALS
              IF (IDUALS(K,INODE).NE.0) THEN
                NCOUNT = NCOUNT + 1
                CO_NODES(NCOUNT,INODE) = IDUALS(K,INODE)
              ENDIF
           ENDDO
           NEDGES(INODE) = NCOUNT
        ENDDO
C
        DO IEL=1, MNE
           INODE = NNEG(2,IEL)
           CO_NODES(NEDGES(INODE)+1,INODE) = NNEG(3,IEL)
           CO_NODES(NEDGES(INODE)+2,INODE) = NNEG(1,IEL)
           NCOUNT = NEDGES(INODE) + 2
           DO K=1, MAXDUALS
              IF (IDUALS(K,INODE).NE.0) THEN
                NCOUNT = NCOUNT + 1
                CO_NODES(NCOUNT,INODE) = IDUALS(K,INODE)
              ENDIF
           ENDDO
           NEDGES(INODE) = NCOUNT
        ENDDO
C
        DO IEL=1, MNE
           INODE = NNEG(3,IEL)
           CO_NODES(NEDGES(INODE)+1,INODE) = NNEG(1,IEL)
           CO_NODES(NEDGES(INODE)+2,INODE) = NNEG(2,IEL)
           NCOUNT = NEDGES(INODE) + 2
           DO K=1, MAXDUALS
              IF (IDUALS(K,INODE).NE.0) THEN
                NCOUNT = NCOUNT + 1
                CO_NODES(NCOUNT,INODE) = IDUALS(K,INODE)
              ENDIF
           ENDDO
           NEDGES(INODE) = NCOUNT
        ENDDO

C
C  REMOVE REDUNDANCY IN NODE LISTS
C
        NEDGETOT = 0           !  This will be twice number of edges
        DO INODE = 1,MNP
           DO J=1, NEDGES(INODE)
              ITVECT(J) = CO_NODES(J,INODE)
           ENDDO
           IF (NEDGES(INODE).GT.1) THEN
             NCOUNT = NEDGES(INODE)
             CALL SORT(NCOUNT,ITVECT)
             JNODE = ITVECT(1)
             CO_NODES(1,INODE) = JNODE
             NCOUNT = 1
             DO J=2, NEDGES(INODE)
                IF (ITVECT(J).NE.JNODE) THEN
                  NCOUNT = NCOUNT + 1
                  JNODE = ITVECT(J)
                  CO_NODES(NCOUNT,INODE) = JNODE
                ENDIF
             ENDDO
           ELSE
              IF ( NPERBC < 0 ) THEN
                 print *, "node = ",INODE," is isolated"
                 CALL EXIT(1)
              ELSE
C.....DW, periodic bcs allow isolated nodes if they are a periodic slave node
                 CO_NODES(1,INODE) = INODE ;
                 NCOUNT = 1 ;
                 IF ( COUNT((IPERCONN(:,2) - INODE) == 0) == 0) THEN
                    print *, "node = ",INODE," is isolated"
                    CALL EXIT(1)
                 END IF
C.....DW
              END IF
           ENDIF
           NEDGES(INODE) = NCOUNT
           NEDGETOT = NEDGETOT + NCOUNT
           if (nedges(inode) == 0) then
             print *, "inode = ", inode, " belongs to no edges"
             CALL EXIT(1)
           endif
        ENDDO
        NEDGETOT = NEDGETOT/2
        print *, "edge count = ",nedgetot

C  check that adjacency matrix is symmetric
C
      SYMMETRIC = .true.
      DO INODE = 1, MNP
      DO J = 1, NEDGES(INODE)
         JNODE = CO_NODES(J,INODE)
         FOUND = .false.
         DO K= 1, NEDGES(JNODE)
            IF (CO_NODES(K,JNODE) == INODE) THEN
              FOUND = .true.
              EXIT
            ENDIF
         ENDDO
         IF (.not. FOUND) THEN
           SYMMETRIC = .false.
      print *, "node ",inode," adjacent to ",jnode," but not visa-versa"
         ENDIF
      ENDDO
      ENDDO
      IF (.not. SYMMETRIC) THEN
         WRITE(*,'(A)') 'bad adjacency matrix: not symmetric!'
         CALL EXIT(1)
      ENDIF
C
C  COMPUTE WEIGHTS OF THE GRAPH VERTICES
C
      DO INODE = 1,MNP
         VWGTS(INODE) = NEDGES(INODE)
         if ( strictBoundaries.eqv..true. ) then
            vwgts(inode) = vwgts(inode) + boundaryWeights(inode)
         endif
      ENDDO
C
C--COMPUTE ADJACENCY LIST OF GRAPH AND ITS EDGE WEIGHTS
C
      XADJ(1) = 1
      ITOT = 0
      DO INODE = 1,MNP
         DO J = 1, NEDGES(INODE)
            ITOT = ITOT + 1
            JNODE = CO_NODES(J,INODE)
            ADJNCY(ITOT) = JNODE
C...... DW
C    Original implementation
            EWGTS(ITOT)  = (VWGTS(JNODE)+VWGTS(INODE))

C            EWGTS(ITOT) = 1 ;
C...... DW
         ENDDO
         XADJ(INODE+1) = ITOT+1
      ENDDO
C
C Dump graph to a file for debugging
C
      IDUMP = 1
      IF (IDUMP.EQ.1) THEN
         OPEN(FILE='metis_graph.txt',UNIT=99)
         WRITE(99,100) MNP, NEDGETOT, 11, 1
         DO INODE=1, MNP
            WRITE(99,200) VWGTS(INODE),
     &     (CO_NODES(J,INODE), EWGTS(XADJ(INODE)+J-1),J=1,NEDGES(INODE))
         ENDDO
         CLOSE(99)
      ENDIF

C
C--CALL K-WAY METIS FOR PARTITIONING (5.1.0)
C
      NCON     = 1 ! number of constraints
      NPARTS = MNPROC
      print *, "Grid Partition Data"

      IERR = METIS_SetDefaultOptions(options)
      if (IERR /= METIS_OK) then
        write(*,*) "METIS_SetDefaultOptions failed with error: ", IERR
        error stop 1
      end if
      options(METIS_OPTION_NUMBERING) = 1     ! Fortran-style numbering

      IERR = METIS_PartGraphKway(nvtxs=MNP,ncon=NCON,xadj=XADJ,adjncy=ADJNCY,
     &          vwgt=VWGTS, nparts=NPARTS, objval=EDGECUT,
     &          part=PROC,options=options)

      if (IERR /= METIS_OK) then
        write(*,*) "METIS_PartGraphKway failed with error: ", IERR
        error stop 1
      end if


      print *, "Total Edges Cut = ",EDGECUT
!      !
!      ! @jasonfleming: Add code to detect boundaries that have less
!      ! than three nodes in a subdomain.
!      do k=1, nbou
!         contiguous = 1
!         subdomainNumber = proc(nbvv(k,1))
!         write(*,'(a,i0,a,i0,a,i0)') 'boundary ',k,' fulldomain node ',
!     &      nbvv(k,1),' subdomain ',subdomainNumber
!         do j=2,nvell(k)
!            write(*,'(a,i0,a,i0,a,i0)') 'boundary ',k,' fulldomain node ',
!     &         nbvv(k,j),' subdomain ',proc(nbvv(k,j))
!            ! check to see if this boundary node is in the same subdomain
!            ! as the previous one on this boundary
!            if ( proc(nbvv(k,j)).eq.subdomainNumber ) then
!               ! same subdomain, increment the number of nodes on this
!               ! boundary in this subdomain
!               contiguous = contiguous + 1
!            else
!               ! this boundary node is on a different subdomain than
!              ! the previous one in this string
!               !
!               ! check to see if there were too few nodes
!               if ( contiguous.lt.3 ) then
!                  print *, 'boundary with less than 3 nodes in a subdomain:'
!                  write(*,'(a,a,i0)') trim(nvellmsg(k)),
!     &              ' in subdomain ',subdomainNumber
!                  do i=1,contiguous
!                     write(*,'(a,i0,a,i0)') 'fulldomain node number ',nbvv(k,j-contiguous),
!     &               ' in subdomain ',proc(nbvv(k,j-contiguous))
!                  end do
!               endif
!               subdomainNumber = proc(nbvv(k,j))
!               contiguous = 1
!            endif
!         end do
!      end do

      write(*,'(a)') 'INFO: Writing mesh partition to partmesh.txt.'
      OPEN(990,FILE='partmesh.txt')
C...... BEG DW, periodic bcs
      IF ( NPERBC > 0 ) THEN
         PROC( IPERCONN(1:NPERBC,2) ) = PROC(IPERCONN(1:NPERBC,1) ) ;
      END IF
C...... END DW
      DO I=1, MNP
         WRITE(990,*) PROC(I)
      ENDDO
      CLOSE(990)

C...... DW, periodic bcs
      IF ( NPERBC > 0 ) THEN
         NNEG = NNEG_TMP ; ! restore the orignial element table
         DEALLOCATE( NNEG_TMP ) ; ! free memory
         DEALLOCATE( PERBC_IDN_MAP ) ;
      END IF
C...... DW

      deallocate(iduals)
      deallocate(itvect,itvect2)
      deallocate(numduals)
      deallocate(xadj,vwgts,nedges)
      deallocate(nedloc)

 100  FORMAT(4I10)
 200  FORMAT(100I10)
C
      RETURN
c-----------------------------------------------------------------------
      END SUBROUTINE METIS
c-----------------------------------------------------------------------
