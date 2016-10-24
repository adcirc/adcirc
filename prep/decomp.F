C----------------------------------------------------------------------------
C
C                           MODULE DECOMP
C
C----------------------------------------------------------------------------
C
C                  For use with ADCPREP Version 1.6 (  5/21/03 )
C
C                     current for ADCIRC v43.03   5/20/2003
C----------------------------------------------------------------------------
C
      SUBROUTINE DECOMP()
      USE pre_global
C
C---------------------------------------------------------------------------C
C                     (  Serial Version 1.1  5/04/99  )                     C
C                                                                           C
C  Decomposes the ADCIRC grid into NPROC subdomains.                        C
C  The Decomposition Variables are defined in the include file adcprep.inc  C
C  This version is compatible with ADCIRC version 34.03                     C
C                                                                           C
C  12/14/98 vjp  Added interface to METIS 4.0                               C
C   3/10/99 vjp  Rewritten to allow Weir-node pairs to both be ghost nodes  C
C   4/05/99 vjp  Fixed bugs in metis interface routine                      C
C                                                                           C
C---------------------------------------------------------------------------C
C
      INTEGER N1, N2, N3, KMIN, VTMAX
      INTEGER I,J,JD,JG,JP,K,M,ITEMP,ITEMP2,IPR,IPR1,ICNT
      INTEGER ITOT,IEL,IELG,ILNODE,ILNODE2,IPROC,IPROC2
      INTEGER IG1,IG2,IG3,IL1,IL2,IL3,PE1,PE2,PE3
      INTEGER I1,DISC,BBN,IBP,NCOUNT
      INTEGER INDX,INDX2, KK
      INTEGER,ALLOCATABLE :: ITVECT(:)
      CHARACTER PE*6
      CHARACTER(6), PARAMETER :: STRINGGLOBAL = "GLOBAL"
C
      VTMAX = 24*MNP
      ALLOCATE ( ITVECT(VTMAX) )

C   STEP 1:
C-- Use Partition of nodes to compute the number of Resident Nodes
C   to be assigned to each processor.
C-- Then construct Local-to-Global and Global-to-Local Node maps
C   for resident nodes: IMAP_NOD_LG(I,PE),  IMAP_NOD_GL(1:2,I)

      DO I=1, NPROC             ! Use METIS 4.0  Partition
         NOD_RES_TOT(I) = 0
      ENDDO
      DO J=1, NNODG
         NCOUNT = NOD_RES_TOT(PROC(J))+1
         IMAP_NOD_GL(1,J) = PROC(J)
         IMAP_NOD_GL(2,J) = NCOUNT
         IMAP_NOD_LG(NCOUNT,PROC(J)) = J
         NOD_RES_TOT(PROC(J)) = NCOUNT
      ENDDO
C     DO I = 1, NNODG
C        print *, I, IMAP_NOD_GL(1,I)
C     ENDDO

C STEP 2:
C  Construct Local-to-Global Element Map: IMAP_EL_LG(:,PE)
C  Add an element to the map if it has an resident node

      DO I = 1,NPROC
         NELP(I) = 0
         DO K = 1,NELG
            N1 = NNEG(1,K)
            N2 = NNEG(2,K)
            N3 = NNEG(3,K)
            PE1 = IMAP_NOD_GL(1,N1) ! Is any vertex a resident node?
            PE2 = IMAP_NOD_GL(1,N2)
            PE3 = IMAP_NOD_GL(1,N3)
            IF ((PE1.EQ.I).OR.(PE2.EQ.I).OR.(PE3.EQ.I)) THEN
              NELP(I) = NELP(I) + 1
              IMAP_EL_LG(NELP(I),I) = K
            ENDIF
         ENDDO
         IF (NELP(I).GT.MNEP) STOP 'NELP(I) > MNEP'
      ENDDO

C STEP 3:
C--Using Local-to-Global Element map
C  Construct Local-to-Global Node map:  IMAP_NOD_LG(I,PE)
C  and reconstruct Global-to-Local map for resident nodes
C
      DO I = 1,NPROC
         ITOT = 0
         DO J = 1,NELP(I)
            IEL = IMAP_EL_LG(J,I)
            DO M=1, 3
               ITOT = ITOT + 1
               ITVECT(ITOT) = NNEG(M,IEL)
            ENDDO
         ENDDO
         ITEMP = ITOT
         IF (ITOT.GT.VTMAX) stop 'step3 decomp'
         CALL SORT(ITEMP,ITVECT) ! Sort and remove multiple occurrences
         ITOT = 1
         IMAP_NOD_LG(1,I) = ITVECT(1)
         IF (IMAP_NOD_GL(1,ITVECT(1)).EQ.I) THEN
           IMAP_NOD_GL(2,ITVECT(1))=1
         ENDIF
         DO J = 2, ITEMP
           IF (ITVECT(J).NE.IMAP_NOD_LG(ITOT,I)) THEN
             ITOT = ITOT + 1
             IMAP_NOD_LG(ITOT,I) = ITVECT(J)
             IF (IMAP_NOD_GL(1,ITVECT(J)).EQ.I) THEN
               IMAP_NOD_GL(2,ITVECT(J))=ITOT
             ENDIF
           ENDIF
         ENDDO
         NNODP(I) = ITOT
         IF (NNODP(I).GT.MNPP) STOP 'NNODP > MNPP'
      ENDDO
c     print *, "Number of Nodes Assigned to PEs"
c     DO I=1, NPROC
c        print *, I-1, NNODP(I)
c        DO J=1, NNODP(I)
c           print *, J,IMAP_NOD_LG(J,I)
c        ENDDO
c     ENDDO

C STEP 4:
C--If there are any global Weir-node pairs, construct
C  Local-to-Global Weir Node maps: WEIRP_LG(:,PE), WEIRDP_LG(:,PE)
C  Rule: if a global Weir node is assigned ( either as a resident or ghost node )
C        then make it and its dual a local Weir-node pair

      IF (NWEIR.GT.0) THEN
        DO I=1, NPROC
           ITOT = 0
!
! Seizo 2008.07.14
        DO J = 1, NWEIR
           INDX = WEIR(J)
           INDX2= WEIRD(J)
           DO K = 1, NNODP(I)
              N1 = IMAP_NOD_LG(K,I)
              IF( (N1 == INDX) .or. (N1 == INDX2) ) THEN
                  ITOT = ITOT + 1
                  ITVECT(ITOT) = J
                  EXIT
              ENDIF
           ENDDO
        ENDDO
!
!Seizo 2008.07.14          DO J = 1,NNODP(I)
!Seizo 2008.07.14             CALL SEARCH(WEIR,NWEIR,IMAP_NOD_LG(J,I),INDX)
!Seizo 2008.07.14             IF (INDX.NE.0) THEN
!Seizo 2008.07.14               ITOT = ITOT+1
!Seizo 2008.07.14               ITVECT(ITOT) = INDX
!Seizo 2008.07.14             ENDIF
!Seizo 2008.07.14             CALL SEARCH(WEIRD,NWEIR,IMAP_NOD_LG(J,I),INDX2)
!Seizo 2008.07.14             IF (INDX2.NE.0) THEN
!Seizo 2008.07.14               ITOT = ITOT+1
!Seizo 2008.07.14               ITVECT(ITOT) = INDX2
!Seizo 2008.07.14             ENDIF
!Seizo 2008.07.14          ENDDO
           NWEIRP(I) = 0
           ITEMP = ITOT
           IF (ITOT.GT.VTMAX) stop 'step4 decomp'
           IF (ITEMP.GT.1) THEN
             CALL SORT(ITEMP,ITVECT)
             ITOT=1
             INDX = ITVECT(1)
             WEIRP_LG(ITOT,I)  = WEIR(INDX)
             WEIRDP_LG(ITOT,I) = WEIRD(INDX)
             DO J = 2,ITEMP
                IF (ITVECT(J).NE.INDX) THEN
                  INDX = ITVECT(J)
                  ITOT = ITOT+1
                  WEIRP_LG(ITOT,I)  = WEIR(INDX)
                  WEIRDP_LG(ITOT,I) = WEIRD(INDX)
                ENDIF
             ENDDO
             NWEIRP(I) = ITOT
           ENDIF
c       DO J = 1, NWEIRP(I)
c          print *, J, WEIRP_LG(J,I),WEIRDP_LG(J,I)
c       ENDDO
c       print *, "decomp: Number of WEIR node-pairs on PE",I-1,
c    &           " = ",NWEIRP(I)
        ENDDO
      ELSE
        DO I=1, NPROC
           NWEIRP(I) = 0
c          print *, "decomp: Number of WEIR node-pairs on PE",I-1,
c    &              " = ",NWEIRP(I)
        ENDDO
      ENDIF

C STEP 5:
C--If there are any global Weir-node pairs,
C  Re-construct Local-to-Global Element Map: IMAP_EL_LG(:,PE)
C  Rule:  Add an element if it has an resident node or
C         has the dual Weir node of a resident or ghost node

      EL_SHARE(:) = -1
      IF (NWEIR.GT.0) THEN
         DO I = 1,NPROC
            NELP(I) = 0
c           print *, "PE = ",I-1
            DO K = 1,NELG
               N1 = NNEG(1,K)
               N2 = NNEG(2,K)
               N3 = NNEG(3,K)
               PE1 = IMAP_NOD_GL(1,N1)   ! Is any vertex a resident node?
               PE2 = IMAP_NOD_GL(1,N2)
               PE3 = IMAP_NOD_GL(1,N3)   ! belong to a Weir-node pair ?
               CALL SEARCH3(WEIRP_LG(1,I),NWEIRP(I),N1,N2,N3,INDX)
               CALL SEARCH3(WEIRDP_LG(1,I),NWEIRP(I),N1,N2,N3,INDX2)
               ! if any vertex is a resident vertex, or is a weir node
               ! paired to a resident vertex
               IF ((PE1.EQ.I).OR.(PE2.EQ.I).OR.(PE3.EQ.I)
     &             .OR.(INDX.NE.0).OR.(INDX2.NE.0)) THEN
c                 print *, K, PE1,PE2,PE3,INDX,INDX2
                  ! increment the number of elements on this subdomain
                  NELP(I) = NELP(I) + 1
                  ! add the fulldomain element number to the mapping;
                  ! make it negative if some other subdomain has already
                  ! listed it as a resident element (the first subdomain
                  ! that has an element that maps to a fulldomain element
                  ! claims it as a resident; the next subdomain(s) that
                  ! have an element that maps to this fulldomain element
                  ! will have that element number listed as negative, 
                  ! or "ghost" in their mapping tables (first come, first
                  ! served)
                  IF (EL_SHARE(K) > -1) THEN
                     IMAP_EL_LG(NELP(I),I) = -K
                  ELSE
                     IMAP_EL_LG(NELP(I),I) = K
                  ENDIF
                  ! in any case, record the fact that a subdomain has already
                  ! claimed this element as a resident 
                  EL_SHARE(K) = I
               ENDIF
            ENDDO
         IF (NELP(I).GT.MNEP) STOP 'NELP(I) > MNEP'
c           print *, "Number of elements on PE",I-1," = ",NELP(I)
c           DO J = 1, NELP(I)
c              print *, J, IMAP_EL_LG(J,I)
c           ENDDO
         ENDDO
      ELSE
         ! jgf: If there are no weirs, then go back through the element
         ! mapping and leave the first subdomain element that maps to a 
         ! particular fulldomain element as a "resident" with a positive
         ! fulldomain element number in IMAP_EL_LG, while giving subsequent 
         ! subdomain elements that map to that same fulldomain element  
         ! a negative value in IMAP_EL_LG (marking them as "ghosts").
         DO I=1,NPROC
            DO J=1,NELP(I)
               K = IMAP_EL_LG(J,I)
               IF (EL_SHARE(K).gt.-1) THEN
                  IMAP_EL_LG(J,I) = -IMAP_EL_LG(J,I)
               ENDIF
               EL_SHARE(K) = I
            END DO
         END DO
      ENDIF

C STEP 6:
C--Using Local-to-Global Element map
C  Construct Local-to-Global Node map:  IMAP_NOD_LG(I,PE)
C  and reconstruct Global-to-Local map for resident nodes
C
      DO I = 1,NPROC
         ITOT = 0
         DO J = 1,NELP(I)
            IEL = abs(IMAP_EL_LG(J,I))
            DO M=1, 3
               ITOT = ITOT + 1
               ITVECT(ITOT) = NNEG(M,IEL)
            ENDDO
         ENDDO
         ITEMP = ITOT
         IF (ITOT.GT.VTMAX) stop 'step6 decomp'
         CALL SORT(ITEMP,ITVECT) ! Sort and remove multiple occurrences
         ITOT = 1
         IMAP_NOD_LG(1,I) = ITVECT(1)
         IF (IMAP_NOD_GL(1,ITVECT(1)).EQ.I) THEN
           IMAP_NOD_GL(2,ITVECT(1))=1
         ENDIF
         DO J = 2, ITEMP
           IF (ITVECT(J).NE.IMAP_NOD_LG(ITOT,I)) THEN
             ITOT = ITOT + 1
             IMAP_NOD_LG(ITOT,I) = ITVECT(J)
             IF (IMAP_NOD_GL(1,ITVECT(J)).EQ.I) THEN
               IMAP_NOD_GL(2,ITVECT(J))=ITOT
             ENDIF
           ENDIF
         ENDDO
         NNODP(I) = ITOT
         IF (NNODP(I).GT.MNPP) STOP 'NNODP > MNPP'
      ENDDO
c     print *, "Number of Nodes Assigned to PEs"
c     DO I=1, NPROC
c        print *, I-1, NNODP(I)
c        DO J=1, NNODP(I)
c           print *, J,IMAP_NOD_LG(J,I)
c        ENDDO
c     ENDDO

C STEP 7:
C--Construct Local Element Connectivity Table for each PE: NNEP(3,I,PE)
C
      ! loop over subdomains
      DO I = 1,NPROC
         ! ITEMP = number of nodes on this subdomain
         ITEMP = NNODP(I)
         ! loop over local nodes on this subdomain
         DO J = 1,NNODP(I)
            ! make array of corresponding fulldomain node numbers
            ITVECT(J) = IMAP_NOD_LG(J,I)
         ENDDO
         ! loop over elements on this subdomain
         DO J = 1,NELP(I)
            ! grab fulldomain element number, whether ghost or resident element
            IELG = abs(IMAP_EL_LG(J,I))
            ! loop over nodes on this element
            DO M = 1,3
               ! get fulldomain node number on this element
               IG1 = NNEG(M,IELG)
               ! look for the fulldomain node number in the list of 
               ! fulldomain nodes in this domain
               CALL LOCATE(ITVECT,ITEMP,IG1,K)
               ! if the fulldomain node number on this element below
               ! the range of fulldomain node numbers on this subdomain...
               IF (K.LE.0) THEN
                  ! if the next index value matches
                  IF (IMAP_NOD_LG(K+1,I).EQ.IG1) THEN
                     ! add the next index value to the local element
                     ! connectivity table
                     NNEP(M,J,I) = K+1
                  ELSE
                     STOP 'ERROR IN IMAP_NOD_LG'
                  ENDIF
               ELSEIF (K.GE.NNODP(I))THEN
                  IF (IMAP_NOD_LG(K,I).EQ.IG1) THEN
                     NNEP(M,J,I) = K
                  ELSE
                     STOP 'ERROR IN IMAP_NOD_LG'
                  ENDIF
               ELSE
                  IF (IMAP_NOD_LG(K,I).EQ.IG1) THEN
                     NNEP(M,J,I) = K
                  ELSEIF (IMAP_NOD_LG(K+1,I).EQ.IG1) THEN
                     NNEP(M,J,I) = K+1
                  ELSE
                     STOP 'ERROR IN IMAP_NOD_LG'
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO

C STEP 8:
C--Compute the number of communicating PEs and their
C  list for each PE:  NUMM_COMM_PE(PE) and COMM_PE_NUM(IPE,PE)
C
      !
      ! ITVECT is a list of the neighboring subdomain numbers  
      ! corresponding to this subdomain's ghost nodes
      ! 
      ! loop over subdomains
      DO I = 1,NPROC
         NUM_COMM_PE(I) = 0  ! zero the num of neighboring subdomains
         ITEMP = 0           ! zero the number of ghost nodes on this subdomain
         ! loop over subdomain nodes
         DO J = 1,NNODP(I)
            ! get corresponding fulldomain node number
            INDX = IMAP_NOD_LG(J,I)
            ! use fulldomain node number to get number of subdomain
            ! where this node is a resident
            IPR = IMAP_NOD_GL(1,INDX)
            ! if the subdomain of residence is not this one
            IF (IPR.NE.I)THEN
               ! increment number of ghost nodes on this subdomain
               ITEMP = ITEMP + 1
               ! add the subdomain number to which the ghost node belongs
               ! to the list of neighboring subdomains
               ITVECT(ITEMP) = IPR
            ENDIF
         ENDDO
         IF (ITEMP.EQ.0) THEN
            NUM_COMM_PE(I) = 0
         ELSE
            IF (ITEMP.GT.VTMAX) stop 'step8 decomp'
            ! sort the list of neighboring subdomains numbers in ascending order
            CALL SORT(ITEMP,ITVECT)
            ! set the first neighboring subdomain
            COMM_PE_NUM(1,I) = ITVECT(1)
            ! reset the subdomain neighbor counter
            ITOT = 1
            ! loop over the subdomain numbers of the ghost nodes 
            ! in this subdomain        
            DO J = 1,ITEMP
               ! if the subdomain number is not a duplicate
               IF (ITVECT(J).NE.COMM_PE_NUM(ITOT,I)) THEN
                  ! increment the number of neighboring subdomains
                  ITOT = ITOT + 1
                  ! record the number of the neighboring subdomain
                  COMM_PE_NUM(ITOT,I) = ITVECT(J)
               ENDIF
            ENDDO
            ! set the total number of neighbor subdomains for this
            ! subdomain
            NUM_COMM_PE(I) = ITOT
            IF (NUM_COMM_PE(I).GT.MNPROC) STOP 'NUM_COMM_PE>MNPROC'
         ENDIF
      ENDDO

Casey 110518: Added this change from Seizo.
!     Send <=> Receive Proc CHECK & CORRECTION !st3 06.19.2009
      ! loop over subdomains
      DO I = 1, NPROC
         ! loop over subdomains neighboring this one
         DO J = 1, NUM_COMM_PE(I)
            ! get neighboring subdomain number
            INDX = COMM_PE_NUM(J,I)
            ! zero the flag indicating that this subdomain is listed
            ! as the neighboring subdomain's neighbor 
            ITOT = 0   
            ! loop over the subdomains neighboring this neighbor 
            DO K = 1, NUM_COMM_PE(INDX)
               ! get the neighbor's neighboring subdomain number
               INDX2 = COMM_PE_NUM(K,INDX)
               ! set a flag if this subdomain is listed as the 
               ! neighbor's neighbor
               IF( I == INDX2 ) ITOT = 1
            ENDDO
            ! if the flag was not set
            IF( ITOT == 0 ) THEN
               ! increment neighbor's number of subdomain neighbors 
               NUM_COMM_PE(INDX) = NUM_COMM_PE(INDX) + 1
               ! add this subdomain to this list of neighbors for
               ! the neighboring subdomain
               COMM_PE_NUM( NUM_COMM_PE(INDX), INDX )  = I
               ! I don't see how this can happen, so I am adding a 
               ! warning message to provide an alert when it does
               ! so that I can examine it more closely. 
               write(6,'(a,i0,a,i0,a)') 'WARNING: Subdomain ',I, 
     &       ' added to the list of subdomain neighbors of subdomain ',
     &         INDX,'.' 
            ENDIF
         ENDDO
      ENDDO

C STEP 9:
C--Construct a Global-to-Local node mapping: IMAP_NOD_GL2(*,J)
C  This is not a function, but is rather a relation
C  It works for both resident and ghost nodes

      DO I = 1,NNODG
         ITOTPROC(I) = 0
      ENDDO
      DO I = 1,NPROC
         DO J = 1,NNODP(I)
            INDX = IMAP_NOD_LG(J,I)
            ITOTPROC(INDX) = ITOTPROC(INDX) + 1
            IF (ITOTPROC(INDX).GT.MNPROC)THEN
              WRITE(6,*)'Some nodes belong to more processors',
     $                  ' than MNPROC'
              STOP
            ENDIF
            ITEMP = (ITOTPROC(INDX)-1)*2 + 1
            IMAP_NOD_GL2(ITEMP,INDX) = I
            IMAP_NOD_GL2(ITEMP+1,INDX) = J
         ENDDO
      ENDDO
c     print *, "Global Nodes assigned to more than one PE"
c     do J=1, NNODG
c        if (ITOTPROC(J).GT.1) print *, J, ITOTPROC(J)
c     enddo

C STEP 10:
C--Print Summary of Decomposition
C
      print *, "Decomposition Data"
      print *, "DOMAIN  RES_NODES  GHOST_NODES  TOT_NODES  ELEMENTS"
      print *, "------  ---------  -----------  ---------  --------"
      WRITE(*,90) STRINGGLOBAL,NNODG,NELG
      DO I=1, NPROC
         PE(1:6) = 'PE0000'
         CALL IWRITE(PE,3,6,I-1)
         WRITE(6,92) PE, NOD_RES_TOT(I),NNODP(I)-NOD_RES_TOT(I),
     &                   NNODP(I),NELP(I)
      ENDDO
 90   FORMAT(1X,A6,25X,I9,2X,I9)
 92   FORMAT(1X,A6,1X,I9,2X,I9,4X,I9,2X,I9)
C
      RETURN
      END


      SUBROUTINE DOMSIZE()
      USE pre_global
C
C---------------------------------------------------------------------------C
C                   (  Serial Version 1.0  12/20/99 vjp )                   C
C                                                                           C
C  Takes dry run through the domain decomp to determine the max number of   C
C  nodes and elements assigned to any subdomain to determine MNPP and MNEP. C
C                                                                           C
C---------------------------------------------------------------------------C
C
      INTEGER N1,N2,N3,VTMAX
      INTEGER I,J,K,M,ITEMP
      INTEGER ITOT,IEL,IELG,ILNODE,ILNODE2,IPROC
      INTEGER IG1,IG2,IG3,IL1,IL2,IL3,PE1,PE2,PE3
      INTEGER INDX,INDX2
      INTEGER RESNODE,NODES,NELEM,ONELEM,NLWEIR
C
      INTEGER,ALLOCATABLE :: ITVECT(:)
      INTEGER,ALLOCATABLE :: NODE_LG(:)
      INTEGER,ALLOCATABLE :: NODE_GL1(:)
      INTEGER,ALLOCATABLE :: NODE_GL2(:)
      INTEGER,ALLOCATABLE :: ELEM_LG(:)
      INTEGER,ALLOCATABLE :: LWEIR_LG(:),LWEIRD_LG(:)
C
      VTMAX = 24*MNP
      ALLOCATE ( ITVECT(VTMAX) )
      ALLOCATE ( NODE_LG(MNP) )
      ALLOCATE ( NODE_GL1(MNP) )
      ALLOCATE ( NODE_GL2(MNP) )
      ALLOCATE ( LWEIR_LG(MNP),LWEIRD_LG(MNP) )
      ALLOCATE ( ELEM_LG(MNE) )
C
      MNPP = 0
      MNEP = 0
C
      DO 1000 IPROC=1, NPROC
C
C   STEP 1:
C-- Use Partition of nodes to compute the number of Resident Nodes
C   to be assigned to each processor.
C-- Then construct Local-to-Global and Global-to-Local Node maps
C   for resident nodes
C
C
      DO J=1, MNP
         NODE_GL1(J) = 0
         NODE_GL2(J) = 0
         LWEIR_LG(J) = 0
         LWEIRD_LG(J) = 0
      ENDDO
C
      DO J=1, MNE
         ELEM_LG(J) = 0
      ENDDO
C
      RESNODE = 0
      NODES   = 0
      DO J=1, NNODG
         IF (IPROC.EQ.PROC(J)) THEN
           RESNODE = RESNODE+1
           NODE_GL1(J) = PROC(J)
           NODE_GL2(J) = RESNODE
           NODE_LG(RESNODE) = J
         ENDIF
      ENDDO
C     DO I = 1, NNODG
C        print *, I, NODE_GL1(I)
C     ENDDO

C STEP 2:
C  Construct Local-to-Global Element Map
C  Add an element to the map if it has an resident node

      NELEM   = 0
      DO K = 1,NELG
         N1 = NNEG(1,K)
         N2 = NNEG(2,K)
         N3 = NNEG(3,K)
         PE1 = NODE_GL1(N1) ! Is any vertex a resident node?
         PE2 = NODE_GL1(N2)
         PE3 = NODE_GL1(N3)
         IF ((PE1.EQ.IPROC).OR.(PE2.EQ.IPROC).OR.(PE3.EQ.IPROC)) THEN
           NELEM = NELEM+1
           ELEM_LG(NELEM) = K
         ENDIF
      ENDDO

C STEP 3:
C--Using Local-to-Global Element map
C  reconstruct Local-to-Global Node map
C  and Global-to-Local map for resident nodes
C
      ITOT = 0
      DO J = 1,NELEM
         IEL = ELEM_LG(J)
         DO M=1, 3
            ITOT = ITOT + 1
            ITVECT(ITOT) = NNEG(M,IEL)
         ENDDO
      ENDDO
      ITEMP = ITOT
      IF (ITOT.GT.VTMAX) stop 'step3 decomp'
      CALL SORT(ITEMP,ITVECT) ! Sort and remove multiple occurrences
      ITOT = 1
      NODE_LG(1) = ITVECT(1)
      IF (NODE_GL1(ITVECT(1)).EQ.IPROC) THEN
        NODE_GL2(ITVECT(1))=1
      ENDIF
      DO J = 2, ITEMP
         IF (ITVECT(J).NE.NODE_LG(ITOT)) THEN
           ITOT = ITOT + 1
           NODE_LG(ITOT) = ITVECT(J)
           IF (NODE_GL1(ITVECT(J)).EQ.IPROC) THEN
             NODE_GL2(ITVECT(J))=ITOT
           ENDIF
         ENDIF
      ENDDO
      NODES = ITOT


C STEP 4:
C--If there are any global Weir-node pairs, construct
C  Local-to-Global Weir Node maps
C  Rule: if a global Weir node is assigned ( as resident or ghost node )
C        then make it and its dual a local Weir-node pair
      IF (NWEIR.GT.0) THEN
        ITOT = 0
! Seizo 2008.07.14
        DO J = 1, NWEIR
           INDX = WEIR(J)
           INDX2= WEIRD(J)
           DO K = 1, NODES
              N1 = NODE_LG(K)
              IF( (N1 == INDX) .or. (N1 == INDX2) ) THEN
                  ITOT = ITOT + 1
                  ITVECT(ITOT) = J
                  EXIT
              ENDIF
           ENDDO
        ENDDO

!Seizo 2008.07.14       DO J = 1,NODES
!Seizo 2008.07.14          CALL SEARCH(WEIR,NWEIR,NODE_LG(J),INDX)
!Seizo 2008.07.14          IF (INDX.NE.0) THEN
!Seizo 2008.07.14            ITOT = ITOT+1
!Seizo 2008.07.14            ITVECT(ITOT) = INDX
!Seizo 2008.07.14          ENDIF
!Seizo 2008.07.14          CALL SEARCH(WEIRD,NWEIR,NODE_LG(J),INDX2)
!Seizo 2008.07.14          IF (INDX2.NE.0) THEN
!Seizo 2008.07.14            ITOT = ITOT+1
!Seizo 2008.07.14            ITVECT(ITOT) = INDX2
!Seizo 2008.07.14          ENDIF
!Seizo 2008.07.14       ENDDO
        NLWEIR  = 0
        ITEMP = ITOT
        IF (ITOT.GT.VTMAX) stop 'step4 decomp'
        IF (ITEMP.GT.1) THEN
          CALL SORT(ITEMP,ITVECT)
          ITOT=1
          INDX = ITVECT(1)
          LWEIR_LG(ITOT)  = WEIR(INDX)
          LWEIRD_LG(ITOT) = WEIRD(INDX)
          DO J = 2,ITEMP
             IF (ITVECT(J).NE.INDX) THEN
               INDX = ITVECT(J)
               ITOT = ITOT+1
               LWEIR_LG(ITOT)  = WEIR(INDX)
               LWEIRD_LG(ITOT) = WEIRD(INDX)
             ENDIF
          ENDDO
          NLWEIR = ITOT
        ENDIF
      ELSE
          NLWEIR = 0
      ENDIF
c     print *, "domsize: Number of WEIR node-pairs on PE",IPROC-1,
c    &         " = ",NLWEIR
      IF (NLWEIR.GT.NWEIR) THEN
        print *, "error in domsize: "
        print *, "local number of weir-pairs exceeds total"
        stop
      ENDIF

C STEP 5:
C--If there are any global Weir-node pairs,
C  Re-construct Local-to-Global Element Map: IMAP_EL_LG(:,PE)
C  Rule:  Add an element if it has an resident node or
C         has the dual Weir node of a resident or ghost node

      ONELEM = NELEM     ! Save NELEM for PEs with no WEIR-pairs
      NELEM  = 0
      IF (NLWEIR.GT.0) THEN
        DO K = 1,NELG
           N1 = NNEG(1,K)
           N2 = NNEG(2,K)
           N3 = NNEG(3,K)
           PE1 = NODE_GL1(N1)   ! Is any vertex a resident node?
           PE2 = NODE_GL1(N2)
           PE3 = NODE_GL1(N3)   ! belong to a Weir-node pair ?
           CALL SEARCH3(LWEIR_LG(1),NLWEIR,N1,N2,N3,INDX)
           CALL SEARCH3(LWEIRD_LG(1),NLWEIR,N1,N2,N3,INDX2)
           IF ((PE1.EQ.IPROC).OR.(PE2.EQ.IPROC).OR.(PE3.EQ.IPROC)
     &          .OR.(INDX.NE.0).OR.(INDX2.NE.0)) THEN
             NELEM = NELEM + 1
             ELEM_LG(NELEM) = K
           ENDIF
        ENDDO
C
      ENDIF

C STEP 6:
C--Using Local-to-Global Element map, reconstruct Local-to-Global Node map
C
      IF (NELEM.EQ.0) NELEM = ONELEM   ! if necessary restore old nelem
      ITOT = 0
      DO J = 1,NELEM
         IEL = ELEM_LG(J)
         DO M=1, 3
            ITOT = ITOT + 1
            ITVECT(ITOT) = NNEG(M,IEL)
         ENDDO
      ENDDO
      ITEMP = ITOT
      IF (ITOT.GT.VTMAX) stop 'step6 decomp'
      CALL SORT(ITEMP,ITVECT) ! Sort and remove multiple occurrences
      ITOT = 1
      NODE_LG(1) = ITVECT(1)
      DO J = 2, ITEMP
         IF (ITVECT(J).NE.NODE_LG(ITOT)) THEN
           ITOT = ITOT + 1
           NODE_LG(ITOT) = ITVECT(J)
         ENDIF
      ENDDO
C
      NODES = ITOT
      IF (NODES.GT.MNPP) MNPP = NODES
c     print *, "Number of nodes on PE",IPROC-1," = ",NODES
      IF (NELEM.GT.MNEP) MNEP = NELEM
c     print *, "Number of elements on PE",IPROC-1," = ",NELEM

1000  CONTINUE
C
      DEALLOCATE ( ITVECT,NODE_LG,NODE_GL1,NODE_GL2,ELEM_LG,
     &             LWEIR_LG,LWEIRD_LG )
C
      print *, " Setting MNPP = ",MNPP
      print *, " Setting MNEP = ",MNEP
C
      RETURN
      END






      SUBROUTINE SORT(N,RA)
      IMPLICIT NONE
      INTEGER N, L, IR, RRA, I, J
      INTEGER RA(N)
C
C---------------------------------------------------------------------------
C  Sorts array RA of length N into ascending order using Heapsort algorithm.
C  N is input; RA is replaced on its output by its sorted rearrangement.
C  Ref: Numerical Recipes
C---------------------------------------------------------------------------
C
      L = N/2 + 1
      IR = N
10    CONTINUE
      IF (L.GT.1)THEN
        L=L-1
        RRA = RA(L)
      ELSE
        RRA=RA(IR)
        RA(IR)=RA(1)
        IR=IR-1
        IF (IR.EQ.1) THEN
          RA(1)=RRA
          RETURN
        ENDIF
      ENDIF
      I=L
      J=L+L
20    IF (J.LE.IR) THEN
        IF (J.LT.IR) THEN
          IF(RA(J).LT.RA(J+1)) J=J+1
        ENDIF
        IF (RRA.LT.RA(J)) THEN
          RA(I)=RA(J)
          I=J
          J=J+J
        ELSE
          J=IR+1
        ENDIF
        GO TO 20
      ENDIF
      RA(I)=RRA
      GO TO 10
      END

      SUBROUTINE LOCATE(XX,N,X,J)
      IMPLICIT NONE
      INTEGER JM,JL,JU,J,N,X,XX(N)
C
C--Given an array XX of length N, and given a value X, returns a value J
C--such that X is between XX(J) and XX(J+1). XX must be monotonic, either
C--increasing or decreasing. J=0 or J=N is returned to indicate that X is out
C--of range.
C--
C--NUMERICAL RECIPES - The Art of Scientific Computing [FORTRAN Version]

C--Initialize lower and upper limits
         JL = 0
         JU = N+1
C
C--If we are not done yet, compute a mid-point, and replace either the lower
C--limit or the upper limit, as appropriate.
C
10       IF(JU-JL.GT.1)THEN
            JM = (JU+JL)/2
            IF((XX(N).GT.XX(1)).EQV.(X.GT.XX(JM)))THEN
              JL = JM
            ELSE
              JU = JM
            ENDIF
C--Repeat until the test condition 10 is satisfied.
         GO TO 10
         ENDIF
C--Then set the output and return.
         J = JL
         RETURN
         END

      SUBROUTINE SEARCH3(MAP,LENG,N1,N2,N3,INDX)
      INTEGER MAP(*),LENG,N1,N2,N3,INDX,IP
cvjp  rewritten 5/3/99
      INDX = 0
      DO I=1,LENG
         IP = MAP(I)
         IF (IP.EQ.N1.OR.IP.EQ.N2.OR.IP.EQ.N3) THEN
           INDX = I
           GOTO 99
         ENDIF
      ENDDO
 99   RETURN
      END

      SUBROUTINE SEARCH(MAP,N,TARGET,INDX)
      INTEGER MAP(*),N,TARGET,INDX,I
C
      INDX = 0
      IF (N.EQ.0) GOTO 99
      DO I=1,N
         IF (MAP(I).EQ.TARGET) THEN
            INDX = I
            GOTO 99
         ENDIF
      ENDDO
99    RETURN
      END
