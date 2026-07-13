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
      SUBROUTINE DECOMP()
      USE pre_global

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
         IF (NELP(I).GT.MNEP) THEN
             WRITE(*,'(A)') 'ERROR: NELP(I) > MNEP'
             CALL EXIT(1)
         ENDIF
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
         IF (ITOT.GT.VTMAX) THEN
             WRITE(*,'(A)') "Error during step3 decomp"
             CALL EXIT(1)
         ENDIF
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
         IF (NNODP(I).GT.MNPP) THEN
             WRITE(*,'(A)') "ERROR:  NNODP > MNPP"
             CALL EXIT(1)
         ENDIF
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
           IF (ITOT.GT.VTMAX) THEN
               WRITE(*,'(A)') 'Error during step4 decomp'
               CALL EXIT(1)
           ENDIF
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
!C                   EL_SHARE(K) = I

!c                   
!c DW: 2023, for NOAA-CMEP
!c           If all three nodes are a ghost node, the element is
!c           taken to be a ghost element.                     
                  IF ( (PE1 .NE. I) .AND. (PE2 .NE. I) 
     &                                .AND. (PE3 .NE. I)  ) THEN
                      IF ( EL_SHARE(K) < 0 ) THEN
                         ! set the element as a ghost element  
                         IMAP_EL_LG(NELP(I),I) = -K ;

                         ! set the element as an unclaimed element 
                         EL_SHARE(K) = -I ;
                      END IF      
                  ELSE
                      EL_SHARE(K) =  I ; 
                  ENDIF        
!c END DW
!c 
               ENDIF
            ENDDO
         IF (NELP(I).GT.MNEP) THEN
             WRITE(*,'(A)') "ERROR: NELP(I) > MNEP"
             CALL EXIT(1)
         ENDIF
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
         IF (ITOT.GT.VTMAX) THEN
             WRITE(*,'(A)') "Error during step6 decomp"
             CALL EXIT(1)
         ENDIF
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
         IF (NNODP(I).GT.MNPP) THEN
             WRITE(*,'(A)') "ERROR: NNODP > MNPP"
             CALL EXIT(1)
         ENDIF
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
                     WRITE(*,'(A)') 'ERROR IN IMAP_NOD_LG'
                     CALL EXIT(1) 
                  ENDIF
               ELSEIF (K.GE.NNODP(I))THEN
                  IF (IMAP_NOD_LG(K,I).EQ.IG1) THEN
                     NNEP(M,J,I) = K
                  ELSE
                     WRITE(*,'(A)') 'ERROR IN IMAP_NOD_LG'
                     CALL EXIT(1)
                  ENDIF
               ELSE
                  IF (IMAP_NOD_LG(K,I).EQ.IG1) THEN
                     NNEP(M,J,I) = K
                  ELSEIF (IMAP_NOD_LG(K+1,I).EQ.IG1) THEN
                     NNEP(M,J,I) = K+1
                  ELSE
                     WRITE(*,'(A)') 'ERROR IN IMAP_NOD_LG'
                     CALL EXIT(1) 
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
            IF (ITEMP.GT.VTMAX) THEN
                WRITE(*,'(A)') "Error during step 8, decomp"
                CALL EXIT(1) 
            ENDIF
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
            IF (NUM_COMM_PE(I).GT.MNPROC) THEN
                WRITE(*,'(A)') "ERROR: NUM_COMM_PE>MNPROC"
                CALL EXIT(1)
            ENDIF
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
              CALL EXIT(1)
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
         WRITE(PE(3:6),'(I4.4)') I-1
         WRITE(6,92) PE, NOD_RES_TOT(I),NNODP(I)-NOD_RES_TOT(I),
     &                   NNODP(I),NELP(I)
      ENDDO
 90   FORMAT(1X,A6,25X,I9,2X,I9)
 92   FORMAT(1X,A6,1X,I9,2X,I9,4X,I9,2X,I9)
C
      RETURN
C----------------------------------------------------------------------
      END SUBROUTINE DECOMP
C----------------------------------------------------------------------

      SUBROUTINE DECOMP_DG()
      USE pre_global

      IMPLICIT NONE
      INTEGER N1, N2, N3, KMIN, VTMAX
      INTEGER I,J,JD,JG,JP,K,M,ITEMP,ITEMP2,IPR,IPR1,ICNT

C.....namo 5/22 add DG parts
      integer LL ! dummy for L which is used in pre_global
      LOGICAL DG_6
Csb-DG1
      INTEGER IPR2,IPR3,M1,M2,M3,A,B,IEDG,INEIGHELG,INEIGHEDG
      INTEGER INEIELG,NELPI
      INTEGER,ALLOCATABLE :: COMM_PE_FLAG(:)
C--
Csb--02/26/06
      INTEGER NM1,NM2,NM3,IK,JJ


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
      DO I = 1, NNODG
        !print *, I, IMAP_NOD_GL(2,I)
      ENDDO

C STEP 2:
C  Construct Local-to-Global Element Map: IMAP_EL_LG(:,PE)
C  Add an element to the map if it has an resident node

      DO I = 1,NPROC
         NELP(I) = 0
         ! Loop through global element list, each element K
         DO K = 1,NELG
            N1 = NNEG(1,K)
            N2 = NNEG(2,K)
            N3 = NNEG(3,K)
            PE1 = IMAP_NOD_GL(1,N1) ! Is any vertex a resident node?
            PE2 = IMAP_NOD_GL(1,N2)
            PE3 = IMAP_NOD_GL(1,N3)
            ! if any node of this element belongs to me
            IF ((PE1.EQ.I).OR.(PE2.EQ.I).OR.(PE3.EQ.I)) THEN
              NELP(I) = NELP(I) + 1
              IMAP_EL_LG(NELP(I),I) = K
            ENDIF
         ENDDO
         !print *, 'NELP(I)', NELP(I)
         IF (NELP(I).GT.MNEP) THEN
             WRITE(*,'(A)') 'ERROR: NELP(I) > MNEP'
             CALL EXIT(1)
         ENDIF
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
         IF (ITOT.GT.VTMAX) THEN
             WRITE(*,'(A)') "Error during step3 decomp"
             CALL EXIT(1)
         ENDIF
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
         IF (NNODP(I).GT.MNPP) THEN
             WRITE(*,'(A)') "ERROR:  NNODP > MNPP"
             CALL EXIT(1)
         ENDIF
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
           IF (ITOT.GT.VTMAX) THEN
               WRITE(*,'(A)') 'Error during step4 decomp'
               CALL EXIT(1)
           ENDIF
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
         IF (NELP(I).GT.MNEP) THEN
             WRITE(*,'(A)') "ERROR: NELP(I) > MNEP"
             CALL EXIT(1)
         ENDIF
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
               ELSE ! namo - mark the first element to claim element K
                  EL_SHARE(K) = I
               ENDIF
            END DO
         END DO
      ENDIF

      !print *, 'NELG', NELG

C     ....namo 5/22 - add DG part

C STEP 5-DG-1:
C--Alloc memory 1
C
      CALL ALLOC_DG1()

      print *, 'Finished alloc_dg1()'
C STEP 5-DG-2
C--Count maximum of the number of the elements associated with a node
      NNDEL(:) = 0
      DO IK=1,MNE
        NNDEL(NNEG(1,IK)) = NNDEL(NNEG(1,IK)) + 1
        NNDEL(NNEG(2,IK)) = NNDEL(NNEG(2,IK)) + 1
        NNDEL(NNEG(3,IK)) = NNDEL(NNEG(3,IK)) + 1
      ENDDO
      MNNDEL = 0
      DO IK=1,MNP
        IF(NNDEL(IK).GT.MNNDEL) MNNDEL = NNDEL(IK)
      ENDDO

C STEP 5-DG-3:
C--Alloc memory 1B
C
      CALL ALLOC_DG1B()

      print *, 'Finished alloc_dg1b'

C STEP 5-DG-4:
C--Make node-to-elements table
      NNDEL(:) = 0

      DO IK=1,MNE
        NM1 = NNEG(1,IK)
        NM2 = NNEG(2,IK)
        NM3 = NNEG(3,IK)

        NNDEL(NM1) = NNDEL(NM1) + 1
        NDEL(NM1,NNDEL(NM1)) = IK

        NNDEL(NM2) = NNDEL(NM2) + 1
        NDEL(NM2,NNDEL(NM2)) = IK

        NNDEL(NM3) = NNDEL(NM3) + 1
        NDEL(NM3,NNDEL(NM3)) = IK
      ENDDO

C STEP 5-DG-5:
C--Make a index of the neighboring element of each edge
C
      IMAP_NEIGHEDG(:,:,:) = 0

      DO I=1,NELG
        DO 2001 A=1,3
          IF(IMAP_NEIGHEDG(1,A,I).NE.0) GOTO 2001

          N1 = NNEG(MOD(A+0,3)+1,I) ! Endnode 1 of edge A
          N2 = NNEG(MOD(A+1,3)+1,I) ! Endnode 2 of edge A

          DO JJ=1,NNDEL(N1)
            J = NDEL(N1,JJ)
            IF(I.NE.J) THEN
              DO B=1,3
                M1 = NNEG(MOD(B+0,3)+1,J) ! Endnode 1 of edge B
                M2 = NNEG(MOD(B+1,3)+1,J) ! Endnode 2 of edge B

                IF(N1.EQ.M2.AND.N2.EQ.M1) THEN
                  IMAP_NEIGHEDG(1,A,I) = J
                  IMAP_NEIGHEDG(2,A,I) = B

                  IMAP_NEIGHEDG(1,B,J) = I
                  IMAP_NEIGHEDG(2,B,J) = A

                  GOTO 2001
                ENDIF
              ENDDO
            ENDIF
          ENDDO
 2001   CONTINUE
      ENDDO

C STEP 5-DG-6:
C--If there are any elements shared by three different sub-domains,
C  Re-construct Local-to-Global Element Map: IMAP_EL_LG(:,PE)
C  Rule: Add non-local elements adjacent to the elements
C        shared by three different sub-domains
C        to the local element list
C
C--namo 5/22 - does reconstruction affect other parts from ADCIRC?
C  leave it as is for now...
C  Should we use negative convention for ghost elements?
      IF (.true.) THEN
      DO I=1,NPROC
        NELPI = NELP(I)
        !print *, NELP(I)
        DO J=1,NELPI
           ! namo - add abs
          IELG = abs(IMAP_EL_LG(J,I))
c$$$          IF (IELG .lt. 0) THEN
c$$$             GOTO 520
c$$$          END IF

          N1 = NNEG(1,IELG)
          N2 = NNEG(2,IELG)
          N3 = NNEG(3,IELG)

          PE1 = IMAP_NOD_GL(1,N1)
          PE2 = IMAP_NOD_GL(1,N2)
          PE3 = IMAP_NOD_GL(1,N3)

C          IF(((PE1.NE.PE2).AND.(PE2.NE.PE3).AND.(PE3.NE.PE1)).AND.
C     &     (PE1.EQ.I.OR.PE2.EQ.I.OR.PE3.EQ.I)) THEN
          IF((PE1.NE.PE2).AND.(PE2.NE.PE3).AND.(PE3.NE.PE1)) THEN
            IF(PE1.EQ.I) THEN
              K = 1
            ELSE IF(PE2.EQ.I) THEN
              K = 2
            ELSE IF(PE3.EQ.I) THEN
              K = 3
            ENDIF

            !IELG = IMAP_EL_LG(J,I)
            INEIELG = IMAP_NEIGHEDG(1,K,IELG)

!Cmm   Do NOT add the neighboring element if the current
! element is on the global boundary.
            if (ineielg /= 0) then

!Cmm   Also do NOT add the neighboring element if it
! is already a member of the current subdomain.
!c                  if (I == 4 .or. i == 8) then
!c                     write(*,*) "NELP(I) = ", NELP(I)
!c                     write(*,*) "NELPI = ", NELPI
!c                  endif
               do LL = 1,NELP(I)
                  !print *, 'LL', IMAP_EL_LG(LL,I), INEIELG
                  if (abs(IMAP_EL_LG(LL,I)) == INEIELG) then
!Cmm   As much as I cringe at using GO TO's,
! I feel one is necessary here
!     since we need to break out of a nested loop
! and give up on element J.
!if (i == 4 .or. i == 8) then
!   write(*,*) "Found a match: L = ", L
!   write(*,*) "INEIELG = ", INEIELG
!endif
                     go to 520
                  end if
               end do

               !print *, NELP(I), MNEP
               IF(NELP(I).GE.MNEP) THEN
                  STOP 'FATAL ERROR (STEP 5-DG-6 (2))'
               ELSE
                  NELP(I) = NELP(I) + 1
                  IMAP_EL_LG(NELP(I),I) = INEIELG
               ENDIF
!c               endif
            ENDIF
         ENDIF
 520     continue
      ENDDO
      ENDDO

      print *, 'Finished DG-6'
      ENDIF
C--
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
         IF (ITOT.GT.VTMAX) THEN
             WRITE(*,'(A)') "Error during step6 decomp"
             CALL EXIT(1)
         ENDIF
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
         IF (NNODP(I).GT.MNPP) THEN
             WRITE(*,'(A)') "ERROR: NNODP > MNPP"
             CALL EXIT(1)
         ENDIF
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
                     WRITE(*,'(A)') 'ERROR IN IMAP_NOD_LG'
                     CALL EXIT(1)
                  ENDIF
               ELSEIF (K.GE.NNODP(I))THEN
                  IF (IMAP_NOD_LG(K,I).EQ.IG1) THEN
                     NNEP(M,J,I) = K
                  ELSE
                     WRITE(*,'(A)') 'ERROR IN IMAP_NOD_LG'
                     CALL EXIT(1)
                  ENDIF
               ELSE
                  IF (IMAP_NOD_LG(K,I).EQ.IG1) THEN
                     NNEP(M,J,I) = K
                  ELSEIF (IMAP_NOD_LG(K+1,I).EQ.IG1) THEN
                     NNEP(M,J,I) = K+1
                  ELSE
                     WRITE(*,'(A)') 'ERROR IN IMAP_NOD_LG'
                     CALL EXIT(1)
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
            IF (ITEMP.GT.VTMAX) THEN
                WRITE(*,'(A)') "Error during step 8, decomp"
                CALL EXIT(1)
            ENDIF
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
            IF (NUM_COMM_PE(I).GT.MNPROC) THEN
                WRITE(*,'(A)') "ERROR: NUM_COMM_PE>MNPROC"
                CALL EXIT(1)
            ENDIF
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
              CALL EXIT(1)
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

C.....namo 5/22 - add DG part
      print *, 'Finished step 9'

!Csb-DG2
!C STEP 9-DG-1:
!C--Compute NIEL_SEND(:) and NIEL_RECV(:) to prepare for ALLOC_DG2
!C
      NIEL_SEND(:) = 0
      NIEL_RECV(:) = 0

      print *, 'Init neil_send/recv'

      DO I=1,NPROC
         print *, I
        DO J=1,NELP(I)
           ! namo - add ABS because of the change to ghost wier nodes in ADCIRC
         IELG = abs(IMAP_EL_LG(J,I))

         !if ((I.eq.2) .and. (J.gt.305)) print *, IELG
         !if ((I.eq.2) .and. (J.eq.313)) print *, 'IELG = ', IELG

         N1 = NNEG(1,IELG)
         N2 = NNEG(2,IELG)
         N3 = NNEG(3,IELG)

         !IF (I .gt. 1) print *, I, ': Finished N1,2,3'

         PE1 = IMAP_NOD_GL(1,N1)
         PE2 = IMAP_NOD_GL(1,N2)
         PE3 = IMAP_NOD_GL(1,N3)

         !if (I.eq.2) print *, J, '/', NELP(I)

! Element J in sub-domain I should receive information from
! a neighboring sub-domain (i.e. element J is a ghost element),
! if two or three nodes of the element belong to a sub-domain
! other than sub-domain I.
! However, there is an exception.  If the elemental nodes belong
! to three different sub-domains and one of the nodes belong to
! sub-domain I, this element is not a ghost element because
! it is immersed in sub-domain I since an adjacent element has
! been added at STEP 5-DG-3.

         IF((PE1.EQ.I.AND.PE2.EQ.I).OR.(PE2.EQ.I.AND.PE3.EQ.I).OR.
     &      (PE3.EQ.I.AND.PE1.EQ.I)) THEN
! A resident element; Do nothing
         ELSE IF((PE1.NE.PE2.AND.PE2.NE.PE3.AND.PE3.NE.PE1).AND.
     &    (PE1.EQ.I.OR.PE2.EQ.I.OR.PE3.EQ.I)) THEN
! An element shared by three different sub-domains; Do nothing
         ELSE IF(PE1.NE.PE2.AND.PE2.NE.PE3.AND.PE3.NE.PE1) THEN
! An element added at STEP 5-DG-3
            ! namo - skip this for now

! Ghost element J can receive information from any of the
! neighboring sub-domains if the elemental nodes belong
! to three different sub-domains and none of the sub-domains
! is sub-domain I.

            NIEL_RECV(I) = NIEL_RECV(I) + 1
            NIEL_SEND(PE1) = NIEL_SEND(PE1) + 1

         ELSE IF((PE1.NE.I.AND.PE2.NE.I).OR.
     &       (PE2.NE.I.AND.PE3.NE.I).OR.
     &       (PE3.NE.I.AND.PE1.NE.I)) THEN
! An element which comes here should have one node which belongs
! to sub-domain I and two nodes which belong to one sub-domain
! different from sub-domain I

! Ghost element J can receive information from a neighboring
! sub-domain to which two of the elemental nodes belong to the
! sub-domain.


            NIEL_RECV(I) = NIEL_RECV(I) + 1

            IF(PE1.EQ.PE2) THEN
               NIEL_SEND(PE1) = NIEL_SEND(PE1) + 1
            ELSE IF(PE2.EQ.PE3) THEN
               NIEL_SEND(PE2) = NIEL_SEND(PE2) + 1
            ELSE IF(PE3.EQ.PE1) THEN
               NIEL_SEND(PE3) = NIEL_SEND(PE3) + 1
            ELSE
               print *, PE1, PE2, PE3
               STOP 'STEP 9-DG-1 You should not see this. (1)'
            ENDIF
         ELSE
            STOP 'STEP 9-DG-1 You should not see this. (2)'
         ENDIF
      END DO

      ENDDO

      print *,' Finished step 9-DG-1'
!C STEP 9-DG-2:
!C--Alloc memory 2
!C
      CALL ALLOC_DG2()

      print *, 'Finished alloc_dg2'

C STEP 9-DG-n
! namo - declare communication using negative map approach

      NIEL_RECV(:) = 0
      NIEL_SEND(:) = 0

      IF (.false.) THEN
      DO I=1,NPROC
       DO J=1,NELP(I)
! if this element is a ghost element, communicate to the
! processor that owns it in EL_SHARE
          IELG = IMAP_EL_LG(J,I)
          IF (IELG.LT.0) THEN
             PE1 = EL_SHARE(abs(IELG))

             NIEL_RECV(I) = NIEL_RECV(I) + 1
             IEL_RECV(1,NIEL_RECV(I),I) = PE1
             IEL_RECV(2,NIEL_RECV(I),I) = J

             NIEL_SEND(PE1) = NIEL_SEND(PE1) + 1
             IEL_SEND(1,NIEL_SEND(PE1),PE1) = I
             IEL_SEND(2,NIEL_SEND(PE1),PE1) = -IELG
          END IF
       END DO
      END DO
      print *, 'NIEL_RECV', NIEL_RECV
      print *, 'NIEL_SEND', NIEL_SEND

      END IF
      IF (.true.) THEN
!C STEP 9-DG-3:
!C--Decompose sub-domain interface edges
!C
         DO I = 1,NPROC
            DO J=1, NELP(I)
! Zero out again
          ! namo - again add abs
         IELG = abs(IMAP_EL_LG(J,I))

         N1 = NNEG(1,IELG)
         N2 = NNEG(2,IELG)
         N3 = NNEG(3,IELG)

         PE1 = IMAP_NOD_GL(1,N1)
         PE2 = IMAP_NOD_GL(1,N2)
         PE3 = IMAP_NOD_GL(1,N3)

! Element J in sub-domain I should receive information from
! a neighboring sub-domain (i.e. element J is a ghost element),
! if two or three nodes of the element belong to a sub-domain
! other than sub-domain I.
! However, there is an exception.  If the elemental nodes belong
! to three different sub-domains and one of the nodes belong to
! sub-domain I, this element is not a ghost element because
! it is immersed in sub-domain I since an adjacent element has
! been added at STEP 5-DG-3.

         IF((PE1.EQ.I.AND.PE2.EQ.I).OR.(PE2.EQ.I.AND.PE3.EQ.I).OR.
     &    (PE3.EQ.I.AND.PE1.EQ.I)) THEN
! A resident element; Do nothing
       ELSE IF((PE1.NE.PE2.AND.PE2.NE.PE3.AND.PE3.NE.PE1).AND.
     &    (PE1.EQ.I.OR.PE2.EQ.I.OR.PE3.EQ.I)) THEN
! An element shared by three different sub-domains; Do nothing
      ELSE IF(PE1.NE.PE2.AND.PE2.NE.PE3.AND.PE3.NE.PE1) THEN
! An element added at STEP 5-DG-3
         ! namo - skip for now

! Ghost element J can receive information from any of the
! neighboring sub-domains if the elemental nodes belong
! to three different sub-domains and none of the sub-domains
! is sub-domain I.

         NIEL_RECV(I) = NIEL_RECV(I) + 1
         IEL_RECV(1,NIEL_RECV(I),I) = PE1
         IEL_RECV(2,NIEL_RECV(I),I) = J

         NIEL_SEND(PE1) = NIEL_SEND(PE1) + 1
         IEL_SEND(1,NIEL_SEND(PE1),PE1) = I
         IEL_SEND(2,NIEL_SEND(PE1),PE1) = IELG

      ELSE IF((PE1.NE.I.AND.PE2.NE.I).OR. (PE2.NE.I.AND.PE3.NE.I).OR.
     &     (PE3.NE.I.AND.PE1.NE.I)) THEN
! An element which comes here should have one node which belongs
! to sub-domain I and two nodes which belong to one sub-domain
! different from sub-domain I

! Ghost element J can receive information from a neighboring
! sub-domain to which two of the elemental nodes belong to the
! sub-domain.

         ! if so, mark that subdomain as the owner, overriding the EL_SHARE

         IF(PE1.EQ.PE2) THEN
            K = PE1
         ELSE IF(PE2.EQ.PE3) THEN
            K = PE2
         ELSE IF(PE3.EQ.PE1) THEN
            K = PE3
         ELSE
            STOP 'STEP 9-DG-1 You should not see this. (1)'
         ENDIF
         EL_SHARE(IELG) = K

         NIEL_RECV(I) = NIEL_RECV(I) + 1
         IEL_RECV(1,NIEL_RECV(I),I) = K
         IEL_RECV(2,NIEL_RECV(I),I) = J

         NIEL_SEND(K) = NIEL_SEND(K) + 1
         IEL_SEND(1,NIEL_SEND(K),K) = I
         IEL_SEND(2,NIEL_SEND(K),K) = IELG

      ELSE
         STOP 'STEP 9-DG-1 You should not see this. (2)'
      ENDIF
      END DO
      ENDDO

      print *, 'NIEL_SEND', NIEL_SEND
      print *, 'NIEL_RECV', NIEL_RECV

      END IF ! if false

!C STEP 9-DG-4:
!C--Replace the temporarily assigned global element IDs with local ones
!C
      DO I=1,NPROC
       DO J=1,NIEL_SEND(I)
         IELG = IEL_SEND(2,J,I)

         DO K=1,NELP(I)
            ! add abs
            IF(abs(IMAP_EL_LG(K,I)) .EQ. IELG) THEN
               IEL_SEND(2,J,I) = K
               cycle
            ENDIF
         enddo
      ENDDO
      ENDDO

!C STEP DG-7:
!C--Compute NUM_COMM_PE_SEND and NUM_COMM_PE_RECV
!C
      ALLOCATE( COMM_PE_FLAG(MNPROC) )

      NUM_COMM_PE_RECV(:) = 0
      NUM_COMM_PE_SEND(:) = 0

      DO I=1,MNPROC
      COMM_PE_FLAG(:) = 0

      DO J=1,NIEL_RECV(I)
          COMM_PE_FLAG(IEL_RECV(1,J,I)) = 1
      ENDDO

      DO J=1,MNPROC
          IF(COMM_PE_FLAG(J).EQ.1) THEN
             NUM_COMM_PE_RECV(I) = NUM_COMM_PE_RECV(I) + 1
          ENDIF
      ENDDO

      COMM_PE_FLAG(:) = 0

      DO J=1,NIEL_SEND(I)
          COMM_PE_FLAG(IEL_SEND(1,J,I)) = 1
      ENDDO

      DO J=1,MNPROC
          IF(COMM_PE_FLAG(J).EQ.1) THEN
             NUM_COMM_PE_SEND(I) = NUM_COMM_PE_SEND(I) + 1
          ENDIF
      ENDDO

      ENDDO

      DEALLOCATE( COMM_PE_FLAG )
      print *, 'passed deallocate'
C STEP 10:
C--Print Summary of Decomposition
C
      print *, "Decomposition Data"
      print *, "DOMAIN  RES_NODES  GHOST_NODES  TOT_NODES  ELEMENTS"
      print *, "------  ---------  -----------  ---------  --------"
      WRITE(*,90) STRINGGLOBAL,NNODG,NELG
      DO I=1, NPROC
         PE(1:6) = 'PE0000'
         WRITE(PE(3:6),'(I4.4)') I-1
         WRITE(6,92) PE, NOD_RES_TOT(I),NNODP(I)-NOD_RES_TOT(I),
     &                   NNODP(I),NELP(I)
      ENDDO
 90   FORMAT(1X,A6,25X,I9,2X,I9)
 92   FORMAT(1X,A6,1X,I9,2X,I9,4X,I9,2X,I9)
C

      print *, 'Finished decomp inside'
      RETURN
C----------------------------------------------------------------------
      END SUBROUTINE DECOMP_DG

C----------------------------------------------------------------------
C----------------------------------------------------------------------
      subroutine checkMappings(partitionOK,numPartitionAttempts)
      use pre_global
      implicit none
      logical, intent(out) :: partitionOK 
      integer, intent(in) :: numPartitionAttempts
      integer i, ilnode, ilnode2, ilnode3, index2, indx, iproc, iproc2
      integer iproc3, itemp, itemp2, j, jd, k
      integer :: hbn ! hanging boundary node
      !
      ! NVELLP(K,PE) Number of Land Boundary Nodes of Segment K on PE
      ! IMAP_NOD_GL2(2(PE-1)+1,I)  = PE assigned to Global Node I
      ! IMAP_NOD_GL2(2(PE-1)+2,I)  = Local Node Number of Global Node I on PE
      ! vjp modified array to drop last dimension to save memory space
      ! LBINDEX_LG(K,I,PE) = Global Index of I-th Node on Land Boundary Segment K on PE
      partitionOK = .true.
      ! 
      ! loop over subdomains
      do iproc = 1,nproc
         nvelp(iproc) = 0
         nvellp(:) = 0   ! initialize number of nodes on each fulldomain boundary
                         ! to zero for this subdomain
         !
         ! loop over full domain boundaries
         do k = 1,nbou
            select case(ibtype(k))
            ! weir boundaries
            case(4,24,64,5,25)
               do i = 1,nvell(k)
                  indx = nbvv(k,i)
                  index2 = ibconnr(k,i)
                  do j = 1,itotproc(indx)
                     itemp = (j-1)*2 + 1
                     iproc2  =  imap_nod_gl2(itemp,indx)
                     ilnode2 =  imap_nod_gl2(itemp+1,indx)
                     if (iproc.eq.iproc2) then
                        do jd = 1, itotproc(index2)
                           itemp2 = (jd-1)*2 + 1
                           iproc3  = imap_nod_gl2(itemp2,index2)
                           ilnode3 = imap_nod_gl2(itemp2+1,index2)
                           if (iproc.eq.iproc3) then
                              nvelp(iproc) = nvelp(iproc) + 1
                              nvellp(k) = nvellp(k) + 1
                              lbindex_lg(k,nvellp(k)) = i
                              nbvvp(k,nvellp(k))   = ilnode2
                              ibconnrp(k,nvellp(k)) = ilnode3
                           endif
                        enddo
                     endif
                  enddo
               enddo
            case default
               ! for each node in this non-weir flux boundary 
               do i = 1,nvell(k)
                  ! find full domain node number
                  indx = nbvv(k,i)
                  ! iterate over the subdomains where this boundary
                  ! node is found
                  do j = 1,itotproc(indx)
                     itemp = (j-1)*2 + 1
                     ! look up subdomain number 
                     iproc2 =  imap_nod_gl2(itemp,indx)
                     ! look up node number in this subdomain                    
                     ilnode =  imap_nod_gl2(itemp+1,indx)
                     ! if the current subdomain matches this one
                     if (iproc.eq.iproc2) then
                        ! increment the total number of flux boundary nodes
                        ! on this subdomain
                        nvelp(iproc) = nvelp(iproc) + 1
                        ! increment number of nodes from this boundary
                        ! that occur in this subdomain
                        nvellp(k) = nvellp(k) + 1
                        ! add fulldomain node number for this boundary 
                        ! and subdomain node array position on this boundary
                        ! to mapping
                        lbindex_lg(k,nvellp(k)) = i
                        ! record subdomain node number for this boundary
                        ! and subdomain boundary position  
                        nbvvp(k,nvellp(k)) = ilnode
                     endif
                  enddo
               enddo
               !
               ! @jasonfleming: it is possible for a single node of an 
               ! island boundary to appear as a ghost node in a nearby
               ! subdomain, and if this single node is the first node
               ! on the fulldomain island boundary, it will also be the 
               ! last node on the fulldomain boundary (because adcirc
               ! requires island boundaries to be closed). As a result,
               ! it can be erroneously translated as a two node subdomain
               ! boundary with the same node number listed twice. 
               if ( nvellp(k).eq.2 ) then
                  if ( nbvvp(k,1).eq.nbvvp(k,2) ) then
                     write(*,'(a)') 
     &           'ERROR: Two node boundary with identical node numbers.'
                     hbn = imap_nod_lg(nbvvp(k,1),iproc)
                     boundaryWeights(hbn) = boundaryWeights(hbn) 
     &                  + 1
                     partitionOK = .false.
                  endif
               endif
C
            end select
C
         enddo  ! end loop over flux specified boundaries
         !
         ! @jasonfleming: it is possible for a single node of an 
         ! land boundary to appear as a ghost node in a nearby
         ! subdomain. 
         do k=1, nbou
            if (nvellp(k).eq.1) then
               write(*,'(a,i0,a,i0,a)') 
     &           'WARNING: The land boundary number ',k,
     &           ' is only one node long in subdomain ',iproc,'.'
               !hbn = imap_nod_lg(nbvvp(k,1),iproc)
               !boundaryWeights(hbn) = boundaryWeights(hbn) + 1
               ! @jasonfleming: large meshes decomposed on many
               ! subdomains may never be able to avoid this without
               ! smarter adjustments of the partition weights, so
               ! the rejection of the partitioning as a result of
               ! this condition is disabled for now. 
               !select case (ibtype(k))
               !case(0,20,30) ! zero flux land boundaries
               !   partitionOK = .false. ! allow discrepancies here ?
               !case default
               !   partitionOK = .false.
               !end select 
            endif
         enddo

      end do ! end loop over subdomains

C----------------------------------------------------------------------
      end subroutine checkMappings
C----------------------------------------------------------------------
      SUBROUTINE DOMSIZE_DG()
      USE pre_global

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
         NODE_GL1(J) = PROC(J) ! namo -sb
         IF (IPROC.EQ.PROC(J)) THEN
           RESNODE = RESNODE+1
           !NODE_GL1(J) = PROC(J)
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
      IF (ITOT.GT.VTMAX) THEN
          WRITE(*,'(A)') "Error during step3 decomp"
          CALL EXIT(1)
      ENDIF
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
        IF (ITOT.GT.VTMAX) THEN
            WRITE(*,'(A)') "Error during step 4 decomp"
            CALL EXIT(1)
        ENDIF
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
        CALL EXIT(1)
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

      ! namo added 2022
      print *, 'Reached step 5 dg 3'
Csb-DG2
C STEP 5-DG-3:
C--If there are any elements shared by three different sub-domains,
C  Re-construct Local-to-Global Element Map: IMAP_EL_LG(:,PE)
C  Rule: Add non-local elements adjacent to the elements
C        shared by three different sub-domains
C        to the local element list
C
      IF (.true.) then
      
      IF (NELEM.eq.0) NELEM = ONELEM
      !NELEM = ONELEM
      ONELEM = NELEM
      NEL_ADDED = 0
      DO J = 1,ONELEM
        K = ELEM_LG(J)

        N1 = NNEG(1,K)
        N2 = NNEG(2,K)
        N3 = NNEG(3,K)

        PE1 = NODE_GL1(N1)
        PE2 = NODE_GL1(N2)
        PE3 = NODE_GL1(N3)

        IF(K.EQ.10352.AND.IPROC.EQ.16) THEN
           print *, 'B=',IPROC,NELEM,K,PE1,PE2,PE3
        ENDIF
        IF(((PE1.NE.PE2).AND.(PE2.NE.PE3).AND.(PE3.NE.PE1)).AND.
     &     (PE1.EQ.IPROC.OR.PE2.EQ.IPROC.OR.PE3.EQ.IPROC)) THEN
C        IF((PE1.NE.PE2).AND.(PE2.NE.PE3).AND.(PE3.NE.PE1)) THEN
          IF(NELEM.GE.MNE) THEN
            STOP 'FATAL ERROR (STEP 5-DG-3 DRY)'
          ELSE
             !print *, 'incrementing nelem'
            NELEM = NELEM + 1
            ELEM_LG(NELEM) = K
            NEL_ADDED = NEL_ADDED + 1
          ENDIF
        ENDIF
      ENDDO
      print *, 'Finished dg dry'
      endif
C--
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
      IF (ITOT.GT.VTMAX) THEN
          WRITE(*,'(A)') "Error during step6 decomp"
          CALL EXIT(1)
      ENDIF
      CALL SORT(ITEMP,ITVECT) ! Sort and remove multiple occurrences
      ITOT = 1
      NODE_LG(1) = ITVECT(1)
      DO J = 2, ITEMP
         IF (ITVECT(J).NE.NODE_LG(ITOT)) THEN
           ITOT = ITOT + 1
           NODE_LG(ITOT) = ITVECT(J)
         ENDIF
      ENDDO

C namo added 2022
Csb
C STEP 6-DG
C--Increase ITOT to take into acount NEL_ADDED
C
      ITOT = ITOT + NEL_ADDED*3
C-----------------------------------------------

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
C----------------------------------------------------------------------
      END SUBROUTINE DOMSIZE_DG

C
C---------------------------------------------------------------------------C
C                   (  Serial Version 1.0  12/20/99 vjp )                   C
C                                                                           C
C  Takes dry run through the domain decomp to determine the max number of   C
C  nodes and elements assigned to any subdomain to determine MNPP and MNEP. C
C                                                                           C
C---------------------------------------------------------------------------C
C
      SUBROUTINE DOMSIZE()
      USE pre_global

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
      IF (ITOT.GT.VTMAX) THEN
          WRITE(*,'(A)') "Error during step3 decomp"
          CALL EXIT(1)
      ENDIF
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
        IF (ITOT.GT.VTMAX) THEN
            WRITE(*,'(A)') "Error during step 4 decomp"
            CALL EXIT(1)  
        ENDIF
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
        CALL EXIT(1)
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
      IF (ITOT.GT.VTMAX) THEN
          WRITE(*,'(A)') "Error during step6 decomp"
          CALL EXIT(1)
      ENDIF
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
C----------------------------------------------------------------------
      END SUBROUTINE DOMSIZE
C----------------------------------------------------------------------
C
C---------------------------------------------------------------------------
C  Sorts array RA of length N into ascending order using Heapsort algorithm.
C  N is input; RA is replaced on its output by its sorted rearrangement.
C  Ref: Numerical Recipes
C---------------------------------------------------------------------------
C
      SUBROUTINE SORT(N,RA)
      IMPLICIT NONE
      INTEGER N, L, IR, RRA, I, J
      INTEGER RA(N)

      L = N/2 + 1
      IR = N
      DO WHILE(.TRUE.)
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
        DO WHILE (J.LE.IR)
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
        ENDDO
        RA(I)=RRA
      ENDDO
C----------------------------------------------------------------------
      END SUBROUTINE SORT
C----------------------------------------------------------------------

C----------------------------------------------------------------------
C--Given an array XX of length N, and given a value X, returns a value J
C--such that X is between XX(J) and XX(J+1). XX must be monotonic, either
C--increasing or decreasing. J=0 or J=N is returned to indicate that X is out
C--of range.
C--
C--NUMERICAL RECIPES - The Art of Scientific Computing [FORTRAN Version]
C----------------------------------------------------------------------
      SUBROUTINE LOCATE(XX,N,X,J)
      IMPLICIT NONE
      INTEGER JM,JL,JU,J,N,X,XX(N)
C

C--Initialize lower and upper limits
         JL = 0
         JU = N+1
C
C--If we are not done yet, compute a mid-point, and replace either the lower
C--limit or the upper limit, as appropriate.
C
         DO WHILE (JU-JL.GT.1)
            JM = (JU+JL)/2
            IF((XX(N).GT.XX(1)).EQV.(X.GT.XX(JM)))THEN
              JL = JM
            ELSE
              JU = JM
            ENDIF
         ENDDO
C--Then set the output and return.
         J = JL
         RETURN
C----------------------------------------------------------------------
         END SUBROUTINE LOCATE
C----------------------------------------------------------------------


C----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE SEARCH3(MAP,LENG,N1,N2,N3,INDX)
      INTEGER MAP(*),LENG,N1,N2,N3,INDX,IP
cvjp  rewritten 5/3/99
      INDX = 0
      DO I=1,LENG
         IP = MAP(I)
         IF (IP.EQ.N1.OR.IP.EQ.N2.OR.IP.EQ.N3) THEN
           INDX = I
           EXIT
         ENDIF
      ENDDO
      RETURN
C----------------------------------------------------------------------
      END SUBROUTINE SEARCH3
C----------------------------------------------------------------------


C----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE SEARCH(MAP,N,TARGET,INDX)
      INTEGER MAP(*),N,TARGET,INDX,I
C
      INDX = 0
      IF (N.EQ.0) RETURN
      DO I=1,N
         IF (MAP(I).EQ.TARGET) THEN
            INDX = I
            EXIT
         ENDIF
      ENDDO
      RETURN
C----------------------------------------------------------------------
      END SUBROUTINE SEARCH
C----------------------------------------------------------------------
