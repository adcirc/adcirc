C******************************************************************************
C (P)ADCIRC VERSION 45.12 03/17/2006                                          *
C  last changes in this file VERSION 45.12                                    *
C                                                                             *
C                                                                             *
C******************************************************************************
C                                                                             *
C                           (P)ADCIRC                                         *
C                                                                             *
C    A (PARALLEL) ADVANCED CIRCULATION MODEL FOR SHELVES, COASTAL SEAS        *
C                         AND ESTUARIES                                       *
C                                                                             *
C                                                                             *
C                          DEVELOPED BY:                                      *
C                                                                             *
C                      DR. R.A. LUETTICH, JR                                  *
C                                                                             *
C             UNIVERSITY OF NORTH CAROLINA AT CHAPEL HILL                     *
C                   INSTITUTE OF MARINE SCIENCES                              *
C                       3431 ARENDELL ST.                                     *
C                   MOREHEAD CITY, NC, 28557                                  *
C                         252-726-6841 EXT. 137                               *
C                   EMAIL  RICK_LUETTICH@UNC.EDU                              *
C                                                                             *
C                        DR. J.J. WESTERINK                                   *
C                                                                             *
C          DEPARTMENT OF CIVIL ENGINEERING AND GEOLOGICAL SCIENCES            *
C                     UNIVERSITY OF NOTRE DAME                                *
C                       NOTRE DAME, IN 46556                                  *
C                         219-631-6475                                        *
C                    EMAIL JJW@PHOTIUS.CE.ND.EDU                              *
C                                                                             *
C                                                                             *
C        MAJOR FUNDING FOR THE DEVELOPMENT OF ADCIRC WAS PROVIDED BY          *
C                                                                             *
C                       DEPARTMENT OF THE ARMY                                *
C                    WATERWAYS EXPERIMENT STATION                             *
C                 COASTAL ENGINEERING RESEARCH CENTER                         *
C                        3909 HALLS FERRY RD                                  *
C                      VICKSBURG, MI 39180-6199                               *
C                                                                             *
C******************************************************************************
C                                                                             *
C          THE ADCIRC SOURCE CODE IS COPYRIGHTED, 1994-2006 BY:               *
C                                                                             *
C                 R.A. LUETTICH, JR AND J.J. WESTERINK                        *
C                                                                             *
C         NO PART OF THIS CODE MAY BE REPRODUCED OR REDISTRIBUTED             *
C               WITHOUT THE WRITTEN PERMISSION OF THE AUTHORS                 *
C                                                                             *
C******************************************************************************
C                                                                             *
C                  Parallelization of ADCIRC 2D was done by                   *
C                    the center for subsurface modeling                       *
C                         The University of Texas                             *
C                             Austin, TX 78712                                *
C                            03/06/98 - 5/19/99                               *
C                                                                             *
C                   Translation of 2D code to Fortran90                       *
C                             Victor J. Parr                                  *
C                              John B. Romo                                   *
C                               8/31/99                                       *
C                                                                             *
C               Parallelization, consolidation of 2D & 3D codes,              *
C             translation of 3D routines to Fortran 90 was done by            *
C                               Tim Campbell                                  *
C                   Naval Research Lab, Stennis Space Center                  *
C                             summer 2002                                     *
C                                                                             *
C******************************************************************************
 
C
      PROGRAM PADCIRC
C
      USE GLOBAL
      USE GLOBAL_3DVS, ONLY : GORho, GORhoOAMB, AMB
      USE HARM
      USE VERSION

#ifdef CMPI
      USE MESSENGER
      USE WRITER, ONLY: WRITER_MAIN, WRITER_INIT
#endif

      IMPLICIT NONE
      INTEGER IE, I, J, N    ! local loop counter
      INTEGER NM1, NM2, NM3  ! local node numbers
      INTEGER NCELE          ! local node code indicator

      REAL(8) AreaEle
 
#ifdef CMPI
      CALL GET_NUMWRITERS()  ! Get number of writer processors from command line arguments
      CALL MSG_INIT()        ! Init MPI and get MPI-rank of this cpu
      ! <ezpp-begin>
      FileFmtVersion = VERSION_NUMBER(FileFmtMajor, FileFmtMinor, FileFmtRev)
      ! File Fmt Version  must match Major and Minor numbers.
      if (MYPROC == 0) print 1999, "PARALLEL", ADC_VERSION

C...  Special instructions for writer processors
      IF(MYPROC >= MNPROC) THEN
C...    02/07/2007 - sb - Writer processors will enter "WRITER_MAIN" subroutine 
C...    and will be terminated before they exit this subroutine.  For this reason,
C...    The "STOP" statement in the next line to "CALL WRITER_MAIN()" will NOT
C...    have any effect, but it is placed there to make it apparent that
C...    writer processors will not execute the rest in this source file (adcirc.F).
        CALL WRITER_MAIN()
        STOP
      ENDIF

      CALL WRITER_INIT()     ! Initialize WRITER module
      CALL MAKE_DIRNAME()    ! Establish Working Directory Name
      CALL READ_INPUT()      ! Establish sizes by reading fort.14 and fort.15
      CALL MSG_TABLE()       ! Read Message-Passing Tables
      CALL MSG_START()       ! Startup message passing
#else
      MNPROC = 1             ! Init number of compute procs for serial
      MNWPROC = 0            ! Init number of writer procs for serial
      MNALLPROC = 1          ! Init number of all procs for serial
      MYPROC = 0             ! Init task id for serial
      FileFmtVersion = VERSION_NUMBER(FileFmtMajor, FileFmtMinor, FileFmtRev)
      ! File Fmt Version  must match Major and Minor numbers.
      print 1999, "",ADC_VERSION
      CALL MAKE_DIRNAME()    ! Establish Working Directory Name
      CALL READ_INPUT()      ! Establish sizes by reading fort.14 and fort.15
#endif         


#ifdef CMACHSUN
      CALL ABRUPT_UNDERFLOW() ! Abrupt underflow for Sun computers
#endif


C...Compute the reciprocal of the number of nodes in the entire domain 

#ifdef CMPI
      CALL ALLNODES(NP_GLOBAL)
      RNP_GLOBAL = 1.0D0/NP_GLOBAL      
#endif

C     jgf45.08 put this here so it can be used in hot start
      IF(NOLIFA.EQ.0) THEN
         IFNLFA=0
      ELSE
         IFNLFA=1
      ENDIF

C...
C...******************** PROGRAM SETUP SECTION ************************
C...
      IF (IHOT.EQ.0) THEN
         CALL COLDSTART()
      ELSE
         CALL HOTSTART()
      ENDIF


C...Determine the number of active elements (MJU), the total number of 
C...elements (NODELE) and the total area of elements (TotalArea) attached to each node

      DO I=1,NP
         MJU(I)=0
         NODELE(I)=0
         TotalArea(I)=0.d0
         NODECODE(I)=NNODECODE(I)
      END DO

      DO IE=1,NE
         NM1=NM(IE,1)
         NM2=NM(IE,2)
         NM3=NM(IE,3)
         NCELE=NODECODE(NM1)*NODECODE(NM2)*NODECODE(NM3)*NOFF(IE) 
         MJU(NM1)=MJU(NM1)+NCELE
         MJU(NM2)=MJU(NM2)+NCELE
         MJU(NM3)=MJU(NM3)+NCELE
         NODELE(NM1)=NODELE(NM1)+1
         NODELE(NM2)=NODELE(NM2)+1
         NODELE(NM3)=NODELE(NM3)+1
         AreaEle=NCEle*Areas(IE)/2.d0
         TotalArea(NM1)=TotalArea(NM1)+AreaEle
         TotalArea(NM2)=TotalArea(NM2)+AreaEle
         TotalArea(NM3)=TotalArea(NM3)+AreaEle
      END DO

      DO I=1,NP
         IF(MJU(I).EQ.0) MJU(I)=1
      END DO


C...
C...************* SET FLAGS AND COEFFICIENTS USED IN TIME STEPPING ***********
C...

C...Set flags for nonlinear terms

      IF(NOLICA.EQ.0) THEN
         IFNLCT=0
      ELSE
         IFNLCT=1
      ENDIF
      IF(NOLICAT.EQ.0) THEN
         IFNLCAT=0
      ELSE
         IFNLCAT=1
      ENDIF

      IFWIND=1
      IF(IM.EQ.1) IFWIND=0

C...Compute several constant coefficients

      GA00=G*A00
      GC00=G*C00
      TADVODT=IFNLCAT/DT
      GB00A00=G*(B00+A00)
      GFAO2=G*IFNLFA/2.D0
      GO3=G/3.D0
      GO2=G/2.D0
      DTO2=DT/2.D0
      DT2=DT*2.D0
      GDTO2=G*DT/2.D0
      SADVDTO3=IFNLCT*DT/3.D0
      GORho=G/RhoWat0
      GORhoOAMB=GORho/AMB

C...
C.....Compute necesary wet/dry coefficients 
C...
      HABSMIN=0.8d0*H0  
      HOFF=1.2d0*H0    

C...
C...************************* BEGIN TIME STEPPING *************************
C...  
      WRITE(16,1112)
      WRITE(16,17931)
      IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,1112)
      IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,17931)
 1112 FORMAT(/,1X,79('_'))
17931 FORMAT(//,1X,'LIMITED RUNTIME INFORMATION SECTION ',//)
      
      NCChange=1 !jgf45.06 set flag to set up GWCE LHS on 1st iteration

      DO 100 ITIME=ITHS+1,NT
C
        CALL TIMESTEP(ITIME)
C
 100    CONTINUE

C...
C...****************** TIME STEPPING LOOP ENDS HERE ********************
C...

C...
C...*************** SOLVE THE HARMONIC ANALYSIS PROBLEM ****************
C...

      IF ((IHARIND.EQ.1).AND.(ITIME.GT.ITHAS)) THEN

C...Compute means and variances for checking the harmonic analysis results
C...Accumulate mean and variance at each node.

        if (CHARMV) then
          IF (FMV.NE.0.) THEN
            DO I=1,NP
              ELAV(I)   = ELAV(I)/NTSTEPS
              XVELAV(I) = XVELAV(I)/NTSTEPS
              YVELAV(I) = YVELAV(I)/NTSTEPS
              ELVA(I)   = ELVA(I)/NTSTEPS   - ELAV(I)*ELAV(I)
              XVELVA(I) = XVELVA(I)/NTSTEPS - XVELAV(I)*XVELAV(I)
              YVELVA(I) = YVELVA(I)/NTSTEPS - YVELAV(I)*YVELAV(I)
              END DO
            TIMEBEG=ITMV*DTDP + (STATIM-REFTIM)*86400.D0
            OPEN(55,FILE=TRIM(LOCALDIR)//'/'//'fort.55')
            WRITE(55,*) NP
            ENDIF
          endif

C......Fill out and decompose the LHS harmonic analaysis matrix

        CALL FULSOL(0)

C......Solve the harmonic analysis problem and write the output

        IF(NHAGE.EQ.1) CALL LSQSOLEG(NP,LOCALDIR,ELAV,ELVA,CHARMV)
C     
        IF(NHAGV.EQ.1) CALL LSQSOLVG(NP,LOCALDIR,
     &                          XVELAV,YVELAV,XVELVA,YVELVA,CHARMV)
C
        IF(NHASE.EQ.1) CALL LSQSOLES(NSTAE,LOCALDIR)
C     
        IF(NHASV.EQ.1) CALL LSQSOLVS(NSTAV,LOCALDIR)
C
        ENDIF

#ifdef CMPI
      CALL MSG_FINI()
#endif
C
C
      STOP
 1999 format(1X, a," ADCIRC Version: ",a)
      END


C******************************************************************************
C                                                                             *
C      Subroutine to generate neighbor tables from a connectivity table.      *
C                                                                             *
C      NOTES                                                                  *
C      a node neighbor table is generated with the node itself is listed as   *
C         neighbor #1 and all other neighbors are sorted and placed in cw     *
C         order from east                                                     *
C      a neighbor element table is generated with:                            *
C         entry 1 = element # defined by neighbors 1,2,3                      *
C         entry 2 = element # defined by neighbors 1,3,4                      *
C         entry 3 = element # defined by neighbors 1,4,5                      *
C          .......                                                            *
C         entry last = element # defined by neighbors 1,nneigh,2              *
C         a zero area means that the defined triangle lies outside the domain *
C                                                                             *
C                                                                             *
C    v1.0   R.L.   6/29/99  used in 3D code                                   *
C    v2.0   R.L.   5/23/02  adapted to provide neighbor el table              *
C******************************************************************************
C                                                                             *
C     -  PARAMETERS WHICH MUST BE SET TO CONTROL THE DIMENSIONING OF ARRAYS   *
C           ARE AS FOLLOWS:                                                   *
C                                                                             *
C          MNP = MAXIMUM NUMBER OF NODAL POINTS                               *
C          MNE = MAXIMUM NUMBER OF ELEMENTS                                   *
C          MNEI= 1+MAXIMUM NUMBER OF NODES CONNECTED TO ANY ONE NODE IN THE   *
C                   FINITE ELEMENT GRID                                       *
C                                                                             *
C******************************************************************************
C                                                                             *
C    VARIABLE DEFINITIONS:                                                    *
C       NE - NUMBER OF ELEMENTS                                               *
C       NP - NUMBER OF NODES                                                  *
C       NM(MNE,3) - NODE NUMBERS ASSOCIATED WITH EACH ELEMENT                 *
C       NNeigh(MNP) NUMBER OF NEIGHBORS FOR EACH NODE                         *
C       NeiTab(MNP,NEIMAX) 2D ARRAY OF NEIGHBORS FOR EACH NODE                *
C       NeiTabEle(MNP,NEIMAX) 2D ARRAY OF NEIGHBOR ELEMENTS FOR EACH NODE     *
C       NEIMIN - 1+MINIMUM NUMBER OF NEIGHBORS FOR ANY NODE                   *
C       NEIMAX - 1+MAXIMUM NUMBER OF NEIGHBORS FOR ANY NODE                   *
C                                                                             *
C******************************************************************************
C

      SUBROUTINE NEIGHB(NE,NP,NM,NNeigh,NeiTab,NeiTabEle,NEIMIN,NEIMAX,
     &     X,Y,NSCREEN,ScreenUnit)

      USE SIZES
      IMPLICIT NONE
      INTEGER :: NP,NE,NEIMIN,NEIMAX,NSCREEN,N,NN,I,J,JJ,K,JLOW
      INTEGER :: ScreenUnit
      INTEGER :: NN1,NN2,NN3,NE1,NE2,NE3
      INTEGER :: NeiTab(MNP,MNEI), NNeigh(MNP), NeiTabEle(MNP,MNEI)
      INTEGER :: NM(MNE,3)
      REAL(8) :: X(MNP),Y(MNP),DELX,DELY,DIST
      REAL(8) :: ANGLELOW,ANGLEMORE,RAD2DEG
      REAL(8), ALLOCATABLE :: ANGLE(:)
      INTEGER, ALLOCATABLE :: NEITEM(:)
      INTEGER, ALLOCATABLE :: NNEIGHELE(:)

C
      ALLOCATE ( ANGLE(MNEI) )
      ALLOCATE ( NEITEM(MNP) )
      ALLOCATE ( NNeighEle(MNP) )
C     
      RAD2DEG=45.0d0/ATAN(1.0d0)
C     
      DO N=1,NP
         NNeigh(N)=0
         NNeighEle(N)=0
         DO NN=1,MNEI
            NeiTab(N,NN)=0
            NeiTabEle(N,NN)=0
         END DO
      END DO

      DO 10 N=1,NE
         NN1 = NM(N,1)
         NN2 = NM(N,2)
         NN3 = NM(N,3)

         NNeighEle(NN1)=NNeighEle(NN1)+1
         NNeighEle(NN2)=NNeighEle(NN2)+1
         NNeighEle(NN3)=NNeighEle(NN3)+1
         NeiTabEle(NN1,NNeighEle(NN1))=N
         NeiTabEle(NN2,NNeighEle(NN2))=N
         NeiTabEle(NN3,NNeighEle(NN3))=N

         DO J=1,NNeigh(NN1)
            IF(NN2.EQ.NeiTab(NN1,J)) GOTO 25
         END DO
         NNeigh(NN1)=NNeigh(NN1)+1
         NNeigh(NN2)=NNeigh(NN2)+1
         IF((NNeigh(NN1).GT.MNEI-1).OR.(NNeigh(NN2).GT.MNEI-1)) GOTO 999
         NeiTab(NN1,NNeigh(NN1))=NN2
         NeiTab(NN2,NNeigh(NN2))=NN1

 25      CONTINUE
         DO J=1,NNeigh(NN1)
            IF(NN3.EQ.NeiTab(NN1,J)) GOTO 35
         END DO
         NNeigh(NN1)=NNeigh(NN1)+1
         NNeigh(NN3)=NNeigh(NN3)+1
         IF((NNeigh(NN1).GT.MNEI-1).OR.(NNeigh(NN3).GT.MNEI-1)) GOTO 999
         NeiTab(NN1,NNeigh(NN1))=NN3
         NeiTab(NN3,NNeigh(NN3))=NN1
         
 35      CONTINUE
         DO J=1,NNeigh(NN2)
            IF(NN3.EQ.NeiTab(NN2,J)) GOTO 10
         END DO
         NNeigh(NN2)=NNeigh(NN2)+1
         NNeigh(NN3)=NNeigh(NN3)+1
         IF((NNeigh(NN2).GT.MNEI-1).OR.(NNeigh(NN3).GT.MNEI-1)) GOTO 999
         NeiTab(NN2,NNeigh(NN2))=NN3
         NeiTab(NN3,NNeigh(NN3))=NN2

 10   CONTINUE

C
C     INSERT NODE ITSELF IN PLACE #1 and SORT other NEIGHBORS by
C     increasing cw angle from East
C
      DO I=1,NP
         DO J=1,NNeigh(I)
            NEITEM(J)=NeiTab(I,J)
            DELX=X(NEITEM(J))-X(I)
            DELY=Y(NEITEM(J))-Y(I)
            DIST=SQRT(DELX*DELX+DELY*DELY)
            IF(DIST.EQ.0.0d0) GOTO 998
            IF(DELY.NE.0.0d0) THEN
               ANGLE(J)=RAD2DEG*ACOS(DELX/DIST)
               IF(DELY.GT.0.0) ANGLE(J)=360.0d0-ANGLE(J)
            ENDIF
            IF(DELY.EQ.0.0d0) THEN
               IF(DELX.GT.0.0d0) ANGLE(J)=0.0d0
               IF(DELX.LT.0.d0) ANGLE(J)=180.0d0
            ENDIF
         END DO
         ANGLEMORE=-1.d0
         DO JJ=1,NNeigh(I)
            ANGLELOW=400.d0
            DO J=1,NNeigh(I)
               IF((ANGLE(J).LT.ANGLELOW).AND.
     &              (ANGLE(J).GT.ANGLEMORE)) THEN
                  ANGLELOW=ANGLE(J)
                  JLOW=J
               ENDIF
            END DO
            NeiTab(I,JJ+1)=NEITEM(JLOW)
            ANGLEMORE=ANGLELOW
         END DO
         NeiTab(I,1)=I
         NNeigh(I)=NNeigh(I)+1
      ENDDO
C
C     MATCH EACH SET OF 3 NODES WITH CORRESPONDING ELEMENT AND REORDER
C     ELEMENTS ACCORDINGLY
C
      DO I=1,NP
         DO K=1,NNeighEle(I)
            NEITEM(K)=NeiTabEle(I,K)
            NeiTabEle(I,K)=0
         END DO
         DO J=2,NNeigh(I)
            NN1=NeiTab(I,1)
            NN3=NeiTab(I,J)
            IF(J.NE.NNeigh(I)) NN2=NeiTab(I,J+1)
            IF(J.EQ.NNeigh(I)) NN2=NeiTab(I,2)
            DO K=1,NNeighEle(I)
               IF(NEITEM(K).NE.0) THEN
                  IF(NM(NEITEM(K),1).EQ.NN1) THEN
                     NE1=NM(NEITEM(K),1)
                     NE2=NM(NEITEM(K),2)
                     NE3=NM(NEITEM(K),3)
                  ENDIF
                  IF(NM(NEITEM(K),2).EQ.NN1) THEN
                     NE1=NM(NEITEM(K),2)
                     NE2=NM(NEITEM(K),3)
                     NE3=NM(NEITEM(K),1)
                  ENDIF
                  IF(NM(NEITEM(K),3).EQ.NN1) THEN
                     NE1=NM(NEITEM(K),3)
                     NE2=NM(NEITEM(K),1)
                     NE3=NM(NEITEM(K),2)
                  ENDIF
                  IF((NE2.EQ.NN2).AND.(NE3.EQ.NN3)) THEN
                     NeiTabEle(I,J-1)=NEITEM(K)
                     NEITEM(K)=0
                  ENDIF
               ENDIF
            END DO
         END DO
      END DO

C
C  DETERMINE THE MAXIMUM AND MINIMUM NUMBER OF NEIGHBORS
C
      NEIMAX = 0
      NEIMIN = 1000
      DO N=1,NP
         IF(NNeigh(N).LT.NEIMIN) NEIMIN=NNeigh(N)
         IF(NNeigh(N).GT.NEIMAX) NEIMAX=NNeigh(N)
      END DO

C
C  WRITE OUT DIAGNOSTIC OUTPUT  
C
C     OPEN(333,file='fort.333')
C     DO N=1,NP
C       WRITE(333,331) (NEIGH(N,J),J=1,NNEIGH(N))
C       WRITE(333,331) N,(NEIGHELE(N,J),J=1,NNEIGH(N)-1)
C       WRITE(333,*) ' '
C331    FORMAT(15(1X,I7))
C       END DO
C     CLOSE (333)

C  Deallocate local work arrays

      DEALLOCATE ( ANGLE )     
      DEALLOCATE ( NEITEM )    
      DEALLOCATE ( NNEIGHELE )

C  DONE

      RETURN
C
C  TERMINATE PROGRAM IF MAXIMUM NUMBER OF NEIGHBORS SET TOO SMALL
C
 999  CONTINUE
      IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) WRITE(ScreenUnit,99311)
      WRITE(16,99311)
99311 FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',
     &      //,1X,'THE DIMENSIONING PARAMETER MNEI IS TOO SMALL',
     &     /,1X,'THERE IS A PROBLEM WITH THE DYNAMIC MEMORY ALLOCATION',
     &     //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
      STOP

 998  CONTINUE
      IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
         WRITE(ScreenUnit,99312) I,NEITEM(J)
      ENDIF
      WRITE(16,99312) I,NEITEM(J)
99312 FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',
     &      //,1X,'NODES ',I7,' AND ',I7,' HAVE THE SAME COORDINATES'
     &     //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
      STOP
      END


C******************************************************************************
C                                                                             *
C    Transform from lon,lat (lamda,phi) coordinates into CPP coordinates.     *
C    Lon,Lat must be in radians.                                              *
C                                                                             *
C******************************************************************************

      SUBROUTINE CPP(X,Y,RLAMBDA,PHI,RLAMBDA0,PHI0)
      IMPLICIT NONE
      REAL*8 X,Y,RLAMBDA,PHI,RLAMBDA0,PHI0,R
      R=6378206.4d0
      X=R*(RLAMBDA-RLAMBDA0)*COS(PHI0)
      Y=PHI*R
      RETURN
      END


C******************************************************************************
C                                                                             *
C    Transform from CPP coordinates to lon,lat (lamda,phi) coordinates        *
C    Lon,Lat is in radians.                                                   *
C                                                                             *
C******************************************************************************

      SUBROUTINE INVCP(XXCP,YYCP,RLAMBDA,PHI,RLAMBDA0,PHI0)
      IMPLICIT NONE
      REAL*8 XXCP,YYCP,RLAMBDA,PHI,RLAMBDA0,PHI0,R
      R=6378206.4d0
      RLAMBDA=RLAMBDA0+XXCP/(R*COS(PHI0))
      PHI=YYCP/R
      RETURN
      END
