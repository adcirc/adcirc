C******************************************************************************
C PADCIRC VERSION 45.12 03/17/2006                                            *
C    last changes in this file prior to VERSION 44.15                         *
C                                                                             *
c                                                                             *
c   added local LOGICAL CHARMV declaration 04/21/2004                         *
c******************************************************************************
c                                                                      *
c   PADCIRC MODULE  ( HARM )                                           *
c                                                                      *
c   HA_SUBS.FOR     V3.01        11/9/95                               *
c                                                                      *
c   Least Square harmonic analysis of timeseries from ADCIRC2DDI_v27   *
c                                                                      *
c    Notes:                                                            *
c    1.)  Both the left hand side matrix and the right hand side       *
c         forcing vectors are continuously updated in time.  This      *
c         eliminates the need to store time series outputs for later   *
c         harmonic analysis.                                           *
c    2.)  The left hand side matrix and the right hand side forcing    *
c         vectors are output in the hotstart file and can be used to   *
c         perform harmonic analysis on an incomplete run.              *
c    3.)  Frequencies should be in rad/sec,times should be in sec.     *
c                                                                      *
c***********************************************************************
c                                                                      *
c    Program Written by:                                               *
c          R.A. Luettich, IMS UNC                                      *
c          J.J. Westerink, CE ND                                       *
c                                                                      *
c    Program Development History:                                      *
c    1.) lsq_stations_v004 by JJW                                      *
c    2.) LSQEX by RL used in 2D binary extr program                    *
c    3.) LSQRL by RL used in 1D test codes                             *
c    4.) LSQ2D v1.00-v2.26 by RL real time Harmonic Analysis for ADCIRC*
c    5.) HA_SUBS v3.01 by RL real time HA for ADCIRC separate          *
c        subroutines for elevation station, velocity station,          *
c        global elevation and global velocity harmonic analysis        *
c                                                                      *
c***********************************************************************
c                                                                      *
c SUBROUTINE LSQUPDLHS updates the LHS matrix                          *
c SUBROUTINE LSQUPDES updates the RHS load vector for elev stations    *
c SUBROUTINE LSQUPDVS updates the RHS load vector for velocity stations*
c SUBROUTINE LSQUPDEG updates the RHS load vector for elevation global *
c SUBROUTINE LSQUPDVG updates the RHS load vector for velocity global  *
c SUBROUTINE FULSOL fills out, decomposes and solves the matricies     *
c SUBROUTINE LSQSOLES solves & writes output for elevation stations    *
c SUBROUTINE LSQSOLVS solves & writes output for velocity stations     *
c SUBROUTINE LSQSOLEG solves & writes output for elevation global      *
c SUBROUTINE LSQSOLVG solves & writes output for velocity global       *
c SUBROUTINE HAHOUT writes HA parameters & LHS matrix to hotstart file *
c SUBROUTINE HAHOUTES writes elev sta RHS load vector to hotstart file *
c SUBROUTINE HAHOUTVS writes vel sta RHS load vector to hotstart file  *
c SUBROUTINE HAHOUTEG writes glob elev RHS load vector to hotstart file*
c SUBROUTINE HAHOUTVG writes glob vel RHS load vector to hotstart file *
c SUBROUTINE HACOLDS initializes HA param & LHS matrix for cold start  *
c SUBROUTINE HACOLDSES initializes elev sta RHS load vec for cold start*
c SUBROUTINE HACOLDSVS initializes vel sta RHS load vec for cold start *
c SUBROUTINE HACOLDSEG initializes glob ele RHS load vec for cold start*
c SUBROUTINE HACOLDSVG initializes glob vel RHS load vec for cold start*
c SUBROUTINE HAHOTS initializes HA params & LHS matrix for a hot start *
c SUBROUTINE HAHOTSES initializes elev sta RHS load vec for hot start  *
c SUBROUTINE HAHOTSVS initializes vel sta RHS load vec for hot start   *
c SUBROUTINE HAHOTSEG initializes glob elev RHS load vec for hot start *
c SUBROUTINE HAHOTSVG initializes glob vel RHS load vec for hot start  *
c                                                                      *
c***********************************************************************
c                                                                      *
c    INPUT FILES:                                                      *
c      - Frequency information is read in by ADCIRC from unit 15.      *
c        This information is passed in common block LSQFREQS.          *
c                                                                      *
c      - If the model is hot start, input is read from UNIT 67 or 68   *
c                                                                      *
c    OUTPUT FILES:                                                     *
C      UNIT 51 : HARMONIC CONSTITUENT ELEVATION VALUES AT SPECIFIED    *
C                  ELEVATION RECORDING STATION COORDINATES (ASCII)     *
C      UNIT 52 : HARMONIC CONSTITUENT VELOCITY VALUES AT SPECIFIED     *
C                  VELOCITY RECORDING STATION COORDINATES  (ASCII)     *
C      UNIT 53 : HARMONIC CONSTITUENT ELEVATIONS AT ALL NODES (ASCII)  *
C      UNIT 54 : HARMONIC CONSTITUENT VELOCITIES AT ALL NODES (ASCII)  *
C      UNIT 55 : COMPARISON BETWEEN THE MEAN AND VARIANCE OF THE TIME  *
C                  SERIES GENERATED BY THE MODEL AND THE MEAN AND      *
C                  VARIANCE OF A TIME SERIES RESYNTHESIZED FROM THE    *
C                  COMPUTED HARMONIC CONSTITUENTS.  THIS GIVES AN      *
C                  INDICATION OF HOW COMPLETE THE HARMONIC ANALYSIS    *
C                  WAS. (ASCII)                                        *
C      UNIT 67 or 68 : HOT START FILES (BINARY)                        *
c                                                                      *
c***********************************************************************
C
      MODULE HARM
C
      USE SIZES
C
      REAL(8),PRIVATE,PARAMETER :: PI=3.141592653589793D0
      INTEGER NFREQ
      CHARACTER*10,ALLOCATABLE ::  NAMEFR(:)
      REAL(SZ),    ALLOCATABLE ::  HAFREQ(:),HAFF(:),HAFACE(:)
C
      INTEGER, SAVE :: NZ, NF, MM, ITUD, ICALL
      REAL(8), SAVE :: TIMEUD
      REAL(SZ),ALLOCATABLE ::  HA(:,:)
      REAL(SZ),ALLOCATABLE ::  HAP(:),HAX(:)
      REAL(SZ),ALLOCATABLE,TARGET :: GLOELV(:,:)
      REAL(SZ),ALLOCATABLE,TARGET ::  GLOULV(:,:),GLOVLV(:,:)
      REAL(SZ),ALLOCATABLE,TARGET ::  STAELV(:,:)
      REAL(SZ),ALLOCATABLE,TARGET ::  STAULV(:,:),STAVLV(:,:)


C-----------------END OF DECLARATIONS---------------------------------------

      CONTAINS


C
C***********************************************************************
C  Allocate arays used by LSQ_HARM.
C
C  vjp 8/99
C***********************************************************************
C
      SUBROUTINE ALLOC_HA()
      IMPLICIT NONE
      ALLOCATE ( HAFREQ(MNHARF),HAFF(MNHARF),HAFACE(MNHARF) )
      ALLOCATE ( NAMEFR(MNHARF) )
C     
      ALLOCATE ( HA(2*MNHARF,2*MNHARF) )
      ALLOCATE ( HAP(2*MNHARF),HAX(2*MNHARF) )
      ALLOCATE ( GLOELV(2*MNHARF,MNP) )
      ALLOCATE ( GLOULV(2*MNHARF,MNP),GLOVLV(2*MNHARF,MNP) )
      ALLOCATE ( STAELV(2*MNHARF,MNSTAE) )
      ALLOCATE ( STAULV(2*MNHARF,MNSTAV),STAVLV(2*MNHARF,MNSTAV) )
      
      RETURN
      END SUBROUTINE
      

c***********************************************************************
c   Subroutine to update the Left Hand Side Matrix                     *
c                                                                      *
c  TIME  - ABSOLUTE MODEL TIME (SEC)                                   *
c  IT    - MODEL TIME STEP                                             *
c  icall - number of times the subroutine has been called              *
c  a     - Left Hand Side Matrix                                       *
c                                                                      *
c                        RL 11/7/95                                    *
c***********************************************************************
c
      SUBROUTINE LSQUPDLHS(TIME,IT)
      IMPLICIT NONE
      INTEGER IT,I,J,I1,I2,J1,J2
      REAL(SZ) TF1,TF2
      REAL(8) TIME
c
      icall = icall + 1
c     
c***** Update the Left Hand Side Matrix
c     Note: this is a symmetric matrix and therefore only store the
c     upper triangular part.  The lower part will be filled out in
c     SUBROUTINE FULSOL prior to the matrix's decomposition

c     Take care of the steady constituent if included in the analysis

      if(nf.eq.1) then
         ha(1,1)=icall
         do j=1,nfreq
            tf1=hafreq(j+nf)*time
            ha(1,2*j)   = ha(1,2*j) + cos(tf1)
            ha(1,2*j+1) = ha(1,2*j+1) + sin(tf1)
         end do
      endif

c   Take care of the other constituents

      do i=1,nfreq
         do j=i,nfreq
            i1=2*i-(1-nf)
            i2=i1+1
            j1=2*j-(1-nf)
            j2=j1+1
            tf1=hafreq(i+nf)*time
            tf2=hafreq(j+nf)*time
            ha(i1,j1) = ha(i1,j1) + cos(tf1)*cos(tf2)
            ha(i1,j2) = ha(i1,j2) + cos(tf1)*sin(tf2)
            ha(i2,j2) = ha(i2,j2) + sin(tf1)*sin(tf2)
            if(i2.le.j1) ha(i2,j1) = ha(i2,j1) + sin(tf1)*cos(tf2)
         end do
      end do

c   Record update time and time step

      TIMEUD = TIME
      ITUD = IT
      
      return
      end subroutine

c***********************************************************************
c   Subroutine to update the Right Hand Side Load Vectors for the      *
c   elevation station harmonic analysis.                               *
c                                                                      *
c  STAE  - STATION ELEVATION VALUES USED TO UPDATE LOAD VECTORS        *
c  NSTAE - NUMBER OF TIDAL ELEVATION RECORDING STATIONS                *
c                                                                      *
c  STAELV - station elevation load vector                              *
c                                                                      *
c                        RL 11/8/95                                    *
c***********************************************************************
c
      SUBROUTINE LSQUPDES(STAE,NSTAE)
      IMPLICIT NONE
      INTEGER NSTAE,N,I,I1,I2
      REAL(SZ) TF1,CTF1,STF1
      REAL(SZ) STAE(MNSTAE)
c
c***** Update the Right Hand Side Load Vectors
c
c   Take care of the steady constituent if included in the analysis

      if(nz.eq.0) then
         do n=1,NSTAE
            STAELV(1,N) = STAELV(1,N) + STAE(N)
         end do
      endif

c   Take care of the other constituents

      do i=1,nfreq
         i1=2*i-nz
         i2=i1+1
         tf1=hafreq(i+nf)*TIMEUD
         ctf1 = cos(tf1)
         stf1 = sin(tf1)
         do n=1,NSTAE
            STAELV(I1,N) = STAELV(I1,N) + STAE(N)*CTF1
            STAELV(I2,N) = STAELV(I2,N) + STAE(N)*STF1
         end do
      end do
C     
      return
      end subroutine
      
c***********************************************************************
c   Subroutine to update the Right Hand Side Load Vectors for the      *
c   velocity station harmonic analysis.                                *
c                                                                      *
c  STAU  - STATION U VELOCITY VALUES USED TO UPDATE LOAD VECTORS       *
c  STAV  - STATION V VELOCITY VALUES USED TO UPDATE LOAD VECTORS       *
c  NSTAV - NUMBER OF TIDAL CURRENT RECORDING STATIONS                  *
c                                                                      *
c  STAULV - station u velocity load vector                             *
c  STAVLV - station v velocity load vector                             *
c                                                                      *
c                        RL 11/8/95                                    *
c***********************************************************************
c
      SUBROUTINE LSQUPDVS(STAU,STAV,NSTAV)
      IMPLICIT NONE
      INTEGER NSTAV,N,I,I1,I2
      REAL(SZ) TF1,CTF1,STF1
      REAL(SZ) STAU(MNSTAV),STAV(MNSTAV)
c
c***** Update the Right Hand Side Load Vectors
c
c     Take care of the steady constituent if included in the analysis

      if(nz.eq.0) then
         do n=1,NSTAV
            STAULV(1,N) = STAULV(1,N) + STAU(N)
            STAVLV(1,N) = STAVLV(1,N) + STAV(N)
         end do
      endif

c     Take care of the other constituents

      do i=1,nfreq
         i1=2*i-nz
         i2=i1+1
         tf1=hafreq(i+nf)*TIMEUD
         ctf1 = cos(tf1)
         stf1 = sin(tf1)
         do n=1,NSTAV
            STAULV(I1,N) = STAULV(I1,N) + STAU(N)*CTF1
            STAVLV(I1,N) = STAVLV(I1,N) + STAV(N)*CTF1
            STAULV(I2,N) = STAULV(I2,N) + STAU(N)*STF1
            STAVLV(I2,N) = STAVLV(I2,N) + STAV(N)*STF1
         end do
      end do
C     
      return
      end subroutine


c***********************************************************************
c   Subroutine to update the Right Hand Side Load Vectors for the      *
c   global elevation harmonic analysis.                                *
c                                                                      *
c  GLOE  - GLOBAL ELEVATION VALUES USED TO UPDATE LOAD VECTORS         *
c  NP    - NUMBER OF POINTS IN GLOBAL GRID                             *
c                                                                      *
c  GLOELV - global elevation load vector                               *
c                                                                      *
c                        RL 11/8/95                                    *
c***********************************************************************
c
      SUBROUTINE LSQUPDEG(GLOE,NP)
      IMPLICIT NONE
      INTEGER I,J,NP,N,I1,I2,IR,IRE,K,JR
      REAL(SZ) TF1,CTF1,STF1
      REAL(SZ) GLOE(MNP)
c     
c*****Update the Right Hand Side Load Vectors
c     
c     Take care of the steady constituent if included in the analysis

      if(nz.eq.0) then
         do n=1,np
            GLOELV(1,N)=GLOELV(1,N)+GLOE(N)
         end do
      endif

c     Take care of the other constituents

      do i=1,nfreq
         i1=2*i-nz
         i2=i1+1
         tf1=hafreq(i+nf)*TIMEUD
         ctf1 = cos(tf1)
         stf1 = sin(tf1)
         do n=1,np
            GLOELV(I1,N)=GLOELV(I1,N)+GLOE(N)*CTF1
            GLOELV(I2,N)=GLOELV(I2,N)+GLOE(N)*STF1
         end do
      end do
C     
      return
      end subroutine


c***********************************************************************
c   Subroutine to update the Right Hand Side Load Vectors for the      *
c   global velocity harmonic analysis.                                 *
c                                                                      *
c  GLOU  - GLOBAL U VELOCITY VALUES USED TO UPDATE LOAD VECTORS        *
c  GLOV  - GLOBAL V VELOCITY VALUES USED TO UPDATE LOAD VECTORS        *
c  NP    - NUMBER OF POINTS IN GLOBAL GRID                             *
c                                                                      *
c  GLOULV - global u velocity load vector                              *
c  GLOVLV - global v velocity load vector                              *
c                                                                      *
c                        RL 11/8/95                                    *
c***********************************************************************
c
      SUBROUTINE LSQUPDVG(GLOU,GLOV,NP)
      IMPLICIT NONE
      INTEGER NP,NPI1,I1,I2,N,I,J,IR,IRE,K,JR
      REAL(SZ) TF1,CTF1,STF1
      REAL(SZ) GLOU(MNP),GLOV(MNP)
c     
c*****Update the Right Hand Side Load Vectors
c     
c     Take care of the steady constituent if included in the analysis

      if(nz.eq.0) then
         do n=1,np
            GLOULV(1,N) = GLOULV(1,N) + GLOU(N)
            GLOVLV(1,N) = GLOVLV(1,N) + GLOV(N)
         end do
      endif

c     Take care of the other constituents

      do i=1,nfreq
         i1=2*i-nz
         i2=i1+1
         tf1=hafreq(i+nf)*TIMEUD
         ctf1 = cos(tf1)
         stf1 = sin(tf1)
         do n=1,np
            GLOULV(I1,N) = GLOULV(I1,N) + GLOU(N)*CTF1
            GLOVLV(I1,N) = GLOVLV(I1,N) + GLOV(N)*CTF1
            GLOULV(I2,N) = GLOULV(I2,N) + GLOU(N)*STF1
            GLOVLV(I2,N) = GLOVLV(I2,N) + GLOV(N)*STF1
         end do
      end do
C     
      return
      end subroutine

c***********************************************************************
c   Subroutine to fill out, decompose and solve the lsq system         *
c   Solves system a*x=b by l*d*l(tr) decomp in full storage mode       *
c                                                                      *
c   NOTE: This routine has been modified so that the filling out and   *
c         decomposition (and only those operations) are done if        *
c         idecom=0.                                                    *
c                                                                      *
c   mm  -  actual dimension of a matrix                                *
c                                                                      *
c                        rl 11/7/95                                    *
c***********************************************************************
c
      subroutine fulsol(idecom)
      implicit none
      integer idecom,i,j,ir,ire,k,jr
      real(sz),allocatable ::  c(:),y(:)

c     
c**** If only want to fill out matrix and decomponse
c     
      if(idecom.eq.0) then
         
c     Set up the lower triangular part of the LHS a matrix
         
         do j=1,mm
            do i=j,mm
               ha(i,j)=ha(j,i)
            end do
         end do
         
c     Decomposition of matrix a

         do 100 ir=1,mm
            ire=ir+1
            do 20 j=ire,mm
 20         ha(ir,j)=ha(ir,j)/ha(ir,ir)
            if(ire.gt.mm) goto 100
            do 40 j=ire,mm
              do 40 k=ire,mm
 40           ha(k,j)=ha(k,j)-ha(k,ir)*ha(ir,j)
            do 50 j=ire,mm
 50           ha(j,ir)=0.0
 100     continue
         return
      endif

c...  solve for y by forward substitution for l*y=p

      allocate ( c(2*MNHARF),y(2*MNHARF) )
c     
      do 120 ir=1,mm
         y(ir)=hap(ir)
         do 110 jr=1,ir-1
 110        y(ir)=y(ir)-ha(jr,ir)*y(jr)
 120     continue

c...  calculate c=d**(-1)*y

         do 130 ir=1,mm
 130        c(ir)=y(ir)/ha(ir,ir)

c...  solve for x by back-substituting for l(tr)*x=c

            ir=mm
 140        continue
            hax(ir)=c(ir)
            do 150 jr=ir+1,mm
 150          hax(ir)=hax(ir)-ha(ir,jr)*hax(jr)
            ir=ir-1
            if(ir.ge.1) goto 140
            return
            end subroutine


c***********************************************************************
c   Subroutine to solve the system and write output for elevation      *
c   stations.                                                          *
c                                                                      *
c   nf=0  if no steady constituent                                     *
c   nf=1 if steady constituent                                         *
c                                                                      *
c  closed unit 51            08/23/05                                  *
c***********************************************************************
c
      SUBROUTINE LSQSOLES(NSTAE,DIRNAME)
      IMPLICIT NONE
      INTEGER NSTAE,N,I,J,K,I1,I2
      REAL(8) CONVRD
      REAL(SZ) EMAG,PHASEE,PHASEDE
      CHARACTER(LEN=*) DIRNAME
C
      convrd=180.d0/pi
c
c**** Open elevation station harmonic output file and write header information
c
      open(51,file=TRIM(LOCALDIR)//'/'//'fort.51')
      write(51,*) nfreq+nf
      do j=1,nfreq+nf
         write(51,3679) hafreq(j),HAFF(j),HAFACE(j),namefr(j)
      end do
 3679 format(1x,e20.10,1x,f10.7,1x,f12.8,1x,a10)
      write(51,*) NSTAE

c
c**** AT each STATION TRANSFER each load vector to p and solve the system
c
      DO N=1,NSTAE
         do k=1,mm
            hap(k)=STAELV(k,n)
         end do
         call fulsol(n)

c        Compute amplitude and phase for each frequency making sure that the
c        phase is between 0 and 360 deg.  Then write output.

         write(51,*) N
         do i=1,nfreq+nf
            if((nf.eq.1).and.(i.eq.1)) then
               emag=hax(i)/haff(i)
               phasee=0.
            else
               i1=2*i-1-nf
               i2=i1+1
               emag=sqrt(hax(i1)*hax(i1)+hax(i2)*hax(i2))/haff(i)
               if((hax(i1).eq.0.).and.(hax(i2).eq.0.)) then
                  phasee=0.
               else
                  phasee = atan2(hax(i2),hax(i1))
               endif
            endif
            phasede=convrd*phasee+haface(i)
            if(phasede.lt.0.) phasede=phasede+360.d0
            if(phasede.ge.360.d0) phasede=phasede-360.d0
               
            write(51,6635) emag,phasede
 6635       format(2x,e16.8,1x,f11.4)
         end do
            
      end do
      
	CLOSE(51)
C     
         return
         end subroutine
      

c***********************************************************************
c   Subroutine to solve the system and write output for velocity       *
c   stations.                                                          *
c                                                                      *
c   nf=0  if no steady constituent                                     *
c   nf=1  if steady constituent                                        *
c                                                                      *
c                        R.L. 11/8/95                                  *
c   closed unit 52            08/23/05                                 *
c***********************************************************************
c
      SUBROUTINE LSQSOLVS(NSTAV,DIRNAME)
      IMPLICIT NONE
      INTEGER NSTAV,I,J,N,K,I1,I2
      REAL(8) CONVRD
      REAL(SZ) UMAG,VMAG,PHASEU,PHASEV,PHASEDU,PHASEDV
      REAL(SZ),ALLOCATABLE :: Y(:)
      CHARACTER(LEN=*) DIRNAME
c
      convrd=180.d0/pi
c
      ALLOCATE ( Y(2*MNHARF) )
c
c**** Open velocity station harmonic output file and write header information
c
      open(52,file=TRIM(DIRNAME)//'/'//'fort.52')
      write(52,*) nfreq+nf
      do j=1,nfreq+nf
         write(52,3679) hafreq(j),HAFF(j),HAFACE(j),namefr(j)
      end do
 3679 format(1x,e20.10,1x,f10.7,1x,f12.8,1x,a10)
      write(52,*) NSTAV
c
c**** AT each STATION, transfer each load vector to p, solve system, 
c**** then write results
c     
      DO N=1,NSTAV
         do k=1,mm
            hap(k) = STAVLV(k,n)
         end do
         call fulsol(n)
         do k=1,mm
            y(k)=hax(k)
         end do
         do k=1,mm
            hap(k) = STAULV(k,n)
         end do
         call fulsol(n)
         
c        Compute amplitude and phase for each frequency making sure that the
c        phase is between 0 and 360 deg.  Then write output.

         write(52,*) N
         do i=1,nfreq+nf
            if((nf.eq.1).and.(i.eq.1)) then
               umag=hax(i)/haff(i)
               vmag=y(i)/haff(i)
               phaseu=0.
               phasev=0.
            else
               i1=2*i-1-nf
               i2=i1+1
               umag=sqrt(hax(i1)*hax(i1)+hax(i2)*hax(i2))/haff(i)
               vmag=sqrt(y(i1)*y(i1)+y(i2)*y(i2))/haff(i)
               if((hax(i1).eq.0.).and.(hax(i2).eq.0.)) then
                  phaseu=0.
               else
                  phaseu = atan2(hax(i2),hax(i1))
               endif
               if((y(i1).eq.0.).and.(y(i2).eq.0.)) then
                  phasev=0.
               else
                  phasev = atan2(y(i2),y(i1))
               endif
            endif
            phasedu=convrd*phaseu+haface(i)
            if(phasedu.lt.0.) phasedu=phasedu+360.d0
            if(phasedu.ge.360.d0) phasedu=phasedu-360.d0
            phasedv=convrd*phasev+haface(i)
            if(phasedv.lt.0.) phasedv=phasedv+360.d0
            if(phasedv.ge.360.d0) phasedv=phasedv-360.d0
            
            write(52,6636) umag,phasedu,vmag,phasedv
 6636       format(2x,e16.8,1x,f11.4,2x,e16.8,1x,f11.4)
         end do
         
      end do

      CLOSE(52)
      
      return
      end subroutine

c***********************************************************************
c   Subroutine to solve the system and write output for elevation      *
c   globally.                                                          *
c                                                                      *
c   nf=0  if no steady constituent                                     *
c   nf=1  if steady constituent                                        *
c                                                                      *
c                        R.L. 11/8/95                                  *
c                                                                      *
c   added local LOGICAL CHARMV declaration 04/09/2004                  *
c   closed unit 53 08/23/05                                            *
c***********************************************************************
c
      SUBROUTINE LSQSOLEG(NP,DIRNAME,ELAV,ELVA,CHARMV)
      IMPLICIT NONE
      LOGICAL CHARMV            !rl 04/09/2004
      integer J,NP,N,K,I,I1,I2,IT,IFR,NEAVMAX,NEAVMIN,
     &  NEVAMAX,NEVAMIN
      REAL(8)  CONVRD 
      REAL(SZ) EAVMAX,EVAMAX,EAVMIN,EVAMIN,EMAGT
      REAL(SZ) PHASEDE,EAV,ESQ,TIME,RSE,FTIME,EAVDIF,EVADIF
      REAL(SZ) ELAV(MNP),ELVA(MNP)
      REAL(SZ),ALLOCATABLE  ::  PHASEE(:),EMAG(:)
C
      INTEGER NTSTEPS,ITMV
      REAL(8) TIMEBEG
      REAL(SZ) DT,FMV
      COMMON /MEANSQ/ TIMEBEG,DT,FMV,NTSTEPS,ITMV
      CHARACTER(LEN=*) DIRNAME
c
      convrd=180.d0/pi
c
      ALLOCATE ( PHASEE(MNHARF),EMAG(MNHARF) )
c
c**** Open velocity station harmonic output file and write header information
c
      open(53,file=TRIM(DIRNAME)//'/'//'fort.53')
      write(53,*) nfreq+nf
      do j=1,nfreq+nf
         write(53,3679) hafreq(j),HAFF(j),HAFACE(j),namefr(j)
      end do
 3679 format(1x,e20.10,1x,f10.7,1x,f12.8,1x,a10)
      write(53,*) NP
C
      if (CHARMV) then
         EAVMAX=-999.
         EVAMAX=-999.
         EAVMIN= 999.
         EVAMIN= 999.
      end if
c
c***** AT each node transfer each load vector to p, solve and write output
c
      DO N=1,NP
         do k=1,mm
            hap(k) = GLOELV(k,n)
         end do
         call fulsol(n)
c
c        Compute amplitude and phase for each frequency making sure that the
c        phase is between 0 and 360 deg.  Then write output.
c
         write(53,*) N
         do i=1,nfreq+nf
            if((nf.eq.1).and.(i.eq.1)) then
               emag(i)=hax(i)
               emagt=emag(i)/haff(i)
               phasee(i)=0.
            else
               i1=2*i-1-nf
               i2=i1+1
               emag(i)=sqrt(hax(i1)*hax(i1)+hax(i2)*hax(i2))
               emagt=emag(i)/haff(i)
               if((hax(i1).eq.0.).and.(hax(i2).eq.0.)) then
                  phasee(i)=0.
               else
                  phasee(i) = atan2(hax(i2),hax(i1))
               endif
            endif
            phasede=convrd*phasee(i)+haface(i)
            if(phasede.lt.0.) phasede=phasede+360.d0
            if(phasede.ge.360.d0) phasede=phasede-360.d0
            write(53,6635) emagt,phasede
 6635       format(2x,e16.8,1x,f11.4)
         end do
         
         if (CHARMV) then
            eav = 0.
            esq = 0.
            do it=1,ntsteps
               TIME=TIMEBEG+DT*IT
               rse=0.
               do ifr=1,nfreq+nf
                  ftime=hafreq(ifr)*time
                  rse=rse+emag(ifr)*cos(ftime-phasee(ifr))
               end do
               eav=eav+rse
               esq=esq+rse*rse
            end do
            
         eav=eav/ntsteps
         esq=esq/ntsteps-eav*eav
         if(elav(n).eq.0.) then
            if(eav.eq.0.) eavdif=1.0d0
            if(eav.ne.0.) eavdif=99d19
         else
            eavdif=eav/elav(n)
         endif
         if(elva(n).eq.0.) then
            if(esq.eq.0.) evadif=1.0d0
            if(esq.ne.0.) evadif=99e19
         else
            evadif=esq/elva(n)
         endif
         write(55,*) n
         write(55,7637) elav(n),eav,eavdif,elva(n),esq,evadif
 7637    format(2x,3(e16.8,1x),2x,3(e16.8,1x))
         
         IF(EAVDIF.GT.EAVMAX) THEN
            EAVMAX=EAVDIF
            NEAVMAX=n
         ENDIF
         IF(EAVDIF.LT.EAVMIN) THEN
            EAVMIN=EAVDIF
            NEAVMIN=n
         ENDIF
         IF(EVADIF.GT.EVAMAX) THEN
            EVAMAX=EVADIF
            NEVAMAX=n
         ENDIF
         IF(EVADIF.LT.EVAMIN) THEN
            EVAMIN=EVADIF
            NEVAMIN=n
         ENDIF
      endif                     ! charmv
      
      end do

      CLOSE(53)
      
      if (charmv) then
c
      WRITE(16,7740)
 7740 FORMAT(///,5X,'THE LARGEST VALUES OF THE RATIO ',
     &              'RESYNTHESIZED ELEV TIME SERIES/RAW TIME SERIES:',/)
      WRITE(16,7741) EAVMAX,NEAVMAX
      WRITE(16,7742) EVAMAX,NEVAMAX
      WRITE(16,7747)
 7747 FORMAT(/,5X,'THE LOWEST VALUES OF THE RATIO ',
     &            'RESYNTHESIZED ELEV TIME SERIES/RAW TIME SERIES:',/)
      WRITE(16,7741) EAVMIN,NEAVMIN
      WRITE(16,7742) EVAMIN,NEVAMIN
 7741 FORMAT(9X,'  AVERAGE ELEVATION RATIO = ',E15.7,' AT NODE ',I8)
 7742 FORMAT(9X,' VARIANCE ELEVATION RATIO = ',E15.7,' AT NODE ',I8)
c     
      endif                     ! charmv
c
      return
      end subroutine

c***********************************************************************
c   Subroutine to solve the system and write output for velocity       *
c   globally.                                                          *
c                                                                      *
c   nf=0  if no steady constituent                                     *
c   nf=1  if steady constituent                                        *
c                                                                      *
c                        R.L. 11/10/95                                 *
c                                                                      *
c   added local LOGICAL CHARMV declaration 04/09/2004                  *
c   closed unit 54 08/23/05                                            *
c***********************************************************************
c
      SUBROUTINE LSQSOLVG(NP,DIRNAME,
     &                    XVELAV,YVELAV,XVELVA,YVELVA,CHARMV)
      IMPLICIT NONE
      LOGICAL CHARMV            !rl 04/09/2004      
      INTEGER NP,I,J,N,K,I1,I2,IT,IFR
      INTEGER NUAVMAX,NUAVMIN,NVAVMAX,NVAVMIN,NUVAMAX,NUVAMIN,
     &  NVVAMAX,NVVAMIN
      REAL(SZ) UAV,VAV,USQ,VSQ,TIME,FTIME,RSU,RSV
      REAL(SZ) UAVMAX,VAVMAX,UVAMAX,VVAMAX,UAVMIN,VAVMIN,
     & UVAMIN,VVAMIN,PHASEDU,PHASEDV,UAVDIF,VAVDIF,UMAGT,VMAGT,
     & UVADIF,VVADIF
      REAL(8) CONVRD
      REAL(SZ) XVELAV(MNP),YVELAV(MNP),XVELVA(MNP),YVELVA(MNP)
      REAL(SZ),ALLOCATABLE :: UMAG(:),VMAG(:),PHASEU(:),PHASEV(:)
      REAL(SZ),ALLOCATABLE :: Y(:)
      INTEGER NTSTEPS,ITMV
      REAL(8) TIMEBEG
      REAL(SZ) DT,FMV
      COMMON /MEANSQ/ TIMEBEG,DT,FMV,NTSTEPS,ITMV
      CHARACTER(LEN=*) DIRNAME
c
      convrd=180.d0/pi
c
      ALLOCATE ( Y(2*MNHARF) )
      ALLOCATE ( UMAG(MNHARF),VMAG(MNHARF) )
      ALLOCATE ( PHASEU(MNHARF),PHASEV(MNHARF) )
c
c**** Open velocity station harmonic output file and write header information
c
      open(54,file=TRIM(DIRNAME)//'/'//'fort.54')
      write(54,*) nfreq+nf
      do j=1,nfreq+nf
         write(54,3679) hafreq(j),HAFF(j),HAFACE(j),namefr(j)
      end do
 3679 format(1x,e20.10,1x,f10.7,1x,f12.8,1x,a10)
      write(54,*) NP
      
      if ( charmv ) then
         UAVMAX=-999.
         VAVMAX=-999.
         UVAMAX=-999.
         VVAMAX=-999.
         UAVMIN= 999.
         VAVMIN= 999.
         UVAMIN= 999.
         VVAMIN= 999.
      endif                     ! charmv
c
c***** AT each node transfer each load vector to p, solve and write output
c
      DO N=1,NP
         do k=1,mm
            hap(k) = GLOVLV(k,n)
         end do
         call fulsol(n)
         do k=1,mm
            y(k)=hax(k)
         end do
         do k=1,mm
            hap(k) = GLOULV(k,n)
         end do
         call fulsol(n)
         write(54,*) n
         do i=1,nfreq+nf
            if((nf.eq.1).and.(i.eq.1)) then
               umag(i)=hax(i)
               umagt=umag(i)/haff(i)
               vmag(i)=y(i)
               vmagt=vmag(i)/haff(i)
               phaseu(i)=0.
               phasev(i)=0.
            else
               i1=2*i-1-nf
               i2=i1+1
               umag(i)=sqrt(hax(i1)*hax(i1)+hax(i2)*hax(i2))
               umagt=umag(i)/haff(i)
               vmag(i)=sqrt(y(i1)*y(i1)+y(i2)*y(i2))
               vmagt=vmag(i)/haff(i)
               if((hax(i1).eq.0.).and.(hax(i2).eq.0.)) then
                  phaseu(i)=0.
               else
                  phaseu(i)=atan2(hax(i2),hax(i1))
               endif
               if((y(i1).eq.0.).and.(y(i2).eq.0.)) then
                  phasev(i)=0.
               else
                  phasev(i)=atan2(y(i2),y(i1))
               endif
            endif
            phasedu=convrd*phaseu(i)+haface(i)
            if(phasedu.lt.0.) phasedu=phasedu+360.d0
            if(phasedu.ge.360.d0) phasedu=phasedu-360.d0
            phasedv=convrd*phasev(i)+haface(i)
            if(phasedv.lt.0.) phasedv=phasedv+360.d0
            if(phasedv.ge.360.d0) phasedv=phasedv-360.d0

            write(54,6636) umagt,phasedu,vmagt,phasedv
 6636       format(2x,e16.8,1x,f11.4,2x,e16.8,1x,f11.4)
         end do

CHARMV...UNCOMMENT THE FOLLOWING LINES TO COMPUTE MEANS AND VARIANCES
CHARMV...FOR CHECKING THE HARMONIC ANALYSIS RESULTS.
CHARMV...Resynthesize the time series to compute the average and variances.
CHARMV...Compare resynthesized values with those computed during time stepping.
         if ( charmv ) then
            uav = 0.
            vav = 0.
            usq = 0.
            vsq = 0.
            do it=1,ntsteps
               TIME=TIMEBEG+DT*IT
               rsu=0.
               rsv=0.
               do ifr=1,nfreq+nf
                  ftime=hafreq(ifr)*time
                  rsu=rsu+umag(ifr)*cos(ftime-phaseu(ifr))
                  rsv=rsv+vmag(ifr)*cos(ftime-phasev(ifr))
               end do
               uav=uav+rsu
               vav=vav+rsv
               usq=usq+rsu*rsu
               vsq=vsq+rsv*rsv
            end do

            uav=uav/ntsteps
            vav=vav/ntsteps
            usq=usq/ntsteps-uav*uav
            vsq=vsq/ntsteps-vav*vav
            if(xvelav(n).eq.0.) then
               if(uav.eq.0.) uavdif=1.0d0
               if(uav.ne.0.) uavdif=99e19
            else
               uavdif=uav/xvelav(n)
            endif
            if(yvelav(n).eq.0.) then
               if(vav.eq.0.) vavdif=1.0d0
               if(vav.ne.0.) vavdif=99e19
            else
               vavdif=vav/yvelav(n)
            endif
            if(xvelva(n).eq.0.) then
               if(usq.eq.0.) uvadif=1.0d0
               if(usq.ne.0.) uvadif=99e19
            else
               uvadif=usq/xvelva(n)
            endif
            if(yvelva(n).eq.0.) then
               if(vsq.eq.0.) vvadif=1.0d0
               if(vsq.ne.0.) vvadif=99e19
            else
               vvadif=vsq/yvelva(n)
            endif
            write(55,*) n
            write(55,7637) xvelav(n),uav,uavdif,xvelva(n),usq,uvadif
            write(55,7637) yvelav(n),vav,vavdif,yvelva(n),vsq,vvadif
 7637       format(2x,3(e16.8,1x),2x,3(e16.8,1x))

            IF(UAVDIF.GT.UAVMAX) THEN
               UAVMAX=UAVDIF
               NUAVMAX=n
            ENDIF
            IF(UAVDIF.LT.UAVMIN) THEN
               UAVMIN=UAVDIF
               NUAVMIN=n
            ENDIF
            IF(VAVDIF.GT.VAVMAX) THEN
               VAVMAX=VAVDIF
               NVAVMAX=n
            ENDIF
            IF(VAVDIF.LT.VAVMIN) THEN
               VAVMIN=VAVDIF
               NVAVMIN=n
            ENDIF
            IF(UVADIF.GT.UVAMAX) THEN
               UVAMAX=UVADIF
               NUVAMAX=n
            ENDIF
            IF(UVADIF.LT.UVAMIN) THEN
               UVAMIN=UVADIF
               NUVAMIN=n
            ENDIF
            IF(VVADIF.GT.VVAMAX) THEN
               VVAMAX=VVADIF
               NVVAMAX=n
            ENDIF
            IF(VVADIF.LT.VVAMIN) THEN
               VVAMIN=VVADIF
               NVVAMIN=n
            ENDIF

         endif                  !  charmv

      end do

      CLOSE(54)

      if ( charmv ) then 
c
         WRITE(16,7740)
 7740    FORMAT(///,5X,'THE LARGEST VALUES OF THE RATIO ',
     &              'RESYNTHESIZED VEL TIME SERIES/RAW TIME SERIES:',/)
         WRITE(16,7743) UAVMAX,NUAVMAX
         WRITE(16,7744) UVAMAX,NUVAMAX
         WRITE(16,7745) VAVMAX,NVAVMAX
         WRITE(16,7746) VVAMAX,NVVAMAX
         WRITE(16,7747)
 7747    FORMAT(//,5X,'THE LOWEST VALUES OF THE RATIO ',
     &             'RESYNTHESIZED VEL TIME SERIES/RAW TIME SERIES:',/)
         WRITE(16,7743) UAVMIN,NUAVMIN
         WRITE(16,7744) UVAMIN,NUVAMIN
         WRITE(16,7745) VAVMIN,NVAVMIN
         WRITE(16,7746) VVAMIN,NVVAMIN
 7743    FORMAT(9X,' AVERAGE U VELOCITY RATIO = ',E15.7,' AT NODE ',I8)
 7744    FORMAT(9X,'VARIANCE U VELOCITY RATIO = ',E15.7,' AT NODE ',I8)
 7745    FORMAT(9X,' AVERAGE V VELOCITY RATIO = ',E15.7,' AT NODE ',I8)
 7746    FORMAT(9X,'VARIANCE V VELOCITY RATIO = ',E15.7,' AT NODE ',I8)
c     
      endif                     ! charmv
c
      return
      end subroutine

C
c***********************************************************************
c   Subroutine to initialize parameters for harmonic analysis with a   *
c   cold start.                                                        *
c                                                                      *
c                        R.L.  11/9/95                                 *
c***********************************************************************
c
      SUBROUTINE HACOLDS(HAFREQ)
      implicit none
      INTEGER I,J
      REAL(SZ) HAFREQ(MNHARF)
C
      if (hafreq(1).eq.0.0) then
         nz=0
         nf=1
      else
         nz=1
         nf=0
      endif
c
      nfreq=nfreq-nf
      mm=2*nfreq+nf
c
      do i=1,mm
         do j=1,mm
            ha(i,j)=0.
         end do
      end do
      icall=0
C     
      return
      end subroutine

c***********************************************************************
c   Subroutine to initialize global elevation load vectors for         *
c   harmonic analysis with a cold start.                               *
c                                                                      *
c                        R.L.  11/9/95                                 *
c***********************************************************************
c
      SUBROUTINE HACOLDSEG(NP)
      implicit none
      INTEGER NP,I,N
c
      do i=1,mm
         do N=1,NP
            GLOELV(I,N)=0.
         end do
      end do
C
      return
      end subroutine

c***********************************************************************
c   Subroutine to initialize elevation station load vectors for        *
c   harmonic analysis with a cold start.                               *
c                                                                      *
c                        R.L.  11/9/95                                 *
c***********************************************************************
c
      SUBROUTINE HACOLDSES(NSTAE)
      implicit none
      INTEGER  NSTAE
      INTEGER I,N
C
      do i=1,mm
         do N=1,NSTAE
            STAELV(I,N)=0.
         end do
      end do
C     
      return
      end subroutine

c***********************************************************************
c   Subroutine to initialize global velocity load vectors for          *
c   harmonic analysis with a cold start.                               *
c                                                                      *
c                        R.L.  11/9/95                                 *
c***********************************************************************
c
      SUBROUTINE HACOLDSVG(NP)
      implicit none
      INTEGER NP,I,N
C
      do i=1,mm
         do N=1,NP
            GLOULV(I,N)=0.
            GLOVLV(I,N)=0.
         end do
      end do
C
      return
      end subroutine

c***********************************************************************
c   Subroutine to initialize elevation station load vectors for        *
c   harmonic analysis with a cold start.                               *
c                                                                      *
c                        R.L.  11/9/95                                 *
c***********************************************************************
c
      SUBROUTINE HACOLDSVS(NSTAV)
      implicit none
      INTEGER NSTAV,I,N
c
      do i=1,mm
         do N=1,NSTAV
            STAULV(I,N)=0.
            STAVLV(I,N)=0.
         enddo
      enddo
C
      return
      end subroutine

c***********************************************************************
c   Subroutine to read in and initialize harmonic analysis for a hot   *
c   start.                                                             *
c                                                                      *
c   Checks are made to ensure agreement between values read in from    *
c   the hotstart file and values read in from the UNIT 15 file.        *
c                                                                      *
c                        R.L. 11/9/95                                  *
c***********************************************************************
c
      SUBROUTINE HAHOTS(NSTAE,NSTAV,NP,ISTAE,ISTAV,IGLOE,IGLOV,
     &   NSCREEN,ScreenUnit,IHOTSTP,IHOT,MYPROC)
      IMPLICIT NONE
      INTEGER NSTAE,NSTAV,NP, MYPROC
      INTEGER ISTAE, ISTAV, IGLOE, IGLOV, IFLAG, I, J
C     
      INTEGER INFREQ, INSTAE, INSTAV, INP, INZ, INF
      INTEGER IISTAE, IISTAV, IIGLOE, IIGLOV
      INTEGER IHOTSTP, IHOT, IMM, IICALL
      INTEGER NSCREEN, ScreenUnit 
      REAL(SZ) FDIFF
C     
      REAL(SZ),ALLOCATABLE ::  IFREQ(:),IFF(:),IFACE(:)
      CHARACTER*10,ALLOCATABLE :: INAMEFR(:)
C     
      CHARACTER*16 FNAME
      CHARACTER*8 FNAM8(2)
      EQUIVALENCE (FNAM8(1),FNAME)
c     
c***** Compute parameter values for checking
c
      if (hafreq(1).eq.0.0) then
         nz=0
         nf=1
      else
         nz=1
         nf=0
      endif
      nfreq=nfreq-nf
      mm=2*nfreq+nf

      ALLOCATE ( IFREQ(MNHARF),IFF(MNHARF),IFACE(MNHARF) )
      ALLOCATE ( INAMEFR(MNHARF) )
c
c***** Read in and check various parameter values
c
      READ(IHOT,REC=IHOTSTP+1) inz
      READ(IHOT,REC=IHOTSTP+2) inf
      READ(IHOT,REC=IHOTSTP+3) imm
      READ(IHOT,REC=IHOTSTP+4) inp
      READ(IHOT,REC=IHOTSTP+5) instae
      READ(IHOT,REC=IHOTSTP+6) instav
      READ(IHOT,REC=IHOTSTP+7) iistae
      READ(IHOT,REC=IHOTSTP+8) iistav
      READ(IHOT,REC=IHOTSTP+9) iigloe
      READ(IHOT,REC=IHOTSTP+10) iiglov
      READ(IHOT,REC=IHOTSTP+11) iicall
      READ(IHOT,REC=IHOTSTP+12) infreq
      IHOTSTP = IHOTSTP+12
c
      iflag=0
      if(nz.ne.inz) iflag=1
      if(nf.ne.inf) iflag=1
      if(mm.ne.imm) iflag=1
      if(np.ne.inp) iflag=1
      if(nstae.ne.instae) iflag=1
      if(nstav.ne.instav) iflag=1
      if(istae.ne.iistae) iflag=1
      if(istav.ne.iistav) iflag=1
      if(igloe.ne.iigloe) iflag=1
      if(iglov.ne.iiglov) iflag=1
      if(nfreq.ne.infreq) iflag=1
c
      do i=1,nfreq+nf
         READ(IHOT,REC=IHOTSTP+1) FNAM8(1)
         READ(IHOT,REC=IHOTSTP+2) FNAM8(2)
         IHOTSTP = IHOTSTP + 2
         INAMEFR(I) = FNAME
         read(IHOT,REC=IHOTSTP+1) ifreq(i)
         read(IHOT,REC=IHOTSTP+2) iff(i)
         read(IHOT,REC=IHOTSTP+3) iface(i)
         IHOTSTP = IHOTSTP + 3

         if(namefr(i).ne.inamefr(i)) iflag=1
         if(abs(hafreq(i)+ifreq(i)).lt.1.0d-30) then
            fdiff=0.
         else
            fdiff=abs(hafreq(i)-ifreq(i))/abs(hafreq(i)+ifreq(i))
         endif
         if(fdiff.ge.1.d-6) iflag=1
         if(abs(HAFF(i)+iFF(i)).lt.1d-30) then
            fdiff=0.
         else
            fdiff=abs(HAFF(i)-iFF(i))/abs(HAFF(i)+iFF(i))
         endif
         if(fdiff.ge.1.d-6) iflag=1
         if(abs(HAFACE(i)+iFACE(i)).lt.1d-30) then
            fdiff=0.
         else
            fdiff=abs(HAFACE(i)-iFACE(i))/abs(HAFACE(i)+iFACE(i))
         endif
         if(fdiff.ge.1.d-6) iflag=1
      end do
      if(iflag.eq.1) goto 999
c
c***** Read in time of most recent H.A. update
c
      READ(IHOT,REC=IHOTSTP+1) TIMEUD
      READ(IHOT,REC=IHOTSTP+2) ITUD
      IHOTSTP = IHOTSTP + 2
c
c***** Read in RHS Matrix
c
      do i=1,mm
         do j=1,mm
            READ(IHOT,REC=IHOTSTP) HA(I,J)
            IHOTSTP = IHOTSTP + 1
         end do
      end do

c
c***** FATAL Error Messages
c
 999  continue
      if(iflag.ne.0) then
         if(NSCREEN.NE.0.AND.MYPROC.EQ.0) write(ScreenUnit,1000)
         write(16,1000)
 1000    FORMAT(////,5x,'***** DISCREPANCY IN HARMONIC ANALYSIS HOT ',
     +        'START FILE *****',/)
      endif

      if(iflag.eq.1) then
         if(nz.ne.inz) then
            if(NSCREEN.NE.0.AND.MYPROC.EQ.0) then
               write(ScreenUnit,2010) inz,nz
            endif
            write(16,2010) inz,nz
 2010       format(5x,'NZ COMPUTED FROM UNIT 14 INPUT = ',I2,
     +           ', NZ READ IN FROM HOT START FILE = ',I2,/)
         endif
         if(nf.ne.inf) then
            if(NSCREEN.NE.0.AND.MYPROC.EQ.0) then
               write(ScreenUnit,2020) inf,nf
            endif
            write(16,2020) inf,nf
 2020       format(5x,'NF COMPUTED FROM UNIT 14 INPUT = ',I2,
     +           ', NF READ IN FROM HOT START FILE = ',I2,/)
         endif
         if(mm.ne.imm) then
            if(NSCREEN.NE.0.AND.MYPROC.EQ.0) then
               write(ScreenUnit,2030) imm,mm
            endif
            write(16,2030) imm,mm
 2030       format(5x,'MM COMPUTED FROM UNIT 14 INPUT = ',I2,
     +           ', MM READ IN FROM HOT START FILE = ',I2,/)
         endif
         if(np.ne.inp) then
            if(NSCREEN.NE.0.AND.MYPROC.EQ.0) then
               write(ScreenUnit,2040) inp,np
            endif
            write(16,2040) inp,np
 2040       format(5x,'NP READ IN FROM UNIT 15 = ',I2,
     +           ', NP READ IN FROM HOT START FILE = ',I2,/)
         endif
         if(nstae.ne.instae) then
            if(NSCREEN.NE.0.AND.MYPROC.EQ.0) then
               write(ScreenUnit,2050) instae,nstae
            endif 
            write(16,2050) instae,nstae
 2050       format(5x,'NSTAE READ IN FROM UNIT 15 = ',I2,
     +           ', NSTAE READ IN FROM HOT START FILE = ',I2,/)
         endif
         if(nstav.ne.instav) then
            if(NSCREEN.NE.0.AND.MYPROC.EQ.0) then
               write(ScreenUnit,2060) instav,nstav
            endif
            write(16,2060) instav,nstav
 2060       format(5x,'NSTAV READ IN FROM UNIT 15 = ',I2,
     +           ', NSTAV READ IN FROM HOT START FILE = ',I2,/)
         endif
         if(istae.ne.iistae) then
            if(NSCREEN.ne.0.AND.MYPROC.EQ.0) then
               write(ScreenUnit,2070) iistae,istae
            endif
            write(16,2070) iistae,istae
 2070       format(5x,'ISTAE READ IN FROM UNIT 15 = ',I2,
     +           ', ISTAE READ IN FROM HOT START FILE = ',I2,/)
         endif
         if(istav.ne.iistav) then
            if(NSCREEN.NE.0.AND.MYPROC.EQ.0) then
               write(ScreenUnit,2080) iistav,istav
            endif
            write(16,2080) iistav,istav
 2080       format(5x,'ISTAV READ IN FROM UNIT 15 = ',I2,
     +           ', ISTAV READ IN FROM HOT START FILE = ',I2,/)
         endif
         if(igloe.ne.iigloe) then
            if(NSCREEN.NE.0.AND.MYPROC.EQ.0) then
               write(ScreenUnit,2090) iigloe,igloe
            endif
            write(16,2090) iigloe,igloe
 2090       format(5x,'IGLOE READ IN FROM UNIT 15 = ',I2,
     +           ', IGLOE READ IN FROM HOT START FILE = ',I2,/)
         endif
         if(iglov.ne.iiglov) then
            if(NSCREEN.NE.0.AND.MYPROC.EQ.0) then
               write(ScreenUnit,2100) iiglov,iglov
            endif
            write(16,2100) iiglov,iglov
 2100       format(5x,'IGLOV READ IN FROM UNIT 15 = ',I2,
     +           ', IGLOV READ IN FROM HOT START FILE = ',I2,/)
         endif
         if(nfreq.ne.infreq) then
            if(NSCREEN.NE.0.AND.MYPROC.EQ.0) then
               write(ScreenUnit,2110) infreq,nfreq
            endif
            write(16,2110) infreq,nfreq
 2110       format(5x,'NFREQ COMPUTED FROM UNIT 15 INPUT = ',I2,
     +           ', NFREQ READ IN FROM HOT START FILE = ',I2,/)
         endif
         do i=1,nfreq+nf
            if(namefr(i).ne.inamefr(i)) then
               if(NSCREEN.NE.0.AND.MYPROC.EQ.0) 
     &              write(ScreenUnit,2120) i,inamefr(i),namefr(i)
               write(16,2120) i,namefr(i),namefr(i)
 2120          format(5x,'FOR CONSTITUENT # ',I3,
     +              ', NAMEFR READ IN FROM UNIT 15 = ',A10,
     +              ', NAMEFR READ IN FROM HOT START FILE = ',A10,/)
            endif
            if(hafreq(i).ne.ifreq(i)) then
               if(NSCREEN.NE.0.AND.MYPROC.EQ.0) 
     &              write(ScreenUnit,2130) i,ifreq(i),hafreq(i)
               write(16,2130) i,ifreq(i),hafreq(i)
 2130          format(5x,'FOR CONSTITUENT # ',I3,
     +              ', FREQ READ IN FROM UNIT 15 = ',D20.10,
     +              ', FREQ READ IN FROM HOT START FILE = ',D20.10,/)
            endif
            if(HAFF(i).ne.iFF(i)) then
               if(NSCREEN.NE.0.AND.MYPROC.EQ.0) 
     &              write(ScreenUnit,2140) i,iff(i),haff(i)
               write(16,2140) i,iff(i),haff(i)
 2140          format(5x,'FOR CONSTITUENT # ',I3,
     +              ', FF READ IN FROM UNIT 15 = ',F10.5,
     +              ', FF READ IN FROM HOT START FILE = ',F10.5,/)
            endif
            if(HAFACE(i).ne.iFACE(i)) then
               if(NSCREEN.NE.0.AND.MYPROC.EQ.0) 
     &              write(ScreenUnit,2150) i,iface(i),haface(i)
               write(16,2150) i,iface(i),haface(i)
 2150          format(5x,'FOR CONSTITUENT # ',I3,
     +              ', FACE READ IN FROM UNIT 15 = ',F10.5,
     +              ', FACE READ IN FROM HOT START FILE = ',F10.5,/)
            endif
         end do
         if(NSCREEN.NE.0.AND.MYPROC.EQ.0) write(ScreenUnit,1010)
         write(16,1010)
 1010    FORMAT(//,5x,'********** RUN TERMINATED **********',/)
         stop
      endif

      return
      end subroutine


c***********************************************************************
c   Subroutine to read in and initialize the global elevation load     *
c   vector for harmonic analysis with a hot start.                     *
c                                                                      *
c                        R.L. 11/9/95                                  *
c***********************************************************************
c
      SUBROUTINE HAHOTSEG(NP,IHOTSTP,IHOT)
      implicit none
      integer np, ihotstp, ihot, n,i
c
c***** Read in Global Elevation LHS load vector
c
      do n=1,np
         do i=1,mm
            READ(IHOT,REC=IHOTSTP) GLOELV(I,N)
            IHOTSTP=IHOTSTP+1
         end do
      end do
C
      return
      end subroutine


c***********************************************************************
c   Subroutine to read in and initialize the elevation station load    *
c   vector for harmonic analysis with a hot start.                     *
c                                                                      *
c                        R.L. 11/9/95                                  *
c***********************************************************************
c
      SUBROUTINE HAHOTSES(NSTAE,IHOTSTP,IHOT)
      implicit none
      INTEGER NSTAE,IHOTSTP,IHOT,I,N
c
c***** Read in Station Elevation LHS load vector
c
      do n=1,NSTAE
         do i=1,mm
            READ(IHOT,REC=IHOTSTP) STAELV(I,N)
            IHOTSTP=IHOTSTP+1
         end do
      end do
C
      return
      end subroutine


c***********************************************************************
c   Subroutine to read in and initialize the global velocity load      *
c   vector for harmonic analysis with a hot start.                     *
c                                                                      *
c                        R.L. 11/9/95                                  *
c***********************************************************************
c
      SUBROUTINE HAHOTSVG(NP,IHOTSTP,IHOT)
      implicit none
      integer np, ihotstp, ihot, n, i
c
c***** Read in Global Velocity LHS load vector
c
      do n=1,np
         do i=1,mm
            READ(IHOT,REC=IHOTSTP+1) GLOULV(I,N)
            READ(IHOT,REC=IHOTSTP+2) GLOVLV(I,N)
            IHOTSTP = IHOTSTP + 2
         end do
      end do
C
      return
      end subroutine


c***********************************************************************
c   Subroutine to read in and initialize the velocity station load     *
c   vector for harmonic analysis with a hot start.                     *
c                                                                      *
c                        R.L. 11/9/95                                  *
c***********************************************************************
c
      SUBROUTINE HAHOTSVS(NSTAV,IHOTSTP,IHOT)
      implicit none
      INTEGER NSTAV,IHOTSTP,IHOT,N,I
c
c***** Read in Station Velocity LHS load vector
c
      do n=1,NSTAV
         do i=1,mm
            READ(IHOT,REC=IHOTSTP+1) STAULV(I,N)
            READ(IHOT,REC=IHOTSTP+2) STAVLV(I,N)
            IHOTSTP = IHOTSTP + 2
         enddo
      enddo
C
      return
      end subroutine


c***********************************************************************
c   Subroutine to write out to the hotstart file (UNITS 67 and 68)     *
c   header information and the LHS matrix for the harmonic analysis    *
c                                                                      *
c                        R.L.  11/8/95                                 *
C     jgf45.07 name changes to variables ISTAE:IGLOV -> NHASE:NHAGV    *
c***********************************************************************
c
      SUBROUTINE HAHOUT(NP,NSTAE,NSTAV,NHASE,NHASV,NHAGE,NHAGV,
     &  IOUNIT,IHOTSTP)
      implicit none
      INTEGER NP,NSTAE,NSTAV,NHASE,AE,NHASV
      INTEGER NHAGE,NHAGV,IOUNIT,IHOTSTP,I,J
      CHARACTER*16 FNAME
      CHARACTER*8 FNAM8(2)
      EQUIVALENCE (FNAM8(1),FNAME)

c
c***** Write Out various parameter values
c
      WRITE(IOUNIT,REC=IHOTSTP+1) NZ
      WRITE(IOUNIT,REC=IHOTSTP+2) NF
      WRITE(IOUNIT,REC=IHOTSTP+3) MM
      WRITE(IOUNIT,REC=IHOTSTP+4) NP
      WRITE(IOUNIT,REC=IHOTSTP+5) NSTAE
      WRITE(IOUNIT,REC=IHOTSTP+6) NSTAV
      WRITE(IOUNIT,REC=IHOTSTP+7) NHASE
      WRITE(IOUNIT,REC=IHOTSTP+8) NHASV
      WRITE(IOUNIT,REC=IHOTSTP+9) NHAGE
      WRITE(IOUNIT,REC=IHOTSTP+10) NHAGV
      WRITE(IOUNIT,REC=IHOTSTP+11) ICALL
      WRITE(IOUNIT,REC=IHOTSTP+12) NFREQ
      IHOTSTP = IHOTSTP+12

      do i=1,nfreq+nf
         FNAME=NAMEFR(I)
         WRITE(IOUNIT,REC=IHOTSTP+1) FNAM8(1)
         WRITE(IOUNIT,REC=IHOTSTP+2) FNAM8(2)
         IHOTSTP=IHOTSTP+2
         WRITE(IOUNIT,REC=IHOTSTP+1) hafreq(i)
         WRITE(IOUNIT,REC=IHOTSTP+2) HAFF(i)
         WRITE(IOUNIT,REC=IHOTSTP+3) HAFACE(i)
         IHOTSTP=IHOTSTP+3
      end do

c
c***** Write Out time of most recent H.A. update
c
      WRITE(IOUNIT,REC=IHOTSTP+1) TIMEUD
      WRITE(IOUNIT,REC=IHOTSTP+2) ITUD
      IHOTSTP=IHOTSTP+2
c
c***** Write Out LHS Matrix
c
      do i=1,mm
         do j=1,mm
            WRITE(IOUNIT,REC=IHOTSTP) HA(I,J)
            IHOTSTP = IHOTSTP + 1
         END DO
      END DO

      return
      end subroutine

c***********************************************************************
c   Subroutine to write global elevation harmonic analysis RHS load    *
c   vector to a hot start file (UNITS 67 and 68)                       *
c                                                                      *
c                        R.L.  11/8/95                                 *
c***********************************************************************
c
      SUBROUTINE HAHOUTEG(NP,IOUNIT,IHOTSTP)
      implicit none
      INTEGER IOUNIT
      INTEGER NP,IHOTSTP,N,I 
c
c***** Write Out Global Elevation RHS load vector
c
      do n=1,np
         do i=1,mm
            WRITE(IOUNIT,REC=IHOTSTP) GLOELV(I,N)
            IHOTSTP=IHOTSTP+1
         end do
      end do
      
      return
      end subroutine

c***********************************************************************
c   Subroutine to write elevation station harmonic analysis RHS load   *
c   vector to a hot start file (UNITS 67 and 68)                       *
c                                                                      *
c                        R.L.  11/8/95                                 *
c***********************************************************************
c
      SUBROUTINE HAHOUTES(NSTAE,IOUNIT,IHOTSTP)
      implicit none
      INTEGER NSTAE,IOUNIT,IHOTSTP,N,I
c
c***** Write Out Station Elevation RHS load vector
c
      do n=1,NSTAE
         do i=1,mm
            WRITE(IOUNIT,REC=IHOTSTP) STAELV(I,N)
            IHOTSTP=IHOTSTP+1
         end do
      end do

      return
      end subroutine

c***********************************************************************
c   Subroutine to write global velocity harmonic analysis RHS load     *
c   vector to a hot start file (UNITS 67 and 68)                       *
c                                                                      *
c                        R.L.  11/8/95                                 *
c***********************************************************************
c
      SUBROUTINE HAHOUTVG(NP,IOUNIT,IHOTSTP)
      implicit none
      INTEGER NP,IOUNIT,IHOTSTP,N,I
c
c***** Write Out Global Velocity RHS load vector
c
      do n=1,np
         do i=1,mm
            WRITE(IOUNIT,REC=IHOTSTP) GLOULV(I,N)
            IHOTSTP=IHOTSTP+1
            WRITE(IOUNIT,REC=IHOTSTP) GLOVLV(I,N)
            IHOTSTP=IHOTSTP+1
         end do
      end do
      
      return
      end subroutine

c***********************************************************************
c   Subroutine to write velocity station harmonic analysis RHS load    *
c   vector to a hot start file (UNITS 67 and 68)                       *
c                                                                      *
c                        R.L.  11/8/95                                 *
c***********************************************************************
c
      SUBROUTINE HAHOUTVS(NSTAV,IOUNIT,IHOTSTP)
      implicit none
      INTEGER NSTAV,IOUNIT,IHOTSTP,N,I
c
c***** Write Out Station Velocity LHS load vector
c
      do N=1,NSTAV
         do i=1,mm
            WRITE(IOUNIT,REC=IHOTSTP) STAULV(I,N)
            IHOTSTP=IHOTSTP+1
            WRITE(IOUNIT,REC=IHOTSTP) STAVLV(I,N)
            IHOTSTP=IHOTSTP+1
         end do
      end do

      return
      end subroutine


      END MODULE HARM
