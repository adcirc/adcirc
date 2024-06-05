!-------------------------------------------------------------------------------!
!
! ADCIRC - The ADvanced CIRCulation model
! Copyright (C) 1994-2023 R.A. Luettich, Jr., J.J. Westerink
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

      MODULE HARM

      USE SIZES, ONLY : MNP, MNHARF, MNSTAE, MNSTAV, LOCALDIR,
     &                   WRITE_LOCAL_HARM_FILES
      USE GLOBAL, ONLY : screenMessage, logMessage, DEBUG, ECHO, INFO,
     &                   WARNING, ERROR, setMessageSource, IFSPROTS,
     &                   unsetMessageSource, allMessage, scratchMessage
      USE ADC_CONSTANTS, ONLY: pi


      ! jgf49.44: parameters describing the harmonic constituents to be
      ! included in the harmonic analysis of model results
      INTEGER :: IHARIND ! =1 if harmonic analysis was requested, 0 otherwise
      LOGICAL :: CHARMV=.false.  ! .true. for timeseries reconstruction
      INTEGER :: NFREQ                      ! number of frequencies
      CHARACTER*10,ALLOCATABLE :: NAMEFR(:) ! "M2", "S2", etc
      REAL(8), ALLOCATABLE :: HAFREQ(:)    ! frequency (rad/s)
      REAL(8), ALLOCATABLE :: HAFF(:)      ! nodal factor
      REAL(8), ALLOCATABLE :: HAFACE(:)    ! equilibrium argument (deg)
C
      INTEGER :: NZ ! nz=0 if there is a steady component, nz=1 otherwise
      INTEGER :: NF ! nf=0 if no steady constituent, nf=1 otherwise
C
      INTEGER :: ITHAS  ! time step upon which harmonic analysis starts
      INTEGER :: ITHAF  ! time step upon which harmonic analysis finishes
      INTEGER :: ITMV   ! time step upon which means and variance calcs start
      INTEGER :: TIMEBEG ! model time (sec) at which means and var calcs start
      INTEGER :: IHABEG ! 1 h.a. time increment past the start of h.a.
      INTEGER :: ICHA   ! spool counter for HA updates
      INTEGER :: NTSTEPS! number of time steps into means and variance calcs
      INTEGER :: ICALL  ! number of times HA update sub has been called
      INTEGER :: ITUD   ! time step of most recent HA update
      REAL(8) :: TIMEUD ! model time of most recent HA update
C
      REAL(8), ALLOCATABLE :: HA(:,:) ! LHS matrix
      INTEGER :: MM ! 2nd dim of HA matrix; 2x num freqs (+ steady)
      REAL(8), ALLOCATABLE :: HAP(:)  ! used in matrix solution
      REAL(8), ALLOCATABLE :: HAX(:)  ! used in matrix solution
C
C     Load Vectors
      REAL(8), ALLOCATABLE, TARGET :: GLOELV(:,:) ! fulldomain elevation
      REAL(8), ALLOCATABLE, TARGET :: GLOULV(:,:) ! fulldomain u velocity
      REAL(8), ALLOCATABLE, TARGET :: GLOVLV(:,:) ! fulldomain v velocity
      REAL(8), ALLOCATABLE, TARGET :: STAELV(:,:) ! station elevation
      REAL(8), ALLOCATABLE, TARGET :: STAULV(:,:) ! station u velocity
      REAL(8), ALLOCATABLE, TARGET :: STAVLV(:,:) ! station v velocity
C
C     Arrays for time series reconstruction.
      REAL(8),ALLOCATABLE,TARGET :: XVELAV(:)
      REAL(8),ALLOCATABLE,TARGET :: YVELAV(:)
      REAL(8),ALLOCATABLE,TARGET :: XVELVA(:)
      REAL(8),ALLOCATABLE,TARGET :: YVELVA(:)
      REAL(8),ALLOCATABLE,TARGET :: ELAV(:)
      REAL(8),ALLOCATABLE,TARGET :: ELVA(:)
C
C     jgf49.44: Parameters that control the spatial locations where
C     harmonic analysis is performed:
      INTEGER :: NHASE ! /= 0 to perform HA at elevation recording stations
      INTEGER :: NHASV ! /= 0 to perform HA at velocity recording stations
      INTEGER :: NHAGE ! /= 0 to do elevation HA at all nodes in the mesh
      INTEGER :: NHAGV ! /= 0 to do velocity HA at all nodes in the mesh
C
C     jgf49.44: Parameters that control the calculation of harmonic
C     constituents:
      REAL(8) :: THAS ! number of days b/f harmonic analysis starts
      REAL(8) :: THAF ! number of days at which harmonic analysis ends
      INTEGER :: NHAINC ! time step increment for including data in HA
      REAL(8) :: FMV  ! fraction of HA period for means and var. calcs (0to1)
C
C     jgf49.44: Harmonic analysis hotstarting:
      REAL(8), ALLOCATABLE, TARGET :: GLOELV_g(:,:)
      REAL(8), ALLOCATABLE, TARGET :: STAELV_g(:,:)
      REAL(8), ALLOCATABLE, TARGET :: GLOULV_g(:,:)
      REAL(8), ALLOCATABLE, TARGET :: GLOVLV_g(:,:)
      REAL(8), ALLOCATABLE, TARGET :: STAULV_g(:,:)
      REAL(8), ALLOCATABLE, TARGET :: STAVLV_g(:,:)
      REAL(8), ALLOCATABLE, TARGET :: ELAV_g(:)
      REAL(8), ALLOCATABLE, TARGET :: ELVA_g(:)
      REAL(8), ALLOCATABLE, TARGET :: XVELAV_g(:)
      REAL(8), ALLOCATABLE, TARGET :: YVELAV_g(:)
      REAL(8), ALLOCATABLE, TARGET :: XVELVA_g(:)
      REAL(8), ALLOCATABLE, TARGET :: YVELVA_g(:)
C
C     jgf49.44: Equivalence variables for reading in hotstart data, either
C     in serial or in parallel.
      REAL(8), POINTER :: STAELV_pt(:,:),STAULV_pt(:,:),STAVLV_pt(:,:)
      REAL(8), POINTER :: GLOELV_pt(:,:),GLOULV_pt(:,:),GLOVLV_pt(:,:)
      REAL(8), POINTER :: ELAV_pt(:), ELVA_pt(:)
      REAL(8), POINTER :: XVELAV_pt(:), YVELAV_pt(:)
      REAL(8), POINTER :: XVELVA_pt(:), YVELVA_pt(:)
C
C     jgf49.44: The following variables are read from the hotstart file and
C     later compared with the values from the fort.15 to ensure that they
C     match.
      INTEGER INHASE, INHASV, INHAGE, INHAGV, IFLAG
      INTEGER INFREQ, INSTAE, INSTAV, INP, INZ, INF, IMM, IICALL
      REAL(8),ALLOCATABLE ::  IFREQ(:),IFF(:),IFACE(:)
      CHARACTER*10,ALLOCATABLE :: INAMEFR(:)
C
      CHARACTER*16 FNAME
      CHARACTER*8 FNAM8(2)
      EQUIVALENCE (FNAM8(1),FNAME)
C
C     jgf49.44: Raised output variables to module level and added
C     full domain counterparts for globalio. These arrays represent
C     magnitudes and phases for each frequency at each station or node.
C     Stations:
      REAL(8), ALLOCATABLE, TARGET :: EMAG(:,:)    ! elevation magnitudes
      REAL(8), ALLOCATABLE, TARGET :: PHASEDE(:,:) ! elevation phases
      REAL(8), ALLOCATABLE, TARGET :: UMAG(:,:)    ! u velocity magnitudes
      REAL(8), ALLOCATABLE, TARGET :: VMAG(:,:)    ! v velocity magnitudes
      REAL(8), ALLOCATABLE, TARGET :: PHASEDU(:,:) ! u velocity phases
      REAL(8), ALLOCATABLE, TARGET :: PHASEDV(:,:) ! v velocity phases
C     Nodes:
      REAL(8), ALLOCATABLE, TARGET :: EMAGT(:,:)    ! elevation magnitudes
      REAL(8), ALLOCATABLE, TARGET :: PHASEDEN(:,:) ! elevation phases
      REAL(8), ALLOCATABLE, TARGET :: UMAGT(:,:)    ! u velocity magnitudes
      REAL(8), ALLOCATABLE, TARGET :: VMAGT(:,:)    ! v velocity magnitudes
      REAL(8), ALLOCATABLE, TARGET :: PHASEDUT(:,:) ! u velocity phases
      REAL(8), ALLOCATABLE, TARGET :: PHASEDVT(:,:) ! v velocity phases
C     Time series reconstruction:
      REAL(8), ALLOCATABLE, TARGET :: EAV(:)
      REAL(8), ALLOCATABLE, TARGET :: ESQ(:)
      REAL(8), ALLOCATABLE, TARGET :: EAVDIF(:)
      REAL(8), ALLOCATABLE, TARGET :: EVADIF(:)
C
      REAL(8), ALLOCATABLE, TARGET :: UAV(:)
      REAL(8), ALLOCATABLE, TARGET :: USQ(:)
      REAL(8), ALLOCATABLE, TARGET :: UAVDIF(:)
      REAL(8), ALLOCATABLE, TARGET :: UVADIF(:)
C
      REAL(8), ALLOCATABLE, TARGET :: VAV(:)
      REAL(8), ALLOCATABLE, TARGET :: VSQ(:)
      REAL(8), ALLOCATABLE, TARGET :: VAVDIF(:)
      REAL(8), ALLOCATABLE, TARGET :: VVADIF(:)

      !jgf51.52.01: Variables used in adcpost.
      REAL(8), ALLOCATABLE,TARGET :: XELAV(:)
      REAL(8), ALLOCATABLE,TARGET :: YELAV(:)
      REAL(8), ALLOCATABLE,TARGET :: XELVA(:)
      REAL(8), ALLOCATABLE,TARGET :: YELVA(:)
     

C-----------------END OF DECLARATIONS---------------------------------------

      CONTAINS
C
C***********************************************************************
C  Allocate arays used by LSQ_HARM.
C
C  vjp 8/99
C***********************************************************************
      SUBROUTINE ALLOC_HA()
      IMPLICIT NONE
C
      call setMessageSource("ALLOC_HA")
#if defined(HARM_TRACE) || defined (ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif
      ALLOCATE ( HAFREQ(MNHARF),HAFF(MNHARF),HAFACE(MNHARF) )
      ALLOCATE ( NAMEFR(MNHARF) )
C
      ALLOCATE ( HA(2*MNHARF,2*MNHARF) )
      ALLOCATE ( HAP(2*MNHARF),HAX(2*MNHARF) )
      ALLOCATE ( GLOELV(2*MNHARF,MNP) )
      ALLOCATE ( GLOULV(2*MNHARF,MNP),GLOVLV(2*MNHARF,MNP) )
      ALLOCATE ( STAELV(2*MNHARF,MNSTAE) )
      ALLOCATE ( STAULV(2*MNHARF,MNSTAV),STAVLV(2*MNHARF,MNSTAV) )
C
      ALLOCATE ( IFREQ(MNHARF),IFF(MNHARF),IFACE(MNHARF) )
      ALLOCATE ( INAMEFR(MNHARF) )
C
      ALLOCATE ( EMAG(MNHARF,MNSTAE), PHASEDE(MNHARF,MNSTAE) )
      ALLOCATE ( UMAG(MNHARF,MNSTAV), PHASEDU(MNHARF,MNSTAV) )
      ALLOCATE ( VMAG(MNHARF,MNSTAV), PHASEDV(MNHARF,MNSTAV) )
      ALLOCATE ( EMAGT(MNHARF,MNP), PHASEDEN(MNHARF,MNP) )
      ALLOCATE ( UMAGT(MNHARF,MNP), PHASEDUT(MNHARF,MNP) )
      ALLOCATE ( VMAGT(MNHARF,MNP), PHASEDVT(MNHARF,MNP) )
C
#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN
C----------------------------------------------------------------------
      END SUBROUTINE ALLOC_HA
C----------------------------------------------------------------------



C----------------------------------------------------------------------
C...
C...Allocate space for harmonic analysis means and variance
C...calculations
C...
C----------------------------------------------------------------------
      SUBROUTINE ALLOC_MAIN14()
      IMPLICIT NONE
C
      call setMessageSource("ALLOC_MAIN14")
#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
      ALLOCATE ( XVELAV(MNP),YVELAV(MNP),XVELVA(MNP),YVELVA(MNP) )
      ALLOCATE ( ELAV(MNP),ELVA(MNP) )
C
      ALLOCATE( EAV(MNP),ESQ(MNP),EAVDIF(MNP),EVADIF(MNP) )
      ALLOCATE( UAV(MNP),USQ(MNP),UAVDIF(MNP),UVADIF(MNP) )
      ALLOCATE( VAV(MNP),VSQ(MNP),VAVDIF(MNP),VVADIF(MNP) )
C

#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN
C----------------------------------------------------------------------
      END SUBROUTINE ALLOC_MAIN14
C----------------------------------------------------------------------



C--------------------------------------------------------------------
C                    S U B R O U T I N E
C     A L L O C A T E   F U L L   D O M A I N   H A   I O   A R R A Y S
C--------------------------------------------------------------------
C     jgf49.44: Allocates memory for the full domain arrays for i/o
C     purposes when executing in parallel.
C--------------------------------------------------------------------
      SUBROUTINE allocateFullDomainHAIOArrays()
      USE GLOBAL, ONLY : NSTAE_G, NSTAV_G, NP_G
      IMPLICIT NONE
C
      call setMessageSource("allocateFullDomainHAIOArrays")
#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      IF (IHARIND.eq.1) THEN
         ALLOCATE(GLOELV_g(2*MNHARF,NP_G))
         ALLOCATE(STAELV_g(2*MNHARF,NSTAE_G))
         ALLOCATE(GLOULV_g(2*MNHARF,NP_G))
         ALLOCATE(GLOVLV_g(2*MNHARF,NP_G))
         ALLOCATE(STAULV_g(2*MNHARF,NSTAV_G))
         ALLOCATE(STAVLV_g(2*MNHARF,NSTAV_G))
      ENDIF
      IF (CHARMV) THEN
         ALLOCATE(ELAV_g(NP_G))
         ALLOCATE(ELVA_g(NP_G))
         ALLOCATE(XVELAV_g(NP_G))
         ALLOCATE(YVELAV_g(NP_G))
         ALLOCATE(XVELVA_g(NP_G))
         ALLOCATE(YVELVA_g(NP_G))
      ENDIF
#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C--------------------------------------------------------------------
      END SUBROUTINE allocateFullDomainHAIOArrays
C--------------------------------------------------------------------

C--------------------------------------------------------------------
C     S U B R O U T I N E   W R I T E   H A 
C--------------------------------------------------------------------
C     jgf50.97: Writes LHS matrix to a file for use in debugging.
C--------------------------------------------------------------------
      SUBROUTINE writeHA(timesec, it)
      USE SIZES, ONLY : globaldir, myproc
      IMPLICIT NONE     
      REAL(8) :: timesec ! time in seconds since cold start
      INTEGER :: it      ! current adcirc time step number
      INTEGER :: i, j    ! loop counters
C
      call setMessageSource("writeHA")
#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      IF (myproc.eq.0) THEN
         OPEN(56,FILE=TRIM(GLOBALDIR)//'/'//'fort.56',
     &      STATUS='UNKNOWN',ACCESS='SEQUENTIAL',ACTION='WRITE',
     &      POSITION='APPEND')
         WRITE(56,2120) timesec, it 
         DO i=1,2*MNHARF
            DO j=1,2*MNHARF
               write(56,2130) i, j, ha(i,j)
            END DO
         END DO
         CLOSE(56)
      ENDIF

 2120 FORMAT(2X,1pE20.10E3,5X,I10)
 2130 FORMAT('ha(',I2,',',I2,')=',1pE20.10E3)

#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C--------------------------------------------------------------------
      END SUBROUTINE writeHA
C--------------------------------------------------------------------


C--------------------------------------------------------------------
C     S U B R O U T I N E  C H E C K   H A R M O N I C   P A R A M E T E R S
C--------------------------------------------------------------------
C     jgf50.41: Checks harmonic analysis settings.
C     WJP 02.20.2018: Allowing for > 1 options for netcdf output
C--------------------------------------------------------------------
      SUBROUTINE checkHarmonicParameters()
      USE GLOBAL, ONLY :  NFOVER, useNetCDF, useNetCDFOutput
      IMPLICIT NONE
C
      call setMessageSource("checkHarmonicParameters")
#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      IF((NHASE.LT.0).OR.(NHASE.GT.5)) THEN
         CALL allMessage(WARNING,"Input value of NHASE is not valid.")
         IF (NFOVER.EQ.1) THEN
            CALL allMessage(WARNING,"NHASE will be set to zero.")
         ELSE
            CALL harmonicTerminate()
         ENDIF
      ELSE IF (NHASE.EQ.1) THEN
         CALL logMessage(INFO,
     &        'STATION ELEVATION HARMONIC ANALYSIS'
     &        //' WILL BE WRITTEN TO UNIT 51 IN ASCII FORMAT.')
      ELSEIF (NHASE > 1) THEN
         CALL logMessage(INFO,
     &        'STATION ELEVATION HARMONIC ANALYSIS'
     &        //' WILL BE WRITTEN TO UNIT 51 IN NETCDF FORMAT.')
         useNetCDF = .true. 
         useNetCDFOutput = .true.
      ENDIF

      IF((NHASV.LT.0).OR.(NHASV.GT.5)) THEN
         CALL allMessage(WARNING,"Input value of NHASV is not valid.")
         IF (NFOVER.EQ.1) THEN
            CALL allMessage(WARNING,"NHASV will be set to zero.")
         ELSE
            CALL harmonicTerminate()
         ENDIF
      ELSE IF (NHASV.EQ.1) THEN
         CALL logMessage(INFO,
     &        'STATION VELOCITY HARMONIC ANALYSIS'
     &        //' WILL BE WRITTEN TO UNIT 52 IN ASCII FORMAT.')
      ELSE IF (NHASV > 1) THEN
         CALL logMessage(INFO,
     &        'STATION VELOCITY HARMONIC ANALYSIS'
     &        //' WILL BE WRITTEN TO UNIT 52 IN NETCDF FORMAT.')
         useNetCDF = .true. 
         useNetCDFOutput = .true.
      ENDIF

      IF((NHAGE.LT.0).OR.(NHAGE.GT.5)) THEN
         CALL allMessage(WARNING,"Input value of NHAGE is not valid.")
         IF (NFOVER.EQ.1) THEN
            CALL allMessage(WARNING,"NHAGE will be set to zero.")
         ELSE
            CALL harmonicTerminate()
         ENDIF
      ELSE IF (NHAGE.EQ.1) THEN
         CALL logMessage(INFO,
     &        'FULL DOMAIN ELEVATION HARMONIC ANALYSIS'
     &        //' WILL BE WRITTEN TO UNIT 53 IN ASCII FORMAT.')
      ELSE IF (NHAGE > 1) THEN
         CALL logMessage(INFO,
     &        'FULL DOMAIN ELEVATION HARMONIC ANALYSIS'
     &        //' WILL BE WRITTEN TO UNIT 53 IN NETCDF FORMAT.')
         useNetCDF = .true. 
         useNetCDFOutput = .true.
      ENDIF

      IF((NHAGV.LT.0).OR.(NHAGV.GT.5)) THEN
         CALL allMessage(WARNING,"Input value of NHAGV is not valid.")
         IF (NFOVER.EQ.1) THEN
            CALL allMessage(WARNING,"NHAGV will be set to zero.")
         ELSE
            CALL harmonicTerminate()
         ENDIF
      ELSE IF (NHAGV.EQ.1) THEN
         CALL logMessage(INFO,
     &        'FULL DOMAIN VELOCITY HARMONIC ANALYSIS'
     &        //' WILL BE WRITTEN TO UNIT 54 ASCII FORMAT.')
      ELSE IF (NHAGV > 1) THEN
         CALL logMessage(INFO,
     &        'FULL DOMAIN VELOCITY HARMONIC ANALYSIS'
     &        //' WILL BE WRITTEN TO UNIT 54 NETCDF FORMAT.')
         useNetCDF = .true. 
         useNetCDFOutput = .true.
      ENDIF

#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C--------------------------------------------------------------------
      END SUBROUTINE checkHarmonicParameters
C--------------------------------------------------------------------


C--------------------------------------------------------------------
C     S U B R O U T I N E  I N I T   H A R M O N I C   P A R A M E T E R S
C--------------------------------------------------------------------
C     jgf49.44: Initializes harmonic analysis arrays and parameters.
C     Based on HACOLDS, written by rl.
C     WJP 02.20.2018: Allowing for > 1 options for netcdf output
C--------------------------------------------------------------------
      SUBROUTINE initHarmonicParameters()
      USE SIZES, ONLY : READ_LOCAL_HOT_START_FILES, MNPROC, myProc,
     &                   WRITE_LOCAL_HARM_FILES, 
     &                   WRITE_LOCAL_HOT_START_FILES
      USE GLOBAL, ONLY : ITHS, NSTAE, NSTAV, NHSTAR
      USE MESH, ONLY : NP
      IMPLICIT NONE
      INTEGER :: i, j, n           ! loop counters
C
      call setMessageSource("initHarmonicParameters")
#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      IF (IHARIND.EQ.0) THEN

#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
         call unsetMessageSource()
         RETURN ! EARLY RETURN if harmonic analysis is not part of this run
      ENDIF
      ! 
      ! jgf51.44: If subdomain harmonic analysis files will be written
      ! in parallel, subdomain hotstart files are also required for
      ! saving the load vectors and rhs of the harmonic analysis problem.
      if ( (mnproc.gt.1).and.
     &     (write_local_harm_files.eqv..true.) ) then
         call logMessage(INFO,'Subdomain harmonic analysis files '
     &                // 'were specified with -m on the command line.')
         if ( (NHSTAR.ne.0).and.
     &        (write_local_hot_start_files.eqv..false.) ) then  
            call logMessage(INFO,'Fulldomain hotstart files will be '
     &                            // 'created by default.')
            call allMessage(ERROR,'Fulldomain hotstart output files '
     &       // 'cannot be used in combination with subdomain harmonic '
     &       // 'analysis output files. '
     &       // 'The -s argument can be used on the ADCIRC command '
     &       // 'line to specify that hotstart files should be written '
     &       // 'for each subdomain; this is required for the writing '
     &       // 'of harmonic analysis output (fort.51 etc) '
     &       // ' to subdomains.') 
             call harmonicTerminate()
         endif
         if ( (NHSTAR.ne.0).and. 
     &        (write_local_hot_start_files.eqv..true.) ) then
            call logMessage(INFO,'Subdomain hotstart files will '
     &       // 'contain subdomain harmonic analysis data.')
         endif
      endif
      
      IHABEG = ITHAS + NHAINC
      ICHA=0
      icall=0
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
C
      ha(:,:)=0.0d0
C
      IF (NHASE.ne.0) THEN
         STAELV(:,:)=0.d0
      ENDIF
      IF (NHASV.ne.0) THEN
         STAULV(:,:)=0.d0
         STAVLV(:,:)=0.d0
      ENDIF
      IF (NHAGE.ne.0) THEN
         GLOELV(:,:)=0.d0
      ENDIF
      IF (NHAGV.ne.0) THEN
         GLOULV(:,:)=0.d0
         GLOVLV(:,:)=0.d0
      ENDIF
      IF (CHARMV.eqv..true.) THEN
         ELAV(:)=0.D0
         XVELAV(:)=0.D0
         YVELAV(:)=0.D0
         ELVA(:)=0.D0
         XVELVA(:)=0.D0
         YVELVA(:)=0.D0
      ENDIF                  !  charmv
C
C     jgf49.44: Initialize all the full domain i/o arrays to zero, if
C     they will be used (i.e., if we are in serial or if we might be
C     reading or writing full domain arrays in parallel).
      if ( (WRITE_LOCAL_HARM_FILES.eqv..false.).and.
     &     (READ_LOCAL_HOT_START_FILES.eqv..false.).and.
     &     (MNPROC.gt.1).and.
     &     (myProc.eq.0) ) then
         GLOELV_g(:,:) = 0.d0
         STAELV_g(:,:) = 0.d0
         GLOULV_g(:,:) = 0.d0
         GLOVLV_g(:,:) = 0.d0
         STAULV_g(:,:) = 0.d0
         STAVLV_g(:,:) = 0.d0
         IF (CHARMV.eqv..true.) THEN
            ELAV_g(:) = 0.d0
            ELVA_g(:) = 0.d0
            XVELAV_g(:) = 0.d0
            YVELAV_g(:) = 0.d0
            XVELVA_g(:) = 0.d0
            YVELVA_g(:) = 0.d0
         ENDIF
      endif
#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C--------------------------------------------------------------------
      end subroutine initHarmonicParameters
C--------------------------------------------------------------------

C--------------------------------------------------------------------
C     S U B R O U T I N E  U P D A T E   H A R M O N I C  A N A L Y S I S
C--------------------------------------------------------------------
C     jgf49.44: Removed from timestep.F and placed here.
C     WJP 02.20.2018: Allowing for > 1 options for netcdf output
C...
C...  IF IHARIND=1 AND THE TIME STEP FALLS WITHIN THE SPECIFIED WINDOW
C...  AND ON THE SPECIFIED INCREMENT, USE MODEL RESULTS TO UPDATE
C...  HARMONIC ANALYSIS MATRIX AND LOAD VECTORS.  NOTE: AN 8 BYTE RECORD
C...  SHOULD BE USED THROUGHOUT THE HARMONIC ANALYSIS SUBROUTINES, EVEN
C...  ON 32 BIT WORKSTATIONS, SINCE IN THAT CASE THE HARMONIC ANALYSIS
C...  IS DONE IN DOUBLE PRECISION.
C...
C--------------------------------------------------------------------
      SUBROUTINE updateHarmonicAnalysis(it, timeh)
      USE GLOBAL, ONLY : ET00, UU00, VV00, ETA2, UU2, VV2,
     &                   STAIE1, STAIE2, STAIE3, STAIV1, STAIV2, STAIV3,
     &                   NSTAE, NSTAV, NNE, NNV
      USE MESH, ONLY : NP, NM, UVECTMP, VVECTMP, DRVMAP2DSPVEC
      IMPLICIT NONE
      INTEGER, intent(in) :: it    ! current ADCIRC timestep
      REAL(8), intent(in) :: timeh ! time (sec) incl. STATIM and REFTIM
C
      INTEGER :: i, j              ! loop counters
      REAL(8) :: EE1, EE2, EE3    ! nodal elevations around an element
      REAL(8) :: U11, U22, U33    ! nodal u velocities around an element
      REAL(8) :: V11, V22, V33    ! nodal v velocities around an element
      INTEGER :: n                 ! loop counter
C
      call setMessageSource("updateHarmonicAnalysis")
#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C...
      IF(IHARIND.EQ.1) THEN
         IF((IT.GT.ITHAS).AND.(IT.LE.ITHAF)) THEN
            IF(ICHA.EQ.NHAINC) ICHA=0
            ICHA=ICHA+1
            IF(ICHA.EQ.NHAINC) THEN
C...
C.....UPDATE THE LHS MATRIX
C...
C.... BEG: DW 2019       
               IF ( (NHASV .EQ. 1) .OR. (NHAGV .EQ. 1) ) THEN 
                  IF ( IFSPROTS .EQ. 1 ) THEN
                     UVECTMP = UU2 ;
                     VVECTMP = VV2 ;
                     
                     CALL DRVMAP2DSPVEC( UU2, VV2,
     &                    UVECTMP, VVECTMP, NP, FWD = .FALSE. )   
                  END IF
               END IF
C.... DW

               CALL LSQUPDLHS(timeh,IT)

C...  IF DESIRED COMPUTE ELEVATION STATION INFORMATION AND UPDATE LOAD
C.....VECTOR
C...
               IF(NHASE.ne.0) THEN
                  DO I=1,NSTAE
                     EE1=ETA2(NM(NNE(I),1))
                     EE2=ETA2(NM(NNE(I),2))
                     EE3=ETA2(NM(NNE(I),3))
                     ET00(I)=EE1*STAIE1(I)+EE2*STAIE2(I)+EE3*STAIE3(I)
                  END DO
                  CALL LSQUPDES(ET00,NSTAE)
               ENDIF
C...  IF DESIRED COMPUTE VELOCITY STATION INFORMATION AND UPDATE LOAD
C.....VECTOR
C...
               IF(NHASV.ne.0) THEN
                  DO I=1,NSTAV
                     U11=UU2(NM(NNV(I),1))
                     U22=UU2(NM(NNV(I),2))
                     U33=UU2(NM(NNV(I),3))
                     V11=VV2(NM(NNV(I),1))
                     V22=VV2(NM(NNV(I),2))
                     V33=VV2(NM(NNV(I),3))
                     UU00(I)=U11*STAIV1(I)+U22*STAIV2(I)+U33*STAIV3(I)
                     VV00(I)=V11*STAIV1(I)+V22*STAIV2(I)+V33*STAIV3(I)
                  END DO
                  CALL LSQUPDVS(UU00,VV00,NSTAV)
               ENDIF
C...
C.....IF DESIRED UPDATE GLOBAL ELEVATION LOAD VECTOR
C...
               IF(NHAGE.ne.0) CALL LSQUPDEG(ETA2,NP)
C...
C.....IF DESIRED UPDATE GLOBAL VELOCITY LOAD VECTOR
C...
               IF(NHAGV.ne.0) CALL LSQUPDVG(UU2,VV2,NP)

C.... BEG: DW 2019       
               IF ( (NHASV .EQ. 1) .OR. (NHAGV .EQ. 1) ) THEN 
                  IF ( IFSPROTS .EQ. 1 ) THEN
                     UU2 = UVECTMP ;
                     VV2 = VVECTMP ;
                  END IF
               END IF
C.... DW

            ENDIF
         ENDIF

C...  LINES TO COMPUTE MEANS AND VARIANCES

         if (CHARMV) then
            IF(IT.GT.ITMV) THEN
               NTSTEPS=NTSTEPS+1
               DO I=1,NP
                  ELAV(I)=ELAV(I)+ETA2(I)
                  XVELAV(I)=XVELAV(I)+UU2(I)
                  YVELAV(I)=YVELAV(I)+VV2(I)
                  ELVA(I)=ELVA(I)+ETA2(I)*ETA2(I)
                  XVELVA(I)=XVELVA(I)+UU2(I)*UU2(I)
                  YVELVA(I)=YVELVA(I)+VV2(I)*VV2(I)
               END DO
            ENDIF
         endif                  !   charmv


      ENDIF



#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C--------------------------------------------------------------------
      end subroutine updateHarmonicAnalysis
C--------------------------------------------------------------------


c***********************************************************************
c   Subroutine to update the Left Hand Side Matrix                     *
c                                                                      *
c  TIMELOC  - ABSOLUTE MODEL TIME (SEC)                                   *
c  IT    - MODEL TIME STEP                                             *
c  icall - number of times the subroutine has been called              *
c  a     - Left Hand Side Matrix                                       *
c                                                                      *
c                        RL 11/7/95                                    *
c***********************************************************************
c
      SUBROUTINE LSQUPDLHS(TIMELOC,IT)
      IMPLICIT NONE
      INTEGER IT,I,J,I1,I2,J1,J2
      REAL(8) TF1,TF2
      REAL(8) TIMELOC
C
      call setMessageSource("LSQUPDLHS")
#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
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
            tf1=hafreq(j+nf)*TIMELOC
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
            tf1=hafreq(i+nf)*TIMELOC
            tf2=hafreq(j+nf)*TIMELOC
            ha(i1,j1) = ha(i1,j1) + cos(tf1)*cos(tf2)
            ha(i1,j2) = ha(i1,j2) + cos(tf1)*sin(tf2)
            ha(i2,j2) = ha(i2,j2) + sin(tf1)*sin(tf2)
            if(i2.le.j1) ha(i2,j1) = ha(i2,j1) + sin(tf1)*cos(tf2)
         end do
      end do

c   Record update time and time step

      TIMEUD = TIMELOC
      ITUD = IT

#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      return
      end subroutine LSQUPDLHS

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
      REAL(8) TF1,CTF1,STF1
      REAL(8) STAE(MNSTAE)
C
      call setMessageSource("LSQUPDES")
#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
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
#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      return
      end subroutine LSQUPDES

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
      REAL(8) TF1,CTF1,STF1
      REAL(8) STAU(MNSTAV),STAV(MNSTAV)
C
      call setMessageSource("LSQUPDVS")
#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
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
#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      return
      end subroutine LSQUPDVS


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
      INTEGER I,NP,N,I1,I2
      REAL(8) TF1,CTF1,STF1
      REAL(8) GLOE(MNP)
C
      call setMessageSource("LSQUPDEG")
#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
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

#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      return
      end subroutine LSQUPDEG


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
      INTEGER NP,I1,I2,N,I,J
      REAL(8) TF1,CTF1,STF1
      REAL(8) GLOU(MNP),GLOV(MNP)
C
      call setMessageSource("LSQUPDVG")
#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
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

#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      return
      end subroutine LSQUPDVG

C--------------------------------------------------------------------
C     S U B R O U T I N E  S O L V E   H A R M O N I C   A N A L Y S I S
C--------------------------------------------------------------------
C     jgf49.44: Solves for amplitudes and phases at the nodes and/or
C     stations as specified. Also finishes time series reconstruction
C     if necessary.
C     WJP 02.20.2018: Allowing for > 1 options for netcdf output
C--------------------------------------------------------------------
      SUBROUTINE solveHarmonicAnalysis(ITIME)
      USE SIZES, ONLY : LOCALDIR
      USE GLOBAL, ONLY : STATIM, REFTIM, DTDP
      USE MESH, ONLY : NP
      IMPLICIT NONE
      INTEGER, intent(in) :: ITIME
      INTEGER :: I
C
      call setMessageSource("solveHarmonicAnalysis")
#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
      IF ((IHARIND.NE.1).OR.(ITIME.LE.ITHAS)) THEN
#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
         call unsetMessageSource()
         RETURN ! EARLY RETURN if we weren't supposed to do harmonic analysis
                ! or it hasn't started yet.
      ENDIF
C
C...Compute means and variances for checking the harmonic analysis results
C...Accumulate mean and variance at each node.
      if ((CHARMV.eqv..true.).and.(FMV.gt.1.0d-3)) then
         DO I=1,NP
            ELAV(I)   = ELAV(I)/NTSTEPS
            XVELAV(I) = XVELAV(I)/NTSTEPS
            YVELAV(I) = YVELAV(I)/NTSTEPS
            ELVA(I)   = ELVA(I)/NTSTEPS   - ELAV(I)*ELAV(I)
            XVELVA(I) = XVELVA(I)/NTSTEPS - XVELAV(I)*XVELAV(I)
            YVELVA(I) = YVELVA(I)/NTSTEPS - YVELAV(I)*YVELAV(I)
         END DO
         TIMEBEG=ITMV*DTDP + (STATIM-REFTIM)*86400.D0
      ENDIF

C......Fill out and decompose the LHS harmonic analaysis matrix

      CALL FULSOL(0)

C......Solve the harmonic analysis problems

      IF(NHAGE.ne.0) CALL LSQSOLEG()
      IF(NHAGV.ne.0) CALL LSQSOLVG()
      IF(NHASE.ne.0) CALL LSQSOLES()
      IF(NHASV.ne.0) CALL LSQSOLVS()

#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C--------------------------------------------------------------------
      END SUBROUTINE solveHarmonicAnalysis
C--------------------------------------------------------------------



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
      real(8),allocatable ::  c(:),y(:)
C
      call setMessageSource("fulsol")
#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
c
c**** If only want to fill out matrix and decompose
c
      if(idecom.eq.0) then

c     Set up the lower triangular part of the LHS a matrix

         do j=1,mm
            do i=j,mm
               ha(i,j)=ha(j,i)
            end do
         end do

c     Decomposition of matrix a

         outer: do ir=1,mm
            ire=ir+1
            do j=ire,mm
               ha(ir,j)=ha(ir,j)/ha(ir,ir)
            enddo
            if (ire.gt.mm) exit outer
            do j=ire,mm
              do k=ire,mm
                 ha(k,j)=ha(k,j)-ha(k,ir)*ha(ir,j)
              enddo
            enddo
            do j=ire,mm
              ha(j,ir)=0.0d0
            enddo
         enddo outer

#if defined(HARM_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG,"Return.")
#endif
         call unsetMessageSource()
         return ! EARLY RETURN if idecom.eq.0
      endif

c...  solve for y by forward substitution for l*y=p

      allocate ( c(2*MNHARF),y(2*MNHARF) )
c
      do ir=1,mm
         y(ir)=hap(ir)
         do jr=1,ir-1
            y(ir)=y(ir)-ha(jr,ir)*y(jr)
         enddo
      enddo

c...  calculate c=d**(-1)*y

      do ir=1,mm
         c(ir)=y(ir)/ha(ir,ir)
      enddo

c...  solve for x by back-substituting for l(tr)*x=c

      ir=mm
      do while(.true.)
         hax(ir)=c(ir)
         do jr=ir+1,mm
            hax(ir)=hax(ir)-ha(ir,jr)*hax(jr)
         enddo
         ir=ir-1
         if(ir.ge.1) cycle
         exit
      enddo

#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      end subroutine fulsol


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
      SUBROUTINE LSQSOLES()
      USE GLOBAL, ONLY : NSTAE
      IMPLICIT NONE
      INTEGER N,I,J,K,I1,I2
      REAL(8) CONVRD
      REAL(8) PHASEE
C
      call setMessageSource("LSQSOLES")
#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
      convrd=180.d0/pi
c
c**** AT each STATION TRANSFER each load vector to p and solve the system
c
      DO N=1,NSTAE
         do k=1,mm
            hap(k)=STAELV(k,n)
         end do
         call fulsol(n)
c
c        Compute amplitude and phase for each frequency making sure that the
c        phase is between 0 and 360 deg.
         do i=1,nfreq+nf
            if((nf.eq.1).and.(i.eq.1)) then
               emag(i,n)=hax(i)/haff(i)
               phasee=0.d0
            else
               i1=2*i-1-nf
               i2=i1+1
               emag(i,n)=sqrt(hax(i1)*hax(i1)+hax(i2)*hax(i2))/haff(i)
               if((hax(i1).eq.0.).and.(hax(i2).eq.0.)) then
                  phasee=0.d0
               else
                  phasee = atan2(hax(i2),hax(i1))
               endif
            endif
            phasede(i,n)=convrd*phasee+haface(i)
            if(phasede(i,n).lt.0.d0) phasede(i,n)=phasede(i,n)+360.d0
            if(phasede(i,n).ge.360.d0) phasede(i,n)=phasede(i,n)-360.d0
         end do
      end do
C

#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      return
      end subroutine lsqsoles


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
      SUBROUTINE LSQSOLVS()
      USE GLOBAL, ONLY : NSTAV
      IMPLICIT NONE
      INTEGER I,J,N,K,I1,I2
      REAL(8) CONVRD
      REAL(8) PHASEU,PHASEV
      REAL(8),ALLOCATABLE :: Y(:)
C
      call setMessageSource("LSQSOLVS")
#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
c
      convrd=180.d0/pi
c
      ALLOCATE ( Y(2*MNHARF) )
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
c        phase is between 0 and 360 deg.

         do i=1,nfreq+nf
            if((nf.eq.1).and.(i.eq.1)) then
               umag(i,n)=hax(i)/haff(i)
               vmag(i,n)=y(i)/haff(i)
               phaseu=0.
               phasev=0.
            else
               i1=2*i-1-nf
               i2=i1+1
               umag(i,n)=sqrt(hax(i1)*hax(i1)+hax(i2)*hax(i2))/haff(i)
               vmag(i,n)=sqrt(y(i1)*y(i1)+y(i2)*y(i2))/haff(i)
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
            phasedu(i,n)=convrd*phaseu+haface(i)
            if(phasedu(i,n).lt.0.) phasedu(i,n)=phasedu(i,n)+360.d0
            if(phasedu(i,n).ge.360.d0) phasedu(i,n)=phasedu(i,n)-360.d0
            phasedv(i,n)=convrd*phasev+haface(i)
            if(phasedv(i,n).lt.0.) phasedv(i,n)=phasedv(i,n)+360.d0
            if(phasedv(i,n).ge.360.d0) phasedv(i,n)=phasedv(i,n)-360.d0
         end do
      end do

#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      return
      end subroutine lsqsolvs

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
      SUBROUTINE LSQSOLEG()
      USE GLOBAL, ONLY : DTDP
      IMPLICIT NONE
      integer J,N,K,I,I1,I2,IT,IFR,NEAVMAX,NEAVMIN,
     &  NEVAMAX,NEVAMIN
      REAL(8)  CONVRD
      REAL(8) EAVMAX,EVAMAX,EAVMIN,EVAMIN
      REAL(8) TIMELOC,RSE,FTIME
      REAL(8),ALLOCATABLE  ::  PHASEE(:),EMAG(:)
C
      call setMessageSource("LSQSOLEG")
#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
      convrd=180.d0/pi
c
      ALLOCATE ( PHASEE(MNHARF),EMAG(MNHARF) )
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
      DO N=1,MNP
         do k=1,mm
            hap(k) = GLOELV(k,n)
         end do
         call fulsol(n)
c
c        Compute amplitude and phase for each frequency making sure that the
c        phase is between 0 and 360 deg.  Then write output.
c
         do i=1,nfreq+nf
            if((nf.eq.1).and.(i.eq.1)) then
               emag(i)=hax(i)
               emagt(i,n)=emag(i)/haff(i)
               phasee(i)=0.
            else
               i1=2*i-1-nf
               i2=i1+1
               emag(i)=sqrt(hax(i1)*hax(i1)+hax(i2)*hax(i2))
               emagt(i,n)=emag(i)/haff(i)
               if((hax(i1).eq.0.).and.(hax(i2).eq.0.)) then
                  phasee(i)=0.
               else
                  phasee(i) = atan2(hax(i2),hax(i1))
               endif
            endif
            phaseden(i,n)=convrd*phasee(i)+haface(i)
            if(phaseden(i,n).lt.0.) phaseden(i,n)=phaseden(i,n)+360.d0
            if(phaseden(i,n).ge.360.d0) phaseden(i,n)=phaseden(i,n)-360.d0
         end do

         if (CHARMV) then
            eav(n) = 0.d0
            esq(n) = 0.d0
            do it=1,ntsteps
               TIMELOC=TIMEBEG+DTDP*IT
               rse=0.d0
               do ifr=1,nfreq+nf
                  ftime=hafreq(ifr)*TIMELOC
                  rse=rse+emag(ifr)*cos(ftime-phasee(ifr))
               end do
               eav(n)=eav(n)+rse
               esq(n)=esq(n)+rse*rse
            end do
C
            eav(n)=eav(n)/ntsteps
            esq(n)=esq(n)/ntsteps-eav(n)*eav(n)
            if(elav(n).eq.0.) then
               if(eav(n).eq.0.) eavdif(n)=1.0d0
               if(eav(n).ne.0.) eavdif(n)=99d19
            else
               eavdif(n)=eav(n)/elav(n)
            endif
            if(elva(n).eq.0.) then
               if(esq(n).eq.0.) evadif(n)=1.0d0
               if(esq(n).ne.0.) evadif(n)=99e19
            else
               evadif(n)=esq(n)/elva(n)
            endif
C
            IF(EAVDIF(n).GT.EAVMAX) THEN
               EAVMAX=EAVDIF(n)
               NEAVMAX=n
            ENDIF
            IF(EAVDIF(n).LT.EAVMIN) THEN
               EAVMIN=EAVDIF(n)
               NEAVMIN=n
            ENDIF
            IF(EVADIF(n).GT.EVAMAX) THEN
               EVAMAX=EVADIF(n)
               NEVAMAX=n
            ENDIF
            IF(EVADIF(n).LT.EVAMIN) THEN
               EVAMIN=EVADIF(n)
               NEVAMIN=n
            ENDIF
         endif                     ! charmv
      end do
C
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

#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      return
      end subroutine lsqsoleg

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
      SUBROUTINE LSQSOLVG()
      USE GLOBAL, ONLY : DTDP
      IMPLICIT NONE
      INTEGER I,J,N,K,I1,I2,IT,IFR
      INTEGER NUAVMAX,NUAVMIN,NVAVMAX,NVAVMIN,NUVAMAX,NUVAMIN,
     &  NVVAMAX,NVVAMIN
      REAL(8) TIMELOC,FTIME,RSU,RSV
      REAL(8) UAVMAX,VAVMAX,UVAMAX,VVAMAX,UAVMIN,VAVMIN,
     & UVAMIN,VVAMIN
      REAL(8) CONVRD
      REAL(8),ALLOCATABLE :: UMAGG(:),VMAGG(:),PHASEUG(:),PHASEVG(:)
      REAL(8),ALLOCATABLE :: Y(:)
C
      call setMessageSource("LSQSOLVG")
#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
c
      convrd=180.d0/pi
c
      ALLOCATE ( UMAGG(MNHARF),VMAGG(MNHARF) )
      ALLOCATE ( PHASEUG(MNHARF),PHASEVG(MNHARF) )
      ALLOCATE ( Y(2*MNHARF) )
c
c
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
      DO N=1,MNP
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
         do i=1,nfreq+nf
            if((nf.eq.1).and.(i.eq.1)) then
               umagg(i)=hax(i)
               umagt(i,n)=umagg(i)/haff(i)
               vmagg(i)=y(i)
               vmagt(i,n)=vmagg(i)/haff(i)
               phaseug(i)=0.d0
               phasevg(i)=0.d0
            else
               i1=2*i-1-nf
               i2=i1+1
               umagg(i)=sqrt(hax(i1)*hax(i1)+hax(i2)*hax(i2))
               umagt(i,n)=umagg(i)/haff(i)
               vmagg(i)=sqrt(y(i1)*y(i1)+y(i2)*y(i2))
               vmagt(i,n)=vmagg(i)/haff(i)
               if((hax(i1).eq.0.).and.(hax(i2).eq.0.)) then
                  phaseug(i)=0.
               else
                  phaseug(i)=atan2(hax(i2),hax(i1))
               endif
               if((y(i1).eq.0.).and.(y(i2).eq.0.)) then
                  phasevg(i)=0.
               else
                  phasevg(i)=atan2(y(i2),y(i1))
               endif
            endif
            phasedut(i,n)=convrd*phaseug(i)+haface(i)
            if(phasedut(i,n).lt.0.) phasedut(i,n)=phasedut(i,n)+360.d0
            if(phasedut(i,n).ge.360.d0) phasedut(i,n)=phasedut(i,n)-360.d0
            phasedvt(i,n)=convrd*phasevg(i)+haface(i)
            if(phasedvt(i,n).lt.0.) phasedvt(i,n)=phasedvt(i,n)+360.d0
            if(phasedvt(i,n).ge.360.d0) phasedvt(i,n)=phasedvt(i,n)-360.d0
 6636       format(2x,e16.8,1x,f11.4,2x,e16.8,1x,f11.4)
         end do

CHARMV...UNCOMMENT THE FOLLOWING LINES TO COMPUTE MEANS AND VARIANCES
CHARMV...FOR CHECKING THE HARMONIC ANALYSIS RESULTS.
CHARMV...Resynthesize the time series to compute the average and variances.
CHARMV...Compare resynthesized values with those computed during time stepping.
         if ( charmv ) then
            uav(n) = 0.d0
            vav(n) = 0.d0
            usq(n) = 0.d0
            vsq(n) = 0.d0
            do it=1,ntsteps
               TIMELOC=TIMEBEG+DTDP*IT
               rsu=0.
               rsv=0.
               do ifr=1,nfreq+nf
                  ftime=hafreq(ifr)*TIMELOC
                  rsu=rsu+umagg(ifr)*cos(ftime-phaseug(ifr))
                  rsv=rsv+vmagg(ifr)*cos(ftime-phasevg(ifr))
               end do
               uav(n)=uav(n)+rsu
               vav(n)=vav(n)+rsv
               usq(n)=usq(n)+rsu*rsu
               vsq(n)=vsq(n)+rsv*rsv
            end do

            uav(n)=uav(n)/ntsteps
            vav(n)=vav(n)/ntsteps
            usq(n)=usq(n)/ntsteps-uav(n)*uav(n)
            vsq(n)=vsq(n)/ntsteps-vav(n)*vav(n)
            if(xvelav(n).eq.0.) then
               if(uav(n).eq.0.) uavdif(n)=1.0d0
               if(uav(n).ne.0.) uavdif(n)=99e19
            else
               uavdif(n)=uav(n)/xvelav(n)
            endif
            if(yvelav(n).eq.0.) then
               if(vav(n).eq.0.) vavdif(n)=1.0d0
               if(vav(n).ne.0.) vavdif(n)=99e19
            else
               vavdif(n)=vav(n)/yvelav(n)
            endif
            if(xvelva(n).eq.0.) then
               if(usq(n).eq.0.) uvadif(n)=1.0d0
               if(usq(n).ne.0.) uvadif(n)=99e19
            else
               uvadif(n)=usq(n)/xvelva(n)
            endif
            if(yvelva(n).eq.0.) then
               if(vsq(n).eq.0.) vvadif(n)=1.0d0
               if(vsq(n).ne.0.) vvadif(n)=99e19
            else
               vvadif(n)=vsq(n)/yvelva(n)
            endif

            IF(UAVDIF(n).GT.UAVMAX) THEN
               UAVMAX=UAVDIF(n)
               NUAVMAX=n
            ENDIF
            IF(UAVDIF(n).LT.UAVMIN) THEN
               UAVMIN=UAVDIF(n)
               NUAVMIN=n
            ENDIF
            IF(VAVDIF(n).GT.VAVMAX) THEN
               VAVMAX=VAVDIF(n)
               NVAVMAX=n
            ENDIF
            IF(VAVDIF(n).LT.VAVMIN) THEN
               VAVMIN=VAVDIF(n)
               NVAVMIN=n
            ENDIF
            IF(UVADIF(n).GT.UVAMAX) THEN
               UVAMAX=UVADIF(n)
               NUVAMAX=n
            ENDIF
            IF(UVADIF(n).LT.UVAMIN) THEN
               UVAMIN=UVADIF(n)
               NUVAMIN=n
            ENDIF
            IF(VVADIF(n).GT.VVAMAX) THEN
               VVAMAX=VVADIF(n)
               NVVAMAX=n
            ENDIF
            IF(VVADIF(n).LT.VVAMIN) THEN
               VVAMIN=VVADIF(n)
               NVVAMIN=n
            ENDIF

         endif                  !  charmv

      end do

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

#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      return
      end subroutine lsqsolvg



C--------------------------------------------------------------------
C     S U B R O U T I N E  R E A D  B I N A R Y  H A  H O T S T A R T
C--------------------------------------------------------------------
C     jgf49.44: Reads harmonic analysis data from a binary hotstart
C     file. Based on HAHOTS, written by rl.
C     WJP 02.20.2018: Allowing for > 1 options for netcdf output
C--------------------------------------------------------------------
      SUBROUTINE readBinaryHAHotstart(lun, counter)
      USE SIZES, ONLY : MNPROC, READ_LOCAL_HOT_START_FILES
      USE GLOBAL, ONLY : NSTAE, NSTAE_G, NSTAV, NSTAV_G,
     &                   ITHS, IMAP_STAE_LG, NP_G, NODES_LG,
     &                   IMAP_STAV_LG
      USE MESH, ONLY : NP
      IMPLICIT NONE
      INTEGER, intent(in) :: lun          ! i/o logical unit number
      INTEGER, intent(inout) :: counter   ! i/o record
      INTEGER :: i, j, n                  ! loop counters
      INTEGER :: num_elev_sta   ! number of elevation stations to read
      INTEGER :: num_vel_sta    ! number of velocity stations to read
      INTEGER :: num_nodes      ! number of nodes to read
      INTEGER :: subdomainStation  ! number of station to map data to
      INTEGER :: fulldomainStation ! number of station to map data from
      INTEGER :: subdomainNode  ! number of node to map data to
      INTEGER :: fulldomainNode ! number of node to map data from
C
      call setMessageSource("readBinaryHAHotstart")
#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
C     Settings for reading in hotstart file in serial, or for reading
C     subdomain hotstart files in parallel.
      IF ((MNPROC.eq.1)
     &   .or.(READ_LOCAL_HOT_START_FILES.eqv..true.)) THEN
         STAELV_pt => STAELV
         num_elev_sta = NSTAE
         STAULV_pt => STAULV
         STAVLV_pt => STAVLV
         num_vel_sta = NSTAV
         GLOELV_pt => GLOELV
         num_nodes = NP
         GLOULV_pt => GLOULV
         GLOVLV_pt => GLOVLV
         ELAV_pt => ELAV
         ELVA_pt => ELVA
         XVELAV_pt => XVELAV
         YVELAV_pt => YVELAV
         XVELVA_pt => XVELVA
         YVELVA_pt => YVELVA
      ELSE
C        ! read full domain hotstart file in parallel
         STAELV_pt => STAELV_G
         num_elev_sta = NSTAE_G
         STAULV_pt => STAULV_G
         STAVLV_pt => STAVLV_G
         num_vel_sta = NSTAV_G
         GLOELV_pt => GLOELV_G
         num_nodes = NP_G
         GLOULV_pt => GLOULV_G
         GLOVLV_pt => GLOVLV_G
         ELAV_pt => ELAV_g
         ELVA_pt => ELVA_g
         XVELAV_pt => XVELAV_g
         YVELAV_pt => YVELAV_g
         XVELVA_pt => XVELVA_g
         YVELVA_pt => YVELVA_g
      ENDIF
c
c***** Read in parameter values
c
      READ(lun,REC=counter+1) inz
      READ(lun,REC=counter+2) inf
      READ(lun,REC=counter+3) imm
      READ(lun,REC=counter+4) inp
      READ(lun,REC=counter+5) instae
      READ(lun,REC=counter+6) instav
      READ(lun,REC=counter+7) inhase
      READ(lun,REC=counter+8) inhasv
      READ(lun,REC=counter+9) inhage
      READ(lun,REC=counter+10) inhagv
      READ(lun,REC=counter+11) iicall
      READ(lun,REC=counter+12) infreq
      counter = counter+12
c
      do i=1,nfreq+nf
         READ(lun,REC=counter+1) FNAM8(1)
         READ(lun,REC=counter+2) FNAM8(2)
         counter = counter + 2
         INAMEFR(I) = FNAME
         read(lun,REC=counter+1) ifreq(i)
         read(lun,REC=counter+2) iff(i)
         read(lun,REC=counter+3) iface(i)
         counter = counter + 3
      end do
c
c***** Read in time of most recent H.A. update
c
      READ(lun,REC=counter+1) TIMEUD
      READ(lun,REC=counter+2) ITUD
      counter = counter + 2
c
c     Read in LHS Matrix
      counter = counter + 1
      
      ! jgf50.96: do not use binaryRead3D(HA, mm, mm, lun, counter)
      ! b/c HA is written with an inner j loop while binaryRead3D is
      ! written with an inner i loop
      do i=1,mm
         do j=1,mm
            READ(lun,REC=counter) HA(I,J)
            counter = counter + 1
         end do
      end do
C
C
c     Read in Station Elevation LHS load vector
      IF (NHASE.ne.0) THEN
         CALL binaryRead3D(STAELV_pt, mm, num_elev_sta, lun, counter)
      ENDIF
c     Read in Station Velocity LHS load vector
      IF (NHASV.ne.0) THEN
         do n=1,num_vel_sta
            do i=1,mm
               READ(lun,REC=counter) STAULV_pt(I,N)
               READ(lun,REC=counter+1) STAVLV_pt(I,N)
               counter = counter + 2
            enddo
         enddo
      ENDIF
c     Read in Global Elevation LHS load vector
      IF (NHAGE.ne.0) THEN
         CALL binaryRead3D(GLOELV_pt, mm, num_nodes, lun, counter)
      ENDIF
c     Read in Global Velocity LHS load vector
      IF (NHAGV.ne.0) THEN
         do n=1,num_nodes
            do i=1,mm
               READ(lun,REC=counter) GLOULV_pt(I,N)
               READ(lun,REC=counter+1) GLOVLV_pt(I,N)
               counter = counter + 2
            end do
         end do
      ENDIF
C..   Read in Means and Squares
      IF (CHARMV.eqv..true.) THEN
         IF ((FMV.NE.0.).AND.(ITHS.GT.ITMV)) THEN
            READ(lun,REC=counter) NTSTEPS ; counter = counter + 1
            IF(NHAGE.ne.0) THEN
               DO I=1,num_nodes
                  READ(lun,REC=counter) ELAV_pt(I)
                  counter = counter + 1
                  READ(lun,REC=counter) ELVA_pt(I)
                  counter = counter + 1
               ENDDO
            ENDIF
            IF (NHAGV.ne.0) THEN
               DO I=1,num_nodes
                  READ(lun,REC=counter) XVELAV_pt(I)
                  counter = counter + 1
                  READ(lun,REC=counter) YVELAV_pt(I)
                  counter = counter + 1
                  READ(lun,REC=counter) XVELVA_pt(I)
                  counter = counter + 1
                  READ(lun,REC=counter) YVELVA_pt(I)
                  counter = counter + 1
               ENDDO
            ENDIF
         ENDIF
      ENDIF   !  charmv
C
C     Map data from fulldomain to subdomain in parallel
      IF ((MNPROC.gt.1)
     &   .and.(READ_LOCAL_HOT_START_FILES.eqv..false.)) THEN
         IF (NHASE.ne.0) THEN
            DO subdomainStation=1,NSTAE
               fulldomainStation = ABS(IMAP_STAE_LG(subdomainStation))
               DO i=1,mm
                  STAELV(i,subdomainStation) =
     &               STAELV_g(i,fulldomainStation)
               END DO
            END DO
         ENDIF
         IF (NHASV.ne.0) THEN
            DO subdomainStation=1,NSTAV
               fulldomainStation = ABS(IMAP_STAV_LG(subdomainStation))
               DO i=1,mm
                  STAULV(i,subdomainStation) =
     &               STAULV_g(i,fulldomainStation)
                  STAVLV(i,subdomainStation) =
     &               STAVLV_g(i,fulldomainStation)
               END DO
            END DO
         ENDIF
         IF (NHAGE.ne.0) THEN
            DO subdomainNode=1,NP
               fulldomainNode = ABS(NODES_LG(subdomainNode))
               DO i=1,mm
                  GLOELV(i,subdomainNode) = GLOELV_g(i,fulldomainNode)
               ENDDO
            ENDDO
         ENDIF
         IF (NHAGV.ne.0) THEN
            DO subdomainNode=1,NP
               fulldomainNode = ABS(NODES_LG(subdomainNode))
               DO i=1,mm
                  GLOULV(i,subdomainNode) = GLOULV_g(i,fulldomainNode)
                  GLOVLV(i,subdomainNode) = GLOVLV_g(i,fulldomainNode)
               ENDDO
            ENDDO
         ENDIF
         
C..      Map Means and Squares
         IF (CHARMV.eqv..true.) THEN
            IF ((FMV.NE.0.).AND.(ITHS.GT.ITMV)) THEN
               IF(NHAGE.ne.0) THEN
                  DO subdomainNode=1,NP
                     fulldomainNode = ABS(NODES_LG(subdomainNode))
                     ELAV(subdomainNode) = ELAV_g(fulldomainNode)
                     ELVA(subdomainNode) = ELVA_g(fulldomainNode)
                  ENDDO
               ENDIF
               IF (NHAGV.ne.0) THEN
                  DO subdomainNode=1,NP
                     fulldomainNode = ABS(NODES_LG(subdomainNode))
                     XVELAV(subdomainNode) = XVELAV_g(fulldomainNode)
                     YVELAV(subdomainNode) = YVELAV_g(fulldomainNode)
                     XVELVA(subdomainNode) = XVELVA_g(fulldomainNode)
                     YVELVA(subdomainNode) = YVELVA_g(fulldomainNode)
                  END DO
               ENDIF
            ENDIF
         ENDIF
      ENDIF ! MNPROC.gt.1 etc

#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C--------------------------------------------------------------------
      end subroutine readBinaryHAHotstart
C--------------------------------------------------------------------

C--------------------------------------------------------------------
C     S U B R O U T I N E  C H E C K   H A  H O T S T A R T
C--------------------------------------------------------------------
C     jgf49.44: Checks harmonic analysis data from a hotstart file
C     by comparing with parameters read from the fort.15.
C     Based on HAHOTS, written by rl.
C--------------------------------------------------------------------
      SUBROUTINE checkHAHotstart()
      USE SIZES, ONLY : READ_LOCAL_HOT_START_FILES, MYPROC
      USE GLOBAL, ONLY : NSTAE_G, NSTAV_G, NSTAE, NSTAV,
     &                   NSCREEN, SCREENUNIT, NP_G
      USE MESH, ONLY : NP
      IMPLICIT NONE
      INTEGER :: i    ! loop counter
      REAL(8) FDIFF  ! difference between frequencies btw hotstart and fort.15
      CHARACTER(len=1024) :: message
C
      call setMessageSource("checkHAHotstart")
#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      call logMessage(INFO,
     & "Comparing harmonic parameters in hotstart file to fort.15...")
C
      iflag=0
      if(nz.ne.inz) iflag=1
      if(nf.ne.inf) iflag=1
      if(mm.ne.imm) iflag=1
      if ((myProc.eq.0).and.(read_local_hot_start_files.eqv..false.)) then
         if(np_g.ne.inp) iflag=1
         if(nstae_g.ne.instae) iflag=1
         if(nstav_g.ne.instav) iflag=1
      else if (READ_LOCAL_HOT_START_FILES.eqv..true.) then
         if(np.ne.inp) iflag=1
         if(nstae.ne.instae) iflag=1
         if(nstav.ne.instav) iflag=1
      endif
      if(nhase.ne.inhase) iflag=1
      if(nhasv.ne.inhasv) iflag=1
      if(nhage.ne.inhage) iflag=1
      if(nhagv.ne.inhagv) iflag=1
      if(nfreq.ne.infreq) iflag=1
c
      do i=1,nfreq+nf
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
C
      if(iflag.eq.0) then
         call logMessage(INFO,
     &      "Harmonic data in hotstart file matches fort.15. PASSED.")
      else
c
c***** FATAL Error Messages
c
         write(message,1000)
         call allMessage(ERROR, message)
 1000    FORMAT('***** DISCREPANCY IN HARMONIC ANALYSIS HOT ',
     +        'START FILE *****')
         if(nz.ne.inz) then
            write(message,2010) nz,inz
            call allMessage(ERROR, message)
 2010       format('NZ COMPUTED FROM UNIT 14 INPUT = ',I2,
     +           ', NZ READ IN FROM HOT START FILE = ',I2)
         endif
         if(nf.ne.inf) then
            write(message,2010) nf, inf
            call allMessage(ERROR, message)
 2020       format('NF COMPUTED FROM UNIT 14 INPUT = ',I2,
     +           ', NF READ IN FROM HOT START FILE = ',I2)
         endif
         if(mm.ne.imm) then
            write(message,2030) mm, imm
            call allMessage(ERROR, message)
 2030       format('MM COMPUTED FROM UNIT 14 INPUT = ',I2,
     +           ', MM READ IN FROM HOT START FILE = ',I2)
         endif
         If (myproc.eq.0) then
            if(np_g.ne.inp) then
               write(message,2040) np_g, inp
               call allMessage(ERROR, message)
 2040          format('NP READ IN FROM UNIT 14 = ',I2,
     +              ', NP READ IN FROM HOT START FILE = ',I2)
            endif
            if(nstae_g.ne.instae) then
               write(message,2050) nstae_g, instae
               call allMessage(ERROR, message)
 2050          format('NSTAE READ IN FROM UNIT 15 = ',I2,
     +           ', NSTAE READ IN FROM HOT START FILE = ',I2)
            endif
            if(nstav_g.ne.instav) then
               write(message,2060) nstav_g, instav
               call allMessage(ERROR, message)
 2060          format('NSTAV READ IN FROM UNIT 15 = ',I2,
     +              ', NSTAV READ IN FROM HOT START FILE = ',I2)
            endif
         endif
         if(nhase.ne.inhase) then
            write(message,2070) NHASE, INHASE
            call allMessage(ERROR, message)
 2070       format('NHASE READ IN FROM UNIT 15 = ',I2,
     +           ', NHASE READ IN FROM HOT START FILE = ',I2)
         endif
         if(nhasv.ne.inhasv) then
            write(message,2080) NHASV, INHASV
            call allMessage(ERROR, message)
 2080       format('NHASV READ IN FROM UNIT 15 = ',I2,
     +           ', NHASV READ IN FROM HOT START FILE = ',I2)
         endif
         if(nhage.ne.inhage) then
            write(message,2090) NHAGE, INHAGE
            call allMessage(ERROR, message)
 2090       format('NHAGE READ IN FROM UNIT 15 = ',I2,
     +           ', NHAGE READ IN FROM HOT START FILE = ',I2)
         endif
         if(nhagv.ne.inhagv) then
            write(message,2100) NHAGV, INHAGV
            call allMessage(ERROR, message)
 2100       format('NHAGV READ IN FROM UNIT 15 = ',I2,
     +           ', NHAGV READ IN FROM HOT START FILE = ',I2)
         endif
         if(nfreq.ne.infreq) then
            write(message,2110) NFREQ,INFREQ
            call allMessage(ERROR, message)
 2110       format('NFREQ COMPUTED FROM UNIT 15 INPUT = ',I2,
     +           ', NFREQ READ IN FROM HOT START FILE = ',I2)
         endif
         do i=1,nfreq+nf
            if(namefr(i).ne.inamefr(i)) then
               write(message,2120) i,namefr(i),inamefr(i)
               call allMessage(ERROR, message)
 2120          format('FOR CONSTITUENT # ',I3,
     +              ', NAMEFR READ IN FROM UNIT 15 = ',A10,
     +              ', NAMEFR READ IN FROM HOT START FILE = ',A10)
            endif
            if(hafreq(i).ne.ifreq(i)) then
               write(message,2130) i,hafreq(i),ifreq(i)
               call allMessage(ERROR, message)
 2130          format('FOR CONSTITUENT # ',I3,
     +              ', FREQ READ IN FROM UNIT 15 = ',D20.10,
     +              ', FREQ READ IN FROM HOT START FILE = ',D20.10)
            endif
            if(HAFF(i).ne.iFF(i)) then
               write(message,2140) i,haff(i),iff(i)
               call allMessage(ERROR, message)
 2140          format('FOR CONSTITUENT # ',I3,
     +              ', FF READ IN FROM UNIT 15 = ',F10.5,
     +              ', FF READ IN FROM HOT START FILE = ',F10.5)
            endif
            if(HAFACE(i).ne.iFACE(i)) then
               write(message,2150) i,haface(i),iface(i)
               call allMessage(ERROR, message)
 2150          format('FOR CONSTITUENT # ',I3,
     +              ', FACE READ IN FROM UNIT 15 = ',F10.5,
     +              ', FACE READ IN FROM HOT START FILE = ',F10.5)
            endif
         end do
         call harmonicTerminate()
      endif

#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C--------------------------------------------------------------------
      end subroutine checkHAHotstart
C--------------------------------------------------------------------

C----------------------------------------------------------------------
C...  jgf50.41: Subroutine to terminate the run cleanly.
C...
C----------------------------------------------------------------------
      SUBROUTINE harmonicTerminate()
#ifdef CMPI
      USE MESSENGER, ONLY : MSG_FINI
#endif
      IMPLICIT NONE
C
      call setMessageSource("harmonicTerminate")
#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
      call allMessage(ERROR,"ADCIRC terminating.")
#ifdef CMPI
      call msg_fini()
#endif
      CALL EXIT(1)
C

#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN
C----------------------------------------------------------------------
      END SUBROUTINE harmonicTerminate
C----------------------------------------------------------------------

C--------------------------------------------------------------------
C     S U B R O U T I N E      B I N A R Y   R E A D   3 D
C--------------------------------------------------------------------
C     jgf49.44: Same as binaryRead2D except it reads data with
C     both horizontal and vertical lengths. Also used for reading
C     harmonic analysis arrays from binary hotstart file.
C--------------------------------------------------------------------
      SUBROUTINE binaryRead3D(array, iLength, jLength, lun,
     &                        counter)
      USE GLOBAL, ONLY : setMessageSource, unsetMessageSource, 
     &                    allMessage, DEBUG
      IMPLICIT NONE
      REAL(8), intent(out), dimension(ilength,jlength) :: array ! data to read from the hotstart file
      INTEGER, intent(in) :: iLength  ! the number of horiz values to read
      INTEGER, intent(in) :: jLength  ! the number of layers
      INTEGER, intent(in) :: lun     ! fortran logical unit number to read from
      INTEGER, intent(inout) :: counter ! i/o position in the file
C
      INTEGER i,j      ! array indices
C
      call setMessageSource("binaryRead3D")
#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
       DO j=1, jLength
          DO i=1, iLength
            READ(lun,REC=counter) array(i,j)
            counter = counter + 1
          END DO
       END DO
C
#if defined(HARM_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C--------------------------------------------------------------------
      END SUBROUTINE binaryRead3D
C--------------------------------------------------------------------


C--------------------------------------------------------------
      END MODULE HARM
C--------------------------------------------------------------
