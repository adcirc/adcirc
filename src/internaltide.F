C---------------------------------------------------------------------C
      MODULE INTERNALTIDE
C---------------------------------------------------------------------C
C     CPB 03/2023: This module was created to make the
C     apply2dinternalwavedrag subroutine in nodalattr.F more readable.
C     For now I have simply moved the calculation of the de-tided
C     velocities here. I plan to move more of the subroutine into this
C     module to hopefully make it a little clearer what is happening in
C     that subroutine.
C---------------------------------------------------------------------C
      USE GLOBAL, ONLY : DEBUG, ECHO, INFO, WARNING, ERROR,
     &          setMessageSource, unsetMessageSource, allMessage,
     &          scratchMessage, logMessage
      IMPLICIT NONE
      ! stored velocity samples for moving average calculation
      REAL(8),ALLOCATABLE,PRIVATE :: UAV(:,:), VAV(:,:)        
      ! averaged velocities so we don't calculate them every single time
      ! step
      REAL(8),ALLOCATABLE :: UBar(:), VBar(:), UTide(:), VTide(:)
      ! weights for filters
      REAL(8),ALLOCATABLE,PRIVATE :: wts(:)
      ! values needed to resynthesize tides
      REAL(8),ALLOCATABLE,PRIVATE :: UTamp(:,:), UTphs(:,:),
     &                                VTamp(:,:), VTphs(:,:)
C---------------------------------------------------------------------C
      CONTAINS
C---------------------------------------------------------------------C
C---------------------------------------------------------------------C
      SUBROUTINE RESYNTH_TIDES(TimeStep)
C---------------------------------------------------------------------C
C     Resynthesizes tides read in from an external harmonic analysis to
C     use in the internal wave drag calculation
C---------------------------------------------------------------------C
      USE MESH, ONLY : NP, DRVMAP2DSPVEC, UVECTMP, VVECTMP
      USE GLOBAL, ONLY: DTDP, AMIGT, FFT, FACET, STatim, RefTim,
     &                  CHOTHS, ITHS, DTDPHS, PERT, IFSPROTS, NP_G
      USE SIZES, ONLY : MNTIF, MYPROC
      USE GL2LOC_MAPPING, ONLY : GATHERGLOBAL_REAL
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: TimeStep    ! current time step
      INTEGER :: NCYC ! for calculating argument to pass to cosine
      REAL(SZ) :: ARGT
      REAL(SZ) :: StaTim_tmp, StaTimHS, TimeH
      LOGICAL,SAVE :: first_call = .TRUE.
      INTEGER :: ii, jj ! for loops
      INTEGER,SAVE :: counter = 1
      REAL(SZ),ALLOCATABLE,DIMENSION(:) :: Udir, Vdir

      IF ( first_call ) THEN
         first_call = .FALSE.
         CALL READ_TIDESFILE()
         ALLOCATE( UTide(NP), VTide(NP) )
      ENDIF
      ! calculate time
      TimeH=TimeStep*DTDP + (StaTim - RefTim)*86400.D0
      StaTim_tmp = StaTim
      IF (CHOTHS.eqv..true.) THEN
         IF ((ITHS+1).EQ.TimeStep) THEN
           StaTimHS=((TimeStep-1)*DTDPHS)/86400.D0
           StaTim_tmp=((TimeStep-1)*DTDP)/86400.D0
         END IF
         TimeH=TimeStep*DTDP + ((StaTimHS - StaTim_tmp) - RefTim)*86400.D0
      END IF
      ! Resynthesize
      UTide = 0d0
      VTide = 0d0
      DO ii = 1,NP
         DO jj = 1,MNTIF
            ! calculate number of cycles in we are
            IF (PERT(jj).EQ.0.) THEN
               NCYC=0
            ELSE
#ifdef IBM
               NCYC=INT(timeh/PERT(jj),KIND(0.0d0))
#else
               NCYC=INT(timeh/PERT(jj))
#endif
            ENDIF
            ARGT = AMIGT(jj)*(TimeH - NCYC*PERT(jj)) +
     &             FACET(jj)
            UTide(ii) = UTide(ii) + 
     &                  FFT(jj)*UTamp(ii,jj)*COS(ARGT-UTphs(ii,jj))
            VTide(ii) = VTide(ii) + 
     &                  FFT(jj)*VTamp(ii,jj)*COS(ARGT-VTphs(ii,jj))
         ENDDO
      ENDDO

      IF (COUNTER.EQ.1) THEN
         IF (MYPROC.EQ.0) THEN
            ALLOCATE( UDir(NP_G), VDir(NP_G) )
         ELSE
            ALLOCATE( UDir(1), VDir(1) )
         ENDIF
         CALL GATHERGLOBAL_REAL(UDir,UTide)
         CALL GATHERGLOBAL_REAL(VDir,VTide)
         IF (MYPROC.EQ.0) THEN
            OPEN(111,ACCESS='SEQUENTIAL',ACTION='WRITE',
     &               FILE='tidevelo_ts0.txt',STATUS='REPLACE')
            WRITE(111,*) NP_G
            WRITE(111,*) COUNTER*DTDP
            DO jj=1,NP_G
               WRITE(111,7638) jj, UDir(jj), VDir(jj)
            ENDDO
            CLOSE(111)
 7638 FORMAT(2x,i8,2x,2E20.5E4)
         ENDIF
         DEALLOCATE( UDir, VDir )
         counter = counter + 1
      ELSE
         COUNTER = COUNTER + 1
      ENDIF
      
      IF (IFSPROTS.EQ.1) THEN
         UVECTMP(1:NP) = UTide(1:NP)
         VVECTMP(1:NP) = VTide(1:NP)
         CALL DRVMAP2DSPVEC( UTide, VTide, UVECTMP, VVECTMP, NP, 
     &                      FWD = .TRUE. )
      ENDIF


C-----------------------------------------------------------------------
      END SUBROUTINE RESYNTH_TIDES
C-----------------------------------------------------------------------
C---------------------------------------------------------------------C
      SUBROUTINE READ_TIDESFILE()
C---------------------------------------------------------------------C
C     This subroutine takes as inputs the current time step velocity,
C     updates the average velocity vectors (if necessary) and outputs
C     the 25 hour lagged average velocity based on a lagged, 25 hour
C     filter with sampling frequency of 12 minutes. It replaces a good
C     chunk of code in the apply2dinternalwavedrag to make that
C     subroutine more readable as well as to eliminate some global
C     variables
C---------------------------------------------------------------------C
      USE MESH, ONLY : NP
      USE GLOBAL, ONLY: TIPOTAG, NP_G
      USE SIZES, ONLY : MYPROC, MNTIF
      USE GL2LOC_MAPPING, ONLY : MAPTOLOCAL_REAL
      USE CONSTANTS, ONLY : DEG2RAD
#ifdef ADCNETCDF
      USE NETCDF
#endif
      IMPLICIT NONE
      CHARACTER(LEN=80) :: TIDESFILE='tidesFile.nc'
      CHARACTER(LEN=10),ALLOCATABLE :: const(:)
      INTEGER :: NC_ID, Temp_id, ua_id, ug_id, va_id, vg_id
      INTEGER :: constlen, num_const, node
      INTEGER :: ii,jj ! looping vars
      REAL(SZ),ALLOCATABLE :: tmp_g(:,:), tmp_l(:)
      LOGICAL :: const_found

      ! determine if this tidesfile will actually work
#ifdef ADCNETCDF
      IF (MYPROC.EQ.0) THEN
         CALL check_err(NF90_OPEN(TIDESFILE,NF90_NOWRITE,NC_ID))
         CALL check_err(NF90_INQ_DIMID(NC_ID,'node',Temp_ID))
         CALL check_err(NF90_INQUIRE_DIMENSION(NC_ID,Temp_ID,len=node))
         IF (node.NE.NP_G) THEN
            WRITE(16,*) "wrong # of nodes yo"
            CALL ITFRICTERMINATE()
         ENDIF
         CALL check_err(NF90_INQ_DIMID(NC_ID,'num_const',Temp_ID))
         CALL check_err(NF90_INQUIRE_DIMENSION(NC_ID,TEMP_ID,
     &                                         len=num_const))
         ALLOCATE( const(num_const) )
         CALL check_err(NF90_INQ_VARID(NC_ID,'const',Temp_ID))
         CALL check_err(NF90_GET_VAR(NC_ID,Temp_ID,const))
         ! get variable id's
         CALL check_err(NF90_INQ_VARID(NC_ID,'u_amp',ua_id))
         CALL check_err(NF90_INQ_VARID(NC_ID,'u_phs',ug_id))
         CALL check_err(NF90_INQ_VARID(NC_ID,'v_amp',va_id))
         CALL check_err(NF90_INQ_VARID(NC_ID,'v_phs',vg_id))
      ENDIF
      ! allocate temporary arrays
      ALLOCATE( tmp_l(NP) )
      IF (MYPROC.EQ.0) THEN
         ALLOCATE( tmp_g(NP_G,4) )
      ELSE
         ALLOCATE( tmp_g(1,4) )
      ENDIF
      ! allocate UTamp, UTphs, VTamp, VTphs
      
      ALLOCATE( UTamp(NP,MNTIF), UTphs(NP,MNTIF), 
     &          VTamp(NP,MNTIF), VTphs(NP,MNTIF) )
      WRITE(16,*) "Allocated utamp"
      ! read in file
      DO ii = 1,MNTIF
         IF (MYPROC.EQ.0) THEN
            ! read in on proc 0
            const_found = .FALSE.
            DO jj = 1,num_const
               IF (TRIM(TIPOTAG(ii)).EQ.TRIM(const(jj))) THEN
                  WRITE(16,*) "TIPOTAG(ii) = ",trim(tipotag(ii))
                  write(16,*) "const(jj) = ",trim(const(jj))
                  const_found = .TRUE.
                  EXIT
               ENDIF
            ENDDO
            IF (const_found.EQV..FALSE.) THEN
               WRITE(16,*) "ERROR: Could not find constituent: ",
     &                     TRIM(TIPOTAG(ii))
               CALL ITFRICTERMINATE()
            ENDIF
            ! u_amp
            CALL check_err(NF90_GET_VAR(NC_ID,ua_id,tmp_g(:,1),
     &                                  start=[jj,1],count=[1,np_g]))
            ! u_phs
            CALL check_err(NF90_GET_VAR(NC_ID,ug_id,tmp_g(:,2),
     &                                  start=[jj,1],count=[1,np_g]))
            ! v_amp
            CALL check_err(NF90_GET_VAR(NC_ID,va_id,tmp_g(:,3),
     &                                  start=[jj,1],count=[1,np_g]))
            ! v_phs
            CALL check_err(NF90_GET_VAR(NC_ID,vg_id,tmp_g(:,4),
     &                                  start=[jj,1],count=[1,np_g]))
         ENDIF
         ! u_amp
         CALL MAPTOLOCAL_REAL(tmp_g(:,1),tmp_l)
         UTamp(:,ii) = tmp_l
         ! u_phs
         CALL MAPTOLOCAL_REAL(tmp_g(:,2),tmp_l)
         UTphs(:,ii) = tmp_l*DEG2RAD
         ! v_amp
         CALL MAPTOLOCAL_REAL(tmp_g(:,3),tmp_l)
         VTamp(:,ii) = tmp_l
         ! v_phs
         CALL MAPTOLOCAL_REAL(tmp_g(:,4),tmp_l)
         VTphs(:,ii) = tmp_l*DEG2RAD
      ENDDO

      DEALLOCATE( tmp_g, tmp_l )

      IF (MYPROC.EQ.0) THEN
         CALL check_err(NF90_CLOSE(NC_ID))
      ENDIF

#endif
C-----------------------------------------------------------------------
      END SUBROUTINE READ_TIDESFILE
C-----------------------------------------------------------------------
      SUBROUTINE UNTIDE(U_i,V_i,TimeStep)
C---------------------------------------------------------------------C
C     This subroutine takes as inputs the current time step velocity,
C     updates the average velocity vectors (if necessary) and outputs
C     the 25 hour lagged average velocity based on a lagged, 25 hour
C     filter with sampling frequency of 12 minutes. It replaces a good
C     chunk of code in the apply2dinternalwavedrag to make that
C     subroutine more readable as well as to eliminate some global
C     variables
C---------------------------------------------------------------------C
      USE MESH, ONLY : NP
      USE GLOBAL, ONLY: DTDP
      IMPLICIT NONE
      REAL(SZ),INTENT(IN),DIMENSION(:) :: U_i, V_i ! current velocity
      INTEGER,INTENT(IN) :: TimeStep    ! current time step
      REAL(SZ) :: filtL = 25d0*3600d0        ! filter length (s)
      REAL(SZ) :: Fs = 12d0*60d0             ! sampling interval (s)
      INTEGER,SAVE :: NS = 1   ! number of samples
      LOGICAL,SAVE :: first_call = .TRUE.
      ! need indices to check if the current and previous time step are
      ! in the same 12 min windo2
      INTEGER :: L, Lm           
      INTEGER :: ii, kk ! for loops
      ! for populating the UAV, VAV matrices at the start of the run
      INTEGER,SAVE :: ISTA

      IF ( first_call ) THEN
         first_call = .FALSE.
         NS = FLOOR(filtL/Fs)
         ista = 1
         ALLOCATE( UAV(NP,NS), VAV(NP,NS) )
         UAV = 0d0
         VAV = 0d0
         ALLOCATE( UBar(NP), VBar(NP) )
         UBar = 0d0
         VBar = 0d0
      ENDIF

      L = FLOOR(TimeStep*DTDP/Fs)
      Lm = FLOOR((TimeStep-1)*DTDP/Fs)
      IF ( L.GT.Lm ) THEN
         IF ( ISTA.GT.NS ) THEN
            DO ii = 1,NP
               DO kk = 1,NS-1
                  UAV(ii,kk) = UAV(ii,kk+1)
                  VAV(ii,kk) = VAV(ii,kk+1)
               ENDDO
               UAV(ii,NS) = U_i(ii)
               VAV(ii,NS) = V_i(ii)
            ENDDO
            DO ii = 1,NP
               UBar(ii) = SUM(UAV(ii,1:NS))/DBLE(NS)
               VBar(ii) = SUM(VAV(ii,1:NS))/DBLE(NS)
            ENDDO
         ELSE
            DO ii = 1,NP
               UAV(ii,ISTA) = U_i(ii)
               VAV(ii,ISTA) = V_i(ii)
               UBar(ii) = SUM(UAV(ii,1:ISTA))/DBLE(ISTA)
               VBar(ii) = SUM(VAV(ii,1:ISTA))/DBLE(ISTA)
            ENDDO
            ISTA = ISTA + 1
         ENDIF
      ENDIF
C-----------------------------------------------------------------------
      END SUBROUTINE UNTIDE
C-----------------------------------------------------------------------
C---------------------------------------------------------------------C
      SUBROUTINE MunkHPFilter(U_i,V_i,TimeStep)
C---------------------------------------------------------------------C
C     This subroutine applies a high-pass filter derived from the
C     so-called Munk "Tide Killer" filter to the velocity field. The 
C     coefficients for the original, low pass filter can be found at: 
C
C        https://www.sonel.org/Filters-for-the-daily-mean-sea.html
C
C     The high-pass filter is derived from the normalized, nonrecursive
C     low-pass filter as:
C
C        W^{HP}_0 = 1-W^{LP}_0
C        W^{HP}_k = -W^{LP}_k   (k not equal to 0)
C
C     The filter is applied as:
C
C        y_n = \sum_{k=-m}^m W_k * x_{n+k}
C
C     Note that the output of the filter is lagged by 25 hours from the
C     current time step (which is fine since we are going to apply this
C     to the internal wave drag and 25 hours is approximately twice the
C     semi-diurnal period).
C
C     Written by: Coleman Blakely 3/2023
C---------------------------------------------------------------------C
      USE MESH, ONLY : NP
      USE GLOBAL, ONLY: DTDP
      IMPLICIT NONE
      REAL(8),INTENT(IN),DIMENSION(:) :: U_i, V_i ! current velocity
      INTEGER,INTENT(IN) :: TimeStep    ! current time step
      REAL(8) :: filtL = 49d0*3600d0        ! filter length (s)
      REAL(8) :: T = 60d0*60d0             ! sampling interval (s)
      INTEGER,SAVE :: NS = 1   ! number of samples
      LOGICAL,SAVE :: first_call = .TRUE.
      ! need indices to check if the current and previous time step are
      ! in the same 1 hour window
      INTEGER :: L, Lm           
      INTEGER :: ii, kk ! for loops
      ! for populating the UAV, VAV matrices at the start of the run
      INTEGER,SAVE :: ISTA

      IF ( first_call ) THEN
         first_call = .FALSE.
         NS = FLOOR(filtL/T)
         ista = 1
         ALLOCATE( UAV(NP,NS), VAV(NP,NS) )
         UAV = 0d0
         VAV = 0d0
         ALLOCATE( UBar(NP), VBar(NP) )
         UBar = 0d0
         VBar = 0d0
         CALL CalcMunkHPWeights()
      ENDIF
      ! check if this and the previous timestep are in the same hour
      L = FLOOR(TimeStep*DTDP/T)
      Lm = FLOOR((TimeStep-1)*DTDP/T)
      IF ( L.GT.Lm ) THEN
         ! check if we are >= 49 hours into the run
         IF ( ISTA.GT.NS ) THEN
            DO ii = 1,NP
               UAV(ii,1:NS-1) = UAV(ii,2:NS)
               VAV(ii,1:NS-1) = VAV(ii,2:NS)
               UAV(ii,NS) = U_i(ii)
               VAV(ii,NS) = V_i(ii)
            ENDDO
            UBar = 0d0
            VBar = 0d0
            DO ii = 1,NP
               DO kk = 1,NS
                  UBar(ii) = UBar(ii) + wts(kk)*UAV(ii,kk)
                  VBar(ii) = VBar(ii) + wts(kk)*VAV(ii,kk)
               ENDDO
            ENDDO
         ELSE
            ! if we do not have a 49 hr record yet just use a lagged
            ! average
            DO ii = 1,NP
               UAV(ii,ISTA) = U_i(ii)
               VAV(ii,ISTA) = V_i(ii)
               UBar(ii) = U_i(ii) - SUM(UAV(ii,1:ISTA))/DBLE(ISTA)
               VBar(ii) = V_i(ii) - SUM(VAV(ii,1:ISTA))/DBLE(ISTA)
            ENDDO
            ISTA = ISTA + 1
         ENDIF
      ENDIF
C---------------------------------------------------------------------C
      END SUBROUTINE MunkHPFilter
C---------------------------------------------------------------------C
C---------------------------------------------------------------------C
      SUBROUTINE MunkLPFilter(U_i,V_i,TimeStep)
C---------------------------------------------------------------------C
C     This subroutine applies a low-pass filter derived from the
C     so-called Munk "Tide Killer" filter to the velocity field. The 
C     coefficients for the original, low pass filter can be found at: 
C
C        https://www.sonel.org/Filters-for-the-daily-mean-sea.html
C
C
C     Written by: Coleman Blakely 3/2023
C---------------------------------------------------------------------C
      USE SIZES, ONLY : MYPROC
      USE MESH, ONLY : NP
      USE GLOBAL, ONLY: DTDP, np_g
      USE GL2LOC_MAPPING, ONLY : GATHERGLOBAL_REAL  
      IMPLICIT NONE
      REAL(8),INTENT(IN),DIMENSION(:) :: U_i, V_i ! current velocity
      INTEGER,INTENT(IN) :: TimeStep    ! current time step
      REAL(8) :: filtL = 49d0*3600d0        ! filter length (s)
      REAL(8) :: T = 60d0*60d0             ! sampling interval (s)
      INTEGER,SAVE :: NS = 1   ! number of samples
      LOGICAL,SAVE :: first_call = .TRUE.
      ! need indices to check if the current and previous time step are
      ! in the same 1 hour window
      INTEGER :: L, Lm           
      INTEGER :: ii, kk ! for loops
      ! for populating the UAV, VAV matrices at the start of the run
      INTEGER,SAVE :: ISTA

      IF ( first_call ) THEN
         first_call = .FALSE.
         NS = FLOOR(filtL/T)
         ista = 1
         ALLOCATE( UAV(NP,NS), VAV(NP,NS) )
         UAV = 0d0
         VAV = 0d0
         ALLOCATE( UBar(NP), VBar(NP) )
         UBar = 0d0
         VBar = 0d0
         CALL CalcMunkLPWeights()
      ENDIF
      ! check if this and the previous timestep are in the same hour
      L = FLOOR(TimeStep*DTDP/T)
      Lm = FLOOR((TimeStep-1)*DTDP/T)
      IF ( L.GT.Lm ) THEN
         ! check if we are >= 49 hours into the run
         IF ( ISTA.GT.NS ) THEN
            DO ii = 1,NP
               DO kk = 1,NS-1
                  UAV(ii,kk) = UAV(ii,kk+1)
                  VAV(ii,kk) = VAV(ii,kk+1)
               ENDDO
               UAV(ii,NS) = U_i(ii)
               VAV(ii,NS) = V_i(ii)
            ENDDO
            UBar = 0d0
            VBar = 0d0
            DO ii = 1,NP
               DO kk = 1,NS
                  UBar(ii) = UBar(ii) + wts(kk)*UAV(ii,kk)
                  VBar(ii) = VBar(ii) + wts(kk)*VAV(ii,kk)
               ENDDO
            ENDDO
         ELSE
            ! if we do not have a 49 hr record yet just use a lagged
            ! average
            DO ii = 1,NP
               UAV(ii,ISTA) = U_i(ii)
               VAV(ii,ISTA) = V_i(ii)
               UBar(ii) = SUM(UAV(ii,1:ISTA))/DBLE(ISTA)
               VBar(ii) = SUM(VAV(ii,1:ISTA))/DBLE(ISTA)
            ENDDO
            ISTA = ISTA + 1
         ENDIF
      ENDIF
C---------------------------------------------------------------------C
      END SUBROUTINE MunkLPFilter
C---------------------------------------------------------------------C
C---------------------------------------------------------------------C
      SUBROUTINE CalcMunkLPWeights()
C---------------------------------------------------------------------C
C     This subroutine sets up the filter weights for use in the
C     high-pass filter derived from the so-called Munk "Tide Killer"
C     low-pass filter. 
C---------------------------------------------------------------------C
      USE SIZES, ONLY : MYPROC
      IMPLICIT NONE
      REAL(8),DIMENSION(25) :: LPwts   ! original weights (one sided)
      INTEGER :: NS = 49                ! length of HP filter needed
      REAL(8) :: K                     ! for normalizing the LP filter
      INTEGER :: ii                     ! for loops
      !
      ! allocate weights
      ALLOCATE( wts(NS) )
      ! define low-pass filter weights (not normalized and one-sided)
      LPwts = (/395287,386839,370094,354118,338603,325633,314959,
     &          300054,278167,251492,234033,219260,208050,195518,
     &          180727,165525,146225,122665,101603,85349,72261,
     &          60772,47028,30073,13307 /)
      ! find sum to normalize low pass filter
      K = LPwts(1)
      DO ii = 2,25
         K = K + 2d0*LPwts(ii)
      END DO
      ! normalize low pass filter weights
      DO ii = 1,25
         LPwts(ii) = LPwts(ii)/K
      END DO
      ! turn one-sided low-pass filter into two-sided high-pass filter
      DO ii = 1,24
         ! -m to -1
         wts(ii) = LPwts(26-ii)
         ! 1 to m
         wts(ii+25) = LPwts(ii+1)
      END DO
      ! 0 (center point)
      wts(25) = LPwts(1)
C---------------------------------------------------------------------C
      END SUBROUTINE CalcMunkLPWeights
C---------------------------------------------------------------------C
C---------------------------------------------------------------------C
      SUBROUTINE CalcMunkHPWeights()
C---------------------------------------------------------------------C
C     This subroutine sets up the filter weights for use in the
C     high-pass filter derived from the so-called Munk "Tide Killer"
C     low-pass filter. 
C---------------------------------------------------------------------C
      IMPLICIT NONE
      REAL(8),DIMENSION(25) :: LPwts   ! original weights (one sided)
      INTEGER :: NS = 49                ! length of HP filter needed
      REAL(8) :: K                     ! for normalizing the LP filter
      INTEGER :: ii                     ! for loops
      !
      ! allocate weights
      ALLOCATE( wts(NS) )
      ! define low-pass filter weights (not normalized and one-sided)
      LPwts = (/395287,386839,370094,354118,338603,325633,314959,
     &          300054,278167,251492,234033,219260,208050,195518,
     &          180727,165525,146225,122665,101603,85349,72261,
     &          60772,47028,30073,13307 /)
      ! find sum to normalize low pass filter
      K = LPwts(1)
      DO ii = 2,25
         K = K + 2d0*LPwts(ii)
      END DO
      ! normalize low pass filter weights
      DO ii = 1,25
         LPwts(ii) = LPwts(ii)/K
      END DO
      ! turn one-sided low-pass filter into two-sided high-pass filter
      DO ii = 1,24
         ! -m to -1
         wts(ii) = -LPwts(26-ii)
         ! 1 to m
         wts(ii+25) = -LPwts(ii+1)
      END DO
      ! 0 (center point)
      wts(25) = 1 - LPwts(1)
C---------------------------------------------------------------------C
      END SUBROUTINE CalcMunkHPWeights
C---------------------------------------------------------------------C
!-----------------------------------------------------------------------
!     S U B R O U T I N E   I T F R I C T E R M I N A T E
!-----------------------------------------------------------------------
      SUBROUTINE ITFRICTERMINATE(NO_MPI_FINALIZE)
      USE SIZES, ONLY : myproc
#ifdef CMPI 
      USE MESSENGER
#endif
      IMPLICIT NONE 
      LOGICAL, OPTIONAL :: NO_MPI_FINALIZE
#ifdef ALL_TRACE
      REAL, ALLOCATABLE :: dummy(:)
#endif
C
      !call setMessageSource("BC3DTerminate")
#ifdef ALL_TRACE
      call allMessage(DEBUG,"Enter.")
#endif

      call allMessage(INFO,"ADCIRC Terminating.")

#ifdef ALL_TRACE
      ! intentionally create a segmentation fault so that we can get
      ! a stack trace to determine the line number of the netcdf call
      ! that went bad ... this assumes that the code was compiled with
      ! debugging symbols, bounds checking, and stack trace turned on.
      write(6,*) 'myproc = ',myProc      
      dummy(1) = 99.9d0
#endif

#ifdef CMPI
      subdomainFatalError = .true.
      IF (PRESENT(NO_MPI_FINALIZE)) THEN
        CALL MSG_FINI(NO_MPI_FINALIZE)
      ELSE
        CALL MSG_FINI()
      ENDIF
#endif
      CALL EXIT(1)
C
#ifdef ALL_TRACE
      call allMessage(DEBUG,"Return.") ! should be unreachable
#endif
      !call unsetMessageSource()
C---------------------------------------------------------------------C
      END SUBROUTINE ITFRICTERMINATE
C---------------------------------------------------------------------C
C-----------------------------------------------------------------------
C     S U B R O U T I N E   C H E C K  _  E R R
C-----------------------------------------------------------------------
C     jgf49.17.02 Checks the return value from netCDF calls; if there
C     was an error, it writes the error message to the screen and to the
C     fort.16 file.
C-----------------------------------------------------------------------
      subroutine check_err(iret)
      USE SIZES, ONLY : myproc
      USE GLOBAL, ONLY : screenUnit
#ifdef CMPI
      USE MESSENGER, ONLY : MSG_FINI
#endif
      USE NETCDF
      IMPLICIT NONE

      INTEGER, intent(in) :: iret

      call setMessageSource("check_err")
#ifdef ALL_TRACE
      call allMessage(DEBUG,"Enter.")
#endif
      if (iret .ne. NF90_NOERR) then
         call allMessage(ERROR,nf90_strerror(iret))
         call ITFRICTerminate()
      endif
#ifdef ALL_TRACE
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      end subroutine check_err
C-----------------------------------------------------------------------
C---------------------------------------------------------------------C
      END MODULE INTERNALTIDE
C---------------------------------------------------------------------C
