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
      USE SIZES, ONLY: SZ
      USE GLOBAL, ONLY : DEBUG, ECHO, INFO, WARNING, ERROR,
     &          setMessageSource, unsetMessageSource, allMessage,
     &          scratchMessage, logMessage      
      IMPLICIT NONE
      ! namelist variables
      CHARACTER(LEN=100) :: filterType, tideVeloFilename
      ! stored velocity samples for moving average calculation
      REAL(SZ),ALLOCATABLE,PRIVATE :: UAV(:,:), VAV(:,:)        
      ! averaged velocities so we don't calculate them every single time
      ! step
      REAL(SZ),ALLOCATABLE :: UBar(:), VBar(:)
      ! tidal velocities
      REAL(SZ),ALLOCATABLE :: UTidal(:), VTidal(:)
      ! weights for filters
      REAL(SZ),ALLOCATABLE,PRIVATE :: wts(:)
      REAL(SZ) :: filtL        ! filter length (s)
      REAL(SZ) :: Fs             ! sampling interval (s)
      INTEGER :: NS   ! number of samples
      ! variables to resynthesize tidal velocities
      REAL(SZ),ALLOCATABLE,PRIVATE :: UTamp(:,:), UTphs(:,:),
     &                                VTamp(:,:), VTphs(:,:),
     &                                UTsin(:,:), UTcos(:,:),
     &                                VTsin(:,:), VTcos(:,:)
C---------------------------------------------------------------------C
      CONTAINS
C---------------------------------------------------------------------C
C---------------------------------------------------------------------C
      SUBROUTINE CalcTidalVelo(U_i,V_i,TimeStep)
C---------------------------------------------------------------------C
C     This subroutine initializes the filters for estimating the tidal
C     velocities and also calls the appropriate subroutine to estimate
C     the tidal velocity based on the internalTide namelist
C---------------------------------------------------------------------C
      USE MESH, ONLY : NP
      USE GLOBAL, ONLY: DTDP
      IMPLICIT NONE
      REAL(SZ),INTENT(IN),DIMENSION(:) :: U_i, V_i ! current velocity
      INTEGER,INTENT(IN) :: TimeStep    ! current time step
      LOGICAL,SAVE :: first_call = .TRUE.

      IF ( first_call ) THEN
         first_call = .FALSE.
         WRITE(16,*) "Inside CalcTidalVelo"
         CALL InitCalcTidalVelo()
      ENDIF
      UTidal = 0d0
      VTidal = 0d0
      SELECT CASE (TRIM(filterType))
         CASE ("MunkHP")               ! NOTE: Default
            CALL MunkHPFilter(U_i,V_i,TimeStep)
         CASE ("LA25")
            CALL UNTIDE_LA25(U_i,V_i,TimeStep)
         CASE ("resynthTide")
            CALL ResynthTideVelo(TimeStep)
      END SELECT


C-----------------------------------------------------------------------
      END SUBROUTINE CalcTidalVelo
C-----------------------------------------------------------------------
C---------------------------------------------------------------------C
      SUBROUTINE InitCalcTidalVelo
C---------------------------------------------------------------------C
C     This subroutine initializes the required variables for tidal
C     filtering based on the namelist variable "filterType". Possible
C     filter types are:
C        1) MunkHP (default): a 49 hour high pass filter which estimates
C           the tidal velocity. It is derived from the so-called Munk
C           "Tide Killer" filter (Groves, 1955). 
C        2) LA25: a 25 hour lagged average filter which estimates the
C           mean flow and subtracts that off of the total velocity to
C           estimate the tidal velocity
C        3) resynthTide: causes ADCIRC to read in a fort.54.nc (global
C           harmonic analysis of velocity) and resynthesize tidal
C           velocities from the constituents used in the fort.15 to
C           force the equilibrium tide.
C---------------------------------------------------------------------C
      USE MESH, ONLY : NP
      IMPLICIT NONE
      
      SELECT CASE (TRIM(filterType))
         CASE ("MunkHP")               ! NOTE: Default
            filtL = 49d0*3600d0        ! filter length (s)
            Fs = 60d0*60d0             ! sampling interval (s)
            NS = FLOOR(filtL/Fs)       ! number of samples
            ALLOCATE( UAV(NP,NS), VAV(NP,NS) )
            UAV = 0d0
            VAV = 0d0
            ALLOCATE( UBar(NP), VBar(NP) )
            UBar = 0d0
            VBar = 0d0
            ALLOCATE( UTidal(NP), VTidal(NP) )
            UTidal = 0d0
            VTidal = 0d0
            CALL CalcMunkWeights()
         CASE ("LA25")
            filtL = 25d0*3600d0        ! filter length (s)
            Fs = 60d0*60d0             ! sampling interval (s)
            NS = FLOOR(filtL/Fs)       ! number of samples
            ALLOCATE( UAV(NP,NS), VAV(NP,NS) )
            UAV = 0d0
            VAV = 0d0
            ALLOCATE( UBar(NP), VBar(NP) )
            UBar = 0d0
            VBar = 0d0
            ALLOCATE( UTidal(NP), VTidal(NP) )
            UTidal = 0d0
            VTidal = 0d0
         CASE ("resynthTide")
            ALLOCATE( UTidal(NP), VTidal(NP) )
            UTidal = 0d0
            VTidal = 0d0
#ifdef ADCNETCDF
            CALL ReadTideVeloFile()
            CALL CALC_SINCOS_TERMS()
            WRITE(16,*) "INFO: Successfully read tidal velocity file "
     &                  //"and calculated sin/cos terms."
#else
            WRITE(16,*) 'ERROR: NetCDF is required to resynthesize'
     &                   //' tidal velocities. Recompile with NetCDF'
     &                   //' and try again.'
            CALL ITFRICTERMINATE()
#endif
         CASE DEFAULT
            WRITE(16,*) 'ERROR: filterType "',filterType,'" not'
     &                  //' recognized.'
            CALL ITFRICTERMINATE()
      END SELECT
C-----------------------------------------------------------------------
      END SUBROUTINE InitCalcTidalVelo
C-----------------------------------------------------------------------
      SUBROUTINE UNTIDE_LA25(U_i,V_i,TimeStep)
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
      ! need indices to check if the current and previous time step are
      ! in the same hour window
      INTEGER :: L, Lm           
      INTEGER :: ii, kk ! for loops
      ! for populating the UAV, VAV matrices at the start of the run
      INTEGER,SAVE :: ISTA = 1

      ! check if this and the previous timestep are in the same hour
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
      DO ii = 1,NP
         UTidal(ii) = U_i(ii) - UBar(ii)
         VTidal(ii) = V_i(ii) - VBar(ii)
      ENDDO
C-----------------------------------------------------------------------
      END SUBROUTINE UNTIDE_LA25
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
      REAL(SZ),INTENT(IN),DIMENSION(:) :: U_i, V_i ! current velocity
      INTEGER,INTENT(IN) :: TimeStep    ! current time step
      ! need indices to check if the current and previous time step are
      ! in the same 1 hour window
      INTEGER :: L, Lm           
      INTEGER :: ii, kk ! for loops
      ! for populating the UAV, VAV matrices at the start of the run
      INTEGER,SAVE :: ISTA = 1

      ! check if this and the previous timestep are in the same hour
      L = FLOOR(TimeStep*DTDP/Fs)
      Lm = FLOOR((TimeStep-1)*DTDP/Fs)
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
               UTidal(ii) = UBar(ii)
               VTidal(ii) = VBar(ii)
            ENDDO
         ELSE
            ! if we do not have a 49 hr record yet just use a lagged
            ! average
            DO ii = 1,NP
               UAV(ii,ISTA) = U_i(ii)
               VAV(ii,ISTA) = V_i(ii)
               UTidal(ii) = U_i(ii) - SUM(UAV(ii,1:ISTA))/DBLE(ISTA)
               VTidal(ii) = V_i(ii) - SUM(VAV(ii,1:ISTA))/DBLE(ISTA)
            ENDDO
            ISTA = ISTA + 1
         ENDIF
      ENDIF
C---------------------------------------------------------------------C
      END SUBROUTINE MunkHPFilter
C---------------------------------------------------------------------C
C---------------------------------------------------------------------C
      SUBROUTINE CalcMunkWeights()
C---------------------------------------------------------------------C
C     This subroutine sets up the filter weights for use in the
C     high-pass filter derived from the so-called Munk "Tide Killer"
C     low-pass filter. 
C---------------------------------------------------------------------C
      IMPLICIT NONE
      REAL(SZ),DIMENSION(25) :: LPwts   ! original weights (one sided)
      INTEGER :: NS = 49                ! length of HP filter needed
      REAL(SZ) :: K                     ! for normalizing the LP filter
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
      END SUBROUTINE CalcMunkWeights
C---------------------------------------------------------------------C
C---------------------------------------------------------------------C
      SUBROUTINE ReadTideVeloFile()
C---------------------------------------------------------------------C
C     This subroutine reads in the tidalvelo.nc file which contains the
C     harmonic decomposition of the zonal and meridional tidal
C     velocities. This file has the same format as the fort.54.nc output
C     file.              
C---------------------------------------------------------------------C
#ifdef ADCNETCDF
      USE NETCDF
#endif
      USE SIZES, ONLY : myproc
      USE GLOBAL, ONLY : deg2rad, TIPOTAG, nodes_lg, np_g, NTIF
      USE MESH, ONLY : NP
#ifdef CMPI
      USE MESSENGER, ONLY : MapToSubdomainRealMPI
#endif      
      IMPLICIT NONE
      ! id for ncfile and temp_id for reading in variables
      INTEGER :: nc_id, temp_id, uamp_id, uphs_id, vamp_id, vphs_id
      ! logical for if we found the right constituent
      LOGICAL :: isfound
      ! looping integers
      INTEGER :: ii, kk
      ! constituent names from fort.112.nc
      CHARACTER(10),ALLOCATABLE :: constnames(:)
      ! number of constituents in fort.112.nc
      INTEGER :: num_const
      ! arrays to store values read in from netcdf
      REAL(SZ),ALLOCATABLE,DIMENSION(:) :: uamp_g, uphs_g, 
     &                                     vamp_g, vphs_g,
     &                                     tmp      
#ifdef ADCNETCDF
      CALL ALLOCATE_TIDAL_DECOMP()
      IF (myproc.eq.0) THEN
         !
         ! Open fort.112.nc
         !
         CALL check_err(nf90_open(tideVeloFileName, nf90_nowrite, 
     &                  nc_id))
         !
         ! Read in constituent names from fort.112.nc
         !
         CALL check_err(nf90_inq_dimid(nc_id, 'num_const', temp_id))
         CALL check_err(nf90_inquire_dimension(nc_id, temp_id, 
     &                  len=num_const))
         ALLOCATE( constnames(num_const) )
         CALL check_err(nf90_inq_varid(nc_id, 'const', temp_id))
         CALL check_err(nf90_get_var(nc_id, temp_id, constnames))
         !
         ! Allocate temporary arrays
         !
         ALLOCATE( uamp_g(np_g), uphs_g(np_g), 
     &             vamp_g(np_g), vphs_g(np_g) )
         !
         ! get id's of variables
         !
         CALL check_err(nf90_inq_varid(nc_id, 'u_amp', uamp_id))
         CALL check_err(nf90_inq_varid(nc_id, 'u_phs', uphs_id))
         CALL check_err(nf90_inq_varid(nc_id, 'v_amp', vamp_id))
         CALL check_err(nf90_inq_varid(nc_id, 'v_phs', vphs_id))
         WRITE(16,*) 'read in global stuff'
      ENDIF
      ALLOCATE( tmp(NP) )
      ! 
      ! read and distribute data
      !
      DO ii = 1, NTIF
         IF (myproc.eq.0) THEN
            ! identify which constituent we want
            isfound = .false.
            DO kk = 1,num_const
               IF ( TRIM(tipotag(ii)).EQ.TRIM(constnames(kk)) ) THEN
                  isfound = .true.
                  exit
               ENDIF
            ENDDO
            IF (.NOT.isfound) THEN
               WRITE(16,*) 'ERROR: fort.112.nc is missing a constituent'
     &                   //' specified in fort.15. check inputs.'
               CALL ITFRICTERMINATE()
            ENDIF
            !
            ! read in this amp and phase
            !
            call check_err(nf90_get_var(nc_id, uamp_id, uamp_g,
     &                     start=[kk,1],count=[1,np_g]))
            call check_err(nf90_get_var(nc_id, uphs_id, uphs_g,
     &                     start=[kk,1],count=[1,np_g]))
            call check_err(nf90_get_var(nc_id, vamp_id, vamp_g,
     &                     start=[kk,1],count=[1,np_g]))
            call check_err(nf90_get_var(nc_id, vphs_id, vphs_g,
     &                     start=[kk,1],count=[1,np_g]))
         ENDIF
         !
         ! distribute
         !
         ! uamp
#ifdef CMPI
         CALL MapToSubdomainRealMPI(np_g,np,tmp,nodes_lg,uamp_g)
         UTAMP(ii,:) = tmp
         ! uphs
         CALL MapToSubdomainRealMPI(np_g,np,tmp,nodes_lg,uphs_g)
         UTPHS(ii,:) = tmp*deg2rad
         ! vamp
         CALL MapToSubdomainRealMPI(np_g,np,tmp,nodes_lg,vamp_g)
         VTAMP(ii,:) = tmp
         ! vphs
         CALL MapToSubdomainRealMPI(np_g,np,tmp,nodes_lg,vphs_g)
         VTPHS(ii,:) = tmp*deg2rad
#endif         
      ENDDO
      IF (myproc.eq.0) THEN
         DEALLOCATE( uamp_g, uphs_g, vamp_g, vphs_g )
      ENDIF
      DEALLOCATE( tmp )      
#else
      WRITE(16,*) 'ERROR: NetCDF is required to resynthesize'
     &                   //' tidal velocities. Recompile with NetCDF'
     &                   //' and try again.'
      CALL ITFRICTERMINATE()
#endif
C---------------------------------------------------------------------C
      END SUBROUTINE ReadTideVeloFile
C---------------------------------------------------------------------C
C---------------------------------------------------------------------C
      SUBROUTINE ALLOCATE_TIDAL_DECOMP()
C---------------------------------------------------------------------C
C     Simple subroutine to allocate amplitude and phases of zonal and
C     meridional tidal velocities. This subroutine is only called in
C     READ_TIDAL_VELO and the amplitudes and phases are deallocated as
C     soon as they can be.
C---------------------------------------------------------------------C
      USE GLOBAL, ONLY : NTIF
      USE MESH, ONLY : NP
      IMPLICIT NONE
C
      ALLOCATE ( UTamp(NTIF,NP), UTphs(NTIF,NP),
     &           VTamp(NTIF,NP), VTphs(NTIF,NP),
     &           UTsin(NTIF,NP), UTcos(NTIF,NP),
     &           VTsin(NTIF,NP), VTcos(NTIF,NP) )
C
      RETURN
      END SUBROUTINE ALLOCATE_TIDAL_DECOMP
C---------------------------------------------------------------------C
C---------------------------------------------------------------------C
      SUBROUTINE CALC_SINCOS_TERMS()
C---------------------------------------------------------------------C
C     Using trigonometric identities, we can simplify the number of
C     sin/cos functions need to be called to resynthesize tidal
C     velocities. This subroutine takes the tidal decomposition read in
C     from the tidalvelo.nc file and calculates the time-invariate term
C     needed to build our tidal velocities  
C---------------------------------------------------------------------C
      USE GLOBAL, ONLY : FACET, FFT, NTIF
      USE MESH, ONLY : NP
      IMPLICIT NONE
      INTEGER :: I, J
      DO I = 1,NTIF
         DO J = 1,NP
            UTsin(I,J) = SIN(FACET(I)-UTphs(I,J))
            UTcos(I,J) = COS(FACET(I)-UTphs(I,J))
            VTsin(I,J) = SIN(FACET(I)-VTphs(I,J))
            VTcos(I,J) = COS(FACET(I)-VTphs(I,J))
            ! multiply by nodal factors and amplitudes since
            ! these also do not change in time.
            UTsin(I,J) = UTsin(I,J)*FFT(I)*UTamp(I,J)
            UTcos(I,J) = UTcos(I,J)*FFT(I)*UTamp(I,J)
            VTsin(I,J) = VTsin(I,J)*FFT(I)*VTamp(I,J)
            VTcos(I,J) = VTcos(I,J)*FFT(I)*VTamp(I,J)
         ENDDO
      ENDDO
      DEALLOCATE( UTamp, UTphs, VTamp, VTphs )      
C---------------------------------------------------------------------C
      END SUBROUTINE CALC_SINCOS_TERMS
C---------------------------------------------------------------------C
C---------------------------------------------------------------------C
      SUBROUTINE ResynthTideVelo(TimeStep)
C---------------------------------------------------------------------C
      USE GLOBAL, ONLY : DTDP, IFSPROTS, AMIGT, PERT, NTIF, StaTim,
     &                   RefTim, DTDPHS, ITHS, CHOTHS
      USE MESH, ONLY : NP, DRVMAP2DSPVEC, UVECTMP, VVECTMP
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: TimeStep
      REAL(SZ) :: timeh, StaTimHS, StaTimTMP
      INTEGER :: I, J, NCYC
      REAL(SZ) :: cosAMIGT, sinAMIGT
      ! calc harmonic time
      TimeH=TimeStep*DTDP + (StaTim - RefTim)*86400.D0
      IF (CHOTHS.eqv..true.) THEN
         IF ((ITHS+1).EQ.TimeStep) THEN
           StaTimHS=((TimeStep-1)*DTDPHS)/86400.D0
           StaTimTMP=((TimeStep-1)*DTDP)/86400.D0
         END IF
         TimeH=TimeStep*DTDP + ((StaTimHS - StaTimTMP) - RefTim)*86400.D0
      END IF

      UTidal = 0d0
      VTidal = 0d0

      DO I = 1,NTIF
         IF (PERT(I).EQ.0) THEN
            NCYC = 0
         ELSE
#ifdef IBM
            NCYC = INT(timeh/PERT(I),KIND(0d0))
#else
            NCYC = INT(timeh/PERT(I))
#endif
         ENDIF
         cosAMIGT = COS(AMIGT(I)*(timeh-NCYC*PERT(I)))
         SINAMIGT = SIN(AMIGT(I)*(timeh-NCYC*PERT(I)))
         DO J = 1,NP
            UTidal(J) = UTidal(J) + cosAMIGT*UTcos(I,J) 
     &                            - sinAMIGT*UTsin(I,J)
            VTidal(J) = VTidal(J) + cosAMIGT*VTcos(I,J)
     &                            - sinAMIGT*VTsin(I,J)
         ENDDO
      ENDDO
      IF (IFSPROTS.EQ.1) THEN
         UVECTMP(1:NP) = UTidal(1:NP)
         VVECTMP(1:NP) = VTidal(1:NP)
         CALL DRVMAP2DSPVEC(UTidal,VTidal,UVECTMP,VVECTMP,NP,
     &        FWD=.TRUE.)           
      ENDIF
C-----------------------------------------------------------------------
      END SUBROUTINE ResynthTideVelo
C-----------------------------------------------------------------------
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
#ifdef ADCNETCDF
      USE NETCDF
#endif
      IMPLICIT NONE
      INTEGER, intent(in) :: iret
      !call setMessageSource("check_err")
#ifdef ALL_TRACE
      call allMessage(DEBUG,"Enter.")
#endif
#ifdef ADCNETCDF
      if (iret .ne. NF90_NOERR) then
         call allMessage(ERROR,nf90_strerror(iret))
         call ITFRICTERMINATE()
      endif
#endif
#ifdef ALL_TRACE
      call allMessage(DEBUG,"Return.")
#endif
      !call unsetMessageSource()
C-----------------------------------------------------------------------
      end subroutine check_err
C-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     S U B R O U T I N E   B C 3 D  T E R M I N A T E
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
      call allMessage(INFO,"ADCIRC Terminating.")
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
!-----------------------------------------------------------------------
      END SUBROUTINE ITFRICTERMINATE
!-----------------------------------------------------------------------  
C---------------------------------------------------------------------C
C---------------------------------------------------------------------C
      END MODULE INTERNALTIDE
C---------------------------------------------------------------------C
