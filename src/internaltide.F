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
      ! stored velocity samples for moving average calculation
      REAL(SZ),ALLOCATABLE,PRIVATE :: UAV(:,:), VAV(:,:)        
      ! averaged velocities so we don't calculate them every single time
      ! step
      REAL(SZ),ALLOCATABLE :: UBar(:), VBar(:)
C---------------------------------------------------------------------C
      CONTAINS
C---------------------------------------------------------------------C
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
      REAL(SZ) :: Fs = 12d0*60d0             ! sampling frequency (s)
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

      END SUBROUTINE UNTIDE
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

      !call setMessageSource("check_err")
#ifdef ALL_TRACE
      call allMessage(DEBUG,"Enter.")
#endif
      if (iret .ne. NF90_NOERR) then
         call allMessage(ERROR,nf90_strerror(iret))
         call ITFRICTERMINATE()
      endif
#ifdef ALL_TRACE
      call allMessage(DEBUG,"Return.")
#endif
      !call unsetMessageSource()
C-----------------------------------------------------------------------
      end subroutine check_err
C-----------------------------------------------------------------------
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
!-----------------------------------------------------------------------
      END SUBROUTINE ITFRICTERMINATE
!-----------------------------------------------------------------------      
C---------------------------------------------------------------------C
      END MODULE INTERNALTIDE
C---------------------------------------------------------------------C
