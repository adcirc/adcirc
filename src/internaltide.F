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
      USE mod_logging, ONLY : DEBUG, ECHO, INFO, WARNING, ERROR,
     &          setMessageSource, unsetMessageSource, allMessage,
     &          scratchMessage, logMessage      
      IMPLICIT NONE
      ! stored velocity samples for moving average calculation
      REAL(8),ALLOCATABLE,PRIVATE :: UAV(:,:), VAV(:,:)        
      ! averaged velocities so we don't calculate them every single time
      ! step
      REAL(8),ALLOCATABLE :: UBar(:), VBar(:)
      ! weights for filters
      REAL(8),ALLOCATABLE,PRIVATE :: wts(:)
C---------------------------------------------------------------------C
      CONTAINS
C---------------------------------------------------------------------C
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
      REAL(8),INTENT(IN),DIMENSION(:) :: U_i, V_i ! current velocity
      INTEGER,INTENT(IN) :: TimeStep    ! current time step
      REAL(8) :: filtL = 25d0*3600d0        ! filter length (s)
      REAL(8) :: Fs = 12d0*60d0             ! sampling interval (s)
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
         CALL CalcMunkWeights()
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
      SUBROUTINE CalcMunkWeights()
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
      END SUBROUTINE CalcMunkWeights
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
      USE mod_logging, ONLY : screenUnit
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
C---------------------------------------------------------------------C
      END MODULE INTERNALTIDE
C---------------------------------------------------------------------C
