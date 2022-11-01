C---------------------------------------------------------------------C
      MODULE INTERNALTIDE
C---------------------------------------------------------------------C
C     CPB 10/2022: This module adds the capability to read in a harmonic
C     decomposition of tidal harmonics and resynthesize them to be used
C     in the internal tide dissipation/generation term. When used, it
C     turns the Cit tensor into a forcing term on the RHS of the
C     momentum equation rather than lumping it with the boundary layer
C     dissipation in the TKM matrix. This functionality is controlled by
C     the itFricControl namelist
C---------------------------------------------------------------------C
      USE SIZES, ONLY: SZ
      USE GLOBAL, ONLY : DEBUG, ECHO, INFO, WARNING, ERROR,
     &          setMessageSource, unsetMessageSource, allMessage,
     &          scratchMessage, logMessage      
      IMPLICIT NONE
      LOGICAL :: resynthTidalVelo
      CHARACTER(LEN=100) :: tidalVeloFilename
      REAL(SZ),ALLOCATABLE,PRIVATE :: UTamp(:,:), UTphs(:,:),
     &                                VTamp(:,:), VTphs(:,:),
     &                                UTsin(:,:), UTcos(:,:),
     &                                VTsin(:,:), VTcos(:,:)
C---------------------------------------------------------------------C
      CONTAINS
C---------------------------------------------------------------------C
      SUBROUTINE READ_TIDAL_VELO()
C---------------------------------------------------------------------C
C     This subroutine reads in the tidalvelo.nc file which contains the
C     harmonic decomposition of the zonal and meridional tidal
C     velocities. This file has the same format as the fort.54.nc output
C     file.              
C---------------------------------------------------------------------C
      USE NETCDF
      USE SIZES, ONLY : myproc
      USE GLOBAL, ONLY : deg2rad, TIPOTAG, nodes_lg, np_g, NTIF
      USE MESH, ONLY : NP
#IFDEF CMPI
      USE MESSENGER, ONLY : MapToSubdomainRealMPI
#ENDIF      
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
      CALL ALLOCATE_TIDAL_DECOMP()
      IF (myproc.eq.0) THEN
         !
         ! Open fort.112.nc
         !
         CALL check_err(nf90_open(tidalVeloFileName, nf90_nowrite, 
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
     &                   //' specified in fort.15. Check inputs.'
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
      END SUBROUTINE READ_TIDAL_VELO
C---------------------------------------------------------------------C
      SUBROUTINE CALC_TIDAL_VELO(Utidal, Vtidal, timeh)
C---------------------------------------------------------------------C
      USE GLOBAL, ONLY : DTDP, IFSPROTS, AMIGT, PERT, NTIF
      USE MESH, ONLY : NP, DRVMAP2DSPVEC, UVECTMP, VVECTMP
      IMPLICIT NONE
      LOGICAL, SAVE :: first_call = .TRUE.
      REAL(SZ), INTENT(OUT), DIMENSION(:) :: Utidal, Vtidal
      REAL(SZ), INTENT(IN) :: timeh
      INTEGER :: I, J, NCYC
      REAL(SZ) :: cosAMIGT, sinAMIGT
      IF (first_call) THEN
         first_call = .FALSE.
         CALL READ_TIDAL_VELO()
         CALL CALC_SINCOS_TERMS()
      ENDIF
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
      END SUBROUTINE CALC_TIDAL_VELO
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
      END SUBROUTINE CALC_SINCOS_TERMS
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
      WRITE(16,*) 'allocated_tidal_decomp'
C      
      RETURN
      END SUBROUTINE ALLOCATE_TIDAL_DECOMP        
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
      END MODULE INTERNALTIDE
C---------------------------------------------------------------------C
