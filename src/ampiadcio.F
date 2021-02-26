      !
      !   1. Closing relevant input files before conducting migration
      !   2. After mirgration, reopening and fast forward input files closed in in step 1 
      !      to the location where they are closed
      !   
      !  DW: July 2020
      !      
      MODULE AMPIIOSUBS
      use CkLunMod, only : CK_LUN
      USE SIZES 
      USE GLOBAL, ONLY: NWS, IHOT, ITHS, ITIME, STATIM, DTDP, DTDPHS,
     &   NDDT, NBFR, NFFR, NWS, NRS, NCICE, RSTIMINC, CICE_TIMINC,
     &   CICE_TIMINC, CloseLogFile, OpenLogFile, PRBCKGRND,
     &   setMessageSource, unsetMessageSource, reopenlogfile

      USE BOUNDARIES, ONLY: NOPE, NFLUXF, NVEL, NETA, LBCODEI 

      USE OWIWIND, ONLY: owiwind_closefiles, NWS12GET, NWS12INIT
      USE MESSENGER, ONLY: MPI_COMM_ADCIRC

      USE GLOBAL_IO, ONLY: NDDT2GET
 
      USE SUBDOMAIN, ONLY: subdomainon, enforcebn
     
      USE OWI_ICE, only: NumSetsIce => NumSets, NCICE1_INIT, NCICE1_GET

      USE RS2, only: rs2init, rs2get
 
#ifdef HAVE_MPI_MOD
      USE MPI
#endif
 
      IMPLICIT NONE

      
      PRIVATE coldstartMeteorologicalReopenInputFile,  
     &     hotstartMeteorologicalReopenInputfile

      CONTAINS
c 
C----------------------------------------------------------------------
C   S U B R O U T I N E    A B O U T _ T O _ M I G R A T E
C----------------------------------------------------------------------
C  Routine called by AMPI before migrating, should close all open files.
C----------------------------------------------------------------------
        SUBROUTINE ABOUT_TO_MIGRATE()
          IMPLICIT NONE

#ifdef AMPI_COMP
          include 'mpif.h'
#endif
c          INTEGER :: RANK = 999999
c!$omp threadprivate(RANK)
c          INTEGER :: IERR, FLAG, OLD_PE

#ifdef AMPI_COMP
c          CALL MPI_COMM_RANK(MPI_COMM_ADCIRC, RANK, IERR)
c          CALL MPI_COMM_GET_ATTR(MPI_COMM_WORLD, AMPI_MY_WTH, OLD_PE, FLAG, IERR)

C          PRINT *,"Rank ",RANK," about to migrate from PE ",OLD_PE

C     TODO: close files temporarily during migration
C          CALL closeLogFile() 

          ! close a metfile
          CALL close_met_files() 

          ! close time dependent elevation and fluxes (if used)
          CALL close_elevbc_riverbc_files( ) 

          ! close time varying bathymetry
          CALL  close_time_varying_bathymetry_file( )

          ! close ice file
          CALL close_ice_files( )

          ! close wave radiation stress file
          CALL close_rs_file()
#endif
          RETURN
        END SUBROUTINE ABOUT_TO_MIGRATE

C----------------------------------------------------------------------
C   S U B R O U T I N E    J U S T _ M I G R A T E D
C----------------------------------------------------------------------
C  Routine called by AMPI after migrating, should re-open all closed files.
C----------------------------------------------------------------------
        SUBROUTINE JUST_MIGRATED()
          IMPLICIT NONE

#ifdef AMPI_COMP
          include 'mpif.h'
#endif
c          INTEGER :: RANK = 999999
c!$omp threadprivate(RANK)
c          INTEGER :: IERR, FLAG, NEW_PE

#ifdef AMPI_COMP
c          CALL MPI_COMM_RANK(MPI_COMM_ADCIRC, RANK, IERR)
c          CALL MPI_COMM_GET_ATTR(MPI_COMM_WORLD, AMPI_MY_WTH, NEW_PE, FLAG, IERR)

C          PRINT *,"Rank ",RANK," just migrated to PE ",NEW_PE

C     TODO: re-open files closed temporarily during migration
C          CALL openLogFile() 
C          CALL reopenLogFile() 


C          PRINT*, "reopen met file "   
          ! open and fast forward to location where the metfile
          ! was closed 
          CALL reopenffwd_met_files( ) 

          ! reopen and fast forward elev & flux input file
          CALL Reopen_elevbc_riverbc_InputFile(  )

          ! reopen and fast forward time varying bathymetry
          CALL Reopen_time_varying_bathymetry( )

          ! reopen and fast forward ice file
          CALL reopenffwd_ice_inputfile( )


          ! reopen and fast forward radiation stress file
          CALL reopenffwd_rs_file( )
   
!     will need to reopen input file and return to the point where
!     it is closed.
      
#endif
          RETURN 
        END SUBROUTINE JUST_MIGRATED



C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        SUBROUTINE ABOUT_TO_MIGRATE_TEST()
          IMPLICIT NONE

          INTEGER :: RANK, IERR, FLAG, OLD_PE

          PRINT *,"Rank ",RANK," about to migrate from PE ",OLD_PE

C     TODO: close files temporarily during migration
          CALL closeLogFile() 

          ! close a metfile
          CALL close_met_files()  

          ! close time dependent elevation and fluxes (if used)
          CALL close_elevbc_riverbc_files( ) 

          ! close time varying bathymetry
          CALL  close_time_varying_bathymetry_file( )

          ! close ice file
          CALL close_ice_files( )

          ! close wave radiation stress file
          CALL close_rs_file()

          RETURN
        END SUBROUTINE ABOUT_TO_MIGRATE_TEST

        SUBROUTINE JUST_MIGRATED_TEST()
          IMPLICIT NONE

          INTEGER :: RANK, IERR, FLAG, NEW_PE

          PRINT *,"Rank ",RANK," just migrated to PE ",NEW_PE

C     TODO: re-open files closed temporarily during migration
          CALL openLogFile()
 
          ! open and fast forward to location where the metfile
          ! was closed 

C          print*, "reopen met file .... "
          CALL reopenffwd_met_files( ) 

          ! reopen and fast forward elev & flux input file
          CALL Reopen_elevbc_riverbc_InputFile(  )

          ! reopen and fast forward time varying bathymetry
          CALL Reopen_time_varying_bathymetry( )

          ! reopen and fast forward ice file
          CALL reopenffwd_ice_inputfile( )


          ! reopen and fast forward radiation stress file
          CALL reopenffwd_rs_file( )
   
!     will need to reopen input file and return to the point where
!     it is closed.
      
          RETURN 
        END SUBROUTINE JUST_MIGRATED_TEST
 
 
C--------------------------------------------------------------------------
C
C CLOSE/REOPEN & FAST FORWARD IF NECESSARY
C
C--------------------------------------------------------------------------

        ! 
        ! Close alll wind input files
        !
        SUBROUTINE close_met_files( )
         IMPLICIT NONE
 
C------------------------MET FORCING---------------------------------------
C
C...  UPDATE THE WIND STRESS AND SURFACE PRESSURE AND READ IN NEW VALUES
C...  FROM UNIT 22.  APPLY Ramp FUNCTION.
C

C        wind, radiation stress or atmospheric pressure forcings are used at all grid nodes and every time step.
        IF( NWS .EQ .1) THEN
           CALL CLOSEFORT22( ) 
        END IF

C     Wind stress and atmospheric pressure are read in at all grid nodes
C     at a time interval that does not equal the model time
C     step. Interpolation in time is used to synchronize the wind and
C     pressure information with the model time step.
         IF(ABS(NWS).EQ.2) THEN
           CALL CLOSEFORT22( )  
         ENDIF

C     Wind velocity in US Navy Fleet Numeric format interpolated in
C     space onto the ADCIRC grid and in time to synchronize the wind and
C     pressure information with the model time step. Garratt's formula
C     is used to compute wind stress from the wind velocity.
         IF(NWS.EQ.3) THEN
           CALL CLOSEFORT22( ) 
         ENDIF

C     Wind velocity and atmospheric pressure are read in (PBL/JAG
C     format) at selected ADCIRC grid nodes. Interpolation in time is
C     used to synchronize the wind and pressure information with the
C     model time step. Garratt's formula is used to compute wind stress
C     from wind velocity.
         IF(ABS(NWS).EQ.4) THEN
          CALL CLOSEFORT22( ) 
         ENDIF

C     Wind velocity and atmospheric pressure are read in at all grid
C     nodes. Interpolation in time is used to synchronize the wind and
C     pressure information with the model time step. Garratt's formula
C     is used to compute wind stress from wind velocity.
         IF(ABS(NWS).EQ.5) THEN
           CALL CLOSEFORT22( ) 
         ENDIF

C     Wind velocity and atmospheric pressure are read in for a
C     rectangular grid (either in Longitude, Latitude or Cartesian
C     coordinates, consistent with the grid coordinates) and
C     interpolated in space onto the ADCIRC grid and in time to
C     synchronize the wind and pressure information with the model time
C     step. Garratt's formula is used to compute wind stress from the
C     wind velocity.
         IF(NWS.EQ.6) THEN
           CALL CLOSEFORT22( ) 
         ENDIF

C     jgf46.01 New option to read in surface wind stress and atmospheric
C     pressure for a rectangular grid (either in Longitude, Latitude or
C     Cartesian coordinates, consistent with the grid coordinates) and
C     interpolate in space onto the ADCIRC grid. Interpolation in time
C     is used to synchronize the wind and pressure information with the
C     model time step.
         IF(ABS(NWS).EQ.7) THEN
          CALL CLOSEFORT22( ) 
         ENDIF

C     jgf46.02 New option to read in hurricane locations and generate
C     generate hurricane winds from the Holland Wind Model.
         IF(ABS(NWS).EQ.8) THEN
          ! Do nothing
         ENDIF
C
C     Wind velocity (10 m) and atmospheric pressure are read in from a
C     sequence of National Weather Service (NWS) Aviation (AVN) model
C     output files. Each AVN file is assumed to contain data on a
C     Gaussian longitude, latitude grid at a single time. Consecutive
C     files in the sequence are separated by N hours in time. Garratt's
C     formula is used to compute wind stress from the wind velocity.
         IF(NWS.EQ.10) THEN
          ! Do nothing
         ENDIF
C
C     Wind velocity (10 m) and atmospheric pressure are read in from a
C     sequence of stripped down National Weather Service (NWS) ETA 29km
C     model output files. Each ETA file is assumed to contain data on an
C     E grid for a single day (8 data sets, one every 3 hours, beginning
C     @ 03:00 and continuing through 24:00 of the given day). The wind
C     data is converted to an east-west, north-south coordinate system
C     inside ADCIRC. Garratt's formula is used to compute wind stress
C     from the wind velocity.
         IF(NWS.EQ.11) THEN
          ! Do nothing
         ENDIF

C...sb46.28sb01 NWS=12 reads in raw OWI files 09/xx/2006
         IF(ABS(NWS).EQ.12) THEN
          ! CALL NWS12GET(WVNX2,WVNY2,PRN2,NP,PRBCKGRND)
          CALL owiwind_closefiles( ) 
         ENDIF

C     RJW ffpl Merged:
!     rjw added nws = 19: asymmetric hurricane winds
         IF(NWS.EQ.19) THEN
           CALL CLOSEFORT22( ) 
         ENDIF
C
C     jie added nws = 20: generalized asymmetric vortex model
         IF(NWS.EQ.20) THEN
           CALL CLOSEFORT22( ) 
         ENDIF    
C      
C     jgf49.1001 Added NWS29 for embedding an asymmetric vortex from
C     NWS19 into an OWI basin scale met field from NWS12 (derived from NAM).
         IF((ABS(NWS).EQ.29).or.(ABS(NWS).eq.30)) THEN
           ! bring in next set of OWI met data
           ! CALL NWS12GET(WVNX2,WVNY2,PRN2,NP,PRBCKGRND)
           CALL owiwind_closefiles( )   
           ! bring in next set of asymmetric vortex met data to separate arrays
           CALL CLOSEFORT22( ) 
         ENDIF

      !
      ! jgf50.38.05: Added NWS=15 for reading HWind data
         IF(ABS(NWS).EQ.15) THEN
           ! Do nothing
         ENDIF

C.... tcm v51.06.02 added for GFDL Met Data
         IF(ABS(NWS).EQ.16) THEN
           ! do nothing
         ENDIF

         RETURN  
      CONTAINS      

        ! close fort.22
        SUBROUTINE CLOSEFORT22( )
          IMPLICIT NONE

          LOGICAL:: iopened
 
          iopened = .FALSE.  
          INQUIRE( UNIT = 22 + CK_LUN, OPENED = iopened )  
          IF ( iopened ) THEN 
            CLOSE( 22 + CK_LUN )  
          END IF

          RETURN  
        END SUBROUTINE CLOSEFORT22 

      END SUBROUTINE close_met_files

      ! reopen a wind files and fast forward to the location 
      ! previously closed 
      SUBROUTINE reopenffwd_met_files( ) 
        IMPLICIT NONE

        REAL (SZ):: cTimeLoc

        SELECT CASE( IHOT )
        CASE (0)
C           print*, "trying to reopen a met file"  
           cTimeLoc = ITIME*DTDP + STATIM*86400.00  
           CALL coldstartMeteorologicalReopenInputFile( nws, cTimeLoc )
        CASE (17,67,68,367,368,567,568) ! non-portable binary or netcdf hotstart file
           cTimeLoc = (ITIME - ITHS)*DTDP + ITHS*DTDPHS +
     &         STATIM*86400.00  
           CALL hotstartMeteorologicalReopenInputfile( nws, cTimeLoc )
        END SELECT

        RETURN   
      END SUBROUTINE reopenffwd_met_files   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C----------------------------------------------------------------------
C                   S U B R O U T I N E 
C    C O L D  S T A R T   M E T E O R O L O G I C A L   F O R C I N G
C----------------------------------------------------------------------
C     Sets up meteorological i/o for a coldstart.
C----------------------------------------------------------------------
      subroutine coldstartMeteorologicalReopenInputFile( nws, ctimeloc )
!
      use global, only: DEBUG, allMessage
      use sizes, only : gblinputdir, inputdir 
      use owiwind, only : nws12init, nws12get, findStormCenter    ! sb added 09/xx/2006
          
!      use nodalAttributes, only : applyDirectionalWindReduction,
!     &   loadCanopyCoef, loadDirEffRLen, applyCanopyCoefficient
      use mesh, only : dp, ics, np
!
     
!      use global, only: prn1, prn2, PRBCKGRND, setMessageSource, 
!     &   unsetMessageSource

      use wind, only: wvnx1, wvnx2, wvny1, wvny2, prn1, prn2, IWTIME, 
     &    IWTIMEP, IWYR,IREFYR, WTIMED, WREFTIM,  nws3get, nws4get,
     &    nws7get, nws6get, WTIMINC 
  
      implicit none
      integer, intent(in) :: nws
      REAL (SZ):: ctimeloc
        
!      ! local 
!      real(sz) :: windx ! wind vel (m/s) in x-dir
!      real(sz) :: windy ! wind vel (m/s) in y-dir
!      real(sz) :: wtratio ! used to time interpolated between met datasets
!      integer :: nhg ! node number temp
!      real(sz) :: bf ! blend factor for multiple sets of wind data
 
      integer :: i, nttime, idf ! node loop counter     
      real (sz):: timeit, offsettime, wt1, wt2
     
      call setMessageSource("coldStartMeteorologicalForcingInput")
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif      

!      WSX2(:) = 0.d0
!      WSY2(:) = 0.D0
!      PR2(:) = PRBCKGRND_MH2O  
!   
!      PRINT*, "global array dim", ubound( wvnx1 ), ubound( wvnx2 ),
!     &  ubound( wvny1 ), ubound( wvny2 ), ubound( prn1 ), ubound( prn2 )

 
      offsettime = STATIM*86400.D0   
      nttime = INT((ctimeloc - offsettime)/dtdp)  

!     print*, "nttime = ", nttime  

      IF(NWS.EQ.1) THEN
         OPEN(22+CK_LUN,FILE=TRIM(INPUTDIR)//'/'//'fort.22')  
         CALL READFORT22DUMMY(  NP )   
         DO idf = 1, nttime
            CALL READFORT22DUMMY(  NP )  
         END DO
      ENDIF

      IF(ABS(NWS).EQ.2) THEN
         OPEN(22+CK_LUN,FILE=TRIM(INPUTDIR)//'/'//'fort.22')

         CALL READFORT22DUMMY( NP )  
         CALL READFORT22DUMMY( NP )  

         WT1 = STATIM*86400.D0 
         WT2 = WT1 + WTIMINC  
         DO idf = 1, nttime 
            TIMEIT = idf*DTDP + STATIM*86400.D0 ! kmd48.33bc - changed the DTDP to DTDPHS
             IF ( TIMEIT .GE. WT2 ) THEN   
               WT1 = WT2 
               WT2 = WT1 + WTIMINC 
               CALL READFORT22DUMMY( NP )  
             END IF 
         END DO
      ENDIF

      IF(NWS.EQ.3) THEN
C....    tcm_v49.04 20091124 -- Changed from local inputdir to global inputdir
C....                           to correspond with having only a global wind file.
!         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         OPEN(22+CK_LUN,FILE=TRIM(GBLINPUTDIR)//'/'//'fort.22')

 2222    CALL NWS3GET()
         IF(IWYR.NE.IREFYR) THEN
            IWTIMEP=IWTIME
            DO I=1,NP
               WVNX1(I)=WVNX2(I)
               WVNY1(I)=WVNY2(I)
            END DO
            GOTO 2222
         ENDIF
         IF(WTIMED.LE.WREFTIM) THEN
            IWTIMEP=IWTIME
            DO I=1,NP
               WVNX1(I)=WVNX2(I)
               WVNY1(I)=WVNY2(I)
            END DO
            GOTO 2222
         ENDIF
         WT2=WTIMED-WREFTIM  !CAST INTO MODEL TIME REFERENCE
         WT1=WT2-WTIMINC

         DO idf = 1, nttime
           TIMEIT = idf*DTDP + STATIM*86400.D0 ! kmd48.33bc - changed the DTDP to DTDPHS
           IF(TIMEIT.GT.WT2) THEN
               WT1=WT2
               WT2=WT2+WTIMINC
               DO I=1,NP
                  WVNX1(I)=WVNX2(I)
                  WVNY1(I)=WVNY2(I)
               END DO
               CALL NWS3GET()
           ENDIF
         END DO
!
      ENDIF

!
      IF(ABS(NWS).EQ.4) THEN
         OPEN(22+CK_LUN,FILE=TRIM(INPUTDIR)//'/'//'fort.22')

         WT1 = STATIM*86400.D0
         WT2 = WT1 + WTIMINC
         CALL NWS4GET()
         wvnx1 = wvnx2 
         wvny1 = wvny2
         prn1 = prn2 
         CALL NWS4GET()
 
         DO idf=1, nttime
            TIMEIT = idf*DTDPHS + STATIM*86400.D0 ! kmd48.33bc - changed the DTDP to DTDPHS
            IF(TIMEIT.GT.WT2) THEN
               WT1=WT2
               WT2=WT2+WTIMINC
               DO I=1,NP
                  WVNX1(I)=WVNX2(I)
                  WVNY1(I)=WVNY2(I)
                  PRN1(I)=PRN2(I)
               END DO
               CALL NWS4GET()
            ENDIF
         END DO
      ENDIF

      IF(ABS(NWS).EQ.5) THEN
         OPEN(22+CK_LUN,FILE=TRIM(INPUTDIR)//'/'//'fort.22')

         CALL READFORT22DUMMY( NP )  
         CALL READFORT22DUMMY( NP )  
         WT1 = STATIM*86400.D0 
         WT2 = WT1 + WTIMINC
  
         DO idf = 1, nttime 
            TIMEIT = idf*DTDP + STATIM*86400.D0 ! kmd48.33bc - changed the DTDP to DTDPHS
             IF ( TIMEIT .GE. WT2 ) THEN   
               WT1 = WT2 
               WT2 = WT1 + WTIMINC 
               CALL READFORT22DUMMY( NP )  
             END IF 
         END DO
      ENDIF

      IF(NWS.EQ.6) THEN
C....    tcm_v49.04 20091124 -- Changed from local inputdir to global inputdir
C....                           to correspond with having only a global wind file.
!         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         OPEN(22+CK_LUN,FILE=TRIM(GBLINPUTDIR)//'/'//'fort.22')
         CALL NWS6GET()
         wvnx1 = wvnx2
         wvny1 = wvny2
         CALL NWS6GET()
         WT1 = STATIM*86400.D0
         WT2 = WT1 + WTIMINC

         DO IDF = 1, NTTIME
           TIMEIT = idf*DTDP + STATIM*86400.D0 ! kmd48.33bc - changed the DTDP to DTDPHS
           IF(TIMEIT .GT. WT2) THEN
              WT1=WT2
              WT2=WT2 + WTIMINC
              DO I=1,NP
                 WVNX1(I) =WVNX2(I)
                 WVNY1(I) =WVNY2(I)
                 PRN1(I)= PRN2(I)
              END DO
!              NWSGGWI=NWSGGWI+1
              CALL NWS6GET()
          ENDIF
        END DO     

      ENDIF

C     jgf46.00 Added option to directly apply surface stress without any
C     other correction factors.
      IF(ABS(NWS).EQ.7) THEN
         OPEN(22+CK_LUN,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         CALL NWS7GET()  !TCMv49.02 Changed (WVNX2,WVNY2,PRN2) to (WVNX1,WVNY1,PRN1)
         wvnx1 = wvnx2
         wvny1 = wvny2

         CALL NWS7GET()
         WT1 = STATIM*86400.D0 
         WT2 = WT1 + WTIMINC
         !
         DO IDF = 1, NTTIME
           TIMEIT = idf*DTDP + STATIM*86400.D0 ! kmd48.33bc - changed the DTDP to DTDPHS
           IF ( TIMEIT .GT. WT2 ) THEN
            WT1=WT2
            WT2=WT2 + WTIMINC
            DO I=1,NP
               WVNX1(I)=WVNX2(I)
               WVNY1(I)=WVNY2(I)
               PRN1(I)=PRN2(I)
            END DO
            CALL NWS7GET()
          ENDIF
        END DO     
        !
      ENDIF

C     jgf46.02 New option to read in hurricane locations and generate
C     generate hurricane winds from the Holland Wind Model.
      IF(ABS(NWS).EQ.8) THEN
         ! no need to do anything
!         CALL HollandGet(WVNX2,WVNY2,PRN2,statim*86400.d0)
      ENDIF

C     RJW Merged:
!    rjw added nws = 19: asymmetric hurricane winds v2.0
      IF(NWS.EQ.19) THEN
        OPEN( 22+CK_LUN, FILE=TRIM(GBLINPUTDIR)//'/'//'fort.22',STATUS ='OLD'  )
!        CALL NWS19GET(WVNX2,WVNY2,PRN2,STATIM*86400.D0)
      ENDIF
      
!     jie added nws = 20: generalized asymmetric vortex model
      IF(NWS.EQ.20) THEN
         OPEN(22+CK_LUN, FILE=TRIM(GBLINPUTDIR)//'/'//'fort.22',STATUS ='OLD')
!        CALL NWS20GET(WVNX2,WVNY2,PRN2,STATIM*86400.D0)
      ENDIF      

C.....sb46.28sb01 NWS=-12 and 12 was added to deal with raw OWI files.  09/xx/2006
      ! jgf49.1001 Added NWS=29, asymmetric vortex wind (NWS19) embedded in
      ! OWI (NWS12) basin scale wind field.
      IF ((ABS(NWS).EQ.12).OR.(ABS(NWS).EQ.29).OR.(ABS(NWS).EQ.30)) THEN
        CALL NWS12INIT(WVNX1,WVNY1,PRN1,NP,PRBCKGRND)
        CALL NWS12GET(WVNX1,WVNY1,PRN1,NP,PRBCKGRND)
        CALL NWS12GET(WVNX2,WVNY2,PRN2,NP,PRBCKGRND)

        WT1 = STATIM*86400 
        WT2 = WT1 + WTIMINC 
            
        DO IDF = 1, NTTIME
           TIMEIT = idf*DTDP + STATIM*86400.D0 ! kmd48.33bc - changed the DTDP to DTDPHS
           IF ( TIMEIT .GT. WT2 ) THEN
            WT1=WT2
            WT2=WT2 + WTIMINC
 
            DO I=1,NP
               WVNX1(I)=WVNX2(I)
               WVNY1(I)=WVNY2(I)
               PRN1(I)=PRN2(I)
            END DO

            CALL NWS12GET( WVNX2, WVNY2, PRN2, NP, PRBCKGRND )
          ENDIF
       END DO     
      ENDIF
      !
      ! jgf49.1001 Added NWS=29, asymmetric vortex wind (NWS19) embedded in
      ! OWI (NWS12) basin scale wind field.
      ! jgf: Added NWS=30 (GAHM+OWI)
      IF ((abs(NWS).EQ.29).or.(abs(NWS).eq.30)) THEN
         ! now blend in the vortex met field
         if (abs(nws).eq.29) then
            OPEN(22+CK_LUN,FILE=TRIM(GBLINPUTDIR)//'/'//'NWS_19_fort.22',
     &           STATUS='OLD')
         endif
         if (abs(nws).eq.30) then
            OPEN(22+CK_LUN,FILE=TRIM(GBLINPUTDIR)//'/'//'NWS_20_fort.22',
     &           STATUS='OLD')
         endif
      ENDIF
      !
      ! jgf51.52.21: Fixed NWS10 and cleaned up interface.
      IF(NWS.EQ.10) THEN
!
      ENDIF

      IF(NWS.EQ.11) THEN
!
      ENDIF
      !
      ! jgf50.38.05: Added NWS=15,-15 for reading HWind data
      IF (ABS(NWS).EQ.15) THEN
!         CALL NWS15INIT(STATIM*86400.0)
      ENDIF
C     
C....tcm v51.06.02: Added NWS=16,-16 for reading ASCII GFDL Met data
      IF (ABS(NWS).EQ.16) THEN
!         CALL INIT_GFDL(STATIM*86400.0)
      ENDIF
      !
      ! jgf52.14: Implementing Casey's fix for uninitialized
      ! wvnx1 and wvny1 arrays  they are not used/needed/initialized
      ! by the parametric vortex meteorological models but are used in
      ! the padcswan_init() subroutine to initialize the SWAN wind 
      ! velocities. Failure to initialize them causes NaNs in the 
      ! solution. 
!      if ( (abs(nws).eq.8).or.(abs(nws).eq.19).or.(abs(nws).eq.20)) then
!         WVNX1 = WVNX2
!         WVNY1 = WVNY2
!      endif

#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()

      CONTAINS

        SUBROUTINE READFORT22DUMMY( NNP )
          IMPLICIT NONE
             
          INTEGER, INTENT(IN):: NNP

          INTEGER:: I, IDUMMY
          REAL (SZ):: RDUMMY1, RDUMMY2, RDUMMY3
                   
          DO I = 1, NNP
             READ(22+CK_LUN,*) IDUMMY, RDUMMY1, RDUMMY2, RDUMMY3
!            READ(22+CK_LUN,*) (NHG,WVNX1(I),WVNY1(I),PRN1(I),I=1,NP)
          END DO

          RETURN  
        END SUBROUTINE READFORT22DUMMY
 

C----------------------------------------------------------------------
      end subroutine coldstartMeteorologicalReopenInputFile
C----------------------------------------------------------------------



C----------------------------------------------------------------------
C                   S U B R O U T I N E 
C    H O T S T A R T   M E T E O R O L O G I C A L   F O R C I N G
C----------------------------------------------------------------------
C     Sets up meteorological i/o for a hotstart.
C----------------------------------------------------------------------
      subroutine hotstartMeteorologicalReopenInputfile( nws, ctimeloc )
      use sizes, only : gblinputdir, inputdir
!     
      use global, only: PRBCKGRND,  setMessageSource, DEBUG,
     &   unsetMessageSource, allMessage

      use wind, only: wvnx1, wvnx2, wvny1, wvny2, IWTIME, 
     &    IWTIMEP, IWYR,IREFYR, WTIMED, WREFTIM,  nws3get, nws4get,
     &    nws7get, nws6get, wtiminc, prn1, prn2 
 
      use owiwind,only : nws12init, nws12get, findStormCenter    ! sb added 09/xx/2006
!      use nodalAttributes, only : applyDirectionalWindReduction,
!     &   loadCanopyCoef, loadDirEffRLen, applyCanopyCoefficient
      use mesh, only : dp, ics, np


      implicit none
      integer, intent(in) :: nws
      real(sz), intent(in) :: ctimeloc
!      real(sz), intent(in) :: rampMete2   ! Ramp for wind and atmospheric pressure @ITHS  

 
      ! local   
      real(sz) :: timeit, wt1, wt2       
!      real(sz) :: windx ! wind vel (m/s) in x-dir
!      real(sz) :: windy ! wind vel (m/s) in y-dir
!      real(sz) :: prdiff
!
!      real(sz) :: wDragCo
!      real(sz) :: windMag
!      real(sz) :: wtratio ! used to time interpolated between met datasets
!      integer :: nhg ! node number temp
!      real(sz) :: bf ! blend factor for multiple sets of wind data 
!      real(sz) :: vortexPRDIFF
!      integer :: it
      integer :: i, j ! node loop counter
!      integer :: j

      integer:: idf, nttime  
      REAL (SZ):: HsTimeLoc     

      call setMessageSource("hotStartMeteorologicalReopen")
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif      

C
C......RESTART WIND AND PRESSURE INFORMATION
       HsTimeLoc=ITHS*DTDPHS + StaTim*86400.D0

     
      nttime = INT((ctimeloc - hstimeloc)/DTDP)  
C     No wind, radiation stress or atmospheric pressure forcings are used.
      IF(NWS.EQ.1) THEN
         OPEN(22+CK_LUN,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         CALL READFORT22DUMMY(  ) 
         DO J=1, (ITHS + nttime)
           CALL READFORT22DUMMY(  )
         END DO
      ENDIF

C     Wind stress and atmospheric pressure are read in at all grid nodes
C     at a time interval that does not equal the model time
C     step. Interpolation in time is used to synchronize the wind and
C     pressure information with the model time step.
      IF(NWS.EQ.2) THEN
         OPEN(22+CK_LUN,FILE=TRIM(INPUTDIR)//'/'//'fort.22')

         WT1 = STATIM*86400.D0
         WT2 = WT1 + WTIMINC
         CALL READFORT22DUMMY(  ) 
         CALL READFORT22DUMMY(  ) 

         CALL NWSGETADVANCE( ITHS + nttime, GETTIMEIT_CS,
     &       READFORT22DUMMY, .TRUE.  )   

       
!         DO IT=1, ITHS + nttime
!           TIMEIT = GETTIMEIT_CS( IT )  
!
!           IF( TIMEIT .GT. WT2 ) THEN
!              WT1=WT2
!              WT2=WT2+WTIMINC
!              CALL READFORT22DUMMY(  NP ) 
!           ENDIF
!         END DO

      ENDIF

      IF(NWS.EQ.-2) THEN
         OPEN(22+CK_LUN,FILE=TRIM(INPUTDIR)//'/'//'fort.22')

         WT1 = HsTimeLoc
         WT2 = WT1 + WTIMINC
         CALL READFORT22DUMMY(  ) 
         CALL READFORT22DUMMY(  ) 

         CALL NWSGETADVANCE( nttime, GETTIMEIT_HS, 
     &       READFORT22DUMMY, .TRUE. )   

!         DO IT=1, nttime 
!            TIMEIT=IT*DTDP + hstimeloc ! kmd48.33bc - changed the DTDP to DTDPHS
!            
!            IF( TIMEIT .GT. WT2 ) THEN
!               WT1=WT2
!               WT2=WT2+WTIMINC
!               DO I=1,NP
!                CALL READFORT22DUMMY(  NP ) 
!              END DO
!            ENDIF
!         END DO

      ENDIF

C     Wind velocity in US Navy Fleet Numeric format interpolated in
C     space onto the ADCIRC grid and in time to synchronize the wind and
C     pressure information with the model time step. Garratt's formula
C     is used to compute wind stress from the wind velocity.
      IF(NWS.EQ.3) THEN
C....    tcm_v49.04 20091124 -- Changed from local inputdir to global inputdir
C....                           to correspond with having only a global wind file.
!         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         OPEN(22+CK_LUN,FILE=TRIM(GBLINPUTDIR)//'/'//'fort.22')
 2223    CALL NWS3GET()
         IF(IWYR.NE.IREFYR) THEN
            IWTIMEP=IWTIME
            DO I=1,NP
               WVNX1(I)=WVNX2(I)
               WVNY1(I)=WVNY2(I)
            END DO
            GOTO 2223
         ENDIF
         IF(WTIMED.LE.WREFTIM) THEN
            IWTIMEP=IWTIME
            DO I=1,NP
               WVNX1(I)=WVNX2(I)
               WVNY1(I)=WVNY2(I)
            END DO
            GOTO 2223
         ENDIF
         WT2=WTIMED-WREFTIM  !CAST INTO MODEL TIME REFERENCE
         WT1=WT2-WTIMINC

         CALL NWSGETADVANCE( ITHS + nttime, GETTIMEIT_CS, NWS3GET  )   

!         DO IT=1,ITHS + nttime
!           TIMEIT = GETTIMEIT_CS( IT )  
!           IF( TIMEIT .GT. WT2 ) THEN
!               WT1=WT2
!               WT2=WT2+WTIMINC
!               DO I=1,NP
!                  WVNX1(I)=WVNX2(I)
!                  WVNY1(I)=WVNY2(I)
!               END DO
!               CALL NWS3GET()
!           ENDIF
!        END DO
      ENDIF

C     Wind velocity and atmospheric pressure are read in (PBL/JAG
C     format) at selected ADCIRC grid nodes. Interpolation in time is
C     used to synchronize the wind and pressure information with the
C     model time step. Garratt's formula is used to compute wind stress
C     from wind velocity.
      IF(NWS.EQ.4) THEN
         OPEN(22+CK_LUN,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         WT1 = STATIM*86400.D0
         WT2 = WT1 + WTIMINC         
         CALL NWS4GET()  !TCM v49.02 Added read for rec. 1
         wvnx1 = wvnx2
         wvny1 = wvny2
         CALL NWS4GET()

         CALL NWSGETADVANCE( ITHS + nttime, GETTIMEIT_CS, NWS4GET  )   
!
!         DO IT=1,ITHS + nttime
!            TIMEIT=GETTIMEIT_CS( IT ) 
!            IF(TIMEIT .GT. WT2) THEN
!               WT1=WT2
!               WT2=WT2+WTIMINC
!               DO I=1,NP
!                  WVNX1(I)=WVNX2(I)
!                  WVNY1(I)=WVNY2(I)
!                  PRN1(I)=PRN2(I)
!               END DO
!               CALL NWS4GET()
!            ENDIF
!         END DO
!
      ENDIF


      IF(NWS.EQ.-4) THEN
         OPEN(22+CK_LUN,FILE=TRIM(INPUTDIR)//'/'//'fort.22')

         WT1 = HsTimeLoc  
         WT2 = WT1 + WTIMINC   
         CALL NWS4GET()
         wvnx1 = wvnx2
         wvny1 = wvny2
         CALL NWS4GET()

         CALL NWSGETADVANCE( nttime, GETTIMEIT_HS, NWS4GET  )   
!
!         DO IT = 1, nttime
!           TIMEIT = GETTIMEIT_HS( RIT= IT )
!           IF(TIMEIT .GT. WT2) THEN
!               WT1=WT2
!               WT2=WT2+WTIMINC
!               DO I=1,NP
!                  WVNX1(I)=WVNX2(I)
!                  WVNY1(I)=WVNY2(I)
!                  PRN1(I)=PRN2(I)
!               END DO
!               CALL NWS4GET()
!            ENDIF
!         END DO
! 
      ENDIF

C     Wind velocity and atmospheric pressure are read in at all grid
C     nodes. Interpolation in time is used to synchronize the wind and
C     pressure information with the model time step. Garratt's formula
C     is used to compute wind stress from wind velocity.
      IF(NWS.EQ.5) THEN
         OPEN(22+CK_LUN,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         WT1 = STATIM*86400.D0
         WT2 = WT1 + WTIMINC
         CALL READFORT22DUMMY(  ) 
         CALL READFORT22DUMMY(  ) 
         
         CALL NWSGETADVANCE( ITHS + nttime, GETTIMEIT_CS, 
     &       READFORT22DUMMY,  .TRUE. )   

!          DO IT=1,ITHS + nttime
!            TIMEIT = GETTIMEIT_CS( IT )  
!            IF(TIMEIT .GT. WT2) THEN
!               WT1=WT2
!               WT2=WT2 + WTIMINC
!               CALL READFORT22DUMMY(  NP ) 
!            ENDIF
!         END DO
      ENDIF

      IF(NWS.EQ.-5) THEN
         OPEN(22+CK_LUN,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
!         WTIME1 = TimeLoc
!         WTIME2 = WTIME1 + WTIMINC
!         READ(22+CK_LUN,*) (NHG,WVNX1(I),WVNY1(I),PRN1(I),I=1,NP)
!         READ(22+CK_LUN,*) (NHG,WVNX2(I),WVNY2(I),PRN2(I),I=1,NP)
         WT1 = HsTimeLoc
         WT2 = WT1 + WTIMINC
         CALL READFORT22DUMMY( ) 
         CALL READFORT22DUMMY( ) 
        
         CALL NWSGETADVANCE( nttime, GETTIMEIT_HS, 
     &       READFORT22DUMMY,  .TRUE. )   
 
!         DO I = 1, nttime
!            TIMEIT = GETTIMEIT_HS( IT )  
!            IF(TIMEIT .GT. WT2) THEN
!               WT1=WT2
!               WT2=WT2 + WTIMINC
!               CALL READFORT22DUMMY(  NP ) 
!            ENDIF
!         END DO

      ENDIF

C     Wind velocity and atmospheric pressure are read in for a
C     rectangular grid (either in Longitude, Latitude or Cartesian
C     coordinates, consistent with the grid coordinates) and
C     interpolated in space onto the ADCIRC grid and in time to
C     synchronize the wind and pressure information with the model time
C     step. Garratt's formula is used to compute wind stress from the
C     wind velocity.
      IF(NWS.EQ.6) THEN
C....    tcm_v49.04 20091124 -- Changed from local inputdir to global inputdir
C....                           to correspond with having only a global wind file.
!         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         OPEN(22+CK_LUN,FILE=TRIM(GBLINPUTDIR)//'/'//'fort.22')
C     The following 3 lines are a hardwire to allow a non standard met
C     file to be read in at time zero in a hot start.  They should be
C     eliminated or commented out for normal operation
c     OPEN(199,FILE=TRIM(INPUTDIR)//'/'//'fort.199')
c     READ(199,*) (NHG,PRN1(I),WVNX1(I),WVNY1(I),I=1,NP)
c     CLOSE(199)
C     The following CALL statement should be uncommented for normal operation
         CALL NWS6GET()
         wvnx1 = wvnx2 
         wvny1 = wvny2 
         prn1 = prn2   
 
         CALL NWS6GET()
         WT1 = HsTimeLoc  
         WT2 = WT1 + WTIMINC 

         CALL NWSGETADVANCE( nttime, GETTIMEIT_HS, NWS6GET )   

!
!         DO IT = 1, nttime
!            TIMEIT = GETTIMEIT_HS( IT )  
!            IF(TIMEIT .GT. WT2) THEN
!              WT1=WT2
!              WT2=WT2 + WTIMINC
!              DO I=1,NP
!                WVNX1(I)=WVNX2(I)
!                WVNY1(I)=WVNY2(I)
!                PRN1(I)=PRN2(I)
!              END DO
!              CALL NWS6GET()
!            ENDIF
!         END DO
!
      ENDIF

C     jgf46.01 New option to read in surface wind stress and atmospheric
C     pressure for a rectangular grid (either in Longitude, Latitude or
C     Cartesian coordinates, consistent with the grid coordinates) and
C     interpolate in space onto the ADCIRC grid. Interpolation in time
C     is used to synchronize the wind and pressure information with the
C     model time step.
C
      IF(NWS.EQ.7) THEN
         OPEN(22+CK_LUN,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         WT1 = STATIM*86400.D0
         CALL NWS7GET()  !TCMv49.02 Changed (WVNX2,WVNY2,PRN2) to (WVNX1,WVNY1,PRN1)
         wvnx1 = wvnx2
         wvny1 = wvny2
         prn1 = prn2 
         CALL NWS7GET()
         WT1 = WT2   
         WT2 = WT1 + WTIMINC   
         !
         CALL NWSGETADVANCE( ITHS + nttime, GETTIMEIT_CS, NWS7GET )   
         
!         DO IT=1, ITHS + nttime
!            TIMEIT=GETTIMEIT_CS( IT ) 
!            IF(TIMEIT .GT. WT2) THEN
!               WT1=WT2
!               WT2=WT2+WTIMINC
!               DO I=1,NP
!                  WVNX1(I)=WVNX2(I)
!                  WVNY1(I)=WVNY2(I)
!                  PRN1(I)=PRN2(I)
!               END DO
!               CALL NWS4GET()
!            ENDIF
!         END DO
! 

      ENDIF

      IF(NWS.EQ.-7) THEN
         OPEN(22+CK_LUN,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
        
         CALL NWS7GET()
         wvnx1 = wvnx2
         wvny1 = wvny2
         prn1 = prn2  
         CALL NWS7GET()

         WT1=HsTimeLoc
         WT2=WT1+WTIMINC

         !
         CALL NWSGETADVANCE( nttime, GETTIMEIT_HS, NWS7GET )     
      ENDIF

C     jgf46.02 New option to read in hurricane locations and generate
C     generate hurricane winds from the Holland Wind Model.
      IF(ABS(NWS).EQ.8) THEN
!        No need to fast forward
      ENDIF

!     rjw added nws = 19: asymmetric hurricane winds v2.0
      IF(ABS(NWS).EQ.19) THEN
!        No need to fast forward
      ENDIF
      
!     jie added nws = 20: generalized asymmetric vortex model
      IF(ABS(NWS).EQ.20) THEN
!       No need to fast forward
      ENDIF            
C
C     Wind velocity (10 m) and atmospheric pressure are read in from a
C     sequence of National Weather Service (NWS) Aviation (AVN) model
C     output files. Each AVN file is assumed to contain data on a
C     Gaussian longitude, latitude grid at a single time. Consecutive
C     files in the sequence are separated by N hours in time. Garratt's
C     formula is used to compute wind stress from the wind velocity.
      IF(NWS.EQ.10) THEN
!        No need to fast forward
      ENDIF

C     Wind velocity (10 m) and atmospheric pressure are read in from a
C     sequence of stripped down National Weather Service (NWS) ETA 29km
C     model output files. Each ETA file is assumed to contain data on an
C     E grid for a single day (8 data sets, one every 3 hours, beginning
C     @ 03:00 and continuing through 24:00 of the given day). The wind
C     data is converted to an east-west, north-south coordinate system
C     inside ADCIRC. Garratt's formula is used to compute wind stress
C     from the wind velocity.
      IF(NWS.EQ.11) THEN
!        No need to fast forward
      ENDIF
C
C.....sb46.28sb01 NWS=12,-12 were added to deal with raw OWI files.  09/xx/2006
      ! jgf: Combined NWS=12 and NWS=-12 and added abs(NWS)=29 and abs(NWS)=30
      IF((abs(NWS).eq.12).or.(abs(NWS).eq.29).or.(abs(NWS).eq.30)) THEN
        WT1 = HsTimeLoc      ! nws=-12, wind data relative to hot start time
        if (nws.gt.0) then
           WT1 = STATIM*86400.D0 ! wind data relative to cold start time
        endif
        WT2 = WT1 + WTIMINC

        CALL NWS12INIT(WVNX1,WVNY1,PRN1,NP,PRBCKGRND)  ! This just initializes some variables and might set WVNX1,etc... if skipping ahead.
        CALL NWS12GET(WVNX1,WVNY1,PRN1,NP,PRBCKGRND)   !TCM v49.02 Added read for rec. 1
        CALL NWS12GET(WVNX2,WVNY2,PRN2,NP,PRBCKGRND)   ! read recond 2

        ! fast forward if wind data are relative to cold start time
        if (nws.gt.0) then
           CALL NWSGETADVANCE( ITHS + nttime, GETTIMEIT_CS, NWS12GETDUMMY )  
        else
           CALL NWSGETADVANCE( nttime, GETTIMEIT_HS, NWS12GETDUMMY )  
        endif
        !          
      ENDIF
C
C     jgf49.1001 Added NWS29, asymmetric vortex winds (NWS19)
C     embedded in an OWI basin scale gridded met field (NWS12)
C     from NAM data
      ! jgf: Added NWS30 (GAHM+OWI)
      IF ((abs(NWS).EQ.29).or.(abs(NWS).eq.30)) THEN
         ! Now blend in the vortex winds into the background OWI wind field.
         if (abs(nws).eq.29) then
            OPEN(22+CK_LUN,FILE=TRIM(GBLINPUTDIR)//'/'//'NWS_19_fort.22',
     &           STATUS='OLD')
         endif
         if (abs(nws).eq.30) then
            OPEN(22+CK_LUN,FILE=TRIM(GBLINPUTDIR)//'/'//'NWS_20_fort.22',
     &           STATUS='OLD')
         endif         
      ENDIF
      !
      ! jgf50.38.05: Added NWS=15 for reading HWind data
      IF(ABS(NWS).EQ.15) THEN
!          CALL NWS15INIT(timeloc)
!          CALL NWS15GET(WVNX2,WVNY2,PRN2,timeloc)
      ENDIF
C
C.....tcm v51.06.02 NWS=16,-16 were added to deal with GFDL Met files.
      IF(NWS.EQ.16) THEN
        !  Appear that there is no need to do anything (?)
!          WT1 = STATIM*86400.D0
!          WT2 = WT1 + WTIMINC
!          CALL INIT_GFDL(wt1)  !this will read in the two time series
!          DO IT=1,ITHS
!             TIMEIT=IT*DTDPHS + STATIM*86400.D0 ! kmd48.33bc - changed the DTDP to DTDPHS
!             IF(TIMEIT.GT.WT2) THEN
!                WT1=WT2
!                WT2=WT2+WTIMINC
!                DO I=1,NP
!                   WVNX1(I)=WVNX2(I)
!                   WVNY1(I)=WVNY2(I)
!                   PRN1(I)=PRN2(I)
!                 END DO
!                 CALL NWS16GET(wt2,WVNX2,WVNY2,PRN2)
!              ENDIF
!          END DO
      ENDIF
C
      IF(NWS.EQ.-16) THEN
        ! Appear that there is not need to dp anything (?)
!          WT1 = TimeLoc
!          WT2 = WT1 + WTIMINC
!          call INIT_GFDL(wt1)  !this reads in the first two time series values
      ENDIF
      !
      ! jgf52.14: Implementing Casey's fix for uninitialized
      ! wvnx1 and wvny1 arrays  they are not used/needed/initialized
      ! by the parametric vortex meteorological models but are used in
      ! the padcswan_init() subroutine to initialize the SWAN wind 
      ! velocities. Failure to initialize them causes NaNs in the 
      ! solution. 
!      if ( (abs(nws).eq.8).or.(abs(nws).eq.19).or.(abs(nws).eq.20)) then 
!         WVNX1 = WVNX2
!         WVNY1 = WVNY2
!      endif
      
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
    
      CONTAINS
        
        !
        SUBROUTINE READFORT22DUMMY(  )
          IMPLICIT NONE 
            
          INTEGER:: I, IDUMMY
          REAL (SZ):: RDUMMY1, RDUMMY2, RDUMMY3
                   
          DO I = 1, NP
             READ(22+CK_LUN,*) IDUMMY, RDUMMY1, RDUMMY2, RDUMMY3
          END DO

          RETURN  
        END SUBROUTINE READFORT22DUMMY

        !
        SUBROUTINE NWS12GETDUMMY( )
          IMPLICIT NONE

          CALL NWS12GET( WVNX2, WVNY2, PRN2, NP, PRBCKGRND )

          RETURN  
        END SUBROUTINE NWS12GETDUMMY

        !
        FUNCTION GETTIMEIT_CS( IT ) RESULT( TIMEIT ) 
          IMPLICIT NONE 
               
          REAL (8):: TIMEIT
          INTEGER:: IT
                
          IF ( IT <= ITHS ) THEN
             TIMEIT=IT*DTDPHS + STATIM*86400.D0 ! kmd48.33bc - changed the DTDP to DTDPHS
          ELSE
             TIMEIT=(IT - ITHS)*DTDP + HsTimeLoc
          END IF

          RETURN  
       END FUNCTION GETTIMEIT_CS

       !
       FUNCTION GETTIMEIT_HS( RIT ) RESULT( TIMEIT ) 
         IMPLICIT NONE

         REAL (8):: TIMEIT 
         INTEGER:: RIT
             
         TIMEIT = RIT*DTDP + HsTimeLoc  

         RETURN 
       END FUNCTION GETTIMEIT_HS  
 
       !
       SUBROUTINE NWSGETADVANCE( NTSTEP, GETTIMEIT, NWSGET, FLAGSKIP )
         IMPLICIT NONE 
        
         INTEGER:: NTSTEP
         LOGICAL, OPTIONAL:: FLAGSKIP 

         INTERFACE 
            FUNCTION GETTIMEIT( IT ) RESULT( TIMEIT ) 
              IMPLICIT NONE 
               
              REAL (8):: TIMEIT
              INTEGER:: IT        
            END FUNCTION GETTIMEIT

             SUBROUTINE NWSGET( )
               IMPLICIT NONE
             END SUBROUTINE NWSGET     
         END INTERFACE

         INTEGER:: IT
         LOGICAL:: ASKIP         

         ASKIP = .FALSE. 
         IF ( PRESENT(FLAGSKIP) ) THEN
           ASKIP = FLAGSKIP   
         END IF

         DO IT=1, NTSTEP
            TIMEIT=REAL(GETTIMEIT( IT ),SZ)   
 
            IF(TIMEIT .GT. WT2) THEN
              WT1=WT2
              WT2=WT2+WTIMINC
              
              IF ( .NOT. ASKIP ) THEN 
                 DO I=1,NP
                   WVNX1(I)=WVNX2(I)
                   WVNY1(I)=WVNY2(I)
                   PRN1(I)=PRN2(I)
                 END DO
              END IF

              CALL NWSGET()
           ENDIF
         END DO
 
         RETURN 
       END SUBROUTINE NWSGETADVANCE

C--------------------END MET FORCING-----------------------------------
      end subroutine hotstartMeteorologicalReopenInputFile
C--------------------END MET FORCING-----------------------------------

C----------------------------------------------------------------------
C----------------------------------------------------------------------
      ! close 
      SUBROUTINE CLOSEFUNIT( FUNIT )
        IMPLICIT NONE

        INTEGER:: FUNIT

        LOGICAL:: iopened
 
        iopened = .FALSE.  
        INQUIRE( UNIT = FUNIT, OPENED = iopened )  
        IF ( iopened ) THEN 
          CLOSE( FUNIT )  
        END IF

        RETURN  
      END SUBROUTINE CLOSEFUNIT

      ! cold start time, 
      FUNCTION CSTIME_CS0( IT, StLocTime ) RESULT( TIMEIT )
        IMPLICIT NONE

        INTEGER:: IT
        REAL (8):: TIMEIT, StLocTime

        TIMEIT=IT*DTDP + StLocTime 
         
        RETURN           
      END FUNCTION CSTIME_CS0

      ! hot start time, IT - relative to the cold strart time
      FUNCTION CSTIME_CS1( IT, StLocTime ) RESULT( TIMEIT )
        IMPLICIT NONE

        INTEGER:: IT
        REAL (SZ):: TIMEIT, StLocTime

        IF ( IT <= ITHS ) THEN
          TIMEIT = IT*DTDP + StLocTime  
        ELSE 
          TIMEIT = (IT - ITHS)*DTDP + ITHS*DTDPHS + StLocTime 
        ENDIF 

        RETURN     
      END FUNCTION CSTIME_CS1

      ! hot start time, IT - relative to the hot start time
      FUNCTION CSTIME_HS( IT, StLocTime ) RESULT( TIMEIT )
        IMPLICIT NONE

        INTEGER:: IT
        REAL (SZ):: TIMEIT, StLocTime

         
        TIMEIT= IT*DTDP + StLocTime 
           
        RETURN     
      END FUNCTION CSTIME_HS
C----------------------------------------------------------------------
C---------------------------------------------------------------------

!C---------------------------------------------------------------------
!C---------------------------------------------------------------------
!C--------- Elevation & River 
!C---------------------------------------------------------------------
!C---------------------------------------------------------------------

!c-------------
      SUBROUTINE close_elevbc_riverbc_files( )
        IMPLICIT NONE

        INTEGER:: FUNIT 

        ! Time dependent elevation boundary
        IF((NOPE.GT.0).AND.(NBFR.EQ.0)) THEN
!         ! NCSU Subdomain Modeling
          if (subdomainOn.eqv..false.) then
            FUNIT = 19 + CK_LUN  
            CALL CLOSEFUNIT( FUNIT )   
          endif
        ENDIF

        ! Time dependent flux bouundary
        IF((NFLUXF.EQ.1).AND.
     &         ((NFFR.EQ.0).OR.(NFFR.EQ.-1))) THEN         
          FUNIT = 20 + CK_LUN   
          CALL CLOSEFUNIT( FUNIT )  
        ENDIF

        RETURN  
      END SUBROUTINE close_elevbc_riverbc_files

!c------------- Reopen and fastforward
      SUBROUTINE Reopen_elevbc_riverbc_InputFile(  )
        use sizes, only : inputdir 

        IMPLICIT NONE

        REAL (SZ):: StartTimeLoc
        INTEGER:: NTTIME

        ! Open boundary
        IF((NOPE.GT.0).AND.(NBFR.EQ.0)) THEN
         !
         if ((subdomainOn.eqv..false.).or.(enforceBN.ne.1)) then ! NCSU Subdomain
            OPEN(19+CK_LUN,FILE=TRIM(INPUTDIR)//'/'//'fort.19')
            !
            StartTimeLoc = STATIM*86400.00
            NTTIME = itime   
            SELECT CASE(IHOT)
            CASE (0) 
              CALL FFWDFORT19( NTTIME, CSTIME_CS0, StartTimeLoc )  
            CASE (17,67,68,367,368,567,568)
              CALL FFWDFORT19( NTTIME, CSTIME_CS1, StartTimeLoc )  
            END SELECT
            !
          endif
        ENDIF

C......FIND PROPER PLACE IN THE APERIODIC NORMAL FLOW BOUNDARY CONDITION
C......FILE IF IT IS REQUIRED
        IF((NFLUXF.EQ.1).AND.(NFFR.EQ.0)) THEN
          OPEN(20+CK_LUN,FILE=TRIM(INPUTDIR)//'/'//'fort.20') 
          !
          StartTimeLoc = STATIM*86400.00
          NTTIME = itime   
          SELECT CASE( IHOT ) 
          CASE (0)
            CALL FFWDFORT20( NTTIME, CSTIME_CS0, StartTimeLoc )
          CASE (17,67,68,367,368,567,568)
            CALL FFWDFORT20( NTTIME, CSTIME_CS1, StartTimeLoc )
          END SELECT
        ENDIF

! kmd - added in a option to start the river information at the hot start time
            ! StartTimeLocc = ITHS*DTDPHS + STATIM*86400  
        IF((NFLUXF.EQ.1).AND.(NFFR.EQ.-1)) THEN

          StartTimeLoc = ITHS*DTDPHS + STATIM*89400.0   
          SELECT CASE( IHOT )
          CASE (0)
            NTTIME = ITIME
            StartTimeLoc = STATIM*86400.D0
            CALL FFWDFORT20( NTTIME, CSTIME_CS0, StartTimeLoc )
          CASE (17,67,68,367,368,567,568)
            StartTimeLoc = ITHS*DTDPHS + STATIM*86400.0  
            NTTIME = itime - ITHS  
            CALL FFWDFORT20( NTTIME, CSTIME_HS, StartTimeLoc )  
          END SELECT
            
        ENDIF

      CONTAINS

        SUBROUTINE READFORT19SNAP( ) 
          IMPLICIT NONE

          INTEGER:: J
          REAL (SZ):: ESBINTMP

          DO J=1, NETA
           READ(19+CK_LUN,*) ESBINTMP
          END DO
          
          RETURN 
        END SUBROUTINE READFORT19SNAP

        SUBROUTINE FFWDFORT19( NTTIME, CSTIME, StartLocTime )
          IMPLICIT NONE
              
          INTEGER, INTENT(IN):: NTTIME
          REAL (SZ):: StartLocTime

          INTERFACE 
            FUNCTION CSTIME( II, StartTimeLoc ) RESULT( TIMEIT )
              IMPLICIT NONE
                   
              INTEGER:: II
              REAL (8):: TIMEIT
              REAL (8):: StartTimeLoc

            END FUNCTION  
          END INTERFACE
 
          INTEGER:: J, IT 
          REAL (SZ):: ET1, ET2, ETIMINC_TMP, ESBINTMP, TIMEIT
    
          READ(19+CK_LUN,*) ETIMINC_TMP
          ET1=StartLocTime 
          ET2=ET1 + ETIMINC_TMP
          CALL READFORT19SNAP( )  
          CALL READFORT19SNAP( )  

          DO IT=1, NTTIME
            ! TIMEIT=(IT - ITHS)*DTDP + ITHS*DTDPHS + STATIM*86400.D0 -- hotstart
            TIMEIT = CSTIME( IT, StartTimeLoc )   
            IF(TIMEIT .GT. ET2) THEN
              ET1 = ET2
              ET2 = ET1 + ETIMINC_TMP
              CALL READFORT19SNAP( )  
            ENDIF
          END DO
                
          RETURN  
        END SUBROUTINE FFWDFORT19

        SUBROUTINE READFORT20SNAP( )
          IMPLICIT NONE

          INTEGER:: J
          REAL (SZ):: QNIN1TMP 

          DO J=1, NVEL
            QNIN1TMP = 0.D0
            IF( (LBCODEI(J).EQ.2) .OR. (LBCODEI(J).EQ.12)
     &          .OR.(LBCODEI(J).EQ.22).OR.(LBCODEI(J).EQ.32)) THEN
               READ(20+CK_LUN,*) QNIN1TMP
            END IF
          END DO
           
          RETURN   
        END SUBROUTINE READFORT20SNAP

        SUBROUTINE FFWDFORT20( NTTIME, CSTIME, StartLocTime )
          IMPLICIT NONE

          INTEGER:: NTTIME
          REAL (SZ):: StartLocTime
          INTERFACE 
            FUNCTION CSTIME( IT, StartTimeLoc ) RESULT( TIMEIT )
              IMPLICIT NONE

              INTEGER:: IT
              REAL (8):: StartTimeLoc
              REAL (8):: TIMEIT

            END FUNCTION
          END INTERFACE 

          INTEGER:: IT 
          REAL (SZ):: TIMEIT, QT1, QT2, FTIMINC_TMP
                           
          READ(20+CK_LUN,*) FTIMINC_TMP
          QT1 = StartLocTime
          QT2 = QT1 + FTIMINC_TMP
          CALL READFORT20SNAP( )  
          CALL READFORT20SNAP( )  

          DO IT=1, NTTIME
            ! TIMEIT=IT*DTDPHS + STATIM*86400.D0 ! kmd48.33bc - changed DTDP to DTDPHS
            TIMEIT = CSTIME( IT, StartTimeLoc )   
            IF(TIMEIT .GT. QT2) THEN
              QT1 = QT2
              QT2 = QT2 + FTIMINC_TMP
              CALL READFORT20SNAP( ) 
            ENDIF
          END DO
 
          RETURN  
        END SUBROUTINE FFWDFORT20

      END SUBROUTINE Reopen_elevbc_riverbc_InputFile

!!------------------------------------------------------------------------------------
!!-------------------------------------------------------------------------------------

      SUBROUTINE close_time_varying_bathymetry_file( )
        IMPLICIT NONE
  
        INTEGER:: FUNIT
        ! Time varying bathmetry
        IF (NDDT.NE.0) then
           FUNIT = 141 + CK_LUN   
           CALL CLOSEFUNIT( FUNIT  )  
        ENDIF  !end test for time varying bathymetry
C
        RETURN  
      END SUBROUTINE close_time_varying_bathymetry_file
 
!!-------------------------------------------------------------------------------------
      SUBROUTINE Reopen_time_varying_bathymetry( )
        use sizes, only : inputdir 
        use mesh, only :  np
        use global, only: BTIMINC, BCHGTIMINC 

        IMPLICIT NONE

        INTEGER:: NTTIME
        REAL (SZ):: StartTimeLoc 
               
C     time varying bathymetry values are read in at all grid nodes
C     at a time interval that does not equal the model time
C     step. Interpolation in time is used to synchronize the 
C     bathymetry information with the model time step.
        IF(Nddt.EQ.1) THEN
          OPEN(141+CK_LUN,FILE=TRIM(INPUTDIR)//'/'//'fort.141')

          NTTIME = itime  
          SELECT CASE( IHOT )
          CASE (0)
            StartTimeLoc = STATIM*86400.D0  
            CALL FFWDFORT141( NTTIME, CSTIME_CS0, StartTimeLoc )
          CASE (17,67,68,367,368,567,568)
            StartTimeLoc = STATIM*86400.D0
            CALL FFWDFORT141( NTTIME, CSTIME_CS1, StartTimeLoc )
          END SELECT
 
        ENDIF

        IF( NDDT .EQ. -1 ) THEN
          OPEN(141+CK_LUN,FILE=TRIM(INPUTDIR)//'/'//'fort.141')
        
          SELECT CASE( IHOT )
          CASE (0)
            NTTIME = itime  
            StartTimeLoc = STATIM*86400.D0  
            CALL FFWDFORT141( NTTIME, CSTIME_CS0, StartTimeLoc )
          CASE (17,67,68,367,368,567,568)
            NTTIME = itime - ITHS 

            StartTimeLoc = ITHS*DTDPHS + STATIM*86400.D0
            CALL FFWDFORT141( NTTIME, CSTIME_HS, StartTimeLoc )
          END SELECT                          
        ENDIF

        IF(NDDT.EQ.2) THEN
          OPEN(141+CK_LUN,FILE=TRIM(INPUTDIR)//'/'//'fort.141')

          NTTIME = itime  
          SELECT CASE( IHOT )
          CASE (0)
            StartTimeLoc = STATIM*86400.D0  
            CALL FFWDFORT141NDDT2( NTTIME, CSTIME_CS0, StartTimeLoc )
          CASE (17,67,68,367,368,567,568)
            StartTimeLoc = STATIM*86400.D0
            CALL FFWDFORT141NDDT2( NTTIME, CSTIME_CS1, StartTimeLoc )
          END SELECT
        END IF

 
        IF(NDDT.EQ.-2) THEN
          OPEN(141+CK_LUN,FILE=TRIM(INPUTDIR)//'/'//'fort.141')

          SELECT CASE( IHOT )
          CASE (0)
            NTTIME = itime  
            StartTimeLoc = STATIM*86400.D0  
            CALL FFWDFORT141NDDT2( NTTIME, CSTIME_CS0, StartTimeLoc )
          CASE (17,67,68,367,368,567,568)
            NTTIME = itime - ITHS  
            StartTimeLoc = ITHS*DTDPHS + STATIM*86400.D0
            CALL FFWDFORT141NDDT2( NTTIME, CSTIME_HS, StartTimeLoc )
          END SELECT
        ENDIF

      CONTAINS

         SUBROUTINE READFORT141SNAP( )
           IMPLICIT NONE

           INTEGER:: I, IDUMMY
           REAL(SZ):: RDUMMY
 
           DO I = 1, NP
              READ(141+CK_LUN,*) IDUMMY, RDUMMY
           END DO

           RETURN
         END SUBROUTINE READFORT141SNAP 
C
         SUBROUTINE FFWDFORT141( NTTIME, CSTIME, StartLocTime )
          IMPLICIT NONE

          INTEGER:: NTTIME
          REAL (SZ):: StartLocTime
          INTERFACE 
            FUNCTION CSTIME( II, StartTimeLoc ) RESULT( TIMEIT )
              IMPLICIT NONE

              INTEGER:: II
              REAL (8):: StartTimeLoc

              REAL (8):: TIMEIT
            END FUNCTION CSTIME
          END INTERFACE

          INTEGER:: IT
          REAL (SZ):: bT1, bT2, bT_END, TIMEIT

!          bT1 = STATIM*86400.D0
          bT1 = StartLocTime 
          bT2 = bT1 + bTIMINC
          BT_END = BT1 + BCHGTIMINC  !ending time for bathymetry changes during the btiminc interval
          CALL READFORT141SNAP( )  

C         DO IT=1,ITHS
          DO IT = 1, NTTIME
C            TIMEIT=IT*DTDPHS + STATIM*86400.D0 ! kmd48.33bc - changed the DTDP to DTDPHS
            TIMEIT = CSTIME( IT, StartTimeLoc )  
            IF(TIMEIT .GT. bT2) THEN
              bT1 = bT2  
              bT2 = bT2+bTIMINC
              BT_END = BT1 + BCHGTIMINC
              CALL READFORT141SNAP( )
            ENDIF  !timeit
          END DO  !IT

          RETURN  
        END SUBROUTINE FFWDFORT141 

        SUBROUTINE FFWDFORT141NDDT2( NTTIME, CSTIME, StartLocTime )
          IMPLICIT NONE   
        
          INTEGER:: NTTIME
          REAL (SZ):: StartLocTime
          INTERFACE 
            FUNCTION CSTIME( II, StartTimeLoc ) RESULT( TIMEIT )
              IMPLICIT NONE

              INTEGER:: II
              REAL (8):: StartTimeLoc

              REAL (8):: TIMEIT
            END FUNCTION CSTIME
          END INTERFACE


          INTEGER:: IT
          REAL (SZ):: bT1, bT2, BT_END, TIMEIT
          REAL (SZ), ALLOCATABLE:: DPTMP(:)

!           bT1 = STATIM*86400.D0
          bT1 = StartLocTime  
          bT2 = bT1 + bTIMINC 
          BT_END = BT1 + BCHGTIMINC  !ending time for bathymetry changes during the btiminc interval
          
          ALLOCATE( DPTMP(NP) ) 
           
!!!     go get first record for only some nodes, all
!!!     other nodes keep their current value
          CALL NDDT2GET( 141 + CK_LUN, DPTMP, -99999.d0 )

C          DO IT=1,ITHS
          DO IT = 1, NTTIME
            ! TIMEIT=IT*DTDPHS + STATIM*86400.D0 ! kmd48.33bc - changed the DTDP to DTDPHS
            TIMEIT = CSTIME( IT, StartTimeLoc ) 
            IF(TIMEIT .GT. bT2) THEN
              bT1 = bT2
              bT2 = bT2 + bTIMINC
              BT_END = BT1 + BCHGTIMINC
              CALL NDDT2GET( 141 + CK_LUN, DPTMP,-99999.d0 )
            ENDIF
          END DO  !IT
          

          DEALLOCATE( DPTMP )  
 
          RETURN  
        END SUBROUTINE FFWDFORT141NDDT2

      END SUBROUTINE Reopen_time_varying_bathymetry

C--------------------------------------------------------------------------
C------------------ Ice
C--------------------------------------------------------------------------
        SUBROUTINE close_ice_files( )
          IMPLICIT NONE

          INTEGER:: FUNIT 

          IF ( NCICE .EQ. 12 ) THEN
             FUNIT = 225 + CK_LUN      
             CALL CLOSEFUNIT( FUNIT )  
             
             IF ( NumSetsIce .gt. 0 ) THEN
               FUNIT =  227 + CK_LUN  
               CALL CLOSEFUNIT( FUNIT )   
             END IF
          ENDIF

          RETURN 
        END SUBROUTINE close_ice_files

        SUBROUTINE reopenffwd_ice_inputfile( )
          use mesh, only : dp, ics, np

          IMPLICIT NONE

          INTEGER:: NTTIME
          REAL (SZ):: StartLocTime

          IF ( NCICE .EQ. 12 ) THEN
            !
            IF ( NWS .GT.  0 ) THEN
              NTTIME = ITIME  
              StartLocTime = STATIM*86400.0  
              SELECT CASE(IHOT)   
              CASE (0)
                  CALL FFWD_OWIICE( NTTIME, CSTIME_CS0, StartLocTime )    
              CASE (17,67,68,367,368,567,568) ! non-portable binary or netcdf hotstart file
                  CALL FFWD_OWIICE( NTTIME, CSTIME_CS1, StartLocTime )   
              END SELECT
            ENDIF
            !
            IF ( NWS .LT. 0 ) THEN
              SELECT CASE(IHOT)   
              CASE (0)
                NTTIME = ITIME   
                StartLocTime = STATIM*86400.0    
                CALL FFWD_OWIICE( NTTIME, CSTIME_CS0, StartLocTime )    
              CASE (17,67,68,367,368,567,568) ! non-portable binary or netcdf hotstart file
                NTTIME = (ITIME - ITHS)  
                StartLocTime = ITHS*DTDPHS +  STATIM*86400.0  
                CALL FFWD_OWIICE( NTTIME, CSTIME_HS, StartLocTime )   
              END SELECT
            END IF
            !
          END IF  

          RETURN
        CONTAINS
 
          SUBROUTINE FFWD_OWIICE( NTTIME, CSTIME, StartLocTime )
           IMPLICIT NONE

           INTEGER:: NTTIME
           REAL (SZ):: StartLocTime
           INTERFACE 
              FUNCTION CSTIME( IT, StartLocTime ) RESULT( TIMEIT )
               IMPLICIT NONE
                
               INTEGER:: IT
               REAL (8):: StartLocTime

               REAL (8):: TIMEIT
              END FUNCTION CSTIME 
           END INTERFACE
           
 
           INTEGER:: IT
           REAL (SZ):: CT1, CT2, TIMEIT
           REAL (SZ), ALLOCATABLE:: CICETMP(:)

           ALLOCATE( CICETMP(NP) )   

!          CT1 = STATIM*86400.D0
           CT1 = StartLocTime  
           CT2 = CT1 + CICE_TIMINC

           ! This just initializes some variables and might set CICE,etc...
           ! if skipping ahead is used
           CALL NCICE1_INIT(CICETMP, NP)  
           CALL NCICE1_GET(CICETMP, NP)  
           CALL NCICE1_GET(CICETMP, NP)  

           !  
           DO IT=1, NTTIME
!             TIMEIT=IT*DTDP + STATIM*86400.D0
             TIMEIT = CSTIME( IT, StartLocTime )  
             IF( TIMEIT .GT. CT2 ) THEN
                  CT1 = CT2
                  CT2 = CT2 + CICE_TIMINC
                  CALL NCICE1_GET( CICETMP, NP)
             ENDIF
           END DO

           DEALLOCATE( CICETMP )  

           RETURN    
          END SUBROUTINE FFWD_OWIICE
             
 
        END SUBROUTINE reopenffwd_ice_inputfile
C------------------------------------------
C--------------------------------------------------------------------------
C------------------ Ice

C--------------------------------------------------------------------------
C--------------------------------------------------------------------------
C  wave radiation stress 
C--------------------------------------------------------------------------
C--------------------------------------------------------------------------


       SUBROUTINE close_rs_file()
         IMPLICIT NONE

          
         INTEGER:: FUNIT

         IF ( NRS .GE. 1 ) THEN
           SELECT CASE( NRS )
           CASE (1,2)
              FUNIT = 23 + CK_LUN   
              CALL CLOSEFUNIT( FUNIT )  
           END SELECT
         ENDIF 

         RETURN 
       END SUBROUTINE close_rs_file 
        
       SUBROUTINE reopenffwd_rs_file( )
         use mesh, only: dp, ics, np
         use wind, only: rsget 
 
         IMPLICIT NONE
   
         INTEGER:: NTTIME
         REAL (SZ):: StartLocTime  
   
         ! jgf51.46: Fixed initialization issue due to mishandling the 
         ! cases where NRS.eq.3 (SWAN coupling) or NRS.eq.4 (STWAVE
         ! coupling). The code was adding the wave radiation stresses
         ! to wind stresses whenever NRS.gt.0, but in the case of 
         ! NRS.eq.3 or NRS.eq.4, the wave radiation stresses were 
         ! uninitialized. This resulted in NaNs getting into the solution
         ! upon hotstarting. I fixed this by initializing RSNX1 and
         ! RSNY1 to 0.0 immediately after the arrays are allocated. 
         IF(NRS.GE.1) THEN ! sb46.28sb03
            SELECT CASE(IHOT)
            CASE (0)
                NTTIME = ITIME  
                StartLocTime = STATIM*86400.0  
                CALL FFWD_RSFILE( NTTIME, CSTIME_CS0, StartLocTime )  
            CASE (17,67,68,367,368,567,568) ! non-portable binary or netcdf hotstart file
                NTTIME = ITIME - ITHS 
                StartLocTime = ITHS*DTDPHS + STATIM*86400.0 
                CALL FFWD_RSFILE( NTTIME, CSTIME_HS, StartLocTime )   
            END SELECT  
         ENDIF

       CONTAINS

           SUBROUTINE FFWD_RSFILE(  NTTIME, CSTIME, StartLocTime  )
             IMPLICIT NONE
  
             INTEGER:: NTTIME
             REAL (SZ):: RT1, RT2, StartLocTime
             INTERFACE 
                FUNCTION CSTIME( IT, StartLocTime ) RESULT( TIMEIT )
                  IMPLICIT NONE

                  INTEGER:: IT
                  REAL (8):: StartLocTime, TIMEIT
                END FUNCTION CSTIME
             END INTERFACE
                      
             INTEGER:: IT
             REAL (SZ):: TIMEIT
             REAL (SZ), ALLOCATABLE:: RSNXTMP(:), RSNYTMP(:) 
           
             ALLOCATE( RSNXTMP(NP), RSNYTMP(NP) )  

!             RT1 = TimeLoc
             RT1 = StartLocTime  
             RT2 = RT1 + RSTIMINC
             SELECT CASE(NRS)
             CASE(1)           
                OPEN(23+CK_LUN,FILE=TRIM(INPUTDIR)//'/'//'fort.23')

                CALL RSGET( RSNXTMP, RSNYTMP )
                CALL RSGET( RSNXTMP, RSNYTMP )
                DO IT = 1, NTTIME
                   TIMEIT = CSTIME( IT, StartLocTime ) 
                   IF ( TIMEIT .GT. RT2 ) THEN 
                      RT1 = RT2  
                      RT2 = RT1 + RSTIMINC  
                      CALL RSGET( RSNXTMP, RSNYTMP )  
                   ENDIF
                END DO    
             CASE(2)
                ! NRS=2 was added. sb46.28sb03 09/xx/2006
                CALL RS2INIT( RSNXTMP ,RSNYTMP, NP)
                CALL RS2GET( RSNXTMP, RSNYTMP, NP)
                CALL RS2GET( RSNXTMP, RSNYTMP, NP)

                DO IT = 1, NTTIME
                   TIMEIT = CSTIME( IT, StartLocTime ) 
                   IF ( TIMEIT .GT. RT2 ) THEN 
                      RT1 = RT2  
                      RT2 = RT1 + RSTIMINC  
                      CALL RS2GET( RSNXTMP, RSNYTMP, NP )  
                   ENDIF 
                END DO 
             CASE(3,4)
               ! do nothing, initialization takes place in timestep.F
             CASE DEFAULT
!               write(scratchMessage,'(a,i0,a)') 
!     &         'The NRS parameter cannot have the value ',NRS,'.' 
!               call allMessage(ERROR,scratchMessage)           
             END SELECT
             !  
            
             RETURN 
           END SUBROUTINE FFWD_RSFILE
 
        
       END SUBROUTINE reopenffwd_rs_file



      END MODULE AMPIIOSUBS
     
