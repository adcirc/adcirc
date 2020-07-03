      ! Sobroutine for 
      !   1. Closing relevant input files before conducting migration
      !   2. After mirgration, reopening and fast forward input files closed in in step 1 
      !      to the location where they are closed         
      MODULE AMPIIOSUBS
      use CkLunMod, only : CK_LUN
      USE SIZES 
      USE GLOBAL, ONLY: NWS

      USE OWIWIND, ONLY: owiwind_closefiles

      IMPLICIT NONE

      CONTAINS

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
          CALL owiwind_closefiles( ) ; 
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
           CALL owiwind_closefiles( ) ; 
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

         RETURN ;
        END SUBROUTINE close_met_files

      CONTAINS      

        ! close fort.22
        SUBROUTINE CLOSEFORT22( )
          IMPLICIT NONE

          LOGICAL:: iopened 
          iopened = .FALSE. ;

          INQUIRE( UNIT = 22 + CK_LUN, OPENED = iopened ) ;
          IF ( iopenned ) THEN 
            CLOSE( 22 + CK_LUN ) ;
          END IF

          RETURN ;
        END SUBROUTINE CLOSEF22 


        ! reopen a wind files and fast forward to the location 
        ! previously closed 
        SUBROUTINE REOPENFFWD22( ) 
          IMPLICIT NONE

          RETURN ; 
        END SUBROUTINE REOPENFFWD22   

      END MODULE AMPIIOSUBS
     
