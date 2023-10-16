C******************************************************************************
C******************************************************************************
C*    Non-ESMF ADCIRC DRIVER
C******************************************************************************
C******************************************************************************
#ifndef AMPI_COMP
      PROGRAM ADCIRC
#else
      SUBROUTINE MPI_MAIN
#endif      
C
      USE ADCIRC_Mod, ONLY : ADCIRC_Init, ADCIRC_Run, ADCIRC_Final
      USE AMPI_MIGRATE_IMPL, ONLY: ABOUT_TO_MIGRATE, JUST_MIGRATED,
     &  migratablelun_setup 
C
#ifdef CSWAN
Casey 090302: Include the following routines for coupling to unstructured SWAN.
      USE Couple2swan, ONLY: PADCSWAN_INIT, PADCSWAN_FINAL
      USE SIZES, ONLY: MNPROC,MYPROC
#endif
C
      IMPLICIT NONE
C
      CALL ADCIRC_Init
C
#ifdef CSWAN
Casey 090302: Allow SWAN to initialize stuff before the start
C             of the time step loop.  This subroutine is inside
C             the 'couple2swan.F' src file.
Casey 110518: Added this IF statement.
      IF(MYPROC.LT.MNPROC)THEN
         CALL PADCSWAN_INIT
      ENDIF
#endif
C
      CALL ADCIRC_Run
C
#ifdef CSWAN
Casey 090302: Let SWAN clean up stuff.
Casey 110518: Added this IF statement.
      IF(MYPROC.LT.MNPROC)THEN
         CALL PADCSWAN_FINAL
      ENDIF
#endif
C
      CALL ADCIRC_Final
C
      CALL EXIT(0)
#ifndef AMPI_COMP      
      END PROGRAM ADCIRC
#else
      END SUBROUTINE MPI_MAIN
#endif      
