C******************************************************************************
C******************************************************************************
C*    Non-ESMF ADCIRC DRIVER
C******************************************************************************
C******************************************************************************
      PROGRAM ADCIRC
C
      USE ADCIRC_Mod, ONLY : ADCIRC_Init, ADCIRC_Run, ADCIRC_Final
C
#ifdef CSWAN
Casey 090302: We need the following routines.
      USE Couple2swan, ONLY: PADCSWAN_INIT, PADCSWAN_FINAL
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
      CALL PADCSWAN_INIT
#endif
C
      CALL ADCIRC_Run
C
#ifdef CSWAN
Casey 090302: Let SWAN clean up stuff.
      CALL PADCSWAN_FINAL
#endif
C
      CALL ADCIRC_Final
C
      STOP
      END PROGRAM ADCIRC
