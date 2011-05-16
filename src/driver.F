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
Casey 110516: Added this IF statement to prevent the writer cores
C             from entering SWAN.
      IF(MYPROC.LT.MNPROC)THEN
         CALL PADCSWAN_INIT
      ENDIF
#endif
C
      CALL ADCIRC_Run
C
#ifdef CSWAN
Casey 090302: Let SWAN clean up stuff.
Casey 110516: Added this IF statement to prevent the writer cores
C             from entering SWAN.
      IF(MYPROC.LT.MNPROC)THEN
         CALL PADCSWAN_FINAL
      ENDIF
#endif
C
      CALL ADCIRC_Final
C
      STOP
      END PROGRAM ADCIRC
