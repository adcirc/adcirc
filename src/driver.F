C******************************************************************************
C******************************************************************************
C*    Non-ESMF ADCIRC DRIVER
C******************************************************************************
C******************************************************************************
      PROGRAM ADCIRC
C
      USE ADCIRC_Mod, ONLY : ADCIRC_Init, ADCIRC_Run, ADCIRC_Final
C
      IMPLICIT NONE
C
      CALL ADCIRC_Init
C
      CALL ADCIRC_Run
C
      CALL ADCIRC_Final
C
      STOP
      END PROGRAM ADCIRC
