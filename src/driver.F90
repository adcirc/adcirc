!******************************************************************************
!******************************************************************************
!*    Non-ESMF ADCIRC DRIVER
!******************************************************************************
!******************************************************************************
    PROGRAM ADCIRC

    USE ADCIRC_Mod, ONLY : ADCIRC_Init, ADCIRC_Run, ADCIRC_Final

#ifdef CSWAN
! Casey 090302: Include the following routines for coupling to unstructured SWAN.
    USE Couple2swan, ONLY: PADCSWAN_INIT, PADCSWAN_FINAL
    USE SIZES, ONLY: MNPROC,MYPROC
#endif

    IMPLICIT NONE

    CALL ADCIRC_Init

#ifdef CSWAN
! Casey 090302: Allow SWAN to initialize stuff before the start
!             of the time step loop.  This subroutine is inside
!             the 'couple2swan.F' src file.
! Casey 110518: Added this IF statement.
    IF(MYPROC < MNPROC)THEN
        CALL PADCSWAN_INIT
    ENDIF
#endif

    CALL ADCIRC_Run

#ifdef CSWAN
! Casey 090302: Let SWAN clean up stuff.
! Casey 110518: Added this IF statement.
    IF(MYPROC < MNPROC)THEN
        CALL PADCSWAN_FINAL
    ENDIF
#endif

    CALL ADCIRC_Final

    STOP
    END PROGRAM ADCIRC
