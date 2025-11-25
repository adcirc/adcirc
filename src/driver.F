!-------------------------------------------------------------------------------!
!
! ADCIRC - The ADvanced CIRCulation model
! Copyright (C) 1994-2025 R.A. Luettich, Jr., J.J. Westerink
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
!******************************************************************************
!******************************************************************************
!*    Non-ESMF ADCIRC DRIVER
!******************************************************************************
!******************************************************************************
      PROGRAM ADCIRC

      USE ADCIRC_Mod, ONLY : ADCIRC_Init, ADCIRC_Run, ADCIRC_Final
      use mod_terminate, only: terminate, ADCIRC_EXIT_SUCCESS

#ifdef CSWAN
! Casey 090302: Include the following routines for coupling to unstructured SWAN.
      USE Couple2swan, ONLY: PADCSWAN_INIT, PADCSWAN_FINAL
      USE SIZES, ONLY: MNPROC,MYPROC
#endif

      IMPLICIT NONE

      CALL ADCIRC_Init()

#ifdef CSWAN
! Casey 090302: Allow SWAN to initialize stuff before the start
!             of the time step loop.  This subroutine is inside
!             the 'couple2swan.F' src file.
! Casey 110518: Added this IF statement.
      IF(MYPROC.LT.MNPROC)THEN
         CALL PADCSWAN_INIT()
      ENDIF
#endif

      CALL ADCIRC_Run()

#ifdef CSWAN
! Casey 090302: Let SWAN clean up stuff.
! Casey 110518: Added this IF statement.
      IF(MYPROC.LT.MNPROC)THEN
         CALL PADCSWAN_FINAL()
      ENDIF
#endif

      CALL ADCIRC_Final()

      call terminate(exit_code=ADCIRC_EXIT_SUCCESS)
      END PROGRAM ADCIRC
