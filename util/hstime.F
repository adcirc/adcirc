C***************************************************************************
C  PROGRAM TO READ A BINARY HOTSTART FILE AND PRINT TIME TO STDOUT         *
C                                                                          *
C  Written by Jason G. Fleming 12 Oct 2006                                 *
C***************************************************************************

      PROGRAM hstime
      IMPLICIT NONE
      INTEGER I,IHOTSTP
      REAL(8) TIMEHSF

      OPEN(99,FILE='fort.67',ACCESS='DIRECT',RECL=8)
      IHOTSTP=3
      READ(99,REC=IHOTSTP) TIMEHSF
      WRITE(*,*) TIMEHSF
      CLOSE(99)           

      END
