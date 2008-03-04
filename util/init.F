C-----------------------------------------------------------------------
C     P R O G R A M   I N I T 
C-----------------------------------------------------------------------
C
C     Program to initialize netcdf output files for a parallel run.      
C
C     Written by Jason G. Fleming. 
C
C-----------------------------------------------------------------------

      PROGRAM INIT
      USE GLOBAL
      USE NETCDFIO

      IMPLICIT NONE
      INTEGER I
      REAL(8) TIMEHSF

      OPEN(99,FILE='fort.67',ACCESS='DIRECT',RECL=8)
      IHOTSTP=3
      READ(99,REC=IHOTSTP) TIMEHSF
      WRITE(*,*) TIMEHSF
      CLOSE(99)           

      END PROGRAM INIT
