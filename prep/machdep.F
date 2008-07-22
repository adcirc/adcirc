C
C--Machine dependent code for adcprep
C

Cifdef LINUX
      SUBROUTINE  IREAD(STRING,INTVAR)
      INTEGER INTVAR
      CHARACTER*4   STRING
c
c-----------------------------------------------------------
c   Wrapper for internal read on LINUX Systems        vjp97
c-----------------------------------------------------------
c
      READ(STRING,'(I4.4)') INTVAR
c
      RETURN
      END

      SUBROUTINE  IWRITE(STRING,I1,I2,INTVAR)
      INTEGER INTVAR,I1,I2
      CHARACTER   STRING*(*)
c
c-----------------------------------------------------------
c   Wrapper for internal write on LINUX Systems        vjp97
c-----------------------------------------------------------
c
      WRITE(STRING(I1:I2),'(I4.4)') INTVAR
c
      RETURN
      END

      SUBROUTINE  IWRITE2(STRING,I1,I2,INTVAR)
      INTEGER INTVAR,I1,I2
      CHARACTER   STRING*(*)
c
c-----------------------------------------------------------
c   Wrapper for internal write on LINUX Systems        vjp97
c-----------------------------------------------------------
c
      WRITE(STRING(I1:I2),'(I3.3)') INTVAR
c
      RETURN
      END

Cendif



