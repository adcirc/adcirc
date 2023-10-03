!-------------------------------------------------------------------------------!
!
! ADCIRC - The ADvanced CIRCulation model
! Copyright (C) 1994-2023 R.A. Luettich, Jr., J.J. Westerink
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
C
C--Machine dependent code for adcprep
C
C... TCM v50.85 -- Added WINDOWS
C#ifdef PC_DIG_FORT
#if defined(PC_DIG_FORT) || defined(WINDOWS)

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

#endif


#ifdef LINUX
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

#endif


#ifdef SGI
      SUBROUTINE  IREAD(STRING,INTVAR)
      INTEGER INTVAR
      CHARACTER*4   STRING
c
c-----------------------------------------------------------
c   Wrapper for internal read on SGI Systems           vjp97
c-----------------------------------------------------------
c
      READ(STRING,'(I4.4)') INTVAR
c
c     print *, "intvar ", INTVAR 
c
      RETURN
      END

      SUBROUTINE  IWRITE(STRING,I1,I2,INTVAR)
      INTEGER INTVAR,I1,I2
      CHARACTER   STRING*(*)
c
c-----------------------------------------------------------
c   Wrapper for internal write on SGI Systems          vjp97
c-----------------------------------------------------------
c
      WRITE(STRING(I1:I2),'(I4.4)') INTVAR
c
c     print *, "string= ", string
c     print *, "I1 = ", I1
c     print *, "I2 = ", I2
c     print *, "INTVAR = ", INTVAR
c
      RETURN
      END

      SUBROUTINE  IWRITE2(STRING,I1,I2,INTVAR)
      INTEGER INTVAR,I1,I2
      CHARACTER   STRING*(*)
c
c-----------------------------------------------------------
c   Wrapper for internal write on SGI Systems          vjp97
c-----------------------------------------------------------
c
      WRITE(STRING(I1:I2),'(I3.3)') INTVAR
c
c     print *, "string= ", string
c     print *, "I1 = ", I1
c     print *, "I2 = ", I2
c     print *, "INTVAR = ", INTVAR
c
      RETURN
      END
#endif
 

#ifdef IBM
      SUBROUTINE  IREAD(STRING,INTVAR)
      INTEGER INTVAR
      CHARACTER*4   STRING
c
c-----------------------------------------------------------
c   Wrapper for internal read on IBM Systems           vjp97
c-----------------------------------------------------------
c
      READ(STRING,'(I4.4)') INTVAR
c
c     print *, "intvar ", INTVAR 
c
      RETURN
      END

      SUBROUTINE  IWRITE(STRING,I1,I2,INTVAR)
      INTEGER INTVAR,I1,I2
      CHARACTER   STRING*(*)
c
c-----------------------------------------------------------
c   Wrapper for internal write on IBM Systems          vjp97
c-----------------------------------------------------------
c
      WRITE(STRING(I1:I2),'(I4.4)') INTVAR
c
c     print *, "string= ", string
c     print *, "I1 = ", I1
c     print *, "I2 = ", I2
c     print *, "INTVAR = ", INTVAR
c
      RETURN
      END

      SUBROUTINE  IWRITE2(STRING,I1,I2,INTVAR)
      INTEGER INTVAR,I1,I2
      CHARACTER   STRING*(*)
c
c-----------------------------------------------------------
c   Wrapper for internal write on IBM Systems          vjp97
c-----------------------------------------------------------
c
      WRITE(STRING(I1:I2),'(I3.3)') INTVAR
c
c     print *, "string= ", string
c     print *, "I1 = ", I1
c     print *, "I2 = ", I2
c     print *, "INTVAR = ", INTVAR
c
      RETURN
      END
#endif

#ifdef MACHSUN
      SUBROUTINE  IREAD(STRING,INTVAR)
      INTEGER INTVAR
      CHARACTER*4   STRING
c
c-----------------------------------------------------------
c   Wrapper for internal read on SUN Systems           vjp97
c-----------------------------------------------------------
c
      READ(STRING,'(I4.4)') INTVAR
c
c     print *, "intvar ", INTVAR 
c
      RETURN
      END

      SUBROUTINE  IWRITE(STRING,I1,I2,INTVAR)
      INTEGER INTVAR,I1,I2
      CHARACTER   STRING*(*)
c
c-----------------------------------------------------------
c   Wrapper for internal write on SUN Systems          vjp97
c-----------------------------------------------------------
c
      WRITE(STRING(I1:I2),'(I4.4)') INTVAR
c
c     print *, "string= ", string
c     print *, "I1 = ", I1
c     print *, "I2 = ", I2
c     print *, "INTVAR = ", INTVAR
c
      RETURN
      END

      SUBROUTINE  IWRITE2(STRING,I1,I2,INTVAR)
      INTEGER INTVAR,I1,I2
      CHARACTER   STRING*(*)
c
c-----------------------------------------------------------
c   Wrapper for internal write on SUN Systems          vjp97
c-----------------------------------------------------------
c
      WRITE(STRING(I1:I2),'(I3.3)') INTVAR
c
c     print *, "string= ", string
c     print *, "I1 = ", I1
c     print *, "I2 = ", I2
c     print *, "INTVAR = ", INTVAR
c
      RETURN
      END
#endif

#ifdef CRAY
      SUBROUTINE  GETARG(NUMARG,STRING)
      INTEGER NUMARG
      CHARACTER(*)   STRING
      INTEGER*8  NUMARG8,ILEN,IERR
c
c-----------------------------------------------------------
c   Wrapper for routine GETARG on all Cray Research Systems
c   vjp97
c-----------------------------------------------------------
c
      NUMARG8 = NUMARG
      CALL PXFGETARG(NUMARG8,STRING,ILEN,IERR)
c
c     print *, "from getarg: ilen = ",ilen
c     print *, "from getarg: ierr = ",ierr
c     print *, "string = ", string
c
      RETURN
      END

      SUBROUTINE  IREAD(STRING,INTVAR)
      INTEGER INTVAR
      CHARACTER  STRING*4
      INTEGER*8  INTVAR8
c
c-----------------------------------------------------------
c   Wrapper for internal read on Cray Research Systems
c   vjp97
c-----------------------------------------------------------
c
      READ(STRING,'(I4.4)') INTVAR8
      INTVAR = INTVAR8
c
c     print *, "intvar8 ", INTVAR 
c
      RETURN
      END

      SUBROUTINE  IWRITE(STRING,IA,IB,INTVAR)
      INTEGER INTVAR,IA,IB
      CHARACTER  STRING*72
      INTEGER*8  INTVAR8, IA8, IB8
c
c-----------------------------------------------------------
c   Wrapper for internal write on Cray Research Systems
c   vjp97
c-----------------------------------------------------------
c
      INTVAR8 = INTVAR
      IA8     = IA
      IB8     = IB
      WRITE(STRING(IA8:IB8),'(I4.4)') INTVAR8
c
c     print *, "string= ", string
c
      RETURN
      END


      SUBROUTINE  IWRITE2(STRING,IA,IB,INTVAR)
      INTEGER INTVAR,IA,IB
      CHARACTER  STRING*72
      INTEGER*8  INTVAR8, IA8, IB8
c
c-----------------------------------------------------------
c   Wrapper for internal write on Cray Research Systems
c   vjp97
c-----------------------------------------------------------
c
      INTVAR8 = INTVAR
      IA8     = IA
      IB8     = IB
      WRITE(STRING(IA8:IB8),'(I3.3)') INTVAR8
c
c     print *, "string= ", string
c
      RETURN
      END
#endif

c     jgf45.06 CRAYX1 code contributed by MEB 04/01/04

#ifdef CRAYX1
      SUBROUTINE  GETARG(NUMARG,STRING)
      INTEGER NUMARG
      CHARACTER(*)   STRING
      INTEGER*8  NUMARG8,ILEN,IERR
c
c-----------------------------------------------------------
c   Wrapper for routine GETARG on all Cray Research Systems
c   vjp97
c-----------------------------------------------------------
c
      NUMARG8 = NUMARG
      CALL PXFGETARG(NUMARG8,STRING,ILEN,IERR)
c
c     print *, "from getarg: ilen = ",ilen
c     print *, "from getarg: ierr = ",ierr
c     print *, "string = ", string
c
      RETURN
      END

      SUBROUTINE  IREAD(STRING,INTVAR)
      INTEGER INTVAR
      CHARACTER  STRING*4
c
c-----------------------------------------------------------
c   Wrapper for internal read on Cray Research Systems
c   vjp97
c-----------------------------------------------------------
c
      READ(STRING,'(I4.4)') INTVAR
c
c     print *, "intvar8 ", INTVAR
c
      RETURN
      END

      SUBROUTINE  IWRITE(STRING,IA,IB,INTVAR)
      INTEGER INTVAR,IA,IB
      CHARACTER  STRING*72
c
c-----------------------------------------------------------
c   Wrapper for internal write on Cray Research Systems
c   vjp97
c-----------------------------------------------------------
c
      WRITE(STRING(IA:IB),'(I4.4)') INTVAR
c
c     print *, "string= ", string
c
      RETURN
      END


      SUBROUTINE  IWRITE2(STRING,IA,IB,INTVAR)
      INTEGER INTVAR,IA,IB
      CHARACTER  STRING*72
c
c-----------------------------------------------------------
c   Wrapper for internal write on Cray Research Systems
c   vjp97
c-----------------------------------------------------------
c
      WRITE(STRING(IA:IB),'(I3.3)') INTVAR
c
c     print *, "string= ", string
c
      RETURN
      END

#endif
