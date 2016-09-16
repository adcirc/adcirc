
!--Machine dependent code for adcprep

!... TCM v50.85 -- Added WINDOWS
! ifdef PC_DIG_FORT
#if defined(PC_DIG_FORT) || defined(WINDOWS)

    SUBROUTINE  IREAD(STRING,INTVAR)
    INTEGER :: INTVAR
    CHARACTER(4) ::   STRING

!-----------------------------------------------------------
!   Wrapper for internal read on LINUX Systems        vjp97
!-----------------------------------------------------------

    READ(STRING,'(I4.4)') INTVAR

    RETURN
    END SUBROUTINE 

    SUBROUTINE  IWRITE(STRING,I1,I2,INTVAR)
    INTEGER :: INTVAR,I1,I2
    CHARACTER   STRING*(*)

!-----------------------------------------------------------
!   Wrapper for internal write on LINUX Systems        vjp97
!-----------------------------------------------------------

    WRITE(STRING(I1:I2),'(I4.4)') INTVAR

    RETURN
    END SUBROUTINE 

    SUBROUTINE  IWRITE2(STRING,I1,I2,INTVAR)
    INTEGER :: INTVAR,I1,I2
    CHARACTER   STRING*(*)

!-----------------------------------------------------------
!   Wrapper for internal write on LINUX Systems        vjp97
!-----------------------------------------------------------

    WRITE(STRING(I1:I2),'(I3.3)') INTVAR

    RETURN
    END SUBROUTINE 

#endif


#ifdef LINUX
    SUBROUTINE  IREAD(STRING,INTVAR)
    INTEGER :: INTVAR
    CHARACTER(4) ::   STRING

!-----------------------------------------------------------
!   Wrapper for internal read on LINUX Systems        vjp97
!-----------------------------------------------------------

    READ(STRING,'(I4.4)') INTVAR

    RETURN
    END SUBROUTINE 

    SUBROUTINE  IWRITE(STRING,I1,I2,INTVAR)
    INTEGER :: INTVAR,I1,I2
    CHARACTER   STRING*(*)

!-----------------------------------------------------------
!   Wrapper for internal write on LINUX Systems        vjp97
!-----------------------------------------------------------

    WRITE(STRING(I1:I2),'(I4.4)') INTVAR

    RETURN
    END SUBROUTINE 

    SUBROUTINE  IWRITE2(STRING,I1,I2,INTVAR)
    INTEGER :: INTVAR,I1,I2
    CHARACTER   STRING*(*)

!-----------------------------------------------------------
!   Wrapper for internal write on LINUX Systems        vjp97
!-----------------------------------------------------------

    WRITE(STRING(I1:I2),'(I3.3)') INTVAR

    RETURN
    END SUBROUTINE 

#endif


#ifdef SGI
    SUBROUTINE  IREAD(STRING,INTVAR)
    INTEGER :: INTVAR
    CHARACTER(4) ::   STRING

!-----------------------------------------------------------
!   Wrapper for internal read on SGI Systems           vjp97
!-----------------------------------------------------------

    READ(STRING,'(I4.4)') INTVAR

!     print *, "intvar ", INTVAR

    RETURN
    END SUBROUTINE 

    SUBROUTINE  IWRITE(STRING,I1,I2,INTVAR)
    INTEGER :: INTVAR,I1,I2
    CHARACTER   STRING*(*)

!-----------------------------------------------------------
!   Wrapper for internal write on SGI Systems          vjp97
!-----------------------------------------------------------

    WRITE(STRING(I1:I2),'(I4.4)') INTVAR

!     print *, "string= ", string
!     print *, "I1 = ", I1
!     print *, "I2 = ", I2
!     print *, "INTVAR = ", INTVAR

    RETURN
    END SUBROUTINE 

    SUBROUTINE  IWRITE2(STRING,I1,I2,INTVAR)
    INTEGER :: INTVAR,I1,I2
    CHARACTER   STRING*(*)

!-----------------------------------------------------------
!   Wrapper for internal write on SGI Systems          vjp97
!-----------------------------------------------------------

    WRITE(STRING(I1:I2),'(I3.3)') INTVAR

!     print *, "string= ", string
!     print *, "I1 = ", I1
!     print *, "I2 = ", I2
!     print *, "INTVAR = ", INTVAR

    RETURN
    END SUBROUTINE 
#endif
     

#ifdef IBM
    SUBROUTINE  IREAD(STRING,INTVAR)
    INTEGER :: INTVAR
    CHARACTER(4) ::   STRING

!-----------------------------------------------------------
!   Wrapper for internal read on IBM Systems           vjp97
!-----------------------------------------------------------

    READ(STRING,'(I4.4)') INTVAR

!     print *, "intvar ", INTVAR

    RETURN
    END SUBROUTINE 

    SUBROUTINE  IWRITE(STRING,I1,I2,INTVAR)
    INTEGER :: INTVAR,I1,I2
    CHARACTER   STRING*(*)

!-----------------------------------------------------------
!   Wrapper for internal write on IBM Systems          vjp97
!-----------------------------------------------------------

    WRITE(STRING(I1:I2),'(I4.4)') INTVAR

!     print *, "string= ", string
!     print *, "I1 = ", I1
!     print *, "I2 = ", I2
!     print *, "INTVAR = ", INTVAR

    RETURN
    END SUBROUTINE 

    SUBROUTINE  IWRITE2(STRING,I1,I2,INTVAR)
    INTEGER :: INTVAR,I1,I2
    CHARACTER   STRING*(*)

!-----------------------------------------------------------
!   Wrapper for internal write on IBM Systems          vjp97
!-----------------------------------------------------------

    WRITE(STRING(I1:I2),'(I3.3)') INTVAR

!     print *, "string= ", string
!     print *, "I1 = ", I1
!     print *, "I2 = ", I2
!     print *, "INTVAR = ", INTVAR

    RETURN
    END SUBROUTINE 
#endif

#ifdef MACHSUN
    SUBROUTINE  IREAD(STRING,INTVAR)
    INTEGER :: INTVAR
    CHARACTER(4) ::   STRING

!-----------------------------------------------------------
!   Wrapper for internal read on SUN Systems           vjp97
!-----------------------------------------------------------

    READ(STRING,'(I4.4)') INTVAR

!     print *, "intvar ", INTVAR

    RETURN
    END SUBROUTINE 

    SUBROUTINE  IWRITE(STRING,I1,I2,INTVAR)
    INTEGER :: INTVAR,I1,I2
    CHARACTER   STRING*(*)

!-----------------------------------------------------------
!   Wrapper for internal write on SUN Systems          vjp97
!-----------------------------------------------------------

    WRITE(STRING(I1:I2),'(I4.4)') INTVAR

!     print *, "string= ", string
!     print *, "I1 = ", I1
!     print *, "I2 = ", I2
!     print *, "INTVAR = ", INTVAR

    RETURN
    END SUBROUTINE 

    SUBROUTINE  IWRITE2(STRING,I1,I2,INTVAR)
    INTEGER :: INTVAR,I1,I2
    CHARACTER   STRING*(*)

!-----------------------------------------------------------
!   Wrapper for internal write on SUN Systems          vjp97
!-----------------------------------------------------------

    WRITE(STRING(I1:I2),'(I3.3)') INTVAR

!     print *, "string= ", string
!     print *, "I1 = ", I1
!     print *, "I2 = ", I2
!     print *, "INTVAR = ", INTVAR

    RETURN
    END SUBROUTINE 
#endif

#ifdef CRAY
    SUBROUTINE  GETARG(NUMARG,STRING)
    INTEGER :: NUMARG
    CHARACTER(*)   STRING
    INTEGER*8 ::  NUMARG8,ILEN,IERR

!-----------------------------------------------------------
!   Wrapper for routine GETARG on all Cray Research Systems
!   vjp97
!-----------------------------------------------------------

    NUMARG8 = NUMARG
    CALL PXFGETARG(NUMARG8,STRING,ILEN,IERR)

!     print *, "from getarg: ilen = ",ilen
!     print *, "from getarg: ierr = ",ierr
!     print *, "string = ", string

    RETURN
    END SUBROUTINE 

    SUBROUTINE  IREAD(STRING,INTVAR)
    INTEGER :: INTVAR
    CHARACTER  STRING*4
    INTEGER*8 ::  INTVAR8

!-----------------------------------------------------------
!   Wrapper for internal read on Cray Research Systems
!   vjp97
!-----------------------------------------------------------

    READ(STRING,'(I4.4)') INTVAR8
    INTVAR = INTVAR8

!     print *, "intvar8 ", INTVAR

    RETURN
    END SUBROUTINE 

    SUBROUTINE  IWRITE(STRING,IA,IB,INTVAR)
    INTEGER :: INTVAR,IA,IB
    CHARACTER  STRING*72
    INTEGER*8 ::  INTVAR8, IA8, IB8

!-----------------------------------------------------------
!   Wrapper for internal write on Cray Research Systems
!   vjp97
!-----------------------------------------------------------

    INTVAR8 = INTVAR
    IA8     = IA
    IB8     = IB
    WRITE(STRING(IA8:IB8),'(I4.4)') INTVAR8

!     print *, "string= ", string

    RETURN
    END SUBROUTINE 


    SUBROUTINE  IWRITE2(STRING,IA,IB,INTVAR)
    INTEGER :: INTVAR,IA,IB
    CHARACTER  STRING*72
    INTEGER*8 ::  INTVAR8, IA8, IB8

!-----------------------------------------------------------
!   Wrapper for internal write on Cray Research Systems
!   vjp97
!-----------------------------------------------------------

    INTVAR8 = INTVAR
    IA8     = IA
    IB8     = IB
    WRITE(STRING(IA8:IB8),'(I3.3)') INTVAR8

!     print *, "string= ", string

    RETURN
    END SUBROUTINE 
#endif

!     jgf45.06 CRAYX1 code contributed by MEB 04/01/04

#ifdef CRAYX1
    SUBROUTINE  GETARG(NUMARG,STRING)
    INTEGER :: NUMARG
    CHARACTER(*)   STRING
    INTEGER*8 ::  NUMARG8,ILEN,IERR

!-----------------------------------------------------------
!   Wrapper for routine GETARG on all Cray Research Systems
!   vjp97
!-----------------------------------------------------------

    NUMARG8 = NUMARG
    CALL PXFGETARG(NUMARG8,STRING,ILEN,IERR)

!     print *, "from getarg: ilen = ",ilen
!     print *, "from getarg: ierr = ",ierr
!     print *, "string = ", string

    RETURN
    END SUBROUTINE 

    SUBROUTINE  IREAD(STRING,INTVAR)
    INTEGER :: INTVAR
    CHARACTER  STRING*4

!-----------------------------------------------------------
!   Wrapper for internal read on Cray Research Systems
!   vjp97
!-----------------------------------------------------------

    READ(STRING,'(I4.4)') INTVAR

!     print *, "intvar8 ", INTVAR

    RETURN
    END SUBROUTINE 

    SUBROUTINE  IWRITE(STRING,IA,IB,INTVAR)
    INTEGER :: INTVAR,IA,IB
    CHARACTER  STRING*72

!-----------------------------------------------------------
!   Wrapper for internal write on Cray Research Systems
!   vjp97
!-----------------------------------------------------------

    WRITE(STRING(IA:IB),'(I4.4)') INTVAR

!     print *, "string= ", string

    RETURN
    END SUBROUTINE 


    SUBROUTINE  IWRITE2(STRING,IA,IB,INTVAR)
    INTEGER :: INTVAR,IA,IB
    CHARACTER  STRING*72

!-----------------------------------------------------------
!   Wrapper for internal write on Cray Research Systems
!   vjp97
!-----------------------------------------------------------

    WRITE(STRING(IA:IB),'(I3.3)') INTVAR

!     print *, "string= ", string

    RETURN
    END SUBROUTINE 

#endif
