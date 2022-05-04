#ifndef CFI_H
#define CFI_H

#if (defined(_CRAYMPP))
#  include <fortran.h>
#  ifndef CRAY_FOR_TYPE
#    define CRAY_FOR_TYPE 1
#  endif
#else
#  define BSD_FOR_TYPE 1
#endif

#ifndef PASTE2
#  define PASTE2(a, b) gEnErIcPaste2(a,b)
#  define gEnErIcPaste2(a, b) a##b
#endif

#if (defined(BSD_FOR_TYPE))
#  define CP2FCD(cp, len, fcd) do{fcd = cp; PASTE2(fcd,_len)=len;}while(0)
#  define FCD2CP(s) s
#  define FCDLEN(s) PASTE2(s,_len)
#  define FDESC(s) char * s; static int PASTE2(s,_len)
#  define FSTRING(a, comma) comma char * a
#  define protoFSTRING(a)  char * a
#  define protoLenFSTRING(a) , int PASTE2(a,_len)
#  if (defined(NO_TRAILING_UNDERSCORE) || defined(_IBMR2))
#    define FORTRAN_NAME(UN,ln) ln
#  else
#    define FORTRAN_NAME(UN, ln) PASTE2(ln,_)
#  endif
#endif

#if (defined(CRAY_FOR_TYPE))
#  define CP2FCD(cp,len,fcd) fcd = _cptofcd(cp,len)
#  define FCD2CP(s) _fcdtocp(s)
#  define FDESC(s) _fcd s
#  define FCDLEN(s) _fcdlen(s)
#  define FSTRING(a,comma) comma _fcd a
#  define protoFSTRING(a) _fcd a
#  define protoLenFSTRING(a)
#  define FORTRAN_NAME(UN,ln) UN
#endif


#endif /* CFI_H */
