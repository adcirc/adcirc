      !=================================================================
      !=================================================================
      !=================================================================
      !      =====                                           =====
      !      =====              MODULE precision             =====
      !      =====                                           =====
      !=================================================================
      !=================================================================
      !=================================================================

      !=================================================================
      ! Define custom precision types for cross-platform portability.
      !
      ! Compilation:
      !    Throw C preprocessor macro switches (DEFINES) to select
      !    precision:
      !       f90 -O -c precision.F -DREAL8 -DINT2
      !       xlf90 -O -c precision.F -WF,-DREAL8,-DINT2
      !
      ! Revision history:
      !    Date         Programmer                 Description of change
      !    ----         ----------                 ---------------------
      !    05/23/06     Craig  Mattocks, UNC-CEP   Wrote original code
      !=================================================================
      MODULE precision

         !--------------------------------------------
         ! Force explicit declaration of all variables
         !--------------------------------------------
         IMPLICIT NONE

         !--------------------------------------------------------------
         ! Retain values of class instance variables between invocations
         !--------------------------------------------------------------
         SAVE

         !----------------------------------------
         ! Define a KIND with 6 digits of accuracy
         ! and exponent range of 10^-37 to 10^37
         !----------------------------------------
         INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(p=6, r=37)

         !-----------------------------------------
         ! Define a KIND with 15 digits of accuracy
         ! and exponent range of 10^-307 to 10^307
         !-----------------------------------------
         INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(p=15, r=307)

         !--------------------------------------------------
         ! Set number of bytes "sz" in REAL(8) declarations
         !--------------------------------------------------
#ifdef REAL8
!        INTEGER, PARAMETER :: sz = 8
         INTEGER, PARAMETER :: sz = dp
#else
!        INTEGER, PARAMETER :: sz = 4
         INTEGER, PARAMETER :: sz = sp
#endif

         !-----------------------------------------------------
         ! Set number of bytes "iz" in INTEGER(iz) declarations
         !-----------------------------------------------------
#ifdef INT8
         INTEGER, PARAMETER :: iz = 8
#elif INT2
         INTEGER, PARAMETER :: iz = 2
#elif INT1
         INTEGER, PARAMETER :: iz = 1
#else
         INTEGER, PARAMETER :: iz = 4
#endif

      END MODULE precision
