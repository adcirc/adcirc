C***********************************************************************
C AMPI virtualization module for LUNs 12/30/2019
C                                                                      *
C***********************************************************************

      module CkLUNMod
      implicit none
C     obtain AMPI rank
C     multiply by LUN factor, default 100000
      integer CK_LUN
!$omp threadprivate(CK_LUN)      
      integer :: CkLUNFactor=100000

      contains 

C-----------------------------------------------------------------------
      subroutine setAMPIzedLUN(thisRank)
      implicit none
      integer, intent(in) :: thisRank

      CK_LUN = thisRank * CkLUNFactor

      RETURN
      
      end subroutine setAMPIzedLUN
C-----------------------------------------------------------------------

      end module CkLUNMod
