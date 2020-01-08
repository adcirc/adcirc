C***********************************************************************
C AMPI virtualization module for LUNs 12/30/2019
C                                                                      *
C***********************************************************************

      module CkLUNMod
      implicit none
      save
C     obtain AMPI rank
C     multiply by LUN factor, default 100000
      integer CK_LUN
!$omp threadprivate(CK_LUN)      
      integer :: CkLUNFactor=100000

      contains 

C-----------------------------------------------------------------------
      integer function setAMPIzedLUN(thisRank)
      implicit none
      integer, intent(in) :: thisRank

      CK_LUN = thisRank * CkLUNFactor
      setAMPIzedLUN = CK_LUN
      RETURN
      
      end function setAMPIzedLUN
C-----------------------------------------------------------------------

      end module CkLUNMod
