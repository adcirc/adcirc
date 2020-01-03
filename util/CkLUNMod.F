C***********************************************************************
C AMPI virtualization module for LUNs 12/30/2019
C                                                                      *
C***********************************************************************

      module CkLUNMod
      implicit none
      save
C     obtain AMPI rank
C     multiply by LUN factor, default 100000
      integer AMPIzedLUN
!$omp threadprivate(AMPIzedLUN)      
      integer :: CkLUNFactor=100000

      contains 

C-----------------------------------------------------------------------
      integer function setAMPIzedLUN(thisRank)
      implicit none
      integer, intent(in) :: thisRank

      AMPIzedLUN = thisRank * CkLUNFactor
      setAMPIzedLUN = AMPIzedLUN
      RETURN
      
      end function setAMPIzedLUN
C-----------------------------------------------------------------------

      end module CkLUNMod
