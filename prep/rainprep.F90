!-----------------------------------------------------------------------
!  MODULE RAINPREP
!-----------------------------------------------------------------------
!> @author Chris Szpilka, University of Oklahoma, cmszpilka@ou.edu
!> 
!> @copyright Dr. R.A. Luettich and Dr. J.J. Westerink
!>
!> @brief This module handles the I/O for adcprep related to adding
!> precipitation source terms to ADCIRC. See src/rain.F90 for
!> more details about these source terms in general, as well as
!> input requirements and expected files.
!>
!> rainType = 1  requires fort.270
!> rainType = 12 requires fort.27
!>                        fort.271
!>                        fort.272 (OPT regional)
!>-----------------------------------------------------------------------
!>  USAGE:
!>
!>----------------------------------------------------------------------- 
      module rainPrep

      implicit none

      ! Namelist default values
      logical :: found_rainControl_nml = .false.
      logical :: useRain = .false. ! default none
      integer :: rainType = 0 ! default none
      integer :: activeRainType = 1 ! default currently wet
      real(8) :: rainTimeInc = 3600 ! input time (sec)
      namelist /rainfallControl/ useRain, rainType, rainTimeInc,       &
     &          activeRainType


      ! --              Public variables
      public :: useRain, rainType, rainTimeInc, activeRainType
      public :: found_rainControl_nml, rainfallControl

      ! --              Public subroutines
      public :: checkRainUsage
      ! --              Private variables
      ! --              Private subroutines
      private :: checkRainFile

!>---------------------end of data declarations--------------------------


      CONTAINS

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                             P U B L I C  
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      subroutine checkRainUsage()

      implicit none
  
      !-------------------------------------------------
      ! Verify that requested rainType is valid and that
      !     files are in the main directory
      !-------------------------------------------------
      select case (abs(rainType))
      case(1)
         write(*,3260) rainType
         write(16,3260) rainType
         call checkRainFile(270)
      case(12)
         write(*,3260) rainType
         write(16,3260) rainType
         call checkRainFile(27)
         call checkRainFile(271)
         call checkRainFile(272)
      case default
         write(*,3261) rainType
         write(16,3261) rainType
         stop   !terminate 
      end select

 3260   format(/,2X, 'rainType = ',I3,  &
     &    /,2X,'Checking for global rain input:',/) 

 3261   format(/,2X,'useRain = T',  & 
     &    /,2X, 'rainType = ',I3,   &
     &    /,2X,'Your selection (a unit 15 parameter) is not valid.', &
     &    /,2X,'ERROR: adcprep terminated.')

      return

!-----------------------------------------------------------------------
      end subroutine checkRainUsage
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                             P R I V A T E  
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!> Called from within checkRainUsage() to verify all input is ready
!> for the run and provide warnings if files are missing.
!-----------------------------------------------------------------------
      subroutine checkRainFile(UnitNumber)
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in) :: UnitNumber     ! i/o unit number to open

      logical found               !.true. if the full domain file exists
      character(len=80) fileName   ! name of full domain file
      character(len=80) defaultName! default name of full domain file
      integer :: dnlen  !length of defaultname (7 or 8)
      integer :: unitwidth !length of unit number
      character (len=20):: unitnumberstr

      found = .false.
      defaultName(:) = ' '  !initialize to all blanks
      defaultName(1:5) = 'fort.'

      unitwidth = ceiling(log(dble(UnitNumber))/log(10.0)) ;        
      dnlen = 5 + unitwidth ;   
      write(unitnumberstr,'(I0)') UnitNumber ;
      defaultName = trim(adjustl(defaultName))// &
     &              trim(adjustl(unitNumberStr))

      ! Search for global file and give error message if not found
      fileName = trim(defaultName(1:dnlen))

      ! Determine if full domain file exists
      inquire(file=fileName,exist=found)

      ! Log messages to user based on file existence
      if ( found ) then
         write(*,1010) FileName
         write(16,1010) FileName
      else
         if (UnitNumber .eq. 272) then
            write(*,1011) FileName
            write(16,1011) FileName
         else
            write(*,1012) FileName
            write(16,1012) FileName
         endif
      end if


 1010 format(2X,'Global file ',A8,' was found!',/)
 1011 format(2X,'WARNING: Regional basin file ',A8,' was NOT found!',   &
     &     /,2X,'adcprep will complete but padcirc may fail at',/,      &
     &       2X,'runtime. Check fort.27 file.')
 1012 format(2X,'WARNING: Global file ',A8,' was NOT found!',           &
     &     /,2X,'adcprep will complete but padcirc will fail at',/,     &
     &       2X,'runtime: add file before running.')

      return
!------------------------------------------------------------------------
      end subroutine checkRainFile
!------------------------------------------------------------------------


!-----------------------------------------------------------------------
      end module rainPrep
