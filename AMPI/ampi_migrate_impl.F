      module AMPI_MIGRATE_IMPL
         use AMPI_LUN_Migratable

         implicit none
#ifdef AMPI_COMP         
         include 'mpif.h'
#endif
      
       contains
       
         subroutine about_to_migrate
           implicit none

           integer :: luncount, ierr
           luncount = AMPI_LUN_close_registered();
         end subroutine about_to_migrate

         subroutine just_migrated
           implicit none
                    
           integer :: ierr
           ierr= AMPI_LUN_reopen_registered();
         end subroutine just_migrated

         subroutine migratablelun_setup()
           implicit none

           integer :: ierr, rank, numranks
           integer(8) :: computed
#ifdef AMPI_COMP           
!       register callbacks
           call ampi_register_about_to_migrate(about_to_migrate, ierr);
           call ampi_register_just_migrated(just_migrated, ierr);
!       create the lun registry
           call AMPI_LUN_create_registry(20000);
#endif
           ! open some files and add them to the registry
         end subroutine migratablelun_setup
       end module AMPI_MIGRATE_IMPL

