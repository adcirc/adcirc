!**********************************************************************
! MODULE: NODAL_ATTRIBUTE_MODULE
! DESCRIPTION: Module containing types and procedures related to nodal
!              attributes in the ADCIRC model.
! WRITTEN BY:  Zachary Cobell
! CONTACT INFO: zcobell@thewaterinstitute.org
!
! This module is use to define the base types for nodal attributes within
! the ADCIRC model. These types can be extended to create specific nodal
! attributes with their own procedures and functions for use within the
! model.
!
!**********************************************************************

MODULE NODAL_ATTRIBUTE_BASE_TYPE_MODULE

  ! Maximum length for character attributes
  integer, parameter :: MAX_NODAL_ATTRIBUTE_CHARACTER_LENGTH = 128

  !**********************************************************************
  ! TYPE: NODAL_ATTRIBUTE_BASE_T
  ! DESCRIPTION: Base type for nodal attributes
  !**********************************************************************
  type, abstract :: NODAL_ATTRIBUTE_BASE_T
    character(len=MAX_NODAL_ATTRIBUTE_CHARACTER_LENGTH) :: attribute_name  !< Name of the nodal attribute
    character(len=MAX_NODAL_ATTRIBUTE_CHARACTER_LENGTH) :: units           !< Units of the nodal attribute
    integer :: n_values                                                    !< Number of values for the nodal attribute
    logical :: loaded = .false.                                            !< Flag indicating whether the nodal attribute has been loaded from file
    logical :: active = .false.                                            !< Flag indicating whether the nodal attribute is active

  contains

    !**********************************************************************
    ! PROCEDURE: init
    ! DESCRIPTION: Initialize the nodal attribute
    !**********************************************************************
    ! procedure(nodal_attribute_t_constructor_interface), pass(self), deferred :: init
    !**********************************************************************

    !**********************************************************************
    ! PROCEDURE: unload
    ! DESCRIPTION: Unload the nodal attribute
    !**********************************************************************
    procedure(nodal_attribute_t_unload_interface), pass(self), deferred :: unload
    !**********************************************************************

    !**********************************************************************
    ! PROCEDURE: read_header_from_file
    ! DESCRIPTION: Read header information from file for the nodal attribute
    !**********************************************************************
    procedure(nodal_attribute_t_read_header_from_file_interface), pass(self), deferred :: read_header_from_file
    !**********************************************************************

    !**********************************************************************
    ! PROCEDURE: read_body_from_file
    ! DESCRIPTION: Read body information from file for the nodal attribute
    !**********************************************************************
    procedure(nodal_attribute_t_read_body_from_file_interface), pass(self), deferred :: read_body_from_file
    !**********************************************************************

  end type NODAL_ATTRIBUTE_BASE_T

  !**********************************************************************
  ! INTERFACE
  ! DESCRIPTION: Interface for nodal attribute procedures
  !**********************************************************************
  interface

    !**********************************************************************
    ! PROCEDURE: nodal_attribute_t_unload_interface
    ! DESCRIPTION: Interface for unloading nodal attributes
    !**********************************************************************
    subroutine nodal_attribute_t_unload_interface(self)
      import NODAL_ATTRIBUTE_BASE_T
      class(NODAL_ATTRIBUTE_BASE_T), intent(inout) :: self          !< Nodal attribute instance
    end subroutine nodal_attribute_t_unload_interface

    !**********************************************************************
    ! PROCEDURE: nodal_attribute_t_read_header_from_file_interface
    ! DESCRIPTION: Interface for reading header information from file
    !**********************************************************************
    subroutine nodal_attribute_t_read_header_from_file_interface(self, io_unit, n_nodes)
      import NODAL_ATTRIBUTE_BASE_T
      class(NODAL_ATTRIBUTE_BASE_T), intent(inout) :: self          !< Nodal attribute instance
      integer, intent(in) :: io_unit                                !< File unit number
      integer, intent(in) :: n_nodes                                !< Number of nodes
    end subroutine nodal_attribute_t_read_header_from_file_interface

    !**********************************************************************
    ! PROCEDURE: nodal_attribute_t_read_body_from_file_interface
    ! DESCRIPTION: Interface for reading body information from file
    !**********************************************************************
    subroutine nodal_attribute_t_read_body_from_file_interface(self, io_unit)
      import NODAL_ATTRIBUTE_BASE_T
      class(NODAL_ATTRIBUTE_BASE_T), intent(inout) :: self          !< Nodal attribute instance
      integer, intent(in) :: io_unit                                !< File unit number
    end subroutine nodal_attribute_t_read_body_from_file_interface
  end interface

  !**********************************************************************
  ! TYPE: NODAL_ATTRIBUTE_VECTOR_REAL_T
  ! DESCRIPTION: Type for nodal attributes with real values
  !**********************************************************************
  type, extends(NODAL_ATTRIBUTE_BASE_T) :: NODAL_ATTRIBUTE_VECTOR_REAL_T
    real(8)              :: default_value                           !< Default value for the nodal attribute
    real(8), allocatable :: values(:)                               !< Values for the nodal attribute

  contains

    !**********************************************************************
    ! PROCEDURE: init
    ! DESCRIPTION: Initialize the nodal attribute
    !**********************************************************************
    procedure :: init => nodal_attribute_vector_real_t_constructor
    !**********************************************************************

    !**********************************************************************
    ! PROCEDURE: unload
    ! DESCRIPTION: Unload the nodal attribute
    !**********************************************************************
    procedure :: unload => nodal_attribute_vector_real_t_unload
    !**********************************************************************

    !**********************************************************************
    ! PROCEDURE: read_header_from_file
    ! DESCRIPTION: Read header information from file for the nodal attribute
    !**********************************************************************
    procedure :: read_header_from_file => nodal_attribute_vector_real_t_read_header_from_file
    !**********************************************************************

    !**********************************************************************
    ! PROCEDURE: read_body_from_file
    ! DESCRIPTION: Read body information from file for the nodal attribute
    !**********************************************************************
    procedure :: read_body_from_file => nodal_attribute_vector_real_t_read_body_from_file
    !**********************************************************************

    !**********************************************************************
    ! PROCEDURE: final
    ! DESCRIPTION: Destructor for the nodal attribute
    !**********************************************************************
    final :: nodal_attribute_vector_real_t_destructor

  end type NODAL_ATTRIBUTE_VECTOR_REAL_T

  !**********************************************************************
  ! TYPE: NODAL_ATTRIBUTE_VECTOR_INT_T
  ! DESCRIPTION: Type for nodal attributes with integer values
  !**********************************************************************
  type, extends(NODAL_ATTRIBUTE_BASE_T) :: NODAL_ATTRIBUTE_VECTOR_INT_T
    integer :: default_value
    integer, allocatable :: values(:)

  contains

    !**********************************************************************
    ! PROCEDURE: init
    ! DESCRIPTION: Initialize the nodal attribute
    !**********************************************************************
    procedure :: init => nodal_attribute_vector_int_t_constructor
    !**********************************************************************

    !**********************************************************************
    ! PROCEDURE: unload
    ! DESCRIPTION: Unload the nodal attribute
    !**********************************************************************
    procedure :: unload => nodal_attribute_vector_int_t_unload
    !**********************************************************************

    !**********************************************************************
    ! PROCEDURE: read_header_from_file
    ! DESCRIPTION: Read header information from file for the nodal attribute
    !**********************************************************************
    procedure :: read_header_from_file => nodal_attribute_vector_int_t_read_header_from_file
    !**********************************************************************

    !**********************************************************************
    ! PROCEDURE: read_body_from_file
    ! DESCRIPTION: Read body information from file for the nodal attribute
    !**********************************************************************
    procedure :: read_body_from_file => nodal_attribute_vector_int_t_read_body_from_file
    !**********************************************************************

    !**********************************************************************
    ! PROCEDURE: final
    ! DESCRIPTION: Destructor for the nodal attribute
    !**********************************************************************
    final :: nodal_attribute_vector_int_t_final

  end type NODAL_ATTRIBUTE_VECTOR_INT_T

  !**********************************************************************
  ! TYPE: NODAL_ATTRIBUTE_VECTOR_LOGICAL_T
  ! DESCRIPTION: Type for nodal attributes with logical values
  !**********************************************************************
  type, extends(NODAL_ATTRIBUTE_BASE_T) :: NODAL_ATTRIBUTE_VECTOR_LOGICAL_T
    logical :: default_value
    logical, allocatable :: values(:)

  contains

    !**********************************************************************
    ! PROCEDURE: init
    ! DESCRIPTION: Initialize the nodal attribute
    !**********************************************************************
    procedure :: init => nodal_attribute_vector_logical_t_constructor
    !**********************************************************************

    !**********************************************************************
    ! PROCEDURE: unload
    ! DESCRIPTION: Unload the nodal attribute
    !**********************************************************************
    procedure :: unload => nodal_attribute_vector_logical_t_unload
    !**********************************************************************

    !**********************************************************************
    ! PROCEDURE: read_header_from_file
    ! DESCRIPTION: Read header information from file for the nodal attribute
    !**********************************************************************
    procedure :: read_header_from_file => nodal_attribute_vector_logical_t_read_header_from_file
    !**********************************************************************

    !**********************************************************************
    ! PROCEDURE: read_body_from_file
    ! DESCRIPTION: Read body information from file for the nodal attribute
    !**********************************************************************
    procedure :: read_body_from_file => nodal_attribute_vector_logical_t_read_body_from_file
    !**********************************************************************

    !**********************************************************************
    ! PROCEDURE: final
    ! DESCRIPTION: Destructor for the nodal attribute
    !**********************************************************************
    final :: nodal_attribute_vector_logical_t_final

  end type NODAL_ATTRIBUTE_VECTOR_LOGICAL_T

  !**********************************************************************
  ! TYPE: NODAL_ATTRIBUTE_MATRIX_REAL_T
  ! DESCRIPTION: Type for nodal attributes with real matrix values
  !**********************************************************************
  type, extends(NODAL_ATTRIBUTE_BASE_T) :: NODAL_ATTRIBUTE_MATRIX_REAL_T
    real(8), allocatable :: default_value(:)                       !< Default value for the nodal attribute
    real(8), allocatable :: values(:,:)                            !< Values for the nodal attribute

  contains

    !**********************************************************************
    ! PROCEDURE: init
    ! DESCRIPTION: Initialize the nodal attribute
    !**********************************************************************
    procedure :: init => nodal_attribute_matrix_real_t_constructor
    !**********************************************************************

    !**********************************************************************
    ! PROCEDURE: unload
    ! DESCRIPTION: Unload the nodal attribute
    !**********************************************************************
    procedure :: unload => nodal_attribute_matrix_real_t_unload
    !**********************************************************************

    !**********************************************************************
    ! PROCEDURE: read_header_from_file
    ! DESCRIPTION: Read header information from file for the nodal attribute
    !**********************************************************************
    procedure :: read_header_from_file => nodal_attribute_matrix_real_t_read_header_from_file
    !**********************************************************************

    !**********************************************************************
    ! PROCEDURE: read_body_from_file
    ! DESCRIPTION: Read body information from file for the nodal attribute
    !**********************************************************************
    procedure :: read_body_from_file => nodal_attribute_matrix_real_t_read_body_from_file
    !**********************************************************************

    !**********************************************************************
    ! PROCEDURE: final
    ! DESCRIPTION: Destructor for the nodal attribute
    !**********************************************************************
    final :: nodal_attribute_matrix_real_t_final

  end type NODAL_ATTRIBUTE_MATRIX_REAL_T

  !**********************************************************************
  ! TYPE: NODAL_ATTRIBUTE_MATRIX_INT_T
  ! DESCRIPTION: Type for nodal attributes with integer matrix values
  !**********************************************************************
  type, extends(NODAL_ATTRIBUTE_BASE_T) :: NODAL_ATTRIBUTE_MATRIX_INT_T
    integer, allocatable :: default_value(:)
    integer, allocatable :: values(:,:)

  contains

    !**********************************************************************
    ! PROCEDURE: init
    ! DESCRIPTION: Initialize the nodal attribute
    !**********************************************************************
    procedure :: init => nodal_attribute_matrix_int_t_constructor
    !**********************************************************************

    !**********************************************************************
    ! PROCEDURE: unload
    ! DESCRIPTION: Unload the nodal attribute
    !**********************************************************************
    procedure :: unload => nodal_attribute_matrix_int_t_unload
    !**********************************************************************

    !**********************************************************************
    ! PROCEDURE: read_header_from_file
    ! DESCRIPTION: Read header information from file for the nodal attribute
    !**********************************************************************
    procedure :: read_header_from_file => nodal_attribute_matrix_int_t_read_header_from_file
    !**********************************************************************

    !**********************************************************************
    ! PROCEDURE: read_body_from_file
    ! DESCRIPTION: Read body information from file for the nodal attribute
    !**********************************************************************
    procedure :: read_body_from_file => nodal_attribute_matrix_int_t_read_body_from_file
    !**********************************************************************

    !**********************************************************************
    ! PROCEDURE: final
    ! DESCRIPTION: Destructor for the nodal attribute
    !**********************************************************************
    final :: nodal_attribute_matrix_int_t_final

  end type NODAL_ATTRIBUTE_MATRIX_INT_T

  !**********************************************************************
  ! TYPE: NODAL_ATTRIBUTE_MATRIX_LOGICAL_T
  ! DESCRIPTION: Type for nodal attributes with logical matrix values
  !**********************************************************************
  type, extends(NODAL_ATTRIBUTE_BASE_T) :: NODAL_ATTRIBUTE_MATRIX_LOGICAL_T
    logical, allocatable :: default_value(:)
    logical, allocatable :: values(:,:)

  contains

    !**********************************************************************
    ! PROCEDURE: init
    ! DESCRIPTION: Initialize the nodal attribute
    !**********************************************************************
    procedure :: init => nodal_attribute_matrix_logical_t_constructor
    !**********************************************************************

    !**********************************************************************
    ! PROCEDURE: unload
    ! DESCRIPTION: Unload the nodal attribute
    !**********************************************************************
    procedure :: unload => nodal_attribute_matrix_logical_t_unload
    !**********************************************************************

    !**********************************************************************
    ! PROCEDURE: read_header_from_file
    ! DESCRIPTION: Read header information from file for the nodal attribute
    !**********************************************************************
    procedure :: read_header_from_file => nodal_attribute_matrix_logical_t_read_header_from_file
    !**********************************************************************

    !**********************************************************************
    ! PROCEDURE: read_body_from_file
    ! DESCRIPTION: Read body information from file for the nodal attribute
    !**********************************************************************
    procedure :: read_body_from_file => nodal_attribute_matrix_logical_t_read_body_from_file
    !**********************************************************************

    !**********************************************************************
    ! PROCEDURE: final
    ! DESCRIPTION: Destructor for the nodal attribute
    !**********************************************************************
    final :: nodal_attribute_matrix_logical_t_final

  end type NODAL_ATTRIBUTE_MATRIX_LOGICAL_T

  !**********************************************************************
  ! IMPLEMENTATION
  ! DESCRIPTION: Implementation of nodal attribute procedures
  !**********************************************************************
contains

  !**********************************************************************
  ! PROCEDURE: nodal_attribute_vector_real_t_constructor
  ! DESCRIPTION: Constructor for the nodal attribute with real vector values
  !**********************************************************************
  subroutine nodal_attribute_vector_real_t_constructor(self, n_nodes, attribute_name, units, n_values, default_value)
    implicit none

    class(NODAL_ATTRIBUTE_VECTOR_REAL_T), intent(inout) :: self   !< Nodal attribute instance
    integer, intent(in) :: n_nodes                                !< Number of nodes
    character(len=*), intent(in) :: attribute_name                !< Name of the nodal attribute
    character(len=*), intent(in) :: units                         !< Units of the nodal attribute
    integer, intent(in) :: n_values                               !< Number of values for the nodal attribute
    real(8), intent(in) :: default_value                          !< Default value for the nodal attribute

    self%attribute_name = trim(adjustl(attribute_name))
    self%units = trim(adjustl(units))
    self%n_values = n_values
    self%active = .true.
    self%default_value = default_value
    allocate(self%values(n_nodes))
    self%values(1:n_nodes) = default_value

  end subroutine nodal_attribute_vector_real_t_constructor

  !**********************************************************************
  ! PROCEDURE: nodal_attribute_vector_int_t_constructor
  ! DESCRIPTION: Constructor for the nodal attribute with integer vector values
  !**********************************************************************
  subroutine nodal_attribute_vector_int_t_constructor(self, n_nodes, attribute_name, units, n_values, default_value)
    implicit none

    class(NODAL_ATTRIBUTE_VECTOR_INT_T), intent(inout) :: self   !< Nodal attribute instance
    integer, intent(in) :: n_nodes                                !< Number of nodes
    character(len=*), intent(in) :: attribute_name                !< Name of the nodal attribute
    character(len=*), intent(in) :: units                         !< Units of the nodal attribute
    integer, intent(in) :: n_values                               !< Number of values for the nodal attribute
    integer, intent(in) :: default_value                          !< Default value for the nodal attribute

    self%attribute_name = trim(adjustl(attribute_name))
    self%units = trim(adjustl(units))
    self%n_values = n_values
    self%active = .true.
    self%default_value = default_value
    allocate(self%values(n_nodes))
    self%values(1:n_nodes) = default_value

  end subroutine nodal_attribute_vector_int_t_constructor

  !**********************************************************************
  ! PROCEDURE: nodal_attribute_vector_logical_t_constructor
  ! DESCRIPTION: Constructor for the nodal attribute with logical vector values
  !**********************************************************************
  subroutine nodal_attribute_vector_logical_t_constructor(self, n_nodes, attribute_name, units, n_values, default_value)
    implicit none

    class(NODAL_ATTRIBUTE_VECTOR_LOGICAL_T), intent(inout) :: self   !< Nodal attribute instance
    integer, intent(in) :: n_nodes                                !< Number of nodes
    character(len=*), intent(in) :: attribute_name                !< Name of the nodal attribute
    character(len=*), intent(in) :: units                         !< Units of the nodal attribute
    integer, intent(in) :: n_values                               !< Number of values for the nodal attribute
    logical, intent(in) :: default_value                          !< Default value for the nodal attribute

    self%attribute_name = trim(adjustl(attribute_name))
    self%units = trim(adjustl(units))
    self%n_values = n_values
    self%active = .true.
    self%default_value = default_value
    allocate(self%values(n_nodes))
    self%values(1:n_nodes) = default_value

  end subroutine nodal_attribute_vector_logical_t_constructor

  !**********************************************************************
  ! PROCEDURE: nodal_attribute_matrix_t_constructor
  ! DESCRIPTION: Constructor for the nodal attribute with real matrix values
  !**********************************************************************
  subroutine nodal_attribute_matrix_real_t_constructor(self, n_nodes, attribute_name, units, n_values, default_value)
    implicit none

    class(NODAL_ATTRIBUTE_MATRIX_REAL_T), intent(inout) :: self   !< Nodal attribute instance
    integer, intent(in) :: n_nodes                                !< Number of nodes
    character(len=*), intent(in) :: attribute_name                !< Name of the nodal attribute
    character(len=*), intent(in) :: units                         !< Units of the nodal attribute
    integer, intent(in) :: n_values                               !< Number of values for the nodal attribute
    real(8), intent(in) :: default_value(n_values)                !< Default value for the nodal attribute
    integer :: i

    self%attribute_name = trim(adjustl(attribute_name))
    self%units = trim(adjustl(units))
    self%n_values = n_values
    self%active = .true.
    allocate(self%default_value(n_values))
    self%default_value(1:n_values) = default_value(1:n_values)
    allocate(self%values(n_nodes, n_values))

    do i = 1, n_nodes
      self%values(i, 1:n_values) = default_value(1:n_values)
    end do

  end subroutine nodal_attribute_matrix_real_t_constructor

  !**********************************************************************
  ! PROCEDURE: nodal_attribute_matrix_int_t_constructor
  ! DESCRIPTION: Constructor for the nodal attribute with integer matrix values
  !**********************************************************************
  subroutine nodal_attribute_matrix_int_t_constructor(self, n_nodes, attribute_name, units, n_values, default_value)
    implicit none

    class(NODAL_ATTRIBUTE_MATRIX_INT_T), intent(inout) :: self    !< Nodal attribute instance
    integer, intent(in) :: n_nodes                                !< Number of nodes
    character(len=*), intent(in) :: attribute_name                !< Name of the nodal attribute
    character(len=*), intent(in) :: units                         !< Units of the nodal attribute
    integer, intent(in) :: n_values                               !< Number of values for the nodal attribute
    integer, intent(in) :: default_value(n_values)                !< Default value for the nodal attribute
    integer :: i

    self%attribute_name = trim(adjustl(attribute_name))
    self%units = trim(adjustl(units))
    self%n_values = n_values
    self%active = .true.
    allocate(self%default_value(n_values))
    self%default_value(1:n_values) = default_value(1:n_values)
    allocate(self%values(n_nodes, n_values))

    do i = 1, n_nodes
      self%values(i, 1:n_values) = default_value(1:n_values)
    end do

  end subroutine nodal_attribute_matrix_int_t_constructor

  !**********************************************************************
  ! PROCEDURE: nodal_attribute_matrix_logical_t_constructor
  ! DESCRIPTION: Constructor for the nodal attribute with logical matrix values
  !**********************************************************************
  subroutine nodal_attribute_matrix_logical_t_constructor(self, n_nodes, attribute_name, units, n_values, default_value)
    implicit none

    class(NODAL_ATTRIBUTE_MATRIX_LOGICAL_T), intent(inout) :: self   !< Nodal attribute instance
    integer, intent(in) :: n_nodes                                   !< Number of nodes
    character(len=*), intent(in) :: attribute_name                   !< Name of the nodal attribute
    character(len=*), intent(in) :: units                            !< Units of the nodal attribute
    integer, intent(in) :: n_values                                  !< Number of values for the nodal attribute
    logical, intent(in) :: default_value(n_values)                   !< Default value for the nodal attribute
    integer :: i

    self%attribute_name = trim(adjustl(attribute_name))
    self%units = trim(adjustl(units))
    self%n_values = n_values
    self%active = .true.
    allocate(self%default_value(n_values))
    self%default_value(1:n_values) = default_value(1:n_values)
    allocate(self%values(n_nodes, n_values))

    do i = 1, n_nodes
      self%values(i, 1:n_values) = default_value(1:n_values)
    end do

  end subroutine nodal_attribute_matrix_logical_t_constructor

  !**********************************************************************
  ! PROCEDURE: nodal_attribute_vector_real_t_destructor
  ! DESCRIPTION: Destructor for the nodal attribute with real vector values
  !**********************************************************************
  subroutine nodal_attribute_vector_real_t_destructor(self)
    implicit none
    type(NODAL_ATTRIBUTE_VECTOR_REAL_T), intent(inout) :: self        !< Nodal attribute instance

    if (allocated(self%values)) deallocate(self%values)

  end subroutine nodal_attribute_vector_real_t_destructor

  !**********************************************************************
  ! PROCEDURE: nodal_attribute_vector_int_t_final
  ! DESCRIPTION: Destructor for the nodal attribute with integer vector values
  !**********************************************************************
  subroutine nodal_attribute_vector_int_t_final(self)
    implicit none
    type(NODAL_ATTRIBUTE_VECTOR_INT_T), intent(inout) :: self        !< Nodal attribute instance

    if (allocated(self%values)) deallocate(self%values)

  end subroutine nodal_attribute_vector_int_t_final

  !**********************************************************************
  ! PROCEDURE: nodal_attribute_vector_logical_t_final
  ! DESCRIPTION: Destructor for the nodal attribute with logical vector values
  !**********************************************************************
  subroutine nodal_attribute_vector_logical_t_final(self)
    implicit none
    type(NODAL_ATTRIBUTE_VECTOR_LOGICAL_T), intent(inout) :: self        !< Nodal attribute instance

    if (allocated(self%values)) deallocate(self%values)

  end subroutine nodal_attribute_vector_logical_t_final

  !**********************************************************************
  ! PROCEDURE: nodal_attribute_matrix_real_t_final
  ! DESCRIPTION: Destructor for the nodal attribute with real matrix values
  !**********************************************************************
  subroutine nodal_attribute_matrix_real_t_final(self)
    implicit none
    type(NODAL_ATTRIBUTE_MATRIX_REAL_T), intent(inout) :: self        !< Nodal attribute instance

    if (allocated(self%default_value)) deallocate(self%default_value)
    if (allocated(self%values)) deallocate(self%values)

  end subroutine nodal_attribute_matrix_real_t_final

  !**********************************************************************
  ! PROCEDURE: nodal_attribute_matrix_int_t_final
  ! DESCRIPTION: Destructor for the nodal attribute with integer matrix values
  !**********************************************************************
  subroutine nodal_attribute_matrix_int_t_final(self)
    implicit none
    type(NODAL_ATTRIBUTE_MATRIX_INT_T), intent(inout) :: self        !< Nodal attribute instance

    if (allocated(self%default_value)) deallocate(self%default_value)
    if (allocated(self%values)) deallocate(self%values)

  end subroutine nodal_attribute_matrix_int_t_final

  !**********************************************************************
  ! PROCEDURE: nodal_attribute_matrix_logical_t_final
  ! DESCRIPTION: Destructor for the nodal attribute with logical matrix values
  !**********************************************************************
  subroutine nodal_attribute_matrix_logical_t_final(self)
    implicit none
    type(NODAL_ATTRIBUTE_MATRIX_LOGICAL_T), intent(inout) :: self        !< Nodal attribute instance

    if (allocated(self%default_value)) deallocate(self%default_value)
    if (allocated(self%values)) deallocate(self%values)

  end subroutine nodal_attribute_matrix_logical_t_final

  !**********************************************************************
  ! PROCEDURE: nodal_attribute_vector_real_t_unload
  ! DESCRIPTION: Unload the nodal attribute with real vector values
  !**********************************************************************
  subroutine nodal_attribute_vector_real_t_unload(self)
    implicit none
    class(NODAL_ATTRIBUTE_VECTOR_REAL_T), intent(inout) :: self        !< Nodal attribute instance

    if (allocated(self%values)) deallocate(self%values)
    self%loaded = .false.

  end subroutine nodal_attribute_vector_real_t_unload

  !**********************************************************************
  ! PROCEDURE: nodal_attribute_vector_int_t_unload
  ! DESCRIPTION: Unload the nodal attribute with integer vector values
  !**********************************************************************
  subroutine nodal_attribute_vector_int_t_unload(self)
    implicit none
    class(NODAL_ATTRIBUTE_VECTOR_INT_T), intent(inout) :: self        !< Nodal attribute instance

    if (allocated(self%values)) deallocate(self%values)
    self%loaded = .false.

  end subroutine nodal_attribute_vector_int_t_unload

  !**********************************************************************
  ! PROCEDURE: nodal_attribute_vector_logical_t_unload
  ! DESCRIPTION: Unload the nodal attribute with logical vector values
  !**********************************************************************
  subroutine nodal_attribute_vector_logical_t_unload(self)
    implicit none
    class(NODAL_ATTRIBUTE_VECTOR_LOGICAL_T), intent(inout) :: self        !< Nodal attribute instance

    if (allocated(self%values)) deallocate(self%values)
    self%loaded = .false.

  end subroutine nodal_attribute_vector_logical_t_unload

  !**********************************************************************
  ! PROCEDURE: nodal_attribute_vector_real_t_read_header_from_file
  ! DESCRIPTION: Read header information from file for the nodal attribute with real vector values
  !**********************************************************************
    subroutine nodal_attribute_vector_real_t_read_header_from_file(self, io_unit, n_nodes)
      implicit none
      class(NODAL_ATTRIBUTE_VECTOR_REAL_T), intent(inout) :: self        !< Nodal attribute instance
      integer, intent(in) :: io_unit                                !< File unit number
      integer, intent(in) :: n_nodes                                !< Number of nodes
      character(len=MAX_NODAL_ATTRIBUTE_CHARACTER_LENGTH) :: attribute_name !< Name of the nodal attribute

      read(io_unit, *) attribute_name
      self%attribute_name = trim(adjustl(attribute_name))

      read(io_unit, *) self%units
      read(io_unit, *) self%n_values
      read(io_unit, *) self%default_value

      if (self%n_values /= 1) then
        write(*,'(A)') '[ERROR]: NODAL_ATTRIBUTE_VECTOR_T only supports a single value'
        stop
      end if

      allocate(self%values(n_nodes))
      self%values(1:n_nodes) = self%default_value

    end subroutine nodal_attribute_vector_real_t_read_header_from_file

    !**********************************************************************
    ! PROCEDURE: nodal_attribute_vector_int_t_read_header_from_file
    ! DESCRIPTION: Read header information from file for the nodal attribute with integer vector values
    !**********************************************************************
    subroutine nodal_attribute_vector_int_t_read_header_from_file(self, io_unit, n_nodes)
      implicit none
      class(NODAL_ATTRIBUTE_VECTOR_INT_T), intent(inout) :: self    !< Nodal attribute instance
      integer, intent(in) :: io_unit                                !< File unit number
      integer, intent(in) :: n_nodes                                !< Number of nodes
      character(len=MAX_NODAL_ATTRIBUTE_CHARACTER_LENGTH) :: attribute_name !< Name of the nodal attribute

      read(io_unit, *) attribute_name
      self%attribute_name = trim(adjustl(attribute_name))

      read(io_unit, *) self%units
      read(io_unit, *) self%n_values
      read(io_unit, *) self%default_value

      if (self%n_values /= 1) then
        write(*,'(A)') '[ERROR]: NODAL_ATTRIBUTE_VECTOR_T only supports a single value'
        stop
      end if

      allocate(self%values(n_nodes))
      self%values(1:n_nodes) = self%default_value

    end subroutine nodal_attribute_vector_int_t_read_header_from_file

    !**********************************************************************
    ! PROCEDURE: nodal_attribute_vector_logical_t_read_header_from_file
    ! DESCRIPTION: Read header information from file for the nodal attribute with logical vector values
    !**********************************************************************
    subroutine nodal_attribute_vector_logical_t_read_header_from_file(self, io_unit, n_nodes)
      implicit none
      class(NODAL_ATTRIBUTE_VECTOR_LOGICAL_T), intent(inout) :: self        !< Nodal attribute instance
      integer, intent(in) :: io_unit                                !< File unit number
      integer, intent(in) :: n_nodes                                !< Number of nodes
      character(len=MAX_NODAL_ATTRIBUTE_CHARACTER_LENGTH) :: attribute_name !< Name of the nodal attribute

      read(io_unit, *) attribute_name
      self%attribute_name = trim(adjustl(attribute_name))

      read(io_unit, *) self%units
      read(io_unit, *) self%n_values
      read(io_unit, *) self%default_value

      if (self%n_values /= 1) then
        write(*,'(A)') '[ERROR]: NODAL_ATTRIBUTE_VECTOR_T only supports a single value'
        stop
      end if

      allocate(self%values(n_nodes))
      self%values(1:n_nodes) = self%default_value

    end subroutine nodal_attribute_vector_logical_t_read_header_from_file

    !**********************************************************************
    ! PROCEDURE: nodal_attribute_vector_real_t_read_body_from_file
    ! DESCRIPTION: Read body information from file for the nodal attribute with real vector values
    !**********************************************************************
    subroutine nodal_attribute_vector_real_t_read_body_from_file(self, io_unit)
      implicit none
      class(NODAL_ATTRIBUTE_VECTOR_REAL_T), intent(inout) :: self     !< Nodal attribute instance
      integer, intent(in) :: io_unit                             !< File unit number
      integer             :: n_non_default                       !< Number of non-default values
      integer             :: i                                   !< Loop index
      integer             :: node_number                         !< Node number
      character(len=MAX_NODAL_ATTRIBUTE_CHARACTER_LENGTH) :: attribute_name !< Name of the nodal attribute

      read(io_unit, *) attribute_name
      if(self%attribute_name /= trim(adjustl(attribute_name)))then
        write(*,'(A)') '[ERROR]: NODAL_ATTRIBUTE_VECTOR_T read from file with incorrect attribute name'
        stop
      end if

      read(io_unit, *) n_non_default
      do i = 1, n_non_default
        read(io_unit, *) node_number, self%values(node_number)
      end do
      self%loaded = .true.

    end subroutine nodal_attribute_vector_real_t_read_body_from_file

    !**********************************************************************
    ! PROCEDURE: nodal_attribute_vector_int_t_read_body_from_file
    ! DESCRIPTION: Read body information from file for the nodal attribute with integer vector values
    !**********************************************************************
    subroutine nodal_attribute_vector_int_t_read_body_from_file(self, io_unit)
      implicit none
      class(NODAL_ATTRIBUTE_VECTOR_INT_T), intent(inout) :: self     !< Nodal attribute instance
      integer, intent(in) :: io_unit                             !< File unit number
      integer             :: n_non_default                       !< Number of non-default values
      integer             :: i                                   !< Loop index
      integer             :: node_number                         !< Node number
      character(len=MAX_NODAL_ATTRIBUTE_CHARACTER_LENGTH) :: attribute_name !< Name of the nodal attribute

      read(io_unit, *) attribute_name
      if(self%attribute_name /= trim(adjustl(attribute_name)))then
        write(*,'(A)') '[ERROR]: NODAL_ATTRIBUTE_VECTOR_T read from file with incorrect attribute name'
        stop
      end if

      read(io_unit, *) n_non_default
      do i = 1, n_non_default
        read(io_unit, *) node_number, self%values(node_number)
      end do
      self%loaded = .true.

    end subroutine nodal_attribute_vector_int_t_read_body_from_file

    !**********************************************************************
    ! PROCEDURE: nodal_attribute_vector_logical_t_read_body_from_file
    ! DESCRIPTION: Read body information from file for the nodal attribute with logical vector values
    !**********************************************************************
    subroutine nodal_attribute_vector_logical_t_read_body_from_file(self, io_unit)
      implicit none
      class(NODAL_ATTRIBUTE_VECTOR_LOGICAL_T), intent(inout) :: self     !< Nodal attribute instance
      integer, intent(in) :: io_unit                             !< File unit number
      integer             :: n_non_default                       !< Number of non-default values
      integer             :: i                                   !< Loop index
      integer             :: node_number                         !< Node number
      character(len=MAX_NODAL_ATTRIBUTE_CHARACTER_LENGTH) :: attribute_name !< Name of the nodal attribute

      read(io_unit, *) attribute_name
      if(self%attribute_name /= trim(adjustl(attribute_name)))then
        write(*,'(A)') '[ERROR]: NODAL_ATTRIBUTE_VECTOR_T read from file with incorrect attribute name'
        stop
      end if

      read(io_unit, *) n_non_default
      do i = 1, n_non_default
        read(io_unit, *) node_number, self%values(node_number)
      end do
      self%loaded = .true.

    end subroutine nodal_attribute_vector_logical_t_read_body_from_file

    !**********************************************************************
    ! PROCEDURE: nodal_attribute_matrix_real_t_unload
    ! DESCRIPTION: Unload the nodal attribute with real matrix values
    !**********************************************************************
    subroutine nodal_attribute_matrix_real_t_unload(self)
      implicit none
      class(NODAL_ATTRIBUTE_MATRIX_REAL_T), intent(inout) :: self        !< Nodal attribute instance

      if(allocated(self%default_value))deallocate(self%default_value)
      if(allocated(self%values))deallocate(self%values)
      self%loaded = .false.

    end subroutine nodal_attribute_matrix_real_t_unload

    !**********************************************************************
    ! PROCEDURE: nodal_attribute_matrix_int_t_unload
    ! DESCRIPTION: Unload the nodal attribute with integer matrix values
    !**********************************************************************
    subroutine nodal_attribute_matrix_int_t_unload(self)
      implicit none
      class(NODAL_ATTRIBUTE_MATRIX_INT_T), intent(inout) :: self        !< Nodal attribute instance

      if(allocated(self%default_value))deallocate(self%default_value)
      if(allocated(self%values))deallocate(self%values)
      self%loaded = .false.

    end subroutine nodal_attribute_matrix_int_t_unload

    !**********************************************************************
    ! PROCEDURE: nodal_attribute_matrix_logical_t_unload
    ! DESCRIPTION: Unload the nodal attribute with logical matrix values
    !**********************************************************************
    subroutine nodal_attribute_matrix_logical_t_unload(self)
      implicit none
      class(NODAL_ATTRIBUTE_MATRIX_LOGICAL_T), intent(inout) :: self        !< Nodal attribute instance

      if(allocated(self%default_value))deallocate(self%default_value)
      if(allocated(self%values))deallocate(self%values)
      self%loaded = .false.

    end subroutine nodal_attribute_matrix_logical_t_unload

    !**********************************************************************
    ! PROCEDURE: nodal_attribute_matrix_real_t_read_header_from_file
    ! DESCRIPTION: Read header information from file for the nodal attribute with real matrix values
    !**********************************************************************
    subroutine nodal_attribute_matrix_real_t_read_header_from_file(self, io_unit, n_nodes)
      implicit none
      class(NODAL_ATTRIBUTE_MATRIX_REAL_T), intent(inout) :: self        !< Nodal attribute instance
      integer, intent(in) :: io_unit                                !< File unit number
      integer, intent(in) :: n_nodes                                !< Number of nodes
      integer :: i                                                  !< Loop index

      read(io_unit, *) self%attribute_name
      read(io_unit, *) self%units
      read(io_unit, *) self%n_values
      allocate(self%default_value(self%n_values))
      read(io_unit, *) (self%default_value(i), i = 1, self%n_values)

      allocate(self%values(n_nodes, self%n_values))
      do i = 1, n_nodes
        self%values(i, 1:self%n_values) = self%default_value(1:self%n_values)
      end do

    end subroutine nodal_attribute_matrix_real_t_read_header_from_file

    !**********************************************************************
    ! PROCEDURE: nodal_attribute_matrix_int_t_read_header_from_file
    ! DESCRIPTION: Read header information from file for the nodal attribute with integer matrix values
    !**********************************************************************
    subroutine nodal_attribute_matrix_int_t_read_header_from_file(self, io_unit, n_nodes)
      implicit none
      class(NODAL_ATTRIBUTE_MATRIX_INT_T), intent(inout) :: self        !< Nodal attribute instance
      integer, intent(in) :: io_unit                                !< File unit number
      integer, intent(in) :: n_nodes                                !< Number of nodes
      integer :: i                                                  !< Loop index

      read(io_unit, *) self%attribute_name
      read(io_unit, *) self%units
      read(io_unit, *) self%n_values
      allocate(self%default_value(self%n_values))
      read(io_unit, *) (self%default_value(i), i = 1, self%n_values)

      allocate(self%values(n_nodes, self%n_values))
      do i = 1, n_nodes
        self%values(i, 1:self%n_values) = self%default_value(1:self%n_values)
      end do

    end subroutine nodal_attribute_matrix_int_t_read_header_from_file

    !**********************************************************************
    ! PROCEDURE: nodal_attribute_matrix_logical_t_read_header_from_file
    ! DESCRIPTION: Read header information from file for the nodal attribute with logical matrix values
    !**********************************************************************
    subroutine nodal_attribute_matrix_logical_t_read_header_from_file(self, io_unit, n_nodes)
      implicit none
      class(NODAL_ATTRIBUTE_MATRIX_LOGICAL_T), intent(inout) :: self        !< Nodal attribute instance
      integer, intent(in) :: io_unit                                !< File unit number
      integer, intent(in) :: n_nodes                                !< Number of nodes
      integer :: i                                                  !< Loop index

      read(io_unit, *) self%attribute_name
      read(io_unit, *) self%units
      read(io_unit, *) self%n_values
      allocate(self%default_value(self%n_values))
      read(io_unit, *) (self%default_value(i), i = 1, self%n_values)

      allocate(self%values(n_nodes, self%n_values))
      do i = 1, n_nodes
        self%values(i, 1:self%n_values) = self%default_value(1:self%n_values)
      end do

    end subroutine nodal_attribute_matrix_logical_t_read_header_from_file

    !**********************************************************************
    ! PROCEDURE: nodal_attribute_matrix_real_t_read_body_from_file
    ! DESCRIPTION: Read body information from file for the nodal attribute with real matrix values
    !**********************************************************************
    subroutine nodal_attribute_matrix_real_t_read_body_from_file(self, io_unit)
      implicit none
      class(NODAL_ATTRIBUTE_MATRIX_REAL_T), intent(inout) :: self        !< Nodal attribute instance
      integer, intent(in) :: io_unit                                !< File unit number
      integer             :: n_non_default                          !< Number of non-default values
      integer             :: i, j                                   !< Loop index
      integer             :: node_number                            !< Node number
      character(len=MAX_NODAL_ATTRIBUTE_CHARACTER_LENGTH) :: attribute_name !< Name of the nodal attribute

      read(io_unit, *) attribute_name
      if(self%attribute_name /= trim(adjustl(attribute_name)))then
        write(*,'(A)') '[ERROR]: NODAL_ATTRIBUTE_MATRIX_T read from file with incorrect attribute name'
        stop
      end if

      read(io_unit, *) n_non_default
      do i = 1, n_non_default
        read(io_unit, *) node_number, (self%values(node_number, j), j = 1, self%n_values)
      end do
      self%loaded = .true.

    end subroutine nodal_attribute_matrix_real_t_read_body_from_file

    !**********************************************************************
    ! PROCEDURE: nodal_attribute_matrix_int_t_read_body_from_file
    ! DESCRIPTION: Read body information from file for the nodal attribute with integer matrix values
    !**********************************************************************
    subroutine nodal_attribute_matrix_int_t_read_body_from_file(self, io_unit)
      implicit none
      class(NODAL_ATTRIBUTE_MATRIX_INT_T), intent(inout) :: self        !< Nodal attribute instance
      integer, intent(in) :: io_unit                                !< File unit number
      integer             :: n_non_default                          !< Number of non-default values
      integer             :: i, j                                   !< Loop index
      integer             :: node_number                            !< Node number
      character(len=MAX_NODAL_ATTRIBUTE_CHARACTER_LENGTH) :: attribute_name !< Name of the nodal attribute

      read(io_unit, *) attribute_name
      if(self%attribute_name /= trim(adjustl(attribute_name)))then
        write(*,'(A)') '[ERROR]: NODAL_ATTRIBUTE_MATRIX_T read from file with incorrect attribute name'
        stop
      end if

      read(io_unit, *) n_non_default
      do i = 1, n_non_default
        read(io_unit, *) node_number, (self%values(node_number, j), j = 1, self%n_values)
      end do
      self%loaded = .true.

    end subroutine nodal_attribute_matrix_int_t_read_body_from_file

    !**********************************************************************
    ! PROCEDURE: nodal_attribute_matrix_logical_t_read_body_from_file
    ! DESCRIPTION: Read body information from file for the nodal attribute with logical matrix values
    !**********************************************************************
    subroutine nodal_attribute_matrix_logical_t_read_body_from_file(self, io_unit)
      implicit none
      class(NODAL_ATTRIBUTE_MATRIX_LOGICAL_T), intent(inout) :: self        !< Nodal attribute instance
      integer, intent(in) :: io_unit                                !< File unit number
      integer             :: n_non_default                          !< Number of non-default values
      integer             :: i, j                                   !< Loop index
      integer             :: node_number                            !< Node number
      character(len=MAX_NODAL_ATTRIBUTE_CHARACTER_LENGTH) :: attribute_name !< Name of the nodal attribute

      read(io_unit, *) attribute_name
      if(self%attribute_name /= trim(adjustl(attribute_name)))then
        write(*,'(A)') '[ERROR]: NODAL_ATTRIBUTE_MATRIX_T read from file with incorrect attribute name'
        stop
      end if

      read(io_unit, *) n_non_default
      do i = 1, n_non_default
        read(io_unit, *) node_number, (self%values(node_number, j), j = 1, self%n_values)
      end do
      self%loaded = .true.

    end subroutine nodal_attribute_matrix_logical_t_read_body_from_file

end module NODAL_ATTRIBUTE_BASE_TYPE_MODULE
