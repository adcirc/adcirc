C
C                      ADCIRC - HASHTABLE MODULE
C
C    ========================================================================
C    |                                                                      |
C    |   This file contains the subroutines required to build a simple      |
C    |   hashtable that maps from integer keys to integer values.           |
C    |   Hash values are calculated using a simple modulus operater and     |
C    |   collisions are handled using a simple linked list.                 |
C    |                                                                      |
C    |   Written by Tristan Dyer, atdyer@ncsu.edu                           |
C    |   North Carolina State University,                                   |
C    |   2017                                                               |
C    |                                                                      |
C    ========================================================================
      module hashtable
      implicit none

      private
      public :: item
      public :: dict
      public :: close_dict
      public :: add_item
      public :: find

C     Define the new item type. Key/value pairs are stored as items, which
C     are also used to implement the linked list.
      type :: item
          integer :: key
          integer :: value
          logical :: empty
          type(item), pointer :: next
      end type item

      contains

C     This subroutine is used to allocate and initialize an empty dictionary.
C       item_array: An unallocated array of items
C       num_items: The number of items to be placed into the dictionary
      subroutine dict ( item_array, num_items )
          implicit none
          type(item), allocatable, intent(inout) :: item_array(:)
          integer, intent(in) :: num_items
          integer :: i

          allocate( item_array(num_items) )
          do i = 1, num_items
              item_array(i)%key = 0
              item_array(i)%value = 0
              item_array(i)%empty = .true.
              nullify(item_array(i)%next)
          enddo

      end subroutine dict

C     Add a key/value pair to the dictionary
C       d: The dictionary to add to
C       key: The key
C       value: The value
      subroutine add_item ( d, key, value )
          implicit none
          type(item), allocatable, target, intent(inout) :: d(:)
          integer, intent(in) :: key
          integer, intent(in) :: value
          integer :: index
          type(item), pointer :: current_item

          ! Calculate index
          index = 1 + modulo(key-1, size(d))

          ! Check for collisions
          if ( d(index)%empty ) then
              d(index)%key = key
              d(index)%value = value
              d(index)%empty = .false.
          else
              current_item => d(index)
              do
                  if ( .not. associated( current_item%next )) then
                      allocate( current_item%next )
                      current_item => current_item%next
                      current_item%key = key
                      current_item%value = value
                      current_item%empty = .false.
                      nullify( current_item%next )
                      exit
                  else
                      current_item => current_item%next
                  endif
              enddo
          endif

      end subroutine add_item

C     Find and return the value associated with a given key
C       d: The dictionary on which to perform the search
C       key: The key to search for
      integer function find ( d, key )
          implicit none
          type(item), allocatable, target, intent(inout) :: d(:)
          integer, intent(in) :: key
          type(item), pointer :: current_item
          integer :: index

          ! Generate hash and calculate index
          index = 1 + modulo(key-1, size(d))

          ! Find the item
          current_item => d(index)
          do
              if ( key == current_item%key ) then
                  find = current_item%value
                  exit
              else
                  if ( associated( current_item%next ) ) then
                      current_item => current_item%next
                  else
                      find = 0
                      exit
                  endif
              endif
          enddo

      end function find

C     Call this subroutine when finished using the dictionary. It
C     ensures that all memory is propery deallocated
C       d: The dictionary to deallocate
      subroutine close_dict ( d )
          implicit none
          type(item), allocatable, target, intent(inout) :: d(:)
          integer :: i
          do i = 1, size(d)
              if ( associated(d(i)%next) ) call remove_item( d(i)%next )
          enddo
          deallocate( d )
      end subroutine close_dict

C     Private subroutine used to recursively deallocate a linked list
      recursive subroutine remove_item ( i )
          implicit none
          type(item), pointer, intent(inout) :: i
          if ( associated(i%next) ) call remove_item( i%next )
          deallocate( i )
      end subroutine remove_item

      end module hashtable
