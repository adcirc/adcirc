!-------------------------------------------------------------------------------!
!
! ADCIRC - The ADvanced CIRCulation model
! Copyright (C) 1994-2025 R.A. Luettich, Jr., J.J. Westerink
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!-------------------------------------------------------------------------------!
!
!                      ADCIRC - HASHTABLE MODULE
!
!    ========================================================================
!    |                                                                      |
!    |   This file contains the subroutines required to build a simple      |
!    |   hashtable that maps from integer keys to integer values.           |
!    |   Hash values are calculated using a simple modulus operater and     |
!    |   collisions are handled using a simple linked list.                 |
!    |                                                                      |
!    |   Written by Tristan Dyer, atdyer@ncsu.edu                           |
!    |   North Carolina State University,                                   |
!    |   2017                                                               |
!    |                                                                      |
!    ========================================================================
module hashtable
   implicit none

   private
   public :: ipair
   public :: dict
   public :: close_dict
   public :: add_ipair
   public :: find

   ! Define the new ipair (integer pair) type. Key/value pairs are stored
   ! as ipairs, which are also used to implement the linked list.
   type :: ipair
      integer :: key
      integer :: value
      logical :: empty
      type(ipair), pointer :: next
   end type ipair

contains

   ! This subroutine is used to allocate and initialize an empty dictionary.
   !   ipair_array: An unallocated array of ipairs
   !   num_ipairs: The number of ipairs to be placed into the dictionary
   subroutine dict(ipair_array, num_ipairs)
      implicit none
      type(ipair), allocatable, intent(inout) :: ipair_array(:)
      integer, intent(in) :: num_ipairs
      integer :: i

      allocate (ipair_array(num_ipairs))
      do i = 1, num_ipairs
         ipair_array(i)%key = 0
         ipair_array(i)%value = 0
         ipair_array(i)%empty = .true.
         nullify (ipair_array(i)%next)
      end do

   end subroutine dict

   ! Add a key/value pair to the dictionary
   !   d: The dictionary to add to
   !   key: The key
   !   value: The value
   subroutine add_ipair(d, key, value)
      implicit none
      type(ipair), allocatable, target, intent(inout) :: d(:)
      integer, intent(in) :: key
      integer, intent(in) :: value
      integer :: index
      type(ipair), pointer :: current_ipair

      ! Calculate index
      index = 1 + modulo(key - 1, size(d))

      ! Check for collisions
      if (d(index)%empty) then
         d(index)%key = key
         d(index)%value = value
         d(index)%empty = .false.
      else
         current_ipair => d(index)
         do
            if (.not. associated(current_ipair%next)) then
               allocate (current_ipair%next)
               current_ipair => current_ipair%next
               current_ipair%key = key
               current_ipair%value = value
               current_ipair%empty = .false.
               nullify (current_ipair%next)
               exit
            else
               current_ipair => current_ipair%next
            end if
         end do
      end if

   end subroutine add_ipair

   ! Find and return the value associated with a given key
   !   d: The dictionary on which to perform the search
   !   key: The key to search for
   integer function find(d, key)
      implicit none
      type(ipair), allocatable, target, intent(inout) :: d(:)
      integer, intent(in) :: key
      type(ipair), pointer :: current_ipair
      integer :: index

      ! Generate hash and calculate index
      index = 1 + modulo(key - 1, size(d))

      ! Find the ipair
      current_ipair => d(index)
      do
         if (key == current_ipair%key) then
            find = current_ipair%value
            exit
         else
            if (associated(current_ipair%next)) then
               current_ipair => current_ipair%next
            else
               find = 0
               exit
            end if
         end if
      end do

   end function find

   ! Call this subroutine when finished using the dictionary. It
   ! ensures that all memory is propery deallocated
   !   d: The dictionary to deallocate
   subroutine close_dict(d)
      implicit none
      type(ipair), allocatable, target, intent(inout) :: d(:)
      integer :: i
      do i = 1, size(d)
         if (associated(d(i)%next)) call remove_ipair(d(i)%next)
      end do
      deallocate (d)
   end subroutine close_dict

   ! Private subroutine used to recursively deallocate a linked list
   recursive subroutine remove_ipair(i)
      implicit none
      type(ipair), pointer, intent(inout) :: i
      if (associated(i%next)) call remove_ipair(i%next)
      deallocate (i)
   end subroutine remove_ipair

end module hashtable
