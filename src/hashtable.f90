module hashtable
implicit none

private
public :: item
public :: dict
public :: close_dict
public :: add_item
public :: find
public :: print_dict

type :: item
    integer :: key
    integer :: value
    logical :: empty
    type(item), pointer :: next
end type item

contains

!!!!!!!!!!!!!!!!!!!!!!
!!!!! Public API !!!!!
!!!!!!!!!!!!!!!!!!!!!!

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

integer function find ( d, key )
    implicit none
    type(item), allocatable, target, intent(inout) :: d(:)
    integer, intent(in) :: key
    type(item), pointer :: current_item
    integer :: index

    ! Generate hash and calculate index
    index = modulo(key, size(d)) + 1

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

subroutine close_dict ( d )
    implicit none
    type(item), allocatable, target, intent(inout) :: d(:)
    integer :: i
    do i = 1, size(d)
        if ( associated(d(i)%next) ) call remove_item( d(i)%next )
    enddo
    deallocate( d )
end subroutine close_dict

subroutine add_item ( d, key, value )
    implicit none
    type(item), allocatable, target, intent(inout) :: d(:)
    integer, intent(in) :: key
    integer, intent(in) :: value
    integer :: index
    type(item), pointer :: current_item

    ! Calculate index
    index = modulo(key, size(d)) + 1

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

subroutine print_dict ( d )
    implicit none
    type(item), allocatable, target, intent(inout) :: d(:)
    type(item), pointer :: next_item
    integer :: i, depth
    do i = 1, size( d )
        depth = 1
        if ( d(i)%empty ) then
            write(*,100) repeat('*', 10), repeat('*', 10), repeat('#', depth) // repeat(' ', 10-depth)
        else
            write(*,200) d(i)%key, d(i)%value, repeat('#', depth) // repeat(' ', 10-depth)
        endif
        next_item => d(i)%next
        do
            depth = depth + 1
            if ( associated(next_item) ) then
                write(*,200) next_item%key, next_item%value, repeat('#', depth) // repeat(' ', 10-depth)
                next_item => next_item%next
            else
                exit
            endif
        enddo
        100 FORMAT (1X, A10, 4X, A10, 4X, A10)
        200 FORMAT (1X, I10, 4X, I10, 4X, A10)
    enddo

end subroutine print_dict



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Private subroutines and functions !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

recursive subroutine remove_item ( i )
    implicit none
    type(item), pointer, intent(inout) :: i
    if ( associated(i%next) ) call remove_item( i%next )
    deallocate( i )
end subroutine remove_item

end module hashtable
