module dict_mod
    use str_mod
    !use Math_util, only: unique_str
    implicit none

    private
    public :: dict, init_dict, set, get

    integer, parameter :: DEFAULT_DICT_ALLOC = 127

    type :: dict
        type(str), allocatable :: keys(:)
        integer, dimension(:), allocatable :: vs(:)
        integer :: size
        contains
        procedure :: set => set
        procedure :: get => get
    end type dict

contains

    subroutine init_dict(d, size)
        type(dict) :: d
        integer, intent(in), optional::size
        if (present(size)) then
            d%size = size
        else
            d%size = DEFAULT_DICT_ALLOC
        endif
        allocate(d%keys(d%size))
        allocate( d%vs(d%size))
        d%size = 0
    end subroutine init_dict

    subroutine set(d, key, v)
        Class(dict), intent(inout) :: d
        type(str), intent(in) :: key
        integer, intent(in) :: v
        integer :: i

        if (.not. allocated(d%keys)) then
            allocate(d%keys(0:0), d%vs(0:0))
            d%size = 0
        end if

        ! Check if key already exists
        do i = 1, d%size
            if (d%keys(i) * key) then
                d%vs(i) = v
                return
            end if
        end do

        ! Key not found, add it to the dictionary
        d%size = d%size + 1
        allocate(d%keys(0:d%size), d%vs(0:d%size))
        d%keys(d%size) = key
        d%vs(d%size) = v
    end subroutine set

    function get(d, key) result(v)
        Class(dict), intent(in) :: d
        type(str), intent(in) :: key
        integer :: v, i

        v = 0

        do i = 1, d%size
            if (d%keys(i) * key) then
                v = d%vs(i)
                return
            end if
        end do
    end function get

!    subroutine label_encoding(Y, result)
!        integer, intent(in) :: Y(:,:)
!        integer, intent(out) :: result(:,:)
!
!        type(dict), dimension(:), allocatable :: idx_list
!        integer :: col, m, d, idx
!        type(str), allocatable :: val
!
!        m = size(Y, 1)
!        d = size(Y, 2)
!        allocate(idx_list(d))
!
!        do col = 1, d
!            call init_dict(idx_list(col))
!
!            do idx = 1, m
!                val = unique_str(Y(:, col))
!                call set(idx_list(col), val, idx)
!            end do
!
!            do idx = 1, m
!                result(idx, col) = get(idx_list(col), unique_str(Y(idx, &
!                    col))%trim())
!            end do
!        end do
!    end subroutine label_encoding



end module dict_mod
