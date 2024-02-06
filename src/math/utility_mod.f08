module utility_mod
    use data_types, only: dp
    use str_mod
    implicit none
    public :: unique_str, mean, std_dev, rand_perm
contains

    function unique_str(arr) result(unique_values_result)
        type(str), dimension(:), intent(in) :: arr
        type(str), allocatable :: unique_values(:)
        type(str), allocatable :: unique_values_result(:)
        integer :: i, j, count

        ! Count unique values
        count = 0
        do i = 1, size(arr)
            if (count == 0) then
                count = 1
                allocate(unique_values(count))
                unique_values(1) = arr(i)
            else
                do j = 1, count
                    if (arr(i) == unique_values(j)) exit
                end do
                if (j > count) then
                    count = count + 1
                    allocate(unique_values(1:count))
                    unique_values(count) = arr(i)
                end if
            end if
        end do

        ! Allocate and copy unique values to the result
        allocate(unique_values_result(1:count))
        unique_values_result(1:count) = unique_values(1:count)

    end function unique_str


    function mean(arr, dim) result(mean_value)
        real(dp), dimension(:), intent(in) :: arr
        integer, intent(in), optional :: dim
        real(dp) :: mean_value
        integer :: i, n

        n = size(arr)
        mean_value = 0.0_dp

        if (present(dim)) then
            ! Calculate mean along specified dimension
            do i = 1, n
                mean_value = mean_value + arr(i)**2
            end do
            mean_value = mean_value / real(n)
        else
            ! Calculate mean of entire array
            do i = 1, n
                mean_value = mean_value + arr(i)
            end do
            mean_value = mean_value / real(n)
        end if

    end function mean

    function std_dev(arr, dim) result(std_dev_value)
        real(dp), dimension(:), intent(in) :: arr
        integer, intent(in), optional :: dim
        real(dp) :: std_dev_value
        real(dp) :: mean_value
        integer :: i, n

        n = size(arr)
        mean_value = mean(arr, dim)
        std_dev_value = 0.0_dp

        if (present(dim)) then
            ! Calculate standard deviation along specified dimension
            do i = 1, n
                std_dev_value = std_dev_value + (arr(i) - mean_value)**2
            end do
            std_dev_value = sqrt(std_dev_value / real(n))
        else
            ! Calculate standard deviation of entire array
            do i = 1, n
                std_dev_value = std_dev_value + (arr(i) - mean_value)**2
            end do
            std_dev_value = sqrt(std_dev_value / real(n))
        end if

    end function std_dev

    !
    !> random permutation on array arr from 1 to n
    subroutine rand_perm(n, arr)
        integer, intent(in) :: n
        integer, dimension(:), intent(out) :: arr
        integer :: i, j, temp
        real(dp) :: r


        ! Initialize array with consecutive integers
        arr = [(i, i=1,n)]

        ! Perform Fisher-Yates shuffle
        do i = n, 2, -1
            call random_number(r)
            j = int(i * r)
            temp = arr(i)
            arr(i) = arr(j)
            arr(j) = temp
        end do

    end subroutine rand_perm

end module utility_mod
