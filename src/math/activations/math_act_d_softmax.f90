module math_act_d_softmax
    use data_types
    use math_act_softmax
    implicit none
    private
    public :: d_softmax
    
    
    interface d_softmax
        module procedure d_softmax_dp
    end interface d_softmax
    

    
    
    contains
    
    function d_softmax_dp(x) result(res)
        real(dp), intent(in)                    :: x(:,:)
        real(dp), dimension(size(x,1),size(x,2)):: res
        integer                                 :: m, d
        real(dp), dimension(:,:), allocatable   :: a
        real(dp), dimension(:,:,:), allocatable :: tensor1, tensor2
        integer                                 :: i, j, k

        ! Get dimensions of x
        m = size(x, 1)
        d = size(x, 2)

        ! Allocate arrays
        allocate(a(m, d), tensor1(m, d, d), tensor2(m, d, d))

        ! Reshape input matrix if necessary
        if (m == 1) then
            a = reshape(x, (/1, d/))
        else
            a = x
        end if
        ! Compute softmax
        a = softmax(a)

        ! Compute tensor1
        do i = 1, m
            do j = 1, d
                do k = 1, d
                    tensor1(i, j, k) = a(i, j) * a(i, k)
                end do
            end do
        end do

        ! Compute tensor2
        do i = 1, m
            do j = 1, d
                do k = 1, d
                    if (j == k) then
                        tensor2(i, j, k) = a(i, j) * (1.0d0 - a(i, j))
                    else
                        tensor2(i, j, k) = -a(i, j) * a(i, k)
                    end if
                end do
            end do
        end do

        ! Compute result
        do i = 1, m
            do j = 1, d
                do k = 1, d
                    res(i, j) = tensor2(i, j, k) - tensor1(i, j, k)
                end do
            end do
        end do

    end function d_softmax_dp

end module math_act_d_softmax