module softmax_mod
    use data_types
    implicit none
    private
    public ::  softmax
    
    
    interface softmax
        module procedure softmax_dp
    end interface softmax
    
    
    contains
    
    function softmax_dp(x) result(res)
        real(dp), dimension(:,:), intent(in)        :: x
        real(dp), dimension(size(x,1), size(x,2))   :: res

        real(dp), dimension(size(x,1), size(x,2))   :: z
        real(dp), dimension(size(x,1))              :: exp_sum
        integer                                     :: i
        
        !TODO implement more efficient
        
        do i = 1, size(x, 1)
            z(i, :) = x(i, :) - maxval(x(i, :))
        end do

        exp_sum = sum(exp(z), dim=2)

        do i = 1, size(x, 1)
            res(i, :) = exp(z(i, :)) / exp_sum(i)
        end do
        
    end function softmax_dp


end module softmax_mod