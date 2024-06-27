module math_act_prelu
    use data_types
    implicit none
    private
    public :: prelu
    
    interface prelu
        module procedure prelu_dp
    end interface prelu
    
    contains
    
    elemental function prelu_dp (x, alpha) result(res)
        real(dp),intent(in) :: x
        real(dp),intent(in) :: alpha
        real(dp)            :: res
        
        res = max(0.0_dp, x) + alpha * min(0.0_dp, x)
    end function prelu_dp
    


end module math_act_prelu