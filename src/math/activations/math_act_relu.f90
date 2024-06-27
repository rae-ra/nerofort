module math_act_relu
    use data_types
    implicit none
    private
    public :: relu
    
    interface relu
        module procedure relu_dp
    end interface relu
    
    contains
    
    elemental function relu_dp (x) result(res)
        real(dp),intent(in) :: x
        real(dp)            :: res
        
        res = max(0.0_dp, x)
    end function relu_dp
    


end module math_act_relu