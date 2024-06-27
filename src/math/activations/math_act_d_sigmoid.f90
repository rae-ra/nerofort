module math_act_d_sigmoid
    use data_types
    use math_act_sigmoid
    implicit none
    private
    public :: d_sigmoid
    
    interface d_sigmoid
        module procedure d_sigmoid_dp
    end interface d_sigmoid
    
    contains
    
    elemental function d_sigmoid_dp (x) result(res)
        real(dp), intent(in) :: x
        real(dp)             :: res
        real(dp)             :: sig_x
        
        sig_x = sigmoid(x)
        
        res = sig_x * (1.0_dp - sig_x)
    end function d_sigmoid_dp


end module math_act_d_sigmoid