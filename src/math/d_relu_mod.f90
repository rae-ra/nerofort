module d_relu_mod
    use data_types
    implicit none
    private
    public :: d_relu
    
    interface d_relu
        module procedure d_relu_dp
    end interface d_relu
    
    contains
    
    elemental function d_relu_dp (x) result(res)
        real(dp),intent(in) :: x
        real(dp)            :: res
        
        res = sign(1.0_dp, x)
    end function d_relu_dp
    


end module d_relu_mod