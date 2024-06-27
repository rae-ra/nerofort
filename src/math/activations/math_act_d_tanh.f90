module math_act_d_tanh
    use data_types
    implicit none
    private
    public d_tanh
    
    
    interface d_tanh
        module procedure d_tanh_dp
    end interface d_tanh
    
    contains
    
    elemental function d_tanh_dp (x) result(res) 
        real(dp), intent(in) :: x
        real(dp)             ::res
        
        res = 1.0 - tanh(x)**2
        
    end function d_tanh_dp

end module math_act_d_tanh