module sigmoid_mod
    use data_types
    implicit none
    private
    public :: sigmoid
    
    interface sigmoid
        module procedure sigmoid_dp
    end interface sigmoid
    
    contains
    
    elemental function sigmoid_dp (x) result(res)
        real(dp),intent(in) :: x
        real(dp)            :: res
        
        res = 1.0 / (1.0 + exp(-x))
    end function sigmoid_dp
    


end module sigmoid_mod