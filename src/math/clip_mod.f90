module clip_mod
    use data_types
    implicit none
    private
    public :: clip
    
    interface clip
        module procedure clip_dp
    end interface clip
    
    contains
    
    
    elemental function clip_dp (x, lower, upper) result(res) 
        real(dp), intent(in) :: x
        real(dp), intent(in) :: lower
        real(dp), intent(in) :: upper
        real(dp)             :: res
        
        res = max(min(x, upper), lower)
        
    end function clip_dp
    


end module clip_mod