module d_prelu_mod
    use data_types
    implicit none
    private
    public :: d_prelu

    interface d_prelu
        module procedure d_prelu_dp
    end interface d_prelu

    contains

    elemental function d_prelu_dp (x, alpha) result(res)
        real(dp),intent(in) :: x
        real(dp),intent(in) :: alpha
        real(dp)            :: res

        res = sign(1.0_dp, x) + (1.0_dp - sign(1.0_dp, x)) * alpha
    end function d_prelu_dp



end module d_prelu_mod