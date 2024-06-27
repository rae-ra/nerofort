module nerofort_learning_rate_decay
    use nerofort_math, only: dp
    implicit none
    private
    public :: constant, time_decay, step_decay, expon_decay

    interface constant
        module procedure constant_dp
    end interface constant

    interface time_decay
        module procedure time_decay_dp
    end interface time_decay

    interface step_decay
        module procedure step_decay_dp
    end interface step_decay

    interface expon_decay
        module procedure expon_decay_dp
    end interface expon_decay

contains

    !TEST needed
    elemental function constant_dp(t, lr_0) result(res)
        integer, intent(in) :: t
        real(dp), intent(in) :: lr_0
        real(dp) :: res

        res = lr_0
    end function constant_dp

    !TEST needed
    elemental function time_decay_dp(t, lr_0, k) result(res)
        integer, intent(in) :: t
        real(dp), intent(in) :: lr_0, k
        real(dp) :: res

        res = lr_0 / (1.0 + k * real(t))
    end function time_decay_dp

    !TEST needed
    elemental function step_decay_dp( t, lr_0, F, D) result(res)
        integer, intent(in) :: t, D
        real(dp), intent(in) :: lr_0, F
        real(dp) :: res
        real(dp) :: mult
        mult = F ** real(floor((1.0 + t) / D))
        res = lr_0 * mult
    end function step_decay_dp

    !TEST needed
    elemental function expon_decay_dp( t, lr_0, k) result(res)
        integer, intent(in) :: t
        real(dp), intent(in) :: lr_0, k
        real(dp) :: res

        res = lr_0 * exp(-k * real(t))
    end function expon_decay_dp

end module nerofort_learning_rate_decay
