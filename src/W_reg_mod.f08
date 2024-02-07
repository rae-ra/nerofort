module W_reg_mod
    use math_util, only: dp
    implicit none
    private
    !public :: make_wreg

    type, public :: W_reg
        character(len=10), allocatable, public :: str
        real(dp), public :: num


    end type W_reg

end module W_reg_mod