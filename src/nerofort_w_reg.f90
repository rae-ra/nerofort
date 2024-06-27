module nerofort_w_reg
    use nerofort_math, only: dp
    implicit none
    private
    !public :: make_wreg

    type, public :: W_reg
        character(len=10), allocatable, public :: str
        real(dp), public :: num


    end type W_reg

end module nerofort_w_reg