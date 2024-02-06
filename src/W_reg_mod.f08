module W_reg_mod
    use math_util, only: dp
    implicit none
    private
    !public :: make_wreg

    type, public :: W_reg
        character(len=10), allocatable, public :: str
        real(dp), public :: num


    end type W_reg

contains
!
!    function make_wreg (str, num) result(wreg)
!        character(len=10), allocatable,intent(in) :: str
!        real(dp), intent(in) :: num
!        type(W_reg) :: wreg
!        wreg%str = str
!        wreg%num = num
!
!    end function make_wreg
end module W_reg_mod