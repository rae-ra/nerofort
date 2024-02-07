module Padding2D_mod
    use Math_UTIL, only: dp
    implicit none

    private
    public :: Padding_init, Padding_get_dimensions, Padding_forward, Padding_backpropagation

    type Padding2D
        character(len=:), allocatable :: ptype
        integer :: pt, pb, pl, pr
        integer :: m, Nc, Nh, Nw
        integer :: output_shape(4)
        real, allocatable :: zeros_r(:,:,:,:), zeros_l(:,:,:,:)
        real, allocatable :: zeros_t(:,:,:,:), zeros_b(:,:,:,:)
    end type Padding2D

    interface Padding_forward
        module procedure :: Padding_forward_3, Padding_forward_4
    end interface Padding_forward

    interface Padding_backpropagation
        module procedure :: Padding_backpropagation_3, &
            Padding_backpropagation_4
    end interface Padding_backpropagation


    contains

    subroutine Padding_init(this, p)
        class(Padding2D), intent(inout) :: this
        character(len=*), intent(in), optional :: p

        if (present(p)) then
            this%ptype = p
        else
            this%ptype = "valid"
        end if
    end subroutine Padding_init

    subroutine Padding_get_dimensions(this, input_shape, kernel_size, s)
        class(Padding2D), intent(inout) :: this
        integer, intent(in) :: input_shape(:), kernel_size(2), s(2)
        integer :: Kh, Kw, sh, sw

        if (size(input_shape) == 4) then
            this%m = input_shape(1)
            this%Nc = input_shape(2)
            this%Nh = input_shape(3)
            this%Nw = input_shape(4)
        elseif (size(input_shape) == 3) then
            this%m = 1
            this%Nc = input_shape(1)
            this%Nh = input_shape(2)
            this%Nw = input_shape(3)
        end if

        Kh = kernel_size(1)
        Kw = kernel_size(2)
        sh = s(1)
        sw = s(2)

        select case(this%ptype)
            case('valid')
                this%pt = 0
                this%pb = 0
                this%pl = 0
                this%pr = 0
            case('same')
                this%pt = (sh - 1) * this%Nh + Kh - sh
                this%pb = this%pt + mod(this%pt + 1, 2)
                this%pl = (sw - 1) * this%Nw + Kw - sw
                this%pr = this%pl + mod(this%pl + 1, 2)
            case default
                this%pt = 0
                this%pb = 0
                this%pl = 0
                this%pr = 0
        end select

        if (size(input_shape) == 4) then
            this%output_shape = [this%m, this%Nc, this%Nh + this%pt + this%pb, &
                  this%Nw + this%pl + this%pr]
        elseif (size(input_shape) == 3) then
            this%output_shape = [this%Nc, this%Nh + this%pt + this%pb, &
                this%Nw + this%pl + this%pr, 1]
        end if
    end subroutine Padding_get_dimensions

    subroutine Padding_forward_4(this, X, kernel_size, s, Xp)
        class(Padding2D), intent(inout) :: this
        real(dp), intent(in) :: X(:,:,:,:)
        integer, intent(in) :: kernel_size(2), s(2)
        real(dp), allocatable, intent(out) :: Xp(:,:,:,:)
        integer :: m, Nc, Nh, Nw, l, k, padded_Nh, padded_Nw

        call Padding_get_dimensions(this, shape(X),kernel_size,s)

        m = this%m
        Nc = this%Nc
        Nh = this%Nh
        Nw = this%Nw

        padded_Nh = Nh + this%pt + this%pb
        padded_Nw = Nw + this%pl + this%pr

        ! Padding top and bottom
        Xp(:, :, 1:this%pt, :) = 0.0
        Xp(:, :, padded_Nh-this%pb+1:padded_Nh, :) = 0.0
        Xp(:, :, this%pt+1:padded_Nh-this%pb, :) = X

        ! Padding left and right
        do l = 1, Nc
            do k = 1, padded_Nh
                Xp(:, l, k, 1:this%pl) = 0.0
                Xp(:, l, k, padded_Nw-this%pr+1:padded_Nw) = 0.0
                Xp(:, l, k, this%pl+1:padded_Nw-this%pr) = &
                    Xp(:, l, k, this%pl+1:padded_Nw-this%pr) &
                                            + X(:, l, k, :)
            end do
        end do

    end subroutine Padding_forward_4

    subroutine Padding_forward_3(this, X, kernel_size, s, Xp)
        class(Padding2D), intent(inout) :: this
        real(dp), intent(in) :: X(:,:,:)
        integer, intent(in) :: kernel_size(2), s(2)
        real(dp), allocatable, intent(out) :: Xp(:,:,:)
        integer :: m, Nc, Nh, Nw, l, k, padded_Nh, padded_Nw

        call Padding_get_dimensions(this, shape(X),kernel_size,s)

        m = this%m
        Nc = this%Nc
        Nh = this%Nh
        Nw = this%Nw

        padded_Nh = Nh + this%pt + this%pb
        padded_Nw = Nw + this%pl + this%pr

        ! Padding top and bottom
        Xp(:, :, 1:this%pt) = 0.0
        Xp(:, :, padded_Nh-this%pb+1:padded_Nh) = 0.0
        Xp(:, :, this%pt+1:padded_Nh-this%pb) = X

        ! Padding left and right
        do l = 1, Nc
            do k = 1, padded_Nh
                Xp( l, k, 1:this%pl) = 0.0
                Xp( l, k, padded_Nw-this%pr+1:padded_Nw) = 0.0
                Xp( l, k, this%pl+1:padded_Nw-this%pr) = &
                    Xp( l, k, this%pl+1:padded_Nw-this%pr) &
                                            + X( l, k, :)
            end do
        end do

    end subroutine Padding_forward_3

    function Padding_backpropagation_4(this, dXp) result (dX)
        class(Padding2D), intent(in) :: this
        real(dp), intent(in) :: dXp(:,:,:,:)
        real(dp), allocatable :: dX(:,:,:,:)
        integer :: m, Nc, Nh, Nw

        m = this%m
        Nc = this%Nc
        Nh = this%Nh
        Nw = this%Nw

        dX = dXp(:, :, this%pt+1:this%pt+Nh, this%pl+1:this%pl+Nw)

    end function Padding_backpropagation_4

    function Padding_backpropagation_3(this, dXp) result (dX)
        class(Padding2D), intent(in) :: this
        real(dp), intent(in) :: dXp(:,:,:)
        real(dp), allocatable :: dX(:,:,:)
        integer :: m, Nc, Nh, Nw

        m = this%m
        Nc = this%Nc
        Nh = this%Nh
        Nw = this%Nw

        dX = dXp( :, this%pt+1:this%pt+Nh, this%pl+1:this%pl+Nw)

    end function Padding_backpropagation_3

    subroutine zeros_array(output_shape, shape_arr, zeros_arr)
        integer, intent(in) :: shape_arr(:)
        integer, intent(out) :: output_shape(:)
        real(dp), allocatable, intent(out) :: zeros_arr(:,:,:,:)

        output_shape = shape_arr
        allocate(zeros_arr(output_shape(1), output_shape(2),output_shape(3), &
            output_shape(4)))

        zeros_arr = 0.0
    end subroutine zeros_array

end module Padding2D_mod
