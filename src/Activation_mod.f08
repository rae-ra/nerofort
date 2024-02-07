module Activation_mod
    use Math_util
    implicit none
    private
    public :: create_activation




    !TODO re-add support for softmax

    !> alpha: slope parameter
    type, public :: Activation
        character(len=10) :: activation_type
        real(dp) :: alpha = 0.2
        real(dp), pointer :: cache_x
        real(dp), pointer :: cache_x_1 (:)
        real(dp), pointer :: cache_x_2 (:,:)
        real(dp), pointer :: cache_x_3 (:,:,:)
        real(dp), pointer :: cache_x_4 (:,:,:,:)
    contains
        procedure :: linear => activ_linear
        procedure :: d_linear => activ_d_linear
        procedure :: print => activation_print
        procedure :: get_activation => get_activation
        procedure :: get_d_activation => get_d_activation
        procedure :: forward => gen_forward_4, gen_forward_3, gen_forward_2, &
            gen_forward_1
        procedure :: backpropagation => gen_back_4, gen_back_3, gen_back_2, &
            gen_back_1
        procedure :: set_alpha => set_alpha
    end type Activation


    interface Activation
        module procedure init_0, init_alpha
    end interface Activation



contains

    function create_activation(activation_type) result(activ_ptr)
        character(len=10), intent(in) :: activation_type
        type(Activation), pointer :: activ_ptr

        activ_ptr = Activation(activation_type)
    end function create_activation

    function init_0(activation_type) result(activ)
        character(len=10), intent(in) :: activation_type

        type(Activation) :: activ

        activ%activation_type = activation_type

        call activ%set_alpha(real(0.2,dp))

    end function init_0

    function init_alpha(activation_type, alpha) result(activ)
        character(len=10), intent(in) :: activation_type
        real(dp) :: alpha

        type(Activation) :: activ

        activ%activation_type = activation_type

        call activ%set_alpha(alpha)

    end function init_alpha


    subroutine set_alpha(this, alpha)
        class(Activation), intent(inout) :: this
        real(dp), intent(in) :: alpha

        this%alpha = alpha

    end subroutine set_alpha


    pure function activ_linear(this, x) result(y)
        class(Activation), intent(in) :: this
        real(dp), intent(in) :: x
        real(dp) :: y

        y = x
    end function activ_linear

    pure function activ_d_linear(this, x) result(dy)
        class(Activation), intent(in) :: this
        real(dp), intent(in) :: x
        real(dp) :: dy

        dy = 1.0
    end function activ_d_linear



    !> not currently used
!    subroutine softmax(a, m, d)
!        real(dp), intent(inout) :: a(:,:)
!        integer, intent(in) :: m, d
!        real(dp) :: exp_sum, max_val
!        integer :: i, j
!
!         Find the maximum value in each row
!        do i = 1, m
!            max_val = maxval(a(i, 1:d))
!             Compute softmax for each element in the row
!            exp_sum = 0.0d0
!            do j = 1, d
!                a(i, j) = exp(a(i, j) - max_val)
!                exp_sum = exp_sum + a(i, j)
!            end do
!             Normalize the row to obtain probabilities
!            do j = 1, d
!                a(i, j) = a(i, j) / exp_sum
!            end do
!        end do
!
!    end subroutine softmax


    subroutine activation_print(this)
        class(Activation), intent(in) :: this
        print *, 'Activation Type: ', this%activation_type
    end subroutine activation_print

    pure function get_activation(this, x) result(y)
        class(Activation), intent(in) :: this
        real(dp), intent(in) :: x
        real(dp) :: y

        select case (this%activation_type)
            case ('linear')
                y = this%linear(x)
            case ('sigmoid')
                y = sigmoid(x)
            case ('tanh')
                y = tanh(x)
            case ('relu')
                y = relu(x)
            case ('prelu')
                y = prelu(x, this%alpha)
            !case ('softmax')
            !    y = softmax(x)

        end select
    end function get_activation

    elemental function get_d_activation(this, x) result(dy)
        class(Activation), intent(in) :: this
        real(dp), intent(in) :: x
        real(dp) :: dy

        select case (this%activation_type)
            case ('linear')
                dy = this%d_linear(x)
            case ('sigmoid')
                dy = d_sigmoid(x)
            case ('tanh')
                dy = d_tanh(x)
            case ('relu')
                dy = d_relu(x)
            case ('prelu')
                dy = d_prelu(x, this%alpha)
            !case ('softmax')
            !    dy = d_softmax(x)

        end select
    end function get_d_activation

    function gen_forward_4(this,X) result(z)
        class(Activation), intent(in) :: this
        real(dp), intent(in) :: X(:,:,:,:)
        real(dp), allocatable :: z(:,:,:,:)

        this%cache_x_4 = X

        z = pure_forward(this, X)
    end function gen_forward_4


    function gen_back_4(this, dz) result(dx)
        class(Activation), intent(in) :: this
        real(dp), intent(in) :: dz (:,:,:,:)
        real(dp), allocatable :: dx (:,:,:,:)

        dx = pure_backpropagation(this,dz, this%cache_x_4)
    end function gen_back_4


    function gen_forward_3(this,X) result(z)
        class(Activation), intent(in) :: this
        real(dp), intent(in) :: X(:,:,:)
        real(dp), allocatable :: z(:,:,:)

        this%cache_x_3 = X

        z = pure_forward(this, X)
    end function gen_forward_3


    function gen_back_3(this, dz) result(dx)
        class(Activation), intent(in) :: this
        real(dp), intent(in) :: dz (:,:,:)
        real(dp), allocatable :: dx (:,:,:)

        dx = pure_backpropagation(this,dz, this%cache_x_3)
    end function gen_back_3


    function gen_forward_2(this,X) result(z)
        class(Activation), intent(in) :: this
        real(dp), intent(in) :: X(:,:)
        real(dp), allocatable :: z(:,:)

        this%cache_x_2 = X

        z = pure_forward(this, X)
    end function gen_forward_2


    function gen_back_2(this, dz) result(dx)
        class(Activation), intent(in) :: this
        real(dp), intent(in) :: dz (:,:)
        real(dp), allocatable :: dx (:,:)

        dx = pure_backpropagation(this, dz, this%cache_x_2)
    end function gen_back_2


    function gen_forward_1(this,X) result(z)
        class(Activation), intent(in) :: this
        real(dp), intent(in) :: X(:)
        real(dp), allocatable :: z(:)

        this%cache_x_1 = X

        z = pure_forward(this, X)
    end function gen_forward_1


    function gen_back_1(this, dz) result(dx)
        class(Activation), intent(in) :: this
        real(dp), intent(in) :: dz (:)
        real(dp), allocatable :: dx (:)

        dx = pure_backpropagation(this, dz, this%cache_x_1)
    end function gen_back_1


    elemental function pure_forward(this, X) result(z)
        class(Activation), intent(in) :: this
        real(dp), intent(in) :: X
        real(dp) :: z

        z = this%get_activation(X)
    end function pure_forward


    elemental function pure_backpropagation(this, dz, cache_x) result(dx)
        class(Activation), intent(in) :: this
        real(dp), intent(in) :: dz
        real(dp), intent(in) :: cache_x
        real(dp) :: dx

        real(dp) :: f_prime

        f_prime = this%get_d_activation(cache_x)

        dx = dz * f_prime
    end function pure_backpropagation

end module Activation_mod