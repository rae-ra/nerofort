module Activation_mod
    implicit none
    private
    public :: create_activation
    
    
    integer, parameter, public :: sp = selected_real_kind(6, 37)
    integer, parameter, public :: dp = selected_real_kind(15, 307)
    integer, parameter, public :: qp = selected_real_kind(33, 4931)
    
    real(dp) :: one_dp = 1.0
    
    real(dp), pointer :: cache_x(:,:)
    
    
    
    !> alpha: slope parameter
    type, public :: Activation
        character(len=10) :: activation_type
        real(dp) :: alpha = 0.2
    contains
        procedure :: linear => activ_linear
        procedure :: d_linear => activ_d_linear
        procedure :: sigmoid => activ_sigmoid
        procedure :: d_sigmoid => activ_d_sigmoid
        procedure :: tanh => activ_tanh
        procedure :: d_tanh => activ_d_tanh
        procedure :: relu => activ_relu
        procedure :: d_relu => activ_d_relu
        procedure :: prelu => activ_prelu
        procedure :: d_prelu => activ_d_prelu
        procedure :: softmax => activ_softmax
        procedure :: d_softmax => activ_d_softmax
        procedure :: print => activation_print
        procedure :: get_activation => get_activation
        procedure :: get_d_activation => get_d_activation
        procedure :: forward => forward
        procedure :: backpropagation => backpropagation
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
        

    function activ_linear(this, x) result(y)
        class(Activation), intent(in) :: this
        real(dp), dimension(:,:), intent(in) :: x
        real(dp), dimension(size(x,1), size(x,2)) :: y

        y = x
    end function activ_linear

    function activ_d_linear(this, x) result(dy)
        class(Activation), intent(in) :: this
        real(dp), dimension(:,:), intent(in) :: x
        real(dp), dimension(size(x,1), size(x,2)) :: dy

        dy = 1.0
    end function activ_d_linear

    function activ_sigmoid(this, x) result(y)
        class(Activation), intent(in) :: this
        real(dp), dimension(:,:), intent(in) :: x
        real(dp), dimension(size(x,1), size(x,2)) :: y

        y = 1.0 / (1.0 + exp(-x))
    end function activ_sigmoid
    
    function activ_d_sigmoid(this, x) result(dy)
        class(Activation), intent(in) :: this
        real(dp), dimension(:,:), intent(in) :: x
        real(dp), dimension(size(x,1), size(x,2)) :: dy
        
        real(dp), dimension(size(x,1), size(x,2)) :: sig_x
        
        sig_x = this%sigmoid(x)
        
        dy = sig_x * (one_dp - sig_x)
    end function activ_d_sigmoid

    function activ_tanh(this, x) result(y)
        class(Activation), intent(in) :: this
        real(dp), dimension(:,:), intent(in) :: x
        real(dp), dimension(size(x,1), size(x,2)) :: y

        y = tanh(x)
    end function activ_tanh
    
    function activ_d_tanh(this, x) result(dy)
        class(Activation), intent(in) :: this
        real(dp), dimension(:,:), intent(in) :: x
        real(dp), dimension(size(x,1), size(x,2)) :: dy

        dy = 1.0 - tanh(x)**2
    end function activ_d_tanh

    function activ_relu(this, x) result(y)
        class(Activation), intent(in) :: this
        real(dp), dimension(:,:), intent(in) :: x
        real(dp), dimension(size(x,1), size(x,2)) :: y

        y = max(0.0, x)
    end function activ_relu
    
    function activ_d_relu(this, x) result(dy)
        class(Activation), intent(in) :: this
        real(dp), dimension(:,:), intent(in) :: x
        real(dp), dimension(size(x,1), size(x,2)) :: dy

        dy = sign(one_dp, x)
    end function activ_d_relu

    function activ_prelu(this, x) result(y)
        class(Activation), intent(in) :: this
        real(dp), dimension(:,:), intent(in) :: x
        real(dp), dimension(size(x,1), size(x,2)) :: y

        y = max(0.0, x) + this%alpha * min(0.0, x)
    end function activ_prelu
    
    function activ_d_prelu(this, x) result(dy)
        class(Activation), intent(in) :: this
        real(dp), dimension(:,:), intent(in) :: x
        real(dp), dimension(size(x,1), size(x,2)) :: dy
        
        
        dy = sign(one_dp, x) + (one_dp - sign(one_dp, x)) * this%alpha
    end function activ_d_prelu

    
    subroutine sub_softmax(x, result)
        real(dp), dimension(:,:), intent(in) :: x
        real(dp), dimension(size(x,1), size(x,2)) :: result

        real(dp), dimension(size(x,1), size(x,2)) :: z
        real(dp), dimension(size(x,1)) :: exp_sum
        integer :: i
        
        !TODO implement more efficient
        
        do i = 1, size(x, 1)
            z(i, :) = x(i, :) - maxval(x(i, :))
        end do

        exp_sum = sum(exp(z), dim=2)

        do i = 1, size(x, 1)
            result(i, :) = exp(z(i, :)) / exp_sum(i)
        end do
        
    end subroutine sub_softmax
    
    function activ_softmax(this, x) result(y)
        class(Activation), intent(in) :: this
        real(dp), dimension(:,:), intent(in) :: x
        real(dp), dimension(size(x,1), size(x,2)) :: y
        
        call sub_softmax(x,y)
        
    end function activ_softmax
    
    

    
    function activ_d_softmax(this, x) result(dy)
        class(Activation), intent(in) :: this
        real(dp), dimension(:,:), intent(in) :: x
        real(dp), dimension(size(x,1), size(x,2)) :: dy

        call d_softmax(x,dy)
        
    end function activ_d_softmax
    
    
    subroutine d_softmax(x, result)
        real(dp), intent(in) :: x(:,:)
        real(dp), intent(out) :: result(:,:)
        integer :: m, d
        real(dp), dimension(:,:), allocatable :: a
        real(dp), dimension(:,:,:), allocatable :: tensor1, tensor2
        integer :: i, j, k

        ! Get dimensions of x
        m = size(x, 1)
        d = size(x, 2)

        ! Allocate arrays
        allocate(a(m, d), tensor1(m, d, d), tensor2(m, d, d))

        ! Reshape input matrix if necessary
        if (m == 1) then
            a = reshape(x, (/1, d/))
        else
            a = x
        end if
        ! Compute softmax
        call sub_softmax(a, a)

        ! Compute tensor1
        do i = 1, m
            do j = 1, d
                do k = 1, d
                    tensor1(i, j, k) = a(i, j) * a(i, k)
                end do
            end do
        end do

        ! Compute tensor2
        do i = 1, m
            do j = 1, d
                do k = 1, d
                    if (j == k) then
                        tensor2(i, j, k) = a(i, j) * (1.0d0 - a(i, j))
                    else
                        tensor2(i, j, k) = -a(i, j) * a(i, k)
                    end if
                end do
            end do
        end do

        ! Compute result
        do i = 1, m
            do j = 1, d
                do k = 1, d
                    result(i, j) = tensor2(i, j, k) - tensor1(i, j, k)
                end do
            end do
        end do

    end subroutine d_softmax
!active_softmax_sub

    !> not currently used
    subroutine softmax(a, m, d)
        real(dp), intent(inout) :: a(:,:)
        integer, intent(in) :: m, d
        real(dp) :: exp_sum, max_val
        integer :: i, j

        ! Find the maximum value in each row
        do i = 1, m
            max_val = maxval(a(i, 1:d))
            ! Compute softmax for each element in the row
            exp_sum = 0.0d0
            do j = 1, d
                a(i, j) = exp(a(i, j) - max_val)
                exp_sum = exp_sum + a(i, j)
            end do
            ! Normalize the row to obtain probabilities
            do j = 1, d
                a(i, j) = a(i, j) / exp_sum
            end do
        end do

    end subroutine softmax


    subroutine activation_print(this)
        class(Activation), intent(in) :: this
        print *, 'Activation Type: ', this%activation_type
    end subroutine activation_print

    function get_activation(this, x) result(y)
        class(Activation), intent(in) :: this
        real(dp), dimension(:,:), intent(in) :: x
        real(dp), dimension(size(x,1), size(x,2)) :: y

        select case (this%activation_type)
            case ('linear')
                y = this%linear(x)
            case ('sigmoid')
                y = this%sigmoid(x)
            case ('tanh')
                y = this%tanh(x)
            case ('relu')
                y = this%relu(x)
            case ('prelu')
                y = this%prelu(x)
            case ('softmax')
                y = activ_softmax(this, x)
            case default
                print *, 'Unknown activation type: ', this%activation_type
                stop
        end select
    end function get_activation

    function get_d_activation(this, x) result(dy)
        class(Activation), intent(in) :: this
        real(dp), dimension(:,:), intent(in) :: x
        real(dp), dimension(size(x,1), size(x,2)) :: dy

        select case (this%activation_type)
            case ('linear')
                dy = this%d_linear(x)
            case ('sigmoid')
                dy = this%d_sigmoid(x)
            case ('tanh')
                dy = this%d_tanh(x)
            case ('relu')
                dy = this%d_relu(x)
            case ('prelu')
                dy = this%d_prelu(x)
            case ('softmax')
                dy = this%d_softmax(x)
            case default
                print *, 'Unknown activation type: ', this%activation_type
                stop
        end select
    end function get_d_activation

    function forward(this, X) result(z)
        class(Activation), intent(in) :: this
        real(dp), dimension(:,:), intent(in) :: X
        real(dp), dimension(size(X,1), size(X,2)) :: z
        
        cache_x = X
        
        z = this%get_activation(X)
    end function forward

    function backpropagation(this, dz) result(dx)
        class(Activation), intent(in) :: this
        real(dp), dimension(:,:), intent(in) :: dz
        real(dp), dimension(size(dz,1), size(dz,2)) :: dx

        real(dp), dimension(size(dz,1), size(dz,2)) :: f_prime
        
        f_prime = this%get_d_activation(cache_x)

        dx = dz * f_prime
    end function backpropagation

end module Activation_mod





    


!program test_activation_circle
!    use class_Circle
!    use Activation_mod
!    implicit none
!
!    type(Circle) :: c           ! Declare a variable of type Circle.
!    type(Activation) :: act     ! Declare a variable of type Activation.
!
!    c = Circle(1.5)             ! Use the implicit constructor, radius = 1.5.
!    call c%print                ! Call the type-bound subroutine for Circle.
!
!    act = Activation('sigmoid') ! Use the constructor with activation type 'sigmoid'.
!    call act%print              ! Call the type-bound subroutine for Activation.
!    
!    
!
!end program test_activation_circle
!