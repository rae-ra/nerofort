module Dense_mod
    use Activation_mod, only: Activation
    use Weights_mod, only: Weights
    use Optimizer_mod, only: Optimizer, get_optimization, optimizer_init
    use Math_util, only :dp
    use w_reg_mod
    implicit none

    private

    type :: Dense
        integer :: neurons
        type(Activation) :: activation
        logical :: use_bias
        character(len=:), allocatable :: weight_initializer_type
        type(W_reg) :: w_reg
        integer :: seed
        integer :: input_dim
        real(dp), allocatable :: W(:,:,:,:)
        real(dp), allocatable :: b(:,:,:,:)
        type(Optimizer) :: optimizer
        real(dp), allocatable :: X(:,:,:,:)
        real(dp), allocatable :: z(:,:,:,:)
        real(dp), allocatable :: dW(:,:,:,:)
        real(dp), allocatable :: db(:,:,:,:)
    contains
        procedure :: init
        procedure :: initialize_parameters
        procedure :: forward
        procedure :: backpropagation
        procedure :: update
    end type Dense

contains

    subroutine init(this, neurons, activation_type, use_bias, &
                     weight_initializer_type, wreg, seed, input_dim)
        class(Dense), intent(inout) :: this
        integer, intent(in) :: neurons
        character(len=*), intent(in), optional :: activation_type
        logical, intent(in), optional :: use_bias
        character(len=*), intent(in), optional :: weight_initializer_type
        type(W_reg), intent(in), optional :: wreg
        integer, intent(in), optional :: seed
        integer, intent(in), optional :: input_dim

        this%neurons = neurons
        if (present(activation_type)) then
            this%activation = Activation(activation_type=activation_type)
        else
            this%activation = Activation("linear")
        end if
        if (present(use_bias)) then
            this%use_bias = use_bias
        else
            this%use_bias = .true.
        end if
        if (present(weight_initializer_type)) then
            this%weight_initializer_type = weight_initializer_type
        else
            this%weight_initializer_type = "he_normal"
        end if
        if (present(wreg)) then
            this%w_reg = wreg
        else
            this%w_reg%str = "L2"
            this%w_reg%num = 0.0_dp

        end if
        if (present(seed)) then
            this%seed = seed
        else
            this%seed = 0
        end if
        if (present(input_dim)) then
            this%input_dim = input_dim
        end if
    end subroutine init

    !TEST needed
    subroutine initialize_parameters(this, hl, optimizer_type)
        class(Dense), intent(inout) :: this
        integer, intent(in) :: hl
        character(len=*), intent(in) :: optimizer_type

        integer :: shape_W(2)
        integer :: shape_b(2)
        type(Weights) :: initializer

        shape_W = [hl, this%neurons]
        shape_b = [this%neurons, 1]
        call initializer%w_init(shape_W, this%weight_initializer_type)
        this%W(:,:,1,1) = initializer%get_initializer()
        allocate(this%b(shape_b(1),shape_b(2),1,1))
        this%b = 0.0_dp
        call optimizer_init(this%optimizer, optimizer_type, shape(this%W), &
            shape(this%b))
    end subroutine initialize_parameters

    !TEST needed
    function forward(this, X) result(a)
        class(Dense), intent(inout) :: this
        real(dp), intent(in) :: X(:,:,:,:)
        real(dp), allocatable :: a(:,:,:,:)

        this%X = X
        this%z(:,:,1,1) = &
            matmul(X(:,:,1,1), this%W(:,:,1,1)) + transpose(this%b(:,:,1,1))
        a = this%activation%forward(this%z)
    end function forward

    !TEST needed
    function backpropagation(this, da) result(dX)
        class(Dense), intent(inout) :: this
        real(dp), intent(in) :: da(:,:,:,:)
        real(dp), allocatable :: dX(:,:,:,:)
        real(dp), allocatable :: dz(:,:,:,:)
        real(dp), allocatable :: dr(:,:,:,:)

        dz = this%activation%backpropagation(da)
        dr = dz
        this%db(:,:,1,1) = reshape(sum(dz, dim=1), [size(dr,2),1])
        this%dW(:,:,1,1) = matmul(transpose(this%X(:,:,1,1)), dr(:,:,1,1))
        dX(:,:,1,1) = matmul(dr(:,:,1,1), transpose(this%W(:,:,1,1)))
    end function backpropagation

    !TEST needed
    subroutine update(this, lr, m, k)
        class(Dense), intent(inout) :: this
        real(dp), intent(in) :: lr
        integer, intent(in) :: m
        integer, intent(in) :: k

        real(dp), allocatable :: dW(:,:,:,:)
        real(dp), allocatable :: db(:,:,:,:)

        call get_optimization(this%optimizer, this%dW, this%db, k, dW, db)
        if (this%w_reg%str == 'L2') then
            dW = dW + this%w_reg%num * this%W
        else if (this%w_reg%str == 'L1') then
            dW = dW + this%w_reg%num * sign(1.0_dp, this%W)
        end if
        this%W = this%W - (dW * lr / m)
        if (this%use_bias) then
            this%b = this%b - (db * lr / m)
        end if
    end subroutine update

end module Dense_mod
