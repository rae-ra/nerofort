module Conv2D_mod
    use Padding2D_mod
    use Activation_mod, only: Activation
    use Weights_mod, only: Weights
    use Optimizer_mod, only: Optimizer, optimizer_init
    use Math_UTIL, only: dp, einsum
    use str_mod
    implicit none

    type :: Conv2D
        type(Padding2D) :: padding
        real(dp), allocatable :: W(:,:) ! weights
        integer :: F
        integer :: input_shape_x(4)
        integer :: kernel_size(2)
        integer :: Kh, Kw
        integer :: s(2)
        integer :: sh, sw
        type(Activation) :: activation
        logical :: use_bias
        character(len=:), allocatable :: weight_init
        character(len=:), allocatable :: kernel_regularizer
        real(dp) :: seed
        integer :: Nc, Nh, Nw
        integer :: Oh, Ow
        integer :: output_shape(4)
        integer :: input_shape(4)
        real(dp), allocatable :: K(:,:,:,:), X(:,:,:,:)
        real(dp), allocatable :: b(:,:,:,:)
        type(Optimizer) :: optimizer
    end type Conv2D

    contains

    subroutine init_Conv2D(this, filters, kernel_size, s, &
        p, activation_type, use_bias, weight_init, &
            kernel_regularizer, seed, input_shape)
        class(Conv2D), intent(inout) :: this
        integer, intent(in) :: filters
        integer, intent(in) :: kernel_size(2)
        integer, intent(in) :: s(2)
        character(len=*), intent(in) :: p
        character(len=*), intent(in), optional :: activation_type
        logical, intent(in), optional :: use_bias
        character(len=*), intent(in), optional :: weight_init
        character(len=*), intent(in), optional :: kernel_regularizer
        real(dp), intent(in), optional :: seed
        integer, intent(in), optional :: input_shape(4)

        !this%padding%init_Padding2D(p)
        call Padding_init(this%padding, p)


        this%F = filters
        if (present(input_shape)) then
            this%input_shape_x = input_shape
        end if
        this%kernel_size = kernel_size
        this%Kh = kernel_size(1)
        this%Kw = kernel_size(2)
        this%s = s
        this%sh = s(1)
        this%sw = s(2)
        if (present(activation_type)) then
            this%activation = Activation(activation_type)
        end if
        if (present(use_bias)) then
            this%use_bias = use_bias
        else
            this%use_bias = .true.
        end if
        if (present(weight_init)) then
            this%weight_init = weight_init
        end if
        if (present(kernel_regularizer)) then
            this%kernel_regularizer = kernel_regularizer
        else
            this%kernel_regularizer = 'L2'
        end if
        if (present(seed)) then
            this%seed = seed
        end if
    end subroutine init_Conv2D

    subroutine get_dimensions(this, input_shape)
        class(Conv2D), intent(inout) :: this
        integer, intent(in) :: input_shape(4)

        integer :: padded_shape(4)
        integer :: m

        !> init padding
        this%input_shape_x = input_shape
        call pad_get_dim(this%padding, this%input_shape_x, &
            this%kernel_size, this%s)

        !> init shapes

        !> input shape
        this%input_shape = this%padding%output_shape

        if (size(input_shape) == 3) then
            this%Nc = padded_shape(1)
            this%Nh = padded_shape(2)
            this%Nw = padded_shape(3)
        elseif (size(input_shape) == 4) then
            m = input_shape(1)
            this%Nc = padded_shape(2)
            this%Nh = padded_shape(3)
            this%Nw = padded_shape(4)
        end if

        !hight & width for  output shape
        this%Oh = (this%Nh - this%Kh) / this%sh + 1
        this%Ow = (this%Nw - this%Kw) / this%sw + 1

        ! output shape
        if (size(input_shape) == 3) then
            this%output_shape = (/this%F, this%Oh, this%Ow, 1/)
        elseif (size(input_shape) == 4) then
            this%output_shape = (/m, this%F, this%Oh, this%Ow/)
        end if
    end subroutine get_dimensions


    !TEST needed
    subroutine initialize_parameters(this, input_shape, optimizer_type)
        class(Conv2D), intent(inout) :: this
        integer, intent(in) :: input_shape(4)
        character(len=*), intent(in) :: optimizer_type
        type(Weights) :: winit

        integer :: shape_b(3)
        integer :: shape_W(4)

        ! init W
        call get_dimensions(this, input_shape)
        shape_b = (/this%F, this%Oh, this%Ow/)
        shape_W = (/this%F, this%Nc, this%Kh, this%Kw/)
        call init(winit, shape_W,&
            this%weight_init)


        this%W = winit%get_initializer()

        ! init b
        allocate(this%b(shape_b(1),shape_b(2),shape_b(3)))
        this%b = 0.0

        !init optimizer
        call init_Optimizer(this%optimizer, optimizer_type, shape_W, shape_b)

    end subroutine initialize_parameters


    !TEST needed
    subroutine dilate2D(this, X, Dr, Xd)
        class(Conv2D), intent(inout) :: this
        real(dp), intent(in) :: X(:,:,:,:)
        integer, intent(in) :: Dr(2)
        real(dp), allocatable, intent(out) :: Xd(:,:,:,:)
        integer :: dh, dw, m, C, H, W, Hd, Wd, ph, pw, i

        dh = Dr(1)
        dw = Dr(2)

        m = size(X,1)
        C = size(X,2)
        H = size(X,3)
        W = size(X,4)

        Hd = H * dh
        Wd = W * dw

        allocate(Xd(m,C,Hd,Wd))
        Xd = 0.0_dp

        do i = 1, m
            do ph = 1, H
                do pw = 1, W
                    Xd(i,:,1+(ph-1)*dh,1+(pw-1)*dw) = X(i,:,ph,pw)
                end do
            end do
        end do

    end subroutine dilate2D


    subroutine prepare_subMatrix(X, Kh, Kw, s, subM)
        integer, intent(in) ::  Kh, Kw, s(2)
        real(dp), intent(in) :: X(:, :, :, :)
        integer :: m, Nc, Nh, Nw, sh, sw
        integer  :: Oh, Ow
        real(dp), intent(out) :: subM(:,:,:,:,:,:) !(m, Nc, Oh, Ow, Kh, Kw)
        integer :: i, j, k, l, n, o, x_shape(4)

        x_shape = shape(X)
        m  = x_shape(1)
        Nc = x_shape(2)
        Nh = x_shape(3)
        Nw = x_shape(4)

        sh = s(1)
        sw = s(2)

        Oh = (Nh - Kh) / sh + 1
        Ow = (Nw - Kw) / sw + 1

        !TODO
            !1 verify if this is correct
            !2 make it without for loops

        do i = 1, m
            do j = 1, Nc
                do k = 1, Oh
                    do l = 1, Ow
                        do n = 1, Kh
                            do o = 1, Kw
                                subM(i, j, k, l, n, o) = &
                                    X(i, j, (k-1)*sh + n, (l-1)*sw + o)
                            end do
                        end do
                    end do
                end do
            end do
        end do


    end subroutine prepare_subMatrix

    subroutine convolve( X, K, s, mode, Xd)
        real(dp), intent(in) :: X(:,:,:,:), K(:,:,:,:)
        integer, intent(in) :: s(2)
        character(len=:), allocatable, intent(in):: mode
        real(dp), allocatable, intent(out) :: Xd(:,:,:,:)
        real(dp), allocatable :: subM(:,:,:,:,:,:)
        integer :: k_shape(4)

        k_shape = shape(K)

        call prepare_subMatrix(X, k_shape(3), k_shape(4), s, subM)

        select case (mode)
            case ("front")
                call einsum('mfkl,mcijkl->fcij', K, subM, Xd)
            case ("back")
                call einsum('mfkl,mcijkl->fcij', K, subM, Xd)
            case ("param")
                call einsum('mfkl,mcijkl->fcij', K, subM, Xd)
        end select


    end subroutine convolve


    subroutine dZ_D_dX(this, dZ, Nh, Nw, dX)
        class(Conv2D), intent(inout) :: this
        real(dp), intent(in) :: dZ(:,:,:,:)
        integer, intent(in) :: Nh, Nw
        real(dp), allocatable, intent(out) :: dX(:,:,:,:)

        integer :: Hd, Wd, ph, pw
        type(Padding2D) :: pad_back
        real(dp), allocatable :: dXp(:,:,:,:), dZ_Dp(:,:,:,:)
        type(str) :: mode

        mode%chars = "back"

        Hd = size(dZ,3)
        Wd = size(dZ,4)

        ph = Nh - Hd + this%Kh - 1
        pw = Nw - Wd + this%Kw - 1


        call Padding_init(pad_back, p = "int" , ph = ph, pw = pw)

        call pad_forward(pad_back, dZ, this%kernel_size, this%s, dZ_Dp)

        ! convolve dZ_Dp with  K rotated by 180

        call convolve(X = dZ_Dp, K = cshift(cshift(this%K, shift = &
            -1, dim = 3), shift = -1, dim = 4), s = (/1,1/),&
                mode = mode%chars, XD = dXp)

        dx = pad_backward(pad_back, dXp)

    end subroutine dZ_D_dX

    subroutine forward(this, X, a)
        class(Conv2D), intent(inout) :: this
        real(dp), intent(in) :: X(:,:,:,:)
        real(dp), allocatable :: Xp(:,:,:,:), Z(:,:,:,:)
        real(dp), intent(out) :: a(:,:,:,:)
        integer :: i
        type(str) :: mode

        mode%chars = "front"

        this%X = X

        call pad_forward(this%padding, X, this%kernel_size, this%s, Xp)

        call convolve(Xp, this%K, this%s, mode=mode%chars, XD = Z)


        Z = Z + this%b



        a = this%activation%forward(Z)

    end subroutine forward

    subroutine backpropagation(this, X, dZ)
        class(Conv2D), intent(inout) :: this
        real(dp), intent(in) :: X(:,:,:,:)
        real(dp), intent(in) :: dZ(:,:,:,:)
        integer :: i

        call dZ_D_dX(this, dZ)

        if (this%use_bias) then
            do i = 1, this%F
                this%db(i,:,:) = sum(dZ(i,:,:,:), dim=3)
            end do
        end if

        call this%activation%apply_derivative(this%X)

    end subroutine backpropagation

    subroutine update(this)
        class(Conv2D), intent(inout) :: this
        real(dp), allocatable :: dK1(:,:,:,:), db1(:,:,:)

        call get_optimization(this%optimizer, this%dK, this%K)
        call this%optimizer%update(this%db, this%b)

        deallocate(this%dK, this%db)
    end subroutine update

    subroutine apply_gradients(this, dZ, lr)
        class(Conv2D), intent(inout) :: this
        real(dp), intent(in) :: dZ(:,:,:,:)
        real(dp), intent(in) :: lr
        integer :: i

        this%dK = dZ

        do i = 1, this%F
            this%dK(i,:,:) = this%dK(i,:,:) * lr
        end do

        call this%update()

    end subroutine apply_gradients

end module Conv2D_mod
