module Pooling2D_mod
    use Padding2D_mod
    use Conv2D_MOD, only: prepare_submatrix
    use Math_UTIL, only: dp
    implicit none

    type :: Pooling2D
        type(Padding2D) :: padding
        integer :: pool_size(2)
        integer :: s(2)
        integer :: Kh, Kw
        integer :: sh, sw
        character(len=6) :: pool_type
    end type Pooling2D

contains

    subroutine init_Pooling2D(this, pool_size, s, p, pool_type)
        class(Pooling2D), intent(inout) :: this
        integer, intent(in) :: pool_size(2), s(2)
        character(len=*), intent(in) :: p, pool_type

        call Padding_init(this%padding, p=p)
        this%pool_size = pool_size
        this%s = s
        this%Kh = pool_size(1)
        this%Kw = pool_size(2)
        this%sh = s(1)
        this%sw = s(2)
        this%pool_type = pool_type
    end subroutine init_Pooling2D

    subroutine get_dimensions(this, input_shape, output_shape)
        class(Pooling2D), intent(inout) :: this
        integer, intent(in) :: input_shape(4)
        integer, intent(out) :: output_shape(4)

        integer :: m, Nc, Nh, Nw, Oh, Ow

        if (size(input_shape) == 4) then
            m = input_shape(1)
            Nc = input_shape(2)
            Nh = input_shape(3)
            Nw = input_shape(4)
        elseif (size(input_shape) == 3) then
            Nc = input_shape(1)
            Nh = input_shape(2)
            Nw = input_shape(3)
        endif

        Oh = (Nh - this%Kh) / this%sh + 1
        Ow = (Nw - this%Kw) / this%sw + 1

        if (size(input_shape) == 4) then
            output_shape = (/m, Nc, Oh, Ow/)
        elseif (size(input_shape) == 3) then
            output_shape = (/Nc, Oh, Ow, 1/)
        endif
    end subroutine get_dimensions

    subroutine max_pool_prep_subMatrix(X, pool_size, s, subM)
        real(dp), intent(in) :: X(:,:,:,:)
        integer, intent(in) :: pool_size(2), s(2)
        real(dp), allocatable, intent(out) :: subM(:,:,:,:,:,:)

        call prepare_submatrix(X, pool_size(1), pool_size(2), s, subM)


    end subroutine max_pool_prep_subMatrix

    subroutine pooling(X, pool_size, s, pool_type, Z)
        real(dp), intent(in) :: X(:,:,:,:)
        integer, intent(in) :: pool_size(2), s(2)
        character(len=3), intent(in) :: pool_type
        real(dp), allocatable, intent(out) :: Z(:,:,:,:)

        real(dp), allocatable :: subM(:,:,:,:,:,:)
        integer :: Oh, Ow

        call max_pool_prep_subMatrix(X, pool_size, s, subM)
        Oh = size(subM, 3)
        Ow = size(subM, 4)

        if (pool_type == 'max') then
            Z = maxval(maxval(subM, dim=5), dim=5)
        elseif (pool_type == 'mean') then
            Z = sum(sum(subM, dim=5), dim=5) / &
                (size(subM, 5) * size(subM, 6))
        else
            stop "Allowed pool types are only 'max' or 'mean'."
        endif

        deallocate(subM)

    end subroutine pooling


    subroutine forward(this, X, Z)
        class(Pooling2D), intent(inout) :: this
        real(dp), intent(in) :: X(:,:,:,:)
        real(dp), allocatable, intent(out) :: Z(:,:,:,:)

        real(dp), allocatable :: Xp(:,:,:,:)

        call pad_forward(this%padding, X, this%pool_size, this%s, Xp)

        call pooling(Xp, this%pool_size, this%s, this%pool_type, Z)

        deallocate(Xp)

    end subroutine forward

end module Pooling2D_mod
