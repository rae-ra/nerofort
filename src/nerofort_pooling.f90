module nerofort_pooling2d
    use nerofort_padding2d
    use nerofort_conv2d, only: prepare_submatrix
    use nerofort_math, only: dp, to_integer
    implicit none

    type :: Pooling2D
        type(padding) :: pool_pad
        integer :: pool_size(2)
        integer :: s(2)
        integer :: Kh, Kw
        integer :: sh, sw
        character(len=6) :: pool_type
        contains
            procedure :: maxpool_backprop
            procedure :: averagepool_backprop
            procedure :: backprop
            procedure :: forward
            procedure :: dZ_dZp
    end type Pooling2D

contains

    subroutine init_Pooling2D(this, pool_size, s, p, pool_type)
        class(Pooling2D), intent(inout) :: this
        integer, intent(in) :: pool_size(2), s(2)
        character(len=*), intent(in) :: p, pool_type

        call Padding_init(this%pool_pad, p=p)
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
        character(len=6), intent(in) :: pool_type
        real(dp), allocatable, intent(out) :: Z(:,:,:,:)

        real(dp), allocatable :: subM(:,:,:,:,:,:)
        integer :: Oh, Ow

        call max_pool_prep_subMatrix(X, pool_size, s, subM)
        Oh = size(subM, 3)
        Ow = size(subM, 4)

        if (pool_type == 'max') then
            Z = maxval(maxval(subM, dim=5), dim=6)
        elseif (pool_type == 'mean') then
            Z = sum(sum(subM, dim=5), dim=6) / &
                (size(subM, 5) * size(subM, 6))
        else
            stop "Allowed pool types are only 'max' or 'mean'."
        endif

        deallocate(subM)
    end subroutine pooling

    function dZ_dZp(this, dZ) result(dZp)
        class(Pooling2D), intent(in) :: this
        real(dp), intent(in) :: dZ(:,:,:,:)
        real(dp), allocatable :: dZp(:,:,:,:)

        integer :: Kh, Kw, sh, sw, jh, jw, L, l1_size, l2_size, i
        integer, dimension(:), allocatable :: mask
        logical, dimension(:), allocatable :: mask_logical
        integer, allocatable :: l1(:), l2(:), r1(:), r2(:) 
        real(dp), allocatable :: ones_4(:,:,:,:)

        Kh = this%pool_size(1)
        Kw = this%pool_size(2)
        sh = this%s(1)
        sw = this%s(2)

        allocate(ones_4(Kh, Kw, 1, 1))
        ones_4 = 1.0_dp

        call kron(dZ, ones_4, dZp)

        jh = Kh - sh
        jw = Kw - sw

        if (jw /= 0) then
            L = size(dZp, 4) - 1
            l1_size = L - sw + 1
            l2_size = L - sw - jw + 1

            allocate(l1(l1_size), l2(l2_size))

            l1 = [(sw + i - 1, i = 1, l1_size)]
            l2 = [(sw + jw + i - 1, i = 1, l2_size)]

            mask_logical = mod([(i, i = 1, l1_size)], jw) /= 0
            mask = to_integer(mask_logical)

            r1 = pack(l1, mask_logical)
            r2 = pack(l2, mask_logical)

            dZp(:, :, :, r1) = dZp(:, :, :, r1) + dZp(:, :, :, r2)

            dZp = pack(dZp(:,:,:), .not. mask, dim=4)

            deallocate(l1, l2, r1, r2)
        endif

        if (jh /= 0) then
            L = size(dZp, 3) - 1
            l1_size = L - sh + 1
            l2_size = L - sh - jh + 1

            allocate(l1(l1_size), l2(l2_size))

            l1 = [(sh + i - 1, i = 1, l1_size)]
            l2 = [(sh + jh + i - 1, i = 1, l2_size)]

            mask = mod([(i, i = 1, l1_size)], jh) /= 0

            r1 = l1(pack(mask, mask))
            r2 = l2(pack(mask, mask))

            dZp(:, :, r1, :) = dZp(:, :, r1, :) + dZp(:, :, r2, :)

            dZp = pack(dZp, .not. mask, dim=3)

            deallocate(l1, l2, r1, r2)
        endif
    end function dZ_dZp

    function averagepool_backprop(this, dZ, X) result(dXp)
        class(Pooling2D), intent(in) :: this
        real(dp), intent(in) :: dZ(:,:,:,:), X(:,:,:,:)
        real(dp), allocatable :: dXp(:,:,:,:), dZp(:,:,:,:), Xp(:,:,:,:)
        type(padding) :: padb
        integer :: m, Nc, Nh, Nw, ph, pw

        call pad_forward(this%pool_pad, X, this%pool_size, this%s, Xp)

        m = size(Xp, 1)
        Nc = size(Xp, 2)
        Nh = size(Xp, 3)
        Nw = size(Xp, 4)

        dZp = this%dZ_dZp(dZ)

        ph = Nh - size(dZp, 3)
        pw = Nw - size(dZp, 4)

        call Padding_init(padb, p='tuple', ph=ph, pw=pw)

        call pad_forward(padb, dZp, this%s, this%pool_size, dXp)

        dXp = dXp / (Nh * Nw)
    end function averagepool_backprop

    function forward(this, X) result(Z)
        class(Pooling2D), intent(in) :: this
        real(dp), intent(in) :: X(:,:,:,:)
        real(dp), allocatable :: Z(:,:,:,:)

        real(dp), allocatable :: Xp(:,:,:,:)

        call pad_forward(this%pool_pad, X, this%pool_size, this%s, Xp)

        call pooling(Xp, this%pool_size, this%s, this%pool_type, Z)

        deallocate(Xp)
    end function forward

    function backprop(this, dZ) result(dX)
        class(Pooling2D), intent(in) :: this
        real(dp), intent(in) :: dZ(:,:,:,:)
        real(dp), allocatable :: dXp(:,:,:,:), dX(:,:,:,:)

        if (this%pool_type == 'max') then
            dXp = this%maxpool_backprop(dZ, this%X)
        elseif (this%pool_type == 'mean') then
            dXp = this%averagepool_backprop(dZ, this%X)
        endif

        call pad_backward(this%pool_pad, dXp, dX)
    end function backprop

    function maxpool_backprop(this, dZ, X) result(dXp)
        class(Pooling2D), intent(in) :: this
        real(dp), intent(in) :: dZ(:,:,:,:), X(:,:,:,:)
        real(dp), allocatable :: dXp(:,:,:,:)
        ! Implementation needed
    end function maxpool_backprop

end module nerofort_pooling2d
