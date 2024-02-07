module Optimizer_mod
    use Math_UTIL, only: dp
    implicit none
    private
    public :: Optimizer, optimizer_init, gd, sgd, adam, rmsprop, get_optimization

    type :: Optimizer
        character(len=5) :: optimizer_type
        real(dp), allocatable :: momentum1(:,:), momentum2(:,:), epsilon(:,:)
        real(dp), allocatable :: vdW(:,:), vdb(:,:), SdW(:,:), Sdb(:,:)
    end type optimizer

    interface optimizer_init
        module procedure optimizer_init_dp
    end interface optimizer_init

    interface gd
        module procedure gd_dp
    end interface gd

    interface sgd
        module procedure sgd_dp
    end interface sgd

    interface rmsprop
        module procedure rmsprop_dp
    end interface rmsprop

    interface adam
        module procedure adam_dp
    end interface adam

    interface get_optimization
        module procedure get_optimization_dp
    end interface get_optimization

contains


    subroutine optimizer_init_dp(self, optimizer_type, shape_W, shape_b, &
        momentum1, momentum2, epsilon)
        type(optimizer), intent(out) :: self
        character(len=*), intent(in), optional :: optimizer_type
        integer, intent(in) :: shape_W(2), shape_b(2)
        real(dp), intent(in), optional :: momentum1, momentum2, epsilon

        if (present(optimizer_type)) then
            self%optimizer_type = optimizer_type
        else
            self%optimizer_type = 'adam'
        end if

        if (present(momentum1)) then
            allocate(self%momentum1(shape_W(1),shape_W(2)))
            self%momentum1 = momentum1
        else
            allocate(self%momentum1(shape_W(1),shape_W(2)))
            self%momentum1 = 0.9
        end if

        if (present(momentum2)) then
            allocate(self%momentum2(shape_W(1),shape_W(2)))
            self%momentum2 = momentum2
        else
            allocate(self%momentum2(shape_W(1),shape_W(2)))
            self%momentum2 = 0.999
        end if

        if (present(epsilon)) then
            allocate(self%epsilon(shape_W(1),shape_W(2)))
            self%epsilon = epsilon
        else
            allocate(self%epsilon(shape_W(1),shape_W(2)))
            self%epsilon = 1.0E-8
        end if

        allocate(self%vdW(shape_W(1),shape_W(2)))
        allocate(self%vdb(shape_b(1),shape_b(2)))
        allocate(self%SdW(shape_W(1),shape_W(2)))
        allocate(self%Sdb(shape_b(1),shape_b(2)))
    end subroutine optimizer_init_dp


    !TEST needed
    subroutine gd_dp(self, dW, db, k, result_dW, result_db)
        type(optimizer), intent(inout) :: self
        real(dp), dimension(:,:), intent(in) :: dW, db
        integer, intent(in) :: k
        real(dp), dimension(:,:), intent(out) :: result_dW, result_db

        result_dW = dW
        result_db = db
    end subroutine gd_dp


    !TEST needed
    subroutine sgd_dp(self, dW, db, k, result_dW, result_db)
        type(optimizer), intent(inout) :: self
        real(dp), dimension(:,:), intent(in) :: dW, db
        integer, intent(in) :: k
        real(dp), dimension(:,:), intent(out) :: result_dW, result_db

        ! Update vdW and vdb
        self%vdW = self%momentum1 * self%vdW + (1.0 - self%momentum1) * dW
        self%vdb = self%momentum1 * self%vdb + (1.0 - self%momentum1) * db

        ! Return updated vdW and vdb
        result_dW = self%vdW
        result_db = self%vdb
    end subroutine sgd_dp

    !TEST needed
    subroutine rmsprop_dp(self, dW, db, k, result_dW, result_db)
        type(optimizer), intent(inout) :: self
        real(dp), dimension(:,:), intent(in) :: dW, db
        integer, intent(in) :: k
        real(dp), dimension(size(dW,1),size(dW,2)) :: den_W
        real(dp), dimension(size(db,1),size(db,2)) :: den_b
        real(dp), dimension(:,:), intent(out) :: result_dW, result_db

        ! Update SdW and Sdb
        self%SdW = self%momentum2 * self%SdW + (1.0 - self%momentum2) * (dW**2)
        self%Sdb = self%momentum2 * self%Sdb + (1.0 - self%momentum2) * (db**2)

        ! Compute den_W and den_b
        den_W = sqrt(self%SdW) + self%epsilon
        den_b = sqrt(self%Sdb) + self%epsilon

        ! Return updated dW and db
        result_dW = dW / den_W
        result_db = db / den_b

    end subroutine rmsprop_dp


    !TEST needed
    subroutine adam_dp(self, dW, db, k, result_dW, result_db)
        type(optimizer), intent(inout)        :: self
        real(dp), dimension(:,:), intent(in)     :: dW, db
        integer, intent(in)                   :: k
        real(dp), dimension(size(self%momentum1,1),size(self%momentum1,2)) :: &
            vdW_h, vdb_h, SdW_h, Sdb_h, den_W, den_b
        real(dp), dimension(:,:), intent(out) :: result_dW, result_db

        ! Update vdW and vdb
        self%vdW = self%momentum1 * self%vdW + (1.0 - self%momentum1) * dW
        self%vdb = self%momentum1 * self%vdb + (1.0 - self%momentum1) * db

        ! Update SdW and Sdb
        self%SdW = self%momentum2 * self%SdW + (1.0 - self%momentum2) * (dW**2)
        self%Sdb = self%momentum2 * self%Sdb + (1.0 - self%momentum2) * (db**2)

        ! Correction
        if (k > 1) then
            vdW_h = self%vdW / (1.0 - (self%momentum1**k))
            vdb_h = self%vdb / (1.0 - (self%momentum1**k))
            SdW_h = self%SdW / (1.0 - (self%momentum2**k))
            Sdb_h = self%Sdb / (1.0 - (self%momentum2**k))
        else
            vdW_h = self%vdW
            vdb_h = self%vdb
            SdW_h = self%SdW
            Sdb_h = self%Sdb
        end if

        ! Compute den_W and den_b
        den_W = sqrt(SdW_h) + self%epsilon
        den_b = sqrt(Sdb_h) + self%epsilon

        ! Return updated vdW_h/den_W and vdb_h/den_b
        result_dw = vdW_h / den_W
        result_db = vdb_h / den_b

        return
    end subroutine adam_dp


    !TEST needed
    subroutine get_optimization_dp(self, dW, db, k, result_dW, result_db)
        type(optimizer), intent(inout) :: self
        real(dp), dimension(:,:), intent(in) :: dW, db
        integer, intent(in) :: k
        real(dp), dimension(:,:), intent(out) :: result_dW, result_db

        select case (self%optimizer_type)

            case ('gd')
                call gd(self, dW, db, k, result_dw, result_db)
            case ('sgd')
                call sgd(self, dW, db, k, result_dw, result_db)
            case ('rmsprop')
                call rmsprop(self, dW, db, k, result_dw, result_db)
            case ('adam')
                call adam(self, dW, db, k, result_dw, result_db)
            case default
                stop 'Invalid optimizer type'
        end select

        return
    end subroutine get_optimization_dp

end module Optimizer_mod
