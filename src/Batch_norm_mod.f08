module Batch_norm_mod
    use Math_UTIL, only: dp
    implicit none

    private


    !> POSSIBLE BUGS WHEN BROADCASTING
    type :: Batch_norm
        real(dp) :: momentum
        real(dp) :: epsilon
        real(dp), allocatable :: gamma(:)
        real(dp), allocatable :: beta(:)
        real(dp), allocatable :: running_mean(:)
        real(dp), allocatable :: running_var(:)
        integer :: m
        integer :: d
        real(dp), allocatable :: mu(:)
        real(dp), allocatable :: var(:)
        real(dp), allocatable :: zmu(:,:)
        real(dp), allocatable :: ivar(:)
        real(dp), allocatable :: zhat(:,:)
        real(dp), allocatable :: dgamma(:)
        real(dp), allocatable :: dbeta(:)
    contains
        procedure :: init
        procedure :: initialize_parameters
        procedure :: forward
        procedure :: backpropagation
        procedure :: update
    end type Batch_norm

contains

    subroutine init(this, momentum, epsilon)
        class(Batch_norm), intent(inout) :: this
        real(dp), intent(in), optional :: momentum
        real(dp), intent(in), optional :: epsilon

        if (present(momentum)) then
            this%momentum = momentum
        else
            this%momentum = 0.9_dp
        end if
        if (present(epsilon)) then
            this%epsilon = epsilon
        else
            this%epsilon = 1.0E-6_dp
        end if
    end subroutine init

    subroutine initialize_parameters(this, d)
        class(Batch_norm), intent(inout) :: this
        integer, intent(in) :: d

        allocate(this%gamma(d))
        this%gamma = 1.0_dp
        allocate(this%beta(d))
        this%beta = 0.0_dp
        allocate(this%running_mean(d))
        this%running_mean = 0.0_dp
        allocate(this%running_var(d))
        this%running_var = 0.0_dp
    end subroutine initialize_parameters

    function spr_d(this, x)
    implicit none

        type(Batch_norm)  ::this
        real(dp), allocatable ::spr_d(:,:)
        real(dp), allocatable ::x(:)

        spr_d = spread(x, 2, this%d)

    end function spr_d

    function spr_m(this, x)
    implicit none

        type(Batch_norm)  ::this
        real(dp), allocatable ::spr_m(:,:)
        real(dp), allocatable ::x(:)

        spr_m = spread(x, 2, this%m)

    end function spr_m

    function forward(this, z, mode) result(q)
        class(Batch_norm), intent(inout) :: this
        real(dp), intent(in) :: z(:,:)
        character(len=*), intent(in), optional :: mode
        real(dp), allocatable :: q(:,:)

        if (present(mode) .and. mode == 'test') then
            q = (z - spr_m(this,this%running_mean)) / &
                sqrt(spr_m(this, this%running_var) + this%epsilon)
            q = spr_m(this,this%gamma) * q + spr_m(this,this%beta)

        else
            this%m = size(z, 1)
            this%d = size(z, 2)
            this%mu = sum(z, dim=1) / this%m
            this%zmu = z - spr_d(this,this%mu)
            this%var = sum(this%zmu**2, dim=1) / this%m

            this%ivar = 1.0_dp / sqrt(this%var + this%epsilon)
            this%zhat = this%zmu * spr_d(this,this%ivar)

            q = spr_m(this,this%gamma) * this%zhat + spr_m(this,this%beta)

            this%running_mean = this%momentum * this%running_mean + (1.0_dp &
                - this%momentum) * this%mu
            this%running_var = this%momentum * this%running_var + (1.0_dp - &
                this%momentum) * this%var
        end if
    end function forward

! def backpropagation(self, dq):
!        self.dgamma = np.sum(dq * self.zhat, axis=0)
!        self.dbeta = np.sum(dq, axis=0)
!        dzhat = dq * self.gamma
!        dvar = np.sum(dzhat * self.zmu * (-.5) * (self.ivar**3), axis=0)
!        dmu = np.sum(dzhat * (-self.ivar), axis=0)
!        dz = dzhat * self.ivar + dvar * (2/self.m) * self.zmu + (1/self.m)*dmu
!        return dz

    function backpropagation(this, dq) result(dz)
        class(Batch_norm), intent(inout) :: this
        real(dp), intent(in) :: dq(:,:)
        real(dp), allocatable :: dz(:,:)

        real(dp), allocatable :: tmp(:,:)
        real(dp), allocatable :: dzhat(:,:)
        real(dp), allocatable :: dvar(:,:)
        real(dp), allocatable :: dmu(:,:)
        real(dp), allocatable :: ivar3(:,:)
        real(dp), allocatable :: gamma_spr(:,:)
        real(dp), allocatable :: temp(:)
        real(dp) :: one_m

        ! Compute dgamma and dbeta
        this%dgamma = sum(dq * this%zhat, dim=1)
        this%dbeta = sum(dq, dim=1)

        gamma_spr = spr_m(this, this%gamma)

        dzhat = dq * gamma_spr

        temp = this%ivar**3
        ivar3 = spr_d(this, temp)
        temp = - sum(dzhat * this%zmu * ivar3  &
            / 2.0_dp, dim=1)

        dvar = spr_d(this, temp)

        temp = sum(dzhat * spr_d(this,this%ivar), dim=1)
        dmu = spr_d(this, temp)
        ! Compute dz
        tmp = dq * gamma_spr * ivar3

        one_m = 1.0_dp/this%m
        dz = dzhat * spr_d(this,this%ivar) + &
            dvar * (2.0_dp * one_m) * &
            this%zmu + dmu * one_m

    end function backpropagation

    subroutine update(this, lr, m, k)
        class(Batch_norm), intent(inout) :: this
        real(dp), intent(in) :: lr
        integer, intent(in) :: m
        integer, intent(in) :: k

        this%gamma = this%gamma - this%dgamma * (lr / real(m))
        this%beta = this%beta - this%dbeta * (lr / real(m))
    end subroutine update

end module Batch_norm_mod
