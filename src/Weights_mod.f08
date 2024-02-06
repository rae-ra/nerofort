module Weights_mod
    use Math_UTIL, only: dp
    implicit none

    private

    ! Declarations
    type, public :: Weights
        character(len=:), allocatable :: initializer_type
        integer :: w_shape(2)
    contains
        procedure, public :: init
        procedure, public :: get_initializer
    end type Weights

contains

    ! Constructor for Weights type
    subroutine init(this, w_shape, initializer_type)
        class(Weights), intent(inout) :: this
        integer, intent(in) :: w_shape(2)
        character(len=*), intent(in), optional :: initializer_type

        this%w_shape = w_shape
        if (present(initializer_type)) then
            this%initializer_type = initializer_type
        else
            this%initializer_type = "he_normal"
        end if

    end subroutine init

    ! Function to select and call the appropriate initializer
    function get_initializer(this) result(w)
        class(Weights), intent(in) :: this
        real(dp), allocatable :: w(:,:)

        call random_seed()

        select case(this%initializer_type)
            case('zeros')
                w = zeros_initializer(this)
            case('ones')
                w = ones_initializer(this)
            case('random_normal')
                w = random_normal_initializer(this)
            case('random_uniform')
                w = random_uniform_initializer(this)
            case('he_normal')
                w = he_initializer(this)
            case('xavier_normal')
                w = xavier_initializer(this)
            case('glorot_normal')
                w = glorot_initializer(this)
            case default
                error stop "valid initializer options are 'zeros',"// &
                     "'ones',"// ""// "'random_normal', 'random_uniform',"// &
                      "'he_normal',"// "'xavier_normal', and 'glorot_normal'"
        end select
    end function get_initializer

    ! Function to initialize with zeros
    function zeros_initializer(this) result(w)
        class(Weights), intent(in) :: this
        real(dp), allocatable :: w(:,:)

        allocate(w(this%w_shape(1), this%w_shape(2)), source=0.0_dp)
    end function zeros_initializer

    ! Function to initialize with ones
    function ones_initializer(this) result(w)
        class(Weights), intent(in) :: this
        real(dp), allocatable :: w(:,:)

        allocate(w(this%w_shape(1), this%w_shape(2)), source=1.0_dp)
    end function ones_initializer

    ! Function to initialize with random normal values
    function random_normal_initializer(this) result(w)
        class(Weights), intent(in) :: this
        real(dp), allocatable :: w(:,:)

        allocate(w(this%w_shape(1), this%w_shape(2)))
        call random_number(w)
    end function random_normal_initializer

    ! Function to initialize with random uniform values
    function random_uniform_initializer(this) result(w)
        class(Weights), intent(in) :: this
        real(dp), allocatable :: w(:,:)

        allocate(w(this%w_shape(1), this%w_shape(2)))
        call random_number(w)
        w = w - 0.5_dp
    end function random_uniform_initializer

    ! Function to initialize with He initializer
    function he_initializer(this) result(w)
        class(Weights), intent(in) :: this
        real(dp), allocatable :: w(:,:)



        w = sqrt(2.0_dp / this%w_shape(1)) * random_normal_initializer(this)

    end function he_initializer

    ! Function to initialize with Xavier initializer
    function xavier_initializer(this) result(w)
        class(Weights), intent(in) :: this
        real(dp), allocatable :: w(:,:)

        w = sqrt(1.0_dp / this%w_shape(1)) * random_normal_initializer(this)

    end function xavier_initializer

    ! Function to initialize with Glorot initializer
    function glorot_initializer(this) result(w)
        class(Weights), intent(in) :: this
        real(dp), allocatable :: w(:,:)

        w = sqrt(2.0_dp / (this%w_shape(1) + this%w_shape(2))) * &
                                    random_normal_initializer(this)

    end function glorot_initializer

end module Weights_mod
