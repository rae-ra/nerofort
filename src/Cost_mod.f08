module Cost_mod
    use Math_Util, only: dp, clip
    implicit none
    
    public :: clip
        

    type :: Cost
        private
        character(len=:), allocatable :: cost_type
    contains
        procedure :: init => Cost_init
        procedure :: mse => Cost_mse
        procedure :: d_mse => Cost_d_mse
        procedure :: cross_entropy => Cost_cross_entropy
        procedure :: d_cross_entropy => Cost_d_cross_entropy
        procedure :: get_cost => Cost_get_cost
        procedure :: get_d_cost => Cost_get_d_cost
    end type Cost

    contains
    

    subroutine Cost_init(this, cost_type)
        class(Cost), intent(inout) :: this
        character(len=*), intent(in), optional :: cost_type

        if (present(cost_type)) then
          this%cost_type = cost_type
        else
          this%cost_type = 'mse'
        end if
        
        if (this%cost_type /= 'mse' .and. this%cost_type /= 'cross-entropy') then
            error stop "Valid cost functions are only 'mse' and 'cross-entropy'"
        end if

    end subroutine Cost_init

    real(dp) function Cost_mse(this, a, y)
        class(Cost), intent(in) :: this
        real(dp), dimension(:, :), intent(in) :: a, y
        
        Cost_mse = 0.5d0 * sum((norm2(a - y, dim=2))**2)

    end function Cost_mse

    function Cost_d_mse(this, a, y) result(part_deriv)
        class(Cost), intent(in) :: this
        real(dp), dimension(:, :), intent(in) :: a, y
        
        real(dp), dimension(size(a,1), size(a,2)) :: part_deriv
        
        part_deriv = a - y

    end function Cost_d_mse
    
    

    real(dp) function Cost_cross_entropy(this, a, y, epsil)
        class(Cost), intent(in) :: this
        real(dp), dimension(:, :), intent(in) :: a, y
        real(dp), intent(in), optional :: epsil
        
        real(dp), dimension(size(a,1), size(a,2)) :: clip_a
        real(dp) :: eps
        
        if (present(epsil)) then
          eps = epsil
        else
          eps = 1.0d-12
        end if
        clip_a = clip(a,eps,1.0d0 - eps)
        Cost_cross_entropy = -sum(y * log(clip_a))

    end function Cost_cross_entropy

    function Cost_d_cross_entropy(this, a, y, epsi) result(part_deriv)
        class(Cost), intent(in) :: this
        real(dp), dimension(:, :), intent(in) :: a, y
        real(dp), intent(in), optional :: epsi
        
        real(dp), dimension(size(a,1), size(a,2)) :: part_deriv
        
        real(dp), dimension(size(a,1), size(a,2)) :: clip_a
        real(dp) :: eps
        if (present(epsi)) then
          eps = epsi
        else
          eps = 1.0d-12
        end if
        clip_a = clip(a,eps,1.0d0 - eps)
        part_deriv = -y / clip_a

    end function Cost_d_cross_entropy

    real(dp) function Cost_get_cost(this, a, y)
        class(Cost), intent(in) :: this
        real(dp), dimension(:, :), intent(in) :: a, y
        
        select case (this%cost_type)
        case ('mse')
          Cost_get_cost = Cost_mse(this, a, y)
        case ('cross-entropy')
          Cost_get_cost = Cost_cross_entropy(this, a, y)
        case default
          error stop "Valid cost functions are only 'mse' and 'cross-entropy'"
        end select

    end function Cost_get_cost

    function Cost_get_d_cost(this, a, y) result(part_deriv)
        class(Cost), intent(in) :: this
        real(dp), dimension(:, :), intent(in) :: a, y
        real(dp), dimension(size(a,1), size(a,2)) :: part_deriv
        
        select case (this%cost_type)
        case ('mse')
          part_deriv = Cost_d_mse(this, a, y)
        case ('cross-entropy')
          part_deriv = Cost_d_cross_entropy(this, a, y)
        case default
          error stop "Valid cost functions are only 'mse' and 'cross-entropy'"
        end select

    end function Cost_get_d_cost

end module Cost_mod
