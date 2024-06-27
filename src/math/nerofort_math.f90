module nerofort_math
    use math_act_d_softmax
    use math_act_softmax
    use data_types
    use clip_mod
    use math_act_d_tanh
    use math_act_sigmoid
    use math_act_d_sigmoid
    use math_act_prelu
    use math_act_d_prelu
    use math_act_relu
    use math_act_d_relu
    use utility_mod
    use einsum_mod
    public :: softmax, d_softmax
    public :: dp, qp, sp
    public :: clip
    public :: d_tanh
    public :: sigmoid, d_sigmoid
    public :: prelu, d_prelu, relu, d_relu
    public :: unique_str, mean, std_dev, rand_perm
    public :: einsum
    public :: to_integer, to_logical
    
contains

    elemental function to_integer(bool) result(num)
        logical, intent(in) :: bool
        integer :: num
        
        if (bool .eqv. .TRUE.) then
            num = 1
        else
            num = 0
        end if
    
    end function to_integer
    
    elemental function to_logical(num) result(bool)
        integer, intent(in) :: num
        logical :: bool
        
        if (num == 1) then
            bool = .TRUE.
        else
            bool = .FALSE.
        end if
    end function
    

end module nerofort_math