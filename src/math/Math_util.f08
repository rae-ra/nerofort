module Math_util
    use d_softmax_mod
    use softmax_mod
    use data_types
    use clip_mod
    use d_tanh_mod
    use sigmoid_mod
    use d_sigmoid_mod
    use prelu_mod
    use d_prelu_mod
    use relu_mod
    use d_relu_mod
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

end module Math_util