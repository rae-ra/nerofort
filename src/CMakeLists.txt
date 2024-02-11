# Set the source files
set(SRC
    data/dict_mod.f90
    data/str_mod.f90
    math/clip_mod.f90
    math/data_types.f90
    math/einsum_mod.f90
    math/Math_util.f90
    math/utility_mod.f90
    math/activations/d_prelu_mod.f90
    math/activations/d_relu_mod.f90
    math/activations/d_sigmoid_mod.f90
    math/activations/d_softmax_mod.f90
    math/activations/d_tanh_mod.f90
    math/activations/prelu_mod.f90
    math/activations/relu_mod.f90
    math/activations/sigmoid_mod.f90
    math/activations/softmax_mod.f90
    Activation_mod.f90
    Batch_norm_mod.f90
    Conv2D_mod.f90
    Cost_mod.f90
    Dense_mod.f90
    Learning_rate_decay_mod.f90
    Main.f90
    Optimizer_mod.f90
    Padding2D_mod.f90
    Pooling2D_mod.f90
    Util.f90
    W_reg_mod.f90
    Weights_mod.f90
)
# Add an executable
add_executable(nerofort ${SRC})
