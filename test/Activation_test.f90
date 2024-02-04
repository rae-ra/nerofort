
    program Activation_test
        use Activation_mod
        
        implicit none

        real(dp), dimension(:, :), allocatable :: input_matrix
        real(dp), dimension(:, :), allocatable  :: output_matrix
        real(dp) :: alpha_value
        type(Activation) :: act 

        ! Initialize input matrix
        allocate(input_matrix(3, 3))
        input_matrix = reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, 7.0_dp, 8.0_dp, 9.0_dp], shape(input_matrix))

        
        
        
        ! Test linear activation
        print *, 'Linear Activation:'
        
        print *, 'Input Matrix:'
        print *, input_matrix
        output_matrix = act%linear(input_matrix)
        print *, 'Output Matrix:'
        print *, output_matrix
        print *, ''

        ! Test sigmoid activation
        print *, 'Sigmoid Activation:'
        output_matrix = act%sigmoid( input_matrix)
        print *, 'Input Matrix:'
        print *, input_matrix
        print *, 'Output Matrix:'
        print *, output_matrix
        print *, ''

        ! Test tanh activation
        print *, 'Tanh Activation:'
        output_matrix = act%tanh(input_matrix)
        print *, 'Input Matrix:'
        print *, input_matrix
        print *, 'Output Matrix:'
        print *, output_matrix
        print *, ''

        ! Test relu activation
        print *, 'ReLU Activation:'
        output_matrix = act%relu(input_matrix)
        print *, 'Input Matrix:'
        print *, input_matrix
        print *, 'Output Matrix:'
        print *, output_matrix
        print *, ''

        ! Test prelu activation with alpha = 0.3
        alpha_value = 0.3_dp
        call act%set_alpha(alpha_value)
        print *, 'PReLU Activation (Alpha = 0.3):'
        output_matrix = act%prelu(input_matrix)
        print *, 'Input Matrix:'
        print *, input_matrix
        print *, 'Output Matrix:'
        print *, output_matrix
        print *, ''

        ! Test softmax activation (assuming input_matrix represents logits)
        print *, 'd_Softmax Activation:'
        output_matrix = act%d_softmax(input_matrix)
        print *, 'Input Matrix (Logits):'
        print *, input_matrix
        print *, 'Output Matrix (Probabilities):'
        print *, output_matrix
        print *, ''

        deallocate(input_matrix)

    end program Activation_test