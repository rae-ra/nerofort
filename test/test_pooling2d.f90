module test_pooling2d
    use nerofort_padding2d
    use nerofort_pooling2d
    use nerofort_activation
    implicit none
    private
    public :: test_averagepool_backprop

contains

    subroutine test_averagepool_backprop()
        implicit none

        ! Declare variables
        type(Pooling2D) :: pool
        real(dp), allocatable :: dZ(:,:,:,:), X(:,:,:,:), dXp(:,:,:,:)
        real(dp), allocatable :: expected_dXp(:,:,:,:)
        integer :: i, j, k, l, m, Nc, Nh, Nw

        ! Initialize the pooling layer with relevant parameters
        call pool%init_pooling2D(pool_size=[2,2], strides=[2,2], &
            pool_pad='valid', pool_type='mean')

        ! Define the input dimensions and allocate arrays
        m = 2
        Nc = 3
        Nh = 4
        Nw = 4

        allocate(X(m, Nc, Nh, Nw))
        allocate(dZ(m, Nc, Nh/2, Nw/2))
        allocate(dXp(m, Nc, Nh, Nw))
        allocate(expected_dXp(m + 2, Nc + 3, Nh, Nw))

        ! Initialize X with some values (example values)
        do i = 1, m
            do j = 1, Nc
                do k = 1, Nh
                    do l = 1, Nw
                        X(i, j, k, l) = real(i + j + k + l, dp)
                    end do
                end do
            end do
        end do

        ! Initialize dZ with some values (example values)
        do i = 1, m
            do j = 1, Nc
                do k = 1, Nh/2
                    do l = 1, Nw/2
                        dZ(i, j, k, l) = real(i + j + k + l, dp)
                    end do
                end do
            end do
        end do

        ! Expected output (manually calculated or from a known good implementation)
        ! Fill in expected_dXp with expected values for the test case

        ! Call the function to test
        call averagepool_backprop(pool, dZ, X, dXp)

        if (all(shape(dXp) - shape(expected_dxp) == 0)) then

            ! Verification step: Compare dXp to expected_dXp
            if (all(abs(dXp - expected_dXp) < 1.0e-6_dp)) then
                print *, "Test passed: averagepool_backprop"
            else
                print *, "Test failed: averagepool_backprop"
                print *, "dXp:"
                print *, dXp
                print *, "Expected dXp:"
                print *, expected_dXp
            end if
        else
            print *, "Test failed: averagepool_backprop"
            print*, "Array bound mismatch."
            print*, "Shape of dXp:", shape(dXp)
            print*, "Shape of expected_dXp:", shape(expected_dXp)
        end if

    end subroutine test_averagepool_backprop

end module test_pooling2d

!program test_all
!    use test_pooling2d
!    implicit none
!
!    call test_averagepool_backprop()
!
!end program test_all