module softmax_mod
    use data_types
    implicit none
    private
    public ::  softmax


    interface softmax
        module procedure softmax_dp_2, softmax_dp_3, &
            softmax_dp_4, softmax_dp_5, softmax_dp_6
    end interface softmax

    contains

    pure function softmax_dp_1(x) result(res)
        real(dp), dimension(:), intent(in)        :: x
        real(dp), dimension(size(x,1))  :: res
        real(dp), dimension(size(x,1)) :: temp

        temp = exp(x - spread(maxval(x, dim=1), 1, size(x,1)))
        res = temp / spread(sum(temp, dim = 1), 1, size(x,1))

    end function softmax_dp_1

    pure function softmax_dp_2(x) result(res)
        real(dp), dimension(:,:), intent(in)        :: x
        real(dp), dimension(size(x,1), size(x,2))  :: res
        real(dp), dimension(size(x,1), size(x,2)) :: temp

        temp = exp(x - spread(maxval(x, dim=2), 2, size(x,2)))
        res = temp / spread(sum(temp, dim = 2), 2, size(x,2))

    end function softmax_dp_2

    pure function softmax_dp_3(x) result(res)
        real(dp), dimension(:,:,:), intent(in)  :: x
        real(dp), dimension(size(x,1), size(x,2), size(x,3))  :: res
        real(dp), dimension(size(x,1), size(x,2), size(x,3)) :: temp

        temp = exp(x - spread(maxval(x, dim=3), 3, size(x,3)))
        res = temp / spread(sum(temp, dim = 3), 3, size(x,3))

    end function softmax_dp_3

    pure function softmax_dp_4(x) result(res)
        real(dp), dimension(:,:,:,:), intent(in)  :: x
        real(dp), dimension(size(x,1), size(x,2), size(x,3), size(x,4)) &
            :: res, temp

        temp = exp(x - spread(maxval(x, dim=4), 4, size(x,4)))
        res = temp / spread(sum(temp, dim = 4), 4, size(x,4))

    end function softmax_dp_4

    pure function softmax_dp_5(x) result(res)
        real(dp), dimension(:,:,:,:,:), intent(in)  :: x
        real(dp), dimension(size(x,1), size(x,2), size(x,3), &
            size(x,4), size(x,5)) :: res, temp

        temp = exp(x - spread(maxval(x, dim=5), 5, size(x,5)))
        res = temp / spread(sum(temp, dim = 5), 5, size(x,5))

    end function softmax_dp_5

    pure function softmax_dp_6(x) result(res)
        real(dp), dimension(:,:,:,:,:,:), intent(in)  :: x
        real(dp), dimension(size(x,1), size(x,2), size(x,3), &
            size(x,4), size(x,5), size(x,6)) :: res, temp

        temp = exp(x - spread(maxval(x, dim=6), 6, size(x,6)))
        res = temp / spread(sum(temp, dim = 6), 6, size(x,6))

    end function softmax_dp_6



end module softmax_mod