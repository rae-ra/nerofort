module nerofort_util
    use nerofort_math, only: dp, unique_str, mean, std_dev, rand_perm
    use dict_mod
    use str_mod
    implicit none


!    interface label_encoding
!        module procedure label_encoding_impl
!    end interface label_encoding

    interface onehot
        module procedure onehot_impl
    end interface onehot

!    interface minmax
!        module procedure minmax_impl
!    end interface minmax

    interface standardize
        module procedure standardize_impl
    end interface standardize

    interface inv_standardize
        module procedure inv_standardize_impl
    end interface inv_standardize

    interface train_test_split
        module procedure train_test_split_impl
    end interface train_test_split

contains


!    subroutine label_encoding_impl(Y, result, idx_list)
!        real(dp), dimension(:, :), intent(in) :: Y
!        real(dp), dimension(:, :), intent(out) :: result
!        type(dict), dimension(:), allocatable :: idx_list
!        integer :: col, m, d, idx
!        character(:), allocatable :: val
!
!        m = size(Y, 1)
!        d = size(Y, 2)
!        allocate(idx_list(d))
!
!        do col = 1, d
!            idx = 0
!            do val = unique_str(Y(:, col))
!                idx = idx + 1
!                idx_list(col)%set(val, idx)
!            end do
!
!            do idx = 1, m
!                result(idx, col) = idx_list(col)%get(Y(idx, col))
!            end do
!        end do
!    end subroutine label_encoding_impl

    subroutine onehot_impl(X, X_onehot, indexes)
        type(str), dimension(:), allocatable, intent(in) :: X
        type(str), allocatable :: uniques(:)
        logical, dimension(:, :), intent(out) :: X_onehot
        type(dict), intent(out) :: indexes

        integer :: m, d, col, idx
        type(str), allocatable :: val_

        m = size(X, 1)
        uniques = unique_str(X)
        d = size(uniques, 1)

        call init_dict(indexes, d)


        do col = 1, d
            val_ = uniques(col)
            call set(indexes, val_, col)
        end do

        do idx = 1, m
            val_ = X(idx)
            X_onehot(idx, indexes%get(val_)) = .true.
        end do

    end subroutine onehot_impl

!    subroutine minmax_impl(X, Z, min_X, max_X)
!        real(dp), dimension(:, :), intent(in) :: X
!        real(dp), dimension(:, :), intent(out) :: Z
!        real(dp), dimension(:), intent(out), optional :: min_X, max_X
!
!        if (present(min_X)) then
!            min_X = min(X, dim=1)
!        end if
!
!        if (present(max_X)) then
!            max_X = max(X, dim=1)
!        end if
!
!        Z = (X - min_X) / (max_X - min_X)
!    end subroutine minmax_impl

    subroutine standardize_impl(X, mu, std, Z)
        real(dp), intent(in) :: X(:,:)
        real(dp), dimension(size(X,1), size(X,2)), intent(out) :: Z
        real(dp), dimension(:), intent(out), optional :: mu, std
        integer :: i

        if (present(mu)) then
            do i = 1, size(X,1)
                mu(i) = mean(X(i,:))
            end do
        end if

        if (present(std)) then
            do i = 1, size(X,1)
                std(i) = std_dev(X(i,:))
            end do
        end if

        do i = 1, size(X,2)
            Z(:,i) = (X(:,i) - mu) / std
        end do
    end subroutine standardize_impl

    subroutine inv_standardize_impl(Z, X, mu, std)
        real(dp), dimension(:, :), intent(in) :: Z
        real(dp), dimension(size(Z,1), size(Z,2)), intent(out) :: X
        real(dp), dimension(:), intent(in) :: mu, std
        integer :: i

        do i = 1, size(Z,2)
            X(:,i) = Z(:,i) * std + mu
        end do
    end subroutine inv_standardize_impl

    subroutine train_test_split_impl(X, y, X_train, X_test, y_train, y_test, test_ratio, seed)
        real(dp), dimension(:, :), intent(in) :: X
        real(dp), dimension(:), intent(in) :: y
        real(dp), dimension(:, :), intent(out) :: X_train, X_test
        real(dp), dimension(:), intent(out) :: y_train, y_test
        real(dp), intent(in) :: test_ratio
        integer, intent(in), optional :: seed

        integer, dimension(:), allocatable :: indices
        integer :: train_size, m

        if (present(seed)) then
            ! todo
        end if

        m = size(X, 1)
        train_size = nint((1.0_dp - test_ratio) * dble(m))

        allocate(indices(m))
        call rand_perm(m, indices)

        X_train = X(indices(1:train_size), :)
        X_test = X(indices(train_size + 1:), :)
        y_train = y(indices(1:train_size))
        y_test = y(indices(train_size + 1:))
    end subroutine train_test_split_impl

end module nerofort_util
