module einsum_mod
    use data_types, only: dp
    implicit none
    private
    public :: einsum

    interface einsum
        module procedure einsum_dp
    end interface einsum

    contains

    subroutine einsum_dp(mode, Ker, Ma, res)
        implicit none

        character(len=*), intent(in) :: mode
        real(dp), intent(in) :: Ker(:,:,:,:)
        real(dp), intent(in) :: Ma(:,:,:,:,:,:)
        real(dp), intent(out) :: res(:,:,:,:)

        integer :: i, j, k, l, m, n, p, q, r, s, t
        integer :: dimK(4), dimM(6), dimRes(4)

        character(len=1) :: charK, charM
        character(len=1), dimension(6) :: chars
        integer, dimension(4) :: posK, posM
        integer, dimension(6) :: posRes

        ! Parse the mode string
        read(mode,'(A)',iostat=i) chars
        if (i /= 0) then
            print *, "Error parsing mode string"
            return
        endif

        ! Determine dimensions of K, M, and result
        dimK = shape(Ker)
        dimM = shape(Ma)
        dimRes = shape(res)

        ! Check if dimensions match Einstein summation convention
        if (size(dimK) /= 4 .or. size(dimM) /= 6) then
            print *, "Invalid dimensions for input arrays"
            return
        endif

        ! Perform contraction
        do i = 1, 4
            do j = 1, 6
                if (chars(j) == 'k') then
                    if (charK == ' ') then
                        charK = 'k'
                        posK(1) = i
                        posM(1) = j
                    else
                        posK(2) = i
                        posM(2) = j
                    endif
                elseif (chars(j) == 'm') then
                    if (charM == ' ') then
                        charM = 'm'
                        posK(3) = i
                        posM(3) = j
                    else
                        posK(4) = i
                        posM(4) = j
                    endif
                elseif (chars(j) == 'f') then
                    posRes(1) = i
                elseif (chars(j) == 'c') then
                    posRes(2) = i
                elseif (chars(j) == 'i') then
                    posRes(3) = i
                elseif (chars(j) == 'j') then
                    posRes(4) = i
                endif
            enddo
        enddo

        ! Perform summation
        do t = 1, dimRes(4)
            do s = 1, dimRes(3)
                do r = 1, dimRes(2)
                    do q = 1, dimRes(1)
                        res(q,r,s,t) = 0.0_dp
                        do p = 1, dimM(6)
                            do n = 1, dimM(5)
                                do m = 1, dimM(4)
                                    do l = 1, dimM(3)
                                        do k = 1, dimM(2)
                                            res(q,r,s,t) = res(q,r,s,t) + &
                                                Ker(posK(1),posK(2),posK(3), &
                                                posK(4)) * Ma(posM(1),posM(2), &
                                                 posM(3),posM(4),l,m,n,p)
                                        enddo
                                    enddo
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo

    end subroutine einsum_dp


end module einsum_mod