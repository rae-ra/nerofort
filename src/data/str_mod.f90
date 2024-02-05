module str_mod
    implicit none

    Type str
        character(len=:), allocatable :: chars
        contains
        procedure :: trim => trim_str
    End Type

    ! Define operator overloads for the str type
    interface operator(+)
        module procedure str_add
    end interface

    interface operator(==)
        module procedure str_equal
    end interface

    interface operator(<)
        module procedure str_less_than
    end interface

    interface operator(>)
        module procedure str_greater_than
    end interface

    interface operator(*)
        module procedure str_equal_trimmed
    end interface

contains

    ! Overload the + operator for str type
    elemental function str_add(a, b) result(result_str)
        type(str), intent(in) :: a, b
        type(str) :: result_str

        ! Perform the desired concatenation
        result_str%chars = trim(a%chars) // trim(b%chars)

    end function str_add

    ! Overload the == operator for str type
    elemental function str_equal(a, b) result(is_equal)
        type(str), intent(in) :: a, b
        logical :: is_equal

        ! Perform the desired comparison
        is_equal = a%chars == b%chars

    end function str_equal

    ! Overload the < operator for str type
    elemental function str_less_than(a, b) result(is_less)
        type(str), intent(in) :: a, b
        logical :: is_less

        ! Perform the desired comparison
        is_less = a%chars < b%chars

    end function str_less_than

    ! Overload the > operator for str type
    elemental function str_greater_than(a, b) result(is_greater)
        type(str), intent(in) :: a, b
        logical :: is_greater

        ! Perform the desired comparison
        is_greater = a%chars > b%chars

    end function str_greater_than

    elemental function trim_str (this) result(res)
        class(str), intent(in) :: this
        type(str) :: res

        res%chars = trim(this%chars)

    end function trim_str

    ! Overload the === operator for str type
    elemental function str_equal_trimmed(a, b) result(is_equal)
        type(str), intent(in) :: a, b
        logical :: is_equal

        ! Perform the desired comparison
        is_equal = trim(a%chars) == trim(b%chars)

    end function str_equal_trimmed

end module str_mod
