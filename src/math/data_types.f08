module data_types
    implicit none
    private
    
    integer, parameter, public :: sp = selected_real_kind(6, 37)
    integer, parameter, public :: dp = selected_real_kind(15, 307)
    integer, parameter, public :: qp = selected_real_kind(33, 4931)



end module data_types