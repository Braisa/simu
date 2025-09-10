module get_kinds
    implicit none
    integer, parameter :: int = SELECTED_INT_KIND(9)
    integer, parameter :: double = SELECTED_REAL_KIND(15, 307)
end module get_kinds