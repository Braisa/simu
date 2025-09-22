module fcc_parameters

    use get_kinds

    implicit none

    integer (kind = int), parameter :: N = 500, k = 4
    real (kind = double), parameter :: L = 10, rc = 5
    real (kind = double), parameter :: l_box = 2.5d00
    real (kind = double), parameter :: E = -680.d00

end module fcc_parameters