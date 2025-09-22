module fcc_parameters

    use get_kinds

    implicit none

    integer (kind = int), parameter :: N = 500, k = 4
    real (kind = double), parameter :: L = 10.d00, rc = 5.d00
    real (kind = double), parameter :: l_box = 2.d00, li_box = l_box**(-1)
    real (kind = double), parameter :: E = -680.d00

end module fcc_parameters