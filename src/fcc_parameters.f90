module fcc_parameters

    use get_kinds

    implicit none

    integer (kind = int), parameter :: N = 500, k = 4
    real (kind = double), parameter :: L = 10.d00, rc = 5.d00, V = L*L*L, D = N/V
    real (kind = double), parameter :: Li = 1/L, rc2 = rc*rc
    real (kind = double), parameter :: l_box = 2.d00
    real (kind = double), parameter :: E = -680.d00

    real (kind = double), parameter :: PI = 3.141592653589793d00
    real (kind = double), parameter :: correction_factor = PI*dble(N)*dble(N)/(V*rc*rc2)
    real (kind = double), parameter :: trc6 = 3.d00*rc2*rc2*rc2
    real (kind = double), parameter :: energy_correction = 8.d00/3.d00 * correction_factor * (1.d00/trc6 - 1.d00)
    real (kind = double), parameter :: rpot_correction = 16.d00 * correction_factor * (1.d00 - 2.d00/trc6)
    real (kind = double), parameter :: r2pot_correction = 16.d00 * correction_factor * (26.d00/trc6 - 7.d00)

end module fcc_parameters