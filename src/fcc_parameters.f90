module fcc_parameters

    use get_kinds

    implicit none

    integer (kind = int), parameter :: N = 500, k = 4
    real (kind = double), parameter :: L = 10.d00, rc = 5.d00, V = L*L*L, D = N/V
    real (kind = double), parameter :: Vi = 1/V, Li = 1/L, rc2 = rc*rc
    real (kind = double), parameter :: l_box = 2.d00
    real (kind = double), parameter :: E = -680.d00

    real (kind = double), parameter :: dt = 0.0001
    real (kind = double), parameter :: dt2 = dt*dt, dt212 = 0.5d00*dt*dt, dt12 = 0.5d00*dt

    real (kind = double), parameter :: PI = 3.141592653589793d00

end module fcc_parameters