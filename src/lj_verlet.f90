subroutine lj_verlet(rx, ry, rz, vx, vy, vz, fx, fy, fz, potential, kinetic)

    use get_kinds
    use fcc_parameters

    implicit none

    ! I/O variables

    real (kind = double), dimension(N), intent(INOUT) :: rx, ry, rz
    real (kind = double), dimension(N), intent(INOUT) :: vx, vy, vz
    real (kind = double), dimension(N), intent(INOUT) :: fx, fy, fz
    real (kind = double), intent(INOUT) :: potential
    real (kind = double), intent(OUT) :: kinetic

    rx = rx + vx*dt + fx*dt212
    ry = ry + vy*dt + fy*dt212
    rz = rz + vz*dt + fz*dt212

    vx = vx + fx*dt12
    vy = vy + fy*dt12
    vz = vz + fz*dt12

    call lj_potential(rx, ry, rz, potential, fx, fy, fz)

    vx = vx + fx*dt12
    vy = vy + fy*dt12
    vz = vz + fz*dt12

    kinetic = 0.5d00 * (sum(vx*vx) + sum(vy*vy) + sum(vz*vz))

end subroutine lj_verlet