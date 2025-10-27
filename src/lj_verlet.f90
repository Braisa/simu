subroutine lj_verlet(rx, ry, rz, vx, vy, vz, ax, ay, az, potential, kinetic)

    use get_kinds
    use fcc_parameters

    implicit none

    ! I/O variables

    real (kind = double), dimension(N), intent(INOUT) :: rx, ry, rz
    real (kind = double), dimension(N), intent(INOUT) :: vx, vy, vz
    real (kind = double), dimension(N), intent(INOUT) :: ax, ay, az
    real (kind = double), intent(INOUT) :: potential
    real (kind = double), intent(OUT) :: kinetic

    rx = rx + vx*dt + ax*dt212
    ry = ry + vy*dt + ay*dt212
    rz = rz + vz*dt + az*dt212

    vx = vx + ax*dt12
    vy = vy + ay*dt12
    vz = vz + az*dt12

    call lj_potential(rx, ry, rz, potential, ax, ay, az)

    vx = vx + ax*dt12
    vy = vy + ay*dt12
    vz = vz + az*dt12

    kinetic = 0.5d00 * (sum(vx*vx) + sum(vy*vy) + sum(vz*vz))

end subroutine lj_verlet