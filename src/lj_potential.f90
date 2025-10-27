subroutine lj_potential(rx, ry, rz, potential, ax, ay, az)

    use get_kinds
    use fcc_parameters

    implicit none

    ! I/O variables

    real (kind = double), dimension(N), intent(IN) :: rx, ry, rz
    real (kind = double), intent(OUT) :: potential
    real (kind = double), dimension(N), intent(OUT) :: ax, ay, az

    ! Auxiliary variables

    real (kind = double) :: rx_i, ry_i, rz_i
    real (kind = double) :: rx_ij, ry_ij, rz_ij
    real (kind = double) :: r2_ij
    real (kind = double) :: r2i_ij, r6i_ij, r12i_ij
    real (kind = double) :: f_mod
    integer (kind = int) :: i, j

    potential = 0.d00
    ax = 0.d00
    ay = 0.d00
    az = 0.d00

    do i = 1, N-1

        rx_i = rx(i)
        ry_i = ry(i)
        rz_i = rz(i)

        do j = i+1, N

            rx_ij = rx_i - rx(j)
            ry_ij = ry_i - ry(j)
            rz_ij = rz_i - rz(j)

            ! Calculate distance with respect to closest image

            rx_ij = rx_ij - L * dnint(rx_ij * Li)
            ry_ij = ry_ij - L * dnint(ry_ij * Li)
            rz_ij = rz_ij - L * dnint(rz_ij * Li)

            r2_ij = rx_ij * rx_ij + ry_ij * ry_ij + rz_ij * rz_ij

            if (r2_ij < rc2) then
                
                r2i_ij = 1 / r2_ij
                r6i_ij = r2i_ij * r2i_ij * r2i_ij
                r12i_ij = r6i_ij * r6i_ij

                potential = potential + r12i_ij - r6i_ij
                f_mod = (2.d00 * r12i_ij - r6i_ij) * r2i_ij
                ax(i) = ax(i) + f_mod * rx_ij
                ax(j) = ax(j) - f_mod * rx_ij
                ay(i) = ay(i) + f_mod * ry_ij
                ay(j) = ay(j) - f_mod * ry_ij
                az(i) = az(i) + f_mod * rz_ij
                az(j) = az(j) - f_mod * rz_ij

            end if
        
        end do
    
    end do

    potential = 4.d00 * potential
    ax = 24.d00 * ax
    ay = 24.d00 * ay
    az = 24.d00 * az

end subroutine lj_potential