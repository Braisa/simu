subroutine lj_potential(rx, ry, rz, potential, fx, fy, fz)

    use get_kinds
    use fcc_parameters

    implicit none

    ! I/O variables

    real (kind = double), dimension(N), intent(IN) :: rx, ry, rz
    real (kind = double), intent(OUT) :: potential
    real (kind = double), dimension(N), intent(OUT) :: fx, fy, fz

    ! Auxiliary variables

    real (kind = double) :: rx_i, ry_i, rz_i
    real (kind = double) :: rx_ij, ry_ij, rz_ij
    real (kind = double) :: r2_ij, r6_ij, r12_ij
    real (kind = double) :: r2i_ij, r6i_ij, r12i_ij
    real (kind = double) :: f_mod
    integer (kind = int) :: i, j


    potential = 0.d00
    fx = 0.d00
    fy = 0.d00
    fz = 0.d00

    do i = 1, N-1

        rx_i = rx(i)
        ry_i = ry(i)
        rz_i = rz(i)

        do j = i+1, N

            rx_ij = rx_i - rx(j)
            ry_ij = ry_i - ry(j)
            rz_ij = rz_i - rz(j)

            rx_ij = rx_ij - l_box * dnint(rx_ij * li_box)
            ry_ij = ry_ij - l_box * dnint(ry_ij * li_box)
            rz_ij = rz_ij - l_box * dnint(rz_ij * li_box)

            r2_ij = rx_ij * rx_ij + ry_ij * ry_ij + rz_ij * rz_ij

            if (r2_ij < rc * rc) then

                r6_ij = r2_ij * r2_ij * r2_ij
                r12_ij = r6_ij * r6_ij
                r2i_ij = r2_ij**(-1)
                r6i_ij = r6_ij**(-1)
                r12i_ij = r12_ij**(-1)

                potential = potential + r12i_ij - r6i_ij
                f_mod = (2.d00 * r12i_ij - r6i_ij) * r2i_ij
                fx(i) = fx(i) + f_mod * rx_ij
                fx(j) = fx(j) - fx(i)
                fy(i) = fy(i) + f_mod * ry_ij
                fy(j) = fy(j) - fy(i)
                fz(i) = fz(i) + f_mod * rz_ij
                fz(j) = fz(j) - fz(i)

            end if
        
        end do
    
    end do

    potential = 4.d00 * potential
    fx = 24.d00 * fx
    fy = 24.d00 * fy
    fz = 24.d00 * fz

end subroutine lj_potential