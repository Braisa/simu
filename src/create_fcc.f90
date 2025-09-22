program create_fcc

    use get_kinds
    use fcc_parameters

    implicit none

    integer (kind = int) :: N_placed, kx, ky, kz
    real (kind = double) :: kx_box, ky_box, kz_box
    real (kind = double) :: half_l_box = l_box * 0.5d00
    integer (kind = int) :: i

    real (kind = double) :: potential, kinetic
    real (kind = double) :: vel_scale
    real (kind = double), dimension(N) :: fx, fy, fz

    real (kind = double), dimension(N) :: rx, ry, rz
    real (kind = double), dimension(N) :: vx, vy, vz

    integer (kind = int) :: random_seed = 1
    real (kind = double) :: random

    ! Place atoms in evenly spaced fcc

    N_placed = 0

    do kx = 0, 4

        kx_box = kx * l_box

        do ky = 0, 4

            ky_box = ky * l_box

            do kz = 0, 4

                kz_box = kz * l_box

                rx(N_placed + 1) = kx_box
                ry(N_placed + 1) = ky_box
                rz(N_placed + 1) = kz_box

                rx(N_placed + 2) = kx_box
                ry(N_placed + 2) = ky_box + half_l_box
                rz(N_placed + 2) = kz_box + half_l_box

                rx(N_placed + 3) = kx_box + half_l_box
                ry(N_placed + 3) = ky_box
                rz(N_placed + 3) = kz_box + half_l_box

                rx(N_placed + 4) = kx_box + half_l_box
                ry(N_placed + 4) = ky_box + half_l_box
                rz(N_placed + 4) = kz_box

                N_placed = N_placed + 4
            
            end do
        
        end do
    
    end do

    print *, 'Placed particles: ', N_placed

    call lj_potential(rx, ry, rz, potential, fx, fy, fz)

    print *, 'Evenly spaced fcc potential: ', potential

    ! Slightly shift particles around evenly spaced positions
    ! Also assign random velocities for each particle, and calculate the resulting kinetic energy

    kinetic = 0.d00

    do i = 1, N

        rx(i) = rx(i) + (2.d00 * random(random_seed) - 1.d00) * 0.1d00 * l_box
        ry(i) = ry(i) + (2.d00 * random(random_seed) - 1.d00) * 0.1d00 * l_box
        rz(i) = rz(i) + (2.d00 * random(random_seed) - 1.d00) * 0.1d00 * l_box

        vx(i) = 2.d00 * random(random_seed) - 1.d00
        vy(i) = 2.d00 * random(random_seed) - 1.d00
        vz(i) = 2.d00 * random(random_seed) - 1.d00

        kinetic = kinetic + 0.5d00 * (vx(i) * vx(i) + vy(i) * vy(i) + vz(i) * vz(i))
    
    end do

    call lj_potential(rx, ry, rz, potential, fx, fy, fz)

    print *, 'Uneven fcc potential: ', potential

    ! Scale velocities so energy is as desired

    vel_scale = dsqrt((E - potential) / kinetic)
    kinetic = 0.d00

    do i = 1, N
        
        vx(i) = vx(i) * vel_scale
        vy(i) = vy(i) * vel_scale
        vz(i) = vz(i) * vel_scale

        kinetic = kinetic + 0.5d00 * (vx(i) * vx(i) + vy(i) * vy(i) + vz(i) * vz(i))

    end do

    print *, 'Kinetic energy: ', kinetic
    print *, 'Total energy: ', potential + kinetic
    print *, 'Desired energy: ', E

end program create_fcc