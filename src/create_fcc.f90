program create_fcc

    use get_kinds
    use fcc_parameters

    implicit none

    ! fcc variables

    integer (kind = int) :: N_placed, kx, ky, kz
    real (kind = double) :: kx_box, ky_box, kz_box
    real (kind = double) :: half_l_box = l_box * 0.5d00
    integer (kind = int) :: i
    real (kind = double) :: position_shift = 0.2d00

    ! Energy and particle variables

    real (kind = double) :: potential, kinetic
    real (kind = double) :: vel_scale
    real (kind = double), dimension(N) :: rx, ry, rz
    real (kind = double), dimension(N) :: vx, vy, vz
    real (kind = double), dimension(N) :: fx, fy, fz

    ! RNG variables

    integer (kind = int) :: random_seed = 1
    real (kind = double) :: random

    ! I/O variables

    character(len = 25) :: data_file, rva_file
    character(len = 3) :: stat
    integer (kind = int) :: io
    logical :: exists

    ! Correction variables

    real (kind = double) :: correction_factor
    real (kind = double) :: trc6
    real (kind = double) :: energy_correction
    real (kind = double) :: rpot_correction
    real (kind = double) :: r2pot_correction

    ! Compute corrections

    correction_factor = PI*dble(N)*dble(N)/(V*rc*rc2)
    trc6 = 3.d00*rc2*rc2*rc2
    energy_correction = 8.d00/3.d00 * correction_factor * (1.d00/trc6 - 1.d00)
    rpot_correction = 16.d00 * correction_factor * (1.d00 - 2.d00/trc6)
    r2pot_correction = 16.d00 * correction_factor * (26.d00/trc6 - 7.d00)

    ! Ask for filenames

    print *, "Data filename?"
    read (*, *) data_file

    print *, "RVA filename?"
    read (*, *) rva_file

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

    print *, "Placed particles: ", N_placed

    call lj_potential(rx, ry, rz, potential, fx, fy, fz)
    potential = potential + energy_correction

    print *, N, V, rc, correction_factor, energy_correction
    print *, "Evenly spaced fcc potential: ", potential

    ! Slightly shift particles around evenly spaced positions

    do i = 1, N

        rx(i) = rx(i) + (2.d00 * random(random_seed) - 1.d00) * position_shift * half_l_box
        ry(i) = ry(i) + (2.d00 * random(random_seed) - 1.d00) * position_shift * half_l_box
        rz(i) = rz(i) + (2.d00 * random(random_seed) - 1.d00) * position_shift * half_l_box
    
    end do

    ! Assign random velocities for each particle, and calculate the resulting kinetic energy
    ! Ensure total momentum is zero by manually assigning the last particle's velocity

    do i = 1, N-1

        vx(i) = 2.d00 * random(random_seed) - 1.d00
        vy(i) = 2.d00 * random(random_seed) - 1.d00
        vz(i) = 2.d00 * random(random_seed) - 1.d00

    end do

    vx(N) = -sum(vx)
    vy(N) = -sum(vy)
    vz(N) = -sum(vz)

    kinetic = 0.5d00 * (sum(vx*vx) + sum(vy*vy) + sum(vz*vz))

    call lj_potential(rx, ry, rz, potential, fx, fy, fz)
    potential = potential + energy_correction

    print *, "Uneven fcc potential: ", potential

    ! Scale velocities so energy is as desired

    vel_scale = dsqrt((E - potential) / kinetic)
    kinetic = 0.d00

    do i = 1, N
        
        vx(i) = vx(i) * vel_scale
        vy(i) = vy(i) * vel_scale
        vz(i) = vz(i) * vel_scale

        kinetic = kinetic + 0.5d00 * (vx(i) * vx(i) + vy(i) * vy(i) + vz(i) * vz(i))

    end do
    
    print *, sum(vx), sum(vy), sum(vz)
    print *, "Kinetic energy: ", kinetic
    print *, "Total energy: ", potential + kinetic
    print *, "Desired energy: ", E

    ! Save data

    inquire(file=data_file, exist=exists)
    if (exists) then
        stat = "old"
    else
        stat = "new"
    end if

    open(newunit=io, file=data_file, status=stat, action="write")

    write(io, *) N, L, Li, rc, rc2
    write(io, *) V, D
    write(io, *) potential + kinetic, potential, kinetic
    write(io, *) data_file
    write(io, *) rva_file

    close(io)

    ! Save RVA

    inquire(file=rva_file, exist=exists)
    if (exists) then
        stat = "old"
    else
        stat = "new"
    end if

    open(newunit=io, file=rva_file, status=stat, action="write", form="unformatted")

    write(io) rx, ry, rz, vx, vy, vz, fx, fy, fz

    close(io)

end program create_fcc