program create_fcc

    use get_kinds
    use fcc_parameters
    use interfaces

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
    real (kind = double), dimension(N) :: ax, ay, az
    real (kind = double) :: px, py, pz

    ! RNG variables

    integer (kind = int) :: random_seed = 3
    real (kind = double) :: random

    ! I/O variables

    character(len = 25) :: data_file, rva_file
    character(len = 3) :: stat
    integer (kind = int) :: io
    logical :: exists
    9000 format (a25)
    9001 format (3(1pe13.6))
    9002 format (i8, 2(1pe13.6))
    9003 format (4(1pe13.6))

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

    call lj_potential(rx, ry, rz, potential, ax, ay, az)
    potential = potential + energy_correction

    print *, N, V, rc, correction_factor, energy_correction
    print *, "Evenly spaced fcc potential: ", potential

    ! Slightly shift particles around evenly spaced positions

    do i = 1, N

        rx(i) = rx(i) + (2.d00 * random(random_seed) - 1.d00) * position_shift * half_l_box
        ry(i) = ry(i) + (2.d00 * random(random_seed) - 1.d00) * position_shift * half_l_box
        rz(i) = rz(i) + (2.d00 * random(random_seed) - 1.d00) * position_shift * half_l_box
    
    end do
    
    call lj_potential(rx, ry, rz, potential, ax, ay, az)
    potential = potential + energy_correction

    print *, "Uneven fcc potential: ", potential

    ! Assign random velocities for each particle, and calculate the resulting kinetic energy
    ! Ensure total momentum is zero by adjusting each particle's velocity

    do i = 1, N

        vx(i) = 2.d00 * random(random_seed) - 1.d00
        vy(i) = 2.d00 * random(random_seed) - 1.d00
        vz(i) = 2.d00 * random(random_seed) - 1.d00

    end do

    px = sum(vx)/N
    py = sum(vy)/N
    pz = sum(vz)/N

    vx = vx - px
    vy = vy - py
    vz = vz - pz

    kinetic = 0.5d00 * (sum(vx*vx) + sum(vy*vy) + sum(vz*vz))

    ! Scale velocities so energy is as desired

    vel_scale = dsqrt((E - potential) / kinetic)

    vx = vx * vel_scale
    vy = vy * vel_scale
    vz = vz * vel_scale

    kinetic = 0.5d00 * (sum(vx*vx) + sum(vy*vy) + sum(vz*vz))
    
    print *, "Total momentum: ", sum(vx), sum(vy), sum(vz)
    print *, "Total forces: ", sum(ax), sum(ay), sum(az)
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

        write(io, 9003) L, Li, rc, rc2
        write(io, 9002) N, V, D
        write(io, 9001) potential + kinetic, potential, kinetic
        write(io, 9000) data_file
        write(io, 9000) rva_file

    close(io)

    ! Save RVA

    inquire(file=rva_file, exist=exists)
    if (exists) then
        stat = "old"
    else
        stat = "new"
    end if

    open(newunit=io, file=rva_file, status=stat, action="write", form="unformatted")

        write(io) rx, ry, rz, vx, vy, vz, ax, ay, az

    close(io)

end program create_fcc