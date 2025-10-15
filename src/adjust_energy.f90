program adjust_energy

    use get_kinds
    use fcc_parameters

    implicit none

    ! Energy and particle variables

    real (kind = double) :: potential, kinetic, energy
    real (kind = double) :: vel_scale
    real (kind = double), dimension(N) :: rx, ry, rz
    real (kind = double), dimension(N) :: vx, vy, vz
    real (kind = double), dimension(N) :: fx, fy, fz
    real (kind = double) :: px, py, pz

    ! I/O variables

    character(len = 25) :: data_file, rva_file
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

    ! Check if file exists

    inquire(file=data_file, exist=exists)
    if (.NOT. exists) then
        print *, "File does not exist."
        stop
    end if

    ! Read files

    open(newunit=io, file=data_file, status="old", action="read")

    read(io, 9003) 
    read(io, 9002) 
    read(io, 9001) energy, potential, kinetic
    read(io, 9000) data_file
    read(io, 9000) rva_file

    close(io)

    open(newunit=io, file=rva_file, status="old", action="read", form="unformatted")

    read(io) rx, ry, rz, vx, vy, vz, fx, fy, fz

    close(io)

    ! Ensure total momentum is still zero

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
    print *, "Kinetic energy: ", kinetic
    print *, "Total energy: ", potential + kinetic
    print *, "Desired energy: ", E
    
end program adjust_energy