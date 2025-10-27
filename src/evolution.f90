program evolution

    use get_kinds
    use fcc_parameters
    use interfaces

    implicit none

    ! Energy and particle variables

    real (kind = double) :: potential, kinetic, energy
    real (kind = double), dimension(N) :: rx, ry, rz
    real (kind = double), dimension(N) :: vx, vy, vz
    real (kind = double), dimension(N) :: ax, ay, az

    ! Evolution variables

    integer (kind = int) :: s, steps
    integer (kind = int) :: logging_interval = 100

    ! I/O variables

    character (len = 25) :: data_file, rva_file, save_rva_file, save_log_file
    character (len = 21) :: save_file
    character (len = 3) :: stat
    integer (kind = int) :: io
    logical :: exists
    9000 format (a25)
    9001 format (3(1pe13.6))
    9002 format (i8, 2(1pe13.6))
    9003 format (4(1pe13.6))
    7000 format (a1)
    7001 format (1pe13.6)

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

    ! Ask for data file and check existence

    print *, "Data filename?"
    read (*, *) data_file

    inquire (file=data_file, exist=exists)
    if (.NOT. exists) then
        print *, "Data file does not exist"
        stop
    end if
    
    ! Read files

    open (newunit=io, file=data_file, status="old", action="read")

        read(io, 9003)
        read(io, 9002)
        read(io, 9001) energy, potential, kinetic
        read(io, 9000) data_file
        read(io, 9000) rva_file

    close(io)

    inquire (file=rva_file, exist=exists)
    if (.NOT. exists) then
        print *, "RVA file does not exist"
        stop
    end if

    open (newunit=io, file=rva_file, status="old", action="read", form="unformatted")

        read(io) rx, ry, rz, vx, vy, vz, ax, ay, az

    close(io)

    ! Ask for save filename

    print *, "Save filename? (without extension, .bin will be added to RVA and .txt to log)"
    read (*, *) save_file

    save_rva_file = TRIM(save_file) // ".bin"
    save_log_file = TRIM(save_file) // ".txt"

    ! Ask for step number

    print *, "Number of steps?"
    read (*, *) steps

    ! Log during evolution

    inquire (file=save_log_file, exist=exists)
    if (exists) then
        stat = "old"
    else
        stat = "new"
    end if

    open (newunit=io, file=save_log_file, status=stat, action="write")

        write(io, 7000) "E"

        ! Execute evolution steps

        do s = 1, steps
            
            call lj_verlet(rx, ry, rz, vx, vy, vz, ax, ay, az, potential, kinetic)

            potential = potential + energy_correction
            energy = potential + kinetic

            if (MOD(s, logging_interval) == 0) then

                write(io, 7001) energy
            
            end if

        end do

    close (io)

    ! Save data

    inquire (file=save_rva_file, exist=exists)
    if (exists) then
        stat = "old"
    else
        stat = "new"
    end if

    open (newunit=io, file=save_rva_file, status=stat, action="write", form="unformatted")

        write(io) rx, ry, rz, vx, vy, vz, ax, ay, az
    
    close(io)

end program evolution