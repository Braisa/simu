program evolution_runs

    use get_kinds
    use fcc_parameters
    use interfaces

    implicit none

    ! Energy and particle variables

    real (kind = double) :: potential, kinetic, energy
    real (kind = double) :: rpot, r2pot, dphi, d2phi, kinetic_inv
    real (kind = double), dimension(N) :: rx, ry, rz
    real (kind = double), dimension(N) :: vx, vy, vz
    real (kind = double), dimension(N) :: ax, ay, az

    ! Evolution variables
    
    integer (kind = int) :: s, r, runs
    integer (kind = int), parameter :: steps = 500000
    integer (kind = int), parameter :: accum_interval = 100, intervals = steps / accum_interval

    ! Averaging variables

    integer (kind = int) :: interval

    real (kind = double), dimension(intervals) :: kinetic_accum, kinetic_inv_accum
    real (kind = double), dimension(intervals) :: dphi_accum, d2phi_accum
    real (kind = double), dimension(intervals) :: dphi_kininv_accum, dphi2_kininv_accum

    real (kind = double) :: kinetic_avg, kinetic_inv_avg
    real (kind = double) :: dphi_avg, d2phi_avg
    real (kind = double) :: dphi_kininv_avg, dphi2_kininv_avg

    real (kind = double) :: kinetic_unc, kinetic_inv_unc
    real (kind = double) :: dphi_unc, d2phi_unc
    real (kind = double) :: dphi_kininv_unc, dphi2_kininv_unc

    ! Magnitude variables

    real (kind = int), parameter :: freedom = 3*N - 3
    real (kind = double), parameter :: freedom_factor = freedom/2.d00 - 1.d00
    real (kind = double) :: temperature, temperature_unc, &
                            volume_heat_capacity, volume_heat_capacity_unc, &
                            pressure, pressure_unc, &
                            energy_expansion, energy_expansion_unc, &
                            gamma, gamma_unc, &
                            adiabatic_compressibility, adiabatic_compressibility_unc, &
                            isothermal_compressibility, isothermal_compressibility_unc, &
                            isobaric_expansion, isobaric_expansion_unc, &
                            adiabatic_expansion, adiabatic_expansion_unc, &
                            pressure_heat_capacity, pressure_heat_capacity_unc

    ! I/O variables

    character (len = 25) :: data_file, rva_file, save_rva_file, magnitude_log_file
    character (len = 18) :: save_file
    character (len = 3) :: stat
    integer (kind = int) :: io, jo
    character (len = 2) :: run_char
    logical :: exists
    9000 format (a25)
    9001 format (3(1pe13.6))
    9002 format (i8, 2(1pe13.6))
    9003 format (4(1pe13.6))
    6000 format (21(a8, 1x))
    6001 format (i8, 1x, 20(SP, 1pe13.6, 1x))

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

    ! INPUT AND INITIALIZATION

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

    ! Ask for run number

    print *, "Number of runs? (<100 to avoid formatting issues)"
    read (*, *) runs

    ! EVOLUTION RUNS

    magnitude_log_file = TRIM(save_file) // "_mg.txt"

    inquire (file=magnitude_log_file, exist=exists)
    if (.NOT. exists) then
        stat = "new"
    else
        stat = "old"
    end if

    open (newunit=io, file=magnitude_log_file, status=stat, action="write")

    write(io, 6000) "RUN", "T", "s_T", "C_V", "s_C_V", "P", "s_P", "alphaE", "s_alphaE", "gamma", "s_gamma", &
                    "kS", "s_kS", "kT", "s_kT", "alphaP", "s_alphaP", "alphaS", "s_alphaS", "C_P", "s_C_P"

    do r = 1, runs

        ! EVOLUTION STEPS FOR EACH RUN

        do s = 1, steps
        
            call lj_verlet_rpot(rx, ry, rz, vx, vy, vz, ax, ay, az, potential, kinetic, rpot, r2pot)

            ! Apply corrections

            potential = potential + energy_correction
            rpot = rpot + rpot_correction
            r2pot = r2pot + r2pot_correction

            ! Compute magnitudes

            energy = potential + kinetic

            if (MOD(s, accum_interval) == 0) then

                interval = s / accum_interval

                kinetic_inv = 1/kinetic
                dphi = Vi/3.d00 * rpot
                d2phi = Vi*Vi/9.d00 * (r2pot - 2.d00 * rpot)

                ! Accumulations
                
                kinetic_accum(interval) = kinetic
                kinetic_inv_accum(interval) = kinetic_inv
                dphi_accum(interval) = dphi
                d2phi_accum(interval) = d2phi
                dphi_kininv_accum(interval) = dphi * kinetic_inv
                dphi2_kininv_accum(interval) = dphi * dphi * kinetic_inv

                write (*, "(a13, i2, a1, i2, a1, a4, i6, a1, i6, a2)", advance="no") &
                      "Progress: Run", r, "/", runs, "|", "Step", s, "/", steps, CHAR(13)
            
            end if

        end do

        ! Compute averages and their uncertainties

        kinetic_avg = SUM(kinetic_accum) / intervals
        kinetic_inv_avg = SUM(kinetic_inv_accum) / intervals
        dphi_avg = SUM(dphi_accum) / intervals
        d2phi_avg = SUM(d2phi_accum) / intervals
        dphi_kininv_avg = SUM(dphi_kininv_accum) / intervals
        dphi2_kininv_avg = SUM(dphi2_kininv_accum) / intervals

        kinetic_unc = DSQRT(SUM((kinetic_avg - kinetic_accum)**2) / intervals / (intervals - 1.d00))
        kinetic_inv_unc = DSQRT(SUM((kinetic_inv_avg - kinetic_inv_accum)**2) / intervals / (intervals - 1.d00))
        dphi_unc = DSQRT(SUM((dphi_avg - dphi_accum)**2) / intervals / (intervals - 1.d00))
        d2phi_unc = DSQRT(SUM((d2phi_avg - d2phi_accum)**2) / intervals / (intervals - 1.d00))
        dphi_kininv_unc = DSQRT(SUM((dphi_kininv_avg - dphi_kininv_accum)**2) / intervals / (intervals - 1.d00))
        dphi2_kininv_unc = DSQRT(SUM((dphi2_kininv_avg - dphi2_kininv_accum)**2) / intervals / (intervals - 1.d00))

        ! Compute direct magnitudes and their uncertainties

        temperature = 2.d00 * kinetic_avg / freedom
        temperature_unc = 2.d00 / freedom * kinetic_unc

        volume_heat_capacity = 1.d00 / (1.d00 + freedom_factor * kinetic_avg * kinetic_inv_avg)
        volume_heat_capacity_unc = volume_heat_capacity**2 * freedom_factor &
                                   * DSQRT((kinetic_inv_avg*kinetic_unc)**2 + (kinetic_avg*kinetic_inv_unc)**2)

        pressure = D * temperature - dphi_avg
        pressure_unc = DSQRT(D**2 * temperature_unc**2 + dphi_unc**2)

        energy_expansion = Vi / (dphi_avg - freedom_factor * kinetic_avg * dphi_kininv_avg)
        energy_expansion_unc = energy_expansion**2 &
                               * DSQRT(dphi_unc**2 + &
                               freedom_factor**2*((dphi_kininv_avg*kinetic_unc)**2 + (kinetic_avg*dphi_kininv_unc)**2))

        gamma = N / volume_heat_capacity + V * freedom_factor * (dphi_avg * kinetic_inv_avg - dphi_kininv_avg)
        gamma_unc = DSQRT((N/volume_heat_capacity**2*volume_heat_capacity_unc)**2 &
                    + (V*freedom_factor)**2*((kinetic_inv_avg*dphi_unc)**2 + (dphi_avg*kinetic_inv_unc)**2 + dphi_kininv_unc**2))

        adiabatic_compressibility = D * temperature * (1.d00 + 2.d00 * gamma - N / volume_heat_capacity) + &
                                    V * d2phi_avg - V * freedom_factor * (dphi2_kininv_avg &
                                    - 2.d00 * dphi_avg * dphi_kininv_avg + dphi_avg * dphi_avg * kinetic_inv_avg)
        adiabatic_compressibility_unc = DSQRT((D*(1.d00 + 2*gamma - N/volume_heat_capacity)*temperature_unc)**2 &
                                        + (2.d00*D*temperature*gamma_unc)**2 + &
                                        (N*D*temperature/volume_heat_capacity**2*volume_heat_capacity_unc)**2 + &
                                        (V*d2phi_unc)**2 + (V*freedom_factor*dphi2_kininv_unc)**2 + &
                                        (2.d00*V*freedom_factor*(dphi_kininv_avg-dphi_avg*kinetic_inv_avg)*dphi_unc)**2 + &
                                        (2.d00*V*freedom_factor*dphi*dphi_kininv_unc)**2 + &
                                        (V*freedom_factor*dphi**2*kinetic_inv_unc)**2)

        adiabatic_compressibility = 1 / adiabatic_compressibility
        adiabatic_compressibility_unc = adiabatic_compressibility**2 * adiabatic_compressibility_unc

        ! Compute indirect magnitudes

        isothermal_compressibility = 1 / adiabatic_compressibility - temperature * volume_heat_capacity * Vi * gamma * gamma
        isothermal_compressibility_unc = DSQRT((adiabatic_compressibility_unc/adiabatic_compressibility**2)**2 + &
                                         (temperature*volume_heat_capacity*Vi*gamma**2)**2 * ((temperature_unc/temperature)**2 + &
                                         (volume_heat_capacity_unc/volume_heat_capacity)**2 + (2*gamma_unc/gamma)**2))

        isothermal_compressibility = 1 / isothermal_compressibility
        isothermal_compressibility_unc = isothermal_compressibility**2 * isothermal_compressibility_unc
        
        isobaric_expansion = volume_heat_capacity * Vi * gamma * isothermal_compressibility
        isobaric_expansion_unc = isobaric_expansion * DSQRT((volume_heat_capacity_unc/volume_heat_capacity)**2 + &
                                 (gamma_unc/gamma)**2 + (isothermal_compressibility_unc/isothermal_compressibility)**2)

        adiabatic_expansion = -1.d00 / gamma / temperature
        adiabatic_expansion_unc = adiabatic_expansion * DSQRT((gamma_unc/gamma)**2 + (temperature_unc/temperature)**2)

        pressure_heat_capacity = volume_heat_capacity * isothermal_compressibility / adiabatic_compressibility
        pressure_heat_capacity_unc = pressure_heat_capacity * DSQRT((volume_heat_capacity_unc/volume_heat_capacity)**2 &
                                     + (isothermal_compressibility_unc/isothermal_compressibility)**2 + &
                                     (adiabatic_compressibility_unc/adiabatic_compressibility)**2)

        write(io, 6001) r, temperature, temperature_unc, &
                        volume_heat_capacity, volume_heat_capacity_unc, &
                        pressure, pressure_unc, &
                        energy_expansion, energy_expansion_unc, &
                        gamma, gamma_unc, &
                        adiabatic_compressibility, adiabatic_compressibility_unc, &
                        isothermal_compressibility, isothermal_compressibility_unc, &
                        isobaric_expansion, isobaric_expansion_unc, &
                        adiabatic_expansion, adiabatic_expansion_unc, &
                        pressure_heat_capacity, pressure_heat_capacity_unc

        ! Store final RVA for the run

        write (run_char, "(i0.2)") r
        save_rva_file = TRIM(save_file) // "_" // TRIM(run_char) // ".bin"

        open (newunit=jo, file=save_rva_file, status="new", action="write", form="unformatted")

            write(jo) rx, ry, rz, vx, vy, vz, ax, ay, az

        close(jo)

    end do

    close(io)

end program evolution_runs