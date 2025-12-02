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
    
    integer (kind = int) :: s, r
    integer (kind = int), parameter :: steps = 500000, runs = 10, log_interval = 100
    
    ! Accumulating variables

    real (kind = double) :: energy_accum, potential_accum
    real (kind = double) :: kinetic_accum, kinetic_inv_accum
    real (kind = double) :: dphi_accum, d2phi_accum
    real (kind = double) :: dphi_kininv_accum, dphi2_kininv_accum

    real (kind = double), dimension(runs) :: energy_avg, potential_avg, kinetic_avg
    real (kind = double) :: kinetic_inv_avg
    real (kind = double) :: dphi_avg, d2phi_avg
    real (kind = double) :: dphi_kininv_avg, dphi2_kininv_avg

    ! Magnitude variables

    real (kind = int), parameter :: freedom = 3*N - 3
    real (kind = double), parameter :: freedom_factor_up = freedom/2.d00 - 1.d00
    real (kind = double), parameter :: freedom_factor_down = 2.d00/freedom - 1.d00

    real (kind = double), dimension(runs) :: temperature, pressure, volume_specific_heat, gammaG, &
                                             adiabatic_compressibility, energy_expansion, alt_energy_expansion

    ! Final magnitude variables

    real (kind = double) :: energy_final, potential_final, kinetic_final, &
                            temperature_final, pressure_final, volume_specific_heat_final, gammaG_final, &
                            adiabatic_compressibility_final, energy_expansion_final, alt_energy_expansion_final

    real (kind = double) :: energy_unc, potential_unc, kinetic_unc, &
                            temperature_unc, pressure_unc, volume_specific_heat_unc, gammaG_unc, &
                            adiabatic_compressibility_unc, energy_expansion_unc, alt_energy_expansion_unc

    real (kind = double) :: isothermal_compressibility, isothermal_compressibility_unc, &
                            isobaric_expansion, isobaric_expansion_unc, &
                            adiabatic_expansion, adiabatic_expansion_unc, &
                            pressure_specific_heat, pressure_specific_heat_unc

    ! I/O variables

    character (len = 25) :: data_file, rva_file, save_rva_file
    character (len = 25) :: magnitude_log_file, final_log_file, run_log_file
    character (len = 15) :: save_file
    character (len = 3) :: stat
    integer (kind = int) :: io, jo, ko
    character (len = 2) :: run_char
    logical :: exists

    ! Data file formats
    9000 format (a25)
    9001 format (3(1pe13.6))
    9002 format (i8, 2(1pe13.6))
    9003 format (4(1pe13.6))
    ! Magnitude log file formats
    6000 format (11(a, 1x))
    6001 format (i0.2, 1x, 10(SP, 1pe13.6, 1x))
    ! Final log file formats
    5000 format (2(a, 1x, 1pe13.6, 1x))

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

    magnitude_log_file = TRIM(save_file) // "_mg.txt"
    final_log_file = TRIM(save_file) // "_fin.txt"

    ! EVOLUTION RUNS

    inquire (file=magnitude_log_file, exist=exists)
    if (.NOT. exists) then
        stat = "new"
    else
        stat = "old"
    end if

    open (newunit=io, file=magnitude_log_file, status=stat, action="write")

    write(io, 6000) "RUN", "E", "V", "K", "T", "P", "cV", "gammaG", "kS", "alphaE", "alt_alphaE"

    do r = 1, runs

        write (run_char, "(i0.2)") r
    
        ! Set accumulators to zero

        energy_accum = 0.d00
        potential_accum = 0.d00
        kinetic_accum = 0.d00
        kinetic_inv_accum = 0.d00
        dphi_accum = 0.d00
        d2phi_accum = 0.d00
        dphi_kininv_accum = 0.d00
        dphi2_kininv_accum = 0.d00

        ! EVOLUTION STEPS FOR EACH RUN

        run_log_file = TRIM(save_file) // "_rv" // run_char // ".bin"

        inquire (file=run_log_file, exist=exists)
        if (.NOT. exists) then
            stat = "new"
        else
            stat = "old"
        end if

        open (newunit=ko, file=run_log_file, status=stat, action="write", form="unformatted")

        do s = 1, steps
        
            call lj_verlet_rpot(rx, ry, rz, vx, vy, vz, ax, ay, az, potential, kinetic, rpot, r2pot)

            ! Apply corrections

            potential = potential + energy_correction
            rpot = rpot + rpot_correction
            r2pot = r2pot + r2pot_correction

            ! Compute magnitudes

            energy = potential + kinetic
            kinetic_inv = 1/kinetic
            dphi = Vi/3.d00 * rpot
            d2phi = Vi*Vi/9.d00 * (r2pot - 2.d00 * rpot)

            ! Accumulations

            energy_accum = energy_accum + energy
            potential_accum = potential_accum + potential
            kinetic_accum = kinetic_accum + kinetic
            kinetic_inv_accum = kinetic_inv_accum + kinetic_inv
            dphi_accum = dphi_accum + dphi
            d2phi_accum = d2phi_accum + d2phi
            dphi_kininv_accum = dphi_kininv_accum + dphi * kinetic_inv
            dphi2_kininv_accum = dphi2_kininv_accum + dphi * dphi * kinetic_inv

            ! Check for logging

            if (MOD(s, log_interval) == 0) then
            
                write(ko) energy, potential, kinetic, rx, ry, rz, vx, vy, vz

                write (*, "(a13, i2, a1, i2, a1, a4, i6, a1, i6, a2)", advance="no") &
                      "Progress: Run", r, "/", runs, "|", "Step", s, "/", steps, CHAR(13)

            end if

        end do

        close(ko)

        ! Compute averages
        
        energy_avg(r) = energy_accum / DBLE(steps)
        potential_avg(r) = potential_accum / DBLE(steps)
        kinetic_avg(r) = kinetic_accum / DBLE(steps)
        kinetic_inv_avg = kinetic_inv_accum / DBLE(steps)
        dphi_avg = dphi_accum / DBLE(steps)
        d2phi_avg = d2phi_accum / DBLE(steps)
        dphi_kininv_avg = dphi_kininv_accum / DBLE(steps)
        dphi2_kininv_avg = dphi2_kininv_accum / DBLE(steps)

        ! Compute direct magnitudes

        temperature(r) = 2.d00/freedom * kinetic_avg(r)
        pressure(r) = D * temperature(r) - dphi_avg
        volume_specific_heat(r) = 1.d00 / (1.d00 + freedom_factor_down * kinetic_avg(r) * kinetic_inv_avg) / N
        gammaG(r) = 1.d00 / volume_specific_heat(r) + V * freedom_factor_up * (dphi_avg * kinetic_inv_avg - dphi_kininv_avg)
        adiabatic_compressibility(r) = D * temperature(r) * (1.d00 + 2.d00 * gammaG(r) - 1.d00 / volume_specific_heat(r)) + &
                                       V * d2phi_avg - V * freedom_factor_up * (dphi2_kininv_avg &
                                       - 2.d00 * dphi_avg * dphi_kininv_avg + dphi_avg * dphi_avg * kinetic_inv_avg)
        adiabatic_compressibility(r) = 1.d00/adiabatic_compressibility(r)
        energy_expansion(r) = -Vi / (dphi_avg + freedom_factor_down * kinetic_avg(r) * dphi_kininv_avg)
        alt_energy_expansion(r) = pressure(r) / (D * volume_specific_heat(r)) - gammaG(r) * temperature(r)
        alt_energy_expansion(r) = 1.d00/alt_energy_expansion(r)

        write(io, 6001) r, energy_avg(r), potential_avg(r), kinetic_avg(r), &
                        temperature(r), pressure(r), volume_specific_heat(r), gammaG(r), &
                        adiabatic_compressibility(r), energy_expansion(r), alt_energy_expansion(r)

        ! Store final RVA for the run

        save_rva_file = TRIM(save_file) // "_" // TRIM(run_char) // ".bin"

        open (newunit=jo, file=save_rva_file, status="new", action="write", form="unformatted")

            write(jo) rx, ry, rz, vx, vy, vz, ax, ay, az
        
        close(jo)

    end do

    close(io)

    ! Compute uncertainties for primary magnitudes

    energy_final = SUM(energy_avg)/runs
    potential_final = SUM(potential_avg)/runs
    kinetic_final = SUM(kinetic_avg)/runs
    temperature_final = SUM(temperature)/runs
    pressure_final = SUM(pressure)/runs
    volume_specific_heat_final = SUM(volume_specific_heat)/runs
    gammaG_final = SUM(gammaG)/runs
    adiabatic_compressibility_final = SUM(adiabatic_compressibility)/runs
    energy_expansion_final = SUM(energy_expansion)/runs
    alt_energy_expansion_final = SUM(alt_energy_expansion)/runs

    energy_unc = DSQRT(SUM((energy_final - energy_avg)**2)/runs/(runs-1.d00))
    potential_unc = DSQRT(SUM((potential_final - potential_avg)**2)/runs/(runs-1.d00))
    kinetic_unc = DSQRT(SUM((kinetic_final - kinetic_avg)**2)/runs/(runs-1.d00))
    temperature_unc = DSQRT(SUM((temperature_final - temperature)**2)/runs/(runs-1.d00))
    pressure_unc = DSQRT(SUM((pressure_final - pressure)**2)/runs/(runs-1.d00))
    volume_specific_heat_unc = DSQRT(SUM((volume_specific_heat_final - volume_specific_heat)**2)/runs/(runs-1.d00))
    gammaG_unc = DSQRT(SUM((gammaG_final - gammaG)**2)/runs/(runs-1.d00))
    adiabatic_compressibility_unc = DSQRT(SUM((adiabatic_compressibility_final - adiabatic_compressibility)**2)/runs/(runs-1.d00))
    energy_expansion_unc = DSQRT(SUM((energy_expansion_final - energy_expansion)**2)/runs/(runs-1.d00))
    alt_energy_expansion_unc = DSQRT(SUM((alt_energy_expansion_final - alt_energy_expansion)**2)/runs/(runs-1.d00))

    ! Compute derivative magnitudes and their uncertainties

    isothermal_compressibility = 1/adiabatic_compressibility_final - &
                                 D * temperature_final * volume_specific_heat_final * gammaG_final**2
    isothermal_compressibility_unc = DSQRT((adiabatic_compressibility_unc/adiabatic_compressibility_final**2)**2 + &
                                     (D*temperature_final*volume_specific_heat_final*gammaG_final**2)**2 * &
                                     ((temperature_unc/temperature_final)**2 + &
                                     (volume_specific_heat_unc/volume_specific_heat_final)**2 + (2*gammaG_unc/gammaG_final)**2))

    isothermal_compressibility = 1/isothermal_compressibility
    isothermal_compressibility_unc = isothermal_compressibility_unc/isothermal_compressibility**2

    isobaric_expansion = D * volume_specific_heat_final * gammaG_final * isothermal_compressibility
    isobaric_expansion_unc = isobaric_expansion * DSQRT((volume_specific_heat_unc/volume_specific_heat_final)**2 + &
                             (gammaG_unc/gammaG_final)**2 + (isothermal_compressibility_unc/isothermal_compressibility)**2)

    adiabatic_expansion = -1.d00 / gammaG_final / temperature_final
    adiabatic_expansion_unc = ABS(adiabatic_expansion) * &
                              DSQRT((gammaG_unc/gammaG_final)**2 + (temperature_unc/temperature_final)**2)

    pressure_specific_heat = volume_specific_heat_final * isothermal_compressibility / adiabatic_compressibility_final
    pressure_specific_heat_unc = pressure_specific_heat * DSQRT((volume_specific_heat_unc/volume_specific_heat_final)**2 &
                                 + (isothermal_compressibility_unc/isothermal_compressibility)**2 + &
                                 (adiabatic_compressibility_unc/adiabatic_compressibility_final)**2)
    
    inquire (file=final_log_file, exist=exists)
    if (.NOT. exists) then
        stat = "new"
    else
        stat = "old"
    end if

    open (newunit=io, file=final_log_file, status=stat, action="write")

        write(io, 5000) "E", energy_final, "s_E", energy_unc
        write(io, 5000) "V", potential_final, "s_V", potential_unc
        write(io, 5000) "K", kinetic_final, "s_K", kinetic_unc
        write(io, 5000) "T", temperature_final, "s_T", temperature_unc
        write(io, 5000) "P", pressure_final, "s_P", pressure_unc
        write(io, 5000) "cV", volume_specific_heat_final, "s_cV", volume_specific_heat_unc
        write(io, 5000) "gammaG", gammaG_final, "s_gammaG", gammaG_unc
        write(io, 5000) "kS", adiabatic_compressibility_final, "s_kS", adiabatic_compressibility_unc
        write(io, 5000) "alphaE", energy_expansion_final, "s_alphaE", energy_expansion_unc
        write(io, 5000) "alt_alphaE", alt_energy_expansion_final, "s_alt_alphaE", alt_energy_expansion_unc
        write(io, 5000) "kT", isothermal_compressibility, "s_kT", isothermal_compressibility_unc
        write(io, 5000) "alphaP", isobaric_expansion, "s_alphaP", isobaric_expansion_unc
        write(io, 5000) "alphaS", adiabatic_expansion, "s_alphaS", adiabatic_expansion_unc
        write(io, 5000) "cP", pressure_specific_heat, "s_cP", pressure_specific_heat_unc

    close(io)

end program evolution_runs