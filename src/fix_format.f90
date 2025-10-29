program fix_format

    use get_kinds
    use fcc_parameters
    
    implicit none

    real (kind = double), dimension(N) :: energy_log, potential_log, kinetic_log
    integer (kind = int) :: i

    character (len = 25) :: txt_file
    logical :: exists
    integer (kind = int) :: io

    6000 format (3(a1))
    6001 format (3(1pe13.6))
    7000 format (3(a1, 1x))
    7001 format (3(SP, 1pe13.6, 1x))

    print *, "Filename to fix format?"
    read (*, *) txt_file

    inquire (file=txt_file, exist=exists)
    if (.NOT. exists) then
        print *, "File does not exist."
        stop
    end if

    open (newunit=io, file=txt_file, status="old", action="read")

        read(io, 6000)
        do i = 1, N
            read(io, 6001) energy_log(i), potential_log(i), kinetic_log(i)
        end do

    close(io)

    open (newunit=io, file=txt_file, status="old", action="write")

        write(io, 7000) "E", "V", "T"
        do i = 1, N
            write(io, 7001) energy_log(i), potential_log(i), kinetic_log(i)
        end do

    close(io)

end program fix_format