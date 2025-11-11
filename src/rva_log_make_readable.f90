program rva_log_make_readable

    use get_kinds
    use fcc_parameters

    implicit none

    ! RVA data

    real (kind = double), dimension(N) :: rx, ry, rz
    real (kind = double), dimension(N) :: vx, vy, vz
    real (kind = double), dimension(N) :: ax, ay, az

    ! I/O variables

    character (len = 25) :: rva_file, save_file
    character (len = 3) :: stat
    integer (kind = int) :: io, jo
    logical :: exists
    8000 format (9(a2, 1x))
    8001 format (9(SP, 1pe13.6, 1x))

    ! Auxiliary

    integer (kind = int) :: M, i, j

    ! Ask for rva file and check existance

    print *, "RVA filename?"
    read (*, *) rva_file

    inquire(file=rva_file, exist=exists)
    if (.NOT. exists) then
        print *, "File does not exist."
        stop
    end if

    open (newunit=jo, file=rva_file, status="old", action="read", form="unformatted")
    
        read(jo) M
        print *, M

    ! Ask for save filename

    print *, "Save filename?"
    read (*, *) save_file

    inquire(file=save_file, exist=exists)
    if (exists) then
        stat = "old"
    else
        stat = "new"
    end if
    
    open (newunit=io, file=save_file, status=stat, action="write")

        write(io, 8000) "rx", "ry", "rz", "vx", "vy", "vz", "ax", "ay", "az"

        do j = 1, 5000

            read(jo) rx, ry, rz, vx, vy, vz, ax, ay, az

            do i = 1, N
                write(io, 8001) rx(i), ry(i), rz(i), vx(i), vy(i), vz(i), ax(i), ay(i), az(i)
            end do

        end do

    close(io)
    
    close(jo)

end program rva_log_make_readable