module interfaces

    use get_kinds
    use fcc_parameters

    implicit none

    interface

        subroutine lj_potential(rx, ry, rz, potential, ax, ay, az)

            use get_kinds
            use fcc_parameters

            implicit none

            real (kind = double), dimension(N), intent(IN) :: rx, ry, rz
            real (kind = double), intent(OUT) :: potential
            real (kind = double), dimension(N), intent(OUT) :: ax, ay, az
        
        end subroutine lj_potential

    end interface

    interface
        
        subroutine lj_verlet(rx, ry, rz, vx, vy, vz, ax, ay, az, potential, kinetic)

            use get_kinds
            use fcc_parameters

            implicit none

            real (kind = double), dimension(N), intent(INOUT) :: rx, ry, rz
            real (kind = double), dimension(N), intent(INOUT) :: vx, vy, vz
            real (kind = double), dimension(N), intent(INOUT) :: ax, ay, az
            real (kind = double), intent(INOUT) :: potential
            real (kind = double), intent(OUT) :: kinetic
        
        end subroutine lj_verlet

    end interface

    interface

        subroutine lj_potential_rpot(rx, ry, rz, potential, ax, ay, az, rpot, r2pot)

            use get_kinds
            use fcc_parameters

            implicit none

            real (kind = double), dimension(N), intent(IN) :: rx, ry, rz
            real (kind = double), intent(OUT) :: potential
            real (kind = double), dimension(N), intent(OUT) :: ax, ay, az
            real (kind = double), intent(OUT) :: rpot, r2pot
        
        end subroutine lj_potential_rpot

    end interface

    interface

        subroutine lj_verlet_rpot(rx, ry, rz, vx, vy, vz, ax, ay, az, potential, kinetic, rpot, r2pot)

            use get_kinds
            use fcc_parameters

            implicit none

            real (kind = double), dimension(N), intent(INOUT) :: rx, ry, rz
            real (kind = double), dimension(N), intent(INOUT) :: vx, vy, vz
            real (kind = double), dimension(N), intent(INOUT) :: ax, ay, az
            real (kind = double), intent(INOUT) :: potential, rpot, r2pot
            real (kind = double), intent(OUT) :: kinetic
        
        end subroutine lj_verlet_rpot
    
    end interface

end module interfaces