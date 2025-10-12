module interfaces

    use get_kinds
    use fcc_parameters

    implicit none

    interface

        subroutine lj_potential(rx, ry, rz, potential, fx, fy, fz)

            use get_kinds
            use fcc_parameters

            implicit none

            real (kind = double), dimension(N), intent(IN) :: rx, ry, rz
            real (kind = double), intent(OUT) :: potential
            real (kind = double), dimension(N), intent(OUT) :: fx, fy, fz
        
        end subroutine lj_potential

    end interface

end module interfaces