program QuantumDynamics

    use plplot
    use wave1d

    implicit none

    integer, parameter :: num_dt = 10000
    real(8), parameter :: sample_length = 1
    integer, parameter :: num_nodes = 200
    integer :: i

    call init_model(sample_length, num_nodes)
    call init_graphics()
    call plot_wave()
    do i = 1, num_dt
        call step()
        if (mod(i, 1) == 0) then
            call plot_wave()
        end if
    end do
    call plspause(.false.)
    call plend()

contains

    subroutine init_graphics()

        call plparseopts(PL_PARSE_FULL)

        ! gnuplot color scheme
        call plscol0(0, 255, 255, 255)  ! white
        call plscol0(1, 255, 0, 0)      ! red
        call plscol0(2, 0, 255, 0)      ! green
        call plscol0(3, 0, 0, 255)      ! blue
        call plscol0(4, 255, 0, 255)    ! magenta
        call plscol0(5, 0, 255, 255)    ! cyan
        call plscol0(6, 255, 255, 0)    ! yellow
        call plscol0(7, 0, 0, 0)        ! black
        call plscol0(8, 255, 76, 0)     ! orange
        call plscol0(9, 128, 128, 128)  ! gray

        call plsdev("xcairo")
        call plinit()

        call plcol0(7)
        call plenv(0d0, sample_length, -0.1d0, 0.1d0, 0, 0)
        call pllab("x", "\psi", "Wave plot")

        call plcol0(1)

    end subroutine

end program
