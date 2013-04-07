program QuantumDynamics

    use plplot
    use Wave1d

    implicit none

    integer, parameter :: tmax = 100
    integer, parameter :: n = 100
    integer :: i

    call init_model(n)
    do i = 1, tmax
        call calc_wave()
        call plot_wave()
    end do

contains

end program
