module wave2d

    implicit none
    private

    real(8), allocatable :: A(:, :), b(:), psi(:)

    public init_model
    public calc_wave
    public plot_wave

contains

    subroutine init_model(n)

        integer, intent(in) :: n
        allocate(A(n**2, n**2), b(n**2), psi(n**2))

        print *, "INIT MODEL"

        psi = 1
        b = 1
        A = 1 
    
    end subroutine

    subroutine calc_wave()

        print *, "calc wave"

    end subroutine

    subroutine plot_wave()

        print *, "plot wave"

    end subroutine

end module
