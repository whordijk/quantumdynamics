module SplitOperator

    use plplot

    implicit none

    include '/usr/include/fftw3.f'

    private

    real(8), allocatable :: potential(:)

    public init_method
    public iterate

contains

    subroutine init_method(input)

        real(8), intent(in) :: input(:)
        
        allocate(potential(size(input)))
        potential = input

    end subroutine

    subroutine iterate(psi, x, b, eps)

        complex(8), intent(inout) :: psi(:)
        complex(8), intent(in) :: b(:)
        real(8), intent(in) :: x(:), eps
        complex(8), parameter :: ii = (0d0, 1d0)
        real(8), parameter :: pi = 4 * atan(1d0)
        complex(8) :: phi(size(psi))
        integer :: plan, n

        n = size(psi)                               
        psi = exp(-ii * potential) * psi
        call dfftw_plan_dft_1d(plan, n, psi, phi, 1, FFTW_ESTIMATE)
        call dfftw_execute_dft(plan, psi, phi)
        call dfftw_destroy_plan(plan)
        phi = exp(ii / 2 * (2 * pi * x / n)**2) * phi
        call dfftw_plan_dft_1d(plan, n, phi, psi, -1, FFTW_ESTIMATE)
        call dfftw_execute_dft(plan, phi, psi)
        call dfftw_destroy_plan(plan)
        psi = psi / n

    end subroutine

end module
