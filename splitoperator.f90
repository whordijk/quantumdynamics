module splitoperator

    use plplot

    implicit none

    include '/usr/include/fftw3.f'

    private

    real(8), parameter :: pi = 4 * atan(1d0)
    complex(8) :: ii = (0d0, 1d0)
    complex(8), allocatable :: psi(:)
    real(8), allocatable :: V(:), x(:)
    real(8) :: L, dx, dt
    integer :: n

    public init_model
    public step
    public plot_wave

contains

    subroutine init_model(sample_length)

        real(8), intent(in) :: sample_length

        L = sample_length
        n = nint(sample_length * 300)
        dx = sample_length / (n - 1)
        dt = dx**2

        allocate(psi(n), x(n), V(n))

        call linspace(n, x)
        call init_potential()
        call init_wave(0.2d0)

    end subroutine

    subroutine step()

        complex(8) :: phi(n)
        integer :: plan
                               
        psi = exp(-ii * dt * V) * psi
        call dfftw_plan_dft_1d(plan, n, psi, phi, 1, FFTW_ESTIMATE)
        call dfftw_execute_dft(plan, psi, phi)
        call dfftw_destroy_plan(plan)
        phi = exp(ii * dt * (2 * pi * x / (L * dx))**2 / 2) * phi
        call dfftw_plan_dft_1d(plan, n, phi, psi, -1, FFTW_ESTIMATE)
        call dfftw_execute_dft(plan, phi, psi)
        call dfftw_destroy_plan(plan)
        psi = psi / n

    end subroutine

    subroutine plot_wave()

        call plclear()
        call plline(x, real(psi))
        call plcol0(6)
        call plline(x, aimag(psi))
        call plcol0(3)
        call plline(x, 2 * pi**2 * real(psi * conjg(psi)))
        call plcol0(7)
        call plline(x, 100000 * V)
        call plcol0(1)
        call plflush()

    end subroutine

    subroutine init_potential()

        !V = (2 / L * x - 1)**60
    
    end subroutine

    subroutine init_wave(p)

        real(8), intent(in) :: p
        real(8), parameter :: k = 50
        real(8) :: arg(n)
        real(8) :: x_0, d

        d = 50
        x_0 = p * L
        arg = x - x_0
        psi = exp(-d * arg**2 / 2) &
            * cmplx(cos(k * arg), sin(k * arg)) &
            / sqrt(2 * pi * d)
        arg = x - (x_0 + L + dx)
        psi = psi + exp(-d * arg**2 / 2) &
            * cmplx(cos(k * arg), sin(k * arg)) &
            / sqrt(2 * pi * d)

    end subroutine

    subroutine linspace(n, x)

        integer, intent(in) :: n
        real(8), intent(out) :: x(n)
        integer :: i

        do i = 1, n
            x(i) = (i - 1) * dx
        end do

    end subroutine

end module
