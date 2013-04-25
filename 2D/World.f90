module World

    use plplot
    use CrankNicolson
    ! use SplitOperator

    implicit none
    private

    real(8), parameter :: pi = 4 * atan(1d0)
    complex(8), parameter :: ii = (0d0, 1d0)
    complex(8), allocatable :: b(:), psi(:, :), psivec(:)
    real(8), allocatable :: potential(:, :), potvec(:), x(:), y(:)
    logical, allocatable :: T(:, :)
    real(8) :: L, W, eps
    integer :: n, m

    public create_world
    public plot_wave
    public step

contains

    subroutine create_world(sample_length, sample_width)

        integer, intent(in) :: sample_length, sample_width

        n = sample_length
        m = sample_width
        L = n
        W = m
        eps = 1d-6

        allocate(b(n * m), psivec(n * m), psi(n, m), x(n), y(m), potvec(n * m), potential(n, m), T(n, m))

        T = .true.
        call linspace(n, x)
        call linspace(m, y)
        call init_potential()
        call init_wave(0d0)
        call init_method(potvec)

    end subroutine

    subroutine step()

        call iterate(psivec, x, y, b, eps)

    end subroutine

    subroutine plot_wave()

        call plclear()
        call plflush()

    end subroutine

    subroutine init_potential()

        ! integer :: i

        potential = 0
        !do i = 1, n
        !    potential = (2 * (x - L / 2) / L)**8
        !end do

        potvec = pack(potential, T)

    end subroutine

    subroutine init_wave(p)

        real(8), intent(in) :: p
        real(8), parameter :: k = 0.35 
        real(8) :: arg(n)
        real(8) :: x_0, d
        integer :: j

        d = L / 30
        x_0 = p * L
        arg = x - x_0
        do j = 1, m
            psi(:, j) = exp(-arg(:)**2 / (2 * d**2)) &
                * exp(ii * k * arg(:))
            arg = x - (x_0 + L + 1)
            psi(:, j) = psi(:, j) + exp(-arg(:)**2 / (2 * d**2)) &
                * exp(ii * k * arg(:))
            psivec = pack(psi, T)
        end do

    end subroutine

    subroutine linspace(n, x)

        integer, intent(in) :: n
        real(8), intent(out) :: x(n)
        integer :: i

        do i = 1, n
            x(i) = i - 1
        end do

    end subroutine
    
end module
