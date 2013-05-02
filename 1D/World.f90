module World

    use plplot
    use TimeStep

    implicit none
    private

    real(8), parameter :: pi = 4 * atan(1d0)
    complex(8), parameter :: ii = (0d0, 1d0)
    complex(8), allocatable :: b(:), psi(:)
    real(8), allocatable :: potential(:), x(:)
    real(8) :: L, eps
    integer :: n

    public create_world
    public plot_wave
    public step

contains

    subroutine create_world(sample_length)

        integer, intent(in) :: sample_length

        n = sample_length
        L = n
        eps = 1d-6

        allocate(b(n), psi(n), x(n), potential(n))

        call linspace(n, x)
        call init_potential()
        call init_wave(0.3d0)
        call init_method(potential)

    end subroutine

    subroutine step()

        call crank_nicolson(psi, b, eps)
        !call split_operator(psi, x)

    end subroutine

    subroutine plot_wave()

        call plclear()
        call plline(x, real(psi))
        call plcol0(6)
        call plline(x, aimag(psi))
        call plcol0(3)
        call plline(x, real(psi * conjg(psi)))
        call plcol0(7)
        call plline(x, 7 * potential)
        call plcol0(1)
        call plflush()

    end subroutine

    subroutine init_potential()

        integer :: i

        potential = 0d0
        do i = 1, n
            !if ((i > 0.4 * n .and. i < (0.4 * n + 2)) .or. & 
                !(i > (0.5 * n  + 2 * pi) .and. i < (0.5 * n + 2 * pi + 2))) then
                !potential(i) = 0.2d0
                potential = (2 * (x - L / 2) / L)**12
            !end if
        end do

    end subroutine

    subroutine init_wave(p)

        real(8), intent(in) :: p
        real(8) :: arg(n)
        real(8) :: x_0, k, d

        k = 300 / L
        d = L / 100
        x_0 = p * L
        arg = x - x_0
        psi = exp(-arg**2 / (2 * d**2)) &
            * exp(ii * k * arg)
        arg = x - (x_0 + L + 1)
        psi = psi + exp(-arg**2 / (2 * d**2)) &
            * exp(ii * k * arg)

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
