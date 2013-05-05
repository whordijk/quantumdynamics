module World

    use plplot
    use TimeStep

    implicit none
    private

    real(8), parameter :: pi = 4 * atan(1d0)
    complex(8), parameter :: ii = (0d0, 1d0)
    complex(8), allocatable :: b(:), psi(:, :), psivec(:)
    real(8), allocatable :: potential(:, :), potvec(:), x(:), y(:), intensity(:, :)
    real(8) :: L, W, eps
    integer :: n, m

    public create_world
    public plot_wave
    public step
    public write_intensity

contains

    subroutine create_world(sample_length, sample_width)

        integer, intent(in) :: sample_length, sample_width

        n = sample_length
        m = sample_width
        L = n
        W = m
        eps = 1d-6

        allocate(b(n * m), psivec(n * m), psi(n, m), x(n), y(m), &
            intensity(n, m), potvec(n * m), potential(n, m))

        intensity = 0
        call linspace(n, x)
        call linspace(m, y)
        call init_potential()
        call init_wave(0.2d0)
        call init_method(potvec, potential)

    end subroutine

    subroutine step()

        call crank_nicolson(psivec, b, m, eps)
        ! call split_operator(psi)

        call calc_intensity()

    end subroutine

    subroutine plot_wave()

        psi = reshape(psivec, (/n, m/))
        call plclear()
        call plcol0(7) 
        call plbox('abc',50d0, 1, 'abc', 10d0, 2)
        call plimage(real(psi), &
            0d0, 1d0 * L, 0d0, 1d0 * W, -1d0, 1d0, 0d0, 1d0 * L, 0d0, 1d0 * W) 
        call plflush()

    end subroutine

    subroutine init_potential()

        integer :: i, j

        potential = 0
        do i = 1, n
            do j = 1, m
                if ((i > 0.3 * L .and. i < (0.3 * L + 1)) .and. &
                    (j < 0.48 * W .or. j > 0.52 * W)) then
                    potential(i, j) = 100d0
                end if
            end do
        end do

        potvec = reshape(potential, (/n * m/))

    end subroutine

    subroutine init_wave(p)

        real(8), intent(in) :: p
        real(8) :: k_x
        real(8) :: arg(n)
        real(8) :: x_0, d
        integer :: j

        k_x = 300 / L
        d = L / 50
        x_0 = p * L
        arg = x - x_0
        psi(:, 1) = exp(-arg**2 / (2 * d**2)) &
            * exp(ii * k_x * arg)
            arg = x - (x_0 + L)
        psi(:, 1) = psi(:, 1) + exp(-arg**2 / (2 * d**2)) &
            * exp(ii * k_x * arg)
        do j = 2, m
            psi(:, j) = psi(:, 1)
        end do
        psivec = reshape(psi, (/n * m/))

    end subroutine

    subroutine linspace(n, x)

        integer, intent(in) :: n
        real(8), intent(out) :: x(n)
        integer :: i

        do i = 1, n
            x(i) = i - 1
        end do

    end subroutine

    subroutine calc_intensity()

        intensity = intensity + real(conjg(psi) * psi)

    end subroutine

    subroutine write_intensity()

        integer :: j

        open (unit = 12 , file = 'intensity.dat' , status = 'replace')
        do j = 1, m
            write (12, *) intensity(nint(0.5 * L), j), j
        end do
        close (unit = 12)

    end subroutine

end module
