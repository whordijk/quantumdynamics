module cranknicolson

    use plplot

    implicit none
    private

    real(8), parameter :: pi = 4 * atan(1d0)
    complex(8), parameter :: ii = (0d0, 1d0)
    complex(8), allocatable :: b(:), psi(:, :), psivec(:) 
    real(8), allocatable :: potential(:, :), potvec(:), x(:), y(:)
    logical, allocatable :: Tmat(:, :)
    real(8) :: L, W, dx, dy, dt, eps
    integer :: n, m

    public init_model
    public step
    public plot_wave

contains

    subroutine init_model(sample_length, sample_width)

        integer, intent(in) :: sample_length, sample_width

        n = sample_length
        m = sample_width
        L = n
        W = m
        dx = 1
        dy = 1
        dt = dx * dy
        eps = 1d-6

        allocate(b(n * m), psi(n, m), psivec(n * m), x(n), y(m), potential(n, m), potvec(n * m), Tmat(n, m))
        Tmat = .true.

        call linspace(n, x)
        call linspace(m, y)
        call init_potential()
        call init_wave(0.2d0)

    end subroutine

    subroutine step()

        call bicgstab_solve(b, psivec, eps)
        b = mult_vec(psivec, -1)

    end subroutine

    subroutine plot_wave()

        psi = unpack(psivec, Tmat, psi)
        call plclear()
        call plmesh(x, y, psi, n, m, opt=DRAW_LINEX)
        !call plimage(real(psi), 0d0, L, 0d0, W, -1d0, 1d0, -1d-3, 1d-3, -1d-3, 1d-3);
        call plflush()

    end subroutine

    subroutine init_potential()

        !integer :: i

        potential = 0
        !do i = 1, n
        !    if (i > 0.4 * L .and. i < 0.5 * L) then
        !        potential(i) = 0.1
        !    end if
        !end do
        potvec = pack(potential, Tmat)

    end subroutine

    subroutine init_wave(p)

        real(8), intent(in) :: p
        real(8), parameter :: k = 0.5 
        real(8) :: arg(n)
        real(8) :: x_0, d
        integer :: j

        d = L / 30
        x_0 = p * L
        arg = x - x_0
        do j = 1, m
            psi(:, j) = exp(-arg(:)**2 / (2 * d**2)) &
                * exp(ii * k * arg(:))
            arg = x - (x_0 + L + dx)
            psi(:, j) = psi(:, j) + exp(-arg(:)**2 / (2 * d**2)) &
                * exp(ii * k * arg(:))
            psivec = pack(psi, Tmat)
            b = mult_vec(psivec, -1)
        end do

    end subroutine

    subroutine linspace(n, x)

        integer, intent(in) :: n
        real(8), intent(out) :: x(n)
        integer :: i

        do i = 1, n
            x(i) = (i - 1) * dx
        end do

    end subroutine
    
        subroutine bicgstab_solve(b, psi, eps)

        complex(8), intent(in) :: b(:)
        real(8), intent(in) :: eps
        complex(8), intent(inout) :: psi(:)
        complex(8), dimension(n * m) :: r, rhat, p, v, s, t
        complex(8) :: rho, rhoold, alpha, beta, omega
        real(8) :: babs, error

        r = b - mult_vec(psi, 1)
        rhat = r
        rho = 1
        alpha = 1
        omega = 1
        v = 0
        p = 0
        babs = sqrt(sum(abs(b)**2))
        error = sqrt(sum(abs(r)**2)) / babs

        do while (error > eps)
            rhoold = rho
            rho = dot_product(rhat, r)
            beta = (rho / rhoold) * (alpha / omega)
            p = r + beta * (p - omega * v)
            v = mult_vec(p, 1)
            alpha = rho / dot_product(rhat, v)
            s = r - alpha * v
            t = mult_vec(s, 1)
            omega = dot_product(t, s) / dot_product(t, t)
            psi = psi + alpha * p + omega * s
            r = s - omega * t
            error = sqrt(sum(abs(r)**2)) / babs
        end do

    end subroutine
    
    function mult_vec(vec, h)
    
        complex(8), intent(inout) :: vec(:)
        integer, intent(in) :: h
        complex(8) :: mult_vec(n * m)

        mult_vec = vec + &
            ii * h * dt / 2 * (-1d0 / 2 * &
                (1 / (dx * dy) * (cshift(vec, -n) + cshift(vec, -1) - 4d0 * vec + cshift(vec, 1) + cshift(vec, n))) &
                + potvec * vec)
    
    end function

end module
