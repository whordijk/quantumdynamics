module cranknicolson

    use plplot

    implicit none
    private

    real(8), parameter :: pi = 4 * atan(1d0)
    complex(8), parameter :: ii = (0d0, 1d0)
    complex(8), allocatable :: A(:, :), b(:), psi(:)
    real(8), allocatable :: potential(:), x(:)
    real(8) :: L, dx, dt, eps
    integer :: n

    public init_model
    public step
    public plot_wave

contains

    subroutine init_model(sample_length)

        integer, intent(in) :: sample_length

        n = sample_length
        L = n
        dx = 1
        dt = dx**2
        eps = 1d-6

        allocate(A(n, n), b(n), psi(n), x(n), potential(n))

        call linspace(n, x)
        call init_potential()
        call init_wave(0.5d0)

    end subroutine

    subroutine step()

        call bicgstab_solve(b, psi, eps)
        b = mult_vec(psi, -1)

    end subroutine

    subroutine plot_wave()

        call plclear()
        call plline(x, real(psi))
        call plcol0(6)
        call plline(x, aimag(psi))
        call plcol0(3)
        call plline(x, real(psi * conjg(psi)))
        call plcol0(7)
        call plline(x, potential)
        call plcol0(1)
        call plflush()

    end subroutine

    subroutine init_potential()

        integer :: i

        potential = 0
        do i = 1, n
            !if (i > 0.5 * L) then
            !    potential(i) = 0.05
            !end if
            potential = (2 * (x - L / 2) / L)**8
        end do

    end subroutine

    subroutine init_wave(p)

        real(8), intent(in) :: p
        real(8), parameter :: k = 0.3 
        real(8) :: arg(n)
        real(8) :: x_0, d

        d = L / 30
        x_0 = p * L
        arg = x - x_0
        psi = exp(-arg**2 / (2 * d**2)) &
            * exp(ii * k * arg)
        arg = x - (x_0 + L + dx)
        psi = psi + exp(-arg**2 / (2 * d**2)) &
            * exp(ii * k * arg)
        b = mult_vec(psi, -1)

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
        complex(8), dimension(size(x)) :: r, rhat, p, v, s, t
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
        complex(8) :: mult_vec(n)

        mult_vec = vec + &
            ii * h * dt / 2 * (-1d0 / 2 * &
                (1 / dx**2 * (cshift(vec, -1) - 2d0 * vec + cshift(vec, 1))) &
                + potential * vec)
    
    end function

end module
