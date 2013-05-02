module TimeStep

    use plplot

    implicit none

    include '/usr/include/fftw3.f'

    private

    real(8), allocatable :: potential(:), potmat(:, :), ksq(:, :)

    public init_method
    public crank_nicolson
    public split_operator

contains

    subroutine init_method(input, inputmat)

        real(8) :: input(:), inputmat(:, :)
        real(8), parameter :: pi = 4 * atan(1d0)
        integer :: i, j, n, m

        allocate(potential(size(input)), potmat(size(inputmat, 1), &
            size(inputmat, 2)), ksq(size(inputmat, 1), size(inputmat, 2)))
        potential = input
        potmat = inputmat
        n = size(inputmat, 1)
        m = size(inputmat, 2)
        do i = 1, n
            do j = 1, m
                ksq(i, j) = 4 * pi**2 * (real(i / n)**2 + real(j / m)**2)
            end do
        end do

    end subroutine

    subroutine crank_nicolson(psi, b, m, eps)

        complex(8), intent(inout) :: psi(:)
        complex(8), intent(inout) :: b(:)
        real(8), intent(in) :: eps
        integer, intent(in) :: m

        b = mult_vec(psi, -1, m)
        call bicgstab_solve(psi, b, m, eps)

    end subroutine

    subroutine bicgstab_solve(psi, b, m, eps)

        complex(8), intent(in) :: b(:)
        real(8), intent(in) :: eps
        integer, intent(in) :: m
        complex(8), intent(inout) :: psi(:)
        complex(8), dimension(size(psi)) :: r, rhat, p, v, s, t
        complex(8) :: rho, rhoold, alpha, beta, omega
        real(8) :: babs, error

        r = b - mult_vec(psi, 1, m)
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
            v = mult_vec(p, 1, m)
            alpha = rho / dot_product(rhat, v)
            s = r - alpha * v
            t = mult_vec(s, 1, m)
            omega = dot_product(t, s) / dot_product(t, t)
            psi = psi + alpha * p + omega * s
            r = s - omega * t
            error = sqrt(sum(abs(r)**2)) / babs
        end do

    end subroutine
    
    function mult_vec(vec, h, m)
    
        complex(8), intent(inout) :: vec(:)
        integer, intent(in) :: h, m
        complex(8), parameter :: ii = (0d0, 1d0)
        complex(8) :: mult_vec(size(vec))

        mult_vec = vec + &
            ii * h / 2 * (-1d0 / 2 * &
                (cshift(vec, -m) + cshift(vec, -1) - 4d0 * vec + cshift(vec, 1) + cshift(vec, m)) &
                + potential * vec)
    
    end function

    subroutine split_operator(psi)

        complex(8), intent(inout) :: psi(:, :)
        complex(8), parameter :: ii = (0d0, 1d0)
        real(8), parameter :: dt = 1
        complex(8) :: phi(size(psi, 1), size(psi, 2))
        integer :: plan, n, m

        n = size(psi, 1)
        m = size(psi, 2)
        psi = exp(-ii * dt * potmat) * psi
        call dfftw_plan_dft_2d(plan, n, m, psi, phi, -1, FFTW_ESTIMATE)
        call dfftw_execute_dft(plan, psi, phi)
        call dfftw_destroy_plan(plan)
        phi = exp(-ii / 2 * dt * ksq) * phi
        call dfftw_plan_dft_2d(plan, n, m, phi, psi, 1, FFTW_ESTIMATE)
        call dfftw_execute_dft(plan, phi, psi)
        call dfftw_destroy_plan(plan)
        psi = psi / n

    end subroutine

end module
