module TimeStep

    use plplot

    implicit none

    include '/usr/include/fftw3.f'

    private

    real(8), allocatable :: potential(:)

    public init_method
    public crank_nicolson
    public split_operator

contains

    subroutine init_method(input)

        real(8) :: input(:)

        allocate(potential(size(input)))
        potential = input

    end subroutine

    subroutine crank_nicolson(psi, b, eps)

        complex(8), intent(inout) :: psi(:)
        complex(8), intent(inout) :: b(:)
        real(8), intent(in) :: eps

        b = mult_vec(psi, -1)
        call bicgstab_solve(psi, b, eps)

    end subroutine

    subroutine bicgstab_solve(psi, b, eps)

        complex(8), intent(in) :: b(:)
        real(8), intent(in) :: eps
        complex(8), intent(inout) :: psi(:)
        complex(8), dimension(size(psi)) :: r, rhat, p, v, s, t
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
        complex(8), parameter :: ii = (0d0, 1d0)
        complex(8) :: mult_vec(size(vec))

        mult_vec = vec + &
            ii * h / 2 * (-1d0 / 2 * &
                (cshift(vec, -1) - 2d0 * vec + cshift(vec, 1)) &
                + potential * vec)
    
    end function

    subroutine split_operator(psi, x)

        complex(8), intent(inout) :: psi(:)
        complex(8), parameter :: ii = (0d0, 1d0)
        real(8), intent(in) :: x(:)
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
