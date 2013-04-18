module BiCGSTAB

    implicit none
    private

    public bicgstab_solve

contains

    subroutine bicgstab_solve(A, b, x, eps)

        complex(8), intent(in) :: A(:, :), b(:)
        real(8), intent(in) :: eps
        complex(8), intent(inout) :: x(:)
        complex(8), dimension(size(x)) :: r, rhat, p, v, s, t
        complex(8) :: rho, rhoold, alpha, beta, omega
        real(8) :: babs, error

        r = b - matmul(A, x)
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
            rho = dot_product(conjg(rhat), r)
            beta = (rho / rhoold) * (alpha / omega)
            p = r + beta * (p - omega * v)
            v = matmul(A, p)
            alpha = rho / dot_product(conjg(rhat), v)
            s = r - alpha * v
            t = matmul(A, s)
            omega = dot_product(conjg(t), s) / dot_product(conjg(t), t)
            x = x + alpha * p + omega * s
            r = s - omega * t
            error = sqrt(sum(abs(r)**2)) / babs
        end do

    end subroutine

end module
