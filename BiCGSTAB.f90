module BiCGSTAB

    implicit none
    private

    public bicgstab_solve

contains

    subroutine bicgstab_solve(A, b, x)

        complex(8), intent(in) :: A(:, :), b(:)
        complex(8), intent(inout) :: x(:)
        complex(8), dimension(size(x)) :: r, rhat, p, v, s, t
        complex(8) :: rho, rhoold, alpha, beta, omega

        r = b - matmul(A, x)
        rhat = r
        rho = 1
        alpha = 1
        omega = 1
        v = 0
        p = 0

        do while (real(sqrt(sum(r**2))) > 1d-3)
            rhoold = rho
            rho = dot_product(rhat, r)
            beta = (rho / rhoold) * (alpha / omega)
            p = r + beta * (p - omega * v)
            v = matmul(A, p)
            alpha = rho / dot_product(rhat, v)
            s = r - alpha * v
            t = matmul(A, s)
            omega = dot_product(t, s) / dot_product(t, t)
            x = x + alpha * p + omega * s
            r = s - omega * t 
        end do

    end subroutine

end module
