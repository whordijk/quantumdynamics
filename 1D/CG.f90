module CG

    implicit none
    private

    public cg_solve

contains

    subroutine cg_solve(A, b, x, eps)

        complex(8), intent(in) :: A(:, :), b(:)
        real(8), intent(in) :: eps
        complex(8), intent(inout) :: x(:)
        complex(8), dimension(size(x)) :: r, rold, p
        real(8) :: alpha, beta
        real(8) :: babs, error

        r = b - matmul(A, x)
        p = r
        babs = sqrt(sum(abs(b)**2))
        error = sqrt(sum(abs(r)**2)) / babs

        do while (error > eps)
            alpha = real(dot_product(conjg(r), r) / dot_product(conjg(p), matmul(A, p)))
            x = x + alpha * p
            rold = r
            r = r - alpha * matmul(A, p)
            beta = real(dot_product(conjg(r), r) / dot_product(conjg(rold), rold))
            p = r + beta * p
            error = sqrt(sum(abs(r)**2)) /  babs
        end do

    end subroutine

end module
