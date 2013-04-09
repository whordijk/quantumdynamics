module wave1d

    use BiCGSTAB
    use plplot

    implicit none
    private

    complex(8), allocatable :: A(:, :), b(:), psi(:)
    real(8), allocatable :: x(:)
    real(8) :: dx, dt

    public init_model
    public step
    public plot_wave

contains

    subroutine init_model(n)

        integer, intent(in) :: n
        real(8) :: p1, p2

        dx = 1d0 / (n - 1)
        dt = 1d-5

        allocate(A(n, n), b(n), psi(n), x(n))

        call init_matrix()
        call linspace(n, x)
        p1 = 0.25
        p2 = p1 + 1 + dx
        psi = exp(-100 * (x - p1)**2 / 2) &
            * cmplx(cos(100 * (x - p1)), sin(100 * (x - p1))) &
            + exp(-100 * (x - p2)**2 / 2) &
                * cmplx(cos(100 * (x - p2)), sin(100 * (x - p2)))
        b = matmul(conjg(A), psi)

    end subroutine

    subroutine step()

        real(8) :: eps

        ! eps = dx**2 * dt**2
        eps = 1d-16
        call bicgstab_solve(A, b, psi, eps)
        b = matmul(conjg(A), psi)        

    end subroutine

    subroutine plot_wave()

        call plclear()
        call plline(x, real(psi)) 
        call plflush()

    end subroutine

    subroutine init_matrix()

        real(8), dimension(size(A, 1), size(A, 2)) :: Hamiltonian, eye, D, V
        integer :: i, n

        n = size(A, 1)
        D = 0
        eye = 0
        V = 0
        do i = 1, n
            D(i, i) = -2
            eye(i, i) = 1
            if (i == 1) then
                D(i + 1, i) = 1
            else if (i == n) then
                D(i - 1, i) = 1
            else
                D(i + 1, i) = 1
                D(i - 1, i) = 1
            end if
        end do
        D(1, n) = 1
        D(n, 1) = 1
        D = 1 / dx**2 * D
        Hamiltonian = -0.5 * D + V
        A = cmplx(eye, dt / 2 * Hamiltonian)

    end subroutine

    subroutine linspace(n, x)

        integer, intent(in) :: n
        real(8), intent(out) :: x(n)
        integer :: i

        do i = 1, n
            x(i) = (i - 1) * dx
        end do

    end subroutine

end module
