module Wave1d

    use BiCGSTAB
    use plplot

    implicit none
    private

    complex(8), allocatable :: A(:, :), b(:), psi(:)
    real(8), allocatable :: x(:)

    public init_model
    public step
    public plot_wave

contains

    subroutine init_model(n)

        integer, intent(in) :: n

        allocate(A(n, n), b(n), psi(n), x(n))

        call init_matrix()
        call linspace(n, x)
        psi = exp(-200 * x**2) * cmplx(cos(200 * x), sin(200 * x))
        b = matmul(conjg(A), psi)

    end subroutine

    subroutine step()

        call bicgstab_solve(A, b, psi)
        b = matmul(conjg(A), psi)        

    end subroutine

    subroutine plot_wave()

        call plclear()
        call plline(x, real(psi)) 
        call plflush()

    end subroutine

    subroutine init_matrix()

        real(8), dimension(size(A, 1), size(A, 2)) :: Hamiltonian, eye, D, V
        real(8) :: h, dt
        integer :: i, n

        n = size(A, 1)
        h = 1d0 / n
        dt = 1d-3
        D = 0
        eye = 0
        V = 0
        do i = 1, n
            D(i, i) = 2
            eye(i, i) = 1
            if (i == 1) then
                D(i + 1, i) = -1
            else if (i == n) then
                D(i - 1, i) = -1
            else
                D(i + 1, i) = -1
                D(i - 1, i) = -1
            end if
        end do
        D = 1 / h**2 * D
        Hamiltonian = -0.5 * D + V
        A = cmplx(eye, dt / 2 * Hamiltonian)

    end subroutine

    subroutine linspace(n, x)

        integer, intent(in) :: n
        real(8), intent(out) :: x(n)
        real(8) :: dx
        integer :: i

        dx = 1d0 / (n - 1)
        do i = 1, n
            x(i) = (i - 1) * dx
        end do

    end subroutine

end module
