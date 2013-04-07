module Wave1d

    use BiCGSTAB

    implicit none
    private

    complex(8), allocatable :: A(:, :), b(:), psi(:)

    public init_model
    public calc_wave
    public plot_wave

contains

    subroutine init_model(n)

        integer, intent(in) :: n
        real(8) :: x(n)
        allocate(A(n, n), b(n), psi(n))

        call init_matrix(A)
        call linspace(n, x)
        psi = exp(-x**2) * cmplx(cos(x), sin(x))
        b = matmul(conjg(A), psi)

        print *, "MODEL INITIALIZED"

    end subroutine

    subroutine calc_wave()

        call bicgstab_solve(A, b, psi)
        b = matmul(conjg(A), psi)        
        print *, "calculated wave"

    end subroutine

    subroutine plot_wave()

        print *, "plot wave"

    end subroutine

    subroutine init_matrix(A)

        complex(8), intent(out) :: A(:, :)
        real(8), dimension(size(A, 1), size(A, 2)) :: Hamiltonian, eye, D, V
        real(8) :: h
        integer :: i, n

        n = size(A, 1)
        h = 1d0 / n
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
        A = cmplx(eye, h / 2 * Hamiltonian)

    end subroutine

    subroutine linspace(n, x)

        integer, intent(in) :: n
        real(8), intent(out) :: x(n)
        real(8) :: dx
        integer :: i

        dx = 1 / (n - 1)
        do i = 1, n
            x(i) = (i - 1) * dx
        end do

    end subroutine

end module
