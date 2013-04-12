module wave1d

    use BiCGSTAB
    use plplot

    implicit none
    private

    real(8), parameter :: pi = 4 * atan(1d0)
    complex(8), allocatable :: A(:, :), b(:), psi(:)
    real(8), allocatable :: potential(:), x(:)
    real(8) :: L, dx, dt, eps
    integer :: n

    public init_model
    public step
    public plot_wave

contains

    subroutine init_model(sample_length)

        real(8), intent(in) :: sample_length

        L = sample_length
        n = nint(sample_length * 100)
        dx = sample_length / (n - 1)
        dt = 5d-4
        eps = dx**2 * dt**2

        allocate(A(n, n), b(n), psi(n), x(n), potential(n))

        call linspace(n, x)
        call init_matrix(.true.) ! argument switches the potential on/off
        call init_wave(0.5d0)

    end subroutine

    subroutine step()

        call bicgstab_solve(A, b, psi, eps)
        b = matmul(conjg(A), psi)        

    end subroutine

    subroutine plot_wave()

        call plclear()
        call plline(x, real(psi))
        call plcol0(6)
        call plline(x, aimag(psi))
        call plcol0(3)
        call plline(x, 2 * pi**2 * real(psi * conjg(psi)))
        call plcol0(7)
        call plline(x, potential)
        call plcol0(1)
        call plflush()

    end subroutine

    subroutine init_matrix(switch)

        logical, intent(in) :: switch
        real(8), dimension(n, n) :: Hamiltonian, eye, D, V
        integer :: i

        D = 0
        eye = 0
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
        
        potential = 0
        V = 0
        if (switch .eqv. .true.) then
            do i = 1, n
                V(i, i) = (2 / L * x(i) - 1)**12
                potential(i) = 0.05 * V(i, i)
                V(i, i) = 10000 * V(i,i)
            end do
        end if
        D(1, n) = 1
        D(n, 1) = 1
        D = 1 / dx**2 * D
        
        Hamiltonian = -0.5 * D + V
        A = cmplx(eye, dt / 2 * Hamiltonian)

    end subroutine

    subroutine init_wave(p)

        real(8), intent(in) :: p
        real(8), parameter :: k = 50
        real(8) :: p1, p2, d

        d = 50
        p1 = p * L
        p2 = p1 + L + dx
        psi = exp(-d * (x - p1)**2 / 2) &
            * cmplx(cos(k * (x - p1)), sin(k * (x - p1))) &
            / sqrt(2 * pi * d) &
            + exp(-d * (x - p2)**2 / 2) &
            * cmplx(cos(k * (x - p2)), sin(k * (x - p2))) &
            / sqrt(2 * pi * d)
        b = matmul(conjg(A), psi)

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
