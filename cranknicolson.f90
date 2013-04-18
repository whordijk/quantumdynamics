module cranknicolson

    use plplot
    use BiCGSTAB

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

        real(8), intent(in) :: sample_length

        L = sample_length
        n = nint(sample_length * 100)
        dx = sample_length / (n - 1)
        dt = dx**2
        eps = dx**2 * dt**3

        allocate(A(n, n), b(n), psi(n), x(n), potential(n))

        call linspace(n, x)
        call init_matrix(.false.) ! argument switches the potential on/off
        call init_wave(0.2d0)    ! fraction of the domain the wave is set

    end subroutine

    subroutine step()

        !complex(8) :: dl(n - 1), du(n - 1)
        !complex(8) :: d(n)
        !integer :: info
        !integer :: i

        !d(1) = A(1, 1)
        !do i = 2, n
        !    dl(i - 1) = A(i, i - 1)
        !    du(i - 1) = A(i - 1, i)
        !    d(i) = A(i, i)
        !end do
        !call zgtsv(n, 1, dl, d, du, b, n, INFO)
        call bicgstab_solve(A, b, psi, eps)
        !psi = b
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
                if (i > 0.375 * n .and. i < 0.425 * n) then
                    V(i, i) = 1
                else if (i > 0.575 * n .and. i < 0.625 * n) then
                    V(i, i) = 1
                end if
                !V(i, i) = (2 / L * x(i) - 1)**80
                potential(i) = 0.01 * V(i, i)
                V(i, i) = 1000 * V(i,i)
            end do
        end if
        ! periodic boundary condtions, need to be turned off for CGTSV method!
        D(1, n) = 1
        D(n, 1) = 1
        D = 1 / dx**2 * D
        
        Hamiltonian = -0.5 * D + V
        A = eye + ii * dt / 2 * Hamiltonian

    end subroutine

    subroutine init_wave(p)

        real(8), intent(in) :: p
        real(8), parameter :: k = 50
        real(8) :: arg(n)
        real(8) :: x_0, d

        d = 50
        x_0 = p * L
        arg = x - x_0
        psi = exp(-d * arg**2 / 2) &
            * exp(ii * k * arg) &
            / sqrt(2 * pi * d)
        arg = x - (x_0 + L + dx)
        psi = psi + exp(-d * arg**2 / 2) &
            * exp(ii * k * arg) &
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
