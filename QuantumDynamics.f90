program QuantumDynamics

    use plplot

    implicit none

    call plparseopts(PL_PARSE_FULL)

    print *, "test"
    call printing()

contains

    subroutine printing()

        print *, "and some more"

    end subroutine

end program
