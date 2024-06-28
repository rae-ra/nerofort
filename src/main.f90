program Main
    implicit none

    character(len=1) :: pathsep
    character(len=10) :: input
    logical :: continue_program

    ! Determine the path separator based on the operating system
    pathsep = "?"
#ifdef WOS
    pathsep = "\"
#endif
#ifdef UOS
    pathsep = "/"
#elif __unix
    pathsep = "/"
#endif
#if __APPLE__ || __MACH__
    pathsep = "/"
#endif

    continue_program = .true.

    ! Main loop to handle user input and execute options
    do while (continue_program)
        print *, "Please enter an option: tests, dummy, exit"
        read *, input

        select case (input)
            case ('tests')
                call run_tests()
            case ('dummy')
                call run_dummy()
            case ('exit')
                continue_program = .false.
            case default
                print *, "Invalid option. Please try again."
        end select
    end do

    contains

    subroutine run_tests()
        use test_pooling2d
        implicit none

        print *, "Running all tests..."
        call test_averagepool_backprop()
    end subroutine run_tests

    subroutine run_dummy()
        implicit none

        character (len=1) :: pathsep
        pathsep = "?"
#ifdef WOS
        pathsep = "\"
#endif
#ifdef UOS
        pathsep = "/"
#elif __unix
        pathsep = "/"
#endif
#if __APPLE__ || __MACH__
        pathsep = "/"
#endif

        print *, "Hey"
        write (*, '( "pathsep is >", A1, "<")' )  pathsep
    end subroutine run_dummy

end program Main
